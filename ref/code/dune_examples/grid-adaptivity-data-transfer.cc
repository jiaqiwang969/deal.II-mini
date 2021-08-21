// File: grid-adaptivity-data-transfer.cc
// Equation: none
// Discretization: First-order Lagrange finite elements
// Grid: Unstructured simplicial grid, implemented with UGGrid
// Solver: none
// Distributed: no
// Discussed in: Chapter 5.9.2

#include "config.h"

#include <iostream>
#include <map>

#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dune/grid/uggrid.hh>
#include <dune/grid/utility/structuredgridfactory.hh>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/exceptions.hh>

using namespace Dune;

template <int dim>
class Sphere
{
    double radius_;
    Dune::FieldVector<double, dim> center_;

public:
    Sphere(const Dune::FieldVector<double, dim> &c, const double &r) : radius_(r), center_(c) {}

    double distanceTo(const Dune::FieldVector<double, dim> &point) const
    {
        return std::abs((center_ - point).two_norm() - radius_);
    }

    void displace(const Dune::FieldVector<double, dim> &increment)
    {
        center_ += increment;
    }
};

// Linear interpolation on a simplex
// { linear_interpolation_begin }
template <int dim>
double interpolate(const std::vector<double> values,
                   FieldVector<double, dim> p)
{
    assert(values.size() == dim + 1);
    double result = values[0];
    for (std::size_t i = 0; i < p.size(); i++)
        result += p[i] * (values[i + 1] - values[0]);
    return result;
}
// { linear_interpolation_end }

// { grid_setup_begin }
int main(int argc, char *argv[])
{
    // Set up MPI if available
    MPIHelper::instance(argc, argv);

    constexpr int dim = 2; // Grid dimension
    using Grid = UGGrid<dim>;

    // Create UGGrid from structured triangle grid
    const std::array<unsigned, dim> n = {8, 25};

    const FieldVector<double, dim> lower = {0, 0};
    const FieldVector<double, dim> upper = {6, 15};

    std::shared_ptr<Grid> grid = StructuredGridFactory<Grid>::createSimplexGrid(lower, upper, n);

    using GridView = Grid::LeafGridView;
    const GridView gridView = grid->leafGridView();
    const auto &indexSet = gridView.indexSet();
    const auto &idSet = grid->localIdSet();
    // { grid_setup_end }

    // Create sphere
    // { sphere_setup_begin }
    Sphere<dim> sphere({3.0, 2.5}, 1.0);

    // Set parameters
    const int steps = 30;
    const FieldVector<double, dim> stepDisplacement = {0, 0.5};

    const double epsilon = 0.4; // Thickness of the refined region
    // around the sphere
    const int levels = 2;
    // { sphere_setup_end }

    // Construct a piecewise linear function on the grid
    // { function_setup_begin }
    auto dataFunction = [](const FieldVector<double, dim> &x)
    {
        return std::sin(x[1]);
    };

    std::vector<double> data(gridView.size(dim));
    for (auto &&v : vertices(gridView))
        data[indexSet.index(v)] = dataFunction(v.geometry().corner(0));
    // { function_setup_end }

    // { coarsening_refinement_begin }
    for (int i = 0; i < steps; ++i)
    {
        std::map<Grid::LocalIdSet::IdType, double> persistentContainer;

        // Coarsen everything
        for (int j = 0; j < levels; ++j)
        {
            for (const auto &element : elements(gridView))
                grid->mark(-1, element); // Mark element for coarsening

            grid->preAdapt();

            for (const auto &vertex : vertices(gridView))
                persistentContainer[idSet.id(vertex)] = data[indexSet.index(vertex)];

            grid->adapt();

            data.resize(gridView.size(dim));

            for (const auto &v : vertices(gridView))
                data[indexSet.index(v)] = persistentContainer[idSet.id(v)];

            grid->postAdapt();
        }
        // { coarsening_end }

        // Refine near the sphere
        // { refinement_begin }
        for (int j = 0; j < levels; ++j)
        {
            // Select elements that are close to the sphere for grid refinement
            for (const auto &element : elements(gridView))
                if (sphere.distanceTo(element.geometry().center()) < epsilon)
                    grid->mark(1, element);

            grid->preAdapt();

            for (const auto &vertex : vertices(gridView))
                persistentContainer[idSet.id(vertex)] = data[indexSet.index(vertex)];

            grid->adapt();
            // { refinement_after_modification }

            // { refinement_data_interpolation_begin }
            data.resize(gridView.size(dim));

            for (const auto &element : elements(gridView))
            {
                if (element.isNew())
                {
                    for (std::size_t k = 0; k < element.subEntities(dim); k++)
                    {
                        auto father = element;
                        auto positionInFather = ReferenceElements<double, dim>::general(element.type())
                                                    .position(k, dim);

                        do
                        {
                            positionInFather = father.geometryInFather().global(positionInFather);
                            father = father.father();
                        } while (father.isNew());

                        // Extract corner values
                        std::vector<double> values(father.subEntities(dim));
                        for (std::size_t l = 0; l < father.subEntities(dim); l++)
                            values[l] = persistentContainer[idSet.subId(father, l, dim)];

                        // Interpolate linearly on the ancestor simplex
                        data[indexSet.subIndex(element, k, dim)] = interpolate(values, positionInFather);
                    }
                }
                else
                    for (std::size_t k = 0; k < element.subEntities(dim); k++)
                        data[indexSet.subIndex(element, k, dim)] = persistentContainer[idSet.subId(element, k, dim)];
            }

            grid->postAdapt();
            // { refinement_data_interpolation_end }
        }
        // { coarsening_refinement_end }

        // Write grid to file
        // { writing_moving_begin }
        VTKWriter<GridView> vtkWriter(gridView);
        vtkWriter.addVertexData(data, "data");
        vtkWriter.write("refined_grid_" + std::to_string(i));

        // Move sphere
        sphere.displace(stepDisplacement);
    }
    // { writing_moving_end }
}