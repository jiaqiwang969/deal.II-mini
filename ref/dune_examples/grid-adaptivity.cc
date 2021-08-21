// File: grid-adaptivity.cc
// Equation: none
// Discretization: none
// Grid: Unstructured simplicial grid, implemented with UGGrid
// Solver: none
// Distributed: no
// Discussed in: Chapter 5.9.1

#include "config.h"

#include <iostream>

#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dune/grid/uggrid.hh>
#include <dune/grid/utility/structuredgridfactory.hh>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/exceptions.hh>

using namespace Dune;

// { sphere_begin }
template <int dim>
class Sphere
{
    double radius_;
    FieldVector<double, dim> center_;

public:
    Sphere(const FieldVector<double, dim> &center, const double &radius)
        : radius_(radius),
          center_(center)
    {
    }

    double distanceTo(const FieldVector<double, dim> &point) const
    {
        return std::abs((center_ - point).two_norm() - radius_);
    }

    void displace(const FieldVector<double, dim> &increment)
    {
        center_ += increment;
    }
};
// { sphere_end }

// { grid_setup_begin }
int main(int argc, char *argv[])
{
    // Set up MPI if available
    MPIHelper::instance(argc, argv);

    constexpr int dim = 2; // Grid and world dimension
    using Grid = UGGrid<dim>;

    // Start with a structured grid
    const std::array<unsigned, dim> n = {8, 25};

    const FieldVector<double, dim> lower = {0, 0};
    const FieldVector<double, dim> upper = {6, 15};

    std::shared_ptr<Grid> grid = StructuredGridFactory<Grid>::createSimplexGrid(lower, upper, n);

    using GridView = Grid::LeafGridView;
    const GridView gridView = grid->leafGridView();
    // { grid_setup_end }

    // Create sphere
    // { sphere_setup_begin }
    Sphere<dim> sphere({3.0, 2.5}, 1.0);

    // Set parameters
    const int steps = 30; // Total number of steps
    const FieldVector<double, dim>
        stepDisplacement = {0, 0.5}; // Sphere displacement per step

    const double epsilon = 0.4; // Thickness of the refined region
    // around the sphere
    const int levels = 3; // Number of refinement levels
    // { sphere_setup_end }

    // { coarsening_refinement_begin }
    for (int i = 0; i < steps; ++i)
    {
        std::cout << "Step " << i << std::endl;

        // Coarsen everything
        for (int k = 0; k < levels - 1; ++k)
        {
            for (const auto &element : elements(gridView))
                grid->mark(-1, element);

            grid->preAdapt();
            grid->adapt();
            grid->postAdapt();
        }

        // Refine near the sphere
        for (int k = 0; k < levels - 1; ++k)
        {
            // Select elements that are close to the sphere for grid refinement
            for (const auto &element : elements(gridView))
                if (sphere.distanceTo(element.geometry().center()) < epsilon)
                    grid->mark(1, element);

            grid->preAdapt();
            grid->adapt();
            grid->postAdapt();
        }
        // { coarsening_refinement_end }

        // Write grid to file
        // { writing_moving_begin }
        VTKWriter<GridView> vtkWriter(gridView);
        vtkWriter.write("refined_grid_" + std::to_string(i));

        // Move sphere
        sphere.displace(stepDisplacement);
    }
    // { writing_moving_end }
}