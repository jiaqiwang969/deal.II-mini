// B.3. Finite Volume Method for the Linear Transport
Equation
// File: getting-started-transport-fv.cc
// Equation: Linear transport equation
// Discretization: Cell-centered finite volumes, with explicit time stepping
// Grid: Uniform quadrilateral grid, implemented with YaspGrid
// Solver: none
// Distributed: no
// Discussed in: Chapter 3.4

#include "config.h"

#include <iostream>
#include <vector>

#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/vtk.hh>

// { using_namespace_dune_begin }
using namespace Dune;
// { using_namespace_dune_end }

// { evolve_signature_begin }
template <class GridView, class Mapper>
void evolve(const GridView &gridView,
            const Mapper &mapper,
            double dt, // Time step size
            std::vector<double> &c,
            const std::function<FieldVector<double, GridView::dimension>(FieldVector<double, GridView::dimension>)> v,
            const std::function<double(FieldVector<double, GridView::dimension>)> inflow)
// { evolve_signature_end }
{
    // { evolve_init_begin }
    // Grid dimension
    constexpr int dim = GridView::dimension;

    // Allocate a temporary vector for the update
    std::vector<double> update(c.size());
    std::fill(update.begin(), update.end(), 0.0);
    // { evolve_init_end }

    // Compute update vector
    // { element_loop_begin }
    for (const auto &element : elements(gridView))
    {
        // Element geometry
        auto geometry = element.geometry();

        // Element volume
        double elementVolume = geometry.volume();

        // Unique element number
        typename Mapper::Index i = mapper.index(element);
        // { element_loop_end }

        // Loop over all intersections
        with neighbors and boundary
            // { intersection_loop_begin }
            for (const auto &intersection : intersections(gridView, element))
        {
            // Geometry of the intersection
            auto intersectionGeometry = intersection.geometry();

            // Center of intersection in global coordinates
            FieldVector<double, dim>
                intersectionCenter = intersectionGeometry.center();

            // Velocity at intersection center vij
            FieldVector<double, dim> velocity = v(intersectionCenter);

            // Center of the intersection in local coordinates
            const auto &intersectionReferenceElement = ReferenceElements<double, dim - 1>::general(intersection.type());
            FieldVector<double, dim - 1> intersectionLocalCenter = intersectionReferenceElement.position(0, 0);

            // Normal vector scaled with intersection area: nij | ij |
            FieldVector<double, dim> integrationOuterNormal = intersection.integrationOuterNormal(intersectionLocalCenter);

            // Compute factor occuring in flux formula: hvij , nij i| ij |
            double intersectionFlow = velocity * integrationOuterNormal;
            // { intersection_loop_initend }

            // { intersection_loop_mainbegin }
            // Outflow contributions
            update[i] -= c[i] * std::max(0.0, intersectionFlow) / elementVolume;

            // Inflow contributions
            if (intersectionFlow <= 0)
            {
                // Handle interior intersection
                if (intersection.neighbor())
                {
                    // Access neighbor
                    auto j = mapper.index(intersection.outside());
                    update[i] -= c[j] * intersectionFlow / elementVolume;
                }

                // Handle boundary intersection
                if (intersection.boundary())
                    update[i] -= inflow(intersectionCenter) * intersectionFlow / elementVolume;
            }
            // { intersection_loopend }
        } // End loop over all intersections
    }     // End loop over the grid elements
    // { element_loop_end }

    // { evolve_laststeps }
    // Update the concentration vector
    for (std::size_t i = 0; i < c.size(); ++i)
        c[i] += dt * update[i];
}
// { evolve_end }

// { main_begin }
int main(int argc, char *argv[])
{
    // Set up MPI, if available
    MPIHelper::instance(argc, argv);
    // { main_signature_end }

    // { create_grid_begin }
    constexpr int dim = 2;
    using Grid = YaspGrid<dim>;
    Grid grid({1.0, 1.0}, // Upper right corner, the lower left one is (0,0)
              {80, 80});  // Number of elements per direction

    using GridView = Grid::LeafGridView;
    GridView gridView = grid.leafGridView();
    // { create_grid_end }

    // Assigns a unique number to each element
    // { create_concentration_begin }
    MultipleCodimMultipleGeomTypeMapper<GridView>
        mapper(gridView, mcmgElementLayout());

    // Allocate a vector for the concentration
    std::vector<double> c(mapper.size());
    // { create_concentration_end }

    // Initial concentration
    // { lambda_initial_concentration_begin }
    auto c0 = [](const FieldVector<double, dim> &x)
    {
        return (x.two_norm() > 0.125 && x.two_norm() < 0.5) ? 1.0 : 0.0;
    };
    // { lambda_initial_concentration_end }

    // { sample_initial_concentration_begin }
    // Iterate over grid elements and evaluate c0 at element centers
    for (const auto &element : elements(gridView))
    {
        // Get element geometry
        auto geometry = element.geometry();

        // Get global coordinate of element center
        auto global = geometry.center();

        // Sample initial concentration c0 at the element center
        c[mapper.index(element)] = c0(global);
    }
    // { sample_initial_concentration_end }

    // Construct VTK writer
    // { construct_vtk_writer_begin }
    auto vtkWriter = std::make_shared<Dune::VTKWriter<GridView>>(gridView);
    VTKSequenceWriter<GridView>
        vtkSequenceWriter(vtkWriter,
                          "getting-started-transport-fv-result"); // File name

    // Write the initial values
    vtkWriter->addCellData(c, "concentration");
    vtkSequenceWriter.write(0.0); // 0.0 is the current time
    // { construct_vtk_writer_end }

    // Now do the time steps
    // { time_loop_begin }
    double t = 0;            // Initial time
    const double tend = 0.6; // Final time
    const double dt = 0.006; // Time step size
    int k = 0;               // Time step counter

    // Inflow boundary values
    auto inflow = [](const FieldVector<double, dim> &x)
    {
        return 0.0;
    };

    // Velocity field
    auto v = [](const FieldVector<double, dim> &x)
    {
        return FieldVector<double, dim>(1.0);
    };

    while (t < tend)
    {
        // Apply finite volume scheme
        evolve(gridView, mapper, dt, c, v, inflow);

        // Augment time and time step counter
        t += dt;
        ++k;

        // Write data. We do not have to call addCellData again!
        vtkSequenceWriter.write(t);

        // Print iteration number, time, and time step size
        std::cout << "k=" << k << " t=" << t << std::endl;
    }
    // { time_loop_end }
}
// { main_end }