// B.9. Dual Mortar Finite Elements
// File: localfunctions-use-dual-mortar-basis.cc
// Equation: none
// Discretization: Dual mortar finite elements
// Grid: Simplicial grid, implemented with UGGrid
// Solver: none
// Distributed: no
// Discussed in: Chapter 8.3.1

#include "config.h"

#include <iostream>

#include <dune/istl/io.hh>
#include <dune/istl/matrix.hh>

#include <dune/localfunctions/lagrange.hh>
#include <dune/localfunctions/dualmortarbasis.hh>

#include <dune/geometry/quadraturerules.hh>

#include <dune/grid/uggrid.hh>
#include <dune/grid/utility/structuredgridfactory.hh

using namespace Dune;

// { use_dual_mortar_basis_begin }
int main(int argc, char *argv[])
{
    // Set up MPI if available
    MPIHelper::instance(argc, argv);

    // Construct a 2x2 structured triangle grid
    constexpr int dim = 2;
    using Grid = UGGrid<dim>;

    auto grid = StructuredGridFactory<Grid>::createSimplexGrid({0.0, 0.0},
                                                               {1.0, 1.0},
                                                               {2, 2});
    auto gridView = grid->leafGridView();
    // { grid_setup_end }

    // Set up Lagrange and dual finite elements for simplex elements
    // { space_setup_begin }
    LagrangeSimplexLocalFiniteElement<double, double, dim, 1> lagrangeFE;
    DualP1LocalFiniteElement<double, double, dim, true> dualFE;

    // Matrix that will store the integrals
    Matrix<double> facetMassMatrix;
    // { space_setup_end }

    // { intersection_loop_begin }
    for (const auto &element : elements(gridView))
    {
        for (const auto &intersection : intersections(gridView, element))
        {
            if (intersection.boundary())
            {
                auto refElement = referenceElement<double, dim>(element.type());
                // { intersection_loop_end }

                // { zero_matrix_begin }
                facetMassMatrix.setSize(lagrangeFE.size(), dualFE.size());
                facetMassMatrix = 0;
                // { zero_matrix_end }

                // { quad_loop_begin }
                constexpr auto intersectionDim = std::decay_t<decltype(intersection)>::mydimension;
                const auto &quad = QuadratureRules<double, intersectionDim>::rule(intersection.type(),
                                                                                  2);

                for (const auto &quadPoint : quad)
                {
                    // { quad_loop_end }
                    auto quadPosInElement = intersection.geometryInInside().global(quadPoint.position());

                    std::vector<FieldVector<double, 1>> lagrangeValues, dualValues;
                    lagrangeFE.localBasis().evaluateFunction(quadPosInElement,
                                                             lagrangeValues);
                    dualFE.localBasis().evaluateFunction(quadPosInElement, dualValues);
                    // { evaluation_end }

                    // { integration_begin }
                    for (std::size_t i = 0; i < lagrangeFE.size(); i++)
                    {
                        LocalKey lagrangeLocalKey = lagrangeFE.localCoefficients().localKey(i);

                        auto lagrangeSubEntities = refElement.subEntities(intersection.indexInInside(),
                                                                          1, // The codimension of facets
                                                                          lagrangeLocalKey.codim());

                        if (lagrangeSubEntities.contains(lagrangeLocalKey.subEntity()))
                        {
                            for (std::size_t j = 0; j < dualFE.size(); j++)
                            {
                                LocalKey dualLocalKey = dualFE.localCoefficients().localKey(j);

                                auto dualSubEntities = refElement.subEntities(intersection.indexInInside(),
                                                                              1, // The codimension of facets
                                                                              dualLocalKey.codim());

                                if (dualSubEntities.contains(dualLocalKey.subEntity()))
                                {
                                    auto integrationElement = intersection.geometry()
                                                                  .integrationElement(quadPoint.position());
                                    facetMassMatrix[i][j] += quadPoint.weight() * integrationElement * lagrangeValues[i] * dualValues[j];
                                }
                            }
                        }
                    }
                } // End of quadrature loop
                // { quad_loop_final }

                // Print the facet mass matrix, observe how it is diagonal!
                // { print_matrix_begin }
                printmatrix(std::cout, facetMassMatrix, "facetMassMatrix", "--");

            } // End of ’if boundary’
        }     // End of intersection loop
    }         //End of element loop
    // { intersection_loop_final }
}