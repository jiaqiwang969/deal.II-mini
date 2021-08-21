// B.11. The Linear Reaction–Diffusion Equation Using
// dune-pdelab
// File: pdelab-linear-reaction-diffusion.cc
// Equation: Linear reaction–diffusion equation
// Discretization: First and second-order Lagrange finite elements, Rannacher–
// Turek finite elements
// Grid: Uniform quadrilateral grid, implemented with YaspGrid
// Solver: CG solver with SSOR preconditioner
// Distributed: no
// Discussed in: Chapter 11.1

#include "config.h"

#include <iostream>
#include <vector>

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>

#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/rannacherturekbasis.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>

#include <dune/pdelab.hh>

using namespace Dune;

// { problem_description_begin }
template <class GridView>
struct ReactionDiffusionProblem
    : public PDELab::ConvectionDiffusionModelProblem<GridView, double>
{
    template <typename Element, typename Coord>
    auto f(const Element &element, const Coord &xi) const
    {
        const Coord center = Coord(0.5);
        auto distanceToCenter = (element.geometry().global(xi) - center).two_norm();
        return (distanceToCenter <= 0.25) ? -10.0 : 10.0;
    }

    template <typename Element, typename Coord>
    auto c(const Element &element, const Coord &xi) const
    {
        return 10.0;
    }
};
// { problem_description_end }

// { solvereactiondiffusionproblem_begin }
template <class Basis>
void solveReactionDiffusionProblem(std::shared_ptr<Basis> basis,
                                   std::string filename,
                                   int averageNumberOfEntriesPerMatrixRow)
{
    // { solvereactiondiffusionproblem_signature_end }
    // Make grid function space
    // { construct_function_space_begin }
    using GridView = typename Basis::GridView;
    using VectorBackend = PDELab::ISTL::VectorBackend<>;

    using GridFunctionSpace = PDELab::Experimental::GridFunctionSpace<Basis,
                                                                      VectorBackend,
                                                                      PDELab::NoConstraints>;

    GridFunctionSpace gridFunctionSpace(basis);
    // { construct_function_space_end }

    // Make grid operator
    // { construct_element_assembler_begin }
    using Problem = ReactionDiffusionProblem<GridView>;
    Problem reactionDiffusionProblem;
    using LocalOperator = PDELab::ConvectionDiffusionFEM<Problem,
                                                         typename GridFunctionSpace::Traits::
                                                             FiniteElementMap>;
    LocalOperator localOperator(reactionDiffusionProblem);
    // { construct_element_assembler_end }

    // { construct_global_assembler_begin }
    using MatrixBackend = PDELab::ISTL::BCRSMatrixBackend<>;
    MatrixBackend matrixBackend(averageNumberOfEntriesPerMatrixRow);
    using GridOperator = PDELab::GridOperator<GridFunctionSpace, // Trial function space
                                              GridFunctionSpace, // Test function space
                                              LocalOperator,     // Element assembler
                                              MatrixBackend,     // Data structure
                                              // for the stiffness matrix
                                              double, // Number type for
                                              // solution vector entries
                                              double, // Number type for
                                              // residual vector entries
                                              double>; // Number type for
    // stiffness matrix entries

    GridOperator gridOperator(gridFunctionSpace,
                              gridFunctionSpace,
                              localOperator,
                              matrixBackend);
    // { construct_global_assembler_end }

    // Select vector data type to hold the iterate
    // { create_initial_iterate_begin }
    using VectorContainer = PDELab::Backend::Vector<GridFunctionSpace, double>;
    VectorContainer u(gridFunctionSpace, 0.0); // Initial value
    // { create_initial_iterate_end }

    // Select a linear solver backend, and solve the linear problem
    // { solving_begin }
    using LinearSolverBackend = PDELab::ISTLBackend_SEQ_CG_SSOR;
    LinearSolverBackend linearSolverBackend(5000, 2);

    using LinearProblemSolver = PDELab::StationaryLinearProblemSolver<GridOperator,
                                                                      LinearSolverBackend,
                                                                      VectorContainer>;

    LinearProblemSolver linearProblemSolver(gridOperator,
                                            linearSolverBackend,
                                            u,
                                            1e-10);
    linearProblemSolver.apply();
    // { solving_end }

    // Output as VTK file
    // { vtk_writing_begin }
    // Make a discrete function from the FE basis and the coefficient vector
    auto uFunction = Functions::makeDiscreteGlobalBasisFunction<double>(
        *basis,
        PDELab::Backend::native(u));

    // We need to subsample, because the dune-grid VTKWriter
    // cannot natively display second-order functions.
    SubsamplingVTKWriter<GridView> vtkWriter(basis->gridView(),
                                             refinementLevels(4));
    vtkWriter.addVertexData(uFunction,
                            VTK::FieldInfo("u",
                                           VTK::FieldInfo::Type::scalar,
                                           1));
    vtkWriter.write(filename);
    // { vtk_writing_end }
}

// { main_begin }
int main(int argc, char *argv[])
{
    // Set up MPI, if available
    MPIHelper::instance(argc, argv);

    // Create a structured grid
    constexpr int dim = 2;
    using Grid = YaspGrid<dim>;
    FieldVector<double, dim> bbox = {1.0, 1.0};
    std::array<int, dim> elements = {8, 8};
    Grid grid(bbox, elements);

    using GridView = Grid::LeafGridView;
    GridView gridView = grid.leafGridView();
    // { grid_setup_end }

    // { simulation_begin }
    // Solve boundary value problem using first-order Lagrange elements
    using Q1Basis = Functions::LagrangeBasis<GridView, 1>;
    auto q1Basis = std::make_shared<Q1Basis>(gridView);

    solveReactionDiffusionProblem(q1Basis,
                                  "pdelab-linear-reaction-diffusion-result-q1",
                                  9); // Average number of matrix entries
    // per row

    // Solve boundary value problem using second-order Lagrange elements
    using Q2Basis = Functions::LagrangeBasis<GridView, 2>;
    auto q2Basis = std::make_shared<Q2Basis>(gridView);

    solveReactionDiffusionProblem(q2Basis,
                                  "pdelab-linear-reaction-diffusion-result-q2",
                                  25); // Average number of matrix entries
    // per row

    // Solve boundary value problem using Rannacher-Turek elements
    using RannacherTurekBasis = Functions::RannacherTurekBasis<GridView>;
    auto rannacherTurekBasis = std::make_shared<RannacherTurekBasis>(gridView);

    solveReactionDiffusionProblem(rannacherTurekBasis,
                                  "pdelab-linear-reaction-diffusion-result-rt",
                                  7); // Average number of matrix entries
    // per row
}
// { main_end }