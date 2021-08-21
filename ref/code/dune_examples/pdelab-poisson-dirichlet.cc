// B.14. The Poisson Equation with Dirichlet Boundary
// Conditions, Using dune-pdelab
// File: pdelab-poisson-dirichlet.cc
// Equation: Poisson equation
// Discretization: First-order Lagrange finite elements
// Grid: Unstructured triangle grid, implemented with UGGrid
// Solver: CG solver with SSOR preconditioner
// Distributed: no
// Discussed in: Chapter 11.3.2

#include "config.h"

#include <memory>

#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/io/file/vtk.hh>
#include <dune/grid/uggrid.hh>

#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>
#include <dune/functions/functionspacebases/interpolate.hh>

#include <dune/pdelab.hh>

using namespace Dune;

// { parameterclass_begin }
template <class GridView, class Range>
class PoissonProblem
    : public PDELab::ConvectionDiffusionModelProblem<GridView, Range>
{
public:
    using Traits = typename PDELab::ConvectionDiffusionModelProblem<GridView, Range>::
        Traits;
    // Source term
    auto f(const typename Traits::ElementType &element,
           const typename Traits::DomainType &xi) const
    {
        return -5.0;
    }

    //! Boundary condition type function
    auto bctype(const typename Traits::IntersectionType &intersection,
                const typename Traits::IntersectionDomainType &xi) const
    {
        auto x = intersection.geometry().global(xi);
        return (x[0] < 1e-8 || x[1] < 1e-8 || (x[0] > 0.4999 && x[1] > 0.4999))
                   ? PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet
                   : PDELab::ConvectionDiffusionBoundaryConditions::Neumann;
    }
};
// { parameterclass_end }

// { main_begin }
int main(int argc, char *argv[])
{
    // Initialize MPI, if available
    MPIHelper::instance(argc, argv);

    constexpr int dim = 2;
    using Grid = UGGrid<dim>;
    std::shared_ptr<Grid> grid = GmshReader<Grid>::read("l-shape.msh");

    grid->globalRefine(2);

    using GridView = Grid::LeafGridView;
    GridView gridView = grid->leafGridView();
    // { grid_setup_end }

    // { make_function_space_begin }
    using Basis = Functions::LagrangeBasis<GridView, 1>;
    auto basis = std::make_shared<Basis>(gridView);

    using VectorBackend = PDELab::ISTL::VectorBackend<>;
    using Constraints = PDELab::ConformingDirichletConstraints;

    using GridFunctionSpace = PDELab::Experimental::GridFunctionSpace<Basis,
                                                                      VectorBackend,
                                                                      Constraints>;

    GridFunctionSpace gridFunctionSpace(basis);
    // { make_function_space_end }

    // Assemble constraints on this space
    // { make_bc_adapter_begin }
    using Problem = PoissonProblem<GridView, double>;
    Problem problem;
    PDELab::ConvectionDiffusionBoundaryConditionAdapter<Problem> bctype(problem);
    // { make_bc_adapter_end }

    // { make_constraints_begin }
    using ConstraintsContainer = GridFunctionSpace::ConstraintsContainer<double>::Type;
    ConstraintsContainer constraintsContainer;
    PDELab::constraints(bctype, gridFunctionSpace, constraintsContainer);
    // { make_constraints_end }

    // Make grid operator
    // { make_assembler_begin }
    using LocalOperator = PDELab::ConvectionDiffusionFEM<Problem,
                                                         GridFunctionSpace::Traits::
                                                             FiniteElementMap>;
    LocalOperator localOperator(problem);

    using MatrixBackend = PDELab::ISTL::BCRSMatrixBackend<>;
    MatrixBackend matrixBackend(7);

    using GridOperator = PDELab::GridOperator<GridFunctionSpace,
                                              GridFunctionSpace,
                                              LocalOperator,
                                              MatrixBackend,
                                              double, double, double,
                                              ConstraintsContainer,  // For trial space
                                              ConstraintsContainer>; // For test space

    GridOperator gridOperator(gridFunctionSpace,
                              constraintsContainer, // Trial space
                              gridFunctionSpace,
                              constraintsContainer, // Test space
                              localOperator,
                              matrixBackend);
    // { make_assembler_end }

    // Make solution vector
    // { make_solution_vector_begin }
    using U = PDELab::Backend::Vector<GridFunctionSpace, double>;
    U u(gridFunctionSpace, 0.0);
    auto g = [](auto x)
    { return (x[0] < 1e-8 || x[1] < 1e-8) ? 0 : 0.5; };
    Functions::interpolate(*basis, PDELab::Backend::native(u), g);
    // { make_solution_vector_end }

    // { linear_solver_begin }
    // Select a linear solver backend
    using LinearSolverBackend = PDELab::ISTLBackend_SEQ_CG_SSOR;
    LinearSolverBackend linearSolverBackend(50, 2);

    // Select linear problem solver
    using LinearProblemSolver = PDELab::StationaryLinearProblemSolver<GridOperator,
                                                                      LinearSolverBackend,
                                                                      U>;
    LinearProblemSolver linearProblemSolver(gridOperator,
                                            linearSolverBackend,
                                            u,
                                            1e-5);
    // Solve linear problem.
    linearProblemSolver.apply();
    // { linear_solver_end }

    // Graphical output
    // { vtk_output_begin }
    VTKWriter<GridView> vtkwriter(gridView);
    auto uFunction = Functions::makeDiscreteGlobalBasisFunction<double>(
        *basis,
        PDELab::Backend::native(u));
    vtkwriter.addVertexData(uFunction,
                            VTK::FieldInfo("u",
                                           VTK::FieldInfo::Type::scalar,
                                           1));
    vtkwriter.write("pdelab-poisson-dirichlet-result");
    // { vtk_output_end }
}