// B.15. Demonstrating the Linear Algebra Backends
// File: pdelab-backends.cc
// Equation: Linear reaction–diffusion equation
// Discretization: Second-order Lagrange finite elements
// Grid: Uniform quadrilateral grid, implemented with YaspGrid
// Solver: CG solver with SSOR preconditioner, CG solver with Jacobi
// preconditioner, damped Richardson iteration
// Distributed: no
// Discussed in: Chapter 11.4

#include "config.h"

#include <iostream>
#include <vector>

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>

#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>

#include <dune/pdelab.hh>

using namespace Dune;

template <class GridView, class Range>
class LinearReactionDiffusionProblem
    : public PDELab::ConvectionDiffusionModelProblem<GridView, Range>
{
public:
    template <typename Element, typename Coord>
    auto f(const Element &element, const Coord &x) const
    {
        auto globalpos = element.geometry().global(x);
        decltype(globalpos) midpoint(0.5);
        globalpos -= midpoint;
        return (globalpos.two_norm() < 0.25) ? -10.0 : 10.0;
    }

    template <typename Element, typename Coord>
    auto c(const Element &element, const Coord &x) const
    {
        return 10.0;
    }
};

// { solve_istl_begin }
template <class GridView>
void solveReactionDiffusionProblemISTL(const GridView &gridView,
                                       std::string filename)
{
    // Make grid function space
    using Basis = Functions::LagrangeBasis<GridView, 2>;
    auto basis = std::make_shared<Basis>(gridView);

    using VectorBackend = PDELab::ISTL::VectorBackend<>;

    using GridFunctionSpace = PDELab::Experimental::GridFunctionSpace<Basis,
                                                                      VectorBackend,
                                                                      PDELab::NoConstraints>;

    GridFunctionSpace gridFunctionSpace(basis);
    // { istl_gfs_end }

    // { istl_gridoperator_begin }
    // Make grid operator
    LinearReactionDiffusionProblem<GridView, double> problem;
    using LocalOperator = PDELab::ConvectionDiffusionFEM<decltype(problem),
                                                         typename GridFunctionSpace::Traits::
                                                             FiniteElementMap>;
    LocalOperator localOperator(problem);

    using MatrixBackend = PDELab::ISTL::BCRSMatrixBackend<>;
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

    MatrixBackend matrixBackend(25);
    GridOperator gridOperator(gridFunctionSpace,
                              gridFunctionSpace,
                              localOperator,
                              matrixBackend);
    // { istl_gridoperator_end }

    // { istl_solver_begin }
    // Select vector data type to hold the iterates
    using VectorContainer = PDELab::Backend::Vector<GridFunctionSpace, double>;
    VectorContainer u(gridFunctionSpace, 0.0); // Initial iterate

    // Select a linear solver backend
    using LinearSolverBackend = PDELab::ISTLBackend_SEQ_CG_SSOR;
    LinearSolverBackend linearSolverBackend(5000, // Maximal number
                                            // of iterations
                                            2); // Verbosity level

    // Solve linear problem
    using LinearProblemSolver = PDELab::StationaryLinearProblemSolver<GridOperator,
                                                                      LinearSolverBackend,
                                                                      VectorContainer>;
    LinearProblemSolver linearProblemSolver(gridOperator,
                                            linearSolverBackend,
                                            u,
                                            1e-10);
    linearProblemSolver.apply();
    // { istl_solver_end }

    // { istl_vtk_begin }
    // Output as VTK file
    SubsamplingVTKWriter<GridView> vtkwriter(gridView, refinementLevels(2));
    auto uFunction = Functions::makeDiscreteGlobalBasisFunction<double>(
        *basis,
        PDELab::Backend::native(u));
    vtkwriter.addVertexData(uFunction,
                            VTK::FieldInfo("u",
                                           VTK::FieldInfo::Type::scalar,
                                           1));
    vtkwriter.write(filename);
}
// { solve_istl_end }

#if HAVE_EIGEN
// { solve_eigen_begin }
template <class GridView>
void solveReactionDiffusionProblemEigen(const GridView &gridView,
                                        std::string filename)
{
    // Make grid function space
    using Basis = Functions::LagrangeBasis<GridView, 2>;
    auto basis = std::make_shared<Basis>(gridView);

    using VectorBackend = PDELab::Eigen::VectorBackend;

    using GridFunctionSpace = PDELab::Experimental::GridFunctionSpace<Basis,
                                                                      VectorBackend,
                                                                      PDELab::NoConstraints>;

    GridFunctionSpace gridFunctionSpace(basis);

    // Make grid operator
    LinearReactionDiffusionProblem<GridView, double> problem;
    using LocalOperator = PDELab::ConvectionDiffusionFEM<decltype(problem),
                                                         typename GridFunctionSpace::Traits::
                                                             FiniteElementMap>;
    LocalOperator localOperator(problem);

    using MatrixBackend = PDELab::Eigen::MatrixBackend<>;
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

    MatrixBackend matrixBackend(25);
    GridOperator gridOperator(gridFunctionSpace,
                              gridFunctionSpace,
                              localOperator,
                              matrixBackend);

    // Select vector data type to hold the iterate
    using VectorContainer = PDELab::Backend::Vector<GridFunctionSpace, double>;
    VectorContainer u(gridFunctionSpace, 0.0); // Initial iterate

    // Select a linear solver backend
    using LinearSolverBackend = PDELab::EigenBackend_CG_Diagonal_Up;
    LinearSolverBackend linearSolverBackend(5000);

    // Solve linear problem
    using LinearProblemSolver = PDELab::StationaryLinearProblemSolver<GridOperator,
                                                                      LinearSolverBackend,
                                                                      VectorContainer>;
    LinearProblemSolver linearProblemSolver(gridOperator,
                                            linearSolverBackend,
                                            u,
                                            1e-10);
    linearProblemSolver.apply();

    // Output as VTK file
    SubsamplingVTKWriter<GridView> vtkwriter(gridView, refinementLevels(2));
    auto uFunction = Functions::makeDiscreteGlobalBasisFunction<double>(
        *basis,
        PDELab::Backend::native(u));
    vtkwriter.addVertexData(uFunction,
                            VTK::FieldInfo("u",
                                           VTK::FieldInfo::Type::scalar,
                                           1));
    vtkwriter.write(filename);
}
// { solve_eigen_end }
#endif

template <class GridView>
void solveReactionDiffusionProblemISTLExplicit(const GridView &gridView,
                                               std::string filename)
{
    // Make grid function space
    using Basis = Functions::LagrangeBasis<GridView, 2>;
    auto basis = std::make_shared<Basis>(gridView);

    using VectorBackend = PDELab::ISTL::VectorBackend<>;

    using GridFunctionSpace = PDELab::Experimental::GridFunctionSpace<Basis,
                                                                      VectorBackend,
                                                                      PDELab::NoConstraints>;

    GridFunctionSpace gridFunctionSpace(basis);

    // Make grid operator
    LinearReactionDiffusionProblem<GridView, double> problem;
    using LocalOperator = PDELab::ConvectionDiffusionFEM<decltype(problem),
                                                         typename GridFunctionSpace::Traits::FiniteElementMap>;
    LocalOperator localOperator(problem);

    using MatrixBackend = PDELab::ISTL::BCRSMatrixBackend<>;
    MatrixBackend matrixBackend(25);

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

    // Vector holding the solver iterates
    // { istl_explicit_make_iterate_begin }
    using VectorContainer = typename GridOperator::Traits::Domain;
    VectorContainer xContainer(gridFunctionSpace, 0.0);
    // { istl_explicit_make_iterate_end }

    // { construct_istl_explicit_begin }
    // Evaluate residual at the zero configuration
    typename GridOperator::Traits::Range rContainer(gridFunctionSpace, 0.0);

    gridOperator.residual(xContainer, rContainer);

    // Compute stiffness matrix
    using MatrixContainer = typename GridOperator::Jacobian;
    MatrixContainer matrixContainer(gridOperator, 0.0);

    gridOperator.jacobian(xContainer, matrixContainer);
    // { construct_istl_explicit_end }

    // { get_native_data_begin }
    auto &x = PDELab::Backend::native(xContainer);
    auto &b = PDELab::Backend::native(rContainer);
    b *= -1.0;
    const auto &stiffnessMatrix = PDELab::Backend::native(matrixContainer);
    // { get_native_data_end }

    // { get_native_types_begin }
    using DomainVector = std::decay_t<decltype(x)>; // Decay from reference
    // to value type
    using RangeVector = std::decay_t<decltype(b)>;
    using Matrix = std::decay_t<decltype(stiffnessMatrix)>;
    // { get_native_types_end }

    // { solve_istl_explicit_begin }
    // Turn the matrix into a linear operator
    MatrixAdapter<Matrix, DomainVector, RangeVector>
        linearOperator(stiffnessMatrix);

    // SSOR as the preconditioner
    SeqSSOR<Matrix, DomainVector, RangeVector> preconditioner(stiffnessMatrix,
                                                              1, // Number of
                                                              // iterations
                                                              1.0); // Dampi

    // Preconditioned conjugate-gradient solver
    CGSolver<DomainVector> cg(linearOperator,
                              preconditioner,
                              1e-4, // Desired residual reduction factor
                              50,   // Maximum number of iterations
                              2);   // Verbosity level

    // Object storing some statistics about the solving process
    InverseOperatorResult statistics;

    // Solve!
    cg.apply(x, b, statistics);
    // { solve_istl_explicit_end }

    // Output as VTK file
    SubsamplingVTKWriter<GridView> vtkwriter(gridView, Dune::refinementLevels(2));
    auto uFunction = Functions::makeDiscreteGlobalBasisFunction<double>(*basis, x);
    vtkwriter.addVertexData(uFunction, VTK::FieldInfo("u", VTK::FieldInfo::Type::scalar, 1));
    vtkwriter.write(filename);
}

// { solve_simple_begin }
template <class GridView>
void solveReactionDiffusionProblemSimple(const GridView &gridView,
                                         std::string filename)
{
    // Make grid function space
    using Basis = Functions::LagrangeBasis<GridView, 2>;
    auto basis = std::make_shared<Basis>(gridView);

    using VectorBackend = PDELab::Simple::VectorBackend<>;

    using GridFunctionSpace = PDELab::Experimental::GridFunctionSpace<Basis,
                                                                      VectorBackend,
                                                                      PDELab::NoConstraints>;

    GridFunctionSpace gridFunctionSpace(basis);
    // { simple_gfs_end }

    // Make grid operator
    // { simple_gridoperator_begin }
    LinearReactionDiffusionProblem<typename Basis::GridView, double> problem;
    using LocalOperator = PDELab::ConvectionDiffusionFEM<decltype(problem),
                                                         typename GridFunctionSpace::Traits::
                                                             FiniteElementMap>;
    LocalOperator localOperator(problem);

    using MatrixBackend = PDELab::Simple::SparseMatrixBackend<>;
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
                              localOperator);
    // { simple_gridoperator_end }

    // { simple_algebraic_problem_begin }
    // Vector for the iterates
    typename GridOperator::Traits::Domain xContainer(gridFunctionSpace, 0.0);

    // Evaluate residual at the zero vector
    typename GridOperator::Traits::Range rContainer(gridFunctionSpace, 0.0);

    gridOperator.residual(xContainer, rContainer);
    rContainer *= -1.0;

    // Compute stiffness matrix
    using MatrixContainer = typename GridOperator::Traits::Jacobian;
    MatrixContainer matrixContainer(gridOperator, 0.0);

    gridOperator.jacobian(xContainer, matrixContainer);
    // { simple_algebraic_problem_end }

    // { simple_solution_begin }
    double omega = 0.2;
    for (int i = 0; i < 200; i++)
    {
        // Damped Richardson iteration
        auto correction = rContainer;
        matrixContainer.usmv(-1.0, xContainer, correction);
        xContainer.axpy(omega, correction);
    }
    // { simple_solution_end }

    // Output as VTK file
    SubsamplingVTKWriter<GridView> vtkwriter(gridView, Dune::refinementLevels(2));
    auto uFunction = Functions::makeDiscreteGlobalBasisFunction<double>(*basis, PDELab::Backend::native(xContainer));
    vtkwriter.addVertexData(uFunction, VTK::FieldInfo("u", VTK::FieldInfo::Type::scalar, 1));
    vtkwriter.write(filename);
}

int main(int argc, char *argv[])
{
    // Set up MPI, if available
    MPIHelper::instance(argc, argv);

    // Create a structured 8x8 grid
    constexpr int dim = 2;
    FieldVector<double, dim> l = {1.0, 1.0};
    std::array<int, dim> n = {4, 4};
    YaspGrid<2> grid(l, n);
    auto gridView = grid.leafGridView();

    // ISTL backend
    solveReactionDiffusionProblemISTL(gridView,
                                      "pdelab-backends-result-istl");

#if HAVE_EIGEN
    // Eigen backend
    solveReactionDiffusionProblemEigen(gridView,
                                       "pdelab-backends-result-eigen");
#else
#warning Skipping the Eigen example, because the Eigen library was not found.
    std::cerr << "Skipping the Eigen example, because the Eigen library was not found." << std::endl;
#endif

    // ISTL backend, but giving the solver explicitly
    solveReactionDiffusionProblemISTLExplicit(gridView,
                                              "pdelab-backends-result-istl-explicit");

    // ’simple’ backend
    solveReactionDiffusionProblemSimple(gridView,
                                        "pdelab-backends-result-simple");
}