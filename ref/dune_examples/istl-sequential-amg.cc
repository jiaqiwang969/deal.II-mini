// B.8. The Sequential AMG Preconditioner
// File: istl-sequential-amg.cc
// Equation: Poisson equation
// Discretization: First-order Lagrange finite elements
// Grid: Unstructured simplicial grid, implemented with UGGrid
// Solver: CG solver with AMG preconditioner
// Distributed: no
// Discussed in: Chapter 7.5.3

#include "config.h"

#include <iostream>

#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/paamg/amg.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/matrixmarket.hh>

using namespace Dune;

// { problem_setup_begin }
int main(int argc, char *argv[])
{
    using Matrix = BCRSMatrix<double>;
    using Vector = BlockVector<double>;
    Matrix A;
    Vector b;
    loadMatrixMarket(A, "getting-started-poisson-fem-4-refinements-matrix.mtx");
    loadMatrixMarket(b, "getting-started-poisson-fem-4-refinements-rhs.mtx");
    // { problem_setup_end }

    // Construct vector to hold the iterates
    // { initial_iterate_begin }
    Vector x(b.size());
    x = b;
    // { initial_iterate_end }

    // Set up the smoother
    // { smoother_setup_begin }
    using Smoother = SeqSSOR<Matrix, Vector, Vector>;

    Amg::SmootherTraits<Smoother>::Arguments smootherArgs;
    smootherArgs.iterations = 3;
    smootherArgs.relaxationFactor = 1;
    // { smoother_setup_end }

    // Set up the coarsening criterion
    // { coarsening_setup_begin }
    using Norm = Amg::FirstDiagonal;
    using Criterion = Amg::CoarsenCriterion<Amg::UnSymmetricCriterion<Matrix, Norm>>;

    Criterion criterion(15, // Maximum number
                        // of multigrid levels
                        2000); // Create coarse levels until
    // problem size gets smaller
    // than this
    criterion.setDefaultValuesIsotropic(2); // Aggregate sizes and shapes
    criterion.setAlpha(.67);                // Connections above this value
    // are "strong"
    criterion.setBeta(1.0e-4); // Connections below this value
    // are treated as zero
    criterion.setGamma(1); // Number of
    // coarse-level iterations
    criterion.setDebugLevel(2); // Print some debugging output
    // { coarsening_setup_end }

    // { solver_setup_begin }
    using LinearOperator = MatrixAdapter<Matrix, Vector, Vector>;
    using AMG = Amg::AMG<LinearOperator, Vector, Smoother>;

    LinearOperator linearOperator(A);
    AMG amg(linearOperator, criterion, smootherArgs);

    CGSolver<Vector> amgCG(linearOperator,
                           amg,  // Preconditioner
                           1e-5, // Desired residual reduction
                           50,   // Maximum number of iterations
                           2);   // Verbosity
    // { solver_setup_end }

    // { solve_begin }
    InverseOperatorResult r;
    amgCG.apply(x, b, r);

    std::cout << "CG+AMG did " << r.iterations << " iterations." << std::endl;
    // { solve_end }
}