// B.16. Error Estimation and Adaptive Grid Refinement
// Using dune-pdelab
// File: pdelab-adaptivity.cc
// Equation: Poisson equation
// Discretization: First-order Lagrange finite elements
// Grid: Unstructured triangle grid, implemented with UGGrid
// Solver: CG solver with SSOR preconditioner
// Distributed: no
// Discussed in: Chapter 11.5.2

#include "config.h"

#include <iostream>
#include <vector>

#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/grid/uggrid.hh>

#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>

#include <dune/pdelab.hh>

using namespace Dune;

template <class GridView, class Range>
class PoissonProblem
    : public PDELab::ConvectionDiffusionModelProblem<GridView, Range>
{
public:
    using Traits = typename PDELab::ConvectionDiffusionModelProblem<GridView, Range>::Traits;

    // Source term
    auto f(const typename Traits::ElementType &element,
           const typename Traits::DomainType &xi) const
    {
        return -5.0;
    }

    //! boundary condition type function
    auto bctype(const typename Traits::IntersectionType &intersection,
                const typename Traits::IntersectionDomainType &xi) const
    {
        auto x = intersection.geometry().global(xi);
        return (x[0] < 1e-8 || x[1] < 1e-8 || (x[0] > 0.4999 && x[1] > 0.4999))
                   ? PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet
                   : PDELab::ConvectionDiffusionBoundaryConditions::Neumann;
    }
};

// { local_estimator_begin }
template <class Problem>
class ResidualErrorEstimator
    : public PDELab::LocalOperatorDefaultFlags
{
    Problem problem_;

public:
    ResidualErrorEstimator(const Problem &problem)
        : problem_(problem)
    {
    }
    // { local_estimator_setup_end }

    // Residual assembly flags
    // { local_estimator_flags_begin }
    static const bool doAlphaVolume = true;
    static const bool doAlphaSkeleton = true;
    static const bool doAlphaBoundary = true;
    // { local_estimator_flags_end }

    // { alpha_volume_begin }
    template <class EntityGeometry,
              class LocalFunctionSpaceU, class Vector,
              class LocalFunctionSpaceV, class Residual>
    void alpha_volume(const EntityGeometry &entityGeometry,
                      const LocalFunctionSpaceU &localFunctionSpaceU,
                      const Vector &x,
                      const LocalFunctionSpaceV &localFunctionSpaceV,
                      Residual &residual) const
    {
        using RangeField = typename LocalFunctionSpaceU::Traits::
            FiniteElement::Traits::LocalBasisType::Traits::RangeType;

        // Element diameter
        auto h_T = diameter(entityGeometry.geometry());

        // Loop over quadrature points
        int intOrder = localFunctionSpaceU.finiteElement().localBasis().order();
        constexpr auto entityDim = EntityGeometry::Entity::mydimension;
        const auto &quadRule = QuadratureRules<double, entityDim>::rule(entityGeometry.geometry()
                                                                            .type(),
                                                                        intOrder)

            for (const auto &quadPoint : quadRule)
        {
            // Laplace of uh is always zero,
            // because we are using first-order elements on simplices
            RangeField laplaceU = 0.0;

            // Evaluate source term function f
            auto f = problem_.f(entityGeometry.entity(), quadPoint.position());

            // Integrate h2 T(f + uh)2
            auto factor = quadPoint.weight() * entityGeometry.geometry()
                                                   .integrationElement(quadPoint.position());
            auto value = std::pow(h_T * (f + laplaceU), 2);
            residual.accumulate(localFunctionSpaceV, 0, value * factor);
        }
    }
    // { alpha_volume_end }

    // { alpha_skeleton_begin }
    template <class IntersectionGeometry,
              class LocalFunctionSpaceU, class Vector,
              class LocalFunctionSpaceV, class Residual>
    void alpha_skeleton(const IntersectionGeometry &intersectionGeometry,
                        const LocalFunctionSpaceU &localFunctionSpaceUIn,
                        const Vector &xIn,
                        const LocalFunctionSpaceV &localFunctionSpaceVIn,
                        const LocalFunctionSpaceU &localFunctionSpaceUOut,
                        const Vector &xOut,
                        const LocalFunctionSpaceV &localFunctionSpaceVOut,
                        Residual &residualIn, Residual &residualOut) const
    {
        // Extract type of shape function gradients
        using TrialFE = typename LocalFunctionSpaceU::Traits::FiniteElementType;
        using LocalBasisU = typename TrialFE::Traits::LocalBasisType;
        using GradientU = typename LocalBasisU::Traits::JacobianType;
        using size_type = typename LocalFunctionSpaceU::Traits::SizeType;

        auto insideGeometry = intersectionGeometry.inside().geometry();
        auto outsideGeometry = intersectionGeometry.outside().geometry();

        auto h_F = diameter(intersectionGeometry.geometry());

        auto geometryInInside = intersectionGeometry.geometryInInside();
        auto geometryInOutside = intersectionGeometry.geometryInOutside();

        std::vector<GradientU> localGradPhiIn(localFunctionSpaceUIn.size());
        std::vector<GradientU> localGradPhiOut(localFunctionSpaceUOut.size());

        // Loop over quadrature points and integrate the jump term
        const int intOrder = 2 * localFunctionSpaceUIn.finiteElement().localBasis().order();
        constexpr auto intersectionDim = IntersectionGeometry::mydimension;
        const auto &quadRule = QuadratureRules<double, intersectionDim>::rule(
            intersectionGeometry.geometry().type(),
            intOrder);

        for (const auto &quadPoint : quadRule)
        {
            // Position of quadrature point in local coordinates of elements
            auto quadPointLocalIn = geometryInInside.global(quadPoint.position());
            auto quadPointLocalOut = geometryInOutside.global(quadPoint.position());

            // Evaluate gradient of basis functions
            localFunctionSpaceUIn.finiteElement().localBasis().evaluateJacobian(quadPointLocalIn, localGradPhiIn);
            localFunctionSpaceUOut.finiteElement().localBasis().evaluateJacobian(quadPointLocalOut, localGradPhiOut);

            // Compute gradient of u
            GradientU localGradUIn(0.0);
            for (size_type i = 0; i < localFunctionSpaceUIn.size(); i++)
                localGradUIn.axpy(xIn(localFunctionSpaceUIn, i), localGradPhiIn[i]);
            GradientU localGradUOut(0.0);
            for (size_type i = 0; i < localFunctionSpaceUOut.size(); i++)
                localGradUOut.axpy(xOut(localFunctionSpaceUOut, i), localGradPhiOut[i]);

            GradientU gradUIn;
            auto jacIn = insideGeometry.jacobianInverseTransposed(quadPointLocalIn);
            jacIn.mv(localGradUIn[0], gradUIn[0]);

            GradientU gradUOut;

            GradientU gradUOut;
            auto jacOut = outsideGeometry.jacobianInverseTransposed(quadPointLocalOut);
            jacOut.mv(localGradUOut[0], gradUOut[0]);

            // Integrate
            const auto n_F = intersectionGeometry.unitOuterNormal(quadPoint.position());
            auto jumpSquared = std::pow((n_F * gradUIn[0]) - (n_F * gradUOut[0]), 2);
            auto factor = quadPoint.weight() * intersectionGeometry.geometry()
                                                   .integrationElement(quadPoint.position());

            // Accumulate indicator
            residualIn.accumulate(localFunctionSpaceVIn,
                                  0,
                                  0.5 * h_F * jumpSquared * factor);
            residualOut.accumulate(localFunctionSpaceVOut,
                                   0,
                                   0.5 * h_F * jumpSquared * factor);
        }
    }
    // { alpha_skeleton_end }

    // { alpha_boundary_begin }
    template <class IntersectionGeometry,
              class LocalFunctionSpaceU, class Vector,
              class LocalFunctionSpaceV, class Residual>
    void alpha_boundary(const IntersectionGeometry &intersectionGeometry,
                        const LocalFunctionSpaceU &localFunctionSpaceUIn,
                        const Vector &xIn,
                        const LocalFunctionSpaceV &localFunctionSpaceVIn,
                        Residual &residualIn) const
    {
        using TrialFE = typename LocalFunctionSpaceU::Traits::FiniteElementType;
        using LocalBasisU = typename TrialFE::Traits::LocalBasisType;
        using GradientU = typename LocalBasisU::Traits::JacobianType;
        using size_type = typename LocalFunctionSpaceU::Traits::SizeType;

        auto insideGeometry = intersectionGeometry.inside().geometry();

        auto geometryInInside = intersectionGeometry.geometryInInside();

        // Intersection diameter
        auto h_F = diameter(intersectionGeometry.geometry());

        std::vector<GradientU> localGradPhi(localFunctionSpaceUIn.size());

        // Loop over quadrature points and integrate the jump term
        const int intOrder = localFunctionSpaceUIn.finiteElement().localBasis().order();
        constexpr auto intersectionDim = IntersectionGeometry::mydimension;
        const auto &quadRule = QuadratureRules<double, intersectionDim>::rule(
            intersectionGeometry.geometry().type(),
            intOrder);

        for (const auto &quadPoint : quadRule)
        {
            // Skip if quadrature point is not on Neumann boundary
            auto bcType = problem_.bctype(intersectionGeometry.intersection(),
                                          quadPoint.position());
            if (bcType != PDELab::ConvectionDiffusionBoundaryConditions::Neumann)
                return;

            // Position of quadrature point in local coordinates of element
            auto quadPointLocalIn = geometryInInside.global(quadPoint.position());

            // Evaluate gradient of trial shape functions
            localFunctionSpaceUIn.finiteElement().localBasis().evaluateJacobian(quadPointLocalIn,
                                                                                localGradPhi);

            // Evaluate gradient of u_h
            GradientU localGradU(0.0);
            for (size_type i = 0; i < localFunctionSpaceUIn.size(); i++)
                localGradU.axpy(xIn(localFunctionSpaceUIn, i), localGradPhi[i]);

            GradientU gradU;
            auto jac = insideGeometry.jacobianInverseTransposed(quadPointLocalIn);
            jac.mv(localGradU[0], gradU[0]);

            // Evaluate Neumann boundary value function
            auto neumann = problem_.j(intersectionGeometry.intersection(),
                                      quadPoint.position());

            // Integrate
            auto factor = quadPoint.weight() * intersectionGeometry.geometry()
                                                   .integrationElement(quadPoint.position());

            const auto n = intersectionGeometry.unitOuterNormal(quadPoint.position());

            residualIn.accumulate(localFunctionSpaceVIn,
                                  0,
                                  h_F * std::pow((neumann - n * gradU[0]), 2) * factor);
        }
    }
    // { alpha_boundary_end }

private:
    template <class Geometry>
    static auto diameter(const Geometry &geometry)
    {
        typename Geometry::ctype h = 0.0;
        for (int i = 0; i < geometry.corners(); i++)
            for (int j = i + 1; j < geometry.corners(); j++)
                h = std::max(h, (geometry.corner(j) - geometry.corner(i)).two_norm());

        return h;
    }
};

// { main_begin }
int main(int argc, char *argv[])
{
    // Initialize MPI, if available
    MPIHelper::instance(argc, argv);

    constexpr int dim = 2;
    using Grid = UGGrid<dim>;
    std::shared_ptr<Grid> grid = GmshReader<Grid>::read("l-shape.msh");
    using GridView = Grid::LeafGridView;
    auto gridView = grid->leafGridView();
    // { grid_setup_end }

    // { define_threshold_begin }
    const double estimatedErrorThreshold = 0.1;
    // { define_threshold_end }

    // { make_function_space_begin }
    // Make grid function space
    using Basis = Functions::LagrangeBasis<GridView, 1>;
    auto basis = std::make_shared<Basis>(gridView);

    using Constraints = PDELab::ConformingDirichletConstraints;
    using VectorBackend = PDELab::ISTL::VectorBackend<>;

    using GridFunctionSpace = PDELab::Experimental::GridFunctionSpace<Basis,
                                                                      VectorBackend,
                                                                      Constraints>;

    GridFunctionSpace gridFunctionSpace(basis);
    // { make_function_space_end }

    // Assemble Dirichlet constraints on this space
    // { make_constraints_begin }
    PoissonProblem<GridView, double> poissonProblem;
    PDELab::
        ConvectionDiffusionBoundaryConditionAdapter<decltype(poissonProblem)>
            bctype(poissonProblem);

    using ConstraintsContainer = GridFunctionSpace::ConstraintsContainer<double>::Type;
    ConstraintsContainer constraintsContainer;
    PDELab::constraints(bctype, gridFunctionSpace, constraintsContainer);
    // { make_constraints_end }

    // Make solution vector
    // { make_solution_vector_begin }
    using U = PDELab::Backend::Vector<GridFunctionSpace, double>;
    U u(gridFunctionSpace, 0.0);
    auto g = [](auto x)
    { return (x[0] < 1e-8 || x[1] < 1e-8) ? 0 : 0.5; };
    Functions::interpolate(*basis, PDELab::Backend::native(u), g);
    // { make_solution_vector_end }

    // { refinement_loop_begin }
    std::size_t i = 0; // Loop variable, for data output
    while (true)       // Loop termination condition
    // is near the middle of the loop body
    {
        std::cout << "Iteration: " << i
                  << ", highest level in grid: " << grid->maxLevel()
                  << ", degrees of freedom: " << gridFunctionSpace.globalSize()
                  << std::endl;
        // { refinement_loop_begin_end }

        // Make grid operator
        // { make_assembler_and_solver_begin }
        using LocalOperator = PDELab::ConvectionDiffusionFEM<decltype(poissonProblem),
                                                             GridFunctionSpace::Traits::
                                                                 FiniteElementMap>;
        LocalOperator localOperator(poissonProblem);

        using MatrixBackend = PDELab::ISTL::BCRSMatrixBackend<>;
        MatrixBackend matrixBackend(7); // Expected average number
        // of matrix entries per row
        using GridOperator = PDELab::GridOperator<GridFunctionSpace,
                                                  GridFunctionSpace,
                                                  LocalOperator,
                                                  MatrixBackend,
                                                  double, double, double,
                                                  ConstraintsContainer,
                                                  ConstraintsContainer>;
        GridOperator gridOperator(gridFunctionSpace, constraintsContainer,
                                  gridFunctionSpace, constraintsContainer,
                                  localOperator,
                                  matrixBackend);

        // Select a linear solver backend
        using LinearSolverBackend = PDELab::ISTLBackend_SEQ_CG_SSOR;
        LinearSolverBackend linearSolverBackend(5000, // Maximum number
                                                // of iterations
                                                1); // Verbosity level

        // Select linear problem solver
        using LinearProblemSolver = PDELab::StationaryLinearProblemSolver<GridOperator,
                                                                          LinearSolverBackend,
                                                                          U>;
        LinearProblemSolver linearProblemSolver(gridOperator,
                                                linearSolverBackend,
                                                u,
                                                1e-10);

        // Solve linear problem
        linearProblemSolver.apply();
        // { make_assembler_and_solver_end }

        // File output
        // { vtk_output_begin }
        VTKWriter<GridView> vtkwriter(gridView);
        auto uFunction = Functions::makeDiscreteGlobalBasisFunction<double>(
            *basis,
            PDELab::Backend::native(u));
        vtkwriter.addVertexData(uFunction,
                                VTK::FieldInfo("u",
                                               VTK::FieldInfo::Type::scalar,
                                               1));
        vtkwriter.write("pdelab-adaptivity-result-" + std::to_string(i));
        // { vtk_output_end }

        // Preparation: Define types for the computation of the error estimate eta
        // { estimation_preparation_begin }
        using P0Basis = Functions::LagrangeBasis<GridView, 0>;
        auto p0Basis = std::make_shared<P0Basis>(gridView);

        using P0GridFunctionSpace = PDELab::Experimental::GridFunctionSpace<P0Basis,
                                                                            VectorBackend,
                                                                            PDELab::NoConstraints>;

        P0GridFunctionSpace p0GridFunctionSpace(p0Basis);
        using U0 = PDELab::Backend::Vector<P0GridFunctionSpace, double>;
        U0 etaSquared(p0GridFunctionSpace, 0.0);
        // { estimation_preparation_end }

        // Compute estimated error eta
        // { estimator_gfs_begin }
        using EstimatorLocalOperator = ResidualErrorEstimator<decltype(poissonProblem)>;
        EstimatorLocalOperator estimatorLocalOperator(poissonProblem);
        using EstimatorGridOperator = PDELab::GridOperator<GridFunctionSpace,
                                                           P0GridFunctionSpace,
                                                           EstimatorLocalOperator,
                                                           MatrixBackend,
                                                           double, double, double>;
        EstimatorGridOperator estimatorGridOperator(gridFunctionSpace,
                                                    p0GridFunctionSpace,
                                                    estimatorLocalOperator,
                                                    matrixBackend);
        // { estimator_gfs_end }

        // { estimator_call_begin }
        estimatorGridOperator.residual(u, etaSquared);
        // { estimator_call_end }

        // { termination_begin }
        // Terminate if desired accuracy has been reached
        double totalEstimatedError = std::sqrt(std::accumulate(etaSquared.begin(),
                                                               etaSquared.end(),
                                                               0.0));
        std::cout << "Total estimated error: " << totalEstimatedError << std::endl;

        if (totalEstimatedError < estimatedErrorThreshold)
            break;
        // { termination_end }

        // Adapt the grid locally...
        // { error_fraction_begin }
        double alpha = 0.4;      // Refinement fraction
        double refineThreshold;  // Refinement threshold
        double beta = 0.0;       // Coarsening fraction
        double coarsenThreshold; // Coarsening threshold
        int verbose = 0;         // No screen output

        PDELab::error_fraction(etaSquared,
                               alpha, beta,
                               refineThreshold, coarsenThreshold,
                               verbose);
        // { error_fraction_end }

        // { mark_begin }
        PDELab::mark_grid(*grid, etaSquared, refineThreshold, coarsenThreshold);
        // { mark_end }
        // { adapt_begin }
        PDELab::adapt_grid(*grid, gridFunctionSpace, u, 2);
        // { adapt_end }

        // { dirichlet_reinterpolation_begin }
        // Reassemble the Dirichlet constraints
        PDELab::constraints(bctype, gridFunctionSpace, constraintsContainer);

        // Interpolate the Dirichlet values function on the entire domain!
        PDELab::Backend::Vector<GridFunctionSpace, double>
            dirichletValues(gridFunctionSpace);
        Functions::interpolate(*basis,
                               PDELab::Backend::native(dirichletValues),
                               g);

        // Copy all Dirichlet degrees of freedom into the actual solution vector
        PDELab::copy_constrained_dofs(constraintsContainer, dirichletValues, u);

        // Increment loop counter
        i++;
    }
    // { dirichlet_reinterpolation_end }
}