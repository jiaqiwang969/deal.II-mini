// B.12. p-Laplace Problem Using dune-pdelab and Finite
// Elements
// File: pdelab-p-laplace.cc
// Equation: p-Laplace equation with linear reaction term
// Discretization: Second-order Lagrange finite elements
// Grid: Uniform quadrilateral grid, implemented with YaspGrid
// Solver: Newton solver with line search
// Distributed: no
// Discussed in: Chapter 11.2.4

#include "config.h"

#include <iostream>
#include <vector>

#include <dune/common/parallel/mpihelper.hh>

#include <dune/geometry/quadraturerules.hh>

#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>
#include <dune/pdelab.hh>
using namespace Dune;
/** Local operator for solving the equation
 *
 * - \Delta_p u + cu = f in \Omega
 * \nabla u \cdot n = 0 on \partial\Omega
 *
 * with conforming finite elements on all types of grids in any dimension
 *
 */
// { localoperator_begin }
template <int dim>
class PLaplaceLocalOperator : public PDELab::LocalOperatorDefaultFlags,
                              public PDELab::FullVolumePattern
{
public:
    // Pattern assembly flags
    static const bool doPatternVolume = true;
    // Residual assembly flags
    static const bool doAlphaVolume = true;
    // { localoperator_header_end }
    // { data_and_constructor_begin }
    // The order of the p-Laplace term
    const double p_;
    // Source term
    std::function<double(const FieldVector<double, dim> &)> f_;
    // Reaction term
    std::function<double(const FieldVector<double, dim> &)> c_;
    PLaplaceLocalOperator(double p,
                          const std::function<double(const FieldVector<double, dim> &)> &f,
                          const std::function<double(const FieldVector<double, dim> &)> &c)
        : p_(p),
          f_(f),
          c_(c)
    {
    }
    // { data_and_constructor_end }
    // Volume contribution of the residual
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
        // Extract types of shape function values and gradients
        using TrialFE = typename LocalFunctionSpaceU::Traits::FiniteElementType;
        using TestFE = typename LocalFunctionSpaceV::Traits::FiniteElementType;
        using LocalBasisU = typename TrialFE::Traits::LocalBasisType;
        using LocalBasisV = typename TestFE::Traits::LocalBasisType;
        using RangeU = typename LocalBasisU::Traits::RangeType;
        using RangeV = typename LocalBasisV::Traits::RangeType;
        using GradientU = typename LocalBasisU::Traits::JacobianType;
        using GradientV = typename LocalBasisV::Traits::JacobianType;
        // { alpha_volume_types_end }
        // { alpha_volume_container_setup }
        std::vector<RangeU> phi(localFunctionSpaceU.size());
        std::vector<RangeV> theta(localFunctionSpaceV.size());
        std::vector<GradientU> localGradPhi(localFunctionSpaceU.size());
        std::vector<GradientV> localGradTheta(localFunctionSpaceV.size());
        std::vector<GradientU> gradPhi(localFunctionSpaceU.size());
        std::vector<GradientV> gradTheta(localFunctionSpaceV.size());
        // { alpha_volume_setup_end }
        // Loop over quadrature points
        // { alpha_volume_quad_begin }
        int intOrder = (localFunctionSpaceU.finiteElement().localBasis().order() - 1) * (p_ - 1) + (localFunctionSpaceV.finiteElement().localBasis().order() - 1);
        const auto &quadRule = QuadratureRules<double, dim>::rule(entityGeometry.entity().type(),
                                                                  intOrder);

        for (const auto &quadPoint : quadRule)
        {
            // { alpha_volume_quad_setup_end }

            // { evaluate_shape_functions_begin }
            // Evaluate basis functions on reference element
            localFunctionSpaceU.finiteElement().localBasis().evaluateFunction(quadPoint.position(), phi);
            localFunctionSpaceV.finiteElement().localBasis().evaluateFunction(quadPoint.position(), theta);

            // Evaluate gradient of basis functions on reference element
            localFunctionSpaceU.finiteElement().localBasis().evaluateJacobian(quadPoint.position(), localGradPhi);
            localFunctionSpaceV.finiteElement().localBasis().evaluateJacobian(quadPoint.position(), localGradTheta);

            // Transform gradients from reference element to grid element
            const auto jacInvTransp = entityGeometry.geometry()
                                          .jacobianInverseTransposed(quadPoint.position());

            for (std::size_t i = 0; i < localFunctionSpaceU.size(); i++)
                jacInvTransp.mv(localGradPhi[i][0], gradPhi[i][0]);

            for (std::size_t i = 0; i < localFunctionSpaceV.size(); i++)
                jacInvTransp.mv(localGradTheta[i][0], gradTheta[i][0]);
            // { evaluate_shape_functions_end }

            // { alpha_volume_compute_u_uh_begin }
            // Compute u at integration point
            RangeU u = 0;
            for (std::size_t i = 0; i < localFunctionSpaceU.size(); ++i)
                u += x(localFunctionSpaceU, i) * phi[i];

            // Compute gradient of u
            GradientU gradU = 0;
            for (std::size_t i = 0; i < localFunctionSpaceU.size(); ++i)
                gradU.axpy(x(localFunctionSpaceU, i), gradPhi[i]);
            // { alpha_volume_compute_u_uh_end }

            // { compute_alpha_volume_term_begin }
            auto globalPos = entityGeometry.geometry().global(quadPoint.position());

            // Integrate |ruh|(p−2)hruh,rii + cuhi − fi
            auto factor = quadPoint.weight() * entityGeometry.geometry()
                                                   .integrationElement(quadPoint.position());

            for (std::size_t i = 0; i < localFunctionSpaceV.size(); ++i)
            {
                auto value = (std::pow(gradU[0].two_norm2(), 0.5 * (p_ - 2)) * (gradU[0] * gradTheta[i][0]) + c_(globalPos) * u * theta[i] - f_(globalPos) * theta[i]) * factor;
                residual.accumulate(localFunctionSpaceV, i, value);
            }
            // { compute_alpha_volume_term_end }
        }
    }

    // Jacobian of volume term
    // { jacobian_volume_begin }
    template <class EntityGeometry,
              class LocalFunctionSpaceU, class Vector,
              class LocalFunctionSpaceV, class Jacobian>
    void jacobian_volume(const EntityGeometry &entityGeometry,
                         const LocalFunctionSpaceU &localFunctionSpaceU,
                         const Vector &x,
                         const LocalFunctionSpaceV &localFunctionSpaceV,
                         Jacobian &jacobian) const
    {
        // { jacobian_volume_signature_end }
        using TrialFE = typename LocalFunctionSpaceU::Traits::FiniteElementType;
        using TestFE = typename LocalFunctionSpaceU::Traits::FiniteElementType;
        using LocalBasisU = typename TrialFE::Traits::LocalBasisType;
        using LocalBasisV = typename TestFE::Traits::LocalBasisType;
        using RangeU = typename LocalBasisU::Traits::RangeType;
        using RangeV = typename LocalBasisV::Traits::RangeType;
        using GradientU = typename LocalBasisU::Traits::JacobianType;
        using GradientV = typename LocalBasisV::Traits::JacobianType;

        std::vector<RangeU> phi(localFunctionSpaceU.size());
        std::vector<RangeV> theta(localFunctionSpaceV.size());

        std::vector<GradientU> localGradPhi(localFunctionSpaceU.size());
        std::vector<GradientV> localGradTheta(localFunctionSpaceV.size());

        std::vector<GradientU> gradPhi(localFunctionSpaceU.size());
        std::vector<GradientV> gradTheta(localFunctionSpaceV.size());
        // { jacobian_volume_setup_end }

        // Loop over quadrature points
        // { jacobian_volume_quad_loop_begin }
        int intOrder = (localFunctionSpaceU.finiteElement().localBasis().order() - 1) * (p_ - 1) * (localFunctionSpaceV.finiteElement().localBasis().order() - 1);
        const auto &quadRule = QuadratureRules<double, dim>::rule(entityGeometry.entity().type(),
                                                                  intOrder);

        for (const auto &quadPoint : quadRule)
        {
            // Evaluate shape functions on reference element
            localFunctionSpaceU.finiteElement().localBasis().evaluateFunction(quadPoint.position(), phi);
            localFunctionSpaceV.finiteElement().localBasis().evaluateFunction(quadPoint.position(), theta);

            // Evaluate gradients of shape functions on reference element
            localFunctionSpaceU.finiteElement().localBasis().evaluateJacobian(quadPoint.position(), localGradPhi);
            localFunctionSpaceV.finiteElement().localBasis().evaluateJacobian(quadPoint.position(), localGradTheta);

            // Transform gradients from reference element to grid element
            const auto jacInvTransp = entityGeometry.geometry()
                                          .jacobianInverseTransposed(quadPoint.position());

            for (std::size_t i = 0; i < localFunctionSpaceU.size(); i++)
                jacInvTransp.mv(localGradPhi[i][0], gradPhi[i][0]);

            for (std::size_t i = 0; i < localFunctionSpaceV.size(); i++)
                jacInvTransp.mv(localGradTheta[i][0], gradTheta[i][0]);
            // { jacobian_volume_evaluate_shape_functions_end }

            // { jacobian_volume_u_uh_begin }
            // Compute gradient of u
            GradientU gradU(0.0);
            for (std::size_t i = 0; i < localFunctionSpaceU.size(); ++i)
                gradU.axpy(x(localFunctionSpaceU, i), gradPhi[i]);
            // { jacobian_volume_u_uh_end }

            // { compute_jacobian_volume_term_begin }
            auto factor = quadPoint.weight() * entityGeometry.geometry()
                                                   .integrationElement(quadPoint.position());
            auto globalPos = entityGeometry.geometry().global(quadPoint.position());

            // Integrate
            for (std::size_t i = 0; i < localFunctionSpaceV.size(); i++)
                for (std::size_t j = 0; j < localFunctionSpaceU.size(); j++)
                {
                    auto value = (p_ - 2) * std::pow(gradU[0].two_norm2(), 0.5 * p_ - 2) * (gradU[0] * gradPhi[j][0]) * (gradU[0] * gradTheta[i][0]);
                    value += std::pow(gradU[0].two_norm2(), 0.5 * p_ - 1) * (gradTheta[i][0] * gradPhi[j][0]);
                    value += c_(globalPos) * theta[i] * phi[j];
                    jacobian.accumulate(localFunctionSpaceV, i,
                                        localFunctionSpaceU, j,
                                        value * factor);
                }
            // { compute_jacobian_volume_term_end }
        }
    }
};

// { main_begin }
int main(int argc, char *argv[])
{
    // Initialize MPI
    MPIHelper::instance(argc, argv);

    constexpr int dim = 2;
    using Grid = YaspGrid<dim>;
    FieldVector<double, dim> upper = {1.0, 1.0};
    std::array<int, dim> elements = {8, 8};
    Grid grid(upper, elements);

    using GridView = Grid::LeafGridView;
    GridView gridView = grid.leafGridView();
    // { grid_setup_end }

    // Make grid function space
    // { make_grid_function_space_begin }
    using Basis = Functions::LagrangeBasis<GridView, 2>;
    auto basis = std::make_shared<Basis>(gridView);

    using VectorBackend = PDELab::ISTL::VectorBackend<>;

    using GridFunctionSpace = PDELab::Experimental::GridFunctionSpace<Basis,
                                                                      VectorBackend,
                                                                      PDELab::NoConstraints>;

    GridFunctionSpace gridFunctionSpace(basis);
    // { make_grid_function_space_end }

    // { make_coefficient_functions_begin }
    auto f = [](const FieldVector<double, dim> &x)
    {
        FieldVector<double, dim> center(0.5);
        auto distanceToCenter = (x - center).two_norm();
        return (distanceToCenter <= 0.25) ? 100 : 0;
    };

    auto c = [](const FieldVector<double, dim> &x)
    {
        return 100.0;
    };

    double p = 5.0;
    // { make_coefficient_functions_end }

    // Make grid operator
    // { make_local_operator_begin }
    using LocalOperator = PLaplaceLocalOperator<dim>;
    LocalOperator localOperator(p, f, c);
    // { make_local_operator_end }

    // { make_grid_operator_begin }
    using MatrixBackend = PDELab::ISTL::BCRSMatrixBackend<>;
    MatrixBackend matrixBackend(25); // 25: Expected average number of entries
    // per matrix row

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
    // { make_grid_operator_end }

    // Select a linear solver backend
    // { solver_preparations_begin }
    using LinearSolverBackend = PDELab::ISTLBackend_SEQ_CG_SSOR;
    LinearSolverBackend linearSolverBackend(5000, // Max. number of iterations
                                            2);   // Verbosity level

    // Make vector of coefficients
    using U = PDELab::Backend::Vector<GridFunctionSpace, double>;
    U u(gridFunctionSpace, 0.0); // Initial iterate
    // { solver_preparations_end }

    // Solve nonlinear problem
    // { newton_begin }
    PDELab::NewtonMethod<GridOperator, LinearSolverBackend>
        newtonMethod(gridOperator,
                     linearSolverBackend);

    newtonMethod.setReduction(1e-10); // Requested reduction
    // of the residual
    newtonMethod.setMinLinearReduction(1e-4); // Minimal reduction of the
    // inner linear solver
    newtonMethod.setVerbosityLevel(2); // Solver verbosity
    newtonMethod.apply(u);
    // { newton_end }

    // Write VTK output file
    // { vtk_writing_begin }
    // Make a discrete function from the finite element basis
    // and the coefficient vector
    auto uFunction = Functions::makeDiscreteGlobalBasisFunction<double>(
        *basis,
        PDELab::Backend::native(u));

    // We need to subsample, because the dune-grid VTKWriter
    // cannot natively display second-order functions.
    SubsamplingVTKWriter<GridView> vtkwriter(gridView, refinementLevels(2));
    vtkwriter.addVertexData(uFunction,
                            VTK::FieldInfo("u",
                                           VTK::FieldInfo::Type::scalar,
                                           1));
    vtkwriter.write("pdelab-p-laplace-result");
    // { vtk_writing_end }
}