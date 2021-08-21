// B.13. The Linear Reaction-Diffusion Equation Using a
// DG Method
// File: pdelab-dg-diffusion.cc
// Equation: Linear reaction–diffusion equation
// Discretization: Symmetric Interior Penalty Galerkin, with second-order
// discontinuous Lagrange finite elements
// Grid: Nonconforming quadrilateral grid, implemented with UGGrid
// Solver: CG solver with ILU preconditioner
// Distributed: no
// Discussed in: Chapter 11.2.6

#include "config.h"

#include <vector>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/fvector.hh>

#include <dune/grid/uggrid.hh>
#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>

#include <dune/functions/functionspacebases/lagrangedgbasis.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>

#include <dune/pdelab.hh>

using namespace Dune;

// A local operator for solving the reaction-diffusion equation with SIPG
// { local_operator_begin }
template <class GridView>
class ReactionDiffusionDG
    : public PDELab::LocalOperatorDefaultFlags,
      public Dune::PDELab::FullVolumePattern,
      public Dune::PDELab::FullSkeletonPattern,
      public Dune::PDELab::FullBoundaryPattern
{
    // { local_operator_signature_end }

    // { diameter_begin }
    template <class Geometry>
    static auto diameter(const Geometry &geometry)
    {
        typename Geometry::ctype h = 0.0;
        for (int i = 0; i < geometry.corners(); i++)
            for (int j = i + 1; j < geometry.corners(); j++)
                h = std::max(h, (geometry.corner(j) - geometry.corner(i)).two_norm());

        return h;
    }
    // { diameter_end }

    // { internal_members_begin }
    double kappa_;

    static constexpr int dim = GridView::dimension;

    std::function<double(const FieldVector<double, dim> &)> c_;
    std::function<double(const FieldVector<double, dim> &)> f_;
    std::function<double(const FieldVector<double, dim> &)> neumann_;
    // { internal_members_end }

public:
    // { flags_begin }
    // Residual assembly flags
    enum
    {
        doAlphaVolume = true
    };
    enum
    {
        doAlphaSkeleton = true
    };
    enum
    {
        doAlphaBoundary = true
    };

    // Pattern assembly flags
    enum
    {
        doPatternVolume = true
    };
    enum
    {
        doPatternSkeleton = true
    };
    enum
    {
        doPatternBoundary = true
    };
    // { flags_end }

    // { constructor_begin }
    ReactionDiffusionDG(std::function<double(const FieldVector<double, dim> &)> c,
                        std::function<double(const FieldVector<double, dim> &)> f,
                        std::function<double(const FieldVector<double, dim> &)>
                            neumann,
                        double kappa)
        : kappa_(kappa),
          c_(c),
          f_(f),
          neumann_(neumann)
    {
    }
    // { constructor_end }

    // Volume integral depending on trial and test functions
    template <class EntityGeometry,
              class LocalFunctionSpaceU, class Vector,
              class LocalFunctionSpaceV, class Residual>
    void alpha_volume(const EntityGeometry &entityGeometry,
                      const LocalFunctionSpaceU &localFunctionSpaceU,
                      const Vector &x,
                      const LocalFunctionSpaceV &localFunctionSpaceV,
                      Residual &residual) const
    {
        using TrialFE = typename LocalFunctionSpaceU::Traits::FiniteElementType;
        using TestFE = typename LocalFunctionSpaceV::Traits::FiniteElementType;
        using RangeU = typename TrialFE::Traits::LocalBasisType::Traits::RangeType;
        using RangeV = typename TestFE::Traits::LocalBasisType::Traits::RangeType;
        using GradientU = typename TrialFE::Traits::LocalBasisType::Traits::JacobianType;
        using GradientV = typename TestFE::Traits::LocalBasisType::Traits::JacobianType;
        using size_type = typename LocalFunctionSpaceU::Traits::SizeType;

        // Containers for shape function values and gradients
        std::vector<RangeU> phi(localFunctionSpaceU.size());
        std::vector<RangeV> theta(localFunctionSpaceV.size());

        std::vector<GradientU> localGradPhi(localFunctionSpaceU.size());
        std::vector<GradientV> localGradTheta(localFunctionSpaceV.size());

        std::vector<GradientU> gradPhi(localFunctionSpaceU.size());
        std::vector<GradientV> gradTheta(localFunctionSpaceV.size());

        // Quadrature loop
        const int quadOrder = 2 * localFunctionSpaceU.finiteElement().localBasis().order();
        constexpr auto entityDim = EntityGeometry::Entity::mydimension;
        const auto &quadRule = QuadratureRules<double, entityDim>::rule(entityGeometry.entity().type(), quadOrder);

        for (const auto &quadPoint : quadRule)
        {
            // Evaluate shape functions
            localFunctionSpaceU.finiteElement().localBasis().evaluateFunction(quadPoint.position(), phi);
            localFunctionSpaceV.finiteElement().localBasis().evaluateFunction(quadPoint.position(), theta);

            // Compute value of u
            RangeU u = 0.0;
            for (size_type i = 0; i < localFunctionSpaceU.size(); i++)
                u += x(localFunctionSpaceU, i) * phi[i];

            // Evaluate gradients of shape functions
            localFunctionSpaceU.finiteElement().localBasis().evaluateJacobian(quadPoint.position(), localGradPhi);
            localFunctionSpaceV.finiteElement().localBasis().evaluateJacobian(quadPoint.position(), localGradTheta);

            // Transform gradients of shape functions to element coordinates
            auto jacInvTransp = entityGeometry.geometry().jacobianInverseTransposed(quadPoint.position());
            for (size_type i = 0; i < localFunctionSpaceU.size(); i++)
                jacInvTransp.mv(localGradPhi[i][0], gradPhi[i][0]);

            for (size_type i = 0; i < localFunctionSpaceV.size(); i++)
                jacInvTransp.mv(localGradTheta[i][0], gradTheta[i][0]);

            // Compute gradient of u
            GradientU gradU = 0.0;
            for (size_type i = 0; i < localFunctionSpaceU.size(); i++)
                gradU.axpy(x(localFunctionSpaceU, i), gradPhi[i]);

            // Evaluate reaction and volume terms
            auto globalQuadPosition = entityGeometry.geometry().global(quadPoint.position());
            auto c = c_(globalQuadPosition);
            auto f = f_(globalQuadPosition);

            // Integrate
            auto factor = quadPoint.weight() * entityGeometry.geometry().integrationElement(quadPoint.position());
            for (size_type i = 0; i < localFunctionSpaceV.size(); i++)
                residual.accumulate(localFunctionSpaceV, i, (gradU[0] * gradTheta[i][0] + c * u * theta[i] - f * theta[i]) * factor);
        }
    }

    // Jacobian of volume contribution of the residual
    template <class EntityGeometry,
              class LocalFunctionSpaceU, class Vector,
              class LocalFunctionSpaceV, class Jacobian>
    void jacobian_volume(const EntityGeometry &entityGeometry,
                         const LocalFunctionSpaceU &localFunctionSpaceU,
                         const Vector &x,
                         const LocalFunctionSpaceV &localFunctionSpaceV,
                         Jacobian &jacobian) const
    {
        using TrialFE = typename LocalFunctionSpaceU::Traits::FiniteElementType;
        using TestFE = typename LocalFunctionSpaceV::Traits::FiniteElementType;
        using RangeU = typename TrialFE::Traits::LocalBasisType::Traits::RangeType;
        using RangeV = typename TestFE::Traits::LocalBasisType::Traits::RangeType;
        using GradientU = typename TrialFE::Traits::LocalBasisType::Traits::JacobianType;
        using GradientV = typename TestFE::Traits::LocalBasisType::Traits::JacobianType;
        using size_type = typename LocalFunctionSpaceU::Traits::SizeType;

        // Containers for shape function values and gradients
        std::vector<RangeU> phi(localFunctionSpaceU.size());
        std::vector<RangeV> theta(localFunctionSpaceV.size());

        std::vector<GradientU> localGradPhi(localFunctionSpaceU.size());
        std::vector<GradientV> localGradTheta(localFunctionSpaceV.size());

        std::vector<GradientU> gradPhi(localFunctionSpaceU.size());
        std::vector<GradientV> gradTheta(localFunctionSpaceV.size());

        // Loop over quadrature points
        const int quadOrder = 2 * localFunctionSpaceU.finiteElement().localBasis().order();
        constexpr auto entityDim = EntityGeometry::Entity::mydimension;
        const auto &quadRule = QuadratureRules<double, entityDim>::rule(entityGeometry.geometry().type(), quadOrder);

        for (const auto &quadPoint : quadRule)
        {
            // Evaluate shape functions
            localFunctionSpaceU.finiteElement().localBasis().evaluateFunction(quadPoint.position(), phi);
            localFunctionSpaceV.finiteElement().localBasis().evaluateFunction(quadPoint.position(), theta);

            // Evaluate gradients of shape functions
            localFunctionSpaceU.finiteElement().localBasis().evaluateJacobian(quadPoint.position(), alGradPhi);
            localFunctionSpaceV.finiteElement().localBasis().evaluateJacobian(quadPoint.position(), alGradTheta);

            // Transform gradients of shape functions to real element
            auto jacInvTransp = entityGeometry.geometry().jacobianInverseTransposed(quadPoint.position());
            for (size_type i = 0; i < localFunctionSpaceU.size(); i++)
                jacInvTransp.mv(localGradPhi[i][0], gradPhi[i][0]);

            for (size_type i = 0; i < localFunctionSpaceV.size(); i++)
                jacInvTransp.mv(localGradTheta[i][0], gradTheta[i][0]);

            // Evaluate reaction term
            auto c = c_(entityGeometry.geometry().global(quadPoint.position()));

            // Integrate
            auto factor = quadPoint.weight() * entityGeometry.geometry().integrationElement(quadPoint.ition());
            for (size_type i = 0; i < localFunctionSpaceV.size(); i++)
                for (size_type j = 0; j < localFunctionSpaceU.size(); j++)
                    jacobian.accumulate(localFunctionSpaceV,
                                        i,
                                        localFunctionSpaceU,
                                        j,
                                        (gradPhi[j][0] * gradTheta[i][0] + c * phi[j] * theta[i]) * factor);
        }
    }

    // Skeleton integral for a given configuration u_h
    // Each intersection is only visited once
    // { alpha_skeleton_interface_begin }
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
        // { alpha_skeleton_interface_end }
        using TrialFE = typename LocalFunctionSpaceU::Traits::FiniteElementType;
        using TestFE = typename LocalFunctionSpaceV::Traits::FiniteElementType;
        using LocalBasisU = typename TrialFE::Traits::LocalBasisType;
        using LocalBasisV = typename TestFE::Traits::LocalBasisType;
        using RangeU = typename LocalBasisU::Traits::RangeType;
        using RangeV = typename LocalBasisV::Traits::RangeType;
        using GradientU = typename LocalBasisU::Traits::JacobianType;
        using GradientV = typename LocalBasisV::Traits::JacobianType;
        using size_type = typename LocalFunctionSpaceV::Traits::SizeType;

        // References to inside and outside elements
        const auto &elementInside = intersectionGeometry.inside();
        const auto &elementOutside = intersectionGeometry.outside();

        // Intersection and element geometries
        auto geo = intersectionGeometry.geometry();
        auto geoIn = elementInside.geometry();
        auto geoOut = elementOutside.geometry();

        // Geometries of intersection in local coordinates
        // of elementInside and elementOutside
        auto geoInInside = intersectionGeometry.geometryInInside();
        auto geoInOutside = intersectionGeometry.geometryInOutside();
        // { alpha_skeleton_variable_setup_end }

        // Intersection diameter
        // { compute_diameter_begin }
        auto h = diameter(geo);
        // { compute_diameter_end }

        // { alpha_skeleton_container_setup_begin }
        // Shape function values
        std::vector<RangeU> phiIn(localFunctionSpaceUIn.size());
        std::vector<RangeV> thetaIn(localFunctionSpaceVIn.size());
        std::vector<RangeU> phiOut(localFunctionSpaceUOut.size());
        std::vector<RangeV> thetaOut(localFunctionSpaceVOut.size());

        // Shape function gradients on the reference element
        std::vector<GradientU> localGradPhiIn(localFunctionSpaceUIn.size());
        std::vector<GradientV> localGradThetaIn(localFunctionSpaceVIn.size());
        std::vector<GradientU> localGradPhiOut(localFunctionSpaceUOut.size());
        std::vector<GradientV> localGradThetaOut(localFunctionSpaceVOut.size());

        // Shape function gradients on the grid element
        std::vector<GradientU> gradPhiIn(localFunctionSpaceUIn.size());
        std::vector<GradientV> gradThetaIn(localFunctionSpaceVIn.size());
        std::vector<GradientU> gradPhiOut(localFunctionSpaceUOut.size());
        std::vector<GradientV> gradThetaOut(localFunctionSpaceVOut.size());
        // { alpha_skeleton_container_setup_end }

        // Loop over quadrature points
        // { alpha_skeleton_quadrature_begin }
        const auto quadOrder = 2 * localFunctionSpaceUIn.finiteElement().localBasis().order();
        constexpr auto intersectionDim = IntersectionGeometry::mydimension;
        const auto &quadRule = QuadratureRules<double, intersectionDim>::
            rule(intersectionGeometry.geometry().type(),
                 quadOrder);

        for (const auto &quadPoint : quadRule)
        {
            // { alpha_skeleton_quadrature_end }
            // Quadrature point position in local coordinates of adjacent elements
            auto quadPointLocalIn = geoInInside.global(quadPoint.position());
            auto quadPointLocalOut = geoInOutside.global(quadPoint.position());

            // Evaluate shape functions
            localFunctionSpaceUIn.finiteElement().localBasis().evaluateFunction(quadPointLocalIn, phiIn);
            localFunctionSpaceUOut.finiteElement().localBasis().evaluateFunction(quadPointLocalOut, phiOut);
            localFunctionSpaceVIn.finiteElement().localBasis().evaluateFunction(quadPointLocalIn, thetaIn);
            localFunctionSpaceVOut.finiteElement().localBasis().evaluateFunction(quadPointLocalOut, thetaOut);

            // Evaluate gradients of shape functions
            localFunctionSpaceUIn.finiteElement().localBasis().evaluateJacobian(quadPointLocalIn,
                                                                                localGradPhiIn);
            localFunctionSpaceUOut.finiteElement().localBasis().evaluateJacobian(quadPointLocalOut,
                                                                                 localGradPhiOut);
            localFunctionSpaceVIn.finiteElement().localBasis().evaluateJacobian(quadPointLocalIn,
                                                                                localGradThetaIn);
            localFunctionSpaceVOut.finiteElement().localBasis().evaluateJacobian(quadPointLocalOut,
                                                                                 localGradThetaOut);
            // { alpha_skeleton_shape_functions_end }

            // Transform gradients of shape functions to grid element coordinates
            // { alpha_skeleton_gradient_transform_begin }
            auto jacInvTranspIn = geoIn.jacobianInverseTransposed(quadPointLocalIn);
            for (size_type i = 0; i < localFunctionSpaceUIn.size(); i++)
                jacInvTranspIn.mv(localGradPhiIn[i][0], gradPhiIn[i][0]);
            for (size_type i = 0; i < localFunctionSpaceVIn.size(); i++)
                jacInvTranspIn.mv(localGradThetaIn[i][0], gradThetaIn[i][0]);

            auto jacInvTranspOut = geoOut.jacobianInverseTransposed(quadPointLocalOut);
            for (size_type i = 0; i < localFunctionSpaceUOut.size(); i++)
                jacInvTranspOut.mv(localGradPhiOut[i][0], gradPhiOut[i][0]);
            for (size_type i = 0; i < localFunctionSpaceVOut.size(); i++)
                jacInvTranspOut.mv(localGradThetaOut[i][0], gradThetaOut[i][0]);
            // { alpha_skeleton_gradient_transform_end }

            // { alpha_skeleton_u_u_h_begin }
            // Compute values of u_h
            RangeU uIn(0.0);
            for (size_type i = 0; i < localFunctionSpaceUIn.size(); i++)
                uIn += xIn(localFunctionSpaceUIn, i) * phiIn[i];
            RangeU uOut(0.0);
            for (size_type i = 0; i < localFunctionSpaceUOut.size(); i++)
                uOut += xOut(localFunctionSpaceUOut, i) * phiOut[i];

            // Compute gradients of u_h
            GradientU gradUIn(0.0);
            for (size_type i = 0; i < localFunctionSpaceUIn.size(); i++)
                gradUIn.axpy(xIn(localFunctionSpaceUIn, i), gradPhiIn[i]);
            GradientU gradUOut(0.0);
            for (size_type i = 0; i < localFunctionSpaceUOut.size(); i++)
                gradUOut.axpy(xOut(localFunctionSpaceUOut, i), gradPhiOut[i]);
            // { alpha_skeleton_u_u_h_end }

            // { alpha_skeleton_compute_values_begin }
            // Unit normal from T_in to T_out
            auto n_F = intersectionGeometry.unitOuterNormal(quadPoint.position());

            // Integration factor
            auto factor = quadPoint.weight() * geo.integrationElement(quadPoint.position());

            // Interior penalty term
            auto interiorPenaltyTerm = -0.5 * (gradUIn[0] * n_F + gradUOut[0] * n_F);
            for (size_type i = 0; i < localFunctionSpaceVIn.size(); i++)
                residualIn.accumulate(localFunctionSpaceVIn,
                                      i,
                                      interiorPenaltyTerm * thetaIn[i] * factor);
            for (size_type i = 0; i < localFunctionSpaceVOut.size(); i++)
                residualOut.accumulate(localFunctionSpaceVOut,
                                       i,
                                       -interiorPenaltyTerm * thetaOut[i] * factor);

            // Symmetric interior penalty term
            for (size_type i = 0; i < localFunctionSpaceVIn.size(); i++)
                residualIn.accumulate(localFunctionSpaceVIn,
                                      i,
                                      -0.5 * (gradThetaIn[i][0] * n_F) * (uIn - uOut) * factor);
            for (size_type i = 0; i < localFunctionSpaceVOut.size(); i++)
                residualOut.accumulate(localFunctionSpaceVOut,
                                       i,
                                       -0.5 * (gradThetaOut[i][0] * n_F) * (uIn - uOut) * factor);

            // Coercivity term
            auto coercivityTerm = (kappa_ / h) * (uIn - uOut);
            for (size_type i = 0; i < localFunctionSpaceVIn.size(); i++)
                residualIn.accumulate(localFunctionSpaceVIn,
                                      i,
                                      coercivityTerm * thetaIn[i] * factor);
            for (size_type i = 0; i < localFunctionSpaceVOut.size(); i++)
                residualOut.accumulate(localFunctionSpaceVOut,
                                       i,
                                       -coercivityTerm * thetaOut[i] * factor);
        }
    }
    // { alpha_skeleton_compute_values_end }

    // { jacobian_skeleton_signature_begin }
    template <class IntersectionGeometry,
              class LocalFunctionSpaceU, class Vector,
              class LocalFunctionSpaceV, class Jacobian>
    void jacobian_skeleton(const IntersectionGeometry &intersectionGeometry,
                           const LocalFunctionSpaceU &localFunctionSpaceUIn,
                           const Vector &xIn,
                           const LocalFunctionSpaceV &localFunctionSpaceVIn,
                           const LocalFunctionSpaceU &localFunctionSpaceUOut,
                           const Vector &xOut,
                           const LocalFunctionSpaceV &localFunctionSpaceVOut,
                           Jacobian &jacobianInIn, Jacobian &jacobianInOut,
                           Jacobian &jacobianOutIn, Jacobian &jacobianOutOut)
        const
    // { jacobian_skeleton_signature_end }
    {
        // { jacobian_skeleton_setup_begin }
        // Define types
        using TrialFE = typename LocalFunctionSpaceU::Traits::FiniteElementType;
        using TestFE = typename LocalFunctionSpaceV::Traits::FiniteElementType;
        using LocalBasisU = typename TrialFE::Traits::LocalBasisType;
        using LocalBasisV = typename TestFE::Traits::LocalBasisType;
        using RangeU = typename LocalBasisU::Traits::RangeType;
        using RangeV = typename LocalBasisV::Traits::RangeType;
        using GradientU = typename LocalBasisU::Traits::JacobianType;
        using GradientV = typename LocalBasisV::Traits::JacobianType;
        using size_type = typename LocalFunctionSpaceV::Traits::SizeType;

        // References to inside and outside elements
        const auto &elementInside = intersectionGeometry.inside();
        const auto &elementOutside = intersectionGeometry.outside();

        // Get geometries
        auto geo = intersectionGeometry.geometry();
        auto geoIn = elementInside.geometry();
        auto geoOut = elementOutside.geometry();

        // Geometries of intersection in local coordinates
        // of elementInside and elementOutside
        auto geoInInside = intersectionGeometry.geometryInInside();
        auto geoInOutside = intersectionGeometry.geometryInOutside();

        // Intersection diameter
        auto h = diameter(geo);

        // Shape function values
        std::vector<RangeU> phiIn(localFunctionSpaceUIn.size());
        std::vector<RangeU> phiOut(localFunctionSpaceUOut.size());
        std::vector<RangeV> thetaIn(localFunctionSpaceVIn.size());
        std::vector<RangeV> thetaOut(localFunctionSpaceVOut.size());

        // Shape function gradients on the reference element
        std::vector<GradientU> localGradPhiIn(localFunctionSpaceUIn.size());
        std::vector<GradientU> localGradPhiOut(localFunctionSpaceUOut.size());
        std::vector<GradientV> localGradThetaIn(localFunctionSpaceVIn.size());
        std::vector<GradientV> localGradThetaOut(localFunctionSpaceVOut.size());

        // Shape function gradients on the grid element
        std::vector<GradientU> gradPhiIn(localFunctionSpaceUIn.size());
        std::vector<GradientU> gradPhiOut(localFunctionSpaceUOut.size());
        std::vector<GradientV> gradThetaIn(localFunctionSpaceVIn.size());
        std::vector<GradientV> gradThetaOut(localFunctionSpaceVOut.size());
        // { jacobian_skeleton_setup_end }

        // Loop over quadrature points
        // { jacobian_skeleton_quad_setup_begin }
        auto quadOrder = 2 * localFunctionSpaceUIn.finiteElement().localBasis().order();
        constexpr auto intersectionDim = IntersectionGeometry::mydimension;
        const auto &quadRule = QuadratureRules<double, intersectionDim>::
            rule(intersectionGeometry.geometry().type(),
                 quadOrder);

        for (const auto &quadPoint : quadRule)
        {
            // Position of quadrature point in local coordinates
            // of inside and outside elements
            auto quadPointLocalIn = geoInInside.global(quadPoint.position());
            auto quadPointLocalOut = geoInOutside.global(quadPoint.position());

            // Evaluate shape functions
            localFunctionSpaceUIn.finiteElement().localBasis().evaluateFunction(quadPointLocalIn, phiIn);
            localFunctionSpaceUOut.finiteElement().localBasis().evaluateFunction(quadPointLocalOut, phiOut);

            localFunctionSpaceVIn.finiteElement().localBasis().evaluateFunction(quadPointLocalIn, thetaIn);
            localFunctionSpaceVOut.finiteElement().localBasis().evaluateFunction(quadPointLocalOut, thetaOut);

            // Evaluate gradients of shape functions
            localFunctionSpaceUIn.finiteElement().localBasis().evaluateJacobian(quadPointLocalIn,
                                                                                localGradPhiIn);
            localFunctionSpaceUOut.finiteElement().localBasis().evaluateJacobian(quadPointLocalOut,
                                                                                 localGradPhiOut);

            localFunctionSpaceVIn.finiteElement().localBasis().evaluateJacobian(quadPointLocalIn,
                                                                                localGradThetaIn);
            localFunctionSpaceVOut.finiteElement().localBasis().evaluateJacobian(quadPointLocalOut,
                                                                                 localGradThetaOut);

            // Transform gradients of shape functions to element coordinates
            auto jacInvTranspIn = geoIn.jacobianInverseTransposed(quadPointLocalIn);
            for (size_type i = 0; i < localFunctionSpaceUIn.size(); i++)
                jacInvTranspIn.mv(localGradPhiIn[i][0], gradPhiIn[i][0]);
            for (size_type i = 0; i < localFunctionSpaceVIn.size(); i++)
                jacInvTranspIn.mv(localGradThetaIn[i][0], gradThetaIn[i][0]);

            auto jacInvTranspOut = geoOut.jacobianInverseTransposed(quadPointLocalOut);
            for (size_type i = 0; i < localFunctionSpaceUOut.size(); i++)
                jacInvTranspOut.mv(localGradPhiOut[i][0], gradPhiOut[i][0]);
            for (size_type i = 0; i < localFunctionSpaceVOut.size(); i++)
                jacInvTranspOut.mv(localGradThetaOut[i][0], gradThetaOut[i][0]);
            // { jacobian_skeleton_quad_setup_end }

            // { jacobian_skeleton_compute_matrices_begin }
            // Unit normal from T_in to T_out
            auto n_F = intersectionGeometry.unitOuterNormal(quadPoint.position());

            // Integration factor
            auto factor = quadPoint.weight() * geo.integrationElement(quadPoint.position());

            // Fill jacobianInIn
            for (size_type i = 0; i < localFunctionSpaceVIn.size(); i++)
            {
                for (size_type j = 0; j < localFunctionSpaceUIn.size(); j++)
                {
                    jacobianInIn.accumulate(localFunctionSpaceVIn,
                                            i,
                                            localFunctionSpaceUIn,
                                            j,
                                            (-0.5 * (gradPhiIn[j][0] * n_F) * thetaIn[i] - 0.5 * phiIn[j] * (gradThetaIn[i][0] * n_F) + (kappa_ / h) * phiIn[j] * thetaIn[i]) * factor);
                }
            }

            // Fill jacobianInOut
            for (size_type i = 0; i < localFunctionSpaceVOut.size(); i++)
            {
                for (size_type j = 0; j < localFunctionSpaceUIn.size(); j++)
                {
                    jacobianInOut.accumulate(localFunctionSpaceVIn,
                                             i,
                                             localFunctionSpaceUOut,
                                             j,
                                             (-0.5 * (n_F * gradPhiOut[j][0]) * thetaIn[i] + 0.5 * phiOut[j] * (n_F * gradThetaIn[i][0]) - (kappa_ / h) * phiOut[j] * thetaIn[i]) * factor);
                }
            }

            // Fill jacobianOutIn
            for (size_type i = 0; i < localFunctionSpaceVIn.size(); i++)
            {
                for (size_type j = 0; j < localFunctionSpaceUOut.size(); j++)
                {
                    jacobianOutIn.accumulate(localFunctionSpaceVOut,
                                             i,
                                             localFunctionSpaceUIn,
                                             j,
                                             (0.5 * (n_F * gradPhiIn[j][0]) * thetaOut[i] - 0.5 * phiIn[j] * (n_F * gradThetaOut[i][0]) - (kappa_ / h) * phiIn[j] * thetaOut[i]) * factor);
                }
            }

            // Fill jacobianOutOut
            for (size_type i = 0; i < localFunctionSpaceVOut.size(); i++)
            {
                for (size_type j = 0; j < localFunctionSpaceUOut.size(); j++)
                {
                    jacobianOutOut.accumulate(localFunctionSpaceVOut,
                                              i,
                                              localFunctionSpaceUOut,
                                              j,
                                              (0.5 * (n_F * gradPhiOut[j][0]) * thetaOut[i] + 0.5 * phiOut[j] * (n_F * gradThetaOut[i][0]) + (kappa_ / h) * phiOut[j] * thetaOut[i]) * factor);
                }
            }
        }
    }
    // { jacobian_skeleton_compute_matrices_end }

    // { alpha_boundary_begin }
    // Boundary integral implementing the Neumann term
    template <class IntersectionGeometry,
              class LocalFunctionSpaceU, class Vector,
              class LocalFunctionSpaceV, class Residual>
    void alpha_boundary(const IntersectionGeometry &intersectionGeometry,
                        const LocalFunctionSpaceU &localFunctionSpaceUIn,
                        const Vector &xIn,
                        const LocalFunctionSpaceV &localFunctionSpaceVIn,
                        Residual &residualIn) const
    {
        using TestFE = typename LocalFunctionSpaceV::Traits::FiniteElementType;
        using RangeV = typename TestFE::Traits::LocalBasisType::Traits::RangeType;
        using size_type = typename LocalFunctionSpaceV::Traits::SizeType;

        // Get geometry of intersection in local coordinates of element_inside
        auto geoInInside = intersectionGeometry.geometryInInside();

        std::vector<RangeV> thetaIn(localFunctionSpaceVIn.size());

        // Loop over quadrature points
        const int quadOrder = 2 * localFunctionSpaceVIn.finiteElement().localBasis().order();
        constexpr auto intersectionDim = IntersectionGeometry::mydimension;
        const auto &quadRule = QuadratureRules<double, intersectionDim>::
            rule(intersectionGeometry.geometry().type(),
                 quadOrder);

        for (const auto &quadPoint : quadRule)
        {
            // Position of quadrature point in local coordinates of elements
            auto quadPointLocalIn = geoInInside.global(quadPoint.position());

            // Evaluate basis functions
            localFunctionSpaceVIn.finiteElement().localBasis().evaluateFunction(quadPointLocalIn, thetaIn);

            // Integration factor
            auto factor = quadPoint.weight() * intersectionGeometry.geometry()
                                                   .integrationElement(quadPoint.position());

            // Evaluate Neumann boundary condition
            auto neumannValue = neumann_(intersectionGeometry.geometry()
                                             .global(quadPoint.position()));

            // Integrate
            for (size_type i = 0; i < localFunctionSpaceVIn.size(); i++)
                residualIn.accumulate(localFunctionSpaceVIn,
                                      i,
                                      -1 * neumannValue * thetaIn[i] * factor);
        }
    }
    // { alpha_boundary_end }

    // { jacobian_boundary_begin }
    template <class IntersectionGeometry,
              class LocalFunctionSpaceU, class Vector,
              class LocalFunctionSpaceV, typename Jacobian>
    void jacobian_boundary(const IntersectionGeometry &intersectionGeometry,
                           const LocalFunctionSpaceU &localFunctionSpaceUIn,
                           const Vector &xIn,
                           const LocalFunctionSpaceV &localFunctionSpaceVIn,
                           Jacobian &jacobianInIn) const
    {
        // Does not do anything---the Jacobian of the boundary term is zero.
    }
    // { jacobian_boundary_end }
};

// { main_begin }
int main(int argc, char *argv[])
{
    // Initialize MPI, if present
    MPIHelper::instance(argc, argv);

    constexpr int dim = 2;
    using Grid = UGGrid<dim>;
    std::shared_ptr<Grid> grid = StructuredGridFactory<Grid>::createCubeGrid({0.0, 0.0},
                                                                             {1.0, 1.0},
                                                                             {4, 4});

    // Nonconforming refinement near the domain center
    grid->setClosureType(Grid::NONE); // UGGrid-specific:
    // turn off green closure

    using Position = FieldVector<Grid::ctype, dim>;
    for (const auto &element : elements(grid->leafGridView()))
        if ((element.geometry().center() - Position(0.5)).two_norm() < 0.25)
            grid->mark(1, element);

    grid->preAdapt();
    grid->adapt();
    grid->postAdapt();

    using GridView = Grid::LeafGridView;
    const GridView gridView = grid->leafGridView();
    // { grid_setup_end }

    // { gfs_setup_begin }
    // Make dune-functions basis
    using DGBasis = Functions::LagrangeDGBasis<GridView, 2>;
    auto dgBasis = std::make_shared<DGBasis>(gridView);

    // Make dune-pdelab grid function space
    using VectorBackend = PDELab::ISTL::VectorBackend<>;
    using GridFunctionSpace = PDELab::Experimental::GridFunctionSpace<DGBasis,
                                                                      VectorBackend,
                                                                      PDELab::NoConstraints>;

    GridFunctionSpace gridFunctionSpace(dgBasis);
    // { gfs_setup_end }

    // { make_parameter_functions_begin }
    // Reaction strength
    auto c = [](const FieldVector<double, dim> &x)
    {
        return 10.0;
    };

    // Source term
    auto f = [](const FieldVector<double, dim> &x)
    {
        const auto center = FieldVector<double, dim>(0.5);
        auto distanceToCenter = (x - center).two_norm();
        return (distanceToCenter <= 0.25) ? -10.0 : 10.0;
    };

    // Neumann boundary data
    auto neumann = [](const FieldVector<double, dim> &x)
    {
        return (x[0] >= 0.999 && x[1] >= 0.5) ? -2 : 0;
    };
    // { make_parameter_functions_end }

    // { assembler_setup_begin }
    // Make local operator
    using LocalOperator = ReactionDiffusionDG<GridView>;
    double kappa = 3;
    LocalOperator localOperator(c, f, neumann, kappa);
    using MatrixBackend = PDELab::ISTL::BCRSMatrixBackend<>;
    MatrixBackend matrixBackend(45); // 45: Expected approximate number
    // of nonzeros per matrix row

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
    // { assembler_setup_end }

    // { construct_and_run_solver_begin }
    // Make a vector of degrees of freedom, and initialize it with 0
    using VectorContainer = typename GridOperator::Traits::Domain;
    VectorContainer u(gridFunctionSpace, 0.0);

    // Make linear solver and solve problem
    using LinearSolverBackend = PDELab::ISTLBackend_SEQ_CG_ILU0;
    LinearSolverBackend linearSolverBackend(50, // Maximum number
                                            // of iterations
                                            2); // Solver verbosity

    using LinearProblemSolver = PDELab::StationaryLinearProblemSolver<GridOperator,
                                                                      LinearSolverBackend,
                                                                      VectorContainer>;

    LinearProblemSolver linearProblemSolver(gridOperator,
                                            linearSolverBackend,
                                            u,
                                            1e-10); // Desired relative
    // residual reduction
    linearProblemSolver.apply();
    // { construct_and_run_solver_end }

    // { output_result_begin }
    // Make a discrete function from the FE basis and the coefficient vector
    auto uFunction = Functions::makeDiscreteGlobalBasisFunction<double>(
        *dgBasis,
        PDELab::Backend::native(u));

    // We need to subsample, because VTKWriter
    // cannot natively display second-order functions.
    SubsamplingVTKWriter<GridView> vtkWriter(gridView,
                                             refinementLevels(5));
    vtkWriter.addVertexData(uFunction,
                            VTK::FieldInfo("u",
                                           VTK::FieldInfo::Type::scalar,
                                           1));
    vtkWriter.write("pdelab-dg-diffusion-result");
    // { output_result_end }
}