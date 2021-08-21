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