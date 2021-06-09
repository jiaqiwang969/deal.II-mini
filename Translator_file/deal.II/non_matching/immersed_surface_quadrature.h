//include/deal.II-translator/non_matching/immersed_surface_quadrature_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

#ifndef dealii_non_matching_immersed_surface_quadrature
#define dealii_non_matching_immersed_surface_quadrature

#include <deal.II/base/config.h>

#include <deal.II/base/point.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/base/tensor.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN
namespace NonMatching
{
  /**
   * 该类定义了一个正交公式，用于对浸没在单元格中的定向曲面 $\hat{S}$ 进行积分。所谓浸入，是指曲面可以以任意方式与单元格相交。正交公式由一组正交点  $\hat{x}_q$  、权重  $w_q$  和归一化表面法线  $\hat{n}_q$  描述。    我们通常要计算实空间中的曲面积分。  一个与实空间中的单元格 $K$ 相交的曲面 $S$ ，可以被映射到与单元格 $\hat{K}$ 相交的曲面 $\hat{S}$ 。  因此，实空间中 $S\cap K$ 上的曲面积分可以根据@f[
   * \int_{S\cap K} f(x) dS =
   * \int_{S\cap K} f(x) |d\bar{S}| =
   * \int_{\hat{S}\cap\hat{K}} f(F_{K}(\hat{x})) \det(J) |\left( J^{-1} \right
   * )^T d\hat{S}|,
   * @f]转换为 $\hat{S} \cap \hat{K}$ 上的曲面积分，其中 $F_K$ 是从参考空间到实空间的映射， $J$ 是其雅各布系数。这种转换是可能的，因为连续表面元素是矢量。  $d\bar{S}, d\hat{S} \in \mathbb{R}^{dim}$ 与 $S$ 和 $\hat{S}$ 的法线平行。所以为了计算实空间的积分，我们需要法线的信息来做转换。    因此，除了存储点和权重之外，这个四分法还存储了每个四分法点的归一化法线。这可以看作是为每个正交点存储了一个离散的曲面元素，@f[
   * \Delta \hat{S}_q \dealcoloneq w_q \hat{n}_q \approx d\hat{S}(\hat{x}_q),
   * @f]。然后，实空间的曲面积分将被近似为@f[
   * \int_{S\cap K} f(x) dS \approx \sum_{q} f \left(F_{K}(\hat{x}_{q})
   * \right) \det(J_q) |\left( J_q^{-1} \right)^T \hat{n}_q| w_q.
   @f]  @image html immersed_surface_quadrature.svg  。
   *
   */
  template <int dim>
  class ImmersedSurfaceQuadrature : public Quadrature<dim>
  {
  public:
    /**
     * 默认构造函数用于初始化没有正交点的正交。
     *
     */
    ImmersedSurfaceQuadrature() = default;

    /**
     * 从点、权重和表面法线的向量构造一个正交公式。这些点、权重和法线应该是关于参考空间的，法线应该被归一化。
     *
     */
    ImmersedSurfaceQuadrature(const std::vector<Point<dim>> &    points,
                              const std::vector<double> &        weights,
                              const std::vector<Tensor<1, dim>> &normals);

    /**
     * 通过一个额外的正交点扩展给定的公式。
     * 该点、权重和法线应该是关于参考空间的，法线应该被归一化。
     * 这个函数的存在是因为沉浸式正交规则的构造可能相当复杂。通常情况下，构造是通过将单元划分为区域并在每个区域上单独构造点来完成的。这可能会使从构造器中创建正交变得很麻烦，因为在创建对象时必须知道所有的正交点。
     * @note 这个函数应该只在构建正交公式时使用。
     *
     */
    void
    push_back(const Point<dim> &    point,
              const double          weight,
              const Tensor<1, dim> &normal);

    /**
     * 返回一个对<tt>i</tt>第1个表面法线的引用。
     *
     */
    const Tensor<1, dim> &
    normal_vector(const unsigned int i) const;

    /**
     * 返回对整个法线%向量的引用。
     *
     */
    const std::vector<Tensor<1, dim>> &
    get_normal_vectors() const;

  protected:
    /**
     * 每个正交点的表面法线的%向量。
     *
     */
    std::vector<Tensor<1, dim>> normals;
  };

} // namespace NonMatching
DEAL_II_NAMESPACE_CLOSE

#endif


