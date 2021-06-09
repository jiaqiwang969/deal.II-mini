//include/deal.II-translator/fe/fe_bernardi_raugel_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2021 by the deal.II authors
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

#ifndef dealii_fe_bernardi_raugel_h
#define dealii_fe_bernardi_raugel_h

#include <deal.II/base/config.h>

#include <deal.II/base/geometry_info.h>
#include <deal.II/base/polynomial.h>
#include <deal.II/base/polynomials_bernardi_raugel.h>
#include <deal.II/base/tensor_product_polynomials.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_poly_tensor.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN

/**
 * Bernardi-Raugel元素。
 * 该类实现了非标准的Bernardi-Raugel（BR）元素，可作为斯托克斯方程的稳定速度/压力对的一部分。BR元素可以看作是
 * $Q_1^d$
 * 元素的丰富版本，在每个边（2D）或面（3D）上增加了气泡函数，或者看作是
 * $Q_2^d$ 元素的缩小版本。它解决了 $Q_1^d\times Q_0$
 * 的组合不是inf-sup稳定的（需要更大的速度空间），以及
 * $Q_2^d\times Q_1$
 * 的组合是稳定的，但是次优的，因为速度空间相对于压力空间太大，无法提供与大量速度未知数的代价相称的额外精度。
 * 该元素在以下论文中被介绍。
 *
 * @code{.bib}
 * @article{BR85,
 * author    = {Christine Bernardi and Genevi{\`e}ve Raugel},
 * title     = {Analysis of some finite elements for the {S}tokes problem},
 * journal   = {Mathematics of Computation},
 * publisher = {American Mathematical Society ({AMS})},
 * volume    = {44},
 * number    = {169},
 * pages     = {71--79},
 * year      = {1985},
 * doi       = {10.1090/s0025-5718-1985-0771031-7},
 * url       = {https://doi.org/10.1090/s0025-5718-1985-0771031-7}
 * }
 * @endcode
 *
 *
 *  <h3>Degrees of freedom</h3>
 * BR1元素在每个顶点上有<i>dim</i>个自由度，在每个面上有1个。形状函数由每个顶点上支持的
 * $(Q_1)^d$
 * 形状函数排序，根据GeometryInfo中元素上的顶点排序增加，然后气泡函数按照PolynomialsBernardiRaugel中给出的排序。
 * 这个元素只有1度（度数 $p=1$
 * ），因为它对斯托克斯问题产生了一个LBB稳定对BR1-P0，其度数比Taylor-Hood元素低。该对元素有时被称为富集的P1-P0元素或缩小的P2-P0元素。
 * 在目前的实现中，该元素不支持悬挂节点或多栅格。
 * 一些数值实验表明，在 step-20
 * 中使用BR1-Q0对混合拉普拉斯方程时，该元素可能以一阶精度收敛。
 *
 *
 */
template <int dim>
class FE_BernardiRaugel : public FE_PolyTensor<dim>
{
public:
  /**
   * 度数为 @p p. 的Bernardi-Raugel元的构造函数
   * 唯一支持的度数为1。  @arg  p: 元素 $p=1$ 的度数为 $BR_1$
   * 。
   *
   */
  FE_BernardiRaugel(const unsigned int p = 1);

  /**
   * 返回一个唯一标识有限元的字符串。该类返回<tt>FE_BR<dim>(degree)</tt>，其中
   * @p dim 和 @p degree 由适当的值代替。
   *
   */
  virtual std::string
  get_name() const override;

  virtual std::unique_ptr<FiniteElement<dim, dim>>
  clone() const override;

  // documentation inherited from the base class
  virtual void
  convert_generalized_support_point_values_to_dof_values(
    const std::vector<Vector<double>> &support_point_values,
    std::vector<double> &              nodal_values) const override;

private:
  /**
   * 仅供内部使用。它的全名是 @p get_dofs_per_object_vector
   * 函数，它创建了 @p dofs_per_object
   * 向量，在构造函数中需要传递给 @p
   * FiniteElementData的构造函数。
   *
   */
  static std::vector<unsigned int>
  get_dpo_vector();

  /**
   * 初始化 FiniteElement<dim>::generalized_support_points 和 FiniteElement<dim>::generalized_face_support_points 字段。从构造函数中调用。更多信息请参见 @ref GlossGeneralizedSupport "广义支持点的词汇表条目"
   * 。
   *
   */
  void
  initialize_support_points();

  /**
   * 初始化包络模式和符号变化模式。
   * @note
   * 这个函数还没有完全充满正确的实现。它需要在未来的版本中统一实现，以便在包含有翻转面的单元格的网格上工作。
   *
   */
  void
  initialize_quad_dof_index_permutation_and_sign_change();
};

DEAL_II_NAMESPACE_CLOSE

#endif


