//include/deal.II-translator/fe/fe_rt_bubbles_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2021 by the deal.II authors
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

#ifndef dealii_fe_raviart_thomas_bubbles_h
#define dealii_fe_raviart_thomas_bubbles_h

#include <deal.II/base/config.h>

#include <deal.II/base/geometry_info.h>
#include <deal.II/base/polynomial.h>
#include <deal.II/base/polynomials_rt_bubbles.h>
#include <deal.II/base/tensor_product_polynomials.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_poly_tensor.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN

/**
 * 该类实现了一个卷曲增强的Raviart-Thomas元素，符合<i>H<sup>div</sup></i>空间。节点函数被定义为Gauss-Lobatto点中的点值。这些元素产生的矢量场在网格单元之间具有连续的法向分量。这个有限元的目的是在使用适当的正交规则时，将节点周围的自由度之间的相互作用本地化，导致一个块对角的质量矩阵（甚至是全张量系数）。
 * 元素是通过对经典的Raviart-Thomas元素的富集与额外的卷曲来定义的，这样就保留了<i>H<sup>div</sup></i>的一致性，而且阶数为k的FE_RT_Bubbles的总自由度数等于阶数为<i>dim</i>的FE_Q副本的自由度数。
 *
 *
 * @note
 * 与Raviart-Thomas不同，这种增强型有限元的最低阶数是1，即
 * $k \ge 1$  。
 * 阶数<i>k</i>的FE_RT_Bubbles的匹配压力空间是阶数<i>k-1</i>的FE_DGQ。在精确积分的情况下，对于矢量变量来说，这一对产生
 * $(k+1)$  -st order of convergence in  $L_2$
 * -norm，对于标量变量来说，产生 $k$  -th order in  $L_2$
 * -norm（与 $BDM_k \times P_{k-1}$  相同）。
 * 对于这个增强的Raviart-Thomas元素，节点值不是相对于某些多项式的单元和面矩，而是Gauss-Lobatto正交点的值。根据一个单元的边（面）的自然排序，边（面）上的节点值首先被评估。内部自由度最后评估。
 * 对于度数为<i>k</i>的RT-Bubbles元素，我们在每个面上选择<i>(k+1)<sup>dim-1</sup></i>个Gauss-Lobatto点。这些点在面的方向上是按字母顺序排列的。在单元格的内部，使用各向异性的Gauss-Lobatto公式计算数值，进行积分。使用这个相同的正交规则组装的质量矩阵是块状对角线的，块状对角线对应于正交点。更多细节见<i><a
 * href="https://arxiv.org/abs/1710.06742">"Higher order multipoint flux mixed
 * finite element methods on quadrilaterals and hexahedra"</a><a
 * href="https://arxiv.org/abs/1710.06742">"Higher order multipoint flux mixed
 * finite element methods on quadrilaterals and hexahedra"</a></i>。
 * <i>2D</i>中的 $k=3$ 和<i>3D</i>中的 $k=2$
 * 度的元素显示在下面的数字中（填充的箭头表示DoF，对于这些DoF需要跨边（<i>3D</i>中的面）的连续性）。
 * <table> <tr> <td align="center">
 @image html rtbubbles.png
 * </td></tr>
 *
 * <tr> <td align="center"> Left
 *
 * - $2D,\,k=3$, right
 *
 * - $3D,\,k=2$.</td></tr> </table>
 * @todo  实施限制矩阵
 *
 *
 */
template <int dim>
class FE_RT_Bubbles : public FE_PolyTensor<dim>
{
public:
  /**
   * 度的RT_Bubbles元素的构造器  @p k.
   *
   */
  FE_RT_Bubbles(const unsigned int k);

  /**
   * 返回一个唯一标识有限元的字符串。该类返回<tt>FE_RT_Bubbles<dim>(degree)</tt>，其中
   * @p dim 和 @p 度被适当的值取代。
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
   * 仅供内部使用。它的全称是 @p get_dofs_per_object_vector
   * 函数，它创建了 @p dofs_per_object
   * 向量，在构造函数中需要传递给 @p
   * FiniteElementData的构造函数。
   *
   */
  static std::vector<unsigned int>
  get_dpo_vector(const unsigned int degree);

  /**
   * 计算用于传递给基类构造函数的 @p restriction_is_additive
   * 字段的向量。
   *
   */
  static std::vector<bool>
  get_ria_vector(const unsigned int degree);

  /**
   * 初始化 FiniteElement<dim>::generalized_support_points 和 FiniteElement<dim>::generalized_face_support_points 字段。从构造函数中调用。    更多信息请参见 @ref GlossGeneralizedSupport "广义支持点的词汇条"
   * 。
   *
   */
  void
  initialize_support_points(const unsigned int rt_degree);

  /**
   * 初始化包络模式和符号变化模式。
   *
   */
  void
  initialize_quad_dof_index_permutation_and_sign_change();
};


DEAL_II_NAMESPACE_CLOSE

#endif


