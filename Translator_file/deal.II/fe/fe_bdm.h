//include/deal.II-translator/fe/fe_bdm_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2021 by the deal.II authors
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

#ifndef dealii_fe_bdm_h
#define dealii_fe_bdm_h

#include <deal.II/base/config.h>

#include <deal.II/base/geometry_info.h>
#include <deal.II/base/polynomial.h>
#include <deal.II/base/polynomials_bdm.h>
#include <deal.II/base/tensor_product_polynomials.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_poly_tensor.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN

/**
 * 布雷兹-道格拉斯-马里尼元素。 <h3>Degrees of freedom</h3>
 * @todo
 * 三维版本表现出一些数值上的不稳定性，特别是对于高阶
 * @todo  限制矩阵丢失。
 * 阶<i>k</i>的FE_BDM的匹配压力空间是阶<i>k-1</i>的FE_DGP元素。
 * 阶数为 @p p
 * 的BDM元素在每个面上都有<i>p+1</i>个自由度。这些都是在每个面上的<i>p+1</i>高斯点中实现的函数值。
 * 此外，对于大于或等于2的阶，我们有额外的<i>p(p-1)</i>，即<i>P<sub>p</sub></i>中的向量值多项式的数量，内部自由度。这些是单元格中<i>p<sup>2</sup></i>高斯点中第一个<i>p(p-1)/2</i>的向量函数值。
 *
 *
 * @ingroup fe
 *
 *
 */
template <int dim>
class FE_BDM : public FE_PolyTensor<dim>
{
public:
  /**
   * 程度为 @p p. 的BDM元素的构造函数
   *
   */
  FE_BDM(const unsigned int p);

  /**
   * 返回一个唯一标识有限元的字符串。该类返回<tt>FE_BDM<dim>(degree)</tt>，其中
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
  get_dpo_vector(const unsigned int degree);

  /**
   * 计算用于传递给基类构造函数的 @p restriction_is_additive
   * 字段的向量。
   *
   */
  static std::vector<bool>
  get_ria_vector(const unsigned int degree);
  /**
   * 初始化 FiniteElement<dim>::generalized_support_points 和 FiniteElement<dim>::generalized_face_support_points 字段。从构造函数中调用。更多信息请参见 @ref GlossGeneralizedSupport "关于广义支持点的词汇条"
   * 。
   *
   */
  void
  initialize_support_points(const unsigned int bdm_degree);
  /**
   * 作为测试函数需要的多项式的面支持点中的值。外侧向量由正交点索引，内部由测试函数索引。测试函数空间是PolynomialsP<dim-1>。
   *
   */
  std::vector<std::vector<double>> test_values_face;
  /**
   * 作为测试函数需要的多项式的内部支持点的值。外侧向量以正交点为索引，内部以测试函数为索引。测试函数空间是PolynomialsP<dim>。
   *
   */
  std::vector<std::vector<double>> test_values_cell;

  /**
   * 初始化置换模式和符号变化模式。
   * @note
   * 这个函数还没有完全填充正确的实现。它需要在未来的版本中统一实现，以便在包含有翻转面的单元格的网格上工作。
   *
   */
  void
  initialize_quad_dof_index_permutation_and_sign_change();
};

DEAL_II_NAMESPACE_CLOSE

#endif


