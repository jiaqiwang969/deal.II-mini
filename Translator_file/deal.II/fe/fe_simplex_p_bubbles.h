//include/deal.II-translator/fe/fe_simplex_p_bubbles_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2020 - 2021 by the deal.II authors
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

#ifndef dealii_fe_fe_p_bubbles_h
#define dealii_fe_fe_p_bubbles_h

#include <deal.II/base/config.h>

#include <deal.II/base/polynomials_barycentric.h>

#include <deal.II/fe/fe_poly.h>

DEAL_II_NAMESPACE_OPEN

/**
 * @brief Enriched version of FE_P that can be used with nodal quadrature.
 * 许多显式时间积分方案需要在每个时间步长解出一个质量矩阵。有各种方法可以绕过这个要求
 *
 * - 例如， step-48 用一个对角线近似来代替质量矩阵，这使得求解步骤变得微不足道。在 step-48 中，也通常用于张量积元素，这是通过用基于有限元节点的低阶正交点（即通过使用形状函数作为插值基础得到的节点正交规则）来计算质量矩阵。
 * 标准的基于单线的有限元的一个主要缺点是它们不能用于节点正交，因为一些正交权重最终要么为零，要么为负，导致质量矩阵的近似无法解决或不稳定。例如：FE_P<2>(2)的形状函数在顶点的支持点的均值为零，所以该元素不能用于质量包络。
 * 这个元素避免了这个问题，它用一个可修正的节点正交规则构造的增强空间取代了FE_P的形状函数。例如，在三角形上增加了一个对应于中心点插值的基础函数（所有其他的基础函数都被更新以保持统一分割的特性）。这就产生了具有正数的形状函数（即有效的节点正交公式）。同样，在三维空间中，FE_P<3>(2)的多项式空间被充实了五个额外的自由度（其中四个在面心有支持点，一个在中心点有支持点），以便能够构建有效的节点正交规则。
 * 由于这个FE空间包括气泡（即只在元素内部非零的额外函数），所以组件基函数的多项式度数高于元素的实际近似度数。例如，用三维的构造参数
 * <code>degree = 2</code>
 * ，多项式实际上是立方的（3度），但逼近的顺序与我们使用二次（2度）有限元的情况相同。
 * 二维二次元最早是在  @cite fried1975finite
 * 中描述的。这里实现的三维二次元在  @cite Geevers_2018
 * 中首次描述。更高程度的元素可以修改为块状，但尚未在本类中实现。
 *
 *
 */
template <int dim, int spacedim = dim>
class FE_SimplexP_Bubbles : public dealii::FE_Poly<dim, spacedim>
{
public:
  /**
   * 构造函数，接受近似度作为参数。多项式空间通常比这个元素的近似空间高一个度数：更多信息请参见这个类的一般文档。
   * @note  对于 <code>degree == 1</code> 这个元素相当于FE_P(1)。
   *
   */
  FE_SimplexP_Bubbles(const unsigned int degree);

  /**
   * @copydoc   dealii::FiniteElement::clone() .
   *
   */
  virtual std::unique_ptr<FiniteElement<dim, spacedim>>
  clone() const override;

  /**
   * 返回一个唯一标识有限元的字符串。该类返回<tt>FE_SimplexP_Bubbles<dim,spacedim>(degree)</tt>，其中
   * @p dim,   @p spacedim, 和 @p degree
   * 由适当的值代替。像往常一样， @p spacedim
   * 在维度为零的情况下被省略了。
   *
   */
  virtual std::string
  get_name() const override;

  /**
   * @copydoc
   * dealii::FiniteElement::convert_generalized_support_point_values_to_dof_values()
   * .
   *
   */
  virtual void
  convert_generalized_support_point_values_to_dof_values(
    const std::vector<Vector<double>> &support_point_values,
    std::vector<double> &              nodal_values) const override;

protected:
  /**
   * 近似的程度（即构造函数参数）。
   *
   */
  unsigned int approximation_degree;
};

DEAL_II_NAMESPACE_CLOSE

#endif


