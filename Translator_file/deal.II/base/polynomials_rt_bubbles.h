//include/deal.II-translator/base/polynomials_rt_bubbles_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2020 by the deal.II authors
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

#ifndef dealii_polynomials_rt_bubbles_h
#define dealii_polynomials_rt_bubbles_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/point.h>
#include <deal.II/base/polynomial.h>
#include <deal.II/base/polynomials_raviart_thomas.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/tensor_polynomials_base.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN

/**
 * 这个类实现了<i>H<sup>div</sup></i>符合要求的、矢量值的增强型Raviart-Thomas多项式。
 * 与经典的Raviart-Thomas空间类似，增强的Raviart-Thomas多项式的构造是这样的：发散是在张量积多项式空间<i>Q<sub>k-1</sub></i>。
 * 这个空间的形式是<i>V<sub>k</sub> = RT<sub>k-1</sub> +
 * B<sub>k</sub></i>，其中<i>B<sub>k</sub></i>定义如下。  <dl>  <dt>
 * 在二维。
 *
 * @f{align*}{
 * B_k^1(E) = \text{span}\left\{x^{a_1-1} y^{a_2}\begin{pmatrix} (a_2+1) x \\
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * -a_1 y \end{pmatrix}\text{ : } a_2=k \right\} \\
 * B_k^2(E) = \text{span}\left\{x^{b_1} y^{b_2-1}\begin{pmatrix}
 *
 * -b_2 x \\
 *   (b_1+1) y \end{pmatrix}\text{ : } b_1=k \right\}
 * @f}
 *
 * <dt> 在三维中。
 * @f{align*}{
 * B_k^1(E) = \text{span}\left\{x^{a_1-1} y^{a_2} z^{a_3}\begin{pmatrix}
 * (a_2+a_3+2) x \\
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * -a_1 y \\
 *
 * -a_1 z \end{pmatrix}\text{ : } a_2=k \text{ or } a_3=k
 * \right\},\\
 * B_k^2(E) = \text{span}\left\{x^{b_1} y^{b_2-1} z^{b_3}\begin{pmatrix}
 *
 * -b_2
 * x \\
 *   (b_1+b_3+2) y \\
 *
 * -b_2 z \end{pmatrix}\text{ : } b_1=k \text{ or } b_3=k
 * \right\},\\
 * B_k^3(E) = \text{span}\left\{x^{c_1}y^{c_2}z^{c_3-1}\begin{pmatrix}
 *
 * -c_3 x
 * \\
 *
 * -c_3y \\ (c_1+c_2+2)z \end{pmatrix}\text{ : } c_1=k \text{ or } c_2=k
 * \right\},
 * @f}
 * </dl>  其中 $0 \le a_1, a_2, a_3 \le k$  。
 *
 *
 * @note
 * 与经典的Raviart-Thomas空间不同，增强空间的最低阶为1，与Brezzi-Douglas-Marini（BDM）多项式空间类似。
 * 空间的总维度<i>dim(V<sub>k</sub>) =
 * d*(k+1)^d</i>，其中<i>d</i>是空间维度。这允许将形状函数与Gauss-Lobatto正交点联系起来，如下图所示。
 * <table> <tr> <td align="center">
 @image html rtbubbles.png
 * </td></tr>
 *
 * <tr> <td align="center"> Left
 *
 * - $2D,\,k=3$, right
 *
 * - $3D,\,k=2$.</td></tr> </table>
 *
 *
 * @ingroup Polynomials
 *
 */

template <int dim>
class PolynomialsRT_Bubbles : public TensorPolynomialsBase<dim>
{
public:
  /**
   * 构造函数。创建给定度数的RT_bubbles多项式的所有基础函数。
   *
   */
  PolynomialsRT_Bubbles(const unsigned int k);

  /**
   * 计算每个RT_bubbles多项式在 @p unit_point.
   * 的值和一、二次导数，向量的大小必须为零或等于<tt>n()</tt>。
   * 在第一种情况下，该函数将不计算这些值。
   * 如果你需要所有张量积多项式的值或导数，那么使用这个函数，而不是使用任何<tt>compute_value</tt>,
   * <tt>compute_grad</tt>或<tt>compute_grad_grad</tt>函数，见下文，在所有张量积多项式上循环。
   *
   */
  void
  evaluate(const Point<dim> &           unit_point,
           std::vector<Tensor<1, dim>> &values,
           std::vector<Tensor<2, dim>> &grads,
           std::vector<Tensor<3, dim>> &grad_grads,
           std::vector<Tensor<4, dim>> &third_derivatives,
           std::vector<Tensor<5, dim>> &fourth_derivatives) const override;

  /**
   * 返回空间的名称，即<tt>RT_Bubbles</tt>。
   *
   */
  std::string
  name() const override;

  /**
   * 返回空间<tt>RT_Bubbles(degree)</tt>中的多项式的数目，而不需要建立PolynomialsRT-Bubbles的对象。这是由FiniteElement类所要求的。
   *
   */
  static unsigned int
  n_polynomials(const unsigned int degree);

  /**
   * @copydoc   TensorPolynomialsBase::clone() .
   *
   */
  virtual std::unique_ptr<TensorPolynomialsBase<dim>>
  clone() const override;

private:
  /**
   * 一个代表空间的Raviart-Thomas部分的对象
   *
   */
  const PolynomialsRaviartThomas<dim> raviart_thomas_space;

  /**
   * 单项式的存储，我们需要从零度到<i>k+1</i>的所有多项式。
   *
   */
  std::vector<Polynomials::Polynomial<double>> monomials;
};


template <int dim>
inline std::string
PolynomialsRT_Bubbles<dim>::name() const
{
  return "RT_bubbles";
}


DEAL_II_NAMESPACE_CLOSE

#endif


