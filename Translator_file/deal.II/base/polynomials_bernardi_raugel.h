//include/deal.II-translator/base/polynomials_bernardi_raugel_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2020 by the deal.II authors
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


#ifndef dealii_polynomials_bernardi_raugel_h
#define dealii_polynomials_bernardi_raugel_h

#include <deal.II/base/config.h>

#include <deal.II/base/geometry_info.h>
#include <deal.II/base/point.h>
#include <deal.II/base/polynomial.h>
#include <deal.II/base/polynomial_space.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/tensor_polynomials_base.h>
#include <deal.II/base/tensor_product_polynomials.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN


/**
 * 这个类实现了Bernardi-Raugel多项式，类似于Christine
 * Bernardi和Geneviève Raugel在1985年发表的<i>Mathematics of
 * Computation</i>论文中的描述。
 * Bernardi-Raugel多项式最初被定义为通过增加气泡函数来丰富Stokes问题的简单网格上的
 * $(P_1)^d$ 元素，产生一个无锁定的有限元，它是 $(P_2)^d$
 * 元素的一个子集。这个实现是对 $(Q_1)^d$
 * 元素的丰富，它是 $(Q_2)^d$
 * 元素的一个子集，用于四边形和六面体网格。 $BR_1$
 * 气泡函数被定义为在面 $e_i$ 的中心有量级1和面
 * $\mathbf{n}_i$
 * 的法线方向，而在所有其他顶点和面有量级0。排序与GeometryInfo中的面的编号一致。矢量
 * $\mathbf{n}_i$
 * 指向正轴方向，不一定是元素的法线，以便在各边上保持一致的方向。
 * <dl>  <dt>二维气泡函数（按顺序）  <dd>   $x=0$  边。
 * $\mathbf{p}_1 = \mathbf{n}_1 (1-x)(y)(1-y)$ $x=1$  边缘。  $\mathbf{p}_2
 * = \mathbf{n}_2 (x)(y)(1-y)$ $y=0$  边缘。  $\mathbf{p}_3 = \mathbf{n}_3
 * (x)(1-x)(1-y)$ $y=1$  边缘。  $\mathbf{p}_4 = \mathbf{n}_4 (x)(1-x)(y)$
 * <dt>三维气泡函数（按顺序）  <dd>   $x=0$  边缘。
 * $\mathbf{p}_1 = \mathbf{n}_1 (1-x)(y)(1-y)(z)(1-z)$ $x=1$  边缘。
 * $\mathbf{p}_2 = \mathbf{n}_2 (x)(y)(1-y)(z)(1-z)$ $y=0$  边缘。
 * $\mathbf{p}_3 = \mathbf{n}_3 (x)(1-x)(1-y)(z)(1-z)$ $y=1$  边缘。
 * $\mathbf{p}_4 = \mathbf{n}_4 (x)(1-x)(y)(z)(1-z)$ $z=0$  边缘。
 * $\mathbf{p}_5 = \mathbf{n}_5 (x)(1-x)(y)(1-y)(1-z)$ $z=1$  边缘。
 * $\mathbf{p}_6 = \mathbf{n}_6 (x)(1-x)(y)(1-y)(z)$ </dl> 那么 $BR_1(E)$
 * 多项式在四边形和六面体上的定义是 $BR_1(E) = Q_1(E) \oplus
 * \mbox{span}\{\mathbf{p}_i, i=1,...,2d\}$  。
 *
 *
 *
 * @ingroup Polynomials
 *
 */
template <int dim>
class PolynomialsBernardiRaugel : public TensorPolynomialsBase<dim>
{
public:
  /**
   * 构造函数。创建给定度数的Bernardi-Raugel多项式的所有基函数。
   * @arg  k
   * Bernardi-Raugel-空间的度数，目前只限于<tt>k=1</tt>的情况。
   *
   */
  PolynomialsBernardiRaugel(const unsigned int k);

  /**
   * 返回空间的名称，即<tt>BernardiRaugel</tt>。
   *
   */
  std::string
  name() const override;

  /**
   * 计算每个Bernardi-Raugel多项式在 @p unit_point.
   * 的值和导数，向量的大小必须是零或者等于<tt>n()</tt>。
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
   * 返回空间<tt>BR(degree)</tt>中的多项式的数量，而不需要建立PolynomialsBernardiRaugel的对象。这是由FiniteElement类所要求的。
   *
   */
  static unsigned int
  n_polynomials(const unsigned int k);

  /**
   * @copydoc   TensorPolynomialsBase::clone() .
   *
   */
  virtual std::unique_ptr<TensorPolynomialsBase<dim>>
  clone() const override;

private:
  /**
   * 一个代表Q函数的多项式空间的对象，通过这些函数与相应的单位ijk向量的外积形成<tt>BR</tt>多项式。
   *
   */
  const AnisotropicPolynomials<dim> polynomial_space_Q;

  /**
   * 一个代表泡沫函数的多项式空间的对象，通过这些函数与相应法线的外积形成<tt>BR</tt>多项式。
   *
   */
  const AnisotropicPolynomials<dim> polynomial_space_bubble;

  /**
   * 一个静态成员函数，用于创建我们用来初始化#polynomial_space_Q成员变量的多项式空间。
   *
   */
  static std::vector<std::vector<Polynomials::Polynomial<double>>>
  create_polynomials_Q();

  /**
   * 一个静态成员函数，用于创建我们用来初始化#polynomial_space_bubble成员变量的多项式空间。
   *
   */
  static std::vector<std::vector<Polynomials::Polynomial<double>>>
  create_polynomials_bubble();
};


template <int dim>
inline std::string
PolynomialsBernardiRaugel<dim>::name() const
{
  return "BernardiRaugel";
}


DEAL_II_NAMESPACE_CLOSE

#endif


