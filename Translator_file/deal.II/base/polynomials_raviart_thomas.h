//include/deal.II-translator/base/polynomials_raviart_thomas_0.txt
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

#ifndef dealii_polynomials_raviart_thomas_h
#define dealii_polynomials_raviart_thomas_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/point.h>
#include <deal.II/base/polynomial.h>
#include <deal.II/base/polynomial_space.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/tensor_polynomials_base.h>
#include <deal.II/base/tensor_product_polynomials.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN

/**
 * 这个类实现了Brezzi和Fortin书中描述的<i>H<sup>div</sup></i>符合要求的、矢量值的Raviart-Thomas多项式。
 * Raviart-Thomas多项式的构造是这样的：发散是在张量积多项式空间<i>Q<sub>k</sub></i>。因此，每个分量的多项式阶数必须在相应的方向上高一阶，产生二维和三维的多项式空间<i>(Q<sub>k+1,k</sub>,
 * Q<sub>k,k+1</sub>)</i>和<i>(Q<sub>k+1,k,k</sub>, Q<sub>k,k+1,k</sub>,
 * Q<sub>k,k,k+1</sub>)</i>，分别。
 *
 *
 * @ingroup Polynomials
 *
 *
 */
template <int dim>
class PolynomialsRaviartThomas : public TensorPolynomialsBase<dim>
{
public:
  /**
   * 构造函数。创建给定度数的Raviart-Thomas多项式的所有基函数。
   * @arg
   * k：Raviart-Thomas-空间的度数，也就是最大张量积多项式空间的度数
   * <i>Q<sub>k</sub></i> 包含。
   *
   */
  PolynomialsRaviartThomas(const unsigned int k);

  /**
   * 计算每个Raviart-Thomas多项式在 @p unit_point.
   * 的值和一、二次导数
   * 向量的大小必须为零或等于<tt>n()</tt>。
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
   * 返回空间的名称，即<tt>RaviartThomas</tt>。
   *
   */
  std::string
  name() const override;

  /**
   * 返回空间<tt>RT(degree)</tt>中的多项式的数目，而不需要建立PolynomialsRaviartThomas的对象。这是由FiniteElement类所要求的。
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
   * 一个代表单个组件的多项式空间的对象。我们可以通过旋转评估点的坐标来重新使用它。
   *
   */
  const AnisotropicPolynomials<dim> polynomial_space;

  /**
   * 一个静态成员函数，用于创建我们用来初始化#polynomial_space成员变量的多项式空间。
   *
   */
  static std::vector<std::vector<Polynomials::Polynomial<double>>>
  create_polynomials(const unsigned int k);
};


template <int dim>
inline std::string
PolynomialsRaviartThomas<dim>::name() const
{
  return "RaviartThomas";
}


DEAL_II_NAMESPACE_CLOSE

#endif


