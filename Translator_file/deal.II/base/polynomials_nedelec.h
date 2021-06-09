//include/deal.II-translator/base/polynomials_nedelec_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2020 by the deal.II authors
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


#ifndef dealii_polynomials_nedelec_h
#define dealii_polynomials_nedelec_h


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
 * 这个类实现了第一个家族<i>H<sup>curl</sup></i>的符合性，向量值的多项式，由J.-C.
 * Nédélec在1980年提出。Nédélec在1980年提出（Numer. Math.
 * 35）。
 * Nédélec多项式的构造是这样的：卷曲是在张量积多项式空间<i>Q<sub>k</sub></i>中。因此，每个分量的多项式阶数必须在相应的两个方向上高一阶，产生二维和三维的多项式空间<i>(Q<sub>k,k+1</sub>,
 * Q<sub>k+1,k</sub>)</i>和<i>(Q<sub>k,k+1,k+1</sub>, Q<sub>k+1,k,k+1</sub>,
 * Q<sub>k+1,k+1,k</sub>)</i>，分别。
 *
 *
 * @ingroup Polynomials
 *
 *
 */
template <int dim>
class PolynomialsNedelec : public TensorPolynomialsBase<dim>
{
public:
  /**
   * 构造函数。创建所有给定度数的Nédélec多项式的基函数。
   * @arg
   * k：Nédélec空间的度数，它是最大的张量积多项式空间的度数
   * <i>Q<sub>k</sub></i> 包含的。
   *
   */
  PolynomialsNedelec(const unsigned int k);

  /**
   * 计算每个Nédélec多项式在 @p unit_point.
   * 处的值和一、二阶导数
   * 向量的大小必须是零或等于<tt>n()</tt>。
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
   * 返回空间的名称，即<tt>Nedelec</tt>。
   *
   */
  std::string
  name() const override;

  /**
   * 返回空间<tt>N(degree)</tt>中的多项式的数量，而不需要建立PolynomialsNedelec的对象。这是由FiniteElement类所要求的。
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
PolynomialsNedelec<dim>::name() const
{
  return "Nedelec";
}


DEAL_II_NAMESPACE_CLOSE

#endif


