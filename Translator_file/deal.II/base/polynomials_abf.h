//include/deal.II-translator/base/polynomials_abf_0.txt
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

#ifndef dealii_polynomials_abf_h
#define dealii_polynomials_abf_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/point.h>
#include <deal.II/base/polynomial.h>
#include <deal.II/base/polynomial_space.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/tensor_polynomials_base.h>
#include <deal.II/base/tensor_product_polynomials.h>
#include <deal.II/base/thread_management.h>

#include <memory>
#include <vector>

DEAL_II_NAMESPACE_OPEN

/**
 * 该类实现了Arnold-Boffi-Falk文章中描述的<i>H<sup>div</sup></i>符合要求的、矢量值的Arnold-Boffi-Falk多项式：四边形H（div）有限元，SIAM
 * J. Numer. 分析。Vol.42, No.6, pp.2429-2451
 *
 * ABF多项式的构造是，发散是在张量积多项式空间<i>Q<sub>k</sub></i>。因此，每个分量的多项式阶数必须在相应的方向上高两阶，产生二维和三维的多项式空间<i>(Q<sub>k+2,k</sub>,
 * Q<sub>k,k+2</sub>)</i>和<i>(Q<sub>k+2,k,k</sub>, Q<sub>k,k+2,k</sub>,
 * Q<sub>k,k,k+2</sub>)</i>。
 *
 *
 * @ingroup Polynomials
 *
 *
 */
template <int dim>
class PolynomialsABF : public TensorPolynomialsBase<dim>
{
public:
  /**
   * 构造函数。创建给定度数的Raviart-Thomas多项式的所有基函数。
   * @arg
   * k：Raviart-Thomas-空间的度数，它是最大的张量积多项式空间的度数
   * <i>Q<sub>k</sub></i> 所包含的。
   *
   */
  PolynomialsABF(const unsigned int k);

  /**
   * 计算每个Raviart-Thomas多项式在 @p unit_point.
   * 处的值和一、二次导数
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
   * 返回空间的名称，即<tt>ABF</tt>。
   *
   */
  std::string
  name() const override;

  /**
   * 返回空间<tt>RT(degree)</tt>中的多项式的数目，而不需要建立PolynomialsABF的对象。这是由FiniteElement类所要求的。
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
   * 一个代表单个组件的多项式空间的对象。我们可以通过旋转评估点的坐标，将其重新用于其他矢量分量。
   *
   */
  const AnisotropicPolynomials<dim> polynomial_space;

  /**
   * 一个突发事件，用于守护以下的抓取数组。
   *
   */
  mutable Threads::Mutex mutex;

  /**
   * 辅助内存。
   *
   */
  mutable std::vector<double> p_values;

  /**
   * 辅助内存。
   *
   */
  mutable std::vector<Tensor<1, dim>> p_grads;

  /**
   * 辅助存储器。
   *
   */
  mutable std::vector<Tensor<2, dim>> p_grad_grads;

  /**
   * 辅助存储器。
   *
   */
  mutable std::vector<Tensor<3, dim>> p_third_derivatives;

  /**
   * 辅助存储器。
   *
   */
  mutable std::vector<Tensor<4, dim>> p_fourth_derivatives;
};


template <int dim>
inline std::string
PolynomialsABF<dim>::name() const
{
  return "ABF";
}


DEAL_II_NAMESPACE_CLOSE

#endif


