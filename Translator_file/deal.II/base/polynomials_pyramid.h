//include/deal.II-translator/base/polynomials_pyramid_0.txt
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


#ifndef dealii_polynomials_pyramid_h
#define dealii_polynomials_pyramid_h

#include <deal.II/base/config.h>

#include <deal.II/base/scalar_polynomials_base.h>

DEAL_II_NAMESPACE_OPEN

/**
 * 一个命名空间，用于提供对单数参考单元实体（即三角形和四面体）的支持的函数和类。
 *
 * @ingroup simplex
 *
 *
 */
/**
 * 定义在金字塔实体上的多项式。该类是FE_PyramidP的基础。
 *
 *
 */
template <int dim>
class ScalarLagrangePolynomialPyramid : public ScalarPolynomialsBase<dim>
{
public:
  /**
   * 使得维度对外界可用。
   *
   */
  static const unsigned int dimension = dim;

  /* 构造函数将多项式 @p degree 作为输入。   
*  @note  目前，只实现了线性多项式（度数=1）。 
*
*/
  ScalarLagrangePolynomialPyramid(const unsigned int degree);

  /**
   * @copydoc   ScalarPolynomialsBase::evaluate() 。
   * @note  目前，只有向量 @p values 和 @p grads 被填充。
   *
   */
  void
  evaluate(const Point<dim> &           unit_point,
           std::vector<double> &        values,
           std::vector<Tensor<1, dim>> &grads,
           std::vector<Tensor<2, dim>> &grad_grads,
           std::vector<Tensor<3, dim>> &third_derivatives,
           std::vector<Tensor<4, dim>> &fourth_derivatives) const override;

  double
  compute_value(const unsigned int i, const Point<dim> &p) const override;

  /**
   * @copydoc   ScalarPolynomialsBase::compute_derivative() 。
   * @note  目前，只对一阶导数实现。
   *
   */
  template <int order>
  Tensor<order, dim>
  compute_derivative(const unsigned int i, const Point<dim> &p) const;

  Tensor<1, dim>
  compute_1st_derivative(const unsigned int i,
                         const Point<dim> & p) const override;

  Tensor<2, dim>
  compute_2nd_derivative(const unsigned int i,
                         const Point<dim> & p) const override;

  /**
   * @copydoc   ScalarPolynomialsBase::compute_3rd_derivative()
   * @note  还没有实现。
   *
   */
  Tensor<3, dim>
  compute_3rd_derivative(const unsigned int i,
                         const Point<dim> & p) const override;

  /**
   * @copydoc   ScalarPolynomialsBase::compute_4th_derivative() .
   * @note  还没有实施。
   *
   */
  Tensor<4, dim>
  compute_4th_derivative(const unsigned int i,
                         const Point<dim> & p) const override;

  /**
   * @copydoc   ScalarPolynomialsBase::compute_grad()
   * ScalarPolynomialsBase::compute_grad() 。
   * @note  还没有实施。
   *
   */
  Tensor<1, dim>
  compute_grad(const unsigned int i, const Point<dim> &p) const override;

  /**
   * @copydoc   ScalarPolynomialsBase::compute_grad_grad()
   * ScalarPolynomialsBase::compute_grad_grad() 。
   * @note  还没有实施。
   *
   */
  Tensor<2, dim>
  compute_grad_grad(const unsigned int i, const Point<dim> &p) const override;

  std::string
  name() const override;

  virtual std::unique_ptr<ScalarPolynomialsBase<dim>>
  clone() const override;
};



template <int dim>
template <int order>
Tensor<order, dim>
ScalarLagrangePolynomialPyramid<dim>::compute_derivative(
  const unsigned int i,
  const Point<dim> & p) const
{
  Tensor<order, dim> der;

  Assert(order == 1, ExcNotImplemented());
  const auto grad = compute_grad(i, p);

  for (unsigned int i = 0; i < dim; i++)
    der[i] = grad[i];

  return der;
}

DEAL_II_NAMESPACE_CLOSE

#endif


