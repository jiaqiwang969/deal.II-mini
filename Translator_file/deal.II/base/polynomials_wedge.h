//include/deal.II-translator/base/polynomials_wedge_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2021 by the deal.II authors
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


#ifndef dealii_base_polynomials_wedge_h
#define dealii_base_polynomials_wedge_h

#include <deal.II/base/config.h>

#include <deal.II/base/polynomials_barycentric.h>
#include <deal.II/base/scalar_polynomials_base.h>

DEAL_II_NAMESPACE_OPEN

/**
 * 定义在楔形实体上的多项式。该类是FE_WedgeP的基础。
 * 多项式是通过一个 BarycentricPolynomials<2>::get_fe_p_basis(degree)
 * 和一个 BarycentricPolynomials<1>::get_fe_p_basis(degree),
 * 的张量乘积创建的，但是，为了更好地匹配FiniteElement的定义，重新进行了列举。
 *
 *
 */
template <int dim>
class ScalarLagrangePolynomialWedge : public ScalarPolynomialsBase<dim>
{
public:
  /**
   * 使得该维度对外部可用。
   *
   */
  static const unsigned int dimension = dim;

  /* 构造函数将多项式 @p degree 作为输入。   
*  @note  目前，只实现了线性（度数=1）和二次多项式（度数=2）。 
*
*/
  ScalarLagrangePolynomialWedge(const unsigned int degree);

  /**
   * @copydoc   ScalarPolynomialsBase::evaluate() .
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

  /**
   * @copydoc   ScalarPolynomialsBase::compute_2nd_derivative()
   * @note  还没有实现。
   *
   */
  Tensor<2, dim>
  compute_2nd_derivative(const unsigned int i,
                         const Point<dim> & p) const override;

  /**
   * @copydoc   ScalarPolynomialsBase::compute_3rd_derivative()
   * @note  还没有实施。
   *
   */
  Tensor<3, dim>
  compute_3rd_derivative(const unsigned int i,
                         const Point<dim> & p) const override;

  /**
   * @copydoc   ScalarPolynomialsBase::compute_4th_derivative()
   * ScalarPolynomialsBase::compute_4th_derivative() 。
   * @note  还没有实施。
   *
   */
  Tensor<4, dim>
  compute_4th_derivative(const unsigned int i,
                         const Point<dim> & p) const override;

  /**
   * @copydoc   ScalarPolynomialsBase::compute_grad() .
   * @note  还没有实施。
   *
   */
  Tensor<1, dim>
  compute_grad(const unsigned int i, const Point<dim> &p) const override;

  /**
   * @copydoc   ScalarPolynomialsBase::compute_grad_grad() .
   * @note  还没有实施。
   *
   */
  Tensor<2, dim>
  compute_grad_grad(const unsigned int i, const Point<dim> &p) const override;

  std::string
  name() const override;

  virtual std::unique_ptr<ScalarPolynomialsBase<dim>>
  clone() const override;

private:
  /**
   * 在一个三角形上定义的标量多项式。
   *
   */
  const BarycentricPolynomials<2> poly_tri;

  /**
   * 在直线上定义的标量多项式。
   *
   */
  const BarycentricPolynomials<1> poly_line;
};



template <int dim>
template <int order>
Tensor<order, dim>
ScalarLagrangePolynomialWedge<dim>::compute_derivative(
  const unsigned int i,
  const Point<dim> & p) const
{
  Tensor<order, dim> der;

  AssertDimension(order, 1);
  const auto grad = compute_grad(i, p);

  for (unsigned int i = 0; i < dim; i++)
    der[i] = grad[i];

  return der;
}

DEAL_II_NAMESPACE_CLOSE

#endif


