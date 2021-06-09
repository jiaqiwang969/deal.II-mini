//include/deal.II-translator/base/polynomials_adini_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2020 by the deal.II authors
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


#ifndef dealii_polynomials_adini_h
#define dealii_polynomials_adini_h

#include <deal.II/base/config.h>

#include <deal.II/base/point.h>
#include <deal.II/base/scalar_polynomials_base.h>
#include <deal.II/base/table.h>
#include <deal.II/base/tensor.h>

DEAL_II_NAMESPACE_OPEN

/**
 * 阿迪尼元素的立方多项式空间
 * 这个空间包括由函数<i>xy<sup>3</sup></i>和<i>x<sup>3</sup>y</i>增强的立方空间<i>P<sub>3</sub></i>。
 * 该空间的基础被选择为与阿迪尼元素的节点函数相匹配。
 * @todo
 * 这个多项式空间只在二维中实现，不计算3阶以上的导数。
 *
 *
 * @ingroup Polynomials
 *
 *
 */
template <int dim>
class PolynomialsAdini : public ScalarPolynomialsBase<dim>
{
public:
  /**
   * 描述空间的多项式的构造函数
   *
   */
  PolynomialsAdini();

  /**
   * 计算每个多项式在<tt>unit_point</tt>的值和一、二次导数。
   * 向量的大小必须等于0或等于n()。在第一种情况下，函数不会计算这些值，也就是说，你要通过调整那些你想要填充的向量的大小来表明你想要计算什么。
   * 如果你需要所有多项式的值或导数，那么使用这个函数，而不是使用任何一个compute_value(),
   * compute_grad()或compute_grad_grad()函数，见下文，在所有多项式上循环。
   *
   */
  void
  evaluate(const Point<dim> &           unit_point,
           std::vector<double> &        values,
           std::vector<Tensor<1, dim>> &grads,
           std::vector<Tensor<2, dim>> &grad_grads,
           std::vector<Tensor<3, dim>> &third_derivatives,
           std::vector<Tensor<4, dim>> &fourth_derivatives) const override;

  /**
   * 计算<tt>i</tt>第1个多项式在<tt>unit_point</tt>的值。
   * 可以考虑用evaluate()代替。
   *
   */
  double
  compute_value(const unsigned int i, const Point<dim> &p) const override;

  /**
   * @copydoc   ScalarPolynomialsBase::compute_1st_derivative() .
   *
   */
  virtual Tensor<1, dim>
  compute_1st_derivative(const unsigned int i,
                         const Point<dim> & p) const override;

  /**
   * @copydoc   ScalarPolynomialsBase::compute_2nd_derivative()
   *
   */
  virtual Tensor<2, dim>
  compute_2nd_derivative(const unsigned int i,
                         const Point<dim> & p) const override;

  /**
   * @copydoc   ScalarPolynomialsBase::compute_3rd_derivative() .
   *
   */
  virtual Tensor<3, dim>
  compute_3rd_derivative(const unsigned int i,
                         const Point<dim> & p) const override;

  /**
   * @copydoc   ScalarPolynomialsBase::compute_4th_derivative()
   * ScalarPolynomialsBase::compute_4th_derivative() 。
   *
   */
  virtual Tensor<4, dim>
  compute_4th_derivative(const unsigned int i,
                         const Point<dim> & p) const override;

  /**
   * 计算<tt>i</tt>次多项式在<tt>unit_point</tt>的梯度。
   * 可以考虑用evaluate()代替。
   *
   */
  Tensor<1, dim>
  compute_grad(const unsigned int i, const Point<dim> &p) const override;

  /**
   * 计算<tt>i</tt>次多项式在<tt>unit_point</tt>的二阶导数(grad_grad)。
   * 可以考虑用evaluate()代替。
   *
   */
  Tensor<2, dim>
  compute_grad_grad(const unsigned int i, const Point<dim> &p) const override;

  /**
   * 返回空间的名称，即<tt>PolynomialsAdini</tt>。
   *
   */
  std::string
  name() const override;

  /**
   * @copydoc   ScalarPolynomialsBase::clone() .
   *
   */
  virtual std::unique_ptr<ScalarPolynomialsBase<dim>>
  clone() const override;

private:
  /**
   * 按照 $1,x,y,x^2,y^2,xy,x^3,y^3,xy^2,x^2y,x^3y,xy^3$
   * 的顺序存储多项式的系数。
   *
   */
  Table<2, double> coef;

  /**
   * 按 $1,x,y,x^2,y^2,xy,x^3,y^3,xy^2,x^2y,x^3y,xy^3$
   * 的顺序存储多项式的x导数的系数。
   *
   */
  Table<2, double> dx;

  /**
   * 按 $1,x,y,x^2,y^2,xy,x^3,y^3,xy^2,x^2y,x^3y,xy^3$
   * 的顺序存储多项式的y阶导数的系数。
   *
   */
  Table<2, double> dy;

  /**
   * 按 $1,x,y,x^2,y^2,xy,x^3,y^3,xy^2,x^2y,x^3y,xy^3$
   * 的顺序存储多项式的第二个x导数的系数。
   *
   */
  Table<2, double> dxx;

  /**
   * 按 $1,x,y,x^2,y^2,xy,x^3,y^3,xy^2,x^2y,x^3y,xy^3$
   * 的顺序存储多项式的第二个y导数的系数。
   *
   */
  Table<2, double> dyy;

  /**
   * 按 $1,x,y,x^2,y^2,xy,x^3,y^3,xy^2,x^2y,x^3y,xy^3$
   * 的顺序存储多项式的第二个混合导数的系数。
   *
   */
  Table<2, double> dxy;
};



template <int dim>
inline Tensor<1, dim>
PolynomialsAdini<dim>::compute_1st_derivative(const unsigned int  /*i*/ ,
                                              const Point<dim> &  /*p*/ ) const
{
  Assert(false, ExcNotImplemented());
  return {};
}



template <int dim>
inline Tensor<2, dim>
PolynomialsAdini<dim>::compute_2nd_derivative(const unsigned int  /*i*/ ,
                                              const Point<dim> &  /*p*/ ) const
{
  Assert(false, ExcNotImplemented());
  return {};
}



template <int dim>
inline Tensor<3, dim>
PolynomialsAdini<dim>::compute_3rd_derivative(const unsigned int  /*i*/ ,
                                              const Point<dim> &  /*p*/ ) const
{
  Assert(false, ExcNotImplemented());
  return {};
}



template <int dim>
inline Tensor<4, dim>
PolynomialsAdini<dim>::compute_4th_derivative(const unsigned int  /*i*/ ,
                                              const Point<dim> &  /*p*/ ) const
{
  Assert(false, ExcNotImplemented());
  return {};
}



template <int dim>
inline std::string
PolynomialsAdini<dim>::name() const
{
  return "PolynomialsAdini";
}



DEAL_II_NAMESPACE_CLOSE

#endif


