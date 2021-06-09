//include/deal.II-translator/base/polynomials_rannacher_turek_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2020 by the deal.II authors
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


#ifndef dealii_polynomials_rannacher_turek_h
#define dealii_polynomials_rannacher_turek_h

#include <deal.II/base/config.h>

#include <deal.II/base/point.h>
#include <deal.II/base/scalar_polynomials_base.h>
#include <deal.II/base/tensor.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN


/**
 * 用于最低阶Rannacher
 * Turek元素的单位平方上的多项式空间的基。
 * 第i个基函数是对应于dof的双基元素，它评估了该函数在第i个面上的平均值。编号可以在GeometryInfo中找到。
 *
 *
 * @ingroup Polynomials
 *
 *
 */
template <int dim>
class PolynomialsRannacherTurek : public ScalarPolynomialsBase<dim>
{
public:
  /**
   * 我们所处的维度。
   *
   */
  static const unsigned int dimension = dim;

  /**
   * 构造函数，检查基础是否在这个维度上实现。
   *
   */
  PolynomialsRannacherTurek();

  /**
   * 基函数 @p i 在 @p p. 的值。
   *
   */
  double
  compute_value(const unsigned int i, const Point<dim> &p) const override;

  /**
   * 在 @p p. 处的基函数 @p i 的<tt>阶</tt>-th
   * 考虑使用evaluate()代替。
   *
   */
  template <int order>
  Tensor<order, dim>
  compute_derivative(const unsigned int i, const Point<dim> &p) const;

  /**
   * @copydoc   ScalarPolynomialsBase::compute_1st_derivative() 。
   *
   */
  virtual Tensor<1, dim>
  compute_1st_derivative(const unsigned int i,
                         const Point<dim> & p) const override;

  /**
   * @copydoc   ScalarPolynomialsBase::compute_2nd_derivative()
   * ScalarPolynomialsBase::compute_2nd_derivative() .
   *
   */
  virtual Tensor<2, dim>
  compute_2nd_derivative(const unsigned int i,
                         const Point<dim> & p) const override;

  /**
   * @copydoc   ScalarPolynomialsBase::compute_3rd_derivative()
   * ScalarPolynomialsBase::compute_3rd_derivative() .
   *
   */
  virtual Tensor<3, dim>
  compute_3rd_derivative(const unsigned int i,
                         const Point<dim> & p) const override;

  /**
   * @copydoc   ScalarPolynomialsBase::compute_4th_derivative()
   * ScalarPolynomialsBase::compute_4th_derivative() .
   *
   */
  virtual Tensor<4, dim>
  compute_4th_derivative(const unsigned int i,
                         const Point<dim> & p) const override;

  /**
   * 基础函数 @p i 在 @p p. 的梯度
   *
   */
  Tensor<1, dim>
  compute_grad(const unsigned int i, const Point<dim> &p) const override;

  /**
   * 基函数 @p i 在 @p p. 处的梯度。
   *
   */
  Tensor<2, dim>
  compute_grad_grad(const unsigned int i, const Point<dim> &p) const override;

  /**
   * 计算 @p unit_point. 处所有基函数的值和导数
   * 矢量的大小必须等于多项式的数量或者是零。大小为零意味着我们不计算向量的条目。
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
   * 返回空间的名称，即<tt>RannacherTurek</tt>。
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
};


namespace internal
{
  namespace PolynomialsRannacherTurekImplementation
  {
    template <int order, int dim>
    inline Tensor<order, dim>
    compute_derivative(const unsigned int, const Point<dim> &)
    {
      Assert(dim == 2, ExcNotImplemented());
      return Tensor<order, dim>();
    }


    template <int order>
    inline Tensor<order, 2>
    compute_derivative(const unsigned int i, const Point<2> &p)
    {
      const unsigned int dim = 2;

      Tensor<order, dim> derivative;
      switch (order)
        {
          case 1:
            {
              Tensor<1, dim> &grad =
                *reinterpret_cast<Tensor<1, dim> *>(&derivative);
              if (i == 0)
                {
                  grad[0] = -2.5 + 3 * p(0);
                  grad[1] = 1.5 - 3 * p(1);
                }
              else if (i == 1)
                {
                  grad[0] = -0.5 + 3.0 * p(0);
                  grad[1] = 1.5 - 3.0 * p(1);
                }
              else if (i == 2)
                {
                  grad[0] = 1.5 - 3.0 * p(0);
                  grad[1] = -2.5 + 3.0 * p(1);
                }
              else if (i == 3)
                {
                  grad[0] = 1.5 - 3.0 * p(0);
                  grad[1] = -0.5 + 3.0 * p(1);
                }
              else
                {
                  Assert(false, ExcNotImplemented());
                }
              return derivative;
            }
          case 2:
            {
              Tensor<2, dim> &grad_grad =
                *reinterpret_cast<Tensor<2, dim> *>(&derivative);
              if (i == 0)
                {
                  grad_grad[0][0] = 3;
                  grad_grad[0][1] = 0;
                  grad_grad[1][0] = 0;
                  grad_grad[1][1] = -3;
                }
              else if (i == 1)
                {
                  grad_grad[0][0] = 3;
                  grad_grad[0][1] = 0;
                  grad_grad[1][0] = 0;
                  grad_grad[1][1] = -3;
                }
              else if (i == 2)
                {
                  grad_grad[0][0] = -3;
                  grad_grad[0][1] = 0;
                  grad_grad[1][0] = 0;
                  grad_grad[1][1] = 3;
                }
              else if (i == 3)
                {
                  grad_grad[0][0] = -3;
                  grad_grad[0][1] = 0;
                  grad_grad[1][0] = 0;
                  grad_grad[1][1] = 3;
                }
              return derivative;
            }
          default:
            {
              // higher derivatives are all zero
              return Tensor<order, dim>();
            }
        }
    }
  } // namespace PolynomialsRannacherTurekImplementation
} // namespace internal



// template functions
template <int dim>
template <int order>
Tensor<order, dim>
PolynomialsRannacherTurek<dim>::compute_derivative(const unsigned int i,
                                                   const Point<dim> & p) const
{
  return internal::PolynomialsRannacherTurekImplementation::compute_derivative<
    order>(i, p);
}



template <int dim>
inline Tensor<1, dim>
PolynomialsRannacherTurek<dim>::compute_1st_derivative(
  const unsigned int i,
  const Point<dim> & p) const
{
  return compute_derivative<1>(i, p);
}



template <int dim>
inline Tensor<2, dim>
PolynomialsRannacherTurek<dim>::compute_2nd_derivative(
  const unsigned int i,
  const Point<dim> & p) const
{
  return compute_derivative<2>(i, p);
}



template <int dim>
inline Tensor<3, dim>
PolynomialsRannacherTurek<dim>::compute_3rd_derivative(
  const unsigned int i,
  const Point<dim> & p) const
{
  return compute_derivative<3>(i, p);
}



template <int dim>
inline Tensor<4, dim>
PolynomialsRannacherTurek<dim>::compute_4th_derivative(
  const unsigned int i,
  const Point<dim> & p) const
{
  return compute_derivative<4>(i, p);
}



template <int dim>
inline std::string
PolynomialsRannacherTurek<dim>::name() const
{
  return "RannacherTurek";
}


DEAL_II_NAMESPACE_CLOSE

#endif


