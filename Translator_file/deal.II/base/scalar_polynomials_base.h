//include/deal.II-translator/base/scalar_polynomials_base_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2021 by the deal.II authors
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

#ifndef dealii_scalar_polynomials_base_h
#define dealii_scalar_polynomials_base_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN

/**
 * 该类为有限元多项式类提供了一个框架，用于从FE_Poly派生的有限元类。这个类型的对象（或者说从这个类派生出来的类型）作为成员变量存储在每个FE_Poly类型的对象中。
 * <h3>Deriving classes</h3>
 * 任何派生类必须为在参考单元上评估的形状函数提供最基本的属性。这包括但不限于实现evaluation()、name()和clone()成员函数。这些函数对于存储派生类中的多项式如何在参考单元上的给定点进行评估的最基本信息是必要的。关于每个函数的更多信息可以在相应函数的文档中找到。
 * 从这个类派生的一些类包括  <ul>   <li>  <tt>PolynomialsAdini</tt>  <li>  <tt>PolynomialsRannacherTurek</tt>  <li>  <tt>PolynomialsP</tt>  <li>  ] <tt>PolynomialSpace</tt>  <li>  <tt>TensorProductPolynomials</tt>  <li>  <tt>TensorProductPolynomialsConst</tt>  <li>  <tt>TensorProductPolynomialsBubbles</tt>  </ul>
 *
 *
 * @ingroup Polynomials
 *
 */
template <int dim>
class ScalarPolynomialsBase
{
public:
  /**
   * 构造函数。这需要空间的度数， @p deg 来自有限元类，
   * @p n, 是空间的多项式数目。
   *
   */
  ScalarPolynomialsBase(const unsigned int deg,
                        const unsigned int n_polynomials);

  /**
   * 移动构造器。
   *
   */
  ScalarPolynomialsBase(ScalarPolynomialsBase<dim> &&) = default; // NOLINT

  /**
   * 复制构造函数。
   *
   */
  ScalarPolynomialsBase(const ScalarPolynomialsBase<dim> &) = default;

  /**
   * 虚拟解构器。确保这个类的指针被正确删除。
   *
   */
  virtual ~ScalarPolynomialsBase() = default;

  /**
   * 计算 @p unit_point.
   * 的多项式的值和导数，向量的大小必须是零或者等于<tt>n()</tt>。
   * 在第一种情况下，函数将不计算这些值。
   * 如果你需要所有多项式的值或导数，那么使用这个函数，而不是使用任何<tt>compute_value</tt>,
   * <tt>compute_grad</tt> 或 <tt>compute_grad_grad</tt>
   * 函数，见下文，在所有张量积多项式上循环。
   *
   */
  virtual void
  evaluate(const Point<dim> &           unit_point,
           std::vector<double> &        values,
           std::vector<Tensor<1, dim>> &grads,
           std::vector<Tensor<2, dim>> &grad_grads,
           std::vector<Tensor<3, dim>> &third_derivatives,
           std::vector<Tensor<4, dim>> &fourth_derivatives) const = 0;

  /**
   * 计算在单位点<tt>p</tt>上的<tt>i</tt>次多项式的值。
   * 可以考虑用evaluate()代替。
   *
   */
  virtual double
  compute_value(const unsigned int i, const Point<dim> &p) const = 0;

  /**
   * 计算<tt>i</tt>次多项式在单位点<tt>p</tt>的<tt>阶</tt>次导数。
   * 可以考虑用evaluate()代替。      @tparam  order 导数的阶数。
   *
   */
  template <int order>
  Tensor<order, dim>
  compute_derivative(const unsigned int i, const Point<dim> &p) const;

  /**
   * 计算<tt>i</tt>第1个多项式在单位点<tt>p</tt>的一阶导数。
   * 可以考虑用evaluate()代替。
   *
   */
  virtual Tensor<1, dim>
  compute_1st_derivative(const unsigned int i, const Point<dim> &p) const = 0;

  /**
   * 计算<tt>i</tt>次多项式在单位点<tt>p</tt>的二阶导数。
   * 可以考虑用evaluate()代替。
   *
   */
  virtual Tensor<2, dim>
  compute_2nd_derivative(const unsigned int i, const Point<dim> &p) const = 0;

  /**
   * 计算<tt>i</tt>次多项式在单位点<tt>p</tt>的三阶导数。
   * 可以考虑用evaluate()代替。
   *
   */
  virtual Tensor<3, dim>
  compute_3rd_derivative(const unsigned int i, const Point<dim> &p) const = 0;

  /**
   * 计算<tt>i</tt>次多项式在单位点<tt>p</tt>的第四导数。
   * 可以考虑用evaluate()代替。
   *
   */
  virtual Tensor<4, dim>
  compute_4th_derivative(const unsigned int i, const Point<dim> &p) const = 0;

  /**
   * 计算<tt>i</tt>次多项式在单位点<tt>p</tt>的梯度。
   * 可以考虑用evaluate()代替。
   *
   */
  virtual Tensor<1, dim>
  compute_grad(const unsigned int  /*i*/, const Point<dim> & /*p*/ ) const = 0;

  /**
   * 计算<tt>i</tt>次多项式在单位点<tt>p</tt>的二阶导数（grad_grad）。
   * 可以考虑用evaluate()代替。
   *
   */
  virtual Tensor<2, dim>
  compute_grad_grad(const unsigned int  /*i*/ ,
                    const Point<dim> &  /*p*/ ) const = 0;

  /**
   * 返回多项式的数量。
   *
   */
  unsigned int
  n() const;

  /**
   * 返回这个类所代表的多项式的最高阶数。如果派生类的值与
   * @p my_degree. 不同，可以重写这个值。
   *
   */
  virtual unsigned int
  degree() const;

  /**
   * 一个虚拟的拷贝构造函数，这个函数返回多项式空间对象的一个拷贝。派生类需要在这个基类中覆盖这里的函数，并返回一个与派生类相同类型的对象。
   * 库中的一些地方，例如FE_Poly的构造函数，需要在不知道其确切类型的情况下制作多项式空间的副本。
   * 他们通过这个函数来实现。
   *
   */
  virtual std::unique_ptr<ScalarPolynomialsBase<dim>>
  clone() const = 0;

  /**
   * 返回空间的名称。
   *
   */
  virtual std::string
  name() const = 0;

  /**
   * 返回这个对象的内存消耗估计值（以字节为单位）。
   *
   */
  virtual std::size_t
  memory_consumption() const;

private:
  /**
   * 此对象所代表的此函数的最高多项式程度。
   *
   */
  const unsigned int polynomial_degree;

  /**
   * 此对象所代表的多项式的数量。
   *
   */
  const unsigned int n_pols;
};



template <int dim>
inline unsigned int
ScalarPolynomialsBase<dim>::n() const
{
  return n_pols;
}



template <int dim>
inline unsigned int
ScalarPolynomialsBase<dim>::degree() const
{
  return polynomial_degree;
}



template <int dim>
template <int order>
inline Tensor<order, dim>
ScalarPolynomialsBase<dim>::compute_derivative(const unsigned int i,
                                               const Point<dim> & p) const
{
  if (order == 1)
    {
      auto derivative = compute_1st_derivative(i, p);
      return *reinterpret_cast<Tensor<order, dim> *>(&derivative);
    }
  if (order == 2)
    {
      auto derivative = compute_2nd_derivative(i, p);
      return *reinterpret_cast<Tensor<order, dim> *>(&derivative);
    }
  if (order == 3)
    {
      auto derivative = compute_3rd_derivative(i, p);
      return *reinterpret_cast<Tensor<order, dim> *>(&derivative);
    }
  if (order == 4)
    {
      auto derivative = compute_4th_derivative(i, p);
      return *reinterpret_cast<Tensor<order, dim> *>(&derivative);
    }
  Assert(false, ExcNotImplemented());
  Tensor<order, dim> empty;
  return empty;
}

DEAL_II_NAMESPACE_CLOSE

#endif


