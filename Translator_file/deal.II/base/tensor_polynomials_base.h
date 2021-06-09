//include/deal.II-translator/base/tensor_polynomials_base_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2020 by the deal.II authors
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

#ifndef dealii_tensor_polynomials_base_h
#define dealii_tensor_polynomials_base_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN

/**
 * 该类为有限元多项式类提供了一个框架，用于从FE_PolyTensor派生的有限元类。这个类型的对象（或者说是从这个类派生出来的类型）被作为成员变量存储在每个FE_PolyTensor类型的对象中。
 * <h3>Deriving classes</h3>
 * 任何派生类必须为在参考单元上评估的形状函数提供最基本的属性。这包括但不限于实现evaluation()、name()和clone()成员函数。这些函数对于存储派生类中的多项式如何在参考单元上的给定点进行评估的最基本信息是必要的。关于每个函数的更多信息可以在相应函数的文档中找到。
 * 从这个类派生的一些类包括  <ul>   <li>  <tt>PolynomialsABF</tt>  <li>  <tt>PolynomialsBDM</tt>  <li>  ] <tt>PolynomialsBernardiRaugel</tt>  <li>  <tt>PolynomialsNedelec</tt>  <li>  <tt>PolynomialsRaviartThomas</tt>  <li>  <tt>PolynomialsRT_Bubbles</tt>  </ul>
 * @ingroup Polynomials
 *
 */
template <int dim>
class TensorPolynomialsBase
{
public:
  /**
   * 构造函数。这需要空间的度数， @p deg 来自有限元类，
   * @p n, 是空间的多项式数目。
   *
   */
  TensorPolynomialsBase(const unsigned int deg,
                        const unsigned int n_polynomials);

  /**
   * 移动构造器。
   *
   */
  TensorPolynomialsBase(TensorPolynomialsBase<dim> &&) = default; // NOLINT

  /**
   * 复制构造函数。
   *
   */
  TensorPolynomialsBase(const TensorPolynomialsBase<dim> &) = default;

  /**
   * 虚拟解构器。确保这个类的指针被正确删除。
   *
   */
  virtual ~TensorPolynomialsBase() = default;

  /**
   * 计算 @p unit_point.
   * 的多项式的值和导数，向量的大小必须是零或等于<tt>n()</tt>。
   * 在第一种情况下，函数将不计算这些值。
   * 如果你需要所有多项式的值或导数，那么使用这个函数，而不是使用任何<tt>compute_value</tt>,
   * <tt>compute_grad</tt> 或 <tt>compute_grad_grad</tt>
   * 函数，见下文，在所有张量积多项式上循环。
   *
   */
  virtual void
  evaluate(const Point<dim> &           unit_point,
           std::vector<Tensor<1, dim>> &values,
           std::vector<Tensor<2, dim>> &grads,
           std::vector<Tensor<3, dim>> &grad_grads,
           std::vector<Tensor<4, dim>> &third_derivatives,
           std::vector<Tensor<5, dim>> &fourth_derivatives) const = 0;

  /**
   * 返回多项式的数量。
   *
   */
  unsigned int
  n() const;

  /**
   * 返回这个类所代表的多项式的最高阶数。如果派生类的值与
   * @p my_degree. 不同，可以重写这个。
   *
   */
  unsigned int
  degree() const;

  /**
   * 一个虚拟的拷贝构造函数，这个函数返回多项式空间对象的一个拷贝。派生类需要在这个基类中覆盖这里的函数，并返回一个与派生类相同类型的对象。
   * 库中的一些地方，例如FE_PolyTensor的构造函数，需要在不知道其确切类型的情况下制作多项式空间的副本。
   * 他们通过这个函数来实现。
   *
   */
  virtual std::unique_ptr<TensorPolynomialsBase<dim>>
  clone() const = 0;

  /**
   * 返回空间的名称。
   *
   */
  virtual std::string
  name() const = 0;

private:
  /**
   * 这个对象所代表的这个函数的最高多项式程度。
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
TensorPolynomialsBase<dim>::n() const
{
  return n_pols;
}



template <int dim>
inline unsigned int
TensorPolynomialsBase<dim>::degree() const
{
  return polynomial_degree;
}



DEAL_II_NAMESPACE_CLOSE

#endif


