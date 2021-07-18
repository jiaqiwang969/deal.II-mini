//include/deal.II-translator/base/function_derivative_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2020 by the deal.II authors
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

#ifndef dealii_function_derivative_h
#define dealii_function_derivative_h

#include <deal.II/base/config.h>

#include <deal.II/base/auto_derivative_function.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/function.h>

DEAL_II_NAMESPACE_OPEN


/**
 * 一个函数对象的派生。
 * 该类的值访问函数返回一个函数相对于构造时提供的方向的导数。如果<tt>b</tt>是向量，则计算出导数<tt>b
 * . grad
 * f</tt>。这个导数是直接计算的，而不是通过计算<tt>f</tt>的梯度和它与<tt>b</tt>的标量乘积。
 * 该导数是通过数值计算的，使用所提供的差分公式之一（见<tt>set_formula</tt>，了解可用方案）。为了获得足够的结果，可能需要对<tt>h</tt>和差分方案进行试验。
 *
 *
 * @ingroup functions
 *
 *
 */
template <int dim>
class FunctionDerivative : public AutoDerivativeFunction<dim>
{
public:
  /**
   * 构造函数。提供了计算导数的函数，微分的方向向量和差分公式的步长<tt>h</tt>。
   *
   */
  FunctionDerivative(const Function<dim> &f,
                     const Point<dim> &   direction,
                     const double         h = 1.e-6);

  /**
   * 构造函数。提供了计算导数的函数和每个正交点的微分方向向量以及差分步长。
   * 这是一个可变速度场的构造函数。最有可能的是，必须为每一组正交点构造一个新的<tt>FunctionDerivative</tt>的对象。
   * 当数值被访问时，正交点的数量必须仍然是相同的。
   *
   */
  FunctionDerivative(const Function<dim> &          f,
                     const std::vector<Point<dim>> &direction,
                     const double                   h = 1.e-6);

  /**
   * 选择差分公式。这在构造函数中被设置为默认值。
   * 现在实现的公式是一阶后向欧拉（<tt>UpwindEuler</tt>），二阶对称欧拉（<tt>Euler</tt>）和一个对称四阶公式（<tt>FourthOrder</tt>）。
   *
   */
  void
  set_formula(typename AutoDerivativeFunction<dim>::DifferenceFormula formula =
                AutoDerivativeFunction<dim>::Euler);
  /**
   * 改变差分公式的基阶大小
   *
   */
  void
  set_h(const double h);

  virtual double
  value(const Point<dim> &p, const unsigned int component = 0) const override;

  virtual void
  vector_value(const Point<dim> &p, Vector<double> &value) const override;

  virtual void
  value_list(const std::vector<Point<dim>> &points,
             std::vector<double> &          values,
             const unsigned int             component = 0) const override;

  /**
   * 返回这个对象的内存消耗估计值，单位是字节。
   * 这不是精确的（但通常会很接近），因为计算树的内存使用量（例如，
   * <tt>std::map</tt>) 很困难。
   *
   */
  virtual std::size_t
  memory_consumption() const override;

private:
  /**
   * 用于微分的函数。
   *
   */
  const Function<dim> &f;

  /**
   * 差分公式的步骤大小。
   *
   */
  double h;

  /**
   * 差分公式。
   *
   */
  typename AutoDerivativeFunction<dim>::DifferenceFormula formula;

  /**
   * 帮助对象。包含公式的增量向量。
   *
   */
  std::vector<Tensor<1, dim>> incr;
};

DEAL_II_NAMESPACE_CLOSE

#endif


