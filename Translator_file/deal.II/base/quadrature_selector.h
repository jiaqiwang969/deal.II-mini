//include/deal.II-translator/base/quadrature_selector_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2020 by the deal.II authors
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


#ifndef dealii_quadrature_selector_h
#define dealii_quadrature_selector_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/quadrature.h>

#include <string>

DEAL_II_NAMESPACE_OPEN

/**
 * 该类实现了以字符串形式传递给其构造函数的正交规则。支持的正交规则有QGauss（所有阶数）、QMidpoint、QMilne、QSimpson、QTrapezoid和QWeddle。
 * 如果你想使用灵活的正交规则，这个类很有用，它可以从一个参数文件中读取（见ParameterHandler）。
 *
 *
 * @ingroup Quadrature
 *
 *
 */
template <int dim>
class QuadratureSelector : public Quadrature<dim>
{
public:
  /**
   * 构造函数。取正交规则的名称（"高斯"、"米尔尼"、"韦德尔
   * "等之一），如果是
   * "高斯"，则取每个坐标方向上的正交点的数量。
   *
   */
  QuadratureSelector(const std::string &s, const unsigned int order = 0);

  /**
   * 该函数以列表形式返回所有可能的正交点名称，并以<tt>|</tt>分隔，这样你可以用它来定义参数文件（详见ParameterHandler）。
   *
   */
  static std::string
  get_quadrature_names();

  /**
   * @addtogroup  Exceptions  
     * @{ 
   *
   */


  /**
   * 异常情况
   *
   */
  DeclException1(ExcInvalidQGaussOrder,
                 int,
                 << "You tried to generate a QGauss object with an invalid "
                 << "number " << arg1
                 << " of quadrature points in each coordinate "
                 << "direction. This number must be greater than or equal "
                 << "to 1.");
  /**
   * 异常情况
   *
   */
  DeclException2(ExcInvalidOrder,
                 std::string,
                 unsigned int,
                 << "You tried to generate a " << arg1
                 << " object; no order is needed for objects of this kind, but "
                 << arg2 << " was given as argument.");
  /**
   * 异常情况
   *
   */
  DeclException1(ExcInvalidQuadrature,
                 std::string,
                 << arg1 << " is not a valid name for a quadrature rule.");
  //@}
private:
  /**
   * 这个静态函数根据给定的字符串名称和适当的顺序（如果名称是
   * "高斯"）创建一个正交对象。它是由构造函数调用的。
   *
   */
  static Quadrature<dim>
  create_quadrature(const std::string &s, const unsigned int order);
};
DEAL_II_NAMESPACE_CLOSE

#endif


