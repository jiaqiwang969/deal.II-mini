//include/deal.II-translator/base/polynomials_integrated_legendre_sz_0.txt
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

#ifndef dealii_polynomials_integrated_legendre_sz_h
#define dealii_polynomials_integrated_legendre_sz_h

#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/point.h>
#include <deal.II/base/polynomial.h>


DEAL_II_NAMESPACE_OPEN

/**
 * 实现Sabine
 * Zaglmayr博士论文中描述的集成Legendre多项式的类。
 * 这个类是以现有的deal.II
 * Legendre类为基础编写的，但调整了系数，使递归公式为Sabine
 * Zaglmayr的博士论文中描述的综合Legendre多项式。该多项式可以从以下方面递归生成。
 *
 *
 *
 * -  $L_{0}(x) =
 *
 * -1$ (添加后可以从0开始递归生成)
 *
 *
 *
 * -  $L_{1}(x) = x$
 *
 * -  $L_{2}(x) = \frac{(x^2
 *
 * - 1)}{2}$
 *
 *
 * -  $(n+1)L_{n+1} = (2n-1)L_{n}
 *
 * - (n-2)L_{n-1}$  .
 * 然而，也可以直接从Legendre多项式中生成它们。 $L_{n} =
 * \frac{l_{n}
 *
 * - l_{n-2}}{2n-1)}$
 *
 *
 */
class IntegratedLegendreSZ : public Polynomials::Polynomial<double>
{
public:
  /**
   * 生成p度的多项式系数的构造器。
   *
   */
  IntegratedLegendreSZ(const unsigned int p);

  /**
   * 返回到给定度数的集成Legendre多项式的完整集合。
   *
   */
  static std::vector<Polynomials::Polynomial<double>>
  generate_complete_basis(const unsigned int degree);

private:
  /**
   * 计算度数为p的多项式的协系数的主函数。
   *
   */
  static const std::vector<double>
  get_coefficients(const unsigned int k);
};

DEAL_II_NAMESPACE_CLOSE

#endif


