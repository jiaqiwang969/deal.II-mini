//include/deal.II-translator/base/polynomials_bernstein_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2013 - 2020 by the deal.II authors
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

#ifndef dealii_polynomials_bernstein_h
#define dealii_polynomials_bernstein_h


#include <deal.II/base/config.h>

#include <deal.II/base/polynomial.h>

#include <fstream>
#include <iostream>


DEAL_II_NAMESPACE_OPEN

/**
 * 该类实现了欲望度的伯恩斯坦基础多项式，如http://www.idav.ucdavis.edu/education/CAGDNotes/Bernstein-Polynomials.pdf
 * "从伯恩斯坦基础转换到功率基础 "一段中所述。
 * 它们被用来创建Bernstein有限元FE_Bernstein。
 *
 *
 * @ingroup Polynomials
 *
 */
template <typename number>
class PolynomialsBernstein : public Polynomials::Polynomial<number>
{
public:
  /**
   * 构建 @p index 。
   *
   * - h伯恩斯坦多项式的度数为 @p degree. 。
   * @f{align*}{
   * B_{\text{index}, \text{degree}} (t)
   * &= \text{binom}(\text{degree}, \text{index})
   *    \cdot t^{\text{index}}
   *    \cdot (1
   *
   * - t)^{\text{degree}
   *
   * - \text{index}} \\
   * &= \sum_{i = \text{index}}^\text{degree}
   *    \cdot (-1)^{i
   *
   * - \text{index}}
   *    \cdot \text{binom}(\text{degree}, i)
   *    \cdot \text{binom}(i, \text{index})
   *    \cdot t^i
   * @f}
   * @param  指数  @param  度数
   *
   */
  PolynomialsBernstein(const unsigned int index, const unsigned int degree);
};


template <typename number>
std::vector<Polynomials::Polynomial<number>>
generate_complete_bernstein_basis(const unsigned int degree)
{
  std::vector<Polynomials::Polynomial<number>> v;
  for (unsigned int i = 0; i < degree + 1; ++i)
    v.push_back(PolynomialsBernstein<number>(i, degree));
  return v;
}

DEAL_II_NAMESPACE_CLOSE

#endif


