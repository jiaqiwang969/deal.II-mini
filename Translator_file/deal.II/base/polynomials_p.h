//include/deal.II-translator/base/polynomials_p_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2020 by the deal.II authors
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

#ifndef dealii_polynomials_P_h
#define dealii_polynomials_P_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/point.h>
#include <deal.II/base/polynomial.h>
#include <deal.II/base/polynomial_space.h>
#include <deal.II/base/tensor.h>

#include <memory>
#include <vector>

DEAL_II_NAMESPACE_OPEN
/**
 * @addtogroup Polynomials   
 * 
     * @{ 
 *
 *
 */

/**
 * 这个类实现了基于单项式 ${1,x,x^2,...}$
 * 的<tt>p</tt>度的多项式空间。即在<tt>d</tt>维度上，它构建了所有形式为
 * $\prod_{i=1}^d x_i^{n_i}$ 的多项式，其中 $\sum_i n_i\leq p$
 * 。基数多项式有一个特定的顺序，例如，在2维中。
 * ${1,x,y,xy,x^2,y^2,x^2y,xy^2,x^3,y^3,...}$  .  $P_k1$
 * 中的单项式的排序与 $P_k2$ 中的单项式的排序相匹配，为
 * $k2>k1$  。
 *
 *
 */
template <int dim>
class PolynomialsP : public PolynomialSpace<dim>
{
public:
  /**
   * 访问此对象的维度，用于检查和自动设置其他类中的维度。
   *
   */
  static const unsigned int dimension = dim;

  /**
   * 构造函数。创建所有  $P_p$  的基函数。  
   * @arg
   * p：多项式空间的度数
   *
   */
  PolynomialsP(const unsigned int p);

  /**
   * 返回多项式空间的度数<tt>p</tt> <tt>P_p</tt>。
   * 注意，这个数字是 <tt>PolynomialSpace::degree()-1</tt>,
   * 比较PolynomialSpace中的定义。
   *
   */
  virtual unsigned int
  degree() const override;

  /**
   * 对于<tt>n</tt>次多项式 $p_n(x,y,z)=x^i y^j z^k$
   * 这个函数给出了x,y,z方向上的度数i,j,k。
   * 在1d和2d中，显然只有i和i,j被返回。
   *
   */
  std::array<unsigned int, dim>
  directional_degrees(unsigned int n) const;

  std::unique_ptr<ScalarPolynomialsBase<dim>>
  clone() const override
  {
    return std::make_unique<PolynomialsP<dim>>(*this);
  }

private:
  /**
   * 填充<tt>index_map</tt>。
   *
   */
  void
  create_polynomial_ordering(std::vector<unsigned int> &index_map) const;

  /**
   * 多项式空间的度数<tt>p</tt>  $P_p$
   * ，即给构造函数的数字<tt>p</tt>。
   *
   */
  const unsigned int p;
};

 /** @} */ 

template <int dim>
inline unsigned int
PolynomialsP<dim>::degree() const
{
  return p;
}


template <int dim>
inline std::array<unsigned int, dim>
PolynomialsP<dim>::directional_degrees(unsigned int n) const
{
  return this->compute_index(n);
}

DEAL_II_NAMESPACE_CLOSE

#endif


