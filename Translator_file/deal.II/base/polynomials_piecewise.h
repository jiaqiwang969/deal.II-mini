//include/deal.II-translator/base/polynomials_piecewise_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2021 by the deal.II authors
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

#ifndef dealii_polynomials_piecewise_h
#define dealii_polynomials_piecewise_h



#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/point.h>
#include <deal.II/base/polynomial.h>
#include <deal.II/base/subscriptor.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN

/**
 * @addtogroup Polynomials   @{
 *
 */

/**
 * 一个命名空间，与描述1d多项式空间有关的类在其中被声明。
 *
 *
 */
namespace Polynomials
{
  /**
   * 单位区间的分片一维多项式的定义。这个空间允许在单位区间的部分上描述插值多项式，类似于在细分元素上定义有限元基函数。这个类的主要目的是允许构造FE_Q_iso_Q1类的形状函数，该类在每个坐标方向上都有一些插值点，但不是用它们来做高阶多项式，而只是选择片状线性形状函数
   *
   * - 实际上，它是一个 $Q_1$ 元素，定义在参考单元的细分上，并在这些子单元中的每一个上进行复制。    这个类不是从ScalarPolynomialsBase基类派生的，因为它实际上不是一个多项式
   *
   * - 它是一个片状多项式。  然而，它与 Polynomials::Polynomial 类是接口兼容的，因此可以作为TensorProductPolynomials的模板参数。
   * @ingroup Polynomials
   *
   */
  template <typename number>
  class PiecewisePolynomial : public Subscriptor
  {
  public:
    /**
     * 拉格朗日多项式的构造函数，该区间是单位区间的一个子集。它使用一个多项式描述，与单位区间相比，子区间的大小、区间的总数（细分）、区间的当前索引以及多项式是否跨越到下一个区间（例如，如果它生活在两个相邻的区间上）。
     * 如果区间数为1，则分片多项式的表现与通常的多项式完全相同。
     *
     */
    PiecewisePolynomial(const Polynomial<number> &coefficients_on_interval,
                        const unsigned int        n_intervals,
                        const unsigned int        interval,
                        const bool                spans_next_interval);

    /**
     * 返回该多项式在给定点的值，评估底层多项式。当超出给定的区间时，该多项式评估为零（当它跨过该区间时，可能是右边的下一个区间）。
     *
     */
    number
    value(const number x) const;

    /**
     * 返回多项式在<tt>x</tt>点的值和导数。 <tt>values[i],
     * i=0,...,values.size()-1</tt>包括<tt>i</tt>的导数。因此，要计算的导数的数量由传递的向量的大小决定。
     * 请注意，所有的导数在单位区间内部的区间边界（假设是精确的算术）评估为零，因为在这种情况下对于片状多项式来说没有唯一的梯度值。这并不总是需要的（例如，当评估元素边界上的梯度跳跃时），但是当没有意义时，用户有责任避免在这些点进行评估。
     *
     */
    void
    value(const number x, std::vector<number> &values) const;

    /**
     * 返回多项式在<tt>x</tt>点的值和导数。 <tt>values[i],
     * i=0,...,n_derivatives</tt>包括<tt>i</tt>的导数。要计算的导数数量由
     * @p n_derivatives 决定， @p values 必须为 @p n_derivatives
     * +1的值提供足够的空间。
     * 请注意，所有导数在单位区间内部的区间边界（假设是精确算术）评估为零，因为在这种情况下，对于片状多项式没有唯一的梯度值。这并不总是需要的（例如，当评估元素边界上的梯度跳跃时），但是当没有意义时，用户有责任避免在这些点进行评估。
     *
     */
    void
    value(const number       x,
          const unsigned int n_derivatives,
          number *           values) const;

    /**
     * 多项式的度数。这是底层基础多项式的度数。
     *
     */
    unsigned int
    degree() const;

    /**
     * 使用[BOOST序列化库](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html)将此对象的数据写入或读出到一个流中，以便进行序列化。
     *
     */
    template <class Archive>
    void
    serialize(Archive &ar, const unsigned int version);

    /**
     * 返回此对象的内存消耗估计值（以字节为单位）。
     *
     */
    virtual std::size_t
    memory_consumption() const;

  protected:
    /**
     * 底层多项式对象，被缩放到一个子区间，并相应地进行连接。
     *
     */
    Polynomial<number> polynomial;

    /**
     * 一个存储单位区间被分成的区间数的变量。
     *
     */
    unsigned int n_intervals;

    /**
     * 一个存储当前多项式在区间范围内的索引的变量。
     *
     */
    unsigned int interval;

    /**
     * 存储如果多项式跨越两个相邻的区间，即子区间中给出的区间和下一个区间。
     *
     */
    bool spans_two_intervals;
  };



  /**
   * 在子区间上给定的度数和区间数的情况下，在单位区间的细分上生成一个完整的拉格朗日基础，并将其划分为更小的区间。
   *
   */
  std::vector<PiecewisePolynomial<double>>
  generate_complete_Lagrange_basis_on_subdivisions(
    const unsigned int n_subdivisions,
    const unsigned int base_degree);

} // namespace Polynomials


 /** @} */ 

 /* -------------------------- inline functions --------------------- */ 

namespace Polynomials
{
  template <typename number>
  inline unsigned int
  PiecewisePolynomial<number>::degree() const
  {
    return polynomial.degree();
  }



  template <typename number>
  inline number
  PiecewisePolynomial<number>::value(const number x) const
  {
    AssertIndexRange(interval, n_intervals);
    number y = x;
    // shift polynomial if necessary
    if (n_intervals > 1)
      {
        const number step = 1. / n_intervals;

        // polynomial spans over two intervals
        if (spans_two_intervals == true)
          {
            const number offset = step * interval;
            if (x < offset)
              return 0;
            else if (x > offset + step + step)
              return 0;
            else if (x < offset + step)
              y = x - offset;
            else
              y = offset + step + step - x;
          }
        else
          {
            const number offset = step * interval;
            if (x < offset || x > offset + step)
              return 0;
            else
              y = x - offset;
          }

        return polynomial.value(y);
      }
    else
      return polynomial.value(x);
  }



  template <typename number>
  template <class Archive>
  inline void
  PiecewisePolynomial<number>::serialize(Archive &ar, const unsigned int)
  {
    // forward to serialization function in the base class.
    ar &static_cast<Subscriptor &>(*this);
    ar &polynomial;
    ar &n_intervals;
    ar &interval;
    ar &spans_two_intervals;
  }

} // namespace Polynomials

DEAL_II_NAMESPACE_CLOSE

#endif


