//include/deal.II-translator/base/incremental_function_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2021 by the deal.II authors
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

#ifndef dealii_incremental_function_h
#define dealii_incremental_function_h


#include <deal.II/base/config.h>

#include <deal.II/base/function.h>
#include <deal.II/base/thread_management.h>

#include <deal.II/lac/vector.h>


DEAL_II_NAMESPACE_OPEN

// Forward declaration
#ifndef DOXYGEN
template <typename number>
class Vector;
#endif

namespace Functions
{
  /**
   * 这个类代表一个增量函数。也就是说，给定一个任意的函数
   * <code>f</code>  ，这个类将返回 <code>f(t)
   *
   * - f(t
   *
   * - delta_t)</code>，其中 <code>f(t)</code> 表示在时间 <code>t</code> 评估的函数，同样地，<code>f(t
   *
   * -delta_t）</code>表示在时间<code>t
   *
   * -delta_t</code>。递减 <code>delta_t</code>
   * 是由set_decrement()方法设置的。这个类的主要应用是将一个给定的Dirichlet边界条件函数转换成增量形式，这是一些非线性求解方案的实现所要求的。
   * @ingroup functions
   *
   */
  template <int dim, typename RangeNumberType = double>
  class IncrementalFunction : public Function<dim, RangeNumberType>
  {
  public:
    /**
     * 将模板参数的值作为一个静态成员常量导出。
     * 这在模板编程的背景下有时很有用。
     *
     */
    static const unsigned int dimension = dim;

    /**
     * 用于表示时间的标量值实数类型。
     *
     */
    using time_type = typename Function<dim, RangeNumberType>::time_type;

    /**
     * 封装给定函数的构造函数  @p base.  。
     * @note 该类存储了对 @p base
     * 的非常量引用，并将在评估过程中调用
     * <code>base.set_time()</code> ，以便在任何任意时间评估 @p
     * base 类。    保证 @p base
     * 的时间状态在该类的每个函数评估后返回到其原始设置。
     *
     */
    IncrementalFunction(Function<dim, RangeNumberType> &base);

    /**
     * 返回函数在给定点的值。
     * 除非只有一个分量（即函数是标量的），否则你应该说明你希望被评估的分量。默认情况下，第一个分量的值被计算出来。
     *
     */
    virtual RangeNumberType
    value(const Point<dim> &p, const unsigned int component = 0) const override;

    /**
     * 返回一个矢量值函数在某一点的所有分量。
     * 在调用这个函数之前，要求 @p values
     * 向量有正确的大小。
     *
     */
    virtual void
    vector_value(const Point<dim> &       p,
                 Vector<RangeNumberType> &values) const override;

    /**
     * 设置时间递减。        希望这个值是正的。
     *
     */
    void
    set_decrement(const time_type delta_t);

  private:
    /**
     * 对被包装的函数的引用。
     *
     */
    Function<dim, RangeNumberType> &base;

    /**
     * 时间递减。
     *
     */
    time_type delta_t;

    /**
     * 一个用于存储数值的辅助向量。
     *
     */
    mutable Vector<RangeNumberType> values_old;

    /**
     * 用于支持多线程背景下的评估的线程互斥。
     *
     */
    mutable Threads::Mutex mutex;
  };

} // namespace Functions


DEAL_II_NAMESPACE_CLOSE

#endif


