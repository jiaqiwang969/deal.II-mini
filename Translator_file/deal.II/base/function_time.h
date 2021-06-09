//include/deal.II-translator/base/function_time_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2020 by the deal.II authors
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

#ifndef dealii_function_time_h
#define dealii_function_time_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>

DEAL_II_NAMESPACE_OPEN
/**
 * 支持与时间有关的函数。该库也是为时间相关问题设计的。为此，函数对象也包含一个存储时间的字段，以及操作它们的函数。与时间无关的问题不应该为其他目的而访问甚至滥用它们，但由于人们通常不会创建数以千计的函数对象，通用性的收益权衡了我们不需要为不与时间有关的问题存储时间值的事实。第二个好处是，派生的标准类如<tt>ZeroFunction</tt>,
 * <tt>ConstantFunction</tt>等也适用于时间相关问题。
 * 对时间的访问通过以下函数进行。  <ul>   <li>  <tt>get_time</tt>: 返回时间变量的现值。   <li>  <tt>set_time</tt>：设置时间值为一个特定的值。   <li>  <tt>advance_time</tt>：将时间增加一定的时间步长。  </ul>  后面两个函数是虚拟的，因此派生类可以执行只需要为每个新时间做一次的计算。例如，如果一个依赖时间的函数有一个因子<tt>sin(t)</tt>，那么在set_time()的派生版本中计算这个因子，将其存储在一个成员变量中并使用该变量，而不是在每次调用<tt>value()</tt>、<tt>value_list</tt>或Function类的其他函数时计算它，可能是个合理的选择。
 * 默认情况下，advance_time()函数会用新的时间调用set_time()函数，所以在大多数情况下，只重载set_time()进行计算就足够了，如上图所示。
 * 这个类的构造函数需要一个时间变量的初始值，默认为零。因为给出了一个默认值，所以如果不需要的话，所有的派生类都不需要为时间变量取一个初始值。
 * @tparam  Number
 * 用于存储时间值的数据类型。几乎在所有情况下，这都是默认的
 * @p double,
 * ，但在某些情况下，人们可能希望将时间存储在不同的（而且总是标量）类型中。一个例子是一个可以存储一个值以及它的不确定性的区间类型。另一个例子是一个允许自动微分的类型（例如，见
 * step-33
 * 中使用的Sacado类型），从而可以生成一个函数的分析性（时间性）导数。
 *
 *
 *
 * @ingroup functions
 *
 */
template <typename Number = double>
class FunctionTime
{
public:
  /**
   * 构造函数。可以为时间变量取一个初始值，默认为零。
   *
   */
  FunctionTime(const Number initial_time = Number(0.0));

  /**
   * 虚拟解构器。
   *
   */
  virtual ~FunctionTime() = default;

  /**
   * 返回时间变量的值。
   *
   */
  Number
  get_time() const;

  /**
   * 设置时间为<tt>new_time</tt>，覆盖旧值。
   *
   */
  virtual void
  set_time(const Number new_time);

  /**
   * 将时间按给定的时间步长<tt>delta_t</tt>推进。
   *
   */
  virtual void
  advance_time(const Number delta_t);

  /**
   * 这个类被初始化的类型，它被用来表示时间。
   *
   */
  using time_type = Number;

private:
  /**
   * 存储现在的时间。
   *
   */
  Number time;
};



 /*----------------------------- Inline functions ----------------------------*/ 

#ifndef DOXYGEN

template <typename Number>
inline Number
FunctionTime<Number>::get_time() const
{
  return time;
}

#endif
DEAL_II_NAMESPACE_CLOSE

#endif


