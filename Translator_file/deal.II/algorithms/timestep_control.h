//include/deal.II-translator/algorithms/timestep_control_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2020 by the deal.II authors
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


#ifndef dealii_time_step_control_h
#define dealii_time_step_control_h

#include <deal.II/base/config.h>

#include <deal.II/base/smartpointer.h>
#include <deal.II/base/subscriptor.h>

#include <deal.II/lac/vector_memory.h>

#include <cstdio>

DEAL_II_NAMESPACE_OPEN

// Forward declaration
#ifndef DOXYGEN
class ParameterHandler;
#endif

namespace Algorithms
{
  /**
   * 时隙方案的控制类。它的主要任务是确定下一个时间步骤的大小和时间间隔中的相应点。  此外，它还控制将解决方案写入文件。    下一个时间步长的确定方法如下。    <ol>   <li>  根据策略，步长暂且加到当前时间上。    <li>  如果得出的时间超过了区间的最终时间，则减少步长，以满足这个时间。    <li>  如果得出的时间低于最终时间，只是步长的一小部分，则增加步长，以满足这个时间。    <li>  产生的步长从当前时间开始使用。    </ol>  变量 @p print_step 可以用来控制分步方案产生的输出量。
   * @note
   * 这个类的许多功能都可以在DiscreteTime中使用，具有更现代的接口和更好的编程保证。考虑使用DiscreteTime而不是TimestepControl。
   *
   */
  class TimestepControl : public Subscriptor
  {
  public:
    /**
     * 设置默认值的构造函数
     *
     */
    TimestepControl(double start      = 0.,
                    double final      = 1.,
                    double tolerance  = 1.e-2,
                    double start_step = 1.e-2,
                    double print_step = -1.,
                    double max_step   = 1.);

    /**
     * 声明参数处理程序的控制参数。
     *
     */
    static void
    declare_parameters(ParameterHandler &param);

    /**
     * 从一个参数处理程序中读取控制参数。
     * 这个函数也会调用restart()，根据刚刚读取的参数，将这个类的所有其他内部参数重置为相应的值。
     *
     */
    void
    parse_parameters(ParameterHandler &param);

    /**
     * 返回时间间隔的左端。
     *
     */
    double
    start() const;
    /**
     * 返回时间间隔的右端。控制机制确保最后的时间步骤在这一点上结束。
     *
     */
    double
    final() const;
    /**
     * 返回控制时间步骤的公差值。
     *
     */
    double
    tolerance() const;
    /**
     * 返回当前时间步长的大小。
     *
     */
    double
    step() const;

    /**
     * 返回当前的时间。
     *
     */
    double
    now() const;

    /**
     * 计算下一步的大小，如果它与当前步长不同，则返回true。用新的步长推进当前时间。
     *
     */
    bool
    advance();

    /**
     * 设置起始值。
     *
     */
    void
    start(double);
    /**
     * 设置最终时间值。
     *
     */
    void
    final(double);
    /**
     * 设置公差
     *
     */
    void
    tolerance(double);

    /**
     * 设置第一个步骤的大小。这可能被时间步进策略所覆盖。
     * @param[in] 步长
     * 第一步的大小，这可能被时间步长策略覆盖。
     *
     */
    void
    start_step(const double step);

    /**
     * 设置最大步长的大小。
     *
     */
    void
    max_step(double);

    /**
     * 设置now()等于start()。将step()和print()初始化为其初始值。
     *
     */
    void
    restart();

    /**
     * 如果这个时间步长应该被写入磁盘，则返回true。
     *
     */
    bool
    print();

  private:
    /**
     * 时间间隔的开始。
     *
     */
    double start_val;

    /**
     * 时间间隔的结束。
     *
     */
    double final_val;

    /**
     * 控制时间步骤的公差值。
     *
     */
    double tolerance_val;

    /**
     * 第一个步骤的大小。
     *
     */
    double start_step_val;

    /**
     * 最大的步长。
     *
     */
    double max_step_val;

    /**
     * 最小步长。
     *
     */
    double min_step_val;

    /**
     * 当前时间步长的大小。如果我们的目标是 @p final_val.
     * ，这可能与 @p step_val, 不同。
     *
     */
    double current_step_val;

    /**
     * 由策略决定的当前时间步长的大小。如果我们的目标是
     * @p final_val. ，这可能与 @p current_step_val, 不同。
     *
     */
    double step_val;

    /**
     * 的当前时间。
     *
     */
    double now_val;

    /**
     * 决定了生成的输出之间的大致时间间隔。
     * 如果是负数，输出将在所有时间步骤中产生。
     *
     */
    double print_step;

    /**
     * 如果当前时间超过此值，则是生成输出的时间。
     *
     */
    double next_print_val;
  };


  inline double
  TimestepControl::start() const
  {
    return start_val;
  }


  inline double
  TimestepControl::final() const
  {
    return final_val;
  }


  inline double
  TimestepControl::step() const
  {
    return current_step_val;
  }


  inline double
  TimestepControl::tolerance() const
  {
    return tolerance_val;
  }


  inline double
  TimestepControl::now() const
  {
    return now_val;
  }


  inline void
  TimestepControl::start(double t)
  {
    start_val = t;
  }


  inline void
  TimestepControl::final(double t)
  {
    final_val = t;
  }


  inline void
  TimestepControl::tolerance(double t)
  {
    tolerance_val = t;
  }


  inline void
  TimestepControl::start_step(const double t)
  {
    start_step_val = t;
  }


  inline void
  TimestepControl::max_step(double t)
  {
    max_step_val = t;
  }


  inline void
  TimestepControl::restart()
  {
    now_val          = start_val;
    step_val         = start_step_val;
    current_step_val = step_val;
    if (print_step > 0.)
      next_print_val = now_val + print_step;
    else
      next_print_val = now_val - 1.;
  }

} // namespace Algorithms

DEAL_II_NAMESPACE_CLOSE

#endif


