//include/deal.II-translator/base/event_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2019 by the deal.II authors
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


#ifndef dealii_event_h
#define dealii_event_h

#include <deal.II/base/config.h>

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

DEAL_II_NAMESPACE_OPEN

namespace Algorithms
{
  /**
   * 这类对象用于通知内部应用程序由外循环引发的变化。它们通过
   * Operator::notify()
   * 被交给应用程序，如何处理它们由实际的应用程序决定。
   * 事件被组织成一个可扩展的二进制枚举器。每个类都可以使用assign()添加自己的事件。一个典型的代码例子是
   * @code
   * class A
   * {
   * static Event event;
   * };
   *
   * Event A::event = Event::assign("Event for A");
   * @endcode
   *
   *
   */
  class Event
  {
  public:
    /**
     * 这个函数注册了一个新的事件类型，并给它分配了一个唯一的标识符。这个函数的结果应该被储存起来以便以后使用。
     *
     */
    static Event
    assign(const std::string &name);

    /**
     * 如果你忘了存储assign的结果，下面是如何在知道名字的情况下检索它。
     *
     */
    //      static Event find(const std::string& name);

    /**
     * 构造函数，生成一个清晰的事件。
     *
     */
    Event();

    /**
     * 清除所有标志
     *
     */
    void
    clear();

    /**
     * 设置所有标志
     *
     */
    void
    all();

    /**
     * 添加其他事件的标志
     *
     */
    Event &
    operator+=(const Event &event);

    /**
     * 清除另一个事件的标志
     *
     */
    Event &
    operator-=(const Event &event);

    /**
     * 测试在其他事件中设置的所有标志是否也在这个事件中设置。
     *
     */
    bool
    test(const Event &event) const;

    /**
     * 如果有事件被设置，返回<tt>true</tt>。
     *
     */
    bool
    any() const;

    /**
     * 列出到一个流的标志。
     *
     */
    template <class OS>
    void
    print(OS &os) const;

    /**
     * 列出所有分配的事件。
     *
     */
    template <class OS>
    static void
    print_assigned(OS &os);

  private:
    /**
     * 有时，必须通过各种手段来采取行动。因此，如果这个值为真，test()总是返回真。
     *
     */
    bool all_true;

    /**
     * 事件的实际列表
     *
     */
    std::vector<bool> flags;

    /**
     * 注册事件的名称
     *
     */
    // TODO: This static field must be guarded by a mutex to be thread-safe!
    static std::vector<std::string> names;
  };

  /**
   * 库操作者使用的事件
   *
   */
  namespace Events
  {
    /**
     * 程序刚刚开始，一切都应该是新的。
     *
     */
    extern const Event initial;

    /**
     * 网格已经改变。
     *
     */
    extern const Event remesh;

    /**
     * 当前的导数导致牛顿方法收敛缓慢。
     *
     */
    extern const Event bad_derivative;

    /**
     * 时间步长方案开始了新的时间步长。
     *
     */
    extern const Event new_time;

    /**
     * 时间步长方案改变了时间步长的大小。
     *
     */
    extern const Event new_timestep_size;
  } // namespace Events


  //----------------------------------------------------------------------//


  inline bool
  Event::any() const
  {
    if (all_true)
      return true;
    return std::find(flags.begin(), flags.end(), true) != flags.end();
  }


  inline bool
  Event::test(const Event &event) const
  {
    // First, test all_true in this
    if (all_true)
      return true;

    const unsigned int n     = flags.size();
    const unsigned int m     = event.flags.size();
    const unsigned int n_min = (n < m) ? n : m;

    // Now, if all_true set in the
    // other, then all must be true
    // in this
    if (event.all_true)
      {
        // Non existing flags are
        // always assumed false
        if (m > n)
          return false;

        // Test all flags separately
        // and return false if one is
        // not set
        return std::find(flags.begin(), flags.end(), false) == flags.end();
      }

    // Finally, compare each flag
    // separately
    for (unsigned int i = 0; i < n_min; ++i)
      if (event.flags[i] && !flags[i])
        return false;
    for (unsigned int i = n_min; i < m; ++i)
      if (event.flags[i])
        return false;
    return true;
  }



  inline Event &
  Event::operator+=(const Event &event)
  {
    all_true |= event.all_true;
    if (all_true)
      return *this;

    if (flags.size() < event.flags.size())
      flags.resize(event.flags.size());
    for (unsigned int i = 0; i < event.flags.size(); ++i)
      flags[i] = flags[i] || event.flags[i];

    return *this;
  }


  inline Event &
  Event::operator-=(const Event &event)
  {
    if (!event.any())
      return *this;

    all_true = false;
    if (event.all_true)
      {
        std::fill(flags.begin(), flags.end(), false);
        return *this;
      }

    if (flags.size() < event.flags.size())
      flags.resize(event.flags.size());
    for (unsigned int i = 0; i < event.flags.size(); ++i)
      if (event.flags[i])
        flags[i] = false;

    return *this;
  }


  template <class OS>
  inline void
  Event::print(OS &os) const
  {
    if (all_true)
      os << " ALL";

    for (unsigned int i = 0; i < flags.size(); ++i)
      if (flags[i])
        os << ' ' << names[i];
  }


  template <class OS>
  inline void
  Event::print_assigned(OS &os)
  {
    for (unsigned int i = 0; i < names.size(); ++i)
      os << i << '\t' << names[i] << std::endl;
  }


  /**
   * 输出事件的移位运算符。调用  Event::print().   @relatesalso
   * 事件
   *
   */
  template <class OS>
  OS &
  operator<<(OS &o, const Event &e)
  {
    e.print(o);
    return o;
  }
} // namespace Algorithms

DEAL_II_NAMESPACE_CLOSE

#endif


