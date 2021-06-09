//include/deal.II-translator/base/conditional_ostream_0.txt
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

#ifndef dealii_conditional_ostream_h
#define dealii_conditional_ostream_h

#include <deal.II/base/config.h>

#include <ostream>

DEAL_II_NAMESPACE_OPEN


/**
 * 一个允许打印到输出流的类，例如 @p std::cout,
 * ，取决于ConditionalOStream对象是否被激活（默认）。这个对象的条件可以通过set_condition()和构造函数来改变。这个类在
 * step-17 、 step-18 、 step-32 、 step-33 和 step-35
 * 的教程程序中使用。
 * 这个类主要在并行计算中有用。通常，你会用 @p std::cout
 * 来打印信息，比如程序目前正在做什么，或者每一步的自由度数量。然而，在并行程序中，这意味着每个MPI进程都要写到屏幕上，这就产生了许多重复的相同文本。为了避免这种情况，我们必须有一个指定的进程，例如MPI进程编号为0的进程，进行输出，并且用一个if-条件来保护每个写入语句。这就变得很麻烦，而且使代码变得很混乱。与其这样做，不如使用本类：其类型的对象就像标准输出流一样，但它们只根据一个条件来打印东西，这个条件可以设置为，例如，<tt>mpi_process==0</tt>，这样，只有一个进程有一个真实的条件，在所有其他进程中，对这个对象的写入只是涅槃消失。
 * 这个类的通常用法是这样的。
 *
 *
 * @code
 * ConditionalOStream pout(std::cout, this_mpi_process==0);
 *
 * // all processes print the following information to standard output
 * std::cout << "Reading parameter file on process "
 *         << this_mpi_process << std::endl;
 *
 * // following is printed by process 0 only
 * pout << "Solving ..." << std::endl;
 * solve();
 * pout << "done" << std::endl;
 * @endcode
 *
 * 这里，"在进程xy上读取参数文件
 * "是由每个进程单独打印的。与此相反，`解决...'和`完成'只打印到标准输出一次，即由进程0打印。
 * 这个类不是从ostream派生的。因此
 *
 * @code
 * system_matrix.print_formatted(pout);
 * @endcode
 * 是 <em> 不是 </em>
 * 的可能。相反，使用is_active()函数来解决这个问题。
 *
 *
 * @code
 * if (pout.is_active())
 * system_matrix.print_formatted(cout);
 * @endcode
 *
 *
 * @ingroup textoutput
 *
 */
class ConditionalOStream
{
public:
  /**
   * 构造函数。设置我们要写入的流，以及实际转发写入的基础条件。默认情况下，一个对象的条件是活动的。
   *
   */
  ConditionalOStream(std::ostream &stream, const bool active = true);

  /**
   * 根据<tt>active</tt>标志将该流的条件设置为active（真）或非active（假）。当且仅当它的条件是激活的，这个类的一个对象就会打印到<tt>cout</tt>。
   *
   */
  void
  set_condition(const bool active);

  /**
   * 返回该对象的条件。
   *
   */
  bool
  is_active() const;

  /**
   * 返回对当前使用的流的引用。
   *
   */
  std::ostream &
  get_stream() const;

  /**
   * 通过这个流输出一个常量的东西。这个函数必须是 @p
   * 常数，这样，这个类型的成员对象也可以从周围类的 @p
   * const 成员函数中使用。
   *
   */
  template <typename T>
  const ConditionalOStream &
  operator<<(const T &t) const;

  /**
   * 处理ostream操纵器。这个函数必须是 @p const
   * ，这样，这个类型的成员对象也可以从周围类的 @p const
   * 成员函数中使用。
   * 注意，编译器希望看到这与上面的一般模板不同的处理方式，因为像
   * @p std::endl
   * 这样的函数实际上是重载的，不能直接绑定到模板类型。
   *
   */
  const ConditionalOStream &
  operator<<(std::ostream &(*p)(std::ostream &)) const;

private:
  /**
   * 对我们要写入的流的引用。
   *
   */
  std::ostream &output_stream;

  /**
   * 存储对象所处的实际条件。
   *
   */
  bool active_flag;
};


// --------------------------- inline and template functions -----------

template <class T>
inline const ConditionalOStream &
ConditionalOStream::operator<<(const T &t) const
{
  if (active_flag == true)
    output_stream << t;

  return *this;
}


inline const ConditionalOStream &
ConditionalOStream::operator<<(std::ostream &(*p)(std::ostream &)) const
{
  if (active_flag == true)
    output_stream << p;

  return *this;
}


inline std::ostream &
ConditionalOStream::get_stream() const
{
  return output_stream;
}


DEAL_II_NAMESPACE_CLOSE

#endif


