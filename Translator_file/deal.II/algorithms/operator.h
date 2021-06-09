//include/deal.II-translator/algorithms/operator_0.txt
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


#ifndef dealii_operator_h
#define dealii_operator_h

#include <deal.II/base/config.h>

#include <deal.II/algorithms/any_data.h>

#include <deal.II/base/event.h>

#include <fstream>

DEAL_II_NAMESPACE_OPEN

/**
 * 包含统一形式的数字算法的命名空间。
 * 这个命名空间中的所有算法类都是从Operator或OutputOperator派生出来的，这取决于它们是否返回一个值。关于如何使用这些类的更多详细信息，请参见这些类的文档。
 *
 *
 */
namespace Algorithms
{
  /**
   * @todo  更新本文档和Operator的文档
   * 本库中所有算法的抽象基类。操作符是一个具有operator()的对象，它将一组命名的向量转换成另一组命名的向量。
   * 此外，一个运算符可以被调用例程通知参数的变化。外层迭代可以通知()运算器一个事件，例如，网格的改变、不同的时间步长或牛顿方法的收敛速度太慢，这将触发矩阵的重新组合或类似的事情。
   * <h3>Usage for nested iterations</h3>
   * 这可能是Operator最突出的用途，一个外部迭代方法调用一个内部求解器，等等。通常，在这样的嵌套系统中，最内层的方法必须使用所有外层迭代的值来计算一个残差。由于这种嵌套的深度和顺序在设计通用工具时很难预测，我们使用AnyData来访问这些向量。通常，<tt>out</tt>中的第一个向量包含调用operator()()时的起始向量，以及函数返回时的解决方案。对象<tt>in</tt>是提供额外的信息，并转发给嵌套迭代的内部Operator对象。
   *
   */
  class OperatorBase : public Subscriptor
  {
  public:
    /**
     * 虚拟析构器。
     *
     */
    virtual ~OperatorBase() override = default;

    /**
     * 实际的操作，在一个派生类中实现。
     *
     */
    virtual void
    operator()(AnyData &out, const AnyData &in) = 0;

    /**
     * 注册一个由外部迭代触发的事件。
     *
     */
    virtual void
    notify(const Event &);

    /**
     * 清除所有的#通知。
     *
     */
    void
    clear_events();

  protected:
    /**
     * 在这里累积事件。如果其中任何一个被设置，终端应用程序的函数solve()必须负责重新组装矩阵。
     *
     */
    Event notifications;
  };

  /**
   * 一个单数运算器基类，目的是在迭代的每个步骤中输出AnyData中的向量。
   *
   */
  template <typename VectorType>
  class OutputOperator : public Subscriptor
  {
  public:
    /**
     * 用无效数据初始化成员变量的构造器。
     *
     */
    OutputOperator();

    /**
     * 复制构造函数被删除，因为这个类的对象不应该被复制。
     *
     */
    OutputOperator(const OutputOperator<VectorType> &) = delete;

    /**
     * 空的虚拟析构器。
     *
     */
    virtual ~OutputOperator() override = default;

    /**
     * 设置数据被写入的流 @p os
     * 。如果没有用这个函数选择流，数据将进入 @p deallog.
     * 。
     *
     */
    void
    initialize_stream(std::ostream &stream);

    /**
     * 设置当前的步骤。
     *
     */
    void
    set_step(const unsigned int step);

    /**
     * 输出AnyData中的所有向量。
     *
     */
    virtual OutputOperator<VectorType> &
    operator<<(const AnyData &vectors);

  protected:
    unsigned int step;

  private:
    std::ostream *os;
  };

  template <typename VectorType>
  inline void
  OutputOperator<VectorType>::set_step(const unsigned int s)
  {
    step = s;
  }


  /**
   * 在OutputOperator中通过移位一个整数值来设置步数。
   * @relatesalso  OutputOperator
   *
   */
  template <typename VectorType>
  inline OutputOperator<VectorType> &
  operator<<(OutputOperator<VectorType> &out, unsigned int step)
  {
    out.set_step(step);
    return out;
  }
} // namespace Algorithms

DEAL_II_NAMESPACE_CLOSE

#endif


