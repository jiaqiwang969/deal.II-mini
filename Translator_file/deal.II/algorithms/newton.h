//include/deal.II-translator/algorithms/newton_0.txt
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


#ifndef dealii_newton_h
#define dealii_newton_h

#include <deal.II/base/config.h>

#include <deal.II/algorithms/any_data.h>
#include <deal.II/algorithms/operator.h>

#include <deal.II/base/smartpointer.h>

#include <deal.II/lac/solver_control.h>

DEAL_II_NAMESPACE_OPEN

// Forward declaration
#ifndef DOXYGEN
class ParameterHandler;
#endif

namespace Algorithms
{
  /**
   * 执行牛顿迭代的运算器类，具有标准步长控制和自适应矩阵生成功能。
   * 这个类执行牛顿迭代，直到由#控制决定的收敛。如果在更新之后，残差的标准值变大，那么步长控制就会被激活，随后更新被除以2，直到残差真正变小（或者达到#n_stepsize_iterations确定的最小缩放系数）。
   * 由于组装矩阵，取决于实现方式，往往是昂贵的，这个方法应用一个自适应的重新组装策略。只有当残差的减少系数超过#阈值时，事件
   * Algorithms::bad_derivative
   * 才会提交给#inverse_derivative。由这个对象来实现相应的重新组合。
   * <h3>Contents of the AnyData objects</h3>
   * 牛顿方法使用的唯一数值是operator()()参数<tt>out</tt>中的第一个向量。它作为牛顿方法的起始向量，并在最后包含了解。所有其他<tt>out</tt>的向量都被牛顿方法及其内部的operator对象所忽略。所有<tt>in</tt>的向量都会被转发给内部的运算符对象，并添加如下的附加信息。
   * 当调用(*#residual)()时，给牛顿迭代的AnyData
   * <tt>in</tt>被一个向量<tt>"Newton
   * iterate"</tt>所预置，这是牛顿迭代的当前值，可用于评估此时的残差。
   * 对于调用(*#inverse_derivative)，向量<tt>"Newton
   * residual"</tt>被插入<tt>"Newton iterate"</tt>之前。
   *
   */
  template <typename VectorType>
  class Newton : public OperatorBase
  {
  public:
    /**
     * 构造函数，分别接收计算残差和解决线性问题的应用程序。
     *
     */
    Newton(OperatorBase &residual, OperatorBase &inverse_derivative);

    /**
     * 声明适用于牛顿方法的参数。
     *
     */
    static void
    declare_parameters(ParameterHandler &param);

    /**
     * 读取ParameterHandler中的参数。
     *
     */
    void
    parse_parameters(ParameterHandler &param);

    /**
     * 初始化用于调试的指针data_out。
     *
     */
    void
    initialize(OutputOperator<VectorType> &output);

    /**
     * 实际的牛顿迭代。初始值在<tt>out(0)</tt>中，它也包含收敛后的结果。<tt>in</tt>中的值不会被牛顿使用，但会被传给对象#residual和#inverse_derivative。
     *
     */
    virtual void
    operator()(AnyData &out, const AnyData &in) override;

    virtual void
    notify(const Event &) override;

    /**
     * 设置允许的最大残差减少量，而不触发下一步的装配。返回之前的值。
     *
     */
    double
    threshold(double new_value);

    /**
     * 牛顿迭代的控制对象。
     *
     */
    ReductionControl control;

  private:
    /**
     * 计算残差的运算器。
     *
     */
    SmartPointer<OperatorBase, Newton<VectorType>> residual;

    /**
     * 将反导数应用于残差的运算器。
     *
     */
    SmartPointer<OperatorBase, Newton<VectorType>> inverse_derivative;

    /**
     * 在debug_vectors为真时处理输出的运算器。
     * 首先调用初始化函数。
     *
     */
    SmartPointer<OutputOperator<VectorType>, Newton<VectorType>> data_out;

    /**
     * 这个标志是由函数assemble()设置的，表示矩阵在启动时必须重新组装。
     *
     */
    bool assemble_now;

    /**
     * 一个用于决定应该进行多少次步长迭代的标志。
     * 默认值为原始值21。        在此输入0以关闭步长控制。
     * @note  由参数文件中的<tt>步长迭代</tt>控制。
     *
     */
    unsigned int n_stepsize_iterations;

    /**
     * 重新组合矩阵的阈值。
     * 如果两个连续残差的商小于这个阈值，系统矩阵就不会在这一步中被组装起来。
     * @note 该参数应调整为内解器的残差增益。
     * 默认值为零，导致在每一个牛顿步骤中重新组装。
     *
     */
    double assemble_threshold;

  public:
    /**
     * 在每一步之后打印残差、更新和更新的解决方案到文件<tt>Newton_NNN</tt>?
     *
     */
    bool debug_vectors;
    /**
     * 将调试输出写到 @p deallog; ，数字越大，输出越多。
     *
     */
    unsigned int debug;
  };
} // namespace Algorithms

DEAL_II_NAMESPACE_CLOSE

#endif


