//include/deal.II-translator/algorithms/theta_timestepping_0.txt
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


#ifndef dealii_theta_timestepping_h
#define dealii_theta_timestepping_h

#include <deal.II/base/config.h>

#include <deal.II/algorithms/operator.h>
#include <deal.II/algorithms/timestep_control.h>

#include <deal.II/base/smartpointer.h>

DEAL_II_NAMESPACE_OPEN

// Forward declaration
#ifndef DOXYGEN
class ParameterHandler;
#endif

namespace Algorithms
{
  /**
   * 一个小结构，收集一个时间步长的大小和当前时间。时间步进方案可以使用它来向实际执行单步的类提供时间步进信息。
   * 什么是 "当前时间
   * "的定义取决于该方案。对于一个显式方案来说，这就是步长开始时的时间。对于隐式方案，它通常是结束时的时间。
   *
   */
  struct TimestepData
  {
    /// The current time
    double time;
    /// The current step size times something
    double step;
  };

  /**
   * 执行theta时间步进方案的应用类。    Theta方案是对隐式和显式欧拉方案、Crank-Nicholson方案以及这些方案的线性组合的抽象。实际方案的选择由参数#theta控制，如下所示。    <ul>   <li>  #theta=0：显式欧拉方案  <li>  #theta=1：隐式欧拉方案  <li>  #theta=½。Crank-Nicholson方案  </ul>  对于固定的#theta，Crank-Nicholson方案是唯一的二阶方案。然而，通过选择大于1/2的#theta可以实现进一步的稳定性，从而引入一个一阶误差项。为了避免收敛顺序的损失，可以使用自适应theta方案，其中<i>#theta=½+c dt</i>。    假设我们要解决方程<i>u' + F(u) = 0</i>，步长为<i>k</i>。 Theta方案的一个步骤可以写成@f[
   * M u_{n+1} + \theta k F(u_{n+1})  = M u_n
   *
   * - (1-\theta)k F(u_n). @f]
   * 这里，<i>M</i>是质量矩阵。我们看到，右手边相当于一个显式欧拉步骤，其步骤大小以弱形式修改（直到M的反转）。左手边对应的是一个隐式欧拉步长的修改步长（右手边给出）。因此，theta方案的实现将使用两个运算器对象，一个用于显式部分，一个用于隐式部分。每个对象都将使用自己的TimestepData来说明修改后的步长（如果问题不是自主的，则使用不同的时间）。请注意，一旦显式部分被计算出来，左手边实际上构成了一个必须被解决的线性或非线性系统。
   * <h3>Usage AnyData</h3>
   * ThetaTimestepping使用AnyData来交流向量和时间步长信息。与外部或内部的操作者对象。它本身并不使用所提供的输入向量，而是将它们转发给显式和隐式算子。
   * <h4>Vector data</h4>
   * 显式运算符#op_explicit在其输入中首先接收向量 "Previous
   * iterate"，这是上一个时间步长后的解值。接下来是提供给
   * ThetaTimestepping::operator()
   * 的所有向量作为输入参数。#op_explicit应该将其结果写入其输出参数的第一个位置，标为
   * "结果"。
   * 隐式运算符#op_implicit在其第一个输入向量中接收#op_explicit的结果，标签为
   * "前次"。它的后面是作为输入参数提供给
   * ThetaTimestepping::operator()
   * 的所有向量。#op_implicit的输出被直接写入给ThetaTimestepping的输出参数中。
   * <h4>Scalar data</h4>
   * 自从引入AnyData后，ThetaTimestepping也能通过AnyData传递当前时间步长信息。
   * 因此，作为#op_explicit和#op_implicit输入的AnyData对象包含两个类型为`const
   * double*`的条目，名为 "Time "和 "Timestep"。请注意，"时间
   * "对于#op_explicit和#op_implicit来说，分别指的是当前步骤开始的时间和结束的时间。
   * <h3>Usage of ThetaTimestepping</h3>
   * ThetaTimestepping的使用比牛顿等更复杂，因为内部运算符通常需要访问TimeStepData。
   * 因此，我们有一个信息的循环依赖，我们包括下面的例子来说明其使用。
   * 首先，我们定义ThetaTimestepping使用的两个算子，并将其称为
   * <code>Implicit</code> and <code>Explicit</code>
   * 。它们都共享Operator的公共接口，另外还提供了要使用的矩阵的存储空间和一个指向TimestepData的指针。注意，我们在这里没有使用SmartPointer，因为TimestepData将在操作者之前被销毁。
   * @code
   * class Explicit : public OperatorBase
   * {
   * public:
   * Explicit(const FullMatrix<double> &matrix);
   * void operator()(AnyData &out, const AnyData &in);
   *
   * private:
   * SmartPointer<const FullMatrix<double>, Explicit> matrix;
   * FullMatrix<double>                               m;
   * };
   *
   * class Implicit : public OperatorBase
   * {
   * public:
   * Implicit(const FullMatrix<double> &matrix);
   * void operator()(AnyData &out, const AnyData &in);
   *
   * private:
   * SmartPointer<const FullMatrix<double>, Implicit> matrix;
   * FullMatrix<double>                               m;
   * };
   * @endcode
   * 这些操作符将在主程序之后实现。但让我们先看看它们是如何被使用的。首先，让我们定义一个用于我们系统的矩阵和一个OutputOperator，以便将每个时间段的数据写到一个文件中。
   * @code
   * int main()
   * {
   * FullMatrix<double> matrix(2);
   * matrix(0, 0) = 0.;
   * matrix(1, 1) = 0.;
   * matrix(0, 1) = 3.14;
   * matrix(1, 0) =
   *
   * -3.14;
   *
   * OutputOperator<Vector<double>> out;
   * out.initialize_stream(std::cout);
   * @endcode
   * 现在我们为步骤的隐式和显式部分以及ThetaTimestepping本身创建对象。我们用输出运算符初始化timestepping，以便能够看到每个步骤中的输出。
   * @code
   * Explicit                          op_explicit(matrix);
   * Implicit                          op_implicit(matrix);
   * ThetaTimestepping<Vector<double>> solver(op_explicit, op_implicit);
   * solver.set_output(out);
   * @endcode
   * 下一步是提供要使用的向量。<tt>value</tt>被填入初始值，也是每个时间步长的解决方案的向量。因为Operator的接口必须能够处理多个向量，所以我们需要将其存储在一个AnyData对象中。
   * 由于我们的问题没有额外的参数，输入的AnyData对象仍然是空的。
   * @code
   * Vector<double> value(2);
   * value(0) = 1.;
   * AnyData indata;
   * AnyData outdata;
   * outdata.add(&value, "value");
   * @endcode
   * 最后，我们准备告诉求解器，我们将从初始时间步长开始，并运行它。
   * @code
   * solver.notify(Events::initial);
   * solver(outdata, indata);
   * }
   * @endcode
   * 除了主函数，我们还需要定义隐式和显式运算符的成员函数。
   * 首先是构造函数，它只是将系统矩阵复制到成员指针中，供以后使用。
   * @code
   * Explicit::Explicit(const FullMatrix<double> &M)
   * : matrix(&M)
   * {
   * m.reinit(M.m(), M.n());
   * }
   * @endcode
   * 现在我们需要研究隐式和显式运算符的应用。我们假设指针
   * <code>matrix</code>
   * 指向主程序中创建的矩阵（构造函数为我们做了这个）。
   * 在这里，我们首先从作为输入的AnyData对象中获得时间步长。然后，如果我们处于第一步或者时间步长发生了变化，我们就填充本地矩阵
   * $m$ ，这样与给定的矩阵 $M$ 一起就变成了\f[ m = I
   *
   * - \Delta t M. \f]
   * 在我们处理完通知后，我们清除它们，这样矩阵就只在必要时生成。
   * @code
   * void Explicit::operator()(AnyData &out, const AnyData &in)
   * {
   * const double timestep =in.read_ptr<double>("Timestep");
   * if (this->notifications.test(Events::initial) ||
   *     this->notifications.test(Events::new_timestep_size))
   *   {
   *     m.equ(-timestep,matrix);
   *     for (unsigned int i = 0; i < m.m(); ++i)
   *       m(i, i) += 1.;
   *   }
   * this->notifications.clear();
   * @endcode
   * 现在我们将输入向量与新矩阵相乘，并存储在输出上。
   * @code
   * m.vmult(*out.entry<Vector<double>>(0),
   *        in.read_ptr<Vector<double>>("Previous iterate"));
   * }
   * @endcode
   * 隐式运算符的代码几乎相同，只是我们改变了时间步长前面的符号，并使用了矩阵的倒数。
   * @code
   * Implicit::Implicit(const FullMatrix<double> &M)
   * : matrix(&M)
   * {
   * m.reinit(M.m(), M.n());
   * }
   *
   * void Implicit::operator()(AnyData &out, const AnyData &in)
   * {
   * const double timestep =in.read_ptr<double>("Timestep");
   * if (this->notifications.test(Events::initial) ||
   *     this->notifications.test(Events::new_timestep_size))
   *   {
   *     m.equ(timestep,matrix);
   *     for (unsigned int i = 0; i < m.m(); ++i)
   *       m(i, i) += 1.;
   *     m.gauss_jordan();
   *   }
   * this->notifications.clear();
   * m.vmult(*out.entry<Vector<double>>(0),
   *        in.read_ptr<Vector<double>>("Previous time"));
   * }
   * @endcode
   *
   *
   */
  template <typename VectorType>
  class ThetaTimestepping : public OperatorBase
  {
  public:
    /**
     * 构造函数，接收存储在#op_explicit和#op_implicit中的两个运算符。关于它们的含义，见这些变量的描述。
     *
     */
    ThetaTimestepping(OperatorBase &op_explicit, OperatorBase &op_implicit);

    /**
     * 时隙方案。          @param
     * in被ThetaTimestepping忽略，但被合并到AnyData对象中，作为运算符#op_explicit和#op_implicit的输入。
     * @param
     * out在其第一个参数中必须包含一个指向VectorType实例的指针，它包含运算符被调用时的初始值。
     * 它包含运算符返回时的最终值。
     *
     */
    virtual void
    operator()(AnyData &out, const AnyData &in) override;

    /**
     * 注册一个由外部迭代触发的事件。
     *
     */
    virtual void
    notify(const Event &) override;

    /**
     * 定义一个运算符，它将在每个步骤中输出结果。请注意，没有这个就不会产生输出。
     *
     */
    void
    set_output(OutputOperator<VectorType> &output);

    /**
     * 在一个参数处理程序中声明参数。
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
     * 在时序计划中的当前时间。
     *
     */
    double
    current_time() const;

    /**
     * 隐性和显性部分之间的权重。
     *
     */
    double
    theta() const;

    /**
     * 设置一个新的权重并返回旧的权重
     *
     */
    double
    theta(double new_theta);

    /**
     * 交给#op_explicit time stepping operator的数据。
     * 这里的时间是当前步骤开始时的时间，时间步长是实际时间步长的（1-#theta）倍。
     *
     */
    const TimestepData &
    explicit_data() const;

    /**
     * 交给#op_implicit时间步长运算器的数据。
     * 这里的时间是当前步骤开始时的时间，时间步长是实际时间步长的#theta倍。
     *
     */
    const TimestepData &
    implicit_data() const;

    /**
     * 允许对控制对象的访问。
     *
     */
    TimestepControl &
    timestep_control();

  private:
    /**
     * 控制时间步长和计算每步新时间的对象。
     *
     */
    TimestepControl control;

    /**
     * 控制参数theta，范围为<tt>[0,1]</tt>。它的默认值是0.5。
     *
     */
    double vtheta;
    /**
     * 如果<tt>true</tt>，使用自适应的#theta。还没有实现。
     *
     */
    bool adaptive;

    /**
     * 方案的显式部分的数据。
     *
     */
    TimestepData d_explicit;

    /**
     * 方案的隐性部分的数据。
     *
     */
    TimestepData d_implicit;


    /**
     * 计算方案的显式部分的运算器。这将在其输入数据中收到当前时间的值，名称为
     * "当前时间解决方案"。它应该从
     * explicit_data()获得当前时间和时间步长。
     * 它的返回值是  $ Mu+cF(u) $  ，其中  $u$
     * 是当前状态向量，  $M$  是质量矩阵，  $F$
     * 是空间算子，  $c$  是调整的时间步长  $(1-\theta) \Delta
     * t$  。
     *
     */
    SmartPointer<OperatorBase, ThetaTimestepping<VectorType>> op_explicit;

    /**
     * 算子解决方案的隐含部分。它将在其输入数据中收到矢量
     * "上一时间"。关于时间步长的信息应该从
     * implicit_data()中获得。
     * 它的返回值是<i>Mu-cF(u)=f</i>的解<i>u</i>，其中<i>f</i>是在输入数据的
     * "前次
     * "条目中找到的对偶空间向量，<i>M</i>是质量矩阵，<i>F</i>是空间中的运算器，<i>c</i>是调整后的时间步长
     * $ \theta \Delta t$
     *
     */
    SmartPointer<OperatorBase, ThetaTimestepping<VectorType>> op_implicit;

    /**
     * 运算器在每个时间步长中写出的输出
     *
     */
    SmartPointer<OutputOperator<VectorType>, ThetaTimestepping<VectorType>>
      output;
  };


  template <typename VectorType>
  inline const TimestepData &
  ThetaTimestepping<VectorType>::explicit_data() const
  {
    return d_explicit;
  }


  template <typename VectorType>
  inline const TimestepData &
  ThetaTimestepping<VectorType>::implicit_data() const
  {
    return d_implicit;
  }


  template <typename VectorType>
  inline TimestepControl &
  ThetaTimestepping<VectorType>::timestep_control()
  {
    return control;
  }

  template <typename VectorType>
  inline void
  ThetaTimestepping<VectorType>::set_output(OutputOperator<VectorType> &out)
  {
    output = &out;
  }


  template <typename VectorType>
  inline double
  ThetaTimestepping<VectorType>::theta() const
  {
    return vtheta;
  }


  template <typename VectorType>
  inline double
  ThetaTimestepping<VectorType>::theta(double new_theta)
  {
    const double tmp = vtheta;
    vtheta           = new_theta;
    return tmp;
  }


  template <typename VectorType>
  inline double
  ThetaTimestepping<VectorType>::current_time() const
  {
    return control.now();
  }
} // namespace Algorithms

DEAL_II_NAMESPACE_CLOSE

#endif


