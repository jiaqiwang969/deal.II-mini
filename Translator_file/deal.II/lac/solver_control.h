//include/deal.II-translator/lac/solver_control_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2020 by the deal.II authors
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

#ifndef dealii_solver_control_h
#define dealii_solver_control_h


#include <deal.II/base/config.h>

#include <deal.II/base/subscriptor.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN

#ifndef DOXYGEN
class ParameterHandler;
#endif

 /*!@addtogroup Solvers */ 
 /*@{*/ 

/**
 * 控制类，用于确定迭代求解器的收敛性。
 * 由迭代方法使用，以确定是否应该继续迭代。为此，虚函数<tt>check()</tt>在每次迭代中被调用，并给出当前的迭代步骤和表示收敛的值（通常是残差）。
 * 迭代结束后，可以用函数last_value()和last_step()来获取迭代的最终状态信息。
 * check()可以在派生类中被替换，以允许进行更复杂的测试。
 *
 *  <h3>State</h3>
 * 检查函数的返回状态是#State类型的，它是本类的一个枚举。它表示解算器所处的状态。
 * State的可能值是  <ul>   <li>  <tt>iterate = 0</tt>: 继续迭代。  <li>   @p success:  目标达成，迭代方法可以成功终止。  <li>   @p failure:  迭代方法应该停止，因为在给定的最大迭代次数内无法实现或至少无法实现收敛。  </ul>
 *
 */
class SolverControl : public Subscriptor
{
public:
  /**
   * 表示求解器可以处于的不同状态的枚举。更多信息请参见该类的一般文档。
   *
   */
  enum State
  {
    /// Continue iteration
    iterate = 0,
    /// Stop iteration, goal reached
    success,
    /// Stop iteration, goal not reached
    failure
  };



  /**
   * 迭代求解器收敛失败时抛出的类，当迭代次数超过极限或残差没有达到期望的极限时，例如崩溃的情况。
   * 最后一次迭代的残差，以及最后一步的迭代数都存储在这个对象中，并且可以在捕捉到这个类别的异常时恢复。
   *
   */

  class NoConvergence : public dealii::ExceptionBase
  {
  public:
    NoConvergence(const unsigned int last_step, const double last_residual)
      : last_step(last_step)
      , last_residual(last_residual)
    {}

    virtual ~NoConvergence() noexcept override = default;

    virtual void
    print_info(std::ostream &out) const override
    {
      out
        << "Iterative method reported convergence failure in step " << last_step
        << ". The residual in the last step was " << last_residual << ".\n\n"
        << "This error message can indicate that you have simply not allowed "
        << "a sufficiently large number of iterations for your iterative solver "
        << "to converge. This often happens when you increase the size of your "
        << "problem. In such cases, the last residual will likely still be very "
        << "small, and you can make the error go away by increasing the allowed "
        << "number of iterations when setting up the SolverControl object that "
        << "determines the maximal number of iterations you allow."
        << "\n\n"
        << "The other situation where this error may occur is when your matrix "
        << "is not invertible (e.g., your matrix has a null-space), or if you "
        << "try to apply the wrong solver to a matrix (e.g., using CG for a "
        << "matrix that is not symmetric or not positive definite). In these "
        << "cases, the residual in the last iteration is likely going to be large."
        << std::endl;
    }

    /**
     * 最后一步的迭代号。
     *
     */
    const unsigned int last_step;

    /**
     * 最后一步的残差。
     *
     */
    const double last_residual;
  };



  /**
   * 构造函数。参数 @p n 和 @p tol
   * 是失败前的最大迭代步数和确定迭代成功的公差。
   * @p log_history
   * 指定是否应将历史记录（即要检查的值和迭代步数）打印到
   * @p  deallog流。 默认为：不打印。同样， @p log_result
   * 指定最终结果是否被记录到 @p deallog. 中，默认为是。
   *
   */
  SolverControl(const unsigned int n           = 100,
                const double       tol         = 1.e-10,
                const bool         log_history = false,
                const bool         log_result  = true);

  /**
   * 需要虚拟析构器，因为这个类中有虚拟函数。
   *
   */
  virtual ~SolverControl() override = default;

  /**
   * 到参数文件的接口。
   *
   */
  static void
  declare_parameters(ParameterHandler &param);

  /**
   * 从文件中读取参数。
   *
   */
  void
  parse_parameters(ParameterHandler &param);

  /**
   * 决定一个迭代的成功或失败。
   * 这个函数获取当前的迭代步骤，以确定是否超过了允许的步骤数，在这种情况下返回
   * @p failure 。如果 @p check_value 低于规定的公差，则返回 @p
   * success.  在所有其他情况下，返回 @p iterate
   * 以建议继续进行迭代程序。
   * 如果残差变成了一个去规范化的值，迭代也会被终止（
   * @p NaN).  <tt>check()</tt>另外保留了 @p step  和  @p check_value.
   * 这些值可以通过<tt>last_value()</tt>和<tt>last_step()</tt>访问。
   * 派生类可以重载这个函数，例如，记录收敛指标（  @p
   * check_value)  或做其他计算。
   *
   */
  virtual State
  check(const unsigned int step, const double check_value);

  /**
   * 返回最后一次检查操作的结果。
   *
   */
  State
  last_check() const;

  /**
   * 返回初始收敛准则。
   *
   */
  double
  initial_value() const;

  /**
   * 返回求解器调用 @p check 的最后一个迭代步骤的收敛值。
   *
   */
  double
  last_value() const;

  /**
   * 最后一个迭代步骤的数量。
   *
   */
  unsigned int
  last_step() const;

  /**
   * 最大的步骤数。
   *
   */
  unsigned int
  max_steps() const;

  /**
   * 改变最大步数。
   *
   */
  unsigned int
  set_max_steps(const unsigned int);

  /**
   * 启用失败检查。如果<tt>residual>failure_residual</tt>与<tt>failure_residual
   * := rel_failure_residual*first_residual</tt>，求解将以 @p ReturnState
   * @p 失败停止。
   *
   */
  void
  set_failure_criterion(const double rel_failure_residual);

  /**
   * 禁用故障检查，并将 @p relative_failure_residual 和 @p
   * failure_residual重置为零。
   *
   */
  void
  clear_failure_criterion();

  /**
   * 宽容。
   *
   */
  double
  tolerance() const;

  /**
   * 改变容忍度。
   *
   */
  double
  set_tolerance(const double);

  /**
   * 使得每一步的残差写成一个向量，以便以后分析。
   *
   */
  void
  enable_history_data();

  /**
   * 提供对收集的残差数据的读取权限。
   *
   */
  const std::vector<double> &
  get_history_data() const;

  /**
   * 所有步骤的平均误差减少。    需要 enable_history_data()
   *
   */
  double
  average_reduction() const;
  /**
   * 最后一步的误差减少；对于静止的迭代，这近似于迭代矩阵的规范。
   * 需要 enable_history_data()
   *
   */
  double
  final_reduction() const;

  /**
   * 任何迭代步骤的误差减少。    需要 enable_history_data()
   *
   */
  double
  step_reduction(unsigned int step) const;

  /**
   * 记录每个迭代步骤。使用 @p log_frequency 来跳过步骤。
   *
   */
  void
  log_history(const bool);

  /**
   * 返回 @p log_history 标志。
   *
   */
  bool
  log_history() const;

  /**
   * 设置记录频率。
   *
   */
  unsigned int
  log_frequency(unsigned int);

  /**
   * 记录开始和结束的步骤。
   *
   */
  void
  log_result(const bool);

  /**
   * 返回 @p log_result 标志。
   *
   */
  bool
  log_result() const;

  /**
   * 如果一个对SolverControl对象的历史数据矢量进行操作的函数被调用，但历史数据的存储没有被enable_history_data()所启用，则会产生这个异常。
   *
   */
  DeclException0(ExcHistoryDataRequired);

protected:
  /**
   * 最大步骤数。
   *
   */
  unsigned int maxsteps;

  /**
   * 要达到的规定的公差。
   *
   */
  double tol;

  /**
   * 最后一次检查操作的结果。
   *
   */
  State lcheck;

  /**
   * 初始值。
   *
   */
  double initial_val;

  /**
   * 收敛标准的最后值。
   *
   */
  double lvalue;

  /**
   * 最后一步。
   *
   */
  unsigned int lstep;

  /**
   * 被 @p set_failure_criterion 设置为 @p true 并启用故障检查。
   *
   */
  bool check_failure;

  /**
   * 存储由 @p set_failure_criterion 设置的 @p rel_failure_residual 。
   *
   */
  double relative_failure_residual;

  /**
   * @p failure_residual 等于第一个残差乘以 @p 由 @p
   * set_failure_criterion 设置的relative_crit（见那里）。
   * 在知道第一个残差之前，它是0。
   *
   */
  double failure_residual;

  /**
   * 对数收敛历史到 @p deallog. 。
   *
   */
  bool m_log_history;

  /**
   * 只记录每第n步。
   *
   */
  unsigned int m_log_frequency;

  /**
   * 记录迭代结果到 @p deallog.
   * 如果是真的，在完成迭代后，关于失败或成功的声明以及
   * @p lstep 和 @p lvalue 被记录下来。
   *
   */
  bool m_log_result;

  /**
   * 对历史数据存储的控制。由enable_history_data()设置。
   *
   */
  bool history_data_enabled;

  /**
   * 存储每个迭代步骤后的结果的向量，用于以后的统计分析。
   * 这个向量的使用由enable_history_data()来启用。
   *
   */
  std::vector<double> history_data;
};


/**
 * @p SolverControl
 * 的特殊化，如果达到了指定的容差，或者初始残差（或解算器类选择的任何标准）减少了一个给定的系数，则返回
 * @p success
 * 。这在你不想精确求解，而是想获得两位数或达到最大迭代数的情况下很有用。比如说。最大的迭代次数是20次，减少系数是1%，公差是0.1%。最初的残差是2.5。如果完成了20次迭代，或者新的残差小于2.5*1%，或者小于0.1%，这个过程将中断。
 *
 *
 */
class ReductionControl : public SolverControl
{
public:
  /**
   * 构造函数。
   * 除了与SolverControl构造函数含义相同的参数外，还提供还原系数。
   *
   */
  ReductionControl(const unsigned int maxiter     = 100,
                   const double       tolerance   = 1.e-10,
                   const double       reduce      = 1.e-2,
                   const bool         log_history = false,
                   const bool         log_result  = true);

  /**
   * 用一个SolverControl对象进行初始化。其结果将通过将 @p
   * reduce 设置为零来模拟SolverControl。
   *
   */
  ReductionControl(const SolverControl &c);

  /**
   * 将一个SolverControl对象赋值给ReductionControl。赋值的结果将通过设置
   * @p reduce 为零来模拟SolverControl。
   *
   */
  ReductionControl &
  operator=(const SolverControl &c);

  /**
   * 由于该类中存在虚拟函数，所以需要虚拟析构器。
   *
   */
  virtual ~ReductionControl() override = default;

  /**
   * 到参数文件的接口。
   *
   */
  static void
  declare_parameters(ParameterHandler &param);

  /**
   * 从文件中读取参数。
   *
   */
  void
  parse_parameters(ParameterHandler &param);

  /**
   * 决定一个迭代的成功或失败。
   * 这个函数调用基类中的函数，但在第一次迭代时将容忍度设置为<tt>还原初始值</tt>。
   *
   */
  virtual State
  check(const unsigned int step, const double check_value) override;

  /**
   * 缩减系数。
   *
   */
  double
  reduction() const;

  /**
   * 改变减少系数。
   *
   */
  double
  set_reduction(const double);

protected:
  /**
   * 希望的减少系数。
   *
   */
  double reduce;

  /**
   * 减少的容忍度。如果达到这个值或者基类表示成功，则停止迭代。
   *
   */
  double reduced_tol;
};

/**
 * @p SolverControl
 * 的特殊化，如果执行了给定的迭代次数，则返回 @p success
 * ，而不考虑实际的残差。这在你不想精确求解，而是想执行固定次数的迭代的情况下很有用，例如在一个内部求解器中。给予该类的参数与SolverControl类的参数完全相同，当达到给定的容差或最大迭代次数之一时，求解器也会同样终止。与SolverControl的唯一区别是，在后一种情况下，求解器会返回成功。
 *
 *
 */
class IterationNumberControl : public SolverControl
{
public:
  /**
   * 构造函数。
   * 提供与SolverControl类的构造函数完全相同的参数。
   *
   */
  IterationNumberControl(const unsigned int maxiter     = 100,
                         const double       tolerance   = 1e-12,
                         const bool         log_history = false,
                         const bool         log_result  = true);

  /**
   * 用一个SolverControl对象进行初始化。结果将模仿SolverControl，将还原目标设置为零。
   *
   */
  IterationNumberControl(const SolverControl &c);

  /**
   * 将一个SolverControl对象分配给ReductionControl。赋值的结果将通过将还原目标设置为零来模拟SolverControl。
   *
   */
  IterationNumberControl &
  operator=(const SolverControl &c);

  /**
   * 需要虚拟析构器，因为这个类中有虚拟函数。
   *
   */
  virtual ~IterationNumberControl() override = default;

  /**
   * 决定一个迭代的成功或失败。这个函数将成功完全建立在是否达到了给定的迭代次数或者检查值正好为零的事实上。
   *
   */
  virtual State
  check(const unsigned int step, const double check_value) override;
};


/**
 * @p SolverControl
 * 的特化，当且仅当一定数量的连续迭代满足指定的公差时，返回
 * SolverControl::State::success
 * 。这在使用非精确Hessian解决非线性问题的情况下很有用。
 * 比如说。要求的连续收敛迭代次数是2，公差是0.2。ConsecutiveControl将只在序列0.5,
 * 0.0005, 1.0, 0.05, 0.01的最后一步返回 SolverControl::State::success
 * 。
 *
 *
 */
class ConsecutiveControl : public SolverControl
{
public:
  /**
   * 构建器。  @p n_consecutive_iterations
   * 是连续迭代的次数，应满足收敛的规定公差。其他参数的含义与SolverControl的构造函数相同。
   *
   */
  ConsecutiveControl(const unsigned int maxiter                  = 100,
                     const double       tolerance                = 1.e-10,
                     const unsigned int n_consecutive_iterations = 2,
                     const bool         log_history              = false,
                     const bool         log_result               = false);

  /**
   * 用一个SolverControl对象进行初始化。结果将通过设置 @p
   * n_consecutive_iterations 为1来模拟SolverControl。
   *
   */
  ConsecutiveControl(const SolverControl &c);

  /**
   * 将一个SolverControl对象分配给ConsecutiveControl。赋值的结果将通过设置
   * @p n_consecutive_iterations 为1来模拟SolverControl。
   *
   */
  ConsecutiveControl &
  operator=(const SolverControl &c);

  /**
   * 由于该类中存在虚拟函数，所以需要虚拟析构器。
   *
   */
  virtual ~ConsecutiveControl() override = default;

  /**
   * 决定一个迭代的成功或失败，见上面的类描述。
   *
   */
  virtual State
  check(const unsigned int step, const double check_value) override;

protected:
  /**
   * 连续迭代的次数，应满足规定的收敛容忍度。
   *
   */
  unsigned int n_consecutive_iterations;

  /**
   * 连续收敛的迭代次数的计数器。
   *
   */
  unsigned int n_converged_iterations;
};

 /*@}*/ 
//---------------------------------------------------------------------------

#ifndef DOXYGEN

inline unsigned int
SolverControl::max_steps() const
{
  return maxsteps;
}



inline unsigned int
SolverControl::set_max_steps(const unsigned int newval)
{
  unsigned int old = maxsteps;
  maxsteps         = newval;
  return old;
}



inline void
SolverControl::set_failure_criterion(const double rel_failure_residual)
{
  relative_failure_residual = rel_failure_residual;
  check_failure             = true;
}



inline void
SolverControl::clear_failure_criterion()
{
  relative_failure_residual = 0;
  failure_residual          = 0;
  check_failure             = false;
}



inline double
SolverControl::tolerance() const
{
  return tol;
}



inline double
SolverControl::set_tolerance(const double t)
{
  double old = tol;
  tol        = t;
  return old;
}


inline void
SolverControl::log_history(const bool newval)
{
  m_log_history = newval;
}



inline bool
SolverControl::log_history() const
{
  return m_log_history;
}


inline void
SolverControl::log_result(const bool newval)
{
  m_log_result = newval;
}


inline bool
SolverControl::log_result() const
{
  return m_log_result;
}


inline double
ReductionControl::reduction() const
{
  return reduce;
}


inline double
ReductionControl::set_reduction(const double t)
{
  double old = reduce;
  reduce     = t;
  return old;
}

#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif


