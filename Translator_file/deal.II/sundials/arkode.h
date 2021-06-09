//include/deal.II-translator/sundials/arkode_0.txt
//-----------------------------------------------------------
//
//    Copyright (C) 2017 - 2021 by the deal.II authors
//
//    This file is part of the deal.II library.
//
//    The deal.II library is free software; you can use it, redistribute
//    it, and/or modify it under the terms of the GNU Lesser General
//    Public License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//    The full text of the license can be found in the file LICENSE.md at
//    the top level directory of deal.II.
//
//---------------------------------------------------------------

#ifndef dealii_sundials_arkode_h
#define dealii_sundials_arkode_h

#include <deal.II/base/config.h>

#include <deal.II/base/mpi.h>

#ifdef DEAL_II_WITH_SUNDIALS

#  include <deal.II/base/conditional_ostream.h>
#  include <deal.II/base/exceptions.h>
#  include <deal.II/base/logstream.h>
#  include <deal.II/base/parameter_handler.h>
#  ifdef DEAL_II_WITH_PETSC
#    include <deal.II/lac/petsc_block_vector.h>
#    include <deal.II/lac/petsc_vector.h>
#  endif
#  include <deal.II/lac/vector.h>
#  include <deal.II/lac/vector_memory.h>

#  include <arkode/arkode.h>
#  if DEAL_II_SUNDIALS_VERSION_LT(4, 0, 0)
#    include <arkode/arkode_impl.h>
#  endif
#  include <nvector/nvector_serial.h>
#  ifdef DEAL_II_WITH_MPI
#    include <nvector/nvector_parallel.h>
#  endif
#  include <deal.II/base/discrete_time.h>

#  include <deal.II/sundials/n_vector.h>
#  include <deal.II/sundials/sunlinsol_wrapper.h>

#  include <boost/signals2.hpp>

#  include <sundials/sundials_linearsolver.h>
#  include <sundials/sundials_math.h>
#  include <sundials/sundials_types.h>

#  include <memory>


DEAL_II_NAMESPACE_OPEN


// Shorthand notation for ARKODE error codes.
#  define AssertARKode(code) Assert(code >= 0, ExcARKodeError(code))

/**
 * 一个命名空间，用于通过SUNDIALS包处理ODE求解器。
 *
 *
 */
namespace SUNDIALS
{
  /**
   * SUNDIALS加性Runge-Kutta方法（ARKode）的接口。
   * ARKode类是对SUNDIALS可变步长、嵌入式、加法Runge-Kutta求解器的封装，它是一个通用的求解器，用于求解具有快速和慢速动态特征的常微分方程系统。
   * 使用嵌套的隐式和显式Runge-Kutta求解器对快速动力学进行隐式处理，对慢速动力学进行显式处理。
   * 直接引自ARKode文档。    ARKode解决 $R^N$
   * 中的ODE初值问题（IVPs）。这些问题应以显式形式提出，如
   * \f[ M\dot y = f_E(t, y) + f_I (t, y), \qquad y(t_0) = y_0. \f] 这里，
   * $t$  是自变量（如时间），因变量由  $y \in R^N$
   * 给出，我们使用符号  $\dot y$  表示  $dy/dt$  。  $M$
   * 是一个用户提供的来自 $R^N \to R^N$ 的非星形算子。
   * 这个算子可能取决于 $t$ ，但不取决于 $y$  。
   * 对于标准的常微分方程系统和使用有限差分或有限体积方法对偏微分方程进行空间半微分所产生的问题，
   * $M$ 通常是身份矩阵， $I$
   * 。对于使用有限元空间半离散化的PDEs， $M$
   * 通常是一个条件良好的质量矩阵。
   * 两个右手边的函数可以描述为。
   *
   *
   *
   *
   *
   *
   * -  $f_E(t, y)$  : 包含系统的 "慢 "时间尺度成分。                 这将使用明确的方法进行整合。
   *
   *
   *
   *
   *
   * -  $f_I(t, y)$  ：包含系统的 "快速 "时间尺度成分。                 这将使用隐式方法进行整合。    ARKode可用于解决刚度、非刚度和多速率问题。  粗略地说，刚度的特点是至少存在一个快速阻尼模式，其时间常数与解本身的时间尺度相比很小。在上面的隐式/显式（ImEx）拆分中，这些刚度成分应该包括在右手边的函数中  $f_I (t, y)$  。    对于多速率问题，用户应提供定义IVP系统的两个函数  $f_E$  和  $f_I$  。    对于非刚性问题，只需提供 $f_E$ ，并假定 $f_I$ 为零，即系统还原为非分叉IVP：\f[
   * M\dot y = f_E(t, y), \qquad y(t_0) = y_0. \f]
   * 在这种情况下，ARK方法还原为经典的显式Runge-Kutta方法（ERK）。对于这些类别的方法，ARKode允许准确度等级
   * $q = \{2, 3, 4, 5, 6, 8\}$  ，嵌入等级  $p = \{1, 2, 3, 4, 5, 7\}$
   * 。这些方法默认为Heun-Euler-2-1-2, Bogacki-Shampine-4-2-3,
   * Zonneveld-5-3-4, Cash-Karp-6-4-5, Verner-8-5-6 和 Fehlberg-13-7-8。
   * 最后，对于刚性（线性或非线性）问题，用户可以只提供
   * $f_I$ ，这意味着 $f_E = 0$
   * ，这样系统就还原为非分叉IVP\f[ M\dot y = f_I(t, y), \qquad
   * y(t_0) = y_0. \f]
   * 与ERK方法类似，在这种情况下，ARK方法还原为经典的对角隐式Runge-Kutta方法（DIRK）。对于这些类别的方法，ARKode允许准确度等级
   * $q = \{2, 3, 4, 5\}$  ，嵌入等级  $p = \{1, 2, 3, 4\}$
   * 。这些默认为SDIRK-2-1-2、ARK-4-2-3（隐式）、SDIRK-5-3-4和ARK-8-4-5（隐式）方法，分别。
   * 对于DIRK和ARK方法，每个阶段都必须求解形式为\f[ G(z_i)
   * \dealcoloneq M z_i
   *
   * - h_n A^I_{i,i} f_I (t^I_{n,i}, z_i)
   *
   * - a_i = 0 \f]的隐式系统 $z_i , i = 1, \ldots, s$
   * ，其中对于ARK方法我们有数据\f[ a_i \dealcoloneq M y_{n-1} +
   * h_n \sum_{j=1}^{i-1} [ A^E_{i,j} f_E(t^E_{n,j}, z_j) + A^I_{i,j} f_I
   * (t^I_{n,j}, z_j)] \f]，对于DIRK方法有\f[ a_i \dealcoloneq M y_{n-1}
   * + h_n \sum_{j=1}^{i-1} A^I_{i,j} f_I (t^I_{n,j}, z_j) \f]。这里
   * $A^I_{i,j}$ 和 $A^E_{i,j}$ 是所选求解器的布彻表。    如果
   * $f_I(t,y)$ 非线性地依赖于 $y$
   * ，那么上面的系统对应于一个非线性方程组；如果 $f_I
   * (t, y)$ 线性地依赖于 $y$
   * ，那么这就是一个线性方程组。通过指定标志`implicit_function_is_linear`，ARKode走了一些捷径，使求解过程更快。
   * 对于这两种类型的系统，ARKode允许选择解决策略。
   * 默认的求解器选择是牛顿方法的变种，\f[ z_i^{m+1} = z_i^m
   * +\delta^{m+1}, \f] 其中 $m$ 是牛顿指数，牛顿更新
   * $\delta^{m+1}$ 需要解决线性牛顿系统\f[ N(z_i^m) \delta^{m+1} =
   *
   * -G(z_i^m), \f] 其中\f[ N \dealcoloneq M
   *
   * - \gamma J, \quad J \dealcoloneq \frac{\partial f_I}{\partial y}, \qquad
   * \gamma\dealcoloneq h_n A^I_{i,i}. \f]
   * 作为牛顿方法的替代方法，ARKode可以对每个阶段 $z_i ,i =
   * 1, \ldots , s$  ] 使用安德森加速固定点迭代\f[ z_i^{m+1} =
   * g(z_i^{m}), m=0,1,\ldots. \f]
   * 与牛顿方法不同，这个选项不需要在每个迭代中求解线性系统，而是选择求解低维最小二乘法来构建非线性更新。
   * 最后，如果用户指定 "implicit_function_is_linear"，即 $f_I(t,
   * y)$ 线性地依赖于 $y$
   * ，并且选择基于牛顿的非线性求解器，那么系统将只使用一次牛顿迭代来求解。请注意，为了使用牛顿求解器，至少应该提供jacobian_times_vector()函数（或者SUNDIALS版本>4.0.0的solve_jacobian_system()）。如果不提供这个函数，则只支持定点迭代，而
   * "implicit_function_is_linear "的设置将被忽略。
   * 最佳解算器（牛顿与定点）在很大程度上取决于问题。
   * 由于定点求解器不需要求解任何线性系统，因此每次迭代的成本可能比牛顿求解器低很多。但是，与牛顿方法相比，这可能是以收敛速度较慢（甚至是发散）为代价的。这些定点求解器确实允许用户指定安德森加速子空间的大小，
   * $m_k$  。虽然所需的求解器内存量与 $m_k N$
   * 成正比增长，但 $m_k$ 的数值越大，收敛速度越快。
   * 即使对于 "小 "值，例如 $1 \leq m_k \leq 5$
   * ，这种改进也可能是显著的，而对于较大的 $m_k$
   * 值，收敛性可能不会改善（甚至恶化）。ARKode使用基于牛顿的迭代作为其默认的求解器，因为它对非常僵硬的问题具有更强的鲁棒性，但强烈建议用户在尝试新问题时也考虑使用定点求解器。
   * 无论是牛顿求解还是定点求解，众所周知，算法的效率和稳健性都与选择一个好的初始猜测密切相关。在ARKode中，任何一种非线性求解方法的初始猜测都是一个预测值
   * $z_i(0)$
   * ，这个预测值是根据之前计算的数据明确计算出来的（例如
   * $y_{n-2}, y_{n-1}$  ，以及 $z_j$ ，其中 $j < i$
   * ）。关于ARKode中实现的具体预测算法的其他信息，在ARKode文档中提供。
   * 用户必须提供以下至少一种（或两种）的实现
   * `std::function`s:  。
   *
   *
   *
   *
   *
   *
   * - implicit_function()
   *
   *
   *
   *
   *
   *
   * - explicit_function() 如果质量矩阵与身份矩阵不同，用户应提供
   *
   *
   *
   *
   *
   * - mass_times_vector() (或solve_mass_system() for SUNDIALS version < 4.0.0) and, optionally,
   *
   *
   *
   *
   *
   * - mass_times_setup() (或setup_mass() 用于SUNDIALS版本< 4.0.0) 如果需要使用牛顿方法，那么用户还应该提供
   *
   *
   *
   *
   *
   * - jacobian_times_vector (或solve_jacobian_system() for SUNDIALS version < 4.0.0)
   *
   *
   *
   *
   *
   * - 可选：jacobian_times_setup() (或者SUNDIALS版本< 4.0.0的setup_jacobian() )
   * @note
   * 尽管SUNDIALS可以提供雅各布系数的差商近似值，但目前不支持通过这个包装器来实现。
   * 仅适用于SUNDIALS版本>4.0.0:
   * SUNDIALS默认的求解器（SPGMR）被用来解决线性系统。如果要对质量矩阵和/或雅各布系数使用自定义的线性求解器，请设置。
   *
   *
   *
   *
   *
   * - solve_mass()和/或
   *
   *
   *
   *
   *
   * - solve_jacobian() 仅适用于SUNDIALS版本>4.0.0: 要使用自定义预处理程序与默认或自定义线性求解器，请设置。
   *
   *
   *
   *
   * - jacobian_preconditioner_solve()和/或mass_preconditioner_solve()，以及，可选择。
   *
   *
   * - jacobian_preconditioner_setup() 和/或 mass_preconditioner_setup() 也可以重写以下函数。默认情况下，它们什么都不做，或者不需要。
   *
   *
   *
   *
   *
   *
   * - solver_should_restart()
   *
   *
   *
   *
   * - get_local_tolerances() 要在固定的步长下产生输出，设置函数
   *
   *
   *
   *
   *
   * - output_step() ARKODE对象的任何其他自定义设置都可以在下面指定
   *
   *
   *
   *
   *
   *
   * - custom_setup() 为了提供一个简单的例子，考虑谐波振荡器问题： \f[
   * \begin{split} u'' & =
   *
   * -k^2 u \\ u (0) & = 0 \\ u'(0) & = k \end{split} \f]
   * 我们用一阶奥德来写它。 \f[ \begin{matrix} y_0' & =  y_1 \\
   * y_1' & =
   *
   * - k^2 y_0 \end{matrix} \f] 也就是 $y' = A y$  其中\f[ A \dealcoloneq
   * \begin{matrix} 0 & 1 \\
   *
   *
   *
   *
   *
   *
   *
   * -k^2 &0 \end{matrix} \f] 和  $y(0)=(0, k)$  。    精确解是  $y_0(t)
   * = \sin(k t)$  ,  $y_1(t) = y_0'(t) = k \cos(k t)$  ,  $y_1'(t) =
   *
   * -k^2 \sin(k t)$  。
   * 一个最小的实现，只使用显式RK方法，由以下代码片断给出。
   * @code
   * using VectorType = Vector<double>;
   *
   * SUNDIALS::ARKode<VectorType> ode;
   *
   * const double kappa = 1.0;
   *
   * ode.explicit_function = [kappa] (double,
   *                                const VectorType &y,
   *                                VectorType &ydot)
   *
   * -> int
   * {
   * ydot[0] = y[1];
   * ydot[1] =
   *
   * -kappa*kappa*y[0];
   * return 0;
   * };
   *
   * Vector<double> y(2);
   * y[1] = kappa;
   *
   * ode.solve_ode(y);
   * @endcode
   *
   */
  template <typename VectorType = Vector<double>>
  class ARKode
  {
  public:
    /**
     * 可以传递给ARKode类的额外参数。
     *
     */
    class AdditionalData
    {
    public:
      /**
       * ARKode的初始化参数。            全局参数。
       * @param  initial_time 初始时间  @param  final_time 最终时间
       * @param  initial_step_size 初始步长  @param  output_period
       * 每次输出之间的期望时间间隔 运行参数。
       * @param  minimum_step_size 最小步长  @param  maximum_order
       * ARK最大顺序  @param  maximum_non_linear_iterations
       * 非线性迭代的最大次数  @param  implicit_function_is_linear
       * 指定问题的隐含部分为线性  @param  ]
       * implicit_function_is_time_independent
       * 指定问题的隐含部分是线性的，与时间无关  @param
       * mass_is_time_independent 指定质量预因子与时间无关  @param
       * anderson_acceleration_subspace
       * 在打包的SUNDIALS求解器中用于安德森加速度的矢量数量。
       * 错误参数。              @param  absolute_tolerance
       * 绝对误差公差  @param  relative_tolerance 相对误差公差
       *
       */
      AdditionalData(
        // Initial parameters
        const double initial_time      = 0.0,
        const double final_time        = 1.0,
        const double initial_step_size = 1e-2,
        const double output_period     = 1e-1,
        // Running parameters
        const double       minimum_step_size                     = 1e-6,
        const unsigned int maximum_order                         = 5,
        const unsigned int maximum_non_linear_iterations         = 10,
        const bool         implicit_function_is_linear           = false,
        const bool         implicit_function_is_time_independent = false,
        const bool         mass_is_time_independent              = false,
        const int          anderson_acceleration_subspace        = 3,
        // Error parameters
        const double absolute_tolerance = 1e-6,
        const double relative_tolerance = 1e-5)
        : initial_time(initial_time)
        , final_time(final_time)
        , initial_step_size(initial_step_size)
        , minimum_step_size(minimum_step_size)
        , absolute_tolerance(absolute_tolerance)
        , relative_tolerance(relative_tolerance)
        , maximum_order(maximum_order)
        , output_period(output_period)
        , maximum_non_linear_iterations(maximum_non_linear_iterations)
        , implicit_function_is_linear(implicit_function_is_linear)
        , implicit_function_is_time_independent(
            implicit_function_is_time_independent)
        , mass_is_time_independent(mass_is_time_independent)
        , anderson_acceleration_subspace(anderson_acceleration_subspace)
      {}

      /**
       * 将所有AdditionalData()参数添加到给定的ParameterHandler对象中。当参数从文件中被解析出来时，内部参数会自动更新。
       * 你在构造时传递的选项被设置为ParameterHandler对象`prm`中的默认值。以后你可以通过使用`prm`解析参数文件来修改它们。每当`prm`的内容被更新时，参数的值就会被更新。
       * 请确保这个类的寿命比`prm`长。如果你破坏了这个类，然后用`prm`解析一个参数文件，将会发生未定义的行为。
       *
       */
      void
      add_parameters(ParameterHandler &prm)
      {
        prm.add_parameter("Initial time", initial_time);
        prm.add_parameter("Final time", final_time);
        prm.add_parameter("Time interval between each output", output_period);

        prm.enter_subsection("Running parameters");
        prm.add_parameter("Initial step size", initial_step_size);
        prm.add_parameter("Minimum step size", minimum_step_size);
        prm.add_parameter("Maximum order of ARK", maximum_order);
        prm.add_parameter("Maximum number of nonlinear iterations",
                          maximum_non_linear_iterations);
        prm.add_parameter("Implicit function is linear",
                          implicit_function_is_linear);
        prm.add_parameter("Implicit function is time independent",
                          implicit_function_is_time_independent);
        prm.add_parameter("Mass is time independent", mass_is_time_independent);
        prm.add_parameter("Anderson-acceleration subspace",
                          anderson_acceleration_subspace);
        prm.leave_subsection();

        prm.enter_subsection("Error control");
        prm.add_parameter("Absolute error tolerance", absolute_tolerance);
        prm.add_parameter("Relative error tolerance", relative_tolerance);
        prm.leave_subsection();
      }

      /**
       * DAE的初始时间。
       *
       */
      double initial_time;

      /**
       * 最终时间。
       *
       */
      double final_time;

      /**
       * 初始步骤大小。
       *
       */
      double initial_step_size;

      /**
       * 最小步长。
       *
       */
      double minimum_step_size;

      /**
       * 自适应时间步进的绝对误差容限。
       *
       */
      double absolute_tolerance;

      /**
       * 自适应时间步进的相对误差容限。
       *
       */
      double relative_tolerance;

      /**
       * ARK的最大顺序。
       *
       */
      unsigned int maximum_order;

      /**
       * 每个输出之间的期望时间周期。实际输出的时间周期可由Arkode调整。
       *
       */
      double output_period;

      /**
       * 在时间推进过程中，牛顿法或定点法的最大迭代次数。
       *
       */
      unsigned int maximum_non_linear_iterations;

      /**
       * 指定问题的隐含部分是否为线性。
       *
       */
      bool implicit_function_is_linear;

      /**
       * 指定问题的隐含部分是否是线性的和时间独立的。
       *
       */
      bool implicit_function_is_time_independent;

      /**
       * 指定质量预因子是否与时间无关。如果没有指定质量，则没有影响。
       *
       */
      bool mass_is_time_independent;

      /**
       * 用于安德森加速度的子空间向量的数量。只有在使用打包的SUNDIALS定点求解器时才有意义。
       *
       */
      int anderson_acceleration_subspace;
    };

    /**
     * 构造器。通过传递一个设置所有求解器参数的AdditionalData()对象，可以对SUNDIALS
     * ARKode求解器进行微调。
     * 在串行情况下，MPI通信器被简单地忽略了。
     * @param  data ARKode配置数据  @param  mpi_comm MPI通信器
     *
     */
    ARKode(const AdditionalData &data     = AdditionalData(),
           const MPI_Comm &      mpi_comm = MPI_COMM_WORLD);

    /**
     * 解构器。
     *
     */
    ~ARKode();

    /**
     * 对初始值问题进行积分。该函数返回最终的计算步骤数。
     * @param  解决方案
     * 在输入时，这个向量包含初始条件。在输出时，它包含最后时间的解决方案。
     *
     */
    unsigned int
    solve_ode(VectorType &solution);

    /**
     * 对初始值问题进行积分。与上面的函数相比，这个函数允许为下一个解指定一个
     * @p intermediate_time 。
     * 这个函数的重复调用必须使用单调增加的 @p
     * intermediate_time. 值，最后一个解的状态与 @p
     * intermediate_time
     * 一起被内部保存，并将作为初始条件在下一次调用中重新使用。
     * 用户可能会发现这个函数在将ARKode集成到他们自己的外部时间循环中时非常有用，特别是当output_step()的限制性太强时。
     * @note   @p intermediate_time 可能大于 AdditionalData::final_time,
     * ，这个函数会忽略。          @param  解决方案
     * 最终解决方案。如果求解器重新启动，无论是因为它是第一次求解还是标志
     * @p reset_solver 被设置，该向量也被用作初始条件。
     * @param  intermediate_time
     * 递增解步骤的时间。必须大于先前调用此函数时使用的最后时间。
     * @param  reset_solver
     * 可选的标志，用于重新创建所有的内部对象，这对于空间适应性方法可能是可取的。如果设置为
     * "true"，在求解ODE之前调用reset()，它将 @p solution
     * 设置为初始条件。这不会*重置以前调用此函数的存储时间。
     *
     */
    unsigned int
    solve_ode_incrementally(VectorType & solution,
                            const double intermediate_time,
                            const bool   reset_solver = false);

    /**
     * 清理内部存储器，以干净的对象开始。当仿真开始时，以及当用户对solver_should_restart()的调用返回true时，该函数被调用。
     * 默认情况下，solver_should_restart()返回false。如果用户需要实现例如空间的局部适应性，可以给solver_should_restart()指定一个不同的函数，执行所有的网格变化，将解转移到新的网格，并返回true。
     * @param  t 新的起始时间  @param  h 新的起始时间步骤
     * @param  y 新的初始解
     *
     */
    void
    reset(const double t, const double h, const VectorType &y);

    /**
     * 提供用户对内部使用的ARKODE内存的访问。
     * 这个功能是为那些希望直接从ARKODE集成器中查询额外信息的用户准备的，关于各种`ARKStepGet...`函数，请参考ARKODE手册。不应调用`ARKStepSet...`函数，因为这可能导致与该ARKode对象执行的各种设置发生冲突。
     * @note
     * 如果需要对ARKODE功能进行自定义设置（无法通过该类的接口实现），应该使用函数custom_setup()。
     * @return
     * 指向ARKODE内存块的指针，可以传递给SUNDIALS函数。
     *
     */
    void *
    get_arkode_memory() const;

    /**
     * 一个用于 "重新启动
     * "给定向量的函数对象。设置这个字段不再有任何作用，所有的辅助向量都会根据solve_ode()中用户提供的向量自动重新启动。
     * @deprecated
     * 这个函数不再使用，可以在用户代码中安全地删除。
     *
     */
    DEAL_II_DEPRECATED
    std::function<void(VectorType &)> reinit_vector;

    /**
     * 一个用户可以提供的函数对象，其目的是计算IVP右手边的显式部分。设置
     * $explicit_f = f_E(t, y)$  。        必须提供 explicit_function()
     * 或 implicit_function()
     * 中的至少一个。根据所提供的函数，显式、隐式或混合RK方法被使用。
     * 这个函数应该返回。
     *
     *
     *
     *
     *
     * - 0: 成功
     *
     *
     *
     *
     * - >0: 可恢复的错误，ARKode将重新尝试解决方案并再次调用此函数。
     *
     *
     *
     *
     *
     *
     * - <0: 无法恢复的错误，计算将被中止，并抛出一个断言。
     *
     */
    std::function<
      int(const double t, const VectorType &y, VectorType &explicit_f)>
      explicit_function;

    /**
     * 一个用户可以提供的函数对象，其目的是计算IVP右手边的隐含部分。设置
     * $implicit_f = f_I(t, y)$  。        必须提供 explicit_function()
     * 或 implicit_function()
     * 中的至少一个。根据所提供的函数，显式、隐式或混合RK方法被使用。
     * 这个函数应该返回。
     *
     *
     *
     *
     *
     * - 0: 成功
     *
     *
     *
     *
     * - >0: 可恢复的错误，ARKode将重新尝试解决方案并再次调用此函数。
     *
     *
     *
     *
     *
     *
     * - <0: 无法恢复的错误，计算将被中止，并抛出一个断言。
     *
     */
    std::function<int(const double t, const VectorType &y, VectorType &res)>
      implicit_function;

#  if DEAL_II_SUNDIALS_VERSION_LT(4, 0, 0)
    /**
     * 一个用户可以提供的函数对象，其目的是为后续调用solve_jacobian_system()准备线性求解器。
     * 确保在调用此函数后，我们知道如何计算系统的解  $A
     * x = b$  ，其中  $A$  是牛顿矩阵的一些近似值，  $M
     *
     * - \gamma \partial f_I/\partial y$
     * 。这个函数是可选的。如果用户不提供，那么solve_jacobian_system()被认为也会在内部进行设置。
     * setup_jacobian()函数可以调用一个用户提供的函数来计算与雅各布矩阵相关的所需数据。或者，它也可以选择检索和使用这些数据的存储值。在这两种情况下，setup_jacobian()也可以根据solve_jacobian_system()的需要对这些数据进行预处理，这可能涉及到调用一个通用函数（比如用于LU因子化）。
     * 这些数据既可以直接使用（在直接线性求解器中），也可以用于预处理器（在预处理的迭代线性求解器中）。setup_jacobian()函数不是在每个求解阶段（甚至是每个时间步长）都被调用，而只是在求解器确定适合执行设置任务时才会频繁调用。这样一来，由setup_jacobian()产生的雅各布相关数据就有望在若干时间步长中被使用。
     * 如果用户使用基于矩阵的雅可比计算，那么这里就应该调用一个装配例程来装配矩阵和雅可比系统的前置条件器。
     * 随后调用（可能不止一次）solve_jacobian_system()可以假设这个函数至少被调用过一次。
     * 注意，这个接口没有假设用户在这个函数中应该做什么。ARKode只假设在调用setup_jacobian()后，有可能调用solve_jacobian_system()，以获得系统
     * $x$ 的解  $J x = b$
     * 。如果没有提供这个函数，那么它就永远不会被调用。
     * 该函数的参数是  @param[in]  t 当前时间  @param[in]  gamma
     * 当前用于雅各布计算的系数  @param[in]  ypred 是当前
     * ARKode 内部步骤的预测  $y$  向量  @param[in]  fpred 是 ypred
     * 处隐含的右手边值，  $f_I (t_n, ypred)$  。
     * @param[in]  convfail
     * 输入标志，用于指示在当前时间步长的非线性方程求解过程中发生的任何问题，线性求解器正在被使用。这个标志可以用来帮助决定线性求解器所保存的雅各布数据是否需要更新。其可能的值是。
     *
     *
     *
     *
     *
     *
     * - ARK_NO_FAILURES：如果这是该步骤的第一次调用，或者在该步骤的前一次尝试中本地错误测试失败（但牛顿迭代收敛了），则该值被传递。
     *
     *
     *
     *
     *
     *
     * - ARK_FAIL_BAD_J：如果(a)之前的牛顿修正器迭代没有收敛，并且线性求解器的设置函数显示其雅各布相关数据不是当前的，或者(b)在之前的牛顿修正器迭代中，线性求解器的求解函数以可恢复的方式失败，并且线性求解器的设置函数显示其雅各布相关数据不是当前的，则传递此值。
     *
     *
     *
     *
     *
     *
     * - ARK_FAIL_OTHER：如果在当前的内部步骤尝试中，即使线性求解器使用的是当前的雅各布相关数据，但之前的牛顿迭代未能收敛，则传递此值。          @param[out]  j_is_current：一个布尔值，由setup_jacobian()填入。    如果调用后雅各布数据是当前的，该值应设置为`true`，如果其雅各布数据不是当前的，应设置为`false`。如果setup_jacobian()调用重新评估雅各布数据（基于convfail和ARKode状态数据），那么它应该无条件地将`j_is_current`设置为`true`，否则会导致无限循环。        这个函数应该返回。
     *
     *
     *
     *
     *
     * - 0: 成功
     *
     *
     *
     *
     * - >0: 可恢复的错误（如果发生这种情况，ARKodeReinit将被调用，然后最后一个函数将被再次尝试
     *
     *
     *
     *
     *
     *
     * - <0: 无法恢复的错误，计算将被中止，并抛出一个断言。
     *
     */
    std::function<int(const int         convfail,
                      const double      t,
                      const double      gamma,
                      const VectorType &ypred,
                      const VectorType &fpred,
                      bool &            j_is_current)>
      setup_jacobian;

    /**
     * 一个用户可以提供的函数对象，其目的是为了解决雅各布线性系统。在setup_jacobian()被调用至少一次后，这个函数将被ARKode调用（可能是多次）。ARKode试图尽最大努力调用setup_jacobian()最少的次数。如果不更新雅各布式就能实现收敛，那么ARKode就不会再次调用setup_jacobian()。相反，如果ARKode内部收敛测试失败，那么ARKode会用更新的向量和系数再次调用setup_jacobian()，这样连续调用solve_jacobian_systems()会导致牛顿过程中更好的收敛。
     * 如果你没有指定solve_jacobian_system()函数，那么就会使用固定点迭代而不是牛顿方法。请注意，这可能不会收敛，或者收敛得很慢。
     * jacobian $J$ 应该是在`t`，`ycur`评估的系统Jacobian\f[ J = M
     *
     * - \gamma \frac{\partial f_I}{\partial y}
     * \f]的（近似）。`fcur`是  $f_I(t,ycur)$  。
     * 对该函数的调用应该在`dst`中存储应用于`src`的 $J^{-1}$
     * 结果，即`J*dst =
     * src`。用户有责任在这个函数中设置适当的求解器和预处理器。
     * 该函数的参数是  @param[in]  t 当前时间  @param[in]  gamma
     * 当前用于雅各布计算的系数  @param[in]  ycur 是当前 ARKode
     * 内部步骤的  $y$  向量  @param[in]  fcur 是 ycur
     * 处隐式右手的当前值，  $f_I (t_n, ypred)$  。
     * 这个函数应该返回。
     *
     *
     *
     *
     *
     * - 0: 成功
     *
     *
     *
     *
     * - >0: 可恢复的错误（如果发生这种情况，ARKodeReinit将被调用，然后最后一个函数将被再次尝试
     *
     *
     *
     *
     *
     *
     * - <0: 无法恢复的错误，计算将被中止，并抛出一个断言。
     *
     */
    std::function<int(const double      t,
                      const double      gamma,
                      const VectorType &ycur,
                      const VectorType &fcur,
                      const VectorType &rhs,
                      VectorType &      dst)>
      solve_jacobian_system;


    /**
     * 一个用户可以提供的函数对象，其目的是设置质量矩阵。ARKode在任何需要更新质量矩阵的时候都会调用这个函数。用户应该计算质量矩阵（或更新所有允许应用质量矩阵的变量）。
     * 这个函数被ARKode调用一次，在任何调用solve_mass_system()之前。
     * ARKode支持质量矩阵可能取决于时间的情况，但不支持质量矩阵取决于解决方案本身的情况。
     * 如果用户没有提供solve_mass_matrix()函数，那么就会使用同一性。如果没有提供setup_mass()函数，那么solve_mass_system()应该自己完成所有工作。
     * 如果用户使用基于矩阵的质量矩阵计算，那么这里就应该调用一个装配程序来装配矩阵和预处理程序。随后调用（可能不止一次）solve_mass_system()可以假设这个函数至少被调用过一次。
     * 请注意，这个接口没有假设用户在这个函数中应该做什么。ARKode只假设在调用setup_mass()后，有可能调用solve_mass_system()，以获得系统
     * $x$ 的解  $M x = b$  。        这个函数应该返回。
     *
     *
     *
     *
     *
     * - 0: 成功
     *
     *
     *
     *
     * - >0: 可恢复的错误（如果发生这种情况，ARKodeReinit将被调用，然后最后一个函数将被再次尝试
     *
     *
     *
     *
     *
     *
     * - <0: 无法恢复的错误，计算将被中止，并抛出一个断言。
     *
     */
    std::function<int(const double t)> setup_mass;

    /**
     * 一个用户可以提供的函数对象，其目的是解决质量矩阵线性系统。在setup_mass()被调用至少一次后，这个函数将被ARKode调用（可能是多次）。ARKode试图尽最大努力调用setup_mass()最少的次数。
     * 对该函数的调用应在`dst`中存储 $M^{-1}$
     * 应用于`src`的结果，即`M*dst =
     * src`。用户有责任在这个函数中设置适当的求解器和预处理器。
     * 这个函数应该返回。
     *
     *
     *
     *
     *
     * - 0: 成功
     *
     *
     *
     *
     * - >0: 可恢复的错误（如果发生这种情况，ARKodeReinit将被调用，然后最后一个函数将被再次尝试
     *
     *
     *
     *
     *
     *
     * - <0: 无法恢复的错误，计算将被中止，并抛出一个断言。
     *
     */
    std::function<int(const VectorType &rhs, VectorType &dst)>
      solve_mass_system;
#  else

    /**
     * 一个用户可以提供的函数对象，其目的是计算质量矩阵与给定矢量`v`的乘积。在mass_times_setup()至少被调用一次后，这个函数将被ARKode调用（可能是多次）。ARKode试图尽最大努力调用mass_times_setup()最少的次数。
     * 对这个函数的调用应该在`Mv`中存储 $M$
     * 应用于`v`的结果。        这个函数应该返回。
     *
     *
     *
     *
     *
     * - 0: 成功
     *
     *
     *
     *
     * - >0: 可恢复的错误，ARKode将重新尝试解决方案并再次调用此函数。
     *
     *
     *
     *
     *
     *
     * - <0: 无法恢复的错误，计算将被中止，并抛出一个断言。
     *
     */
    std::function<int(double t, const VectorType &v, VectorType &Mv)>
      mass_times_vector;

    /**
     * 一个用户可以提供的函数对象，其目的是设置质量矩阵。ARKode在任何需要更新质量矩阵的时候都会调用这个函数。用户应该计算质量矩阵（或者更新所有允许应用质量矩阵的变量）。
     * 这个函数被保证至少被ARKode调用一次，在任何调用mass_times_vector()之前。
     * ARKode支持质量矩阵可能取决于时间的情况，但不支持质量矩阵取决于解决方案本身的情况。
     * 如果用户没有提供mass_times_vector()函数，那么就会使用同一性。如果没有提供mass_times_setup()函数，那么mass_times_vector()应该自己完成所有工作。
     * 如果用户使用基于矩阵的质量矩阵计算，那么这就是应该调用汇编例程来汇编矩阵的正确位置。对mass_times_vector()的后续调用（可能不止一次）可以假定这个函数至少被调用过一次。
     * @note
     * 这个接口没有假设用户在这个函数中应该做什么。ARKode只假设在调用mass_times_setup()后，有可能调用mass_times_vector()。
     * @param  t 当前的评估时间 这个函数应该返回。
     *
     *
     *
     *
     *
     * - 0: 成功
     *
     *
     *
     *
     * - >0: 可恢复的错误，ARKode将重新尝试解决方案并再次调用此函数。
     *
     *
     *
     *
     *
     *
     * - <0: 无法恢复的错误，计算将被中止，并抛出一个断言。
     *
     */
    std::function<int(const double t)> mass_times_setup;

    /**
     * 一个用户可以提供的函数对象，其目的是计算雅各布矩阵与给定矢量`v`的乘积。这里的Jacobian指的是
     * $J=\frac{\partial f_I}{\partial y}$
     * ，即用户指定的implicit_function的Jacobian。
     * 对这个函数的调用应该在`Jv`中存储 $J$
     * 应用于`v`的结果。        该函数的参数是  @param[in]  v
     * 要与Jacobian相乘的向量  @param[out]  Jv
     * 要用乘积J*v填充的向量  @param[in]  t 当前的时间
     * @param[in]  y 当前ARKode内部步骤的当前  $y$  向量
     * @param[in]  fy 在y的隐式右手边的当前值，  $f_I (t_n, y)$
     * 。        这个函数应该返回。
     *
     *
     *
     *
     *
     * - 0: 成功
     *
     *
     *
     *
     * - >0: 可恢复的错误，ARKode将重新尝试解决方案并再次调用此函数。
     *
     *
     *
     *
     *
     *
     * - <0: 无法恢复的错误，计算将被中止，并抛出一个断言。
     *
     */
    std::function<int(const VectorType &v,
                      VectorType &      Jv,
                      double            t,
                      const VectorType &y,
                      const VectorType &fy)>
      jacobian_times_vector;

    /**
     * 一个用户可以提供的函数对象，其目的是为jacobian_times_vector()的应用设置所有必要的数据。ARKode在任何需要更新雅各布的时候都会调用这个函数。
     * 用户应该计算Jacobian（或者更新所有允许应用Jacobian的变量）。这个函数保证至少被ARKode调用一次，在任何调用jacobian_times_vector()之前。
     * 如果没有提供jacobian_times_setup()函数，那么jacobian_times_vector()应该自己完成所有工作。
     * 如果用户使用基于矩阵的雅各布式计算，那么这里就应该调用一个汇编例程来组装矩阵。随后对jacobian_times_vector()的调用（可能不止一次）可以假定这个函数至少被调用过一次。
     * @note
     * 这个接口没有假设用户在这个函数中应该做什么。ARKode只假设在调用jacobian_times_setup()后，有可能调用jacobian_times_vector()。
     * @param  t 当前的时间  @param  y 当前的ARKode内部求解向量
     * $y$   @param  fy 在当前时间  $t$  和状态  $y$
     * 评估的隐式右手函数，即  $f_I(y,t)$
     * 这个函数应该返回。
     *
     *
     *
     *
     *
     * - 0: 成功
     *
     *
     *
     *
     * - >0: 可恢复的错误，ARKode将重新尝试解决方案并再次调用此函数。
     *
     *
     *
     *
     *
     *
     * - <0: 无法恢复的错误，计算将被中止，并抛出一个断言。
     *
     */
    std::function<int(realtype t, const VectorType &y, const VectorType &fy)>
      jacobian_times_setup;

    /**
     * 一个用户可以提供的LinearSolveFunction对象，其目的是为了解决线性化系统
     * $Ax=b$ ，其中 $A = M-\gamma J$
     * 是非线性残差的雅各布系数。质量矩阵 $M$ 和雅各布
     * $J$
     * 的应用通过函数mass_times_vector()和jacobian_times_vector()知道，而
     * $\gamma$ 是SUNDIALS提供的一个因子。矩阵-向量乘积 $Ax$
     * 在提供的SundialsOperator中被编码。如果通过jacobian_preconditioner_solve()设置了一个预处理程序，它被编码为SundialsPreconditioner。如果没有以这种方式提供预处理，则预处理是身份矩阵，即没有预处理。用户可以在这个函数对象中自由地使用一个没有通过SUNDIALS提供的自定义预处理。
     * 如果你没有指定solve_linearized_system()函数，那么就会使用SUNDIALS打包的SPGMR求解器，并采用默认设置。
     * 关于函数类型的更多细节请参考LinearSolveFunction。
     *
     */
    LinearSolveFunction<VectorType> solve_linearized_system;

    /**
     * 一个用户可以提供的LinearSolveFunction对象，其目的是为了解决质量系统
     * $Mx=b$  。矩阵-向量乘积  $Mx$
     * 在提供的SundialsOperator中被编码。如果通过mass_preconditioner_solve()设置了一个预处理程序，它被编码在SundialsPreconditioner中。如果没有以这种方式提供预调节器，预调节器就是身份矩阵，也就是说，没有预调节器。用户可以在这个函数对象中自由地使用一个没有通过SUNDIALS提供的自定义预处理。
     * 如果在mass_times_vector()中使用和应用非同一质量矩阵，用户必须指定这个函数。
     * 关于该函数类型的更多细节，请参考LinearSolveFunction。
     *
     */
    LinearSolveFunction<VectorType> solve_mass;


    /**
     * 用户可以提供一个函数对象，用于向SUNDIALS内置的求解器传递一个预处理程序，或者在用户自己的线性求解中应用一个自定义的预处理程序，该程序在solve_linearized_system()中指定。
     * 这个函数应该计算预处理方程 $Pz=r$
     * 的解，并将其存储在 @p z. 中。在这个方程 $P$
     * 中应该近似于非线性系统的雅各布 $M-\gamma J$ 。
     * @param[in]  t 当前的时间  @param[in]  y
     * 当前ARKode内部步骤的 $y$  矢量  @param[in]  fy
     * 在y处的隐式右手边的当前值，  $f_I (t_n, y)$  。
     * @param[in]  r 先决条件方程的右手边  @param[out]  z
     * 应用先决条件的解决方案，即解决  $Pz=r$   @param[in]
     * gamma 先决条件方程的值  $\gamma$   @param[in]  ] tol
     * 系统应该被解决的公差  @param[in]  lr
     * 一个输入标志，表示预处理程序的求解是使用左预处理程序（lr
     * = 1）还是右预处理程序（lr =
     * 2）。只有在与SUNDIALS打包的求解器一起使用时才相关。如果与自定义的solven_mass()函数一起使用，这个函数将被设置为零。
     * 这个函数应该返回。
     *
     *
     *
     *
     *
     * - 0: 成功
     *
     *
     *
     *
     * - >0: 可恢复的错误，ARKode将重新尝试解决方案并再次调用此函数。
     *
     *
     *
     *
     *
     *
     * - <0: 无法恢复的错误，计算将被中止，并抛出一个断言。
     *
     */
    std::function<int(double            t,
                      const VectorType &y,
                      const VectorType &fy,
                      const VectorType &r,
                      VectorType &      z,
                      double            gamma,
                      double            tol,
                      int               lr)>
      jacobian_preconditioner_solve;

    /**
     * 一个用户可以提供的函数对象，用于设置jacobian_preconditioner_solve()中指定的预处理程序。
     * 这个函数应该准备好预处理方程的解  $Pz=r$
     * 。在这个方程中 $P$
     * 应该近似于非线性系统的雅各布系数 $M-\gamma J$ 。
     * 如果没有提供jacobian_preconditioner_setup()函数，那么jacobian_preconditioner_solve()应该自己完成所有工作。
     * @note
     * 这个接口没有假设用户在这个函数中应该做什么。ARKode只假设在调用jacobian_preconditioner_setup()后，有可能调用jacobian_preconditioner_solve()。
     * @param[in]  t 当前的时间  @param[in]  y
     * 当前ARKode内部步骤的 $y$  向量  @param[in]  fy
     * 在y处的隐式右手边的当前值，  $f_I (t_n, y)$  。
     * @param[in]  jok
     * 一个输入标志，表示是否需要更新雅各布相关数据。jok参数提供了在预调节器求解函数中重复使用雅各布数据。当jok
     * = SUNFALSE时，雅各布相关数据应从头开始重新计算。
     * 当jok=SUNTRUE时，如果之前调用该函数时保存的雅各布数据可以被重新使用（用当前的gamma值）。jok
     * = SUNTRUE的调用只能发生在jok = SUNFALSE的调用之后。
     * @param[out]  jcur
     * 在输出时，如果重新计算了Jacobian数据，则应设置为SUNTRUE；如果没有重新计算Jacobian数据，但仍重新使用了保存的数据，则应设置为SUNFALSE。
     * @param[in]  gamma  $M-\gamma J$  中的值  $\gamma$
     * 。预处理程序应该近似于这个矩阵的逆值。
     * 这个函数应该返回。
     *
     *
     *
     *
     *
     * - 0: 成功
     *
     *
     *
     *
     * - >0: 可恢复的错误，ARKode将重新尝试解决方案并再次调用此函数。
     *
     *
     *
     *
     *
     *
     * - <0: 无法恢复的错误，计算将被中止，并抛出一个断言。
     *
     */
    std::function<int(double            t,
                      const VectorType &y,
                      const VectorType &fy,
                      int               jok,
                      int &             jcur,
                      double            gamma)>
      jacobian_preconditioner_setup;

    /**
     * 一个用户可以提供的函数对象，用于向SUNDIALS的内置求解器传递一个前置条件，或者在用户自己的线性求解中应用一个自定义的前置条件，在solve_mass()中指定。
     * 这个函数应该计算出预处理方程 $Pz=r$
     * 的解，并将其存储在 @p z. 中，在这个方程中 $P$
     * 应该近似于质量矩阵 $M$  。          @param[in]  t
     * 当前时间  @param[in]  r 先决条件方程的右侧  @param[out]  z
     * 应用先决条件的解决方案，即解决  $Pz=r$   @param[in]
     * gamma 先决条件方程中的数值  $\gamma$   @param[in]  ] tol
     * 系统应该被解决的公差  @param[in]  lr
     * 一个输入标志，表示预处理程序的求解是使用左预处理程序（lr
     * = 1）还是右预处理程序（lr =
     * 2）。只有在与SUNDIALS打包的求解器一起使用时才相关。如果与自定义的solven_mass()函数一起使用，这个函数将被设置为零。
     * 这个函数应该返回。
     *
     *
     *
     *
     *
     * - 0: 成功
     *
     *
     *
     *
     * - >0: 可恢复的错误，ARKode将重新尝试解决方案并再次调用此函数。
     *
     *
     *
     *
     *
     *
     * - <0: 无法恢复的错误，计算将被中止，并抛出一个断言。
     *
     */
    std::function<
      int(double t, const VectorType &r, VectorType &z, double tol, int lr)>
      mass_preconditioner_solve;

    /**
     * 一个用户可以提供的函数对象，用于设置在mass_preconditioner_solve()中指定的预处理程序。
     * 该函数应准备好预调节器方程的解  $Pz=r$
     * 。在这个方程中  $P$  应该近似于质量矩阵  $M$  。
     * 如果没有提供mass_preconditioner_setup()函数，那么mass_preconditioner_solve()应该自己完成所有工作。
     * @note
     * 本接口没有假设用户在这个函数中应该做什么。ARKode只假设在调用mass_preconditioner_setup()后，有可能调用mass_preconditioner_solve()。
     * @param[in]  t 当前时间 这个函数应该返回。
     *
     *
     *
     *
     *
     * - 0: 成功
     *
     *
     *
     *
     * - >0: 可恢复的错误，ARKode将重新尝试解决方案并再次调用此函数。
     *
     *
     *
     *
     *
     *
     * - <0: 无法恢复的错误，计算将被中止，并抛出一个断言。
     *
     */
    std::function<int(double t)> mass_preconditioner_setup;
#  endif

    /**
     * 一个用户可以提供的函数对象，其目的是对解决方案进行后处理。ARKode以固定的时间增量（每隔`output_period`秒）调用该函数，它被传递给解决方案的多项式插值，使用当前ARK顺序和（内部存储的）先前计算的解决方案步骤进行计算。
     * @note
     * ARKode内部计算的时间步长很可能比`output_period`步长大得多，因此通过简单地执行所有中间插值连续多次调用这个函数。这个函数被调用的次数和实际计算的时间步数之间没有关系。
     *
     */
    std::function<void(const double       t,
                       const VectorType & sol,
                       const unsigned int step_number)>
      output_step;

    /**
     * 一个用户可以提供的函数对象，其目的是评估求解器是否应该被重新启动（例如因为自由度的数量已经改变）。
     * 这个函数应该执行`sol'中所有必要的操作，以确保结果的向量是一致的，并且最终大小正确。
     * 例如，我们可以决定在时间t有必要进行局部细化。这个函数应该返回true，并改变`sol`的尺寸以反映新的尺寸。由于ARKode不知道新的维度，内部重置是必要的。
     * 默认实现只是返回
     * "false"，也就是说，在演化过程中不进行重启。
     *
     */
    std::function<bool(const double t, VectorType &sol)> solver_should_restart;

    /**
     * 一个用户可以提供的函数对象，目的是返回一个向量，其组成部分是ARKode用来计算向量法线的权重。这个函数的实现是可选的，它只在实现时使用。
     *
     */
    std::function<VectorType &()> get_local_tolerances;

    /**
     * 一个用户可以提供的函数对象，其目的是对提供的 @p
     * arkode_mem
     * 对象进行自定义设置。关于有效的选项请参考SUNDIALS文档。
     * 例如，下面的代码附加了两个文件，用于内部ARKODE实现的诊断和错误输出。
     * @code
     *    ode.custom_setup = [&](voidarkode_mem) {
     *      ARKStepSetErrFile(arkode_mem, errfile);
     *      ARKStepSetDiagnostics(arkode_mem, diagnostics_file);
     *    };
     * @endcode
     *
     * @note
     * 这个函数将在所有其他设置的末尾被调用，就在实际的时间演化开始或用solve_ode()继续之前。当求解器重新启动时也会调用这个函数，见solver_should_restart()。请查阅SUNDIALS手册，看看这时哪些选项仍然可用。
     * @param  arkode_mem
     * 指向ARKODE内存块的指针，可用于自定义调用`ARKStepSet...`方法。
     *
     */
    std::function<void(void *arkode_mem)> custom_setup;

  private:
    /**
     * 当一个具有给定名称的函数没有实现时，抛出一个异常。
     *
     */
    DeclException1(ExcFunctionNotProvided,
                   std::string,
                   << "Please provide an implementation for the function \""
                   << arg1 << "\"");

    /**
     * 反复调用ARKode的内部例程。
     *
     */
    int
    do_evolve_time(VectorType &          solution,
                   dealii::DiscreteTime &time,
                   const bool            do_reset);

#  if DEAL_II_SUNDIALS_VERSION_GTE(4, 0, 0)

    /**
     * 根据用户指定的函数，在ARKODE内存对象中设置（非）线性求解器和前置条件器。
     * @param solution 解的向量，它被用作创建新向量的模板。
     *
     */
    void
    setup_system_solver(const VectorType &solution);

    /**
     * 根据用户指定的函数，为ARKODE内存对象中的非同质性质量矩阵设置求解器和预处理器。
     * @param solution 被用作创建新向量模板的解向量。
     *
     */
    void
    setup_mass_solver(const VectorType &solution);

#  endif

    /**
     * 这个函数在构造时被执行，以设置上述 std::function
     * ，如果它们没有被实现，则触发断言。
     *
     */
    void
    set_functions_to_trigger_an_assert();

    /**
     * ARKode配置数据。
     *
     */
    AdditionalData data;

    /**
     * ARKode内存对象。
     *
     */
    void *arkode_mem;

    /**
     * MPI通信器。SUNDIALS求解器可以愉快地并行运行。注意，如果库的编译没有MPI支持，MPI_Comm被别名为int。
     *
     */
    MPI_Comm communicator;

    /**
     * 最后一次调用solve_ode()时的最后时间。
     *
     */
    double last_end_time;

#  if DEAL_II_SUNDIALS_VERSION_GTE(4, 0, 0)
    std::unique_ptr<internal::LinearSolverWrapper<VectorType>> linear_solver;
    std::unique_ptr<internal::LinearSolverWrapper<VectorType>> mass_solver;
#  endif

#  ifdef DEAL_II_WITH_PETSC
#    ifdef PETSC_USE_COMPLEX
    static_assert(!std::is_same<VectorType, PETScWrappers::MPI::Vector>::value,
                  "Sundials does not support complex scalar types, "
                  "but PETSc is configured to use a complex scalar type!");

    static_assert(
      !std::is_same<VectorType, PETScWrappers::MPI::BlockVector>::value,
      "Sundials does not support complex scalar types, "
      "but PETSc is configured to use a complex scalar type!");
#    endif // PETSC_USE_COMPLEX
#  endif   // DEAL_II_WITH_PETSC
  };


  /**
   * 处理 ARKode 异常。
   *
   */
  DeclException1(ExcARKodeError,
                 int,
                 << "One of the SUNDIALS ARKode internal functions "
                 << " returned a negative error code: " << arg1
                 << ". Please consult SUNDIALS manual.");

} // namespace SUNDIALS

DEAL_II_NAMESPACE_CLOSE
#endif


#endif


