//include/deal.II-translator/sundials/ida_0.txt
//-----------------------------------------------------------
//
//    Copyright (C) 2017 - 2020 by the deal.II authors
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

#ifndef dealii_sundials_ida_h
#define dealii_sundials_ida_h

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

#  ifdef DEAL_II_SUNDIALS_WITH_IDAS
#    include <idas/idas.h>
#  else
#    include <ida/ida.h>
#  endif

#  include <sundials/sundials_config.h>
#  if DEAL_II_SUNDIALS_VERSION_LT(3, 0, 0)
#    include <ida/ida_spbcgs.h>
#    include <ida/ida_spgmr.h>
#    include <ida/ida_sptfqmr.h>
#  endif
#  include <deal.II/sundials/sunlinsol_wrapper.h>

#  include <boost/signals2.hpp>

#  include <nvector/nvector_serial.h>
#  include <sundials/sundials_math.h>
#  include <sundials/sundials_types.h>

#  include <memory>


DEAL_II_NAMESPACE_OPEN

// Shorthand notation for IDA error codes.
#  define AssertIDA(code) Assert(code >= 0, ExcIDAError(code))

namespace SUNDIALS
{
  /**
   * SUNDIALS隐式微分代数（IDA）求解器的接口。
   * IDA类是SUNDIALS隐式微分代数求解器的一个封装器，它是一个通用的微分代数方程（DAE）系统的求解器。
   * 用户必须提供以下内容的实现  std::functions:  。
   *
   *
   *
   *
   *
   * - reinit_vector;
   * - 残留物。
   *
   * - setup_jacobian;
   *
   *
   *
   *
   *
   * - solve_jacobian_system/solve_with_jacobian; 函数`solve_jacobian_system`应该在SUNDIALS < 4.0.0中实现。对于后来的版本，你应该使用`solve_with_jacobian`来利用更好的非线性算法。    也可以选择提供以下函数。默认情况下，它们什么都不做，或者不需要。如果你在调用构造函数时需要一个未实现的函数，将抛出一个断言。
   *
   *
   *
   *
   *
   *
   * - solver_should_restart;
   *
   *
   * - 差异化的成分。
   *
   *
   *
   *
   * - get_local_tolerances; 要输出步骤，请将一个函数连接到信号上
   *
   *
   *
   *
   *
   * - output_step; 引自SUNDIALS文档。      考虑一个一般形式的微分代数方程组 \f[
   * \begin{cases} F(t,y,\dot y) = 0\, , \\ y(t_0) = y_0\, , \\ \dot y (t_0) =
   * \dot y_0\, . \end{cases} \f] ，其中  $y,\dot y$  是  $\mathbb{R}^n$
   * 中的向量，  $t$
   * 通常是时间（但也可以是一个参数量），以及
   * $F:\mathbb{R}\times\mathbb{R}^n\times
   * \mathbb{R}^n\rightarrow\mathbb{R}^n$  。
   * 这样的问题是用牛顿迭代和直线搜索全局策略来解决的。IDA中使用的积分方法是变阶、变系数的BDF（后向微分公式），为固定导程系数。方法的阶数从1到5，阶数为
   * $q$ 的BDF由多步公式\f[ \sum_{i=0}^q \alpha_{n,i}\,y_{n-i}=h_n\,\dot
   * y_n\, , \label{eq:bdf} \f]给出，其中 $y_n$ 和 $\dot y_n$ 分别是
   * $y(t_n)$ 和 $\dot y(t_n)$ 的计算近似值，步长为
   * $h_n=t_n-t_{n-1}$  。系数 $\alpha_{n,i}$ 由阶数 $q$
   * 和步长的历史唯一地决定。将BDF方法应用于DAE系统的结果是在每个时间步长上需要解决一个非线性代数系统。
   * \f[ G(y_n)\equiv F\left(t_n,y_n,\dfrac{1}{h_n}\sum_{i=0}^q
   * \alpha_{n,i}\,y_{n-i}\right)=0\, . \f]
   * 牛顿方法导致了一个形式为\f[
   * J[y_{n(m+1)}-y_{n(m)}]=-G(y_{n(m)})\, , \f]的线性系统，其中
   * $y_{n(m)}$ 是 $m$ 对 $y_n$ 的第1次近似， $J$
   * 是系统雅各布\f[ J=\dfrac{\partial G}{\partial y} = \dfrac{\partial
   * F}{\partial y} + \alpha \dfrac{\partial F}{\partial \dot y}\, , \f]和
   * $\alpha = \alpha_{n,0}/h_n$
   * 的近似。值得一提的是，每当步长或方法阶数改变时，标量
   * $\alpha$ 就会改变。
   * 为了提供一个简单的例子，考虑下面的谐波振荡器问题：
   * \f[ \begin{split} u'' & =
   *
   * -k^2 u \\ u (0) & = 0 \\ u'(0) & = k \end{split} \f]
   * 我们用一阶奥德来写它。 \f[ \begin{matrix} y_0' &
   *
   * -y_1      & = 0 \\ y_1' & + k^2 y_0 & = 0 \end{matrix} \f] 即 $F(y', y,
   * t) = y' + A y = 0 $  其中\f[ \begin{matrix} 0 &
   *
   * -1 \\ k^2 &0 \end{matrix} \f] 和 $y(0)=(0, k)$  ,  $y'(0) = (k, 0)$  。
   * 精确解是  $y_0(t) = \sin(k t)$  ,  $y_1(t) = y_0'(t) = k \cos(k t)$
   * ,  $y_1'(t) =
   *
   * -k^2 \sin(k t)$  。    要集合的雅各布系数如下。   $J =
   * \alpha I + A$  .     这可以通过以下代码片段实现。
   * @code
   * using VectorType = Vector<double>;
   *
   * VectorType y(2);
   * VectorType y_dot(2);
   *
   * double kappa = 1.0;
   *
   * FullMatrix<double> A(2,2);
   * A(0,1) =
   *
   * -1.0;
   * A(1,0) = kappa*kappa;
   *
   * FullMatrix<double> J(2,2);
   * FullMatrix<double> Jinv(2,2);
   *
   * IDA time_stepper;
   *
   * time_stepper.reinit_vector = [&] (VectorType&v)
   * {
   * v.reinit(2);
   * };
   *
   * time_stepper.residual = [&](const double t,
   *                           const VectorType &y,
   *                           const VectorType &y_dot,
   *                           VectorType &res)
   *
   * ->int
   * {
   * res = y_dot;
   * A.vmult_add(res, y);
   * return 0;
   * };
   *
   * time_stepper.setup_jacobian = [&](const double ,
   *                                 const VectorType &,
   *                                 const VectorType &,
   *                                 const double alpha)
   *
   * ->int
   * {
   * J = A;
   *
   * J(0,0) = alpha;
   * J(1,1) = alpha;
   *
   * Jinv.invert(J);
   * return 0;
   * };
   *
   * time_stepper.solve_jacobian_system = [&](const VectorType &src,
   *                                        VectorType &dst)
   *
   * ->int
   * {
   * Jinv.vmult(dst,src);
   * return 0;
   * };
   *
   * y[1] = kappa;
   * y_dot[0] = kappa;
   * time_stepper.solve_dae(y,y_dot);
   * @endcode
   *
   *
   */
  template <typename VectorType = Vector<double>>
  class IDA
  {
  public:
    /**
     * 可以传递给IDA类的额外参数。
     *
     */
    class AdditionalData
    {
    public:
      /**
       * IDA是一个微分代数求解器。因此，它对一阶导数也需要初始条件。如果你没有提供一致的初始条件，（即F(y_dot(0),
       * y(0),
       * 0)=0的条件），你可以要求SUNDIALS为你计算初始条件，在`initial_time`（`ic_type`）和复位后（`reset_type`）指定初始条件的InitialConditionCorrection。
       *
       */
      enum InitialConditionCorrection
      {
        /**
         * 不要试图使初始条件一致。
         *
         */
        none = 0,

        /**
         * 计算y的代数分量和y_dot的微分分量，给定y的微分分量。该选项要求用户在函数get_differential_components中指定微分和代数分量。
         *
         */
        use_y_diff = 1,

        /**
         * 计算y的所有分量，给定y_dot。
         *
         */
        use_y_dot = 2
      };

      /**
       * IDA的初始化参数。            全局参数。
       * @param  initial_time 初始时间  @param  final_time 最终时间
       * @param  initial_step_size 初始步长  @param  output_period
       * 每次输出的时间间隔 运行参数。              @param
       * minimum_step_size 最小步长  @param  maximum_order 最大BDF阶数
       * @param  maximum_non_linear_iterations 最大非线性迭代次数
       * @param  ls_norm_factor
       * 从积分器容限到线性求解器容限迭代的转换系数
       * 错误参数。              @param  absolute_tolerance
       * 绝对误差公差  @param  relative_tolerance 相对误差公差
       * @param  ignore_algebraic_terms_for_errors
       * 忽略误差计算的代数项。              @param  ic_type
       * 初始条件修正类型  @param  reset_type
       * 重启后的初始条件修正类型  @param
       * maximum_non_linear_iterations_ic 初始条件Newton最大迭代次数
       *
       */
      AdditionalData( // Initial parameters
        const double initial_time      = 0.0,
        const double final_time        = 1.0,
        const double initial_step_size = 1e-2,
        const double output_period     = 1e-1,
        // Running parameters
        const double       minimum_step_size             = 1e-6,
        const unsigned int maximum_order                 = 5,
        const unsigned int maximum_non_linear_iterations = 10,
        const double       ls_norm_factor                = 0,
        // Error parameters
        const double absolute_tolerance                = 1e-6,
        const double relative_tolerance                = 1e-5,
        const bool   ignore_algebraic_terms_for_errors = true,
        // Initial conditions parameters
        const InitialConditionCorrection &ic_type    = use_y_diff,
        const InitialConditionCorrection &reset_type = use_y_diff,
        const unsigned int                maximum_non_linear_iterations_ic = 5)
        : initial_time(initial_time)
        , final_time(final_time)
        , initial_step_size(initial_step_size)
        , minimum_step_size(minimum_step_size)
        , absolute_tolerance(absolute_tolerance)
        , relative_tolerance(relative_tolerance)
        , maximum_order(maximum_order)
        , output_period(output_period)
        , ignore_algebraic_terms_for_errors(ignore_algebraic_terms_for_errors)
        , ic_type(ic_type)
        , reset_type(reset_type)
        , maximum_non_linear_iterations_ic(maximum_non_linear_iterations_ic)
        , maximum_non_linear_iterations(maximum_non_linear_iterations)
        , ls_norm_factor(ls_norm_factor)
      {}

      /**
       * 将所有AdditionalData()参数添加到给定的ParameterHandler对象中。当参数从文件中被解析出来时，内部参数会自动更新。
       * 声明了以下参数。
       * @code
       * set Final time                        = 1.000000
       * set Initial time                      = 0.000000
       * set Time interval between each output = 0.2
       * subsection Error control
       * set Absolute error tolerance                      = 0.000001
       * set Ignore algebraic terms for error computations = true
       * set Relative error tolerance                      = 0.00001
       * set Use local tolerances                          = false
       * end
       * subsection Initial condition correction parameters
       * set Correction type at initial time        = none
       * set Correction type after restart          = none
       * set Maximum number of nonlinear iterations = 5
       * end
       * subsection Running parameters
       * set Initial step size                      = 0.1
       * set Maximum number of nonlinear iterations = 10
       * set Maximum order of BDF                   = 5
       * set Minimum step size                      = 0.000001
       * end
       * @endcode
       * 这些参数与你在构建时可以传递的选项是一一对应的。
       * 你在构建时传递的选项在ParameterHandler对象`prm`中被设置为默认值。之后你可以通过使用`prm`解析参数文件来修改它们。每当`prm`的内容被更新时，参数的值就会被更新。
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
        prm.add_parameter("Maximum order of BDF", maximum_order);
        prm.add_parameter("Maximum number of nonlinear iterations",
                          maximum_non_linear_iterations);
        prm.leave_subsection();

        prm.enter_subsection("Error control");
        prm.add_parameter("Absolute error tolerance", absolute_tolerance);
        prm.add_parameter("Relative error tolerance", relative_tolerance);
        prm.add_parameter(
          "Ignore algebraic terms for error computations",
          ignore_algebraic_terms_for_errors,
          "Indicate whether or not to suppress algebraic variables "
          "in the local error test.");
        prm.leave_subsection();

        prm.enter_subsection("Initial condition correction parameters");
        static std::string ic_type_str = "use_y_diff";
        prm.add_parameter(
          "Correction type at initial time",
          ic_type_str,
          "This is one of the following three options for the "
          "initial condition calculation. \n"
          " none: do not try to make initial conditions consistent. \n"
          " use_y_diff: compute the algebraic components of y and differential\n"
          "    components of y_dot, given the differential components of y. \n"
          "    This option requires that the user specifies differential and \n"
          "    algebraic components in the function get_differential_components.\n"
          " use_y_dot: compute all components of y, given y_dot.",
          Patterns::Selection("none|use_y_diff|use_y_dot"));
        prm.add_action("Correction type at initial time",
                       [&](const std::string &value) {
                         if (value == "use_y_diff")
                           ic_type = use_y_diff;
                         else if (value == "use_y_dot")
                           ic_type = use_y_dot;
                         else if (value == "none")
                           ic_type = none;
                         else
                           AssertThrow(false, ExcInternalError());
                       });

        static std::string reset_type_str = "use_y_diff";
        prm.add_parameter(
          "Correction type after restart",
          reset_type_str,
          "This is one of the following three options for the "
          "initial condition calculation. \n"
          " none: do not try to make initial conditions consistent. \n"
          " use_y_diff: compute the algebraic components of y and differential\n"
          "    components of y_dot, given the differential components of y. \n"
          "    This option requires that the user specifies differential and \n"
          "    algebraic components in the function get_differential_components.\n"
          " use_y_dot: compute all components of y, given y_dot.",
          Patterns::Selection("none|use_y_diff|use_y_dot"));
        prm.add_action("Correction type after restart",
                       [&](const std::string &value) {
                         if (value == "use_y_diff")
                           reset_type = use_y_diff;
                         else if (value == "use_y_dot")
                           reset_type = use_y_dot;
                         else if (value == "none")
                           reset_type = none;
                         else
                           AssertThrow(false, ExcInternalError());
                       });
        prm.add_parameter("Maximum number of nonlinear iterations",
                          maximum_non_linear_iterations_ic);
        prm.add_parameter(
          "Factor to use when converting from the integrator tolerance to the linear solver tolerance",
          ls_norm_factor);
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
       * BDF的最大顺序。
       *
       */
      unsigned int maximum_order;

      /**
       * 每个输出之间的时间周期。
       *
       */
      double output_period;

      /**
       * 忽略错误的代数项。
       *
       */
      bool ignore_algebraic_terms_for_errors;

      /**
       * 对初始条件的修正类型。
       * 如果你没有提供一致的初始条件，（即 $F(y_dot(0),
       * y(0), 0) = 0$
       * 的条件），你可以在构造时使用`ic_type`参数，要求SUNDIALS为你计算初始条件。
       * 请注意，原则上你可以使用这个功能来解决稳态问题，将y_dot设置为零，并要求计算出满足
       * $F(0, y(0), 0) = 0$ 的 $y(0)$
       * ，然而IDA内部使用的非线性求解器可能对有几百万个未知数的复杂问题不够强大。
       *
       */
      InitialConditionCorrection ic_type;

      /**
       * 解算器重启后使用的初始条件修正类型。
       * 如果你在重启后没有一致的初始条件，（即F(y_dot(t_restart),
       * y(t_restart), t_restart) =
       * 0的条件），你可以要求SUNDIALS在构造时使用`reset_type`参数为你计算新的初始条件。
       *
       */
      InitialConditionCorrection reset_type;

      /**
       * IC计算中牛顿法的最大迭代次数。
       *
       */
      unsigned maximum_non_linear_iterations_ic;

      /**
       * 在时间推进过程中，牛顿方法的最大迭代次数。
       *
       */
      unsigned int maximum_non_linear_iterations;

      /**
       * 从积分器公差转换到线性求解器公差时使用的系数。
       *
       */
      double ls_norm_factor;
    };

    /**
     * 构造器。通过传递一个设定所有求解器参数的AdditionalData()对象，可以对SUNDIALS
     * IDA求解器进行微调。
     * IDA是一个微分代数求解器。因此，它对一阶导数也需要初始条件。如果你没有提供一致的初始条件，（即F(y_dot(0),
     * y(0),
     * 0)=0的条件），你可以在构造时使用`ic_type`参数要求SUNDIALS为你计算初始条件。
     * 你有三个选择
     *
     *
     *
     *
     *
     * - 无：不要试图使初始条件一致。
     *
     *
     *
     * - use_y_diff: 计算y的代数成分和y_dot的微分成分，给定y的微分成分。这个选项要求用户在函数get_differential_components中指定微分和代数成分。
     *
     *
     *
     *
     *
     *
     * - use_y_dot: 计算y的所有分量，给定y_dot。        默认情况下，该类假设所有分量都是微分，并且你想解决一个标准的颂歌。在这种情况下，初始分量类型被设置为`use_y_diff'，因此在时间t=`初始时间'的`y_dot'是通过解决变量`y_dot'的非线性问题 $F(y_dot,
     * y(t0), t0) = 0$ 计算的。
     * 请注意，牛顿求解器被用于这一计算。牛顿求解器的参数可以通过作用于`ic_alpha`和`ic_max_iter`进行调整。
     * 如果你在某个时候重置求解器，你可能想在重置后为初始条件选择一个不同的计算。比如说，你完善了一个网格，在将解转移到新的网格后，初始条件不再一致了。那么你可以选择如何使其一致，使用与`reset_type`中初始条件相同的三个选项。
     * 在串行情况下，MPI通信器被简单地忽略了。
     * @param  data IDA配置数据  @param  mpi_comm MPI通信器
     *
     */
    IDA(const AdditionalData &data     = AdditionalData(),
        const MPI_Comm &      mpi_comm = MPI_COMM_WORLD);

    /**
     * 解构器。
     *
     */
    ~IDA();

    /**
     * 对微分代数方程进行积分。该函数返回最终的计算步骤数。
     *
     */
    unsigned int
    solve_dae(VectorType &solution, VectorType &solution_dot);

    /**
     * 清理内部内存，以干净的对象开始。这个函数在模拟开始时和用户对solver_should_restart()的调用返回true时被调用。
     * 默认情况下，solver_should_restart()返回false。如果用户需要实现例如空间的局部适应性，可以给solver_should_restart()指定一个不同的函数，执行所有的网格变化，将解和解点转移到新的网格，并返回true。
     * 在reset()过程中，y和yp都会被检查是否一致，根据指定的ic_type（如果t==initial_time）或reset_type（如果t>initial_time），yp、y或两者都被修改以获得一致的初始数据集。
     * @param[in]  t 新的起始时间  @param[in]  h
     * 新的（暂定）起始时间步骤  @param[in,out]  y
     * 新的（暂定）初始解  @param[in,out]  yp
     * 新的（暂定）初始解_dot
     *
     */
    void
    reset(const double t, const double h, VectorType &y, VectorType &yp);

    /**
     * 重新设置向量，使其具有正确的大小和MPI通信器等。
     *
     */
    std::function<void(VectorType &)> reinit_vector;

    /**
     * 计算残差。返回  $F(t, y, \dot y)$  。
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
     * - >0: 可恢复的错误（如果发生这种情况，将调用IDAReinit，然后再次尝试最后一个函数
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
                      const VectorType &y,
                      const VectorType &y_dot,
                      VectorType &      res)>
      residual;

    /**
     * 计算雅各布系数。IDA在任何需要更新雅各布的时候都会调用这个函数。用户应该计算Jacobian（或者更新所有允许应用Jacobian的变量）。IDA在调用solve_jacobian_system()（适用于SUNDIALS
     * < 4.0.0）或solve_with_jacobian()（适用于SUNDIALS >=
     * 4.0.0）之前，会调用该函数一次。        雅各邦 $J$
     * 应该是\f[ J=\dfrac{\partial G}{\partial y} = \dfrac{\partial
     * F}{\partial y} + \alpha \dfrac{\partial F}{\partial \dot y}.
     * \f]的一个（可能是不精确的）计算，如果用户使用基于矩阵的雅各邦计算，那么在这里应该调用一个装配例程来装配雅各邦系统的矩阵和预处理程序。
     * 随后调用（可能不止一次）solve_jacobian_system()或solve_with_jacobian()可以假设这个函数至少被调用过一次。
     * 请注意，这个接口没有假设用户在这个函数中应该做什么。IDA只假设在调用setup_jacobian()后，有可能调用solve_jacobian_system()或solve_with_jacobian()来获得系统的解
     * $x$  。        这个函数应该返回。
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
     * - >0: 可恢复的错误（如果发生这种情况，将调用IDAReinit，然后再次尝试最后一个函数
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
                      const VectorType &y,
                      const VectorType &y_dot,
                      const double      alpha)>
      setup_jacobian;

    /**
     * 解决Jacobian线性系统。这个函数将在setup_jacobian()被调用至少一次之后被IDA调用（可能是多次）。IDA试图尽最大努力调用setup_jacobian()最少的次数。如果不更新雅各布式就能实现收敛，那么IDA就不会再次调用setup_jacobian()。相反，如果IDA内部收敛测试失败，那么IDA会用更新的向量和系数再次调用setup_jacobian()，这样连续调用solve_jacobian_systems()会导致牛顿过程中更好的收敛。
     * 雅可比 $J$ 应该是系统雅可比\f[ J=\dfrac{\partial G}{\partial
     * y} = \dfrac{\partial F}{\partial y} + \alpha \dfrac{\partial
     * F}{\partial \dot y}.
     * \f]的（近似值）。对该函数的调用应该在`dst`中存储
     * $J^{-1}$ 应用于`src`的结果，即`J*dst =
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
     * - >0: 可恢复的错误（如果发生这种情况，将调用IDAReinit，然后再次尝试最后一个函数
     *
     *
     *
     *
     *
     *
     * - <0: 无法恢复的错误，计算将被中止，并抛出一个断言。          @warning  从SUNDIALS 4.1开始，SUNDIALS提供了指定分辨率的公差的可能性。从公差的一部分只提供`rhs`，`dst`需要返回。
     *
     */
    DEAL_II_DEPRECATED
    std::function<int(const VectorType &rhs, VectorType &dst)>
      solve_jacobian_system;

    /**
     * 解决雅各布线性系统，直到指定的公差。这个函数将在setup_jacobian()被调用至少一次之后被IDA调用（可能是多次）。IDA试图尽最大努力调用setup_jacobian()的最少次数。如果不更新雅各布式就能实现收敛，那么IDA就不会再次调用setup_jacobian()。相反，如果IDA内部收敛测试失败，那么IDA会用更新的向量和系数再次调用setup_jacobian()，这样连续调用solve_with_jacobian()会导致牛顿过程中更好的收敛。
     * 雅各邦 $J$ 应该是系统雅各邦\f[ J=\dfrac{\partial G}{\partial
     * y} = \dfrac{\partial F}{\partial y} + \alpha \dfrac{\partial
     * F}{\partial \dot y}.
     * \f]的（近似值），函数的参数是。          @param[in]  rhs
     * 要解决的系统右侧。      @param[out]  dst  $J^{-1} src$
     * 的解。      @param[in]  tolerance 解决线性方程组的公差。
     * 对该函数的调用应该在`dst`中存储 $J^{-1}$
     * 应用于`src`的结果，即线性系统`J*dst = src`的解。
     * 用户有责任在这个函数中或在`setup_jacobian()`函数中设置适当的求解器和预处理器。例如，后者是
     * step-77
     * 程序所做的。所有昂贵的操作都发生在`setup_jacobian()`中，因为该函数被调用的频率远低于当前函数）。)
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
     * - >0: 可恢复的错误（如果发生这种情况，将调用IDAReinit，然后再次尝试最后一个函数）。
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
      int(const VectorType &rhs, VectorType &dst, const double tolerance)>
      solve_with_jacobian;

    /**
     * 处理解决方案。这个函数由IDA在固定的时间步数上调用，每隔`output_period`秒，它被传递给解决方案的多项式插值及其时间导数，使用当前的BDF顺序和（内部存储的）先前计算的解决方案步数计算。
     * 请注意，IDA内部计算的时间步长很可能比`output_period`步长大得多，因此通过简单地执行所有中间插值，连续多次调用这个函数。这个函数被调用的次数和实际计算的时间步数之间没有关系。
     *
     */
    std::function<void(const double       t,
                       const VectorType & sol,
                       const VectorType & sol_dot,
                       const unsigned int step_number)>
      output_step;

    /**
     * 评估求解器是否应该被重新启动（例如因为自由度的数量发生了变化）。
     * 这个函数应该执行所有在`sol`和`sol_dot`中需要的操作，以确保得到的向量是一致的，并且最终大小正确。
     * 例如，我们可以决定在时间t有必要进行局部细化。这个函数应该返回true，并改变sol和sol_dot的尺寸以反映新的尺寸。由于IDA不知道新的维度，所以内部重置是必要的。
     * 默认实现只是返回
     * "false"，也就是说，在演化过程中不进行重启。
     *
     */
    std::function<bool(const double t, VectorType &sol, VectorType &sol_dot)>
      solver_should_restart;

    /**
     * 返回一个包含微分成分的索引集。
     * 这个函数的实现是可选的。默认是返回一个完整的索引集。如果你的方程也是代数的（即它包含代数约束，或拉格朗日乘数），你应该覆盖这个函数，以便只返回系统的微分成分。
     * 当并行运行时，每个进程都会独立地调用这个函数，同步将在初始化设置结束时发生，以沟通哪些组件是本地的。确保你只返回本地拥有的（或本地相关的）组件，以减少进程间的通信。
     *
     */
    std::function<IndexSet()> differential_components;

    /**
     * 返回一个向量，其成分是IDA用来计算向量法线的权重。这个函数的实现是可选的。如果用户没有提供实现，则假设权重为所有的1。
     *
     */
    std::function<VectorType &()> get_local_tolerances;

    /**
     * 处理IDA的异常。
     *
     */
    DeclException1(ExcIDAError,
                   int,
                   << "One of the SUNDIALS IDA internal functions "
                   << " returned a negative error code: " << arg1
                   << ". Please consult SUNDIALS manual.");


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
     * 这个函数在构造时被执行，以设置上面的 std::function
     * ，如果它们没有被实现，则触发一个断言。
     *
     */
    void
    set_functions_to_trigger_an_assert();

    /**
     * IDA配置数据。
     *
     */
    const AdditionalData data;

    /**
     * IDA内存对象。
     *
     */
    void *ida_mem;

    /**
     * MPI通信器。SUNDIALS解算器可以愉快地并行运行。注意，如果库的编译没有MPI支持，MPI_Comm被别名为int。
     *
     */
    MPI_Comm communicator;

    /**
     * 向量的内存池。
     *
     */
    GrowingVectorMemory<VectorType> mem;

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
} // namespace SUNDIALS

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_SUNDIALS

#endif


