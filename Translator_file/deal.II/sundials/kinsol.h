//include/deal.II-translator/sundials/kinsol_0.txt
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

#ifndef dealii_sundials_kinsol_h
#define dealii_sundials_kinsol_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_SUNDIALS


#  include <deal.II/base/conditional_ostream.h>
#  include <deal.II/base/exceptions.h>
#  include <deal.II/base/logstream.h>
#  include <deal.II/base/mpi.h>
#  include <deal.II/base/parameter_handler.h>

#  include <deal.II/lac/vector.h>
#  include <deal.II/lac/vector_memory.h>

#  include <boost/signals2.hpp>

#  include <kinsol/kinsol.h>
#  if DEAL_II_SUNDIALS_VERSION_LT(4, 1, 0)
#    include <kinsol/kinsol_impl.h>
#  endif
#  include <nvector/nvector_serial.h>
#  include <sundials/sundials_math.h>
#  include <sundials/sundials_types.h>

#  include <memory>


DEAL_II_NAMESPACE_OPEN

// Shorthand notation for KINSOL error codes.
#  define AssertKINSOL(code) Assert(code >= 0, ExcKINSOLError(code))

namespace SUNDIALS
{
  /**
   * SUNDIALS的非线性求解器（KINSOL）的接口。
   * KINSOL是一个用于残差形式 $F(u) = 0$ 或定点形式 $G(u) = u$
   * 的非线性代数系统的求解器，其中 $u$
   * 是一个向量，我们将假设它在 ${\mathbb R}^n$ 或 ${\mathbb
   * C}^n$
   * 中，但它也可能有一个块结构，可以分布在并行计算中；函数
   * $F$ 和 $G$  满足 $F,G:{\mathbb R}^N \to{\mathbb R}^N$  或
   * $F,G:{\mathbb C}^N \to{\mathbb C}^N$
   * 。它包括牛顿-克雷洛夫求解器以及皮卡德和定点求解器，这两种求解器都可以用安德森加速器进行加速。KINSOL是基于Brown和Saad以前的Fortran软件包NKSOL。一个使用KINSOL的例子可以在
   * step-77 的教程程序中找到。
   * KINSOL的牛顿求解器采用了不精确的牛顿方法。由于这个求解器主要用于大型系统，用户需要提供他们自己的求解器函数。
   * 在最高级别，KINSOL实现了以下迭代方案。
   *
   *
   *
   *
   *
   *
   * - 设置 $u_0$ =一个初始猜测
   *
   *
   *
   *
   *
   *
   * - 对于 $n = 0, 1, 2, \ldots$ 直到收敛做。
   *
   *
   *
   *
   * - 解决  $J(u_n)\delta_n =
   *
   * -F(u_n)$
   *
   *
   *
   *
   *
   * - 设置 $u_{n+1} = u_n + \lambda \delta_n, 0 < \lambda \leq 1$
   *
   *
   *
   *
   * - 检验收敛性 这里， $u_n$ 是 $n$ 到 $u$ 的第1次迭代，而 $J(u) = \nabla_u F(u)$ 是系统的雅各布系数。在迭代过程的每个阶段，步骤 $\delta_n$ 的标量倍数被添加到 $u_n$ ，以产生一个新的迭代， $u_{n+1}$ 。在继续迭代之前，要进行收敛测试。    除非用户另有规定，否则KINSOL会尽可能少地更新雅各布信息，以平衡矩阵操作的高成本和其他成本。具体来说，这些更新在以下情况下发生。
   *
   *
   *
   *
   * - 问题被初始化了。
   *
   *
   *
   *
   *
   * -  $\|\lambda \delta_{n-1} \|_{D_u,\infty} \geq 1.5$ （仅有不确切的牛顿，见下文对 $\| \cdot \|_{D_u,\infty}$ 的定义）。
   *
   *
   *
   *
   *
   * - 自上次更新以来，已经过了指定数量的非线性迭代。
   *
   *
   *
   *
   * - 线性求解器在过期的雅各布信息下失败，可以恢复。
   *
   *
   *
   *
   * - 全局策略因过期的雅各布信息而失败，或
   *
   *
   *
   *
   *
   * -  $\|\lambda \delta_{n} \|_{D_u,\infty} \leq $ 容忍*与过时的雅各布信息。    KINSOL允许通过可选的求解器输入来改变上述策略。用户可以禁用初始雅各布信息评估，或者改变强制更新雅各布信息的非线性迭代次数的默认值。    为了解决条件不良的非线性系统的情况，KINSOL允许为解向量和残差向量规定比例系数。为了使用缩放，用户可以提供函数get_solution_scaling()，该函数返回值 $D_u$ ，它是缩放矩阵的对角线元素，当 $u_n$ 接近解时， $D_u u_n$ 的所有分量都大致相同。 ] 接近于一个解决方案，而get_function_scaling()，提供的值 $D_F$ ，是缩放矩阵的对角线元素，当 $u_n$ 不*太接近于一个解决方案时， $D_F F$ 的所有分量大致相同。    当为解向量提供缩放值时，如果用户没有通过solve_jacobian_system()函数提供雅各布信息的默认差分商数近似值，这些值将自动并入计算的扰动。    有两种将计算的步骤 $\delta_n$ 应用于先前计算的解向量的方法被实现。第一种也是最简单的是标准的牛顿策略，它以常数 $\lambda$ 始终设置为1来应用更新。另一种方法是全局策略，它试图以最有效的方式使用 $\delta_n$ 所暗示的方向来促进非线性问题的收敛。这种技术在第二个策略中实现，称为Linesearch。这个方案同时采用了J.E.Dennis和R.B.B.Dennis所给出的Goldstein-Armijo linesearch算法的 $\alpha$ 和 $\beta$ 的条件。E. Dennis和R. B. Schnabel。"Numerical Methods for Unconstrained Optimization and Nonlinear Equations."。SIAM, Philadelphia, 1996.*，其中 $\lambda$ 的选择是为了保证 $F$ 相对于步长的充分下降，以及相对于 $F$ 的初始下降率的最小步长。该算法的一个特性是，牛顿全步趋向于在接近解的地方进行。    在KINSOL中实现的基本定点迭代方案是这样的：。
   *
   *
   *
   *
   *
   * - 设置 $u_0 =$ 一个初始猜测
   *
   *
   *
   *
   *
   *
   * - 对于 $n = 0, 1, 2, \dots$ 直到收敛做。
   *
   *
   *
   *
   *
   * - 设置 $u_{n+1} = G(u_n)$
   *
   *
   *
   *
   * - 收敛性测试 在迭代过程的每个阶段，函数  $G$  被应用于当前迭代，产生一个新的迭代，  $u_{n+1}$  。在继续迭代之前要进行收敛性测试。    对于KINSOL中实现的Picard迭代，我们考虑非线性函数 $F$ 的特殊形式，如 $F(u) = Lu
   *
   * - N(u)$ ，其中 $L$ 是一个常数非辛格矩阵， $N$
   * 是（一般）非线性的。    那么固定点函数 $G$ 定义为
   * $G(u) = u
   *
   * - L^{-1}F(u)$  。
   * 在每个迭代中，Picard步骤被计算出来，然后加入到 $u_n$
   * 中，产生新的迭代。接下来，非线性残差函数在新迭代处被评估，并检查收敛性。使用安德森的方法可以大大加快Picard和固定点方法的速度。
   * 用户必须提供以下内容的实现 std::functions:  。
   *
   *
   *
   *
   *
   * - reinit_vector；并且只有一个
   *
   *
   * - 剩余的；或
   * - iteration_function; 指定residual()允许用户使用牛顿和皮卡德策略（即 $F(u)=0$ 将被解决），而指定iteration_function()，将使用固定点迭代（即 $G(u)=u$ 将被解决）。    如果希望使用牛顿或皮卡德方法，那么用户还应该提供
   *
   *
   *
   *
   *
   * - solve_jacobian_system 或 solve_with_jacobian；和可选的
   *
   *
   *
   *
   * - setup_jacobian; 固定点迭代不需要解决任何线性系统。    另外，以下函数也可以改写，以便在收敛检查时为解和残差评估提供额外的缩放因子。
   *
   *
   *
   *
   * - get_solution_scaling;
   *
   *
   * - get_function_scaling;
   *
   */
  template <typename VectorType = Vector<double>>
  class KINSOL
  {
  public:
    /**
     * 可以传递给KINSOL类的额外参数。
     *
     */
    class AdditionalData
    {
    public:
      /**
       * KINSOL的求解策略。KINSOL包括一个Newton-Krylov求解器（包括局部和全局）以及Picard和固定点求解器。
       *
       */
      enum SolutionStrategy
      {
        /**
         * 标准牛顿迭代。
         *
         */
        newton = KIN_NONE,
        /**
         * 牛顿迭代与线搜索。
         *
         */
        linesearch = KIN_LINESEARCH,
        /**
         * 固定点迭代。
         *
         */
        fixed_point = KIN_FP,
        /**
         * 皮卡德迭代。
         *
         */
        picard = KIN_PICARD,
      };

      /**
       * KINSOL的初始化参数。            全局参数。
       * @param  strategy 解决策略  @param  maximum_non_linear_iterations
       * 非线性迭代的最大次数  @param  function_tolerance
       * %函数规范的停止公差  @param  step_tolerance
       * 缩放的步骤停止公差 牛顿参数。              @param
       * no_init_setup 无初始矩阵设置  @param  maximum_setup_calls
       * 无矩阵设置的最大迭代次数  @param  maximum_newton_step
       * 牛顿步骤的最大允许比例长度  @param  dq_relative_error
       * 不同商数计算的相对误差 行搜索参数。
       * @param  maximum_beta_failures 最大β条件失败次数
       * 固定点和Picard参数。              @param
       * anderson_subspace_size 安德森加速子空间大小
       *
       */
      AdditionalData(const SolutionStrategy &strategy = linesearch,
                     const unsigned int maximum_non_linear_iterations = 200,
                     const double       function_tolerance            = 0.0,
                     const double       step_tolerance                = 0.0,
                     const bool         no_init_setup                 = false,
                     const unsigned int maximum_setup_calls           = 0,
                     const double       maximum_newton_step           = 0.0,
                     const double       dq_relative_error             = 0.0,
                     const unsigned int maximum_beta_failures         = 0,
                     const unsigned int anderson_subspace_size        = 0);

      /**
       * 将所有AdditionalData()参数添加到给定的ParameterHandler对象中。当参数从文件中被解析出来时，内部参数会自动更新。
       * 以下参数被声明。
       * @code
       * set Function norm stopping tolerance       = 0
       * set Maximum number of nonlinear iterations = 200
       * set Scaled step stopping tolerance         = 0
       * set Solution strategy                      = linesearch
       * subsection Fixed point and Picard parameters
       * set Anderson acceleration subspace size = 5
       * end
       * subsection Linesearch parameters
       * set Maximum number of beta-condition failures = 0
       * end
       * subsection Newton parameters
       * set Maximum allowable scaled length of the Newton step = 0
       * set Maximum iterations without matrix setup            = 0
       * set No initial matrix setup                            = false
       * set Relative error for different quotient computation  = 0
       * end
       * @endcode
       * 这些参数与你在构建时可以传递的选项是一一对应的。
       * 你在构建时传递的选项在ParameterHandler对象`prm`中被设置为默认值。之后你可以通过使用`prm`解析参数文件来修改它们。每当`prm`的内容被更新时，参数的值就会被更新。
       * 请确保这个类的寿命比`prm`长。如果你破坏了这个类，然后用`prm`解析一个参数文件，将会发生未定义的行为。
       *
       */
      void
      add_parameters(ParameterHandler &prm);

      /**
       * 要使用的解决策略。如果你选择 SolutionStrategy::newton
       * 或 SolutionStrategy::linesearch,
       * ，你必须同时提供函数residual（）。如果你选择
       * SolutionStrategy::picard 或 SolutionStrategy::fixed_point,
       * ，你还必须提供函数iteration_function()。
       *
       */
      SolutionStrategy strategy;

      /**
       * 允许的最大非线性迭代次数。
       *
       */
      unsigned int maximum_non_linear_iterations;

      /**
       * 一个标量，用来作为系统函数的最大尺度的停止容忍度
       * $F(u)$  或  $G(u)$  。
       * 如果设置为零，将使用KINSOL提供的默认值。
       *
       */
      double function_tolerance;

      /**
       * 一个标量，用作最小缩放步长的停止公差。
       * 如果设置为零，将使用KINSOL提供的默认值。
       *
       */
      double step_tolerance;

      /**
       * 是否应该对预处理程序或雅各布设置函数进行初始调用。
       * 当解决一连串的问题时，调用这个函数是很有用的，在这个过程中，一个问题的最终预处理或雅各布值将被用于下一个问题的初始化。
       *
       */
      bool no_init_setup;

      /**
       * 在调用setup_jacobian()函数之间可以执行的最大非线性迭代次数。
       * 如果设置为零，将使用KINSOL提供的默认值，在实践中，这通常意味着KINSOL将在以后的迭代中重新使用在一次迭代中计算的雅各布矩阵。
       *
       */
      unsigned int maximum_setup_calls;

      /**
       * 牛顿步长的最大允许比例长度。
       * 如果设置为零，将使用KINSOL提供的默认值。
       *
       */
      double maximum_newton_step;

      /**
       * 计算 $F(u)$
       * 的相对误差，当用户没有提供solve_jacobian_system_matrix()函数时，该误差用于对Jacobian矩阵的差分商数近似。
       * 如果设置为零，将使用KINSOL提供的默认值。
       *
       */
      double dq_relative_error;

      /**
       * 线路搜索算法中β条件失败的最大数量。只在
       * strategy==SolutionStrategy::linesearch. 时使用。
       *
       */
      unsigned int maximum_beta_failures;

      /**
       * 安德森加速与Picard或定点迭代一起使用的子空间的大小。
       * 如果你设置为0，则不使用加速。
       *
       */
      unsigned int anderson_subspace_size;
    };

    /**
     * 构造器。通过传递一个设置所有求解器参数的AdditionalData()对象，可以对SUNDIALS
     * KINSOL求解器进行微调。          @param  data KINSOL配置数据
     * @param  mpi_comm MPI通信器
     *
     */
    KINSOL(const AdditionalData &data     = AdditionalData(),
           const MPI_Comm &      mpi_comm = MPI_COMM_WORLD);

    /**
     * 解构器。
     *
     */
    ~KINSOL();

    /**
     * 解决非线性系统。返回为收敛而采取的非线性步骤的数量。KINSOL使用`initial_guess_and_solution`的内容作为初始猜测，并将最终解存储在同一个向量中。
     *
     */
    unsigned int
    solve(VectorType &initial_guess_and_solution);

    /**
     * 一个用户需要提供的函数对象，其目的是将给定的向量重新itize为正确的大小、块结构（如果使用块向量）和MPI通信器（如果向量使用MPI分布在多个处理器上），以及任何其他必要的属性。
     *
     */
    std::function<void(VectorType &)> reinit_vector;

    /**
     * 一个用户应该提供的函数对象，其目的是计算残差`dst
     * = F(src)`。该函数仅在选择了 SolutionStrategy::newton 或
     * SolutionStrategy::linesearch 策略时使用。
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
     * - >0: 可恢复的错误（KINSOL将尝试改变其内部参数并尝试一个新的解决步骤
     *
     *
     *
     *
     *
     * - <0: 无法恢复的错误，计算将被中止，并抛出一个断言。
     *
     */
    std::function<int(const VectorType &src, VectorType &dst)> residual;

    /**
     * 一个用户应该提供的函数对象，其目的是为了计算固定点和Picard迭代的迭代函数
     * $G(u)$ 。该函数仅在选择 SolutionStrategy::fixed_point 或
     * SolutionStrategy::picard 策略时使用。
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
     * - >0: 可恢复的错误（KINSOL将尝试改变其内部参数并尝试一个新的解决步骤
     *
     *
     *
     *
     *
     * - <0: 无法恢复的错误；计算将被中止，并抛出一个断言。
     *
     */
    std::function<int(const VectorType &src, VectorType &dst)>
      iteration_function;

    /**
     * 一个用户可以提供的函数对象，其目的是为后续调用solve_jacobian_system()准备线性求解器。
     * setup_jacobian()的工作是为后续调用solve_jacobian_system()准备线性求解器，在解决线性系统
     * $Ax = b$
     * 。这个系统的确切性质取决于已经选择的SolutionStrategy。
     * 在strategy =  SolutionStrategy::newton  或
     * SolutionStrategy::linesearch,  的情况下  $A$  是雅各布 $J =
     * \partial
     * F/\partial u$  。如果策略= SolutionStrategy::picard,   $A$
     * 是近似的雅各布矩阵  $L$  。如果策略=
     * SolutionStrategy::fixed_point,
     * ，则不会出现线性系统，这个函数也不会被调用。
     * setup_jacobian()函数可以调用一个用户提供的函数，或者线性求解器模块中的一个函数，以计算线性求解器所需要的雅各布相关数据。它还可以根据solve_jacobian_system()的需要对这些数据进行预处理，这可能涉及到调用一个通用函数（比如用于LU因子化），或者更普遍的是，从组装的雅各布式建立预处理程序。在任何情况下，这样生成的数据都可以在解决线性系统的时候被使用。
     * 这个函数的意义在于，setup_jacobian()函数不是在每次牛顿迭代时都被调用，而是在求解器确定适合执行设置任务时才被调用。这样，由setup_jacobian()函数生成的雅各布相关数据有望在若干次牛顿迭代中被使用。KINSOL自己决定何时重新生成雅各布系数和相关信息（比如为雅各布系数计算的预处理程序）是有益的，从而尽可能地节省重新生成雅各布系数矩阵和其预处理程序的努力。
     * @param  current_u  $u$   @param  current_f  $F(u)$  或  $G(u)$
     * 的当前值 这个函数应该返回。
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
     * - >0: 可恢复的错误（KINSOL将尝试改变其内部参数并尝试一个新的解决步骤
     *
     *
     *
     *
     *
     * - <0: 无法恢复的错误，计算将被中止，并抛出一个断言。
     *
     */
    std::function<int(const VectorType &current_u, const VectorType &current_f)>
      setup_jacobian;

    /**
     * @deprecated
     * 4.0之后的SUNDIALS版本不再提供这个回调的所有必要信息（见下文）。使用下面描述的`solve_with_jacobian`回调。
     * 一个用户可以提供的函数对象，其目的是用雅各布矩阵来解决一个线性系统。这个函数将在setup_jacobian()被调用至少一次后被KINSOL调用（可能是多次）。KINSOL试图尽最大努力调用setup_jacobian()最少的次数。如果不更新雅各布式就能达到收敛，那么KINSOL就不会再次调用setup_jacobian()。相反，如果KINSOL内部收敛测试失败，那么KINSOL就会用更新的向量和系数再次调用setup_jacobian()，以便连续调用solve_jacobian_systems()导致牛顿过程中更好的收敛。
     * 如果你没有指定`solve_jacobian_system`或`solve_with_jacobian`函数，那么只能使用固定点迭代策略。注意，这可能不会收敛，或者收敛得很慢。
     * 对这个函数的调用应该在`dst`中存储 $J^{-1}$
     * 应用于`rhs`的结果，即，`J*dst =
     * rhs`。用户有责任在这个函数中（或在上面的`setup_jacobian`回调中）设置适当的求解器和预处理器。
     * 该函数的参数是。          @param[in]  ycur
     * 当前KINSOL内部步骤的 $y$ 向量。在上面的文档中，这个
     * $y$ 向量一般用 $u$ 表示。      @param[in]  fcur
     * `ycur`处隐含的右手边的当前值，  $f_I (t_n, ypred)$  。
     * @param[in]  rhs 要解决的系统右手边  @param[out]  dst  $J^{-1}
     * src$  的解决方案 这个函数应该返回。
     *
     *
     *
     *
     * - 0: 成功
     *
     *
     *
     *
     * - >0: 可恢复的错误（KINSOL将尝试改变其内部参数并尝试一个新的解决步骤
     *
     *
     *
     *
     *
     * - <0: 无法恢复的错误，计算将被中止，并抛出一个断言。          @warning  从SUNDIALS 4.1开始，SUNDIALS不再提供`ycur'和`fcur'变量。
     *
     * - 只提供`rhs`，并且需要返回`dst`。因此在这种情况下，前两个参数将是空的向量。在实践中，这意味着我们不能再在这个函数中计算当前迭代的雅各布矩阵。相反，这必须发生在上面的`setup_jacobian`函数中，该函数接收这些信息。      如果雅各布矩阵与当前*迭代对应是很重要的（而不是一个重复的迭代）。
     * 迭代（而不是重新使用在之前的迭代中计算过的雅各布矩阵，因此对应于之前的*迭代），那么你还必须将
     * AdditionalData::maximum_newton_step
     * 变量设置为1，表示雅各布矩阵应该在每次迭代中重新计算。
     *
     */
    DEAL_II_DEPRECATED
    std::function<int(const VectorType &ycur,
                      const VectorType &fcur,
                      const VectorType &rhs,
                      VectorType &      dst)>
      solve_jacobian_system;

    /**
     * 一个用户可以提供的函数对象，其目的是用雅各布矩阵求解一个线性系统。这个函数将在setup_jacobian()被调用至少一次之后被KINSOL调用（可能是多次）。KINSOL试图尽最大努力调用setup_jacobian()最少的次数。如果不更新雅各布式就能达到收敛，那么KINSOL就不会再次调用setup_jacobian()。相反，如果KINSOL内部收敛测试失败，那么KINSOL就会用更新的向量和系数再次调用setup_jacobian()，以便连续调用solve_jacobian_systems()导致牛顿过程中更好的收敛。
     * 如果你没有指定`solve_with_jacobian`函数，那么只能使用固定点迭代策略。注意，这可能不会收敛，或者收敛得很慢。
     * 对这个函数的调用应该在`dst`中存储 $J^{-1}$
     * 应用于`rhs`的结果，即，`J*dst =
     * rhs`。用户有责任在这个函数中（或在上面的`setup_jacobian`回调中）设置适当的求解器和预处理器。附在这个回调上的函数还提供了对线性求解器的容忍度，表明不需要用雅各布矩阵精确地求解线性系统，而只需要达到KINSOL将随着时间推移而适应的容忍度。
     * 该函数的参数是。          @param[in]  rhs
     * 要解决的系统右侧。      @param[out]  dst  $J^{-1} src$
     * 的解。      @param[in]  tolerance
     * 用来解决线性方程组的公差。
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
     * - >0: 可恢复的错误（KINSOL将尝试改变其内部参数并尝试一个新的解决步骤
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
     * 一个用户可以提供的函数对象，目的是返回一个向量，其组成部分是KINSOL用来计算解决方案的向量规范的权重。这个函数的实现是可选的，只有在实现时才会使用。
     * 这个比例因子的目的是为了解决解的不同成分在数值上有很大差异的问题。
     *
     * - 通常是因为它们有不同的物理单位和代表不同的东西。例如，如果要解决一个非线性斯托克斯问题，解的向量有对应于速度的成分和对应于压力的其他成分。这些都有不同的物理单位，根据人们选择的单位，它们可能有大致可比的数字大小，也可能没有。仅举一个例子，在模拟地球内部的流动时，速度可能是每年10厘米，压力大约是100GPa。如果用国际单位来表示，这相当于大约 $0.000,000,003=3 \times 10^{-9}$ 米/秒的速度，和大约 $10^9 \text{kg}/\text{m}/\text{s}^2$ 的压力，也就是说，差别很大。在这种情况下，计算解决方案类型向量的 $l_2$ 准则（例如，前一个解决方案和当前解决方案之间的差异）是没有意义的，因为准则将被速度分量或压力分量所支配。该函数返回的缩放向量旨在为解的每个分量提供一个缩放因子，该因子通常被选为 "典型速度 "或 "典型压力 "的倒数，这样，当一个向量分量乘以相应的缩放向量分量时，可以得到一个数量级为1的数字（即典型速度/压力的1倍的合理小数）。KINSOL手册中对此有如下说明。"用户应提供数值 $D_u$ ，这些数值是缩放矩阵的对角线元素，当 $U$ 接近解时， $D_u U$ 的所有分量的大小大致相同"。        如果没有向KINSOL对象提供任何函数，那么这将被解释为隐含地表示所有这些缩放因子都应被视为一个。
     *
     */
    std::function<VectorType &()> get_solution_scaling;

    /**
     * 一个用户可以提供的函数对象，其目的是返回一个向量，其组成部分是KINSOL用来计算远离解的函数评估的向量规范的权重。这个函数的实现是可选的，只有在实现时才会使用。
     * 这个函数的要点和它返回的缩放向量与上面讨论的`get_solution_scaling'相似，只是在计算规范时，它是对函数
     * $F(U)$ 的分量进行缩放，而不是 $U$
     * 的分量。如上所述，如果没有提供函数，那么这就相当于使用一个缩放向量，其分量都等于1。
     *
     */
    std::function<VectorType &()> get_function_scaling;

    /**
     * 处理KINSOL异常。
     *
     */
    DeclException1(ExcKINSOLError,
                   int,
                   << "One of the SUNDIALS KINSOL internal functions "
                   << "returned a negative error code: " << arg1
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
     * KINSOL配置数据。
     *
     */
    AdditionalData data;

    /**
     * KINSOL内存对象。
     *
     */
    void *kinsol_mem;

    /**
     * 向量的内存池。
     *
     */
    GrowingVectorMemory<VectorType> mem;
  };

} // namespace SUNDIALS


DEAL_II_NAMESPACE_CLOSE

#endif

#endif


