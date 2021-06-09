//include/deal.II-translator/sundials/sunlinsol_wrapper_0.txt
//-----------------------------------------------------------
//
//    Copyright (C) 2021 by the deal.II authors
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
//-----------------------------------------------------------

#ifndef dealii_sundials_sunlinsol_wrapper_h
#define dealii_sundials_sunlinsol_wrapper_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_SUNDIALS
#  if DEAL_II_SUNDIALS_VERSION_GTE(4, 0, 0)

#    include <sundials/sundials_linearsolver.h>

#    include <functional>
#    include <memory>

DEAL_II_NAMESPACE_OPEN

namespace SUNDIALS
{
#    ifndef DOXYGEN
  // forward declarations
  namespace internal
  {
    template <typename VectorType>
    struct LinearSolverContent;
  }
#    endif

  /**
   * 一个包裹SUNDIALS功能的线性运算符。
   *
   */
  template <typename VectorType>
  struct SundialsOperator
  {
    /**
     * 将此线性运算符应用于 @p src 并将结果存储在 @p dst.
     * 中。
     *
     */
    void
    vmult(VectorType &dst, const VectorType &src) const;

    /**
     * 构造函数。          @param  A_data  @p a_times_fn 所需的数据
     * @param  a_times_fn 一个指向计算A*v的函数的指针。
     *
     */
    SundialsOperator(void *A_data, ATimesFn a_times_fn);

  private:
    /**
     * 评估a_times_fn的必要数据。
     *
     */
    void *A_data;

    /**
     * 由SUNDIALS声明的%函数指针，用于评估矩阵向量乘积。
     *
     */
    ATimesFn a_times_fn;
  };



  /**
   * 由SUNDIALS指定的包裹预处理功能的线性运算器。vmult()函数解决预处理方程
   * $Px=b$ ，即计算 $x=P^{-1}b$  。
   *
   */
  template <typename VectorType>
  struct SundialsPreconditioner
  {
    /**
     * 应用包裹的预处理程序，即解出  $Px=b$  ，其中  $x$  是
     * @p dst  矢量，  $b$  是  @p src  矢量。          @param  dst
     * 前处理程序应用的结果向量  @param  src
     * 前处理程序应用的目标向量
     *
     */
    void
    vmult(VectorType &dst, const VectorType &src) const;

    /**
     * 构造函数。          @param  P_data  @p p_solve_fn 所需的数据
     * @param  p_solve_fn 计算A*v的函数指针  @param  tol
     * 迭代求解器用来判断收敛的公差
     *
     */
    SundialsPreconditioner(void *P_data, PSolveFn p_solve_fn, double tol);

  private:
    /**
     * 调用p_solve_fn的必要数据
     *
     */
    void *P_data;

    /**
     * %函数指针，用于计算预调节器应用的函数。
     *
     */
    PSolveFn p_solve_fn;

    /**
     * 潜在的容忍度，用于预处理程序方程的内部求解。
     *
     */
    double tol;
  };

  /**
   * 与SUNDIALS线性求解器接口的函数对象类型
   * 该函数类型封装了求解的动作  $P^{-1}Ax=P^{-1}b$  。
   * LinearOperator  @p op  封装了矩阵向量乘积  $Ax$
   * ，LinearOperator  @p prec  封装了前置条件器的应用  $P^{-1}z$
   * 。
   * 用户可以指定这种类型的函数对象，将自定义的线性求解器例程附加到SUNDIALS。两个LinearOperators
   * @p op  和  @p prec
   * 是由SUNDIALS根据用户的设置在内部建立的。参数的解释如下。
   * @param[in]  op 一个应用矩阵向量乘积的LinearOperator
   * @param[in]  prec 一个应用预处理的LinearOperator  @param[out]  x
   * 输出的解向量  @param[in]  b 右手边  @param[in]  tol
   * 迭代求解器的公差 这个函数应该返回。
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
  template <typename VectorType>
  using LinearSolveFunction =
    std::function<int(SundialsOperator<VectorType> &      op,
                      SundialsPreconditioner<VectorType> &prec,
                      VectorType &                        x,
                      const VectorType &                  b,
                      double                              tol)>;

  namespace internal
  {
    /*!     为SUNDIALS的线性求解器接口附加包装函数。我们假装用户提供的线性求解器是无矩阵的，尽管它可以是基于矩阵的。这样，SUNDIALS就不需要理解我们的矩阵类型。   
*
*/
    template <typename VectorType>
    class LinearSolverWrapper
    {
    public:
      explicit LinearSolverWrapper(LinearSolveFunction<VectorType> lsolve);

      ~LinearSolverWrapper();

      /**
       * 隐式转换为SUNLinearSolver。
       *
       */
      operator SUNLinearSolver();

    private:
      SUNLinearSolver                                  sun_linear_solver;
      std::unique_ptr<LinearSolverContent<VectorType>> content;
    };
  } // namespace internal
} // namespace SUNDIALS

DEAL_II_NAMESPACE_CLOSE

#  endif
#endif
#endif

