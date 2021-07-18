//include/deal.II-translator/base/time_stepping_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2014 - 2020 by the deal.II authors
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

#ifndef dealii_time_stepping_h
#define dealii_time_stepping_h


#include <deal.II/base/config.h>

#include <deal.II/base/signaling_nan.h>

#include <functional>
#include <vector>

DEAL_II_NAMESPACE_OPEN

/**
 * 包含时间步进方法的命名空间。
 *
 *
 */

namespace TimeStepping
{
  /**
   * 下列Runge-Kutta方法可用。
   *
   *
   *
   *
   *
   *
   * - 明确的方法（见 ExplicitRungeKutta::initialize): ）。
   *
   *
   *
   *
   *
   *
   * - FORWARD_EULER (第一顺序)
   *
   *
   *
   * - RK_THIRD_ORDER (三阶 Runge-Kutta)
   *
   * - SSP_THIRD_ORDER (三阶SSP Runge-Kutta)
   *
   * - RK_CLASSIC_FOURTH_ORDER (经典四阶Runge-Kutta)
   *
   *
   * - 低存储（显式）Runge-Kutta方法
   *
   *
   * - LOW_STORAGE_RK_STAGE3_ORDER3 (三段式和三阶式)
   *
   *
   *
   * - LOW_STORAGE_RK_STAGE5_ORDER4 (五级四阶)
   *
   *
   *
   * - LOW_STORAGE_RK_STAGE7_ORDER4 (七级四阶)
   *
   *
   *
   * - LOW_STORAGE_RK_STAGE9_ORDER5 (九级五阶)
   *
   *
   * - 隐式方法（见 ImplicitRungeKutta::initialize): ）。
   *
   *
   *
   *
   *
   *
   * - BACKWARD_EULER (一阶)
   *
   *
   *
   * - IMPLICIT_MIDPOINT (二阶)
   *
   * - CRANK_NICOLSON (二阶)
   *
   *
   *
   *
   *
   * - SDIRK_TWO_STAGES (二阶)
   *
   *
   *
   * - 嵌入式显式方法（见 EmbeddedExplicitRungeKutta::initialize): ）。
   *
   *
   *
   *
   *
   * - HEUN_EULER (二阶)
   *
   *
   *
   * - BOGACKI_SHAMPINE (三阶)
   *
   * - DOPRI (Dormand-Prince method, fifth order; this is the method used by ode45 in MATLAB)
   *
   *
   *
   * - FEHLBERG (第五阶)
   * - CASH_KARP (第五阶)
   *
   */
  enum runge_kutta_method
  {
    /**
     * 正向欧拉法，一阶。
     *
     */
    FORWARD_EULER,
    /**
     * 三阶Runge-Kutta方法。
     *
     */
    RK_THIRD_ORDER,
    /**
     * 三阶强稳定（SSP）Runge-Kutta方法（SSP时间离散化在文献中也称为总变差递减（TVD）方法，见
     * @cite gottlieb2001strong ）。
     *
     */
    SSP_THIRD_ORDER,
    /**
     * 经典的四阶Runge-Kutta方法。
     *
     */
    RK_CLASSIC_FOURTH_ORDER,
    /**
     * Kennedy等人的三阶方案  @cite KennedyCarpenterLewis2000
     * 。它的稳定区域明显小于高阶方案，但由于只有三个阶段，它在每个阶段的工作方面非常有竞争力。
     *
     */
    LOW_STORAGE_RK_STAGE3_ORDER3,
    /**
     * 四阶的五级方案，在Kennedy等人的论文中定义  @cite
     * KennedyCarpenterLewis2000  。
     *
     */
    LOW_STORAGE_RK_STAGE5_ORDER4,
    /**
     * Tselios和Simos的论文中定义的四阶七级方案  @cite
     * TseliosSimos2007  。
     *
     */
    LOW_STORAGE_RK_STAGE7_ORDER4,
    /**
     * 肯尼迪等人的论文中定义的五阶九段式方案  @cite
     * KennedyCarpenterLewis2000  。
     *
     */
    LOW_STORAGE_RK_STAGE9_ORDER5,
    /**
     * 后退欧拉法，一阶。
     *
     */
    BACKWARD_EULER,
    /**
     * 隐式中点法，二阶。
     *
     */
    IMPLICIT_MIDPOINT,
    /**
     * Crank-Nicolson方法，二阶。
     *
     */
    CRANK_NICOLSON,
    /**
     * 两阶段SDIRK方法（"单对角隐式Runge-Kutta
     * "的简称），二阶。
     *
     */
    SDIRK_TWO_STAGES,
    /**
     * Heun方法（改进的欧拉方法），二阶。
     *
     */
    HEUN_EULER,
    /**
     * Bogacki-Shampine方法，三阶。
     *
     */
    BOGACKI_SHAMPINE,
    /**
     * Dormand-Prince方法，五阶；这是MATLAB中ODE45所使用的方法。
     *
     */
    DOPRI,
    /**
     * Fehlberg方法，五阶。
     *
     */
    FEHLBERG,
    /**
     * Cash-Karp方法，五阶。
     *
     */
    CASH_KARP,
    /**
     * 无效。
     *
     */
    invalid
  };



  /**
   * 使用嵌入式方法时退出evolve_one_time_step的原因。  delta_t,
   * min_delta_t, max_delta_t.
   *
   */
  enum embedded_runge_kutta_time_step
  {
    /**
     * 时间步长在有效范围内。
     *
     */
    DELTA_T,
    /**
     * 时间步长被增加到可接受的最小时间步长。
     *
     */
    MIN_DELTA_T,
    /**
     * 时间步长被减少到可接受的最大时间步长。
     *
     */
    MAX_DELTA_T
  };



  /**
   * 时间步长方法的抽象类。这些方法假定方程的形式。  $
   * \frac{\partial y}{\partial t} = f(t,y) $  .
   *
   */
  template <typename VectorType>
  class TimeStepping
  {
  public:
    /**
     * 虚拟解构器。
     *
     */
    virtual ~TimeStepping() = default;

    /**
     * 纯粹的虚拟函数。这个函数用于从时间 @p t推进到t+  @p
     * delta_t.   @p F 是一个应该被整合的函数 $ f(t,y) $
     * 的向量，输入参数是时间t和向量y，输出是此时的f值。
     * @p J_inverse
     * 是一个计算与隐式问题相关的雅各布式的逆函数的向量。输入参数是时间，
     * $ \tau $
     * 和一个矢量。输出是函数在这一点上的值。该函数返回时间步数结束时的时间。
     *
     */
    virtual double
    evolve_one_time_step(
      std::vector<std::function<VectorType(const double, const VectorType &)>>
        &                                                             F,
      std::vector<std::function<
        VectorType(const double, const double, const VectorType &)>> &J_inverse,
      double                                                          t,
      double                                                          delta_t,
      VectorType &                                                    y) = 0;

    /**
     * 用来存储信息的空结构。
     *
     */
    struct Status
    {};

    /**
     * 纯粹的虚拟函数，返回状态。
     *
     */
    virtual const Status &
    get_status() const = 0;
  };



  /**
   * Runge-Kutta方法的基类
   *
   */
  template <typename VectorType>
  class RungeKutta : public TimeStepping<VectorType>
  {
  public:
    /**
     * 虚拟解构器。
     *
     */
    virtual ~RungeKutta() override = default;

    /**
     * 用于初始化Runge-Kutta方法的纯虚拟方法。
     *
     */
    virtual void
    initialize(const runge_kutta_method method) = 0;

    /**
     * 该函数用于从时间 @p t 推进到t+ @p delta_t.   @p F
     * 是一个应被积分的函数 $ f(t,y) $
     * 的向量，输入参数是时间t和向量y，输出是f在此点的值。
     * @p J_inverse
     * 是一个计算与隐式问题相关的雅各布式的逆函数的向量。输入参数是时间，
     * $ \tau $  ，和一个向量。
     * 输出是函数在这一点上的值。该函数返回时间步数结束时的时间。当使用Runge-Kutta方法时，
     * @p F 和@J_inverse只能包含一个元素。
     *
     */
    double
    evolve_one_time_step(
      std::vector<std::function<VectorType(const double, const VectorType &)>>
        &                                                             F,
      std::vector<std::function<
        VectorType(const double, const double, const VectorType &)>> &J_inverse,
      double                                                          t,
      double                                                          delta_t,
      VectorType &y) override;

    /**
     * 纯粹的虚拟函数。这个函数用于从时间 @p t推进到t+ @p
     * delta_t.   @p f 是应该被积分的函数 $ f(t,y) $
     * ，输入参数是时间t和向量y，输出是f在这一点上的值。
     * @p id_minus_tau_J_inverse 是一个计算 $ inv(I-\tau J)$
     * 的函数，其中 $ I $ 是身份矩阵， $ \tau $ 是给定的， $
     * J $ 是雅各布 $ \frac{\partial J}{\partial y} $
     * 。输入参数是时间， $ \tau $
     * ，和一个矢量。输出是函数在这一点上的值。
     * evolve_one_time_step返回时间步数结束时的时间。
     *
     */
    virtual double
    evolve_one_time_step(
      const std::function<VectorType(const double, const VectorType &)> &f,
      const std::function<
        VectorType(const double, const double, const VectorType &)>
        &         id_minus_tau_J_inverse,
      double      t,
      double      delta_t,
      VectorType &y) = 0;

  protected:
    /**
     * Runge-Kutta方法的阶段数。
     *
     */
    unsigned int n_stages;

    /**
     * Butcher tableau系数。
     *
     */
    std::vector<double> b;

    /**
     * 屠夫 tableau 系数。
     *
     */
    std::vector<double> c;

    /**
     * 屠夫 tableau 系数。
     *
     */
    std::vector<std::vector<double>> a;
  };



  /**
   * ExplicitRungeKutta派生自RungeKutta并实现显式方法。
   *
   */
  template <typename VectorType>
  class ExplicitRungeKutta : public RungeKutta<VectorType>
  {
  public:
    using RungeKutta<VectorType>::evolve_one_time_step;

    /**
     * 默认构造函数。这个构造函数创建了一个对象，在它被使用之前，你要为它调用
     * <code>initialize(runge_kutta_method)</code> 。
     *
     */
    ExplicitRungeKutta() = default;

    /**
     * 构造函数。这个函数调用initialize(runge_kutta_method)。
     *
     */
    ExplicitRungeKutta(const runge_kutta_method method);

    /**
     * 初始化显式Runge-Kutta方法。
     *
     */
    void
    initialize(const runge_kutta_method method) override;

    /**
     * 该函数用于从时间 @p t 推进到t+ @p delta_t.   @p f
     * 是应该被积分的函数 $ f(t,y) $
     * ，输入参数是时间t和矢量y，输出是此时的f值。  @p
     * id_minus_tau_J_inverse 是一个计算 $ inv(I-\tau J)$
     * 的函数，其中 $ I $ 是身份矩阵， $ \tau $ 是给定的， $
     * J $ 是雅各布 $ \frac{\partial J}{\partial y} $
     * 。输入参数是时间， $ \tau $
     * ，和一个矢量。输出是函数在这一点上的值。
     * evolve_one_time_step返回时间步数结束时的时间。
     *
     */
    double
    evolve_one_time_step(
      const std::function<VectorType(const double, const VectorType &)> &f,
      const std::function<
        VectorType(const double, const double, const VectorType &)>
        &         id_minus_tau_J_inverse,
      double      t,
      double      delta_t,
      VectorType &y) override;

    /**
     * 这个函数用于从时间 @p t 推进到t+ @p delta_t.
     * 这个函数类似于从RungeKutta导出的函数，但不需要id_minus_tau_J_inverse，因为它不用于显式方法。
     * evolve_one_time_step 返回时间步长结束时的时间。
     *
     */
    double
    evolve_one_time_step(
      const std::function<VectorType(const double, const VectorType &)> &f,
      double                                                             t,
      double      delta_t,
      VectorType &y);

    /**
     * 该结构存储了所使用方法的名称。
     *
     */
    struct Status : public TimeStepping<VectorType>::Status
    {
      Status()
        : method(invalid)
      {}

      runge_kutta_method method;
    };

    /**
     * 返回当前对象的状态。
     *
     */
    const Status &
    get_status() const override;

  private:
    /**
     * 计算所需的不同阶段。
     *
     */
    void
    compute_stages(
      const std::function<VectorType(const double, const VectorType &)> &f,
      const double                                                       t,
      const double             delta_t,
      const VectorType &       y,
      std::vector<VectorType> &f_stages) const;

    /**
     * 对象的状态结构。
     *
     */
    Status status;
  };



  /**
   * LowStorageRungeKutta类派生于RungeKutta，实现了一类特殊的显式方法。低存储方法的主要优点是降低了内存消耗和减少了内存访问。
   *
   */
  template <typename VectorType>
  class LowStorageRungeKutta : public RungeKutta<VectorType>
  {
  public:
    using RungeKutta<VectorType>::evolve_one_time_step;

    /**
     * 默认构造函数。这个构造函数创建了一个对象，在使用该对象之前，你要先调用
     * <code>initialize(runge_kutta_method)</code> 。
     *
     */
    LowStorageRungeKutta() = default;

    /**
     * 构造函数。这个函数调用initialize(runge_kutta_method)。
     *
     */
    LowStorageRungeKutta(const runge_kutta_method method);

    /**
     * 初始化显式Runge-Kutta方法。
     *
     */
    void
    initialize(const runge_kutta_method method) override;

    /**
     * 该函数用于从时间 @p t 推进到t+ @p delta_t.   @p f
     * 是应该被积分的函数 $ f(t,y) $
     * ，输入参数是时间t和向量y，输出是此时的f值。  @p
     * id_minus_tau_J_inverse 是一个计算 $ inv(I-\tau J)$
     * 的函数，其中 $ I $ 是身份矩阵， $ \tau $ 是给定的， $
     * J $ 是雅各布 $ \frac{\partial J}{\partial y} $
     * 。输入参数是时间， $ \tau $
     * ，和一个矢量。输出是函数在这一点上的值。
     * evolve_one_time_step返回时间步数结束时的时间。
     *
     */
    double
    evolve_one_time_step(
      const std::function<VectorType(const double, const VectorType &)> &f,
      const std::function<
        VectorType(const double, const double, const VectorType &)>
        &         id_minus_tau_J_inverse,
      double      t,
      double      delta_t,
      VectorType &y) override;

    /**
     * 这个函数用于从时间 @p t 推进到t+ @p delta_t.
     * 这个函数类似于从RungeKutta导出的函数，但是不需要id_minus_tau_J_inverse，因为它不用于显式方法。
     * evolve_one_time_step返回时间步长结束时的时间。注意，vec_ki保存微分算子的评估，vec_ri保存微分算子应用的右手边。
     *
     */
    double
    evolve_one_time_step(
      const std::function<VectorType(const double, const VectorType &)> &f,
      double                                                             t,
      double      delta_t,
      VectorType &solution,
      VectorType &vec_ri,
      VectorType &vec_ki);

    /**
     * 获取该方案的系数。    注意这里的向量 @p a
     * 不是传统意义上的布彻表的定义，而只是其中的一个子对角线。更多细节可以在
     * step-67 和其中的参考文献中找到。
     *
     */
    void
    get_coefficients(std::vector<double> &a,
                     std::vector<double> &b,
                     std::vector<double> &c) const;

    /**
     * 该结构存储所使用的方法的名称。
     *
     */
    struct Status : public TimeStepping<VectorType>::Status
    {
      Status()
        : method(invalid)
      {}

      runge_kutta_method method;
    };

    /**
     * 返回当前对象的状态。
     *
     */
    const Status &
    get_status() const override;

  private:
    /**
     * 计算一个阶段的低存储量rk。
     *
     */
    void
    compute_one_stage(
      const std::function<VectorType(const double, const VectorType &)> &f,
      const double                                                       t,
      const double      factor_solution,
      const double      factor_ai,
      const VectorType &corrent_ri,
      VectorType &      vec_ki,
      VectorType &      solution,
      VectorType &      next_ri) const;

    /**
     * 对象的状态结构。
     *
     */
    Status status;
  };



  /**
   * 该类派生于RungeKutta并实现了隐式方法。
   * 这个类只对对角线隐式Runge-Kutta（DIRK）方法起作用。
   *
   */
  template <typename VectorType>
  class ImplicitRungeKutta : public RungeKutta<VectorType>
  {
  public:
    using RungeKutta<VectorType>::evolve_one_time_step;

    /**
     * 默认构造函数。初始化(runge_kutta_method)和set_newton_solver_parameters(unsigned
     * int,double)需要在使用该对象之前被调用。
     *
     */
    ImplicitRungeKutta() = default;

    /**
     * 构造函数。该函数调用initialize(runge_kutta_method)并初始化牛顿求解器的最大迭代次数和容忍度。
     *
     */
    ImplicitRungeKutta(const runge_kutta_method method,
                       const unsigned int       max_it    = 100,
                       const double             tolerance = 1e-6);

    /**
     * 初始化隐式Runge-Kutta方法。
     *
     */
    void
    initialize(const runge_kutta_method method) override;

    /**
     * 该函数用于从时间 @p t 推进到t+ @p delta_t.   @p f
     * 是应该被积分的函数 $ f(t,y) $
     * ，输入参数是时间t和矢量y，输出是此时的f值。  @p
     * id_minus_tau_J_inverse 是一个计算 $ (I-\tau J)^{-1}$
     * 的函数，其中 $ I $ 是身份矩阵， $ \tau $ 是给定的， $
     * J $ 是雅各布 $ \frac{\partial J}{\partial y} $
     * 。这个函数收到的输入参数是时间， $ \tau $
     * ，和一个向量。    输出是函数在这一点上的值。
     * evolve_one_time_step返回时间步数结束时的时间。
     *
     */
    double
    evolve_one_time_step(
      const std::function<VectorType(const double, const VectorType &)> &f,
      const std::function<
        VectorType(const double, const double, const VectorType &)>
        &         id_minus_tau_J_inverse,
      double      t,
      double      delta_t,
      VectorType &y) override;

    /**
     * 设置牛顿求解器使用的最大迭代次数和公差。
     *
     */
    void
    set_newton_solver_parameters(const unsigned int max_it,
                                 const double       tolerance);

    /**
     * 存储方法名称、牛顿迭代次数和退出牛顿求解器时的残差准则的结构。
     *
     */
    struct Status : public TimeStepping<VectorType>::Status
    {
      Status()
        : method(invalid)
        , n_iterations(numbers::invalid_unsigned_int)
        , norm_residual(numbers::signaling_nan<double>())
      {}

      runge_kutta_method method;
      unsigned int       n_iterations;
      double             norm_residual;
    };

    /**
     * 返回当前对象的状态。
     *
     */
    const Status &
    get_status() const override;

  private:
    /**
     * 计算所需的不同阶段。
     *
     */
    void
    compute_stages(
      const std::function<VectorType(const double, const VectorType &)> &f,
      const std::function<
        VectorType(const double, const double, const VectorType &)>
        &                      id_minus_tau_J_inverse,
      double                   t,
      double                   delta_t,
      VectorType &             y,
      std::vector<VectorType> &f_stages);

    /**
     * 用于隐含阶段的牛顿求解器。
     *
     */
    void
    newton_solve(
      const std::function<void(const VectorType &, VectorType &)> &get_residual,
      const std::function<VectorType(const VectorType &)>
        &         id_minus_tau_J_inverse,
      VectorType &y);

    /**
     * 计算牛顿求解器所需的残差。
     *
     */
    void
    compute_residual(
      const std::function<VectorType(const double, const VectorType &)> &f,
      double                                                             t,
      double            delta_t,
      const VectorType &new_y,
      const VectorType &y,
      VectorType &      tendency,
      VectorType &      residual) const;

    /**
     * 当使用SDIRK时，不需要计算各阶段的线性组合。因此，当这个标志为真时，线性组合被跳过。
     *
     */
    bool skip_linear_combi;

    /**
     * 牛顿求解器的最大迭代次数。
     *
     */
    unsigned int max_it;

    /**
     * 牛顿求解器的容忍度。
     *
     */
    double tolerance;

    /**
     * 对象的状态结构。
     *
     */
    Status status;
  };



  /**
   * 该类派生于RungeKutta，实现了嵌入式显式方法。
   *
   */
  template <typename VectorType>
  class EmbeddedExplicitRungeKutta : public RungeKutta<VectorType>
  {
  public:
    using RungeKutta<VectorType>::evolve_one_time_step;

    /**
     * 默认构造函数。初始化(runge_kutta_method)和set_time_adaptation_parameters(double,
     * double, double, double, double,
     * double)需要在使用该对象之前被调用。
     *
     */
    EmbeddedExplicitRungeKutta() = default;

    /**
     * 构造函数。这个函数调用initialize(runge_kutta_method)并初始化时间适应所需的参数。
     *
     */
    EmbeddedExplicitRungeKutta(const runge_kutta_method method,
                               const double             coarsen_param = 1.2,
                               const double             refine_param  = 0.8,
                               const double             min_delta     = 1e-14,
                               const double             max_delta     = 1e100,
                               const double             refine_tol    = 1e-8,
                               const double             coarsen_tol   = 1e-12);

    /**
     * 解构器。
     *
     */
    ~EmbeddedExplicitRungeKutta() override
    {
      free_memory();
    }

    /**
     * 如果有必要，删除对象分配的内存。
     *
     */
    void
    free_memory();

    /**
     * 初始化嵌入式显式Runge-Kutta方法。
     *
     */
    void
    initialize(const runge_kutta_method method) override;

    /**
     * 该函数用于从时间 @p t 推进到t+ @p delta_t.   @p f
     * 是应该被积分的函数 $ f(t,y) $
     * ，输入参数是时间t和矢量y，输出是此时的f值。  @p
     * id_minus_tau_J_inverse 是一个计算 $ inv(I-\tau J)$
     * 的函数，其中 $ I $ 是身份矩阵， $ \tau $ 是给定的， $
     * J $ 是雅各布系数 $ \frac{\partial J}{\partial y} $
     * 。输入参数是时间， $ \tau $
     * ，和一个矢量。输出是函数在这一点上的值。
     * evolve_one_time_step返回时间步数结束时的时间。
     *
     */
    double
    evolve_one_time_step(
      const std::function<VectorType(const double, const VectorType &)> &f,
      const std::function<
        VectorType(const double, const double, const VectorType &)>
        &         id_minus_tau_J_inverse,
      double      t,
      double      delta_t,
      VectorType &y) override;

    /**
     * 这个函数用于从时间 @p t 推进到t+ @p delta_t.
     * 这个函数类似于从TimeStepping导出的函数，但不需要id_minus_tau_J_inverse，因为它不用于显式方法。
     * evolve_one_time_step返回时间步数结束时的时间。
     *
     */
    double
    evolve_one_time_step(
      const std::function<VectorType(const double, const VectorType &)> &f,
      double                                                             t,
      double      delta_t,
      VectorType &y);

    /**
     * 设置时间适应的必要参数。
     *
     */
    void
    set_time_adaptation_parameters(const double coarsen_param,
                                   const double refine_param,
                                   const double min_delta,
                                   const double max_delta,
                                   const double refine_tol,
                                   const double coarsen_tol);

    /**
     * 存储方法名称的结构，退出evolve_one_time_step的原因，n_iterations里面的迭代次数，对下一个时间步长的猜测，以及对误差规范的估计。
     *
     */
    struct Status : public TimeStepping<VectorType>::Status
    {
      runge_kutta_method             method;
      embedded_runge_kutta_time_step exit_delta_t;
      unsigned int                   n_iterations;
      double                         delta_t_guess;
      double                         error_norm;
    };

    /**
     * 返回当前对象的状态。
     *
     */
    const Status &
    get_status() const override;

  private:
    /**
     * 计算所需的不同阶段。
     *
     */
    void
    compute_stages(
      const std::function<VectorType(const double, const VectorType &)> &f,
      const double                                                       t,
      const double             delta_t,
      const VectorType &       y,
      std::vector<VectorType> &f_stages);

    /**
     * 这个参数是时间步长可以粗化时，时间步长所乘的系数（>1）。
     *
     */
    double coarsen_param;

    /**
     * 该参数是当时间步长必须细化时，时间步长乘以的系数（<1）。
     *
     */
    double refine_param;

    /**
     * 允许的最小的时间步长。
     *
     */
    double min_delta_t;

    /**
     * 允许的最大的时间步长。
     *
     */
    double max_delta_t;

    /**
     * 细化容忍度：如果误差估计值大于refine_tol，则对时间步长进行细化。
     *
     */
    double refine_tol;

    /**
     * 粗化公差：如果误差估计值小于coarse_tol，则对时间步长进行粗化。
     *
     */
    double coarsen_tol;

    /**
     * 如果该标志为真，最后阶段与第一阶段相同，可以保存f的一次评估。
     *
     */
    bool last_same_as_first;

    /**
     * 屠夫 tableau 系数。
     *
     */
    std::vector<double> b1;

    /**
     * 屠夫 tableau 系数。
     *
     */
    std::vector<double> b2;

    /**
     * 如果last_same_as_first标志被设置为
     * "true"，最后一个阶段将被保存并作为下一个时间步骤的第一阶段重新使用。
     *
     */
    VectorType *last_stage;

    /**
     * 对象的状态结构。
     *
     */
    Status status;
  };
} // namespace TimeStepping

DEAL_II_NAMESPACE_CLOSE

#endif


