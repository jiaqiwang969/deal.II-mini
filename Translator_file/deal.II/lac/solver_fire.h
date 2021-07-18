//include/deal.II-translator/lac/solver_fire_0.txt
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

#ifndef dealii_solver_fire_h
#define dealii_solver_fire_h


#include <deal.II/base/config.h>

#include <deal.II/base/logstream.h>

#include <deal.II/lac/diagonal_matrix.h>
#include <deal.II/lac/solver.h>

#include <functional>


DEAL_II_NAMESPACE_OPEN


 /*!@addtogroup Solvers */ 
 /*@{*/ 

/**
 * FIRE（快速惯性放松引擎）用于最小化（潜在的非线性）目标函数
 * $E(\mathbf x)$  ，  $\mathbf x$  是一个  $n$  变量的向量（  $n$
 * 是目标函数的变量数）。像所有其他求解器类一样，只要满足一定的要求，它可以在任何类型的向量和矩阵上工作（关于使用该类的矩阵和向量的要求，请参见求解器基类的文档）。求解向量的类型必须作为模板参数传递，默认为
 * dealii::Vector<double>.  。 FIRE是Bitzek等人在<a
 * href="https://doi.org/10.1103/PhysRevLett.97.170201">Structural Relaxation
 * Made
 * Simple</a>中描述的一种阻尼动力学方法，通常用于寻找计算材料科学中原子系统的稳定平衡配置。从给定的原子系统初始配置开始，该算法依靠惯性来获得势能最小的（最近的）配置。
 * 符号。
 *
 *
 *
 *
 * - 未知变量的全局向量。  $\mathbf x$  .
 *
 *
 *
 *
 * - 目标函数。                      $E(\mathbf x)$  .
 *
 *
 *
 *
 * - 未知数的变化率。              $\mathbf v$  .
 *
 *
 *
 *
 * - 目标函数的梯度与未知数的关系。                 $\mathbf g = \nabla E(\mathbf x)$  .
 *
 *
 *
 *
 * - 质量矩阵。                             $\mathbf M$  .
 *
 *
 *
 *
 *
 * - 未知数的初始猜测。               $\mathbf x_0$  .
 *
 *
 *
 *
 *
 * - 时间步骤。                               $\Delta t$  .
 * 给出 $\Delta t$ 、 $\alpha = \alpha_0$ 、 $\epsilon$ 、 $\mathbf x =
 * \mathbf x_0$ 和 $\mathbf v= \mathbf 0$
 * 的初始值以及给定的质量矩阵 $\mathbf M$
 * ，FIRE算法如下：1.计算 $\mathbf g = \nabla E(\mathbf x)$
 * 并检查收敛情况（ $\mathbf g \cdot \mathbf g < \epsilon^2 $ ）。2.
 * 使用简单的（正向）欧拉积分步骤更新 $\mathbf x$ 和 $V$
 * ，<BR>  $\mathbf x = \mathbf x + \Delta t \mathbf v$  ，<BR>  $\mathbf v
 * = \mathbf v + \Delta t \mathbf M^{-1} \cdot \mathbf g$  。3. 计算 $p =
 * \mathbf g \cdot \mathbf v$  。4. 设置 $\mathbf v = (1-\alpha) \mathbf v
 * + \alpha \frac{|\mathbf v|}{|\mathbf g|} \mathbf g$  。5. 如果 $p<0$
 * 和自 $p$
 * 最后一次为负数以来的步数大于特定值，则增加时间步数
 * $\Delta t$ 并减少 $\alpha$  。6. 如果  $p>0$
 * ，则减少时间步长，冻结系统，即  $\mathbf v = \mathbf 0$
 * 并重置  $\alpha = \alpha_0$  。7. 返回到1。 也见Eidel等人的<a
 * href="http://onlinelibrary.wiley.com/doi/10.1002/pamm.201110246/full">
 * Energy-Minimization in Atomic-to-Continuum Scale-Bridging Methods
 * </a>，2011年。
 *
 *
 */
template <typename VectorType = Vector<double>>
class SolverFIRE : public SolverBase<VectorType>
{
public:
  /**
   * 标准化的数据结构，用于向求解器输送额外的数据。
   *
   */
  struct AdditionalData
  {
    /**
     * 构造函数。默认情况下，将（正向）欧拉积分步骤的初始时间步长设置为0.1，最大时间步长设置为1，任何变量的最大允许变化（每次迭代）为1。
     *
     */
    explicit AdditionalData(const double initial_timestep    = 0.1,
                            const double maximum_timestep    = 1,
                            const double maximum_linfty_norm = 1);

    /**
     * (正向)欧拉积分步骤的初始时间步长。
     *
     */
    const double initial_timestep;

    /**
     * (正向)欧拉积分步骤的最大时间步长。
     *
     */
    const double maximum_timestep;

    /**
     * 允许目标函数的任何变量的最大变化。
     *
     */
    const double maximum_linfty_norm;
  };

  /**
   * 构造器。
   *
   */
  SolverFIRE(SolverControl &           solver_control,
             VectorMemory<VectorType> &vector_memory,
             const AdditionalData &    data = AdditionalData());

  /**
   * 构造函数。使用一个GrowingVectorMemory类型的对象作为默认分配内存。
   *
   */
  SolverFIRE(SolverControl &       solver_control,
             const AdditionalData &data = AdditionalData());

  /**
   * 获得一组变量 @p x ，使多态函数包装器 @p compute,
   * 描述的目标函数最小化，并具有给定的预处理器 @p
   * inverse_mass_matrix 和初始 @p x 值。  函数 @p compute
   * 返回目标函数的值，并根据第二个参数--变量的状态，作为第一个参数传入时，更新目标函数的梯度（相对于变量）。
   *
   */
  template <typename PreconditionerType = DiagonalMatrix<VectorType>>
  void
  solve(const std::function<double(VectorType &, const VectorType &)> &compute,
        VectorType &                                                   x,
        const PreconditionerType &inverse_mass_matrix);

  /**
   * 当 $E(\mathbf x) = \frac{1}{2} \mathbf x^{T} \mathbf A \mathbf x
   *
   * - \mathbf x^{T} \mathbf b$ 时，求解x，使 $E(\mathbf x)$
   * 的<EM>特殊情况</EM>最小。
   *
   */
  template <typename MatrixType, typename PreconditionerType>
  void
  solve(const MatrixType &        A,
        VectorType &              x,
        const VectorType &        b,
        const PreconditionerType &preconditioner);

protected:
  /**
   * 派生类的接口。这个函数在每个步骤中获得当前迭代 @p
   * x （变量）， @p v （x的时间导数）和 @p g （梯度）。
   * 它可以用于收敛历史的图形输出。
   *
   */
  virtual void
  print_vectors(const unsigned int,
                const VectorType &x,
                const VectorType &v,
                const VectorType &g) const;

  /**
   * 解算器的附加数据。
   *
   */
  const AdditionalData additional_data;
};

 /*@}*/ 

 /*------------------------- Implementation ----------------------------*/ 

#ifndef DOXYGEN

template <typename VectorType>
SolverFIRE<VectorType>::AdditionalData::AdditionalData(
  const double initial_timestep,
  const double maximum_timestep,
  const double maximum_linfty_norm)
  : initial_timestep(initial_timestep)
  , maximum_timestep(maximum_timestep)
  , maximum_linfty_norm(maximum_linfty_norm)
{
  AssertThrow(initial_timestep > 0. && maximum_timestep > 0. &&
                maximum_linfty_norm > 0.,
              ExcMessage("Expected positive values for initial_timestep, "
                         "maximum_timestep and maximum_linfty_norm but one "
                         "or more of the these values are not positive."));
}



template <typename VectorType>
SolverFIRE<VectorType>::SolverFIRE(SolverControl &           solver_control,
                                   VectorMemory<VectorType> &vector_memory,
                                   const AdditionalData &    data)
  : SolverBase<VectorType>(solver_control, vector_memory)
  , additional_data(data)
{}



template <typename VectorType>
SolverFIRE<VectorType>::SolverFIRE(SolverControl &       solver_control,
                                   const AdditionalData &data)
  : SolverBase<VectorType>(solver_control)
  , additional_data(data)
{}



template <typename VectorType>
template <typename PreconditionerType>
void
SolverFIRE<VectorType>::solve(
  const std::function<double(VectorType &, const VectorType &)> &compute,
  VectorType &                                                   x,
  const PreconditionerType &inverse_mass_matrix)
{
  LogStream::Prefix prefix("FIRE");

  // FIRE algorithm constants
  const double DELAYSTEP       = 5;
  const double TIMESTEP_GROW   = 1.1;
  const double TIMESTEP_SHRINK = 0.5;
  const double ALPHA_0         = 0.1;
  const double ALPHA_SHRINK    = 0.99;

  using real_type = typename VectorType::real_type;

  typename VectorMemory<VectorType>::Pointer v(this->memory);
  typename VectorMemory<VectorType>::Pointer g(this->memory);

  // Set velocities to zero but not gradients
  // as we are going to compute them soon.
  v->reinit(x, false);
  g->reinit(x, true);

  // Refer to v and g with some readable names.
  VectorType &velocities = *v;
  VectorType &gradients  = *g;

  // Update gradients for the new x.
  compute(gradients, x);

  unsigned int iter = 0;

  SolverControl::State conv = SolverControl::iterate;
  conv = this->iteration_status(iter, gradients * gradients, x);
  if (conv != SolverControl::iterate)
    return;

  // Refer to additional data members with some readable names.
  const auto &maximum_timestep = additional_data.maximum_timestep;
  double      timestep         = additional_data.initial_timestep;

  // First scaling factor.
  double alpha = ALPHA_0;

  unsigned int previous_iter_with_positive_v_dot_g = 0;

  while (conv == SolverControl::iterate)
    {
      ++iter;
      // Euler integration step.
      x.add(timestep, velocities);                     // x += dt     * v
      inverse_mass_matrix.vmult(gradients, gradients); // g  = M^{-1} * g
      velocities.add(-timestep, gradients);            // v -= dt     * h

      // Compute gradients for the new x.
      compute(gradients, x);

      const real_type gradient_norm_squared = gradients * gradients;
      conv = this->iteration_status(iter, gradient_norm_squared, x);
      if (conv != SolverControl::iterate)
        break;

      // v_dot_g = V * G
      const real_type v_dot_g = velocities * gradients;

      if (v_dot_g < 0.)
        {
          const real_type velocities_norm_squared = velocities * velocities;

          // Check if we divide by zero in DEBUG mode.
          Assert(gradient_norm_squared > 0., ExcInternalError());

          // beta = - alpha |V|/|G|
          const real_type beta =
            -alpha * std::sqrt(velocities_norm_squared / gradient_norm_squared);

          // V = (1-alpha) V + beta G.
          velocities.sadd(1. - alpha, beta, gradients);

          if (iter - previous_iter_with_positive_v_dot_g > DELAYSTEP)
            {
              // Increase timestep and decrease alpha.
              timestep = std::min(timestep * TIMESTEP_GROW, maximum_timestep);
              alpha *= ALPHA_SHRINK;
            }
        }
      else
        {
          // Decrease timestep, reset alpha and set V = 0.
          previous_iter_with_positive_v_dot_g = iter;
          timestep *= TIMESTEP_SHRINK;
          alpha      = ALPHA_0;
          velocities = 0.;
        }

      real_type vmax = velocities.linfty_norm();

      // Change timestep if any dof would move more than maximum_linfty_norm.
      if (vmax > 0.)
        {
          const double minimal_timestep =
            additional_data.maximum_linfty_norm / vmax;
          if (minimal_timestep < timestep)
            timestep = minimal_timestep;
        }

      print_vectors(iter, x, velocities, gradients);

    } // While we need to iterate.

  // In the case of failure: throw exception.
  if (conv != SolverControl::success)
    AssertThrow(false,
                SolverControl::NoConvergence(iter, gradients * gradients));
}



template <typename VectorType>
template <typename MatrixType, typename PreconditionerType>
void
SolverFIRE<VectorType>::solve(const MatrixType &        A,
                              VectorType &              x,
                              const VectorType &        b,
                              const PreconditionerType &preconditioner)
{
  std::function<double(VectorType &, const VectorType &)> compute_func =
    [&](VectorType &g, const VectorType &x) -> double {
    // Residual of the quadratic form $ \frac{1}{2} xAx - xb $.
    // G = b - Ax
    A.residual(g, x, b);

    // Gradient G = Ax -b.
    g *= -1.;

    // The quadratic form $\frac{1}{2} xAx - xb $.
    return 0.5 * A.matrix_norm_square(x) - x * b;
  };

  this->solve(compute_func, x, preconditioner);
}



template <typename VectorType>
void
SolverFIRE<VectorType>::print_vectors(const unsigned int,
                                      const VectorType &,
                                      const VectorType &,
                                      const VectorType &) const
{}



#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif


