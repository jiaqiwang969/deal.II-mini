//include/deal.II-translator/lac/solver_bicgstab_0.txt
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

#ifndef dealii_solver_bicgstab_h
#define dealii_solver_bicgstab_h


#include <deal.II/base/config.h>

#include <deal.II/base/logstream.h>
#include <deal.II/base/signaling_nan.h>
#include <deal.II/base/subscriptor.h>

#include <deal.II/lac/solver.h>
#include <deal.II/lac/solver_control.h>

#include <cmath>

DEAL_II_NAMESPACE_OPEN

 /*!@addtogroup Solvers */ 
 /*@{*/ 

namespace internal
{
  /**
   * 包含SolverBicgstab类所使用的非参数非模板值的类。
   *
   */
  class SolverBicgstabData
  {
  protected:
    /**
     * 辅助值。
     *
     */
    double alpha;
    /**
     * 辅助值。
     *
     */
    double beta;
    /**
     * 辅助值。
     *
     */
    double omega;
    /**
     * 辅助值。
     *
     */
    double rho;
    /**
     * 辅助值。
     *
     */
    double rhobar;

    /**
     * 当前迭代步骤。
     *
     */
    unsigned int step;

    /**
     * 残差。
     *
     */
    double res;

    /**
     * 默认构造函数。这是受保护的，所以只有SolverBicgstab可以创建实例。
     *
     */
    SolverBicgstabData();
  };
} // namespace internal

/**
 * van der Vorst的Bicgstab算法。
 * 关于使用该类时对矩阵和向量的要求，请参见Solver基类的文档。
 * 像所有其他求解器类一样，这个类有一个叫做 @p
 * AdditionalData的局部结构，用来向求解器传递额外的参数，如阻尼参数或临时向量的数量。我们使用这个额外的结构，而不是直接将这些值传递给构造函数，因为这使得
 * @p SolverSelector
 * 和其他类的使用更加容易，并保证即使某个求解器的额外参数的数量或类型发生变化，这些也能继续工作。
 * Bicgstab方法有两个额外的参数：第一个是布尔值，决定是在每一步中计算实际的残差（
 * @p true)  还是使用计算的正交残差的长度（  @p false).
 * 注意，计算残差会在每一步中引起第三次矩阵-向量-乘法，尽管没有额外的预处理。这样做的原因是，在迭代过程中计算的正交残差的大小可能比真正的残差要大几个数量级。这是由于与条件不良的矩阵有关的数值不稳定性造成的。由于这种不稳定性导致了不好的停止标准，这个参数的默认值是
 * @p true.
 * 。每当用户知道估计的残差工作得很合理时，为了提高求解器的性能，这个标志应该设置为
 * @p false 。
 * 第二个参数是分解准则的大小。很难找到一个普遍的好准则，所以如果事情对你不适用，可以尝试改变这个值。
 *
 *  <h3>Observing the progress of linear solver iterations</h3>
 * 这个类的solve()函数使用Solver基类中描述的机制来确定收敛性。这个机制也可以用来观察迭代的进度。
 *
 *
 */
template <typename VectorType = Vector<double>>
class SolverBicgstab : public SolverBase<VectorType>,
                       protected internal::SolverBicgstabData
{
public:
  /**
   * 计算残差有两种可能：一种是使用计算值 @p tau.
   * 进行估计
   * 另一种是使用另一个矩阵向量乘法进行精确计算。这增加了算法的成本，所以只要问题允许，就应该将其设置为假。
   * Bicgstab很容易发生故障，所以我们需要一个参数来告诉我们哪些数字被认为是零。
   *
   */
  struct AdditionalData
  {
    /**
     * 构造函数。
     * 默认是进行精确的残差计算，分解参数是VectorType的value_type可表示的最小有限值。
     *
     */
    explicit AdditionalData(
      const bool   exact_residual = true,
      const double breakdown =
        std::numeric_limits<typename VectorType::value_type>::min())
      : exact_residual(exact_residual)
      , breakdown(breakdown)
    {}
    /**
     * 残差精确计算的标志。
     *
     */
    bool exact_residual;
    /**
     * 破解阈值。
     *
     */
    double breakdown;
  };

  /**
   * 构造函数。
   *
   */
  SolverBicgstab(SolverControl &           cn,
                 VectorMemory<VectorType> &mem,
                 const AdditionalData &    data = AdditionalData());

  /**
   * 构造器。使用一个GrowingVectorMemory类型的对象作为默认分配内存。
   *
   */
  SolverBicgstab(SolverControl &       cn,
                 const AdditionalData &data = AdditionalData());

  /**
   * 虚拟解构器。
   *
   */
  virtual ~SolverBicgstab() override = default;

  /**
   * 只解决原始问题。
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
   * 传递给solve()的解决方案矢量的指针。
   *
   */
  VectorType *Vx;

  /**
   * 辅助向量。
   *
   */
  typename VectorMemory<VectorType>::Pointer Vr;

  /**
   * 辅助向量。
   *
   */
  typename VectorMemory<VectorType>::Pointer Vrbar;

  /**
   * 辅助向量。
   *
   */
  typename VectorMemory<VectorType>::Pointer Vp;

  /**
   * 辅助向量。
   *
   */
  typename VectorMemory<VectorType>::Pointer Vy;

  /**
   * 辅助向量。
   *
   */
  typename VectorMemory<VectorType>::Pointer Vz;

  /**
   * 辅助向量。
   *
   */
  typename VectorMemory<VectorType>::Pointer Vt;

  /**
   * 辅助向量。
   *
   */
  typename VectorMemory<VectorType>::Pointer Vv;

  /**
   * 指向传递给solve()的右手边向量的指针。
   *
   */
  const VectorType *Vb;

  /**
   * 停止准则的计算。
   *
   */
  template <typename MatrixType>
  double
  criterion(const MatrixType &A, const VectorType &x, const VectorType &b);

  /**
   * 派生类的接口。
   * 这个函数在每一步中获得当前迭代向量、残差和更新向量。它可以用于收敛历史的图形输出。
   *
   */
  virtual void
  print_vectors(const unsigned int step,
                const VectorType & x,
                const VectorType & r,
                const VectorType & d) const;

  /**
   * 附加参数。
   *
   */
  AdditionalData additional_data;

private:
  /**
   * iterate()函数返回的一个结构，代表它发现在迭代过程中发生的事情。
   *
   */
  struct IterationResult
  {
    bool                 breakdown;
    SolverControl::State state;
    unsigned int         last_step;
    double               last_residual;

    IterationResult(const bool                 breakdown,
                    const SolverControl::State state,
                    const unsigned int         last_step,
                    const double               last_residual);
  };

  /**
   * 迭代循环本身。该函数返回一个结构，表示在这个函数中发生了什么。
   *
   */
  template <typename MatrixType, typename PreconditionerType>
  IterationResult
  iterate(const MatrixType &A, const PreconditionerType &preconditioner);
};

 /*@}*/ 
 /*-------------------------Inline functions -------------------------------*/ 

#ifndef DOXYGEN


template <typename VectorType>
SolverBicgstab<VectorType>::IterationResult::IterationResult(
  const bool                 breakdown,
  const SolverControl::State state,
  const unsigned int         last_step,
  const double               last_residual)
  : breakdown(breakdown)
  , state(state)
  , last_step(last_step)
  , last_residual(last_residual)
{}



template <typename VectorType>
SolverBicgstab<VectorType>::SolverBicgstab(SolverControl &           cn,
                                           VectorMemory<VectorType> &mem,
                                           const AdditionalData &    data)
  : SolverBase<VectorType>(cn, mem)
  , Vx(nullptr)
  , Vb(nullptr)
  , additional_data(data)
{}



template <typename VectorType>
SolverBicgstab<VectorType>::SolverBicgstab(SolverControl &       cn,
                                           const AdditionalData &data)
  : SolverBase<VectorType>(cn)
  , Vx(nullptr)
  , Vb(nullptr)
  , additional_data(data)
{}



template <typename VectorType>
template <typename MatrixType>
double
SolverBicgstab<VectorType>::criterion(const MatrixType &A,
                                      const VectorType &x,
                                      const VectorType &b)
{
  A.vmult(*Vt, x);
  Vt->add(-1., b);
  res = Vt->l2_norm();

  return res;
}



template <typename VectorType>
void
SolverBicgstab<VectorType>::print_vectors(const unsigned int,
                                          const VectorType &,
                                          const VectorType &,
                                          const VectorType &) const
{}



template <typename VectorType>
template <typename MatrixType, typename PreconditionerType>
typename SolverBicgstab<VectorType>::IterationResult
SolverBicgstab<VectorType>::iterate(const MatrixType &        A,
                                    const PreconditionerType &preconditioner)
{
  A.vmult(*Vr, *Vx);
  Vr->sadd(-1., 1., *Vb);
  res = Vr->l2_norm();

  SolverControl::State state = this->iteration_status(step, res, *Vx);
  if (state == SolverControl::State::success)
    return IterationResult(false, state, step, res);

  alpha = omega = rho = 1.;

  VectorType &r    = *Vr;
  VectorType &rbar = *Vrbar;
  VectorType &p    = *Vp;
  VectorType &y    = *Vy;
  VectorType &z    = *Vz;
  VectorType &t    = *Vt;
  VectorType &v    = *Vv;

  rbar         = r;
  bool startup = true;

  do
    {
      ++step;

      rhobar = r * rbar;
      if (std::fabs(rhobar) < additional_data.breakdown)
        {
          return IterationResult(true, state, step, res);
        }
      beta = rhobar * alpha / (rho * omega);
      rho  = rhobar;
      if (startup == true)
        {
          p       = r;
          startup = false;
        }
      else
        {
          p.sadd(beta, 1., r);
          p.add(-beta * omega, v);
        }

      preconditioner.vmult(y, p);
      A.vmult(v, y);
      rhobar = rbar * v;
      if (std::fabs(rhobar) < additional_data.breakdown)
        {
          return IterationResult(true, state, step, res);
        }

      alpha = rho / rhobar;

      res = std::sqrt(r.add_and_dot(-alpha, v, r));

      // check for early success, see the lac/bicgstab_early testcase as to
      // why this is necessary
      //
      // note: the vector *Vx we pass to the iteration_status signal here is
      // only the current approximation, not the one we will return with, which
      // will be x=*Vx + alpha*y
      if (this->iteration_status(step, res, *Vx) == SolverControl::success)
        {
          Vx->add(alpha, y);
          print_vectors(step, *Vx, r, y);
          return IterationResult(false, SolverControl::success, step, res);
        }

      preconditioner.vmult(z, r);
      A.vmult(t, z);
      rhobar         = t * r;
      auto t_squared = t * t;
      if (t_squared < additional_data.breakdown)
        {
          return IterationResult(true, state, step, res);
        }
      omega = rhobar / (t * t);
      Vx->add(alpha, y, omega, z);

      if (additional_data.exact_residual)
        {
          r.add(-omega, t);
          res = criterion(A, *Vx, *Vb);
        }
      else
        res = std::sqrt(r.add_and_dot(-omega, t, r));

      state = this->iteration_status(step, res, *Vx);
      print_vectors(step, *Vx, r, y);
    }
  while (state == SolverControl::iterate);

  return IterationResult(false, state, step, res);
}



template <typename VectorType>
template <typename MatrixType, typename PreconditionerType>
void
SolverBicgstab<VectorType>::solve(const MatrixType &        A,
                                  VectorType &              x,
                                  const VectorType &        b,
                                  const PreconditionerType &preconditioner)
{
  LogStream::Prefix prefix("Bicgstab");

  // Allocate temporary memory.
  Vr    = typename VectorMemory<VectorType>::Pointer(this->memory);
  Vrbar = typename VectorMemory<VectorType>::Pointer(this->memory);
  Vp    = typename VectorMemory<VectorType>::Pointer(this->memory);
  Vy    = typename VectorMemory<VectorType>::Pointer(this->memory);
  Vz    = typename VectorMemory<VectorType>::Pointer(this->memory);
  Vt    = typename VectorMemory<VectorType>::Pointer(this->memory);
  Vv    = typename VectorMemory<VectorType>::Pointer(this->memory);

  Vr->reinit(x, true);
  Vrbar->reinit(x, true);
  Vp->reinit(x, true);
  Vy->reinit(x, true);
  Vz->reinit(x, true);
  Vt->reinit(x, true);
  Vv->reinit(x, true);

  Vx = &x;
  Vb = &b;

  step = 0;

  IterationResult state(false, SolverControl::failure, 0, 0);
  do
    {
      state = iterate(A, preconditioner);
    }
  while (state.state == SolverControl::iterate);


  // Release the temporary memory again.
  Vr.reset();
  Vrbar.reset();
  Vp.reset();
  Vy.reset();
  Vz.reset();
  Vt.reset();
  Vv.reset();

  // In case of failure: throw exception
  AssertThrow(state.state == SolverControl::success,
              SolverControl::NoConvergence(state.last_step,
                                           state.last_residual));
  // Otherwise exit as normal
}

#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif


