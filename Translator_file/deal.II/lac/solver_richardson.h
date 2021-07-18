//include/deal.II-translator/lac/solver_richardson_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2020 by the deal.II authors
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

#ifndef dealii_solver_richardson_h
#define dealii_solver_richardson_h


#include <deal.II/base/config.h>

#include <deal.II/base/logstream.h>
#include <deal.II/base/signaling_nan.h>

#include <deal.II/lac/solver.h>
#include <deal.II/lac/solver_control.h>

DEAL_II_NAMESPACE_OPEN

 /*!@addtogroup Solvers */ 
 /*@{*/ 

/**
 * 预设条件的Richardson迭代法的实现。停止的标准是残差的规范。
 * 关于使用该类时对矩阵和向量的要求，请参见求解器基类的文档。
 * 像所有其他求解器类一样，该类有一个名为 @p
 * AdditionalData的局部结构，用于向求解器传递额外的参数，如阻尼参数或临时向量的数量。我们使用这个额外的结构，而不是直接将这些值传递给构造函数，因为这使得
 * @p SolverSelector
 * 和其他类的使用更加容易，并保证即使某个求解器的额外参数的数量或类型发生变化，这些也能继续工作。
 * 对于Richardson方法，附加数据是阻尼参数，它是 @p
 * AdditionalData
 * 结构的唯一内容。默认情况下，该结构的构造函数将其设置为1。
 *
 *  <h3>Observing the progress of linear solver iterations</h3>
 * 该类的solve()函数使用Solver基类中描述的机制来确定收敛性。这个机制也可以用来观察迭代的进度。
 *
 *
 */
template <class VectorType = Vector<double>>
class SolverRichardson : public SolverBase<VectorType>
{
public:
  /**
   * 标准化的数据结构，用于向求解器输送额外的数据。
   *
   */
  struct AdditionalData
  {
    /**
     * 构造函数。默认情况下，将阻尼参数设置为1。
     *
     */
    explicit AdditionalData(const double omega                       = 1,
                            const bool   use_preconditioned_residual = false);

    /**
     * 松弛参数。
     *
     */
    double omega;

    /**
     * 停止准则的参数。
     *
     */
    bool use_preconditioned_residual;
  };

  /**
   * 构造函数。
   *
   */
  SolverRichardson(SolverControl &           cn,
                   VectorMemory<VectorType> &mem,
                   const AdditionalData &    data = AdditionalData());

  /**
   * 构造函数。使用一个GrowingVectorMemory类型的对象作为默认分配内存。
   *
   */
  SolverRichardson(SolverControl &       cn,
                   const AdditionalData &data = AdditionalData());

  /**
   * 虚拟解构器。
   *
   */
  virtual ~SolverRichardson() override = default;

  /**
   * 求解x的线性系统 $Ax=b$ 。
   *
   */
  template <typename MatrixType, typename PreconditionerType>
  void
  solve(const MatrixType &        A,
        VectorType &              x,
        const VectorType &        b,
        const PreconditionerType &preconditioner);

  /**
   * 对 $x$ 求解 $A^Tx=b$  。
   *
   */
  template <typename MatrixType, typename PreconditionerType>
  void
  Tsolve(const MatrixType &        A,
         VectorType &              x,
         const VectorType &        b,
         const PreconditionerType &preconditioner);

  /**
   * 设置阻尼系数。默认为1.，即没有阻尼。
   *
   */
  void
  set_omega(const double om = 1.);

  /**
   * 派生类的接口。这个函数得到当前的迭代向量，残差和每一步的更新向量。它可以用于收敛历史的图形输出。
   *
   */
  virtual void
  print_vectors(const unsigned int step,
                const VectorType & x,
                const VectorType & r,
                const VectorType & d) const;

protected:
  /**
   * 实现残差的常数计算。
   * 根据给求解器的标志，该函数的默认实现使用实际残差，
   * @p r, 或预处理残差， @p d.  。
   *
   */
  virtual typename VectorType::value_type
  criterion(const VectorType &r, const VectorType &d) const;

  /**
   * 控制参数。
   *
   */
  AdditionalData additional_data;
};

 /*@}*/ 
 /*----------------- Implementation of the Richardson Method ------------------*/ 

#ifndef DOXYGEN

template <class VectorType>
inline SolverRichardson<VectorType>::AdditionalData::AdditionalData(
  const double omega,
  const bool   use_preconditioned_residual)
  : omega(omega)
  , use_preconditioned_residual(use_preconditioned_residual)
{}


template <class VectorType>
SolverRichardson<VectorType>::SolverRichardson(SolverControl &           cn,
                                               VectorMemory<VectorType> &mem,
                                               const AdditionalData &    data)
  : SolverBase<VectorType>(cn, mem)
  , additional_data(data)
{}



template <class VectorType>
SolverRichardson<VectorType>::SolverRichardson(SolverControl &       cn,
                                               const AdditionalData &data)
  : SolverBase<VectorType>(cn)
  , additional_data(data)
{}



template <class VectorType>
template <typename MatrixType, typename PreconditionerType>
void
SolverRichardson<VectorType>::solve(const MatrixType &        A,
                                    VectorType &              x,
                                    const VectorType &        b,
                                    const PreconditionerType &preconditioner)
{
  SolverControl::State conv = SolverControl::iterate;

  double last_criterion = -std::numeric_limits<double>::max();

  unsigned int iter = 0;

  // Memory allocation.
  // 'Vr' holds the residual, 'Vd' the preconditioned residual
  typename VectorMemory<VectorType>::Pointer Vr(this->memory);
  typename VectorMemory<VectorType>::Pointer Vd(this->memory);

  VectorType &r = *Vr;
  r.reinit(x);

  VectorType &d = *Vd;
  d.reinit(x);

  LogStream::Prefix prefix("Richardson");

  // Main loop
  while (conv == SolverControl::iterate)
    {
      // Do not use residual,
      // but do it in 2 steps
      A.vmult(r, x);
      r.sadd(-1., 1., b);
      preconditioner.vmult(d, r);

      // get the required norm of the (possibly preconditioned)
      // residual
      last_criterion = criterion(r, d);
      conv           = this->iteration_status(iter, last_criterion, x);
      if (conv != SolverControl::iterate)
        break;

      x.add(additional_data.omega, d);
      print_vectors(iter, x, r, d);

      ++iter;
    }

  // in case of failure: throw exception
  if (conv != SolverControl::success)
    AssertThrow(false, SolverControl::NoConvergence(iter, last_criterion));
  // otherwise exit as normal
}



template <class VectorType>
template <typename MatrixType, typename PreconditionerType>
void
SolverRichardson<VectorType>::Tsolve(const MatrixType &        A,
                                     VectorType &              x,
                                     const VectorType &        b,
                                     const PreconditionerType &preconditioner)
{
  SolverControl::State conv           = SolverControl::iterate;
  double               last_criterion = -std::numeric_limits<double>::max();

  unsigned int iter = 0;

  // Memory allocation.
  // 'Vr' holds the residual, 'Vd' the preconditioned residual
  typename VectorMemory<VectorType>::Pointer Vr(this->memory);
  typename VectorMemory<VectorType>::Pointer Vd(this->memory);

  VectorType &r = *Vr;
  r.reinit(x);

  VectorType &d = *Vd;
  d.reinit(x);

  LogStream::Prefix prefix("RichardsonT");

  // Main loop
  while (conv == SolverControl::iterate)
    {
      // Do not use Tresidual,
      // but do it in 2 steps
      A.Tvmult(r, x);
      r.sadd(-1., 1., b);
      preconditioner.Tvmult(d, r);

      last_criterion = criterion(r, d);
      conv           = this->iteration_status(iter, last_criterion, x);
      if (conv != SolverControl::iterate)
        break;

      x.add(additional_data.omega, d);
      print_vectors(iter, x, r, d);

      ++iter;
    }

  // in case of failure: throw exception
  if (conv != SolverControl::success)
    AssertThrow(false, SolverControl::NoConvergence(iter, last_criterion));

  // otherwise exit as normal
}


template <class VectorType>
void
SolverRichardson<VectorType>::print_vectors(const unsigned int,
                                            const VectorType &,
                                            const VectorType &,
                                            const VectorType &) const
{}



template <class VectorType>
inline typename VectorType::value_type
SolverRichardson<VectorType>::criterion(const VectorType &r,
                                        const VectorType &d) const
{
  if (!additional_data.use_preconditioned_residual)
    return r.l2_norm();
  else
    return d.l2_norm();
}


template <class VectorType>
inline void
SolverRichardson<VectorType>::set_omega(const double om)
{
  additional_data.omega = om;
}

#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif


