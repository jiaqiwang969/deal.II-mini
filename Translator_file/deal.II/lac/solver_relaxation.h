//include/deal.II-translator/lac/solver_relaxation_0.txt
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

#ifndef dealii_solver_relaxation_h
#define dealii_solver_relaxation_h


#include <deal.II/base/config.h>

#include <deal.II/base/logstream.h>
#include <deal.II/base/subscriptor.h>

#include <deal.II/lac/solver.h>
#include <deal.II/lac/solver_control.h>

DEAL_II_NAMESPACE_OPEN

/**
 * 实现一个基于松弛方法的迭代求解器。停止的标准是残差的准则。
 * 关于使用该类时对矩阵和向量的要求，请参见Solver基类的文档。
 * 像所有其他求解器类一样，该类有一个名为 @p
 * AdditionalData的局部结构，用于向求解器传递额外的参数，如阻尼参数或临时向量的数量。我们使用这个额外的结构，而不是直接将这些值传递给构造函数，因为这使得
 * @p SolverSelector
 * 和其他类的使用更加容易，并保证即使某个求解器的额外参数的数量或类型发生变化，这些也能继续工作。这个类的AdditionalData目前不包含任何数据。
 *
 *  <h3>Observing the progress of linear solver iterations</h3>
 * 该类的solve()函数使用Solver基类中描述的机制来确定收敛性。这个机制也可以用来观察迭代的进度。
 *
 *
 *
 * @ingroup Solvers
 *
 *
 */
template <typename VectorType = Vector<double>>
class SolverRelaxation : public SolverBase<VectorType>
{
public:
  /**
   * 标准化的数据结构，用于向求解器输送额外数据。这里没有松弛方法的数据。
   *
   */
  struct AdditionalData
  {};

  /**
   * 构造函数。
   *
   */
  SolverRelaxation(SolverControl &       cn,
                   const AdditionalData &data = AdditionalData());

  /**
   * 使用松弛方法  $x_{k+1} = R(x_k,b)$  解决系统  $Ax = b$
   * 。矩阵<i>A</i>本身只用于计算残差。
   *
   */
  template <typename MatrixType, class RelaxationType>
  void
  solve(const MatrixType &    A,
        VectorType &          x,
        const VectorType &    b,
        const RelaxationType &R);
};

//----------------------------------------------------------------------//

template <class VectorType>
SolverRelaxation<VectorType>::SolverRelaxation(SolverControl &cn,
                                               const AdditionalData &)
  : SolverBase<VectorType>(cn)
{}



template <class VectorType>
template <typename MatrixType, class RelaxationType>
void
SolverRelaxation<VectorType>::solve(const MatrixType &    A,
                                    VectorType &          x,
                                    const VectorType &    b,
                                    const RelaxationType &R)
{
  GrowingVectorMemory<VectorType> mem;
  SolverControl::State            conv = SolverControl::iterate;

  // Memory allocation
  typename VectorMemory<VectorType>::Pointer Vr(mem);
  VectorType &                               r = *Vr;
  r.reinit(x);
  typename VectorMemory<VectorType>::Pointer Vd(mem);
  VectorType &                               d = *Vd;
  d.reinit(x);

  LogStream::Prefix prefix("Relaxation");

  int iter = 0;
  // Main loop
  for (; conv == SolverControl::iterate; iter++)
    {
      // Compute residual
      A.vmult(r, x);
      r.sadd(-1., 1., b);

      // The required norm of the
      // (preconditioned)
      // residual is computed in
      // criterion() and stored
      // in res.
      conv = this->iteration_status(iter, r.l2_norm(), x);
      if (conv != SolverControl::iterate)
        break;
      R.step(x, b);
    }

  // in case of failure: throw exception
  AssertThrow(conv == SolverControl::success,
              SolverControl::NoConvergence(iter, r.l2_norm()));
  // otherwise exit as normal
}


DEAL_II_NAMESPACE_CLOSE

#endif


