//include/deal.II-translator/lac/solver_minres_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2020 by the deal.II authors
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

#ifndef dealii_solver_minres_h
#define dealii_solver_minres_h


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

/**
 * 对称矩阵的最小残差法。
 * 为了使用这个类，对矩阵和向量的要求，见求解器基类的文档。
 * 像所有其他求解器类一样，这个类有一个名为 @p
 * AdditionalData的局部结构，用来向求解器传递额外的参数，如阻尼参数或临时向量的数量。我们使用这个额外的结构，而不是直接将这些值传递给构造函数，因为这使得
 * @p SolverSelector
 * 和其他类的使用更加容易，并保证即使某个求解器的额外参数的数量或类型发生变化，这些也能继续工作。
 * 然而，由于MinRes方法不需要额外的数据，相应的结构是空的，不提供任何功能。构造函数有一个默认参数，所以你可以在没有附加参数的情况下调用它。
 * 前提条件必须是正定和对称的。 该算法取自Astrid
 * Battermann的硕士论文，并做了一些修改。全文可在http://scholar.lib.vt.edu/theses/public/etd-12164379662151/etd-title.html。
 *
 *  <h3>Observing the progress of linear solver iterations</h3>
 * 这个类的solve()函数使用Solver基类中描述的机制来确定收敛性。这个机制也可以用来观察迭代的进度。
 *
 *
 */
template <class VectorType = Vector<double>>
class SolverMinRes : public SolverBase<VectorType>
{
public:
  /**
   * 标准化的数据结构，用于向求解器输送额外的数据。这个求解器还不需要额外的数据。
   *
   */
  struct AdditionalData
  {};

  /**
   * 构造函数。
   *
   */
  SolverMinRes(SolverControl &           cn,
               VectorMemory<VectorType> &mem,
               const AdditionalData &    data = AdditionalData());

  /**
   * 构造函数。使用一个GrowingVectorMemory类型的对象作为默认分配内存。
   *
   */
  SolverMinRes(SolverControl &       cn,
               const AdditionalData &data = AdditionalData());

  /**
   * 虚拟解构器。
   *
   */
  virtual ~SolverMinRes() override = default;

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
   * @addtogroup  Exceptions  
     * @{ 
   *
   */

  /**
   * 异常情况
   *
   */
  DeclException0(ExcPreconditionerNotDefinite);
  //@}

protected:
  /**
   * 实现计算残差的准则。
   *
   */
  virtual double
  criterion();

  /**
   * 派生类的接口。这个函数得到当前的迭代向量，残差和每一步的更新向量。它可以用于收敛历史的图形输出。
   *
   */
  virtual void
  print_vectors(const unsigned int step,
                const VectorType & x,
                const VectorType & r,
                const VectorType & d) const;

  /**
   * 在迭代循环中，残差向量的平方被存储在这个变量中。函数
   * @p criterion
   * 使用这个变量来计算收敛值，在这个类别中，它是残差向量的规范，因此是
   * @p res2 值的平方根。
   *
   */
  double res2;
};

 /*@}*/ 
 /*------------------------- Implementation ----------------------------*/ 

#ifndef DOXYGEN

template <class VectorType>
SolverMinRes<VectorType>::SolverMinRes(SolverControl &           cn,
                                       VectorMemory<VectorType> &mem,
                                       const AdditionalData &)
  : SolverBase<VectorType>(cn, mem)
  , res2(numbers::signaling_nan<double>())
{}



template <class VectorType>
SolverMinRes<VectorType>::SolverMinRes(SolverControl &cn,
                                       const AdditionalData &)
  : SolverBase<VectorType>(cn)
  , res2(numbers::signaling_nan<double>())
{}



template <class VectorType>
double
SolverMinRes<VectorType>::criterion()
{
  return res2;
}


template <class VectorType>
void
SolverMinRes<VectorType>::print_vectors(const unsigned int,
                                        const VectorType &,
                                        const VectorType &,
                                        const VectorType &) const
{}



template <class VectorType>
template <typename MatrixType, typename PreconditionerType>
void
SolverMinRes<VectorType>::solve(const MatrixType &        A,
                                VectorType &              x,
                                const VectorType &        b,
                                const PreconditionerType &preconditioner)
{
  LogStream::Prefix prefix("minres");

  // Memory allocation
  typename VectorMemory<VectorType>::Pointer Vu0(this->memory);
  typename VectorMemory<VectorType>::Pointer Vu1(this->memory);
  typename VectorMemory<VectorType>::Pointer Vu2(this->memory);

  typename VectorMemory<VectorType>::Pointer Vm0(this->memory);
  typename VectorMemory<VectorType>::Pointer Vm1(this->memory);
  typename VectorMemory<VectorType>::Pointer Vm2(this->memory);

  typename VectorMemory<VectorType>::Pointer Vv(this->memory);

  // define some aliases for simpler access
  using vecptr     = VectorType *;
  vecptr      u[3] = {Vu0.get(), Vu1.get(), Vu2.get()};
  vecptr      m[3] = {Vm0.get(), Vm1.get(), Vm2.get()};
  VectorType &v    = *Vv;

  // resize the vectors, but do not set the values since they'd be overwritten
  // soon anyway.
  u[0]->reinit(b, true);
  u[1]->reinit(b, true);
  u[2]->reinit(b, true);
  m[0]->reinit(b, true);
  m[1]->reinit(b, true);
  m[2]->reinit(b, true);
  v.reinit(b, true);

  // some values needed
  double delta[3] = {0, 0, 0};
  double f[2]     = {0, 0};
  double e[2]     = {0, 0};

  double r_l2 = 0;
  double r0   = 0;
  double tau  = 0;
  double c    = 0;
  double s    = 0;
  double d_   = 0;

  // The iteration step.
  unsigned int j = 1;


  // Start of the solution process
  A.vmult(*m[0], x);
  *u[1] = b;
  *u[1] -= *m[0];
  // Precondition is applied.
  // The preconditioner has to be
  // positive definite and symmetric

  // M v = u[1]
  preconditioner.vmult(v, *u[1]);

  delta[1] = v * (*u[1]);
  // Preconditioner positive
  Assert(delta[1] >= 0, ExcPreconditionerNotDefinite());

  r0   = std::sqrt(delta[1]);
  r_l2 = r0;


  u[0]->reinit(b);
  delta[0] = 1.;
  m[0]->reinit(b);
  m[1]->reinit(b);
  m[2]->reinit(b);

  SolverControl::State conv = this->iteration_status(0, r_l2, x);
  while (conv == SolverControl::iterate)
    {
      if (delta[1] != 0)
        v *= 1. / std::sqrt(delta[1]);
      else
        v.reinit(b);

      A.vmult(*u[2], v);
      u[2]->add(-std::sqrt(delta[1] / delta[0]), *u[0]);

      const double gamma = *u[2] * v;
      u[2]->add(-gamma / std::sqrt(delta[1]), *u[1]);
      *m[0] = v;

      // precondition: solve M v = u[2]
      // Preconditioner has to be positive
      // definite and symmetric.
      preconditioner.vmult(v, *u[2]);

      delta[2] = v * (*u[2]);

      Assert(delta[2] >= 0, ExcPreconditionerNotDefinite());

      if (j == 1)
        {
          d_   = gamma;
          e[1] = std::sqrt(delta[2]);
        }
      if (j > 1)
        {
          d_   = s * e[0] - c * gamma;
          e[0] = c * e[0] + s * gamma;
          f[1] = s * std::sqrt(delta[2]);
          e[1] = -c * std::sqrt(delta[2]);
        }

      const double d = std::sqrt(d_ * d_ + delta[2]);

      if (j > 1)
        tau *= s / c;
      c = d_ / d;
      tau *= c;

      s = std::sqrt(delta[2]) / d;

      if (j == 1)
        tau = r0 * c;

      m[0]->add(-e[0], *m[1]);
      if (j > 1)
        m[0]->add(-f[0], *m[2]);
      *m[0] *= 1. / d;
      x.add(tau, *m[0]);
      r_l2 *= std::fabs(s);

      conv = this->iteration_status(j, r_l2, x);

      // next iteration step
      ++j;
      // All vectors have to be shifted
      // one iteration step.
      // This should be changed one time.
      swap(*m[2], *m[1]);
      swap(*m[1], *m[0]);

      // likewise, but reverse direction:
      //   u[0] = u[1];
      //   u[1] = u[2];
      swap(*u[0], *u[1]);
      swap(*u[1], *u[2]);

      // these are scalars, so need
      // to bother
      f[0]     = f[1];
      e[0]     = e[1];
      delta[0] = delta[1];
      delta[1] = delta[2];
    }

  // in case of failure: throw exception
  AssertThrow(conv == SolverControl::success,
              SolverControl::NoConvergence(j, r_l2));

  // otherwise exit as normal
}

#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif


