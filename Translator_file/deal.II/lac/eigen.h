//include/deal.II-translator/lac/eigen_0.txt
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

#ifndef dealii_eigen_h
#define dealii_eigen_h


#include <deal.II/base/config.h>

#include <deal.II/lac/linear_operator.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/solver_minres.h>
#include <deal.II/lac/vector_memory.h>

#include <cmath>

DEAL_II_NAMESPACE_OPEN


 /*!@addtogroup Solvers */ 
 /*@{*/ 

/**
 * 用于特征值计算的功率法（von Mises）。
 * 这种方法通过将一个矩阵的增幂应用于一个向量来确定该矩阵的最大特征值。如果有一个绝对值占优势的特征值
 * $l$ ，迭代的向量将与它的特征空间和 $Ax = lx$ 对齐。
 * 一个移动参数允许移动频谱，所以也可以计算最小的特征值。
 * 这种方法的收敛性已知是很慢的。
 *
 *
 */
template <typename VectorType = Vector<double>>
class EigenPower : private SolverBase<VectorType>
{
public:
  /**
   * 声明容器大小的类型。
   *
   */
  using size_type = types::global_dof_index;

  /**
   * 标准化的数据结构，用于向求解器输送额外的数据。
   *
   */
  struct AdditionalData
  {
    /**
     * 移位参数。这个参数允许转移频谱以计算不同的特征值。
     *
     */
    double shift;
    /**
     * 构造函数。设置移位参数。
     *
     */
    AdditionalData(const double shift = 0.)
      : shift(shift)
    {}
  };

  /**
   * 构造函数。
   *
   */
  EigenPower(SolverControl &           cn,
             VectorMemory<VectorType> &mem,
             const AdditionalData &    data = AdditionalData());


  /**
   * 功率法。  @p x
   * 是功率法的（不一定是归一化的，但不为零）起始向量。迭代后，
   * @p value 是近似的特征值， @p x
   * 是相应的特征向量，相对于l2-norm进行归一化。
   *
   */
  template <typename MatrixType>
  void
  solve(double &value, const MatrixType &A, VectorType &x);

protected:
  /**
   * 移位参数。
   *
   */
  AdditionalData additional_data;
};

/**
 * 用于特征值计算的逆向迭代（Wieland）。
 * 该类实现了Wieland反迭代的自适应版本。
 * 对于停止标准有两种选择：默认情况下，计算残差 $A x
 *
 * - l x$
 * 的准则。由于这对于非对称矩阵的非琐碎约旦块来说可能不会收敛为零，因此可以用检查连续的特征值之差来代替。使用
 * AdditionalData::use_residual 来切换这个选项。
 * 通常，进入这个方法的初始猜测在每一步之后都会被更新，用新的特征值的近似值来代替它。使用参数
 * AdditionalData::relaxation
 * 在0和1之间，这个更新可以被阻尼。如果放松参数为0，则不进行更新。这种阻尼允许较慢地适应转移值，以确保方法收敛到最接近初始猜测的特征值。这可以通过参数
 * AdditionalData::start_adaption,
 * 来帮助实现，该参数表示应该调整移位值的第一个迭代步骤。
 *
 *
 */
template <typename VectorType = Vector<double>>
class EigenInverse : private SolverBase<VectorType>
{
public:
  /**
   * 声明容器大小的类型。
   *
   */
  using size_type = types::global_dof_index;

  /**
   * 标准化的数据结构，用于向求解器输送额外的数据。
   *
   */
  struct AdditionalData
  {
    /**
     * 更新移位值的阻尼。
     *
     */
    double relaxation;

    /**
     * 自适应移位参数的起始步骤。
     *
     */
    unsigned int start_adaption;
    /**
     * 停止标准的标志。
     *
     */
    bool use_residual;
    /**
     * 构造函数。
     *
     */
    AdditionalData(double       relaxation     = 1.,
                   unsigned int start_adaption = 6,
                   bool         use_residual   = true)
      : relaxation(relaxation)
      , start_adaption(start_adaption)
      , use_residual(use_residual)
    {}
  };

  /**
   * 构造函数。
   *
   */
  EigenInverse(SolverControl &           cn,
               VectorMemory<VectorType> &mem,
               const AdditionalData &    data = AdditionalData());

  /**
   * 逆向方法。  @p value 是特征值的起始猜测， @p x
   * 是功率法的（不一定是归一化的，但不为零）起始向量。迭代后，
   * @p value 是近似的特征值， @p x
   * 是相应的特征向量，相对于l2-norm进行归一化。
   *
   */
  template <typename MatrixType>
  void
  solve(double &value, const MatrixType &A, VectorType &x);

protected:
  /**
   * 执行的标志。
   *
   */
  AdditionalData additional_data;
};

 /*@}*/ 
//---------------------------------------------------------------------------


template <class VectorType>
EigenPower<VectorType>::EigenPower(SolverControl &           cn,
                                   VectorMemory<VectorType> &mem,
                                   const AdditionalData &    data)
  : SolverBase<VectorType>(cn, mem)
  , additional_data(data)
{}



template <class VectorType>
template <typename MatrixType>
void
EigenPower<VectorType>::solve(double &value, const MatrixType &A, VectorType &x)
{
  SolverControl::State conv = SolverControl::iterate;

  LogStream::Prefix prefix("Power method");

  typename VectorMemory<VectorType>::Pointer Vy(this->memory);
  VectorType &                               y = *Vy;
  y.reinit(x);
  typename VectorMemory<VectorType>::Pointer Vr(this->memory);
  VectorType &                               r = *Vr;
  r.reinit(x);

  double length     = x.l2_norm();
  double old_length = 0.;
  x *= 1. / length;

  A.vmult(y, x);

  // Main loop
  int iter = 0;
  for (; conv == SolverControl::iterate; iter++)
    {
      y.add(additional_data.shift, x);

      // Compute absolute value of eigenvalue
      old_length = length;
      length     = y.l2_norm();

      // do a little trick to compute the sign
      // with not too much effect of round-off errors.
      double    entry  = 0.;
      size_type i      = 0;
      double    thresh = length / x.size();
      do
        {
          Assert(i < x.size(), ExcInternalError());
          entry = y(i++);
        }
      while (std::fabs(entry) < thresh);

      --i;

      // Compute unshifted eigenvalue
      value = (entry * x(i) < 0.) ? -length : length;
      value -= additional_data.shift;

      // Update normalized eigenvector
      x.equ(1 / length, y);

      // Compute residual
      A.vmult(y, x);

      // Check the change of the eigenvalue
      // Brrr, this is not really a good criterion
      conv = this->iteration_status(iter,
                                    std::fabs(1. / length - 1. / old_length),
                                    x);
    }

  // in case of failure: throw exception
  AssertThrow(conv == SolverControl::success,
              SolverControl::NoConvergence(
                iter, std::fabs(1. / length - 1. / old_length)));

  // otherwise exit as normal
}

//---------------------------------------------------------------------------

template <class VectorType>
EigenInverse<VectorType>::EigenInverse(SolverControl &           cn,
                                       VectorMemory<VectorType> &mem,
                                       const AdditionalData &    data)
  : SolverBase<VectorType>(cn, mem)
  , additional_data(data)
{}



template <class VectorType>
template <typename MatrixType>
void
EigenInverse<VectorType>::solve(double &          value,
                                const MatrixType &A,
                                VectorType &      x)
{
  LogStream::Prefix prefix("Wielandt");

  SolverControl::State conv = SolverControl::iterate;

  // Prepare matrix for solver
  auto   A_op          = linear_operator(A);
  double current_shift = -value;
  auto   A_s           = A_op + current_shift * identity_operator(A_op);

  // Define solver
  ReductionControl        inner_control(5000, 1.e-16, 1.e-5, false, false);
  PreconditionIdentity    prec;
  SolverGMRES<VectorType> solver(inner_control, this->memory);

  // Next step for recomputing the shift
  unsigned int goal = additional_data.start_adaption;

  // Auxiliary vector
  typename VectorMemory<VectorType>::Pointer Vy(this->memory);
  VectorType &                               y = *Vy;
  y.reinit(x);
  typename VectorMemory<VectorType>::Pointer Vr(this->memory);
  VectorType &                               r = *Vr;
  r.reinit(x);

  double length    = x.l2_norm();
  double old_value = value;

  x *= 1. / length;

  // Main loop
  double    res  = -std::numeric_limits<double>::max();
  size_type iter = 0;
  for (; conv == SolverControl::iterate; iter++)
    {
      solver.solve(A_s, y, x, prec);

      // Compute absolute value of eigenvalue
      length = y.l2_norm();

      // do a little trick to compute the sign
      // with not too much effect of round-off errors.
      double    entry  = 0.;
      size_type i      = 0;
      double    thresh = length / x.size();
      do
        {
          Assert(i < x.size(), ExcInternalError());
          entry = y(i++);
        }
      while (std::fabs(entry) < thresh);

      --i;

      // Compute unshifted eigenvalue
      value = (entry * x(i) < 0. ? -1. : 1.) / length - current_shift;

      if (iter == goal)
        {
          const auto & relaxation = additional_data.relaxation;
          const double new_shift =
            relaxation * (-value) + (1. - relaxation) * current_shift;

          A_s           = A_op + new_shift * identity_operator(A_op);
          current_shift = new_shift;

          ++goal;
        }

      // Update normalized eigenvector
      x.equ(1. / length, y);
      // Compute residual
      if (additional_data.use_residual)
        {
          y.equ(value, x);
          A.vmult(r, x);
          r.sadd(-1., value, x);
          res = r.l2_norm();
          // Check the residual
          conv = this->iteration_status(iter, res, x);
        }
      else
        {
          res  = std::fabs(1. / value - 1. / old_value);
          conv = this->iteration_status(iter, res, x);
        }
      old_value = value;
    }

  // in case of failure: throw
  // exception
  AssertThrow(conv == SolverControl::success,
              SolverControl::NoConvergence(iter, res));
  // otherwise exit as normal
}

DEAL_II_NAMESPACE_CLOSE

#endif


