//include/deal.II-translator/lac/solver_qmrs_0.txt
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

#ifndef dealii_solver_qmrs_h
#define dealii_solver_qmrs_h

#include <deal.II/base/config.h>

#include <deal.II/base/logstream.h>
#include <deal.II/base/subscriptor.h>

#include <deal.II/lac/solver.h>
#include <deal.II/lac/solver_control.h>

#include <cmath>

DEAL_II_NAMESPACE_OPEN

 /*!@addtogroup Solvers */ 
 /*@{*/ 

/**
 * <h3>Quasi-minimal method for symmetric matrices (SQMR)</h3>
 * SQMR（对称准最小残差）方法应该是用来解决具有对称的、不一定是确定的预处理的对称不定线性系统。它是原始准最小残差法（QMR）的变体，产生相同的迭代解。这个版本的SQMR是由Freund/Nachtigal各自给出的对称QMR-from-BiCG算法改编的。A
 * new Krylov-subspace method for symmetric indefinite linear systems, NASA
 * STI/Recon Technical Report N, 95 (1994) and Freund/Nachtigal: Software for
 * simplified Lanczos and QMR algorithms, Appl. Num. Math. 19 (1995), pp.
 * 319-341，并提供了右和左（但不是分裂）预处理。
 *
 *  <h3>Trade off of stability to simplicity</h3>
 * 请注意，给定算法所基于的QMR实现是由经典的BiCG导出的。可以证明（Freund/Szeto:
 * A transpos-free quasi-minimal residual squared algorithm for non-Hermitian
 * linear systems, Advances in Computer Methods for Partial Differential
 * Equations VII (IMACS, New Brunswick, NJ, 1992)
 * pp.258-264)，QMR迭代结果可以通过一个额外的向量和一些标量更新从BiCG迭代中生成。因此，BiCG可能出现的故障（准确地说是除以0）显然会转移到这个简单的无先验算法上。
 * 与经典的QMR或BiCGStab相比，该算法的成本很低，每次迭代只使用一个与系统矩阵的矩阵向量乘积和一个预处理程序的应用。
 * 用于衡量收敛性的残差只是通过一个上界来近似计算。如果该值低于AdditionalData结构中规定的阈值，那么当前QMR迭代的精确残差将通过与系统矩阵的另一次乘法来计算。根据经验（根据Freund和Nachtigal），这种技术对于阈值是求解容忍度的10倍是有用的，在这种情况下，只在完整迭代的最后一到两步使用。
 * 为了使用这个类，对矩阵和向量的要求，见求解器基类的文档。
 * 像所有其他求解器类一样，该类有一个名为 @p
 * AdditionalData的局部结构，用于向求解器传递额外的参数，如阻尼参数或临时向量的数量。我们使用这个额外的结构，而不是直接将这些值传递给构造函数，因为这使得
 * @p SolverSelector
 * 和其他类的使用更加容易，并保证即使某个求解器的额外参数的数量或类型发生变化，这些也能继续工作。
 *
 *  <h3>Observing the progress of linear solver iterations</h3>
 * 这个类的solve()函数使用Solver基类中描述的机制来确定收敛性。这个机制也可以用来观察迭代的进度。
 *
 *
 */
template <typename VectorType = Vector<double>>
class SolverQMRS : public SolverBase<VectorType>
{
public:
  /**
   * 标准化的数据结构，用于向求解器输送额外的数据。
   * 用户能够在右预处理和左预处理之间进行切换，也就是使用相应的参数分别求解<i>P<sup>-1</sup>A</i>和<i>AP<sup>-1</sup></i>的系统。请注意，左预处理意味着采用预处理的（BiCG-）残差，否则采用未预处理的残差。默认情况是采用右侧的应用。
   * @p solver_tolerance
   * 阈值用于定义所说的边界，在该边界以下的残差被精确计算。更多信息请参见类文件。默认值是1e-9，也就是默认解算精度乘以10。
   * SQMR容易受到分解（除以零）的影响，所以我们需要一个参数告诉我们哪些数字被认为是零。适当的分解标准非常不清楚，所以这里可能需要进行实验。
   * 甚至有可能在除以小数的情况下也能实现收敛。甚至在有些情况下，接受这样的除法是有利的，因为廉价的迭代成本使该算法成为所有可用的不定式迭代求解器中最快的。尽管如此，默认的分解阈值是1e-16。
   *
   */
  struct AdditionalData
  {
    /**
     * 构造器。        默认是右预处理， @p solver_tolerance
     * 选择为1e-9， @p breakdown_threshold 设置为1e-16。
     *
     */
    explicit AdditionalData(const bool   left_preconditioning = false,
                            const double solver_tolerance     = 1.e-9,
                            const bool   breakdown_testing    = true,
                            const double breakdown_threshold  = 1.e-16)
      : left_preconditioning(left_preconditioning)
      , solver_tolerance(solver_tolerance)
      , breakdown_testing(breakdown_testing)
      , breakdown_threshold(breakdown_threshold)
    {}

    /**
     * 标志着使用左置条件的版本。
     *
     */
    bool left_preconditioning;

    /**
     * 准确计算当前残差的阈值。
     *
     */
    double solver_tolerance;

    /**
     * 击穿测试的标志。
     *
     */
    bool breakdown_testing;

    /**
     * 崩溃阈值。测量到这个界限的标度被用于划分。
     *
     */
    double breakdown_threshold;
  };

  /**
   * 构造器。
   *
   */
  SolverQMRS(SolverControl &           cn,
             VectorMemory<VectorType> &mem,
             const AdditionalData &    data = AdditionalData());

  /**
   * 构造函数。使用一个GrowingVectorMemory类型的对象作为默认分配内存。
   *
   */
  SolverQMRS(SolverControl &cn, const AdditionalData &data = AdditionalData());

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
   * 派生类的接口。这个函数在每一步中获得当前的迭代向量、残差和更新向量。它可以用于收敛历史的图形输出。
   *
   */
  virtual void
  print_vectors(const unsigned int step,
                const VectorType & x,
                const VectorType & r,
                const VectorType & d) const;

protected:
  /**
   * 附加参数。
   *
   */
  AdditionalData additional_data;

private:
  /**
   * 由iterate()函数返回的一个结构，代表它发现在迭代过程中发生的事情。
   *
   */
  struct IterationResult
  {
    SolverControl::State state;
    double               last_residual;

    IterationResult(const SolverControl::State state,
                    const double               last_residual);
  };

  /**
   * 迭代循环本身。该函数返回一个结构，表示在这个函数中发生了什么。
   *
   */
  template <typename MatrixType, typename PreconditionerType>
  IterationResult
  iterate(const MatrixType &        A,
          VectorType &              x,
          const VectorType &        b,
          const PreconditionerType &preconditioner,
          VectorType &              r,
          VectorType &              u,
          VectorType &              q,
          VectorType &              t,
          VectorType &              d);

  /**
   * 当前迭代的编号（在重启过程中累积）。
   *
   */
  unsigned int step;
};

 /*@}*/ 
 /*------------------------- Implementation ----------------------------*/ 

#ifndef DOXYGEN


template <class VectorType>
SolverQMRS<VectorType>::IterationResult::IterationResult(
  const SolverControl::State state,
  const double               last_residual)
  : state(state)
  , last_residual(last_residual)
{}



template <class VectorType>
SolverQMRS<VectorType>::SolverQMRS(SolverControl &           cn,
                                   VectorMemory<VectorType> &mem,
                                   const AdditionalData &    data)
  : SolverBase<VectorType>(cn, mem)
  , additional_data(data)
  , step(0)
{}

template <class VectorType>
SolverQMRS<VectorType>::SolverQMRS(SolverControl &       cn,
                                   const AdditionalData &data)
  : SolverBase<VectorType>(cn)
  , additional_data(data)
  , step(0)
{}

template <class VectorType>
void
SolverQMRS<VectorType>::print_vectors(const unsigned int,
                                      const VectorType &,
                                      const VectorType &,
                                      const VectorType &) const
{}

template <class VectorType>
template <typename MatrixType, typename PreconditionerType>
void
SolverQMRS<VectorType>::solve(const MatrixType &        A,
                              VectorType &              x,
                              const VectorType &        b,
                              const PreconditionerType &preconditioner)
{
  LogStream::Prefix prefix("SQMR");


  // temporary vectors, allocated trough the @p VectorMemory object at the
  // start of the actual solution process and deallocated at the end.
  typename VectorMemory<VectorType>::Pointer Vr(this->memory);
  typename VectorMemory<VectorType>::Pointer Vu(this->memory);
  typename VectorMemory<VectorType>::Pointer Vq(this->memory);
  typename VectorMemory<VectorType>::Pointer Vt(this->memory);
  typename VectorMemory<VectorType>::Pointer Vd(this->memory);


  // resize the vectors, but do not set
  // the values since they'd be overwritten
  // soon anyway.
  Vr->reinit(x, true);
  Vu->reinit(x, true);
  Vq->reinit(x, true);
  Vt->reinit(x, true);
  Vd->reinit(x, true);

  step = 0;

  IterationResult state(SolverControl::failure, 0);

  do
    {
      if (step > 0)
        deallog << "Restart step " << step << std::endl;
      state = iterate(A, x, b, preconditioner, *Vr, *Vu, *Vq, *Vt, *Vd);
    }
  while (state.state == SolverControl::iterate);


  // in case of failure: throw exception
  AssertThrow(state.state == SolverControl::success,
              SolverControl::NoConvergence(step, state.last_residual));
  // otherwise exit as normal
}

template <class VectorType>
template <typename MatrixType, typename PreconditionerType>
typename SolverQMRS<VectorType>::IterationResult
SolverQMRS<VectorType>::iterate(const MatrixType &        A,
                                VectorType &              x,
                                const VectorType &        b,
                                const PreconditionerType &preconditioner,
                                VectorType &              r,
                                VectorType &              u,
                                VectorType &              q,
                                VectorType &              t,
                                VectorType &              d)
{
  SolverControl::State state = SolverControl::iterate;

  int it = 0;

  double tau, rho, theta = 0;
  double res;

  // Compute the start residual
  A.vmult(r, x);
  r.sadd(-1., 1., b);

  // Doing the initial preconditioning
  if (additional_data.left_preconditioning)
    {
      // Left preconditioning
      preconditioner.vmult(t, r);
      q = t;
    }
  else
    {
      // Right preconditioning
      t = r;
      preconditioner.vmult(q, t);
    }

  tau = t.norm_sqr();
  res = std::sqrt(tau);

  if (this->iteration_status(step, res, x) == SolverControl::success)
    return IterationResult(SolverControl::success, res);

  rho = q * r;

  while (state == SolverControl::iterate)
    {
      step++;
      it++;
      //--------------------------------------------------------------
      // Step 1: apply the system matrix and compute one inner product
      //--------------------------------------------------------------
      A.vmult(t, q);
      const double sigma = q * t;

      // Check the breakdown criterion
      if (additional_data.breakdown_testing == true &&
          std::fabs(sigma) < additional_data.breakdown_threshold)
        return IterationResult(SolverControl::iterate, res);
      // Update the residual
      const double alpha = rho / sigma;
      r.add(-alpha, t);

      //--------------------------------------------------------------
      // Step 2: update the solution vector
      //--------------------------------------------------------------
      const double theta_old = theta;

      // Apply the preconditioner
      if (additional_data.left_preconditioning)
        {
          // Left Preconditioning
          preconditioner.vmult(t, r);
        }
      else
        {
          // Right Preconditioning
          t = r;
        }

      // Double updates
      theta            = t * t / tau;
      const double psi = 1. / (1. + theta);
      tau *= theta * psi;

      // Actual update of the solution vector
      d.sadd(psi * theta_old, psi * alpha, q);
      x += d;

      print_vectors(step, x, r, d);

      // Check for convergence
      // Compute a simple and cheap upper bound of the norm of the residual
      // vector b-Ax
      res = std::sqrt((it + 1) * tau);
      // If res lies close enough, within the desired tolerance, calculate the
      // exact residual
      if (res < additional_data.solver_tolerance)
        {
          A.vmult(u, x);
          u.sadd(-1., 1., b);
          res = u.l2_norm();
        }
      state = this->iteration_status(step, res, x);
      if ((state == SolverControl::success) ||
          (state == SolverControl::failure))
        return IterationResult(state, res);

      //--------------------------------------------------------------
      // Step 3: check breakdown criterion and update the vectors
      //--------------------------------------------------------------
      if (additional_data.breakdown_testing == true &&
          std::fabs(sigma) < additional_data.breakdown_threshold)
        return IterationResult(SolverControl::iterate, res);

      const double rho_old = rho;

      // Applying the preconditioner
      if (additional_data.left_preconditioning)
        {
          // Left preconditioning
          u = t;
        }
      else
        {
          // Right preconditioning
          preconditioner.vmult(u, t);
        }

      // Double and vector updates
      rho               = u * r;
      const double beta = rho / rho_old;
      q.sadd(beta, 1., u);
    }
  return IterationResult(SolverControl::success, res);
}

#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif


