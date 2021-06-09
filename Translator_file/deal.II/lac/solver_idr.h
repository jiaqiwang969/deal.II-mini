//include/deal.II-translator/lac/solver_idr_0.txt
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

#ifndef dealii_solver_idr_h
#define dealii_solver_idr_h


#include <deal.II/base/config.h>

#include <deal.II/base/logstream.h>
#include <deal.II/base/signaling_nan.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/solver.h>
#include <deal.II/lac/solver_control.h>

#include <cmath>
#include <random>

DEAL_II_NAMESPACE_OPEN

 /*!@addtogroup Solvers */ 
 /*@{*/ 

namespace internal
{
  /**
   * 一个用于IDR(s)求解器的辅助类的命名空间。
   *
   */
  namespace SolverIDRImplementation
  {
    /**
     * 保存临时向量的类，其大小取决于求解器的参数s。
     *
     */
    template <typename VectorType>
    class TmpVectors
    {
    public:
      /**
       * 构造函数。准备一个长度为 @p s_param. 的 @p VectorType
       * 阵列。
       *
       */
      TmpVectors(const unsigned int s_param, VectorMemory<VectorType> &vmem);

      /**
       * 销毁器。删除所有分配的向量。
       *
       */
      ~TmpVectors() = default;

      /**
       * 获取向量编号  @p i.
       * 如果这个向量之前未被使用，会发生错误。
       *
       */
      VectorType &operator[](const unsigned int i) const;

      /**
       * 获取向量编号  @p i.  必要时分配它。
       * 如果必须分配一个向量，则使用 @p temp
       * 将其重新设定为适当的尺寸。
       *
       */
      VectorType &
      operator()(const unsigned int i, const VectorType &temp);

    private:
      /**
       * 获取向量的池子。
       *
       */
      VectorMemory<VectorType> &mem;

      /**
       * 用于存储向量的字段。
       *
       */
      std::vector<typename VectorMemory<VectorType>::Pointer> data;
    };
  } // namespace SolverIDRImplementation
} // namespace internal

/**
 * 该类实现了IDR(s)方法，用于求解非对称、不确定的线性系统，在<a
 * href="https://epubs.siam.org/doi/abs/10.1137/070685804"> IDR(s): A Family
 * of Simple and Fast Algorithms for Solving Large Nonsymmetric Systems of
 * Linear Equations by Martin B. van Gijzen and Peter Sonneveld
 * </a>中开发。这里的实现是<a
 * href="https://dl.acm.org/citation.cfm?id=2049667"> Algorithm 913: An
 * Elegant IDR(s) Variant that Efficiently Exploits Biorthogonality Properties
 * by Martin B. van Gijzen and Peter
 * Sonneveld</a>中的预处理版本。本地结构 @p AdditionalData
 * 的参数s的值可以是大于等于1的任何整数。对于
 * <code>s=1</code> ，该方法的收敛性与BiCGStab相似。
 *
 *
 * @note  IDR(s)的每一次迭代都需要 <code>s+1</code>
 * 的预处理步骤和矩阵-向量积。在这个实现中，残差被更新，并在外迭代内部的每个步骤后检查收敛性。如果用户启用历史数据，这些步骤的残差将被存储，因此每次迭代将有多个值。
 *
 *
 */
template <class VectorType = Vector<double>>
class SolverIDR : public SolverBase<VectorType>
{
public:
  /**
   * 用于存储求解器所需的额外数据的结构。
   *
   */
  struct AdditionalData
  {
    /**
     * 构造函数。默认情况下，使用一个IDR(2)方法。
     *
     */
    explicit AdditionalData(const unsigned int s = 2)
      : s(s)
    {}

    const unsigned int s;
  };

  /**
   * 构造函数。
   *
   */
  SolverIDR(SolverControl &           cn,
            VectorMemory<VectorType> &mem,
            const AdditionalData &    data = AdditionalData());

  /**
   * 构造函数。使用一个GrowingVectorMemory类型的对象作为默认分配内存。
   *
   */
  explicit SolverIDR(SolverControl &       cn,
                     const AdditionalData &data = AdditionalData());

  /**
   * 虚拟解构器。
   *
   */
  virtual ~SolverIDR() override = default;

  /**
   * 求解x的线性系统 <code>Ax=b</code> 。
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
   * 派生类的接口。这个函数在每一步中获得当前的迭代向量，残差和更新向量。它可以用于收敛历史的图形输出。
   *
   */
  virtual void
  print_vectors(const unsigned int step,
                const VectorType & x,
                const VectorType & r,
                const VectorType & d) const;

private:
  /**
   * 附加求解器参数。
   *
   */
  AdditionalData additional_data;
};

 /*@}*/ 
 /*------------------------- Implementation ----------------------------*/ 

#ifndef DOXYGEN


namespace internal
{
  namespace SolverIDRImplementation
  {
    template <class VectorType>
    inline TmpVectors<VectorType>::TmpVectors(const unsigned int        s_param,
                                              VectorMemory<VectorType> &vmem)
      : mem(vmem)
      , data(s_param)
    {}



    template <class VectorType>
    inline VectorType &TmpVectors<VectorType>::
                       operator[](const unsigned int i) const
    {
      AssertIndexRange(i, data.size());

      Assert(data[i] != nullptr, ExcNotInitialized());
      return *data[i];
    }



    template <class VectorType>
    inline VectorType &
    TmpVectors<VectorType>::operator()(const unsigned int i,
                                       const VectorType & temp)
    {
      AssertIndexRange(i, data.size());
      if (data[i] == nullptr)
        {
          data[i] = std::move(typename VectorMemory<VectorType>::Pointer(mem));
          data[i]->reinit(temp);
        }
      return *data[i];
    }
  } // namespace SolverIDRImplementation
} // namespace internal



template <class VectorType>
SolverIDR<VectorType>::SolverIDR(SolverControl &           cn,
                                 VectorMemory<VectorType> &mem,
                                 const AdditionalData &    data)
  : SolverBase<VectorType>(cn, mem)
  , additional_data(data)
{}



template <class VectorType>
SolverIDR<VectorType>::SolverIDR(SolverControl &cn, const AdditionalData &data)
  : SolverBase<VectorType>(cn)
  , additional_data(data)
{}



template <class VectorType>
void
SolverIDR<VectorType>::print_vectors(const unsigned int,
                                     const VectorType &,
                                     const VectorType &,
                                     const VectorType &) const
{}



template <class VectorType>
template <typename MatrixType, typename PreconditionerType>
void
SolverIDR<VectorType>::solve(const MatrixType &        A,
                             VectorType &              x,
                             const VectorType &        b,
                             const PreconditionerType &preconditioner)
{
  LogStream::Prefix prefix("IDR(s)");

  SolverControl::State iteration_state = SolverControl::iterate;
  unsigned int         step            = 0;

  const unsigned int s = additional_data.s;

  // Define temporary vectors which do not do not depend on s
  typename VectorMemory<VectorType>::Pointer r_pointer(this->memory);
  typename VectorMemory<VectorType>::Pointer v_pointer(this->memory);
  typename VectorMemory<VectorType>::Pointer vhat_pointer(this->memory);
  typename VectorMemory<VectorType>::Pointer uhat_pointer(this->memory);
  typename VectorMemory<VectorType>::Pointer ghat_pointer(this->memory);

  VectorType &r    = *r_pointer;
  VectorType &v    = *v_pointer;
  VectorType &vhat = *vhat_pointer;
  VectorType &uhat = *uhat_pointer;
  VectorType &ghat = *ghat_pointer;

  r.reinit(x, true);
  v.reinit(x, true);
  vhat.reinit(x, true);
  uhat.reinit(x, true);
  ghat.reinit(x, true);

  // Initial residual
  A.vmult(r, x);
  r.sadd(-1.0, 1.0, b);

  // Check for convergent initial guess
  double res      = r.l2_norm();
  iteration_state = this->iteration_status(step, res, x);
  if (iteration_state == SolverControl::success)
    return;

  // Initialize sets of vectors/matrices whose size dependent on s
  internal::SolverIDRImplementation::TmpVectors<VectorType> G(s, this->memory);
  internal::SolverIDRImplementation::TmpVectors<VectorType> U(s, this->memory);
  internal::SolverIDRImplementation::TmpVectors<VectorType> Q(s, this->memory);
  FullMatrix<double>                                        M(s, s);

  // Random number generator for vector entries of
  // Q (normal distribution, mean=0 sigma=1)
  std::mt19937               rng;
  std::normal_distribution<> normal_distribution(0.0, 1.0);
  for (unsigned int i = 0; i < s; ++i)
    {
      VectorType &tmp_g = G(i, x);
      VectorType &tmp_u = U(i, x);
      tmp_g             = 0;
      tmp_u             = 0;

      // Compute random set of s orthonormalized vectors Q
      // Note: the first vector is chosen to be the initial
      // residual to match BiCGStab (as is done in comparisons
      // with BiCGStab in the papers listed in the documentation
      // of this function)
      VectorType &tmp_q = Q(i, x);
      if (i != 0)
        {
          for (auto indx : tmp_q.locally_owned_elements())
            tmp_q(indx) = normal_distribution(rng);
          tmp_q.compress(VectorOperation::insert);
        }
      else
        tmp_q = r;

      for (unsigned int j = 0; j < i; ++j)
        {
          v = Q[j];
          v *= (v * tmp_q) / (tmp_q * tmp_q);
          tmp_q.add(-1.0, v);
        }

      if (i != 0)
        tmp_q *= 1.0 / tmp_q.l2_norm();

      M(i, i) = 1.;
    }

  double omega = 1.;

  bool early_exit = false;

  // Outer iteration
  while (iteration_state == SolverControl::iterate)
    {
      ++step;

      // Compute phi
      Vector<double> phi(s);
      for (unsigned int i = 0; i < s; ++i)
        phi(i) = Q[i] * r;

      // Inner iteration over s
      for (unsigned int k = 0; k < s; ++k)
        {
          // Solve M(k:s)*gamma = phi(k:s)
          Vector<double> gamma(s - k);
          {
            Vector<double>            phik(s - k);
            FullMatrix<double>        Mk(s - k, s - k);
            std::vector<unsigned int> indices;
            unsigned int              j = 0;
            for (unsigned int i = k; i < s; ++i, ++j)
              {
                indices.push_back(i);
                phik(j) = phi(i);
              }
            Mk.extract_submatrix_from(M, indices, indices);

            FullMatrix<double> Mk_inv(s - k, s - k);
            Mk_inv.invert(Mk);
            Mk_inv.vmult(gamma, phik);
          }

          {
            v = r;

            unsigned int j = 0;
            for (unsigned int i = k; i < s; ++i, ++j)
              v.add(-1.0 * gamma(j), G[i]);
            preconditioner.vmult(vhat, v);

            uhat = vhat;
            uhat *= omega;
            j = 0;
            for (unsigned int i = k; i < s; ++i, ++j)
              uhat.add(gamma(j), U[i]);
            A.vmult(ghat, uhat);
          }

          // Update G and U
          // Orthogonalize ghat to Q0,..,Q_{k-1}
          // and update uhat
          for (unsigned int i = 0; i < k; ++i)
            {
              double alpha = (Q[i] * ghat) / M(i, i);
              ghat.add(-alpha, G[i]);
              uhat.add(-alpha, U[i]);
            }
          G[k] = ghat;
          U[k] = uhat;

          // Update kth column of M
          for (unsigned int i = k; i < s; ++i)
            M(i, k) = Q[i] * G[k];

          // Orthogonalize r to Q0,...,Qk,
          // update x
          {
            double beta = phi(k) / M(k, k);
            r.add(-1.0 * beta, G[k]);
            x.add(beta, U[k]);

            print_vectors(step, x, r, U[k]);

            // Check for early convergence. If so, store
            // information in early_exit so that outer iteration
            // is broken before recomputing the residual
            res             = r.l2_norm();
            iteration_state = this->iteration_status(step, res, x);
            if (iteration_state != SolverControl::iterate)
              {
                early_exit = true;
                break;
              }

            // Update phi
            if (k + 1 < s)
              {
                for (unsigned int i = 0; i < k + 1; ++i)
                  phi(i) = 0.0;
                for (unsigned int i = k + 1; i < s; ++i)
                  phi(i) -= beta * M(i, k);
              }
          }
        }
      if (early_exit == true)
        break;

      // Update r and x
      preconditioner.vmult(vhat, r);
      A.vmult(v, vhat);

      omega = (v * r) / (v * v);

      r.add(-1.0 * omega, v);
      x.add(omega, vhat);

      print_vectors(step, x, r, vhat);

      // Check for convergence
      res             = r.l2_norm();
      iteration_state = this->iteration_status(step, res, x);
      if (iteration_state != SolverControl::iterate)
        break;
    }

  if (iteration_state != SolverControl::success)
    AssertThrow(false, SolverControl::NoConvergence(step, res));
}


#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif


