//include/deal.II-translator/lac/utilities_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2020 by the deal.II authors
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

#ifndef dealii_lac_utilities_h
#define dealii_lac_utilities_h

#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/signaling_nan.h>

#include <deal.II/lac/lapack_support.h>
#include <deal.II/lac/vector_memory.h>

#include <array>
#include <complex>

DEAL_II_NAMESPACE_OPEN

namespace Utilities
{
  /**
   * 一个线性代数实用程序的集合。
   *
   */
  namespace LinearAlgebra
  {
    /**
     * 返回一个连续的吉文斯旋转矩阵的元素和输入矢量的规范。
     * 即对于给定的一对 @p x 和 @p y, ，返回 $c$  ， $s$ 和
     * $\sqrt{x^2+y^2}$ ，以便\f[ \begin{bmatrix} c  & s \\
     *
     *
     *
     *
     *
     *
     *
     *
     *
     *
     *
     * -s & c \end{bmatrix} \begin{bmatrix} x \\ y \end{bmatrix} =
     * \begin{bmatrix} \sqrt{x^2+y^2} \\ 0 \end{bmatrix} \f] 。
     * @note  该函数只对实值数实现。
     *
     */
    template <typename NumberType>
    std::array<NumberType, 3>
    givens_rotation(const NumberType &x, const NumberType &y);

    /**
     * 返回一个双曲旋转矩阵的元素。
     * 即对于给定的一对  @p x  和  @p y,  返回  $c$  ,  $s$  和
     * $r$  这样 \f[ \begin{bmatrix} c  &
     *
     * -s \\
     *
     *
     *
     *
     *
     *
     *
     *
     *
     *
     *
     * -s & c \end{bmatrix} \begin{bmatrix} x \\ y \end{bmatrix} =
     * \begin{bmatrix} r \\ 0 \end{bmatrix} \f] 只有在  $|x|>|g|$
     * 时存在实值解，否则函数将抛出一个错误。
     * @note  该函数只对实值数实现。
     *
     */
    template <typename NumberType>
    std::array<NumberType, 3>
    hyperbolic_rotation(const NumberType &x, const NumberType &y);

    /**
     * 通过一个 @p k 估计 @p H 的最大特征值的上界。
     *
     * - tep Lanczos过程从初始向量 @p v0. 开始， @p k 的典型值低于10。这个估计器计算一个k步Lanczos分解  $H V_k=V_k T_k+f_k e_k^T$  ，其中  $V_k$  包含k个Lanczos基，  $V_k^TV_k=I_k$  ，  $T_k$  是三对角Lanczos矩阵，  $f_k$  是一个残差向量  $f_k^TV_k=0$  ，  $e_k$  是  $R^k$  的第k个典范基。返回的值是  $ ||T_k||_2 + ||f_k||_2$  。    如果 @p eigenvalues 不是 <code>nullptr</code> ， $T_k$ 的特征值将被写在那里。          @p vector_memory  用于为临时向量分配内存。    OperatorType必须提供 <code>vmult</code> 与VectorType的操作。        这个函数实现了以下的算法
     * @code{.bib}
     * @article{Zhou2006,
     * Title   = {Self-consistent-field Calculations Using Chebyshev-filtered
     *            Subspace Iteration},
     * Author  = {Zhou, Yunkai and Saad, Yousef and Tiago, Murilo L. and
     *            Chelikowsky, James R.},
     * Journal = {Journal of Computational Physics},
     * Year    = {2006},
     * Volume  = {219},
     * Pages   = {172--184},
     * }
     * @endcode
     *
     * @note  这个函数使用Lapack例程来计算  $T_k$
     * 的最大特征值。
     * @note  这个函数提供了一个替代的估计，即用
     * SolverCG<VectorType>::connect_eigenvalues_slot().
     * 从SolverCG的几个步骤中得到的估计。
     *
     */
    template <typename OperatorType, typename VectorType>
    double
    lanczos_largest_eigenvalue(const OperatorType &      H,
                               const VectorType &        v0,
                               const unsigned int        k,
                               VectorMemory<VectorType> &vector_memory,
                               std::vector<double> *     eigenvalues = nullptr);

    /**
     * 将算子 @p H 的切比雪夫多项式应用于 @p x.
     * 对于一个具有完整特征对 $H \psi_i = \lambda_i \psi_i$
     * 的非缺陷算子 $H$ ，多项式滤波器 $p$ 的作用由 $p(H)x
     * =\sum_i a_i p(\lambda_i) \psi_i$ 给出，其中 $x=: \sum_i a_i
     * \psi_i$
     * 。因此，通过适当地选择多项式滤波器，可以改变  $x$
     * 中包含的特征模式。
     * 这个函数使用切比雪夫第一类多项式。下面是一个度数为
     * $n=8$ 的多项式 $T_n(x)$ 在 $-1.2$ 处归一的例子。  <table>
     * <tr> <td align="center">
           @image html chebyshev8.png
     * </td> </tr> </table>  通过引入从 @p unwanted_spectrum 到 $[-1,1]$
     * 的线性映射 $L$ ，我们可以在 @p x.
     * 甩出相应的模式，多项式的度数 $n$ 越高，它在 $[-1,1]$
     * 之外的增长越快。为了避免数字溢出，我们在  @p tau.
     * 处将多项式滤波归一，因此，滤波后的算子是  $p(H) =
     * T_n(L(H))/T_n(L(\tau))$  。
     * 切比雪夫滤波器的作用只需要对 <code>vmult()</code> 的 @p
     * H
     * 进行评估，并基于切比雪夫多项式的递归方程，其度数为
     * $n$  :  $T_{n}(x) = 2x T_{n-1}(x)
     *
     * - T_{n-2}(x)$  与  $T_0(x)=1$  和  $T_1(x)=x$  。          @p
     * vector_memory 被用来为临时对象分配内存。
     * 这个函数实现了 $\sigma_1$ 中的算法（对 $\sigma_1$
     * 的符号做了一个小的修正）。
     * @code{.bib}
     * @article{Zhou2014,
     * Title   = {Chebyshev-filtered subspace iteration method free of sparse
     *            diagonalization for solving the Kohn--Sham equation},
     * Author  = {Zhou, Yunkai and Chelikowsky, James R and Saad, Yousef},
     * Journal = {Journal of Computational Physics},
     * Year    = {2014},
     * Volume  = {274},
     * Pages   = {770--782},
     * }
     * @endcode
     *
     * @note  如果 @p tau 等于
     * <code>std::numeric_limits<double>::infinity()</code>
     * ，将不进行规范化处理。
     *
     */
    template <typename OperatorType, typename VectorType>
    void
    chebyshev_filter(VectorType &                    x,
                     const OperatorType &            H,
                     const unsigned int              n,
                     const std::pair<double, double> unwanted_spectrum,
                     const double                    tau,
                     VectorMemory<VectorType> &      vector_memory);

  } // namespace LinearAlgebra

} // namespace Utilities


 /*------------------------- Implementation ----------------------------*/ 

#ifndef DOXYGEN

namespace internal
{
  namespace UtilitiesImplementation
  {
    // We want to avoid including our own LAPACK wrapper header in any external
    // headers to avoid possible conflicts with other packages that may define
    // their own such header. At the same time we want to be able to call some
    // LAPACK functions from the template functions below. To resolve both
    // problems define some extra wrappers here that can be in the header:
    template <typename Number>
    void
    call_stev(const char            jobz,
              const types::blas_int n,
              Number *              d,
              Number *              e,
              Number *              z,
              const types::blas_int ldz,
              Number *              work,
              types::blas_int *     info);
  } // namespace UtilitiesImplementation
} // namespace internal

namespace Utilities
{
  namespace LinearAlgebra
  {
    template <typename NumberType>
    std::array<std::complex<NumberType>, 3>
    hyperbolic_rotation(const std::complex<NumberType> &  /*f*/ ,
                        const std::complex<NumberType> &  /*g*/ )
    {
      AssertThrow(false, ExcNotImplemented());
      std::array<NumberType, 3> res;
      return res;
    }



    template <typename NumberType>
    std::array<NumberType, 3>
    hyperbolic_rotation(const NumberType &f, const NumberType &g)
    {
      Assert(f != 0, ExcDivideByZero());
      const NumberType tau = g / f;
      AssertThrow(std::abs(tau) < 1.,
                  ExcMessage(
                    "real-valued Hyperbolic rotation does not exist for (" +
                    std::to_string(f) + "," + std::to_string(g) + ")"));
      const NumberType u =
        std::copysign(std::sqrt((1. - tau) * (1. + tau)),
                      f); // <-- more stable than std::sqrt(1.-tau*tau)
      std::array<NumberType, 3> csr;
      csr[0] = 1. / u;       // c
      csr[1] = csr[0] * tau; // s
      csr[2] = f * u;        // r
      return csr;
    }



    template <typename NumberType>
    std::array<std::complex<NumberType>, 3>
    givens_rotation(const std::complex<NumberType> &  /*f*/ ,
                    const std::complex<NumberType> &  /*g*/ )
    {
      AssertThrow(false, ExcNotImplemented());
      std::array<NumberType, 3> res;
      return res;
    }



    template <typename NumberType>
    std::array<NumberType, 3>
    givens_rotation(const NumberType &f, const NumberType &g)
    {
      std::array<NumberType, 3> res;
      // naive calculation for "r" may overflow or underflow:
      // c =  x / \sqrt{x^2+y^2}
      // s = -y / \sqrt{x^2+y^2}

      // See Golub 2013, Matrix computations, Chapter 5.1.8
      // Algorithm 5.1.3
      // and
      // Anderson (2000),
      // Discontinuous Plane Rotations and the Symmetric Eigenvalue Problem.
      // LAPACK Working Note 150, University of Tennessee, UT-CS-00-454,
      // December 4, 2000.
      // Algorithm 4
      // We implement the latter below:
      if (g == NumberType())
        {
          res[0] = std::copysign(1., f);
          res[1] = NumberType();
          res[2] = std::abs(f);
        }
      else if (f == NumberType())
        {
          res[0] = NumberType();
          res[1] = std::copysign(1., g);
          res[2] = std::abs(g);
        }
      else if (std::abs(f) > std::abs(g))
        {
          const NumberType tau = g / f;
          const NumberType u   = std::copysign(std::sqrt(1. + tau * tau), f);
          res[0]               = 1. / u;       // c
          res[1]               = res[0] * tau; // s
          res[2]               = f * u;        // r
        }
      else
        {
          const NumberType tau = f / g;
          const NumberType u   = std::copysign(std::sqrt(1. + tau * tau), g);
          res[1]               = 1. / u;       // s
          res[0]               = res[1] * tau; // c
          res[2]               = g * u;        // r
        }

      return res;
    }



    template <typename OperatorType, typename VectorType>
    double
    lanczos_largest_eigenvalue(const OperatorType &      H,
                               const VectorType &        v0_,
                               const unsigned int        k,
                               VectorMemory<VectorType> &vector_memory,
                               std::vector<double> *     eigenvalues)
    {
      // Do k-step Lanczos:

      typename VectorMemory<VectorType>::Pointer v(vector_memory);
      typename VectorMemory<VectorType>::Pointer v0(vector_memory);
      typename VectorMemory<VectorType>::Pointer f(vector_memory);

      v->reinit(v0_);
      v0->reinit(v0_);
      f->reinit(v0_);

      // two vectors to store diagonal and subdiagonal of the Lanczos
      // matrix
      std::vector<double> diagonal;
      std::vector<double> subdiagonal;

      // 1. Normalize input vector
      (*v)     = v0_;
      double a = v->l2_norm();
      Assert(a != 0, ExcDivideByZero());
      (*v) *= 1. / a;

      // 2. Compute f = Hv; a = f*v; f <- f - av; T(0,0)=a;
      H.vmult(*f, *v);
      a = (*f) * (*v);
      f->add(-a, *v);
      diagonal.push_back(a);

      // 3. Loop over steps
      for (unsigned int i = 1; i < k; ++i)
        {
          // 4. L2 norm of f
          const double b = f->l2_norm();
          Assert(b != 0, ExcDivideByZero());
          // 5. v0 <- v; v <- f/b
          *v0 = *v;
          *v  = *f;
          (*v) *= 1. / b;
          // 6. f = Hv; f <- f - b v0;
          H.vmult(*f, *v);
          f->add(-b, *v0);
          // 7. a = f*v; f <- f - a v;
          a = (*f) * (*v);
          f->add(-a, *v);
          // 8. T(i,i-1) = T(i-1,i) = b;  T(i,i) = a;
          diagonal.push_back(a);
          subdiagonal.push_back(b);
        }

      Assert(diagonal.size() == k, ExcInternalError());
      Assert(subdiagonal.size() == k - 1, ExcInternalError());

      // Use Lapack dstev to get ||T||_2 norm, i.e. the largest eigenvalue
      // of T
      const types::blas_int n = k;
      std::vector<double>   Z;       // unused for eigenvalues-only ("N") job
      const types::blas_int ldz = 1; // ^^   (>=1)
      std::vector<double>   work;    // ^^
      types::blas_int       info;
      // call lapack_templates.h wrapper:
      internal::UtilitiesImplementation::call_stev('N',
                                                   n,
                                                   diagonal.data(),
                                                   subdiagonal.data(),
                                                   Z.data(),
                                                   ldz,
                                                   work.data(),
                                                   &info);

      Assert(info == 0, LAPACKSupport::ExcErrorCode("dstev", info));

      if (eigenvalues != nullptr)
        {
          eigenvalues->resize(diagonal.size());
          std::copy(diagonal.begin(), diagonal.end(), eigenvalues->begin());
        }

      // note that the largest eigenvalue of T is below the largest
      // eigenvalue of the operator.
      // return ||T||_2 + ||f||_2, although it is not guaranteed to be an upper
      // bound.
      return diagonal[k - 1] + f->l2_norm();
    }


    template <typename OperatorType, typename VectorType>
    void
    chebyshev_filter(VectorType &                    x,
                     const OperatorType &            op,
                     const unsigned int              degree,
                     const std::pair<double, double> unwanted_spectrum,
                     const double                    a_L,
                     VectorMemory<VectorType> &      vector_memory)
    {
      const double a = unwanted_spectrum.first;
      const double b = unwanted_spectrum.second;
      Assert(degree > 0, ExcMessage("Only positive degrees make sense."));

      const bool scale = (a_L < std::numeric_limits<double>::infinity());
      Assert(
        a < b,
        ExcMessage(
          "Lower bound of the unwanted spectrum should be smaller than the upper bound."));

      Assert(a_L <= a || a_L >= b || !scale,
             ExcMessage(
               "Scaling point should be outside of the unwanted spectrum."));

      // Setup auxiliary vectors:
      typename VectorMemory<VectorType>::Pointer p_y(vector_memory);
      typename VectorMemory<VectorType>::Pointer p_yn(vector_memory);

      p_y->reinit(x);
      p_yn->reinit(x);

      // convenience to avoid pointers
      VectorType &y  = *p_y;
      VectorType &yn = *p_yn;

      // Below is an implementation of
      // Algorithm 3.2 in Zhou et al, Journal of Computational Physics 274
      // (2014) 770-782 with **a bugfix for sigma1**. Here is the original
      // algorithm verbatim:
      //
      // [Y]=chebyshev_filter_scaled(X, m, a, b, aL).
      // e=(b-a)/2; c=(a+b)/2; σ=e/(c-aL); τ=2/σ;
      // Y=(H∗X-c∗X)∗(σ/e);
      // for i=2 to m do
      //   σnew =1/(τ - σ);
      //   Yt =(H∗Y - c∗Y)∗(2∗σnew/e)-(σ∗σnew)∗X;
      //   X =Y; Y =Yt; σ =σnew;

      const double e     = (b - a) / 2.;
      const double c     = (a + b) / 2.;
      const double alpha = 1. / e;
      const double beta  = -c / e;

      const double sigma1 =
        e / (a_L - c); // BUGFIX which is relevant for odd degrees
      double       sigma = scale ? sigma1 : 1.;
      const double tau   = 2. / sigma;
      op.vmult(y, x);
      y.sadd(alpha * sigma, beta * sigma, x);

      for (unsigned int i = 2; i <= degree; ++i)
        {
          const double sigma_new = scale ? 1. / (tau - sigma) : 1.;
          op.vmult(yn, y);
          yn.sadd(2. * alpha * sigma_new, 2. * beta * sigma_new, y);
          yn.add(-sigma * sigma_new, x);
          x.swap(y);
          y.swap(yn);
          sigma = sigma_new;
        }

      x.swap(y);
    }

  } // namespace LinearAlgebra
} // namespace Utilities

#endif



DEAL_II_NAMESPACE_CLOSE


#endif


