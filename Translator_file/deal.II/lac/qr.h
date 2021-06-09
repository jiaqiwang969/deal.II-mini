//include/deal.II-translator/lac/qr_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2020 by the deal.II authors
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

#ifndef dealii_qr_h
#define dealii_qr_h

#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>

#include <deal.II/lac/lapack_full_matrix.h>
#include <deal.II/lac/utilities.h>

#include <boost/signals2/signal.hpp>

#include <memory>

DEAL_II_NAMESPACE_OPEN

/**
 * 用于薄型QR实现的基类。
 * 这个类和从它派生出来的类是为了一次建立 $Q$ 和 $R$
 * 矩阵的行/列，即通过从一个空的 $0\times 0$ 矩阵增长到
 * $N\times N$ ，其中 $N$ 是增加的列向量数量。
 * 因此，与每个向量有相同行数的矩阵（即 $Q$
 * 矩阵）被存储为`VectorType`的向量集合。
 *
 *
 */
template <typename VectorType>
class BaseQR
{
  /**
   * R矩阵的数字类型。
   *
   */
  using Number = typename VectorType::value_type;

protected:
  /**
   * 默认的私有构造函数。
   *
   */
  BaseQR();

public:
  /**
   * 解构器。
   *
   */
  virtual ~BaseQR() = default;

  /**
   * 将 @p column 追加到QR因式分解中。
   * 如果结果是成功的，即列是线性独立的，则返回
   * <code>true</code> 。否则 @p column 被拒绝，返回值为
   * <code>false</code>  。
   *
   */
  virtual bool
  append_column(const VectorType &column) = 0;

  /**
   * 删除一列 @p k 并更新QR分解。
   *
   */
  virtual void
  remove_column(const unsigned int k = 0) = 0;

  /**
   * 返回子空间的大小。
   *
   */
  unsigned int
  size() const;

  /**
   * 返回当前的上三角矩阵R。
   *
   */
  const LAPACKFullMatrix<Number> &
  get_R() const;

  /**
   * 解决  $Rx=y$  . 矢量 @p x 和 @p y
   * 应该与子空间的当前大小一致。  如果 @p transpose 是
   * <code>true</code> ，则 $R^Tx=y$ 被解决。
   *
   */
  void
  solve(Vector<Number> &      x,
        const Vector<Number> &y,
        const bool            transpose = false) const;

  /**
   * 设置  $y = Qx$  。 $x$ 的大小应与R矩阵的大小一致。
   *
   */
  virtual void
  multiply_with_Q(VectorType &y, const Vector<Number> &x) const = 0;

  /**
   * 设置  $y = Q^Tx$  .  $x$  的大小应与列向量的大小一致。
   *
   */
  virtual void
  multiply_with_QT(Vector<Number> &y, const VectorType &x) const = 0;

  /**
   * 设置  $y = QRx$  .  $x$ 的大小应与R矩阵的大小一致。
   *
   */
  virtual void
  multiply_with_A(VectorType &y, const Vector<Number> &x) const = 0;

  /**
   * 设置  $y = R^T Q^Tx$  .  $x$
   * 的大小应与列向量的大小一致。
   *
   */
  virtual void
  multiply_with_AT(Vector<Number> &y, const VectorType &x) const = 0;

  /**
   * 连接一个槽，以便在进行吉文斯旋转时检索到一个通知。
   * 该函数需要两个指数， @p i 和 @p j,
   * 描述旋转的平面，以及三组数字 @p csr
   * （余弦、正弦和半径，见
   * Utilities::LinearAlgebra::givens_rotation()) ，代表旋转矩阵。
   *
   */
  boost::signals2::connection
  connect_givens_slot(
    const std::function<void(const unsigned int           i,
                             const unsigned int           j,
                             const std::array<Number, 3> &csr)> &slot);

protected:
  /**
   * 计算 $y=Hx$ ，其中 $H$
   * 是由该对象存储的列向量形成的矩阵。
   *
   */
  void
  multiply_with_cols(VectorType &y, const Vector<Number> &x) const;

  /**
   * 与存储在该对象中的转置列相乘。
   *
   */
  void
  multiply_with_colsT(Vector<Number> &y, const VectorType &x) const;

  /**
   * 存储列的唯一指针的向量。
   *
   */
  std::vector<std::unique_ptr<VectorType>> columns;

  /**
   * 用于存储R的矩阵。
   *
   */
  LAPACKFullMatrix<Number> R;

  /**
   * 当前大小（Q中的列数）。
   *
   */
  unsigned int current_size;

  /**
   * 当吉文斯旋转在`(i,j)`面进行时，用于检索通知的信号。
   *
   */
  boost::signals2::signal<void(const unsigned int i,
                               const unsigned int j,
                               const std::array<Number, 3> &)>
    givens_signal;
};

// clang-format off
/**
 * 一个用于计算和存储由一组列向量代表的矩阵的QR分解的类。
 * 该类被设计用来更新矩阵 $A$
 * 的给定（可能是空的）QR因式分解（通过提供其列而逐步构建），由于向
 * $A$
 * 添加了一个新的列向量。这相当于通过Gram-Schmidt程序构建一个正态基。该类还提供了当第一列被移除时的更新功能。
 * `VectorType`模板参数可以是并行和串行矢量，只需要有基本操作，如加法、标量乘积等。它还需要有一个复制构造器。
 * 见第335-338页的6.5.2-6.5.3节。
 *
 * @code{.bib}
 * @Book{Golub2013,
 * title     = {Matrix computations},
 * publisher = {Johns Hopkins University Press},
 * year      = {2013},
 * author    = {Golub, Gene H and Van Loan, Charles F},
 * edition   = {4},
 * }
 * @endcode
 * 以及
 *
 * @code{.bib}
 * @article{Daniel1976,
 * author   = {Daniel, James W and Gragg, Walter Bill and Kaufman, Linda and Stewart, Gilbert W},
 * title    = {{Reorthogonalization and stable algorithms for updating the Gram-Schmidt QR factorization}},
 * journal  = {Mathematics of Computation},
 * year     = {1976},
 * volume   = {30},
 * number   = {136},
 * pages    = {772--795},
 * }
 * @Article{Reichel1990,
 * author     = {Reichel, L. and Gragg, W. B.},
 * title      = {{Algorithm 686: FORTRAN Subroutines for Updating the QR Decomposition}},
 * journal    = {ACM Trans. Math. Softw.},
 * year       = {1990},
 * volume     = {16},
 * number     = {4},
 * pages      = {369--377},
 * month      = dec,
 * issn       = {0098-3500},
 * acmid      = {98291},
 * address    = {New York, NY, USA},
 * doi        = {10.1145/98267.98291},
 * issue_date = {Dec. 1990},
 * numpages   = {9},
 * publisher  = {ACM},
 * url        = {http://doi.acm.org/10.1145/98267.98291},
 * }
 * @endcode
 *
 *
 */
// clang-format on
template <typename VectorType>
class QR : public BaseQR<VectorType>
{
public:
  /**
   * R矩阵的数字类型。
   *
   */
  using Number = typename VectorType::value_type;

  /**
   * 默认构造函数。
   *
   */
  QR();

  /**
   * 解构器。
   *
   */
  virtual ~QR() = default;

  /**
   * @copydoc   BaseQR::append_column
   * @note  目前这个函数总是返回  <code>true</code>  。
   *
   */
  virtual bool
  append_column(const VectorType &column);

  /**
   * 删除第一列并更新QR分解。    从给定的QR分解 $QR= A =
   * [a_1\,\dots a_n], \quad a_i \in {\mathbb R}^m$
   * 开始，我们的目标是计算 $\tilde Q \tilde R= \tilde A =
   * [a_2\,\dots a_n], \quad a_i \in {\mathbb R}^m$ 的因式分解。
   * 标准方法是将 $R$ 划分为\f[ R = \begin{bmatrix} r_{11} & w^T \\ 0
   * & R_{33} \end{bmatrix} \f]，然后得出\f[ Q^T \tilde A =
   * \begin{bmatrix} 0 & w^T \\ 0 & R_{33} \end{bmatrix}
   * \f]是上海森堡，其中不需要的对角线元素可以通过一连串的吉文斯旋转来归零。
   * 请注意， $\tilde R^T \tilde R = \tilde A^T \tilde A$
   * ，其中RHS包括在 $A^T A = R^T R$  中。因此 $\tilde R$
   * 可以通过Cholesky分解得到。
   *
   */
  virtual void
  remove_column(const unsigned int k = 0);

  virtual void
  multiply_with_Q(VectorType &y, const Vector<Number> &x) const;

  virtual void
  multiply_with_QT(Vector<Number> &y, const VectorType &x) const;

  virtual void
  multiply_with_A(VectorType &y, const Vector<Number> &x) const;

  virtual void
  multiply_with_AT(Vector<Number> &y, const VectorType &x) const;

private:
  /**
   * 对 @p Q 和 @p R 进行`(i,j)`面的基文旋转，使
   * <code>R(k,k)</code> 被归零。    见Golub
   * 2013的第5.1.9章，矩阵计算。
   *
   */
  void
  apply_givens_rotation(const unsigned int i, const unsigned int k);

  /**
   * 对Q进行吉文斯旋转所需的临时向量。
   *
   */
  VectorType tmp;
};



/**
 * 获得 $R$ 因式分解的三角 $A=QR$ 矩阵以及矩阵 $A$
 * 本身的一个类。正交矩阵 $Q$
 * 没有被明确存储，该类的名称。与 $Q$ 的乘法可以表示为
 * $Q=A R^{-1}$ ，而与 $Q^T$ 的乘法则由 $Q^T=R^{-T}A^T$ 给出。
 * 该类旨在更新一个给定的（可能是空的）QR因式分解，由于增加了一个新的列向量。这相当于通过Gram-Schmidt程序构建一个正态基础。当列被移除时，该类也提供更新功能。
 * `VectorType`模板参数可以是并行和串行矢量，只需要有基本的操作，如加法、标量积等。它还需要有一个复制构造器。
 *
 *
 */
template <typename VectorType>
class ImplicitQR : public BaseQR<VectorType>
{
public:
  /**
   * R矩阵的数字类型。
   *
   */
  using Number = typename VectorType::value_type;

  /**
   * 默认构造函数。
   *
   */
  ImplicitQR();

  /**
   * 解构器。
   *
   */
  virtual ~ImplicitQR() = default;

  virtual bool
  append_column(const VectorType &column);

  /**
   * 删除列并更新QR分解。    从给定的QR分解 $QR= A =
   * [a_1\,\dots a_n], \quad a_i \in R^m$ 开始，我们旨在计算 $\tilde
   * Q \tilde R= \tilde A = [a_2\,\dots a_n], \quad a_i \in R^m$
   * 的因式分解。    请注意， $\tilde R^T \tilde R = \tilde A^T
   * \tilde A$  ，其中RHS包括在 $A^T A = R^T R$  中。因此 $\tilde R$
   * 可以通过Cholesky分解得到。
   *
   */
  virtual void
  remove_column(const unsigned int k = 0);

  virtual void
  multiply_with_Q(VectorType &y, const Vector<Number> &x) const;

  virtual void
  multiply_with_QT(Vector<Number> &y, const VectorType &x) const;

  virtual void
  multiply_with_A(VectorType &y, const Vector<Number> &x) const;

  virtual void
  multiply_with_AT(Vector<Number> &y, const VectorType &x) const;

  /**
   * 连接一个槽，以实现在添加一列时对线性依赖的自定义检查。
   * 这里， @p u 是待R矩阵的最后一列， @p rho 是其对角线，
   * @p col_norm_sqr 是该列的 $l2$ 准则的平方。
   * 如果新列是线性独立的，该函数应该返回 <code>true</code>
   * 。
   *
   */
  boost::signals2::connection
  connect_append_column_slot(
    const std::function<bool(const Vector<Number> &u,
                             const Number &        rho2,
                             const Number &        col_norm_sqr)> &slot);

private:
  /**
   * 在"(i,k) "平面上应用给定的旋转，将 $R(k,k)$ 归零。
   *
   */
  void
  apply_givens_rotation(const unsigned int i, const unsigned int k);

  /**
   * 用来决定新列是否是线性依赖的信号。    这里， @p u
   * 是要成为R矩阵的最后一列， @p rho 是其对角线， @p
   * col_norm_sqr 是该列的 $l2$ 准则的平方。
   * 如果新列是线性独立的，该函数应该返回 <code>true</code>
   * 。
   *
   */
  boost::signals2::signal<bool(const Vector<Number> &u,
                               const Number &        rho,
                               const Number &        col_norm_sqr)>
    column_signal;
};

// -------------------  inline and template functions ----------------
#ifndef DOXYGEN

namespace internal
{
  namespace QRImplementation
  {
    // We want to avoid including our own LAPACK wrapper header in any external
    // headers to avoid possible conflicts with other packages that may define
    // their own such header. At the same time we want to be able to call some
    // LAPACK functions from the template functions below. To resolve both
    // problems define some extra wrappers here that can be in the header:
    template <typename Number>
    void
    call_trmv(const char            uplo,
              const char            trans,
              const char            diag,
              const types::blas_int n,
              const Number *        a,
              const types::blas_int lda,
              Number *              x,
              const types::blas_int incx);

    template <typename Number>
    void
    call_trtrs(const char            uplo,
               const char            trans,
               const char            diag,
               const types::blas_int n,
               const types::blas_int nrhs,
               const Number *        a,
               const types::blas_int lda,
               Number *              b,
               const types::blas_int ldb,
               types::blas_int *     info);
  } // namespace QRImplementation
} // namespace internal



template <typename VectorType>
BaseQR<VectorType>::BaseQR()
  : current_size(0)
{
  R.set_property(LAPACKSupport::upper_triangular);
}



template <typename VectorType>
unsigned int
BaseQR<VectorType>::size() const
{
  return current_size;
}



template <typename VectorType>
const LAPACKFullMatrix<typename BaseQR<VectorType>::Number> &
BaseQR<VectorType>::get_R() const
{
  return R;
}



template <typename VectorType>
void
BaseQR<VectorType>::solve(Vector<Number> &      x,
                          const Vector<Number> &y,
                          const bool            transpose) const
{
  Assert(x.size() == this->current_size,
         ExcDimensionMismatch(x.size(), this->current_size));
  Assert(y.size() == this->current_size,
         ExcDimensionMismatch(y.size(), this->current_size));

  // copy if the two vectors are not the same
  if (&x != &y)
    x = y;

  const int lda   = this->current_size;
  const int ldb   = this->current_size;
  const int N     = this->current_size;
  const int n_rhs = 1;
  int       info  = 0;
  internal::QRImplementation::call_trtrs('U',
                                         transpose ? 'T' : 'N',
                                         'N',
                                         N,
                                         n_rhs,
                                         &this->R(0, 0),
                                         lda,
                                         &x(0),
                                         ldb,
                                         &info);
}



template <typename VectorType>
void
BaseQR<VectorType>::multiply_with_cols(VectorType &          y,
                                       const Vector<Number> &x) const
{
  Assert(x.size() == this->current_size,
         ExcDimensionMismatch(x.size(), this->current_size));

  y = 0.;
  for (unsigned int j = 0; j < this->current_size; ++j)
    y.add(x[j], *this->columns[j]);
}



template <typename VectorType>
void
BaseQR<VectorType>::multiply_with_colsT(Vector<Number> &  y,
                                        const VectorType &x) const
{
  Assert(y.size() == this->current_size,
         ExcDimensionMismatch(y.size(), this->current_size));

  for (unsigned int j = 0; j < this->current_size; ++j)
    y[j] = (*this->columns[j]) * x;
}



template <class VectorType>
boost::signals2::connection
BaseQR<VectorType>::connect_givens_slot(
  const std::function<void(const unsigned int i,
                           const unsigned int j,
                           const std::array<Number, 3> &)> &slot)
{
  return givens_signal.connect(slot);
}



template <class VectorType>
boost::signals2::connection
ImplicitQR<VectorType>::connect_append_column_slot(
  const std::function<bool(const Vector<Number> &u,
                           const Number &        rho,
                           const Number &        col_norm_sqr)> &slot)
{
  return column_signal.connect(slot);
}



template <typename VectorType>
ImplicitQR<VectorType>::ImplicitQR()
  : BaseQR<VectorType>()
{}



template <typename VectorType>
bool
ImplicitQR<VectorType>::append_column(const VectorType &column)
{
  if (this->current_size == 0)
    {
      this->R.grow_or_shrink(this->current_size + 1);
      this->columns.push_back(std::make_unique<VectorType>(column));
      this->R(0, 0) = column.l2_norm();
      ++this->current_size;
    }
  else
    {
      // first get scalar products with A^T
      Vector<Number> u(this->current_size);
      this->multiply_with_AT(u, column);

      // now solve R^T x = (A^T * column)
      const int lda   = this->current_size;
      const int ldb   = this->current_size;
      const int N     = this->current_size;
      const int n_rhs = 1;
      int       info  = 0;
      internal::QRImplementation::call_trtrs(
        'U', 'T', 'N', N, n_rhs, &this->R(0, 0), lda, &u(0), ldb, &info);

      // finally get the diagonal element:
      // rho2 = |column|^2 - |u|^2
      const Number column_norm_sqr = column.norm_sqr();
      const Number rho2            = column_norm_sqr - u.norm_sqr();
      const bool   linearly_independent =
        column_signal.empty() ? rho2 > 0 :
                                column_signal(u, rho2, column_norm_sqr).get();

      // bail out if it turns out to be linearly dependent
      if (!linearly_independent)
        return false;

      // at this point we update is successful and we can enlarge R
      // and store the column:
      this->columns.push_back(std::make_unique<VectorType>(column));
      this->R.grow_or_shrink(this->current_size + 1);
      this->R(this->current_size, this->current_size) = std::sqrt(rho2);
      for (unsigned int i = 0; i < this->current_size; ++i)
        this->R(i, this->current_size) = u(i);

      this->current_size++;
    }

  return true;
}



template <typename VectorType>
void
ImplicitQR<VectorType>::apply_givens_rotation(const unsigned int i,
                                              const unsigned int k)
{
  AssertIndexRange(i, k);
  AssertIndexRange(k, this->current_size);
  const std::array<Number, 3> csr =
    dealii::Utilities::LinearAlgebra::givens_rotation<Number>(this->R(i, k),
                                                              this->R(k, k));

  // first, set k'th column:
  this->R(i, k) = csr[2];
  this->R(k, k) = 0.;
  // now do the rest:
  for (unsigned int j = 0; j < this->R.n(); ++j)
    if (j != k)
      {
        const Number t = this->R(i, j);
        this->R(i, j)  = csr[0] * this->R(i, j) + csr[1] * this->R(k, j);
        this->R(k, j)  = -csr[1] * t + csr[0] * this->R(k, j);
      }

  if (!this->givens_signal.empty())
    this->givens_signal(i, k, csr);
}



template <typename VectorType>
void
ImplicitQR<VectorType>::remove_column(const unsigned int k)
{
  // before actually removing a column from Q and resizing R,
  // apply givens rotations to bring H into upper triangular form:
  for (unsigned int j = k + 1; j < this->R.n(); ++j)
    {
      const unsigned int i = j - 1;
      apply_givens_rotation(i, j);
    }

  // remove last row and k-th column
  --this->current_size;
  this->R.remove_row_and_column(this->current_size, k);

  // Finally remove the column from A
  this->columns.erase(this->columns.begin() + k);
}



template <typename VectorType>
void
ImplicitQR<VectorType>::multiply_with_Q(VectorType &          y,
                                        const Vector<Number> &x) const
{
  // A = QR
  // A R^{-1} = Q
  Vector<Number> x1 = x;
  BaseQR<VectorType>::solve(x1, x1, false);
  multiply_with_A(y, x1);
}



template <typename VectorType>
void
ImplicitQR<VectorType>::multiply_with_QT(Vector<Number> &  y,
                                         const VectorType &x) const
{
  // A = QR
  // A^T = R^T Q^T
  // {R^T}^{-1} A^T = Q^T
  multiply_with_AT(y, x);
  BaseQR<VectorType>::solve(y, y, true);
}



template <typename VectorType>
void
ImplicitQR<VectorType>::multiply_with_A(VectorType &          y,
                                        const Vector<Number> &x) const
{
  BaseQR<VectorType>::multiply_with_cols(y, x);
}



template <typename VectorType>
void
ImplicitQR<VectorType>::multiply_with_AT(Vector<Number> &  y,
                                         const VectorType &x) const
{
  BaseQR<VectorType>::multiply_with_colsT(y, x);
}



template <typename VectorType>
QR<VectorType>::QR()
  : BaseQR<VectorType>()
{}



template <typename VectorType>
bool
QR<VectorType>::append_column(const VectorType &column)
{
  // resize R:
  this->R.grow_or_shrink(this->current_size + 1);
  this->columns.push_back(std::make_unique<VectorType>(column));

  // now a Gram-Schmidt part: orthonormalize the new column
  // against everything we have so far:
  auto &last_col = *this->columns.back();
  for (unsigned int i = 0; i < this->current_size; ++i)
    {
      const auto &i_col              = *this->columns[i];
      this->R(i, this->current_size) = i_col * last_col;
      last_col.add(-this->R(i, this->current_size), i_col);
    }

  this->R(this->current_size, this->current_size) = last_col.l2_norm();

  Assert(this->R(this->current_size, this->current_size) > 0.,
         ExcDivideByZero());
  last_col *= 1. / this->R(this->current_size, this->current_size);

  ++this->current_size;
  return true;
}



template <typename VectorType>
void
QR<VectorType>::apply_givens_rotation(const unsigned int i,
                                      const unsigned int k)
{
  AssertIndexRange(i, k);
  AssertIndexRange(k, this->current_size);
  const std::array<Number, 3> csr =
    dealii::Utilities::LinearAlgebra::givens_rotation<Number>(this->R(i, k),
                                                              this->R(k, k));

  // first, set k'th column:
  this->R(i, k) = csr[2];
  this->R(k, k) = 0.;
  // now do the rest:
  for (unsigned int j = 0; j < this->R.n(); ++j)
    if (j != k)
      {
        const Number t = this->R(i, j);
        this->R(i, j)  = csr[0] * this->R(i, j) + csr[1] * this->R(k, j);
        this->R(k, j)  = -csr[1] * t + csr[0] * this->R(k, j);
      }

  // now adjust i,k columns due to multiplication with the
  // transpose Givens matrix from right:
  auto &col_i = *this->columns[i];
  auto &col_k = *this->columns[k];
  // save column i:
  tmp = col_i;
  col_i.sadd(csr[0], csr[1], col_k);
  col_k.sadd(csr[0], -csr[1], tmp);

  if (!this->givens_signal.empty())
    this->givens_signal(i, k, csr);
}



template <typename VectorType>
void
QR<VectorType>::remove_column(const unsigned int k)
{
  AssertIndexRange(k, this->current_size);
  Assert(this->current_size > 0,
         ExcMessage("Can not remove a column if QR is empty"));
  // apply a sequence of Givens rotations
  // see section 6.5 "Updating matrix factorizations" in Golub 2013, Matrix
  // computations

  // So we want to have QR for \tilde A \in R^{m*(n-1)}
  // if we remove the column k, we end up with upper Hessenberg matrix
  //      x x x x x
  //        x x x x
  // H =      x x x
  //          x x x
  //            x x
  //              x
  // where k = 2 (3rd column), m = 7, n = 6
  //
  // before actually removing a column from Q and resizing R,
  // apply givens rotations to bring H into upper triangular form:
  for (unsigned int j = k + 1; j < this->R.n(); ++j)
    {
      const unsigned int i = j - 1;
      apply_givens_rotation(i, j);
    }

  // now we can throw away the column from Q and adjust R
  // since we do thin-QR, after Givens rotations we need to throw
  // away the last column:
  const unsigned int size_minus_1 = this->columns.size() - 1;
  this->columns.erase(this->columns.begin() + size_minus_1);

  // remove last row and k-th column
  --this->current_size;
  this->R.remove_row_and_column(this->current_size, k);
}



template <typename VectorType>
void
QR<VectorType>::multiply_with_Q(VectorType &y, const Vector<Number> &x) const
{
  BaseQR<VectorType>::multiply_with_cols(y, x);
}



template <typename VectorType>
void
QR<VectorType>::multiply_with_QT(Vector<Number> &y, const VectorType &x) const
{
  BaseQR<VectorType>::multiply_with_colsT(y, x);
}



template <typename VectorType>
void
QR<VectorType>::multiply_with_A(VectorType &y, const Vector<Number> &x) const
{
  Vector<Number> x1   = x;
  const int      N    = this->current_size;
  const int      lda  = N;
  const int      incx = 1;
  internal::QRImplementation::call_trmv(
    'U', 'N', 'N', N, &this->R(0, 0), lda, &x1[0], incx);

  multiply_with_Q(y, x1);
}



template <typename VectorType>
void
QR<VectorType>::multiply_with_AT(Vector<Number> &y, const VectorType &x) const
{
  multiply_with_QT(y, x);

  const int N    = this->current_size;
  const int lda  = N;
  const int incx = 1;
  internal::QRImplementation::call_trmv(
    'U', 'T', 'N', N, &this->R(0, 0), lda, &y[0], incx);
}

#endif // no DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif


