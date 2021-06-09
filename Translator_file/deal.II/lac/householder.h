//include/deal.II-translator/lac/householder_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2020 by the deal.II authors
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

#ifndef dealii_householder_h
#define dealii_householder_h


#include <deal.II/base/config.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector_memory.h>

#include <cmath>
#include <vector>

DEAL_II_NAMESPACE_OPEN


// forward declarations
#ifndef DOXYGEN
template <typename number>
class Vector;
#endif

/*!   @addtogroup  Matrix2  
     * @{ 

 
*
*/


/**
 * 一个完整矩阵的QR分解。
 * 该类通过Householder算法计算给定矩阵的QR分解。然后，函数
 * least_squares() 可以用来计算给定向量  $x$  的最小化
 * $\|Ax-b\|$  的向量。 $A$
 * 的QR分解对这一目的很有用，因为最小化器由方程
 * $x=(A^TA)^{-1}A^Tb=(R^TQ^TQR)^{-1}R^TQ^Tb$
 * 给出，这很容易计算，因为 $Q$ 是一个正交矩阵，因此
 * $Q^TQ=I$  。因此，
 * $x=(R^TR)^{-1}R^TQ^Tb=R^{-1}R^{-T}R^TQ^Tb=R^{-1}Q^Tb$  . 此外， $R$
 * 是三角形的，所以将 $R^{-1}$
 * 应用于一个向量只涉及到后向或前向解。
 *
 *  <h3>Implementation details</h3> 该类实际上没有明确地将 $Q$ 和
 * $R$ 的因子存储为矩阵。它确实存储了 $R$ ，但是 $Q$
 * 因子被存储为形式为 $Q_i = I-v_i v_i^T$
 * 的Householder反射的乘积，其中向量 $v_i$
 * 是为了使它们可以存储在底层矩阵对象的下三角部分，而
 * $R$ 则存储在上三角部分。 $v_i$ 向量和 $R$
 * 矩阵现在发生了冲突，因为它们都想使用矩阵的对角线条目，但我们当然只能在这些位置存储一个。因此，
 * $(v_i)_i$ 的条目被单独存储在`对角线`成员变量中。
 *
 *
 * @note
 * 这个模板的实例化提供给<tt>  @<float@>  和  @<double@></tt>;
 * 其他的可以在应用程序中生成（见手册中的 @ref
 * Instantiations 部分）。
 *
 *
 */
template <typename number>
class Householder
{
public:
  /**
   * 声明容器尺寸类型的类型。
   *
   */
  using size_type = types::global_dof_index;

  /**
   * 创建一个空对象。
   *
   */
  Householder() = default;

  /**
   * 创建一个持有矩阵的QR分解的对象  $A$  。
   *
   */
  template <typename number2>
  Householder(const FullMatrix<number2> &A);

  /**
   * 计算给定矩阵的QR分解  $A$  。
   * 这将覆盖任何先前计算的QR分解。
   *
   */
  template <typename number2>
  void
  initialize(const FullMatrix<number2> &A);

  /**
   * 解决右手边<tt>src</tt>的最小二乘法问题。返回的标量值是近似误差的欧几里得准则。
   * @arg   @c  dst包含返回的最小二乘问题的解。      @arg   @c
   * src包含最小二乘法问题的右手边<i>b</i>。它将在算法过程中被改变，在返回时无法使用。
   *
   */
  template <typename number2>
  double
  least_squares(Vector<number2> &dst, const Vector<number2> &src) const;

  /**
   * 这个函数与前一个函数的作用相同，但对BlockVectors而言。
   *
   */
  template <typename number2>
  double
  least_squares(BlockVector<number2> &      dst,
                const BlockVector<number2> &src) const;

  /**
   * 一个对least_squares()的包装器，实现了标准的MatrixType接口。
   *
   */
  template <class VectorType>
  void
  vmult(VectorType &dst, const VectorType &src) const;

  /**
   * 一个对least_squares()的封装器，实现了与转置矩阵的乘法。
   *
   */
  template <class VectorType>
  void
  Tvmult(VectorType &dst, const VectorType &src) const;


private:
  /**
   * 存储正交变换的对角线元素。更多信息请参见类文件。
   *
   */
  std::vector<number> diagonal;

  /**
   * 内部用于Householder变换的存储。
   *
   */
  FullMatrix<double> storage;
};

 /*@}*/ 

#ifndef DOXYGEN
 /*-------------------------Inline functions -------------------------------*/ 

// QR-transformation cf. Stoer 1 4.8.2 (p. 191)

template <typename number>
template <typename number2>
void
Householder<number>::initialize(const FullMatrix<number2> &M)
{
  const size_type m = M.n_rows(), n = M.n_cols();
  storage.reinit(m, n);
  storage.fill(M);
  Assert(!storage.empty(), typename FullMatrix<number2>::ExcEmptyMatrix());
  diagonal.resize(m);

  // m > n, src.n() = m
  Assert(storage.n_cols() <= storage.n_rows(),
         ExcDimensionMismatch(storage.n_cols(), storage.n_rows()));

  for (size_type j = 0; j < n; ++j)
    {
      number2   sigma = 0;
      size_type i;
      // sigma = ||v||^2
      for (i = j; i < m; ++i)
        sigma += storage(i, j) * storage(i, j);
      // We are ready if the column is
      // empty. Are we?
      if (std::fabs(sigma) < 1.e-15)
        return;

      number2 s = (storage(j, j) < 0) ? std::sqrt(sigma) : -std::sqrt(sigma);
      //
      number2 beta = std::sqrt(1. / (sigma - s * storage(j, j)));

      // Make column j the Householder
      // vector, store first entry in
      // diagonal
      diagonal[j]   = beta * (storage(j, j) - s);
      storage(j, j) = s;

      for (i = j + 1; i < m; ++i)
        storage(i, j) *= beta;


      // For all subsequent columns do
      // the Householder reflection
      for (size_type k = j + 1; k < n; ++k)
        {
          number2 sum = diagonal[j] * storage(j, k);
          for (i = j + 1; i < m; ++i)
            sum += storage(i, j) * storage(i, k);

          storage(j, k) -= sum * this->diagonal[j];
          for (i = j + 1; i < m; ++i)
            storage(i, k) -= sum * storage(i, j);
        }
    }
}



template <typename number>
template <typename number2>
Householder<number>::Householder(const FullMatrix<number2> &M)
{
  initialize(M);
}



template <typename number>
template <typename number2>
double
Householder<number>::least_squares(Vector<number2> &      dst,
                                   const Vector<number2> &src) const
{
  Assert(!storage.empty(), typename FullMatrix<number2>::ExcEmptyMatrix());
  AssertDimension(dst.size(), storage.n());
  AssertDimension(src.size(), storage.m());

  const size_type m = storage.m(), n = storage.n();

  GrowingVectorMemory<Vector<number2>>            mem;
  typename VectorMemory<Vector<number2>>::Pointer aux(mem);
  aux->reinit(src, true);
  *aux = src;
  // m > n, m = src.n, n = dst.n

  // Multiply Q_n ... Q_2 Q_1 src
  // Where Q_i = I - v_i v_i^T
  for (size_type j = 0; j < n; ++j)
    {
      // sum = v_i^T dst
      number2 sum = diagonal[j] * (*aux)(j);
      for (size_type i = j + 1; i < m; ++i)
        sum += static_cast<number2>(storage(i, j)) * (*aux)(i);
      // dst -= v * sum
      (*aux)(j) -= sum * diagonal[j];
      for (size_type i = j + 1; i < m; ++i)
        (*aux)(i) -= sum * static_cast<number2>(storage(i, j));
    }
  // Compute norm of residual
  number2 sum = 0.;
  for (size_type i = n; i < m; ++i)
    sum += (*aux)(i) * (*aux)(i);
  AssertIsFinite(sum);

  // Compute solution
  storage.backward(dst, *aux);

  return std::sqrt(sum);
}



template <typename number>
template <typename number2>
double
Householder<number>::least_squares(BlockVector<number2> &      dst,
                                   const BlockVector<number2> &src) const
{
  Assert(!storage.empty(), typename FullMatrix<number2>::ExcEmptyMatrix());
  AssertDimension(dst.size(), storage.n());
  AssertDimension(src.size(), storage.m());

  const size_type m = storage.m(), n = storage.n();

  GrowingVectorMemory<BlockVector<number2>>            mem;
  typename VectorMemory<BlockVector<number2>>::Pointer aux(mem);
  aux->reinit(src, true);
  *aux = src;
  // m > n, m = src.n, n = dst.n

  // Multiply Q_n ... Q_2 Q_1 src
  // Where Q_i = I-v_i v_i^T
  for (size_type j = 0; j < n; ++j)
    {
      // sum = v_i^T dst
      number2 sum = diagonal[j] * (*aux)(j);
      for (size_type i = j + 1; i < m; ++i)
        sum += storage(i, j) * (*aux)(i);
      // dst -= v * sum
      (*aux)(j) -= sum * diagonal[j];
      for (size_type i = j + 1; i < m; ++i)
        (*aux)(i) -= sum * storage(i, j);
    }
  // Compute norm of residual
  number2 sum = 0.;
  for (size_type i = n; i < m; ++i)
    sum += (*aux)(i) * (*aux)(i);
  AssertIsFinite(sum);

  // backward works for
  // Vectors only, so copy
  // them before
  Vector<number2> v_dst, v_aux;
  v_dst = dst;
  v_aux = *aux;
  // Compute solution
  storage.backward(v_dst, v_aux);
  // copy the result back
  // to the BlockVector
  dst = v_dst;

  return std::sqrt(sum);
}


template <typename number>
template <class VectorType>
void
Householder<number>::vmult(VectorType &dst, const VectorType &src) const
{
  least_squares(dst, src);
}


template <typename number>
template <class VectorType>
void
Householder<number>::Tvmult(VectorType &, const VectorType &) const
{
  Assert(false, ExcNotImplemented());
}



#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif


