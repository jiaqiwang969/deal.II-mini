//include/deal.II-translator/lac/tensor_product_matrix_0.txt
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

#ifndef dealii_tensor_product_matrix_h
#define dealii_tensor_product_matrix_h


#include <deal.II/base/config.h>

#include <deal.II/base/array_view.h>
#include <deal.II/base/thread_management.h>

#include <deal.II/lac/lapack_full_matrix.h>

#include <deal.II/matrix_free/tensor_product_kernels.h>

DEAL_II_NAMESPACE_OPEN

// Forward declarations
#ifndef DOXYGEN
template <typename>
class Vector;
template <typename>
class FullMatrix;
#endif

/**
 * 这是一个抽象的基类，用于一个特殊的矩阵类，即TensorProductMatrixSymmetricSum。
 * 首先，这个基类就像一个容器，存储一维质量矩阵和一维导数矩阵，以及每个张量方向的广义特征值和特征向量。对于这些矩阵和相应的广义特征问题的详细定义，我们可以参考TensorProductMatrixSymmetricSum的主要文档。
 *
 *
 * @note
 * 这个基类没有计算所给质量和导数矩阵的特征值和特征向量的功能。初始化数据成员的责任完全在于派生类。
 * 其次，它用张量乘积矩阵（vmult()）和它的逆（apply_inverse()）实现了矩阵-向量乘积，如TensorProductMatrixSymmetricSum的主要文档中所述。
 *
 *
 * @note
 * 该类使用一个临时数组来存储中间结果，它是一个类成员。一个mutex被用来保护对这个数组的访问，并确保正确的结果。如果几个线程并行运行这个类的实例，建议每个线程持有自己的矩阵版本。
 * @tparam  dim
 * 问题的尺寸。目前，一维、二维和三维代码已经实现。
 * @tparam  Number 底层数组元素的算术类型。
 * @tparam  n_rows_1d
 * 编译时一维矩阵的行数（只有在每个维度的行数和列数重合时才有效）。默认情况下，在
 *
 * - ，这意味着行数是在运行时通过传递给reinit()函数的矩阵来确定。
 *
 *
 */
template <int dim, typename Number, int n_rows_1d = -1>
class TensorProductMatrixSymmetricSumBase
{
public:
  /**
   * 矩阵条目的类型。这个别名类似于标准库容器中的<tt>value_type</tt>。
   *
   */
  using value_type = Number;

  /**
   * 一维矩阵的静态行数。详情见模板参数<tt>n_rows_1d</tt>的描述。
   *
   */
  static constexpr int n_rows_1d_static = n_rows_1d;

  /**
   * 返回由一维矩阵的克朗克积产生的张量积矩阵的行数，在TensorProductMatrixSymmetricSum的主文档中有所描述。
   *
   */
  unsigned int
  m() const;

  /**
   * 返回由一维矩阵的克朗克积产生的张量积矩阵的列数，这在TensorProductMatrixSymmetricSum的主要文档中描述。
   *
   */
  unsigned int
  n() const;

  /**
   * 实现与底层矩阵的矩阵-向量乘积，在TensorProductMatrixSymmetricSum的主要文档中有所描述。
   * 这个函数在ArrayView上操作，以允许检查与 @p dst 和 @p src.
   * 有关的阵列界限。
   *
   */
  void
  vmult(const ArrayView<Number> &dst, const ArrayView<const Number> &src) const;

  /**
   * 如TensorProductMatrixSymmetricSum的主文档中所述，实现了与底层矩阵的矩阵-向量乘积。
   * 这个函数在ArrayView上操作，以允许检查与 @p dst 和 @p src.
   * 有关的阵列界限。
   *
   */
  void
  apply_inverse(const ArrayView<Number> &      dst,
                const ArrayView<const Number> &src) const;

protected:
  /**
   * 默认构造函数。
   *
   */
  TensorProductMatrixSymmetricSumBase() = default;

  /**
   * 一个包含每个张量方向的质量矩阵的数组。
   *
   */
  std::array<Table<2, Number>, dim> mass_matrix;

  /**
   * 一个包含每个张量方向的导数矩阵的数组。
   *
   */
  std::array<Table<2, Number>, dim> derivative_matrix;

  /**
   * 一个存储每个张量方向的广义特征值的数组。
   *
   */
  std::array<AlignedVector<Number>, dim> eigenvalues;

  /**
   * 一个存储每个张量方向的广义特征向量的数组。
   *
   */
  std::array<Table<2, Number>, dim> eigenvectors;

private:
  /**
   * 一个用于临时数据的数组。
   *
   */
  mutable AlignedVector<Number> tmp_array;

  /**
   * 一个保护对数组访问的突变器  @p tmp_array.  。
   *
   */
  mutable Threads::Mutex mutex;
};



/**
 * 这是一个特殊的矩阵类，定义为一维矩阵的张量积（或克朗克积），其类型为
 *
 * @f{align*}{
 * L &= A_1 \otimes M_0 + M_1 \otimes A_0
 * @f}
 * 在二维和
 *
 * @f{align*}{
 * L &= A_2 \otimes M_1 \otimes M_0 + M_2 \otimes A_1 \otimes M_0 + M_2 \otimes
 * M_1 \otimes A_0
 * @f}
 * 在三维中。典型的应用设置是拉普拉斯 $L$
 * 在笛卡尔（轴对齐）几何上的离散化，它可以准确地由一维质量矩阵
 * $M$ 和一维拉普拉斯矩阵 $A$
 * 在每个张量方向的克朗克或张量乘积来表示（由于对称性
 * $M$ 和 $A$
 * 在每个维度都相同）。结果类的维度是一维矩阵的乘积。
 * 这个类实现了两个基本操作，即通常的向量乘法和逆向操作。对于这两种操作，可以应用快速的张量技术，在
 * $\text{size}(M)^{d+1}$
 * 的算术运算中实现运算符的评估，大大低于
 * $\text{size}(M)^{2d}$ 的天真正向变换和 $\text{size}(M)^{3d}$
 * 的设置逆运算。 有趣的是，由于1964年<a
 * href="http://dl.acm.org/citation.cfm?id=2716130">R. E. Lynch, J. R. Rice,
 * D. H. Thomas, Direct solution of partial difference equations by tensor
 * product methods, Numerische Mathematik 6,
 * 185-199</a>的一篇文章，矩阵 $L$
 * 的精确逆可以通过张量积找到。
 *
 * @f{align*}{
 * L^{-1} &= S_1 \otimes S_0 (\Lambda_1 \otimes I + I \otimes \Lambda_0)^{-1}
 * S_1^\mathrm T \otimes S_0^\mathrm T,
 * @f}
 * 其中 $S_d$
 * 是给定张量方向上广义特征值问题的特征向量矩阵 $d$  。
 *
 * @f{align*}{
 * A_d s  &= \lambda M_d s, d = 0, \quad \ldots,\mathrm{dim},
 * @f}
 * 而 $\Lambda_d$ 是代表广义特征值的对角线矩阵  $\lambda$
 * 。请注意，向量 $s$ 是这样的：它们同时对角化 $A_d$ 和
 * $M_d$ ，即 $S_d^{\mathrm T} A_d S_d = \Lambda_d$ 和 $S_d^{\mathrm T} M_d
 * S_d = I$  。这种矩阵反演的方法被称为快速对角线化法。
 * 这个类需要LAPACK支持。
 * 请注意，这个类允许两种使用模式。第一种是对矩阵维度使用运行时常数的使用模式，通过设置可选的模板参数<tt>n_rows_1d</tt>来实现。
 *
 * - . 第二种是更快的使用模式，允许将模板参数设置为编译时常量，特别是对于小尺寸的矩阵，给出明显更快的代码。
 * @tparam  dim
 * 问题的尺寸。目前，已经实现了一维、二维和三维代码。
 * @tparam  Number
 * 基础数组元素的算术类型。请注意，底层的LAPACK实现只支持浮点数和双数，所以目前通用类只支持这两种类型。然而，存在一个针对矢量化类型VectorizedArray<float>和VectorizedArray<double>的模板特化。这对于为每个矢量化通道执行LAPACK计算是必要的，也就是说，对于支持的浮点数和双数。
 * @tparam  n_rows_1d
 * 编译时一维矩阵的行数（只有在每个维度的行数和列数重合时才有效）。默认情况下，在
 *
 * - ，这意味着行数是在运行时通过传递给reinit()函数的矩阵来确定。
 *
 *
 */
template <int dim, typename Number, int n_rows_1d = -1>
class TensorProductMatrixSymmetricSum
  : public TensorProductMatrixSymmetricSumBase<dim, Number, n_rows_1d>
{
public:
  /**
   * 默认构造函数。
   *
   */
  TensorProductMatrixSymmetricSum() = default;

  /**
   * 相当于空构造函数的构造函数，并立即调用 reinit(
   * std::array<Table<2,Number>,  dim>&,const  std::array<Table<2,Number>,
   * dim>&) 。
   *
   */
  TensorProductMatrixSymmetricSum(
    const std::array<Table<2, Number>, dim> &mass_matrix,
    const std::array<Table<2, Number>, dim> &derivative_matrix);

  /**
   * 等同于空构造函数的构造函数，并立即调用 reinit(const
   * std::array<FullMatrix<Number>,dim>&,const
   * std::array<FullMatrix<Number>,dim>&). )
   *
   */
  TensorProductMatrixSymmetricSum(
    const std::array<FullMatrix<Number>, dim> &mass_matrix,
    const std::array<FullMatrix<Number>, dim> &derivative_matrix);

  /**
   * 相当于空构造函数并立即调用 reinit(const
   * Table<2,Number>&,const Table<2,Number>&) 的构造函数。
   *
   */
  TensorProductMatrixSymmetricSum(const Table<2, Number> &mass_matrix,
                                  const Table<2, Number> &derivative_matrix);

  /**
   * 通过将一维质量矩阵 @p mass_matrix 和一维导数矩阵 @p
   * derivative_matrix
   * 的数组分别复制到它的基类对应物中，并通过在
   * TensorProductMatrixSymmetricSumBase::eigenvalues 和
   * TensorProductMatrixSymmetricSumBase::eigenvectors,
   * 中分别组装关于广义特征值和特征向量来初始化张量积矩阵。
   * 注意，目前的实现要求每个 $M_{d}$
   * 都是对称的和正定的，每个 $A_{d}$
   * 都是对称的和可逆的，但不一定是正定的。
   *
   */
  void
  reinit(const std::array<Table<2, Number>, dim> &mass_matrix,
         const std::array<Table<2, Number>, dim> &derivative_matrix);

  /**
   * 这个函数等同于之前的 reinit()，只是 @p mass_matrix 和 @p
   * derivative_matrix 中的一维矩阵分别以FullMatrix的形式传递。
   *
   */
  void
  reinit(const std::array<FullMatrix<Number>, dim> &mass_matrix,
         const std::array<FullMatrix<Number>, dim> &derivative_matrix);

  /**
   * 这个函数等同于第一个reinit()，只是我们为每个张量方向考虑相同的一维质量矩阵
   * @p mass_matrix 和相同的一维导数矩阵 @p derivative_matrix 。
   *
   */
  void
  reinit(const Table<2, Number> &mass_matrix,
         const Table<2, Number> &derivative_matrix);

private:
  /**
   * 一个基于完美转发的所有reinit()函数的通用实现，允许传递lvalue以及rvalue参数。
   * @tparam  MatrixArray 必须可以转换为
   * TensorProductMatrixSymmetricSumBase::mass_matrix 和
   * TensorProductMatrixSymmetricSumBase::derivative_matrix.
   * 的基本类型。
   *
   */
  template <typename MatrixArray>
  void
  reinit_impl(MatrixArray &&mass_matrix, MatrixArray &&derivative_matrix);
};



/**
 * 这是VectorizedArray<Number>的模板特化，是算术模板。详细说明见通用TensorProductMatrixSymmetricSum类的主文件。
 *
 *
 */
template <int dim, typename Number, int n_rows_1d>
class TensorProductMatrixSymmetricSum<dim, VectorizedArray<Number>, n_rows_1d>
  : public TensorProductMatrixSymmetricSumBase<dim,
                                               VectorizedArray<Number>,
                                               n_rows_1d>
{
public:
  /**
   * 默认构造函数。
   *
   */
  TensorProductMatrixSymmetricSum() = default;

  /**
   * 等同于空构造函数的构造函数，并立即调用
   * reinit(constitution  std::array<Table<2,VectorizedArray<Number>  >,
   * dim>&,const  std::array<Table<2,VectorizedArray<Number>  >, dim>&) 。
   *
   */
  TensorProductMatrixSymmetricSum(
    const std::array<Table<2, VectorizedArray<Number>>, dim> &mass_matrix,
    const std::array<Table<2, VectorizedArray<Number>>, dim>
      &derivative_matrix);

  /**
   * 相当于空构造函数并立即调用
   * reinit(Table<2,VectorizedArray<Number> >&,const
   * Table<2,VectorizedArray<Number> >&) 的构造函数。
   *
   */
  TensorProductMatrixSymmetricSum(
    const Table<2, VectorizedArray<Number>> &mass_matrix,
    const Table<2, VectorizedArray<Number>> &derivative_matrix);

  /**
   * 通过将一维质量矩阵 @p mass_matrix 和一维导数矩阵 @p
   * derivative_matrix
   * 的数组分别复制到它的基类对应物中，并通过在
   * TensorProductMatrixSymmetricSumBase::eigenvalues 和
   * TensorProductMatrixSymmetricSumBase::eigenvectors,
   * 中分别组装关于广义特征值和特征向量来初始化张量积矩阵。
   * 请注意，目前的实现要求每个 $M_{d}$
   * 都是对称的和正定的，每个 $A_{d}$
   * 都是对称的和可逆的，但不一定是正定的。
   *
   */
  void
  reinit(const std::array<Table<2, VectorizedArray<Number>>, dim> &mass_matrix,
         const std::array<Table<2, VectorizedArray<Number>>, dim>
           &derivative_matrix);

  /**
   * 这个函数等同于之前的reinit()，只是我们为每个张量方向考虑相同的一维质量矩阵
   * @p mass_matrix 和相同的一维导数矩阵 @p derivative_matrix 。
   *
   */
  void
  reinit(const Table<2, VectorizedArray<Number>> &mass_matrix,
         const Table<2, VectorizedArray<Number>> &derivative_matrix);

private:
  /**
   * 一个基于完美转发的所有reinit()函数的通用实现，允许传递lvalue以及rvalue参数。
   * @tparam  MatrixArray 必须可以转换为
   * TensorProductMatrixSymmetricSumBase::mass_matrix 和
   * TensorProductMatrixSymmetricSumBase::derivative_matrix.
   * 的基本类型。
   *
   */
  template <typename MatrixArray>
  void
  reinit_impl(MatrixArray &&mass_matrix, MatrixArray &&derivative_matrix);
};


 /*----------------------- Inline functions ----------------------------------*/ 

#ifndef DOXYGEN

namespace internal
{
  namespace TensorProductMatrix
  {
    /**
     * 计算实际广义对称特征问题的广义特征值和特征向量
     * $A v = \lambda M v$
     * 。由于我们是在纯指针上操作，我们需要事先知道矩阵的大小。注意，特征值和特征向量的数据数组也必须被初始化为适当的大小。(不可能检查数组的界限)
     *
     */
    template <typename Number>
    void
    spectral_assembly(const Number *     mass_matrix,
                      const Number *     derivative_matrix,
                      const unsigned int n_rows,
                      const unsigned int n_cols,
                      Number *           eigenvalues,
                      Number *           eigenvectors)
    {
      Assert(n_rows == n_cols, ExcNotImplemented());

      auto &&transpose_fill_nm = [](Number *           out,
                                    const Number *     in,
                                    const unsigned int n,
                                    const unsigned int m) {
        for (unsigned int mm = 0; mm < m; ++mm)
          for (unsigned int nn = 0; nn < n; ++nn)
            out[mm + nn * m] = *(in++);
      };

      std::vector<dealii::Vector<Number>> eigenvecs(n_rows);
      LAPACKFullMatrix<Number>            mass_copy(n_rows, n_cols);
      LAPACKFullMatrix<Number>            deriv_copy(n_rows, n_cols);

      transpose_fill_nm(&(mass_copy(0, 0)), mass_matrix, n_rows, n_cols);
      transpose_fill_nm(&(deriv_copy(0, 0)), derivative_matrix, n_rows, n_cols);

      deriv_copy.compute_generalized_eigenvalues_symmetric(mass_copy,
                                                           eigenvecs);
      AssertDimension(eigenvecs.size(), n_rows);
      for (unsigned int i = 0; i < n_rows; ++i)
        for (unsigned int j = 0; j < n_cols; ++j, ++eigenvectors)
          *eigenvectors = eigenvecs[j][i];

      for (unsigned int i = 0; i < n_rows; ++i, ++eigenvalues)
        *eigenvalues = deriv_copy.eigenvalue(i).real();
    }
  } // namespace TensorProductMatrix
} // namespace internal


template <int dim, typename Number, int n_rows_1d>
inline unsigned int
TensorProductMatrixSymmetricSumBase<dim, Number, n_rows_1d>::m() const
{
  unsigned int m = mass_matrix[0].n_rows();
  for (unsigned int d = 1; d < dim; ++d)
    m *= mass_matrix[d].n_rows();
  return m;
}



template <int dim, typename Number, int n_rows_1d>
inline unsigned int
TensorProductMatrixSymmetricSumBase<dim, Number, n_rows_1d>::n() const
{
  unsigned int n = mass_matrix[0].n_cols();
  for (unsigned int d = 1; d < dim; ++d)
    n *= mass_matrix[d].n_cols();
  return n;
}



template <int dim, typename Number, int n_rows_1d>
inline void
TensorProductMatrixSymmetricSumBase<dim, Number, n_rows_1d>::vmult(
  const ArrayView<Number> &      dst_view,
  const ArrayView<const Number> &src_view) const
{
  AssertDimension(dst_view.size(), this->m());
  AssertDimension(src_view.size(), this->n());
  std::lock_guard<std::mutex> lock(this->mutex);
  const unsigned int          n = Utilities::fixed_power<dim>(
    n_rows_1d > 0 ? n_rows_1d : eigenvalues[0].size());
  tmp_array.resize_fast(n * 2);
  constexpr int kernel_size = n_rows_1d > 0 ? n_rows_1d : 0;
  internal::EvaluatorTensorProduct<internal::evaluate_general,
                                   dim,
                                   kernel_size,
                                   kernel_size,
                                   Number>
                eval(AlignedVector<Number>{},
         AlignedVector<Number>{},
         AlignedVector<Number>{},
         mass_matrix[0].n_rows(),
         mass_matrix[0].n_rows());
  Number *      t   = tmp_array.begin();
  const Number *src = src_view.begin();
  Number *      dst = dst_view.data();

  if (dim == 1)
    {
      const Number *A = &derivative_matrix[0](0, 0);
      eval.template apply<0, false, false>(A, src, dst);
    }

  else if (dim == 2)
    {
      const Number *A0 = &derivative_matrix[0](0, 0);
      const Number *M0 = &mass_matrix[0](0, 0);
      const Number *A1 = &derivative_matrix[1](0, 0);
      const Number *M1 = &mass_matrix[1](0, 0);
      eval.template apply<0, false, false>(M0, src, t);
      eval.template apply<1, false, false>(A1, t, dst);
      eval.template apply<0, false, false>(A0, src, t);
      eval.template apply<1, false, true>(M1, t, dst);
    }

  else if (dim == 3)
    {
      const Number *A0 = &derivative_matrix[0](0, 0);
      const Number *M0 = &mass_matrix[0](0, 0);
      const Number *A1 = &derivative_matrix[1](0, 0);
      const Number *M1 = &mass_matrix[1](0, 0);
      const Number *A2 = &derivative_matrix[2](0, 0);
      const Number *M2 = &mass_matrix[2](0, 0);
      eval.template apply<0, false, false>(M0, src, t + n);
      eval.template apply<1, false, false>(M1, t + n, t);
      eval.template apply<2, false, false>(A2, t, dst);
      eval.template apply<1, false, false>(A1, t + n, t);
      eval.template apply<0, false, false>(A0, src, t + n);
      eval.template apply<1, false, true>(M1, t + n, t);
      eval.template apply<2, false, true>(M2, t, dst);
    }

  else
    AssertThrow(false, ExcNotImplemented());
}



template <int dim, typename Number, int n_rows_1d>
inline void
TensorProductMatrixSymmetricSumBase<dim, Number, n_rows_1d>::apply_inverse(
  const ArrayView<Number> &      dst_view,
  const ArrayView<const Number> &src_view) const
{
  AssertDimension(dst_view.size(), this->n());
  AssertDimension(src_view.size(), this->m());
  std::lock_guard<std::mutex> lock(this->mutex);
  const unsigned int n = n_rows_1d > 0 ? n_rows_1d : eigenvalues[0].size();
  tmp_array.resize_fast(Utilities::fixed_power<dim>(n));
  constexpr int kernel_size = n_rows_1d > 0 ? n_rows_1d : 0;
  internal::EvaluatorTensorProduct<internal::evaluate_general,
                                   dim,
                                   kernel_size,
                                   kernel_size,
                                   Number>
                eval(AlignedVector<Number>(),
         AlignedVector<Number>(),
         AlignedVector<Number>(),
         mass_matrix[0].n_rows(),
         mass_matrix[0].n_rows());
  Number *      t   = tmp_array.begin();
  const Number *src = src_view.data();
  Number *      dst = dst_view.data();

  // NOTE: dof_to_quad has to be interpreted as 'dof to eigenvalue index'
  //       --> apply<.,true,.> (S,src,dst) calculates dst = S^T * src,
  //       --> apply<.,false,.> (S,src,dst) calculates dst = S * src,
  //       while the eigenvectors are stored column-wise in S, i.e.
  //       rows correspond to dofs whereas columns to eigenvalue indices!
  if (dim == 1)
    {
      const Number *S = &eigenvectors[0](0, 0);
      eval.template apply<0, true, false>(S, src, t);
      for (unsigned int i = 0; i < n; ++i)
        t[i] /= eigenvalues[0][i];
      eval.template apply<0, false, false>(S, t, dst);
    }

  else if (dim == 2)
    {
      const Number *S0 = &(eigenvectors[0](0, 0));
      const Number *S1 = &(eigenvectors[1](0, 0));
      eval.template apply<0, true, false>(S0, src, t);
      eval.template apply<1, true, false>(S1, t, dst);
      for (unsigned int i1 = 0, c = 0; i1 < n; ++i1)
        for (unsigned int i0 = 0; i0 < n; ++i0, ++c)
          dst[c] /= (eigenvalues[1][i1] + eigenvalues[0][i0]);
      eval.template apply<0, false, false>(S0, dst, t);
      eval.template apply<1, false, false>(S1, t, dst);
    }

  else if (dim == 3)
    {
      const Number *S0 = &eigenvectors[0](0, 0);
      const Number *S1 = &eigenvectors[1](0, 0);
      const Number *S2 = &eigenvectors[2](0, 0);
      eval.template apply<0, true, false>(S0, src, t);
      eval.template apply<1, true, false>(S1, t, dst);
      eval.template apply<2, true, false>(S2, dst, t);
      for (unsigned int i2 = 0, c = 0; i2 < n; ++i2)
        for (unsigned int i1 = 0; i1 < n; ++i1)
          for (unsigned int i0 = 0; i0 < n; ++i0, ++c)
            t[c] /=
              (eigenvalues[2][i2] + eigenvalues[1][i1] + eigenvalues[0][i0]);
      eval.template apply<0, false, false>(S0, t, dst);
      eval.template apply<1, false, false>(S1, dst, t);
      eval.template apply<2, false, false>(S2, t, dst);
    }

  else
    Assert(false, ExcNotImplemented());
}


//---------------------- TensorProductMatrixSymmetricSum ----------------------

template <int dim, typename Number, int n_rows_1d>
inline TensorProductMatrixSymmetricSum<dim, Number, n_rows_1d>::
  TensorProductMatrixSymmetricSum(
    const std::array<Table<2, Number>, dim> &mass_matrix,
    const std::array<Table<2, Number>, dim> &derivative_matrix)
{
  reinit(mass_matrix, derivative_matrix);
}



template <int dim, typename Number, int n_rows_1d>
inline TensorProductMatrixSymmetricSum<dim, Number, n_rows_1d>::
  TensorProductMatrixSymmetricSum(
    const std::array<FullMatrix<Number>, dim> &mass_matrix,
    const std::array<FullMatrix<Number>, dim> &derivative_matrix)
{
  reinit(mass_matrix, derivative_matrix);
}



template <int dim, typename Number, int n_rows_1d>
inline TensorProductMatrixSymmetricSum<dim, Number, n_rows_1d>::
  TensorProductMatrixSymmetricSum(const Table<2, Number> &mass_matrix,
                                  const Table<2, Number> &derivative_matrix)
{
  reinit(mass_matrix, derivative_matrix);
}



template <int dim, typename Number, int n_rows_1d>
template <typename MatrixArray>
inline void
TensorProductMatrixSymmetricSum<dim, Number, n_rows_1d>::reinit_impl(
  MatrixArray &&mass_matrices_,
  MatrixArray &&derivative_matrices_)
{
  auto &&mass_matrices       = std::forward<MatrixArray>(mass_matrices_);
  auto &&derivative_matrices = std::forward<MatrixArray>(derivative_matrices_);
  this->mass_matrix          = mass_matrices;
  this->derivative_matrix    = derivative_matrices;

  for (int dir = 0; dir < dim; ++dir)
    {
      Assert(n_rows_1d == -1 ||
               (n_rows_1d > 0 && static_cast<unsigned int>(n_rows_1d) ==
                                   mass_matrices[dir].n_rows()),
             ExcDimensionMismatch(n_rows_1d, mass_matrices[dir].n_rows()));
      AssertDimension(mass_matrices[dir].n_rows(), mass_matrices[dir].n_cols());
      AssertDimension(mass_matrices[dir].n_rows(),
                      derivative_matrices[dir].n_rows());
      AssertDimension(mass_matrices[dir].n_rows(),
                      derivative_matrices[dir].n_cols());

      this->eigenvectors[dir].reinit(mass_matrices[dir].n_cols(),
                                     mass_matrices[dir].n_rows());
      this->eigenvalues[dir].resize(mass_matrices[dir].n_cols());
      internal::TensorProductMatrix::spectral_assembly<Number>(
        &(mass_matrices[dir](0, 0)),
        &(derivative_matrices[dir](0, 0)),
        mass_matrices[dir].n_rows(),
        mass_matrices[dir].n_cols(),
        this->eigenvalues[dir].begin(),
        &(this->eigenvectors[dir](0, 0)));
    }
}



template <int dim, typename Number, int n_rows_1d>
inline void
TensorProductMatrixSymmetricSum<dim, Number, n_rows_1d>::reinit(
  const std::array<Table<2, Number>, dim> &mass_matrix,
  const std::array<Table<2, Number>, dim> &derivative_matrix)
{
  reinit_impl(mass_matrix, derivative_matrix);
}



template <int dim, typename Number, int n_rows_1d>
inline void
TensorProductMatrixSymmetricSum<dim, Number, n_rows_1d>::reinit(
  const std::array<FullMatrix<Number>, dim> &mass_matrix,
  const std::array<FullMatrix<Number>, dim> &derivative_matrix)
{
  std::array<Table<2, Number>, dim> mass_copy;
  std::array<Table<2, Number>, dim> deriv_copy;

  std::transform(mass_matrix.cbegin(),
                 mass_matrix.cend(),
                 mass_copy.begin(),
                 [](const FullMatrix<Number> &m) -> Table<2, Number> {
                   return m;
                 });
  std::transform(derivative_matrix.cbegin(),
                 derivative_matrix.cend(),
                 deriv_copy.begin(),
                 [](const FullMatrix<Number> &m) -> Table<2, Number> {
                   return m;
                 });

  reinit_impl(std::move(mass_copy), std::move(deriv_copy));
}



template <int dim, typename Number, int n_rows_1d>
inline void
TensorProductMatrixSymmetricSum<dim, Number, n_rows_1d>::reinit(
  const Table<2, Number> &mass_matrix,
  const Table<2, Number> &derivative_matrix)
{
  std::array<Table<2, Number>, dim> mass_matrices;
  std::array<Table<2, Number>, dim> derivative_matrices;

  std::fill(mass_matrices.begin(), mass_matrices.end(), mass_matrix);
  std::fill(derivative_matrices.begin(),
            derivative_matrices.end(),
            derivative_matrix);

  reinit_impl(std::move(mass_matrices), std::move(derivative_matrices));
}



//------------- vectorized spec.: TensorProductMatrixSymmetricSum -------------

template <int dim, typename Number, int n_rows_1d>
inline TensorProductMatrixSymmetricSum<dim,
                                       VectorizedArray<Number>,
                                       n_rows_1d>::
  TensorProductMatrixSymmetricSum(
    const std::array<Table<2, VectorizedArray<Number>>, dim> &mass_matrix,
    const std::array<Table<2, VectorizedArray<Number>>, dim> &derivative_matrix)
{
  reinit(mass_matrix, derivative_matrix);
}



template <int dim, typename Number, int n_rows_1d>
inline TensorProductMatrixSymmetricSum<dim,
                                       VectorizedArray<Number>,
                                       n_rows_1d>::
  TensorProductMatrixSymmetricSum(
    const Table<2, VectorizedArray<Number>> &mass_matrix,
    const Table<2, VectorizedArray<Number>> &derivative_matrix)
{
  reinit(mass_matrix, derivative_matrix);
}



template <int dim, typename Number, int n_rows_1d>
template <typename MatrixArray>
inline void
TensorProductMatrixSymmetricSum<dim, VectorizedArray<Number>, n_rows_1d>::
  reinit_impl(MatrixArray &&mass_matrices_, MatrixArray &&derivative_matrices_)
{
  auto &&mass_matrix       = std::forward<MatrixArray>(mass_matrices_);
  auto &&derivative_matrix = std::forward<MatrixArray>(derivative_matrices_);
  this->mass_matrix        = mass_matrix;
  this->derivative_matrix  = derivative_matrix;

  constexpr unsigned int macro_size = VectorizedArray<Number>::size();
  std::size_t            n_rows_max = (n_rows_1d > 0) ? n_rows_1d : 0;
  if (n_rows_1d == -1)
    for (unsigned int d = 0; d < dim; ++d)
      n_rows_max = std::max(n_rows_max, mass_matrix[d].n_rows());
  const std::size_t nm_flat_size_max = n_rows_max * n_rows_max * macro_size;
  const std::size_t n_flat_size_max  = n_rows_max * macro_size;

  std::vector<Number> mass_matrix_flat;
  std::vector<Number> deriv_matrix_flat;
  std::vector<Number> eigenvalues_flat;
  std::vector<Number> eigenvectors_flat;
  mass_matrix_flat.resize(nm_flat_size_max);
  deriv_matrix_flat.resize(nm_flat_size_max);
  eigenvalues_flat.resize(n_flat_size_max);
  eigenvectors_flat.resize(nm_flat_size_max);
  std::array<unsigned int, macro_size> offsets_nm;
  std::array<unsigned int, macro_size> offsets_n;
  for (int dir = 0; dir < dim; ++dir)
    {
      Assert(n_rows_1d == -1 ||
               (n_rows_1d > 0 && static_cast<unsigned int>(n_rows_1d) ==
                                   mass_matrix[dir].n_rows()),
             ExcDimensionMismatch(n_rows_1d, mass_matrix[dir].n_rows()));
      AssertDimension(mass_matrix[dir].n_rows(), mass_matrix[dir].n_cols());
      AssertDimension(mass_matrix[dir].n_rows(),
                      derivative_matrix[dir].n_rows());
      AssertDimension(mass_matrix[dir].n_rows(),
                      derivative_matrix[dir].n_cols());

      const unsigned int n_rows = mass_matrix[dir].n_rows();
      const unsigned int n_cols = mass_matrix[dir].n_cols();
      const unsigned int nm     = n_rows * n_cols;
      for (unsigned int vv = 0; vv < macro_size; ++vv)
        offsets_nm[vv] = nm * vv;

      vectorized_transpose_and_store(false,
                                     nm,
                                     &(mass_matrix[dir](0, 0)),
                                     offsets_nm.cbegin(),
                                     mass_matrix_flat.data());
      vectorized_transpose_and_store(false,
                                     nm,
                                     &(derivative_matrix[dir](0, 0)),
                                     offsets_nm.cbegin(),
                                     deriv_matrix_flat.data());

      const Number *mass_cbegin    = mass_matrix_flat.data();
      const Number *deriv_cbegin   = deriv_matrix_flat.data();
      Number *      eigenvec_begin = eigenvectors_flat.data();
      Number *      eigenval_begin = eigenvalues_flat.data();
      for (unsigned int lane = 0; lane < macro_size; ++lane)
        internal::TensorProductMatrix::spectral_assembly<Number>(
          mass_cbegin + nm * lane,
          deriv_cbegin + nm * lane,
          n_rows,
          n_cols,
          eigenval_begin + n_rows * lane,
          eigenvec_begin + nm * lane);

      this->eigenvalues[dir].resize(n_rows);
      this->eigenvectors[dir].reinit(n_rows, n_cols);
      for (unsigned int vv = 0; vv < macro_size; ++vv)
        offsets_n[vv] = n_rows * vv;
      vectorized_load_and_transpose(n_rows,
                                    eigenvalues_flat.data(),
                                    offsets_n.cbegin(),
                                    this->eigenvalues[dir].begin());
      vectorized_load_and_transpose(nm,
                                    eigenvectors_flat.data(),
                                    offsets_nm.cbegin(),
                                    &(this->eigenvectors[dir](0, 0)));
    }
}



template <int dim, typename Number, int n_rows_1d>
inline void
TensorProductMatrixSymmetricSum<dim, VectorizedArray<Number>, n_rows_1d>::
  reinit(
    const std::array<Table<2, VectorizedArray<Number>>, dim> &mass_matrix,
    const std::array<Table<2, VectorizedArray<Number>>, dim> &derivative_matrix)
{
  reinit_impl(mass_matrix, derivative_matrix);
}



template <int dim, typename Number, int n_rows_1d>
inline void
TensorProductMatrixSymmetricSum<dim, VectorizedArray<Number>, n_rows_1d>::
  reinit(const Table<2, VectorizedArray<Number>> &mass_matrix,
         const Table<2, VectorizedArray<Number>> &derivative_matrix)
{
  std::array<Table<2, VectorizedArray<Number>>, dim> mass_matrices;
  std::array<Table<2, VectorizedArray<Number>>, dim> derivative_matrices;

  std::fill(mass_matrices.begin(), mass_matrices.end(), mass_matrix);
  std::fill(derivative_matrices.begin(),
            derivative_matrices.end(),
            derivative_matrix);

  reinit_impl(std::move(mass_matrices), std::move(derivative_matrices));
}



#endif

DEAL_II_NAMESPACE_CLOSE

#endif


