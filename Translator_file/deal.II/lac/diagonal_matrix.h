//include/deal.II-translator/lac/diagonal_matrix_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2020 by the deal.II authors
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

#ifndef dealii_diagonal_matrix_h
#define dealii_diagonal_matrix_h


#include <deal.II/base/config.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/vector_operation.h>

DEAL_II_NAMESPACE_OPEN

/**
 * 这个类表示一个基于大小为<i>n</i>的向量的<i>n x
 * n</i>对角矩阵。矩阵-向量的乘积由 @p  VectorType::scale,
 * 实现，所以模板向量类需要提供一个 @p scale() 方法。
 * 当把这个类和 ConstraintsMatrix::distribute_local_to_global(),
 * 一起使用时，底层向量需要提供对装配过程中单元格引用的所有条目的写入权限。这意味着该类也需要访问除调用进程外的其他进程所拥有的幽灵条目。在实践中，这需要对向量进行如下初始化
 *
 * @code
 * DiagonalMatrix<LinearAlgebra::distributed::Vector<double> > diagonal_matrix;
 * LinearAlgebra::distributed::Vector<double> &diagonal_vector =
 * diagonal_matrix.get_vector();
 * diagonal_vector.reinit(locally_owned_dofs,
 *                      locally_relevant_dofs,
 *                      mpi_communicator);
 * @endcode
 *
 *
 *
 */
template <typename VectorType = Vector<double>>
class DiagonalMatrix : public Subscriptor
{
public:
  using value_type = typename VectorType::value_type;
  using size_type  = typename VectorType::size_type;

  /**
   * 默认构造函数。该对象仍然需要重新初始化才能使用。
   *
   */
  DiagonalMatrix() = default;

  /**
   * 构造函数将此对象初始化为大小为`n x
   * n`的对角线矩阵，其中`n`是矢量的大小，其对角线条目等于
   * @p vec. 的元素。
   *
   */
  explicit DiagonalMatrix(const VectorType &vec);

  /**
   * 通过复制向量 @p vec.
   * 的内容，用一个给定的向量进行初始化
   *
   */
  void
  reinit(const VectorType &vec);

  /**
   * 压缩数据结构，并允许产生的矩阵用于所有其他操作，如矩阵-向量乘积。这是一个集体操作，即在并行使用时需要在所有处理器上运行。
   *
   */
  void
  compress(VectorOperation::values operation);

  /**
   * 返回一个对底层向量的引用，以便对矩阵对角线上的条目进行操作。
   *
   */
  VectorType &
  get_vector();

  /**
   * 清除此对象的内容并重置为默认构造函数的状态。
   *
   */
  void
  clear();

  /**
   * 返回一个对底层向量的只读引用。
   *
   */
  const VectorType &
  get_vector() const;

  /**
   * 这个矩阵的行数。这个数字对应于底层向量的大小。
   *
   */
  size_type
  m() const;

  /**
   * 此矩阵的列数。这个数字与基础向量的大小相对应。
   *
   */
  size_type
  n() const;

  /**
   * 对一个值的只读访问。由于矩阵存储的原因，这被限制在<i>i==j</i>的情况下。
   * 如果代表对角线的向量是用MPI分布的，实际上可能不是所有的指数<i>i</i>都能被访问。请参考方法
   * <code>get_vector().locally_owned_elements()</code>
   * ，了解实际可访问的条目。
   *
   */
  value_type
  operator()(const size_type i, const size_type j) const;

  /**
   * 对一个值的读写访问。由于矩阵存储的原因，这被限制在<i>i==j</i>的情况下。
   * 如果代表对角线的向量是用MPI分布的，那么实际上可能不是所有的指数<i>i</i>都可以访问。请参考方法
   * <code>get_vector().locally_owned_elements()</code>
   * ，了解实际可访问的条目。
   *
   */
  value_type &
  operator()(const size_type i, const size_type j);

  /**
   * 在给定的全局矩阵行中，在col_indices指定的列中添加一个由<tt>values</tt>给出的数值数组。由于该矩阵的存储方式，条目只被添加到矩阵的对角线上。所有其他的条目都被忽略，并且不抛出异常。
   * 这个函数是为了与deal.II中的其他矩阵类保持一致的接口，可以在
   * AffineConstraints::distribute_local_to_global
   * 中使用，得到与组装成稀疏矩阵时完全相同的对角线。
   *
   */
  template <typename number2>
  void
  add(const size_type  row,
      const size_type  n_cols,
      const size_type *col_indices,
      const number2 *  values,
      const bool       elide_zero_values      = true,
      const bool       col_indices_are_sorted = false);

  /**
   * 为元素(i,j)加值。
   * 由于这个矩阵的存储方式，条目只被添加到矩阵的对角线上。所有其他的条目都被忽略，并且不抛出异常。
   *
   */
  void
  add(const size_type i, const size_type j, const value_type value);

  /**
   * 与给定的矩阵进行矩阵-向量乘法。
   *
   */
  void
  vmult(VectorType &dst, const VectorType &src) const;

  /**
   * 对给定的矩阵进行转置矩阵-向量乘法。因为这代表了一个对角线矩阵，与vmult()完全相同。
   *
   */
  void
  Tvmult(VectorType &dst, const VectorType &src) const;

  /**
   * 将矩阵-向量乘法的结果添加到目标向量dst中。需要创建一个临时向量，这使得性能比
   * @p vmult(). 的要慢。
   *
   */
  void
  vmult_add(VectorType &dst, const VectorType &src) const;

  /**
   * 将转置矩阵向量乘法的结果添加到目标向量dst中。需要创建一个临时向量，这使得性能比
   * @p Tvmult(). 慢。
   *
   */
  void
  Tvmult_add(VectorType &dst, const VectorType &src) const;

  /**
   * 初始化向量 @p dst ，使其具有与本类的 @p diagonal
   * 成员相同的大小和分区。
   * 这是linear_operator()所要求的接口的一部分。
   *
   */
  void
  initialize_dof_vector(VectorType &dst) const;

  /**
   * 返回这个对象的内存消耗。
   *
   */
  std::size_t
  memory_consumption() const;

private:
  /**
   * 存储的向量。
   *
   */
  VectorType diagonal;
};

 /* ---------------------------------- Inline functions ------------------- */ 

#ifndef DOXYGEN

template <typename VectorType>
DiagonalMatrix<VectorType>::DiagonalMatrix(const VectorType &vec)
  : diagonal(vec)
{}



template <typename VectorType>
void
DiagonalMatrix<VectorType>::clear()
{
  diagonal.reinit(0);
}



template <typename VectorType>
std::size_t
DiagonalMatrix<VectorType>::memory_consumption() const
{
  return diagonal.memory_consumption();
}



template <typename VectorType>
void
DiagonalMatrix<VectorType>::reinit(const VectorType &vec)
{
  diagonal = vec;
}



template <typename VectorType>
void
DiagonalMatrix<VectorType>::initialize_dof_vector(VectorType &dst) const
{
  dst.reinit(diagonal);
}



template <typename VectorType>
void
DiagonalMatrix<VectorType>::compress(VectorOperation::values operation)
{
  diagonal.compress(operation);
}



template <typename VectorType>
VectorType &
DiagonalMatrix<VectorType>::get_vector()
{
  return diagonal;
}



template <typename VectorType>
const VectorType &
DiagonalMatrix<VectorType>::get_vector() const
{
  return diagonal;
}



template <typename VectorType>
typename VectorType::size_type
DiagonalMatrix<VectorType>::m() const
{
  return diagonal.size();
}



template <typename VectorType>
typename VectorType::size_type
DiagonalMatrix<VectorType>::n() const
{
  return diagonal.size();
}



template <typename VectorType>
typename VectorType::value_type
DiagonalMatrix<VectorType>::operator()(const size_type i,
                                       const size_type j) const
{
  Assert(i == j, ExcIndexRange(j, i, i + 1));
  (void)j;
  return diagonal(i);
}



template <typename VectorType>
typename VectorType::value_type &
DiagonalMatrix<VectorType>::operator()(const size_type i, const size_type j)
{
  Assert(i == j, ExcIndexRange(j, i, i + 1));
  (void)j;
  return diagonal(i);
}



template <typename VectorType>
template <typename number2>
void
DiagonalMatrix<VectorType>::add(const size_type  row,
                                const size_type  n_cols,
                                const size_type *col_indices,
                                const number2 *  values,
                                const bool,
                                const bool)
{
  for (size_type i = 0; i < n_cols; ++i)
    if (col_indices[i] == row)
      diagonal(row) += values[i];
}



template <typename VectorType>
void
DiagonalMatrix<VectorType>::add(const size_type  i,
                                const size_type  j,
                                const value_type value)
{
  if (i == j)
    diagonal(i) += value;
}



template <typename VectorType>
void
DiagonalMatrix<VectorType>::vmult(VectorType &dst, const VectorType &src) const
{
  dst = src;
  dst.scale(diagonal);
}



template <typename VectorType>
void
DiagonalMatrix<VectorType>::Tvmult(VectorType &dst, const VectorType &src) const
{
  vmult(dst, src);
}



template <typename VectorType>
void
DiagonalMatrix<VectorType>::vmult_add(VectorType &      dst,
                                      const VectorType &src) const
{
  VectorType tmp(src);
  tmp.scale(diagonal);
  dst += tmp;
}



template <typename VectorType>
void
DiagonalMatrix<VectorType>::Tvmult_add(VectorType &      dst,
                                       const VectorType &src) const
{
  vmult_add(dst, src);
}


#endif

DEAL_II_NAMESPACE_CLOSE

#endif


