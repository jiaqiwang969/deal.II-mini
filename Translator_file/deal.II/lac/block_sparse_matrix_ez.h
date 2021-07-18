//include/deal.II-translator/lac/block_sparse_matrix_ez_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2002 - 2020 by the deal.II authors
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

#ifndef dealii_block_sparse_matrix_ez_h
#define dealii_block_sparse_matrix_ez_h


// TODO: Derive BlockSparseMatrixEZ from BlockMatrixBase, like all the
// other block matrices as well; this would allow to instantiate a few
// functions with this template argument as well (in particular
// AffineConstraints::distribute_local_to_global)

#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/base/table.h>

#include <deal.II/lac/block_indices.h>
#include <deal.II/lac/sparse_matrix_ez.h>

DEAL_II_NAMESPACE_OPEN

// Forward declaration
#ifndef DOXYGEN
template <typename Number>
class BlockVector;
#endif

/*!   @addtogroup  Matrix1  
     * @{  

* 
*
*/


/**
 * 一个由SparseMatrixEZ类型的块组成的块矩阵。
 * 与其他块对象一样，这个矩阵可以像SparseMatrixEZ一样使用，当涉及到对条目的访问时。然后，还有一些函数用于与BlockVector相乘，以及对单个块的访问。
 * @see   @ref GlossBlockLA  "块（线性代数）"
 *
 *
 */
template <typename Number>
class BlockSparseMatrixEZ : public Subscriptor
{
public:
  /**
   * 声明容器大小的类型。
   *
   */
  using size_type = types::global_dof_index;

  /**
   * 默认构造函数。其结果是一个尺寸为零的空对象。
   *
   */
  BlockSparseMatrixEZ() = default;

  /**
   * 构造函数以给定的块的行和列的数量设置一个对象。块本身仍有零维度。
   *
   */
  BlockSparseMatrixEZ(const unsigned int block_rows,
                      const unsigned int block_cols);

  /**
   * 复制构造函数。这对于一些容器类来说是需要的。它创建一个具有相同数量的块行和块列的对象。因为它调用了SparseMatrixEZ的复制构造函数，所以块s必须是空的。
   *
   */
  BlockSparseMatrixEZ(const BlockSparseMatrixEZ<Number> &);

  /**
   * 拷贝操作符。和复制构造函数一样，它只能被调用于空块的对象。
   *
   */
  BlockSparseMatrixEZ &
  operator=(const BlockSparseMatrixEZ<Number> &);

  /**
   * 这个操作符将一个标量分配给一个矩阵。因为这通常没有什么意义（我们应该把所有的矩阵条目都设置为这个值吗？仅仅是稀疏模式的非零条目？），这个操作只允许在实际要分配的值为零时进行。这个操作符的存在只是为了允许明显的符号<tt>matrix=0</tt>，它将矩阵的所有元素设置为零，但保留之前使用的稀疏模式。
   *
   */
  BlockSparseMatrixEZ &
  operator=(const double d);


  /**
   * 将矩阵设置为零维，并释放内存。
   *
   */
  void
  clear();

  /**
   * 初始化为给定的块数。
   * 在此操作之后，矩阵将具有所提供的块的尺寸。每个区块的尺寸都是零，必须随后进行初始化。在设置块的尺寸后，必须调用collect_sizes()来更新内部数据结构。
   *
   */
  void
  reinit(const unsigned int n_block_rows, const unsigned int n_block_cols);
  /**
   * 这个函数收集了子对象的尺寸，并将其存储在内部数组中，以便能够将矩阵的全局索引转为子对象的索引。在你改变了子对象的大小之后，你必须*每次都调用这个函数。
   *
   */
  void
  collect_sizes();

  /**
   * 访问具有给定坐标的块。
   *
   */
  SparseMatrixEZ<Number> &
  block(const unsigned int row, const unsigned int column);


  /**
   * 访问具有给定坐标的块。对于常量对象的版本。
   *
   */
  const SparseMatrixEZ<Number> &
  block(const unsigned int row, const unsigned int column) const;

  /**
   * 返回一列中的块的数量。
   *
   */
  unsigned int
  n_block_rows() const;

  /**
   * 返回一行中的块数。
   *
   */
  unsigned int
  n_block_cols() const;

  /**
   * 返回该对象是否为空。如果没有分配内存，它就是空的，这与两个维度都是零是一样的。这个函数只是对所有子矩阵的各自调用的串联。
   *
   */
  bool
  empty() const;

  /**
   * 返回该矩阵的行数，相当于共域（或范围）空间的维数。它是这个矩阵的子矩阵块上的行数之和。回顾一下，该矩阵的大小为m()乘以n()。
   *
   */
  size_type
  m() const;

  /**
   * 返回该矩阵的列数，等于域空间的维度。它是该矩阵的子矩阵块的列数之和。回顾一下，该矩阵的大小为m()乘以n()。
   *
   */
  size_type
  n() const;

  /**
   * 将元素<tt>(i,j)</tt>设为 @p value.
   * 如果该条目不存在或<tt>值</tt>不是有限数，则抛出一个错误。尽管如此，它仍然允许在不存在的字段中存储零值。
   *
   */
  void
  set(const size_type i, const size_type j, const Number value);

  /**
   * 在元素<tt>(i,j)</tt>上添加 @p value 。
   * 如果该条目不存在或者<tt>value</tt>不是一个有限的数字，则抛出一个错误。尽管如此，它仍然允许在不存在的字段中存储零值。
   *
   */
  void
  add(const size_type i, const size_type j, const Number value);


  /**
   * 矩阵-向量乘法：让 $dst = M*src$ 与 $M$ 是这个矩阵。
   *
   */
  template <typename somenumber>
  void
  vmult(BlockVector<somenumber> &dst, const BlockVector<somenumber> &src) const;

  /**
   * 矩阵-向量乘法：让 $dst = M^T*src$ 与 $M$
   * 为这个矩阵。这个函数与vmult()的作用相同，但需要转置的矩阵。
   *
   */
  template <typename somenumber>
  void
  Tvmult(BlockVector<somenumber> &      dst,
         const BlockVector<somenumber> &src) const;

  /**
   * 加法 矩阵-向量乘法。在 $dst$ 上添加 $M*src$ ， $M$
   * 为该矩阵。
   *
   */
  template <typename somenumber>
  void
  vmult_add(BlockVector<somenumber> &      dst,
            const BlockVector<somenumber> &src) const;

  /**
   * 添加矩阵-向量乘法。将 $M^T*src$ 加到 $dst$ ， $M$
   * 是这个矩阵。这个函数与vmult_add()的作用相同，但需要转置的矩阵。
   *
   */
  template <typename somenumber>
  void
  Tvmult_add(BlockVector<somenumber> &      dst,
             const BlockVector<somenumber> &src) const;


  /**
   * 打印统计数据。如果 @p full 是 @p true,
   * ，则打印所有现有行长和分配行长的直方图。否则，只显示已分配和已使用条目的关系。
   *
   */
  template <class StreamType>
  void
  print_statistics(StreamType &s, bool full = false);

private:
  /**
   * 对象存储和管理行索引到子对象索引的转换。
   *
   */
  BlockIndices row_indices;

  /**
   * 储存和管理列索引到子对象索引的转换的对象。
   *
   */
  BlockIndices column_indices;

  /**
   * 实际的矩阵
   *
   */
  Table<2, SparseMatrixEZ<Number>> blocks;
};

 /*@}*/ 
 /*----------------------------------------------------------------------*/ 


template <typename Number>
inline unsigned int
BlockSparseMatrixEZ<Number>::n_block_rows() const
{
  return row_indices.size();
}



template <typename Number>
inline unsigned int
BlockSparseMatrixEZ<Number>::n_block_cols() const
{
  return column_indices.size();
}



template <typename Number>
inline SparseMatrixEZ<Number> &
BlockSparseMatrixEZ<Number>::block(const unsigned int row,
                                   const unsigned int column)
{
  AssertIndexRange(row, n_block_rows());
  AssertIndexRange(column, n_block_cols());

  return blocks[row][column];
}



template <typename Number>
inline const SparseMatrixEZ<Number> &
BlockSparseMatrixEZ<Number>::block(const unsigned int row,
                                   const unsigned int column) const
{
  AssertIndexRange(row, n_block_rows());
  AssertIndexRange(column, n_block_cols());

  return blocks[row][column];
}



template <typename Number>
inline typename BlockSparseMatrixEZ<Number>::size_type
BlockSparseMatrixEZ<Number>::m() const
{
  return row_indices.total_size();
}



template <typename Number>
inline typename BlockSparseMatrixEZ<Number>::size_type
BlockSparseMatrixEZ<Number>::n() const
{
  return column_indices.total_size();
}



template <typename Number>
inline void
BlockSparseMatrixEZ<Number>::set(const size_type i,
                                 const size_type j,
                                 const Number    value)
{
  AssertIsFinite(value);

  const std::pair<size_type, size_type> row_index =
                                          row_indices.global_to_local(i),
                                        col_index =
                                          column_indices.global_to_local(j);
  block(row_index.first, col_index.first)
    .set(row_index.second, col_index.second, value);
}



template <typename Number>
inline void
BlockSparseMatrixEZ<Number>::add(const size_type i,
                                 const size_type j,
                                 const Number    value)
{
  AssertIsFinite(value);

  const std::pair<unsigned int, size_type> row_index =
                                             row_indices.global_to_local(i),
                                           col_index =
                                             column_indices.global_to_local(j);
  block(row_index.first, col_index.first)
    .add(row_index.second, col_index.second, value);
}


template <typename Number>
template <typename somenumber>
void
BlockSparseMatrixEZ<Number>::vmult(BlockVector<somenumber> &      dst,
                                   const BlockVector<somenumber> &src) const
{
  Assert(dst.n_blocks() == n_block_rows(),
         ExcDimensionMismatch(dst.n_blocks(), n_block_rows()));
  Assert(src.n_blocks() == n_block_cols(),
         ExcDimensionMismatch(src.n_blocks(), n_block_cols()));

  dst = 0.;

  for (unsigned int row = 0; row < n_block_rows(); ++row)
    for (unsigned int col = 0; col < n_block_cols(); ++col)
      block(row, col).vmult_add(dst.block(row), src.block(col));
}



template <typename Number>
template <typename somenumber>
void
BlockSparseMatrixEZ<Number>::vmult_add(BlockVector<somenumber> &      dst,
                                       const BlockVector<somenumber> &src) const
{
  Assert(dst.n_blocks() == n_block_rows(),
         ExcDimensionMismatch(dst.n_blocks(), n_block_rows()));
  Assert(src.n_blocks() == n_block_cols(),
         ExcDimensionMismatch(src.n_blocks(), n_block_cols()));

  for (unsigned int row = 0; row < n_block_rows(); ++row)
    for (unsigned int col = 0; col < n_block_cols(); ++col)
      block(row, col).vmult_add(dst.block(row), src.block(col));
}



template <typename Number>
template <typename somenumber>
void
BlockSparseMatrixEZ<Number>::Tvmult(BlockVector<somenumber> &      dst,
                                    const BlockVector<somenumber> &src) const
{
  Assert(dst.n_blocks() == n_block_cols(),
         ExcDimensionMismatch(dst.n_blocks(), n_block_cols()));
  Assert(src.n_blocks() == n_block_rows(),
         ExcDimensionMismatch(src.n_blocks(), n_block_rows()));

  dst = 0.;

  for (unsigned int row = 0; row < n_block_rows(); ++row)
    for (unsigned int col = 0; col < n_block_cols(); ++col)
      block(row, col).Tvmult_add(dst.block(col), src.block(row));
}



template <typename Number>
template <typename somenumber>
void
BlockSparseMatrixEZ<Number>::Tvmult_add(
  BlockVector<somenumber> &      dst,
  const BlockVector<somenumber> &src) const
{
  Assert(dst.n_blocks() == n_block_cols(),
         ExcDimensionMismatch(dst.n_blocks(), n_block_cols()));
  Assert(src.n_blocks() == n_block_rows(),
         ExcDimensionMismatch(src.n_blocks(), n_block_rows()));

  for (unsigned int row = 0; row < n_block_rows(); ++row)
    for (unsigned int col = 0; col < n_block_cols(); ++col)
      block(row, col).Tvmult_add(dst.block(col), src.block(row));
}


template <typename number>
template <class StreamType>
inline void
BlockSparseMatrixEZ<number>::print_statistics(StreamType &out, bool full)
{
  size_type              used_total      = 0;
  size_type              allocated_total = 0;
  size_type              reserved_total  = 0;
  std::vector<size_type> used_by_line_total;

  size_type              used;
  size_type              allocated;
  size_type              reserved;
  std::vector<size_type> used_by_line;

  for (size_type i = 0; i < n_block_rows(); ++i)
    for (size_type j = 0; j < n_block_cols(); ++j)
      {
        used_by_line.clear();
        out << "block:\t" << i << '\t' << j << std::endl;
        block(i, j).compute_statistics(
          used, allocated, reserved, used_by_line, full);

        out << "used:" << used << std::endl
            << "allocated:" << allocated << std::endl
            << "reserved:" << reserved << std::endl;

        used_total += used;
        allocated_total += allocated;
        reserved_total += reserved;

        if (full)
          {
            used_by_line_total.resize(used_by_line.size());
            for (size_type i = 0; i < used_by_line.size(); ++i)
              if (used_by_line[i] != 0)
                {
                  out << "row-entries\t" << i << "\trows\t" << used_by_line[i]
                      << std::endl;
                  used_by_line_total[i] += used_by_line[i];
                }
          }
      }
  out << "Total" << std::endl
      << "used:" << used_total << std::endl
      << "allocated:" << allocated_total << std::endl
      << "reserved:" << reserved_total << std::endl;
  for (size_type i = 0; i < used_by_line_total.size(); ++i)
    if (used_by_line_total[i] != 0)
      {
        out << "row-entries\t" << i << "\trows\t" << used_by_line_total[i]
            << std::endl;
      }
}


DEAL_II_NAMESPACE_CLOSE

#endif // dealii_block_sparse_matrix_ez_h


