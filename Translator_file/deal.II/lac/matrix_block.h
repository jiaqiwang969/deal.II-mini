//include/deal.II-translator/lac/matrix_block_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2007 - 2020 by the deal.II authors
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

#ifndef dealii_matrix_block_h
#define dealii_matrix_block_h

#include <deal.II/base/config.h>

#include <deal.II/algorithms/any_data.h>

#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/mg_level_object.h>
#include <deal.II/base/smartpointer.h>

#include <deal.II/lac/block_indices.h>
#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>

#include <memory>


DEAL_II_NAMESPACE_OPEN

// Forward declarations
#ifndef DOXYGEN
template <typename MatrixType>
class MatrixBlock;
#endif

namespace internal
{
  template <typename MatrixType>
  void
  reinit(MatrixBlock<MatrixType> &v, const BlockSparsityPattern &p);

  template <typename number>
  void
  reinit(MatrixBlock<dealii::SparseMatrix<number>> &v,
         const BlockSparsityPattern &               p);
} // namespace internal

/**
 * 一个围绕矩阵对象的封装器，将坐标也存储在一个块状矩阵中。
 * 这个类是BlockMatrixBase的一个替代品，如果你只想生成系统的一个块，而不是整个系统。使用这个类的add()函数，可以使用用于块状矩阵的标准装配函数，但只输入其中一个块，仍然可以避免涉及的索引计算。
 * 这个类的原因是，在一个块系统中的不同块，我们可能需要不同数量的矩阵。例如，Oseen系统的预处理程序可以建立为一个块系统，其中压力块的形式为<b>M</b><sup>-1</sup><b>FA</b><sup>-1</sup>，<b>M</b>是压力质量矩阵，<b>A</b>是压力拉普拉斯，<b>F</b>是应用于压力空间的平流扩散算符。由于其他区块只需要一个矩阵，使用BlockSparseMatrix或类似的方法将是对内存的浪费。
 * 虽然add()函数使MatrixBlock看起来像一个用于组装的块状矩阵，但vmult()、Tvmult()、vmult_add()和Tvmult_add()函数使它的行为像一个MatrixType，当它应用于一个矢量时。这种行为允许我们在向量中存储MatrixBlock对象，例如在MGLevelObject中，而不用先提取#matrix。
 * MatrixBlock在使用BlockMatrixArray时很方便。一旦MatrixBlock被正确地初始化和填充，它就可以在最简单的情况下被用来作为。
 *
 * @code
 * MatrixBlockVector<SparseMatrix<double> > > blocks;
 *
 * ...
 *
 * BlockMatrixArray matrix (n_blocks, n_blocks);
 *
 * for (size_type i=0;i<blocks.size;++i)
 * matrix.enter(blocks.block(i).row, blocks.block(i).column,
 * blocks.matrix(i));
 * @endcode
 *
 * 在这里，我们的收获并不大，只是我们不需要在块系统中设置空块。
 *
 *
 * @note
 * 这个类期望，系统的行和列的BlockIndices对象是相等的。如果它们不相等，一些函数会抛出ExcNotImplemented。
 * @todo  压力舒尔补码的乘积预处理的例子。
 *
 *
 * @ingroup Matrix2
 *
 * @ingroup vector_valued
 * @see   @ref GlossBlockLA  "块（线性代数）"
 *
 */
template <typename MatrixType>
class MatrixBlock : public Subscriptor
{
public:
  /**
   * 申报容器尺寸的类型。
   *
   */
  using size_type = types::global_dof_index;

  /**
   * 为矩阵条目声明一个类型。
   *
   */
  using value_type = typename MatrixType::value_type;

  /**
   * 构造函数渲染一个未初始化的对象。
   *
   */
  MatrixBlock();

  /**
   * 复制构造函数。
   *
   */
  MatrixBlock(const MatrixBlock<MatrixType> &M) = default;

  /**
   * 赋值运算符。
   *
   */
  MatrixBlock<MatrixType> &
  operator=(const MatrixBlock<MatrixType> &) = default;

  /**
   * 构造函数设置块坐标，但不初始化矩阵。
   *
   */

  MatrixBlock(size_type i, size_type j);

  /**
   * 为一个新的BlockSparsityPattern重新初始化矩阵。这将调整#matrix以及#row_indices和#column_indices。
   * @note  稀疏模式的行和列块结构必须相等。
   *
   */
  void
  reinit(const BlockSparsityPattern &sparsity);

  operator MatrixType &();
  operator const MatrixType &() const;

  /**
   * 向元素添加<tt>value</tt>（<i>i,j</i>）。如果该条目不存在或在不同的块中，则抛出一个错误。
   *
   */
  void
  add(const size_type                       i,
      const size_type                       j,
      const typename MatrixType::value_type value);

  /**
   * 将FullMatrix中的所有元素添加到由<tt>indices</tt>给出的稀疏矩阵位置。这个函数假设一个二次元稀疏矩阵和一个二次元full_matrix。
   * 全局位置被转换为该块中的位置，如果全局索引没有指向#row和#column所指的块，则抛出ExcBlockIndexMismatch。
   * @todo  <tt>elide_zero_values</tt>目前被忽略。
   * 可选参数<tt>elide_zero_values</tt>可以用来指定是无论如何都要添加零值，还是要过滤掉这些数据，只添加非零值。默认值是<tt>true</tt>，也就是说，零值不会被添加到矩阵中。
   *
   */
  template <typename number>
  void
  add(const std::vector<size_type> &indices,
      const FullMatrix<number> &    full_matrix,
      const bool                    elide_zero_values = true);

  /**
   * 将FullMatrix中的所有元素添加到分别由<tt>row_indices</tt>和<tt>col_indices</tt>给出的全局位置。全局位置被转换为该块中的位置，如果全局索引没有指向#row和#column所指的块，则抛出ExcBlockIndexMismatch。
   * @todo  <tt>elide_zero_values</tt>目前被忽略。
   * 可选参数<tt>elide_zero_values</tt>可以用来指定是无论如何都要添加零值，还是要过滤掉这些数据，只添加非零值。默认值是<tt>true</tt>，也就是说，零值不会被添加到矩阵中。
   *
   */
  template <typename number>
  void
  add(const std::vector<size_type> &row_indices,
      const std::vector<size_type> &col_indices,
      const FullMatrix<number> &    full_matrix,
      const bool                    elide_zero_values = true);

  /**
   * 将矩阵指定行中的几个元素与<tt>col_indices</tt>给出的列索引设置为相应的值。这是为添加完整矩阵的做实际工作的函数。全局位置<tt>row_index</tt>和<tt>col_indices</tt>被转换为该块中的位置，如果全局索引没有指向#row和#column所指的块，则抛出ExcBlockIndexMismatch。
   * @todo  <tt>elide_zero_values</tt>目前被忽略。
   * 可选参数<tt>elide_zero_values</tt>可以用来指定是无论如何都要添加零值，还是要过滤掉这些数据，只添加非零值。默认值是<tt>true</tt>，也就是说，零值不会被添加到矩阵中。
   *
   */
  template <typename number>
  void
  add(const size_type               row_index,
      const std::vector<size_type> &col_indices,
      const std::vector<number> &   values,
      const bool                    elide_zero_values = true);

  /**
   * 在给定的全局矩阵行中，在稀疏矩阵中由col_indices指定的列中添加一个由<tt>values</tt>给出的数值阵列。
   * 可选的参数<tt>elide_zero_values</tt>可以用来指定是无论如何都要添加零值，还是要过滤掉这些数据，只添加非零值。默认值是<tt>true</tt>，也就是说，零值不会被添加到矩阵中。
   *
   */
  template <typename number>
  void
  add(const size_type  row,
      const size_type  n_cols,
      const size_type *col_indices,
      const number *   values,
      const bool       elide_zero_values      = true,
      const bool       col_indices_are_sorted = false);

  /**
   * 矩阵-向量-乘法，转发到MatrixType中的相同函数。没有进行索引计算，因此，向量需要有与#matrix匹配的大小。
   *
   */
  template <class VectorType>
  void
  vmult(VectorType &w, const VectorType &v) const;

  /**
   * 矩阵-向量-乘法，转发到MatrixType中的相同函数。没有做索引计算，因此，向量需要有与#matrix匹配的大小。
   *
   */
  template <class VectorType>
  void
  vmult_add(VectorType &w, const VectorType &v) const;

  /**
   * 矩阵-向量-乘法，转发到MatrixType中的相同函数。没有进行索引计算，因此，向量需要有与#matrix匹配的大小。
   *
   */
  template <class VectorType>
  void
  Tvmult(VectorType &w, const VectorType &v) const;

  /**
   * 矩阵-向量-乘法，转发到MatrixType中的相同函数。没有做索引计算，因此，向量需要有与#matrix匹配的大小。
   *
   */
  template <class VectorType>
  void
  Tvmult_add(VectorType &w, const VectorType &v) const;

  /**
   * 这个对象使用的内存。
   *
   */
  std::size_t
  memory_consumption() const;

  /**
   * 通过使用BlockIndices从索引中计算出的块编号与存储在此对象中的块坐标不匹配。
   *
   */
  DeclException2(ExcBlockIndexMismatch,
                 size_type,
                 size_type,
                 << "Block index " << arg1 << " does not match " << arg2);

  /**
   * 行坐标。 这是该数据成员矩阵在全局矩阵上的位置。
   *
   */
  size_type row;
  /**
   * 列坐标。 这是该数据成员矩阵在全局矩阵上的位置。
   *
   */
  size_type column;

  /**
   * 矩阵本身
   *
   */
  MatrixType matrix;

private:
  /**
   * 整个系统的行BlockIndices。使用row()，这可以让我们找到这个块的第一行自由度的索引。
   *
   */
  BlockIndices row_indices;
  /**
   * 整个系统的列BlockIndices。使用column()，我们可以找到该块的第一个列自由度的索引。
   *
   */
  BlockIndices column_indices;

  template <class OTHER_MatrixType>
  friend void
  dealii::internal::reinit(MatrixBlock<OTHER_MatrixType> &,
                           const BlockSparsityPattern &);

  template <typename number>
  friend void
  internal::reinit(MatrixBlock<dealii::SparseMatrix<number>> &v,
                   const BlockSparsityPattern &               p);
};


/**
 * 一个MatrixBlock的向量，使用共享指针来实现，以便于复制和重新排列。每个矩阵块都可以通过名称来识别。
 * @relatesalso  MatrixBlock
 *
 * @ingroup vector_valued
 *
 *
 */
template <typename MatrixType>
class MatrixBlockVector : private AnyData
{
public:
  /**
   * 声明容器尺寸的类型。
   *
   */
  using size_type = types::global_dof_index;

  /**
   * 存储对象的类型。
   *
   */
  using value_type = MatrixBlock<MatrixType>;

  /**
   * 用于存储对象的指针类型。我们使用一个分片的指针，这样，当不再使用时，它们会被自动删除。
   *
   */
  using ptr_type = std::shared_ptr<value_type>;

  /**
   * 在块系统中的<tt>(row,column)</tt>位置添加一个新的矩阵块。
   *
   */
  void
  add(size_type row, size_type column, const std::string &name);

  /**
   * 对于使用SparsityPattern的矩阵，这个函数用块系统的正确模式重新初始化向量中的每个矩阵。
   *
   */
  void
  reinit(const BlockSparsityPattern &sparsity);

  /**
   * 清除对象。
   * 由于通常只需要清除单个矩阵，而不需要清除块本身，所以有一个可选的参数。如果缺少这个参数或者
   * @p false,
   * ，所有的矩阵将被清空，但是这个对象的大小和块的位置将不会改变。如果
   * @p really_clean是 @p true,
   * ，那么该对象在最后将不包含任何块。
   *
   */
  void
  clear(bool really_clean = false);

  /**
   * 这个对象所使用的内存。
   *
   */
  std::size_t
  memory_consumption() const;

  /**
   * 访问位置<i>i</i>的块的常数引用。
   *
   */
  const value_type &
  block(size_type i) const;

  /**
   * 访问位置<i>i</i>的块的引用。
   *
   */
  value_type &
  block(size_type i);

  /**
   * 访问位置<i>i</i>的矩阵，进行读和写访问。
   *
   */
  MatrixType &
  matrix(size_type i);

  /**
   * 从私有基类中导入函数
   *
   */
  using AnyData::name;
  using AnyData::size;
  using AnyData::subscribe;
  using AnyData::unsubscribe;
};


/**
 * 一个MGLevelObject<MatrixBlock>的向量，使用共享指针来实现，以便于复制和重新排列。每个矩阵块都可以通过名称来识别。
 * @relatesalso  MatrixBlock
 *
 * @ingroup vector_valued
 *
 *
 */
template <typename MatrixType>
class MGMatrixBlockVector : public Subscriptor
{
public:
  /**
   * 声明容器尺寸的类型。
   *
   */
  using size_type = types::global_dof_index;

  /**
   * 存储对象的类型。
   *
   */
  using value_type = MGLevelObject<MatrixBlock<MatrixType>>;
  /**
   * 构造函数，决定哪些矩阵应该被存储。
   * 如果<tt>edge_matrices</tt>为真，则分配面的自由度离散的边缘矩阵的对象。
   * 如果<tt>edge_flux_matrices</tt>为真，则分配细化边上的危险通量对象。
   *
   */
  MGMatrixBlockVector(const bool edge_matrices      = false,
                      const bool edge_flux_matrices = false);

  /**
   * 块的数量。
   *
   */
  unsigned int
  size() const;

  /**
   * 在块系统中<tt>(行,列)</tt>的位置添加一个新的矩阵块。第三个参数允许给矩阵一个名称，以便以后识别。
   *
   */
  void
  add(size_type row, size_type column, const std::string &name);

  /**
   * 对于使用SparsityPattern的矩阵，这个函数用块系统的正确模式重新初始化向量中的每个矩阵。
   * 这个函数重新初始化了水平矩阵。
   *
   */
  void
  reinit_matrix(const MGLevelObject<BlockSparsityPattern> &sparsity);
  /**
   * 对于使用SparsityPattern的矩阵，此函数用块系统的正确模式重新初始化向量中的每个矩阵。
   * 这个函数对细化边上的自由度的矩阵进行重新初始化。
   *
   */
  void
  reinit_edge(const MGLevelObject<BlockSparsityPattern> &sparsity);
  /**
   * 对于使用SparsityPattern的矩阵，此函数用块系统的正确模式重新初始化向量中的每个矩阵。
   * 这个函数在细化边上重新初始化通量矩阵。
   *
   */
  void
  reinit_edge_flux(const MGLevelObject<BlockSparsityPattern> &sparsity);

  /**
   * 清除对象。
   * 因为通常只需要清除单个矩阵，而不需要清除块本身，所以有一个可选的参数。如果缺少这个参数或者
   * @p false,
   * ，所有的矩阵将被清空，但是这个对象的大小和块的位置将不会改变。如果
   * @p really_clean是 @p true,
   * ，那么该对象在最后将不包含任何块。
   *
   */
  void
  clear(bool really_clean = false);

  /**
   * 访问位于<i>i</i>位置的矩阵块的常量引用。
   *
   */
  const value_type &
  block(size_type i) const;

  /**
   * 访问位置为<i>i</i>的矩阵块的引用。
   *
   */
  value_type &
  block(size_type i);

  /**
   * 访问位置<i>i</i>的边缘矩阵块的一个常量引用。
   *
   */
  const value_type &
  block_in(size_type i) const;

  /**
   * 访问位置<i>i</i>的边缘矩阵块的一个引用。
   *
   */
  value_type &
  block_in(size_type i);

  /**
   * 访问位置<i>i</i>的边缘矩阵块的一个常量引用。
   *
   */
  const value_type &
  block_out(size_type i) const;

  /**
   * 访问位置<i>i</i>的边缘矩阵块的一个引用。
   *
   */
  value_type &
  block_out(size_type i);

  /**
   * 访问位置<i>i</i>的边缘通量矩阵块的一个常数参考。
   *
   */
  const value_type &
  block_up(size_type i) const;

  /**
   * 访问位置<i>i</i>的边缘通量矩阵块的引用。
   *
   */
  value_type &
  block_up(size_type i);

  /**
   * 访问位置<i>i</i>的边缘通量矩阵块的一个常数参考。
   *
   */
  const value_type &
  block_down(size_type i) const;

  /**
   * 访问位置<i>i</i>的边缘通量矩阵块的引用。
   *
   */
  value_type &
  block_down(size_type i);

  /**
   * 这个对象所使用的内存。
   *
   */
  std::size_t
  memory_consumption() const;

private:
  /// Clear one of the matrix objects
  void
  clear_object(AnyData &);

  /// Flag for storing matrices_in and matrices_out
  const bool edge_matrices;

  /// Flag for storing flux_matrices_up and flux_matrices_down
  const bool edge_flux_matrices;

  /// The level matrices
  AnyData matrices;
  /// The matrix from the interior of a level to the refinement edge
  AnyData matrices_in;
  /// The matrix from the refinement edge to the interior of a level
  AnyData matrices_out;
  /// The DG flux from a level to the lower level
  AnyData flux_matrices_down;
  /// The DG flux from the lower level to a level
  AnyData flux_matrices_up;
};


//----------------------------------------------------------------------//

namespace internal
{
  template <typename MatrixType>
  void
  reinit(MatrixBlock<MatrixType> &v, const BlockSparsityPattern &p)
  {
    v.row_indices    = p.get_row_indices();
    v.column_indices = p.get_column_indices();
  }


  template <typename number>
  void
  reinit(MatrixBlock<dealii::SparseMatrix<number>> &v,
         const BlockSparsityPattern &               p)
  {
    v.row_indices    = p.get_row_indices();
    v.column_indices = p.get_column_indices();
    v.matrix.reinit(p.block(v.row, v.column));
  }
} // namespace internal


template <typename MatrixType>
inline MatrixBlock<MatrixType>::MatrixBlock()
  : row(numbers::invalid_size_type)
  , column(numbers::invalid_size_type)
{}


template <typename MatrixType>
inline MatrixBlock<MatrixType>::MatrixBlock(size_type i, size_type j)
  : row(i)
  , column(j)
{}


template <typename MatrixType>
inline void
MatrixBlock<MatrixType>::reinit(const BlockSparsityPattern &sparsity)
{
  internal::reinit(*this, sparsity);
}


template <typename MatrixType>
inline MatrixBlock<MatrixType>::operator MatrixType &()
{
  return matrix;
}


template <typename MatrixType>
inline MatrixBlock<MatrixType>::operator const MatrixType &() const
{
  return matrix;
}


template <typename MatrixType>
inline void
MatrixBlock<MatrixType>::add(const size_type                       gi,
                             const size_type                       gj,
                             const typename MatrixType::value_type value)
{
  Assert(row_indices.size() != 0, ExcNotInitialized());
  Assert(column_indices.size() != 0, ExcNotInitialized());

  const std::pair<unsigned int, size_type> bi = row_indices.global_to_local(gi);
  const std::pair<unsigned int, size_type> bj =
    column_indices.global_to_local(gj);

  Assert(bi.first == row, ExcBlockIndexMismatch(bi.first, row));
  Assert(bj.first == column, ExcBlockIndexMismatch(bj.first, column));

  matrix.add(bi.second, bj.second, value);
}


template <typename MatrixType>
template <typename number>
inline void
MatrixBlock<MatrixType>::add(const std::vector<size_type> &r_indices,
                             const std::vector<size_type> &c_indices,
                             const FullMatrix<number> &    values,
                             const bool                    elide_zero_values)
{
  Assert(row_indices.size() != 0, ExcNotInitialized());
  Assert(column_indices.size() != 0, ExcNotInitialized());

  AssertDimension(r_indices.size(), values.m());
  AssertDimension(c_indices.size(), values.n());

  for (size_type i = 0; i < row_indices.size(); ++i)
    add(r_indices[i],
        c_indices.size(),
        c_indices.data(),
        &values(i, 0),
        elide_zero_values);
}


template <typename MatrixType>
template <typename number>
inline void
MatrixBlock<MatrixType>::add(const size_type  b_row,
                             const size_type  n_cols,
                             const size_type *col_indices,
                             const number *   values,
                             const bool,
                             const bool)
{
  Assert(row_indices.size() != 0, ExcNotInitialized());
  Assert(column_indices.size() != 0, ExcNotInitialized());

  const std::pair<unsigned int, size_type> bi =
    row_indices.global_to_local(b_row);

  // In debug mode, we check whether
  // all indices are in the correct
  // block.

  // Actually, for the time being, we
  // leave it at this. While it may
  // not be the most efficient way,
  // it is at least thread safe.
  //#ifdef DEBUG
  Assert(bi.first == row, ExcBlockIndexMismatch(bi.first, row));

  for (size_type j = 0; j < n_cols; ++j)
    {
      const std::pair<unsigned int, size_type> bj =
        column_indices.global_to_local(col_indices[j]);
      Assert(bj.first == column, ExcBlockIndexMismatch(bj.first, column));

      matrix.add(bi.second, bj.second, values[j]);
    }
  //#endif
}


template <typename MatrixType>
template <typename number>
inline void
MatrixBlock<MatrixType>::add(const std::vector<size_type> &indices,
                             const FullMatrix<number> &    values,
                             const bool                    elide_zero_values)
{
  Assert(row_indices.size() != 0, ExcNotInitialized());
  Assert(column_indices.size() != 0, ExcNotInitialized());

  AssertDimension(indices.size(), values.m());
  Assert(values.n() == values.m(), ExcNotQuadratic());

  for (size_type i = 0; i < indices.size(); ++i)
    add(indices[i],
        indices.size(),
        indices.data(),
        &values(i, 0),
        elide_zero_values);
}



template <typename MatrixType>
template <typename number>
inline void
MatrixBlock<MatrixType>::add(const size_type               row,
                             const std::vector<size_type> &col_indices,
                             const std::vector<number> &   values,
                             const bool                    elide_zero_values)
{
  Assert(row_indices.size() != 0, ExcNotInitialized());
  Assert(column_indices.size() != 0, ExcNotInitialized());

  AssertDimension(col_indices.size(), values.size());
  add(row,
      col_indices.size(),
      col_indices.data(),
      values.data(),
      elide_zero_values);
}


template <typename MatrixType>
template <class VectorType>
inline void
MatrixBlock<MatrixType>::vmult(VectorType &w, const VectorType &v) const
{
  matrix.vmult(w, v);
}


template <typename MatrixType>
template <class VectorType>
inline void
MatrixBlock<MatrixType>::vmult_add(VectorType &w, const VectorType &v) const
{
  matrix.vmult_add(w, v);
}


template <typename MatrixType>
template <class VectorType>
inline void
MatrixBlock<MatrixType>::Tvmult(VectorType &w, const VectorType &v) const
{
  matrix.Tvmult(w, v);
}


template <typename MatrixType>
template <class VectorType>
inline void
MatrixBlock<MatrixType>::Tvmult_add(VectorType &w, const VectorType &v) const
{
  matrix.Tvmult_add(w, v);
}


template <typename MatrixType>
inline std::size_t
MatrixBlock<MatrixType>::memory_consumption() const
{
  return (sizeof(*this) + MemoryConsumption::memory_consumption(matrix) -
          sizeof(matrix));
}

//----------------------------------------------------------------------//

template <typename MatrixType>
inline void
MatrixBlockVector<MatrixType>::add(size_type          row,
                                   size_type          column,
                                   const std::string &name)
{
  ptr_type p(new value_type(row, column));
  AnyData::add(p, name);
}


template <typename MatrixType>
inline void
MatrixBlockVector<MatrixType>::reinit(const BlockSparsityPattern &sparsity)
{
  for (size_type i = 0; i < this->size(); ++i)
    {
      block(i).reinit(sparsity);
    }
}


template <typename MatrixType>
inline void
MatrixBlockVector<MatrixType>::clear(bool really_clean)
{
  if (really_clean)
    {
      Assert(false, ExcNotImplemented());
    }
  else
    {
      for (size_type i = 0; i < this->size(); ++i)
        matrix(i).clear();
    }
}



template <typename MatrixType>
inline const MatrixBlock<MatrixType> &
MatrixBlockVector<MatrixType>::block(size_type i) const
{
  return *this->read<ptr_type>(i);
}


template <typename MatrixType>
inline MatrixBlock<MatrixType> &
MatrixBlockVector<MatrixType>::block(size_type i)
{
  return *this->entry<ptr_type>(i);
}


template <typename MatrixType>
inline MatrixType &
MatrixBlockVector<MatrixType>::matrix(size_type i)
{
  return this->entry<ptr_type>(i)->matrix;
}



//----------------------------------------------------------------------//

template <typename MatrixType>
inline MGMatrixBlockVector<MatrixType>::MGMatrixBlockVector(const bool e,
                                                            const bool f)
  : edge_matrices(e)
  , edge_flux_matrices(f)
{}


template <typename MatrixType>
inline unsigned int
MGMatrixBlockVector<MatrixType>::size() const
{
  return matrices.size();
}


template <typename MatrixType>
inline void
MGMatrixBlockVector<MatrixType>::add(size_type          row,
                                     size_type          column,
                                     const std::string &name)
{
  MGLevelObject<MatrixBlock<MatrixType>> p(0, 1);
  p[0].row    = row;
  p[0].column = column;

  matrices.add(p, name);
  if (edge_matrices)
    {
      matrices_in.add(p, name);
      matrices_out.add(p, name);
    }
  if (edge_flux_matrices)
    {
      flux_matrices_up.add(p, name);
      flux_matrices_down.add(p, name);
    }
}


template <typename MatrixType>
inline const MGLevelObject<MatrixBlock<MatrixType>> &
MGMatrixBlockVector<MatrixType>::block(size_type i) const
{
  return *matrices.read<const MGLevelObject<MatrixType> *>(i);
}


template <typename MatrixType>
inline MGLevelObject<MatrixBlock<MatrixType>> &
MGMatrixBlockVector<MatrixType>::block(size_type i)
{
  return *matrices.entry<MGLevelObject<MatrixType> *>(i);
}


template <typename MatrixType>
inline const MGLevelObject<MatrixBlock<MatrixType>> &
MGMatrixBlockVector<MatrixType>::block_in(size_type i) const
{
  return *matrices_in.read<const MGLevelObject<MatrixType> *>(i);
}


template <typename MatrixType>
inline MGLevelObject<MatrixBlock<MatrixType>> &
MGMatrixBlockVector<MatrixType>::block_in(size_type i)
{
  return *matrices_in.entry<MGLevelObject<MatrixType> *>(i);
}


template <typename MatrixType>
inline const MGLevelObject<MatrixBlock<MatrixType>> &
MGMatrixBlockVector<MatrixType>::block_out(size_type i) const
{
  return *matrices_out.read<const MGLevelObject<MatrixType> *>(i);
}


template <typename MatrixType>
inline MGLevelObject<MatrixBlock<MatrixType>> &
MGMatrixBlockVector<MatrixType>::block_out(size_type i)
{
  return *matrices_out.entry<MGLevelObject<MatrixType> *>(i);
}


template <typename MatrixType>
inline const MGLevelObject<MatrixBlock<MatrixType>> &
MGMatrixBlockVector<MatrixType>::block_up(size_type i) const
{
  return *flux_matrices_up.read<const MGLevelObject<MatrixType> *>(i);
}


template <typename MatrixType>
inline MGLevelObject<MatrixBlock<MatrixType>> &
MGMatrixBlockVector<MatrixType>::block_up(size_type i)
{
  return *flux_matrices_up.entry<MGLevelObject<MatrixType> *>(i);
}


template <typename MatrixType>
inline const MGLevelObject<MatrixBlock<MatrixType>> &
MGMatrixBlockVector<MatrixType>::block_down(size_type i) const
{
  return *flux_matrices_down.read<const MGLevelObject<MatrixType> *>(i);
}


template <typename MatrixType>
inline MGLevelObject<MatrixBlock<MatrixType>> &
MGMatrixBlockVector<MatrixType>::block_down(size_type i)
{
  return *flux_matrices_down.entry<MGLevelObject<MatrixType> *>(i);
}


template <typename MatrixType>
inline void
MGMatrixBlockVector<MatrixType>::reinit_matrix(
  const MGLevelObject<BlockSparsityPattern> &sparsity)
{
  for (size_type i = 0; i < this->size(); ++i)
    {
      MGLevelObject<MatrixBlock<MatrixType>> &o   = block(i);
      const size_type                         row = o[o.min_level()].row;
      const size_type                         col = o[o.min_level()].column;

      o.resize(sparsity.min_level(), sparsity.max_level());
      for (size_type level = o.min_level(); level <= o.max_level(); ++level)
        {
          o[level].row    = row;
          o[level].column = col;
          internal::reinit(o[level], sparsity[level]);
        }
    }
}


template <typename MatrixType>
inline void
MGMatrixBlockVector<MatrixType>::reinit_edge(
  const MGLevelObject<BlockSparsityPattern> &sparsity)
{
  for (size_type i = 0; i < this->size(); ++i)
    {
      MGLevelObject<MatrixBlock<MatrixType>> &o   = block(i);
      const size_type                         row = o[o.min_level()].row;
      const size_type                         col = o[o.min_level()].column;

      block_in(i).resize(sparsity.min_level(), sparsity.max_level());
      block_out(i).resize(sparsity.min_level(), sparsity.max_level());
      for (size_type level = o.min_level(); level <= o.max_level(); ++level)
        {
          block_in(i)[level].row    = row;
          block_in(i)[level].column = col;
          internal::reinit(block_in(i)[level], sparsity[level]);
          block_out(i)[level].row    = row;
          block_out(i)[level].column = col;
          internal::reinit(block_out(i)[level], sparsity[level]);
        }
    }
}


template <typename MatrixType>
inline void
MGMatrixBlockVector<MatrixType>::reinit_edge_flux(
  const MGLevelObject<BlockSparsityPattern> &sparsity)
{
  for (size_type i = 0; i < this->size(); ++i)
    {
      MGLevelObject<MatrixBlock<MatrixType>> &o   = block(i);
      const size_type                         row = o[o.min_level()].row;
      const size_type                         col = o[o.min_level()].column;

      block_up(i).resize(sparsity.min_level(), sparsity.max_level());
      block_down(i).resize(sparsity.min_level(), sparsity.max_level());
      for (size_type level = o.min_level(); level <= o.max_level(); ++level)
        {
          block_up(i)[level].row    = row;
          block_up(i)[level].column = col;
          internal::reinit(block_up(i)[level], sparsity[level]);
          block_down(i)[level].row    = row;
          block_down(i)[level].column = col;
          internal::reinit(block_down(i)[level], sparsity[level]);
        }
    }
}


template <typename MatrixType>
inline void
MGMatrixBlockVector<MatrixType>::clear_object(AnyData &mo)
{
  for (size_type i = 0; i < mo.size(); ++i)
    {
      MGLevelObject<MatrixBlock<MatrixType>> &o =
        mo.entry<MGLevelObject<MatrixType> *>(i);
      for (size_type level = o.min_level(); level <= o.max_level(); ++level)
        o[level].matrix.clear();
    }
}


template <typename MatrixType>
inline void
MGMatrixBlockVector<MatrixType>::clear(bool really_clean)
{
  if (really_clean)
    {
      Assert(false, ExcNotImplemented());
    }
  else
    {
      clear_object(matrices);
      clear_object(matrices_in);
      clear_object(matrices_out);
      clear_object(flux_matrices_up);
      clear_object(flux_matrices_down);
    }
}



DEAL_II_NAMESPACE_CLOSE

#endif


