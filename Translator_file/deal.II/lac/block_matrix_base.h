//include/deal.II-translator/lac/block_matrix_base_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2020 by the deal.II authors
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

#ifndef dealii_block_matrix_base_h
#define dealii_block_matrix_base_h


#include <deal.II/base/config.h>

#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/base/table.h>
#include <deal.II/base/thread_management.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/block_indices.h>
#include <deal.II/lac/exceptions.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/matrix_iterator.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/vector_operation.h>

#include <cmath>

DEAL_II_NAMESPACE_OPEN


// Forward declaration
#ifndef DOXYGEN
template <typename>
class MatrixIterator;
#endif


/*!   @addtogroup  Matrix1  
     * @{ 

* 
*
*/

/**
 * 实现块状矩阵中迭代器的名称空间。
 *
 *
 */
namespace BlockMatrixIterators
{
  /**
   * 块状矩阵访问器的基类，实现对矩阵的步进。
   *
   */
  template <class BlockMatrixType>
  class AccessorBase
  {
  public:
    /**
     * 声明容器大小的类型。
     *
     */
    using size_type = types::global_dof_index;

    /**
     * Typedef我们指向的矩阵的值类型。
     *
     */
    using value_type = typename BlockMatrixType::value_type;

    /**
     * 将数据字段初始化为默认值。
     *
     */
    AccessorBase();

    /**
     * 该对象所代表的元素的块行。
     *
     */
    unsigned int
    block_row() const;

    /**
     * 此对象所代表的元素的块列。
     *
     */
    unsigned int
    block_column() const;

  protected:
    /**
     * 我们目前所指向的块行。
     *
     */
    unsigned int row_block;

    /**
     * 我们现在所指向的块列。
     *
     */
    unsigned int col_block;

    // Let the iterator class be a friend.
    template <typename>
    friend class MatrixIterator;
  };



  /**
   * 块矩阵中的访问器类。
   *
   */
  template <class BlockMatrixType, bool Constness>
  class Accessor;


  /**
   * 非恒定矩阵的块状矩阵访问器。
   *
   */
  template <class BlockMatrixType>
  class Accessor<BlockMatrixType, false> : public AccessorBase<BlockMatrixType>
  {
  public:
    /**
     * 声明容器大小的类型。
     *
     */
    using size_type = types::global_dof_index;

    /**
     * 在这个访问器中使用的矩阵的类型。
     *
     */
    using MatrixType = BlockMatrixType;

    /**
     * Typedef我们指向的矩阵的值类型。
     *
     */
    using value_type = typename BlockMatrixType::value_type;

    /**
     * 构造器。因为我们只使用访问器进行读取访问，所以一个常量矩阵指针就足够了。
     * 将迭代器放在矩阵给定行的开头，如果 @p row
     * 等于矩阵的总行数，则创建末端指针。
     *
     */
    Accessor(BlockMatrixType *m, const size_type row, const size_type col);

    /**
     * 这个对象所代表的元素的行数。
     *
     */
    size_type
    row() const;

    /**
     * 这个对象所代表的元素的列号。
     *
     */
    size_type
    column() const;

    /**
     * 当前位置的条目值。
     *
     */
    value_type
    value() const;

    /**
     * 设置新值。
     *
     */
    void
    set_value(value_type newval) const;

  protected:
    /**
     * 访问的矩阵。
     *
     */
    BlockMatrixType *matrix;

    /**
     * 底层矩阵类的迭代器。
     *
     */
    typename BlockMatrixType::BlockType::iterator base_iterator;

    /**
     * 向前移动一个元素。
     *
     */
    void
    advance();

    /**
     * 将这个访问器与另一个访问器进行比较，以确定是否相等。
     *
     */
    bool
    operator==(const Accessor &a) const;

    template <typename>
    friend class MatrixIterator;
    friend class Accessor<BlockMatrixType, true>;
  };

  /**
   * 用于常数矩阵的块状矩阵访问器，实现对矩阵的步进。
   *
   */
  template <class BlockMatrixType>
  class Accessor<BlockMatrixType, true> : public AccessorBase<BlockMatrixType>
  {
  public:
    /**
     * 声明容器大小的类型。
     *
     */
    using size_type = types::global_dof_index;

    /**
     * 该访问器中使用的矩阵的类型。
     *
     */
    using MatrixType = const BlockMatrixType;

    /**
     * Typedef我们指向的矩阵的值类型。
     *
     */
    using value_type = typename BlockMatrixType::value_type;

    /**
     * 构造器。因为我们只使用访问器进行读取访问，所以一个常量矩阵指针就足够了。
     * 将迭代器放在矩阵给定行的开头，如果 @p row
     * 等于矩阵的总行数，则创建末端指针。
     *
     */
    Accessor(const BlockMatrixType *m,
             const size_type        row,
             const size_type        col);

    /**
     * 从非常量访问器初始化常量访问器。
     *
     */
    Accessor(const Accessor<BlockMatrixType, false> &);

    /**
     * 这个对象所代表的元素的行数。
     *
     */
    size_type
    row() const;

    /**
     * 这个对象所代表的元素的列号。
     *
     */
    size_type
    column() const;

    /**
     * 当前位置的条目值。
     *
     */
    value_type
    value() const;

  protected:
    /**
     * 访问的矩阵。
     *
     */
    const BlockMatrixType *matrix;

    /**
     * 底层矩阵类的迭代器。
     *
     */
    typename BlockMatrixType::BlockType::const_iterator base_iterator;

    /**
     * 向前移动一个元素。
     *
     */
    void
    advance();

    /**
     * 将这个访问器与另一个访问器进行比较，以确定是否相等。
     *
     */
    bool
    operator==(const Accessor &a) const;

    // Let the iterator class be a friend.
    template <typename>
    friend class dealii::MatrixIterator;
  };
} // namespace BlockMatrixIterators



/**
 * 封锁的矩阵类。这种类型的对象的行为几乎和通常的矩阵对象一样，大部分的功能都在两个类中实现。主要的区别是，这个对象所代表的矩阵是由一个矩阵数组（例如SparseMatrix<number>类型）组成的，对这个对象的元素的所有访问都被转发到对基础矩阵的访问。这个矩阵的各个块的实际类型是模板参数的类型，例如可以是通常的SparseMatrix或
 * PETScWrappers::SparseMatrix.  。
 * 除了通常的矩阵访问和线性代数函数外，还有一些函数block()，允许访问矩阵的不同块。例如，当你想实现舒尔补码方法或块状预处理时，这可能会有帮助，因为每个块都属于你目前正在离散化的方程的一个特定部分。
 * 请注意，块和行的数量是由使用的稀疏模式对象隐含决定的。
 * 当一个微分方程系统的解属于不同类别的变量时，这种类型的对象经常被使用。例如，斯托克斯或纳维尔-斯托克斯方程的解有
 * @p dim
 * 个速度分量和一个压力分量。在这种情况下，将线性方程组视为2x2块的系统可能是有意义的，人们可以在这个2x2块结构的基础上构建预调节器或求解器。在这种情况下，这个类可以帮助你，因为它允许把矩阵看作是一个大矩阵，或者看作是一些单独的块。
 *
 *  <h3>Inheriting from this class</h3>
 * 由于这个类只是将其调用转发给子对象（如果需要的话，在调整了表示哪个子对象的索引后），这个类完全独立于子对象的实际类型。然而，设置块状矩阵和销毁块状矩阵的函数必须在派生类中实现。这些函数也必须填充这个基类所提供的数据成员，因为它们在这个类中只是被动地使用。
 *
 * 大多数函数都需要一个向量或块向量参数。一般来说，只有当这个矩阵的各个块实现了各自在有关向量类型上操作的函数时，这些函数才能被成功编译。例如，如果你有一个在deal.II
 * SparseMatrix对象上的块状稀疏矩阵，那么你很可能无法用
 * PETScWrappers::SparseMatrix
 * 对象上的块状向量形成矩阵-向量乘法。如果你试图这样做，你可能会得到一系列的编译器错误。
 *
 *
 * @note
 * 这个模板的实例化提供给<tt>  @<float@>  和  @<double@></tt>;
 * 其他可以在应用程序中生成（见手册中 @ref Instantiations
 * 一节）。
 * @see   @ref GlossBlockLA  "块（线性代数）"
 *
 *
 */
template <typename MatrixType>
class BlockMatrixBase : public Subscriptor
{
public:
  /**
   * Typedef底层矩阵的类型。
   *
   */
  using BlockType = MatrixType;

  /**
   * 矩阵条目的类型。这些类似于标准库容器中的别名。
   *
   */
  using value_type      = typename BlockType::value_type;
  using real_type       = typename numbers::NumberTraits<value_type>::real_type;
  using pointer         = value_type *;
  using const_pointer   = const value_type *;
  using reference       = value_type &;
  using const_reference = const value_type &;
  using size_type       = types::global_dof_index;

  using iterator =
    MatrixIterator<BlockMatrixIterators::Accessor<BlockMatrixBase, false>>;

  using const_iterator =
    MatrixIterator<BlockMatrixIterators::Accessor<BlockMatrixBase, true>>;


  /**
   * 默认构造函数。
   *
   */
  BlockMatrixBase() = default;

  /**
   * 解构器。
   *
   */
  ~BlockMatrixBase() override;

  /**
   * 将作为参数的矩阵复制到当前对象中。
   * 复制矩阵是一个昂贵的操作，我们不希望通过编译器生成的代码意外发生
   * <code>operator=</code>
   * 。（例如，如果不小心声明了一个当前类型为<i>by
   * value</i>而不是<i>by
   * reference</i>的函数参数，就会发生这种情况）。复制矩阵的功能是在这个成员函数中实现的。因此，该类型对象的所有复制操作都需要一个明确的函数调用。
   * 源矩阵可以是一个任意类型的矩阵，只要其数据类型可以转换为该矩阵的数据类型。
   * 该函数返回一个对<tt>this</tt>的引用。
   *
   */
  template <class BlockMatrixType>
  BlockMatrixBase &
  copy_from(const BlockMatrixType &source);

  /**
   * 访问具有给定坐标的块。
   *
   */
  BlockType &
  block(const unsigned int row, const unsigned int column);


  /**
   * 访问具有给定坐标的区块。常量对象的版本。
   *
   */
  const BlockType &
  block(const unsigned int row, const unsigned int column) const;

  /**
   * 返回共域（或范围）空间的维数。注意，矩阵的维度是
   * $m \times n$  。
   *
   */
  size_type
  m() const;

  /**
   * 返回域空间的维度。请注意，矩阵的维度是 $m \times n$  .
   *
   */
  size_type
  n() const;


  /**
   * 返回一列中的块数。如果目前没有与此矩阵相关的稀疏模式，则返回0。
   *
   */
  unsigned int
  n_block_rows() const;

  /**
   * 返回一个行中的块数。如果目前没有与该矩阵相关的稀疏模式，则返回0。
   *
   */
  unsigned int
  n_block_cols() const;

  /**
   * 将元素<tt>(i,j)/tt>设置为<tt>值</tt>。如果该条目不存在或者<tt>value</tt>不是一个有限的数字，则抛出一个错误。尽管如此，它仍然允许在不存在的字段中存储零值。
   *
   */
  void
  set(const size_type i, const size_type j, const value_type value);

  /**
   * 将FullMatrix中给出的所有元素设置到<tt>indices</tt>给出的稀疏矩阵位置。换句话说，这个函数将<tt>full_matrix</tt>中的元素写入调用的矩阵中，对矩阵的行和列都使用<tt>indices</tt>指定的本地到全球的索引。这个函数假设一个二次稀疏矩阵和一个二次全矩阵，这是FE计算中的通常情况。
   * 可选参数<tt>elide_zero_values</tt>可以用来指定是无论如何都要设置零值，还是要过滤掉零值（如果存在的话，不改变相应元素中的先前内容）。默认值是<tt>false</tt>，也就是说，即使是零值也要处理。
   *
   */
  template <typename number>
  void
  set(const std::vector<size_type> &indices,
      const FullMatrix<number> &    full_matrix,
      const bool                    elide_zero_values = false);

  /**
   * 与之前的函数相同，但现在包括了使用矩形full_matrices的可能性，以及在行和列上分别使用不同的本地到全球索引。
   *
   */
  template <typename number>
  void
  set(const std::vector<size_type> &row_indices,
      const std::vector<size_type> &col_indices,
      const FullMatrix<number> &    full_matrix,
      const bool                    elide_zero_values = false);

  /**
   * 将矩阵的指定行中的几个元素与<tt>col_indices</tt>给出的列索引设置为相应的值。
   * 可选的参数<tt>elide_zero_values</tt>可以用来指定是无论如何都要设置零值，还是要过滤掉零值（如果存在的话，不改变相应元素中的先前内容）。默认值是<tt>false</tt>，也就是说，即使是零值也要处理。
   *
   */
  template <typename number>
  void
  set(const size_type               row,
      const std::vector<size_type> &col_indices,
      const std::vector<number> &   values,
      const bool                    elide_zero_values = false);

  /**
   * 将几个元素设置为由<tt>values</tt>给出的值，在给定的行和col_indices给出的列中设置为稀疏矩阵。
   * 可选的参数<tt>elide_zero_values</tt>可以用来指定是无论如何都要插入零值还是要过滤掉它们。默认值是<tt>false</tt>，也就是说，即使是零值也要插入/替换。
   *
   */
  template <typename number>
  void
  set(const size_type  row,
      const size_type  n_cols,
      const size_type *col_indices,
      const number *   values,
      const bool       elide_zero_values = false);

  /**
   * 向元素添加<tt>value</tt>（<i>i,j</i>）。
   * 如果该条目不存在或者<tt>value</tt>不是一个有限的数字，则抛出一个错误。尽管如此，它仍然允许在不存在的字段中存储零值。
   *
   */
  void
  add(const size_type i, const size_type j, const value_type value);

  /**
   * 将FullMatrix<double>中给出的所有元素添加到由<tt>indices</tt>给出的稀疏矩阵位置。换句话说，这个函数将<tt>full_matrix</tt>中的元素添加到调用矩阵的相应条目中，使用<tt>indices</tt>为矩阵的行和列指定的本地到全球索引。这个函数假定一个二次稀疏矩阵和一个二次全矩阵，这是FE计算中通常的情况。
   * 可选参数<tt>elide_zero_values</tt>可以用来指定是无论如何都要添加零值，还是要过滤掉这些零值，只添加非零数据。默认值是<tt>true</tt>，也就是说，零值不会被添加到矩阵中。
   *
   */
  template <typename number>
  void
  add(const std::vector<size_type> &indices,
      const FullMatrix<number> &    full_matrix,
      const bool                    elide_zero_values = true);

  /**
   * 与之前的函数相同，但现在包括了使用矩形full_matrices的可能性，以及在行和列上分别使用不同的本地到全球索引。
   *
   */
  template <typename number>
  void
  add(const std::vector<size_type> &row_indices,
      const std::vector<size_type> &col_indices,
      const FullMatrix<number> &    full_matrix,
      const bool                    elide_zero_values = true);

  /**
   * 将矩阵的指定行中的几个元素与<tt>col_indices</tt>给出的列索引设置为相应的值。
   * 可选参数<tt>elide_zero_values</tt>可以用来指定是无论如何都要添加零值，还是要过滤掉这些零值，只添加非零数据。默认值是<tt>true</tt>，也就是说，零值不会被添加到矩阵中。
   *
   */
  template <typename number>
  void
  add(const size_type               row,
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
   * 将<tt>matrix</tt>按<tt>factor</tt>的比例添加到这个矩阵中，也就是说，矩阵<tt>factor*matrix</tt>被添加到<tt>this</tt>。如果调用矩阵的稀疏性模式不包含输入矩阵的稀疏性模式中的所有元素，这个函数将抛出一个异常。
   * 然而，根据MatrixType，可能会出现额外的限制。
   * 一些稀疏矩阵格式要求<tt>matrix</tt>是基于与调用矩阵相同的稀疏模式的。
   *
   */
  void
  add(const value_type factor, const BlockMatrixBase<MatrixType> &matrix);

  /**
   * 返回条目(i,j)的值。
   * 这可能是一个昂贵的操作，你应该始终注意在哪里调用这个函数。
   * 为了避免滥用，如果想要的元素在矩阵中不存在，这个函数会抛出一个异常。
   *
   */
  value_type
  operator()(const size_type i, const size_type j) const;

  /**
   * 这个函数主要像operator()()，它返回矩阵条目<tt>(i,j)</tt>的值。唯一的区别是，如果这个条目不存在于稀疏模式中，那么就不会引发异常，而是返回0。虽然这在某些情况下可能很方便，但请注意，由于没有使用矩阵的稀疏性，写出的算法与最优解相比很简单，很慢。
   *
   */
  value_type
  el(const size_type i, const size_type j) const;

  /**
   * 返回第<i>i</i>行中的主对角线元素。如果矩阵不是二次方的，以及矩阵的对角线块不是二次方的，这个函数会抛出一个错误。
   * 这个函数比operator()()快得多，因为对于二次矩阵来说，对角线条目可能是每行中第一个被存储的，因此访问时不需要搜索正确的列号。
   *
   */
  value_type
  diag_element(const size_type i) const;

  /**
   * 在矩阵的所有子块上调用compress()函数。      参见 @ref GlossCompress "压缩分布式对象 "
   * 以获得更多信息。
   *
   */
  void
  compress(::dealii::VectorOperation::values operation);

  /**
   * 将整个矩阵乘以一个固定系数。
   *
   */
  BlockMatrixBase &
  operator*=(const value_type factor);

  /**
   * 将整个矩阵除以一个固定系数。
   *
   */
  BlockMatrixBase &
  operator/=(const value_type factor);

  /**
   * 加法 矩阵-向量乘法。在 $dst$ 上添加 $M*src$ ， $M$
   * 为该矩阵。
   *
   */
  template <class BlockVectorType>
  void
  vmult_add(BlockVectorType &dst, const BlockVectorType &src) const;

  /**
   * 添加矩阵-向量乘法。将<i>M<sup>T</sup>src</i>加到<i>dst</i>上，<i>M</i>是这个矩阵。这个函数与vmult_add()的作用相同，但需要转置的矩阵。
   *
   */
  template <class BlockVectorType>
  void
  Tvmult_add(BlockVectorType &dst, const BlockVectorType &src) const;

  /**
   * 返回向量<i>v</i>相对于该矩阵诱导的规范，即<i>v<sup>T</sup>Mv)</i>的规范。这很有用，例如在有限元背景下，一个函数的<i>L<sup>T</sup></i>规范等于相对于代表有限元函数节点值的向量矩阵的矩阵规范。请注意，尽管函数的名称可能暗示了一些不同的东西，但由于历史原因，返回的不是法线，而是它的平方，正如上面所定义的标量积。
   * 很明显，对于这个操作，矩阵需要是平方的。
   *
   */
  template <class BlockVectorType>
  value_type
  matrix_norm_square(const BlockVectorType &v) const;

  /**
   * 返回矩阵的frobenius
   * norm，即矩阵中所有条目的平方和的平方根。
   *
   */
  real_type
  frobenius_norm() const;

  /**
   * 计算矩阵标量乘积  $\left(u,Mv\right)$  。
   *
   */
  template <class BlockVectorType>
  value_type
  matrix_scalar_product(const BlockVectorType &u,
                        const BlockVectorType &v) const;

  /**
   * 计算残差<i>r=b-Ax</i>。将残差写进<tt>dst</tt>。
   *
   */
  template <class BlockVectorType>
  value_type
  residual(BlockVectorType &      dst,
           const BlockVectorType &x,
           const BlockVectorType &b) const;

  /**
   * 打印矩阵到给定的数据流中，使用格式<tt>(line,col)
   * value</tt>，即每行有一个矩阵的非零条目。可选的标志是根据底层稀疏矩阵的类型，以不同的风格输出稀疏模式。
   *
   */
  void
  print(std::ostream &out, const bool alternative_output = false) const;

  /**
   * 迭代器从第一个条目开始。
   *
   */
  iterator
  begin();

  /**
   * 最后的迭代器。
   *
   */
  iterator
  end();

  /**
   * 从<tt>r</tt>行的第一个条目开始的迭代器。
   *
   */
  iterator
  begin(const size_type r);

  /**
   * 行<tt>r</tt>的最终迭代器。
   *
   */
  iterator
  end(const size_type r);
  /**
   * 从第一条开始的迭代器。
   *
   */
  const_iterator
  begin() const;

  /**
   * 最后的迭代器。
   *
   */
  const_iterator
  end() const;

  /**
   * 从<tt>r</tt>行的第一个条目开始的迭代器。
   *
   */
  const_iterator
  begin(const size_type r) const;

  /**
   * 行<tt>r</tt>的最终迭代器。
   *
   */
  const_iterator
  end(const size_type r) const;

  /**
   * 返回对行的底层BlockIndices数据的引用。
   *
   */
  const BlockIndices &
  get_row_indices() const;

  /**
   * 返回对列的基本BlockIndices数据的引用。
   *
   */
  const BlockIndices &
  get_column_indices() const;

  /**
   * 确定此对象的内存消耗（以字节为单位）的估计值。注意，如果是在基于MPI的程序中调用，则只返回当前处理器上保留的内存。
   *
   */
  std::size_t
  memory_consumption() const;

  /**
   * @addtogroup  Exceptions  
     * @{ 
   *
   */

  /**
   * 异常情况
   *
   */
  DeclException4(ExcIncompatibleRowNumbers,
                 int,
                 int,
                 int,
                 int,
                 << "The blocks [" << arg1 << ',' << arg2 << "] and [" << arg3
                 << ',' << arg4 << "] have differing row numbers.");
  /**
   * 异常情况
   *
   */
  DeclException4(ExcIncompatibleColNumbers,
                 int,
                 int,
                 int,
                 int,
                 << "The blocks [" << arg1 << ',' << arg2 << "] and [" << arg3
                 << ',' << arg4 << "] have differing column numbers.");
  //@}
protected:
  /**
   * 释放所有内存并返回到与调用默认构造函数后相同的状态。它也忘记了它之前绑定的稀疏模式。
   * 这个函数对所有的子矩阵进行清除，然后将这个对象重置为完全没有块。
   * 这个函数是受保护的，因为它可能需要释放额外的结构。如果足够的话，一个派生类可以再次将其公开。
   *
   */
  void
  clear();

  /**
   * 行和列的索引数组。
   *
   */
  BlockIndices row_block_indices;
  BlockIndices column_block_indices;

  /**
   * 子矩阵的数组。
   *
   */
  Table<2, SmartPointer<BlockType, BlockMatrixBase<MatrixType>>> sub_objects;

  /**
   * 这个函数收集了子对象的大小，并将其存储在内部数组中，以便能够将矩阵的全局索引转为子对象的索引。在你改变了子对象的大小之后，你必须*每次都调用这个函数。
   * 每当子对象的大小发生变化， @p X_block_indices
   * 数组需要更新时，派生类应该调用这个函数。
   * 注意，这个函数不是公开的，因为不是所有的派生类都需要导出其接口。例如，对于通常的deal.II
   * SparseMatrix类，每当调用reinit()时，大小都是隐式确定的，而且单个块不能被调整大小。因此，对于该类，这个函数不必是公共的。另一方面，对于PETSc类来说，没有相关的稀疏模式对象来决定块的大小，对于这些类来说，这个函数必须是公开可用的。因此，这些类导出了这个函数。
   *
   */
  void
  collect_sizes();

  /**
   * 矩阵-向量乘法：让 $dst = M*src$ 与 $M$ 是这个矩阵。
   * 由于在vmult/Tvmult函数的块和非块版本之间衍生模板参数的问题，实际的函数是在派生类中实现的，实现者将调用转发给这里提供的实现者，在一个独特的名称下，模板参数可以被编译器衍生。
   *
   */
  template <class BlockVectorType>
  void
  vmult_block_block(BlockVectorType &dst, const BlockVectorType &src) const;

  /**
   * 矩阵-向量乘法。就像之前的函数一样，但只适用于矩阵只有一个块列的情况。
   * 由于在vmult/Tvmult函数的块和非块版本之间派生模板参数的问题，实际的函数是在派生类中实现的，实现者将调用转发给这里提供的实现者，在一个独特的名称下，模板参数可以被编译器派生。
   *
   */
  template <class BlockVectorType, class VectorType>
  void
  vmult_block_nonblock(BlockVectorType &dst, const VectorType &src) const;

  /**
   * 矩阵-向量乘法。就像之前的函数一样，但只适用于矩阵只有一个块行的情况。
   * 由于在vmult/Tvmult函数的块和非块版本之间派生模板参数的问题，实际的函数是在派生类中实现的，实现者将调用转发给这里提供的实现者，在一个独特的名字下，模板参数可以被编译器派生。
   *
   */
  template <class BlockVectorType, class VectorType>
  void
  vmult_nonblock_block(VectorType &dst, const BlockVectorType &src) const;

  /**
   * 矩阵-向量乘法。就像前面的函数，但只适用于矩阵只有一个块的情况。
   * 由于在vmult/Tvmult函数的块和非块版本之间派生模板参数的问题，实际的函数是在派生类中实现的，实现者将调用转发给这里提供的实现者，其模板参数可以由编译器派生出来的唯一名称。
   *
   */
  template <class VectorType>
  void
  vmult_nonblock_nonblock(VectorType &dst, const VectorType &src) const;

  /**
   * 矩阵-向量乘法：让 $dst = M^T*src$ 与 $M$
   * 是这个矩阵。这个函数与vmult()的作用相同，但需要转置的矩阵。
   * 由于vmult/Tvmult函数的块和非块版本之间的模板参数派生问题，实际的函数是在派生类中实现的，实现者将调用转发给这里提供的实现者，其模板参数可以由编译器派生出来的唯一名称。
   *
   */
  template <class BlockVectorType>
  void
  Tvmult_block_block(BlockVectorType &dst, const BlockVectorType &src) const;

  /**
   * 矩阵-向量乘法。就像前面的函数，但只适用于矩阵只有一个块行的情况。
   * 由于在vmult/Tvmult函数的块和非块版本之间派生模板参数的问题，实际的函数是在派生类中实现的，实现者将调用转发给这里提供的实现者，在一个唯一的名称下，模板参数可以被编译器派生。
   *
   */
  template <class BlockVectorType, class VectorType>
  void
  Tvmult_block_nonblock(BlockVectorType &dst, const VectorType &src) const;

  /**
   * 矩阵-向量乘法。就像前面的函数一样，但只适用于矩阵只有一个块列的情况。
   * 由于在vmult/Tvmult函数的块和非块版本之间派生模板参数的问题，实际的函数是在派生类中实现的，实现者将调用转发给这里提供的实现者，在一个独特的名称下，模板参数可以被编译器派生。
   *
   */
  template <class BlockVectorType, class VectorType>
  void
  Tvmult_nonblock_block(VectorType &dst, const BlockVectorType &src) const;

  /**
   * 矩阵-向量乘法。就像前面的函数，但只适用于矩阵只有一个块的情况。
   * 由于在vmult/Tvmult函数的块和非块版本之间派生模板参数的问题，实际的函数是在派生类中实现的，实现者将调用转发给这里提供的实现者，其模板参数可以由编译器派生的唯一名称。
   *
   */
  template <class VectorType>
  void
  Tvmult_nonblock_nonblock(VectorType &dst, const VectorType &src) const;


protected:
  /**
   * 一些矩阵类型，特别是PETSc，需要同步设置和添加操作。这必须对BlockMatrix中的所有矩阵进行操作。本例程通过通知所有块来准备添加元素。在添加元素之前，所有的内部例程都会调用它。
   *
   */
  void
  prepare_add_operation();

  /**
   * 通知所有区块，让它们为设置元素做准备，见prepare_add_operation()。
   *
   */
  void
  prepare_set_operation();


private:
  /**
   * 一个包含set()和add()函数使用的一些字段的结构，用于对输入字段进行预排序。由于人们可以合理地期望从多个线程同时调用set()和add()，只要所触及的矩阵索引是不相干的，这些临时数据字段需要由一个mutex来保护；因此该结构包含这样一个mutex作为成员变量。
   *
   */
  struct TemporaryData
  {
    /**
     * 临时向量，用于计算在做集体添加或设置时写入各个块的元素。
     *
     */
    std::vector<size_type> counter_within_block;

    /**
     * 临时向量，用于在每个稀疏矩阵上将局部数据写入全局数据时，每个块上的列索引。
     *
     */
    std::vector<std::vector<size_type>> column_indices;

    /**
     * 用于存储本地值的临时向量（在将本地值写入全局时需要对它们进行重新排序）。
     *
     */
    std::vector<std::vector<value_type>> column_values;

    /**
     * 一个突变变量，用于保护对该结构的成员变量的访问。
     *
     */
    std::mutex mutex;

    /**
     * 拷贝操作符。需要这样做是因为这个类的默认拷贝操作符被删除了（因为
     * std::mutex
     * 是不可拷贝的），因此包围类的默认拷贝操作符也被删除。
     * 这里的实现只是没有做任何事情
     *
     * - TemporaryData对象只是在开始使用时被调整大小的抓取对象，所以实际上没有必要复制任何东西。
     *
     */
    TemporaryData &
    operator=(const TemporaryData &)
    {
      return *this;
    }
  };

  /**
   * 一组可以被add()和set()函数使用的从头开始的数组，这些函数在使用前采取数据的指针来预排序索引。来自多个线程的访问通过作为该结构一部分的mutex变量进行同步。
   *
   */
  TemporaryData temporary_data;

  // Make the iterator class a friend. We have to work around a compiler bug
  // here again.
  template <typename, bool>
  friend class BlockMatrixIterators::Accessor;

  template <typename>
  friend class MatrixIterator;
};


 /*@}*/ 

#ifndef DOXYGEN
 /* ------------------------- Template functions ---------------------- */ 


namespace BlockMatrixIterators
{
  template <class BlockMatrixType>
  inline AccessorBase<BlockMatrixType>::AccessorBase()
    : row_block(0)
    , col_block(0)
  {}


  template <class BlockMatrixType>
  inline unsigned int
  AccessorBase<BlockMatrixType>::block_row() const
  {
    Assert(row_block != numbers::invalid_unsigned_int, ExcIteratorPastEnd());

    return row_block;
  }


  template <class BlockMatrixType>
  inline unsigned int
  AccessorBase<BlockMatrixType>::block_column() const
  {
    Assert(col_block != numbers::invalid_unsigned_int, ExcIteratorPastEnd());

    return col_block;
  }


  template <class BlockMatrixType>
  inline Accessor<BlockMatrixType, true>::Accessor(
    const BlockMatrixType *matrix,
    const size_type        row,
    const size_type        col)
    : matrix(matrix)
    , base_iterator(matrix->block(0, 0).begin())
  {
    (void)col;
    Assert(col == 0, ExcNotImplemented());

    // check if this is a regular row or
    // the end of the matrix
    if (row < matrix->m())
      {
        const std::pair<unsigned int, size_type> indices =
          matrix->row_block_indices.global_to_local(row);

        // find the first block that does
        // have an entry in this row
        for (unsigned int bc = 0; bc < matrix->n_block_cols(); ++bc)
          {
            base_iterator =
              matrix->block(indices.first, bc).begin(indices.second);
            if (base_iterator !=
                matrix->block(indices.first, bc).end(indices.second))
              {
                this->row_block = indices.first;
                this->col_block = bc;
                return;
              }
          }

        // hm, there is no block that has
        // an entry in this column. we need
        // to take the next entry then,
        // which may be the first entry of
        // the next row, or recursively the
        // next row, or so on
        *this = Accessor(matrix, row + 1, 0);
      }
    else
      {
        // we were asked to create the end
        // iterator for this matrix
        this->row_block = numbers::invalid_unsigned_int;
        this->col_block = numbers::invalid_unsigned_int;
      }
  }


  //   template <class BlockMatrixType>
  //   inline
  //   Accessor<BlockMatrixType, true>::Accessor (const
  //   Accessor<BlockMatrixType, true>& other)
  //                :
  //                matrix(other.matrix),
  //                base_iterator(other.base_iterator)
  //   {
  //     this->row_block = other.row_block;
  //     this->col_block = other.col_block;
  //   }


  template <class BlockMatrixType>
  inline Accessor<BlockMatrixType, true>::Accessor(
    const Accessor<BlockMatrixType, false> &other)
    : matrix(other.matrix)
    , base_iterator(other.base_iterator)
  {
    this->row_block = other.row_block;
    this->col_block = other.col_block;
  }


  template <class BlockMatrixType>
  inline typename Accessor<BlockMatrixType, true>::size_type
  Accessor<BlockMatrixType, true>::row() const
  {
    Assert(this->row_block != numbers::invalid_unsigned_int,
           ExcIteratorPastEnd());

    return (matrix->row_block_indices.local_to_global(this->row_block, 0) +
            base_iterator->row());
  }


  template <class BlockMatrixType>
  inline typename Accessor<BlockMatrixType, true>::size_type
  Accessor<BlockMatrixType, true>::column() const
  {
    Assert(this->col_block != numbers::invalid_unsigned_int,
           ExcIteratorPastEnd());

    return (matrix->column_block_indices.local_to_global(this->col_block, 0) +
            base_iterator->column());
  }


  template <class BlockMatrixType>
  inline typename Accessor<BlockMatrixType, true>::value_type
  Accessor<BlockMatrixType, true>::value() const
  {
    Assert(this->row_block != numbers::invalid_unsigned_int,
           ExcIteratorPastEnd());
    Assert(this->col_block != numbers::invalid_unsigned_int,
           ExcIteratorPastEnd());

    return base_iterator->value();
  }



  template <class BlockMatrixType>
  inline void
  Accessor<BlockMatrixType, true>::advance()
  {
    Assert(this->row_block != numbers::invalid_unsigned_int,
           ExcIteratorPastEnd());
    Assert(this->col_block != numbers::invalid_unsigned_int,
           ExcIteratorPastEnd());

    // Remember current row inside block
    size_type local_row = base_iterator->row();

    // Advance one element inside the
    // current block
    ++base_iterator;

    // while we hit the end of the row of a
    // block (which may happen multiple
    // times if rows inside a block are
    // empty), we have to jump to the next
    // block and take the
    while (base_iterator ==
           matrix->block(this->row_block, this->col_block).end(local_row))
      {
        // jump to next block in this block
        // row, if possible, otherwise go
        // to next row
        if (this->col_block < matrix->n_block_cols() - 1)
          {
            ++this->col_block;
            base_iterator =
              matrix->block(this->row_block, this->col_block).begin(local_row);
          }
        else
          {
            // jump back to next row in
            // first block column
            this->col_block = 0;
            ++local_row;

            // see if this has brought us
            // past the number of rows in
            // this block. if so see
            // whether we've just fallen
            // off the end of the whole
            // matrix
            if (local_row ==
                matrix->block(this->row_block, this->col_block).m())
              {
                local_row = 0;
                ++this->row_block;
                if (this->row_block == matrix->n_block_rows())
                  {
                    this->row_block = numbers::invalid_unsigned_int;
                    this->col_block = numbers::invalid_unsigned_int;
                    return;
                  }
              }

            base_iterator =
              matrix->block(this->row_block, this->col_block).begin(local_row);
          }
      }
  }


  template <class BlockMatrixType>
  inline bool
  Accessor<BlockMatrixType, true>::operator==(const Accessor &a) const
  {
    if (matrix != a.matrix)
      return false;

    if (this->row_block == a.row_block && this->col_block == a.col_block)
      // end iterators do not necessarily
      // have to have the same
      // base_iterator representation, but
      // valid iterators have to
      return (((this->row_block == numbers::invalid_unsigned_int) &&
               (this->col_block == numbers::invalid_unsigned_int)) ||
              (base_iterator == a.base_iterator));

    return false;
  }

  //----------------------------------------------------------------------//


  template <class BlockMatrixType>
  inline Accessor<BlockMatrixType, false>::Accessor(BlockMatrixType *matrix,
                                                    const size_type  row,
                                                    const size_type  col)
    : matrix(matrix)
    , base_iterator(matrix->block(0, 0).begin())
  {
    (void)col;
    Assert(col == 0, ExcNotImplemented());
    // check if this is a regular row or
    // the end of the matrix
    if (row < matrix->m())
      {
        const std::pair<unsigned int, size_type> indices =
          matrix->row_block_indices.global_to_local(row);

        // find the first block that does
        // have an entry in this row
        for (size_type bc = 0; bc < matrix->n_block_cols(); ++bc)
          {
            base_iterator =
              matrix->block(indices.first, bc).begin(indices.second);
            if (base_iterator !=
                matrix->block(indices.first, bc).end(indices.second))
              {
                this->row_block = indices.first;
                this->col_block = bc;
                return;
              }
          }

        // hm, there is no block that has
        // an entry in this column. we need
        // to take the next entry then,
        // which may be the first entry of
        // the next row, or recursively the
        // next row, or so on
        *this = Accessor(matrix, row + 1, 0);
      }
    else
      {
        // we were asked to create the end
        // iterator for this matrix
        this->row_block = numbers::invalid_size_type;
        this->col_block = numbers::invalid_size_type;
      }
  }


  template <class BlockMatrixType>
  inline typename Accessor<BlockMatrixType, false>::size_type
  Accessor<BlockMatrixType, false>::row() const
  {
    Assert(this->row_block != numbers::invalid_size_type, ExcIteratorPastEnd());

    return (matrix->row_block_indices.local_to_global(this->row_block, 0) +
            base_iterator->row());
  }


  template <class BlockMatrixType>
  inline typename Accessor<BlockMatrixType, false>::size_type
  Accessor<BlockMatrixType, false>::column() const
  {
    Assert(this->col_block != numbers::invalid_size_type, ExcIteratorPastEnd());

    return (matrix->column_block_indices.local_to_global(this->col_block, 0) +
            base_iterator->column());
  }


  template <class BlockMatrixType>
  inline typename Accessor<BlockMatrixType, false>::value_type
  Accessor<BlockMatrixType, false>::value() const
  {
    Assert(this->row_block != numbers::invalid_size_type, ExcIteratorPastEnd());
    Assert(this->col_block != numbers::invalid_size_type, ExcIteratorPastEnd());

    return base_iterator->value();
  }



  template <class BlockMatrixType>
  inline void
  Accessor<BlockMatrixType, false>::set_value(
    typename Accessor<BlockMatrixType, false>::value_type newval) const
  {
    Assert(this->row_block != numbers::invalid_size_type, ExcIteratorPastEnd());
    Assert(this->col_block != numbers::invalid_size_type, ExcIteratorPastEnd());

    base_iterator->value() = newval;
  }



  template <class BlockMatrixType>
  inline void
  Accessor<BlockMatrixType, false>::advance()
  {
    Assert(this->row_block != numbers::invalid_size_type, ExcIteratorPastEnd());
    Assert(this->col_block != numbers::invalid_size_type, ExcIteratorPastEnd());

    // Remember current row inside block
    size_type local_row = base_iterator->row();

    // Advance one element inside the
    // current block
    ++base_iterator;

    // while we hit the end of the row of a
    // block (which may happen multiple
    // times if rows inside a block are
    // empty), we have to jump to the next
    // block and take the
    while (base_iterator ==
           matrix->block(this->row_block, this->col_block).end(local_row))
      {
        // jump to next block in this block
        // row, if possible, otherwise go
        // to next row
        if (this->col_block < matrix->n_block_cols() - 1)
          {
            ++this->col_block;
            base_iterator =
              matrix->block(this->row_block, this->col_block).begin(local_row);
          }
        else
          {
            // jump back to next row in
            // first block column
            this->col_block = 0;
            ++local_row;

            // see if this has brought us
            // past the number of rows in
            // this block. if so see
            // whether we've just fallen
            // off the end of the whole
            // matrix
            if (local_row ==
                matrix->block(this->row_block, this->col_block).m())
              {
                local_row = 0;
                ++this->row_block;
                if (this->row_block == matrix->n_block_rows())
                  {
                    this->row_block = numbers::invalid_size_type;
                    this->col_block = numbers::invalid_size_type;
                    return;
                  }
              }

            base_iterator =
              matrix->block(this->row_block, this->col_block).begin(local_row);
          }
      }
  }



  template <class BlockMatrixType>
  inline bool
  Accessor<BlockMatrixType, false>::operator==(const Accessor &a) const
  {
    if (matrix != a.matrix)
      return false;

    if (this->row_block == a.row_block && this->col_block == a.col_block)
      // end iterators do not necessarily
      // have to have the same
      // base_iterator representation, but
      // valid iterators have to
      return (((this->row_block == numbers::invalid_size_type) &&
               (this->col_block == numbers::invalid_size_type)) ||
              (base_iterator == a.base_iterator));

    return false;
  }
} // namespace BlockMatrixIterators


//---------------------------------------------------------------------------

template <typename MatrixType>
inline BlockMatrixBase<MatrixType>::~BlockMatrixBase()
{
  try
    {
      clear();
    }
  catch (...)
    {}
}


template <class MatrixType>
template <class BlockMatrixType>
inline BlockMatrixBase<MatrixType> &
BlockMatrixBase<MatrixType>::copy_from(const BlockMatrixType &source)
{
  for (unsigned int r = 0; r < n_block_rows(); ++r)
    for (unsigned int c = 0; c < n_block_cols(); ++c)
      block(r, c).copy_from(source.block(r, c));

  return *this;
}


template <class MatrixType>
std::size_t
BlockMatrixBase<MatrixType>::memory_consumption() const
{
  std::size_t mem =
    MemoryConsumption::memory_consumption(row_block_indices) +
    MemoryConsumption::memory_consumption(column_block_indices) +
    MemoryConsumption::memory_consumption(sub_objects) +
    MemoryConsumption::memory_consumption(temporary_data.counter_within_block) +
    MemoryConsumption::memory_consumption(temporary_data.column_indices) +
    MemoryConsumption::memory_consumption(temporary_data.column_values) +
    sizeof(temporary_data.mutex);

  for (unsigned int r = 0; r < n_block_rows(); ++r)
    for (unsigned int c = 0; c < n_block_cols(); ++c)
      {
        MatrixType *p = this->sub_objects[r][c];
        mem += MemoryConsumption::memory_consumption(*p);
      }

  return mem;
}



template <class MatrixType>
inline void
BlockMatrixBase<MatrixType>::clear()
{
  for (unsigned int r = 0; r < n_block_rows(); ++r)
    for (unsigned int c = 0; c < n_block_cols(); ++c)
      {
        MatrixType *p           = this->sub_objects[r][c];
        this->sub_objects[r][c] = nullptr;
        delete p;
      }
  sub_objects.reinit(0, 0);

  // reset block indices to empty
  row_block_indices = column_block_indices = BlockIndices();
}



template <class MatrixType>
inline typename BlockMatrixBase<MatrixType>::BlockType &
BlockMatrixBase<MatrixType>::block(const unsigned int row,
                                   const unsigned int column)
{
  AssertIndexRange(row, n_block_rows());
  AssertIndexRange(column, n_block_cols());

  return *sub_objects[row][column];
}



template <class MatrixType>
inline const typename BlockMatrixBase<MatrixType>::BlockType &
BlockMatrixBase<MatrixType>::block(const unsigned int row,
                                   const unsigned int column) const
{
  AssertIndexRange(row, n_block_rows());
  AssertIndexRange(column, n_block_cols());

  return *sub_objects[row][column];
}


template <class MatrixType>
inline typename BlockMatrixBase<MatrixType>::size_type
BlockMatrixBase<MatrixType>::m() const
{
  return row_block_indices.total_size();
}



template <class MatrixType>
inline typename BlockMatrixBase<MatrixType>::size_type
BlockMatrixBase<MatrixType>::n() const
{
  return column_block_indices.total_size();
}



template <class MatrixType>
inline unsigned int
BlockMatrixBase<MatrixType>::n_block_cols() const
{
  return column_block_indices.size();
}



template <class MatrixType>
inline unsigned int
BlockMatrixBase<MatrixType>::n_block_rows() const
{
  return row_block_indices.size();
}



// Write the single set manually,
// since the other function has a lot
// of overhead in that case.
template <class MatrixType>
inline void
BlockMatrixBase<MatrixType>::set(const size_type  i,
                                 const size_type  j,
                                 const value_type value)
{
  prepare_set_operation();

  AssertIsFinite(value);

  const std::pair<unsigned int, size_type>
    row_index = row_block_indices.global_to_local(i),
    col_index = column_block_indices.global_to_local(j);
  block(row_index.first, col_index.first)
    .set(row_index.second, col_index.second, value);
}



template <class MatrixType>
template <typename number>
inline void
BlockMatrixBase<MatrixType>::set(const std::vector<size_type> &row_indices,
                                 const std::vector<size_type> &col_indices,
                                 const FullMatrix<number> &    values,
                                 const bool elide_zero_values)
{
  Assert(row_indices.size() == values.m(),
         ExcDimensionMismatch(row_indices.size(), values.m()));
  Assert(col_indices.size() == values.n(),
         ExcDimensionMismatch(col_indices.size(), values.n()));

  for (size_type i = 0; i < row_indices.size(); ++i)
    set(row_indices[i],
        col_indices.size(),
        col_indices.data(),
        &values(i, 0),
        elide_zero_values);
}



template <class MatrixType>
template <typename number>
inline void
BlockMatrixBase<MatrixType>::set(const std::vector<size_type> &indices,
                                 const FullMatrix<number> &    values,
                                 const bool elide_zero_values)
{
  Assert(indices.size() == values.m(),
         ExcDimensionMismatch(indices.size(), values.m()));
  Assert(values.n() == values.m(), ExcNotQuadratic());

  for (size_type i = 0; i < indices.size(); ++i)
    set(indices[i],
        indices.size(),
        indices.data(),
        &values(i, 0),
        elide_zero_values);
}



template <class MatrixType>
template <typename number>
inline void
BlockMatrixBase<MatrixType>::set(const size_type               row,
                                 const std::vector<size_type> &col_indices,
                                 const std::vector<number> &   values,
                                 const bool elide_zero_values)
{
  Assert(col_indices.size() == values.size(),
         ExcDimensionMismatch(col_indices.size(), values.size()));

  set(row,
      col_indices.size(),
      col_indices.data(),
      values.data(),
      elide_zero_values);
}



// This is a very messy function, since
// we need to calculate to each position
// the location in the global array.
template <class MatrixType>
template <typename number>
inline void
BlockMatrixBase<MatrixType>::set(const size_type  row,
                                 const size_type  n_cols,
                                 const size_type *col_indices,
                                 const number *   values,
                                 const bool       elide_zero_values)
{
  prepare_set_operation();

  // lock access to the temporary data structure to
  // allow multiple threads to call this function concurrently
  std::lock_guard<std::mutex> lock(temporary_data.mutex);

  // Resize scratch arrays
  if (temporary_data.column_indices.size() < this->n_block_cols())
    {
      temporary_data.column_indices.resize(this->n_block_cols());
      temporary_data.column_values.resize(this->n_block_cols());
      temporary_data.counter_within_block.resize(this->n_block_cols());
    }

  // Resize sub-arrays to n_cols. This
  // is a bit wasteful, but we resize
  // only a few times (then the maximum
  // row length won't increase that
  // much any more). At least we know
  // that all arrays are going to be of
  // the same size, so we can check
  // whether the size of one is large
  // enough before actually going
  // through all of them.
  if (temporary_data.column_indices[0].size() < n_cols)
    {
      for (unsigned int i = 0; i < this->n_block_cols(); ++i)
        {
          temporary_data.column_indices[i].resize(n_cols);
          temporary_data.column_values[i].resize(n_cols);
        }
    }

  // Reset the number of added elements
  // in each block to zero.
  for (unsigned int i = 0; i < this->n_block_cols(); ++i)
    temporary_data.counter_within_block[i] = 0;

  // Go through the column indices to
  // find out which portions of the
  // values should be set in which
  // block of the matrix. We need to
  // touch all the data, since we can't
  // be sure that the data of one block
  // is stored contiguously (in fact,
  // indices will be intermixed when it
  // comes from an element matrix).
  for (size_type j = 0; j < n_cols; ++j)
    {
      number value = values[j];

      if (value == number() && elide_zero_values == true)
        continue;

      const std::pair<unsigned int, size_type> col_index =
        this->column_block_indices.global_to_local(col_indices[j]);

      const size_type local_index =
        temporary_data.counter_within_block[col_index.first]++;

      temporary_data.column_indices[col_index.first][local_index] =
        col_index.second;
      temporary_data.column_values[col_index.first][local_index] = value;
    }

#  ifdef DEBUG
  // If in debug mode, do a check whether
  // the right length has been obtained.
  size_type length = 0;
  for (unsigned int i = 0; i < this->n_block_cols(); ++i)
    length += temporary_data.counter_within_block[i];
  Assert(length <= n_cols, ExcInternalError());
#  endif

  // Now we found out about where the
  // individual columns should start and
  // where we should start reading out
  // data. Now let's write the data into
  // the individual blocks!
  const std::pair<unsigned int, size_type> row_index =
    this->row_block_indices.global_to_local(row);
  for (unsigned int block_col = 0; block_col < n_block_cols(); ++block_col)
    {
      if (temporary_data.counter_within_block[block_col] == 0)
        continue;

      block(row_index.first, block_col)
        .set(row_index.second,
             temporary_data.counter_within_block[block_col],
             temporary_data.column_indices[block_col].data(),
             temporary_data.column_values[block_col].data(),
             false);
    }
}



template <class MatrixType>
inline void
BlockMatrixBase<MatrixType>::add(const size_type  i,
                                 const size_type  j,
                                 const value_type value)
{
  AssertIsFinite(value);

  prepare_add_operation();

  // save some cycles for zero additions, but
  // only if it is safe for the matrix we are
  // working with
  using MatrixTraits = typename MatrixType::Traits;
  if ((MatrixTraits::zero_addition_can_be_elided == true) &&
      (value == value_type()))
    return;

  const std::pair<unsigned int, size_type>
    row_index = row_block_indices.global_to_local(i),
    col_index = column_block_indices.global_to_local(j);
  block(row_index.first, col_index.first)
    .add(row_index.second, col_index.second, value);
}



template <class MatrixType>
template <typename number>
inline void
BlockMatrixBase<MatrixType>::add(const std::vector<size_type> &row_indices,
                                 const std::vector<size_type> &col_indices,
                                 const FullMatrix<number> &    values,
                                 const bool elide_zero_values)
{
  Assert(row_indices.size() == values.m(),
         ExcDimensionMismatch(row_indices.size(), values.m()));
  Assert(col_indices.size() == values.n(),
         ExcDimensionMismatch(col_indices.size(), values.n()));

  for (size_type i = 0; i < row_indices.size(); ++i)
    add(row_indices[i],
        col_indices.size(),
        col_indices.data(),
        &values(i, 0),
        elide_zero_values);
}



template <class MatrixType>
template <typename number>
inline void
BlockMatrixBase<MatrixType>::add(const std::vector<size_type> &indices,
                                 const FullMatrix<number> &    values,
                                 const bool elide_zero_values)
{
  Assert(indices.size() == values.m(),
         ExcDimensionMismatch(indices.size(), values.m()));
  Assert(values.n() == values.m(), ExcNotQuadratic());

  for (size_type i = 0; i < indices.size(); ++i)
    add(indices[i],
        indices.size(),
        indices.data(),
        &values(i, 0),
        elide_zero_values);
}



template <class MatrixType>
template <typename number>
inline void
BlockMatrixBase<MatrixType>::add(const size_type               row,
                                 const std::vector<size_type> &col_indices,
                                 const std::vector<number> &   values,
                                 const bool elide_zero_values)
{
  Assert(col_indices.size() == values.size(),
         ExcDimensionMismatch(col_indices.size(), values.size()));

  add(row,
      col_indices.size(),
      col_indices.data(),
      values.data(),
      elide_zero_values);
}



// This is a very messy function, since
// we need to calculate to each position
// the location in the global array.
template <class MatrixType>
template <typename number>
inline void
BlockMatrixBase<MatrixType>::add(const size_type  row,
                                 const size_type  n_cols,
                                 const size_type *col_indices,
                                 const number *   values,
                                 const bool       elide_zero_values,
                                 const bool       col_indices_are_sorted)
{
  prepare_add_operation();

  // TODO: Look over this to find out
  // whether we can do that more
  // efficiently.
  if (col_indices_are_sorted == true)
    {
#  ifdef DEBUG
      // check whether indices really are
      // sorted.
      size_type before = col_indices[0];
      for (size_type i = 1; i < n_cols; ++i)
        if (col_indices[i] <= before)
          Assert(false,
                 ExcMessage("Flag col_indices_are_sorted is set, but "
                            "indices appear to not be sorted.")) else before =
            col_indices[i];
#  endif
      const std::pair<unsigned int, size_type> row_index =
        this->row_block_indices.global_to_local(row);

      if (this->n_block_cols() > 1)
        {
          const size_type *first_block =
            Utilities::lower_bound(col_indices,
                                   col_indices + n_cols,
                                   this->column_block_indices.block_start(1));

          const size_type n_zero_block_indices = first_block - col_indices;
          block(row_index.first, 0)
            .add(row_index.second,
                 n_zero_block_indices,
                 col_indices,
                 values,
                 elide_zero_values,
                 col_indices_are_sorted);

          if (n_zero_block_indices < n_cols)
            this->add(row,
                      n_cols - n_zero_block_indices,
                      first_block,
                      values + n_zero_block_indices,
                      elide_zero_values,
                      false);
        }
      else
        {
          block(row_index.first, 0)
            .add(row_index.second,
                 n_cols,
                 col_indices,
                 values,
                 elide_zero_values,
                 col_indices_are_sorted);
        }

      return;
    }

  // Lock scratch arrays, then resize them
  std::lock_guard<std::mutex> lock(temporary_data.mutex);

  if (temporary_data.column_indices.size() < this->n_block_cols())
    {
      temporary_data.column_indices.resize(this->n_block_cols());
      temporary_data.column_values.resize(this->n_block_cols());
      temporary_data.counter_within_block.resize(this->n_block_cols());
    }

  // Resize sub-arrays to n_cols. This
  // is a bit wasteful, but we resize
  // only a few times (then the maximum
  // row length won't increase that
  // much any more). At least we know
  // that all arrays are going to be of
  // the same size, so we can check
  // whether the size of one is large
  // enough before actually going
  // through all of them.
  if (temporary_data.column_indices[0].size() < n_cols)
    {
      for (unsigned int i = 0; i < this->n_block_cols(); ++i)
        {
          temporary_data.column_indices[i].resize(n_cols);
          temporary_data.column_values[i].resize(n_cols);
        }
    }

  // Reset the number of added elements
  // in each block to zero.
  for (unsigned int i = 0; i < this->n_block_cols(); ++i)
    temporary_data.counter_within_block[i] = 0;

  // Go through the column indices to
  // find out which portions of the
  // values should be written into
  // which block of the matrix. We need
  // to touch all the data, since we
  // can't be sure that the data of one
  // block is stored contiguously (in
  // fact, data will be intermixed when
  // it comes from an element matrix).
  for (size_type j = 0; j < n_cols; ++j)
    {
      number value = values[j];

      if (value == number() && elide_zero_values == true)
        continue;

      const std::pair<unsigned int, size_type> col_index =
        this->column_block_indices.global_to_local(col_indices[j]);

      const size_type local_index =
        temporary_data.counter_within_block[col_index.first]++;

      temporary_data.column_indices[col_index.first][local_index] =
        col_index.second;
      temporary_data.column_values[col_index.first][local_index] = value;
    }

#  ifdef DEBUG
  // If in debug mode, do a check whether
  // the right length has been obtained.
  size_type length = 0;
  for (unsigned int i = 0; i < this->n_block_cols(); ++i)
    length += temporary_data.counter_within_block[i];
  Assert(length <= n_cols, ExcInternalError());
#  endif

  // Now we found out about where the
  // individual columns should start and
  // where we should start reading out
  // data. Now let's write the data into
  // the individual blocks!
  const std::pair<unsigned int, size_type> row_index =
    this->row_block_indices.global_to_local(row);
  for (unsigned int block_col = 0; block_col < n_block_cols(); ++block_col)
    {
      if (temporary_data.counter_within_block[block_col] == 0)
        continue;

      block(row_index.first, block_col)
        .add(row_index.second,
             temporary_data.counter_within_block[block_col],
             temporary_data.column_indices[block_col].data(),
             temporary_data.column_values[block_col].data(),
             false,
             col_indices_are_sorted);
    }
}



template <class MatrixType>
inline void
BlockMatrixBase<MatrixType>::add(const value_type                   factor,
                                 const BlockMatrixBase<MatrixType> &matrix)
{
  AssertIsFinite(factor);

  prepare_add_operation();

  // save some cycles for zero additions, but
  // only if it is safe for the matrix we are
  // working with
  using MatrixTraits = typename MatrixType::Traits;
  if ((MatrixTraits::zero_addition_can_be_elided == true) && (factor == 0))
    return;

  for (unsigned int row = 0; row < n_block_rows(); ++row)
    for (unsigned int col = 0; col < n_block_cols(); ++col)
      // This function should throw if the sparsity
      // patterns of the two blocks differ
      block(row, col).add(factor, matrix.block(row, col));
}



template <class MatrixType>
inline typename BlockMatrixBase<MatrixType>::value_type
BlockMatrixBase<MatrixType>::operator()(const size_type i,
                                        const size_type j) const
{
  const std::pair<unsigned int, size_type>
    row_index = row_block_indices.global_to_local(i),
    col_index = column_block_indices.global_to_local(j);
  return block(row_index.first, col_index.first)(row_index.second,
                                                 col_index.second);
}



template <class MatrixType>
inline typename BlockMatrixBase<MatrixType>::value_type
BlockMatrixBase<MatrixType>::el(const size_type i, const size_type j) const
{
  const std::pair<unsigned int, size_type>
    row_index = row_block_indices.global_to_local(i),
    col_index = column_block_indices.global_to_local(j);
  return block(row_index.first, col_index.first)
    .el(row_index.second, col_index.second);
}



template <class MatrixType>
inline typename BlockMatrixBase<MatrixType>::value_type
BlockMatrixBase<MatrixType>::diag_element(const size_type i) const
{
  Assert(n_block_rows() == n_block_cols(), ExcNotQuadratic());

  const std::pair<unsigned int, size_type> index =
    row_block_indices.global_to_local(i);
  return block(index.first, index.first).diag_element(index.second);
}



template <class MatrixType>
inline void
BlockMatrixBase<MatrixType>::compress(
  ::dealii::VectorOperation::values operation)
{
  for (unsigned int r = 0; r < n_block_rows(); ++r)
    for (unsigned int c = 0; c < n_block_cols(); ++c)
      block(r, c).compress(operation);
}



template <class MatrixType>
inline BlockMatrixBase<MatrixType> &
BlockMatrixBase<MatrixType>::operator*=(const value_type factor)
{
  Assert(n_block_cols() != 0, ExcNotInitialized());
  Assert(n_block_rows() != 0, ExcNotInitialized());

  for (unsigned int r = 0; r < n_block_rows(); ++r)
    for (unsigned int c = 0; c < n_block_cols(); ++c)
      block(r, c) *= factor;

  return *this;
}



template <class MatrixType>
inline BlockMatrixBase<MatrixType> &
BlockMatrixBase<MatrixType>::operator/=(const value_type factor)
{
  Assert(n_block_cols() != 0, ExcNotInitialized());
  Assert(n_block_rows() != 0, ExcNotInitialized());
  Assert(factor != 0, ExcDivideByZero());

  const value_type factor_inv = 1. / factor;

  for (unsigned int r = 0; r < n_block_rows(); ++r)
    for (unsigned int c = 0; c < n_block_cols(); ++c)
      block(r, c) *= factor_inv;

  return *this;
}



template <class MatrixType>
const BlockIndices &
BlockMatrixBase<MatrixType>::get_row_indices() const
{
  return this->row_block_indices;
}



template <class MatrixType>
const BlockIndices &
BlockMatrixBase<MatrixType>::get_column_indices() const
{
  return this->column_block_indices;
}



template <class MatrixType>
template <class BlockVectorType>
void
BlockMatrixBase<MatrixType>::vmult_block_block(BlockVectorType &      dst,
                                               const BlockVectorType &src) const
{
  Assert(dst.n_blocks() == n_block_rows(),
         ExcDimensionMismatch(dst.n_blocks(), n_block_rows()));
  Assert(src.n_blocks() == n_block_cols(),
         ExcDimensionMismatch(src.n_blocks(), n_block_cols()));

  for (size_type row = 0; row < n_block_rows(); ++row)
    {
      block(row, 0).vmult(dst.block(row), src.block(0));
      for (size_type col = 1; col < n_block_cols(); ++col)
        block(row, col).vmult_add(dst.block(row), src.block(col));
    };
}



template <class MatrixType>
template <class BlockVectorType, class VectorType>
void
BlockMatrixBase<MatrixType>::vmult_nonblock_block(
  VectorType &           dst,
  const BlockVectorType &src) const
{
  Assert(n_block_rows() == 1, ExcDimensionMismatch(1, n_block_rows()));
  Assert(src.n_blocks() == n_block_cols(),
         ExcDimensionMismatch(src.n_blocks(), n_block_cols()));

  block(0, 0).vmult(dst, src.block(0));
  for (size_type col = 1; col < n_block_cols(); ++col)
    block(0, col).vmult_add(dst, src.block(col));
}



template <class MatrixType>
template <class BlockVectorType, class VectorType>
void
BlockMatrixBase<MatrixType>::vmult_block_nonblock(BlockVectorType & dst,
                                                  const VectorType &src) const
{
  Assert(dst.n_blocks() == n_block_rows(),
         ExcDimensionMismatch(dst.n_blocks(), n_block_rows()));
  Assert(1 == n_block_cols(), ExcDimensionMismatch(1, n_block_cols()));

  for (size_type row = 0; row < n_block_rows(); ++row)
    block(row, 0).vmult(dst.block(row), src);
}



template <class MatrixType>
template <class VectorType>
void
BlockMatrixBase<MatrixType>::vmult_nonblock_nonblock(
  VectorType &      dst,
  const VectorType &src) const
{
  Assert(1 == n_block_rows(), ExcDimensionMismatch(1, n_block_rows()));
  Assert(1 == n_block_cols(), ExcDimensionMismatch(1, n_block_cols()));

  block(0, 0).vmult(dst, src);
}



template <class MatrixType>
template <class BlockVectorType>
void
BlockMatrixBase<MatrixType>::vmult_add(BlockVectorType &      dst,
                                       const BlockVectorType &src) const
{
  Assert(dst.n_blocks() == n_block_rows(),
         ExcDimensionMismatch(dst.n_blocks(), n_block_rows()));
  Assert(src.n_blocks() == n_block_cols(),
         ExcDimensionMismatch(src.n_blocks(), n_block_cols()));

  for (unsigned int row = 0; row < n_block_rows(); ++row)
    for (unsigned int col = 0; col < n_block_cols(); ++col)
      block(row, col).vmult_add(dst.block(row), src.block(col));
}



template <class MatrixType>
template <class BlockVectorType>
void
BlockMatrixBase<MatrixType>::Tvmult_block_block(
  BlockVectorType &      dst,
  const BlockVectorType &src) const
{
  Assert(dst.n_blocks() == n_block_cols(),
         ExcDimensionMismatch(dst.n_blocks(), n_block_cols()));
  Assert(src.n_blocks() == n_block_rows(),
         ExcDimensionMismatch(src.n_blocks(), n_block_rows()));

  dst = 0.;

  for (unsigned int row = 0; row < n_block_rows(); ++row)
    {
      for (unsigned int col = 0; col < n_block_cols(); ++col)
        block(row, col).Tvmult_add(dst.block(col), src.block(row));
    };
}



template <class MatrixType>
template <class BlockVectorType, class VectorType>
void
BlockMatrixBase<MatrixType>::Tvmult_block_nonblock(BlockVectorType & dst,
                                                   const VectorType &src) const
{
  Assert(dst.n_blocks() == n_block_cols(),
         ExcDimensionMismatch(dst.n_blocks(), n_block_cols()));
  Assert(1 == n_block_rows(), ExcDimensionMismatch(1, n_block_rows()));

  dst = 0.;

  for (unsigned int col = 0; col < n_block_cols(); ++col)
    block(0, col).Tvmult_add(dst.block(col), src);
}



template <class MatrixType>
template <class BlockVectorType, class VectorType>
void
BlockMatrixBase<MatrixType>::Tvmult_nonblock_block(
  VectorType &           dst,
  const BlockVectorType &src) const
{
  Assert(1 == n_block_cols(), ExcDimensionMismatch(1, n_block_cols()));
  Assert(src.n_blocks() == n_block_rows(),
         ExcDimensionMismatch(src.n_blocks(), n_block_rows()));

  block(0, 0).Tvmult(dst, src.block(0));

  for (size_type row = 1; row < n_block_rows(); ++row)
    block(row, 0).Tvmult_add(dst, src.block(row));
}



template <class MatrixType>
template <class VectorType>
void
BlockMatrixBase<MatrixType>::Tvmult_nonblock_nonblock(
  VectorType &      dst,
  const VectorType &src) const
{
  Assert(1 == n_block_cols(), ExcDimensionMismatch(1, n_block_cols()));
  Assert(1 == n_block_rows(), ExcDimensionMismatch(1, n_block_rows()));

  block(0, 0).Tvmult(dst, src);
}



template <class MatrixType>
template <class BlockVectorType>
void
BlockMatrixBase<MatrixType>::Tvmult_add(BlockVectorType &      dst,
                                        const BlockVectorType &src) const
{
  Assert(dst.n_blocks() == n_block_cols(),
         ExcDimensionMismatch(dst.n_blocks(), n_block_cols()));
  Assert(src.n_blocks() == n_block_rows(),
         ExcDimensionMismatch(src.n_blocks(), n_block_rows()));

  for (unsigned int row = 0; row < n_block_rows(); ++row)
    for (unsigned int col = 0; col < n_block_cols(); ++col)
      block(row, col).Tvmult_add(dst.block(col), src.block(row));
}



template <class MatrixType>
template <class BlockVectorType>
typename BlockMatrixBase<MatrixType>::value_type
BlockMatrixBase<MatrixType>::matrix_norm_square(const BlockVectorType &v) const
{
  Assert(n_block_rows() == n_block_cols(), ExcNotQuadratic());
  Assert(v.n_blocks() == n_block_rows(),
         ExcDimensionMismatch(v.n_blocks(), n_block_rows()));

  value_type norm_sqr = 0;
  for (unsigned int row = 0; row < n_block_rows(); ++row)
    for (unsigned int col = 0; col < n_block_cols(); ++col)
      if (row == col)
        norm_sqr += block(row, col).matrix_norm_square(v.block(row));
      else
        norm_sqr +=
          block(row, col).matrix_scalar_product(v.block(row), v.block(col));
  return norm_sqr;
}



template <class MatrixType>
typename BlockMatrixBase<MatrixType>::real_type
BlockMatrixBase<MatrixType>::frobenius_norm() const
{
  value_type norm_sqr = 0;

  // For each block, get the Frobenius norm, and add the square to the
  // accumulator for the full matrix
  for (unsigned int row = 0; row < n_block_rows(); ++row)
    {
      for (unsigned int col = 0; col < n_block_cols(); ++col)
        {
          const value_type block_norm = block(row, col).frobenius_norm();
          norm_sqr += block_norm * block_norm;
        }
    }

  return std::sqrt(norm_sqr);
}



template <class MatrixType>
template <class BlockVectorType>
typename BlockMatrixBase<MatrixType>::value_type
BlockMatrixBase<MatrixType>::matrix_scalar_product(
  const BlockVectorType &u,
  const BlockVectorType &v) const
{
  Assert(u.n_blocks() == n_block_rows(),
         ExcDimensionMismatch(u.n_blocks(), n_block_rows()));
  Assert(v.n_blocks() == n_block_cols(),
         ExcDimensionMismatch(v.n_blocks(), n_block_cols()));

  value_type result = 0;
  for (unsigned int row = 0; row < n_block_rows(); ++row)
    for (unsigned int col = 0; col < n_block_cols(); ++col)
      result +=
        block(row, col).matrix_scalar_product(u.block(row), v.block(col));
  return result;
}



template <class MatrixType>
template <class BlockVectorType>
typename BlockMatrixBase<MatrixType>::value_type
BlockMatrixBase<MatrixType>::residual(BlockVectorType &      dst,
                                      const BlockVectorType &x,
                                      const BlockVectorType &b) const
{
  Assert(dst.n_blocks() == n_block_rows(),
         ExcDimensionMismatch(dst.n_blocks(), n_block_rows()));
  Assert(b.n_blocks() == n_block_rows(),
         ExcDimensionMismatch(b.n_blocks(), n_block_rows()));
  Assert(x.n_blocks() == n_block_cols(),
         ExcDimensionMismatch(x.n_blocks(), n_block_cols()));
  // in block notation, the residual is
  // r_i = b_i - \sum_j A_ij x_j.
  // this can be written as
  // r_i = b_i - A_i0 x_0 - \sum_{j>0} A_ij x_j.
  //
  // for the first two terms, we can
  // call the residual function of
  // A_i0. for the other terms, we
  // use vmult_add. however, we want
  // to subtract, so in order to
  // avoid a temporary vector, we
  // perform a sign change of the
  // first two term before, and after
  // adding up
  for (unsigned int row = 0; row < n_block_rows(); ++row)
    {
      block(row, 0).residual(dst.block(row), x.block(0), b.block(row));

      for (size_type i = 0; i < dst.block(row).size(); ++i)
        dst.block(row)(i) = -dst.block(row)(i);

      for (unsigned int col = 1; col < n_block_cols(); ++col)
        block(row, col).vmult_add(dst.block(row), x.block(col));

      for (size_type i = 0; i < dst.block(row).size(); ++i)
        dst.block(row)(i) = -dst.block(row)(i);
    };

  value_type res = 0;
  for (size_type row = 0; row < n_block_rows(); ++row)
    res += dst.block(row).norm_sqr();
  return std::sqrt(res);
}



template <class MatrixType>
inline void
BlockMatrixBase<MatrixType>::print(std::ostream &out,
                                   const bool    alternative_output) const
{
  for (unsigned int row = 0; row < n_block_rows(); ++row)
    for (unsigned int col = 0; col < n_block_cols(); ++col)
      {
        if (!alternative_output)
          out << "Block (" << row << ", " << col << ")" << std::endl;

        block(row, col).print(out, alternative_output);
      }
}



template <class MatrixType>
inline typename BlockMatrixBase<MatrixType>::const_iterator
BlockMatrixBase<MatrixType>::begin() const
{
  return const_iterator(this, 0);
}



template <class MatrixType>
inline typename BlockMatrixBase<MatrixType>::const_iterator
BlockMatrixBase<MatrixType>::end() const
{
  return const_iterator(this, m());
}



template <class MatrixType>
inline typename BlockMatrixBase<MatrixType>::const_iterator
BlockMatrixBase<MatrixType>::begin(const size_type r) const
{
  AssertIndexRange(r, m());
  return const_iterator(this, r);
}



template <class MatrixType>
inline typename BlockMatrixBase<MatrixType>::const_iterator
BlockMatrixBase<MatrixType>::end(const size_type r) const
{
  AssertIndexRange(r, m());
  return const_iterator(this, r + 1);
}



template <class MatrixType>
inline typename BlockMatrixBase<MatrixType>::iterator
BlockMatrixBase<MatrixType>::begin()
{
  return iterator(this, 0);
}



template <class MatrixType>
inline typename BlockMatrixBase<MatrixType>::iterator
BlockMatrixBase<MatrixType>::end()
{
  return iterator(this, m());
}



template <class MatrixType>
inline typename BlockMatrixBase<MatrixType>::iterator
BlockMatrixBase<MatrixType>::begin(const size_type r)
{
  AssertIndexRange(r, m());
  return iterator(this, r);
}



template <class MatrixType>
inline typename BlockMatrixBase<MatrixType>::iterator
BlockMatrixBase<MatrixType>::end(const size_type r)
{
  AssertIndexRange(r, m());
  return iterator(this, r + 1);
}



template <class MatrixType>
void
BlockMatrixBase<MatrixType>::collect_sizes()
{
  std::vector<size_type> row_sizes(this->n_block_rows());
  std::vector<size_type> col_sizes(this->n_block_cols());

  // first find out the row sizes
  // from the first block column
  for (unsigned int r = 0; r < this->n_block_rows(); ++r)
    row_sizes[r] = sub_objects[r][0]->m();
  // then check that the following
  // block columns have the same
  // sizes
  for (unsigned int c = 1; c < this->n_block_cols(); ++c)
    for (unsigned int r = 0; r < this->n_block_rows(); ++r)
      Assert(row_sizes[r] == sub_objects[r][c]->m(),
             ExcIncompatibleRowNumbers(r, 0, r, c));

  // finally initialize the row
  // indices with this array
  this->row_block_indices.reinit(row_sizes);


  // then do the same with the columns
  for (unsigned int c = 0; c < this->n_block_cols(); ++c)
    col_sizes[c] = sub_objects[0][c]->n();
  for (unsigned int r = 1; r < this->n_block_rows(); ++r)
    for (unsigned int c = 0; c < this->n_block_cols(); ++c)
      Assert(col_sizes[c] == sub_objects[r][c]->n(),
             ExcIncompatibleRowNumbers(0, c, r, c));

  // finally initialize the row
  // indices with this array
  this->column_block_indices.reinit(col_sizes);
}



template <class MatrixType>
void
BlockMatrixBase<MatrixType>::prepare_add_operation()
{
  for (unsigned int row = 0; row < n_block_rows(); ++row)
    for (unsigned int col = 0; col < n_block_cols(); ++col)
      block(row, col).prepare_add();
}



template <class MatrixType>
void
BlockMatrixBase<MatrixType>::prepare_set_operation()
{
  for (unsigned int row = 0; row < n_block_rows(); ++row)
    for (unsigned int col = 0; col < n_block_cols(); ++col)
      block(row, col).prepare_set();
}

#endif // DOXYGEN


DEAL_II_NAMESPACE_CLOSE

#endif // dealii_block_matrix_base_h


