//include/deal.II-translator/lac/block_sparsity_pattern_0.txt
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

#ifndef dealii_block_sparsity_pattern_h
#define dealii_block_sparsity_pattern_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/base/table.h>

#include <deal.II/lac/block_indices.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/trilinos_sparsity_pattern.h>

DEAL_II_NAMESPACE_OPEN

// Forward declarations
#ifndef DOXYGEN
template <typename number>
class BlockSparseMatrix;
class BlockDynamicSparsityPattern;
#endif

/*!   @addtogroup Sparsity  
     * @{  

* 
*
*/


/**
 * 这是块版本的稀疏模式和动态稀疏模式类的基类。它的功能不多，只是管理一个疏散模式对象的数组，并将工作委托给它们。它与SparsityPattern和DynamicSparsityPattern的接口基本相同，只是将对其成员函数的调用转换为对成员疏散模式的各自成员函数的调用。
 * SparsityPattern和DynamicSparsityPattern类和这个类之间最大的区别是，大多数情况下，矩阵有不同的属性，你会想在构成矩阵的块上工作，而不是整个矩阵。你可以使用<tt>block(row,col)</tt>函数访问不同的块。
 * 注意：如果它的一个子对象的大小发生变化，这个对象不会被自动通知。因此，在你初始化子对象的大小后，你将不得不调用这个类的<tt>collect_sizes()</tt>函数!
 * 注意，当然，一个（块）行中的所有子矩阵必须有相同的行数，而一个（块）列中的所有子矩阵必须有相同的列数。
 * 一般来说，你不会想使用这个类，而是使用其中的一个派生类。
 * @todo  正确处理底层SparsityPattern的对角线元素的优化。
 * @see   @ref GlossBlockLA  "块（线性代数）"
 *
 *
 */
template <typename SparsityPatternType>
class BlockSparsityPatternBase : public Subscriptor
{
public:
  /**
   * 申报容器大小的类型。
   *
   */
  using size_type = types::global_dof_index;

  /**
   * 定义一个值，用来表示 @p
   * colnums数组中的某个值是未使用的，即不代表某个列号索引。
   * 这个值只是SparsityPattern类中各自数值的别名。
   *
   */
  static const size_type invalid_entry = SparsityPattern::invalid_entry;

  /**
   * 初始化矩阵为空，也就是没有分配内存。如果你想让这样的对象在其他类中作为成员变量，这很有用。你可以通过调用reinit()函数使该结构可用。
   *
   */
  BlockSparsityPatternBase();

  /**
   * 用给定的块的行和列的数量初始化矩阵。
   * 块本身仍然是空的，你必须在给它们分配尺寸后调用
   * collect_sizes()。
   *
   */
  BlockSparsityPatternBase(const size_type n_block_rows,
                           const size_type n_block_columns);

  /**
   * 复制构造函数。这个构造函数只允许在要复制的稀疏模式是空的情况下被调用，也就是说，目前没有分配的块。这与SparsityPattern的原因相同，详见那里。
   *
   */
  BlockSparsityPatternBase(const BlockSparsityPatternBase &bsp);

  /**
   * 解构器。
   *
   */
  ~BlockSparsityPatternBase() override;

  /**
   * 通过设置块的行数和列数，调整矩阵的大小。这将删除所有的块，并用未初始化的块代替它们，也就是说，那些块的大小还没有被设置。你必须通过调用块本身的reinit()函数来做到这一点。不要忘了在这之后对这个对象调用collect_sizes()。
   * 你必须自己设置块的大小的原因是，大小可能是变化的，每行的最大元素数可能是变化的，等等。在这里不复制SparsityPattern类的接口是比较简单的，而是让用户调用他们想要的任何函数。
   *
   */
  void
  reinit(const size_type n_block_rows, const size_type n_block_columns);

  /**
   * 复制操作符。对于这一点，与复制构造函数的情况相同：它被声明、定义并可以被调用，但后者只针对空对象。
   *
   */
  BlockSparsityPatternBase &
  operator=(const BlockSparsityPatternBase &);

  /**
   * 这个函数收集子对象的大小，并将其存储在内部数组中，以便能够将矩阵的全局索引转为子对象的索引。在你改变了子对象的大小之后，你必须*每次都调用这个函数。
   *
   */
  void
  collect_sizes();

  /**
   * 访问具有给定坐标的块。
   *
   */
  SparsityPatternType &
  block(const size_type row, const size_type column);


  /**
   * 访问具有给定坐标的块。对于常量对象的版本。
   *
   */
  const SparsityPatternType &
  block(const size_type row, const size_type column) const;

  /**
   * 授权访问描述行指数分布到各个块的对象。
   *
   */
  const BlockIndices &
  get_row_indices() const;

  /**
   * 允许访问描述各个块的列索引分布的对象。
   *
   */
  const BlockIndices &
  get_column_indices() const;

  /**
   * 这个函数压缩了这个对象所代表的稀疏结构。它只是为所有子对象调用
   * @p compress 。
   *
   */
  void
  compress();

  /**
   * 返回一列中的块的数量。
   *
   */
  size_type
  n_block_rows() const;

  /**
   * 返回一行中的块数。
   *
   */
  size_type
  n_block_cols() const;

  /**
   * 返回该对象是否为空。如果没有分配内存，它就是空的，这与两个维度都是零是一样的。这个函数只是对所有子矩阵的相应调用的串联。
   *
   */
  bool
  empty() const;

  /**
   * 返回每行的最大条目数。它返回每行的最大条目数，在一行的所有块上累积，以及所有行的最大条目数。
   *
   */
  size_type
  max_entries_per_row() const;

  /**
   * 向矩阵添加一个非零条目。这个函数只能用于非压缩的稀疏模式。
   * 如果该条目已经存在，则不会发生任何坏事。
   * 这个函数只是找出<tt>(i,j)</tt>属于哪个块，然后转发到该块。
   *
   */
  void
  add(const size_type i, const size_type j);

  /**
   * 向指定的矩阵行添加几个非零条目。
   * 这个函数只可用于非压缩的稀疏模式。
   * 如果有些条目已经存在，则不会发生什么坏事。
   * 这个函数只是找出迭代器范围内<tt>col</tt>的块<tt>(row,col)</tt>，然后转发到这些块。
   *
   */
  template <typename ForwardIterator>
  void
  add_entries(const size_type row,
              ForwardIterator begin,
              ForwardIterator end,
              const bool      indices_are_sorted = false);

  /**
   * 返回该矩阵的行数，等于图像空间的维度。它是子矩阵的（块-）行数之和。
   *
   */
  size_type
  n_rows() const;

  /**
   * 返回该矩阵的列数，相当于范围空间的维度。它是子矩阵的（块）列的列数之和。
   *
   */
  size_type
  n_cols() const;

  /**
   * 检查某个位置的值是否可能是非零。
   *
   */
  bool
  exists(const size_type i, const size_type j) const;

  /**
   * 某一行的条目数，与构成这一行的所有块相加。
   *
   */
  unsigned int
  row_length(const size_type row) const;

  /**
   * 返回这个矩阵的非零元素的数量。实际上，它返回的是稀疏模式中的条目数；如果任何一个条目碰巧是零，无论如何都会被计算在内。
   * 这个函数只有在矩阵结构被压缩的情况下才能被调用。否则就没有太大意义了。
   * 在当前情况下，它是子对象返回的数值之和。
   *
   */
  size_type
  n_nonzero_elements() const;

  /**
   * 打印矩阵的稀疏度。输出包括每行一行，格式为<tt>[i,j1,j2,j3,...]/tt>。<i>i</i>是行号，<i>jn</i>是该行中分配的列。
   *
   */
  void
  print(std::ostream &out) const;

  /**
   * 以<tt>gnuplot</tt>理解的格式打印矩阵的稀疏度，可以用来以图形方式绘制稀疏度模式。这与通常的稀疏度模式实现的功能相同，见
   * SparsityPatternBase::print_gnuplot().  。
   *
   */
  void
  print_gnuplot(std::ostream &out) const;

  /**
   * 以<tt>svg</tt>格式打印矩阵的稀疏度。这与通常的稀疏度模式实现的功能相同，见
   * SparsityPatternBase::print_svg().  。
   *
   */
  void
  print_svg(std::ostream &out) const;

  /**
   * @addtogroup Exceptions  
   * 
     * @{ 
   *
   */

  /**
   * 例外情况
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
   * 块的行数。
   *
   */
  size_type rows;

  /**
   * 块的列数。
   *
   */
  size_type columns;

  /**
   * 稀疏模式的阵列。
   *
   */
  Table<2,
        SmartPointer<SparsityPatternType,
                     BlockSparsityPatternBase<SparsityPatternType>>>
    sub_objects;

  /**
   * 储存和管理行索引到子对象索引的转换的对象。
   *
   */
  BlockIndices row_indices;

  /**
   * 储存和管理列索引到子对象索引的转换的对象。
   *
   */
  BlockIndices column_indices;

private:
  /**
   * 临时向量，用于计算在做集体添加或设置时写入各个块的元素。
   *
   */
  std::vector<size_type> counter_within_block;

  /**
   * 临时向量，用于在每个稀疏矩阵上将局部数据写入全局数据时，每个块上的列索引。
   *
   */
  std::vector<std::vector<size_type>> block_column_indices;

  // Make the block sparse matrix a friend, so that it can use our
  // #row_indices and #column_indices objects.
  template <typename number>
  friend class BlockSparseMatrix;
};



/**
 * 这个类扩展了基类，实现了一个可以被块状稀疏矩阵对象使用的稀疏模式阵列。它只增加了一些额外的成员函数，但主要的接口源于基类，更多信息请看那里。
 * 这个类是  @ref Sparsity  的 "静态 "
 * 类型的一个例子。
 *
 *
 */
class BlockSparsityPattern : public BlockSparsityPatternBase<SparsityPattern>
{
public:
  /**
   * 初始化矩阵为空，也就是没有分配内存。如果你想让这样的对象作为其他类中的成员变量，这很有用。你可以通过调用reinit()函数使该结构可用。
   *
   */
  BlockSparsityPattern() = default;

  /**
   * 用给定的块的行和列的数量初始化矩阵。
   * 块本身仍然是空的，你必须在给它们分配尺寸后调用
   * collect_sizes()。
   *
   */
  BlockSparsityPattern(const size_type n_rows, const size_type n_columns);

  /**
   * 转发到 BlockSparsityPatternBase::reinit(). 。
   *
   */
  void
  reinit(const size_type n_block_rows, const size_type n_block_columns);

  /**
   * 用两个BlockIndices初始化模式，用于矩阵行和列的块结构，以及一个行长向量。
   * 行长向量应该是由DoFTools产生的格式。
   * 另外，还有一个简化版本，每个内向量的长度为1。然后，相应的条目被用作最大的行长。
   * 对于对角线块，内部SparsityPattern用优化的对角线初始化，而对于非对角线块则不这样做。
   *
   */
  void
  reinit(const BlockIndices &                          row_indices,
         const BlockIndices &                          col_indices,
         const std::vector<std::vector<unsigned int>> &row_lengths);


  /**
   * 返回该结构是否被压缩，即所有子矩阵是否被压缩。
   *
   */
  bool
  is_compressed() const;

  /**
   * 确定该对象的内存消耗（以字节为单位）的估计值。
   *
   */
  std::size_t
  memory_consumption() const;

  /**
   * 从BlockDynamicSparsityPattern类型的对象中复制数据，即把这个对象的大小调整为给定参数的大小，并把每个子对象的内容复制过去。这个对象以前的内容会丢失。
   *
   */
  void
  copy_from(const BlockDynamicSparsityPattern &dsp);
};



/**
 * 这个类扩展了基类，实现了一个压缩稀疏模式的数组，可以用来初始化BlockSparsityPattern类型的对象。它没有增加额外的成员函数，而是作为一个
 * @p alias
 * 来介绍这个类的名称，而不要求用户指定基类的模板名称。关于这个类的接口的信息请参考基类。各个区块是基于DynamicSparsityPattern类的。
 * 这个类是  @ref Sparsity  的 "动态 "
 * 类型的一个例子。 <h3>Example</h3>
 * 这个类的用法与DynamicSparsityPattern非常相似，但由于使用块索引会引起一些额外的复杂情况，我们举一个简短的例子。
 * 在DoFHandler <tt>dof</tt>和AffineConstraints
 * <tt>constraints</tt>已经设置了系统元素后，我们必须计算每个矩阵块中的自由度。
 *
 *
 * @code
 * const std::vector<unsigned int> dofs_per_block =
 * DoFTools::count_dofs_per_fe_block(dof);
 * @endcode
 *
 * 现在，我们准备设置BlockDynamicSparsityPattern。
 *
 *
 * @code
 * BlockDynamicSparsityPattern dsp(fe.n_blocks(), fe.n_blocks());
 * for (unsigned int i = 0; i < fe.n_blocks(); ++i)
 * for (unsigned int j = 0; j < fe.n_blocks(); ++j)
 *   dsp.block(i, j).reinit(dofs_per_block[i], dofs_per_block[j]);
 * dsp.collect_sizes();
 * @endcode
 *
 * 它被填充，就像它是一个正常的模式一样
 *
 *
 * @code
 * DoFTools::make_sparsity_pattern(dof, dsp);
 * constraints.condense(dsp);
 * @endcode
 *
 * 最后，它被复制到一个正常的BlockSparsityPattern中供以后使用。
 *
 *
 * @code
 * BlockSparsityPattern sparsity;
 * sparsity.copy_from(dsp);
 * @endcode
 *
 *
 */

class BlockDynamicSparsityPattern
  : public BlockSparsityPatternBase<DynamicSparsityPattern>
{
public:
  /**
   * 初始化矩阵为空，也就是没有分配内存。如果你想把这样的对象作为其他类中的成员变量，这是很有用的。你可以通过调用reinit()函数使该结构可用。
   *
   */
  BlockDynamicSparsityPattern() = default;

  /**
   * 用给定的块的行和列的数量初始化矩阵。
   * 块本身仍然是空的，你必须在给它们分配尺寸后调用
   * collect_sizes()。
   *
   */
  BlockDynamicSparsityPattern(const size_type n_rows,
                              const size_type n_columns);

  /**
   * 用两个BlockIndices初始化模式，用于矩阵行和列的块结构。这个函数相当于用两个索引向量的长度调用前面的构造函数，然后输入索引值。
   *
   */
  BlockDynamicSparsityPattern(const std::vector<size_type> &row_block_sizes,
                              const std::vector<size_type> &col_block_sizes);

  /**
   * 用对称块初始化模式。向量中IndexSets的数量决定了块的行和列的数量。每个块的大小由各自IndexSet的size()决定。
   * 每个块只存储由IndexSet中的值给出的行，这对分布式内存并行计算很有用，通常对应于本地拥有的DoF。
   *
   */
  BlockDynamicSparsityPattern(const std::vector<IndexSet> &partitioning);

  /**
   * 用两个BlockIndices初始化模式，用于矩阵行和列的块结构。
   *
   */
  BlockDynamicSparsityPattern(const BlockIndices &row_indices,
                              const BlockIndices &col_indices);


  /**
   * 将模式调整为由参数定义的尺寸的矩阵的张量乘积。
   * 矩阵将有与两个参数中的条目一样多的块行和块列。在位置（<i>i,j</i>）的块将有<tt>row_block_sizes[i]</tt>乘以<tt>col_block_sizes[j]</tt>的尺寸。
   *
   */
  void
  reinit(const std::vector<size_type> &row_block_sizes,
         const std::vector<size_type> &col_block_sizes);

  /**
   * 用每个IndexSet的size()决定的对称块来调整模式的大小。详见接受一个IndexSets向量的构造函数。
   *
   */
  void
  reinit(const std::vector<IndexSet> &partitioning);

  /**
   * 将矩阵调整为由参数定义的尺寸的矩阵的张量积。两个BlockIndices对象必须被初始化，之后的稀疏模式将具有相同的块结构。
   *
   */
  void
  reinit(const BlockIndices &row_indices, const BlockIndices &col_indices);

  /**
   * 访问列号字段。返回第 @p index 行中第 @p row.
   * 个条目的列号。
   *
   */
  size_type
  column_number(const size_type row, const unsigned int index) const;

  /**
   * 也允许使用基类的 reinit 函数。
   *
   */
  using BlockSparsityPatternBase<DynamicSparsityPattern>::reinit;
};

 /*@}*/ 


#ifdef DEAL_II_WITH_TRILINOS


namespace TrilinosWrappers
{
  /*!   @addtogroup  TrilinosWrappers  
     * @{   
*
*/

  /**
   * 这个类扩展了基类，实现了一个Trilinos稀疏模式的数组，可以用来初始化Trilinos块状稀疏矩阵，可以分布在不同的处理器中。除了建立在 TrilinosWrappers::SparsityPattern 而不是 dealii::SparsityPattern. 的基础上，它的使用方式与 dealii::BlockSparsityPattern 相同。该类具有 @ref Sparsity 的 "动态 "
   * 类型的属性（即如果分配的元素太少，它可以扩展内存），但其他方面更像基本的deal.II
   * SparsityPattern（即在使用该模式之前需要调用compress()方法）。
   * 这个类在  step-32  中使用。
   *
   */
  class BlockSparsityPattern
    : public dealii::BlockSparsityPatternBase<SparsityPattern>
  {
  public:
    /**
     * 初始化矩阵为空，也就是没有分配内存。如果你想让这样的对象作为其他类的成员变量，这很有用。
     * 你可以通过调用reinit()函数使该结构可用。
     *
     */
    BlockSparsityPattern() = default;

    /**
     * 用给定的块的行和列的数量初始化矩阵。
     * 块本身仍然是空的，你必须在给它们分配尺寸后调用collect_sizes()。
     *
     */
    BlockSparsityPattern(const size_type n_rows, const size_type n_columns);

    /**
     * 用两个BlockIndices初始化模式，用于矩阵行和列的块结构。这个函数相当于用两个索引向量的长度调用前面的构造函数，然后输入索引值。
     *
     */
    BlockSparsityPattern(const std::vector<size_type> &row_block_sizes,
                         const std::vector<size_type> &col_block_sizes);

    /**
     * 用一个索引集数组初始化模式，该数组指定了矩阵的行和列（所以最终的矩阵将是一个方形矩阵），其中IndexSets的size()指定了块的大小，每个IndexSet中的值表示每个块中要保存的行。
     *
     */
    BlockSparsityPattern(const std::vector<IndexSet> &parallel_partitioning,
                         const MPI_Comm &communicator = MPI_COMM_WORLD);

    /**
     * 用两个指定矩阵的行和列的索引集数组初始化该模式，其中IndexSets的size()指定了块的大小，每个IndexSet中的值表示每个块中要保存的行。额外的索引集
     * writable_rows 用来设置所有我们允许本地写入的行。
     * 这个构造函数用于创建允许多个线程同时向矩阵中写入的矩阵（当然是向不同的行写入），更多细节请参见方法
     * TrilinosWrappers::SparsityPattern::reinit
     * 带有三个索引集参数的方法。
     *
     */
    BlockSparsityPattern(
      const std::vector<IndexSet> &row_parallel_partitioning,
      const std::vector<IndexSet> &column_parallel_partitioning,
      const std::vector<IndexSet> &writeable_rows,
      const MPI_Comm &             communicator = MPI_COMM_WORLD);

    /**
     * 调整矩阵的大小，使其成为一个由参数定义的尺寸的矩阵的张量积。
     * 矩阵将有与两个参数中的条目一样多的块行和块列。在位置（<i>i,j</i>）的块将有<tt>row_block_sizes[i]</tt>乘以<tt>col_block_sizes[j]</tt>的尺寸。
     *
     */
    void
    reinit(const std::vector<size_type> &row_block_sizes,
           const std::vector<size_type> &col_block_sizes);

    /**
     * 将矩阵调整为矩阵的平方张量积。详见接受一个IndexSets向量的构造函数。
     *
     */
    void
    reinit(const std::vector<IndexSet> &parallel_partitioning,
           const MPI_Comm &             communicator = MPI_COMM_WORLD);

    /**
     * 将矩阵调整为一个矩形块状矩阵。这种方法允许行和列是不同的，在外部块结构和块内都是如此。
     *
     */
    void
    reinit(const std::vector<IndexSet> &row_parallel_partitioning,
           const std::vector<IndexSet> &column_parallel_partitioning,
           const MPI_Comm &             communicator = MPI_COMM_WORLD);

    /**
     * 将矩阵调整为一个矩形块矩阵，进一步明确指定每个块中的可写行。这个方法用于创建允许几个线程同时向矩阵中写入的矩阵（当然是向不同的行写入），更多细节请参见方法
     * TrilinosWrappers::SparsityPattern::reinit
     * 带有三个索引集参数的方法。
     *
     */
    void
    reinit(const std::vector<IndexSet> &row_parallel_partitioning,
           const std::vector<IndexSet> &column_parallel_partitioning,
           const std::vector<IndexSet> &writeable_rows,
           const MPI_Comm &             communicator = MPI_COMM_WORLD);

    /**
     * 也允许使用基类的reinit函数。
     *
     */
    using BlockSparsityPatternBase<SparsityPattern>::reinit;
  };

   /*@}*/ 

}  /* namespace TrilinosWrappers */ 

#endif

 /*--------------------- Template functions ----------------------------------*/ 



template <typename SparsityPatternType>
inline SparsityPatternType &
BlockSparsityPatternBase<SparsityPatternType>::block(const size_type row,
                                                     const size_type column)
{
  AssertIndexRange(row, rows);
  AssertIndexRange(column, columns);
  return *sub_objects[row][column];
}



template <typename SparsityPatternType>
inline const SparsityPatternType &
BlockSparsityPatternBase<SparsityPatternType>::block(
  const size_type row,
  const size_type column) const
{
  AssertIndexRange(row, rows);
  AssertIndexRange(column, columns);
  return *sub_objects[row][column];
}



template <typename SparsityPatternType>
inline const BlockIndices &
BlockSparsityPatternBase<SparsityPatternType>::get_row_indices() const
{
  return row_indices;
}



template <typename SparsityPatternType>
inline const BlockIndices &
BlockSparsityPatternBase<SparsityPatternType>::get_column_indices() const
{
  return column_indices;
}



template <typename SparsityPatternType>
inline void
BlockSparsityPatternBase<SparsityPatternType>::add(const size_type i,
                                                   const size_type j)
{
  // if you get an error here, are
  // you sure you called
  // <tt>collect_sizes()</tt> before?
  const std::pair<size_type, size_type> row_index =
                                          row_indices.global_to_local(i),
                                        col_index =
                                          column_indices.global_to_local(j);
  sub_objects[row_index.first][col_index.first]->add(row_index.second,
                                                     col_index.second);
}



template <typename SparsityPatternType>
template <typename ForwardIterator>
void
BlockSparsityPatternBase<SparsityPatternType>::add_entries(
  const size_type row,
  ForwardIterator begin,
  ForwardIterator end,
  const bool      indices_are_sorted)
{
  // Resize scratch arrays
  if (block_column_indices.size() < this->n_block_cols())
    {
      block_column_indices.resize(this->n_block_cols());
      counter_within_block.resize(this->n_block_cols());
    }

  const size_type n_cols = static_cast<size_type>(end - begin);

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
  if (block_column_indices[0].size() < n_cols)
    for (size_type i = 0; i < this->n_block_cols(); ++i)
      block_column_indices[i].resize(n_cols);

  // Reset the number of added elements
  // in each block to zero.
  for (size_type i = 0; i < this->n_block_cols(); ++i)
    counter_within_block[i] = 0;

  // Go through the column indices to
  // find out which portions of the
  // values should be set in which
  // block of the matrix. We need to
  // touch all the data, since we can't
  // be sure that the data of one block
  // is stored contiguously (in fact,
  // indices will be intermixed when it
  // comes from an element matrix).
  for (ForwardIterator it = begin; it != end; ++it)
    {
      const size_type col = *it;

      const std::pair<size_type, size_type> col_index =
        this->column_indices.global_to_local(col);

      const size_type local_index = counter_within_block[col_index.first]++;

      block_column_indices[col_index.first][local_index] = col_index.second;
    }

#ifdef DEBUG
  // If in debug mode, do a check whether
  // the right length has been obtained.
  size_type length = 0;
  for (size_type i = 0; i < this->n_block_cols(); ++i)
    length += counter_within_block[i];
  Assert(length == n_cols, ExcInternalError());
#endif

  // Now we found out about where the
  // individual columns should start and
  // where we should start reading out
  // data. Now let's write the data into
  // the individual blocks!
  const std::pair<size_type, size_type> row_index =
    this->row_indices.global_to_local(row);
  for (size_type block_col = 0; block_col < n_block_cols(); ++block_col)
    {
      if (counter_within_block[block_col] == 0)
        continue;
      sub_objects[row_index.first][block_col]->add_entries(
        row_index.second,
        block_column_indices[block_col].begin(),
        block_column_indices[block_col].begin() +
          counter_within_block[block_col],
        indices_are_sorted);
    }
}



template <typename SparsityPatternType>
inline bool
BlockSparsityPatternBase<SparsityPatternType>::exists(const size_type i,
                                                      const size_type j) const
{
  // if you get an error here, are
  // you sure you called
  // <tt>collect_sizes()</tt> before?
  const std::pair<size_type, size_type> row_index =
                                          row_indices.global_to_local(i),
                                        col_index =
                                          column_indices.global_to_local(j);
  return sub_objects[row_index.first][col_index.first]->exists(
    row_index.second, col_index.second);
}



template <typename SparsityPatternType>
inline unsigned int
BlockSparsityPatternBase<SparsityPatternType>::row_length(
  const size_type row) const
{
  const std::pair<size_type, size_type> row_index =
    row_indices.global_to_local(row);

  unsigned int c = 0;

  for (size_type b = 0; b < rows; ++b)
    c += sub_objects[row_index.first][b]->row_length(row_index.second);

  return c;
}



template <typename SparsityPatternType>
inline typename BlockSparsityPatternBase<SparsityPatternType>::size_type
BlockSparsityPatternBase<SparsityPatternType>::n_block_cols() const
{
  return columns;
}



template <typename SparsityPatternType>
inline typename BlockSparsityPatternBase<SparsityPatternType>::size_type
BlockSparsityPatternBase<SparsityPatternType>::n_block_rows() const
{
  return rows;
}


inline BlockDynamicSparsityPattern::size_type
BlockDynamicSparsityPattern::column_number(const size_type    row,
                                           const unsigned int index) const
{
  // .first= ith block, .second = jth row in that block
  const std::pair<size_type, size_type> row_index =
    row_indices.global_to_local(row);

  AssertIndexRange(index, row_length(row));

  size_type c             = 0;
  size_type block_columns = 0; // sum of n_cols for all blocks to the left
  for (unsigned int b = 0; b < columns; ++b)
    {
      unsigned int rowlen =
        sub_objects[row_index.first][b]->row_length(row_index.second);
      if (index < c + rowlen)
        return block_columns +
               sub_objects[row_index.first][b]->column_number(row_index.second,
                                                              index - c);
      c += rowlen;
      block_columns += sub_objects[row_index.first][b]->n_cols();
    }

  Assert(false, ExcInternalError());
  return 0;
}


inline void
BlockSparsityPattern::reinit(const size_type n_block_rows,
                             const size_type n_block_columns)
{
  BlockSparsityPatternBase<SparsityPattern>::reinit(n_block_rows,
                                                    n_block_columns);
}


DEAL_II_NAMESPACE_CLOSE

#endif


