//include/deal.II-translator/lac/sparsity_pattern_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2021 by the deal.II authors
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

#ifndef dealii_sparsity_pattern_h
#define dealii_sparsity_pattern_h


#include <deal.II/base/config.h>

#include <deal.II/base/array_view.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/linear_index_iterator.h>
#include <deal.II/base/subscriptor.h>

// boost::serialization::make_array used to be in array.hpp, but was
// moved to a different file in BOOST 1.64
#include <boost/version.hpp>
#if BOOST_VERSION >= 106400
#  include <boost/serialization/array_wrapper.hpp>
#else
#  include <boost/serialization/array.hpp>
#endif
#include <boost/serialization/split_member.hpp>

#include <algorithm>
#include <iostream>
#include <memory>
#include <vector>

DEAL_II_NAMESPACE_OPEN

// Forward declarations
#ifndef DOXYGEN
class SparsityPattern;
class SparsityPatternBase;
class DynamicSparsityPattern;
class ChunkSparsityPattern;
template <typename number>
class FullMatrix;
template <typename number>
class SparseMatrix;
template <typename number>
class SparseLUDecomposition;
template <typename number>
class SparseILU;

namespace ChunkSparsityPatternIterators
{
  class Accessor;
}
#endif

/*!   @addtogroup Sparsity  
     * @{  

* 
*
*/

namespace internals
{
  namespace SparsityPatternTools
  {
    /**
     * 申报容器大小的类型。
     *
     */
    using size_type = types::global_dof_index;

    /**
     * 在copy_from()函数中，如果内部迭代器类型指向普通的无符号整数，则用辅助函数从一个被解读的迭代器中获取列索引。
     *
     */
    size_type
    get_column_index_from_iterator(const size_type i);

    /**
     * 在copy_from()函数中，如果内部迭代器类型指向无符号整数和其他值的对，则用辅助函数从一个被解除引用的迭代器中获取列索引。
     *
     */
    template <typename value>
    size_type
    get_column_index_from_iterator(const std::pair<size_type, value> &i);

    /**
     * 同样，但有时需要用于某些类型的容器，使对中的第一个元素成为常数（如
     * <tt>std::map</tt>).  ）。
     *
     */
    template <typename value>
    size_type
    get_column_index_from_iterator(const std::pair<const size_type, value> &i);

  } // namespace SparsityPatternTools
} // namespace internals


/**
 * SparsityPattern类型的对象上的迭代器。
 *
 *
 */
namespace SparsityPatternIterators
{
  // forward declaration
  class Iterator;

  /**
   * 为容器大小声明类型。
   *
   */
  using size_type = types::global_dof_index;

  /**
   * 进入稀疏性模式的迭代器的访问器类。这个类也是进入稀疏矩阵的常量和非常量访问器类的基类。
   * 请注意，这个类只允许对元素进行读取访问，提供它们的行号和列号（或者是完整的稀疏模式中的索引）。它不允许修改稀疏模式本身。
   *
   */
  class Accessor
  {
  public:
    /**
     * SparsityPattern的大小类型。
     *
     */
    using size_type = SparsityPatternIterators::size_type;

    /**
     * 构造函数。
     *
     */
    Accessor(const SparsityPatternBase *matrix, const std::size_t linear_index);

    /**
     * 构造器。为给定的稀疏性模式构造结束访问器。
     *
     */
    Accessor(const SparsityPatternBase *matrix);

    /**
     * 默认构造器创建一个假访问器。这个构造函数在这里只是为了能够在STL容器中存储访问器，例如
     * `std::vector`.  。
     *
     */
    Accessor();

    /**
     * 这个对象所代表的元素的行号。这个函数只能为is_valid_entry()为真的条目调用。
     *
     */
    size_type
    row() const;

    /**
     * 这个对象所代表的元素在当前行中的索引。
     * 这个函数只能对is_valid_entry()为true的条目进行调用。
     *
     */
    size_type
    index() const;

    /**
     * 这个函数返回当前迭代器指向的整个稀疏模式中的第多少个条目。虽然疏散模式中的条目顺序通常并不重要，但这个函数允许使用线性索引来索引疏散模式的条目。
     * 这个函数只能为is_valid_entry()为真的条目调用。
     *
     */
    size_type
    global_index() const;

    /**
     * 这个对象所代表的元素的列号。这个函数只能为is_valid_entry()为真的条目调用。
     *
     */
    size_type
    column() const;

    /**
     * 返回这个迭代器所指向的稀疏模式条目是否有效。注意，在压缩稀疏模式后，所有条目都是有效的。然而，在压缩之前，稀疏模式分配了一些内存来使用，同时仍在增加新的非零条目；如果你在稀疏模式生命周期的这个阶段创建迭代器，你将迭代那些无效的元素。
     * 如果是这样的话，那么这个函数将返回false。
     *
     */
    bool
    is_valid_entry() const;

    /**
     * 比较。真，如果两个迭代器都指向同一个矩阵位置。
     *
     */
    bool
    operator==(const Accessor &) const;

    /**
     * 比较运算符。如果第一行数字较小，或者行数字相等且第一个索引较小，则结果为真。
     * 这个函数只有在两个迭代器都指向同一个稀疏模式时才有效。
     *
     */
    bool
    operator<(const Accessor &) const;

  protected:
    DeclExceptionMsg(DummyAccessor,
                     "The instance of this class was initialized"
                     " without SparsityPattern object, which"
                     " means that it is a dummy accessor that can"
                     " not do any operations.");

    /**
     * 我们所操作的稀疏模式被访问。
     *
     */
    const SparsityPatternBase *container;

    /**
     * 全局稀疏度模式中的索引。这个索引代表了迭代器/访问器在SparsityPattern类的数组中所指向的位置，该数组存储了列号。它也是稀疏矩阵的数值数组内的索引，该数组存储了该站点的相应数值。
     *
     */
    std::size_t linear_index;

    /**
     * 将访问器移动到矩阵中的下一个非零条目。
     *
     */
    void
    advance();

    // Grant access to iterator class.
    friend class LinearIndexIterator<Iterator, Accessor>;

    // Grant access to accessor class of ChunkSparsityPattern.
    friend class ChunkSparsityPatternIterators::Accessor;
  };



  /**
   * 一个迭代器类，用于在稀疏模式的元素上行走。
   * 这些迭代器的典型用途是迭代稀疏模式的元素（或者，因为它们也是迭代相关矩阵的元素的基础，所以是迭代稀疏矩阵的元素），或者迭代单个行的元素。不能保证行的元素实际上是按照列数单调增加的顺序来遍历的。更多信息请参见SparsityPattern类的文档。
   * @note
   * 该类直接对SparsityPatternBase类的内部数据结构进行操作。因此，有些操作很便宜，有些则不然。特别是，访问指向的稀疏模式条目的列索引是很便宜的。另一方面，访问行索引是很昂贵的（对于一个有
   * $N$ 行的矩阵，这需要 $O(\log(N))$
   * 次操作）。因此，当你设计使用这些迭代器的算法时，通常的做法是不一次性循环疏散模式的<i>all</i>个元素，而是在所有行上有一个外循环，在这个循环中迭代这个行的元素。这样，你只需要解除对迭代器的引用以获得列索引，而通过使用循环索引可以避免对行索引的（昂贵）查找。
   *
   */
  class Iterator : public LinearIndexIterator<Iterator, Accessor>
  {
  public:
    /**
     * 尺寸类型。
     *
     */
    using size_type = types::global_dof_index;

    /**
     * 存储指针的类型。
     *
     */
    using container_pointer_type = SparsityPatternBase *;

    /**
     * 构造函数。为给定的全局索引（即从第2行开始计算的给定元素的索引）创建一个进入稀疏模式
     * @p sp 的迭代器。
     *
     */
    Iterator(const SparsityPatternBase *sp, const std::size_t linear_index);

    /**
     * 构造函数。为给定的访问器创建一个进入稀疏模式 @p
     * sp 的迭代器。
     *
     */
    Iterator(const Accessor &accessor);
  };
} // namespace SparsityPatternIterators



/**
 * 一个可以存储矩阵中哪些元素是非零的类（或者，事实上，<i>may</i>是非零的），我们必须为其分配内存以存储其值。这个类是 "静态 "
 * 类型的稀疏格局的一个例子（见 @ref Sparsity  ）。它使用<a
 * href="https://en.wikipedia.org/wiki/Sparse_matrix">compressed row storage
 * (CSR)</a>格式来存储数据，并被用作派生的SparsityPattern类和SparseMatrix类的基础。
 * SparsityPatternBase的元素，对应于SparseMatrix对象可以存储非零条目的地方，被逐行存储。每行中的非零元素的排序（即增加列索引的顺序）取决于派生类。
 *
 *
 */
class SparsityPatternBase : public Subscriptor
{
public:
  /**
   * 声明容器大小的类型。
   *
   */
  using size_type = types::global_dof_index;

  /**
   * Typedef一个迭代器类，允许在一个稀疏模式的所有非零元素上行走。
   *
   */
  using const_iterator = SparsityPatternIterators::Iterator;

  /**
   * Typedef一个迭代器类，允许在一个稀疏模式的所有非零元素上行走。
   * 由于该迭代器不允许修改稀疏模式，该类型与 @p
   * const_iterator. 的类型相同。
   *
   */
  using iterator = SparsityPatternIterators::Iterator;

  /**
   * 定义一个值，用来表示#colnums数组中的某个值是未使用的，即不代表某个列号索引。
   * 具有这种无效值的索引被用来使用add()成员函数向稀疏模式插入新条目，并在调用compress()时被删除。
   * 你不应该认为这里声明的变量有一定的价值。这里给出的初始化只是为了让编译器进行一些优化，但变量的实际值可能会随着时间的推移而改变。
   *
   */
  static const size_type invalid_entry = numbers::invalid_size_type;

  /**
   * @name  构造和初始化 构造器，析构器，初始化、复制和填充对象的函数。
   *
   */
  // @{
  /**
   * 初始化矩阵是空的，也就是没有分配内存。如果你想让这样的对象作为其他类中的成员变量，这是很有用的。你可以通过调用reinit()函数使该结构可用。
   *
   */
  SparsityPatternBase();

  /**
   * 解构器。
   *
   */
  ~SparsityPatternBase() override = default;

  /**
   * 为一个新的矩阵重新分配内存并设置数据结构，该矩阵有
   * @p m 行和 @p n 列，每行最多有 @p max_per_row 个非零条目。
   * 这个函数只是将其操作映射到另一个reinit()函数。
   *
   */
  void
  reinit(const size_type m, const size_type n, const unsigned int max_per_row);

  /**
   * 为大小为 @p m 乘以 @p n.
   * 的矩阵重新分配内存，每行的条目数取自数组 @p
   * row_lengths ，它必须给出每行的这个数字 $i=1\ldots m$  。
   * 如果<tt>m*n==0</tt>所有的内存被释放，导致对象的完全重新初始化。如果是非零，只有当新的大小扩展到旧的大小时，才会分配新的内存。这样做是为了节省时间和避免堆的碎片化。
   *
   */
  void
  reinit(const size_type                  m,
         const size_type                  n,
         const std::vector<unsigned int> &row_lengths);

  /**
   * 和上面一样，但用ArrayView参数代替。
   * 派生类负责实现这个函数。
   *
   */
  virtual void
  reinit(const size_type                      m,
         const size_type                      n,
         const ArrayView<const unsigned int> &row_lengths) = 0;

  /**
   * 通过添加转置对象的稀疏模式使稀疏模式对称化。
   * 如果稀疏度模式不代表二次矩阵，这个函数会抛出一个异常。
   *
   */
  void
  symmetrize();

  /**
   * 给矩阵添加一个非零条目。
   * 这个函数只能用于非压缩的稀疏度模式。
   * 如果该条目已经存在，则不会发生任何坏事。
   *
   */
  void
  add(const size_type i, const size_type j);

  // @}

  /**
   * @name  迭代器
   *
   */
  // @{

  /**
   * 迭代器从矩阵的第一个条目开始。由此产生的迭代器可以用来走过稀疏模式的所有非零条目。
   * 访问元素的顺序取决于派生类实现的存储方案。
   *
   */
  iterator
  begin() const;

  /**
   * 最终迭代器。
   *
   */
  iterator
  end() const;

  /**
   * 迭代器从行<tt>r</tt>的第一个条目开始。
   * 注意，如果给定的行是空的，即不包含任何非零条目，那么这个函数返回的迭代器就等于<tt>end(r)</tt>。还要注意的是，在这种情况下，迭代器可能不能被解除引用。
   * 元素被访问的顺序取决于派生类实现的存储方案。
   *
   */
  iterator
  begin(const size_type r) const;

  /**
   * 行<tt>r</tt>的最终迭代器。它指向超过行 @p r,
   * 末尾的第一个元素或超过整个稀疏模式的末尾。
   * 请注意，结束迭代器不一定是可被解除引用的。特别是如果它是一个矩阵的最后一行的结束迭代器，情况更是如此。
   *
   */
  iterator
  end(const size_type r) const;


  // @}

  /**
   * @name  查询信息
   *
   */
  // @{

  /**
   * 测试两个SparsityPatterns是否相等。
   *
   */
  bool
  operator==(const SparsityPatternBase &) const;

  /**
   * 返回该对象是否为空。如果没有分配内存，它就是空的，这与两个维度都是零是一样的。
   *
   */
  bool
  empty() const;

  /**
   * 检查在某一位置的值是否可能是非零。
   *
   */
  bool
  exists(const size_type i, const size_type j) const;

  /**
   * 返回每行的最大条目数。在压缩之前，这等于给构造函数的数字，而在压缩之后，它等于用户实际分配的最大条目数。
   *
   */
  size_type
  max_entries_per_row() const;

  /**
   * 计算这个结构所代表的矩阵的带宽。该带宽是 $|i-j|$
   * 的最大值，其中索引对 $(i,j)$
   * 代表矩阵的一个非零条目。因此， $n\times m$
   * 矩阵的最大带宽为 $\max\{n-1,m-1\}$
   * ，对角线矩阵的带宽为0，如果带宽为 $q$
   * ，则每行最多有 $2*q+1$
   * 个条目。返回的数量有时在文献中被称为 "半带宽"。
   *
   */
  size_type
  bandwidth() const;

  /**
   * 返回这个矩阵的非零元素的数量。实际上，它返回的是稀疏模式中的条目数；如果任何一个条目碰巧是零，无论如何都会被计算在内。
   * 这个函数只有在矩阵结构被压缩的情况下才能被调用。否则就没有太大意义了。
   *
   */
  std::size_t
  n_nonzero_elements() const;

  /**
   * 返回该结构是否被压缩。
   *
   */
  bool
  is_compressed() const;

  /**
   * 返回该矩阵的行数，相当于图像空间的维度。
   *
   */
  size_type
  n_rows() const;

  /**
   * 返回该矩阵的列数，相当于范围空间的维度。
   *
   */
  size_type
  n_cols() const;

  /**
   * 特定行中的条目数。
   *
   */
  unsigned int
  row_length(const size_type row) const;

  /**
   * 确定此对象的内存消耗（以字节为单位）的估计值。见MemoryConsumption。
   *
   */
  std::size_t
  memory_consumption() const;

  // @}

  /**
   * @name  访问条目
   *
   */
  // @{

  /**
   * 访问列号字段。
   * 返回<tt>index</tt>中的第<tt>row</tt>条目的列号。注意，如果对角线元素被优化，每行的第一个元素就是对角线元素，即<tt>column_number(row,0)==row</tt>。
   * 如果稀疏模式已经被压缩，那么（除了对角线元素），条目按列排序，即：<tt>column_number(row,i)</tt>
   * <tt></tt> <tt>column_number(row,i+1)</tt>。
   *
   */
  size_type
  column_number(const size_type row, const unsigned int index) const;

  /**
   * 一个全局矩阵条目在其行中的索引。
   * 这个函数类似于operator()，但它计算的索引不是关于总域的，而只是关于行的<tt>j</tt>。
   *
   */
  size_type
  row_position(const size_type i, const size_type j) const;

  /**
   * 这是operator()()的逆向操作：给定一个全局索引，找出它所属的矩阵条目的行和列。返回值是由行和列索引组成的一对。
   * 这个函数只有在稀疏模式是封闭的情况下才能被调用。那么全局索引必须在0和n_nonzero_elements()之间。
   * 如果<tt>N</tt>是这个矩阵的行数，那么这个函数的复杂性是<i>log(N)</i>。
   *
   */
  std::pair<size_type, size_type>
  matrix_position(const std::size_t global_index) const;

  // @}

  /**
   * @name  输入/输出
   *
   */
  // @{

  /**
   * 打印矩阵的稀疏度。输出包括每行一行，格式为<tt>[i,j1,j2,j3,...]/tt>。<i>i</i>是行号，<i>jn</i>是这一行中分配的列。
   *
   */
  void
  print(std::ostream &out) const;

  /**
   * 以<tt>gnuplot</tt>能理解的格式打印矩阵的稀疏度，该格式可用于以图形方式绘制稀疏度模式。该格式由成对的<tt>i
   * j</tt>非零元素组成，每个元素代表该矩阵的一个条目，输出文件中每行一个。指数从零开始计算，和平常一样。由于稀疏模式的打印方式与矩阵的显示方式相同，我们打印的是列索引的负数，这意味着<tt>(0,0)</tt>元素位于左上角而不是左下角。
   * 在gnuplot中通过将数据样式设置为点或点来打印稀疏模式，并使用<tt>plot</tt>命令。
   *
   */
  void
  print_gnuplot(std::ostream &out) const;

  /**
   * 在一个.svg文件中打印出矩阵的稀疏度，可以在网络浏览器中打开。该.svg文件包含与矩阵中的条目相对应的方块。矩阵中包含非零值的条目对应的是一个红色的方块，而矩阵中零值的条目对应的是一个白色方块。
   *
   */
  void
  print_svg(std::ostream &out) const;

  /**
   * 使用[BOOST序列化库](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html)将此对象的数据写入一个流中，以便进行序列化。
   *
   */
  template <class Archive>
  void
  save(Archive &ar, const unsigned int version) const;

  /**
   * 使用[BOOST序列化库](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html)从一个流中读取此对象的数据，以达到序列化的目的。
   *
   */
  template <class Archive>
  void
  load(Archive &ar, const unsigned int version);

#ifdef DOXYGEN
  /**
   * 使用[BOOST序列化库](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html)从流中写入和读取此对象的数据，以达到序列化的目的。
   *
   */
  template <class Archive>
  void
  serialize(Archive &archive, const unsigned int version);
#else
  // This macro defines the serialize() method that is compatible with
  // the templated save() and load() method that have been implemented.
  BOOST_SERIALIZATION_SPLIT_MEMBER()
#endif

  // @}

  /**
   * @addtogroup  Exceptions  
     * @{ 
   *
   */

  /**
   * 只有在设置了SparsityPattern并且调用了compress()之后才允许进行该操作。
   *
   */
  DeclExceptionMsg(
    ExcNotCompressed,
    "The operation you attempted is only allowed after the SparsityPattern "
    "has been set up and compress() was called.");

  /**
   * 你试图向某一行添加一个元素，但没有剩余的空间。
   *
   */
  DeclException2(ExcNotEnoughSpace,
                 int,
                 int,
                 << "Upon entering a new entry to row " << arg1
                 << ": there was no free entry any more. " << std::endl
                 << "(Maximum number of entries for this row: " << arg2
                 << "; maybe the matrix is already compressed?)");

  /**
   * 这个操作改变了SparsityPattern的结构，在调用了compress()之后就不可能了。
   *
   */
  DeclExceptionMsg(
    ExcMatrixIsCompressed,
    "The operation you attempted changes the structure of the SparsityPattern "
    "and is not possible after compress() has been called.");

  // @}


protected:
  /**
   * 可以存储在#rowstart数组中的最大行数。
   * 因为只有在当前数组太小的情况下才会重新分配该数组，而不是当该矩阵结构的大小缩小时，#max_dim可能大于#rows，在这种情况下，#rowstart的元素比使用的要多。
   *
   */
  size_type max_dim;

  /**
   * 该稀疏结构应表示的行数。
   *
   */
  size_type rows;

  /**
   * 该稀疏结构应代表的列数。
   *
   */
  size_type cols;

  /**
   * 实际分配的数组#colnums的大小。这里，与#rowstart数组的情况相同，即它可能大于数组的实际使用部分。
   *
   */
  std::size_t max_vec_len;

  /**
   * 每行的最大元素数。这个值被设置为给reinit()函数（或构造函数）的值，或者在调用更灵活的构造函数或reinit版本时，设置为从向量中计算出的最大行长度。在调用compress()后，它的值或多或少没有意义。
   *
   */
  unsigned int max_row_length;

  /**
   * 数组，为每一行保存属于该行的#colnums中的第一个元素。请注意，数组的大小比行数大一个，因为最后一个元素是用于<tt>row</tt>=#rows，即超过最后使用的行。#rowstart[#rows]}的值等于#colnums中超过终点的元素的索引；这样，我们就能够写出类似<tt>for
   * (i=rowstart[k]; i<rowstart[k+1];
   * ++i)</tt>的循环，也是为了最后一行。
   * 请注意，分配的内存的实际大小可能比使用的区域要大。被分配的实际元素数被存储在#max_dim中。
   *
   */
  std::unique_ptr<std::size_t[]> rowstart;

  /**
   * 列号的数组。在这个数组中，我们为每个非零元素存储其列号。第<i>r</i>行的元素的列号被存储在#rowstart[<i>r</i>]...#rowstart[<i>r+1</i>]的索引范围内。因此要找出一个给定的元素（<i>r,c</i>）是否存在，我们必须检查列号<i>c</i>是否存在于这个数组的上述范围内。如果它存在，比如在这个数组中的<i>p</i>位置，那么稀疏矩阵中相应元素的值也将在该类值数组的<i>p</i>位置。
   * 在开始时，这个数组的所有元素都被设置为 @p 。
   *
   * - 表示无效的（未使用的）列号（不过如果要求优化存储，对角线元素会被预设）。现在，如果非零元素被添加，那么在行的各自范围内的一个列号在另一个之后被设置为添加元素的列号。当压缩被调用时，未使用的元素（由列号 @p 表示
   *
   * - ）通过复制后续行的列号来消除，每行内的列号（对角线元素可能除外）被排序，这样就可以通过二进制搜索来寻找一个元素是否存在并确定其位置。
   *
   */
  std::unique_ptr<size_type[]> colnums;

  /**
   * 存储是否为这个对象调用了compress()函数。
   *
   */
  bool compressed;

  // Make all sparse matrices friends of this class.
  template <typename number>
  friend class SparseMatrix;
  template <typename number>
  friend class SparseLUDecomposition;
  template <typename number>
  friend class SparseILU;
  template <typename number>
  friend class ChunkSparseMatrix;

  friend class ChunkSparsityPattern;
  friend class DynamicSparsityPattern;

  // Also give access to internal details to the iterator/accessor classes.
  friend class SparsityPatternIterators::Iterator;
  friend class SparsityPatternIterators::Accessor;
  friend class ChunkSparsityPatternIterators::Accessor;
};

/**
 * 这个类以<a
 * href="https://en.wikipedia.org/wiki/Sparse_matrix">compressed row storage
 * (CSR)</a>的格式来存储数据的稀疏性模式，并作为SparseMatrix类的基础。
 * SparsityPattern的元素，对应于SparseMatrix对象可以存储非零条目的地方，被逐行存储。在每一行中，元素通常是按照列索引递增的顺序从左到右存储的；这一规则的例外情况是，如果矩阵是方形的（n_rows()
 * ==
 * n_columns()），那么对角线条目会被存储为每一行的第一个元素，以使应用雅可比或SSOR预处理程序等操作更快。因此，如果你在迭代器的帮助下遍历SparsityPattern的一行元素（使用
 * SparsityPattern::begin 和 SparsityPattern::end)
 * ，你会发现只要矩阵是方形的，每一行的元素就不是按列索引排序的（第一项是对角线，其次是按列索引排序的其他条目）。
 * 虽然这个类构成了SparseMatrix对象的存储格式的基础，因此在建立线性系统中起着核心作用，但由于它的信息存储方式，很少被直接建立。相反，人们通常会先通过一个中间格式，例如参见
 * step-2  教程以及文档模块  @ref Sparsity  。
 * 你可以使用begin()、end()、begin(row)和end(row)对模式中的条目进行迭代。这些函数返回一个类型为
 * SparsityPatternIterators::Iterator.
 * 的迭代器，当取消引用一个迭代器 @p it, 时，你可以访问
 * SparsityPatternIterators::Accessor,
 * 中的成员函数，如<tt>it->column()</tt>和<tt>it->row()</tt>。
 *
 *
 */
class SparsityPattern : public SparsityPatternBase
{
public:
  /**
   * 声明容器尺寸的类型。
   *
   */
  using size_type = SparsityPatternBase::size_type;

  /**
   * 类型化一个迭代器类，允许在一个稀疏模式的所有非零元素上行走。
   *
   */
  using const_iterator = SparsityPatternBase::const_iterator;

  /**
   * Typedef一个迭代器类，允许在一个稀疏模式的所有非零元素上行走。
   * 由于该迭代器不允许修改稀疏模式，该类型与 @p
   * const_iterator. 的类型相同。
   *
   */
  using iterator = SparsityPatternBase::iterator;

  /**
   * 由于这个类只需要实现一个reinit()函数，我们需要将所有的基础reinit()函数带入范围，以便编译器能够找到它们。
   *
   */
  using SparsityPatternBase::reinit;

  /**
   * @name  构造和设置 构造器、析构器、初始化、复制和填充对象的函数。
   *
   */
  // @{
  /**
   * 初始化矩阵是空的，也就是没有分配内存。如果你想让这样的对象作为其他类中的成员变量，这是很有用的。你可以通过调用reinit()函数使该结构可用。
   *
   */
  SparsityPattern();

  /**
   * 复制构造函数。只有当要复制的矩阵结构为空时，才允许调用该构造函数。这样做是为了防止非自愿的复制对象的临时性，这可能会使用大量的计算时间。然而，如果想把SparsityPattern放在一个容器中，例如，写诸如<tt>v.push_back(SparsityPattern());</tt>这样的语句，<tt>v</tt>一个
   * std::vector 的SparsityPattern对象，则需要复制构造函数。
   * 通常，使用显式关键字来禁止不需要的暂存器就足够了，但这对
   * <tt>std::vector</tt>s.
   * 不起作用，因为无论如何复制这样的结构是没有用的，因为多个矩阵可以使用相同的稀疏度结构，所以只允许对空对象进行复制，如上所述。
   *
   */
  SparsityPattern(const SparsityPattern &);

  /**
   * 初始化一个大小为<tt>m x n</tt>的矩形图案。      @param[in]
   * m 行的数量。    @param[in]  n 列的数量。    @param[in]
   * max_per_row 每行的最大非零条目数。
   *
   */
  SparsityPattern(const size_type    m,
                  const size_type    n,
                  const unsigned int max_per_row);


  /**
   * 初始化一个大小为<tt>m x n</tt>的矩形图案。      @param[in]
   * m 行的数量。    @param[in]  n 列的数量。    @param[in]
   * row_lengths 每行的非零条目的可能数量。
   * 这个向量的每一行必须有一个条目。
   *
   */
  SparsityPattern(const size_type                  m,
                  const size_type                  n,
                  const std::vector<unsigned int> &row_lengths);

  /**
   * 初始化一个维度为<tt>m</tt>的二次元模式，每行最多有<tt>max_per_row</tt>非零条目。
   * 这个构造函数自动启用对角线元素的优化存储。为了避免这种情况，请使用分别取行号和列号的构造函数。
   *
   */
  SparsityPattern(const size_type m, const unsigned int max_per_row);

  /**
   * 初始化一个大小为<tt>m x m</tt>的二次方格。      @param[in]
   * m 行和列的数量。    @param[in]  row_lengths
   * 每行的最大非零条目数。
   * 这个向量的每一行必须有一个条目。
   *
   */
  SparsityPattern(const size_type                  m,
                  const std::vector<unsigned int> &row_lengths);

  /**
   * 制作一个带有额外对角线的副本。
   * 这样构造的对象是为了应用ILU(n)方法或其他不完全分解。
   * 因此，在原始入口结构之外，在主对角线的两边为<tt>extra_off_diagonals</tt>侧对角线提供空间。
   * <tt>max_per_row</tt>是该结构每行容纳非零元素的最大数量。假设这个数字足够大，以容纳<tt>original</tt>中的元素以及由这个构造函数创建的新的非对角线元素。你通常希望给出与你给<tt>original</tt>相同的数字，再加上对角线的数量乘以2。然而，如果你希望根据其他标准而不是在边对角线上为分解增加更多的非零条目，你可以给一个更大的值。
   * 这个函数要求<tt>original</tt>指的是一个二次方矩阵结构。
   * 它必须被压缩。这个函数完成后，矩阵结构不会被压缩。
   *
   */
  SparsityPattern(const SparsityPattern &original,
                  const unsigned int     max_per_row,
                  const size_type        extra_off_diagonals);

  /**
   * 解构器。
   *
   */
  ~SparsityPattern() override = default;

  /**
   * 复制操作符。对于这一点，与复制构造函数的情况相同：它被声明、定义并可以被调用，但后者只针对空对象。
   *
   */
  SparsityPattern &
  operator=(const SparsityPattern &);

  /**
   * 为一个大小为 @p m 乘以 @p n.
   * 的矩阵重新分配内存，每一行的条目数取自ArrayView  @p
   * row_lengths ，它必须给出每一行的这个数字 $i=0\ldots m-1$
   * 。
   *
   */
  virtual void
  reinit(const size_type                      m,
         const size_type                      n,
         const ArrayView<const unsigned int> &row_lengths) override;

  /**
   * 这个函数压缩了这个对象所代表的稀疏结构。
   * 它通过消除未使用的条目和对剩余的条目进行排序来实现，以便通过使用二进制搜索算法更快地访问。一个特殊的排序方案被用于二次矩阵的对角线条目，它总是每行的第一个条目。
   * 不再需要的内存会被释放。
   * SparseMatrix对象需要对其初始化的SparsityPattern对象进行压缩，以减少内存需求。
   *
   */
  void
  compress();


  /**
   * 如果你事先确切地知道将形成矩阵稀疏模式的条目，这个函数可以用来替代reinit()、对add()的后续调用和对close()的最后调用。
   * 前两个参数决定了矩阵的大小。对于最后两个，请注意稀疏矩阵可以由一连串的行来描述，每一个行都由一连串的列索引和值来表示。在这里，begin()和end()参数指定了进入一个容器的迭代器（正向迭代器类型），一个代表一行。因此，begin()和end()之间的距离应该等于n_rows()。这些迭代器可以是
   * <tt>std::vector</tt>,   <tt>std::list</tt>,
   * 指向C风格数组的迭代器，或者任何其他满足正向迭代器要求的迭代器。这些迭代器所指向的对象（即我们在对这些迭代器之一应用<tt>operator*</tt>或<tt>operator-></tt>后得到的东西）必须是一个容器本身，它提供了<tt>begin</tt>和<tt>end</tt>函数，指定了描述一行内容的迭代器的范围。解除这些内部迭代器必须产生一对作为列索引的无符号整数和一个任意类型的值（如果我们想用一个这样的对象来描述一个稀疏矩阵，就会使用这样的类型），或者只是一个无符号整数（如果我们只想描述一个稀疏的模式）。该函数能够通过一些模板魔法来确定我们在解读内部迭代器后得到的是无符号整数还是一对整数。
   * 虽然外迭代器的顺序表示矩阵的不同行，但表示列的内迭代器的顺序并不重要，因为无论如何它们在这个函数的内部是被排序的。
   * 由于这一切听起来非常复杂，请考虑下面的示例代码，它可能被用来填补一个稀疏模式。
   * @code
   * std::vector<std::vector<unsigned int> > column_indices (n_rows);
   * for (unsigned int row=0; row<n_rows; ++row)
   *       // generate necessary columns in this row
   * fill_row (column_indices[row]);
   *
   * sparsity.copy_from (n_rows, n_cols,
   *                   column_indices.begin(),
   *                   column_indices.end());
   * @endcode
   * 注意这个例子是有效的，因为被解读的迭代器产生的容器有<tt>begin</tt>和<tt>end</tt>函数（即
   * <tt>std::vector</tt>s),
   * ，被解读的内部迭代器产生无符号整数作为列索引。请注意，我们可以用
   * <tt>std::list</tt>, 替换两个 <tt>std::vector</tt>
   * 的出现，也可以用 <tt>std::set</tt> 替换内部的。
   * 另一个例子如下，我们初始化整个矩阵，而不仅仅是一个稀疏模式。
   * @code
   * std::vector<std::map<unsigned int,double> > entries (n_rows);
   * for (unsigned int row=0; row<n_rows; ++row)
   *       // generate necessary pairs of columns
   *       // and corresponding values in this row
   * fill_row (entries[row]);
   *
   * sparsity.copy_from (n_rows, n_cols,
   *                   column_indices.begin(),
   *                   column_indices.end());
   * matrix.reinit (sparsity);
   * matrix.copy_from (column_indices.begin(),
   *                 column_indices.end());
   * @endcode
   * 这个例子是可行的，因为解读内部类型的迭代器会产生一对无符号整数和一个值，我们把其中的第一个作为列索引。如前所述，外部的
   * <tt>std::vector</tt> 可以用 <tt>std::list</tt>, 代替，内部的
   * <tt>std::map<unsigned  int,double></tt>可以用
   * <tt>std::vector<std::pair<unsigned
   * int,double></tt>代替，或者用一个列表或一组这样的对，因为它们都返回指向这种对的迭代器。
   *
   */
  template <typename ForwardIterator>
  void
  copy_from(const size_type       n_rows,
            const size_type       n_cols,
            const ForwardIterator begin,
            const ForwardIterator end);

  /**
   * 从一个DynamicSparsityPattern中复制数据。这个对象以前的内容会丢失，之后的稀疏模式会处于压缩模式。
   *
   */
  void
  copy_from(const DynamicSparsityPattern &dsp);

  /**
   * 从一个SparsityPattern复制数据。这个对象以前的内容会丢失，而之后的稀疏模式处于压缩模式。
   *
   */
  void
  copy_from(const SparsityPattern &sp);

  /**
   * 取一个完整的矩阵并使用其非零条目为这个对象生成一个稀疏矩阵条目模式。
   * 这个对象以前的内容会丢失，之后的稀疏模式处于压缩模式。
   * 一旦你用这个函数建立了一个稀疏模式，你可能想给它附加一个SparseMatrix对象。然后可以使用以FullMatrix对象为参数的
   * SparseMatrix::copy_from()
   * 版本，将原`matrix`对象复制到这个SparseMatrix对象中。通过这个过程，你可以将一个FullMatrix转换为一个SparseMatrix。
   *
   */
  template <typename number>
  void
  copy_from(const FullMatrix<number> &matrix);

  /**
   * 在指定的矩阵行中添加几个非零条目。
   * 这个函数只能用于非压缩的稀疏模式。
   * 如果其中一些条目已经存在，则不会发生任何坏事。
   *
   */
  template <typename ForwardIterator>
  void
  add_entries(const size_type row,
              ForwardIterator begin,
              ForwardIterator end,
              const bool      indices_are_sorted = false);

  // @}


  /**
   * @name  查询信息
   *
   */
  // @{
  /**
   * 测试两个SparsityPatterns是否相等。
   *
   */
  bool
  operator==(const SparsityPattern &) const;

  /**
   * 返回这个对象是否只存储那些显式添加的条目，或者稀疏模式是否包含在构建它时通过其他方式（隐式）添加的元素。对于当前的类，当且仅当它是正方形时，结果是假的，因为此时它无条件地存储对角线条目，无论它们是否被显式添加。
   * 这个函数的主要作用是在几种稀疏模式可以作为模板参数传递的情况下描述当前类。
   *
   */
  bool
  stores_only_added_elements() const;

  /**
   * 确定此对象的内存消耗（以字节为单位）的估计值。参见MemoryConsumption。
   *
   */
  std::size_t
  memory_consumption() const;

  // @}
  /**
   * @name  访问条目
   *
   */
  // @{
  /**
   * 返回行号为<tt>i</tt>、列号为<tt>j</tt>的矩阵元素的索引。如果矩阵元素不是非零，则返回
   * SparsityPattern::invalid_entry.  这个函数通常由
   * SparseMatrix::operator()().
   * 它可能只被用于压缩稀疏模式，因为在这种情况下，搜索条目是否存在可以用二进制排序算法相当快地完成，因为列号被排序了。
   * 如果<tt>m</tt>是<tt>row</tt>中的条目数，那么如果稀疏模式被压缩，这个函数的复杂度是<i>log(m)</i>
   * 。
   * @note
   * 这个函数并不便宜，因为它必须在给定行<tt>i</tt>的所有元素中进行搜索，以寻找索引<tt>j</tt>是否存在。因此，在你想循环查看这个稀疏模式（或与之相关的稀疏矩阵）的所有非零元素或单一行的非零元素的情况下，它比必要的要昂贵。在这种情况下，使用遍历稀疏模式或稀疏矩阵的元素的迭代器会更有效。
   *
   */
  size_type
  operator()(const size_type i, const size_type j) const;

  // @}
  /**
   * @name  输入/输出
   *
   */
  // @{

  /**
   * 将此对象的数据全部写到文件中。这是以二进制模式进行的，所以输出的数据既不能被人类阅读，也不能（可能）被其他使用不同操作系统或数字格式的计算机阅读。
   * 这个函数的目的是，如果你的内存不足，想在不同的程序之间进行交流，或者允许对象在程序的不同运行中持续存在，你可以把矩阵和稀疏模式换掉。
   *
   */
  void
  block_write(std::ostream &out) const;

  /**
   * 从文件中读取先前由block_write()写入的数据。
   * 这是用上述函数的逆运算来完成的，所以它的速度相当快，因为除了前面的几个数字，比特流是不被解释的。
   * 在这个操作中，对象被调整了大小，所有以前的内容都被丢失。
   * 一个原始形式的错误检查被执行，它将识别最直白的尝试，将一些数据解释为按位数存储到文件的向量，但不会超过。
   *
   */
  void
  block_read(std::istream &in);

  /**
   * 为了序列化的目的，将此对象的数据写入一个流中
   *
   */
  template <class Archive>
  void
  save(Archive &ar, const unsigned int version) const;

  /**
   * 为了序列化的目的，从一个流中读取此对象的数据
   *
   */
  template <class Archive>
  void
  load(Archive &ar, const unsigned int version);

#ifdef DOXYGEN
  /**
   * 为了序列化的目的，从一个流中写和读这个对象的数据。
   *
   */
  template <class Archive>
  void
  serialize(Archive &archive, const unsigned int version);
#else
  // This macro defines the serialize() method that is compatible with
  // the templated save() and load() method that have been implemented.
  BOOST_SERIALIZATION_SPLIT_MEMBER()
#endif

  // @}

  /**
   * @addtogroup  Exceptions  @{
   *
   */
  /**
   * 异常情况
   *
   */
  DeclException2(ExcIteratorRange,
                 int,
                 int,
                 << "The iterators denote a range of " << arg1
                 << " elements, but the given number of rows was " << arg2);
  /**
   * 异常情况
   *
   */
  DeclException1(ExcInvalidNumberOfPartitions,
                 int,
                 << "The number of partitions you gave is " << arg1
                 << ", but must be greater than zero.");
  //@}
private:
  /**
   * 是否启用对角线的特殊处理？
   *
   */
  bool store_diagonal_first_in_row;

  // Make all sparse matrices friends of this class.
  template <typename number>
  friend class SparseMatrix;
  template <typename number>
  friend class SparseLUDecomposition;
  template <typename number>
  friend class SparseILU;
  template <typename number>
  friend class ChunkSparseMatrix;

  friend class ChunkSparsityPattern;
  friend class DynamicSparsityPattern;

  // Also give access to internal details to the iterator/accessor classes.
  friend class SparsityPatternIterators::Iterator;
  friend class SparsityPatternIterators::Accessor;
  friend class ChunkSparsityPatternIterators::Accessor;
};


 /*@}*/ 
 /*---------------------- Inline functions -----------------------------------*/ 

#ifndef DOXYGEN


namespace SparsityPatternIterators
{
  inline Accessor::Accessor(const SparsityPatternBase *sparsity_pattern,
                            const std::size_t          i)
    : container(sparsity_pattern)
    , linear_index(i)
  {}



  inline Accessor::Accessor(const SparsityPatternBase *sparsity_pattern)
    : container(sparsity_pattern)
    , linear_index(container->rowstart[container->rows])
  {}



  inline Accessor::Accessor()
    : container(nullptr)
    , linear_index(numbers::invalid_size_type)
  {}



  inline bool
  Accessor::is_valid_entry() const
  {
    Assert(container != nullptr, DummyAccessor());
    return (linear_index < container->rowstart[container->rows] &&
            container->colnums[linear_index] != SparsityPattern::invalid_entry);
  }



  inline size_type
  Accessor::row() const
  {
    Assert(is_valid_entry() == true, ExcInvalidIterator());

    const std::size_t *insert_point =
      std::upper_bound(container->rowstart.get(),
                       container->rowstart.get() + container->rows + 1,
                       linear_index);
    return insert_point - container->rowstart.get() - 1;
  }



  inline size_type
  Accessor::column() const
  {
    Assert(is_valid_entry() == true, ExcInvalidIterator());

    return (container->colnums[linear_index]);
  }



  inline size_type
  Accessor::index() const
  {
    Assert(is_valid_entry() == true, ExcInvalidIterator());

    return linear_index - container->rowstart[row()];
  }



  inline size_type
  Accessor::global_index() const
  {
    Assert(is_valid_entry() == true, ExcInvalidIterator());

    return linear_index;
  }



  inline bool
  Accessor::operator==(const Accessor &other) const
  {
    return (container == other.container && linear_index == other.linear_index);
  }



  inline bool
  Accessor::operator<(const Accessor &other) const
  {
    Assert(container != nullptr, DummyAccessor());
    Assert(other.container != nullptr, DummyAccessor());
    Assert(container == other.container, ExcInternalError());

    return linear_index < other.linear_index;
  }



  inline void
  Accessor::advance()
  {
    Assert(container != nullptr, DummyAccessor());
    Assert(linear_index < container->rowstart[container->rows],
           ExcIteratorPastEnd());
    ++linear_index;
  }


  inline Iterator::Iterator(const SparsityPatternBase *sp,
                            const std::size_t          linear_index)
    : LinearIndexIterator<Iterator, Accessor>(Accessor(sp, linear_index))
  {}


  inline Iterator::Iterator(const Accessor &accessor)
    : LinearIndexIterator<Iterator, Accessor>(accessor)
  {}


} // namespace SparsityPatternIterators



inline SparsityPatternBase::iterator
SparsityPatternBase::begin() const
{
  if (n_rows() > 0)
    return {this, rowstart[0]};
  else
    return end();
}



inline SparsityPatternBase::iterator
SparsityPatternBase::end() const
{
  if (n_rows() > 0)
    return {this, rowstart[rows]};
  else
    return {nullptr, 0};
}



inline SparsityPatternBase::iterator
SparsityPatternBase::begin(const size_type r) const
{
  AssertIndexRange(r, n_rows());

  return {this, rowstart[r]};
}



inline SparsityPatternBase::iterator
SparsityPatternBase::end(const size_type r) const
{
  AssertIndexRange(r, n_rows());

  return {this, rowstart[r + 1]};
}



inline SparsityPatternBase::size_type
SparsityPatternBase::n_rows() const
{
  return rows;
}



inline SparsityPatternBase::size_type
SparsityPatternBase::n_cols() const
{
  return cols;
}



inline bool
SparsityPatternBase::is_compressed() const
{
  return compressed;
}



inline bool
SparsityPattern::stores_only_added_elements() const
{
  return (store_diagonal_first_in_row == false);
}



inline unsigned int
SparsityPatternBase::row_length(const size_type row) const
{
  AssertIndexRange(row, rows);
  return rowstart[row + 1] - rowstart[row];
}



inline SparsityPattern::size_type
SparsityPatternBase::column_number(const size_type    row,
                                   const unsigned int index) const
{
  AssertIndexRange(row, rows);
  AssertIndexRange(index, row_length(row));

  return colnums[rowstart[row] + index];
}



inline std::size_t
SparsityPatternBase::n_nonzero_elements() const
{
  Assert(compressed, ExcNotCompressed());

  if ((rowstart != nullptr) && (colnums != nullptr))
    return rowstart[rows] - rowstart[0];
  else
    // the object is empty or has zero size
    return 0;
}



template <class Archive>
inline void
SparsityPatternBase::save(Archive &ar, const unsigned int) const
{
  // forward to serialization function in the base class.
  ar &boost::serialization::base_object<const Subscriptor>(*this);

  ar &max_dim &rows &cols &max_vec_len &max_row_length &compressed;

  ar &boost::serialization::make_array(rowstart.get(), max_dim + 1);
  ar &boost::serialization::make_array(colnums.get(), max_vec_len);
}



template <class Archive>
inline void
SparsityPatternBase::load(Archive &ar, const unsigned int)
{
  // forward to serialization function in the base class.
  ar &boost::serialization::base_object<Subscriptor>(*this);

  ar &max_dim &rows &cols &max_vec_len &max_row_length &compressed;

  rowstart = std::make_unique<std::size_t[]>(max_dim + 1);
  colnums  = std::make_unique<size_type[]>(max_vec_len);

  ar &boost::serialization::make_array(rowstart.get(), max_dim + 1);
  ar &boost::serialization::make_array(colnums.get(), max_vec_len);
}



template <class Archive>
inline void
SparsityPattern::save(Archive &ar, const unsigned int) const
{
  // forward to serialization function in the base class.
  ar &boost::serialization::base_object<const SparsityPatternBase>(*this);
  ar &store_diagonal_first_in_row;
}



template <class Archive>
inline void
SparsityPattern::load(Archive &ar, const unsigned int)
{
  // forward to serialization function in the base class.
  ar &boost::serialization::base_object<SparsityPatternBase>(*this);
  ar &store_diagonal_first_in_row;
}



inline bool
SparsityPatternBase::operator==(const SparsityPatternBase &sp2) const
{
  // it isn't quite necessary to compare *all* member variables. by only
  // comparing the essential ones, we can say that two sparsity patterns are
  // equal even if one is compressed and the other is not (in which case some
  // of the member variables are not yet set correctly)
  if (rows != sp2.rows || cols != sp2.cols || compressed != sp2.compressed)
    return false;

  for (size_type i = 0; i < rows + 1; ++i)
    if (rowstart[i] != sp2.rowstart[i])
      return false;

  for (size_type i = 0; i < rowstart[rows]; ++i)
    if (colnums[i] != sp2.colnums[i])
      return false;

  return true;
}



inline bool
SparsityPattern::operator==(const SparsityPattern &sp2) const
{
  return (static_cast<const SparsityPatternBase &>(*this) == sp2) &&
         (store_diagonal_first_in_row == sp2.store_diagonal_first_in_row);
}



namespace internal
{
  namespace SparsityPatternTools
  {
    /**
     * 声明容器大小的类型。
     *
     */
    using size_type = types::global_dof_index;

    inline size_type
    get_column_index_from_iterator(const size_type i)
    {
      return i;
    }



    template <typename value>
    inline size_type
    get_column_index_from_iterator(const std::pair<size_type, value> &i)
    {
      return i.first;
    }



    template <typename value>
    inline size_type
    get_column_index_from_iterator(const std::pair<const size_type, value> &i)
    {
      return i.first;
    }
  } // namespace SparsityPatternTools
} // namespace internal



template <typename ForwardIterator>
void
SparsityPattern::copy_from(const size_type       n_rows,
                           const size_type       n_cols,
                           const ForwardIterator begin,
                           const ForwardIterator end)
{
  Assert(static_cast<size_type>(std::distance(begin, end)) == n_rows,
         ExcIteratorRange(std::distance(begin, end), n_rows));

  // first determine row lengths for each row. if the matrix is quadratic,
  // then we might have to add an additional entry for the diagonal, if that
  // is not yet present. as we have to call compress anyway later on, don't
  // bother to check whether that diagonal entry is in a certain row or not
  const bool                is_square = (n_rows == n_cols);
  std::vector<unsigned int> row_lengths;
  row_lengths.reserve(n_rows);
  for (ForwardIterator i = begin; i != end; ++i)
    row_lengths.push_back(std::distance(i->begin(), i->end()) +
                          (is_square ? 1 : 0));
  reinit(n_rows, n_cols, row_lengths);

  // now enter all the elements into the matrix. note that if the matrix is
  // quadratic, then we already have the diagonal element preallocated
  //
  // for use in the inner loop, we define an alias to the type of the inner
  // iterators
  size_type row = 0;
  using inner_iterator =
    typename std::iterator_traits<ForwardIterator>::value_type::const_iterator;
  for (ForwardIterator i = begin; i != end; ++i, ++row)
    {
      size_type *          cols = &colnums[rowstart[row]] + (is_square ? 1 : 0);
      const inner_iterator end_of_row = i->end();
      for (inner_iterator j = i->begin(); j != end_of_row; ++j)
        {
          const size_type col =
            internal::SparsityPatternTools::get_column_index_from_iterator(*j);
          AssertIndexRange(col, n_cols);

          if ((col != row) || !is_square)
            *cols++ = col;
        }
    }

  // finally compress everything. this also sorts the entries within each row
  compress();
}


#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif


