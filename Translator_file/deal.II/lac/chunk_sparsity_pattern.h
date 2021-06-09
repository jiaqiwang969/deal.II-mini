//include/deal.II-translator/lac/chunk_sparsity_pattern_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2020 by the deal.II authors
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

#ifndef dealii_chunk_sparsity_pattern_h
#define dealii_chunk_sparsity_pattern_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/subscriptor.h>

#include <deal.II/lac/sparsity_pattern.h>

#include <iostream>
#include <vector>

DEAL_II_NAMESPACE_OPEN


// Forward declaration
#ifndef DOXYGEN
template <typename>
class ChunkSparseMatrix;
#endif

/*!   @addtogroup Sparsity  
     * @{  

* 
*
*/



/**
 * 稀疏性模式的迭代器
 *
 *
 */
namespace ChunkSparsityPatternIterators
{
  // forward declaration
  class Iterator;

  /**
   * 进入稀疏模式的迭代器的访问器类。这个类也是进入稀疏矩阵的常量和非常量访问器类的基类。
   * 请注意，这个类只允许对元素进行读取访问，提供它们的行和列号。它不允许修改稀疏模式本身。
   *
   */
  class Accessor
  {
  public:
    /**
     * 声明容器大小的类型。
     *
     */
    using size_type = types::global_dof_index;

    /**
     * 构造函数。
     *
     */
    Accessor(const ChunkSparsityPattern *matrix, const size_type row);

    /**
     * 构造器。为给定的稀疏模式构造结束访问器。
     *
     */
    Accessor(const ChunkSparsityPattern *matrix);

    /**
     * 这个对象所代表的元素的行号。这个函数只能为is_valid_entry()为真的条目调用。
     *
     */
    size_type
    row() const;

    /**
     * 返回全局索引，该索引来自减少的稀疏模式。
     *
     */
    std::size_t
    reduced_index() const;

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
    /**
     * 我们所操作的稀疏模式被访问。
     *
     */
    const ChunkSparsityPattern *sparsity_pattern;

    /**
     * (减少的)稀疏模式的访问器。
     *
     */
    SparsityPatternIterators::Accessor reduced_accessor;

    /**
     * 当前块的行数。
     *
     */
    size_type chunk_row;

    /**
     * 当前块的列数。
     *
     */
    size_type chunk_col;

    /**
     * 将访问器移到矩阵的下一个非零条目。
     *
     */
    void
    advance();

    // Grant access to iterator class.
    friend class Iterator;
  };



  /**
   * 遍历稀疏模式的元素的迭代器。
   *
   */
  class Iterator
  {
  public:
    /**
     * 声明容器大小的类型。
     *
     */
    using size_type = types::global_dof_index;

    /**
     * 构造函数。为给定的行和其中的索引创建一个进入稀疏模式
     * @p sp 的迭代器。
     *
     */
    Iterator(const ChunkSparsityPattern *sp, const size_type row);

    /**
     * 前缀增量。
     *
     */
    Iterator &
    operator++();

    /**
     * 后缀增量。
     *
     */
    Iterator
    operator++(int);

    /**
     * 撤消运算符。
     *
     */
    const Accessor &operator*() const;

    /**
     * 解除引用操作符。
     *
     */
    const Accessor *operator->() const;

    /**
     * 比较。真，如果两个迭代器都指向同一个矩阵位置。
     *
     */
    bool
    operator==(const Iterator &) const;

    /**
     * <tt>==</tt>的倒数。
     *
     */
    bool
    operator!=(const Iterator &) const;

    /**
     * 比较运算符。如果第一个行号较小，或者行号相等且第一个索引较小，则结果为真。
     * 这个函数只有在两个迭代器都指向同一个矩阵时才有效。
     *
     */
    bool
    operator<(const Iterator &) const;

  private:
    /**
     * 存储一个访问器类的对象。
     *
     */
    Accessor accessor;
  };
} // namespace ChunkSparsityPatternIterators



/**
 * 代表稀疏矩阵的稀疏模式的结构。这个类是  @ref Sparsity  的 "静态 "
 * 类型的一个例子。它使用压缩行存储（CSR）格式来存储数据。
 * 这个类的使用在  step-51  中得到了证明。
 *
 *
 */
class ChunkSparsityPattern : public Subscriptor
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
  using const_iterator = ChunkSparsityPatternIterators::Iterator;

  /**
   * Typedef一个迭代器类，允许在一个稀疏模式的所有非零元素上行走。
   * 由于该迭代器不允许修改稀疏模式，该类型与 @p
   * const_iterator. 的类型相同。
   *
   */
  using iterator = ChunkSparsityPatternIterators::Iterator;

  /**
   * 定义一个值，用来表示colnums数组中的某个值是未使用的，即不代表某个列号索引。
   * 具有这个无效值的索引被用来使用add()成员函数向稀疏模式插入新条目，并在调用compress()时被移除。
   * 你不应该认为这里声明的变量有一定的价值。这里给出的初始化只是为了让编译器进行一些优化，但变量的实际值可能会随着时间的推移而改变。
   *
   */
  static const size_type invalid_entry = SparsityPattern::invalid_entry;

  /**
   * 初始化矩阵为空，也就是没有分配内存。如果你想让这样的对象作为其他类中的成员变量，这很有用。你可以通过调用reinit()函数使该结构可用。
   *
   */
  ChunkSparsityPattern();

  /**
   * 复制构造函数。只有当要复制的矩阵结构为空时，才允许调用该构造函数。这样做是为了防止非自愿的复制对象的临时性，这可能会使用大量的计算时间。然而，如果想把ChunkSparsityPattern放在一个容器中，例如，写诸如<tt>v.push_back
   * (ChunkSparsityPattern());</tt>这样的语句，就需要复制构造函数，而<tt>v</tt>是ChunkSparsityPattern对象的向量。
   * 通常情况下，使用显式关键字来禁止不需要的暂存器就足够了，但这对
   * <tt>std::vector</tt>.
   * 不起作用，因为无论如何复制这样的结构是没有用的，因为多个矩阵可以使用同一个稀疏度结构，所以只允许对空对象进行复制，如上所述。
   *
   */
  ChunkSparsityPattern(const ChunkSparsityPattern &);

  /**
   * 初始化一个矩形的矩阵。      @arg  m行数  @arg  n列数  @arg
   * max_per_row 每行非零条目的最大数量
   *
   */
  ChunkSparsityPattern(const size_type m,
                       const size_type n,
                       const size_type max_chunks_per_row,
                       const size_type chunk_size);

  /**
   * 初始化一个矩形的矩阵。      @arg  m行数  @arg  n列数  @arg
   * row_lengths 每行可能的非零条目数。
   * 这个向量的每一行必须有一个条目。
   *
   */
  ChunkSparsityPattern(const size_type               m,
                       const size_type               n,
                       const std::vector<size_type> &row_lengths,
                       const size_type               chunk_size);

  /**
   * 初始化一个维度为<tt>n</tt>的二次矩阵，每行最多有<tt>max_per_row</tt>个非零条目。
   * 这个构造函数自动启用对角线元素的优化存储。为了避免这种情况，请使用分别取行号和列号的构造函数。
   *
   */
  ChunkSparsityPattern(const size_type n,
                       const size_type max_per_row,
                       const size_type chunk_size);

  /**
   * 初始化一个四维矩阵。      @arg  m 行数和列数  @arg
   * row_lengths 每一行的非零条目的可能数量。
   * 这个向量的每一行必须有一个条目。
   *
   */
  ChunkSparsityPattern(const size_type               m,
                       const std::vector<size_type> &row_lengths,
                       const size_type               chunk_size);

  /**
   * 解构器。
   *
   */
  ~ChunkSparsityPattern() override = default;

  /**
   * 复制操作符。对于这一点，与复制构造函数的情况相同：它被声明、定义并可以被调用，但后者只针对空对象。
   *
   */
  ChunkSparsityPattern &
  operator=(const ChunkSparsityPattern &);

  /**
   * 为一个新的矩阵重新分配内存并设置数据结构，该矩阵具有<tt>m
   * </tt>行和<tt>n</tt>列，每行最多具有<tt>max_per_row</tt>非零条目。
   * 这个函数只是将其操作映射到另一个<tt>reinit</tt>函数。
   *
   */
  void
  reinit(const size_type m,
         const size_type n,
         const size_type max_per_row,
         const size_type chunk_size);

  /**
   * 为一个大小为<tt>m x
   * n</tt>的矩阵重新分配内存。每一行的条目数取自数组<tt>row_lengths</tt>，它必须给出每一行的这个数字<tt>i=1...m</tt>。
   * 如果<tt>m*n==0</tt>所有的内存被释放，导致对象完全重新初始化。如果是非零，只有当新的大小扩展到旧的大小时才会分配新的内存。这样做是为了节省时间和避免堆的碎片化。
   * 如果行的数量等于列的数量，那么对角线元素将首先存储在每一行中，以便在SparseMatrix的放松方法中进行优化访问。
   *
   */
  void
  reinit(const size_type               m,
         const size_type               n,
         const std::vector<size_type> &row_lengths,
         const size_type               chunk_size);

  /**
   * 和上面一样，但用ArrayView参数代替。
   *
   */
  void
  reinit(const size_type                   m,
         const size_type                   n,
         const ArrayView<const size_type> &row_lengths,
         const size_type                   chunk_size);

  /**
   * 这个函数压缩这个对象所代表的稀疏结构。
   * 它通过消除未使用的条目并对剩余的条目进行排序，以便通过使用二进制搜索算法更快地访问。一个特殊的排序方案被用于二次矩阵的对角线条目，它总是每行的第一个条目。
   * 不再需要的内存会被释放。
   * SparseMatrix对象需要对其初始化的ChunkSparsityPattern对象进行压缩，以减少内存需求。
   *
   */
  void
  compress();

  /**
   * 如果你事先确切地知道将形成矩阵稀疏模式的条目，这个函数可以用来替代reinit()、对add()的后续调用和对close()的最后调用。
   * 前两个参数决定了矩阵的大小。对于最后两个，请注意稀疏矩阵可以由一连串的行来描述，每一个行都由一连串的列索引和值来表示。在这里，begin()和end()参数指定了进入一个容器的迭代器（正向迭代器类型），一个代表一行。因此，begin()和end()之间的距离应该等于n_rows()。这些迭代器可以是
   * <tt>std::vector</tt>,   <tt>std::list</tt>,
   * 指向C风格数组的迭代器，或者任何其他满足正向迭代器要求的迭代器。这些迭代器所指向的对象（即我们在对这些迭代器之一应用<tt>operator*</tt>或<tt>operator-></tt>之后得到的东西）必须是一个容器本身，它提供了<tt>begin</tt>和<tt>end</tt>函数，指定了描述一行内容的迭代器的范围。解除这些内部迭代器必须产生一对作为列索引的整数和一个任意类型的值（如果我们想用一个这样的对象来描述一个稀疏矩阵，就会使用这样的类型），或者只是一个整数（如果我们只想描述一个稀疏模式）。该函数能够通过一些模板魔法来确定我们在解引用内部迭代器后得到的是一个整数还是一对。
   * 虽然外迭代器的顺序表示矩阵的不同行，但表示列的内迭代器的顺序并不重要，因为无论如何它们在这个函数中都是被排序的。
   * 由于这一切听起来非常复杂，请考虑下面的示例代码，它可能被用来填补一个稀疏模式。
   * @code
   * std::vector<std::vector<size_type> > column_indices (n_rows);
   * for (size_type row=0; row<n_rows; ++row)
   *       // generate necessary columns in this row
   * fill_row (column_indices[row]);
   *
   * sparsity.copy_from (n_rows, n_cols,
   *                   column_indices.begin(),
   *                   column_indices.end());
   * @endcode
   * 请注意，这个例子是有效的，因为被解读的迭代器产生了具有<tt>begin</tt>和<tt>end</tt>函数的容器（即
   * <tt>std::vector</tt>s),
   * ，被解读的内部迭代器产生整数作为列索引。请注意，我们可以用
   * <tt>std::list</tt>, 替换两个 <tt>std::vector</tt>
   * 的出现，也可以用 <tt>std::set</tt> 替换内部的。
   * 另一个例子如下，我们初始化整个矩阵，而不仅仅是一个稀疏模式。
   * @code
   * std::vector<std::map<size_type,double> > entries (n_rows);
   * for (size_type row=0; row<n_rows; ++row)
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
   * 这个例子是可行的，因为解读内部类型的迭代器会产生一对整数和一个值，我们把其中的第一个值作为列索引。如前所述，外部的
   * <tt>std::vector</tt> 可以用 <tt>std::list</tt>, 代替，内部的
   * <tt>std::map<size_type,double></tt> 可以用
   * <tt>std::vector<std::pair<size_type,double>
   * 代替></tt>，或者一个列表或一组这样的对，因为它们都返回指向这样的对的迭代器。
   *
   */
  template <typename ForwardIterator>
  void
  copy_from(const size_type       n_rows,
            const size_type       n_cols,
            const ForwardIterator begin,
            const ForwardIterator end,
            const size_type       chunk_size);

  /**
   * 从一个DynamicSparsityPattern类型的对象中复制数据。这个对象以前的内容会丢失，之后的稀疏模式会处于压缩模式。
   *
   */
  template <typename SparsityPatternType>
  void
  copy_from(const SparsityPatternType &dsp, const size_type chunk_size);

  /**
   * 取一个完整的矩阵并使用其非零条目为这个对象生成一个稀疏矩阵条目模式。
   * 这个对象以前的内容会丢失，之后的稀疏模式处于压缩模式。
   *
   */
  template <typename number>
  void
  copy_from(const FullMatrix<number> &matrix, const size_type chunk_size);

  /**
   * 设置大块稀疏模式的稀疏模式，由<tt>chunk_size*chunksize</tt>块指定的大块稀疏模式给出。注意，疏散模式的最终行数<tt>m</tt>大约是<tt>sparsity_pattern_for_chunks.n_rows()
   * chunk_size</tt>(modulo padding elements in the last
   * chunk)，同样，列数<tt>n</tt>。
   * 这是一个特殊的初始化选项，以防你可以从一开始就知道大块的位置，而不需要使用<tt>make_sparsity_pattern</tt>调用来生成稀疏模式。这绕过了对块的搜索，但当然需要小心处理，以便给出一个正确的稀疏模式。
   * 这个对象以前的内容会丢失，而且之后的稀疏度模式处于压缩模式。
   *
   */
  template <typename Sparsity>
  void
  create_from(const size_type m,
              const size_type n,
              const Sparsity &sparsity_pattern_for_chunks,
              const size_type chunk_size,
              const bool      optimize_diagonal = true);

  /**
   * 返回该对象是否为空。如果没有分配内存，它就是空的，这与两个维度都是零是一样的。
   *
   */
  bool
  empty() const;

  /**
   * 返回构建此对象时作为参数给出的分块大小。
   *
   */
  size_type
  get_chunk_size() const;

  /**
   * 返回每行的最大条目数。在压缩之前，这等于给构造函数的数字，而在压缩之后，它等于用户实际分配的最大条目数。
   *
   */
  size_type
  max_entries_per_row() const;

  /**
   * 向矩阵添加一个非零条目。这个函数只能用于非压缩的稀疏模式。
   * 如果该条目已经存在，则不会发生任何坏事。
   *
   */
  void
  add(const size_type i, const size_type j);

  /**
   * 通过添加转置对象的稀疏度模式使稀疏度模式对称化。
   * 如果稀疏度模式不代表二次矩阵，这个函数会抛出一个异常。
   *
   */
  void
  symmetrize();

  /**
   * 返回该矩阵的行数，等于图像空间的维度。
   *
   */
  inline size_type
  n_rows() const;

  /**
   * 返回该矩阵的列数，相当于范围空间的维度。
   *
   */
  inline size_type
  n_cols() const;

  /**
   * 检查某个位置的数值是否可能为非零。
   *
   */
  bool
  exists(const size_type i, const size_type j) const;

  /**
   * 特定行中的条目数。
   *
   */
  size_type
  row_length(const size_type row) const;

  /**
   * 计算这个结构所代表的矩阵的带宽。带宽是 $|i-j|$
   * 的最大值，其中索引对 $(i,j)$
   * 代表矩阵的一个非零条目。因此， $n\times m$
   * 矩阵的最大带宽是 $\max\{n-1,m-1\}$  。
   *
   */
  size_type
  bandwidth() const;

  /**
   * 返回该矩阵的非零元素的数量。实际上，它返回的是稀疏模式中的条目数；如果任何一个条目碰巧是零，无论如何都会被计算在内。
   * 这个函数只有在矩阵结构被压缩的情况下才能被调用。否则就没有太大意义了。
   *
   */
  size_type
  n_nonzero_elements() const;

  /**
   * 返回该结构是否被压缩。
   *
   */
  bool
  is_compressed() const;

  /**
   * 返回这个对象是否只存储那些显式添加的条目，或者稀疏模式是否包含在构建时通过其他方式（隐式）添加的元素。对于当前的类，当且仅当它是正方形时，结果为真，因为此时它无条件地存储对角线条目，无论它们是否被显式添加。
   * 这个函数的主要作用是在几种稀疏模式可以作为模板参数传递的情况下描述当前类。
   *
   */
  bool
  stores_only_added_elements() const;

  /**
   * 迭代器从矩阵的第一个条目开始。由此产生的迭代器可以用来走过稀疏模式的所有非零条目。
   *
   */
  iterator
  begin() const;

  /**
   * 最后的迭代器。
   *
   */
  iterator
  end() const;

  /**
   * 迭代器从<tt>r</tt>行的第一个条目开始。
   * 注意，如果给定的行是空的，即不包含任何非零条目，那么这个函数返回的迭代器就等于<tt>end(r)</tt>。还要注意的是，在这种情况下，迭代器可能不能被解除引用。
   *
   */
  iterator
  begin(const size_type r) const;

  /**
   * 行<tt>r</tt>的最终迭代器。它指向过了 @p r,
   * 行末尾的第一个元素，或者过了整个稀疏模式的末尾。
   * 请注意，结束迭代器不一定是可被解除引用的。特别是如果它是一个矩阵的最后一行的结束迭代器，情况更是如此。
   *
   */
  iterator
  end(const size_type r) const;

  /**
   * 将此对象的数据全部写到文件中。这是以二进制模式进行的，所以输出的数据既不能被人类阅读，也不能（可能）被其他使用不同操作系统的数字格式的计算机阅读。
   * 这个函数的目的是，如果你的内存不足，想在不同的程序之间进行交流，或者允许对象在程序的不同运行中持续存在，你可以把矩阵和稀疏模式换出来。
   *
   */
  void
  block_write(std::ostream &out) const;

  /**
   * 从文件中读取先前由block_write()写入的数据。
   * 这是用上述函数的逆运算来完成的，所以它的速度相当快，因为除了前面的几个数字，比特流是不被解释的。
   * 在这个操作中，对象被调整了大小，所有以前的内容都被丢失。
   * 一个原始形式的错误检查被执行，它将识别最直白的尝试，将一些数据解释为一个向量，以位方式存储到文件中，但不会更多。
   *
   */
  void
  block_read(std::istream &in);

  /**
   * 打印矩阵的稀疏度。输出包括每行一行，格式为<tt>[i,j1,j2,j3,...]</tt>。<i>i</i>是行号，<i>jn</i>是该行中分配的列。
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
   * 确定此对象的内存消耗（以字节为单位）的估计值。见MemoryConsumption。
   *
   */
  std::size_t
  memory_consumption() const;

  /**
   * @addtogroup  Exceptions  @{
   *
   */
  /**
   * 异常情况
   *
   */
  DeclException1(ExcInvalidNumber,
                 size_type,
                 << "The provided number is invalid here: " << arg1);
  /**
   * 异常情况
   *
   */
  DeclException2(ExcInvalidIndex,
                 size_type,
                 size_type,
                 << "The given index " << arg1 << " should be less than "
                 << arg2 << ".");
  /**
   * 异常情况
   *
   */
  DeclException2(ExcNotEnoughSpace,
                 size_type,
                 size_type,
                 << "Upon entering a new entry to row " << arg1
                 << ": there was no free entry any more. " << std::endl
                 << "(Maximum number of entries for this row: " << arg2
                 << "; maybe the matrix is already compressed?)");
  /**
   * 只有在设置了SparsityPattern并且调用了compress()之后才允许进行该操作。
   *
   */
  DeclExceptionMsg(
    ExcNotCompressed,
    "The operation you attempted is only allowed after the SparsityPattern "
    "has been set up and compress() was called.");
  /**
   * 该操作改变了SparsityPattern的结构，在调用了compress()后不可能进行。
   *
   */
  DeclExceptionMsg(
    ExcMatrixIsCompressed,
    "The operation you attempted changes the structure of the SparsityPattern "
    "and is not possible after compress() has been called.");
  /**
   * 异常情况
   *
   */
  DeclException0(ExcEmptyObject);
  /**
   * 异常情况
   *
   */
  DeclException2(ExcIteratorRange,
                 size_type,
                 size_type,
                 << "The iterators denote a range of " << arg1
                 << " elements, but the given number of rows was " << arg2);
  /**
   * 异常情况
   *
   */
  DeclException0(ExcMETISNotInstalled);
  /**
   * 异常情况
   *
   */
  DeclException1(ExcInvalidNumberOfPartitions,
                 size_type,
                 << "The number of partitions you gave is " << arg1
                 << ", but must be greater than zero.");
  /**
   * 异常情况
   *
   */
  DeclException2(ExcInvalidArraySize,
                 size_type,
                 size_type,
                 << "The array has size " << arg1 << " but should have size "
                 << arg2);
  //@}
private:
  /**
   * 该稀疏结构应代表的行数。
   *
   */
  size_type rows;

  /**
   * 该稀疏结构应代表的列数。
   *
   */
  size_type cols;

  /**
   * 块的大小。
   *
   */
  size_type chunk_size;

  /**
   * 减少的稀疏性模式。我们只存储存在的块，每个块是大小为chunk_size
   * by chunk_size的矩阵中的块。
   *
   */
  SparsityPattern sparsity_pattern;

  // Make all the chunk sparse matrix kinds friends.
  template <typename>
  friend class ChunkSparseMatrix;

  // Make the accessor class a friend.
  friend class ChunkSparsityPatternIterators::Accessor;
};


 /*@}*/ 
 /*---------------------- Inline functions -----------------------------------*/ 

#ifndef DOXYGEN

namespace ChunkSparsityPatternIterators
{
  inline Accessor::Accessor(const ChunkSparsityPattern *sparsity_pattern,
                            const size_type             row)
    : sparsity_pattern(sparsity_pattern)
    , reduced_accessor(row == sparsity_pattern->n_rows() ?
                         *sparsity_pattern->sparsity_pattern.end() :
                         *sparsity_pattern->sparsity_pattern.begin(
                           row / sparsity_pattern->get_chunk_size()))
    , chunk_row(row == sparsity_pattern->n_rows() ?
                  0 :
                  row % sparsity_pattern->get_chunk_size())
    , chunk_col(0)
  {}



  inline Accessor::Accessor(const ChunkSparsityPattern *sparsity_pattern)
    : sparsity_pattern(sparsity_pattern)
    , reduced_accessor(*sparsity_pattern->sparsity_pattern.end())
    , chunk_row(0)
    , chunk_col(0)
  {}



  inline bool
  Accessor::is_valid_entry() const
  {
    return reduced_accessor.is_valid_entry() &&
           sparsity_pattern->get_chunk_size() * reduced_accessor.row() +
               chunk_row <
             sparsity_pattern->n_rows() &&
           sparsity_pattern->get_chunk_size() * reduced_accessor.column() +
               chunk_col <
             sparsity_pattern->n_cols();
  }



  inline Accessor::size_type
  Accessor::row() const
  {
    Assert(is_valid_entry() == true, ExcInvalidIterator());

    return sparsity_pattern->get_chunk_size() * reduced_accessor.row() +
           chunk_row;
  }



  inline Accessor::size_type
  Accessor::column() const
  {
    Assert(is_valid_entry() == true, ExcInvalidIterator());

    return sparsity_pattern->get_chunk_size() * reduced_accessor.column() +
           chunk_col;
  }



  inline std::size_t
  Accessor::reduced_index() const
  {
    Assert(is_valid_entry() == true, ExcInvalidIterator());

    return reduced_accessor.linear_index;
  }



  inline bool
  Accessor::operator==(const Accessor &other) const
  {
    // no need to check for equality of sparsity patterns as this is done in
    // the reduced case already and every ChunkSparsityPattern has its own
    // reduced sparsity pattern
    return (reduced_accessor == other.reduced_accessor &&
            chunk_row == other.chunk_row && chunk_col == other.chunk_col);
  }



  inline bool
  Accessor::operator<(const Accessor &other) const
  {
    Assert(sparsity_pattern == other.sparsity_pattern, ExcInternalError());

    if (chunk_row != other.chunk_row)
      {
        if (reduced_accessor.linear_index ==
            reduced_accessor.container->n_nonzero_elements())
          return false;
        if (other.reduced_accessor.linear_index ==
            reduced_accessor.container->n_nonzero_elements())
          return true;

        const auto global_row = sparsity_pattern->get_chunk_size() *
                                  reduced_accessor.row() +
                                chunk_row,
                   other_global_row = sparsity_pattern->get_chunk_size() *
                                        other.reduced_accessor.row() +
                                      other.chunk_row;
        if (global_row < other_global_row)
          return true;
        else if (global_row > other_global_row)
          return false;
      }

    return (
      reduced_accessor.linear_index < other.reduced_accessor.linear_index ||
      (reduced_accessor.linear_index == other.reduced_accessor.linear_index &&
       chunk_col < other.chunk_col));
  }


  inline void
  Accessor::advance()
  {
    const auto chunk_size = sparsity_pattern->get_chunk_size();
    Assert(chunk_row < chunk_size && chunk_col < chunk_size,
           ExcIteratorPastEnd());
    Assert(reduced_accessor.row() * chunk_size + chunk_row <
               sparsity_pattern->n_rows() &&
             reduced_accessor.column() * chunk_size + chunk_col <
               sparsity_pattern->n_cols(),
           ExcIteratorPastEnd());
    if (chunk_size == 1)
      {
        reduced_accessor.advance();
        return;
      }

    ++chunk_col;

    // end of chunk
    if (chunk_col == chunk_size ||
        reduced_accessor.column() * chunk_size + chunk_col ==
          sparsity_pattern->n_cols())
      {
        const auto reduced_row = reduced_accessor.row();
        // end of row
        if (reduced_accessor.linear_index + 1 ==
            reduced_accessor.container->rowstart[reduced_row + 1])
          {
            ++chunk_row;

            chunk_col = 0;

            // end of chunk rows or end of matrix
            if (chunk_row == chunk_size ||
                (reduced_row * chunk_size + chunk_row ==
                 sparsity_pattern->n_rows()))
              {
                chunk_row = 0;
                reduced_accessor.advance();
              }
            // go back to the beginning of the same reduced row but with
            // chunk_row increased by one
            else
              reduced_accessor.linear_index =
                reduced_accessor.container->rowstart[reduced_row];
          }
        // advance within chunk
        else
          {
            reduced_accessor.advance();
            chunk_col = 0;
          }
      }
  }



  inline Iterator::Iterator(const ChunkSparsityPattern *sparsity_pattern,
                            const size_type             row)
    : accessor(sparsity_pattern, row)
  {}



  inline Iterator &
  Iterator::operator++()
  {
    accessor.advance();
    return *this;
  }



  inline Iterator
  Iterator::operator++(int)
  {
    const Iterator iter = *this;
    accessor.advance();
    return iter;
  }



  inline const Accessor &Iterator::operator*() const
  {
    return accessor;
  }



  inline const Accessor *Iterator::operator->() const
  {
    return &accessor;
  }


  inline bool
  Iterator::operator==(const Iterator &other) const
  {
    return (accessor == other.accessor);
  }



  inline bool
  Iterator::operator!=(const Iterator &other) const
  {
    return !(accessor == other.accessor);
  }


  inline bool
  Iterator::operator<(const Iterator &other) const
  {
    return accessor < other.accessor;
  }

} // namespace ChunkSparsityPatternIterators



inline ChunkSparsityPattern::iterator
ChunkSparsityPattern::begin() const
{
  return {this, 0};
}


inline ChunkSparsityPattern::iterator
ChunkSparsityPattern::end() const
{
  return {this, n_rows()};
}



inline ChunkSparsityPattern::iterator
ChunkSparsityPattern::begin(const size_type r) const
{
  AssertIndexRange(r, n_rows());
  return {this, r};
}



inline ChunkSparsityPattern::iterator
ChunkSparsityPattern::end(const size_type r) const
{
  AssertIndexRange(r, n_rows());
  return {this, r + 1};
}



inline ChunkSparsityPattern::size_type
ChunkSparsityPattern::n_rows() const
{
  return rows;
}


inline ChunkSparsityPattern::size_type
ChunkSparsityPattern::n_cols() const
{
  return cols;
}



inline ChunkSparsityPattern::size_type
ChunkSparsityPattern::get_chunk_size() const
{
  return chunk_size;
}



inline bool
ChunkSparsityPattern::is_compressed() const
{
  return sparsity_pattern.compressed;
}



template <typename ForwardIterator>
void
ChunkSparsityPattern::copy_from(const size_type       n_rows,
                                const size_type       n_cols,
                                const ForwardIterator begin,
                                const ForwardIterator end,
                                const size_type       chunk_size)
{
  Assert(static_cast<size_type>(std::distance(begin, end)) == n_rows,
         ExcIteratorRange(std::distance(begin, end), n_rows));

  // first determine row lengths for each row. if the matrix is quadratic,
  // then we might have to add an additional entry for the diagonal, if that
  // is not yet present. as we have to call compress anyway later on, don't
  // bother to check whether that diagonal entry is in a certain row or not
  const bool             is_square = (n_rows == n_cols);
  std::vector<size_type> row_lengths;
  row_lengths.reserve(n_rows);
  for (ForwardIterator i = begin; i != end; ++i)
    row_lengths.push_back(std::distance(i->begin(), i->end()) +
                          (is_square ? 1 : 0));
  reinit(n_rows, n_cols, row_lengths, chunk_size);

  // now enter all the elements into the matrix
  size_type row = 0;
  using inner_iterator =
    typename std::iterator_traits<ForwardIterator>::value_type::const_iterator;
  for (ForwardIterator i = begin; i != end; ++i, ++row)
    {
      const inner_iterator end_of_row = i->end();
      for (inner_iterator j = i->begin(); j != end_of_row; ++j)
        {
          const size_type col =
            internal::SparsityPatternTools::get_column_index_from_iterator(*j);
          Assert(col < n_cols, ExcInvalidIndex(col, n_cols));

          add(row, col);
        }
    }

  // finally compress everything. this also sorts the entries within each row
  compress();
}


#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif


