//include/deal.II-translator/lac/dynamic_sparsity_pattern_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2021 by the deal.II authors
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

#ifndef dealii_dynamic_sparsity_pattern_h
#define dealii_dynamic_sparsity_pattern_h


#include <deal.II/base/config.h>

#include <deal.II/base/index_set.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/exceptions.h>

#include <algorithm>
#include <iostream>
#include <vector>

DEAL_II_NAMESPACE_OPEN

// Forward declaration
#ifndef DOXYGEN
class DynamicSparsityPattern;
#endif

/*!   @addtogroup Sparsity  
     * @{  

* 
*
*/


/**
 * DynamicSparsityPattern类型对象的迭代器。
 *
 *
 */
namespace DynamicSparsityPatternIterators
{
  // forward declaration
  class Iterator;

  /**
   * 声明容器大小的类型。
   *
   */
  using size_type = types::global_dof_index;

  /**
   * 迭代器进入DynamicSparsityPattern类型的对象的访问器类。
   * 请注意，这个类只允许对元素进行读取访问，提供它们的行和列号（或者是完整的稀疏模式中的索引）。它不允许修改稀疏度模式本身。
   *
   */
  class Accessor
  {
  public:
    /**
     * 构造函数。
     *
     */
    Accessor(const DynamicSparsityPattern *sparsity_pattern,
             const size_type               row,
             const unsigned int            index_within_row);

    /**
     * 构造函数。为给定的稀疏模式构造结束访问器。
     *
     */
    Accessor(const DynamicSparsityPattern *sparsity_pattern);

    /**
     * 默认构造器创建一个假访问器。这个构造函数在这里只是为了能够在STL容器中存储访问器，例如
     * `std::vector`.  。
     *
     */
    Accessor();

    /**
     * 这个对象所代表的元素的行数。
     *
     */
    size_type
    row() const;

    /**
     * 这个对象所代表的元素在当前行中的索引。
     *
     */
    size_type
    index() const;

    /**
     * 这个对象所代表的元素的列号。
     *
     */
    size_type
    column() const;

    /**
     * 比较。如果两个迭代器都指向同一个矩阵位置，则为真。
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
                     " without DynamicSparsityPattern object, which"
                     " means that it is a dummy accessor that can"
                     " not do any operations.");

    /**
     * 我们所操作的稀疏模式被访问。
     *
     */
    const DynamicSparsityPattern *sparsity_pattern;

    /**
     * 我们当前指向的行。
     *
     */
    size_type current_row;

    /**
     * 一个指向当前行内的元素的指针，我们当前指向的元素。
     *
     */
    std::vector<size_type>::const_iterator current_entry;

    /**
     * 一个指向当前行结束的指针。我们存储这个指针是为了与行末迭代器进行比较，这样做比较便宜，否则它需要进行IndexSet转换，从行的索引到DynamicSparsityPattern的'lines'数组中的索引。
     *
     */
    std::vector<size_type>::const_iterator end_of_row;

    /**
     * 将访问器移动到矩阵中的下一个非零条目。
     *
     */
    void
    advance();

    // Grant access to iterator class.
    friend class Iterator;
  };



  /**
   * 一个迭代器类，用于在稀疏性模式的元素上行走。
   * 这些迭代器的典型用途是迭代稀疏模式的元素（或者，因为它们也是迭代相关矩阵的元素的基础，所以是迭代稀疏矩阵的元素），或者迭代单个行的元素。不能保证行的元素实际上是按照列数单调增加的顺序来遍历的。更多信息请参见SparsityPattern类的文档。
   * @note
   * 该类直接对DynamicSparsityPattern类的内部数据结构进行操作。因此，有些操作很便宜，有些则不然。特别是，访问指向的稀疏模式条目的列索引很便宜。另一方面，计算两个迭代器之间的距离是很昂贵的。因此，当你设计使用这些迭代器的算法时，通常的做法是不一次性循环疏散模式的<i>all</i>元素，而是在所有行上有一个外循环，并在这个循环中迭代这个行的元素。这样，你只需要解除对迭代器的引用以获得列索引，而通过使用循环索引可以避免对行索引的（昂贵）查找。
   *
   */
  class Iterator
  {
  public:
    /**
     * 构造函数。为给定的全局索引（即从第2行开始计算的给定元素的索引）创建一个进入稀疏模式
     * @p sp 的迭代器。
     *
     */
    Iterator(const DynamicSparsityPattern *sp,
             const size_type               row,
             const unsigned int            index_within_row);

    /**
     * 构造函数。创建一个无效的（结束）迭代器进入稀疏模式
     * @p sp.  。
     *
     */
    Iterator(const DynamicSparsityPattern *sp);

    /**
     * 默认的构造函数，创建一个无效的迭代器。这个构造函数在这里只是为了能够在STL容器中存储迭代器，如
     * `std::vector`. 。
     *
     */
    Iterator() = default;

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

    /**
     * 返回当前迭代器与参数之间的距离。这个距离是通过对当前迭代器应用operator++多少次才能得到参数（对于正的返回值），或operator--（对于负的返回值）。
     *
     */
    int
    operator-(const Iterator &p) const;

  private:
    /**
     * 存储一个访问器类的对象。
     *
     */
    Accessor accessor;
  };
} // namespace DynamicSparsityPatternIterators


/**
 * 这个类作为SparsityPattern类的一个中间形式。从界面上看，它主要代表一个SparsityPattern对象，该对象在任何时候都保持压缩。然而，由于最终的稀疏度模式在构建时并不为人所知，因此在任何时候都保持压缩的模式只能以增加内存或使用时的运行时间消耗为代价。这个类的主要目的是为了避免一些内存瓶颈，所以我们选择实现它的内存保守性。不过，所选择的数据格式太不适合用于实际的矩阵。因此，在实际的矩阵中使用它之前，有必要先把这个对象的数据复制到SparsityPattern类型的对象上。
 * 另一个观点是，这个类不需要预先分配一定的内存，而是根据需要增长。 在
 * @ref Sparsity
 * 模块的文档中可以找到关于稀疏性模式的广泛描述。
 * 这个类是  @ref Sparsity  的 "动态 "
 * 类型的一个例子。它以这样或那样的方式被用于大多数教程程序中。
 * <h3>Interface</h3>
 * 因为这个类是作为SparsityPattern类的中间替代物，所以它的接口大部分是相同的，在必要的地方有小的改动。特别是add()函数，以及查询稀疏模式属性的函数都是一样的。
 *
 *  <h3>Usage</h3> 这个类的用法在 step-2 （无约束）和 step-6
 * （有AffineConstraints）中解释，通常看起来如下。
 *
 * @code
 * DynamicSparsityPattern dynamic_pattern (dof_handler.n_dofs());
 * DoFTools::make_sparsity_pattern (dof_handler,
 *                                dynamic_pattern,
 *                                constraints);
 * SparsityPattern sp;
 * sp.copy_from (dynamic_pattern);
 * @endcode
 *
 *
 *
 */
class DynamicSparsityPattern : public Subscriptor
{
public:
  /**
   * 声明容器尺寸的类型。
   *
   */
  using size_type = types::global_dof_index;

  /**
   * Typedef一个迭代器类，允许在一个稀疏模式的所有非零元素上行走。
   * 由于迭代器不允许修改稀疏模式，这个类型与 @p
   * const_iterator. 的类型相同。
   *
   */
  using iterator = DynamicSparsityPatternIterators::Iterator;

  /**
   * 迭代器类的类型定义，允许在稀疏模式的所有非零元素上行走。
   *
   */
  using const_iterator = DynamicSparsityPatternIterators::Iterator;

  /**
   * 初始化为一个空对象。如果你想让这样的对象作为其他类的成员变量，这很有用。你可以通过调用reinit()函数使该结构可用。
   *
   */
  DynamicSparsityPattern();

  /**
   * 复制构造函数。这个构造函数只允许在要复制的稀疏结构为空时被调用。这样做是为了防止临时性的对象的非自愿复制，这可能会使用大量的计算时间。
   * 然而，如果你想把一个DynamicSparsityPattern放在一个容器中，例如写这样的语句，如<tt>v.push_back(DynamicSparsityPattern());</tt>，
   * @p v是一个 @p DynamicSparsityPattern
   * 对象的向量，就需要拷贝构造函数。
   *
   */
  DynamicSparsityPattern(const DynamicSparsityPattern &);

  /**
   * 初始化一个具有 @p m 行和 @p n 列的矩形稀疏模式。 @p
   * rowset 限制了对这个集合的行的元素的存储。
   * 添加这个集合以外的元素没有影响。默认参数保留所有条目。
   *
   */
  DynamicSparsityPattern(const size_type m,
                         const size_type n,
                         const IndexSet &rowset = IndexSet());

  /**
   * 使用给定的索引集创建一个方形SparsityPattern。总大小由
   * @p indexset 的大小给出，并且只有与 @p indexset
   * 中的索引相对应的行被存储在当前处理器上。
   *
   */
  DynamicSparsityPattern(const IndexSet &indexset);

  /**
   * 初始化一个尺寸为 @p n. 的方形图案。
   *
   */
  DynamicSparsityPattern(const size_type n);

  /**
   * 复制操作。对于这一点，与复制构造函数的情况相同：它被声明、定义并可以被调用，但后者只针对空对象。
   *
   */
  DynamicSparsityPattern &
  operator=(const DynamicSparsityPattern &);

  /**
   * 为具有 @p m 行和 @p n
   * 列的新稀疏模式重新分配内存并设置数据结构。 @p rowset
   * 限制了对这个集合的行的元素的存储。
   * 添加这个集合以外的元素没有影响。默认参数保留所有条目。
   *
   */
  void
  reinit(const size_type m,
         const size_type n,
         const IndexSet &rowset = IndexSet());

  /**
   * 由于这个对象无论如何都是保持压缩的，所以这个函数什么都不做，但声明这个函数是为了使这个类的接口与SparsityPattern类的接口一样。
   *
   */
  void
  compress();

  /**
   * 返回该对象是否为空。如果没有分配内存，它就是空的，这与两个维度都是零是一样的。
   *
   */
  bool
  empty() const;

  /**
   * 返回每行的最大条目数。注意，这个数字可能会随着条目的增加而改变。
   *
   */
  size_type
  max_entries_per_row() const;

  /**
   * 添加一个非零的条目。如果该条目已经存在，这个调用不做任何事情。
   *
   */
  void
  add(const size_type i, const size_type j);

  /**
   * 在指定的行中添加几个非零的条目。已经存在的条目会被忽略。
   *
   */
  template <typename ForwardIterator>
  void
  add_entries(const size_type row,
              ForwardIterator begin,
              ForwardIterator end,
              const bool      indices_are_unique_and_sorted = false);

  /**
   * 检查某个位置的值是否可能为非零。
   *
   */
  bool
  exists(const size_type i, const size_type j) const;

  /**
   * 返回该稀疏模式的视图。  也就是说，对于 @p rows
   * 中的所有行，提取非空列。
   * 得到的稀疏模式的行数将等于`rows.n_elements()`。
   *
   */
  DynamicSparsityPattern
  get_view(const IndexSet &rows) const;

  /**
   * 通过添加转置对象的稀疏度模式使稀疏度模式对称。
   * 如果稀疏模式不代表一个正方形矩阵，这个函数会抛出一个异常。
   *
   */
  void
  symmetrize();

  /**
   * 构建并在此对象中存储对应于 @p left 和 @p right
   * 稀疏模式的乘积的稀疏模式。
   *
   */
  template <typename SparsityPatternTypeLeft, typename SparsityPatternTypeRight>
  void
  compute_mmult_pattern(const SparsityPatternTypeLeft & left,
                        const SparsityPatternTypeRight &right);

  /**
   * 构造并在此对象中存储对应于转置 @p left 和非转置 @p
   * right 稀疏模式乘积的稀疏模式。
   *
   */
  template <typename SparsityPatternTypeLeft, typename SparsityPatternTypeRight>
  void
  compute_Tmmult_pattern(const SparsityPatternTypeLeft & left,
                         const SparsityPatternTypeRight &right);

  /**
   * 打印稀疏性模式。输出包括每行一行，格式为<tt>[i,j1,j2,j3,...]/tt>。<i>i</i>是行号，<i>jn</i>是这一行中分配的列。
   *
   */
  void
  print(std::ostream &out) const;

  /**
   * 以 @p gnuplot
   * 理解的格式打印稀疏模式，该格式可用于以图形方式绘制稀疏模式。该格式由成对的<tt>i
   * j</tt>非零元素组成，每个代表一个条目，在输出文件中每行一个。指数从零开始计算，如常。由于稀疏模式的打印方式与矩阵的显示方式相同，我们打印的是列索引的负数，这意味着<tt>(0,0)</tt>元素位于左上角而不是左下角。
   * 在gnuplot中通过将数据样式设置为点或点来打印稀疏模式，并使用
   * @p plot 命令。
   *
   */
  void
  print_gnuplot(std::ostream &out) const;

  /**
   * 返回行数，相当于图像空间的尺寸。
   *
   */
  size_type
  n_rows() const;

  /**
   * 返回列数，相当于范围空间的尺寸。
   *
   */
  size_type
  n_cols() const;

  /**
   * 特定行中的条目数。只有当给定的行是我们要存储的行的索引集的成员时，才能调用这个函数。
   *
   */
  size_type
  row_length(const size_type row) const;

  /**
   * 清除存储在某一特定行中的所有条目。
   *
   */
  void
  clear_row(const size_type row);

  /**
   * 访问列号字段。 返回 @p 中的第 @p row.
   * 个索引条目的列号。
   *
   */
  size_type
  column_number(const size_type row, const size_type index) const;

  /**
   * 返回 @p row. 行中 @p col 列的索引
   * 如果该列在此稀疏模式中不存在，返回值将是
   * 'numbers::invalid_size_type'.  。
   *
   */
  size_type
  column_index(const size_type row, const size_type col) const;

  /**
   * @name  迭代器
   *
   */
  // @{

  /**
   * 迭代器从矩阵的第一个条目开始。产生的迭代器可以用来走过疏散模式的所有非零条目。
   * 注意这个类的一般文档中关于元素被访问的顺序的讨论。
   * @note
   * 如果稀疏模式已经用IndexSet初始化，表示要存储哪些行，那么迭代器将简单地跳过没有存储的行。换句话说，它们看起来像空行，但是在迭代这些行的时候不会产生异常。
   *
   */
  iterator
  begin() const;

  /**
   * 最终的迭代器。
   *
   */
  iterator
  end() const;

  /**
   * 迭代器从<tt>r</tt>行的第一个条目开始。
   * 注意，如果给定的行是空的，即不包含任何非零条目，那么这个函数返回的迭代器就等于<tt>end(r)</tt>。还要注意的是，在这种情况下，迭代器可能是不能被解除引用的。
   * 还要注意这个类的一般文档中关于元素被访问的顺序的讨论。
   * @note
   * 如果稀疏模式已经用IndexSet初始化，表示要存储哪些行，那么迭代器将简单地跳过没有存储的行。换句话说，它们看起来像空行，但是在迭代这些行的时候不会产生异常。
   *
   */
  iterator
  begin(const size_type r) const;

  /**
   * 行<tt>r</tt>的最终迭代器。它指向超过 @p r,
   * 行末尾的第一个元素，或者超过整个稀疏模式的末尾。
   * 请注意，结束迭代器不一定是可被解除引用的。特别是如果它是一个矩阵的最后一行的结束迭代器，情况更是如此。
   *
   */
  iterator
  end(const size_type r) const;

  // @}

  /**
   * 计算这个结构所代表的矩阵的带宽。带宽是 $|i-j|$
   * 的最大值，其中索引对 $(i,j)$
   * 代表矩阵的一个非零条目。
   *
   */
  size_type
  bandwidth() const;

  /**
   * 返回通过该稀疏模式分配的非零元素的数量。
   *
   */
  size_type
  n_nonzero_elements() const;

  /**
   * 返回设定哪些行在当前处理器上是活动的IndexSet。它对应于在构造函数或reinit函数中给这个类的IndexSet。
   *
   */
  const IndexSet &
  row_index_set() const;

  /**
   * 返回包含所有列的条目的IndexSet，其中至少有一个元素存在于这个稀疏模式中。
   * @note 在并行情况下，这只考虑本地存储的行。
   *
   */
  IndexSet
  nonempty_cols() const;

  /**
   * 返回IndexSet，该索引集包含所有行的条目，其中至少有一个元素存在于该稀疏模式中。
   * @note 在并行情况下，这只考虑本地存储的行。
   *
   */
  IndexSet
  nonempty_rows() const;

  /**
   * 返回这个对象是否只存储那些显式添加的条目，或者稀疏模式是否包含在构建它时通过其他方式（隐式）添加的元素。对于当前的类，结果总是为真。
   * 这个函数的主要作用是在几种稀疏模式可以作为模板参数传递的情况下描述当前类。
   *
   */
  static bool
  stores_only_added_elements();

  /**
   * 确定这个对象的内存消耗（以字节为单位）的估计值。
   *
   */
  size_type
  memory_consumption() const;

private:
  /**
   * 一个标志，存储到目前为止是否有任何条目被添加。
   *
   */
  bool have_entries;

  /**
   * 该稀疏结构应代表的行数。
   *
   */
  size_type rows;

  /**
   * 这个稀疏结构所代表的列的数量。
   *
   */
  size_type cols;

  /**
   * 一个包含有效行的集合。
   *
   */

  IndexSet rowset;


  /**
   * 为每一行存储一些数据，描述这一行的哪些条目是非零的。数据被分类存储在
   * @p entries   std::vector.
   * 每行的向量在插入时是动态增长的，每次将其内存翻倍。
   *
   */
  struct Line
  {
  public:
    /**
     * 该行的列索引的存储。这个数组总是保持排序的。
     *
     */
    std::vector<size_type> entries;

    /**
     * 将给定的列号添加到这一行。
     *
     */
    void
    add(const size_type col_num);

    /**
     * 将迭代器范围指定的列添加到这一行。
     *
     */
    template <typename ForwardIterator>
    void
    add_entries(ForwardIterator begin,
                ForwardIterator end,
                const bool      indices_are_sorted);

    /**
     * 估计内存消耗。
     *
     */
    size_type
    memory_consumption() const;
  };


  /**
   * 实际数据：为每行存储非零条目的集合。
   *
   */
  std::vector<Line> lines;

  // make the accessor class a friend
  friend class DynamicSparsityPatternIterators::Accessor;
};

 /*@}*/ 
 /*---------------------- Inline functions -----------------------------------*/ 


namespace DynamicSparsityPatternIterators
{
  inline Accessor::Accessor(const DynamicSparsityPattern *sparsity_pattern,
                            const size_type               row,
                            const unsigned int            index_within_row)
    : sparsity_pattern(sparsity_pattern)
    , current_row(row)
    , current_entry(
        ((sparsity_pattern->rowset.size() == 0) ?
           sparsity_pattern->lines[current_row].entries.begin() :
           sparsity_pattern
             ->lines[sparsity_pattern->rowset.index_within_set(current_row)]
             .entries.begin()) +
        index_within_row)
    , end_of_row(
        (sparsity_pattern->rowset.size() == 0) ?
          sparsity_pattern->lines[current_row].entries.end() :
          sparsity_pattern
            ->lines[sparsity_pattern->rowset.index_within_set(current_row)]
            .entries.end())
  {
    AssertIndexRange(current_row, sparsity_pattern->n_rows());
    Assert((sparsity_pattern->rowset.size() == 0) ||
             sparsity_pattern->rowset.is_element(current_row),
           ExcMessage("You can't create an iterator into a "
                      "DynamicSparsityPattern's row that is not "
                      "actually stored by that sparsity pattern "
                      "based on the IndexSet argument to it."));
    AssertIndexRange(
      index_within_row,
      ((sparsity_pattern->rowset.size() == 0) ?
         sparsity_pattern->lines[current_row].entries.size() :
         sparsity_pattern
           ->lines[sparsity_pattern->rowset.index_within_set(current_row)]
           .entries.size()));
  }


  inline Accessor::Accessor(const DynamicSparsityPattern *sparsity_pattern)
    : sparsity_pattern(sparsity_pattern)
    , current_row(numbers::invalid_size_type)
    , current_entry()
    , end_of_row()
  {}



  inline Accessor::Accessor()
    : sparsity_pattern(nullptr)
    , current_row(numbers::invalid_size_type)
    , current_entry()
    , end_of_row()
  {}


  inline size_type
  Accessor::row() const
  {
    Assert(sparsity_pattern != nullptr, DummyAccessor());
    Assert(current_row < sparsity_pattern->n_rows(), ExcInternalError());

    return current_row;
  }


  inline size_type
  Accessor::column() const
  {
    Assert(sparsity_pattern != nullptr, DummyAccessor());
    Assert(current_row < sparsity_pattern->n_rows(), ExcInternalError());

    return *current_entry;
  }


  inline size_type
  Accessor::index() const
  {
    Assert(sparsity_pattern != nullptr, DummyAccessor());
    Assert(current_row < sparsity_pattern->n_rows(), ExcInternalError());

    return (current_entry -
            ((sparsity_pattern->rowset.size() == 0) ?
               sparsity_pattern->lines[current_row].entries.begin() :
               sparsity_pattern
                 ->lines[sparsity_pattern->rowset.index_within_set(current_row)]
                 .entries.begin()));
  }



  inline bool
  Accessor::operator==(const Accessor &other) const
  {
    Assert(sparsity_pattern != nullptr, DummyAccessor());
    Assert(other.sparsity_pattern != nullptr, DummyAccessor());
    // compare the sparsity pattern the iterator points into, the
    // current row, and the location within this row. ignore the
    // latter if the row is past-the-end because in that case the
    // current_entry field may not point to a deterministic location
    return (sparsity_pattern == other.sparsity_pattern &&
            current_row == other.current_row &&
            ((current_row == numbers::invalid_size_type) ||
             (current_entry == other.current_entry)));
  }



  inline bool
  Accessor::operator<(const Accessor &other) const
  {
    Assert(sparsity_pattern != nullptr, DummyAccessor());
    Assert(other.sparsity_pattern != nullptr, DummyAccessor());
    Assert(sparsity_pattern == other.sparsity_pattern, ExcInternalError());

    // if *this is past-the-end, then it is less than no one
    if (current_row == numbers::invalid_size_type)
      return (false);
    // now *this should be an valid value
    Assert(current_row < sparsity_pattern->n_rows(), ExcInternalError());

    // if other is past-the-end
    if (other.current_row == numbers::invalid_size_type)
      return (true);
    // now other should be an valid value
    Assert(other.current_row < sparsity_pattern->n_rows(), ExcInternalError());

    // both iterators are not one-past-the-end
    return ((current_row < other.current_row) ||
            ((current_row == other.current_row) &&
             (current_entry < other.current_entry)));
  }


  inline void
  Accessor::advance()
  {
    Assert(sparsity_pattern != nullptr, DummyAccessor());
    Assert(current_row < sparsity_pattern->n_rows(), ExcInternalError());

    // move to the next element in this row
    ++current_entry;

    // if this moves us beyond the end of the row, go to the next row
    // if possible, or set the iterator to an invalid state if not.
    //
    // going to the next row is a bit complicated because we may have
    // to skip over empty rows, and because we also have to avoid rows
    // that aren't listed in a possibly passed IndexSet argument of
    // the sparsity pattern. consequently, rather than trying to
    // duplicate code here, just call the begin() function of the
    // sparsity pattern itself
    if (current_entry == end_of_row)
      {
        if (current_row + 1 < sparsity_pattern->n_rows())
          *this = *sparsity_pattern->begin(current_row + 1);
        else
          *this = Accessor(sparsity_pattern); // invalid object
      }
  }



  inline Iterator::Iterator(const DynamicSparsityPattern *sparsity_pattern,
                            const size_type               row,
                            const unsigned int            index_within_row)
    : accessor(sparsity_pattern, row, index_within_row)
  {}



  inline Iterator::Iterator(const DynamicSparsityPattern *sparsity_pattern)
    : accessor(sparsity_pattern)
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
    return !(*this == other);
  }


  inline bool
  Iterator::operator<(const Iterator &other) const
  {
    return accessor < other.accessor;
  }


  inline int
  Iterator::operator-(const Iterator &other) const
  {
    (void)other;
    Assert(accessor.sparsity_pattern == other.accessor.sparsity_pattern,
           ExcInternalError());
    Assert(false, ExcNotImplemented());

    return 0;
  }
} // namespace DynamicSparsityPatternIterators


inline void
DynamicSparsityPattern::Line::add(const size_type j)
{
  // first check the last element (or if line is still empty)
  if ((entries.size() == 0) || (entries.back() < j))
    {
      entries.push_back(j);
      return;
    }

  // do a binary search to find the place where to insert:
  std::vector<size_type>::iterator it =
    Utilities::lower_bound(entries.begin(), entries.end(), j);

  // If this entry is a duplicate, exit immediately
  if (*it == j)
    return;

  // Insert at the right place in the vector. Vector grows automatically to
  // fit elements. Always doubles its size.
  entries.insert(it, j);
}



inline DynamicSparsityPattern::size_type
DynamicSparsityPattern::n_rows() const
{
  return rows;
}



inline types::global_dof_index
DynamicSparsityPattern::n_cols() const
{
  return cols;
}



inline void
DynamicSparsityPattern::add(const size_type i, const size_type j)
{
  AssertIndexRange(i, rows);
  AssertIndexRange(j, cols);

  if (rowset.size() > 0 && !rowset.is_element(i))
    return;

  have_entries = true;

  const size_type rowindex =
    rowset.size() == 0 ? i : rowset.index_within_set(i);
  lines[rowindex].add(j);
}



template <typename ForwardIterator>
inline void
DynamicSparsityPattern::add_entries(const size_type row,
                                    ForwardIterator begin,
                                    ForwardIterator end,
                                    const bool      indices_are_sorted)
{
  AssertIndexRange(row, rows);

  if (rowset.size() > 0 && !rowset.is_element(row))
    return;

  if (!have_entries && begin < end)
    have_entries = true;

  const size_type rowindex =
    rowset.size() == 0 ? row : rowset.index_within_set(row);
  lines[rowindex].add_entries(begin, end, indices_are_sorted);
}



inline types::global_dof_index
DynamicSparsityPattern::row_length(const size_type row) const
{
  AssertIndexRange(row, n_rows());

  if (!have_entries)
    return 0;

  if (rowset.size() > 0 && !rowset.is_element(row))
    return 0;

  const size_type rowindex =
    rowset.size() == 0 ? row : rowset.index_within_set(row);
  return lines[rowindex].entries.size();
}



inline types::global_dof_index
DynamicSparsityPattern::column_number(const size_type row,
                                      const size_type index) const
{
  AssertIndexRange(row, n_rows());
  Assert(rowset.size() == 0 || rowset.is_element(row), ExcInternalError());

  const size_type local_row =
    rowset.size() ? rowset.index_within_set(row) : row;
  AssertIndexRange(index, lines[local_row].entries.size());
  return lines[local_row].entries[index];
}



inline DynamicSparsityPattern::iterator
DynamicSparsityPattern::begin() const
{
  if (n_rows() > 0)
    return begin(0);
  else
    return end();
}


inline DynamicSparsityPattern::iterator
DynamicSparsityPattern::end() const
{
  return {this};
}



inline DynamicSparsityPattern::iterator
DynamicSparsityPattern::begin(const size_type r) const
{
  AssertIndexRange(r, n_rows());

  if (!have_entries)
    return {this};

  if (rowset.size() > 0)
    {
      // We have an IndexSet that describes the locally owned set. For
      // performance reasons we need to make sure that we don't do a
      // linear search over 0..n_rows(). Instead, find the first entry
      // >= row r in the locally owned set (this is done in log
      // n_ranges time inside at()). From there, we move forward until
      // we find a non-empty row. By iterating over the IndexSet instead
      // of incrementing the row index, we potentially skip over entries
      // not in the rowset.
      IndexSet::ElementIterator it = rowset.at(r);
      if (it == rowset.end())
        return end(); // we don't own any row between r and the end

      // Instead of using row_length(*it)==0 in the while loop below,
      // which involves an expensive index_within_set() call, we
      // look at the lines vector directly. This works, because we are
      // walking over this vector entry by entry anyways.
      size_type rowindex = rowset.index_within_set(*it);

      while (it != rowset.end() && lines[rowindex].entries.size() == 0)
        {
          ++it;
          ++rowindex;
        }

      if (it == rowset.end())
        return end();
      else
        return {this, *it, 0};
    }

  // Without an index set we have to do a linear search starting at
  // row r until we find a non-empty one. We will check the lines vector
  // directly instead of going through the slower row_length() function
  size_type row = r;

  while (row < n_rows() && lines[row].entries.size() == 0)
    {
      ++row;
    }

  if (row == n_rows())
    return {this};
  else
    return {this, row, 0};
}



inline DynamicSparsityPattern::iterator
DynamicSparsityPattern::end(const size_type r) const
{
  AssertIndexRange(r, n_rows());

  const size_type row = r + 1;
  if (row == n_rows())
    return {this};
  else
    return begin(row);
}



inline const IndexSet &
DynamicSparsityPattern::row_index_set() const
{
  return rowset;
}



inline bool
DynamicSparsityPattern::stores_only_added_elements()
{
  return true;
}


DEAL_II_NAMESPACE_CLOSE

#endif


