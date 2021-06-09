//include/deal.II-translator/base/index_set_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2021 by the deal.II authors
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

#ifndef dealii_index_set_h
#define dealii_index_set_h

#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/thread_management.h>
#include <deal.II/base/utilities.h>

#include <boost/serialization/vector.hpp>

#include <algorithm>
#include <iterator>
#include <vector>


#ifdef DEAL_II_WITH_TRILINOS
#  include <Epetra_Map.h>
#  ifdef DEAL_II_TRILINOS_WITH_TPETRA
#    include <Tpetra_Map.hpp>
#  endif
#endif

#if defined(DEAL_II_WITH_MPI) || defined(DEAL_II_WITH_PETSC)
#  include <mpi.h>
#else
using MPI_Comm = int;
#  ifndef MPI_COMM_WORLD
#    define MPI_COMM_WORLD 0
#  endif
#endif

DEAL_II_NAMESPACE_OPEN

/**
 * 一个代表较大集合中的一个索引子集的类。例如，它可以用来表示在
 * $[0,\mathrm{dof\_handler.n\_dofs()})$
 * 范围内属于某个特定子域的自由度集合，或者在分布式并行计算中存储在某个特定处理器上的所有自由度中的自由度。
 * 这个类可以表示半开范围的指数集合，也可以表示单个元素。为实用起见，它还存储了这些指数可以承担的总体范围。换句话说，你需要指定索引空间的大小
 * $[0,\text{size})$ ，这个类的对象是其中的一个子集。
 * 有两种方法可以在IndexSets上进行迭代。首先，begin()和end()允许对集合中的单个索引进行迭代。第二，begin_interval()和end_interval()允许对上述的半开范围进行迭代。
 * 这个类中使用的数据结构以及原理可以在 @ref distributed_paper "分布式计算论文 "
 * 中找到。
 *
 *
 */
class IndexSet
{
public:
  // forward declarations:
  class ElementIterator;
  class IntervalIterator;

  /**
   * @p size_type 是用于存储索引集的大小和单个条目的类型。
   *
   */
  using size_type = types::global_dof_index;

  /**
   * 我们可以把 IndexSet 看作是 size()
   * 的容器，其中容器的元素是 bool
   * 值，这些元素要么是假的，要么是真的，取决于某个特定的索引是否是
   * IndexSet 的一个元素。换句话说，IndexSet
   * 有点像一个矢量，其中我们存储的元素是布尔值。另一方面，
   * @p bool
   * 的缺点是它不是一个数字类型，例如，它不允许与 @p
   * double.
   * 相乘。换句话说，在一个允许使用其他向量的地方，我们不能轻易使用一个布尔的向量。因此，我们声明这样一个向量的元素的类型是有符号的整数。这是因为在C++语言中，booleans可以隐含地转换为整数。换句话说，将向量元素的类型声明为有符号整数只是一个小谎言，但它是一个有用的谎言。
   *
   */
  using value_type = signed int;


  /**
   * 默认构造函数。
   *
   */
  IndexSet();

  /**
   * 构造函数，也设置了索引范围的整体大小。
   *
   */
  explicit IndexSet(const size_type size);

  /**
   * 复制构造函数。
   *
   */
  IndexSet(const IndexSet &) = default;

  /**
   * 拷贝赋值操作符。
   *
   */
  IndexSet &
  operator=(const IndexSet &) = default;

  /**
   * 移动构造函数。通过转移输入集的内部数据来创建一个新的IndexSet。
   *
   */
  IndexSet(IndexSet &&is) noexcept;

  /**
   * 移动赋值运算符。将输入集的内部数据转移到当前集中。
   *
   */
  IndexSet &
  operator=(IndexSet &&is) noexcept;

#ifdef DEAL_II_WITH_TRILINOS
  /**
   * 来自Trilinos Epetra_BlockMap的构造函数。
   *
   */
  explicit IndexSet(const Epetra_BlockMap &map);
#endif

  /**
   * 从这个索引集中删除所有的索引。然而，索引集仍保留其大小。
   *
   */
  void
  clear();

  /**
   * 设置这个对象所操作的索引的最大尺寸。
   * 这个函数只有在索引集还不包含任何元素时才能被调用。
   * 这可以通过调用clear()来实现，比如说。
   *
   */
  void
  set_size(const size_type size);

  /**
   * 返回这个索引集是其子集的索引空间的大小。
   * 注意，这个结果不等于这个索引集内的索引数。后者的信息是由n_elements()返回的。
   *
   */
  size_type
  size() const;

  /**
   * 将半开范围 $[\text{begin},\text{end})$
   * 添加到该类所代表的索引集合中。    @param[in]  开始
   * 要添加的范围的第一个元素。    @param[in]  结束
   * 要添加的范围的过末元素。
   *
   */
  void
  add_range(const size_type begin, const size_type end);

  /**
   * 将一个单独的索引添加到索引集合中。
   *
   */
  void
  add_index(const size_type index);

  /**
   * 添加一整组索引，通过取消引用迭代器范围的每个元素来描述
   * <code>[begin,end)</code>  。      @param[in]  开始
   * 迭代器到要添加的索引范围的第一个元素  @param[in]
   * 结束 要添加的元素范围的过去-结束迭代器。  @pre
   * 需要满足条件 <code>begin@<=end</code> 。
   *
   */
  template <typename ForwardIterator>
  void
  add_indices(const ForwardIterator &begin, const ForwardIterator &end);

  /**
   * 将给定的索引集 @p other
   * 添加到当前的索引集，构建this和 @p other. 的联合体 如果
   * @p offset 参数为非零，那么 @p other
   * 中的每个索引在被添加到当前索引集之前都会被 @p offset
   * 移位。例如，这允许从其他几个应该代表不同范围的索引集中构建一个索引集（例如，当从一个向量的各个块的非零元素集中构建一个块向量的非零条目集时）。
   * 如果 @p other
   * 索引集的任何索引（可能是移位的）位于当前对象所代表的
   * <code>[0,size())</code>
   * 范围之外，这个函数将产生一个异常。
   *
   */
  void
  add_indices(const IndexSet &other, const size_type offset = 0);

  /**
   * 返回指定索引是否是索引集的一个元素。
   *
   */
  bool
  is_element(const size_type index) const;

  /**
   * 返回这个对象所存储的索引集是否定义了一个连续的范围。如果根本没有存储任何索引，这也是真的。
   *
   */
  bool
  is_contiguous() const;

  /**
   * 返回这个对象所存储的索引集是否不包含任何元素。
   * 这类似于，但比检查  <code>n_elements() == 0</code>  要快。
   *
   */
  bool
  is_empty() const;

  /**
   * 返回IndexSets是否相对于MPI进程号是升序的，并且是1:1，也就是说，每个索引正好包含在一个IndexSet中（在不同进程上存储的索引中），每个进程存储连续的索引子集，并且进程
   * $p+1$  上的索引集开始于比进程  $p$
   * 上存储的最后一个索引大的索引。
   * 在只有一个MPI进程的情况下，这仅仅意味着索引集已经完成。
   *
   */
  bool
  is_ascending_and_one_to_one(const MPI_Comm &communicator) const;

  /**
   * 返回存储在这个索引集中的元素的数量。
   *
   */
  size_type
  n_elements() const;

  /**
   * 返回存储在这个索引集中的编号为 @p local_index
   * 的局部索引的全局索引。  @p local_index
   * 显然需要小于n_elements()。
   *
   */
  size_type
  nth_index_in_set(const size_type local_index) const;

  /**
   * 返回这个集合的第多少个元素（按升序计数）  @p
   * global_index是。  @p global_index 需要小于size()。如果索引 @p
   * global_index
   * 实际上不是这个索引集的成员，即如果is_element(global_index)是假的，这个函数返回
   * numbers::invalid_dof_index 。
   *
   */
  size_type
  index_within_set(const size_type global_index) const;

  /**
   * 每个索引集可以表示为若干个连续的索引区间的联合，如果有必要，区间可以只由单个元素组成，以表示索引集的孤立成员。
   * 这个函数返回代表所考虑的索引集所需的最小数量的这种区间。
   *
   */
  unsigned int
  n_intervals() const;

  /**
   * 该函数返回最大区间开始的局部索引。
   * 换句话说，返回值是nth_index_in_set(x)，其中x是IndexSet中最大连续范围的第一个索引。因此，返回值等于集合中在最大范围之前的元素的数量。
   * 这个调用假定IndexSet是不空的。
   *
   */
  size_type
  largest_range_starting_index() const;

  /**
   * 通过合并具有连续范围的单个元素来压缩内部表示，等等。这个函数没有任何外部影响。
   *
   */
  void
  compress() const;

  /**
   * 对索引集的平等进行比较。这个操作只允许在两个集合的大小相同的情况下进行（当然它们不一定要有相同数量的索引）。
   *
   */
  bool
  operator==(const IndexSet &is) const;

  /**
   * 索引集不平等的比较。这个操作只允许在两个集合的大小相同的情况下进行（当然它们不一定要有相同的索引数）。
   *
   */
  bool
  operator!=(const IndexSet &is) const;

  /**
   * 返回当前索引集和所给参数的交集，即一组索引是两个索引集的元素。这两个索引集必须具有相同的大小（当然，它们不一定具有相同的索引数）。
   *
   */
  IndexSet operator&(const IndexSet &is) const;

  /**
   * 这个命令接收一个区间<tt>[begin,
   * end)</tt>，并返回当前索引集与该区间的交集，移到<tt>[0,
   * end-begin)</tt>范围。
   * 换句话说，这个操作的结果是当前对象所代表的集合与区间<tt>[begin,
   * end)</tt>的交集，如<i>within the interval <tt>[begin,
   * end)</tt></i>所见，将交集操作的结果向左移动<tt>begin</tt>。这对应于一个<i>view</i>的概念。区间<tt>[begin,
   * end)</tt>是一个<i>window</i>，通过它我们可以看到当前对象所代表的集合。
   *
   */
  IndexSet
  get_view(const size_type begin, const size_type end) const;

  /**
   * 将此对象所代表的集合指数分割成由 @p n_indices_per_block
   * 结构给出的块。其条目之和必须与当前对象的全局大小一致。
   *
   */
  std::vector<IndexSet>
  split_by_block(
    const std::vector<types::global_dof_index> &n_indices_per_block) const;

  /**
   * 从这个集合中移除 @p other
   * 中包含的所有元素。换句话说，如果 $x$ 是当前对象，
   * $o$ 是参数，那么我们计算 $x \leftarrow x \backslash o$  。
   *
   */
  void
  subtract_set(const IndexSet &other);

  /**
   * 返回一个新的IndexSet，全局大小等于`this->size()*other.size()`，对于这个IndexSet的每个元素`n`，包含
   * @p other IndexSet的半开放范围`[n*other.size(),
   * (n+1)*other.size()]中的条目。
   * 这个名字来自于这样的观点：我们从一个IndexSet开始，与另一个有`other.size()`元素的IndexSet进行张量乘积；这将产生一个大小为`this->size()`乘以`other.size()`的矩阵，在这个IndexSet包含一个索引的行和
   * @p other
   * IndexSet包含一个索引的列中有1。然后，这个矩阵被再次
   * "解卷"，逐一浏览每一行，并以连续的顺序对矩阵的条目重新进行索引。矩阵中的
   * "1
   * "对应于重新索引的IndexSet中的一个条目，这个函数将返回该条目。
   *
   */
  IndexSet
  tensor_product(const IndexSet &other) const;

  /**
   * 删除并返回最后一个范围的最后一个元素。
   * 如果IndexSet是空的，这个函数会抛出一个异常。
   *
   */
  size_type
  pop_back();

  /**
   * 删除并返回第一个区域的第一个元素。
   * 如果IndexSet是空的，这个函数会抛出一个异常。
   *
   */
  size_type
  pop_front();

  /**
   * 用这个IndexSet中包含的所有索引填充给定的向量。
   *
   */
  void
  fill_index_vector(std::vector<size_type> &indices) const;

  /**
   * 用零或一的元素填充给定的向量，提供这个索引集的二进制表示。假设给定的向量已经有了正确的大小。
   * 给定的参数用整数值0和1填充，使用
   * <code>vector.operator[]</code>
   * 。因此，只要允许将零和一的整数转换为向量的元素，任何具有这种操作符的对象都可以被使用。具体来说，对于类Vector、BlockVector，以及
   * std::vector@<bool@>,   std::vector@<int@>, 和 std::vector@<double@>.
   * 都是这样的情况。
   *
   */
  template <typename VectorType>
  void
  fill_binary_vector(VectorType &vector) const;

  /**
   * 将此IndexSet的文本表示法输出到给定的流中。用于测试。
   *
   */
  template <class StreamType>
  void
  print(StreamType &out) const;

  /**
   * 将IndexSet写成基于文本的文件格式，可以用read()函数再次读入。
   *
   */
  void
  write(std::ostream &out) const;

  /**
   * 从由write()函数写入的流 @p in
   * 给出的基于文本的表示法中构建IndexSet。
   *
   */
  void
  read(std::istream &in);

  /**
   * 将 IndexSet 写入一个二进制的、紧凑的表示法，可以使用
   * block_read() 函数再次读入。
   *
   */
  void
  block_write(std::ostream &out) const;

  /**
   * 从由write_block()函数写入的流 @p in
   * 给出的二进制表示法中构造IndexSet。
   *
   */
  void
  block_read(std::istream &in);

#ifdef DEAL_II_WITH_TRILINOS
  /**
   * 给定一个MPI通信器，创建一个Trilinos
   * map对象，该对象表示向量元素或矩阵行的分布，我们将在其中本地存储那些我们在当前索引集中存储索引的元素或行，而所有其他元素/行则在其他MPI进程之一上的其他地方。
   * 最后一个参数只有在通信器是一个并行的通信器，将计算分布在多个处理器上时才起作用。在这种情况下，如果最后一个参数是假的，那么就假定在所有处理器上调用这个函数的索引集是互斥的，但一起列举每个索引时正好是一次。换句话说，如果你在两个处理器上调用这个函数，那么这个函数被调用的索引集必须一起拥有从0到size()-1的所有可能的索引，并且没有索引必须出现在两个索引集中。例如，这相当于我们要把向量的元素分成独特的子集，存储在不同的处理器上的情况
   *
   * - 任何元素都不应被一个以上的处理器所拥有，但每个元素必须被一个处理器所拥有。    另一方面，如果第二个参数为真，那么索引集可以是重叠的，它们也不需要横跨整个索引集。如果我们想创建的向量不仅包含本地拥有的索引，而且还包含对应于位于幽灵单元上的自由度的元素，那么这是一个有用的操作。这个方法的另一个应用是选择一个向量的元素子集，例如，只提取某些解的成分。
   *
   */
  Epetra_Map
  make_trilinos_map(const MPI_Comm &communicator = MPI_COMM_WORLD,
                    const bool      overlapping  = false) const;

#  ifdef DEAL_II_TRILINOS_WITH_TPETRA
  Tpetra::Map<int, types::global_dof_index>
  make_tpetra_map(const MPI_Comm &communicator = MPI_COMM_WORLD,
                  const bool      overlapping  = false) const;
#  endif
#endif


  /**
   * 确定这个对象的内存消耗（以字节为单位）的估计值。
   *
   */
  std::size_t
  memory_consumption() const;

  DeclException1(ExcIndexNotPresent,
                 size_type,
                 << "The global index " << arg1
                 << " is not an element of this set.");

  /**
   * 使用[BOOST序列化库](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html)将此对象的数据写入或读出到一个流中，以便进行序列化。
   *
   */
  template <class Archive>
  void
  serialize(Archive &ar, const unsigned int version);


  /**
   * @name  迭代器  
     * @{ 
   *
   */

  /**
   * 解除对IntervalIterator的引用将返回对该类型对象的引用。它允许访问被迭代的IndexSet的一个连续的区间
   * $[a,b[$ （也称为范围）。
   *
   */
  class IntervalAccessor
  {
  public:
    /**
     * 给出一个IndexSet和要指向的范围的索引 @p range_idx
     * ，构造一个有效的访问器。
     *
     */
    IntervalAccessor(const IndexSet *idxset, const size_type range_idx);

    /**
     * 为IndexSet构造一个无效的访问器。
     *
     */
    explicit IntervalAccessor(const IndexSet *idxset);

    /**
     * 这个区间内的元素数量。
     *
     */
    size_type
    n_elements() const;

    /**
     * 如果为真，我们就指向IndexSet中的一个有效区间。
     *
     */
    bool
    is_valid() const;

    /**
     * 返回一个迭代器，指向这个区间的第一个索引。
     *
     */
    ElementIterator
    begin() const;

    /**
     * 返回一个直接指向该区间最后一个索引之后的迭代器。
     *
     */
    ElementIterator
    end() const;

    /**
     * 返回这个区间的最后一个索引的索引。
     *
     */
    size_type
    last() const;

  private:
    /**
     * 私有的复制构造函数。
     *
     */
    IntervalAccessor(const IntervalAccessor &other);
    /**
     * 私有的复制操作符。
     *
     */
    IntervalAccessor &
    operator=(const IntervalAccessor &other);

    /**
     * 平等测试，由IntervalIterator使用。
     *
     */
    bool
    operator==(const IntervalAccessor &other) const;
    /**
     * 小于运算符，由IntervalIterator使用。
     *
     */
    bool
    operator<(const IntervalAccessor &other) const;
    /**
     * 推进这个访问器，指向 @p index_set中的下一个区间。
     *
     */
    void
    advance();
    /**
     * 对IndexSet的引用。
     *
     */
    const IndexSet *index_set;

    /**
     * 进入index_set. ranges[]的索引。如果无效，则设置为
     * numbers::invalid_dof_index ，或者结束迭代器。
     *
     */
    size_type range_idx;

    friend class IntervalIterator;
  };

  /**
   * 表示指向由 IndexSet::begin_interval(). 返回的连续区间 $[a,b[$
   * 的迭代器的类。
   *
   */
  class IntervalIterator
  {
  public:
    /**
     * 构建一个有效的迭代器，指向索引为 @p
     * range_idx的区间。
     *
     */
    IntervalIterator(const IndexSet *idxset, const size_type range_idx);

    /**
     * 构建一个无效的迭代器（作为end()使用）。
     *
     */
    explicit IntervalIterator(const IndexSet *idxset);

    /**
     * 构建一个空的迭代器。
     *
     */
    IntervalIterator();

    /**
     * 从 @p other 迭代器复制构造函数。
     *
     */
    IntervalIterator(const IntervalIterator &other) = default;

    /**
     * 对另一个迭代器进行赋值。
     *
     */
    IntervalIterator &
    operator=(const IntervalIterator &other) = default;

    /**
     * 前缀增量。
     *
     */
    IntervalIterator &
    operator++();

    /**
     * 后缀增量。
     *
     */
    IntervalIterator
    operator++(int);

    /**
     * 撤消运算符，返回一个IntervalAccessor。
     *
     */
    const IntervalAccessor &operator*() const;

    /**
     * 去引用操作符，返回一个指向IntervalAccessor的指针。
     *
     */
    const IntervalAccessor *operator->() const;

    /**
     * 比较。
     *
     */
    bool
    operator==(const IntervalIterator &) const;

    /**
     * <tt>==</tt>的倒数。
     *
     */
    bool
    operator!=(const IntervalIterator &) const;

    /**
     * 比较运算符。
     *
     */
    bool
    operator<(const IntervalIterator &) const;

    /**
     * 返回当前迭代器与参数之间的距离。这个距离是通过对当前迭代器应用operator++多少次才能得到参数（对于一个正的返回值），或operator--（对于一个负的返回值）。
     *
     */
    int
    operator-(const IntervalIterator &p) const;

    /**
     * 将该类标记为前向迭代器，并声明一些别名，这些别名是迭代器的标准，被算法用来查询它们所工作的迭代器的具体内容。
     *
     */
    using iterator_category = std::forward_iterator_tag;
    using value_type        = IntervalAccessor;
    using difference_type   = std::ptrdiff_t;
    using pointer           = IntervalAccessor *;
    using reference         = IntervalAccessor &;

  private:
    /**
     * 访问器，包含了我们所指向的IndexSet和间隔的内容。
     *
     */
    IntervalAccessor accessor;
  };

  /**
   * 表示指向IndexSet中单个元素的迭代器的类，如
   * IndexSet::begin(). 返回的那样。
   *
   */
  class ElementIterator
  {
  public:
    /**
     * 在区间 @p range_idx 内构造一个指向全局索引 @p index
     * 的迭代器。
     *
     */
    ElementIterator(const IndexSet *idxset,
                    const size_type range_idx,
                    const size_type index);

    /**
     * 构建一个指向IndexSet末端的迭代器。
     *
     */
    explicit ElementIterator(const IndexSet *idxset);

    /**
     * 去引用操作符。返回值是IndexSet中元素的索引。
     *
     */
    size_type operator*() const;

    /**
     * 这个迭代器是否指向一个现有的元素？
     *
     */
    bool
    is_valid() const;

    /**
     * 前缀增量。
     *
     */
    ElementIterator &
    operator++();

    /**
     * 后缀增量。
     *
     */
    ElementIterator
    operator++(int);

    /**
     * 比较。
     *
     */
    bool
    operator==(const ElementIterator &) const;

    /**
     * <tt>==</tt>的倒数。
     *
     */
    bool
    operator!=(const ElementIterator &) const;

    /**
     * 比较运算符。
     *
     */
    bool
    operator<(const ElementIterator &) const;

    /**
     * 返回当前迭代器与参数之间的距离。在表达式
     * <code>it_left-it_right</code> 中，距离是由对右操作数 @p
     * it_right应用operator++多少次才能得到左操作数 @p it_left
     * （对于正返回值），或者对 @p it_left
     * 应用多少次才能得到 @p it_right
     * （对于负返回值）而给出。
     *
     */
    std::ptrdiff_t
    operator-(const ElementIterator &p) const;

    /**
     * 将该类标记为前向迭代器，并声明一些别名，这些别名是迭代器的标准，被算法用来查询它们所工作的迭代器的具体信息。
     *
     */
    using iterator_category = std::forward_iterator_tag;
    using value_type        = size_type;
    using difference_type   = std::ptrdiff_t;
    using pointer           = size_type *;
    using reference         = size_type &;

  private:
    /**
     * 将迭代器提前一个。
     *
     */
    void
    advance();

    /**
     * 父索引集。
     *
     */
    const IndexSet *index_set;
    /**
     * 索引到index_set. ranges。
     *
     */
    size_type range_idx;
    /**
     * 这个迭代器所指向的全局索引。
     *
     */
    size_type idx;
  };

  /**
   * 返回一个迭代器，指向这个IndexSet中包含的第一个索引。
   *
   */
  ElementIterator
  begin() const;

  /**
   * 返回一个指向全局索引 @p global_index
   * 的元素的迭代器，如果该索引不在集合中，则指向下一个更大的元素。这等同于
   * @code
   * auto p = begin();
   * while (*p<global_index)
   * ++p;
   * return p;
   * @endcode
   * 如果在这个IndexSet中没有位于 @p global_index,
   * 或后面的元素，这个方法将返回end()。
   *
   */
  ElementIterator
  at(const size_type global_index) const;

  /**
   * 返回一个迭代器，该迭代器指向这个IndexSet中包含的最后一个索引之后。
   *
   */
  ElementIterator
  end() const;

  /**
   * 返回一个迭代器，该迭代器指向此IndexSet的第一个区间。
   *
   */
  IntervalIterator
  begin_intervals() const;

  /**
   * 返回一个指向此IndexSet的最后一个区间后的迭代器。
   *
   */
  IntervalIterator
  end_intervals() const;

  /**
   * @}
   *
   */

private:
  /**
   * 一个表示半开索引范围的类型  <code>[begin,end)</code>  。
   * nth_index_in_set表示当前范围的第一个元素在这个IndexSet中的第几个索引。只有当
   * IndexSet::compress()
   * 在最后一次插入后被调用，这个信息才是准确的。
   *
   */
  struct Range
  {
    size_type begin;
    size_type end;

    size_type nth_index_in_set;

    /**
     * 默认构造函数。由于默认构造的区间没有有用的选择，这个构造函数简单地创建了类似于一个无效区间的东西。我们需要这个构造函数来进行序列化，但是无效区间在使用之前应该用从存档中读取的东西来填充，所以我们应该希望永远不会在野外看到一个无效区间。
     *
     */
    Range();

    /**
     * 构造函数。用给定的指数创建一个半开放的区间。
     * @param  i1 区间的左端点。      @param  i2
     * 大于指定范围最后一个索引的第一个索引。
     *
     */
    Range(const size_type i1, const size_type i2);

    friend inline bool
    operator<(const Range &range_1, const Range &range_2)
    {
      return (
        (range_1.begin < range_2.begin) ||
        ((range_1.begin == range_2.begin) && (range_1.end < range_2.end)));
    }

    static bool
    end_compare(const IndexSet::Range &x, const IndexSet::Range &y)
    {
      return x.end < y.end;
    }

    static bool
    nth_index_compare(const IndexSet::Range &x, const IndexSet::Range &y)
    {
      return (x.nth_index_in_set + (x.end - x.begin) <
              y.nth_index_in_set + (y.end - y.begin));
    }

    friend inline bool
    operator==(const Range &range_1, const Range &range_2)
    {
      return ((range_1.begin == range_2.begin) && (range_1.end == range_2.end));
    }

    static std::size_t
    memory_consumption()
    {
      return sizeof(Range);
    }

    /**
     * 为了使用[BOOST序列化库](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html)进行序列化，将此对象的数据写入或读出到一个流中。
     *
     */
    template <class Archive>
    void
    serialize(Archive &ar, const unsigned int version);
  };

  /**
   * 一组连续的索引范围，构成了这个索引集的（一部分）。这个变量总是保持排序。
   * 这个变量被标记为
   * "可变"，因此它可以被compress()改变，当然这并不改变这个索引集的外部表示。
   *
   */
  mutable std::vector<Range> ranges;

  /**
   * 如果compress()在索引集的最后一次改变后被调用，则为真。
   * 该变量被标记为
   * "可变"，因此它可以被compress()改变，当然这并不改变这个索引集的外部表示。
   *
   */
  mutable bool is_compressed;

  /**
   * 索引范围的总体大小。这个索引集的元素必须有一个比这个值小的数字。
   *
   */
  size_type index_space_size;

  /**
   * 这个整数缓存了 @p ranges.
   * 中最大范围的索引，这给了<tt>O(1)</tt>对具有最多元素的范围的访问，而一般的访问成本<tt>O(log(n_ranges))</tt>。最大的范围需要用于
   * @p is_element(),   @p index_within_set(),   @p nth_index_in_set.
   * 的方法。在许多应用中，最大的范围包含大多数元素（本地拥有的范围），而只有少数其他元素（幽灵）。
   *
   */
  mutable size_type largest_range;

  /**
   * 一个mutex，用于同步do_compress()函数的操作，该函数通过compress()从许多'const'函数调用。
   *
   */
  mutable Threads::Mutex compress_mutex;

  /**
   * 实际执行compress()操作。
   *
   */
  void
  do_compress() const;
};


/**
 * 创建并返回一个大小为 $N$
 * 的索引集，其中包含这个范围内的每一个索引。从本质上讲，这个函数返回一个由以下方式创建的索引集
 *
 * @code
 * IndexSet is (N);
 * is.add_range(0, N);
 * @endcode
 * 这个函数的存在是为了让人们能够在一个步骤中创建和初始化完整的索引集，或者说是为了让人们能够写出这样的代码
 *
 * @code
 * if (my_index_set == complete_index_set(my_index_set.size())
 *   ...
 * @endcode
 * @relatesalso  IndexSet
 *
 *
 */
inline IndexSet
complete_index_set(const IndexSet::size_type N)
{
  IndexSet is(N);
  is.add_range(0, N);
  is.compress();
  return is;
}

 /* ------------------ inline functions ------------------ */ 


 /* IntervalAccessor */ 

inline IndexSet::IntervalAccessor::IntervalAccessor(
  const IndexSet *          idxset,
  const IndexSet::size_type range_idx)
  : index_set(idxset)
  , range_idx(range_idx)
{
  Assert(range_idx < idxset->n_intervals(),
         ExcInternalError("Invalid range index"));
}



inline IndexSet::IntervalAccessor::IntervalAccessor(const IndexSet *idxset)
  : index_set(idxset)
  , range_idx(numbers::invalid_dof_index)
{}



inline IndexSet::IntervalAccessor::IntervalAccessor(
  const IndexSet::IntervalAccessor &other)
  : index_set(other.index_set)
  , range_idx(other.range_idx)
{
  Assert(range_idx == numbers::invalid_dof_index || is_valid(),
         ExcMessage("invalid iterator"));
}



inline IndexSet::size_type
IndexSet::IntervalAccessor::n_elements() const
{
  Assert(is_valid(), ExcMessage("invalid iterator"));
  return index_set->ranges[range_idx].end - index_set->ranges[range_idx].begin;
}



inline bool
IndexSet::IntervalAccessor::is_valid() const
{
  return index_set != nullptr && range_idx < index_set->n_intervals();
}



inline IndexSet::ElementIterator
IndexSet::IntervalAccessor::begin() const
{
  Assert(is_valid(), ExcMessage("invalid iterator"));
  return {index_set, range_idx, index_set->ranges[range_idx].begin};
}



inline IndexSet::ElementIterator
IndexSet::IntervalAccessor::end() const
{
  Assert(is_valid(), ExcMessage("invalid iterator"));

  // point to first index in next interval unless we are the last interval.
  if (range_idx < index_set->ranges.size() - 1)
    return {index_set, range_idx + 1, index_set->ranges[range_idx + 1].begin};
  else
    return index_set->end();
}



inline IndexSet::size_type
IndexSet::IntervalAccessor::last() const
{
  Assert(is_valid(), ExcMessage("invalid iterator"));

  return index_set->ranges[range_idx].end - 1;
}



inline IndexSet::IntervalAccessor &
IndexSet::IntervalAccessor::operator=(const IndexSet::IntervalAccessor &other)
{
  index_set = other.index_set;
  range_idx = other.range_idx;
  Assert(range_idx == numbers::invalid_dof_index || is_valid(),
         ExcMessage("invalid iterator"));
  return *this;
}



inline bool
IndexSet::IntervalAccessor::
operator==(const IndexSet::IntervalAccessor &other) const
{
  Assert(index_set == other.index_set,
         ExcMessage(
           "Can not compare accessors pointing to different IndexSets"));
  return range_idx == other.range_idx;
}



inline bool
IndexSet::IntervalAccessor::
operator<(const IndexSet::IntervalAccessor &other) const
{
  Assert(index_set == other.index_set,
         ExcMessage(
           "Can not compare accessors pointing to different IndexSets"));
  return range_idx < other.range_idx;
}



inline void
IndexSet::IntervalAccessor::advance()
{
  Assert(
    is_valid(),
    ExcMessage(
      "Impossible to advance an IndexSet::IntervalIterator that is invalid"));
  ++range_idx;

  // set ourselves to invalid if we walk off the end
  if (range_idx >= index_set->ranges.size())
    range_idx = numbers::invalid_dof_index;
}


 /* IntervalIterator */ 

inline IndexSet::IntervalIterator::IntervalIterator(
  const IndexSet *          idxset,
  const IndexSet::size_type range_idx)
  : accessor(idxset, range_idx)
{}



inline IndexSet::IntervalIterator::IntervalIterator()
  : accessor(nullptr)
{}



inline IndexSet::IntervalIterator::IntervalIterator(const IndexSet *idxset)
  : accessor(idxset)
{}



inline IndexSet::IntervalIterator &
IndexSet::IntervalIterator::operator++()
{
  accessor.advance();
  return *this;
}



inline IndexSet::IntervalIterator
IndexSet::IntervalIterator::operator++(int)
{
  const IndexSet::IntervalIterator iter = *this;
  accessor.advance();
  return iter;
}



inline const IndexSet::IntervalAccessor &IndexSet::IntervalIterator::
                                         operator*() const
{
  return accessor;
}



inline const IndexSet::IntervalAccessor *IndexSet::IntervalIterator::
                                         operator->() const
{
  return &accessor;
}



inline bool
IndexSet::IntervalIterator::
operator==(const IndexSet::IntervalIterator &other) const
{
  return accessor == other.accessor;
}



inline bool
IndexSet::IntervalIterator::
operator!=(const IndexSet::IntervalIterator &other) const
{
  return !(*this == other);
}



inline bool
IndexSet::IntervalIterator::
operator<(const IndexSet::IntervalIterator &other) const
{
  return accessor < other.accessor;
}



inline int
IndexSet::IntervalIterator::
operator-(const IndexSet::IntervalIterator &other) const
{
  Assert(accessor.index_set == other.accessor.index_set,
         ExcMessage(
           "Can not compare iterators belonging to different IndexSets"));

  const size_type lhs = (accessor.range_idx == numbers::invalid_dof_index) ?
                          accessor.index_set->ranges.size() :
                          accessor.range_idx;
  const size_type rhs =
    (other.accessor.range_idx == numbers::invalid_dof_index) ?
      accessor.index_set->ranges.size() :
      other.accessor.range_idx;

  if (lhs > rhs)
    return static_cast<int>(lhs - rhs);
  else
    return -static_cast<int>(rhs - lhs);
}



 /* ElementIterator */ 

inline IndexSet::ElementIterator::ElementIterator(
  const IndexSet *          idxset,
  const IndexSet::size_type range_idx,
  const IndexSet::size_type index)
  : index_set(idxset)
  , range_idx(range_idx)
  , idx(index)
{
  Assert(range_idx < index_set->ranges.size(),
         ExcMessage(
           "Invalid range index for IndexSet::ElementIterator constructor."));
  Assert(
    idx >= index_set->ranges[range_idx].begin &&
      idx < index_set->ranges[range_idx].end,
    ExcInternalError(
      "Invalid index argument for IndexSet::ElementIterator constructor."));
}



inline IndexSet::ElementIterator::ElementIterator(const IndexSet *idxset)
  : index_set(idxset)
  , range_idx(numbers::invalid_dof_index)
  , idx(numbers::invalid_dof_index)
{}



inline bool
IndexSet::ElementIterator::is_valid() const
{
  Assert((range_idx == numbers::invalid_dof_index &&
          idx == numbers::invalid_dof_index) ||
           (range_idx < index_set->ranges.size() &&
            idx < index_set->ranges[range_idx].end),
         ExcInternalError("Invalid ElementIterator state."));

  return (range_idx < index_set->ranges.size() &&
          idx < index_set->ranges[range_idx].end);
}



inline IndexSet::size_type IndexSet::ElementIterator::operator*() const
{
  Assert(
    is_valid(),
    ExcMessage(
      "Impossible to dereference an IndexSet::ElementIterator that is invalid"));
  return idx;
}



inline bool
IndexSet::ElementIterator::
operator==(const IndexSet::ElementIterator &other) const
{
  Assert(index_set == other.index_set,
         ExcMessage(
           "Can not compare iterators belonging to different IndexSets"));
  return range_idx == other.range_idx && idx == other.idx;
}



inline void
IndexSet::ElementIterator::advance()
{
  Assert(
    is_valid(),
    ExcMessage(
      "Impossible to advance an IndexSet::ElementIterator that is invalid"));
  if (idx < index_set->ranges[range_idx].end)
    ++idx;
  // end of this range?
  if (idx == index_set->ranges[range_idx].end)
    {
      // point to first element in next interval if possible
      if (range_idx < index_set->ranges.size() - 1)
        {
          ++range_idx;
          idx = index_set->ranges[range_idx].begin;
        }
      else
        {
          // we just fell off the end, set to invalid:
          range_idx = numbers::invalid_dof_index;
          idx       = numbers::invalid_dof_index;
        }
    }
}



inline IndexSet::ElementIterator &
IndexSet::ElementIterator::operator++()
{
  advance();
  return *this;
}



inline IndexSet::ElementIterator
IndexSet::ElementIterator::operator++(int)
{
  const IndexSet::ElementIterator it = *this;
  advance();
  return it;
}



inline bool
IndexSet::ElementIterator::
operator!=(const IndexSet::ElementIterator &other) const
{
  return !(*this == other);
}



inline bool
IndexSet::ElementIterator::
operator<(const IndexSet::ElementIterator &other) const
{
  Assert(index_set == other.index_set,
         ExcMessage(
           "Can not compare iterators belonging to different IndexSets"));
  return range_idx < other.range_idx ||
         (range_idx == other.range_idx && idx < other.idx);
}



inline std::ptrdiff_t
IndexSet::ElementIterator::
operator-(const IndexSet::ElementIterator &other) const
{
  Assert(index_set == other.index_set,
         ExcMessage(
           "Can not compare iterators belonging to different IndexSets"));
  if (*this == other)
    return 0;
  if (!(*this < other))
    return -(other - *this);

  // only other can be equal to end() because of the checks above.
  Assert(is_valid(), ExcInternalError());

  // Note: we now compute how far advance *this in "*this < other" to get other,
  // so we need to return -c at the end.

  // first finish the current range:
  std::ptrdiff_t c = index_set->ranges[range_idx].end - idx;

  // now walk in steps of ranges (need to start one behind our current one):
  for (size_type range = range_idx + 1;
       range < index_set->ranges.size() && range <= other.range_idx;
       ++range)
    c += index_set->ranges[range].end - index_set->ranges[range].begin;

  Assert(
    other.range_idx < index_set->ranges.size() ||
      other.range_idx == numbers::invalid_dof_index,
    ExcMessage(
      "Inconsistent iterator state. Did you invalidate iterators by modifying the IndexSet?"));

  // We might have walked too far because we went until the end of
  // other.range_idx, so walk backwards to other.idx:
  if (other.range_idx != numbers::invalid_dof_index)
    c -= index_set->ranges[other.range_idx].end - other.idx;

  return -c;
}


 /* Range */ 

inline IndexSet::Range::Range()
  : begin(numbers::invalid_dof_index)
  , end(numbers::invalid_dof_index)
  , nth_index_in_set(numbers::invalid_dof_index)
{}



inline IndexSet::Range::Range(const size_type i1, const size_type i2)
  : begin(i1)
  , end(i2)
  , nth_index_in_set(numbers::invalid_dof_index)
{}



 /* IndexSet itself */ 

inline IndexSet::IndexSet()
  : is_compressed(true)
  , index_space_size(0)
  , largest_range(numbers::invalid_unsigned_int)
{}



inline IndexSet::IndexSet(const size_type size)
  : is_compressed(true)
  , index_space_size(size)
  , largest_range(numbers::invalid_unsigned_int)
{}



inline IndexSet::IndexSet(IndexSet &&is) noexcept
  : ranges(std::move(is.ranges))
  , is_compressed(is.is_compressed)
  , index_space_size(is.index_space_size)
  , largest_range(is.largest_range)
{
  is.ranges.clear();
  is.is_compressed    = true;
  is.index_space_size = 0;
  is.largest_range    = numbers::invalid_unsigned_int;

  compress();
}



inline IndexSet &
IndexSet::operator=(IndexSet &&is) noexcept
{
  ranges           = std::move(is.ranges);
  is_compressed    = is.is_compressed;
  index_space_size = is.index_space_size;
  largest_range    = is.largest_range;

  is.ranges.clear();
  is.is_compressed    = true;
  is.index_space_size = 0;
  is.largest_range    = numbers::invalid_unsigned_int;

  compress();

  return *this;
}



inline IndexSet::ElementIterator
IndexSet::begin() const
{
  compress();
  if (ranges.size() > 0)
    return {this, 0, ranges[0].begin};
  else
    return end();
}



inline IndexSet::ElementIterator
IndexSet::at(const size_type global_index) const
{
  compress();
  AssertIndexRange(global_index, size());

  if (ranges.empty())
    return end();

  std::vector<Range>::const_iterator main_range =
    ranges.begin() + largest_range;

  Range r(global_index, global_index + 1);
  // This optimization makes the bounds for lower_bound smaller by checking
  // the largest range first.
  std::vector<Range>::const_iterator range_begin, range_end;
  if (global_index < main_range->begin)
    {
      range_begin = ranges.begin();
      range_end   = main_range;
    }
  else
    {
      range_begin = main_range;
      range_end   = ranges.end();
    }

  // This will give us the first range p=[a,b[ with b>=global_index using
  // a binary search
  const std::vector<Range>::const_iterator p =
    Utilities::lower_bound(range_begin, range_end, r, Range::end_compare);

  // We couldn't find a range, which means we have no range that contains
  // global_index and also no range behind it, meaning we need to return end().
  if (p == ranges.end())
    return end();

  // Finally, we can have two cases: Either global_index is not in [a,b[,
  // which means we need to return an iterator to a because global_index, ...,
  // a-1 is not in the IndexSet (if branch). Alternatively, global_index is in
  // [a,b[ and we will return an iterator pointing directly at global_index
  // (else branch).
  if (global_index < p->begin)
    return {this, static_cast<size_type>(p - ranges.begin()), p->begin};
  else
    return {this, static_cast<size_type>(p - ranges.begin()), global_index};
}



inline IndexSet::ElementIterator
IndexSet::end() const
{
  compress();
  return IndexSet::ElementIterator(this);
}



inline IndexSet::IntervalIterator
IndexSet::begin_intervals() const
{
  compress();
  if (ranges.size() > 0)
    return IndexSet::IntervalIterator(this, 0);
  else
    return end_intervals();
}



inline IndexSet::IntervalIterator
IndexSet::end_intervals() const
{
  compress();
  return IndexSet::IntervalIterator(this);
}



inline void
IndexSet::clear()
{
  // reset so that there are no indices in the set any more; however,
  // as documented, the index set retains its size
  ranges.clear();
  is_compressed = true;
  largest_range = numbers::invalid_unsigned_int;
}



inline void
IndexSet::set_size(const size_type sz)
{
  Assert(ranges.empty(),
         ExcMessage("This function can only be called if the current "
                    "object does not yet contain any elements."));
  index_space_size = sz;
  is_compressed    = true;
}



inline IndexSet::size_type
IndexSet::size() const
{
  return index_space_size;
}



inline void
IndexSet::compress() const
{
  if (is_compressed == true)
    return;

  do_compress();
}



inline void
IndexSet::add_index(const size_type index)
{
  AssertIndexRange(index, index_space_size);

  const Range new_range(index, index + 1);
  if (ranges.size() == 0 || index > ranges.back().end)
    ranges.push_back(new_range);
  else if (index == ranges.back().end)
    ranges.back().end++;
  else
    ranges.insert(Utilities::lower_bound(ranges.begin(),
                                         ranges.end(),
                                         new_range),
                  new_range);
  is_compressed = false;
}



inline void
IndexSet::add_range(const size_type begin, const size_type end)
{
  Assert((begin < index_space_size) ||
           ((begin == index_space_size) && (end == index_space_size)),
         ExcIndexRangeType<size_type>(begin, 0, index_space_size));
  Assert(end <= index_space_size,
         ExcIndexRangeType<size_type>(end, 0, index_space_size + 1));
  AssertIndexRange(begin, end + 1);

  if (begin != end)
    {
      const Range new_range(begin, end);

      // the new index might be larger than the last index present in the
      // ranges. Then we can skip the binary search
      if (ranges.size() == 0 || begin > ranges.back().end)
        ranges.push_back(new_range);
      else
        ranges.insert(Utilities::lower_bound(ranges.begin(),
                                             ranges.end(),
                                             new_range),
                      new_range);
      is_compressed = false;
    }
}



template <typename ForwardIterator>
inline void
IndexSet::add_indices(const ForwardIterator &begin, const ForwardIterator &end)
{
  if (begin == end)
    return;

  // identify ranges in the given iterator range by checking whether some
  // indices happen to be consecutive. to avoid quadratic complexity when
  // calling add_range many times (as add_range() going into the middle of an
  // already existing range must shift entries around), we first collect a
  // vector of ranges.
  std::vector<std::pair<size_type, size_type>> tmp_ranges;
  bool                                         ranges_are_sorted = true;
  for (ForwardIterator p = begin; p != end;)
    {
      const size_type begin_index = *p;
      size_type       end_index   = begin_index + 1;
      ForwardIterator q           = p;
      ++q;
      while ((q != end) && (*q == end_index))
        {
          ++end_index;
          ++q;
        }

      tmp_ranges.emplace_back(begin_index, end_index);
      p = q;

      // if the starting index of the next go-around of the for loop is less
      // than the end index of the one just identified, then we will have at
      // least one pair of ranges that are not sorted, and consequently the
      // whole collection of ranges is not sorted.
      if (p != end && *p < end_index)
        ranges_are_sorted = false;
    }

  if (!ranges_are_sorted)
    std::sort(tmp_ranges.begin(), tmp_ranges.end());

  // if we have many ranges, we first construct a temporary index set (where
  // we add ranges in a consecutive way, so fast), otherwise, we work with
  // add_range(). the number 9 is chosen heuristically given the fact that
  // there are typically up to 8 independent ranges when adding the degrees of
  // freedom on a 3D cell or 9 when adding degrees of freedom of faces. if
  // doing cell-by-cell additions, we want to avoid repeated calls to
  // IndexSet::compress() which gets called upon merging two index sets, so we
  // want to be in the other branch then.
  if (tmp_ranges.size() > 9)
    {
      IndexSet tmp_set(size());
      tmp_set.ranges.reserve(tmp_ranges.size());
      for (const auto &i : tmp_ranges)
        tmp_set.add_range(i.first, i.second);
      this->add_indices(tmp_set);
    }
  else
    for (const auto &i : tmp_ranges)
      add_range(i.first, i.second);
}



inline bool
IndexSet::is_element(const size_type index) const
{
  if (ranges.empty() == false)
    {
      compress();

      // fast check whether the index is in the largest range
      Assert(largest_range < ranges.size(), ExcInternalError());
      if (index >= ranges[largest_range].begin &&
          index < ranges[largest_range].end)
        return true;

      // get the element after which we would have to insert a range that
      // consists of all elements from this element to the end of the index
      // range plus one. after this call we know that if p!=end() then
      // p->begin<=index unless there is no such range at all
      //
      // if the searched for element is an element of this range, then we're
      // done. otherwise, the element can't be in one of the following ranges
      // because otherwise p would be a different iterator
      //
      // since we already know the position relative to the largest range (we
      // called compress!), we can perform the binary search on ranges with
      // lower/higher number compared to the largest range
      std::vector<Range>::const_iterator p = std::upper_bound(
        ranges.begin() +
          (index < ranges[largest_range].begin ? 0 : largest_range + 1),
        index < ranges[largest_range].begin ? ranges.begin() + largest_range :
                                              ranges.end(),
        Range(index, size() + 1));

      if (p == ranges.begin())
        return ((index >= p->begin) && (index < p->end));

      Assert((p == ranges.end()) || (p->begin > index), ExcInternalError());

      // now move to that previous range
      --p;
      Assert(p->begin <= index, ExcInternalError());

      return (p->end > index);
    }

  // didn't find this index, so it's not in the set
  return false;
}



inline bool
IndexSet::is_contiguous() const
{
  compress();
  return (ranges.size() <= 1);
}



inline bool
IndexSet::is_empty() const
{
  return ranges.empty();
}



inline IndexSet::size_type
IndexSet::n_elements() const
{
  // make sure we have non-overlapping ranges
  compress();

  size_type v = 0;
  if (!ranges.empty())
    {
      Range &r = ranges.back();
      v        = r.nth_index_in_set + r.end - r.begin;
    }

#ifdef DEBUG
  size_type s = 0;
  for (const auto &range : ranges)
    s += (range.end - range.begin);
  Assert(s == v, ExcInternalError());
#endif

  return v;
}



inline unsigned int
IndexSet::n_intervals() const
{
  compress();
  return ranges.size();
}



inline IndexSet::size_type
IndexSet::largest_range_starting_index() const
{
  Assert(ranges.empty() == false, ExcMessage("IndexSet cannot be empty."));

  compress();
  const std::vector<Range>::const_iterator main_range =
    ranges.begin() + largest_range;

  return main_range->nth_index_in_set;
}



inline IndexSet::size_type
IndexSet::nth_index_in_set(const size_type n) const
{
  AssertIndexRange(n, n_elements());

  compress();

  // first check whether the index is in the largest range
  Assert(largest_range < ranges.size(), ExcInternalError());
  std::vector<Range>::const_iterator main_range =
    ranges.begin() + largest_range;
  if (n >= main_range->nth_index_in_set &&
      n < main_range->nth_index_in_set + (main_range->end - main_range->begin))
    return main_range->begin + (n - main_range->nth_index_in_set);

  // find out which chunk the local index n belongs to by using a binary
  // search. the comparator is based on the end of the ranges. Use the
  // position relative to main_range to subdivide the ranges
  Range r(n, n + 1);
  r.nth_index_in_set = n;
  std::vector<Range>::const_iterator range_begin, range_end;
  if (n < main_range->nth_index_in_set)
    {
      range_begin = ranges.begin();
      range_end   = main_range;
    }
  else
    {
      range_begin = main_range + 1;
      range_end   = ranges.end();
    }

  const std::vector<Range>::const_iterator p =
    Utilities::lower_bound(range_begin, range_end, r, Range::nth_index_compare);

  Assert(p != ranges.end(), ExcInternalError());
  return p->begin + (n - p->nth_index_in_set);
}



inline IndexSet::size_type
IndexSet::index_within_set(const size_type n) const
{
  // to make this call thread-safe, compress() must not be called through this
  // function
  Assert(is_compressed == true, ExcMessage("IndexSet must be compressed."));
  AssertIndexRange(n, size());

  // return immediately if the index set is empty
  if (is_empty())
    return numbers::invalid_dof_index;

  // check whether the index is in the largest range. use the result to
  // perform a one-sided binary search afterward
  Assert(largest_range < ranges.size(), ExcInternalError());
  std::vector<Range>::const_iterator main_range =
    ranges.begin() + largest_range;
  if (n >= main_range->begin && n < main_range->end)
    return (n - main_range->begin) + main_range->nth_index_in_set;

  Range                              r(n, n);
  std::vector<Range>::const_iterator range_begin, range_end;
  if (n < main_range->begin)
    {
      range_begin = ranges.begin();
      range_end   = main_range;
    }
  else
    {
      range_begin = main_range + 1;
      range_end   = ranges.end();
    }

  std::vector<Range>::const_iterator p =
    Utilities::lower_bound(range_begin, range_end, r, Range::end_compare);

  // if n is not in this set
  if (p == range_end || p->end == n || p->begin > n)
    return numbers::invalid_dof_index;

  Assert(p != ranges.end(), ExcInternalError());
  Assert(p->begin <= n, ExcInternalError());
  Assert(n < p->end, ExcInternalError());
  return (n - p->begin) + p->nth_index_in_set;
}



inline bool
IndexSet::operator==(const IndexSet &is) const
{
  Assert(size() == is.size(), ExcDimensionMismatch(size(), is.size()));

  compress();
  is.compress();

  return ranges == is.ranges;
}



inline bool
IndexSet::operator!=(const IndexSet &is) const
{
  Assert(size() == is.size(), ExcDimensionMismatch(size(), is.size()));

  compress();
  is.compress();

  return ranges != is.ranges;
}



template <typename Vector>
void
IndexSet::fill_binary_vector(Vector &vector) const
{
  Assert(vector.size() == size(), ExcDimensionMismatch(vector.size(), size()));

  compress();
  // first fill all elements of the vector with zeroes.
  std::fill(vector.begin(), vector.end(), 0);

  // then write ones into the elements whose indices are contained in the
  // index set
  for (const auto &range : ranges)
    for (size_type i = range.begin; i < range.end; ++i)
      vector[i] = 1;
}



template <class StreamType>
inline void
IndexSet::print(StreamType &out) const
{
  compress();
  out << "{";
  std::vector<Range>::const_iterator p;
  for (p = ranges.begin(); p != ranges.end(); ++p)
    {
      if (p->end - p->begin == 1)
        out << p->begin;
      else
        out << "[" << p->begin << "," << p->end - 1 << "]";

      if (p != --ranges.end())
        out << ", ";
    }
  out << "}" << std::endl;
}



template <class Archive>
inline void
IndexSet::Range::serialize(Archive &ar, const unsigned int)
{
  ar &begin &end &nth_index_in_set;
}



template <class Archive>
inline void
IndexSet::serialize(Archive &ar, const unsigned int)
{
  ar &ranges &is_compressed &index_space_size &largest_range;
}

DEAL_II_NAMESPACE_CLOSE

#endif


