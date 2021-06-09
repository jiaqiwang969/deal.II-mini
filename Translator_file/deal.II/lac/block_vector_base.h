//include/deal.II-translator/lac/block_vector_base_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2021 by the deal.II authors
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

#ifndef dealii_block_vector_base_h
#define dealii_block_vector_base_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/numbers.h>
#include <deal.II/base/subscriptor.h>

#include <deal.II/lac/block_indices.h>
#include <deal.II/lac/exceptions.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/vector_operation.h>

#include <cmath>
#include <cstddef>
#include <iterator>
#include <type_traits>
#include <vector>

DEAL_II_NAMESPACE_OPEN


/*! @addtogroup Vectors
 *@{
 */

// Forward declaration
#ifndef DOXYGEN
template <typename>
class BlockVectorBase;
#endif

/**
 * 一个可以用来确定某个类型是否是块向量类型的类。比如说。
 *
 * @code
 * IsBlockVector<Vector<double> >::value
 * @endcode
 * 的值是false，而
 *
 * @code
 * IsBlockVector<BlockVector<double> >::value
 * @endcode
 * 是真。这在模板环境中有时很有用，我们可能想根据模板类型是表示普通的还是块状向量类型来做不同的事情。
 *
 *
 */
template <typename VectorType>
struct IsBlockVector
{
private:
  /**
   * 如果该类派生自BlockVectorBase，重载返回true，这就是块向量的作用。
   *
   */
  template <typename T>
  static std::true_type
  check_for_block_vector(const BlockVectorBase<T> *);

  /**
   * 对所有其他不是块状矩阵的潜在向量类型进行捕捉。
   *
   */
  static std::false_type
  check_for_block_vector(...);

public:
  /**
   * 一个可静态计算的值，表示该类的模板参数是否是块向量（事实上，该类型是否来自BlockVectorBase<T>）。
   *
   */
  static const bool value =
    std::is_same<decltype(check_for_block_vector(std::declval<VectorType *>())),
                 std::true_type>::value;
};


// instantiation of the static member
template <typename VectorType>
const bool IsBlockVector<VectorType>::value;



namespace internal
{
  /**
   * 实现块向量中迭代器的命名空间。
   *
   */
  namespace BlockVectorIterators
  {
    /**
     * 用于块向量的通用随机访问迭代器类。由于我们不希望有两个非常量迭代器和常量迭代器的类，我们采取第二个模板参数，表示我们指向的向量是否是常量对象。第一个模板参数总是使用中的块向量的数字类型。
     * 这个类满足了C++标准中定义的随机访问迭代器的所有要求。对这些迭代器的操作在块向量中的元素数量上是常数。然而，它们有时与向量中的块数成线性关系，但由于这在应用程序中很少动态变化，所以这是一个常数，我们再次得到该迭代器满足了随机访问迭代器的要求。
     *
     */
    template <class BlockVectorType, bool Constness>
    class Iterator
    {
    public:
      /**
       * 声明容器大小的类型。
       *
       */
      using size_type = types::global_dof_index;

      /**
       * 这个迭代器所指向的数字的类型。根据第二个模板参数的值，这要么是一个常数，要么是非常数。
       *
       */
      using value_type =
        typename std::conditional<Constness,
                                  const typename BlockVectorType::value_type,
                                  typename BlockVectorType::value_type>::type;

      /**
       * 声明一些别名，这些别名是迭代器的标准，被算法用来查询它们所工作的迭代器的具体内容。
       *
       */
      using iterator_category = std::random_access_iterator_tag;
      using difference_type   = std::ptrdiff_t;
      using reference         = typename BlockVectorType::reference;
      using pointer           = value_type *;

      using dereference_type = typename std::conditional<
        Constness,
        value_type,
        typename BlockVectorType::BlockType::reference>::type;

      /**
       * Typedef块向量的类型（根据第二个模板参数的不同，其constness也不同）。
       *
       */
      using BlockVector = typename std::
        conditional<Constness, const BlockVectorType, BlockVectorType>::type;

      /**
       * 从我们指向的向量和指向的元素的全局索引构造一个迭代器。
       * 根据这个类的<tt>Constness</tt>模板参数的值，这个构造函数的第一个参数要么是一个const，要么是非const的引用。
       *
       */
      Iterator(BlockVector &parent, const size_type global_index);

      /**
       * 从不同常量的迭代器复制构造函数。
       * @note
       * 从一个常数迭代器构造一个非常数迭代器是没有意义的。试图这样做将导致一个编译时错误（通过
       * <code>static_assert</code>  ）。
       *
       */
      Iterator(const Iterator<BlockVectorType, !Constness> &c);


      /**
       * 从具有相同常量的迭代器复制构造器。
       *
       */
      Iterator(const Iterator &c);

    private:
      /**
       * 本类内部使用的构造函数。参数与各自成员变量的值完全匹配。
       *
       */
      Iterator(BlockVector &   parent,
               const size_type global_index,
               const size_type current_block,
               const size_type index_within_block,
               const size_type next_break_forward,
               const size_type next_break_backward);

    public:
      /**
       * 复制操作符。
       *
       */
      Iterator &
      operator=(const Iterator &c);

      /**
       * 去引用操作符。如果模板参数<tt>Constness</tt>是<tt>true</tt>，那么就不可能向结果写入，使之成为一个const_iterator。
       *
       */
      dereference_type operator*() const;

      /**
       * 随机访问操作符，允许访问相对于当前指向的元素的任意元素。
       *
       */
      dereference_type operator[](const difference_type d) const;

      /**
       * 前缀递增操作符。这个操作符将迭代器推进到下一个元素，并返回一个对<tt>*this</tt>的引用。
       *
       */
      Iterator &
      operator++();

      /**
       * 后缀增量操作符。这个操作符将迭代器推进到下一个元素，并返回这个迭代器的旧值的副本。
       *
       */
      Iterator
      operator++(int);

      /**
       * 前缀递减运算符。这个操作符将迭代器缩回到上一个元素，并返回一个对<tt>*this</tt>的引用。
       *
       */
      Iterator &
      operator--();

      /**
       * 后缀递减运算符。该操作符将迭代器缩回到前一个元素，并返回该迭代器的旧值的副本。
       *
       */
      Iterator
      operator--(int);

      /**
       * 迭代器的平等比较。该操作符检查所指向的向量是否相同，如果不相同，则抛出一个异常。
       *
       */
      template <bool OtherConstness>
      bool
      operator==(const Iterator<BlockVectorType, OtherConstness> &i) const;

      /**
       * 比较迭代器的不平等。该操作符检查指向的向量是否相同，如果不相同，则抛出一个异常。
       *
       */
      template <bool OtherConstness>
      bool
      operator!=(const Iterator<BlockVectorType, OtherConstness> &i) const;

      /**
       * 检查这个迭代器是否指向给定参数所指向的元素之前的一个元素。这个操作符检查所指向的向量是否相同，如果不相同，则抛出一个异常。
       *
       */
      template <bool OtherConstness>
      bool
      operator<(const Iterator<BlockVectorType, OtherConstness> &i) const;

      /**
       * 与上面的比较运算符相同。
       *
       */
      template <bool OtherConstness>
      bool
      operator<=(const Iterator<BlockVectorType, OtherConstness> &i) const;

      /**
       * 与上面的比较运算符相同。
       *
       */
      template <bool OtherConstness>
      bool
      operator>(const Iterator<BlockVectorType, OtherConstness> &i) const;

      /**
       * 与上面的比较运算符相同。
       *
       */
      template <bool OtherConstness>
      bool
      operator>=(const Iterator<BlockVectorType, OtherConstness> &i) const;

      /**
       * 返回两个迭代器之间的距离，单位是元素。
       *
       */
      template <bool OtherConstness>
      difference_type
      operator-(const Iterator<BlockVectorType, OtherConstness> &i) const;

      /**
       * 返回一个迭代器，该迭代器在当前迭代器的前面有一定数量的元素。
       *
       */
      Iterator
      operator+(const difference_type &d) const;

      /**
       * 返回一个迭代器，该迭代器是在当前迭代器后面的给定数量的元素。
       *
       */
      Iterator
      operator-(const difference_type &d) const;

      /**
       * 一次性向前移动迭代器<tt>d</tt>元素，并返回结果。
       *
       */
      Iterator &
      operator+=(const difference_type &d);

      /**
       * 一次性向后移动迭代器<tt>d</tt>元素，并返回结果。
       *
       */
      Iterator &
      operator-=(const difference_type &d);

      /**
       * @addtogroup  Exceptions 
     * @{ 
       *
       */

      /**
       * 当人们对属于两个不同块向量的迭代器进行算术比较时产生的异常。
       *
       */
      DeclExceptionMsg(ExcPointerToDifferentVectors,
                       "Your program tried to compare iterators pointing to "
                       "different block vectors. There is no reasonable way "
                       "to do this.");

      //@}
    private:
      /**
       * 指向该迭代器所指向的块向量对象的指针。
       * 根据这个类的<tt>Constness</tt>模板参数的值，这是一个<tt>const</tt>或非<tt>const</tt>的指针。
       *
       */
      BlockVector *parent;

      /**
       * 我们目前所指向的元素的全局索引。
       *
       */
      size_type global_index;

      /**
       * 当前的块和当前指向的元素在该块中的索引。
       *
       */
      unsigned int current_block;
      size_type    index_within_block;

      /**
       * 全局元素地址的索引，当向前和向后移动时，我们必须转到另一个块。这些索引作为缓存保存，因为这比总是询问父对象要有效得多。
       *
       */
      size_type next_break_forward;
      size_type next_break_backward;

      /**
       * 向前移动一个元素。
       *
       */
      void
      move_forward();

      /**
       * 向后移动一个元素。
       *
       */
      void
      move_backward();

      // Mark all other instances of this template as friends.
      template <typename, bool>
      friend class Iterator;
    };
  } // namespace BlockVectorIterators
} // namespace internal


/**
 * 一个由几个块组成的向量，每个块代表它自己的一个向量。
 * BlockVector是一个给定类型的向量的集合（例如，deal.II
 * Vector对象， PETScWrappers::MPI::Vector
 * 对象，等等）。里面的每个向量都可以有不同的大小。
 * BlockVector的功能包括一个Vector所能做的一切，加上通过<tt>block(i)</tt>访问BlockVector内部的单个Vector。它也有一个完整的随机访问迭代器，就像其他Vector类或标准C++库模板
 * <tt>std::vector</tt>.
 * 一样。因此，所有在迭代器上工作的算法也能与这个类的对象工作。
 * 虽然这个基类通过将对其成员函数的调用分派给每个独立块上的相应函数来实现大部分功能，但这个类实际上并没有分配一些内存或改变向量的大小。为此，派生类的构造函数、赋值运算符和
 * reinit()
 * 函数负责。这个类只处理与块向量的实际向量类型无关的共同部分。
 *
 *  <h3>Accessing individual blocks, and resizing vectors</h3>
 * 除了将这个对象作为一个整体使用外，你还可以使用block()函数，将每个块作为一个向量单独使用。
 * 有一个注意事项：如果你改变了一个或几个块的大小，你必须调用块向量的函数
 * collect_sizes() 来更新其内部结构。
 * @attention  警告。如果你在没有调用 collect_sizes()
 * 的情况下改变单个块的大小，结果可能是不可预测的。由于性能的原因，调试版在此不检查一致性!
 * @see   @ref GlossBlockLA  "块（线性代数）"
 *
 *
 */
template <class VectorType>
class BlockVectorBase : public Subscriptor
{
public:
  /**
   * Typedef底层向量的类型。
   *
   */
  using BlockType = VectorType;

  /* 声明所有容器中使用的标准类型。这些类型与<tt>C++</tt>标准库 <tt>std::vector<...></tt> 类中的类型平行。这包括迭代器类型。 
*
*/
  using value_type    = typename BlockType::value_type;
  using pointer       = value_type *;
  using const_pointer = const value_type *;
  using iterator =
    dealii::internal::BlockVectorIterators::Iterator<BlockVectorBase, false>;
  using const_iterator =
    dealii::internal::BlockVectorIterators::Iterator<BlockVectorBase, true>;
  using reference       = typename BlockType::reference;
  using const_reference = typename BlockType::const_reference;
  using size_type       = types::global_dof_index;

  /**
   * 声明一个类型，该类型具有持有实值数字的功能，其精度与该类的模板参数相同。如果这个类的模板参数是一个实数数据类型，那么real_type就等于模板参数。
   * 如果模板参数是一个 std::complex
   * 类型，那么real_type等于复数的基础类型。
   * 这个别名被用来表示规范的返回类型。
   *
   */
  using real_type = typename BlockType::real_type;

  /**
   * 默认构造函数。
   *
   */
  BlockVectorBase() = default;

  /**
   * 复制构造函数。
   *
   */
  BlockVectorBase(const BlockVectorBase &  /*V*/ ) = default;

  /**
   * 移动构造函数。如果底层的 <code>VectorType</code>
   * 可以移动构造，那么参数向量的每个块都被移动到当前对象中，否则它们被复制。
   *
   */
  BlockVectorBase(BlockVectorBase &&  /*V*/ ) noexcept = default;

  /**
   * 在调整向量大小后更新内部结构。每当你重新引用一个块状向量的一个块时，内部数据结构就会被破坏。
   * 因此，你应该在所有块得到新的大小后再调用这个函数。
   *
   */
  void
  collect_sizes();

  /**
   * 在矩阵的所有子块上调用compress()函数。    只有在使用基于MPI的向量时才需要调用这个功能，为了兼容，在其他对象中也存在。    更多信息见 @ref GlossCompress  "压缩分布式对象"
   * 。
   *
   */
  void
  compress(::dealii::VectorOperation::values operation);

  /**
   * 对单个块的访问。
   *
   */
  BlockType &
  block(const unsigned int i);

  /**
   * 对单个区块的只读访问。
   *
   */
  const BlockType &
  block(const unsigned int i) const;

  /**
   * 返回一个描述块和全局索引之间映射的对象上的引用。这个函数的使用已被高度废弃，它应该在下一个版本中消失。
   *
   */
  const BlockIndices &
  get_block_indices() const;

  /**
   * 块的数量。
   *
   */
  unsigned int
  n_blocks() const;

  /**
   * 返回向量的尺寸。这是所有组件的尺寸之和。
   *
   */
  std::size_t
  size() const;

  /**
   * 返回向量的局部尺寸。这是所有组件的本地尺寸（即存储在当前处理器上的值）的总和。
   *
   */
  std::size_t
  locally_owned_size() const;

  /**
   * 返回一个索引集，描述这个向量的哪些元素是由当前处理器拥有的。请注意，这个索引集不包括这个向量可能作为鬼魂元素存储在本地，但实际上是由另一个处理器拥有的元素。因此，如果这是一个分布式向量，在不同处理器上返回的索引集将形成不相交的集合，加起来就是完整的索引集。很明显，如果一个向量只在一个处理器上创建，那么结果将满足
   * @code
   * vec.locally_owned_elements() == complete_index_set (vec.size())
   * @endcode
   * 对于块状向量，这个函数返回各个块的本地拥有的元素的并集，以它们各自的索引偏移量进行移位。
   *
   */
  IndexSet
  locally_owned_elements() const;

  /**
   * 返回一个指向第一个元素的迭代器。
   *
   */
  iterator
  begin();

  /**
   * 返回一个指向恒定块向量的第一个元素的迭代器。
   *
   */
  const_iterator
  begin() const;

  /**
   * 返回一个指向结束后的元素的迭代器。
   *
   */
  iterator
  end();

  /**
   * 返回一个迭代器，指向一个恒定块向量的结束后的元素。
   *
   */
  const_iterator
  end() const;

  /**
   * 访问组件，返回U(i)。
   *
   */
  value_type
  operator()(const size_type i) const;

  /**
   * 访问组件，返回U(i)作为一个可写的引用。
   *
   */
  reference
  operator()(const size_type i);

  /**
   * 访问组件，返回U(i)。    与operator()完全相同。
   *
   */
  value_type operator[](const size_type i) const;

  /**
   * 访问组件，返回U(i)作为一个可写的引用。
   * 与operator()完全相同。
   *
   */
  reference operator[](const size_type i);

  /**
   * 与通过operator()获取向量的单个元素不同，这个函数允许一次性获取整个元素集。要读取的元素的索引在第一个参数中说明，相应的值在第二个参数中返回。
   * 如果当前的向量被称为 @p v,
   * ，那么这个函数就等同于代码
   * @code
   * for (unsigned int i=0; i<indices.size(); ++i)
   *   values[i] = v[indices[i]];
   * @endcode
   * @pre  @p indices 和 @p values 数组的大小必须是一致的。
   *
   */
  template <typename OtherNumber>
  void
  extract_subvector_to(const std::vector<size_type> &indices,
                       std::vector<OtherNumber> &    values) const;

  /**
   * 这个函数不是通过operator()获得向量的单个元素，而是允许一次获得整个元素集。与前一个函数不同的是，这个函数通过取消引用前两个参数提供的迭代器范围内的所有元素来获得元素的索引，并将向量的值放入通过取消引用从第三个参数指向的位置开始的迭代器范围获得的内存位置。
   * 如果当前的向量被称为 @p v,
   * ，那么这个函数就等同于代码
   * @code
   * ForwardIterator indices_p = indices_begin;
   * OutputIterator  values_p  = values_begin;
   * while (indices_p != indices_end)
   * {
   *  values_p = v[*indices_p];
   *   ++indices_p;
   *   ++values_p;
   * }
   * @endcode
   * @pre  必须能够写进从 @p values_begin
   * 开始的内存位置，因为在 @p indices_begin 和 @p indices_end.
   * 之间有许多迭代器。
   *
   */
  template <typename ForwardIterator, typename OutputIterator>
  void
  extract_subvector_to(ForwardIterator       indices_begin,
                       const ForwardIterator indices_end,
                       OutputIterator        values_begin) const;

  /**
   * 复制操作：用给定的标量值填充向量的所有组件。
   *
   */
  BlockVectorBase &
  operator=(const value_type s);

  /**
   * 对相同类型的参数进行复制操作。
   *
   */
  BlockVectorBase &
  operator=(const BlockVectorBase &V);

  /**
   * 移动赋值运算符。如果`VectorType`是可移动的，将给定参数向量的每个块移动到当前对象中，否则复制它们。
   *
   */
  BlockVectorBase &
  operator=(BlockVectorBase &&  /*V*/ ) = default; // NOLINT

  /**
   * 对不同类型的模板参数进行复制操作。
   *
   */
  template <class VectorType2>
  BlockVectorBase &
  operator=(const BlockVectorBase<VectorType2> &V);

  /**
   * 从非块向量到块向量的复制操作。
   *
   */
  BlockVectorBase &
  operator=(const VectorType &v);

  /**
   * 检查两个块向量类型是否相等。只有当两个向量已经具有相同的块结构时，才允许该操作。
   *
   */
  template <class VectorType2>
  bool
  operator==(const BlockVectorBase<VectorType2> &v) const;

  /**
   * $U = U V$  : 标量乘积。
   *
   */
  value_type operator*(const BlockVectorBase &V) const;

  /**
   * 返回 $l_2$ -norm的平方。
   *
   */
  real_type
  norm_sqr() const;

  /**
   * 返回这个向量的元素的平均值。
   *
   */
  value_type
  mean_value() const;

  /**
   * 返回该向量的 $l_1$ -norm，即绝对值之和。
   *
   */
  real_type
  l1_norm() const;

  /**
   * 返回向量的 $l_2$ -Norm，即元素的平方根之和。
   *
   */
  real_type
  l2_norm() const;

  /**
   * 返回该向量元素的最大绝对值，也就是向量的 $l_\infty$
   * -norm。
   *
   */
  real_type
  linfty_norm() const;

  /**
   * 执行一个矢量加法和随后的内积的组合操作，返回内积的值。换句话说，这个函数的结果与用户调用的
   * @code
   * this->add(a, V);
   * return_value =this W;
   * @endcode
   * 这个函数存在的原因是这个操作比在deal.II的向量类（Vector<Number>和
   * LinearAlgebra::distributed::Vector<double>).
   * ）上分别调用这两个函数涉及的内存转移要少。这个方法只需要加载三个向量，
   * @p this,   @p V,   @p W,
   * 而调用单独的方法意味着要加载两次调用向量 @p this
   * 。由于大多数向量操作都有内存传输限制，这就使时间减少了25\%（如果
   * @p W 等于 @p this).
   * ，则减少50\%）对于复值向量，第二步中的标量乘法被实现为
   * $\left<v,w\right>=\sum_i v_i \bar{w_i}$  。
   *
   */
  value_type
  add_and_dot(const value_type       a,
              const BlockVectorBase &V,
              const BlockVectorBase &W);

  /**
   * 如果给定的全局索引在这个处理器的局部范围内，返回真。询问相应的块。
   *
   */
  bool
  in_local_range(const size_type global_index) const;

  /**
   * 返回向量是否只包含值为0的元素。这个函数主要用于内部一致性检查，在非调试模式下应该很少使用，因为它需要花费相当多的时间。
   *
   */
  bool
  all_zero() const;

  /**
   * 如果向量没有负的条目，即所有条目都是零或正的，则返回
   * @p true
   * 。例如，这个函数用于检查细化指标是否真的都是正的（或零）。
   *
   */
  bool
  is_non_negative() const;

  /**
   * 加法运算符。 快速等同于<tt>U.add(1, V)</tt>。
   *
   */
  BlockVectorBase &
  operator+=(const BlockVectorBase &V);

  /**
   * 减法运算符。 快速等同于<tt>U.add(-1, V)</tt>。
   *
   */
  BlockVectorBase &
  operator-=(const BlockVectorBase &V);


  /**
   * 一个集体添加操作。这个函数将存储在 @p values
   * 中的一整套数值添加到 @p indices. 指定的向量成分中。
   *
   */
  template <typename Number>
  void
  add(const std::vector<size_type> &indices, const std::vector<Number> &values);

  /**
   * 这是第二次集体添加操作。作为区别，这个函数需要一个deal.II的数值向量。
   *
   */
  template <typename Number>
  void
  add(const std::vector<size_type> &indices, const Vector<Number> &values);

  /**
   * 取一个<tt>n_elements</tt>连续存储的地址，并将其添加到向量中。处理上述其他两个<tt>add()</tt>函数未涵盖的所有情况。
   *
   */
  template <typename Number>
  void
  add(const size_type  n_elements,
      const size_type *indices,
      const Number *   values);

  /**
   * $U(0-DIM)+=s$  .
   * 在所有组件上增加<tt>s</tt>。注意，<tt>s</tt>是一个标量而不是一个矢量。
   *
   */
  void
  add(const value_type s);

  /**
   * U+=a*V。缩放向量的简单相加。
   *
   */
  void
  add(const value_type a, const BlockVectorBase &V);

  /**
   * U+=a*V+b*W。缩放向量的多次加法。
   *
   */
  void
  add(const value_type       a,
      const BlockVectorBase &V,
      const value_type       b,
      const BlockVectorBase &W);

  /**
   * U=s*U+V。缩放和简单的向量相加。
   *
   */
  void
  sadd(const value_type s, const BlockVectorBase &V);

  /**
   * U=s*U+a*V。缩放和简单的加法。
   *
   */
  void
  sadd(const value_type s, const value_type a, const BlockVectorBase &V);

  /**
   * U=s*U+a*V+b*W。缩放和多重加法。
   *
   */
  void
  sadd(const value_type       s,
       const value_type       a,
       const BlockVectorBase &V,
       const value_type       b,
       const BlockVectorBase &W);

  /**
   * U=s*U+a*V+b*W+c*X。缩放和多重加法。
   *
   */
  void
  sadd(const value_type       s,
       const value_type       a,
       const BlockVectorBase &V,
       const value_type       b,
       const BlockVectorBase &W,
       const value_type       c,
       const BlockVectorBase &X);

  /**
   * 将向量的每个元素按一个常数进行缩放。
   *
   */
  BlockVectorBase &
  operator*=(const value_type factor);

  /**
   * 用给定值的倒数来缩放向量的每个元素。
   *
   */
  BlockVectorBase &
  operator/=(const value_type factor);

  /**
   * 将该向量的每个元素乘以<tt>v</tt>的相应元素。
   *
   */
  template <class BlockVector2>
  void
  scale(const BlockVector2 &v);

  /**
   * U=a*V。赋值。
   *
   */
  template <class BlockVector2>
  void
  equ(const value_type a, const BlockVector2 &V);

  /**
   * 通过调用 <code>update_ghost_values</code>
   * 更新每个区块的鬼魂值。
   *
   */
  void
  update_ghost_values() const;

  /**
   * 确定这个对象的内存消耗（以字节为单位）的估计值。
   *
   */
  std::size_t
  memory_consumption() const;

protected:
  /**
   * 指向组件阵列的指针。
   *
   */
  std::vector<VectorType> components;

  /**
   * 管理全局索引和不同块内索引之间转换的对象。
   *
   */
  BlockIndices block_indices;

  // Make the iterator class a friend.
  template <typename N, bool C>
  friend class dealii::internal::BlockVectorIterators::Iterator;

  template <typename>
  friend class BlockVectorBase;
};


 /*@}*/ 

 /*----------------------- Inline functions ----------------------------------*/ 


#ifndef DOXYGEN
namespace internal
{
  namespace BlockVectorIterators
  {
    template <class BlockVectorType, bool Constness>
    inline Iterator<BlockVectorType, Constness>::Iterator(
      const Iterator<BlockVectorType, Constness> &c)
      : parent(c.parent)
      , global_index(c.global_index)
      , current_block(c.current_block)
      , index_within_block(c.index_within_block)
      , next_break_forward(c.next_break_forward)
      , next_break_backward(c.next_break_backward)
    {}



    template <class BlockVectorType, bool Constness>
    inline Iterator<BlockVectorType, Constness>::Iterator(
      const Iterator<BlockVectorType, !Constness> &c)
      : parent(c.parent)
      , global_index(c.global_index)
      , current_block(c.current_block)
      , index_within_block(c.index_within_block)
      , next_break_forward(c.next_break_forward)
      , next_break_backward(c.next_break_backward)
    {
      // Only permit copy-constructing const iterators from non-const
      // iterators, and not vice versa (i.e., Constness must always be
      // true).
      static_assert(Constness == true,
                    "Constructing a non-const iterator from a const iterator "
                    "does not make sense.");
    }



    template <class BlockVectorType, bool Constness>
    inline Iterator<BlockVectorType, Constness>::Iterator(
      BlockVector &   parent,
      const size_type global_index,
      const size_type current_block,
      const size_type index_within_block,
      const size_type next_break_forward,
      const size_type next_break_backward)
      : parent(&parent)
      , global_index(global_index)
      , current_block(current_block)
      , index_within_block(index_within_block)
      , next_break_forward(next_break_forward)
      , next_break_backward(next_break_backward)
    {}



    template <class BlockVectorType, bool Constness>
    inline Iterator<BlockVectorType, Constness> &
    Iterator<BlockVectorType, Constness>::operator=(const Iterator &c)
    {
      parent              = c.parent;
      global_index        = c.global_index;
      index_within_block  = c.index_within_block;
      current_block       = c.current_block;
      next_break_forward  = c.next_break_forward;
      next_break_backward = c.next_break_backward;

      return *this;
    }



    template <class BlockVectorType, bool Constness>
    inline typename Iterator<BlockVectorType, Constness>::dereference_type
      Iterator<BlockVectorType, Constness>::operator*() const
    {
      return parent->block(current_block)(index_within_block);
    }



    template <class BlockVectorType, bool Constness>
    inline typename Iterator<BlockVectorType, Constness>::dereference_type
      Iterator<BlockVectorType, Constness>::
      operator[](const difference_type d) const
    {
      // if the index pointed to is
      // still within the block we
      // currently point into, then we
      // can save the computation of
      // the block
      if ((global_index + d >= next_break_backward) &&
          (global_index + d <= next_break_forward))
        return parent->block(current_block)(index_within_block + d);

      // if the index is not within the
      // block of the block vector into
      // which we presently point, then
      // there is no way: we have to
      // search for the block. this can
      // be done through the parent
      // class as well.
      return (*parent)(global_index + d);
    }



    template <class BlockVectorType, bool Constness>
    inline Iterator<BlockVectorType, Constness> &
    Iterator<BlockVectorType, Constness>::operator++()
    {
      move_forward();
      return *this;
    }



    template <class BlockVectorType, bool Constness>
    inline Iterator<BlockVectorType, Constness>
    Iterator<BlockVectorType, Constness>::operator++(int)
    {
      const Iterator old_value = *this;
      move_forward();
      return old_value;
    }



    template <class BlockVectorType, bool Constness>
    inline Iterator<BlockVectorType, Constness> &
    Iterator<BlockVectorType, Constness>::operator--()
    {
      move_backward();
      return *this;
    }



    template <class BlockVectorType, bool Constness>
    inline Iterator<BlockVectorType, Constness>
    Iterator<BlockVectorType, Constness>::operator--(int)
    {
      const Iterator old_value = *this;
      move_backward();
      return old_value;
    }



    template <class BlockVectorType, bool Constness>
    template <bool OtherConstness>
    inline bool
    Iterator<BlockVectorType, Constness>::
    operator==(const Iterator<BlockVectorType, OtherConstness> &i) const
    {
      Assert(parent == i.parent, ExcPointerToDifferentVectors());

      return (global_index == i.global_index);
    }



    template <class BlockVectorType, bool Constness>
    template <bool OtherConstness>
    inline bool
    Iterator<BlockVectorType, Constness>::
    operator!=(const Iterator<BlockVectorType, OtherConstness> &i) const
    {
      Assert(parent == i.parent, ExcPointerToDifferentVectors());

      return (global_index != i.global_index);
    }



    template <class BlockVectorType, bool Constness>
    template <bool OtherConstness>
    inline bool
    Iterator<BlockVectorType, Constness>::
    operator<(const Iterator<BlockVectorType, OtherConstness> &i) const
    {
      Assert(parent == i.parent, ExcPointerToDifferentVectors());

      return (global_index < i.global_index);
    }



    template <class BlockVectorType, bool Constness>
    template <bool OtherConstness>
    inline bool
    Iterator<BlockVectorType, Constness>::
    operator<=(const Iterator<BlockVectorType, OtherConstness> &i) const
    {
      Assert(parent == i.parent, ExcPointerToDifferentVectors());

      return (global_index <= i.global_index);
    }



    template <class BlockVectorType, bool Constness>
    template <bool OtherConstness>
    inline bool
    Iterator<BlockVectorType, Constness>::
    operator>(const Iterator<BlockVectorType, OtherConstness> &i) const
    {
      Assert(parent == i.parent, ExcPointerToDifferentVectors());

      return (global_index > i.global_index);
    }



    template <class BlockVectorType, bool Constness>
    template <bool OtherConstness>
    inline bool
    Iterator<BlockVectorType, Constness>::
    operator>=(const Iterator<BlockVectorType, OtherConstness> &i) const
    {
      Assert(parent == i.parent, ExcPointerToDifferentVectors());

      return (global_index >= i.global_index);
    }



    template <class BlockVectorType, bool Constness>
    template <bool OtherConstness>
    inline typename Iterator<BlockVectorType, Constness>::difference_type
    Iterator<BlockVectorType, Constness>::
    operator-(const Iterator<BlockVectorType, OtherConstness> &i) const
    {
      Assert(parent == i.parent, ExcPointerToDifferentVectors());

      return (static_cast<signed int>(global_index) -
              static_cast<signed int>(i.global_index));
    }



    template <class BlockVectorType, bool Constness>
    inline Iterator<BlockVectorType, Constness>
    Iterator<BlockVectorType, Constness>::
    operator+(const difference_type &d) const
    {
      // if the index pointed to is
      // still within the block we
      // currently point into, then we
      // can save the computation of
      // the block
      if ((global_index + d >= next_break_backward) &&
          (global_index + d <= next_break_forward))
        return Iterator(*parent,
                        global_index + d,
                        current_block,
                        index_within_block + d,
                        next_break_forward,
                        next_break_backward);
      else
        // outside present block, so
        // have to seek new block
        // anyway
        return Iterator(*parent, global_index + d);
    }



    template <class BlockVectorType, bool Constness>
    inline Iterator<BlockVectorType, Constness>
    Iterator<BlockVectorType, Constness>::
    operator-(const difference_type &d) const
    {
      // if the index pointed to is
      // still within the block we
      // currently point into, then we
      // can save the computation of
      // the block
      if ((global_index - d >= next_break_backward) &&
          (global_index - d <= next_break_forward))
        return Iterator(*parent,
                        global_index - d,
                        current_block,
                        index_within_block - d,
                        next_break_forward,
                        next_break_backward);
      else
        // outside present block, so
        // have to seek new block
        // anyway
        return Iterator(*parent, global_index - d);
    }



    template <class BlockVectorType, bool Constness>
    inline Iterator<BlockVectorType, Constness> &
    Iterator<BlockVectorType, Constness>::operator+=(const difference_type &d)
    {
      // if the index pointed to is
      // still within the block we
      // currently point into, then we
      // can save the computation of
      // the block
      if ((global_index + d >= next_break_backward) &&
          (global_index + d <= next_break_forward))
        {
          global_index += d;
          index_within_block += d;
        }
      else
        // outside present block, so
        // have to seek new block
        // anyway
        *this = Iterator(*parent, global_index + d);

      return *this;
    }



    template <class BlockVectorType, bool Constness>
    inline Iterator<BlockVectorType, Constness> &
    Iterator<BlockVectorType, Constness>::operator-=(const difference_type &d)
    {
      // if the index pointed to is
      // still within the block we
      // currently point into, then we
      // can save the computation of
      // the block
      if ((global_index - d >= next_break_backward) &&
          (global_index - d <= next_break_forward))
        {
          global_index -= d;
          index_within_block -= d;
        }
      else
        // outside present block, so
        // have to seek new block
        // anyway
        *this = Iterator(*parent, global_index - d);

      return *this;
    }


    template <class BlockVectorType, bool Constness>
    Iterator<BlockVectorType, Constness>::Iterator(BlockVector &   parent,
                                                   const size_type global_index)
      : parent(&parent)
      , global_index(global_index)
    {
      // find which block we are
      // in. for this, take into
      // account that it happens at
      // times that people want to
      // initialize iterators
      // past-the-end
      if (global_index < parent.size())
        {
          const std::pair<size_type, size_type> indices =
            parent.block_indices.global_to_local(global_index);
          current_block      = indices.first;
          index_within_block = indices.second;

          next_break_backward =
            parent.block_indices.local_to_global(current_block, 0);
          next_break_forward =
            (parent.block_indices.local_to_global(current_block, 0) +
             parent.block_indices.block_size(current_block) - 1);
        }
      else
        // past the end. only have one
        // value for this
        {
          this->global_index  = parent.size();
          current_block       = parent.n_blocks();
          index_within_block  = 0;
          next_break_backward = global_index;
          next_break_forward  = numbers::invalid_size_type;
        };
    }



    template <class BlockVectorType, bool Constness>
    void
    Iterator<BlockVectorType, Constness>::move_forward()
    {
      if (global_index != next_break_forward)
        ++index_within_block;
      else
        {
          // ok, we traverse a boundary
          // between blocks:
          index_within_block = 0;
          ++current_block;

          // break backwards is now old
          // break forward
          next_break_backward = next_break_forward + 1;

          // compute new break forward
          if (current_block < parent->block_indices.size())
            next_break_forward +=
              parent->block_indices.block_size(current_block);
          else
            // if we are beyond the end,
            // then move the next
            // boundary arbitrarily far
            // away
            next_break_forward = numbers::invalid_size_type;
        };

      ++global_index;
    }



    template <class BlockVectorType, bool Constness>
    void
    Iterator<BlockVectorType, Constness>::move_backward()
    {
      if (global_index != next_break_backward)
        --index_within_block;
      else if (current_block != 0)
        {
          // ok, we traverse a boundary
          // between blocks:
          --current_block;
          index_within_block =
            parent->block_indices.block_size(current_block) - 1;

          // break forwards is now old
          // break backward
          next_break_forward = next_break_backward - 1;

          // compute new break forward
          next_break_backward -=
            parent->block_indices.block_size(current_block);
        }
      else
        // current block was 0, we now
        // get into unspecified terrain
        {
          --current_block;
          index_within_block  = numbers::invalid_size_type;
          next_break_forward  = 0;
          next_break_backward = 0;
        };

      --global_index;
    }


  } // namespace BlockVectorIterators

} // namespace internal



template <class VectorType>
inline std::size_t
BlockVectorBase<VectorType>::size() const
{
  return block_indices.total_size();
}



template <class VectorType>
inline std::size_t
BlockVectorBase<VectorType>::locally_owned_size() const
{
  std::size_t local_size = 0;
  for (unsigned int b = 0; b < n_blocks(); ++b)
    local_size += block(b).locally_owned_size();
  return local_size;
}



template <class VectorType>
inline IndexSet
BlockVectorBase<VectorType>::locally_owned_elements() const
{
  IndexSet is(size());

  // copy index sets from blocks into the global one, shifted
  // by the appropriate amount for each block
  for (unsigned int b = 0; b < n_blocks(); ++b)
    {
      IndexSet x = block(b).locally_owned_elements();
      is.add_indices(x, block_indices.block_start(b));
    }

  is.compress();

  return is;
}



template <class VectorType>
inline unsigned int
BlockVectorBase<VectorType>::n_blocks() const
{
  return block_indices.size();
}


template <class VectorType>
inline typename BlockVectorBase<VectorType>::BlockType &
BlockVectorBase<VectorType>::block(const unsigned int i)
{
  AssertIndexRange(i, n_blocks());

  return components[i];
}



template <class VectorType>
inline const typename BlockVectorBase<VectorType>::BlockType &
BlockVectorBase<VectorType>::block(const unsigned int i) const
{
  AssertIndexRange(i, n_blocks());

  return components[i];
}



template <class VectorType>
inline const BlockIndices &
BlockVectorBase<VectorType>::get_block_indices() const
{
  return block_indices;
}


template <class VectorType>
inline void
BlockVectorBase<VectorType>::collect_sizes()
{
  std::vector<size_type> sizes(n_blocks());

  for (size_type i = 0; i < n_blocks(); ++i)
    sizes[i] = block(i).size();

  block_indices.reinit(sizes);
}



template <class VectorType>
inline void
BlockVectorBase<VectorType>::compress(
  ::dealii::VectorOperation::values operation)
{
  for (unsigned int i = 0; i < n_blocks(); ++i)
    block(i).compress(operation);
}



template <class VectorType>
inline typename BlockVectorBase<VectorType>::iterator
BlockVectorBase<VectorType>::begin()
{
  return iterator(*this, 0U);
}



template <class VectorType>
inline typename BlockVectorBase<VectorType>::const_iterator
BlockVectorBase<VectorType>::begin() const
{
  return const_iterator(*this, 0U);
}


template <class VectorType>
inline typename BlockVectorBase<VectorType>::iterator
BlockVectorBase<VectorType>::end()
{
  return iterator(*this, size());
}



template <class VectorType>
inline typename BlockVectorBase<VectorType>::const_iterator
BlockVectorBase<VectorType>::end() const
{
  return const_iterator(*this, size());
}


template <class VectorType>
inline bool
BlockVectorBase<VectorType>::in_local_range(const size_type global_index) const
{
  const std::pair<size_type, size_type> local_index =
    block_indices.global_to_local(global_index);

  return components[local_index.first].in_local_range(global_index);
}


template <class VectorType>
bool
BlockVectorBase<VectorType>::all_zero() const
{
  for (size_type i = 0; i < n_blocks(); ++i)
    if (components[i].all_zero() == false)
      return false;

  return true;
}



template <class VectorType>
bool
BlockVectorBase<VectorType>::is_non_negative() const
{
  for (size_type i = 0; i < n_blocks(); ++i)
    if (components[i].is_non_negative() == false)
      return false;

  return true;
}



template <class VectorType>
typename BlockVectorBase<VectorType>::value_type BlockVectorBase<VectorType>::
                                                 operator*(const BlockVectorBase<VectorType> &v) const
{
  Assert(n_blocks() == v.n_blocks(),
         ExcDimensionMismatch(n_blocks(), v.n_blocks()));

  value_type sum = 0.;
  for (size_type i = 0; i < n_blocks(); ++i)
    sum += components[i] * v.components[i];

  return sum;
}


template <class VectorType>
typename BlockVectorBase<VectorType>::real_type
BlockVectorBase<VectorType>::norm_sqr() const
{
  real_type sum = 0.;
  for (size_type i = 0; i < n_blocks(); ++i)
    sum += components[i].norm_sqr();

  return sum;
}



template <class VectorType>
typename BlockVectorBase<VectorType>::value_type
BlockVectorBase<VectorType>::mean_value() const
{
  value_type sum = 0.;
  // need to do static_cast as otherwise it won't work with
  // value_type=complex<T>
  for (size_type i = 0; i < n_blocks(); ++i)
    sum += components[i].mean_value() *
           (typename numbers::NumberTraits<value_type>::real_type(
             components[i].size()));

  return sum / (typename numbers::NumberTraits<value_type>::real_type(size()));
}



template <class VectorType>
typename BlockVectorBase<VectorType>::real_type
BlockVectorBase<VectorType>::l1_norm() const
{
  real_type sum = 0.;
  for (size_type i = 0; i < n_blocks(); ++i)
    sum += components[i].l1_norm();

  return sum;
}



template <class VectorType>
typename BlockVectorBase<VectorType>::real_type
BlockVectorBase<VectorType>::l2_norm() const
{
  return std::sqrt(norm_sqr());
}



template <class VectorType>
typename BlockVectorBase<VectorType>::real_type
BlockVectorBase<VectorType>::linfty_norm() const
{
  real_type sum = 0.;
  for (size_type i = 0; i < n_blocks(); ++i)
    {
      value_type newval = components[i].linfty_norm();
      if (sum < newval)
        sum = newval;
    }
  return sum;
}



template <class VectorType>
typename BlockVectorBase<VectorType>::value_type
BlockVectorBase<VectorType>::add_and_dot(
  const typename BlockVectorBase<VectorType>::value_type a,
  const BlockVectorBase<VectorType> &                    V,
  const BlockVectorBase<VectorType> &                    W)
{
  AssertDimension(n_blocks(), V.n_blocks());
  AssertDimension(n_blocks(), W.n_blocks());

  value_type sum = 0.;
  for (size_type i = 0; i < n_blocks(); ++i)
    sum += components[i].add_and_dot(a, V.components[i], W.components[i]);

  return sum;
}



template <class VectorType>
BlockVectorBase<VectorType> &
BlockVectorBase<VectorType>::operator+=(const BlockVectorBase<VectorType> &v)
{
  Assert(n_blocks() == v.n_blocks(),
         ExcDimensionMismatch(n_blocks(), v.n_blocks()));

  for (size_type i = 0; i < n_blocks(); ++i)
    {
      components[i] += v.components[i];
    }

  return *this;
}



template <class VectorType>
BlockVectorBase<VectorType> &
BlockVectorBase<VectorType>::operator-=(const BlockVectorBase<VectorType> &v)
{
  Assert(n_blocks() == v.n_blocks(),
         ExcDimensionMismatch(n_blocks(), v.n_blocks()));

  for (size_type i = 0; i < n_blocks(); ++i)
    {
      components[i] -= v.components[i];
    }
  return *this;
}



template <class VectorType>
template <typename Number>
inline void
BlockVectorBase<VectorType>::add(const std::vector<size_type> &indices,
                                 const std::vector<Number> &   values)
{
  Assert(indices.size() == values.size(),
         ExcDimensionMismatch(indices.size(), values.size()));
  add(indices.size(), indices.data(), values.data());
}



template <class VectorType>
template <typename Number>
inline void
BlockVectorBase<VectorType>::add(const std::vector<size_type> &indices,
                                 const Vector<Number> &        values)
{
  Assert(indices.size() == values.size(),
         ExcDimensionMismatch(indices.size(), values.size()));
  const size_type n_indices = indices.size();
  for (size_type i = 0; i < n_indices; ++i)
    (*this)(indices[i]) += values(i);
}



template <class VectorType>
template <typename Number>
inline void
BlockVectorBase<VectorType>::add(const size_type  n_indices,
                                 const size_type *indices,
                                 const Number *   values)
{
  for (size_type i = 0; i < n_indices; ++i)
    (*this)(indices[i]) += values[i];
}



template <class VectorType>
void
BlockVectorBase<VectorType>::add(const value_type a)
{
  AssertIsFinite(a);

  for (size_type i = 0; i < n_blocks(); ++i)
    {
      components[i].add(a);
    }
}



template <class VectorType>
void
BlockVectorBase<VectorType>::add(const value_type                   a,
                                 const BlockVectorBase<VectorType> &v)
{
  AssertIsFinite(a);

  Assert(n_blocks() == v.n_blocks(),
         ExcDimensionMismatch(n_blocks(), v.n_blocks()));

  for (size_type i = 0; i < n_blocks(); ++i)
    {
      components[i].add(a, v.components[i]);
    }
}



template <class VectorType>
void
BlockVectorBase<VectorType>::add(const value_type                   a,
                                 const BlockVectorBase<VectorType> &v,
                                 const value_type                   b,
                                 const BlockVectorBase<VectorType> &w)
{
  AssertIsFinite(a);
  AssertIsFinite(b);

  Assert(n_blocks() == v.n_blocks(),
         ExcDimensionMismatch(n_blocks(), v.n_blocks()));
  Assert(n_blocks() == w.n_blocks(),
         ExcDimensionMismatch(n_blocks(), w.n_blocks()));


  for (size_type i = 0; i < n_blocks(); ++i)
    {
      components[i].add(a, v.components[i], b, w.components[i]);
    }
}



template <class VectorType>
void
BlockVectorBase<VectorType>::sadd(const value_type                   x,
                                  const BlockVectorBase<VectorType> &v)
{
  AssertIsFinite(x);

  Assert(n_blocks() == v.n_blocks(),
         ExcDimensionMismatch(n_blocks(), v.n_blocks()));

  for (size_type i = 0; i < n_blocks(); ++i)
    {
      components[i].sadd(x, v.components[i]);
    }
}



template <class VectorType>
void
BlockVectorBase<VectorType>::sadd(const value_type                   x,
                                  const value_type                   a,
                                  const BlockVectorBase<VectorType> &v)
{
  AssertIsFinite(x);
  AssertIsFinite(a);

  Assert(n_blocks() == v.n_blocks(),
         ExcDimensionMismatch(n_blocks(), v.n_blocks()));

  for (size_type i = 0; i < n_blocks(); ++i)
    {
      components[i].sadd(x, a, v.components[i]);
    }
}



template <class VectorType>
void
BlockVectorBase<VectorType>::sadd(const value_type                   x,
                                  const value_type                   a,
                                  const BlockVectorBase<VectorType> &v,
                                  const value_type                   b,
                                  const BlockVectorBase<VectorType> &w)
{
  AssertIsFinite(x);
  AssertIsFinite(a);
  AssertIsFinite(b);

  Assert(n_blocks() == v.n_blocks(),
         ExcDimensionMismatch(n_blocks(), v.n_blocks()));
  Assert(n_blocks() == w.n_blocks(),
         ExcDimensionMismatch(n_blocks(), w.n_blocks()));

  for (size_type i = 0; i < n_blocks(); ++i)
    {
      components[i].sadd(x, a, v.components[i], b, w.components[i]);
    }
}



template <class VectorType>
void
BlockVectorBase<VectorType>::sadd(const value_type                   x,
                                  const value_type                   a,
                                  const BlockVectorBase<VectorType> &v,
                                  const value_type                   b,
                                  const BlockVectorBase<VectorType> &w,
                                  const value_type                   c,
                                  const BlockVectorBase<VectorType> &y)
{
  AssertIsFinite(x);
  AssertIsFinite(a);
  AssertIsFinite(b);
  AssertIsFinite(c);

  Assert(n_blocks() == v.n_blocks(),
         ExcDimensionMismatch(n_blocks(), v.n_blocks()));
  Assert(n_blocks() == w.n_blocks(),
         ExcDimensionMismatch(n_blocks(), w.n_blocks()));
  Assert(n_blocks() == y.n_blocks(),
         ExcDimensionMismatch(n_blocks(), y.n_blocks()));

  for (size_type i = 0; i < n_blocks(); ++i)
    {
      components[i].sadd(
        x, a, v.components[i], b, w.components[i], c, y.components[i]);
    }
}



template <class VectorType>
template <class BlockVector2>
void
BlockVectorBase<VectorType>::scale(const BlockVector2 &v)
{
  Assert(n_blocks() == v.n_blocks(),
         ExcDimensionMismatch(n_blocks(), v.n_blocks()));
  for (size_type i = 0; i < n_blocks(); ++i)
    components[i].scale(v.block(i));
}



template <class VectorType>
std::size_t
BlockVectorBase<VectorType>::memory_consumption() const
{
  return (MemoryConsumption::memory_consumption(this->block_indices) +
          MemoryConsumption::memory_consumption(this->components));
}



template <class VectorType>
template <class BlockVector2>
void
BlockVectorBase<VectorType>::equ(const value_type a, const BlockVector2 &v)
{
  AssertIsFinite(a);

  Assert(n_blocks() == v.n_blocks(),
         ExcDimensionMismatch(n_blocks(), v.n_blocks()));

  for (size_type i = 0; i < n_blocks(); ++i)
    components[i].equ(a, v.components[i]);
}



template <class VectorType>
void
BlockVectorBase<VectorType>::update_ghost_values() const
{
  for (size_type i = 0; i < n_blocks(); ++i)
    block(i).update_ghost_values();
}



template <class VectorType>
BlockVectorBase<VectorType> &
BlockVectorBase<VectorType>::operator=(const value_type s)
{
  AssertIsFinite(s);

  for (size_type i = 0; i < n_blocks(); ++i)
    components[i] = s;

  return *this;
}


template <class VectorType>
BlockVectorBase<VectorType> &
BlockVectorBase<VectorType>::operator=(const BlockVectorBase<VectorType> &v)
{
  AssertDimension(n_blocks(), v.n_blocks());

  for (size_type i = 0; i < n_blocks(); ++i)
    components[i] = v.components[i];

  return *this;
}


template <class VectorType>
template <class VectorType2>
BlockVectorBase<VectorType> &
BlockVectorBase<VectorType>::operator=(const BlockVectorBase<VectorType2> &v)
{
  AssertDimension(n_blocks(), v.n_blocks());

  for (size_type i = 0; i < n_blocks(); ++i)
    components[i] = v.components[i];

  return *this;
}



template <class VectorType>
BlockVectorBase<VectorType> &
BlockVectorBase<VectorType>::operator=(const VectorType &v)
{
  Assert(size() == v.size(), ExcDimensionMismatch(size(), v.size()));

  size_type index_v = 0;
  for (size_type b = 0; b < n_blocks(); ++b)
    for (size_type i = 0; i < block(b).size(); ++i, ++index_v)
      block(b)(i) = v(index_v);

  return *this;
}



template <class VectorType>
template <class VectorType2>
inline bool
BlockVectorBase<VectorType>::
operator==(const BlockVectorBase<VectorType2> &v) const
{
  Assert(block_indices == v.block_indices, ExcDifferentBlockIndices());

  for (size_type i = 0; i < n_blocks(); ++i)
    if (!(components[i] == v.components[i]))
      return false;

  return true;
}



template <class VectorType>
inline BlockVectorBase<VectorType> &
BlockVectorBase<VectorType>::operator*=(const value_type factor)
{
  AssertIsFinite(factor);

  for (size_type i = 0; i < n_blocks(); ++i)
    components[i] *= factor;

  return *this;
}



template <class VectorType>
inline BlockVectorBase<VectorType> &
BlockVectorBase<VectorType>::operator/=(const value_type factor)
{
  AssertIsFinite(factor);
  Assert(factor != 0., ExcDivideByZero());

  const value_type inverse_factor = value_type(1.) / factor;

  for (size_type i = 0; i < n_blocks(); ++i)
    components[i] *= inverse_factor;

  return *this;
}


template <class VectorType>
inline typename BlockVectorBase<VectorType>::value_type
BlockVectorBase<VectorType>::operator()(const size_type i) const
{
  const std::pair<unsigned int, size_type> local_index =
    block_indices.global_to_local(i);
  return components[local_index.first](local_index.second);
}



template <class VectorType>
inline typename BlockVectorBase<VectorType>::reference
BlockVectorBase<VectorType>::operator()(const size_type i)
{
  const std::pair<unsigned int, size_type> local_index =
    block_indices.global_to_local(i);
  return components[local_index.first](local_index.second);
}



template <class VectorType>
inline typename BlockVectorBase<VectorType>::value_type
  BlockVectorBase<VectorType>::operator[](const size_type i) const
{
  return operator()(i);
}



template <class VectorType>
inline typename BlockVectorBase<VectorType>::reference
  BlockVectorBase<VectorType>::operator[](const size_type i)
{
  return operator()(i);
}



template <typename VectorType>
template <typename OtherNumber>
inline void
BlockVectorBase<VectorType>::extract_subvector_to(
  const std::vector<size_type> &indices,
  std::vector<OtherNumber> &    values) const
{
  for (size_type i = 0; i < indices.size(); ++i)
    values[i] = operator()(indices[i]);
}



template <typename VectorType>
template <typename ForwardIterator, typename OutputIterator>
inline void
BlockVectorBase<VectorType>::extract_subvector_to(
  ForwardIterator       indices_begin,
  const ForwardIterator indices_end,
  OutputIterator        values_begin) const
{
  while (indices_begin != indices_end)
    {
      *values_begin = operator()(*indices_begin);
      indices_begin++;
      values_begin++;
    }
}

#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif


