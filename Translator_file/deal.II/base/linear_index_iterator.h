//include/deal.II-translator/base/linear_index_iterator_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2020 by the deal.II authors
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

#ifndef dealii_linear_index_iterator_h
#define dealii_linear_index_iterator_h

#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>


DEAL_II_NAMESPACE_OPEN
/**
 * deal.II中的许多类，如FullMatrix、TransposeTable和SparseMatrix，将它们的数据存储在连续的缓冲区中（当然，这些缓冲区的元素所代表的
 * <em> 解释 </em>
 * 可能很复杂）。例如，FullMatrix和TransposeTable分别以行大和列大的顺序存储数据，而对于SparseMatrix，从缓冲区位置到矩阵条目的映射
 * $\mathbf{A}(i, j)$
 * 更为复杂。然而，在任何情况下，元素的连续排列都可以实现随机访问迭代。
 * LinearIndexIterator提供了为这些类编写迭代器所需的大部分功能。LinearIndexIterator本质上是
 * <code>boost::iterator_facade</code> 的一个简化版本，它假定
 * <code>AccessorType</code>
 * 提供了某些成员（下文有记录），这些成员完全描述了迭代器的状态。这个类的预期用途是让容器定义自己的访问器类，然后使用奇怪的重复出现的模板模式（CRTP）技术来定义它们的迭代器。例如，这里是一个使用
 * LinearIndexIterator 来定义自己的迭代器类的容器。
 *
 *
 * @code
 * template <typename T>
 * class Container
 * {
 * protected:
 * // forward declaration for friendship
 * template <bool Constness>
 * class Iterator;
 *
 * template <bool Constness>
 * class Accessor
 * {
 * public:
 *   // const iterators store a const pointer
 *   using container_pointer_type
 *     = typename std::conditional<Constness,
 *                                 const Container<T>*,
 *                                 Container<T>*>::type;
 *
 *   // This alias is assumed to exist.
 *   using size_type = std::size_t;
 *
 *   // constructor.
 *   Accessor(const container_pointer_type container,
 *            const std::ptrdiff_t index);
 *
 *   // constructor.
 *   Accessor();
 *
 *   // get a constant reference to the current value.
 *   const T& value() const;
 *
 * protected:
 *   container_pointer_type container;
 *   std::ptrdiff_t linear_index;
 *
 *   // LinearIndexIterator needs access to linear_index and container.
 *   friend class LinearIndexIterator<Iterator<Constness>,
 *                                    Accessor<Constness>>;
 * };
 *
 * template <bool Constness>
 * class Iterator : public LinearIndexIterator<Iterator<Constness>,
 *                                             Accessor<Constness>>
 * {
 *   // Constructor.
 *   Iterator(Container<T> const container, const std::ptrdiff_t index);
 *
 *   // implement additional constructors here, but all state should be
 *   // contained in the Accessor, which is a member of the base class.
 * };
 *
 * public:
 * using size_type = std::size_t;
 * using const_iterator = Iterator<true>;
 * using iterator = Iterator<false>;
 *
 * iterator begin ();
 * iterator end ();
 *
 * const_iterator begin () const;
 * const_iterator end () const;
 * };
 * @endcode
 *
 * @tparam  DerivedIterator
 * 如上例所示，具体的迭代器类应该使用该类的CRTP技术：这为迭代器提供了模板式的比较和算术运算符。这是必要的，例如，
 * LinearIndexIterator::operator++() 要返回正确的类型。
 * @tparam  AccessorType LinearIndexIterator假定 <code>AccessorType</code>  模板参数有以下成员，这些成员完全描述了迭代器的当前状态。  <ol>   <li>  一个名为 <code>container</code> 的指针指向原始容器（例如，相关的稀疏矩阵）。对于 <code>const</code> 迭代器来说，这应该是一个 <code>const</code> 的指针。 </li>   <li>  一个名为 <code>linear_index</code> 的数组索引，用于存储容器存储缓冲区中的当前位置。  <code>linear_index</code>  不需要是一个整数：它可以是一个实现了  <code>operator+=</code>, <code>operator&lt;</code>  , 和  <code>operator==</code>  的类类型（可转换为容器的正确索引类型）。例如，可以通过实现  <code>operator+=</code>  和  <code>operator-</code>  的乘法因子来实现strided iterator。 </li>   </ol>  此外， <code>AccessorType</code>  应该声明相关的LinearIndexIterator实例化是一个 <code>friend</code> ，并定义一个 <code>size_type</code>  类型。
 *
 *
 * @note  TransposeTable使用这个模板来实现其迭代器。
 *
 *
 */
template <class DerivedIterator, class AccessorType>
class LinearIndexIterator
{
public:
  /**
   * 迭代器类别。
   *
   */
  using iterator_category = std::random_access_iterator_tag;

  /**
   * 当你解除对当前种类的迭代器的定义时，你得到的类型的别名。
   *
   */
  using value_type = AccessorType;

  /**
   * 差异类型。
   *
   */
  using difference_type = std::ptrdiff_t;

  /**
   * 引用类型。
   *
   */
  using reference = const value_type &;

  /**
   * 指针类型。
   *
   */
  using pointer = const value_type *;

  /**
   * 底层容器使用的尺寸类型。
   *
   */
  using size_type = typename value_type::size_type;

  /**
   * 复制操作符。
   *
   */
  DerivedIterator &
  operator=(const DerivedIterator &it);

  /**
   * 前缀增量。
   *
   */
  DerivedIterator &
  operator++();

  /**
   * 后缀增量。
   *
   */
  DerivedIterator
  operator++(int);

  /**
   * 前缀递减。
   *
   */
  DerivedIterator &
  operator--();

  /**
   * 后缀递减。
   *
   */
  DerivedIterator
  operator--(int);

  /**
   * 返回一个比当前迭代器提前 @p n 个条目的迭代器。
   *
   */
  DerivedIterator
  operator+(const difference_type n) const;

  /**
   * 返回一个比当前迭代器晚 @p n 个条目的迭代器。
   *
   */
  DerivedIterator
  operator-(const difference_type n) const;

  /**
   * 将迭代器的位置增加 @p n.  。
   *
   */
  DerivedIterator &
  operator+=(const difference_type n);

  /**
   * 递减迭代器的位置 @p n. 。
   *
   */
  DerivedIterator &
  operator-=(const difference_type n);

  /**
   * 返回当前迭代器与参数之间的距离。这个距离是通过对当前迭代器应用operator++()的次数来获得参数（对于正的返回值），或者operator--()（对于负的返回值）。
   *
   */
  difference_type
  operator-(const DerivedIterator &p) const;

  /**
   * 解除引用操作符。
   *
   */
  reference operator*() const;

  /**
   * 解除引用操作符。
   *
   */
  pointer operator->() const;

  /**
   * 比较运算符。如果两个迭代器都指向同一个容器中的同一个条目，则返回
   * <code>true</code> 。
   *
   */
  template <typename OtherIterator>
  friend typename std::enable_if<
    std::is_convertible<OtherIterator, DerivedIterator>::value,
    bool>::type
  operator==(const LinearIndexIterator &left, const OtherIterator &right)
  {
    const auto &right_2 = static_cast<const DerivedIterator &>(right);
    return left.accessor == right_2.accessor;
  }

  /**
   * 与operator==()相对应。
   *
   */
  template <typename OtherIterator>
  friend typename std::enable_if<
    std::is_convertible<OtherIterator, DerivedIterator>::value,
    bool>::type
  operator!=(const LinearIndexIterator &left, const OtherIterator &right)
  {
    return !(left == right);
  }

  /**
   * 比较运算符：使用与operator<()相同的排序，但也检查是否相等。
   * 这个函数只有在两个迭代器都指向同一个容器时才有效。
   *
   */
  bool
  operator<=(const DerivedIterator &) const;

  /**
   * 比较运算符：使用与operator>()相同的排序，但也检查是否相等。
   * 这个函数只有在两个迭代器都指向同一个容器时才有效。
   *
   */
  bool
  operator>=(const DerivedIterator &) const;

  /**
   * 比较运算符。如果第一个行号较小，或者行号相等且第一个索引较小，结果为真。
   * 这个函数只有在两个迭代器都指向同一个容器时才有效。
   *
   */
  bool
  operator<(const DerivedIterator &) const;

  /**
   * 比较运算符。与operator<()的工作方式相同，只是反过来了。
   *
   */
  bool
  operator>(const DerivedIterator &) const;

protected:
  /* 继承类应该有一个默认的构造函数。 
*
*/
  LinearIndexIterator() = default; // NOLINT

  /**
   * 复制一个访问器的构造函数。
   *
   */
  LinearIndexIterator(const AccessorType accessor);

protected:
  /**
   * 存储一个访问器类的对象。
   *
   */
  AccessorType accessor;
};



template <class DerivedIterator, class AccessorType>
inline DerivedIterator &
LinearIndexIterator<DerivedIterator, AccessorType>::
operator=(const DerivedIterator &it)
{
  accessor.container    = it.container;
  accessor.linear_index = it.linear_index;
  return static_cast<DerivedIterator &>(*this);
}



template <class DerivedIterator, class AccessorType>
inline DerivedIterator &
LinearIndexIterator<DerivedIterator, AccessorType>::operator++()
{
  return operator+=(1);
}



template <class DerivedIterator, class AccessorType>
inline DerivedIterator
LinearIndexIterator<DerivedIterator, AccessorType>::operator++(int)
{
  const DerivedIterator copy(this->accessor);
                        operator+=(1);
  return copy;
}



template <class DerivedIterator, class AccessorType>
inline DerivedIterator &
LinearIndexIterator<DerivedIterator, AccessorType>::operator--()
{
  return operator+=(-1);
}



template <class DerivedIterator, class AccessorType>
inline DerivedIterator
LinearIndexIterator<DerivedIterator, AccessorType>::operator--(int)
{
  const DerivedIterator copy(this->accessor);
                        operator+=(-1);
  return copy;
}



template <class DerivedIterator, class AccessorType>
inline DerivedIterator
LinearIndexIterator<DerivedIterator, AccessorType>::
operator+(const difference_type n) const
{
  DerivedIterator copy(this->accessor);
  copy += n;
  return copy;
}



template <class DerivedIterator, class AccessorType>
inline DerivedIterator
LinearIndexIterator<DerivedIterator, AccessorType>::
operator-(const difference_type n) const
{
  DerivedIterator copy(this->accessor);
  copy += -n;
  return copy;
}



template <class DerivedIterator, class AccessorType>
inline DerivedIterator &
LinearIndexIterator<DerivedIterator, AccessorType>::
operator+=(const difference_type n)
{
  accessor.linear_index += n;
  return static_cast<DerivedIterator &>(*this);
}



template <class DerivedIterator, class AccessorType>
inline DerivedIterator &
LinearIndexIterator<DerivedIterator, AccessorType>::
operator-=(const difference_type n)
{
  return operator+=(-n);
}



template <class DerivedIterator, class AccessorType>
inline
  typename LinearIndexIterator<DerivedIterator, AccessorType>::difference_type
  LinearIndexIterator<DerivedIterator, AccessorType>::
  operator-(const DerivedIterator &other) const
{
  Assert(this->accessor.container == other.accessor.container,
         ExcMessage(
           "Only iterators pointing to the same container can be compared."));
  return this->accessor.linear_index - other.accessor.linear_index;
}



template <class DerivedIterator, class AccessorType>
inline typename LinearIndexIterator<DerivedIterator, AccessorType>::reference
  LinearIndexIterator<DerivedIterator, AccessorType>::operator*() const
{
  return accessor;
}



template <class DerivedIterator, class AccessorType>
inline typename LinearIndexIterator<DerivedIterator, AccessorType>::pointer
  LinearIndexIterator<DerivedIterator, AccessorType>::operator->() const
{
  return &accessor;
}



template <class DerivedIterator, class AccessorType>
inline bool
LinearIndexIterator<DerivedIterator, AccessorType>::
operator<=(const DerivedIterator &other) const
{
  return (*this == other) || (*this < other);
}



template <class DerivedIterator, class AccessorType>
inline bool
LinearIndexIterator<DerivedIterator, AccessorType>::
operator>=(const DerivedIterator &other) const
{
  return !(*this < other);
}



template <class DerivedIterator, class AccessorType>
inline bool
LinearIndexIterator<DerivedIterator, AccessorType>::
operator<(const DerivedIterator &other) const
{
  Assert(this->accessor.container == other.accessor.container,
         ExcMessage(
           "Only iterators pointing to the same container can be compared."));
  return this->accessor.linear_index < other.accessor.linear_index;
}



template <class DerivedIterator, class AccessorType>
inline bool
LinearIndexIterator<DerivedIterator, AccessorType>::
operator>(const DerivedIterator &other) const
{
  return other < static_cast<const DerivedIterator &>(*this);
}



template <class DerivedIterator, class AccessorType>
inline LinearIndexIterator<DerivedIterator, AccessorType>::LinearIndexIterator(
  const AccessorType accessor)
  : accessor(accessor)
{}


DEAL_II_NAMESPACE_CLOSE

#endif


