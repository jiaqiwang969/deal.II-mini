//include/deal.II-translator/base/iterator_range_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2014 - 2020 by the deal.II authors
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

#ifndef dealii_iterator_range_h
#define dealii_iterator_range_h


#include <deal.II/base/config.h>

#include <deal.II/base/template_constraints.h>

#include <iterator>


DEAL_II_NAMESPACE_OPEN


// Forward declaration
#ifndef DOXYGEN
template <typename Iterator>
class IteratorOverIterators;
#endif


/**
 * 一个用来表示迭代器集合的类，可以用迭代器的范围来表示，其特征是有一个开始和一个结束的迭代器。正如C++中常见的那样，这些范围被指定为由一个开始迭代器和一个过去的结束迭代器定义的半开区间。
 * 这个类的目的是让诸如Triangulation和DoFHandler这样的类能够使用当前类型的对象从诸如 Triangulation::cell_iterators() 这样的函数中返回单元格迭代器的范围，然后这样的对象可以被用于C++11支持的基于范围的for循环，也见 @ref CPP11  "C++11标准"
 * 。
 * 例如，如果目标是在每个活动单元上设置用户标志，这样的循环可以是这样的。
 *
 * @code
 * Triangulation<dim> triangulation;
 * ...
 * for (auto &cell : triangulation.active_cell_iterators())
 *   cell->set_user_flag();
 * @endcode
 * 换句话说， <code>cell</code> 对象是迭代器，而由
 * Triangulation::active_cell_iterators()
 * 和类似函数返回的范围对象在概念上被认为是<i>collections
 * of iterators</i>。
 * 当然，这个类也可以用来表示其他迭代器范围，使用不同种类的迭代器进入其他容器。
 *
 *  <h3>Class design: Motivation</h3> 非正式地，C++11标准描述<a
 * href="http://en.wikipedia.org/wiki/C%2B%2B11#Range-based_for_loop">range-
 * based for loops</a>的方式如下。一个<i>range-based for
 * loop</i>的形式是
 *
 * @code
 * Container c;
 * for (auto v : c)
 *   statement;
 * @endcode
 * 其中 <code>c</code>
 * 是一个容器或集合，等同于以下的循环。
 *
 * @code
 * Container c;
 * for (auto tmp=c.begin(); tmp!=c.end(); ++tmp)
 *   {
 *     auto v =tmp;
 *     statement;
 *   }
 * @endcode
 * （精确的定义可以在这里找到：https://en.cppreference.com/w/cpp/language/range-for
 * 。换句话说，编译器引入了一个临时变量，<i>iterates</i>超过了容器或集合的元素，而出现在基于范围的for循环中的原始变量
 * <code>v</code> 代表了这些迭代器的<i>dereferenced</i>状态
 *
 * - 即集合的<i>elements</i>。
 * 在单元格上的循环中，我们通常希望保留循环变量是一个迭代器而不是一个值的事实。这是因为在deal.II中，我们从未实际使用过单元格迭代器的<i>dereferenced
 * state</i>：从概念上讲，它代表一个单元格，技术上它是由CellAccessor和DoFCellAccessor等类实现的，但这些类从未被明确使用。因此，我们所希望的是，像
 * Triangulation::active_cell_iterators()
 * 这样的调用返回一个对象，该对象代表一个<code>{begin,
 * begin+1, ..., end-1}</code>的<i>collection of
 * iterators</i>。这很方便地表示为半开区间
 * <code>[begin,end)</code>
 * 。然后，基于范围的for循环中的循环变量将依次接受这些迭代器中的每一个。
 *
 *  <h3>Class design: Implementation</h3>
 * 为了表示上述所需的语义，这个类存储了给定模板类型的半开区间的迭代器
 * <code>[b,e)</code>
 * 。其次，该类需要提供begin()和end()函数，如果你<i>dereference</i>
 * IteratorRange::begin(), 的结果，你会得到 <code>b</code>
 * 迭代器。此外，你必须能够递增由 IteratorRange::begin()
 * 返回的对象，以便 <code>*(++begin()) == b+1</code>
 * 。换句话说， IteratorRange::begin()
 * 必须返回一个迭代器，该迭代器在被解除引用时返回一个模板类型的迭代器
 * <code>Iterator</code>  :
 * 它是一个迭代器上的迭代器，其意义与你有一个指针进入一个指针数组相同。
 * 这是以 IteratorRange::IteratorOverIterators 类的形式实现的。
 *
 *
 * @ingroup CPP11
 *
 *
 */
template <typename Iterator>
class IteratorRange
{
public:
  /**
   * 迭代器类型的类型定义，该类型在其他迭代器上进行迭代。
   *
   */
  using IteratorOverIterators = dealii::IteratorOverIterators<Iterator>;


  /**
   * 本类所代表的迭代器类型的类型定义。
   *
   */
  using iterator = Iterator;

  /**
   * 默认构造函数。创建一个由两个默认构造的迭代器代表的范围。这个范围可能是空的（取决于迭代器的类型）。
   *
   */
  IteratorRange();

  /**
   * 构造函数。给出开始和结束迭代器，构造一个范围。
   * @param[in]  begin 指向该范围的第一个元素的迭代器
   * @param[in]  end
   * 指向该范围所代表的最后一个元素的迭代器。
   *
   */
  IteratorRange(const iterator begin, const iterator end);

  /**
   * 返回指向这个范围的第一个元素的迭代器。
   *
   */
  IteratorOverIterators
  begin();

  /**
   * 返回指向此范围的第一个元素的迭代器。
   *
   */
  IteratorOverIterators
  begin() const;

  /**
   * 返回指向此范围最后一个元素之后的元素的迭代器。
   *
   */
  IteratorOverIterators
  end() const;

  /**
   * 返回指向该范围最后一个元素之后的元素的迭代器。
   *
   */
  IteratorOverIterators
  end();

private:
  /**
   * 表征该范围的开始和结束的迭代器。
   *
   */
  const IteratorOverIterators it_begin;
  const IteratorOverIterators it_end;
};



/**
 * 一个实现迭代器上的迭代器的语义的类，在IteratorRange类的设计部分中讨论过。
 *
 *
 */
template <typename Iterator>
class IteratorOverIterators
{
public:
  /**
   * 对集合中的元素进行类型化定义，给它们一个更明显的名称。
   *
   */
  using BaseIterator = Iterator;

  /**
   * 构造函数。以这样的方式初始化这个迭代器-over-iterator，使其指向给定的参数。
   * @param  iterator 这个对象应该指向的一个迭代器。
   *
   */
  explicit IteratorOverIterators(const BaseIterator &iterator);

  /**
   * 去引用操作符。    @return  当前指向的集合中的迭代器。
   *
   */
  const BaseIterator &operator*() const;

  /**
   * 解除引用操作符。    @return
   * 当前指向的集合中的迭代器。
   *
   */
  const BaseIterator *operator->() const;

  /**
   * 前缀增量运算符。将当前的迭代器移动到集合的下一个元素，并返回新的值。
   *
   */
  IteratorOverIterators &
  operator++();

  /**
   * 后缀增量运算符。将当前的迭代器移动到集合的下一个元素，但返回迭代器的前一个值。
   *
   */
  IteratorOverIterators
  operator++(int);

  /**
   * 比较运算符  @param  i_o_i 另一个迭代器的迭代器。
   * @return
   * 返回当前迭代器是否指向一个与参数所代表的迭代器不同的对象。
   *
   */
  bool
  operator!=(const IteratorOverIterators &i_o_i) const;

  /**
   * 隐式转换操作符。      @warning
   * 当你调用这个转换操作符时（即，你将这个迭代器-over-iterators转换为我们当前指向的迭代器），你获得一个对这个底层迭代器的`const`引用。你真正能对这个结果做的唯一事情是取消引用本身：它可能指向一些有用的东西，但由于你不知道指向的对象在哪里，你不应该增加或减少你从这个操作符得到的迭代器。因此，返回的迭代器被标记为
   * "const"，因为这将防止你做任何其他事情，而不是取消引用它。
   *
   */
  operator const BaseIterator &() const;

  /**
   * 将该类标记为前向迭代器，并声明一些别名，这些别名是迭代器的标准，被算法用来查询它们所工作的迭代器的具体内容。
   *
   */
  using iterator_category = std::forward_iterator_tag;
  using value_type        = Iterator;
  using difference_type   = typename Iterator::difference_type;
  using pointer           = Iterator *;
  using reference         = Iterator &;

private:
  /**
   * 这个迭代器当前指向的对象。
   *
   */
  BaseIterator element_of_iterator_collection;
};



/**
 * 给出开始和结束的迭代器，创建一个IteratorRange类型的对象。
 *
 *
 */
template <typename BaseIterator>
IteratorRange<BaseIterator>
make_iterator_range(const BaseIterator &                         begin,
                    const typename identity<BaseIterator>::type &end)
{
  IteratorRange<BaseIterator> ir(begin, end);
  return ir;
}


// ------------------- template member functions


template <typename Iterator>
inline IteratorOverIterators<Iterator>::IteratorOverIterators(
  const BaseIterator &iterator)
  : element_of_iterator_collection(iterator)
{}



template <typename Iterator>
inline const typename IteratorOverIterators<Iterator>::BaseIterator &
  IteratorOverIterators<Iterator>::operator*() const
{
  return element_of_iterator_collection;
}



template <typename Iterator>
inline const typename IteratorOverIterators<Iterator>::BaseIterator *
  IteratorOverIterators<Iterator>::operator->() const
{
  return &element_of_iterator_collection;
}



template <typename Iterator>
inline IteratorOverIterators<Iterator> &
IteratorOverIterators<Iterator>::operator++()
{
  ++element_of_iterator_collection;
  return *this;
}



template <typename Iterator>
inline IteratorOverIterators<Iterator>
IteratorOverIterators<Iterator>::operator++(int)
{
  const IteratorOverIterators old_value = *this;
  ++element_of_iterator_collection;
  return *old_value;
}



template <typename Iterator>
inline bool
IteratorOverIterators<Iterator>::
operator!=(const IteratorOverIterators &i_o_i) const
{
  return element_of_iterator_collection != i_o_i.element_of_iterator_collection;
}



template <typename Iterator>
inline IteratorOverIterators<Iterator>::operator const BaseIterator &() const
{
  return element_of_iterator_collection;
}



template <typename Iterator>
inline IteratorRange<Iterator>::IteratorRange()
  : it_begin()
  , it_end()
{}



template <typename Iterator>
inline IteratorRange<Iterator>::IteratorRange(const iterator b,
                                              const iterator e)
  : it_begin(b)
  , it_end(e)
{}


template <typename Iterator>
inline typename IteratorRange<Iterator>::IteratorOverIterators
IteratorRange<Iterator>::begin()
{
  return it_begin;
}


template <typename Iterator>
inline typename IteratorRange<Iterator>::IteratorOverIterators
IteratorRange<Iterator>::begin() const
{
  return it_begin;
}


template <typename Iterator>
inline typename IteratorRange<Iterator>::IteratorOverIterators
IteratorRange<Iterator>::end()
{
  return it_end;
}


template <typename Iterator>
inline typename IteratorRange<Iterator>::IteratorOverIterators
IteratorRange<Iterator>::end() const
{
  return it_end;
}


DEAL_II_NAMESPACE_CLOSE

#endif


