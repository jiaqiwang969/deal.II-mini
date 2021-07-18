//include/deal.II-translator/base/synchronous_iterator_0.txt
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

#ifndef dealii_synchronous_iterator_h
#define dealii_synchronous_iterator_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>

#include <iterator>
#include <tuple>

DEAL_II_NAMESPACE_OPEN

/**
 * 一个表示一组迭代器的类，每个迭代器都同时递增一个。这通常用于像
 * <code>std::transform(a.begin(), a.end(), b.begin(), functor);</code>
 * 这样的调用，我们有同步的迭代器在容器中行进
 * <code>a,b</code>
 * 。如果这种类型的对象代表一个范围的结束，只有第一个元素被考虑（我们只有
 * <code>a.end()</code>  ，没有  <code>b.end()</code>  ）。在  step-35
 * 中给出了一个关于如何使用这个类的例子。
 * 当前类的模板参数应是 <code>std::tuple</code>
 * 类型，参数等于迭代器类型。 可以使用
 * <code>std::get<X>(synchronous_iterator.iterators)</code>
 * 访问各个迭代器，其中X是对应于所需迭代器的数字。
 * 这个类型，以及与之相关的辅助函数，被用作Threading
 * Building Blocks的blocked_range类型的Value概念。
 *
 *
 */
template <typename Iterators>
struct SynchronousIterators
{
  /**
   * 构造函数。
   *
   */
  SynchronousIterators(const Iterators &i);

  /**
   * 解除引用常量操作符。返回一个对当前类所代表的迭代器的常量引用。
   *
   */
  const Iterators &operator*() const;

  /**
   * 解除引用操作符。返回一个对当前类所代表的迭代器的引用。
   *
   */
  Iterators &operator*();

private:
  /**
   * 当前类所代表的迭代器的存储。
   *
   */
  Iterators iterators;
};



template <typename Iterators>
inline SynchronousIterators<Iterators>::SynchronousIterators(const Iterators &i)
  : iterators(i)
{}



template <typename Iterators>
inline const Iterators &SynchronousIterators<Iterators>::operator*() const
{
  return iterators;
}



template <typename Iterators>
inline Iterators &SynchronousIterators<Iterators>::operator*()
{
  return iterators;
}



/**
 * 返回第一个参数的第一个元素是否小于第二个参数的第一个元素。由于被比较的对象同时前进所有的元素，所以比较第一个元素就足够了。
 * @relatesalso  SynchronousIterators
 *
 *
 */
template <typename Iterators>
inline bool
operator<(const SynchronousIterators<Iterators> &a,
          const SynchronousIterators<Iterators> &b)
{
  return std::get<0>(*a) < std::get<0>(*b);
}



/**
 * 返回第一个和第二个参数之间的距离。由于被比较的对象在同一时间向前推进所有的元素，所以对第一个元素进行差分就足够了。
 * @relatesalso  SynchronousIterators
 *
 *
 */
template <typename Iterators>
inline std::size_t
operator-(const SynchronousIterators<Iterators> &a,
          const SynchronousIterators<Iterators> &b)
{
  Assert(std::distance(std::get<0>(*b), std::get<0>(*a)) >= 0,
         ExcInternalError());
  return std::distance(std::get<0>(*b), std::get<0>(*a));
}


/**
 * 通过 $n$ 推进一个迭代器的元组。
 * @relatesalso  SynchronousIterators
 *
 *
 */
template <typename I1, typename I2>
inline void
advance(std::tuple<I1, I2> &t, const unsigned int n)
{
  std::advance(std::get<0>(t), n);
  std::advance(std::get<1>(t), n);
}

/**
 * 以 $n$ 推进一个迭代器的元组。
 * @relatesalso  SynchronousIterators
 *
 *
 */
template <typename I1, typename I2, typename I3>
inline void
advance(std::tuple<I1, I2, I3> &t, const unsigned int n)
{
  std::advance(std::get<0>(t), n);
  std::advance(std::get<1>(t), n);
  std::advance(std::get<2>(t), n);
}

/**
 * 以 $n$ 推进一个迭代器的元组。
 * @relatesalso  SynchronousIterators
 *
 *
 */
template <typename I1, typename I2, typename I3, typename I4>
inline void
advance(std::tuple<I1, I2, I3, I4> &t, const unsigned int n)
{
  std::advance(std::get<0>(t), n);
  std::advance(std::get<1>(t), n);
  std::advance(std::get<2>(t), n);
  std::advance(std::get<3>(t), n);
}



/**
 * 将一个迭代器的元组提前1。
 * @relatesalso  SynchronousIterators
 *
 *
 */
template <typename I1, typename I2>
inline void
advance_by_one(std::tuple<I1, I2> &t)
{
  ++std::get<0>(t);
  ++std::get<1>(t);
}

/**
 * 将一个迭代器的元组前进1。
 * @relatesalso  SynchronousIterators
 *
 *
 */
template <typename I1, typename I2, typename I3>
inline void
advance_by_one(std::tuple<I1, I2, I3> &t)
{
  ++std::get<0>(t);
  ++std::get<1>(t);
  ++std::get<2>(t);
}

/**
 * 将一个迭代器的元组推进1。
 * @relatesalso  SynchronousIterators
 *
 *
 */
template <typename I1, typename I2, typename I3, typename I4>
inline void
advance_by_one(std::tuple<I1, I2, I3, I4> &t)
{
  ++std::get<0>(t);
  ++std::get<1>(t);
  ++std::get<2>(t);
  ++std::get<3>(t);
}



/**
 * 将这个迭代器的元素向前推进  $n$  。
 * @relatesalso  SynchronousIterators
 *
 *
 */
template <typename Iterators>
inline SynchronousIterators<Iterators>
operator+(const SynchronousIterators<Iterators> &a, const std::size_t n)
{
  SynchronousIterators<Iterators> x(a);
  dealii::advance(*x, n);
  return x;
}

/**
 * 将这个迭代器的元素前进1。
 * @relatesalso  SynchronousIterators
 *
 *
 */
template <typename Iterators>
inline SynchronousIterators<Iterators>
operator++(SynchronousIterators<Iterators> &a)
{
  dealii::advance_by_one(*a);
  return a;
}


/**
 * 比较同步迭代器的不平等。因为它们是同步进行的，所以只比较第一个元素就可以了。
 * @relatesalso  SynchronousIterators
 *
 *
 */
template <typename Iterators>
inline bool
operator!=(const SynchronousIterators<Iterators> &a,
           const SynchronousIterators<Iterators> &b)
{
  return (std::get<0>(*a) != std::get<0>(*b));
}

DEAL_II_NAMESPACE_CLOSE

#endif


