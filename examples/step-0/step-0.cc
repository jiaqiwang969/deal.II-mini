/* ---------------------------------------------------------------------
 *
 * Copyright (C) 1999 - 2019 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE.md at
 * the top level directory of deal.II.
 *
 * ---------------------------------------------------------------------

 */

// @sect3{Include files}


// 这是stl的vector容器的头文件。
#include <vector>
// 这是stl的algorithm头文件。
#include <algorithm>
// 这是stl的afunctional头文件。
#include <functional>
// 这是输入输出的头文件。
#include <iostream>

// count_if 核心算法，摘抄自STL标准模版库。这里用test
// namespace包裹，可以有效的防止相同名称的函数产生误调用的行为。
namespace test
{
  // 模版
  template <typename _InputIterator, typename _Predicate>
  typename std::iterator_traits<_InputIterator>::difference_type
  __count_if(_InputIterator __first, _InputIterator __last, _Predicate __pred)
  {
    typename std::iterator_traits<_InputIterator>::difference_type __n = 0;
    for (; __first != __last; ++__first)
      if (__pred(__first)) // 返回真假
        ++__n;
    return __n;
  }

  /**
   *  @brief Count the elements of a sequence for which a predicate is true.
   *  @ingroup non_mutating_algorithms
   *  @param  __first  An input iterator.
   *  @param  __last   An input iterator.
   *  @param  __pred   A predicate.
   *  @return   The number of iterators @c i in the range @p [__first,__last)
   *  for which @p __pred(*i) is true.
   */
  template <typename _InputIterator, typename _Predicate>
  inline typename std::iterator_traits<_InputIterator>::difference_type
  count_if(_InputIterator __first, _InputIterator __last, _Predicate __pred)
  {
    // concept requirements
    __glibcxx_function_requires(_InputIteratorConcept<_InputIterator>)
      __glibcxx_function_requires(
        _UnaryPredicateConcept<
          _Predicate,
          typename std::iterator_traits<_InputIterator>::value_type>)
        __glibcxx_requires_valid_range(__first, __last);

    return test::__count_if(__first,
                            __last,
                            __gnu_cxx::__ops::__pred_iter(__pred));
  }



} // namespace test

// 该函数的目的是实现计算出vector中大于40元素的个数。
int main()
{
  // `ia` 定义为基本的数组。
  int ia[6] = {27, 210, 12, 47, 109, 83};
  // std::vector是一个容器。众所周知,
  // 常用的数据结构不外乎array（数组）、list（链表）、 tree (树）、hash(表)
  // $\cdots$等等。根捨“数据在容器中的排列”特性,这些数据结构分为序列式(sequence)和关联式(associative)两种。

  // 所谓序列式容器, 其中的元素都可序(ordered)，但未必有序(sorted)。

  // 所谓关联式容器，观念上类似关联式数据库（实际上则简单许多）：每笔数据(每个元素）都有一个键值（key）和一个实值（value)。当元素被插人到关联式容器中时,容器内部结构（可能是
  // RB-tree, 也可能是 hash-tab $1 \mathrm{e}$)
  // 便依照其键值大小,以某种特定规则将这个元素放置于适当位置。关联式容器没有所谓头尾(只有最大元素和最小元素），所以不会有所谓push_back(),
  // push_front()，pop_back()、pop_front(), begin(), end() 这样的操作行为。

  // 这里，我们在vector容器里面放入 `int`
  // 元素类型。第二个vector的参数是分配器，用来分配内存，每一次分配 `int`
  // 大小的内存。通常情况下，自动匹配设置为默认值，因此可以省略。

  std::vector<int, std::allocator<int>> vi(ia, ia + 6);
  // 接下来，对vector进行一定的算法操作。
  // cout_if
  // 算法用来计算一定数据范围下，例如vector的头和尾，符合给定某一个算法条件的个数。
  // 接下来，设置初值，其中一种办法是找到 `ia`
  // vi.begin()和vi.end()输出的类型为迭代器，一种范化的指针。
  // 迭代器是一种行为类似指针的对象,
  // 而指针的各种行为中最常见也最重要的便是内容提领（dereference)
  // 和成员访问（member access ) , 因此, 迭代器最重要的编程工作就是对 operator*
  // 和 operator-> 进行重载（overloading ) 工作。
  // 可以说，迭代器是一种智能指针。
  // not1 和 bind2nd 为 adapter 适配器。 bind2nd 表示绑定第二参数为40。
  // 整个cout_if算法的第三个参数
  // 代表一个条件或动作，Predicate，它会传回真或者假。
  // less 为仿函数，或者叫函数对象。
  std::cout << test::count_if(vi.begin(),
                              vi.end(),
                              std::not1(std::bind2nd(std::less<int>(), 40)))
            << std::endl;
  return 0;
}
// 整个程序包含了配置器allocator、迭代器interator、仿函数functor、容器container、算法algorithm、配接器adpater六大组成成分。
