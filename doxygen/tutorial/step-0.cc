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


#include <vector>
#include <algorithm>
#include <functional>
#include <iostream>

namespace test
{
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

int main()
{
  int ia[6] = {27, 210, 12, 47, 109, 83};




  std::vector<int, std::allocator<int>> vi(ia, ia + 6);
  std::cout << test::count_if(vi.begin(),
                              vi.end(),
                              std::not1(std::bind2nd(std::less<int>(), 40)))
            << std::endl;
  return 0;
}
