//include/deal.II-translator/base/ndarray_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2021 by the deal.II authors
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

#ifndef dealii_ndarray_h
#define dealii_ndarray_h

#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>

#include <array>

DEAL_II_NAMESPACE_OPEN

#ifndef DOXYGEN
namespace internal
{
  namespace ndarray
  {
    // clang-format off
    /**
     * 一个变量模板辅助类，用于递归地 "展开
     * "ndarray的大小信息。这在一个例子中得到了最好的解释。
     * @code
     *  HelperArray<double, 1, 2, 3, 4>::type
     * == std::array<HelperArray<double, 2, 3, 4>::type, 1>
     * == std::array<std::array<HelperArray<double, 3, 4>::type, 2>, 1>
     * == std::array<std::array<std::array<HelperArray<double, 4>::type, 3>, 2>, 1>
     * == std::array<std::array<std::array<std::array<HelperArray<double>::type, 4>, 3>, 2>, 1>
     * == std::array<std::array<std::array<std::array<double, 4>, 3>, 2>, 1>
     * @endcode
     *
     *
     */
    template <typename T, std::size_t... Ns>
    struct HelperArray;
    // clang-format on

    /**
     * 递归地定义HelperArray<T, N, ...Ns>的类型别名
     * "type"，在HelperArray<T,  Ns...>::type 周围包裹一个 std::array
     * 。
     *
     */
    template <typename T, std::size_t N, std::size_t... Ns>
    struct HelperArray<T, N, Ns...>
    {
      using type = std::array<typename HelperArray<T, Ns...>::type, N>;
    };

    /**
     * 一旦没有 std::size_t
     * 模板参数，就结束递归，并简单地将类型别名设置为T类型。
     *
     */
    template <typename T>
    struct HelperArray<T>
    {
      using type = T;
    };
  } // namespace ndarray
} // namespace internal
#endif // DOXYGEN

/**
 * 用于方便地定义多维<a
 * href="https://en.cppreference.com/w/cpp/container/array">std::array</a>的（变量模板）类型别名
 * 我们试图用类型别名解决的问题如下。假设你想创建一个多维的双数数组，例如，等级为3，大小为2，3，4的第一，中间和最后的索引。那么使用C风格的数组，你可以简单地写道
 *
 * @code
 * double my_array[2][3][4] = { ... };
 * @endcode
 * 现在，有很多很好的理由可以说明为什么不鼓励使用C风格的数组（从与STL函数的不兼容到需要笨拙的包装器，以及在比较相等时的意外，等等），如果你想做同样的事情，使用更现代（和鼓励）的
 * `std::array` 类，那么你必须声明
 *
 * @code
 * std::array<std::array<std::array<double, 4>, 3>, 2> = { ... };
 * @endcode
 * `std::array`
 * 的重复看起来很别扭，更糟糕的是，索引范围已经颠倒了：最左边的索引范围是[0,2]，中间的索引范围是[0,3]，最右边的索引范围是[0,4)。我们通过提供一个ndarray类来解决这个问题，它允许你通过简单的书写来声明上述堆叠的
 * `std::array` 类型。
 *
 * @code
 * dealii::ndarray<double, 2, 3, 4> my_array = { ... };
 * @endcode
 *
 *
 *
 * @note   dealii::ndarray 只是以<a
 * href="https://en.cppreference.com/w/cpp/language/type_alias">type
 * alias</a>（"使用
 * "声明）为形式的语法糖。它不是一个deal.II特定的类，而只是一个帮助工具，用于干净地定义由
 * "堆叠" `std::array` 类实现的多维数组。
 *
 *
 */
template <typename T, std::size_t... Ns>
using ndarray = typename internal::ndarray::HelperArray<T, Ns...>::type;

DEAL_II_NAMESPACE_CLOSE

#endif


