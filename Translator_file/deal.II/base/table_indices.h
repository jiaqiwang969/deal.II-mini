//include/deal.II-translator/base/table_indices_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2021 by the deal.II authors
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

#ifndef dealii_table_indices_h
#define dealii_table_indices_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/template_constraints.h>

#include <algorithm>
#include <iterator>
#include <ostream>


DEAL_II_NAMESPACE_OPEN


/**
 * 一个代表固定大小的索引阵列的类。
 * 它被用于像TableBase和SymmetricTensor类这样的张量对象中，表示指数的嵌套选择。
 * @tparam  N 每个对象中存储的指数的数量。
 *
 *
 * @ingroup data
 *
 *
 */
template <int N>
class TableIndices
{
public:
  static_assert(N > 0,
                "TableIndices objects need to represent at least one index.");


  /**
   * 默认构造函数。这个构造函数将所有的指数设置为零。
   *
   */
  constexpr TableIndices() = default;

  /**
   * 构造函数。通过给定的参数 @p indices
   * 初始化此对象存储的指数，如果模板参数 @p N
   * 与参数的数量不同，这个构造函数将导致编译器错误。
   *
   */
  template <typename... T>
  constexpr TableIndices(const T... indices);

  /**
   * 只读访问<tt>i</tt>第1个索引的值。
   *
   */
  constexpr std::size_t operator[](const unsigned int i) const;

  /**
   * 写入访问<tt>i</tt>th索引的值。
   *
   */
  constexpr std::size_t &operator[](const unsigned int i);

  /**
   * 比较两个索引字段是否相等。
   *
   */
  constexpr bool
  operator==(const TableIndices<N> &other) const;

  /**
   * 比较两个索引字段的不等式。
   *
   */
  constexpr bool
  operator!=(const TableIndices<N> &other) const;

  /**
   * 按升序对索引进行排序。虽然这个操作对表对象不是很有用，但它被用于SymmetricTensor类。
   *
   */
  DEAL_II_CONSTEXPR void
  sort();

  /**
   * 使用[BOOST序列化库](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html)将此对象的数据写入或读出到一个流中，以便进行序列化。
   *
   */
  template <class Archive>
  void
  serialize(Archive &ar, const unsigned int version);

protected:
  /**
   * 在一个数组中存储索引。
   *
   */
  std::size_t indices[N]{};
};



 /* --------------------- Template and inline functions ---------------- */ 

template <int N>
template <typename... T>
constexpr TableIndices<N>::TableIndices(const T... args)
  : indices{static_cast<std::size_t>(args)...}
{
  static_assert(internal::TemplateConstraints::all_true<
                  std::is_integral<T>::value...>::value,
                "Not all of the parameters have integral type!");
  static_assert(sizeof...(T) == N, "Wrong number of constructor arguments!");
}


template <int N>
constexpr inline std::size_t TableIndices<N>::
                             operator[](const unsigned int i) const
{
  AssertIndexRange(i, N);
  return indices[i];
}


template <int N>
constexpr inline std::size_t &TableIndices<N>::operator[](const unsigned int i)
{
  AssertIndexRange(i, N);
  return indices[i];
}


template <int N>
constexpr bool
TableIndices<N>::operator==(const TableIndices<N> &other) const
{
  return std::equal(std::begin(indices),
                    std::end(indices),
                    std::begin(other.indices));
}


template <int N>
constexpr bool
TableIndices<N>::operator!=(const TableIndices<N> &other) const
{
  return !(*this == other);
}


template <int N>
DEAL_II_CONSTEXPR inline void
TableIndices<N>::sort()
{
  std::sort(std::begin(indices), std::end(indices));
}


template <int N>
template <class Archive>
inline void
TableIndices<N>::serialize(Archive &ar, const unsigned int)
{
  ar &indices;
}


/**
 * TableIndices对象的输出运算符；像这样以列表形式报告它们。
 * <code>[i1,i2,...]</code>  .
 * @relatesalso  TableIndices
 *
 *
 */
template <int N>
std::ostream &
operator<<(std::ostream &out, const TableIndices<N> &indices)
{
  out << '[';
  for (unsigned int i = 0; i < N; ++i)
    {
      out << indices[i];
      if (i + 1 != N)
        out << ',';
    }
  out << ']';

  return out;
}


DEAL_II_NAMESPACE_CLOSE

#endif


