//include/deal.II-translator/numerics/history_0.txt
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

#ifndef dealii_storage_h
#define dealii_storage_h

#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>

#include <deque>
#include <memory>
#include <type_traits>

DEAL_II_NAMESPACE_OPEN

/**
 * 一个辅助类，用于存储一个有限大小的`T'类型的对象集合。如果元素的数量超过了指定的容器的最大尺寸，最古老的元素就会被移除。此外，元素的随机访问和移除也被实现。索引是相对于最后添加的元素进行的。
 * 为了优化容器，使其适用于需要内存的对象（即线性代数向量），删除一个元素并不释放内存。相反，该元素被保存在一个单独的缓存中，这样后续的添加就不需要重新分配内存了。
 * 该类的主要用途是在求解器中存储向量的历史。也就是说，如果在迭代
 * $k$ 时，我们存储了 $m$ 以前迭代 $\{k-1,k-2,...,k-m\}$
 * 的向量，那么新元素的加入将使对象包含迭代
 * $\{k,k-1,k-2,...,k-m+1\}$ 的元素。
 *
 *
 */
template <typename T>
class FiniteSizeHistory
{
public:
  static_assert(
    std::is_default_constructible<T>::value,
    "This class requires that the elements of type T are default constructible.");

  /**
   * 构造函数。      @param  max_elements
   * 储存在历史中的最大元素数。
   *
   */
  FiniteSizeHistory(const std::size_t max_elements = 0);

  /**
   * 通过复制添加新对象。
   * 如果达到最大的元素数，最老的元素会被删除。
   *
   */
  void
  add(const T &element);

  /**
   * 删除一个从最后添加的元素开始计算的索引 @p index,
   * 的元素。  `index==0`因此对应于删除新闻组元素。
   *
   */
  void
  remove(const std::size_t index);

  /**
   * 读/写访问一个索引为 @p index,
   * 的元素，从最后添加的元素算起。
   * `index==0`因此对应于newset元素。
   *
   */
  T &operator[](const std::size_t index);

  /**
   * 读取索引为 @p index,
   * 的元素，从最后添加的元素开始计算。
   * `index==0`因此对应于newset元素。
   *
   */
  const T &operator[](const std::size_t index) const;

  /**
   * 返回历史记录的当前大小。
   *
   */
  std::size_t
  size() const;

  /**
   * 返回历史记录的最大允许大小。
   *
   */
  std::size_t
  max_size() const;

  /**
   * 清除内容，包括缓存。
   *
   */
  void
  clear();

private:
  /**
   * 要存储的最大元素数。
   *
   */
  std::size_t max_n_elements;

  /**
   * 一个存储数据的deque。
   *
   */
  std::deque<std::unique_ptr<T>> data;

  /**
   * 一个用于缓存数据的deque，特别是在移除后再添加的情况下。
   *
   */
  std::deque<std::unique_ptr<T>> cache;
};



// -------------------  inline and template functions ----------------
#ifndef DOXYGEN



template <typename T>
FiniteSizeHistory<T>::FiniteSizeHistory(const std::size_t max_elements)
  : max_n_elements(max_elements)
{}



template <typename T>
void
FiniteSizeHistory<T>::remove(const std::size_t ind)
{
  AssertIndexRange(ind, data.size());
  auto el = std::move(data[ind]);
  data.erase(data.begin() + ind);

  cache.push_back(std::move(el));

  // whatever we do, we shall not store more than the maximum number of
  // elements
  Assert(data.size() + cache.size() <= max_n_elements, ExcInternalError());
}



template <typename T>
void
FiniteSizeHistory<T>::add(const T &element)
{
  std::unique_ptr<T> new_el;
  if (data.size() < max_n_elements)
    // have not reached the maximum number of elements yet
    {
      if (cache.size() == 0)
        // nothing is cached, just copy a given element
        {
          new_el = std::make_unique<T>(element);
        }
      else
        // something is cached, take one element and copy
        // the user provided one there.
        {
          new_el    = std::move(cache.back());
          (*new_el) = element;

          cache.pop_back(); // removes a pointer that is now a nullptr anyway
        }
    }
  else
    // we reached the maximum number of elements and
    // thus have to re-order/cycle elements currently stored
    {
      new_el    = std::move(data.back());
      (*new_el) = element;

      data.pop_back(); // removes a pointer that is now a nullptr anyway
    }

  // finally insert the new one where appropriate
  data.push_front(std::move(new_el));

  // whatever we do, we shall not store more than the maximum number of
  // elements
  Assert(data.size() + cache.size() <= max_n_elements, ExcInternalError());
}



template <typename T>
T &FiniteSizeHistory<T>::operator[](const std::size_t ind)
{
  AssertIndexRange(ind, data.size());
  return *data[ind];
}



template <typename T>
const T &FiniteSizeHistory<T>::operator[](const std::size_t ind) const
{
  AssertIndexRange(ind, data.size());
  return *data[ind];
}



template <typename T>
std::size_t
FiniteSizeHistory<T>::size() const
{
  return data.size();
}



template <typename T>
std::size_t
FiniteSizeHistory<T>::max_size() const
{
  return max_n_elements;
}



template <typename T>
void
FiniteSizeHistory<T>::clear()
{
  data.clear();
  cache.clear();
}

#endif // Doxygen

DEAL_II_NAMESPACE_CLOSE

#endif // dealii_storage_h


