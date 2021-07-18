//include/deal.II-translator/hp/collection_0.txt
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

#ifndef dealii_hp_collection_h
#define dealii_hp_collection_h

#include <deal.II/base/config.h>

#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/subscriptor.h>

#include <memory>
#include <vector>

DEAL_II_NAMESPACE_OPEN

namespace hp
{
  /**
   * 这个类实现了一个对象的集合。
   * 它实现了doxygen文档中描述的 @ref hpcollection
   * 模块中的概念。
   * @ingroup hp hpcollection
   *
   */
  template <typename T>
  class Collection : public Subscriptor
  {
  public:
    /**
     * 默认构造函数。导致一个空的集合，以后可以用push_back()来填充。
     *
     */
    Collection() = default;

    /**
     * 添加一个新的对象。
     *
     */
    void
    push_back(const std::shared_ptr<const T> &new_entry);

    /**
     * 返回用户为活动FE索引指定的对象，该索引是作为参数提供给该方法的。
     * @pre   @p index 必须在0和集合的元素数之间。
     *
     */
    const T &operator[](const unsigned int index) const;

    /**
     * 返回存储在这个容器中的对象的数量。
     *
     */
    unsigned int
    size() const;

    /**
     * 确定这个对象的内存消耗（以字节为单位）的估计值。
     *
     */
    std::size_t
    memory_consumption() const;

  private:
    /**
     * 真正的容器，它存储了指向不同对象的指针。
     *
     */
    std::vector<std::shared_ptr<const T>> entries;
  };


   /* --------------- inline functions ------------------- */ 



  template <typename T>
  std::size_t
  Collection<T>::memory_consumption() const
  {
    return (sizeof(*this) + MemoryConsumption::memory_consumption(entries));
  }



  template <typename T>
  void
  Collection<T>::push_back(const std::shared_ptr<const T> &new_entry)
  {
    entries.push_back(new_entry);
  }



  template <typename T>
  inline unsigned int
  Collection<T>::size() const
  {
    return entries.size();
  }



  template <typename T>
  inline const T &Collection<T>::operator[](const unsigned int index) const
  {
    AssertIndexRange(index, entries.size());
    return *entries[index];
  }

} // namespace hp


DEAL_II_NAMESPACE_CLOSE

#endif


