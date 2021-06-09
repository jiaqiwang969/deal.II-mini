//include/deal.II-translator/base/memory_consumption_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2020 by the deal.II authors
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

#ifndef dealii_memory_consumption_h
#define dealii_memory_consumption_h


#include <deal.II/base/config.h>

#include <array>
#include <complex>
#include <cstddef>
#include <cstring>
#include <memory>
#include <string>
#include <type_traits>
#include <vector>

DEAL_II_NAMESPACE_OPEN

/**
 * 这个命名空间提供了帮助确定对象所使用的内存量的函数。其目的不一定是要给出到最后一位的内存用量（
 * <tt>std::map</tt>
 * 对象使用的内存是多少？），而是要帮助寻找内存瓶颈。
 * 这个命名空间有一个成员函数memory_consumption()和大量的特殊化。根据函数的参数类型，有几种操作模式。
 * <ol>   <li>  如果参数是一个基本的C++数据类型（如<tt>bool</tt>, <tt>float</tt>, <tt>double</tt>或任何一个整数类型），那么 memory_consumption() 只是返回其参数的<tt>sizeof</tt>。该库还提供了对一个 <tt>std::string</tt>. * <li> 所占用的内存量的估计。
 * <li>
 * 对于既不是标准类型，也不是向量的对象，memory_consumption()将简单地调用同名的成员函数。这取决于数据类型的实现，以提供对所使用的内存量的良好估计。在这个函数中，对类的化合物使用
 * MemoryConsumption::memory_consumption()
 * 有助于获得这种估计。deal.II库中的大多数类都有这样一个成员函数。
 * <li>
 * 对于对象的向量和C++数组，memory_consumption()对所有条目进行递归调用，并将结果加入对象本身的大小。存在一些针对标准数据类型的优化特殊化。
 * <li>
 * 对于普通指针的向量，memory_consumption(T*)返回指针向量的大小，忽略对象的大小。
 * </ol>
 * <h3>Extending this namespace</h3>
 * 这个名字空间中的函数和它所提供的功能依赖于这样的假设：在这个名字空间中，有一个函数<tt>memory_consumption(T)</tt>确定<tt>T</tt>类型的对象所使用的内存量，或者<tt>T</tt>类有一个这个名字的成员函数。虽然后者对deal.II中几乎所有的类都是如此，但我们只为最常见的数据类型实现了第一类函数，如基本类型、字符串、C++向量、C-style数组和C++对。因此，这些函数并不包括，例如，C++地图、列表等。如果你需要这样的函数，请随时实现它们，并将它们发送给我们，以便纳入。
 *
 *
 * @ingroup memory  Wells
 *
 *
 */
namespace MemoryConsumption
{
  /**
   * 计算一个基本类型的内存消耗。参见EnableIfScalar，讨论如何实现这一限制（SFINAE）。
   *
   */
  template <typename T>
  inline
    typename std::enable_if<std::is_fundamental<T>::value, std::size_t>::type
    memory_consumption(const T &t);

  /**
   * 估计一个对象的内存消耗。如果对<tt>T</tt>类型没有进一步的模板专门化（超过这个），那么这个函数返回成员函数<tt>t.memory_consumption()</tt>的值。
   *
   */
  template <typename T>
  inline typename std::enable_if<!(std::is_fundamental<T>::value ||
                                   std::is_pointer<T>::value),
                                 std::size_t>::type
  memory_consumption(const T &t);

  /**
   * 确定一个C风格的字符串所消耗的内存量。返回的值不包括指针的大小。这个函数只测量到（并包括）NUL字节；底层缓冲区可能更大。
   *
   */
  inline std::size_t
  memory_consumption(const char *string);

  /**
   * 确定一个 <tt>std::complex</tt> 变量所消耗的内存字节数。
   *
   */
  template <typename T>
  inline std::size_t
  memory_consumption(const std::complex<T> &);

  /**
   * 确定一个<tt>VectorizedArray</tt>变量所消耗的内存字节数。
   *
   */
  template <typename T, std::size_t width>
  inline std::size_t
  memory_consumption(const VectorizedArray<T, width> &);

  /**
   * 确定一个 <tt>std::string</tt>
   * 变量所消耗的内存的估计字节数。
   *
   */
  inline std::size_t
  memory_consumption(const std::string &s);

  /**
   * 通过调用每个条目的memory_consumption()，确定一个<tt>T</tt>类型的元素的
   * <tt>std::vector</tt> 所消耗的内存字节量。
   * 这个函数在向量的所有条目上循环，并对每个<tt>v[i]</tt>使用memory_consumption()确定其大小。如果条目的大小是恒定的，可能有另一个全局函数memory_consumption()用于这个数据类型，或者有一个该类的成员函数的名字返回一个恒定的值，编译器会解开这个循环，这样操作就很快了。如果数据元素的大小是可变的，例如它们自己做内存分配，那么这个操作必然会更加昂贵。
   * 使用该算法，特别是所有元素的循环，也可以计算向量的向量、字符串的向量等的内存消耗，其中单个元素的大小可能有很大不同。
   * 请注意，这个算法也考虑到了被这个向量分配但目前没有使用的元素的大小。
   * 对于最常用的向量，有一些特殊的函数可以不通过循环来计算其大小。这也适用于bools向量的特殊情况。
   *
   */
  template <typename T>
  inline std::size_t
  memory_consumption(const std::vector<T> &v);

  /**
   * 通过对每个条目调用memory_consumption()，确定一个<tt>N</tt>类型<tt>元素的
   * <tt>std::array</tt> 所消耗的内存字节数。
   * 这个函数循环遍历数组的所有条目，并对每个<tt>v[i]</tt>使用memory_consumption()确定其大小。如果条目的大小是恒定的，可能有另一个全局函数memory_consumption()用于该数据类型，或者有一个该类的成员函数的名字返回一个恒定的值，编译器会解开这个循环，这样操作就会很快。如果数据元素的大小是可变的，例如它们自己做内存分配，那么这个操作必然会更加昂贵。
   * 使用该算法，特别是所有元素的循环，也可以计算向量数组、字符串数组等的内存消耗，其中单个元素的大小可能有很大不同。
   *
   */
  template <typename T, std::size_t N>
  inline std::size_t
  memory_consumption(const std::array<T, N> &v);

  /**
   * 估计一个C型数组所占用的内存量（以字节为单位）。
   * 由于在这个库中，我们通常不在这种数组中存储简单的数据元素，如<tt>double</tt>s（而是使用
   * <tt>std::vector</tt>s 或deal.II
   * <tt>Vector</tt>对象），我们不提供像 <tt>std::vector</tt>
   * 数组的特殊化，而是始终使用所有元素的循环。
   *
   */
  template <typename T, int N>
  inline std::size_t
  memory_consumption(const T (&v)[N]);

  /**
   * 确定一个向量的内存消耗的特殊化，这里是针对一个<tt>bool</tt>s的向量。
   * 这是一个特殊情况，因为bool不是一个一个地存储，而是作为一个位域。
   *
   */
  inline std::size_t
  memory_consumption(const std::vector<bool> &v);

  /**
   * 确定一对数值所消耗的内存字节数的估计值。
   *
   */
  template <typename A, typename B>
  inline std::size_t
  memory_consumption(const std::pair<A, B> &p);

  /**
   * 计算一个指针的内存消耗。
   * @note
   * 这个函数对于C风格的字符串是重载的；关于这种情况，请看该函数的文档。
   * @note
   * 这个函数返回指针的大小，而不是指向的对象的大小。
   *
   */
  template <typename T>
  inline std::size_t
  memory_consumption(const T *const);

  /**
   * 返回一个共享指针所使用的内存量。
   * @note
   * 这将返回指针的大小，而不是所指向的对象的大小。
   *
   */
  template <typename T>
  inline std::size_t
  memory_consumption(const std::shared_ptr<T> &);

  /**
   * 返回一个 std::unique_ptr 对象所使用的内存量。
   * @note  这将返回指针的大小，而不是指向的对象的大小。
   *
   */
  template <typename T>
  inline std::size_t
  memory_consumption(const std::unique_ptr<T> &);
} // namespace MemoryConsumption



// now comes the implementation of these functions

namespace MemoryConsumption
{
  template <typename T>
  inline
    typename std::enable_if<std::is_fundamental<T>::value, std::size_t>::type
    memory_consumption(const T &)
  {
    return sizeof(T);
  }



  inline std::size_t
  memory_consumption(const char *string)
  {
    if (string == nullptr)
      {
        return 0;
      }
    else
      {
        return sizeof(char) * (strlen(string)  /*Remember the NUL*/  + 1);
      }
  }



  template <typename T>
  inline std::size_t
  memory_consumption(const std::complex<T> &)
  {
    return sizeof(std::complex<T>);
  }



  template <typename T, std::size_t width>
  inline std::size_t
  memory_consumption(const VectorizedArray<T, width> &)
  {
    return sizeof(VectorizedArray<T, width>);
  }



  inline std::size_t
  memory_consumption(const std::string &s)
  {
    return sizeof(s) + s.length();
  }



  template <typename T>
  std::size_t
  memory_consumption(const std::vector<T> &v)
  {
    // shortcut for types that do not allocate memory themselves
    if (std::is_fundamental<T>::value || std::is_pointer<T>::value)
      {
        return v.capacity() * sizeof(T) + sizeof(v);
      }
    else
      {
        std::size_t mem = sizeof(std::vector<T>);
        for (unsigned int i = 0; i < v.size(); ++i)
          {
            mem += memory_consumption(v[i]);
          }
        mem += (v.capacity() - v.size()) * sizeof(T);
        return mem;
      }
  }



  template <typename T, std::size_t N>
  std::size_t
  memory_consumption(const std::array<T, N> &v)
  {
    // shortcut for types that do not allocate memory themselves
    if (std::is_fundamental<T>::value || std::is_pointer<T>::value)
      {
        return sizeof(v);
      }
    else
      {
        std::size_t mem = 0;
        for (std::size_t i = 0; i != N; ++i)
          mem += memory_consumption(v[i]);
        return mem;
      }
  }



  template <typename T, int N>
  std::size_t
  memory_consumption(const T (&v)[N])
  {
    std::size_t mem = 0;
    for (unsigned int i = 0; i < N; ++i)
      mem += memory_consumption(v[i]);
    return mem;
  }



  inline std::size_t
  memory_consumption(const std::vector<bool> &v)
  {
    return v.capacity() / 8 + sizeof(v);
  }



  template <typename A, typename B>
  inline std::size_t
  memory_consumption(const std::pair<A, B> &p)
  {
    return (memory_consumption(p.first) + memory_consumption(p.second));
  }



  template <typename T>
  inline std::size_t
  memory_consumption(const T *const)
  {
    return sizeof(T *);
  }



  template <typename T>
  inline std::size_t
  memory_consumption(const std::shared_ptr<T> &)
  {
    return sizeof(std::shared_ptr<T>);
  }



  template <typename T>
  inline std::size_t
  memory_consumption(const std::unique_ptr<T> &)
  {
    return sizeof(std::unique_ptr<T>);
  }



  template <typename T>
  inline typename std::enable_if<!(std::is_fundamental<T>::value ||
                                   std::is_pointer<T>::value),
                                 std::size_t>::type
  memory_consumption(const T &t)
  {
    return t.memory_consumption();
  }
} // namespace MemoryConsumption

DEAL_II_NAMESPACE_CLOSE

#endif


