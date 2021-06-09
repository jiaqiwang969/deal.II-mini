//include/deal.II-translator/algorithms/any_data_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2018 by the deal.II authors
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

#ifndef dealii_any_data_h
#define dealii_any_data_h

#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/subscriptor.h>

#include <boost/any.hpp>

#include <algorithm>
#include <typeinfo>
#include <vector>

DEAL_II_NAMESPACE_OPEN

/**
 * 存储任何数量的由标识符字符串访问的任何类型的数据。
 * @todo
 * GK：废除了通过索引访问AnyData的做法，改为使用地图。
 *
 *
 */
class AnyData : public Subscriptor
{
public:
  /// Default constructor for empty object
  AnyData() = default;

  /// Number of stored data objects.
  unsigned int
  size() const;

  /// Add a new data object
  template <typename type>
  void
  add(type entry, const std::string &name);

  /**
   * @brief Merge the data of another AnyData to the end of this object.
   *
   */
  void
  merge(const AnyData &other);

  /**
   * @brief Access to stored data object by name.
   * 找到具有给定名称的对象，尝试将其转换为<tt>type</tt>并返回。如果名字不存在或者转换失败，这个函数会抛出一个异常。如果不希望出现这样的异常，请使用try_read()代替。
   *
   */
  template <typename type>
  type
  entry(const std::string &name);

  /**
   * @brief Read-only access to stored data object by name.
   * 找到具有给定名称的对象，尝试将其转换为<tt>type</tt>并返回。如果名字不存在或者转换失败，这个函数会抛出一个异常。如果不希望出现这样的异常，请使用try_read()代替。
   *
   */
  template <typename type>
  const type
  entry(const std::string &name) const;

  /**
   * @brief Dedicated read only access by name.
   * 对于一个常量对象，这个函数等于
   * entry()。对于一个非常量对象，它强制对数据进行只读访问。特别是，如果没有找到该对象或无法转换类型，它会抛出一个异常。
   * 如果不希望出现这样的异常，请使用try_read()代替。
   * @warning  不要对作为指针的存储对象使用这个函数。
   * 使用read_ptr()代替!
   *
   */
  template <typename type>
  const type
  read(const std::string &name) const;

  /**
   * @brief Dedicated read only access by name for pointer data.
   * 如果存储的数据对象是一个指向常量对象的指针，访问的逻辑就变得相当复杂。也就是说，标准的读取函数可能会失败，这取决于它是一个常量指针还是一个普通指针。
   * 这个函数修复了这个逻辑，并确定该对象不会因为意外而变得可变。
   *
   */
  template <typename type>
  const type *
  read_ptr(const std::string &name) const;

  /**
   * 执行与 read_ptr()
   * 相同的操作，但如果指针不存在，则不抛出异常。而是返回一个空指针。
   *
   */
  template <typename type>
  const type *
  try_read_ptr(const std::string &name) const;

  /**
   * @brief Dedicated read only access by name without exceptions.
   * 这个函数试图在列表中找到名字并返回一个指向相关对象的指针。如果没有找到名字或者对象不能被转换为返回类型，则返回一个空指针。
   *
   */
  template <typename type>
  const type *
  try_read(const std::string &name) const;

  /**
   * 通过索引访问存储的数据对象。
   *
   */
  template <typename type>
  type
  entry(const unsigned int i);

  /// Read-only access to stored data object by index.
  template <typename type>
  const type
  entry(const unsigned int i) const;

  /// Dedicated read only access.
  template <typename type>
  const type
  read(const unsigned int i) const;

  /// Dedicated read only access to pointer object.
  template <typename type>
  const type *
  read_ptr(const unsigned int i) const;

  /// Dedicated read only access to pointer object without exception.
  template <typename type>
  const type *
  try_read_ptr(const unsigned int i) const;

  /// Dedicated read only access without exception.
  template <typename type>
  const type *
  try_read(const unsigned int i) const;

  /// Name of object at index.
  const std::string &
  name(const unsigned int i) const;

  /**
   * @brief Find index of a named object
   * 尝试找到该对象并返回其在列表中的索引。如果没有找到该对象，则抛出一个异常。
   *
   */
  unsigned int
  find(const std::string &name) const;

  /**
   * @brief Try to find index of a named object
   * 试图找到该对象并返回其在列表中的索引。如果没有找到名称，则返回
   * numbers::invalid_unsigned_int 。
   *
   */
  unsigned int
  try_find(const std::string &name) const;

  /// Find out if object is of a certain type
  template <typename type>
  bool
  is_type(const unsigned int i) const;

  /// List the contents to a stream
  template <class StreamType>
  void
  list(StreamType &os) const;

  /// An entry with this name does not exist in the AnyData object.
  DeclException1(ExcNameNotFound,
                 std::string,
                 << "No entry with the name " << arg1 << " exists.");

  /// The requested type and the stored type are different
  DeclException2(ExcTypeMismatch,
                 std::string,
                 std::string,
                 << "The requested type " << arg1 << " and the stored type "
                 << arg2 << " must coincide.");

  /**
   * 异常，表明一个函数期望一个向量有一个特定的名字，但我们在该位置存储了一个不同的名字。
   *
   */
  DeclException2(ExcNameMismatch,
                 int,
                 std::string,
                 << "Name at position " << arg1 << " is not equal to " << arg2
                 << ".");

private:
  /// The stored data
  std::vector<boost::any> data;
  /// The names of the stored data
  std::vector<std::string> names;
};


unsigned int inline AnyData::size() const
{
  AssertDimension(data.size(), names.size());
  return data.size();
}


template <typename type>
inline type
AnyData::entry(const unsigned int i)
{
  AssertIndexRange(i, size());
  type *p = boost::any_cast<type>(&data[i]);
  Assert(p != nullptr,
         ExcTypeMismatch(typeid(type).name(), data[i].type().name()));
  return *p;
}


template <typename type>
inline const type
AnyData::entry(const unsigned int i) const
{
  AssertIndexRange(i, size());
  const type *p = boost::any_cast<type>(&data[i]);
  if (p == nullptr)
    p = boost::any_cast<const type>(&data[i]);
  Assert(p != nullptr,
         ExcTypeMismatch(typeid(type).name(), data[i].type().name()));
  return *p;
}


template <typename type>
inline const type
AnyData::read(const unsigned int i) const
{
  AssertIndexRange(i, size());
  const type *p = boost::any_cast<type>(&data[i]);
  if (p == nullptr)
    p = boost::any_cast<const type>(&data[i]);
  Assert(p != nullptr,
         ExcTypeMismatch(typeid(type).name(), data[i].type().name()));
  return *p;
}


template <typename type>
inline const type *
AnyData::read_ptr(const unsigned int i) const
{
  AssertIndexRange(i, size());
  const type *const *p = boost::any_cast<type *>(&data[i]);
  if (p == nullptr)
    p = boost::any_cast<const type *>(&data[i]);
  Assert(p != nullptr,
         ExcTypeMismatch(typeid(type *).name(), data[i].type().name()));
  return *p;
}


template <typename type>
inline const type *
AnyData::try_read_ptr(const unsigned int i) const
{
  AssertIndexRange(i, size());
  const type *const *p = boost::any_cast<type *>(&data[i]);
  if (p == nullptr)
    p = boost::any_cast<const type *>(&data[i]);
  if (p == nullptr)
    return nullptr;
  return *p;
}


template <typename type>
inline const type *
AnyData::try_read(const unsigned int i) const
{
  AssertIndexRange(i, size());
  const type *p = boost::any_cast<type>(&data[i]);
  if (p == 0)
    p = boost::any_cast<const type>(&data[i]);
  return p;
}


inline const std::string &
AnyData::name(const unsigned int i) const
{
  AssertIndexRange(i, size());
  return names[i];
}


inline unsigned int
AnyData::try_find(const std::string &n) const
{
  std::vector<std::string>::const_iterator it =
    std::find(names.begin(), names.end(), n);

  if (it == names.end())
    return numbers::invalid_unsigned_int;

  return it - names.begin();
}


inline unsigned int
AnyData::find(const std::string &n) const
{
  const unsigned int i = try_find(n);
  Assert(i != numbers::invalid_unsigned_int, ExcNameNotFound(n));

  return i;
}


template <typename type>
inline bool
AnyData::is_type(const unsigned int i) const
{
  return data[i].type() == typeid(type);
}


template <typename type>
inline type
AnyData::entry(const std::string &n)
{
  const unsigned int i = find(n);
  type *             p = boost::any_cast<type>(&data[i]);
  Assert(p != 0, ExcTypeMismatch(typeid(type).name(), data[i].type().name()));
  return *p;
}


template <typename type>
inline const type
AnyData::entry(const std::string &n) const
{
  const unsigned int i = find(n);
  const type *       p = boost::any_cast<type>(&data[i]);
  Assert(p != nullptr,
         ExcTypeMismatch(typeid(type).name(), data[i].type().name()));
  return *p;
}


template <typename type>
inline const type
AnyData::read(const std::string &n) const
{
  const unsigned int i = find(n);
  const type *       p = boost::any_cast<type>(&data[i]);
  Assert(p != 0, ExcTypeMismatch(typeid(type).name(), data[i].type().name()));
  return *p;
}


template <typename type>
inline const type *
AnyData::read_ptr(const std::string &n) const
{
  const unsigned int i = find(n);
  const type *const *p = boost::any_cast<type *>(&data[i]);
  if (p == nullptr)
    p = boost::any_cast<const type *>(&data[i]);
  Assert(p != nullptr,
         ExcTypeMismatch(typeid(type).name(), data[i].type().name()));
  return *p;
}


template <typename type>
inline const type *
AnyData::try_read_ptr(const std::string &n) const
{
  const unsigned int i = try_find(n);
  if (i == numbers::invalid_unsigned_int)
    return 0;

  const type *const *p = boost::any_cast<type *>(&data[i]);
  if (p == 0)
    p = boost::any_cast<const type *>(&data[i]);
  return *p;
}


template <typename type>
inline const type *
AnyData::try_read(const std::string &n) const
{
  // Try to find name
  std::vector<std::string>::const_iterator it =
    std::find(names.begin(), names.end(), n);
  // Return null pointer if not found
  if (it == names.end())
    return nullptr;

  // Compute index and return casted pointer
  unsigned int i = it - names.begin();
  const type * p = boost::any_cast<type>(&data[i]);
  return p;
}


template <typename type>
inline void
AnyData::add(type ent, const std::string &n)
{
  boost::any e = ent;
  data.push_back(e);
  names.push_back(n);
}


inline void
AnyData::merge(const AnyData &other)
{
  for (unsigned int i = 0; i < other.size(); ++i)
    {
      names.push_back(other.names[i]);
      data.push_back(other.data[i]);
    }
}


template <class StreamType>
inline void
AnyData::list(StreamType &os) const
{
  for (unsigned int i = 0; i < names.size(); ++i)
    {
      os << i << '\t' << names[i] << '\t' << data[i].type().name() << std::endl;
    }
}


//----------------------------------------------------------------------//



DEAL_II_NAMESPACE_CLOSE

#endif


