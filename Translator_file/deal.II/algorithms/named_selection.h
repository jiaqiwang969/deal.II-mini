//include/deal.II-translator/algorithms/named_selection_0.txt
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

#ifndef dealii_named_selection_h
#define dealii_named_selection_h

#include <deal.II/base/config.h>

#include <deal.II/algorithms/any_data.h>

#include <string>

DEAL_II_NAMESPACE_OPEN

/**
 * 从AnyData中选择与所附名称对应的数据。
 * 给出一个要搜索的名字列表（由add()提供），该类的对象提供一个所选数据的索引列表。
 *
 *
 */
class NamedSelection
{
public:
  /**
   * 在initialize()中提供的 @p data
   * 中添加一个新的名字来进行搜索。
   * @note 名字将被添加到当前列表的末尾。
   *
   */
  void
  add(const std::string &name);


  /**
   * 创建指向AnyData对象的索引向量。
   *
   */
  void
  initialize(const AnyData &data);


  /**
   * 这个对象中的名字的数量。无论之前是否调用过initialize()，都可以使用这个函数。
   *
   */
  unsigned int
  size() const;


  /**
   * 返回提供给上次initialize()的AnyData对象中的相应索引。如果initialize()之前没有被调用，则是一个错误。
   * 指数的顺序与调用add()的顺序相同。
   *
   */
  unsigned int
  operator()(unsigned int i) const;


private:
  /**
   * 选定的名称。
   *
   */
  std::vector<std::string> names;

  /**
   * 由initialize()生成并由operator()访问的索引图。
   *
   */
  std::vector<unsigned int> indices;
};


inline unsigned int
NamedSelection::size() const
{
  return names.size();
}


inline void
NamedSelection::add(const std::string &s)
{
  names.push_back(s);
}


inline unsigned int
NamedSelection::operator()(unsigned int i) const
{
  Assert(indices.size() == names.size(), ExcNotInitialized());

  AssertIndexRange(i, size());

  return indices[i];
}

DEAL_II_NAMESPACE_CLOSE

#endif


