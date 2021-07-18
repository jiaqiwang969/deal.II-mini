//include/deal.II-translator/grid/tria_iterator_base_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2018 by the deal.II authors
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

#ifndef dealii_tria_iterator_base_h
#define dealii_tria_iterator_base_h


#include <deal.II/base/config.h>

DEAL_II_NAMESPACE_OPEN

/**
 * 名称空间，其中声明了一个枚举，表示迭代器可以处于的状态。
 *
 *
 * @ingroup Iterators
 *
 */
namespace IteratorState
{
  /**
   * 迭代器可以处于的三种状态：有效、过期和无效。
   *
   */
  enum IteratorStates
  {
    /// Iterator points to a valid object
    valid,
    /// Iterator reached end of container
    past_the_end,
    /// Iterator is invalid, probably due to an error
    invalid
  };
} // namespace IteratorState



DEAL_II_NAMESPACE_CLOSE

#endif


