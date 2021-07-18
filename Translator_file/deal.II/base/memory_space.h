//include/deal.II-translator/base/memory_space_0.txt
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


#ifndef dealii_memory_space_h
#define dealii_memory_space_h

#include <deal.II/base/config.h>

DEAL_II_NAMESPACE_OPEN

 /**
  */
namespace MemorySpace
{
  /**
   * 描述主机内存空间的结构。
   *
   */
  struct Host
  {};



  /**
   * 描述CUDA内存空间的结构。
   *
   */
  struct CUDA
  {};

} // namespace MemorySpace

DEAL_II_NAMESPACE_CLOSE

#endif


