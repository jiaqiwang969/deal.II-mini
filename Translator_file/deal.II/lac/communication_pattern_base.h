//include/deal.II-translator/lac/communication_pattern_base_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2021 by the deal.II authors
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

#ifndef dealii_lac_communication_pattern_base_h
#define dealii_lac_communication_pattern_base_h

#include <deal.II/base/communication_pattern_base.h>

DEAL_II_NAMESPACE_OPEN

DEAL_II_WARNING(
  "This file is deprecated. Use deal.II/base/communication_pattern_base.h instead!")

namespace LinearAlgebra
{
  /**
   * Utilities::MPI::CommunicationPatternBase. 的别名
   * 该类最初是在LinearAlgebra命名空间中定义的，但现在被用于更普遍的目的。
   *
   */
  using CommunicationPatternBase DEAL_II_DEPRECATED =
    Utilities::MPI::CommunicationPatternBase;
} // end of namespace LinearAlgebra

DEAL_II_NAMESPACE_CLOSE

#endif


