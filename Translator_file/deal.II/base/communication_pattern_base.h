//include/deal.II-translator/base/communication_pattern_base_0.txt
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

#ifndef dealii_base_communication_pattern_base_h
#define dealii_base_communication_pattern_base_h

#include <deal.II/base/config.h>

#include <deal.II/base/mpi.h>

DEAL_II_NAMESPACE_OPEN

// Forward declaration
#ifndef DOXYGEN
class IndexSet;
#endif

namespace Utilities
{
  namespace MPI
  {
    /**
     * 通信模式（CommunicationPattern）是一个抽象类，用于定义一个可以重复调用的通信计划，以有效地获得非处理器元素。其想法是将通信模式与需要通信的数据解耦。目标是为不同的容器重复使用相同的通信模式。这与SparseMatrix和SparsityPattern的工作方式类似。
     *
     */
    class CommunicationPatternBase
    {
    public:
      /**
       * 解构器。
       *
       */
      virtual ~CommunicationPatternBase() = default;

      /**
       * 重新初始化通信模式。第一个参数`vector_space_vector_index_set`是与一个VectorSpaceVector对象相关的索引集。第二个参数`read_write_vector_index_set`是与一个ReadWriteVector对象相关的索引集。
       *
       */
      virtual void
      reinit(const IndexSet &vector_space_vector_index_set,
             const IndexSet &read_write_vector_index_set,
             const MPI_Comm &communicator) = 0;

      /**
       * 返回一个对底层MPI通信器的常数引用。
       *
       */
      virtual const MPI_Comm &
      get_mpi_communicator() const = 0;
    };

  } // namespace MPI

} // namespace Utilities

DEAL_II_NAMESPACE_CLOSE

#endif


