//include/deal.II-translator/lac/trilinos_epetra_communication_pattern_0.txt
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

#ifndef dealii_trilinos_epetra_communication_pattern_h
#define dealii_trilinos_epetra_communication_pattern_h


#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_TRILINOS

#  ifdef DEAL_II_WITH_MPI

#    include <deal.II/base/communication_pattern_base.h>

#    include <Epetra_Import.h>

#    include <memory>

DEAL_II_NAMESPACE_OPEN

namespace LinearAlgebra
{
  namespace EpetraWrappers
  {
    /**
     * 该类实现了对Trilinos
     * Epetra_Import对象的包装器，用于需要
     * Utilities::MPI::CommunicationPatternBase 对象的地方。
     *
     */
    class CommunicationPattern : public Utilities::MPI::CommunicationPatternBase
    {
    public:
      /**
       * 初始化通信模式。第一个参数 @p
       * vector_space_vector_index_set是与一个VectorSpaceVector对象相关的索引集。第二个参数
       * @p
       * read_write_vector_index_set是与一个ReadWriteVector对象相关的索引集。
       *
       */
      CommunicationPattern(const IndexSet &vector_space_vector_index_set,
                           const IndexSet &read_write_vector_index_set,
                           const MPI_Comm &communicator);

      /**
       * 重新初始化该对象。
       *
       */
      virtual void
      reinit(const IndexSet &vector_space_vector_index_set,
             const IndexSet &read_write_vector_index_set,
             const MPI_Comm &communicator) override;

      /**
       * 返回底层的MPI通信器。
       *
       */
      virtual const MPI_Comm &
      get_mpi_communicator() const override;

      /**
       * 返回底层的Epetra_Import对象。
       *
       */
      const Epetra_Import &
      get_epetra_import() const;

    private:
      /**
       * 用于MPI通信器的共享指针。
       *
       */
      std::shared_ptr<const MPI_Comm> comm;

      /**
       * 使用的Epetra_Import对象的共享指针。
       *
       */
      std::unique_ptr<Epetra_Import> importer;
    };
  } // end of namespace EpetraWrappers
} // end of namespace LinearAlgebra

DEAL_II_NAMESPACE_CLOSE

#  endif

#endif

#endif


