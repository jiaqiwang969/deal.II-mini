//include/deal.II-translator/base/mpi_noncontiguous_partitioner_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2020 - 2021 by the deal.II authors
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

#ifndef dealii_mpi_noncontiguous_partitioner_h
#define dealii_mpi_noncontiguous_partitioner_h

#include <deal.II/base/config.h>

#include <deal.II/base/communication_pattern_base.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/mpi_compute_index_owner_internal.h>
#include <deal.II/base/mpi_tags.h>

#include <deal.II/lac/vector_space_vector.h>


DEAL_II_NAMESPACE_OPEN

namespace Utilities
{
  namespace MPI
  {
    /**
     * 一个灵活的Partitioner类，它不对基础索引集的顺序进行限制。
     *
     */
    class NoncontiguousPartitioner
      : public Utilities::MPI::CommunicationPatternBase
    {
    public:
      /**
       * 默认构造函数。需要调用reinit()函数之一来创建一个有效的对象。
       *
       */
      NoncontiguousPartitioner() = default;

      /**
       * 构造函数。根据MPI通信器 @p indexset_locally_owned 和 @p
       * indexset_ghost 的IndexSets参数设置点对点通信模式。
       *
       */
      NoncontiguousPartitioner(const IndexSet &indexset_locally_owned,
                               const IndexSet &indexset_ghost,
                               const MPI_Comm &communicator);

      /**
       * 构造函数。与上述相同，但对于索引为 @p
       * indices_locally_owned 和 @p indices_ghost.
       * 的向量，这允许索引不被排序，并且在update_values()、update_values_start()和update_values_finish()期间自动在向量的正确位置读写值。允许包括值为
       * numbers::invalid_dof_index
       * 的条目，这些条目不占索引交换的一部分，但作为填充存在于数据向量中。
       *
       */
      NoncontiguousPartitioner(
        const std::vector<types::global_dof_index> &indices_locally_owned,
        const std::vector<types::global_dof_index> &indices_ghost,
        const MPI_Comm &                            communicator);

      /**
       * 根据预先计算的通信模式，用来自 @p locally_owned_array.
       * @pre 的值填充向量 @p ghost_array
       * ，向量只需提供一个方法begin()，允许访问其原始数据。
       * @pre
       * 两个向量的大小必须至少与传递给构造函数或reinit()函数的索引集的条目数一样大。
       * @note  这个函数依次调用了 update_values_start() 和
       * update_values_finish()
       * 方法。用户可以单独调用这两个函数，从而使通信和计算重叠。
       *
       */
      template <typename Number>
      void
      export_to_ghosted_array(
        const ArrayView<const Number> &locally_owned_array,
        const ArrayView<Number> &      ghost_array) const;

      /**
       * 同上，但接口与
       * Utilities::MPI::Partitioner::export_to_ghosted_array_start 和
       * Utilities::MPI::Partitioner::export_to_ghosted_array_finish. 类似
       * 在这个函数中，用户可以提供要使用的临时数据结构。
       * @pre   @p temporary_storage
       * 向量的大小必须至少是temporary_storage_size。其原因是这个向量被用作发送和接收数据的缓冲区。
       * @note 任何小于10的值都是 @p communication_channel.
       * 的有效值。
       *
       */
      template <typename Number>
      void
      export_to_ghosted_array(
        const unsigned int             communication_channel,
        const ArrayView<const Number> &locally_owned_array,
        const ArrayView<Number> &      temporary_storage,
        const ArrayView<Number> &      ghost_array,
        std::vector<MPI_Request> &     requests) const;

      /**
       * 开始更新：数据被打包，非阻塞的发送和接收开始。
       * @note  与函数
       * Utilities::MPI::Partitioner::export_to_ghosted_array_start,
       * 相比，用户不传递对目标向量的引用，因为数据被接收到缓冲区的指定部分
       * @p temporary_storage.
       * 这允许对接收的数据进行填充和其他后处理。
       * @pre  向量的要求大小与上面的函数相同。
       * @note  任何小于10的值都是 @p communication_channel.
       * 的有效值。
       *
       */
      template <typename Number>
      void
      export_to_ghosted_array_start(
        const unsigned int             communication_channel,
        const ArrayView<const Number> &locally_owned_array,
        const ArrayView<Number> &      temporary_storage,
        std::vector<MPI_Request> &     requests) const;

      /**
       * 完成更新。该方法一直等待，直到所有数据都被发送和接收。一旦收到来自任何进程的数据，它将被处理并放置在向量
       * @p dst. 的正确位置。
       * @note  与函数
       * Utilities::MPI::Partitioner::export_to_ghosted_array_finish,
       * 相比，用户还必须传递对缓冲区 @p temporary_storage,
       * 的引用，因为数据已经被接收到缓冲区而不是目标向量中。
       * @pre  向量的要求大小与上面的函数相同。
       *
       */
      template <typename Number>
      void
      export_to_ghosted_array_finish(
        const ArrayView<const Number> &temporary_storage,
        const ArrayView<Number> &      ghost_array,
        std::vector<MPI_Request> &     requests) const;

      /**
       * 返回本进程发送数据的进程数和本进程接收数据的进程数。
       *
       */
      std::pair<unsigned int, unsigned int>
      n_targets() const;

      /**
       * 返回export_to_ghosted_array()函数所需的临时存储的大小，如果临时存储是由用户代码处理的。
       *
       */
      unsigned int
      temporary_storage_size() const;

      /**
       * 以Byte为单位返回内存消耗。
       *
       */
      types::global_dof_index
      memory_consumption();

      /**
       * 返回底层通信器。
       *
       */
      const MPI_Comm &
      get_mpi_communicator() const override;

      /**
       * 初始化内部数据结构。
       *
       */
      void
      reinit(const IndexSet &indexset_locally_owned,
             const IndexSet &indexset_ghost,
             const MPI_Comm &communicator) override;

      /**
       * 初始化内部数据结构。
       *
       */
      void
      reinit(const std::vector<types::global_dof_index> &indices_locally_owned,
             const std::vector<types::global_dof_index> &indices_ghost,
             const MPI_Comm &                            communicator);

    private:
      /**
       * MPI通信器。
       *
       */
      MPI_Comm communicator;

      /**
       * 这个进程发送数据的行列。
       *
       */
      std::vector<unsigned int> send_ranks;

      /**
       * send_buffer内每个进程的偏移量。
       * @note 与`send_indices`一起构成一个CRS数据结构。
       *
       */
      std::vector<types::global_dof_index> send_ptr;

      /**
       * send_buffer中每个条目在目标向量中的本地索引。
       * @note  与`send_ptr`一起构成一个CRS数据结构。
       *
       */
      std::vector<types::global_dof_index> send_indices;

      /**
       * 这个进程接收数据的行列。
       *
       */
      std::vector<unsigned int> recv_ranks;

      /**
       * recv_buffer中每个进程的偏移量。
       * @note 与`recv_indices`一起构成一个CRS数据结构。
       *
       */
      std::vector<types::global_dof_index> recv_ptr;

      /**
       * recv_buffer中每个条目在目标向量中的本地索引。
       * @note 与`recv_ptr`一起构成一个CRS数据结构。
       *
       */
      std::vector<types::global_dof_index> recv_indices;

      /**
       * 包含按等级排序的数值的缓冲区，用于发送和接收。
       * @note  只有在用户不提供外部的情况下才会分配。
       * @note
       * 在这个地方我们不知道要发送的数据的类型。所以我们使用一个大小为1字节的任意类型。在相关的函数中，该类型被转换为要求的类型。
       *
       */
      mutable std::vector<uint8_t> buffers;

      /**
       * 发送和接收的MPI请求。
       * @note 只有在用户不从外部提供的情况下才分配。
       *
       */
      mutable std::vector<MPI_Request> requests;
    };

  } // namespace MPI
} // namespace Utilities

DEAL_II_NAMESPACE_CLOSE

#endif


