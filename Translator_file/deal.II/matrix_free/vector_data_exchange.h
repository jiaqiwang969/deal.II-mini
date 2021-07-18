//include/deal.II-translator/matrix_free/vector_data_exchange_0.txt
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


#ifndef dealii_matrix_free_vector_data_exchange_h
#define dealii_matrix_free_vector_data_exchange_h


#include <deal.II/base/config.h>

#include <deal.II/base/partitioner.h>


DEAL_II_NAMESPACE_OPEN

namespace internal
{
  namespace MatrixFreeFunctions
  {
    /**
     * 包含MatrixFree中用于进程间数据交换（即用于update_ghost_values和compress）的类的命名空间。
     *
     */
    namespace VectorDataExchange
    {
      /**
       * MatrixFree所需的接口。
       *
       */
      class Base
      {
      public:
        virtual ~Base() = default;

        virtual unsigned int
        locally_owned_size() const = 0;

        virtual unsigned int
        n_ghost_indices() const = 0;

        virtual unsigned int
        n_import_indices() const = 0;

        virtual unsigned int
        n_import_sm_procs() const = 0;

        virtual types::global_dof_index
        size() const = 0;

        virtual void
        export_to_ghosted_array_start(
          const unsigned int                          communication_channel,
          const ArrayView<const double> &             locally_owned_array,
          const std::vector<ArrayView<const double>> &shared_arrays,
          const ArrayView<double> &                   ghost_array,
          const ArrayView<double> &                   temporary_storage,
          std::vector<MPI_Request> &                  requests) const = 0;

        virtual void
        export_to_ghosted_array_finish(
          const ArrayView<const double> &             locally_owned_array,
          const std::vector<ArrayView<const double>> &shared_arrays,
          const ArrayView<double> &                   ghost_array,
          std::vector<MPI_Request> &                  requests) const = 0;

        virtual void
        import_from_ghosted_array_start(
          const VectorOperation::values               vector_operation,
          const unsigned int                          communication_channel,
          const ArrayView<const double> &             locally_owned_array,
          const std::vector<ArrayView<const double>> &shared_arrays,
          const ArrayView<double> &                   ghost_array,
          const ArrayView<double> &                   temporary_storage,
          std::vector<MPI_Request> &                  requests) const = 0;

        virtual void
        import_from_ghosted_array_finish(
          const VectorOperation::values               vector_operation,
          const ArrayView<double> &                   locally_owned_storage,
          const std::vector<ArrayView<const double>> &shared_arrays,
          const ArrayView<double> &                   ghost_array,
          const ArrayView<const double> &             temporary_storage,
          std::vector<MPI_Request> &                  requests) const = 0;

        virtual void
        reset_ghost_values(const ArrayView<double> &ghost_array) const = 0;

        virtual void
        export_to_ghosted_array_start(
          const unsigned int                         communication_channel,
          const ArrayView<const float> &             locally_owned_array,
          const std::vector<ArrayView<const float>> &shared_arrays,
          const ArrayView<float> &                   ghost_array,
          const ArrayView<float> &                   temporary_storage,
          std::vector<MPI_Request> &                 requests) const = 0;

        virtual void
        export_to_ghosted_array_finish(
          const ArrayView<const float> &             locally_owned_array,
          const std::vector<ArrayView<const float>> &shared_arrays,
          const ArrayView<float> &                   ghost_array,
          std::vector<MPI_Request> &                 requests) const = 0;

        virtual void
        import_from_ghosted_array_start(
          const VectorOperation::values              vector_operation,
          const unsigned int                         communication_channel,
          const ArrayView<const float> &             locally_owned_array,
          const std::vector<ArrayView<const float>> &shared_arrays,
          const ArrayView<float> &                   ghost_array,
          const ArrayView<float> &                   temporary_storage,
          std::vector<MPI_Request> &                 requests) const = 0;

        virtual void
        import_from_ghosted_array_finish(
          const VectorOperation::values              vector_operation,
          const ArrayView<float> &                   locally_owned_storage,
          const std::vector<ArrayView<const float>> &shared_arrays,
          const ArrayView<float> &                   ghost_array,
          const ArrayView<const float> &             temporary_storage,
          std::vector<MPI_Request> &                 requests) const = 0;

        virtual void
        reset_ghost_values(const ArrayView<float> &ghost_array) const = 0;
      };


      /**
       * 简单地将任务委托给一个 Utilities::MPI::Partitioner.
       * 的类。
       *
       */
      class PartitionerWrapper : public Base
      {
      public:
        PartitionerWrapper(
          const std::shared_ptr<const Utilities::MPI::Partitioner>
            &partitioner);

        virtual ~PartitionerWrapper() = default;

        unsigned int
        locally_owned_size() const override;

        unsigned int
        n_ghost_indices() const override;

        unsigned int
        n_import_indices() const override;

        unsigned int
        n_import_sm_procs() const override;

        types::global_dof_index
        size() const override;

        void
        export_to_ghosted_array_start(
          const unsigned int                          communication_channel,
          const ArrayView<const double> &             locally_owned_array,
          const std::vector<ArrayView<const double>> &shared_arrays,
          const ArrayView<double> &                   ghost_array,
          const ArrayView<double> &                   temporary_storage,
          std::vector<MPI_Request> &                  requests) const override;

        void
        export_to_ghosted_array_finish(
          const ArrayView<const double> &             locally_owned_array,
          const std::vector<ArrayView<const double>> &shared_arrays,
          const ArrayView<double> &                   ghost_array,
          std::vector<MPI_Request> &                  requests) const override;

        void
        import_from_ghosted_array_start(
          const VectorOperation::values               vector_operation,
          const unsigned int                          communication_channel,
          const ArrayView<const double> &             locally_owned_array,
          const std::vector<ArrayView<const double>> &shared_arrays,
          const ArrayView<double> &                   ghost_array,
          const ArrayView<double> &                   temporary_storage,
          std::vector<MPI_Request> &                  requests) const override;

        void
        import_from_ghosted_array_finish(
          const VectorOperation::values               vector_operation,
          const ArrayView<double> &                   locally_owned_storage,
          const std::vector<ArrayView<const double>> &shared_arrays,
          const ArrayView<double> &                   ghost_array,
          const ArrayView<const double> &             temporary_storage,
          std::vector<MPI_Request> &                  requests) const override;

        void
        reset_ghost_values(const ArrayView<double> &ghost_array) const override;

        void
        export_to_ghosted_array_start(
          const unsigned int                         communication_channel,
          const ArrayView<const float> &             locally_owned_array,
          const std::vector<ArrayView<const float>> &shared_arrays,
          const ArrayView<float> &                   ghost_array,
          const ArrayView<float> &                   temporary_storage,
          std::vector<MPI_Request> &                 requests) const override;

        void
        export_to_ghosted_array_finish(
          const ArrayView<const float> &             locally_owned_array,
          const std::vector<ArrayView<const float>> &shared_arrays,
          const ArrayView<float> &                   ghost_array,
          std::vector<MPI_Request> &                 requests) const override;

        void
        import_from_ghosted_array_start(
          const VectorOperation::values              vector_operation,
          const unsigned int                         communication_channel,
          const ArrayView<const float> &             locally_owned_array,
          const std::vector<ArrayView<const float>> &shared_arrays,
          const ArrayView<float> &                   ghost_array,
          const ArrayView<float> &                   temporary_storage,
          std::vector<MPI_Request> &                 requests) const override;

        void
        import_from_ghosted_array_finish(
          const VectorOperation::values              vector_operation,
          const ArrayView<float> &                   locally_owned_storage,
          const std::vector<ArrayView<const float>> &shared_arrays,
          const ArrayView<float> &                   ghost_array,
          const ArrayView<const float> &             temporary_storage,
          std::vector<MPI_Request> &                 requests) const override;

        void
        reset_ghost_values(const ArrayView<float> &ghost_array) const override;

      private:
        template <typename Number>
        void
        reset_ghost_values_impl(const ArrayView<Number> &ghost_array) const;

        const std::shared_ptr<const Utilities::MPI::Partitioner> partitioner;
      };



      /**
       * 与上述类似，但使用分区器中的内部数据结构，以便识别处于同一共享内存区域的自由度的索引。
       *
       */
      class Full : public Base
      {
      public:
        Full(
          const std::shared_ptr<const Utilities::MPI::Partitioner> &partitioner,
          const MPI_Comm &communicator_sm);

        unsigned int
        locally_owned_size() const override;

        unsigned int
        n_ghost_indices() const override;

        unsigned int
        n_import_indices() const override;

        virtual unsigned int
        n_import_sm_procs() const override;

        virtual types::global_dof_index
        size() const override;

        const MPI_Comm &
        get_sm_mpi_communicator() const;

        void
        export_to_ghosted_array_start(
          const unsigned int                          communication_channel,
          const ArrayView<const double> &             locally_owned_array,
          const std::vector<ArrayView<const double>> &shared_arrays,
          const ArrayView<double> &                   ghost_array,
          const ArrayView<double> &                   temporary_storage,
          std::vector<MPI_Request> &                  requests) const override;

        void
        export_to_ghosted_array_finish(
          const ArrayView<const double> &             locally_owned_array,
          const std::vector<ArrayView<const double>> &shared_arrays,
          const ArrayView<double> &                   ghost_array,
          std::vector<MPI_Request> &                  requests) const override;

        void
        import_from_ghosted_array_start(
          const VectorOperation::values               vector_operation,
          const unsigned int                          communication_channel,
          const ArrayView<const double> &             locally_owned_array,
          const std::vector<ArrayView<const double>> &shared_arrays,
          const ArrayView<double> &                   ghost_array,
          const ArrayView<double> &                   temporary_storage,
          std::vector<MPI_Request> &                  requests) const override;

        void
        import_from_ghosted_array_finish(
          const VectorOperation::values               vector_operation,
          const ArrayView<double> &                   locally_owned_storage,
          const std::vector<ArrayView<const double>> &shared_arrays,
          const ArrayView<double> &                   ghost_array,
          const ArrayView<const double> &             temporary_storage,
          std::vector<MPI_Request> &                  requests) const override;

        void
        reset_ghost_values(const ArrayView<double> &ghost_array) const override;

        void
        export_to_ghosted_array_start(
          const unsigned int                         communication_channel,
          const ArrayView<const float> &             locally_owned_array,
          const std::vector<ArrayView<const float>> &shared_arrays,
          const ArrayView<float> &                   ghost_array,
          const ArrayView<float> &                   temporary_storage,
          std::vector<MPI_Request> &                 requests) const override;

        void
        export_to_ghosted_array_finish(
          const ArrayView<const float> &             locally_owned_array,
          const std::vector<ArrayView<const float>> &shared_arrays,
          const ArrayView<float> &                   ghost_array,
          std::vector<MPI_Request> &                 requests) const override;

        void
        import_from_ghosted_array_start(
          const VectorOperation::values              vector_operation,
          const unsigned int                         communication_channel,
          const ArrayView<const float> &             locally_owned_array,
          const std::vector<ArrayView<const float>> &shared_arrays,
          const ArrayView<float> &                   ghost_array,
          const ArrayView<float> &                   temporary_storage,
          std::vector<MPI_Request> &                 requests) const override;

        void
        import_from_ghosted_array_finish(
          const VectorOperation::values              vector_operation,
          const ArrayView<float> &                   locally_owned_storage,
          const std::vector<ArrayView<const float>> &shared_arrays,
          const ArrayView<float> &                   ghost_array,
          const ArrayView<const float> &             temporary_storage,
          std::vector<MPI_Request> &                 requests) const override;

        void
        reset_ghost_values(const ArrayView<float> &ghost_array) const override;

      private:
        template <typename Number>
        void
        export_to_ghosted_array_start_impl(
          const unsigned int                          communication_channel,
          const ArrayView<const Number> &             locally_owned_array,
          const std::vector<ArrayView<const Number>> &shared_arrays,
          const ArrayView<Number> &                   ghost_array,
          const ArrayView<Number> &                   temporary_storage,
          std::vector<MPI_Request> &                  requests) const;

        template <typename Number>
        void
        export_to_ghosted_array_finish_impl(
          const ArrayView<const Number> &             locally_owned_array,
          const std::vector<ArrayView<const Number>> &shared_arrays,
          const ArrayView<Number> &                   ghost_array,
          std::vector<MPI_Request> &                  requests) const;

        template <typename Number>
        void
        import_from_ghosted_array_start_impl(
          const VectorOperation::values               vector_operation,
          const unsigned int                          communication_channel,
          const ArrayView<const Number> &             locally_owned_array,
          const std::vector<ArrayView<const Number>> &shared_arrays,
          const ArrayView<Number> &                   ghost_array,
          const ArrayView<Number> &                   temporary_storage,
          std::vector<MPI_Request> &                  requests) const;

        template <typename Number>
        void
        import_from_ghosted_array_finish_impl(
          const VectorOperation::values               vector_operation,
          const ArrayView<Number> &                   locally_owned_storage,
          const std::vector<ArrayView<const Number>> &shared_arrays,
          const ArrayView<Number> &                   ghost_array,
          const ArrayView<const Number> &             temporary_storage,
          std::vector<MPI_Request> &                  requests) const;

        template <typename Number>
        void
        reset_ghost_values_impl(const ArrayView<Number> &ghost_array) const;

      private:
        /**
         * 全局通信器。
         *
         */
        const MPI_Comm comm;

        /**
         * 共享内存子通信器。
         *
         */
        const MPI_Comm comm_sm;

        /**
         * 本地拥有的向量项的数量。
         *
         */
        const unsigned int n_local_elements;

        /**
         * 幽灵向量条目的数量。
         *
         */
        const unsigned int n_ghost_elements;

        /**
         * 全局向量条目的数量。
         *
         */
        const types::global_dof_index n_global_elements;

        /**
         * 一个变量按等级缓存更大的索引集中的鬼魂索引数。
         *
         */
        std::vector<unsigned int> n_ghost_indices_in_larger_set_by_remote_rank;

        /**
         * 为一个IndexSet出现的索引集，该索引集是一个较大的索引集的子集，以压缩的方式出现在每个等级中。
         *
         */
        std::pair<std::vector<unsigned int>,
                  std::vector<std::pair<unsigned int, unsigned int>>>
          ghost_indices_subset_data;

        /**
         * 一个包含信息的数组，我的鬼魂索引属于哪个处理器，在哪个偏移量，以及这些索引的数量。
         *
         */
        std::vector<std::array<unsigned int, 3>> ghost_targets_data;

        /**
         * 向我们发送幽灵数据的处理器的集合和数据字段的长度。
         * @note 结构为ghost_targets_data。
         *
         */
        std::vector<std::array<unsigned int, 3>> import_targets_data;

        /**
         * 一个数组，用于缓存每个MPI等级的导入索引中的块数。其长度为
         * import_indices_data.size()+1。
         * 我们在compress()过程中从远程进程中导入的（本地）索引集，即属于本地范围的其他人的幽灵。
         *
         */
        std::pair<std::vector<unsigned int>,
                  std::vector<std::pair<unsigned int, unsigned int>>>
          import_indices_data;

        /**
         * 共享内存行列，在export_to_ghosted_array_finish()过程中，数据是从这个行列复制过来的。
         *
         */
        std::vector<unsigned int> sm_ghost_ranks;

        /**
         * 在export_to_ghosted_array_finish()过程中从哪里复制数据的索引。
         *
         */
        std::pair<std::vector<unsigned int>,
                  std::vector<std::pair<unsigned int, unsigned int>>>
          sm_export_data;

        /**
         * 在export_to_ghosted_array_finish()过程中，将数据复制到哪里的索引。
         *
         */
        std::pair<std::vector<unsigned int>,
                  std::vector<std::pair<unsigned int, unsigned int>>>
          sm_export_data_this;

        /**
         * 在import_from_ghosted_array_finish()过程中，从哪里复制数据的共享内存等级。
         *
         */
        std::vector<unsigned int> sm_import_ranks;

        /**
         * 在import_from_ghosted_array_finish()过程中从哪里复制数据的索引。
         *
         */
        std::pair<std::vector<unsigned int>,
                  std::vector<std::pair<unsigned int, unsigned int>>>
          sm_import_data;

        /**
         * 在import_from_ghosted_array_finish()过程中，将数据复制到哪里的索引。
         *
         */
        std::pair<std::vector<unsigned int>,
                  std::vector<std::pair<unsigned int, unsigned int>>>
          sm_import_data_this;
      };

    } // namespace VectorDataExchange
  }   // end of namespace MatrixFreeFunctions
} // end of namespace internal

DEAL_II_NAMESPACE_CLOSE

#endif


