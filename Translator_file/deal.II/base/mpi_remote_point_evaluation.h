//include/deal.II-translator/base/mpi_remote_point_evaluation_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2021 by the deal.II authors
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

#ifndef dealii_mpi_mpi_remote_point_evaluation_h
#define dealii_mpi_mpi_remote_point_evaluation_h

#include <deal.II/base/config.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_tools_cache.h>

DEAL_II_NAMESPACE_OPEN


namespace Utilities
{
  namespace MPI
  {
    /**
     * 用于访问非匹配网格上的值的帮助类。
     * @note
     * 字段的名称是在考虑到evaluation_and_process()方法的情况下选择的。在这里，数量是在指定的任意定位点（甚至在MPI宇宙中的远程进程上）逐个单元计算的，这些值被发送给请求进程，请求进程接收结果并根据各点的情况对结果进行求助。
     *
     */
    template <int dim, int spacedim = dim>
    class RemotePointEvaluation
    {
    public:
      /**
       * 构造器。              @param  tolerance
       * 在reinit()过程中，用于确定传递给类的点周围所有单元的单元坐标的公差。根据问题的不同，可能需要调整公差，以便能够确定一个单元。
       * 浮点运算意味着，一般来说，一个点不会完全位于一个顶点、边缘或面。
       * @param  enforce_unique_mapping
       * 强制执行唯一映射，即点和单元的（一对一）关系。
       * @param  rtree_level 构建边界框时使用的RTree级别。
       *
       */
      RemotePointEvaluation(const double       tolerance              = 1e-6,
                            const bool         enforce_unique_mapping = false,
                            const unsigned int rtree_level            = 0);

      /**
       * 销毁器。
       *
       */
      ~RemotePointEvaluation();

      /**
       * 根据点列表 @p points 和网格描述（ @p tria 和 @p
       * 的映射）设置内部数据结构和通信模式。
       * @warning
       * 这是一个集体调用，需要由通信器中的所有处理器来执行。
       *
       */
      void
      reinit(const std::vector<Point<spacedim>> &points,
             const Triangulation<dim, spacedim> &tria,
             const Mapping<dim, spacedim> &      mapping);

      /**
       * 定位在一个单元中的点的数据。
       *
       */
      struct CellData
      {
        /**
         * 单元的级别和索引。
         *
         */
        std::vector<std::pair<int, int>> cells;

        /**
         * 指向与单元格相关的（参考）点的开始和结束的指针。
         *
         */
        std::vector<unsigned int> reference_point_ptrs;

        /**
         * 区间[0,1]^dim中的参考点。
         *
         */
        std::vector<Point<dim>> reference_point_values;
      };

      /**
       * 在给定的点和三角结构中评估函数 @p evaluation_function
       * 。结果存储在 @p output. 中。
       * @note
       * 如果点到单元格的映射不是一对一的关系（is_map_unique()==false），结果需要借助get_point_ptrs()来处理。如果一个点与一个被多个单元共享的几何实体（例如，顶点）重合，或者一个点在计算域之外，就会出现这种情况。
       * @warning
       * 这是一个集体调用，需要由通信器中的所有处理器执行。
       *
       */
      template <typename T>
      void
      evaluate_and_process(
        std::vector<T> &output,
        std::vector<T> &buffer,
        const std::function<void(const ArrayView<T> &, const CellData &)>
          &evaluation_function) const;

      /**
       * 这个方法是evaluate_and_process()方法的逆过程。它使 @p
       * input, 提供的各点数据在 @p evaluation_function.  @warning
       * 函数中可用。这是一个集体调用，需要由通信器中的所有处理器执行。
       *
       */
      template <typename T>
      void
      process_and_evaluate(
        const std::vector<T> &input,
        std::vector<T> &      buffer,
        const std::function<void(const ArrayView<const T> &, const CellData &)>
          &evaluation_function) const;

      /**
       * 返回一个类似CRS的数据结构，以确定结果对应的一个点的位置和数量。
       *
       */
      const std::vector<unsigned int> &
      get_point_ptrs() const;

      /**
       * 如果点和单元格有一对一的关系，则返回。如果一个点不为任何单元所拥有（该点在域外）或多个单元拥有该点（该点位于相邻单元共享的几何实体上），则不是这种情况。
       *
       */
      bool
      is_map_unique() const;

      /**
       * 返回在reinit()过程中使用的三角测量对象。
       *
       */
      const Triangulation<dim, spacedim> &
      get_triangulation() const;

      /**
       * 返回reinit()过程中使用的Mapping对象。
       *
       */
      const Mapping<dim, spacedim> &
      get_mapping() const;

      /**
       * 返回内部数据结构是否已经设置好，如果是，它们是否仍然有效（并且没有因为三角结构的变化而失效）。
       *
       */
      bool
      is_ready() const;

    private:
      /**
       * 在确定一个点的周围单元时要使用的公差。
       *
       */
      const double tolerance;

      /**
       * 强制执行唯一映射，即点和单元的（一对一）关系。
       *
       */
      const bool enforce_unique_mapping;

      /**
       * 在构建边界框的过程中，要使用RTree级别。
       *
       */
      const unsigned int rtree_level;

      /**
       * 存储三角测量信号的状态。
       *
       */
      boost::signals2::connection tria_signal;

      /**
       * 指示是否调用过reinit()函数的标志，如果是的话，三角剖分此后没有被修改（可能会使通信模式无效）。
       *
       */
      bool ready_flag;

      /**
       * 对在reinit()过程中使用的Triangulation对象的引用。
       *
       */
      SmartPointer<const Triangulation<dim, spacedim>> tria;

      /**
       * 对reinit()过程中使用的Mapping对象的引用。
       *
       */
      SmartPointer<const Mapping<dim, spacedim>> mapping;

      /**
       * 点和单元格的（一对一）关系。
       *
       */
      bool unique_mapping;

      /**
       * 因为对于每个点来说，可以有多个或没有结果，所以这个向量中的指针以类似CRS的方式表示与一个点相关的第一个和最后一个条目。
       *
       */
      std::vector<unsigned int> point_ptrs;

      /**
       * 在recv缓冲区内的互斥索引。
       *
       */
      std::vector<unsigned int> recv_permutation;

      /**
       * 一个接收缓冲区内的范围指针，由recv_ranks指定的等级来填充。
       *
       */
      std::vector<unsigned int> recv_ptrs;

      /**
       * 接收数据的等级。
       *
       */
      std::vector<unsigned int> recv_ranks;

      /**
       * 根据单元格排序的点数据，以便对每个单元格只需进行一次评估（包括自由度的读取）。
       *
       */
      CellData cell_data;

      /**
       * 发送缓冲区内的互换索引。
       *
       */
      std::vector<unsigned int> send_permutation;

      /**
       * 要发送的等级。
       *
       */
      std::vector<unsigned int> send_ranks;

      /**
       * 发送缓冲区内的范围指针，将被发送至send_ranks指定的等级。
       *
       */
      std::vector<unsigned int> send_ptrs;
    };


    template <int dim, int spacedim>
    template <typename T>
    void
    RemotePointEvaluation<dim, spacedim>::evaluate_and_process(
      std::vector<T> &output,
      std::vector<T> &buffer,
      const std::function<void(const ArrayView<T> &, const CellData &)>
        &evaluation_function) const
    {
#ifndef DEAL_II_WITH_MPI
      Assert(false, ExcNeedsMPI());
      (void)output;
      (void)buffer;
      (void)evaluation_function;
#else
      static CollectiveMutex      mutex;
      CollectiveMutex::ScopedLock lock(mutex, tria->get_communicator());

      output.resize(point_ptrs.back());
      buffer.resize(send_permutation.size() * 2);
      ArrayView<T> buffer_1(buffer.data(), buffer.size() / 2);
      ArrayView<T> buffer_2(buffer.data() + buffer.size() / 2,
                            buffer.size() / 2);

      // evaluate functions at points
      evaluation_function(buffer_1, cell_data);

      // sort for communication
      for (unsigned int i = 0; i < send_permutation.size(); ++i)
        buffer_2[send_permutation[i]] = buffer_1[i];

      // process remote quadrature points and send them away
      std::map<unsigned int, std::vector<char>> temp_map;

      std::vector<MPI_Request> requests;
      requests.reserve(send_ranks.size());

      const unsigned int my_rank =
        Utilities::MPI::this_mpi_process(tria->get_communicator());

      std::map<unsigned int, std::vector<T>> temp_recv_map;

      for (unsigned int i = 0; i < send_ranks.size(); ++i)
        {
          if (send_ranks[i] == my_rank)
            {
              // process locally-owned values
              temp_recv_map[my_rank] =
                std::vector<T>(buffer_2.begin() + send_ptrs[i],
                               buffer_2.begin() + send_ptrs[i + 1]);
              continue;
            }

          temp_map[send_ranks[i]] =
            Utilities::pack(std::vector<T>(buffer_2.begin() + send_ptrs[i],
                                           buffer_2.begin() + send_ptrs[i + 1]),
                            false);

          auto &buffer = temp_map[send_ranks[i]];

          requests.push_back(MPI_Request());

          const int ierr = MPI_Isend(buffer.data(),
                                     buffer.size(),
                                     MPI_CHAR,
                                     send_ranks[i],
                                     internal::Tags::remote_point_evaluation,
                                     tria->get_communicator(),
                                     &requests.back());
          AssertThrowMPI(ierr);
        }

      for (const auto recv_rank : recv_ranks)
        {
          if (recv_rank == my_rank)
            continue;

          MPI_Status status;
          int        ierr = MPI_Probe(MPI_ANY_SOURCE,
                               internal::Tags::remote_point_evaluation,
                               tria->get_communicator(),
                               &status);
          AssertThrowMPI(ierr);

          int message_length;
          ierr = MPI_Get_count(&status, MPI_CHAR, &message_length);
          AssertThrowMPI(ierr);

          std::vector<char> buffer(message_length);

          ierr = MPI_Recv(buffer.data(),
                          buffer.size(),
                          MPI_CHAR,
                          status.MPI_SOURCE,
                          internal::Tags::remote_point_evaluation,
                          tria->get_communicator(),
                          MPI_STATUS_IGNORE);
          AssertThrowMPI(ierr);

          temp_recv_map[status.MPI_SOURCE] =
            Utilities::unpack<std::vector<T>>(buffer, false);
        }

      // make sure all messages have been sent
      const int ierr =
        MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
      AssertThrowMPI(ierr);

      // copy received data into output vector
      auto it = recv_permutation.begin();
      for (const auto &j : temp_recv_map)
        for (const auto &i : j.second)
          {
            output[*it] = i;
            it++;
          }
#endif
    }


    template <int dim, int spacedim>
    template <typename T>
    void
    RemotePointEvaluation<dim, spacedim>::process_and_evaluate(
      const std::vector<T> &input,
      std::vector<T> &      buffer,
      const std::function<void(const ArrayView<const T> &, const CellData &)>
        &evaluation_function) const
    {
#ifndef DEAL_II_WITH_MPI
      Assert(false, ExcNeedsMPI());
      (void)input;
      (void)buffer;
      (void)evaluation_function;
#else
      static CollectiveMutex      mutex;
      CollectiveMutex::ScopedLock lock(mutex, tria->get_communicator());

      const auto &ptr = this->get_point_ptrs();

      std::map<unsigned int, std::vector<T>> temp_recv_map;

      for (unsigned int i = 0; i < recv_ranks.size(); ++i)
        temp_recv_map[recv_ranks[i]].resize(recv_ptrs[i + 1] - recv_ptrs[i]);

      const unsigned int my_rank =
        Utilities::MPI::this_mpi_process(tria->get_communicator());

#  ifdef DEBUG
      {
        unsigned int i = 0;

        for (auto &j : temp_recv_map)
          i += j.second.size();

        AssertDimension(recv_permutation.size(), i);
      }
#  endif

      {
        // duplicate data to be able to sort it more easily in the next step
        std::vector<T> buffer_(ptr.back());
        for (unsigned int i = 0, c = 0; i < ptr.size() - 1; ++i)
          {
            const auto n_entries = ptr[i + 1] - ptr[i];

            for (unsigned int j = 0; j < n_entries; ++j, ++c)
              buffer_[c] = input[i];
          }

        // sort data according to the ranks
        auto it = recv_permutation.begin();
        for (auto &j : temp_recv_map)
          for (auto &i : j.second)
            {
              i = buffer_[*it];
              it++;
            }
      }

      // buffer.resize(point_ptrs.back());
      buffer.resize(send_permutation.size() * 2);
      ArrayView<T> buffer_1(buffer.data(), buffer.size() / 2);
      ArrayView<T> buffer_2(buffer.data() + buffer.size() / 2,
                            buffer.size() / 2);

      // process remote quadrature points and send them away
      std::map<unsigned int, std::vector<char>> temp_map;

      std::vector<MPI_Request> requests;
      requests.reserve(recv_ranks.size());

      for (const auto recv_rank : recv_ranks)
        {
          if (recv_rank == my_rank)
            continue;

          temp_map[recv_rank] =
            Utilities::pack(temp_recv_map[recv_rank], false);

          auto &buffer_send = temp_map[recv_rank];

          requests.push_back(MPI_Request());

          const int ierr = MPI_Isend(buffer_send.data(),
                                     buffer_send.size(),
                                     MPI_CHAR,
                                     recv_rank,
                                     internal::Tags::remote_point_evaluation,
                                     tria->get_communicator(),
                                     &requests.back());
          AssertThrowMPI(ierr);
        }

      for (unsigned int i = 0; i < send_ranks.size(); ++i)
        {
          if (send_ranks[i] == my_rank)
            {
              const auto &buffer_send = temp_recv_map[send_ranks[i]];
              // process locally-owned values
              const unsigned int j = std::distance(send_ranks.begin(),
                                                   std::find(send_ranks.begin(),
                                                             send_ranks.end(),
                                                             my_rank));

              AssertDimension(buffer_send.size(),
                              send_ptrs[j + 1] - send_ptrs[j]);

              for (unsigned int i = send_ptrs[j], c = 0; i < send_ptrs[j + 1];
                   ++i, ++c)
                buffer_1[i] = buffer_send[c];

              continue;
            }

          MPI_Status status;
          int        ierr = MPI_Probe(MPI_ANY_SOURCE,
                               internal::Tags::remote_point_evaluation,
                               tria->get_communicator(),
                               &status);
          AssertThrowMPI(ierr);

          int message_length;
          ierr = MPI_Get_count(&status, MPI_CHAR, &message_length);
          AssertThrowMPI(ierr);

          std::vector<char> recv_buffer(message_length);

          ierr = MPI_Recv(recv_buffer.data(),
                          recv_buffer.size(),
                          MPI_CHAR,
                          status.MPI_SOURCE,
                          internal::Tags::remote_point_evaluation,
                          tria->get_communicator(),
                          MPI_STATUS_IGNORE);
          AssertThrowMPI(ierr);

          const auto recv_buffer_unpacked =
            Utilities::unpack<std::vector<T>>(recv_buffer, false);

          auto ptr =
            std::find(send_ranks.begin(), send_ranks.end(), status.MPI_SOURCE);

          Assert(ptr != send_ranks.end(), ExcNotImplemented());

          const unsigned int j = std::distance(send_ranks.begin(), ptr);

          AssertDimension(recv_buffer_unpacked.size(),
                          send_ptrs[j + 1] - send_ptrs[j]);

          for (unsigned int i = send_ptrs[j], c = 0; i < send_ptrs[j + 1];
               ++i, ++c)
            {
              AssertIndexRange(i, buffer_1.size());
              AssertIndexRange(c, recv_buffer_unpacked.size());
              buffer_1[i] = recv_buffer_unpacked[c];
            }
        }

      const int ierr =
        MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
      AssertThrowMPI(ierr);

      // sort for easy access during function call
      for (unsigned int i = 0; i < send_permutation.size(); ++i)
        buffer_2[i] = buffer_1[send_permutation[i]];

      // evaluate function at points
      evaluation_function(buffer_2, cell_data);
#endif
    }

  } // end of namespace MPI
} // end of namespace Utilities


DEAL_II_NAMESPACE_CLOSE

#endif


