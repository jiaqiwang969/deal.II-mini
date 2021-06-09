//include/deal.II-translator/base/mpi_compute_index_owner_internal_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2019 - 2021 by the deal.II authors
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

#ifndef dealii_base_mpi_compute_index_owner_internal_h
#define dealii_base_mpi_compute_index_owner_internal_h

#include <deal.II/base/config.h>

#include <deal.II/base/mpi.h>
#include <deal.II/base/mpi_consensus_algorithms.h>

DEAL_II_NAMESPACE_OPEN

namespace Utilities
{
  namespace MPI
  {
    namespace internal
    {
      /**
       * 一个用于 Utilities::MPI::compute_index_owner() 和
       * Utilities::MPI::Partitioner::set_ghost_indices().
       * 的内部命名空间。
       *
       */
      namespace ComputeIndexOwner
      {
        /**
         * ConsensusAlgorithms::Process
         * 的特殊化，用于设置字典，即使在IndexSet空间中有不属于任何进程的范围。
         * @note  仅供内部使用。
         *
         */
        class DictionaryPayLoad
          : public ConsensusAlgorithms::Process<
              std::pair<types::global_dof_index, types::global_dof_index>,
              unsigned int>
        {
        public:
          /**
           * 构造函数。
           *
           */
          DictionaryPayLoad(
            const std::map<unsigned int,
                           std::vector<std::pair<types::global_dof_index,
                                                 types::global_dof_index>>>
              &                        buffers,
            std::vector<unsigned int> &actually_owning_ranks,
            const std::pair<types::global_dof_index, types::global_dof_index>
              &                        local_range,
            std::vector<unsigned int> &actually_owning_rank_list)
            : buffers(buffers)
            , actually_owning_ranks(actually_owning_ranks)
            , local_range(local_range)
            , actually_owning_rank_list(actually_owning_rank_list)
          {}

          /**
           * Utilities::MPI::ConsensusAlgorithms::Process::compute_targets().
           * 的实现
           *
           */
          virtual std::vector<unsigned int>
          compute_targets() override
          {
            std::vector<unsigned int> targets;
            for (const auto &rank_pair : buffers)
              targets.push_back(rank_pair.first);

            return targets;
          }

          /**
           * Utilities::MPI::ConsensusAlgorithms::Process::create_request().
           * 的实现
           *
           */
          virtual void
          create_request(const unsigned int other_rank,
                         std::vector<std::pair<types::global_dof_index,
                                               types::global_dof_index>>
                           &send_buffer) override
          {
            send_buffer = this->buffers.at(other_rank);
          }

          /**
           * 实现
           * Utilities::MPI::ConsensusAlgorithms::Process::answer_request(). 
           */
          virtual void
          answer_request(
            const unsigned int                                     other_rank,
            const std::vector<std::pair<types::global_dof_index,
                                        types::global_dof_index>> &buffer_recv,
            std::vector<unsigned int> &request_buffer) override
          {
            (void)request_buffer; // not needed


            // process message: loop over all intervals
            for (auto interval : buffer_recv)
              {
#ifdef DEBUG
                for (types::global_dof_index i = interval.first;
                     i < interval.second;
                     i++)
                  Assert(actually_owning_ranks[i - local_range.first] ==
                           numbers::invalid_unsigned_int,
                         ExcInternalError());
                Assert(interval.first >= local_range.first &&
                         interval.first < local_range.second,
                       ExcInternalError());
                Assert(interval.second > local_range.first &&
                         interval.second <= local_range.second,
                       ExcInternalError());
#endif
                std::fill(actually_owning_ranks.data() + interval.first -
                            local_range.first,
                          actually_owning_ranks.data() + interval.second -
                            local_range.first,
                          other_rank);
              }
            actually_owning_rank_list.push_back(other_rank);
          }

        private:
          const std::map<unsigned int,
                         std::vector<std::pair<types::global_dof_index,
                                               types::global_dof_index>>>
            &buffers;

          std::vector<unsigned int> &actually_owning_ranks;

          const std::pair<types::global_dof_index, types::global_dof_index>
            &local_range;

          std::vector<unsigned int> &actually_owning_rank_list;
        };





        /**
         * 具有基本分区的字典类，以所有MPI等级都知道的固定大小的单一区间为单位，用于两阶段索引查找。
         *
         */
        struct Dictionary
        {
          /**
           * 区间的最小颗粒大小。
           * 我们选择限制两阶段查找的区间的最小尺寸是考虑到以下两个相互冲突的目标。一方面，我们不希望字典中的区间变得太短。对于不均匀的未知数分布（有些等级有几千个未知数，有些则没有），查询的DoFs为
           *
           * ->
           * 字典就涉及到从一个MPI行列向许多其他持有字典间隔的MPI行列发送，导致一些行列必须发送的消息数量过高。另外，较少的较长的区间通常在查找时更有效率。另一方面，范围大小过大会导致相反的效果，即在查找DoFs中，许多消息会进入一个特定的字典所有者中
           *
           * ->
           * 字典。在目前的设置下，在每个MPI等级有1个DoF的情况下，我们最多得到64条消息进入一个MPI等级，这是很合理的低。同时，不均匀的分布可以用最多64条信息来处理，最高可达4096的系数。
           *
           */
          static const unsigned int range_minimum_grain_size = 64;

          /**
           * 一个向量，其条目数与当前进程的字典中的道夫一样多，每个条目包含该道夫的所有者在IndexSet
           * `owned_indices`中的等级。这在索引查找中被查询到，所以我们保留一个扩展的列表。
           *
           */
          std::vector<unsigned int> actually_owning_ranks;

          /**
           * 一个排序的向量，包含出现在`actually_owning_ranks`中的MPI等级。
           *
           */
          std::vector<unsigned int> actually_owning_rank_list;

          /**
           * 用于索引空间分割的每个MPI等级上的字典中的未知数。为了简化索引查询，不需要额外的通信，这个数字在所有MPI等级上都是一样的。
           *
           */
          types::global_dof_index dofs_per_process;

          /**
           * 在字典中表示的全局索引空间的局部范围，由`dofs_per_process`、当前MPI等级和range_minimum_grain_size计算得出。
           *
           */
          std::pair<types::global_dof_index, types::global_dof_index>
            local_range;

          /**
           * 实际大小，计算为dofs_per_process的最小值和索引空间的可能终点。相当于`local_range.second
           *
           * - local_range.first`.
           *
           */
          types::global_dof_index locally_owned_size;

          /**
           * 索引空间的全局大小。
           *
           */
          types::global_dof_index size;

          /**
           * `owned_indices`索引集分布的等级数量。
           *
           */
          unsigned int n_dict_procs_in_owned_indices;

          /**
           * 一个stride，用于在MPI行列中更均匀地分配工作，以防止颗粒大小迫使我们的范围少于我们的处理器。
           *
           */
          unsigned int stride_small_size;

          /**
           * 通过计算全局大小的分区来设置字典，并将本地拥有的范围的等级信息发送给字典部分的所有者。
           *
           */
          void
          reinit(const IndexSet &owned_indices, const MPI_Comm &comm)
          {
            // 1) set up the partition
            this->partition(owned_indices, comm);

#ifdef DEAL_II_WITH_MPI
            unsigned int my_rank = this_mpi_process(comm);

            types::global_dof_index dic_local_received = 0;
            std::map<unsigned int,
                     std::vector<std::pair<types::global_dof_index,
                                           types::global_dof_index>>>
              buffers;

            std::fill(actually_owning_ranks.begin(),
                      actually_owning_ranks.end(),
                      numbers::invalid_subdomain_id);

            // 2) collect relevant processes and process local dict entries
            for (auto interval = owned_indices.begin_intervals();
                 interval != owned_indices.end_intervals();
                 ++interval)
              {
                // Due to the granularity of the dictionary, the interval
                // might be split into several ranges of processor owner
                // ranks. Here, we process the interval by breaking into
                // smaller pieces in terms of the dictionary number.
                std::pair<types::global_dof_index, types::global_dof_index>
                                   index_range(*interval->begin(), interval->last() + 1);
                const unsigned int owner_last =
                  dof_to_dict_rank(interval->last());
                unsigned int owner_first = numbers::invalid_unsigned_int;
                while (owner_first != owner_last)
                  {
                    Assert(index_range.first < index_range.second,
                           ExcInternalError());

                    owner_first = dof_to_dict_rank(index_range.first);

                    // this explicitly picks up the formula of
                    // dof_to_dict_rank, so the two places must be in sync
                    types::global_dof_index next_index =
                      std::min(get_index_offset(owner_first + 1),
                               index_range.second);

                    Assert(next_index > index_range.first, ExcInternalError());

#  ifdef DEBUG
                    // make sure that the owner is the same on the current
                    // interval
                    for (types::global_dof_index i = index_range.first + 1;
                         i < next_index;
                         ++i)
                      AssertDimension(owner_first, dof_to_dict_rank(i));
#  endif

                    // add the interval, either to the local range or into a
                    // buffer to be sent to another processor
                    if (owner_first == my_rank)
                      {
                        std::fill(actually_owning_ranks.data() +
                                    index_range.first - local_range.first,
                                  actually_owning_ranks.data() + next_index -
                                    local_range.first,
                                  my_rank);
                        dic_local_received += next_index - index_range.first;
                        if (actually_owning_rank_list.empty())
                          actually_owning_rank_list.push_back(my_rank);
                      }
                    else
                      buffers[owner_first].emplace_back(index_range.first,
                                                        next_index);

                    index_range.first = next_index;
                  }
              }

            n_dict_procs_in_owned_indices = buffers.size();
            std::vector<MPI_Request> request;

            // Check if index set space is partitioned globally without gaps.
            if (Utilities::MPI::sum(owned_indices.n_elements(), comm) ==
                owned_indices.size())
              {
                // no gaps: setup is simple! Processes send their locally owned
                // indices to the dictionary. The dictionary stores the sending
                // rank for each index. The dictionary knows exactly
                // when it is set up when all indices it is responsible for
                // have been processed.

                request.reserve(n_dict_procs_in_owned_indices);

                // protect the following communication steps using a mutex:
                static CollectiveMutex      mutex;
                CollectiveMutex::ScopedLock lock(mutex, comm);

                const int mpi_tag =
                  Utilities::MPI::internal::Tags::dictionary_reinit;


                // 3) send messages with local dofs to the right dict process
                for (const auto &rank_pair : buffers)
                  {
                    request.push_back(MPI_Request());
                    const int ierr = MPI_Isend(rank_pair.second.data(),
                                               rank_pair.second.size() * 2,
                                               DEAL_II_DOF_INDEX_MPI_TYPE,
                                               rank_pair.first,
                                               mpi_tag,
                                               comm,
                                               &request.back());
                    AssertThrowMPI(ierr);
                  }

                // 4) receive messages until all dofs in dict are processed
                while (this->locally_owned_size != dic_local_received)
                  {
                    // wait for an incoming message
                    MPI_Status status;
                    int        ierr =
                      MPI_Probe(MPI_ANY_SOURCE, mpi_tag, comm, &status);
                    AssertThrowMPI(ierr);

                    // retrieve size of incoming message
                    int number_amount;
                    ierr = MPI_Get_count(&status,
                                         DEAL_II_DOF_INDEX_MPI_TYPE,
                                         &number_amount);
                    AssertThrowMPI(ierr);

                    const auto other_rank = status.MPI_SOURCE;
                    actually_owning_rank_list.push_back(other_rank);

                    // receive message
                    Assert(number_amount % 2 == 0, ExcInternalError());
                    std::vector<std::pair<types::global_dof_index,
                                          types::global_dof_index>>
                      buffer(number_amount / 2);
                    ierr = MPI_Recv(buffer.data(),
                                    number_amount,
                                    DEAL_II_DOF_INDEX_MPI_TYPE,
                                    status.MPI_SOURCE,
                                    status.MPI_TAG,
                                    comm,
                                    MPI_STATUS_IGNORE);
                    AssertThrowMPI(ierr);
                    // process message: loop over all intervals
                    for (auto interval : buffer)
                      {
#  ifdef DEBUG
                        for (types::global_dof_index i = interval.first;
                             i < interval.second;
                             i++)
                          Assert(actually_owning_ranks[i - local_range.first] ==
                                   numbers::invalid_unsigned_int,
                                 ExcInternalError());
                        Assert(interval.first >= local_range.first &&
                                 interval.first < local_range.second,
                               ExcInternalError());
                        Assert(interval.second > local_range.first &&
                                 interval.second <= local_range.second,
                               ExcInternalError());
#  endif

                        std::fill(actually_owning_ranks.data() +
                                    interval.first - local_range.first,
                                  actually_owning_ranks.data() +
                                    interval.second - local_range.first,
                                  other_rank);
                        dic_local_received += interval.second - interval.first;
                      }
                  }
              }
            else
              {
                // with gap: use a ConsensusAlgorithm to determine when all
                // dictionaries have been set up.

                // 3/4) use a ConsensusAlgorithm to send messages with local
                // dofs to the right dict process
                DictionaryPayLoad temp(buffers,
                                       actually_owning_ranks,
                                       local_range,
                                       actually_owning_rank_list);

                ConsensusAlgorithms::Selector<
                  std::pair<types::global_dof_index, types::global_dof_index>,
                  unsigned int>
                  consensus_algo(temp, comm);
                consensus_algo.run();
              }

            std::sort(actually_owning_rank_list.begin(),
                      actually_owning_rank_list.end());

            for (unsigned int i = 1; i < actually_owning_rank_list.size(); ++i)
              Assert(actually_owning_rank_list[i] >
                       actually_owning_rank_list[i - 1],
                     ExcInternalError());

            // 5) make sure that all messages have been sent
            if (request.size() > 0)
              {
                const int ierr = MPI_Waitall(request.size(),
                                             request.data(),
                                             MPI_STATUSES_IGNORE);
                AssertThrowMPI(ierr);
              }

#else
            (void)owned_indices;
            (void)comm;
#endif
          }

          /**
           * 使用`dofs_per_process`将全局dof索引转换为字典中的MPI等级。我们乘以`stride_small_size`，以确保在MPI等级上的平衡，由于晶粒大小。
           *
           */
          unsigned int
          dof_to_dict_rank(const types::global_dof_index i)
          {
            // note: this formula is also explicitly used in
            // get_index_offset(), so keep the two in sync
            return (i / dofs_per_process) * stride_small_size;
          }

          /**
           * 给出一个任意处理器的MPI等级ID，返回该处理器的本地范围开始的索引偏移。
           *
           */
          types::global_dof_index
          get_index_offset(const unsigned int rank)
          {
            return std::min(dofs_per_process *
                              static_cast<types::global_dof_index>(
                                (rank + stride_small_size - 1) /
                                stride_small_size),
                            size);
          }

          /**
           * 给出来自`actually_owning_ranks`的自有索引中的等级，这将返回`actually_owning_rank_list`中的等级索引。
           *
           */
          unsigned int
          get_owning_rank_index(const unsigned int rank_in_owned_indices,
                                const unsigned int guess = 0)
          {
            AssertIndexRange(guess, actually_owning_rank_list.size());
            if (actually_owning_rank_list[guess] == rank_in_owned_indices)
              return guess;
            else
              {
                auto it = std::lower_bound(actually_owning_rank_list.begin(),
                                           actually_owning_rank_list.end(),
                                           rank_in_owned_indices);
                Assert(it != actually_owning_rank_list.end(),
                       ExcInternalError());
                Assert(*it == rank_in_owned_indices, ExcInternalError());
                return it - actually_owning_rank_list.begin();
              }
          }

        private:
          /**
           * 从索引空间的全局大小和等级的数量计算分区。
           *
           */
          void
          partition(const IndexSet &owned_indices, const MPI_Comm &comm)
          {
#ifdef DEAL_II_WITH_MPI
            const unsigned int n_procs = n_mpi_processes(comm);
            const unsigned int my_rank = this_mpi_process(comm);

            size = owned_indices.size();

            Assert(size > 0, ExcNotImplemented());

            dofs_per_process = (size + n_procs - 1) / n_procs;
            if (dofs_per_process < range_minimum_grain_size)
              {
                dofs_per_process  = range_minimum_grain_size;
                stride_small_size = dofs_per_process * n_procs / size;
              }
            else
              stride_small_size = 1;
            local_range.first  = get_index_offset(my_rank);
            local_range.second = get_index_offset(my_rank + 1);

            locally_owned_size = local_range.second - local_range.first;

            actually_owning_ranks = {};
            actually_owning_ranks.resize(locally_owned_size,
                                         numbers::invalid_unsigned_int);
#else
            (void)owned_indices;
            (void)comm;
#endif
          }
        };



        /**
         * 在 Utilities::MPI::compute_index_owner() 和
         * Utilities::MPI::Partitioner::set_ghost_indices() 的背景下，对
         * ConsensusAlgorithms::Process
         * 进行专业化处理，并增加有效载荷。
         *
         */
        class ConsensusAlgorithmsPayload
          : public ConsensusAlgorithms::Process<
              std::pair<types::global_dof_index, types::global_dof_index>,
              unsigned int>
        {
        public:
          /**
           * 构造函数。
           *
           */
          ConsensusAlgorithmsPayload(const IndexSet &owned_indices,
                                     const IndexSet &indices_to_look_up,
                                     const MPI_Comm &comm,
                                     std::vector<unsigned int> &owning_ranks,
                                     const bool track_index_requests = false)
            : owned_indices(owned_indices)
            , indices_to_look_up(indices_to_look_up)
            , comm(comm)
            , my_rank(this_mpi_process(comm))
            , n_procs(n_mpi_processes(comm))
            , track_index_requests(track_index_requests)
            , owning_ranks(owning_ranks)
          {
            dict.reinit(owned_indices, comm);
            requesters.resize(dict.actually_owning_rank_list.size());
          }

          /**
           * 描述本地拥有的空间的索引空间。
           *
           */
          const IndexSet &owned_indices;

          /**
           * 在一个给定等级上是 "幽灵 "的指数，应该从
           * owned_indices中根据其所有者等级来查找。
           *
           */
          const IndexSet &indices_to_look_up;

          /**
           * 底层的MPI通信器。
           *
           */
          const MPI_Comm comm;

          /**
           * 目前的MPI等级。
           *
           */
          const unsigned int my_rank;

          /**
           * 参与MPI通信器 "comm "的行列总数。
           *
           */
          const unsigned int n_procs;

          /**
           * 控制鬼魂所有者的起源是否也应该被存储。如果是，它将被添加到`requesters`中，并可以通过`get_requesters()`查询。
           *
           */
          const bool track_index_requests;

          /**
           * 索引所有者计算的结果。对`indices_to_look_up`中包含的每个索引，这个向量包含`owned_indices`中所有者的MPI等级。
           *
           */
          std::vector<unsigned int> &owning_ranks;

          /**
           * 追踪请求的来源。该数据结构的布局如下。最外层的向量有与
           * Dictionary::actually_owning_rank_list
           * 一样多的条目，代表我们应该从现在的字典条目中送回给所有者的信息。然后，第二个向量收集了一个已请求数据的MPI等级列表，使用第一对条目中的等级和索引范围的列表作为第二条目。
           *
           */
          std::vector<std::vector<
            std::pair<unsigned int,
                      std::vector<std::pair<unsigned int, unsigned int>>>>>
            requesters;

          /**
           * 处理请求的字典。
           *
           */
          Dictionary dict;

          /**
           * 用于收集要查询的索引的数组，按字典中的等级排序。
           *
           */
          std::map<unsigned int, std::vector<types::global_dof_index>>
            indices_to_look_up_by_dict_rank;

          /**
           * 存储从进程中传入数据的索引的字段。
           *
           */
          std::map<unsigned int, std::vector<unsigned int>> recv_indices;

          /**
           * 实现
           * Utilities::MPI::ConsensusAlgorithms::Process::answer_request(),
           * 在request_buffer中添加特定索引的所有者（并跟踪谁请求了一个特定的索引，以防该信息也被需要）。
           *
           */
          virtual void
          answer_request(
            const unsigned int                                     other_rank,
            const std::vector<std::pair<types::global_dof_index,
                                        types::global_dof_index>> &buffer_recv,
            std::vector<unsigned int> &request_buffer) override
          {
            unsigned int owner_index = 0;
            for (const auto &interval : buffer_recv)
              for (auto i = interval.first; i < interval.second; ++i)
                {
                  const unsigned int actual_owner =
                    dict.actually_owning_ranks[i - dict.local_range.first];
                  request_buffer.push_back(actual_owner);

                  if (track_index_requests)
                    append_index_origin(i, owner_index, other_rank);
                }
          }

          /**
           * Utilities::MPI::ConsensusAlgorithms::Process::compute_targets().
           * 的实现
           *
           */
          virtual std::vector<unsigned int>
          compute_targets() override
          {
            std::vector<unsigned int> targets;

            // 1) collect relevant processes and process local dict entries
            {
              unsigned int index       = 0;
              unsigned int owner_index = 0;
              for (auto i : indices_to_look_up)
                {
                  unsigned int other_rank = dict.dof_to_dict_rank(i);
                  if (other_rank == my_rank)
                    {
                      owning_ranks[index] =
                        dict.actually_owning_ranks[i - dict.local_range.first];
                      if (track_index_requests)
                        append_index_origin(i, owner_index, my_rank);
                    }
                  else if (targets.empty() || targets.back() != other_rank)
                    targets.push_back(other_rank);
                  index++;
                }
            }


            for (auto i : targets)
              {
                recv_indices[i]                    = {};
                indices_to_look_up_by_dict_rank[i] = {};
              }

            // 3) collect indices for each process
            {
              unsigned int index = 0;
              for (auto i : indices_to_look_up)
                {
                  unsigned int other_rank = dict.dof_to_dict_rank(i);
                  if (other_rank != my_rank)
                    {
                      recv_indices[other_rank].push_back(index);
                      indices_to_look_up_by_dict_rank[other_rank].push_back(i);
                    }
                  index++;
                }
            }

            Assert(targets.size() == recv_indices.size() &&
                     targets.size() == indices_to_look_up_by_dict_rank.size(),
                   ExcMessage("Size does not match!"));

            return targets;
          }

          /**
           * Utilities::MPI::ConsensusAlgorithms::Process::create_request().
           * 的实施
           *
           */
          virtual void
          create_request(const unsigned int other_rank,
                         std::vector<std::pair<types::global_dof_index,
                                               types::global_dof_index>>
                           &send_buffer) override
          {
            // create index set and compress data to be sent
            auto &   indices_i = indices_to_look_up_by_dict_rank[other_rank];
            IndexSet is(dict.size);
            is.add_indices(indices_i.begin(), indices_i.end());
            is.compress();

            for (auto interval = is.begin_intervals();
                 interval != is.end_intervals();
                 ++interval)
              send_buffer.emplace_back(*interval->begin(),
                                       interval->last() + 1);
          }

          /**
           * 执行
           * Utilities::MPI::ConsensusAlgorithms::Process::prepare_buffer_for_answer().
           * 。
           *
           */
          virtual void
          prepare_buffer_for_answer(
            const unsigned int         other_rank,
            std::vector<unsigned int> &recv_buffer) override
          {
            recv_buffer.resize(recv_indices[other_rank].size());
          }

          /**
           * 执行
           * Utilities::MPI::ConsensusAlgorithms::Process::read_answer(). 。
           *
           */
          virtual void
          read_answer(const unsigned int               other_rank,
                      const std::vector<unsigned int> &recv_buffer) override
          {
            Assert(recv_indices[other_rank].size() == recv_buffer.size(),
                   ExcMessage("Sizes do not match!"));

            for (unsigned int j = 0; j < recv_indices[other_rank].size(); j++)
              owning_ranks[recv_indices[other_rank][j]] = recv_buffer[j];
          }

          /**
           * 通过将共识算法运行期间积累的字典所有者方面的信息送回原始IndexSet中的所有者，来解决请求的来源。这需要一些点对点的通信。
           * @return
           * 从当前等级请求的处理器和相关指数范围的地图
           *
           */
          std::map<unsigned int, IndexSet>
          get_requesters()
          {
            Assert(track_index_requests,
                   ExcMessage("Must enable index range tracking in "
                              "constructor of ConsensusAlgorithmProcess"));

            std::map<unsigned int, dealii::IndexSet> requested_indices;

#ifdef DEAL_II_WITH_MPI

            static CollectiveMutex      mutex;
            CollectiveMutex::ScopedLock lock(mutex, comm);

            const int mpi_tag = Utilities::MPI::internal::Tags::
              consensus_algorithm_payload_get_requesters;

            // reserve enough slots for the requests ahead; depending on
            // whether the owning rank is one of the requesters or not, we
            // might have one less requests to execute, so fill the requests
            // on demand.
            std::vector<MPI_Request> send_requests;
            send_requests.reserve(requesters.size());

            // We use an integer vector for the data exchange. Since we send
            // data associated to intervals with different requesters, we will
            // need to send (a) the MPI rank of the requester, (b) the number
            // of intervals directed to this requester, and (c) a list of
            // intervals, i.e., two integers per interval. The number of items
            // sent in total can be deduced both via the MPI status message at
            // the receiver site as well as be counting the buckets from
            // different requesters.
            std::vector<std::vector<unsigned int>> send_data(requesters.size());
            for (unsigned int i = 0; i < requesters.size(); ++i)
              {
                // special code for our own indices
                if (dict.actually_owning_rank_list[i] == my_rank)
                  {
                    for (const auto &j : requesters[i])
                      {
                        const types::global_dof_index index_offset =
                          dict.get_index_offset(my_rank);
                        IndexSet &my_index_set = requested_indices[j.first];
                        my_index_set.set_size(owned_indices.size());
                        for (const auto &interval : j.second)
                          my_index_set.add_range(index_offset + interval.first,
                                                 index_offset +
                                                   interval.second);
                      }
                  }
                else
                  {
                    for (const auto &j : requesters[i])
                      {
                        send_data[i].push_back(j.first);
                        send_data[i].push_back(j.second.size());
                        for (const auto &interval : j.second)
                          {
                            send_data[i].push_back(interval.first);
                            send_data[i].push_back(interval.second);
                          }
                      }
                    send_requests.push_back(MPI_Request());
                    const int ierr =
                      MPI_Isend(send_data[i].data(),
                                send_data[i].size(),
                                MPI_UNSIGNED,
                                dict.actually_owning_rank_list[i],
                                mpi_tag,
                                comm,
                                &send_requests.back());
                    AssertThrowMPI(ierr);
                  }
              }

            // receive the data
            for (unsigned int c = 0; c < dict.n_dict_procs_in_owned_indices;
                 ++c)
              {
                // wait for an incoming message
                MPI_Status status;
                int ierr = MPI_Probe(MPI_ANY_SOURCE, mpi_tag, comm, &status);
                AssertThrowMPI(ierr);

                // retrieve size of incoming message
                int number_amount;
                ierr = MPI_Get_count(&status, MPI_UNSIGNED, &number_amount);
                AssertThrowMPI(ierr);

                // receive message
                Assert(number_amount % 2 == 0, ExcInternalError());
                std::vector<std::pair<unsigned int, unsigned int>> buffer(
                  number_amount / 2);
                ierr = MPI_Recv(buffer.data(),
                                number_amount,
                                MPI_UNSIGNED,
                                status.MPI_SOURCE,
                                status.MPI_TAG,
                                comm,
                                &status);
                AssertThrowMPI(ierr);

                // unpack the message and translate the dictionary-local
                // indices coming via MPI to the global index range
                const types::global_dof_index index_offset =
                  dict.get_index_offset(status.MPI_SOURCE);
                unsigned int offset = 0;
                while (offset < buffer.size())
                  {
                    AssertIndexRange(offset + buffer[offset].second,
                                     buffer.size());

                    IndexSet my_index_set(owned_indices.size());
                    for (unsigned int i = offset + 1;
                         i < offset + buffer[offset].second + 1;
                         ++i)
                      my_index_set.add_range(index_offset + buffer[i].first,
                                             index_offset + buffer[i].second);

                    // the underlying index set is able to merge ranges coming
                    // from different ranks due to the partitioning in the
                    // dictionary
                    IndexSet &index_set =
                      requested_indices[buffer[offset].first];
                    if (index_set.size() == 0)
                      index_set.set_size(owned_indices.size());
                    index_set.add_indices(my_index_set);

                    offset += buffer[offset].second + 1;
                  }
                AssertDimension(offset, buffer.size());
              }

            if (send_requests.size() > 0)
              {
                const auto ierr = MPI_Waitall(send_requests.size(),
                                              send_requests.data(),
                                              MPI_STATUSES_IGNORE);
                AssertThrowMPI(ierr);
              }


#  ifdef DEBUG
            for (const auto &it : requested_indices)
              {
                IndexSet copy_set = it.second;
                copy_set.subtract_set(owned_indices);
                Assert(copy_set.n_elements() == 0,
                       ExcInternalError(
                         "The indices requested from the current "
                         "MPI rank should be locally owned here!"));
              }
#  endif

#endif // DEAL_II_WITH_MPI

            return requested_indices;
          }

        private:
          /**
           * 在 "requesters
           * "字段中存储索引请求。我们首先找出被请求的索引的所有者（使用`owner_index`中的猜测，因为我们通常可能在同一等级上连续查找几次，这避免了
           * Dictionary::get_owning_rank_index(). 中的二进制搜索
           * 一旦我们知道所有者的等级，我们就用请求的等级的向量条目。在这里，我们利用了请求被逐级处理的事实，所以我们可以简单地在向量的末尾看看是否已经有一些数据被存储。最后，我们建立范围，再次利用索引列表被排序的特点，因此我们只需要在最后追加。
           *
           */
          void
          append_index_origin(const types::global_dof_index index,
                              unsigned int &                owner_index,
                              const unsigned int            rank_of_request)
          {
            // remember who requested which index. We want to use an
            // std::vector with simple addressing, via a good guess from the
            // preceding index, rather than std::map, because this is an inner
            // loop and it avoids the map lookup in every iteration
            const unsigned int rank_of_owner =
              dict.actually_owning_ranks[index - dict.local_range.first];
            owner_index =
              dict.get_owning_rank_index(rank_of_owner, owner_index);
            if (requesters[owner_index].empty() ||
                requesters[owner_index].back().first != rank_of_request)
              requesters[owner_index].emplace_back(
                rank_of_request,
                std::vector<std::pair<unsigned int, unsigned int>>());
            if (requesters[owner_index].back().second.empty() ||
                requesters[owner_index].back().second.back().second !=
                  index - dict.local_range.first)
              requesters[owner_index].back().second.emplace_back(
                index - dict.local_range.first,
                index - dict.local_range.first + 1);
            else
              ++requesters[owner_index].back().second.back().second;
          }
        };

      } // namespace ComputeIndexOwner
    }   // namespace internal
  }     // namespace MPI
} // namespace Utilities

DEAL_II_NAMESPACE_CLOSE

#endif


