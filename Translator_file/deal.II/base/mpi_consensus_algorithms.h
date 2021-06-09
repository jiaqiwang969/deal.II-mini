//include/deal.II-translator/base/mpi_consensus_algorithms_0.txt
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

#ifndef dealii_mpi_consensus_algorithm_h
#define dealii_mpi_consensus_algorithm_h

#include <deal.II/base/config.h>

#include <deal.II/base/mpi.h>
#include <deal.II/base/mpi.templates.h>

DEAL_II_NAMESPACE_OPEN


namespace Utilities
{
  namespace MPI
  {
    /**
     * 一个为动态稀疏通信模式设计的共识算法命名空间。
     * @ingroup MPI
     *
     */
    namespace ConsensusAlgorithms
    {
      /**
       * 一个能够使用接口类的接口。实现的主要功能是返回此进程想要的数据的进程等级列表，并处理由ConsensusAlgorithm类发送/接收的消息的可选有效载荷。
       * 有两种消息。
       *
       *
       *
       *
       * - 发送/请求消息。一个由数据请求组成的消息，应该由另一个进程来回答。该消息被接收等级视为请求消息。
       *
       *
       *
       *
       *
       * - recv消息。对发送/请求信息的回答。              @tparam  T1 要发送的向量元素的类型  @tparam  T2 要接收的向量元素的类型
       * @note
       * 由于消息的有效载荷是可选的，用户必须自己处理缓冲区。ConsensusAlgorithm类1）只提供对空向量（大小为0）的引用，要发送的数据可以被插入或读出，2）盲目地交流这些向量。
       *
       */
      template <typename T1, typename T2>
      class Process
      {
      public:
        /**
         * 销毁器。
         *
         */
        virtual ~Process() = default;

        /**
         * @return  这个进程要发送请求的行列的向量。
         * @note
         * 这是唯一必须实现的方法，因为信息的有效载荷是可选的。
         *
         */
        virtual std::vector<unsigned int>
        compute_targets() = 0;

        /**
         * 在对进程的请求中加入指定等级的有效载荷。
         * @param[in]  进程的其他等级  @param[out]  send_buffer
         * 要发送的请求的一部分数据（可选）。
         * @note
         * 缓冲区是空的。在使用它之前，你必须设置它的大小。
         *
         */
        virtual void
        create_request(const unsigned int other_rank,
                       std::vector<T1> &  send_buffer);

        /**
         * 准备好缓冲区，在这个缓冲区中保存着对指定等级的进程的请求的回答的有效载荷。最明显的任务是调整缓冲区的大小，因为当函数被调用时它是空的。
         * @param[in]  进程的其他等级  @param[out]  recv_buffer
         * 要发送请求部分的数据（可选）。
         *
         */
        virtual void
        prepare_buffer_for_answer(const unsigned int other_rank,
                                  std::vector<T2> &  recv_buffer);

        /**
         * 准备一个缓冲区，在这个缓冲区中保存对具有指定等级的进程的请求的回答的有效载荷。
         * @param[in]  进程的其他等级  @param[in]  buffer_recv
         * 收到的有效载荷（可选）  @param[out]  request_buffer
         * 将作为请求的一部分发送的有效载荷（可选）。
         * @note
         * request_buffer是空的。在使用它之前，你必须设置它的大小。
         *
         */
        virtual void
        answer_request(const unsigned int     other_rank,
                       const std::vector<T1> &buffer_recv,
                       std::vector<T2> &      request_buffer);

        /**
         * 将请求的答案的有效载荷处理给具有指定等级的进程。
         * @param[in]  进程的其他等级  @param[in]  recv_buffer
         * 将发送请求的一部分数据（可选）。
         *
         */
        virtual void
        read_answer(const unsigned int     other_rank,
                    const std::vector<T2> &recv_buffer);
      };



      /**
       * 一个基类，用于实现提出通信模式的任务，以动态稀疏的方式从其他进程中检索数据的算法。在计算机科学中，这通常被称为<a
       * href="https://en.wikipedia.org/wiki/Consensus_algorithm">consensus
       * problem</a>。            动态稀疏在这里意味着。
       *
       *
       *
       *
       *
       * - 当这个函数被调用时，其他进程还不知道它们必须回答请求。
       *
       *
       *
       * - 每个进程只需要与MPI通信器的一小部分进程进行通信。            当然，用户必须提供。
       *
       *
       *
       *
       * - 一个交流者。
       *
       * - 对于每个等级，这个进程应该与之通信的进程的等级列表。
       *
       *
       *
       *
       *
       * - 打包/解包要发送/接收的数据的功能。            这个基类只介绍了实现这些目标的基本接口，而派生类实现了不同的算法来实际计算这种通信模式。      本段上面列表中的最后两个特征是在派生类中实现的  ConsensusAlgorithm::Process.   @tparam  T1 要发送的向量元素的类型。        @tparam  T2 要接收的向量元素的类型。
       *
       */
      template <typename T1, typename T2>
      class Interface
      {
      public:
        Interface(Process<T1, T2> &process, const MPI_Comm &comm);

        /**
         * 解构器。
         *
         */
        virtual ~Interface() = default;

        /**
         * 运行共识算法并返回请求的进程。
         *
         */
        virtual std::vector<unsigned int>
        run() = 0;

      protected:
        /**
         * 对用户提供的进程的引用。
         *
         */
        Process<T1, T2> &process;

        /**
         * MPI通信器。
         *
         */
        const MPI_Comm &comm;

        /**
         * 如果工作支持MPI，则为缓存。
         *
         */
        const bool job_supports_mpi;

        /**
         * 这个进程的等级。
         *
         */
        const unsigned int my_rank;

        /**
         * 通信器中的进程数。
         *
         */
        const unsigned int n_procs;
      };


      /**
       * 该类实现了 ConsensusAlgorithms::Interface
       * 基类的具体算法，仅使用点对点通信和单一的IBarrier。
       * @note  该类紧跟  @cite hoefler2010scalable
       * 。由于那里显示的算法没有考虑有效载荷，这里对算法进行了修改，同步发送（Issend）被等效的Isend/Irecv所取代，其中Irecv接收请求（带有效载荷）的答案。
       * @tparam  T1 要发送的向量元素的类型。        @tparam  T2
       * 要接收的向量元素的类型。
       *
       */
      template <typename T1, typename T2>
      class NBX : public Interface<T1, T2>
      {
      public:
        /**
         * 构造函数。                  @param  process
         * 共识算法期间要运行的进程。          @param  comm
         * MPI通信器
         *
         */
        NBX(Process<T1, T2> &process, const MPI_Comm &comm);

        /**
         * 解构器。
         *
         */
        virtual ~NBX() = default;

        /**
         * @copydoc   Interface::run() .
         *
         */
        virtual std::vector<unsigned int>
        run() override;

      private:
#ifdef DEAL_II_WITH_MPI
        /**
         * 此进程要发送请求的进程列表。
         *
         */
        std::vector<unsigned int> targets;

        /**
         * 用于发送请求的缓冲区。
         *
         */
        std::vector<std::vector<T1>> send_buffers;

        /**
         * 用于发送请求的请求。
         *
         */
        std::vector<MPI_Request> send_requests;

        /**
         * 用于接收请求的答复的缓冲区。
         *
         */
        std::vector<std::vector<T2>> recv_buffers;


        /**
         * 用于接收请求的答复的请求。
         *
         */
        std::vector<MPI_Request> recv_requests;

        /**
         * 用于发送请求的答案的缓冲区。
         *
         */
        std::vector<std::unique_ptr<std::vector<T2>>> request_buffers;

        /**
         * 用于发送对请求的回答的请求。
         *
         */
        std::vector<std::unique_ptr<MPI_Request>> request_requests;

        // request for barrier
        MPI_Request barrier_request;
#endif

        /**
         * 向本进程发出请求的进程列表。
         *
         */
        std::set<unsigned int> requesting_processes;

        /**
         * 检查这个等级是否已经收到了所有的请求答案。
         *
         */
        bool
        check_own_state();

        /**
         * 向所有其他等级发出信号，表明这个等级已经通过进入IBarrier收到了所有的请求答案。
         *
         */
        void
        signal_finish();

        /**
         * 检查所有等级是否已经收到所有的请求答案，即所有等级都已到达
         * "障碍"。
         *
         */
        bool
        check_global_state();

        /**
         * 已收到另一个等级的请求信息：处理请求并发送答复。
         *
         */
        void
        answer_requests();

        /**
         * 开始通过ISend发送所有的请求，并对收到的回答信息发布IRecvs。
         *
         */
        void
        start_communication();

        /**
         * 在所有等级都收到所有答案后，可以释放MPI数据结构，并对收到的答案进行处理。
         *
         */
        void
        clean_up_and_end_communication();
      };

      /**
       * 这个类实现了 ConsensusAlgorithms::Interface
       * 基类的一个具体算法，使用了两步方法。
       * 在第一步中，源等级被确定，在第二步中，进行静态稀疏数据交换。
       * @note
       * 与NBX不同，该类将同一任务分成两个不同的步骤。在第一步中，确定所有要向该进程发送请求的进程。在第二步中，进行数据交换。然而，由于
       *
       * - 在第二步中
       *
       * - 现在很清楚有多少个请求要被回答，也就是说，当这个进程可以停止等待请求时，就不需要IBarrier了。
       * @note  函数
       * Utilities::MPI::compute_point_to_point_communication_pattern()
       * 用于确定源进程，它实现了  @cite hoefler2010scalable
       * 的PEX-算法。              @tparam  T1
       * 要发送的向量的元素的类型。        @tparam  T2
       * 要接收的向量元素的类型。
       *
       */
      template <typename T1, typename T2>
      class PEX : public Interface<T1, T2>
      {
      public:
        /**
         * 构造函数。                  @param  process
         * 共识算法期间要运行的进程。          @param  comm
         * MPI通信器
         *
         */
        PEX(Process<T1, T2> &process, const MPI_Comm &comm);

        /**
         * 解构器。
         *
         */
        virtual ~PEX() = default;

        /**
         * @copydoc   Interface::run() .
         *
         */
        virtual std::vector<unsigned int>
        run() override;

      private:
#ifdef DEAL_II_WITH_MPI
        /**
         * 此进程要发送请求的进程等级列表。
         *
         */
        std::vector<unsigned int> targets;

        /**
         * 希望向该进程发送请求的进程的等级列表。
         *
         */
        std::vector<unsigned int> sources;

        // data structures to send and receive requests

        /**
         * 用于发送请求的缓冲区。
         *
         */
        std::vector<std::vector<T1>> send_buffers;

        /**
         * 用于接收请求的答案的缓冲区。
         *
         */
        std::vector<std::vector<T2>> recv_buffers;

        /**
         * 用于发送请求和接收请求答案的请求。
         *
         */
        std::vector<MPI_Request> send_and_recv_buffers;

        /**
         * 用于发送请求的答案的缓冲区。
         *
         */
        std::vector<std::vector<T2>> requests_buffers;

        /**
         * 用于发送请求的答案的请求。
         *
         */
        std::vector<MPI_Request> requests_answers;
#endif
        /**
         * 向本进程发出请求的进程列表。
         *
         */
        std::set<unsigned int> requesting_processes;

        /**
         * 已收到来自另一等级的第1个请求信息：处理该请求并发送一个回答。
         *
         */
        void
        answer_requests(int index);

        /**
         * 开始通过ISend发送所有请求，并对收到的应答信息发布IRecvs。
         *
         */
        unsigned int
        start_communication();

        /**
         * 在所有的答案被交换后，可以释放MPI数据结构，并处理收到的答案。
         *
         */
        void
        clean_up_and_end_communication();
      };

      /**
       * 上述类的串行回退，以允许独立于是否使用MPI进行编程。
       *
       */
      template <typename T1, typename T2>
      class Serial : public Interface<T1, T2>
      {
      public:
        /**
         * 构造函数。                  @param  进程
         * 在共识算法过程中运行的进程。          @param  comm
         * MPI通信器（忽略）。
         *
         */
        Serial(Process<T1, T2> &process, const MPI_Comm &comm);

        /**
         * @copydoc   Interface::run() .
         *
         */
        virtual std::vector<unsigned int>
        run() override;
      };

      /**
       * 一个类，根据MPI通信器中的进程数量，将其任务委托给其他
       * ConsensusAlgorithms::Interface
       * 实现。对于少量的进程它使用PEX，对于大量的进程使用NBX。阈值取决于程序是以调试模式还是发布模式编译的。
       * @tparam  T1 要发送的向量元素的类型。        @tparam  T2
       * 要接收的向量元素的类型。
       *
       */
      template <typename T1, typename T2>
      class Selector : public Interface<T1, T2>
      {
      public:
        /**
         * 构造函数。                  @param  process
         * 共识算法期间要运行的进程。          @param  comm
         * MPI通信器。
         *
         */
        Selector(Process<T1, T2> &process, const MPI_Comm &comm);

        /**
         * 解构器。
         *
         */
        virtual ~Selector() = default;

        /**
         * @copydoc   Interface::run() .
         * @note  函数调用被委托给另一个
         * ConsensusAlgorithms::Interface 实现。
         *
         */
        virtual std::vector<unsigned int>
        run() override;

      private:
        // Pointer to the actual ConsensusAlgorithms::Interface implementation.
        std::shared_ptr<Interface<T1, T2>> consensus_algo;
      };

      /**
       * 该类使用用户提供的函数包装器实现了
       * Utilities::MPI::ConsensusAlgorithms::Process, 。
       * 这个类的优点是，用户不必编写自己的实现，而是可以直接注册lambda函数。
       *
       */
      template <typename T1, typename T2>
      class AnonymousProcess : public Process<T1, T2>
      {
      public:
        /**
         * 注册为实现Process的接口而应该被调用的函数。
         * @param  function_compute_targets在`compute_targets'期间调用。
         * @param  function_create_request在`create_request`时被调用。
         * @param  在 "answer_request "中调用 function_answer_request。
         * @param  在 "prepare_buffer_for_answer "中调用
         * function_prepare_buffer_for_answer。          @param  在
         * "read_answer "中调用 function_read_answer。
         *
         */
        AnonymousProcess(
          const std::function<std::vector<unsigned int>()>
            &function_compute_targets,
          const std::function<void(const unsigned int, std::vector<T1> &)>
            &function_create_request =
              [](const unsigned int, std::vector<T1> &) {},
          const std::function<void(const unsigned int,
                                   const std::vector<T1> &,
                                   std::vector<T2> &)>
            &function_answer_request = [](const unsigned int,
                                          const std::vector<T1> &,
                                          std::vector<T2> &) {},
          const std::function<void(const unsigned int, std::vector<T2> &)>
            &function_prepare_buffer_for_answer =
              [](const unsigned int, std::vector<T2> &) {},
          const std::function<void(const unsigned int, const std::vector<T2> &)>
            &function_read_answer =
              [](const unsigned int, const std::vector<T2> &) {});

        /**
         * @copydoc   Process::compute_targets() .
         *
         */
        std::vector<unsigned int>
        compute_targets() override;

        /**
         * @copydoc   Process::create_request() .
         *
         */
        void
        create_request(const unsigned int other_rank,
                       std::vector<T1> &  send_buffer) override;

        /**
         * @copydoc   Process::answer_request()   Process::answer_request()
         * 。
         *
         */
        void
        answer_request(const unsigned int     other_rank,
                       const std::vector<T1> &buffer_recv,
                       std::vector<T2> &      request_buffer) override;

        /**
         * @copydoc   Process::prepare_buffer_for_answer()
         *
         */
        void
        prepare_buffer_for_answer(const unsigned int other_rank,
                                  std::vector<T2> &  recv_buffer) override;

        /**
         * @copydoc   Process::read_answer()
         *
         */
        void
        read_answer(const unsigned int     other_rank,
                    const std::vector<T2> &recv_buffer) override;

      private:
        const std::function<std::vector<unsigned int>()>
          function_compute_targets;
        const std::function<void(const int, std::vector<T1> &)>
          function_create_request;
        const std::function<
          void(const unsigned int, const std::vector<T1> &, std::vector<T2> &)>
          function_answer_request;
        const std::function<void(const int, std::vector<T2> &)>
          function_prepare_buffer_for_answer;
        const std::function<void(const int, const std::vector<T2> &)>
          function_read_answer;
      };



      template <typename T1, typename T2>
      AnonymousProcess<T1, T2>::AnonymousProcess(
        const std::function<std::vector<unsigned int>()>
          &function_compute_targets,
        const std::function<void(const unsigned int, std::vector<T1> &)>
          &                                           function_create_request,
        const std::function<void(const unsigned int,
                                 const std::vector<T1> &,
                                 std::vector<T2> &)> &function_answer_request,
        const std::function<void(const unsigned int, std::vector<T2> &)>
          &function_prepare_buffer_for_answer,
        const std::function<void(const unsigned int, const std::vector<T2> &)>
          &function_read_answer)
        : function_compute_targets(function_compute_targets)
        , function_create_request(function_create_request)
        , function_answer_request(function_answer_request)
        , function_prepare_buffer_for_answer(function_prepare_buffer_for_answer)
        , function_read_answer(function_read_answer)
      {}



      template <typename T1, typename T2>
      std::vector<unsigned int>
      AnonymousProcess<T1, T2>::compute_targets()
      {
        return function_compute_targets();
      }



      template <typename T1, typename T2>
      void
      AnonymousProcess<T1, T2>::create_request(const unsigned int other_rank,
                                               std::vector<T1> &  send_buffer)
      {
        function_create_request(other_rank, send_buffer);
      }



      template <typename T1, typename T2>
      void
      AnonymousProcess<T1, T2>::answer_request(
        const unsigned int     other_rank,
        const std::vector<T1> &buffer_recv,
        std::vector<T2> &      request_buffer)
      {
        function_answer_request(other_rank, buffer_recv, request_buffer);
      }



      template <typename T1, typename T2>
      void
      AnonymousProcess<T1, T2>::prepare_buffer_for_answer(
        const unsigned int other_rank,
        std::vector<T2> &  recv_buffer)
      {
        function_prepare_buffer_for_answer(other_rank, recv_buffer);
      }



      template <typename T1, typename T2>
      void
      AnonymousProcess<T1, T2>::read_answer(const unsigned int     other_rank,
                                            const std::vector<T2> &recv_buffer)
      {
        function_read_answer(other_rank, recv_buffer);
      }



    } // namespace ConsensusAlgorithms
  }   // end of namespace MPI
} // end of namespace Utilities


DEAL_II_NAMESPACE_CLOSE

#endif


