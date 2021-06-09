//include/deal.II-translator/base/mpi_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2021 by the deal.II authors
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

#ifndef dealii_mpi_h
#define dealii_mpi_h

#include <deal.II/base/config.h>

#include <deal.II/base/array_view.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/mpi_tags.h>
#include <deal.II/base/numbers.h>

#include <boost/signals2.hpp>

#include <map>
#include <numeric>
#include <set>
#include <vector>

#if !defined(DEAL_II_WITH_MPI) && !defined(DEAL_II_WITH_PETSC)
// without MPI, we would still like to use
// some constructs with MPI data
// types. Therefore, create some dummies
using MPI_Comm     = int;
using MPI_Request  = int;
using MPI_Datatype = int;
using MPI_Op       = int;
#  ifndef MPI_COMM_WORLD
#    define MPI_COMM_WORLD 0
#  endif
#  ifndef MPI_COMM_SELF
#    define MPI_COMM_SELF 0
#  endif
#  ifndef MPI_REQUEST_NULL
#    define MPI_REQUEST_NULL 0
#  endif
#  ifndef MPI_MIN
#    define MPI_MIN 0
#  endif
#  ifndef MPI_MAX
#    define MPI_MAX 0
#  endif
#  ifndef MPI_SUM
#    define MPI_SUM 0
#  endif
#  ifndef MPI_LOR
#    define MPI_LOR 0
#  endif
#endif



/**
 * 帮助性宏，用于从一些MPI_*的指针参数中移除const。
 * 函数的指针参数中删除const。
 * 这是需要的，因为像MPI_Allgather()这样的函数的输入参数在OpenMPI
 * 1.6.5中没有标记为const。如果使用MPI
 * 3或更新的版本，这个宏是一个NOOP，而我们在其他情况下做以下工作。
 * 1.从 @p expr 的类型中移除 2.从结果类型中移除const
 * 3.添加到结果类型中 4.将给定的表达式 @p expr
 * const_cast到这个新类型中。
 *
 *
 */
#ifdef DEAL_II_WITH_MPI
#  if DEAL_II_MPI_VERSION_GTE(3, 0)

#    define DEAL_II_MPI_CONST_CAST(expr) (expr)

#  else

#    include <type_traits>

#    define DEAL_II_MPI_CONST_CAST(expr)     \
      const_cast<typename std::remove_const< \
        typename std::remove_pointer<decltype(expr)>::type>::type *>(expr)

#  endif
#endif



DEAL_II_NAMESPACE_OPEN


// Forward type declarations to allow MPI sums over tensorial types
#ifndef DOXYGEN
template <int rank, int dim, typename Number>
class Tensor;
template <int rank, int dim, typename Number>
class SymmetricTensor;
template <typename Number>
class SparseMatrix;
class IndexSet;
#endif

namespace Utilities
{
  /**
   * 给出元素总数 @p total_size, ，为整个 @p n_partitions.
   * 的元素创建一个均匀分布的1:1分区，本地大小将等于 @p
   * total_size
   * 除以分区的数量，再加上余下的部分在第一个进程中分配。每个进程将存储一个连续的索引子集，进程p+1上的索引集从比进程p上存储的最后一个索引大的索引开始。例如，一个
   * @p total_size 为11的3个进程将产生索引集{ [0,4), [4,8), [8,11)]
   * }，这个函数将返回 @p my_partition_id  的索引集。
   *
   */
  IndexSet
  create_evenly_distributed_partitioning(const unsigned int my_partition_id,
                                         const unsigned int n_partitions,
                                         const IndexSet::size_type total_size);

  /**
   * 一个命名空间，用于抽象使用消息传递接口（MPI）的某些操作，或者在deal.II被配置为完全不使用MPI的情况下提供后备操作的实用函数。
   * @ingroup utilities
   *
   */
  namespace MPI
  {
    /**
     * 返回给定的 @ref GlossMPICommunicator "communicator "
     * 对象中存在的MPI进程的数量。如果这是一个顺序作业（即，程序根本没有使用MPI，或者使用了MPI但只启动了一个MPI进程），那么通信器必然只涉及一个进程，函数返回1。
     *
     */
    unsigned int
    n_mpi_processes(const MPI_Comm &mpi_communicator);

    /**
     *
     */
    unsigned int
    this_mpi_process(const MPI_Comm &mpi_communicator);

    /**
     * 返回一个行列向量（在 @p comm_large) 指定的进程子集的
     * @p comm_small. 内）。
     *
     */
    const std::vector<unsigned int>
    mpi_processes_within_communicator(const MPI_Comm &comm_large,
                                      const MPI_Comm &comm_small);

    /**
     * 考虑一个非结构化的通信模式，MPI宇宙中的每个进程都想向其他进程的一个子集发送一些数据。要做到这一点，其他处理器需要知道从谁那里期待消息。这个函数可以计算这个信息。          @param  mpi_comm 一个 @ref GlossMPICommunicator  "通信器"
     * ，描述要相互通信的处理器。          @param  destinations
     * 当前进程想要发送信息的处理器列表。这个列表不需要以任何方式进行排序。如果它包含重复的条目，那就意味着有多条信息是要发给某个目的地的。
     * @return
     * 已表示要向当前处理器发送东西的处理器的列表。由此产生的列表没有被排序。
     * 如果处理器在其目的地列表中多次输入同一个目的地，它可能包含重复的条目。
     *
     */
    std::vector<unsigned int>
    compute_point_to_point_communication_pattern(
      const MPI_Comm &                 mpi_comm,
      const std::vector<unsigned int> &destinations);

    /**
     * compute_point_to_point_communication_pattern()的简化版本（为了提高效率），它只计算MPI宇宙中期待通信的进程数。          @param  mpi_comm 一个 @ref GlossMPICommunicator "通信器"
     * ，描述要相互通信的处理器。          @param  destinations
     * 当前进程想要发送信息的处理器的列表。这个列表不需要以任何方式进行排序。如果它包含重复的条目，那就意味着有多条信息是要发给某个目的地的。
     * @return 想要向当前处理器发送东西的处理器的数量。
     *
     */
    unsigned int
    compute_n_point_to_point_communications(
      const MPI_Comm &                 mpi_comm,
      const std::vector<unsigned int> &destinations);

    /**
     * 给定一个 @ref GlossMPICommunicator "通信器"
     * ，生成一个新的通信器，该通信器包含相同的处理器集合，但有一个不同的、唯一的标识。
     * 这个功能可以用来确保不同的对象，如分布式矩阵，都有唯一的通信器，它们可以在上面进行交互而不互相干扰。
     * 当不再需要时，这里创建的通信器需要用free_communicator()来销毁。
     * 这个函数等同于调用  <code>MPI_Comm_dup(mpi_communicator,
     * &return_value);</code>  。
     *
     */
    MPI_Comm
    duplicate_communicator(const MPI_Comm &mpi_communicator);

    /**
     * 释放给定的 @ref GlossMPICommunicator "通信器"
     * @p mpi_communicator
     * ，该通信器是使用diplicate_communicator()复制的。
     * 参数是通过引用传递的，并且将被无效化并设置为MPI空手柄。这个函数等同于调用
     * <code>MPI_Comm_free(&mpi_communicator);</code>  。
     *
     */
    void
    free_communicator(MPI_Comm &mpi_communicator);

    /**
     * 帮助类，用于自动复制和释放MPI  @ref GlossMPICommunicator  "communicator"
     * 。
     * 这个类使用diplicate_communicator()复制构造函数中给出的通信器，并在这个对象被销毁时通过调用free_communicator()自动释放它。你可以使用operator*来访问被包裹的通信器。
     * 这个类的存在是为了轻松地允许复制通信器，而不必担心使用后何时以及如何释放它。
     *
     */
    class DuplicatedCommunicator
    {
    public:
      /**
       * 创建一个给定的 @p communicator. 的复制品。
       *
       */
      explicit DuplicatedCommunicator(const MPI_Comm &communicator)
        : comm(duplicate_communicator(communicator))
      {}

      /**
       * 不允许制作副本。
       *
       */
      DuplicatedCommunicator(const DuplicatedCommunicator &) = delete;

      /**
       * 解构器会自动释放通讯器。
       *
       */
      ~DuplicatedCommunicator()
      {
        free_communicator(comm);
      }

      /**
       * 访问存储的通信器。
       *
       */
      const MPI_Comm &operator*() const
      {
        return comm;
      }


      /**
       * 不允许对这个类进行赋值。
       *
       */
      DuplicatedCommunicator &
      operator=(const DuplicatedCommunicator &) = delete;

    private:
      /**
       * 当然的通讯器。
       *
       */
      MPI_Comm comm;
    };

    /**
     * 这个类代表一个mutex，在使用MPI的并行计算中为一组处理器守卫一个关键部分。
     * lock()命令会等待，直到通信器中的所有MPI等级都使用unlock()释放了之前的锁。
     * 一个典型的用法是使用锁防护来守护一个关键部分。
     * @code
     * {
     * static CollectiveMutex      mutex;
     * CollectiveMutex::ScopedLock lock(mutex, comm);
     * // [ critical code to be guarded]
     * }
     * @endcode
     * 这里，关键代码将在所有处理器上完成，然后才能再次获得mutex（例如通过第二次执行上面的块。关键代码块通常涉及MPI通信，如果没有锁，会产生不正确的结果。例如，如果代码包含具有MPI_ANY_SOURCE的非阻塞接收，数据包在迭代之间可能会被混淆。
     * 请注意，在调用同一关键区域之间，mutex需要是同一个实例。虽然不是必须的，但这可以通过使实例静态化来实现（就像上面的例子）。该变量也可以是一个全局变量，或执行函数所属对象的成员变量。
     *
     */
    class CollectiveMutex
    {
    public:
      /**
       * 这个辅助类为CollectiveMutex提供了一个范围内的锁。
       * 详见CollectiveMutex的类文档。
       *
       */
      class ScopedLock
      {
      public:
        /**
         * 构造函数。阻塞直到它能获得锁。
         *
         */
        explicit ScopedLock(CollectiveMutex &mutex, const MPI_Comm &comm)
          : mutex(mutex)
          , comm(comm)
        {
          mutex.lock(comm);
        }

        /**
         * 销毁器。释放锁。
         *
         */
        ~ScopedLock()
        {
          mutex.unlock(comm);
        }

      private:
        /**
         * 对mutex的引用。
         *
         */
        CollectiveMutex &mutex;
        /**
         * 通信器。
         *
         */
        const MPI_Comm comm;
      };

      /**
       * 该类的构造函数。
       *
       */
      explicit CollectiveMutex();

      /**
       * 销毁mutex。假设当前没有持有锁。
       *
       */
      ~CollectiveMutex();

      /**
       * 获取mutex，如果有必要，等待我们可以这样做。
       * 这是一个集体调用，需要由通信器中的所有处理器来执行。
       *
       */
      void
      lock(const MPI_Comm &comm);

      /**
       * 释放锁。
       * 这是一个集体调用，需要由通信器中的所有处理器来执行。
       *
       */
      void
      unlock(const MPI_Comm &comm);

    private:
      /**
       * 保持跟踪，如果我们现在有这个锁。
       *
       */
      bool locked;

      /**
       * 追踪非阻塞屏障的请求。
       *
       */
      MPI_Request request;
    };



    /**
     * 如果 @p comm
     * 是一个内部通信器，这个函数返回一个新的通信器 @p
     * newcomm ，其通信组由 @p group
     * 参数定义。该函数只对实际想要创建通信器的进程组进行集合，即在
     * @p group
     * 参数中被命名的进程。如果一个给定进程的多个线程同时执行create_group()操作，用户必须通过提供不同的
     * @p tag 或 @p comm 参数来区分这些操作。
     * 这个函数是在MPI-3.0标准中引入的。如果可用，则使用所提供的MPI实现中的相应函数。
     * 否则，该实现遵循以下出版物中描述的实现。
     * @code{.bib}
     * @inproceedings{dinan2011noncollective,
     * title        = {Noncollective communicator creation in MPI},
     * author       = {Dinan, James and Krishnamoorthy, Sriram and Balaji,
     *                 Pavan and Hammond, Jeff R and Krishnan, Manojkumar and
     *                 Tipparaju, Vinod and Vishnu, Abhinav},
     * booktitle    = {European MPI Users' Group Meeting},
     * pages        = {282--291},
     * year         = {2011},
     * organization = {Springer}
     * }
     * @endcode
     *
     *
     */
#ifdef DEAL_II_WITH_MPI
    int
    create_group(const MPI_Comm & comm,
                 const MPI_Group &group,
                 const int        tag,
                 MPI_Comm *       new_comm);
#endif

    /**
     * 考虑到本地拥有的元素数量 @p locally_owned_size,
     * ，在整个MPI通信器中创建一个1:1的元素分区 @p comm.
     * ，元素的总大小是整个MPI通信器中 @p locally_owned_size
     * 的总和。
     * 每个进程将存储连续的索引子集，进程p+1上的索引集从比进程p上存储的最后一个索引大的一个索引开始。
     *
     */
    std::vector<IndexSet>
    create_ascending_partitioning(const MPI_Comm &          comm,
                                  const IndexSet::size_type locally_owned_size);

    /**
     * 给定元素总数  @p total_size,
     * 在MPI通信器上创建一个均匀分布的1:1的元素分区  @p
     * comm.  使用  @p comm
     * 来确定分区数量和处理器ID，以调用上述  @p
     * create_evenly_distributed_partitioning()  函数。
     *
     */
    IndexSet
    create_evenly_distributed_partitioning(
      const MPI_Comm &          comm,
      const IndexSet::size_type total_size);

#ifdef DEAL_II_WITH_MPI
    /**
     * 计算整个MPI通信器 @p comm
     * 的平均值和标准偏差，提供的数值为一个范围`[begin,end)`。
     * 平均值计算为 $\bar x=\frac 1N \sum x_k$ ，其中 $x_k$
     * 是所有处理器上的`begin'和`end'迭代器所指向的元素（即，每个处理器的`[begin,end]范围指向整体元素数量的一个子集）。标准偏差的计算方法是
     * $\sigma=\sqrt{\frac {1}{N-1} \sum |x_k
     *
     * -\bar x|^2}$  ，这被称为无偏的样本方差。          @tparam
     * Number指定了存储均值的类型。
     * 标准偏差被存储为相应的实数类型。
     * 例如，这允许从整数输入值计算统计数据。
     *
     */
    template <class Iterator, typename Number = long double>
    std::pair<Number, typename numbers::NumberTraits<Number>::real_type>
    mean_and_standard_deviation(const Iterator  begin,
                                const Iterator  end,
                                const MPI_Comm &comm);
#endif

    /**
     * 返回所有处理器上的数值之和  @p t.  这个函数是在 @ref GlossMPICommunicator "通信器 "
     * 中给出的所有处理器上的集体。
     * 如果deal.II没有被配置为使用MPI，这个函数只是返回  @p
     * t.  这个函数对应于  <code>MPI_Allreduce</code>
     * 函数，即所有处理器都收到这个操作的结果。
     * @note
     * 有时，并非所有处理器都需要一个结果，在这种情况下，人们会调用
     * <code>MPI_Reduce</code> 函数而不是 <code>MPI_Allreduce</code>
     * 函数。后者的费用最多是前者的两倍，所以如果你关心性能，可能值得调查一下你的算法是否确实到处需要结果。
     * @note  这个函数只对某些模板参数实现  <code>T</code>,
     * namely <code>float, double, int, unsigned int</code>  。
     *
     */
    template <typename T>
    T
    sum(const T &t, const MPI_Comm &mpi_communicator);

    /**
     * 和前面的函数一样，但是对T类型的数组的元素进行求和。换句话说，结果数组的第i个元素是每个处理器的输入数组的第i个条目之和。T和U必须衰减到相同的类型，例如，它们的区别只是其中一个有const类型限定符，另一个没有。
     * 输入和输出数组可以是相同的。
     *
     */
    template <typename T, typename U>
    void
    sum(const T &values, const MPI_Comm &mpi_communicator, U &sums);

    /**
     * 和前面的函数一样，但是要对ArrayView参数指定的数组元素进行求和。
     * 换句话说，结果数组的第i个元素是每个处理器的输入数组的第i个条目之和。
     * 输入和输出数组可以是相同的。
     *
     */
    template <typename T>
    void
    sum(const ArrayView<const T> &values,
        const MPI_Comm &          mpi_communicator,
        const ArrayView<T> &      sums);

    /**
     * 对一个对称张量的条目进行MPI求和。          @relatesalso
     * SymmetricTensor
     *
     */
    template <int rank, int dim, typename Number>
    SymmetricTensor<rank, dim, Number>
    sum(const SymmetricTensor<rank, dim, Number> &local,
        const MPI_Comm &                          mpi_communicator);

    /**
     * 对一个张量的条目进行MPI求和。          @relatesalso
     * 张量
     *
     */
    template <int rank, int dim, typename Number>
    Tensor<rank, dim, Number>
    sum(const Tensor<rank, dim, Number> &local,
        const MPI_Comm &                 mpi_communicator);

    /**
     * 对稀疏矩阵的条目进行MPI求和。
     * @note   @p local  和  @p global
     * 应该具有相同的稀疏模式，而且对所有MPI进程都应该是相同的。
     * @relatesalso  稀疏矩阵
     *
     */
    template <typename Number>
    void
    sum(const SparseMatrix<Number> &local,
        const MPI_Comm &            mpi_communicator,
        SparseMatrix<Number> &      global);

    /**
     * 返回所有处理器上的最大值  @p t.  这个函数是在 @ref GlossMPICommunicator  "通信器 "
     * 中给出的所有处理器上的集合。
     * 如果deal.II没有被配置为使用MPI，这个函数只是返回 @p
     * t. 的值，这个函数对应于 <code>MPI_Allreduce</code>
     * 函数，即所有处理器都收到这个操作的结果。
     * @note
     * 有时，并非所有处理器都需要一个结果，在这种情况下，人们会调用
     * <code>MPI_Reduce</code> 函数而不是 <code>MPI_Allreduce</code>
     * 函数。后者的费用最多是前者的两倍，所以如果你关心性能，可能值得调查一下你的算法是否确实到处需要结果。
     * @note  这个函数只对某些模板参数实现  <code>T</code>,
     * namely <code>float, double, int, unsigned int</code>  。
     *
     */
    template <typename T>
    T
    max(const T &t, const MPI_Comm &mpi_communicator);

    /**
     * 和前面的函数一样，但是在一个T类型的数组的元素上取最大值。换句话说，结果数组的第i个元素是每个处理器的输入数组的第i个条目的最大值。T和U必须衰减到相同的类型，例如，它们的区别只是其中一个有const类型限定符，另一个没有。
     * 输入和输出向量可以是相同的。
     *
     */
    template <typename T, typename U>
    void
    max(const T &values, const MPI_Comm &mpi_communicator, U &maxima);

    /**
     * 和前面的函数一样，但在ArrayView参数指定的数组元素上取最大值。
     * 换句话说，结果数组的第i个元素是每个处理器的输入数组的第i个条目上的最大值。
     * 输入和输出数组可以是相同的。
     *
     */
    template <typename T>
    void
    max(const ArrayView<const T> &values,
        const MPI_Comm &          mpi_communicator,
        const ArrayView<T> &      maxima);

    /**
     * 返回所有处理器上的最小值  @p t.  这个函数是在 @ref GlossMPICommunicator "通信器 "
     * 中给出的所有处理器上的集合。
     * 如果deal.II没有被配置为使用MPI，这个函数只是返回  @p
     * t.  这个函数对应于  <code>MPI_Allreduce</code>
     * 函数，即所有处理器都收到这个操作的结果。
     * @note
     * 有时，并非所有处理器都需要一个结果，在这种情况下，人们会调用
     * <code>MPI_Reduce</code> 函数而不是 <code>MPI_Allreduce</code>
     * 函数。后者的费用最多是前者的两倍，所以如果你关心性能，可能值得调查一下你的算法是否确实到处需要结果。
     * @note  这个函数只对某些模板参数实现  <code>T</code>,
     * namely <code>float, double, int, unsigned int</code>  。
     *
     */
    template <typename T>
    T
    min(const T &t, const MPI_Comm &mpi_communicator);

    /**
     * 和前面的函数一样，但在一个T类型的数组的元素上取最小值。换句话说，结果数组的第i个元素是每个处理器的输入数组的第i个条目的最小值。T和U必须衰减到相同的类型，例如，它们的区别只是其中一个有const类型限定符，另一个没有。
     * 输入和输出数组可以是相同的。
     *
     */
    template <typename T, typename U>
    void
    min(const T &values, const MPI_Comm &mpi_communicator, U &minima);

    /**
     * 和前面的函数一样，但是在ArrayView参数指定的数组元素中取最小值。
     * 换句话说，结果数组的第i个元素是每个处理器的输入数组的第i个条目上的最小值。
     * 输入和输出数组可以是相同的。
     *
     */
    template <typename T>
    void
    min(const ArrayView<const T> &values,
        const MPI_Comm &          mpi_communicator,
        const ArrayView<T> &      minima);

    /**
     * 在数值为 @p t. 的所有处理器上执行<i>logical or</i>的操作。 <i>logical or</i>操作符`||`如果任一或所有操作数为`真'，则返回布尔值`真'，否则返回`假`。如果提供的值 @p t 在其相关的数据类型`T`中对应于`0`，它将被解释为`false`，否则为`true`。数据类型`T`必须是`integral`类型，即`bool`、`char`、`short`、`int`、`long`，或它们的任何变化。        这个函数是在 @ref GlossMPICommunicator "通信器 "
     * 中给出的所有处理器上的集体。
     * 如果deal.II没有被配置为使用MPI，这个函数只是返回 @p
     * value. 的值，这个函数对应于 <code>MPI_Allreduce</code>
     * 函数，即所有处理器都收到这个操作的结果。
     * @note
     * 有时，并非所有处理器都需要一个结果，在这种情况下，人们会调用
     * <code>MPI_Reduce</code> 函数而不是 <code>MPI_Allreduce</code>
     * 函数。后者的费用最多是前者的两倍，所以如果你关心性能，可能值得调查一下你的算法是否确实到处需要结果。
     *
     */
    template <typename T>
    T
    logical_or(const T &t, const MPI_Comm &mpi_communicator);

    /**
     * 和前面的函数一样，但是对数组中的每个元素执行<i>logical
     * or</i>操作。换句话说，结果数组的第i个元素是对每个处理器的输入数组的第i个条目应用<i>logical
     * or</i>操作的结果。T和U必须衰减到相同的类型，例如，它们只是因其中一个有const类型限定符而另一个没有而不同。
     * 输入和输出数组可以是相同的。
     * @note
     * 根据你的标准库，这个函数可能无法与数据类型`bool`的
     * `std::vector`
     * 的特殊化一起工作。在这种情况下，请使用一个不同的容器或数据类型。
     *
     */
    template <typename T, typename U>
    void
    logical_or(const T &values, const MPI_Comm &mpi_communicator, U &results);

    /**
     * 和前面的函数一样，但是对ArrayView参数指定的数组中的每个元素执行<i>logical
     * or</i>操作。
     * 换句话说，结果数组的第i个元素是对每个处理器的输入数组的第i个条目应用<i>logical
     * or</i>操作的结果。
     * 输入和输出数组可以是相同的。
     *
     */
    template <typename T>
    void
    logical_or(const ArrayView<const T> &values,
               const MPI_Comm &          mpi_communicator,
               const ArrayView<T> &      results);

    /**
     * @note
     * 这个结构没有构造函数，因为MPI要求它是一个POD类型。
     *
     */
    struct MinMaxAvg
    {
      /**
       * 参与调用min_max_avg()的处理器贡献的所有值的总和。
       *
       */
      double sum;

      /**
       * 参与调用min_max_avg()的处理器贡献的所有数值的最小值。
       *
       */
      double min;

      /**
       * 参与调用min_max_avg()的处理器贡献的所有数值的最大值。
       *
       */
      double max;

      /**
       *
       */
      unsigned int min_index;

      /**
       *
       */
      unsigned int max_index;

      /**
       * 参与调用min_max_avg()的处理器所贡献数值的平均值。
       *
       */
      double avg;
    };

    /**
     * 在给定的MPI  @ref GlossMPICommunicator  "communicator"
     * @p mpi_communicator.
     * 的集体操作中，返回总和、平均值、最小值、最大值、最小和最大值的处理器ID，并将返回结果。该结果在所有机器上都可用。
     * @note
     * 有时，并非所有处理器都需要结果，在这种情况下，人们会调用
     * <code>MPI_Reduce</code> 函数而不是 <code>MPI_Allreduce</code>
     * 函数。后者最多只有两倍的费用，所以如果你关心性能，可能值得调查一下你的算法是否确实到处需要结果。
     *
     */
    MinMaxAvg
    min_max_avg(const double my_value, const MPI_Comm &mpi_communicator);

    /**
     * 与上述相同，但在给定的MPI  @ref GlossMPICommunicator  "communicator"
     * @p mpi_communicator
     * 上对向量的每个条目返回总和、平均数、最小值、最大值、最小和最大值的进程ID作为集体操作。
     * @note  该函数执行一次缩减扫频。          @pre
     * 输入向量的大小在所有进程中必须是相同的。
     *
     */
    std::vector<MinMaxAvg>
    min_max_avg(const std::vector<double> &my_value,
                const MPI_Comm &           mpi_communicator);


    /**
     * 和上面一样，但是在给定的MPI  @ref GlossMPICommunicator  "communicator"
     * @p mpi_communicator  上对ArrayView的每个条目返回sum, average,
     * minimum, maximum, process id of minimum and
     * maximum作为一个集体操作。
     * @note  该函数执行一次缩减扫频。          @pre
     * 输入ArrayView的大小在所有进程中必须是相同的，并且输入和输出ArrayVew必须有相同的大小。
     *
     */
    void
    min_max_avg(const ArrayView<const double> &my_values,
                const ArrayView<MinMaxAvg> &   result,
                const MPI_Comm &               mpi_communicator);


    /**
     * 一个用于在程序开始时初始化MPI系统并在程序结束时再次关闭的类。它还允许你控制每个MPI进程中使用的线程数量。
     * 如果deal.II配置了PETSc，PETSc将在开始时通过`PetscInitialize`初始化（该类的构造函数），在结束时通过`PetscFinalize`去初始化（即在该类的析构函数中）。这对SLEPc来说也是如此。
     * 如果deal.II配置了p4est，该库也将在开始时被初始化，并在结束时被解除初始化（通过调用sc_init(),
     * p4est_init(), 和sc_finalize()）。
     * 如果一个程序使用MPI，通常只需在  <code>main()</code>
     * 的开头创建一个这种类型的对象。然后，这个类的构造函数带着给定的参数运行
     * <code>MPI_Init()</code>
     * ，同时初始化上面提到的其他库。在程序结束时，编译器将调用这个对象的析构器，反过来调用
     * <code>MPI_Finalize</code> 来关闭MPI系统。        这个类在
     * step-17 、 step-18 、 step-40 、 step-32
     * 和其他一些地方使用。
     * @note
     * 该类通过`MPI_COMM_WORLD`通信器执行MPI子系统以及上面列出的依赖库的初始化。这意味着你必须在<i>all</i>MPI进程上创建一个MPI_InitFinalize对象，无论你是否打算在特定的处理器上使用deal.II。在大多数使用情况下，人们当然希望使用基本相同的程序在所有MPI进程上工作，因此这不是一个问题。但是如果你计划只在MPI进程的一个子集上运行基于deal.II的工作，使用@
     * ref GlossMPICommunicator "MPI
     * communicator"，它是`MPI_COMM_WORLD`的一个子集（例如，在客户-服务器设置中，只有一个子集的进程负责有限元通信，其余进程做其他事情），那么你仍然需要在程序开始时在所有MPI进程中创建这个对象，因为它在初始化时使用`MPI_COMM_WORLD`。
     *
     */
    class MPI_InitFinalize
    {
    public:
      /**
       * 初始化MPI（如果deal.II被配置为使用它，则初始化PETSc），并将deal.II使用的线程数（通过底层的线程构件库）设置为给定参数。
       * @param[in,out]  argc
       * 对传递给main的'argc'参数的引用。这个参数用于初始化MPI（可能还有PETSc），因为它们从命令行读取参数。
       * @param[in,out]  argv 对传递给main的'argv'参数的引用。
       * @param[in]  max_num_threads
       * 这个MPI进程应该利用的最大线程数。如果这个参数被设置为
       * numbers::invalid_unsigned_int
       * （默认值），那么线程的数量将以如下方式自动确定：在这个MPI进程上运行的线程数量是以你的节点上所有的核心都被使用的方式设置的。换句话说，如果你在每个节点上启动了一个MPI进程，设置这个参数就相当于把它设置为这个MPI进程所运行的节点上的核心数。如果你在每个节点上启动的MPI进程与每个节点上的核数一样多，那么这就相当于将1作为参数传递。另一方面，例如，如果你在每个16核节点上启动4个MPI进程，那么这个选项将为每个节点启动4个工作线程。如果你在一个8核节点上启动3个进程，那么它们将分别启动3、3和2个线程。
       * @note  这个函数用 @p max_num_threads
       * 或者按照上面的讨论，用等于分配给这个MPI进程的核数的线程数调用
       * MultithreadInfo::set_thread_limit() 。然而，
       * MultithreadInfo::set_thread_limit()
       * 反过来也评估了环境变量DEAL_II_NUM_THREADS。最后，工作线程只能在当前MPI进程可以访问的核上创建；一些MPI实现将每个进程可以访问的核数量限制在一个或一个子集上，以确保更好的缓存行为。因此，真正被创建的线程数将是这里传递的参数、环境变量（如果设置了）和线程可访问的核心数的最小值。
       * @note   MultithreadInfo::set_thread_limit()
       * 只有在创建任何线程之前调用它才能发挥作用。因此，对它的调用最安全的地方是在
       * <code>main()</code>  的开头。
       * 因此，这延伸到了当前的类：创建这种类型的对象的最佳位置也是在
       * <code>main()</code>  的顶部或接近顶部。
       *
       */
      MPI_InitFinalize(
        int &              argc,
        char **&           argv,
        const unsigned int max_num_threads = numbers::invalid_unsigned_int);

      /**
       * 解构器。在该类拥有MPI进程的情况下调用<tt>MPI_Finalize()</tt>。
       *
       */
      ~MPI_InitFinalize();

      /**
       * 注册一个MPI_Request的引用，在调用`MPI_Finalize'之前，我们需要对其调用`MPI_Wait'。
       * 当MPI_Finalize被调用时，该对象 @p request
       * 需要存在，这意味着该请求通常是静态分配的。否则，你需要在请求超出范围之前调用unregister_request()。注意，一个请求已经被等待（并因此被重置为MPI_REQUEST_NULL）是可以接受的。
       * 在同一个实例中多次调用这个函数是可以接受的（就像下面的例子中所做的）。
       * 通常情况下，这个函数被CollectiveMutex使用，而不是直接使用，但它也可以像这样直接使用。
       * @code
       * void my_fancy_communication()
       * {
       * static MPI_Request request = MPI_REQUEST_NULL;
       * MPI_InitFinalize::register_request(request);
       * MPI_Wait(&request, MPI_STATUS_IGNORE);
       * // [some algorithm that is not safe to be executed twice in a row.]
       * MPI_IBarrier(comm, &request);
       * }
       * @endcode
       *
       *
       */
      static void
      register_request(MPI_Request &request);

      /**
       * 取消先前使用register_request()添加的请求的注册。
       *
       */
      static void
      unregister_request(MPI_Request &request);

      /**
       * 一个具有 boost::signal 对象的结构，用于注册MPI
       * init或finalize之后的回调运行。
       * 关于信号的文档，见http://www.boost.org/doc/libs/release/libs/signals2
       * 。
       *
       */
      struct Signals
      {
        /**
         * 在我们用 <code>MPI_Init()</code>
         * 初始化MPI上下文后立即触发的信号。
         *
         */
        boost::signals2::signal<void()> at_mpi_init;

        /**
         * 一个在我们用  <code>MPI_Finalize()</code>
         * 关闭MPI上下文之前触发的信号。它可用于在调用
         * <code>MPI_Finalize()</code>
         * 之前取消静态分配的MPI资源，这些资源需要被取消分配。
         *
         */
        boost::signals2::signal<void()> at_mpi_finalize;
      };

      static Signals signals;

    private:
      /**
       * 在最终确定之前对MPI_Wait的请求
       *
       */
      static std::set<MPI_Request *> requests;
    };

    /**
     * 返回(i)deal.II是否已被编译为支持MPI（例如用
     * <code>CXX=mpiCC</code> 编译），如果是，是否(ii)
     * <code>MPI_Init()</code> 已被调用（例如使用
     * Utilities::MPI::MPI_InitFinalize
     * 类）。换句话说，该结果表明当前作业是否在MPI下运行。
     * @note
     * 该函数没有考虑到一个MPI作业是否实际运行在一个以上的处理器上，或者实际上是一个恰好在MPI下运行的单节点作业。
     *
     */
    bool
    job_supports_mpi();

    /**
     * 发起一个某某通信，并在处理器之间交换任意对象（类T应该是可序列化的，使用
     * boost::serialize) 。          @param[in]  comm MPI通信器。
     * @param[in]  objects_to_send
     * 从意在接收数据的进程的等级（无符号int）和要发送的对象的映射（类型`T`必须是可序列化的，这个函数才能正常工作）。如果这个映射包含一个键值等于当前进程等级的条目（即一个向进程发送数据给自己的指令），那么这个数据项就被简单地复制到返回的对象中。
     * @return
     * 从发送数据的进程的等级（无符号int）和收到的对象的映射。
     *
     */
    template <typename T>
    std::map<unsigned int, T>
    some_to_some(const MPI_Comm &                 comm,
                 const std::map<unsigned int, T> &objects_to_send);

    /**
     * 经典MPI_Allgather函数的泛化，它接受任意的数据类型T，只要
     * boost::serialize 接受T作为参数。          @param[in]  comm
     * MPI通信器。      @param[in]  object_to_send
     * 一个要发送给所有其他进程的对象  @return
     * 一个对象的向量，其大小等于MPI通信器中的进程数。每个条目包含从处理器收到的对象，在通信器中具有相应的等级。
     *
     */
    template <typename T>
    std::vector<T>
    all_gather(const MPI_Comm &comm, const T &object_to_send);

    /**
     * 经典MPI_Gather函数的泛化，它接受任意的数据类型T，只要
     * boost::serialize 接受T作为参数。          @param[in]  comm
     * MPI通信器。      @param[in]  object_to_send
     * 一个要发送给根进程的对象  @param[in]  root_process
     * 进程，它接收来自所有进程的对象。默认情况下，等级为0的进程是根进程。
     * @return   @p root_process
     * 接收一个对象的向量，其大小等于MPI通信器中的进程数。每个条目包含从通信器中具有相应等级的处理器接收的对象。所有其他进程收到一个空的向量。
     *
     */
    template <typename T>
    std::vector<T>
    gather(const MPI_Comm &   comm,
           const T &          object_to_send,
           const unsigned int root_process = 0);

    /**
     * 从进程 @p root_process 发送一个对象 @p object_to_send
     * 到所有其他进程。
     * 经典的`MPI_Bcast`函数的泛化，接受任意的数据类型`T`，只要
     * Utilities::pack() （它反过来使用 `boost::serialize`, ，详见
     * Utilities::pack() ）接受`T`作为参数。          @param[in]  comm
     * MPI通信器。      @param[in]  object_to_send
     * 一个要发送给所有进程的对象。      @param[in]
     * root_process
     * 向所有进程发送对象的进程。默认情况下，等级为0的进程是根进程。
     * @return  在根进程上，返回一份  @p object_to_send.
     * 在其他每个进程上，返回一份由  @p root_process.
     * 发送的对象的副本。
     *
     */
    template <typename T>
    T
    broadcast(const MPI_Comm &   comm,
              const T &          object_to_send,
              const unsigned int root_process = 0);

    /**
     * 一个通过用户指定的二进制操作 @p combiner 在 @p
     * root_process. 上结合所有进程的值 @p local_value
     * 的函数，因此，这个函数类似于MPI_Reduce（和
     * Utilities::MPI::min/max()):
     * 然而，一方面由于用户指定的二进制操作，它对内置类型较慢，但另一方面，可以处理一般对象类型，包括存储可变数量数据的对象。
     * 与all_reduce相反，结果将只在单一等级上可用。在所有其他进程中，返回值是未定义的。
     *
     */
    template <typename T>
    T
    reduce(const T &                                     local_value,
           const MPI_Comm &                              comm,
           const std::function<T(const T &, const T &)> &combiner,
           const unsigned int                            root_process = 0);

    /**
     * 一个通过用户指定的二进制操作 @p combiner
     * 将所有进程的值 @p local_value
     * 结合起来并将结果分配给所有进程的函数。因此，这个函数类似于MPI_Allreduce（如果它是由一个全局还原和一个广播步骤实现的），但由于用户指定的二进制操作，也可以处理一般的对象类型，包括那些存储可变数量的数据的对象。
     *
     */
    template <typename T>
    T
    all_reduce(const T &                                     local_value,
               const MPI_Comm &                              comm,
               const std::function<T(const T &, const T &)> &combiner);

    /**
     * 给定一个分割的索引集空间，根据分割的索引集，计算第二个索引集的每个元素的自有MPI进程等级。这个函数的一个自然用法是为每个鬼魂的自由度计算拥有该索引的进程的MPI等级。
     * 人们可能会想："但是我们知道一个鬼魂自由度属于哪个等级，基于它所在单元的子域ID"。但是这种启发式方法对于不同子域ID的幽灵单元之间的接口上的DoF，或者幽灵单元和人工单元之间的接口上的DoF是失败的。此外，这个函数可以实现完全抽象的信息交换，而不需要借助于网格中的邻居。
     * 传给这个函数的第一个参数， @p owned_indices,
     * 必须在所有进程之间唯一地划分一个索引空间。
     * 否则，对这个参数没有任何限制。特别是，没有必要将索引空间划分为连续的子集。此外，对第二个索引集
     * @p indices_to_look_up
     * 没有任何限制，只要其大小与第一个索引集相符。它可以在每个进程上任意独立选择。在第二个索引集也包含本地拥有的索引的情况下，这些索引将被正确处理，并为这些条目返回本进程的等级。
     * @note
     * 这是一个集体操作：给定通信器内的所有进程都必须调用这个函数。由于这个函数不使用MPI_Alltoall或MPI_Allgather，而是使用非阻塞的点对点通信来代替，并且只使用一个非阻塞的屏障，所以它大大减少了内存消耗。这个函数适合于具有>100k
     * MPI行列的大规模模拟。          @param[in]  owned_indices
     * 指数集，包含本进程本地拥有的指数。      @param[in]
     * indices_to_look_up
     * 包含用户对拥有的进程的等级感兴趣的指数的索引集。
     * @param[in]  comm MPI通信器。          @return
     * 包含索引集中每个条目的MPI进程等级的列表  @p
     * indices_to_look_up.  顺序与ElementIterator中的顺序相吻合。
     *
     */
    std::vector<unsigned int>
    compute_index_owner(const IndexSet &owned_indices,
                        const IndexSet &indices_to_look_up,
                        const MPI_Comm &comm);

    /**
     * 计算MPI通信器中所有进程的输入向量 @p vec 的联合 @p
     * comm.  。
     * @note
     * 这是一个集体操作。其结果将在所有进程中可用。
     *
     */
    template <typename T>
    std::vector<T>
    compute_set_union(const std::vector<T> &vec, const MPI_Comm &comm);

    /**
     * 与上述相同，但对 std::set. 而言。
     *
     */
    template <typename T>
    std::set<T>
    compute_set_union(const std::set<T> &set, const MPI_Comm &comm);

#ifndef DOXYGEN
    // declaration for an internal function that lives in mpi.templates.h
    namespace internal
    {
      template <typename T>
      void
      all_reduce(const MPI_Op &            mpi_op,
                 const ArrayView<const T> &values,
                 const MPI_Comm &          mpi_communicator,
                 const ArrayView<T> &      output);
    }

    // Since these depend on N they must live in the header file
    template <typename T, unsigned int N>
    void
    sum(const T (&values)[N], const MPI_Comm &mpi_communicator, T (&sums)[N])
    {
      internal::all_reduce(MPI_SUM,
                           ArrayView<const T>(values, N),
                           mpi_communicator,
                           ArrayView<T>(sums, N));
    }

    template <typename T, unsigned int N>
    void
    max(const T (&values)[N], const MPI_Comm &mpi_communicator, T (&maxima)[N])
    {
      internal::all_reduce(MPI_MAX,
                           ArrayView<const T>(values, N),
                           mpi_communicator,
                           ArrayView<T>(maxima, N));
    }

    template <typename T, unsigned int N>
    void
    min(const T (&values)[N], const MPI_Comm &mpi_communicator, T (&minima)[N])
    {
      internal::all_reduce(MPI_MIN,
                           ArrayView<const T>(values, N),
                           mpi_communicator,
                           ArrayView<T>(minima, N));
    }

    template <typename T, unsigned int N>
    void
    logical_or(const T (&values)[N],
               const MPI_Comm &mpi_communicator,
               T (&results)[N])
    {
      static_assert(std::is_integral<T>::value,
                    "The MPI_LOR operation only allows integral data types.");

      internal::all_reduce(MPI_LOR,
                           ArrayView<const T>(values, N),
                           mpi_communicator,
                           ArrayView<T>(results, N));
    }

    template <typename T>
    std::map<unsigned int, T>
    some_to_some(const MPI_Comm &                 comm,
                 const std::map<unsigned int, T> &objects_to_send)
    {
#  ifndef DEAL_II_WITH_MPI
      (void)comm;
      Assert(objects_to_send.size() < 2,
             ExcMessage("Cannot send to more than one processor."));
      Assert(objects_to_send.find(0) != objects_to_send.end() ||
               objects_to_send.size() == 0,
             ExcMessage("Can only send to myself or to nobody."));
      return objects_to_send;
#  else
      const auto my_proc = this_mpi_process(comm);

      std::map<unsigned int, T> received_objects;

      std::vector<unsigned int> send_to;
      send_to.reserve(objects_to_send.size());
      for (const auto &m : objects_to_send)
        if (m.first == my_proc)
          received_objects[my_proc] = m.second;
        else
          send_to.emplace_back(m.first);

      const unsigned int n_point_point_communications =
        Utilities::MPI::compute_n_point_to_point_communications(comm, send_to);

      // Protect the following communication:
      static CollectiveMutex      mutex;
      CollectiveMutex::ScopedLock lock(mutex, comm);

      // If we have something to send, or we expect something from other
      // processors, we need to visit one of the two scopes below. Otherwise,
      // no other action is required by this mpi process, and we can safely
      // return.
      if (send_to.size() == 0 && n_point_point_communications == 0)
        return received_objects;

      const int mpi_tag =
        internal::Tags::compute_point_to_point_communication_pattern;

      // Sending buffers
      std::vector<std::vector<char>> buffers_to_send(send_to.size());
      std::vector<MPI_Request>       buffer_send_requests(send_to.size());
      {
        unsigned int i = 0;
        for (const auto &rank_obj : objects_to_send)
          if (rank_obj.first != my_proc)
            {
              const auto &rank   = rank_obj.first;
              buffers_to_send[i] = Utilities::pack(rank_obj.second,
                                                    /*allow_compression=*/ false);
              const int ierr     = MPI_Isend(buffers_to_send[i].data(),
                                         buffers_to_send[i].size(),
                                         MPI_CHAR,
                                         rank,
                                         mpi_tag,
                                         comm,
                                         &buffer_send_requests[i]);
              AssertThrowMPI(ierr);
              ++i;
            }
      }

      // Fill the output map
      {
        std::vector<char> buffer;
        // We do this on a first come/first served basis
        for (unsigned int i = 0; i < n_point_point_communications; ++i)
          {
            // Probe what's going on. Take data from the first available sender
            MPI_Status status;
            int        ierr = MPI_Probe(MPI_ANY_SOURCE, mpi_tag, comm, &status);
            AssertThrowMPI(ierr);

            // Length of the message
            int len;
            ierr = MPI_Get_count(&status, MPI_CHAR, &len);
            AssertThrowMPI(ierr);
            buffer.resize(len);

            // Source rank
            const unsigned int rank = status.MPI_SOURCE;

            // Actually receive the message
            ierr = MPI_Recv(buffer.data(),
                            len,
                            MPI_CHAR,
                            status.MPI_SOURCE,
                            status.MPI_TAG,
                            comm,
                            MPI_STATUS_IGNORE);
            AssertThrowMPI(ierr);
            Assert(received_objects.find(rank) == received_objects.end(),
                   ExcInternalError(
                     "I should not receive again from this rank"));
            received_objects[rank] =
              Utilities::unpack<T>(buffer,
                                    /*allow_compression=*/ false);
          }
      }

      // Wait to have sent all objects.
      const int ierr = MPI_Waitall(send_to.size(),
                                   buffer_send_requests.data(),
                                   MPI_STATUSES_IGNORE);
      AssertThrowMPI(ierr);

      return received_objects;
#  endif // deal.II with MPI
    }

    template <typename T>
    std::vector<T>
    all_gather(const MPI_Comm &comm, const T &object)
    {
      if (job_supports_mpi() == false)
        return {object};

#  ifndef DEAL_II_WITH_MPI
      (void)comm;
      std::vector<T> v(1, object);
      return v;
#  else
      const auto n_procs = dealii::Utilities::MPI::n_mpi_processes(comm);

      std::vector<char> buffer = Utilities::pack(object);

      int n_local_data = buffer.size();

      // Vector to store the size of loc_data_array for every process
      std::vector<int> size_all_data(n_procs, 0);

      // Exchanging the size of each buffer
      int ierr = MPI_Allgather(
        &n_local_data, 1, MPI_INT, size_all_data.data(), 1, MPI_INT, comm);
      AssertThrowMPI(ierr);

      // Now computing the displacement, relative to recvbuf,
      // at which to store the incoming buffer
      std::vector<int> rdispls(n_procs);
      rdispls[0] = 0;
      for (unsigned int i = 1; i < n_procs; ++i)
        rdispls[i] = rdispls[i - 1] + size_all_data[i - 1];

      // Step 3: exchange the buffer:
      std::vector<char> received_unrolled_buffer(rdispls.back() +
                                                 size_all_data.back());

      ierr = MPI_Allgatherv(buffer.data(),
                            n_local_data,
                            MPI_CHAR,
                            received_unrolled_buffer.data(),
                            size_all_data.data(),
                            rdispls.data(),
                            MPI_CHAR,
                            comm);
      AssertThrowMPI(ierr);

      std::vector<T> received_objects(n_procs);
      for (unsigned int i = 0; i < n_procs; ++i)
        {
          std::vector<char> local_buffer(received_unrolled_buffer.begin() +
                                           rdispls[i],
                                         received_unrolled_buffer.begin() +
                                           rdispls[i] + size_all_data[i]);
          received_objects[i] = Utilities::unpack<T>(local_buffer);
        }

      return received_objects;
#  endif
    }

    template <typename T>
    std::vector<T>
    gather(const MPI_Comm &   comm,
           const T &          object_to_send,
           const unsigned int root_process)
    {
#  ifndef DEAL_II_WITH_MPI
      (void)comm;
      (void)root_process;
      std::vector<T> v(1, object_to_send);
      return v;
#  else
      const auto n_procs = dealii::Utilities::MPI::n_mpi_processes(comm);
      const auto my_rank = dealii::Utilities::MPI::this_mpi_process(comm);

      AssertIndexRange(root_process, n_procs);

      std::vector<char> buffer       = Utilities::pack(object_to_send);
      int               n_local_data = buffer.size();

      // Vector to store the size of loc_data_array for every process
      // only the root process needs to allocate memory for that purpose
      std::vector<int> size_all_data;
      if (my_rank == root_process)
        size_all_data.resize(n_procs, 0);

      // Exchanging the size of each buffer
      int ierr = MPI_Gather(&n_local_data,
                            1,
                            MPI_INT,
                            size_all_data.data(),
                            1,
                            MPI_INT,
                            root_process,
                            comm);
      AssertThrowMPI(ierr);

      // Now computing the displacement, relative to recvbuf,
      // at which to store the incoming buffer; only for root
      std::vector<int> rdispls;
      if (my_rank == root_process)
        {
          rdispls.resize(n_procs, 0);
          for (unsigned int i = 1; i < n_procs; ++i)
            rdispls[i] = rdispls[i - 1] + size_all_data[i - 1];
        }
      // exchange the buffer:
      std::vector<char> received_unrolled_buffer;
      if (my_rank == root_process)
        received_unrolled_buffer.resize(rdispls.back() + size_all_data.back());

      ierr = MPI_Gatherv(buffer.data(),
                         n_local_data,
                         MPI_CHAR,
                         received_unrolled_buffer.data(),
                         size_all_data.data(),
                         rdispls.data(),
                         MPI_CHAR,
                         root_process,
                         comm);
      AssertThrowMPI(ierr);

      std::vector<T> received_objects;

      if (my_rank == root_process)
        {
          received_objects.resize(n_procs);

          for (unsigned int i = 0; i < n_procs; ++i)
            {
              const std::vector<char> local_buffer(
                received_unrolled_buffer.begin() + rdispls[i],
                received_unrolled_buffer.begin() + rdispls[i] +
                  size_all_data[i]);
              received_objects[i] = Utilities::unpack<T>(local_buffer);
            }
        }
      return received_objects;
#  endif
    }



    template <typename T>
    T
    broadcast(const MPI_Comm &   comm,
              const T &          object_to_send,
              const unsigned int root_process)
    {
#  ifndef DEAL_II_WITH_MPI
      (void)comm;
      (void)root_process;
      return object_to_send;
#  else
      const auto n_procs = dealii::Utilities::MPI::n_mpi_processes(comm);
      AssertIndexRange(root_process, n_procs);
      (void)n_procs;

      std::vector<char> buffer;
      unsigned int      buffer_size = numbers::invalid_unsigned_int;

      // On the root process, pack the data and determine what the
      // buffer size needs to be.
      if (this_mpi_process(comm) == root_process)
        {
          buffer      = Utilities::pack(object_to_send, false);
          buffer_size = buffer.size();
        }

      // Exchange the size of buffer
      int ierr = MPI_Bcast(&buffer_size, 1, MPI_UNSIGNED, root_process, comm);
      AssertThrowMPI(ierr);

      // If not on the root process, correctly size the buffer to
      // receive the data, then do exactly that.
      if (this_mpi_process(comm) != root_process)
        buffer.resize(buffer_size);

      ierr =
        MPI_Bcast(buffer.data(), buffer_size, MPI_CHAR, root_process, comm);
      AssertThrowMPI(ierr);

      if (Utilities::MPI::this_mpi_process(comm) == root_process)
        return object_to_send;
      else
        return Utilities::unpack<T>(buffer, false);
#  endif
    }


#  ifdef DEAL_II_WITH_MPI
    template <class Iterator, typename Number>
    std::pair<Number, typename numbers::NumberTraits<Number>::real_type>
    mean_and_standard_deviation(const Iterator  begin,
                                const Iterator  end,
                                const MPI_Comm &comm)
    {
      // below we do simple and straight-forward implementation. More elaborate
      // options are:
      // http://dx.doi.org/10.1145/2807591.2807644 section 3.1.2
      // https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Welford's_online_algorithm
      // https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Online
      using Std        = typename numbers::NumberTraits<Number>::real_type;
      const Number sum = std::accumulate(begin, end, Number(0.));

      const auto size = Utilities::MPI::sum(std::distance(begin, end), comm);
      Assert(size > 0, ExcDivideByZero());
      const Number mean =
        Utilities::MPI::sum(sum, comm) / static_cast<Std>(size);
      Std sq_sum = 0.;
      std::for_each(begin, end, [&mean, &sq_sum](const Number &v) {
        sq_sum += numbers::NumberTraits<Number>::abs_square(v - mean);
      });
      sq_sum = Utilities::MPI::sum(sq_sum, comm);
      return std::make_pair(mean,
                            std::sqrt(sq_sum / static_cast<Std>(size - 1)));
    }
#  endif

#endif
  } // end of namespace MPI
} // end of namespace Utilities


DEAL_II_NAMESPACE_CLOSE

#endif


