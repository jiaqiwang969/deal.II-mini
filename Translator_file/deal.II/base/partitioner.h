//include/deal.II-translator/base/partitioner_0.txt
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

#ifndef dealii_partitioner_h
#define dealii_partitioner_h

#include <deal.II/base/config.h>

#include <deal.II/base/array_view.h>
#include <deal.II/base/communication_pattern_base.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/memory_space.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/types.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/vector_operation.h>

#include <limits>


DEAL_II_NAMESPACE_OPEN

namespace Utilities
{
  namespace MPI
  {
    /**
     这个类定义了一个模型，用于在使用MPI的处理器之间分割一个向量（或者，事实上，任何线性数据结构）。        分区器将全局向量大小和本地拥有的范围作为一个半开放的区间[  @p lower,   @p upper) 存储在每个进程上。    此外，它还包括一个点对点通信模式的结构。它允许通过IndexSet包含幽灵指数（即当前处理器需要访问的指数，但被另一个进程拥有）。此外，它还存储了属于当前处理器的其他处理器的幽灵索引（见import_targets()），这些是其他处理器可能需要信息的索引。从某种意义上说，这些导入指数构成了幽灵指数的对偶。这些信息在构建分区器时被收集一次，这就避免了在交换数据时的后续全局通信步骤。        下图给出了一个索引空间 $[0,74)$ 被分割成四个部分的例子，这四个部分分别由一个MPI进程拥有。      @image html partitioner.png  第一行（粗黑线上方）显示哪个进程拥有哪些元素。在它下面，接下来的四行显示了每个处理器想要了解整个数组的哪些元素
     *
     * - 这通常是本地拥有的元素的超集，区别在于所谓的 "幽灵元素"。        为了理解图中的其余部分（以及这个类），请记住，在MPI中，你不能只是向另一个进程索取数据。    这并不完全正确：在较新的MPI标准中有这样的机制，但作为一般规则，这是真的）。相反，如果你需要信息，你需要向另一个进程发送一个消息，另一个进程需要期待该消息并以自己的消息作出适当的回应。因此，在实践中，如果每个进程都已经*知道它将被问到什么，并在适当的时候，只需发送该数据，就会更容易和更快。余下的几行信息就设定了这种方案。        为此，请注意，进程0想知道它自己不拥有的五个元素：那些指数为20、21的元素（由进程1拥有）；以及40、41、43（由进程2拥有）。类似的信息可以从以下几行得到。为了满足这种知识需求，因此，如果进程1存储了在适当的时候，它将不得不把元素20、21发送给进程0，这将是相当有用的。同样，如果你通过粗黑线下面的第2-4行，你会看到进程0应该知道它将需要把元素1、2、13、18、19发送到进程1；18、19发送到进程2；元素1、2、13发送到进程3。这些被称为 "导入索引"，因为其他进程会想要导入它们。与其将这些索引作为一个集合来存储，不如使用半开放的索引集来代替，因此上面列出的导入索引形成了以下集合：`[1,3)`, `[13,14)`, `[18,20)`, `[1,3)`, `[13,14]`。这就是上图中进口指数的显示方式。我们现在只需要知道这些半开的集合中哪些是要发给谁的。这是在上面一行完成的，我们在这里列出了 "导入目标"（即导入操作的目标进程）。进程1将收到5个元素（由前三个半开放目标索引集组成），进程2将收到2个索引（第四个半开放区间），进程3将收到3个索引（其余两个半开放区间）。这些信息被编码为`{1,5}`、`{2,2}`、`{3,3}`这几对作为导入目标。对于进程1、2、3要发送的内容也可以做类似的考虑。        最后，在接收信息时，知道每个进程将从谁那里接收多少个索引是很有用的，因为这样就可以预先分配合适大小的缓冲区。这在最后一行的 "幽灵目标 "中列出。进程0将从进程1接收两个元素（即索引20、21的元素），从进程2接收三个元素（即索引40、41、43的元素）。    这被编码为一对`{1,2}`和`{2,3}`。同样，对于进程1、2、3应该期待什么，然后在其各自的列中显示什么，也可以做类似的考虑。        这个类的主要目的是设置这些数据结构，只知道哪个进程拥有哪些元素，以及每个进程需要了解哪些额外的幽灵元素。            <h4>Local and global numbering</h4> 分区器包括一个将全局索引转换为局部索引和局部索引转换为全局索引的机制。在内部，这个类使用如下的约定来存储向量元素。本地范围与本地索引[0, locally_owned_size()]相关联，而鬼魂索引连续存储在[local_owned_size(), locally_owned_size() + n_ghost_indices()]。鬼魂索引根据其全局索引进行排序。            <h4>Parallel data exchange</h4> 这个类也处理对象分区数组的幽灵数据交换问题
     *
     * - 即，如果一个单独的类在每个进程上存储本地拥有的元素，那么这个类就可以方便地导入那些本地需要的、但存储在其他地方的元素，作为ghost。使用这个类的一个例子是 LinearAlgebra::distributed::Vector 类。        数据交换通过四个函数发生。      <ul>   <li>  export_to_ghosted_array_start()用于启动一个导出操作，将数据从本地拥有的数据域，以数组形式传递，发送到其他进程的ghost数据数组（根据存储在本类的ghost索引）。    这个调用启动了非阻塞的MPI发送通信例程，但并不等待例程的完成。因此，用户可能不会向底层数组的相应位置写入数据，因为这些数据可能仍然被MPI需要。      <li>  export_to_ghosted_array_finish()最终完成在export_to_ghosted_array_start()中开始的MPI数据交换，并发出信号，表示数组中的数据可以用于进一步处理或适当地修改。      <li>  import_from_ghosted_array_start()用于启动一个导入操作，根据存储在本类中的ghost索引，将数据从作为数组传递的ghost数据域发送到本地拥有的数组中。可以通过一个 VectorOperation::values 标志来决定如何将鬼域中的数据与所有者的数据结合起来，因为两者都与同一个数据条目有关。在汇编中，这通常是一个加到操作，即需要把所有进程对本地拥有的元素的贡献加起来。这个调用启动了非阻塞的MPI通信例程，但并不等待例程的完成。因此，用户可能不会向底层数组的各个位置写入数据，因为MPI可能仍然需要这些数据。      <li>  import_from_ghosted_array_finish()最终完成在import_from_ghosted_array_start()中开始的MPI数据交换，并发出信号，表示数组中的数据可以用于进一步处理或适当地修改。      </ul>  MPI通信例程是点对点通信模式。            <h4>Sending only selected ghost data</h4> 这个分区器类对一组固定的鬼魂索引进行操作，并且必须始终与它所代表的分区数组内的鬼魂索引兼容。在某些情况下，人们只想把存在于一个向量中的一些鬼魂索引送来送去，但不需要创建一个具有合适索引集的向量副本
     *
     * - 想想例如局部时间步进，其中一个向量的不同区域可能在一个时间步进片的不同阶段被交换。该类通过以下模型支持这种情况。一个向量首先以完整的重影索引集被创建。然后，创建第二个Partitioner实例，用更紧密的索引集来设置鬼魂索引，但指定更大的索引集作为set_ghost_indices()调用的第二个参数。当数据被交换时，export_to_ghosted_array_start()和import_f_ghosted_array_start()会检测到这种情况，并且只发送选定的索引，从完整的ghost条目阵列中提取。
     *
     */
    class Partitioner : public Utilities::MPI::CommunicationPatternBase
    {
    public:
      /**
       * 默认构造函数。
       *
       */
      Partitioner();

      /**
       * 带有大小参数的构造函数。创建一个MPI_COMM_SELF结构，其中没有真正的并行布局。
       *
       */
      Partitioner(const unsigned int size);

      /**
       * 构造函数，接收本地拥有的自由度数量 @p local_size
       * 和幽灵自由度数量 @p ghost_size.
       * 本地索引范围以升序和一对一的方式转化为全局索引，即进程
       * $p$ 的索引正好位于进程 $p-1$ 和 $p+1$
       * 的索引之间，分别。
       * @note  将 @p ghost_size
       * 变量设置为适当的值，在向量的内存分配中为幽灵数据提供内存空间，并允许通过local_element()访问它。然而，在这种情况下，相关的全局索引必须由外部来处理。
       *
       */
      Partitioner(const types::global_dof_index local_size,
                  const types::global_dof_index ghost_size,
                  const MPI_Comm &              communicator);

      /**
       * 带有索引集参数的构造函数。这个构造函数基于给定的通信器创建一个分布式布局，一个IndexSet描述本地拥有的范围，另一个用于描述由其他处理器拥有的、但我们需要有读或写权限的幽灵指数。
       *
       */
      Partitioner(const IndexSet &locally_owned_indices,
                  const IndexSet &ghost_indices_in,
                  const MPI_Comm &communicator_in);

      /**
       * 带有一个索引集参数的构造函数。这个构造函数基于给定的通信器创建一个分布式布局，以及一个描述本地拥有范围的IndexSet。它允许在以后的时间里设置鬼魂索引。除此以外，它与其他具有两个索引集的构造函数类似。
       *
       */
      Partitioner(const IndexSet &locally_owned_indices,
                  const MPI_Comm &communicator_in);

      /**
       * 重新初始化通信模式。第一个参数`vector_space_vector_index_set`是与一个VectorSpaceVector对象相关的索引集。第二个参数`read_write_vector_index_set`是与一个ReadWriteVector对象相关的索引集。
       *
       */
      virtual void
      reinit(const IndexSet &vector_space_vector_index_set,
             const IndexSet &read_write_vector_index_set,
             const MPI_Comm &communicator) override;

      /**
       * 设置本地拥有的索引。在构造函数中使用。
       *
       */
      void
      set_owned_indices(const IndexSet &locally_owned_indices);

      /**
       * 在构造函数被调用后，设置幽灵索引。
       * 可选的参数 @p larger_ghost_index_set
       * 允许定义一个间接寻址到一个更大的鬼魂索引集。如果一个分布式向量是基于那个更大的鬼魂索引集，但根据
       * @p ghost_indices.
       * ，只有一个更紧密的子集应该被传达，那么这种设置就很有用。
       *
       */
      void
      set_ghost_indices(const IndexSet &ghost_indices,
                        const IndexSet &larger_ghost_index_set = IndexSet());

      /**
       * 返回全局大小。
       *
       */
      types::global_dof_index
      size() const;

      /**
       * 返回本地拥有的索引数，即local_range().second减去local_range().first。
       * 当所有进程的索引数相加时，返回的数字需要与索引总数相加
       * @deprecated
       * 使用命名更明确的函数local_owned_size()代替。
       *
       */
      DEAL_II_DEPRECATED
      unsigned int
      local_size() const;

      /**
       * 返回本地拥有索引的数量，即local_range().second减去local_range().first。
       * 当所有进程的索引数相加时，返回的数字需要与索引总数相加。
       *
       */
      unsigned int
      locally_owned_size() const;

      /**
       * 返回一个本地范围的IndexSet表示。这个类只支持连续的局部范围，所以IndexSet实际上只由一个单一的数据范围组成，并等同于local_range()的结果。
       *
       */
      const IndexSet &
      locally_owned_range() const;

      /**
       * 返回本地范围。返回的数据对由第一个元素的索引和本地最后一个元素之后的一个元素的索引组成。
       *
       */
      std::pair<types::global_dof_index, types::global_dof_index>
      local_range() const;

      /**
       * 如果给定的全局索引在这个处理器的本地范围内，返回true。
       *
       */
      bool
      in_local_range(const types::global_dof_index global_index) const;

      /**
       * 返回与给定全局索引对应的本地索引。如果给定的全局索引既不是本地拥有的，也不是ghost，则会抛出一个异常。
       * 注意，对于本地拥有的索引，返回的本地索引将在0和local_owned_size()之间。
       *
       * - 1之间，而鬼魂的本地索引在local_owned_size()和local_owned_size()+ n_ghost_indices()之间。
       *
       * - 1.
       *
       */
      unsigned int
      global_to_local(const types::global_dof_index global_index) const;

      /**
       * 返回对应于给定的本地索引的全局索引。
       * 注意，本地拥有的索引的本地索引必须在0和local_owned_size()之间。
       *
       * - 1，鬼魂的本地索引必须在local_owned_size()和local_owned_size() + n_ghost_indices()之间。
       *
       * - 1.
       *
       */
      types::global_dof_index
      local_to_global(const unsigned int local_index) const;

      /**
       * 返回给定的全局索引是否是当前处理器上的一个幽灵索引。对于本地拥有的索引和根本不存在的索引返回false。
       *
       */
      bool
      is_ghost_entry(const types::global_dof_index global_index) const;

      /**
       * 返回所有幽灵索引的IndexSet表示。
       *
       */
      const IndexSet &
      ghost_indices() const;

      /**
       * 返回鬼魂索引的数量。与ghost_indices().n_elements()相同，但为了更简单的访问而进行缓存。
       *
       */
      unsigned int
      n_ghost_indices() const;

      /**
       * 如果分区器是为了将ghost索引定义为一个更大的ghost集合中的索引子集而建立的，这个函数会返回该集合中的范围的编号。与IndexSet中的结构类似，但为迭代而定制。
       * 如果分区器没有考虑到第二组鬼魂索引，这个子集就简单地定义为半开区间
       * <code>[0, n_ghost_indices())</code>  .
       *
       */
      const std::vector<std::pair<unsigned int, unsigned int>> &
      ghost_indices_within_larger_ghost_set() const;

      /**
       * 返回一个处理器的列表（第一项）和该处理器拥有的鬼魂自由度的数量（第二项）。后者在所有处理器上的总和等于n_ghost_indices（）。
       *
       */
      const std::vector<std::pair<unsigned int, unsigned int>> &
      ghost_targets() const;

      /**
       * 返回一个我们在compress()过程中导入的本地指数的范围向量，即属于本地范围的其他人的鬼点子。
       * 与IndexSet中的结构类似，但为迭代而量身定做，一些指数可能会被重复。返回的对包括第一个元素的索引和一个范围内的最后一个元素的索引。
       *
       */
      const std::vector<std::pair<unsigned int, unsigned int>> &
      import_indices() const;

      /**
       * 导入指数的数量，即在其他处理器上是幽灵的指数，我们将从这些指数中接收数据。
       *
       */
      unsigned int
      n_import_indices() const;

      /**
       * 返回一个处理器的列表（第一项）和在compress()操作中从它那里导入的自由度的数量（第二项），所有的处理器的数据都是从那里获得的，即在其他处理器上是幽灵的本地拥有的指数。
       * @note
       * 返回的向量只包含那些第二项为非零的处理器ID。
       *
       */
      const std::vector<std::pair<unsigned int, unsigned int>> &
      import_targets() const;

      /**
       * 检查给定的分区器是否与当前的分区器兼容。如果两个分区器有相同的本地大小和相同的ghost索引，那么它们就是兼容的。它们不一定需要对应于基于这些分区器对象而存储的相同数据。这只是一个局部操作，即如果只有一些处理器决定分区不兼容，只有这些处理器会返回
       * @p false, ，而其他处理器会返回 @p true.  。
       *
       */
      bool
      is_compatible(const Partitioner &part) const;

      /**
       * 检查给定的分区器是否与当前的分区器兼容。如果两个分区器具有相同的本地大小和相同的ghost索引，那么它们就是兼容的。它们不一定需要对应于基于这些分区器对象而存储的相同数据。与is_compatible()相反，该方法检查所有处理器之间的兼容性，如果分区器在所有处理器上是相同的，该方法只返回
       * @p true
       * 。换句话说，它在所有参与的进程上对is_compatible()返回的结果进行全局
       * "和 "操作。
       * 这个方法执行全局通信，所以要确保只在所有处理器调用它的次数相同的情况下使用它。
       *
       */
      bool
      is_globally_compatible(const Partitioner &part) const;

      /**
       * 返回调用处理器的MPI ID。缓存，以便有简单的访问。
       *
       */
      unsigned int
      this_mpi_process() const;

      /**
       * 返回参与给定分区器的MPI处理器的总数。缓存以获得简单的访问。
       *
       */
      unsigned int
      n_mpi_processes() const;

      /**
       * 返回分区器对象背后的MPI通信器。
       *
       */
      virtual const MPI_Comm &
      get_mpi_communicator() const override;

      /**
       * 返回是否已经明确添加了作为 @p
       * ghost_indices参数的幽灵索引。如果reinit()调用或构造函数提供了该参数，则仅为真。
       *
       */
      bool
      ghost_indices_initialized() const;

#ifdef DEAL_II_WITH_MPI
      /**
       * 开始将本地拥有的数组中的数据导出到该类的ghost索引所描述的范围。
       * @param  communication_channel
       * 为MPI_Isend和MPI_Irecv调用设置一个偏移量，避免与其他正在进行的不同条目上的export_to_ghosted_array_start()调用产生干扰。通常是在块向量的块内处理。任何小于200的值都是有效值。
       * @param  locally_owned_array
       * 从中提取数据并发送给远程处理器上的幽灵条目的数据阵列。
       * @param  temporary_storage
       * 一个长度为n_import_indices()的临时存储数组，用于存放要发送的
       * @p
       * locally_owned_array中的打包数据。注意，在各自的export_to_ghosted_array_finish()调用之前，这个数组不能被触及，因为该模型使用非阻塞通信。
       * @param  ghost_array
       * 将接收导出数据的数组，也就是远程处理器发送给调用进程的条目。它的大小必须是n_ghost_indices()或者等于作为set_ghost_indices()第二个参数的较大索引集中的ghost索引数。如果只发送选定的索引，我们不保证没有被设置的条目。其中一些可能被用来组织传输，随后被重置为零，所以要确保不在计算中使用它们。
       * @param  requests
       * 正在进行的非阻塞通信的MPI请求列表，将在export_to_ghosted_array_finish()调用中最终完成。
       * 这个功能在
       * LinearAlgebra::distributed::Vector::update_ghost_values().
       * 中使用。
       *
       */
      template <typename Number, typename MemorySpaceType = MemorySpace::Host>
      void
      export_to_ghosted_array_start(
        const unsigned int                              communication_channel,
        const ArrayView<const Number, MemorySpaceType> &locally_owned_array,
        const ArrayView<Number, MemorySpaceType> &      temporary_storage,
        const ArrayView<Number, MemorySpaceType> &      ghost_array,
        std::vector<MPI_Request> &                      requests) const;

      /**
       * 完成将本地拥有的数组中的数据导出到该类的ghost索引所描述的范围。
       * @param  ghost_array 将接收 @p export_to_ghosted_array_start().
       * 中开始的导出数据的数组
       * 这必须是传递给该函数的同一个数组，否则行为就无法定义。
       * @param  requests
       * 正在进行的非阻塞通信的MPI请求列表，这些请求是在export_to_ghosted_array_start()调用中启动的。这必须是传递给该函数的同一个数组，否则MPI可能会抛出一个错误。
       * 这个功能在
       * LinearAlgebra::distributed::Vector::update_ghost_values().
       * 中使用。
       *
       */
      template <typename Number, typename MemorySpaceType = MemorySpace::Host>
      void
      export_to_ghosted_array_finish(
        const ArrayView<Number, MemorySpaceType> &ghost_array,
        std::vector<MPI_Request> &                requests) const;

      /**
       * 开始导入由该类的鬼魂索引组成的数组上的数据，该数组后来用import_from_ghosted_array_finish()累积成一个本地拥有的数组。
       * @param  vector_operation
       * 定义了发送给所有者的数据应该如何与现有的条目相结合，例如，添加到。
       * @param  communication_channel
       * 设置MPI_Isend和MPI_Irecv调用的偏移量，避免与不同条目上其他正在进行的import_from_ghosted_array_start()调用产生干扰。通常是在块向量的块内处理。
       * 任何小于200的值都是有效值。              @param
       * ghost_array
       * 发送给向量中相应索引的远程所有者的幽灵数据数组。它的大小必须是n_ghost_indices()，或者等于作为set_ghost_indices()第二个参数的较大索引集中的鬼魂索引数。这个函数或随后的import_from_ghosted_array_finish()函数（顺序取决于实现）将把
       * @p  ghost_array后面的所有数据项设置为0。
       * @param  temporary_storage
       * 一个长度为n_import_indices()的临时存储数组，用于存放MPI通信的打包数据，这些数据以后将被写入本地拥有的数组中。注意，在各自的import_from_ghosted_array_finish()调用之前，这个数组不得被触及，因为该模型使用非阻塞通信。
       * @param  requests
       * 正在进行的非阻塞通信的MPI请求列表，将在export_to_ghosted_array_finish()调用中最终完成。
       * 这个功能在 LinearAlgebra::distributed::Vector::compress().
       * 中使用。
       *
       */
      template <typename Number, typename MemorySpaceType = MemorySpace::Host>
      void
      import_from_ghosted_array_start(
        const VectorOperation::values             vector_operation,
        const unsigned int                        communication_channel,
        const ArrayView<Number, MemorySpaceType> &ghost_array,
        const ArrayView<Number, MemorySpaceType> &temporary_storage,
        std::vector<MPI_Request> &                requests) const;

      /**
       * 完成从由该类的ghost索引的数组中导入数据到指定的本地拥有的数组中，根据给定的输入来组合结果
       * @p vector_operation.   @param  vector_operation
       * 定义发送给所有者的数据应该如何与现有的条目相结合，例如，加入到。
       * @param  temporary_storage
       * 给予import_from_ghosted_array_start()调用的相同数组，包含MPI通信的打包数据。在这样的函数中，根据
       * @p vector_operation.   @param  ghost_array
       * 在一个向量中的相应索引，被发送到远程所有者的鬼魂数据阵列。它的大小必须是n_ghost_indices()，或者等于作为set_ghost_indices()第二个参数的较大索引集中的鬼魂索引数。在import_from_ghosted_array_start()调用中尚未完成的情况下，该函数将把
       * @p ghost_array 后面的所有数据项设置为零。
       * @param  locally_owned_storage
       * 由远程进程发送至调用进程的结果数据将被累积到的数据阵列。
       * @param  requests
       * 在import_to_ghosted_array_finish()调用中已经启动的用于持续非阻塞通信的MPI请求列表。这必须是传递给该函数的同一个数组，否则MPI可能会抛出一个错误。
       * 这个功能在 LinearAlgebra::distributed::Vector::compress().
       * 中使用。
       *
       */
      template <typename Number, typename MemorySpaceType = MemorySpace::Host>
      void
      import_from_ghosted_array_finish(
        const VectorOperation::values                   vector_operation,
        const ArrayView<const Number, MemorySpaceType> &temporary_storage,
        const ArrayView<Number, MemorySpaceType> &      locally_owned_storage,
        const ArrayView<Number, MemorySpaceType> &      ghost_array,
        std::vector<MPI_Request> &                      requests) const;
#endif

      /**
       * 计算这个结构的内存消耗。
       *
       */
      std::size_t
      memory_consumption() const;

      /**
       * 异常情况
       *
       */
      DeclException2(ExcIndexNotPresent,
                     types::global_dof_index,
                     unsigned int,
                     << "Global index " << arg1
                     << " neither owned nor ghost on proc " << arg2 << ".");

      /**
       * 异常情况
       *
       */
      DeclException3(ExcGhostIndexArrayHasWrongSize,
                     unsigned int,
                     unsigned int,
                     unsigned int,
                     << "The size of the ghost index array (" << arg1
                     << ") must either equal the number of ghost in the "
                     << "partitioner (" << arg2
                     << ") or be equal in size to a more comprehensive index"
                     << "set which contains " << arg3
                     << " elements for this partitioner.");

    private:
      /**
       * 从import_indices_data初始化import_indices_plain_dev。这个函数只在使用CUDA感知的MPI时使用。
       *
       */
      void
      initialize_import_indices_plain_dev() const;

      /**
       * 所有处理器上的向量的全局大小
       *
       */
      types::global_dof_index global_size;

      /**
       * 本地存储的向量的范围。
       *
       */
      IndexSet locally_owned_range_data;

      /**
       * 本地存储的向量的范围。由于性能原因，从local_owned_range中提取。
       *
       */
      std::pair<types::global_dof_index, types::global_dof_index>
        local_range_data;

      /**
       * 我们需要有读取权限但不属于本地的索引集。
       *
       */
      IndexSet ghost_indices_data;

      /**
       * 一个缓存鬼魂索引数量的变量。使用 @p
       * ghost_indices.n_elements() 来计算这个会很昂贵。
       *
       */
      unsigned int n_ghost_indices_data;

      /**
       * 一个数组，包含我的幽灵索引属于哪个处理器以及这些索引的数量。
       *
       */
      std::vector<std::pair<unsigned int, unsigned int>> ghost_targets_data;

      /**
       * 我们在compress()过程中导入的（本地）索引集，即属于本地范围的其他人的ghost。与IndexSet中的结构类似，但为迭代而量身定做，有些指数可能是重复的。
       *
       */
      std::vector<std::pair<unsigned int, unsigned int>> import_indices_data;

      /**
       * 我们在compress()过程中导入的（本地）指数集，也就是属于本地范围的其他人的鬼魂。存储的数据与import_indices_data中的数据相同，但数据是以纯数组展开的。这个变量只在使用CUDA感知的MPI时使用。
       *
       */
      // The variable is mutable to enable lazy initialization in
      // export_to_ghosted_array_start(). This way partitioner does not have to
      // be templated on the MemorySpaceType.
      mutable std::vector<
        std::pair<std::unique_ptr<unsigned int[], void (*)(unsigned int *)>,
                  unsigned int>>
        import_indices_plain_dev;

      /**
       * 一个缓存鬼魂索引数量的变量。如果通过迭代导入索引并将其累积起来来计算，会很昂贵。
       *
       */
      unsigned int n_import_indices_data;

      /**
       * 向我们发送幽灵数据的处理器集合和数据字段的长度。
       *
       */
      std::vector<std::pair<unsigned int, unsigned int>> import_targets_data;

      /**
       * 一个数组，用于缓存每个MPI等级的导入索引中的块数。其长度为
       * import_indices_data.size()+1。
       *
       */
      std::vector<unsigned int> import_indices_chunks_by_rank_data;

      /**
       * 一个变量，用于缓存由set_ghost_indices()的可选参数给出的较大的索引集中的鬼魂索引数量。
       *
       */
      unsigned int n_ghost_indices_in_larger_set;

      /**
       * 一个数组，用于缓存每个MPI等级的导入索引中的块数。其长度为ghost_indices_subset_data.size()+1。
       *
       */
      std::vector<unsigned int> ghost_indices_subset_chunks_by_rank_data;

      /**
       * 出现在作为较大集合的子集的IndexSet的索引集合。与IndexSet中的所有鬼魂索引的结构相似，但为迭代而定制。
       *
       */
      std::vector<std::pair<unsigned int, unsigned int>>
        ghost_indices_subset_data;

      /**
       * MPI网络中当前处理器的ID
       *
       */
      unsigned int my_pid;

      /**
       * 问题中活跃的处理器的总数量
       *
       */
      unsigned int n_procs;

      /**
       * 问题中涉及的MPI通信器
       *
       */
      MPI_Comm communicator;

      /**
       * 一个变量，存储鬼魂指数是否被明确设置。
       *
       */
      bool have_ghost_indices;
    };



     /*--------------------- Inline functions --------------------------------*/ 

#ifndef DOXYGEN

    inline types::global_dof_index
    Partitioner::size() const
    {
      return global_size;
    }



    inline const IndexSet &
    Partitioner::locally_owned_range() const
    {
      return locally_owned_range_data;
    }



    inline std::pair<types::global_dof_index, types::global_dof_index>
    Partitioner::local_range() const
    {
      return local_range_data;
    }



    inline unsigned int
    Partitioner::local_size() const
    {
      return locally_owned_size();
    }



    inline unsigned int
    Partitioner::locally_owned_size() const
    {
      types::global_dof_index size =
        local_range_data.second - local_range_data.first;
      Assert(size <= std::numeric_limits<unsigned int>::max(),
             ExcNotImplemented());
      return static_cast<unsigned int>(size);
    }



    inline bool
    Partitioner::in_local_range(
      const types::global_dof_index global_index) const
    {
      return (local_range_data.first <= global_index &&
              global_index < local_range_data.second);
    }



    inline bool
    Partitioner::is_ghost_entry(
      const types::global_dof_index global_index) const
    {
      // if the index is in the global range, it is trivially not a ghost
      if (in_local_range(global_index) == true)
        return false;
      else
        return ghost_indices().is_element(global_index);
    }



    inline unsigned int
    Partitioner::global_to_local(
      const types::global_dof_index global_index) const
    {
      Assert(in_local_range(global_index) || is_ghost_entry(global_index),
             ExcIndexNotPresent(global_index, my_pid));
      if (in_local_range(global_index))
        return static_cast<unsigned int>(global_index - local_range_data.first);
      else if (is_ghost_entry(global_index))
        return (locally_owned_size() +
                static_cast<unsigned int>(
                  ghost_indices_data.index_within_set(global_index)));
      else
        // should only end up here in optimized mode, when we use this large
        // number to trigger a segfault when using this method for array
        // access
        return numbers::invalid_unsigned_int;
    }



    inline types::global_dof_index
    Partitioner::local_to_global(const unsigned int local_index) const
    {
      AssertIndexRange(local_index,
                       locally_owned_size() + n_ghost_indices_data);
      if (local_index < locally_owned_size())
        return local_range_data.first + types::global_dof_index(local_index);
      else
        return ghost_indices_data.nth_index_in_set(local_index -
                                                   locally_owned_size());
    }



    inline const IndexSet &
    Partitioner::ghost_indices() const
    {
      return ghost_indices_data;
    }



    inline unsigned int
    Partitioner::n_ghost_indices() const
    {
      return n_ghost_indices_data;
    }



    inline const std::vector<std::pair<unsigned int, unsigned int>> &
    Partitioner::ghost_indices_within_larger_ghost_set() const
    {
      return ghost_indices_subset_data;
    }



    inline const std::vector<std::pair<unsigned int, unsigned int>> &
    Partitioner::ghost_targets() const
    {
      return ghost_targets_data;
    }


    inline const std::vector<std::pair<unsigned int, unsigned int>> &
    Partitioner::import_indices() const
    {
      return import_indices_data;
    }



    inline unsigned int
    Partitioner::n_import_indices() const
    {
      return n_import_indices_data;
    }



    inline const std::vector<std::pair<unsigned int, unsigned int>> &
    Partitioner::import_targets() const
    {
      return import_targets_data;
    }



    inline unsigned int
    Partitioner::this_mpi_process() const
    {
      // return the id from the variable stored in this class instead of
      // Utilities::MPI::this_mpi_process() in order to make this query also
      // work when MPI is not initialized.
      return my_pid;
    }



    inline unsigned int
    Partitioner::n_mpi_processes() const
    {
      // return the number of MPI processes from the variable stored in this
      // class instead of Utilities::MPI::n_mpi_processes() in order to make
      // this query also work when MPI is not initialized.
      return n_procs;
    }



    inline const MPI_Comm &
    Partitioner::get_mpi_communicator() const
    {
      return communicator;
    }



    inline bool
    Partitioner::ghost_indices_initialized() const
    {
      return have_ghost_indices;
    }

#endif // ifndef DOXYGEN

  } // end of namespace MPI

} // end of namespace Utilities


DEAL_II_NAMESPACE_CLOSE

#endif


