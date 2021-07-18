//include/deal.II-translator/lac/la_parallel_vector_0.txt
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

#ifndef dealii_la_parallel_vector_h
#define dealii_la_parallel_vector_h

#include <deal.II/base/config.h>

#include <deal.II/base/communication_pattern_base.h>
#include <deal.II/base/memory_space.h>
#include <deal.II/base/memory_space_data.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/numbers.h>
#include <deal.II/base/partitioner.h>
#include <deal.II/base/thread_management.h>

#include <deal.II/lac/vector_operation.h>
#include <deal.II/lac/vector_space_vector.h>
#include <deal.II/lac/vector_type_traits.h>

#include <iomanip>
#include <memory>

DEAL_II_NAMESPACE_OPEN

// Forward declarations
#ifndef DOXYGEN
namespace LinearAlgebra
{
  /**
   * 一个用于向量的并行实现的命名空间。
   *
   */
  namespace distributed
  {
    template <typename>
    class BlockVector;
  }

  template <typename>
  class ReadWriteVector;
} // namespace LinearAlgebra

#  ifdef DEAL_II_WITH_PETSC
namespace PETScWrappers
{
  namespace MPI
  {
    class Vector;
  }
} // namespace PETScWrappers
#  endif

#  ifdef DEAL_II_WITH_TRILINOS
namespace TrilinosWrappers
{
  namespace MPI
  {
    class Vector;
  }
} // namespace TrilinosWrappers
#  endif
#endif

namespace LinearAlgebra
{
  namespace distributed
  {
    /*!   @addtogroup Vectors 
     * @{     
*
*/

    /**
     * 并行向量类的实现。这个类的设计与deal.II中的标准 ::dealii::Vector 类相似，不同的是，存储是用MPI分布的。        矢量被设计为以下平行分区的方案。      <ul>   <li>  MPI并行化中各个进程（本地拥有的部分）所持有的索引形成一个连续的范围  <code>[my_first_index,my_last_index)</code>  。      <li>  驻留在其他处理器任意位置的幽灵指数是允许的。一般来说，如果鬼魂索引被聚在一起，效率会更高，因为它们被存储为一组间隔。在调用函数<code>reinit (local_owned, ghost_indices, communicator)</code>时，鬼魂索引的通信模式被确定，并保留到分区被改变。    这允许指数的有效并行通信。特别是，它存储了通信模式，而不是在每次通信时都要再次计算它。关于鬼魂向量的更多信息，也可参见 @ref GlossGhostedVector "关于有鬼魂元素的向量的词汇表条目"。      <li>  除了通常的全局访问操作符()，还可以用函数 @p  local_element()访问本地索引空间的向量条目。本地拥有的索引被放在首位，[0，local_owned_size()]，然后所有的鬼魂索引在它们之后连续排列，[local_owned_size(), locally_owned_size()+n_ghost_entries()]。      </ul>  与并行功能有关的函数。      <ul>   <li>  函数 <code>compress()</code> 穿过与鬼魂索引相关的数据，并将其传达给所有者进程，然后它可以将其添加到正确的位置。例如，这可以在运行了一个涉及到填充这个向量的鬼魂的汇编程序之后使用。注意， @p compress() 的 @p insert 模式并不设置包含在ghost条目中的元素，而是简单地丢弃它们，假设拥有的处理器已经将它们设置为所需的值（也见 @ref GlossCompress "关于压缩的词汇条"）。      <li>   <code>update_ghost_values()</code>  函数将数据从拥有处理器导入到ghost索引中，以便提供对与ghost相关的数据的读取访问。      <li>  可以将上述函数分成两个阶段，第一个阶段启动通信，第二个阶段完成通信。    这些函数可以用来将通信与代码中其他部分的计算重叠起来。      <li>  当然，还原操作（如规范）利用了集体的全对全MPI通信。      </ul>  这个向量相对于鬼魂元素来说可以采取两种不同的状态。      <ul>   <li>  在创建后，只要调用zero_out_ghost_values()(或 <code>operator= (0.)</code> )，该向量只允许向鬼魂元素写入，而不允许从鬼魂元素中读出。      <li>  在调用update_ghost_values()后，向量不允许写进ghost元素，只允许读出ghost元素。这是为了避免在修改一些向量条目后调用compress()时出现不希望出现的ghost数据假象。ghost条目的当前状态（读模式或写模式）可以通过has_ghost_elements()方法查询，当ghost元素被更新时，该方法准确地返回 <code>true</code> ，否则返回 <code>false</code> ，而不考虑向量布局中ghost条目的实际数量（要获得该信息，请使用n_ghost_entries()代替）。      </ul>  这个向量使用类 dealii::Vector<Number> 的设施来实现对向量的本地范围的操作。特别是，它还继承了线程并行性，如果程序使用多个线程，它可以将大多数向量-向量操作分割成更小的块。这在与MPI一起工作时可能需要也可能不需要。        <h4>Limitations regarding the vector size</h4> 这个向量类是基于两种不同的数字类型进行索引的。    所谓的全局索引类型是对向量的整体大小进行编码。    它的类型是 types::global_dof_index. ，最大的可能值是 <code>2^32-1</code> 或大约40亿，如果在配置deal.II时禁用64位整数（默认情况），或者如果启用64位整数，则为 <code>2^64-1</code> or approximately <code>10^19</code> （进一步信息见术语表中的 @ref GlobalDoFIndex 条目）。        第二个相关的索引类型是在一个MPI等级内使用的本地索引。与全局索引相反，该实现无条件地假定为32位无符号整数。换句话说，要实际使用一个超过40亿条目的向量，你需要使用一个以上等级的MPI（一般来说，这是一个安全的假设，因为40亿条目对于浮点数来说至少要消耗16GB的内存，对于双数来说至少要消耗32GB的内存），并启用64位索引。如果有超过40亿的局部元素存在，实现会尝试检测，这将触发一个异常并中止代码。然而，请注意，对溢出的检测是很棘手的，在某些情况下，检测机制可能会失败。因此，强烈建议不要依赖这个类来自动检测不支持的情况。        <h4>CUDA support</h4> 这个向量类支持两个不同的内存空间。主机和CUDA。默认情况下，内存空间是Host，所有的数据都分配在CPU上。当内存空间为CUDA时，所有的数据都被分配在GPU上。    对向量的操作是在选择的内存空间上进行的。    从主机来看，当使用CUDA内存空间时，有两种方法可以访问向量的元素。      <ul>   <li>  使用get_values（）。
     * @code
     * Vector<double, MemorySpace::CUDA> vector(local_range, comm);
     * double* vector_dev = vector.get_values();
     * std::vector<double> vector_host(local_range.n_elements(), 1.);
     * Utilities::CUDA::copy_to_dev(vector_host, vector_dev);
     * @endcode
     * <li>  使用import()。
     * @code
     * Vector<double, MemorySpace::CUDA> vector(local_range, comm);
     * ReadWriteVector<double> rw_vector(local_range);
     * for (auto & val : rw_vector)
     * val = 1.;
     * vector.import(rw_vector, VectorOperations::insert);
     * @endcode
     * </ul>  import方法要安全得多，必要时将执行MPI通信。因为可能会执行MPI通信，所以需要在所有处理器上调用import。
     * @note  默认情况下，所有的等级将尝试访问设备0。
     * 如果你每个节点有一个等级，每个节点有一个GPU，这就很好。如果你在一个节点上有多个GPU，我们需要每个进程访问一个不同的GPU。如果每个节点有相同数量的GPU，这可以按以下方式进行。
     * <code> int n_devices = 0; cudaGetDeviceCount(&n_devices); int device_id
     * = my_rank % n_devices; cudaSetDevice(device_id); </code> <h4>MPI-3
     * shared-memory support</h4>
     * 在主机模式下，该类允许使用MPI-3共享内存功能，提供一个单独的MPI通信器，由同一共享内存域中的进程组成。通过调用`vector.shared_vector_data();`，用户可以只读访问共享内存通信器中的进程的本地所有值和幽灵值（
     * @p comm_sm  in reinit() ）。
     * 要做到这一点，你必须调用这个类的构造函数或其中一个reinit()函数的`comm_sm`参数的非默认值，该参数对应于一个由同一共享内存域的所有进程组成的通信器。这种通信器可以用下面的代码片断来创建。
     * @code
     * MPI_Comm comm_sm;
     * MPI_Comm_split_type(comm, MPI_COMM_TYPE_SHARED, rank, MPI_INFO_NULL,
     *                     &comm_sm);
     * @endcode
     * @see  CUDAWrappers
     *
     */
    template <typename Number, typename MemorySpace = MemorySpace::Host>
    class Vector : public ::dealii::LinearAlgebra::VectorSpaceVector<Number>,
                   public Subscriptor
    {
    public:
      using memory_space    = MemorySpace;
      using value_type      = Number;
      using pointer         = value_type *;
      using const_pointer   = const value_type *;
      using iterator        = value_type *;
      using const_iterator  = const value_type *;
      using reference       = value_type &;
      using const_reference = const value_type &;
      using size_type       = types::global_dof_index;
      using real_type       = typename numbers::NumberTraits<Number>::real_type;

      static_assert(
        std::is_same<MemorySpace, ::dealii::MemorySpace::Host>::value ||
          std::is_same<MemorySpace, ::dealii::MemorySpace::CUDA>::value,
        "MemorySpace should be Host or CUDA");

      /**
       * @name  1: 基本对象处理
       *
       */
      //@{
      /**
       * 空的构造器。
       *
       */
      Vector();

      /**
       * 复制构造函数。使用 @p in_vector. 的并行分区
       * 应该注意的是，这个构造函数会自动将ghost值设置为零。如果需要一个重影向量，请在构造后直接调用
       * @p update_ghost_values() 。
       *
       */
      Vector(const Vector<Number, MemorySpace> &in_vector);

      /**
       * 构建一个给定的全局大小的平行向量，没有任何实际的平行分布。
       *
       */
      Vector(const size_type size);

      /**
       * 构建一个平行向量。本地范围由 @p local_owned_set指定（注意，这必须是一个连续的区间，多个区间是不可能的）。IndexSet  @p ghost_indices  指定了鬼魂索引，即人们可能需要从中读取数据或积累数据的索引。允许鬼魂索引集也包含本地范围，但不需要。            这个函数涉及到全局通信，所以对于一个给定的布局，它应该只被调用一次。使用带有Vector<Number>参数的构造函数来创建具有相同平行布局的额外向量。              @see   @ref GlossGhostedVector  "有幽灵元素的向量"
       *
       */
      Vector(const IndexSet &local_range,
             const IndexSet &ghost_indices,
             const MPI_Comm &communicator);

      /**
       * 与上述构造函数相同，但没有任何鬼魂索引。
       *
       */
      Vector(const IndexSet &local_range, const MPI_Comm &communicator);

      /**
       * 基于 @p
       * 分区器中描述的并行分区来创建向量。输入参数是一个共享指针，它只存储一次分区器数据，并在具有相同布局的几个向量之间共享它。
       *
       */
      Vector(
        const std::shared_ptr<const Utilities::MPI::Partitioner> &partitioner);

      /**
       * 解构器。
       *
       */
      virtual ~Vector() override;

      /**
       * 将向量的全局大小设置为 @p size
       * ，而没有任何实际的并行分布。
       *
       */
      void
      reinit(const size_type size, const bool omit_zeroing_entries = false);

      /**
       * 使用输入向量 @p in_vector
       * 的并行布局，并为该向量分配内存。当需要创建多个具有相同布局的向量时，推荐使用初始化函数。
       * 如果标志 @p omit_zeroing_entries
       * 被设置为false，内存将被初始化为0，否则内存将不被动用（用户在使用前必须确保用合理的数据填充）。
       *
       */
      template <typename Number2>
      void
      reinit(const Vector<Number2, MemorySpace> &in_vector,
             const bool                          omit_zeroing_entries = false);

      /**
       * 初始化向量。本地范围由 @p local_owned_set指定（注意，这必须是一个连续的区间，多个区间是不可能的）。IndexSet  @p ghost_indices  指定了鬼魂索引，即人们可能需要从中读取数据或积累数据的索引。允许鬼魂索引集也包含本地范围，但不需要。            这个函数涉及到全局通信，所以对于一个给定的布局，它应该只被调用一次。使用带有Vector<Number>参数的 @p reinit 函数来创建具有相同平行布局的额外向量。              @see   @ref GlossGhostedVector  "有幽灵元素的向量"
       *
       */
      void
      reinit(const IndexSet &local_range,
             const IndexSet &ghost_indices,
             const MPI_Comm &communicator);

      /**
       * 与上述相同，但没有鬼魂条目。
       *
       */
      void
      reinit(const IndexSet &local_range, const MPI_Comm &communicator);

      /**
       * 初始化给与 @p partitioner.
       * 中描述的并行分区的向量，输入参数是一个共享指针，它只存储一次分区器数据，并在具有相同布局的几个向量之间共享。
       * 可选的参数 @p comm_sm,
       * 由同一共享内存域上的进程组成，允许用户对共享内存通信器中合并的进程的本地拥有值和幽灵值进行只读访问。关于这个参数的更多信息，请参见该类的一般文档。
       *
       */
      void
      reinit(
        const std::shared_ptr<const Utilities::MPI::Partitioner> &partitioner,
        const MPI_Comm &comm_sm = MPI_COMM_SELF);

      /**
       * 用 @p local_size 本地拥有和 @p ghost_size
       * 幽灵自由度初始化向量。            可选的参数 @p
       * comm_sm,
       * 由同一共享内存域的进程组成，允许用户对共享内存通信器中的进程的本地拥有值和幽灵值进行只读访问。关于这个参数的更多信息，请参见该类的一般文档。
       * @note
       * 在创建的底层分区器中，本地索引范围以升序和一对一的方式转化为全局索引，即进程
       * $p$ 的索引正好位于进程 $p-1$ 和 $p+1$
       * 的索引之间，分别。将 @p ghost_size
       * 变量设置为一个适当的值，在向量的内存分配中为幽灵数据提供内存空间，并允许通过local_element()访问它。然而，在这种情况下，相关的全局索引必须由外部来处理。
       *
       */
      void
      reinit(const types::global_dof_index local_size,
             const types::global_dof_index ghost_size,
             const MPI_Comm &              comm,
             const MPI_Comm &              comm_sm = MPI_COMM_SELF);

      /**
       * 交换这个向量和另一个向量的内容  @p v.
       * 人们可以通过一个临时变量和复制数据元素来完成这个操作，但是这个函数明显更有效率，因为它只交换两个向量的数据指针，因此不需要分配临时存储和移动数据。
       * 这个函数类似于所有C++标准容器的 @p swap
       * 函数。此外，还有一个全局函数<tt>swap(u,v)</tt>，它简单地调用<tt>u.swap(v)</tt>，同样与标准函数相类似。
       *
       */
      void
      swap(Vector<Number, MemorySpace> &v);

      /**
       * 将向量分配给输入向量 @p in_vector,
       * 的并行分区，并复制所有数据。
       * 如果输入向量或调用向量中的一个（在赋值运算符的左边）在此操作之前设置了鬼魂元素，调用向量将设置鬼魂值。否则，它将处于写模式。如果输入向量根本没有任何鬼魂元素，该向量也将更新其鬼魂值，类似于Trilinos和PETSc向量的各自设置。
       *
       */
      Vector<Number, MemorySpace> &
      operator=(const Vector<Number, MemorySpace> &in_vector);

      /**
       * 将向量分配给输入向量 @p in_vector,
       * 的并行分区，并复制所有数据。
       * 如果输入向量或调用向量中的一个（在赋值运算符的左边）在此操作之前设置了鬼魂元素，调用向量将设置鬼魂值。否则，它将处于写模式。如果输入向量根本没有任何鬼魂元素，该向量也将更新其鬼魂值，类似于Trilinos和PETSc向量的各自设置。
       *
       */
      template <typename Number2>
      Vector<Number, MemorySpace> &
      operator=(const Vector<Number2, MemorySpace> &in_vector);

      //@}

      /**
       * @name  2：并行数据交换
       *
       */
      //@{
      /**
       * 这个函数将积累在数据缓冲器中的鬼魂索引的数据复制到拥有的处理器中。关于参数 @p operation, 的含义，见术语表中的 @ref GlossCompress "压缩分布式向量和矩阵 "
       * 条目。            这个函数有四个变体。如果用参数 @p
       * 调用 VectorOperation::add
       * ，则将所有积累在幽灵元素中的数据添加到拥有处理器的相应元素中，并在之后清除幽灵阵列。如果用参数
       * @p   VectorOperation::insert,
       * 调用，则执行设置操作。因为在一个有鬼魂元素的向量中设置元素是不明确的（因为可以同时设置鬼魂位置和拥有位置上的元素），这个操作假设所有的数据都在拥有处理器上正确设置。因此在调用
       * compress(VectorOperation::insert),
       * 时，所有的ghost条目都被简单地清零（使用zero_ghost_values()）。
       * 在调试模式下，会对数据集在处理器之间是否真的一致进行检查，也就是说，每当发现一个非零的ghost元素，就会与拥有处理器上的值进行比较，如果这些元素不一致，就会抛出一个异常。
       * 如果与 VectorOperation::min 或 VectorOperation::max,
       * 一起调用，将设置所有处理器上的元素的最小值或最大值。
       * @note
       * 这个向量类有一组固定的鬼魂条目附加到本地表示。因此，所有的ghost条目都被认为是有效的，并将根据给定的VectorOperation无条件地进行交换。确保用给定的VectorOperation的中性元素初始化所有的ghost条目，或者触摸所有的ghost条目。对于
       * VectorOperation::add 和 VectorOperation::insert,
       * 来说，中性元素是零 `+inf`对于 VectorOperation::min,
       * 来说，`-inf`对于 VectorOperation::max.
       * 来说，如果所有的值被初始化为低于零的值，并且随后用
       * VectorOperation::max
       * 调用compress两次，第二次计算后的最大值将是零。
       *
       */
      virtual void
      compress(::dealii::VectorOperation::values operation) override;

      /**
       * 用存储在所属处理器各自位置上的值来填充ghost索引的数据域。在从鬼魂中读取数据之前需要这个函数。该函数是 @p const ，即使ghost数据被改变。需要这样做是为了让具有 @p 常量向量的函数在不创建临时性的情况下进行数据交换。            调用这个方法后，禁止对向量的ghost元素进行写入访问，并抛出一个异常。在这种状态下，只允许对鬼魂元素进行读取访问。请注意，随后对这个向量的所有操作，比如全局向量加法等，也会在操作后通过调用这个方法来更新ghost值。然而，像规范或内积这样的全局缩减操作将总是忽略幽灵元素，以避免多次计算幽灵数据。为了允许再次写入ghost元素，可以调用zero_out_ghost_values()。              @see   @ref GlossGhostedVector  "有鬼魂元素的向量"
       *
       */
      void
      update_ghost_values() const;

      /**
       * 为 @p compress()
       * 函数启动非阻塞式通信。该函数不等待传输完成，以便在所有数据到达之前的时间内进行其他计算。
       * 在数据实际交换之前，该函数必须接着调用 @p
       * compress_finish(). 。如果在调用 @p
       * compress_finish()之前，该函数被调用超过一个向量，则必须为每个此类调用指定一个唯一的通信通道，以避免几个具有相同ID的消息破坏该操作。任何小于100的通信通道都是有效值（特别是，
       * $[100, 200)$ 范围保留给
       * LinearAlgebra::distributed::BlockVector). 。
       *
       */
      void
      compress_start(
        const unsigned int                communication_channel = 0,
        ::dealii::VectorOperation::values operation = VectorOperation::add);

      /**
       * 对于所有在compress_start中启动的请求，等待通信结束。一旦完成，将数据（取决于标志操作）添加或设置到拥有处理器的相应位置，并清除幽灵数据字段中的内容。
       * 这个参数的含义与compress()中的相同。
       * 这个函数应该在调用compress_start后，每个向量精确地调用一次，否则结果就无法定义。特别是，在调用compress_finished之前，对同一个向量再次调用compress_start是不好定义的。然而，没有警告来防止这种情况。
       * 必须在对 @p compress_start 函数的调用之后。
       * 当MemorySpace是CUDA并且MPI不支持CUDA时，调用compress_start后设备上改变的数据将丢失。
       *
       */
      void
      compress_finish(::dealii::VectorOperation::values operation);

      /**
       * 为 @p update_ghost_values()
       * 函数启动非阻塞式通信。这个函数不等待传输完成，以便在所有数据到达之前的时间内允许进行其他计算。
       * 在数据实际交换之前，该函数必须紧接着调用 @p
       * update_ghost_values_finish(). 。如果在调用 @p
       * update_ghost_values_finish()之前，该函数被调用超过一个向量，则必须为每个此类调用指定一个唯一的通信通道，以避免几个相同ID的消息破坏该操作。
       * 任何小于100的通信通道都是有效值（特别是， $[100,
       * 200)$ 范围保留给 LinearAlgebra::distributed::BlockVector). 。
       *
       */
      void
      update_ghost_values_start(
        const unsigned int communication_channel = 0) const;


      /**
       * 对于所有在update_ghost_values_start中启动的请求，等待通信结束。
       * 必须在调用 @p update_ghost_values_start
       * 函数之后，才能从ghost索引中读取数据。
       *
       */
      void
      update_ghost_values_finish() const;

      /**
       * 这个方法将ghost
       * dofs上的条目清零，但不触及本地拥有的DoF。
       * 调用此方法后，禁止对向量的ghost元素进行读取访问，并抛出一个异常。在这种状态下，只允许对鬼魂元素进行写入访问。
       * @deprecated  使用zero_out_ghost_values()代替。
       *
       */
      DEAL_II_DEPRECATED void
      zero_out_ghosts() const;

      /**
       * 这个方法将ghost
       * dofs上的条目清零，但不涉及本地拥有的DoF。
       * 调用该方法后，禁止对向量中的ghost元素进行读取访问，并且会抛出一个异常。在这种状态下，只允许对鬼魂元素进行写入访问。
       *
       */
      void
      zero_out_ghost_values() const;

      /**
       * 返回向量当前是否处于可以读取ghost值的状态。这与其他并行向量的功能相同。如果这个方法返回false，这只意味着对ghost元素的读访问被禁止，而写访问仍然是可能的（对那些在初始化时指定为ghost的条目），而不是说根本就没有ghost元素。              @see   @ref GlossGhostedVector  "有幽灵元素的向量"
       *
       */
      bool
      has_ghost_elements() const;

      /**
       * 该方法将本地拥有的范围内的数据从另一个分布式向量
       * @p src
       * 复制到调用向量中。与同样包括鬼魂条目的operator=相反，这个操作忽略了鬼魂范围。唯一的前提是调用向量和给定向量
       * @p src
       * 的本地范围在所有处理器上是相同的。明确允许两个向量有不同的重影元素，这些元素可能是也可能不是相互关联的。
       * 由于没有进行数据交换，确保 @p src
       * 和调用向量都没有待处理的通信，以便获得正确的结果。
       *
       */
      template <typename Number2>
      void
      copy_locally_owned_data_from(const Vector<Number2, MemorySpace> &src);

      /**
       * 导入分布式向量中存在的所有元素  @p src.
       * VectorOperation::values   @p operation  用于决定  @p V
       * 中的元素是否应该添加到当前向量中或替换当前元素。这个函数的主要目的是将数据从一个内存空间，如CUDA，获取到另一个内存空间，如Host。
       * @note
       * 两个分布式向量的分区器需要相同，因为没有进行MPI通信。
       *
       */
      template <typename MemorySpace2>
      void
      import(const Vector<Number, MemorySpace2> &src,
             VectorOperation::values             operation);

      //@}

      /**
       * @name  3：VectorSpaceVector的实现
       *
       */
      //@{

      /**
       * 将维度改为向量V的维度，V的元素不被复制。
       *
       */
      virtual void
      reinit(const VectorSpaceVector<Number> &V,
             const bool omit_zeroing_entries = false) override;

      /**
       * 将整个向量乘以一个固定的因子。
       *
       */
      virtual Vector<Number, MemorySpace> &
      operator*=(const Number factor) override;

      /**
       * 用整个向量除以一个固定的因子。
       *
       */
      virtual Vector<Number, MemorySpace> &
      operator/=(const Number factor) override;

      /**
       * 将向量 @p V 添加到现在的向量中。
       *
       */
      virtual Vector<Number, MemorySpace> &
      operator+=(const VectorSpaceVector<Number> &V) override;

      /**
       * 从现在的向量中减去向量 @p V 。
       *
       */
      virtual Vector<Number, MemorySpace> &
      operator-=(const VectorSpaceVector<Number> &V) override;

      /**
       * 从输入向量 @p V.  VectorOperation::values  @p operation
       * 中导入所有存在于向量IndexSet中的元素，用于决定 @p
       * V
       * 中的元素是否应该被添加到当前向量中，或者替换当前元素。如果多次使用同一通信模式，可以使用最后一个参数。这可以用来提高性能。
       * @note
       * 如果MemorySpace是CUDA，ReadWriteVector中的数据将被移动到设备上。
       *
       */
      virtual void
      import(const LinearAlgebra::ReadWriteVector<Number> &V,
             VectorOperation::values                       operation,
             std::shared_ptr<const Utilities::MPI::CommunicationPatternBase>
               communication_pattern = {}) override;

      /**
       * 返回两个向量的标量乘积。
       *
       */
      virtual Number
      operator*(const VectorSpaceVector<Number> &V) const override;

      /**
       * 将 @p a 加到所有分量上。注意 @p a
       * 是一个标量而不是一个向量。
       *
       */
      virtual void
      add(const Number a) override;

      /**
       * 向量的倍数的简单加法，即<tt>*this += a*V</tt>。
       *
       */
      virtual void
      add(const Number a, const VectorSpaceVector<Number> &V) override;

      /**
       * 缩放向量的多重加法，即：<tt>*this += a*V+b*W</tt>。
       *
       */
      virtual void
      add(const Number                     a,
          const VectorSpaceVector<Number> &V,
          const Number                     b,
          const VectorSpaceVector<Number> &W) override;

      /**
       * 一个集体加法操作。这个函数将存储在 @p values
       * 中的一整套数值添加到 @p indices.
       * 指定的向量成分中。
       *
       */
      virtual void
      add(const std::vector<size_type> &indices,
          const std::vector<Number> &   values);

      /**
       * 缩放和简单的向量倍数相加，即<tt>*this =
       * s*(*this)+a*V</tt>。
       *
       */
      virtual void
      sadd(const Number                     s,
           const Number                     a,
           const VectorSpaceVector<Number> &V) override;

      /**
       * 用参数中的相应元素来缩放这个向量的每个元素。这个函数主要是为了模拟对角线缩放矩阵的乘法（和立即重新分配）。
       *
       */
      virtual void
      scale(const VectorSpaceVector<Number> &scaling_factors) override;

      /**
       * 赋值 <tt>*this = a*V</tt>.
       *
       */
      virtual void
      equ(const Number a, const VectorSpaceVector<Number> &V) override;

      /**
       * 返回向量的l<sub>1</sub>规范（即所有处理器中所有条目的绝对值之和）。
       *
       */
      virtual real_type
      l1_norm() const override;

      /**
       * 返回向量的 $l_2$
       * 准则（即所有处理器中所有条目的平方之和的平方根）。
       *
       */
      virtual real_type
      l2_norm() const override;

      /**
       * 返回向量的 $l_2$ 规范的平方。
       *
       */
      real_type
      norm_sqr() const;

      /**
       * 返回向量的最大规范（即所有条目和所有处理器中的最大绝对值）。
       *
       */
      virtual real_type
      linfty_norm() const override;

      /**
       * 执行一个向量加法和后续内积的组合操作，返回内积的值。换句话说，这个函数的结果与用户调用的
       * @code
       * this->add(a, V);
       * return_value =this W;
       * @endcode
       * 这个函数存在的原因是这个操作比单独调用这两个函数涉及的内存转移要少。这个方法只需要加载三个向量，
       * @p this,   @p V,   @p W,
       * ，而调用单独的方法意味着要加载调用向量 @p this
       * 两次。由于大多数向量操作都有内存传输限制，这就使时间减少了25\%（如果
       * @p W 等于 @p this).
       * ，则减少50\%）对于复值向量，第二步中的标量乘法被实现为
       * $\left<v,w\right>=\sum_i v_i \bar{w_i}$  。
       *
       */
      virtual Number
      add_and_dot(const Number                     a,
                  const VectorSpaceVector<Number> &V,
                  const VectorSpaceVector<Number> &W) override;

      /**
       * 返回向量的全局大小，等于所有处理器中本地拥有的索引数之和。
       *
       */
      virtual size_type
      size() const override;

      /**
       * 返回一个索引集，描述这个向量的哪些元素为当前处理器所拥有。因此，如果这是一个分布式向量，在不同处理器上返回的索引集将形成不相交的集合，加起来就是完整的索引集。很明显，如果一个向量只在一个处理器上创建，那么结果将满足
       * @code
       * vec.locally_owned_elements() == complete_index_set(vec.size())
       * @endcode
       *
       *
       */
      virtual dealii::IndexSet
      locally_owned_elements() const override;

      /**
       * 将向量打印到输出流  @p out.  。
       *
       */
      virtual void
      print(std::ostream &     out,
            const unsigned int precision  = 3,
            const bool         scientific = true,
            const bool         across     = true) const override;

      /**
       * 返回这个类的内存消耗，单位是字节。
       *
       */
      virtual std::size_t
      memory_consumption() const override;
      //@}

      /**
       * @name  4：不包括在VectorSpaceVector中的其他向量操作
       *
       */
      //@{

      /**
       * 将向量的所有元素设置为标量  @p s.
       * 如果标量为零，也将鬼魂元素设置为零，否则它们保持不变。
       *
       */
      virtual Vector<Number, MemorySpace> &
      operator=(const Number s) override;

      /**
       * 这是一个集体添加操作，将存储在 @p values
       * 中的整组值添加到 @p indices. 指定的向量成分中。
       *
       */
      template <typename OtherNumber>
      void
      add(const std::vector<size_type> &       indices,
          const ::dealii::Vector<OtherNumber> &values);

      /**
       * 取一个连续存储n_elements的地址，将它们添加到向量中。
       *
       */
      template <typename OtherNumber>
      void
      add(const size_type    n_elements,
          const size_type *  indices,
          const OtherNumber *values);

      /**
       * 缩放和简单的向量添加，即<tt>*this = s*(*this)+V</tt>。
       *
       */
      void
      sadd(const Number s, const Vector<Number, MemorySpace> &V);

      //@}


      /**
       * @name  5：条目访问和本地数据表示法
       *
       */
      //@{

      /**
       * 返回向量的本地大小，即本地拥有的索引数。
       * @deprecated  使用local_owned_size()代替。
       *
       */
      DEAL_II_DEPRECATED
      size_type
      local_size() const;

      /**
       * 返回向量的本地大小，即本地拥有的索引数。
       *
       */
      size_type
      locally_owned_size() const;

      /**
       * 如果给定的全局索引在这个处理器的本地范围内，则返回true。
       *
       */
      bool
      in_local_range(const size_type global_index) const;

      /**
       * 使 @p Vector
       * 类有点像C++标准库中的<tt>vector<>/tt>类，返回该向量的<i>locally
       * owned</i>元素的开始和结束的迭代器。
       * 它认为，end()
       *
       * - begin() == local_owned_size()。
       * @note 对于CUDA内存空间，迭代器指向设备上的内存。
       *
       */
      iterator
      begin();

      /**
       * 返回常数迭代器，指向向量中本地拥有的元素的起点。
       * @note
       * 对于CUDA内存空间，该迭代器指向设备上的内存。
       *
       */
      const_iterator
      begin() const;

      /**
       * 返回一个迭代器，指向本地拥有的条目数组结束后的元素。
       * @note
       * 对于CUDA内存空间，该迭代器指向设备上的内存。
       *
       */
      iterator
      end();

      /**
       * 返回一个恒定的迭代器，指向本地拥有的条目数组末端之后的元素。
       * @note
       * 对于CUDA内存空间，该迭代器指向设备上的内存。
       *
       */
      const_iterator
      end() const;

      /**
       * 对 @p
       * global_index对应位置的数据进行读取访问。该索引必须在向量的本地范围内，或者在构造时被指定为一个幽灵索引。
       * 性能。<tt>O(1)</tt>对于代表连续范围的本地拥有的元素，<tt>O(log(n<sub>ranges</sub>))</tt>对于幽灵元素（相当快，但是比local_element()慢）。
       *
       */
      Number
      operator()(const size_type global_index) const;

      /**
       * 读取和写入对应于 @p
       * global_index位置的数据。该索引必须在向量的本地范围内，或者在构造时被指定为一个幽灵索引。
       * 性能。<tt>O(1)</tt>对于代表连续范围的本地拥有的元素，<tt>O(log(n<sub>ranges</sub>))</tt>对于幽灵元素（相当快，但是比local_element()慢）。
       *
       */
      Number &
      operator()(const size_type global_index);

      /**
       * 读取访问对应于 @p
       * global_index位置的数据。该索引必须在向量的本地范围内，或者在构造时被指定为一个鬼魂索引。
       * 这个函数与operator()的作用相同。
       *
       */
      Number operator[](const size_type global_index) const;
      /**
       * 读取和写入对应于 @p
       * global_index位置的数据。该索引必须在向量的本地范围内，或者在构造时被指定为一个幽灵索引。
       * 这个函数与operator()的作用相同。
       *
       */
      Number &operator[](const size_type global_index);

      /**
       * 读取访问由 @p local_index. 指定的数据字段
       * 本地拥有的指数可以用指数
       * <code>[0,locally_owned_size)</code>
       * 访问，而鬼魂指数可以用指数
       * <code>[locally_owned_size,locally_owned_size+ n_ghost_entries]</code>
       * 访问。            性能。直接数组访问（快速）。
       *
       */
      Number
      local_element(const size_type local_index) const;

      /**
       * 对 @p local_index. 指定的数据字段进行读写访问
       * 可以用索引 <code>[0,locally_owned_size())</code>
       * 访问本地拥有的索引，用索引 <code>[locally_owned_size(),
       * locally_owned_size()+n_ghosts]</code> 访问幽灵索引。
       * 性能。直接数组访问（快速）。
       *
       */
      Number &
      local_element(const size_type local_index);

      /**
       * 返回指向底层原始数组的指针。
       * @note  对于CUDA内存空间，指针指向设备上的内存。
       *
       */
      Number *
      get_values() const;

      /**
       * 与通过operator()获取向量的单个元素不同，这个函数允许一次性获取一整组元素。要读取的元素的索引在第一个参数中说明，相应的值在第二个参数中返回。
       * 如果当前的向量被称为 @p v,
       * ，那么这个函数就等同于代码
       * @code
       * for (unsigned int i=0; i<indices.size(); ++i)
       *   values[i] = v[indices[i]];
       * @endcode
       * @pre   @p indices 和 @p values 数组的大小必须是一致的。
       * @note  这个函数在CUDA内存空间没有实现。
       *
       */
      template <typename OtherNumber>
      void
      extract_subvector_to(const std::vector<size_type> &indices,
                           std::vector<OtherNumber> &    values) const;

      /**
       * 这个函数不是通过operator()来获取向量的单个元素，而是允许一次性获取一整组元素。与前一个函数不同的是，这个函数通过取消引用前两个参数提供的迭代器范围内的所有元素来获得元素的索引，并将向量的值放入通过取消引用从第三个参数指向的位置开始的迭代器范围获得的内存位置。
       * 如果当前的向量被称为 @p v,
       * ，那么这个函数就相当于代码中的
       * @code
       * ForwardIterator indices_p = indices_begin;
       * OutputIterator  values_p  = values_begin;
       * while (indices_p != indices_end)
       * {
       *  values_p = v[*indices_p];
       *   ++indices_p;
       *   ++values_p;
       * }
       * @endcode
       * @pre  必须能够写进从 @p values_begin
       * 开始的尽可能多的内存位置，因为有 @p indices_begin 和
       * @p indices_end. 之间的迭代器。
       *
       */
      template <typename ForwardIterator, typename OutputIterator>
      void
      extract_subvector_to(ForwardIterator       indices_begin,
                           const ForwardIterator indices_end,
                           OutputIterator        values_begin) const;
      /**
       * 返回该向量是否只包含值为0的元素。
       * 这是一个集体操作。这个函数很昂贵，因为有可能所有的元素都要被检查。
       *
       */
      virtual bool
      all_zero() const override;

      /**
       * 计算向量中所有条目的平均值。
       *
       */
      virtual Number
      mean_value() const override;

      /**
       * $l_p$
       * --向量的正负值。元素绝对值的p次方之和的p次根。
       *
       */
      real_type
      lp_norm(const real_type p) const;
      //@}

      /**
       * @name  6：混合东西
       *
       */
      //@{

      /**
       * 返回一个对与此向量一起使用的MPI通信器对象的引用。
       *
       */
      const MPI_Comm &
      get_mpi_communicator() const;

      /**
       * 返回描述向量的并行布局的MPI分区器。这个对象可以用来用各自的reinit()调用初始化另一个向量，用于有关并行通信的额外查询，或者分区器的兼容性。
       *
       */
      const std::shared_ptr<const Utilities::MPI::Partitioner> &
      get_partitioner() const;

      /**
       * 检查给定的分区器是否与用于该向量的分区器兼容。如果两个分区器具有相同的局部大小和相同的ghost索引，那么它们就是兼容的。它们不一定需要是共享指针的同一个数据域。
       * 这只是一个局部操作，即如果只有一些处理器决定分区不兼容，只有这些处理器会返回
       * @p false, ，而其他处理器会返回 @p true.  。
       *
       */
      bool
      partitioners_are_compatible(
        const Utilities::MPI::Partitioner &part) const;

      /**
       * 检查给定的分区器是否与该向量使用的分区器兼容。如果两个分区器具有相同的局部大小和相同的ghost索引，那么它们就是兼容的。它们不一定需要是同一个数据域。与partitioners_are_compatible()相反，该方法检查所有处理器之间的兼容性，如果分区器在所有处理器上是相同的，该方法只返回
       * @p true 。
       * 这个方法执行全局通信，所以要确保只在所有处理器调用它的次数相同的情况下使用它。
       *
       */
      bool
      partitioners_are_globally_compatible(
        const Utilities::MPI::Partitioner &part) const;

      /**
       * 将这个向量的幽灵状态改为 @p ghosted. 。
       *
       */
      void
      set_ghost_state(const bool ghosted) const;

      /**
       * 获取指向同一共享内存域的其他进程的起始值的指针。
       *
       */
      const std::vector<ArrayView<const Number>> &
      shared_vector_data() const;

      //@}

      /**
       * 试图在两个不兼容的向量类型之间执行操作。
       * @ingroup Exceptions
       *
       */
      DeclException0(ExcVectorTypeNotCompatible);

      /**
       * 试图执行一个在设备上没有实现的操作。
       * @ingroup Exceptions
       *
       */
      DeclException0(ExcNotAllowedForCuda);

      /**
       * 异常情况
       *
       */
      DeclException3(ExcNonMatchingElements,
                     Number,
                     Number,
                     unsigned int,
                     << "Called compress(VectorOperation::insert), but"
                     << " the element received from a remote processor, value "
                     << std::setprecision(16) << arg1
                     << ", does not match with the value "
                     << std::setprecision(16) << arg2
                     << " on the owner processor " << arg3);

      /**
       * 异常情况
       *
       */
      DeclException4(
        ExcAccessToNonLocalElement,
        size_type,
        size_type,
        size_type,
        size_type,
        << "You tried to access element " << arg1
        << " of a distributed vector, but this element is not "
        << "stored on the current processor. Note: The range of "
        << "locally owned elements is [" << arg2 << "," << arg3
        << "], and there are " << arg4 << " ghost elements "
        << "that this vector can access."
        << "\n\n"
        << "A common source for this kind of problem is that you "
        << "are passing a 'fully distributed' vector into a function "
        << "that needs read access to vector elements that correspond "
        << "to degrees of freedom on ghost cells (or at least to "
        << "'locally active' degrees of freedom that are not also "
        << "'locally owned'). You need to pass a vector that has these "
        << "elements as ghost entries.");

    private:
      /**
       * 简单的向量的倍数加法，即<tt>*this +=
       * a*V</tt>，无需MPI通信。
       *
       */
      void
      add_local(const Number a, const VectorSpaceVector<Number> &V);

      /**
       * 一个向量的倍数的缩放和简单加法，即：<tt>*this =
       * s*(*this)+a*V</tt>，无需MPI通信。
       *
       */
      void
      sadd_local(const Number                     s,
                 const Number                     a,
                 const VectorSpaceVector<Number> &V);

      /**
       * 两个向量的内积的局部部分。
       *
       */
      template <typename Number2>
      Number
      inner_product_local(const Vector<Number2, MemorySpace> &V) const;

      /**
       * norm_sqr()的本地部分。
       *
       */
      real_type
      norm_sqr_local() const;

      /**
       * Mean_value()的局部部分。
       *
       */
      Number
      mean_value_local() const;

      /**
       * l1_norm()的本地部分。
       *
       */
      real_type
      l1_norm_local() const;

      /**
       * lp_norm()的本地部分。
       *
       */
      real_type
      lp_norm_local(const real_type p) const;

      /**
       * linfty_norm()的本地部分。
       *
       */
      real_type
      linfty_norm_local() const;

      /**
       * 两个向量的内积之后的加法的局部部分。这同样适用于复值向量的add_and_dot()函数。
       *
       */
      Number
      add_and_dot_local(const Number                       a,
                        const Vector<Number, MemorySpace> &V,
                        const Vector<Number, MemorySpace> &W);

      /**
       * 共享指针，用于存储并行分区的信息。这个信息可以在具有相同分区的几个向量之间共享。
       *
       */
      std::shared_ptr<const Utilities::MPI::Partitioner> partitioner;

      /**
       * 目前分配在val数组中的大小。
       *
       */
      size_type allocated_size;

      /**
       * 储存该向量本地元素的底层数据结构。
       *
       */
      mutable ::dealii::MemorySpace::MemorySpaceData<Number, MemorySpace> data;

      /**
       * 对于带有TBB的并行循环，这个成员变量存储循环的亲和力信息。
       *
       */
      mutable std::shared_ptr<::dealii::parallel::internal::TBBPartitioner>
        thread_loop_partitioner;

      /**
       * 临时存储，用于保存在compress()中发送到这个处理器的数据，或者在update_ghost_values()中从这个处理器发送的数据。
       *
       */
      mutable ::dealii::MemorySpace::MemorySpaceData<Number, MemorySpace>
        import_data;

      /**
       * 存储向量当前是否允许读取ghost元素。请注意，这是为了确保一致的鬼魂数据，并不表明向量是否真的可以存储鬼魂元素。特别是在组装一个向量时，我们不允许读取元素，只允许写入元素。
       *
       */
      mutable bool vector_is_ghosted;

#ifdef DEAL_II_WITH_MPI
      /**
       * 一个收集所有来自compress()操作的请求的向量。
       * 这个类使用持久的MPI通信器，也就是说，在连续调用一个给定的函数时，通信通道被存储起来。这减少了设置MPI机器所涉及的开销，但它并没有消除在数据真正被发送之前发布接收操作的需要。
       *
       */
      std::vector<MPI_Request> compress_requests;

      /**
       * 一个向量，收集来自update_ghost_values()操作的所有请求。这个类使用持久的MPI通信器。
       *
       */
      mutable std::vector<MPI_Request> update_ghost_values_requests;
#endif

      /**
       * 一个锁，确保压缩()和update_ghost_values()函数在与多个线程一起使用时也能给出合理的结果。
       *
       */
      mutable std::mutex mutex;

      /**
       * 用于共享内存领域的通信器。关于 "comm_sm
       * "的用途的更多信息，请参见该类的一般文档。
       *
       */
      MPI_Comm comm_sm;

      /**
       * 一个帮助函数，用于清除compress_requests和update_ghost_values_requests域。在reinit()函数中使用。
       *
       */
      void
      clear_mpi_requests();

      /**
       * 一个帮助函数，用于调整val数组的大小。
       *
       */
      void
      resize_val(const size_type new_allocated_size,
                 const MPI_Comm &comm_sm = MPI_COMM_SELF);

      // Make all other vector types friends.
      template <typename Number2, typename MemorySpace2>
      friend class Vector;

      // Make BlockVector type friends.
      template <typename Number2>
      friend class BlockVector;
    };
     /*@}*/ 


     /*-------------------- Inline functions ---------------------------------*/ 

#ifndef DOXYGEN

    namespace internal
    {
      template <typename Number, typename MemorySpace>
      struct Policy
      {
        static inline typename Vector<Number, MemorySpace>::iterator
        begin(::dealii::MemorySpace::MemorySpaceData<Number, MemorySpace> &)
        {
          return nullptr;
        }

        static inline typename Vector<Number, MemorySpace>::const_iterator
        begin(
          const ::dealii::MemorySpace::MemorySpaceData<Number, MemorySpace> &)
        {
          return nullptr;
        }

        static inline Number *
        get_values(
          ::dealii::MemorySpace::MemorySpaceData<Number, MemorySpace> &)
        {
          return nullptr;
        }
      };



      template <typename Number>
      struct Policy<Number, ::dealii::MemorySpace::Host>
      {
        static inline
          typename Vector<Number, ::dealii::MemorySpace::Host>::iterator
          begin(::dealii::MemorySpace::
                  MemorySpaceData<Number, ::dealii::MemorySpace::Host> &data)
        {
          return data.values.get();
        }

        static inline
          typename Vector<Number, ::dealii::MemorySpace::Host>::const_iterator
          begin(const ::dealii::MemorySpace::
                  MemorySpaceData<Number, ::dealii::MemorySpace::Host> &data)
        {
          return data.values.get();
        }

        static inline Number *
        get_values(::dealii::MemorySpace::
                     MemorySpaceData<Number, ::dealii::MemorySpace::Host> &data)
        {
          return data.values.get();
        }
      };



      template <typename Number>
      struct Policy<Number, ::dealii::MemorySpace::CUDA>
      {
        static inline
          typename Vector<Number, ::dealii::MemorySpace::CUDA>::iterator
          begin(::dealii::MemorySpace::
                  MemorySpaceData<Number, ::dealii::MemorySpace::CUDA> &data)
        {
          return data.values_dev.get();
        }

        static inline
          typename Vector<Number, ::dealii::MemorySpace::CUDA>::const_iterator
          begin(const ::dealii::MemorySpace::
                  MemorySpaceData<Number, ::dealii::MemorySpace::CUDA> &data)
        {
          return data.values_dev.get();
        }

        static inline Number *
        get_values(::dealii::MemorySpace::
                     MemorySpaceData<Number, ::dealii::MemorySpace::CUDA> &data)
        {
          return data.values_dev.get();
        }
      };
    } // namespace internal


    template <typename Number, typename MemorySpace>
    inline bool
    Vector<Number, MemorySpace>::has_ghost_elements() const
    {
      return vector_is_ghosted;
    }



    template <typename Number, typename MemorySpace>
    inline typename Vector<Number, MemorySpace>::size_type
    Vector<Number, MemorySpace>::size() const
    {
      return partitioner->size();
    }



    template <typename Number, typename MemorySpace>
    inline typename Vector<Number, MemorySpace>::size_type
    Vector<Number, MemorySpace>::local_size() const
    {
      return locally_owned_size();
    }



    template <typename Number, typename MemorySpace>
    inline typename Vector<Number, MemorySpace>::size_type
    Vector<Number, MemorySpace>::locally_owned_size() const
    {
      return partitioner->locally_owned_size();
    }



    template <typename Number, typename MemorySpace>
    inline bool
    Vector<Number, MemorySpace>::in_local_range(
      const size_type global_index) const
    {
      return partitioner->in_local_range(global_index);
    }



    template <typename Number, typename MemorySpace>
    inline IndexSet
    Vector<Number, MemorySpace>::locally_owned_elements() const
    {
      IndexSet is(size());

      is.add_range(partitioner->local_range().first,
                   partitioner->local_range().second);

      return is;
    }



    template <typename Number, typename MemorySpace>
    inline typename Vector<Number, MemorySpace>::iterator
    Vector<Number, MemorySpace>::begin()
    {
      return internal::Policy<Number, MemorySpace>::begin(data);
    }



    template <typename Number, typename MemorySpace>
    inline typename Vector<Number, MemorySpace>::const_iterator
    Vector<Number, MemorySpace>::begin() const
    {
      return internal::Policy<Number, MemorySpace>::begin(data);
    }



    template <typename Number, typename MemorySpace>
    inline typename Vector<Number, MemorySpace>::iterator
    Vector<Number, MemorySpace>::end()
    {
      return internal::Policy<Number, MemorySpace>::begin(data) +
             partitioner->locally_owned_size();
    }



    template <typename Number, typename MemorySpace>
    inline typename Vector<Number, MemorySpace>::const_iterator
    Vector<Number, MemorySpace>::end() const
    {
      return internal::Policy<Number, MemorySpace>::begin(data) +
             partitioner->locally_owned_size();
    }



    template <typename Number, typename MemorySpace>
    const std::vector<ArrayView<const Number>> &
    Vector<Number, MemorySpace>::shared_vector_data() const
    {
      return data.values_sm;
    }



    template <typename Number, typename MemorySpace>
    inline Number
    Vector<Number, MemorySpace>::operator()(const size_type global_index) const
    {
      Assert((std::is_same<MemorySpace, ::dealii::MemorySpace::Host>::value),
             ExcMessage(
               "This function is only implemented for the Host memory space"));
      Assert(
        partitioner->in_local_range(global_index) ||
          partitioner->ghost_indices().is_element(global_index),
        ExcAccessToNonLocalElement(global_index,
                                   partitioner->local_range().first,
                                   partitioner->local_range().second - 1,
                                   partitioner->ghost_indices().n_elements()));
      // do not allow reading a vector which is not in ghost mode
      Assert(partitioner->in_local_range(global_index) ||
               vector_is_ghosted == true,
             ExcMessage("You tried to read a ghost element of this vector, "
                        "but it has not imported its ghost values."));
      return data.values[partitioner->global_to_local(global_index)];
    }



    template <typename Number, typename MemorySpace>
    inline Number &
    Vector<Number, MemorySpace>::operator()(const size_type global_index)
    {
      Assert((std::is_same<MemorySpace, ::dealii::MemorySpace::Host>::value),
             ExcMessage(
               "This function is only implemented for the Host memory space"));
      Assert(
        partitioner->in_local_range(global_index) ||
          partitioner->ghost_indices().is_element(global_index),
        ExcAccessToNonLocalElement(global_index,
                                   partitioner->local_range().first,
                                   partitioner->local_range().second - 1,
                                   partitioner->ghost_indices().n_elements()));
      // we would like to prevent reading ghosts from a vector that does not
      // have them imported, but this is not possible because we might be in a
      // part of the code where the vector has enabled ghosts but is non-const
      // (then, the compiler picks this method according to the C++ rule book
      // even if a human would pick the const method when this subsequent use
      // is just a read)
      return data.values[partitioner->global_to_local(global_index)];
    }



    template <typename Number, typename MemorySpace>
    inline Number Vector<Number, MemorySpace>::
                  operator[](const size_type global_index) const
    {
      return operator()(global_index);
    }



    template <typename Number, typename MemorySpace>
    inline Number &Vector<Number, MemorySpace>::
                   operator[](const size_type global_index)
    {
      return operator()(global_index);
    }



    template <typename Number, typename MemorySpace>
    inline Number
    Vector<Number, MemorySpace>::local_element(
      const size_type local_index) const
    {
      Assert((std::is_same<MemorySpace, ::dealii::MemorySpace::Host>::value),
             ExcMessage(
               "This function is only implemented for the Host memory space"));
      AssertIndexRange(local_index,
                       partitioner->locally_owned_size() +
                         partitioner->n_ghost_indices());
      // do not allow reading a vector which is not in ghost mode
      Assert(local_index < locally_owned_size() || vector_is_ghosted == true,
             ExcMessage("You tried to read a ghost element of this vector, "
                        "but it has not imported its ghost values."));

      return data.values[local_index];
    }



    template <typename Number, typename MemorySpace>
    inline Number &
    Vector<Number, MemorySpace>::local_element(const size_type local_index)
    {
      Assert((std::is_same<MemorySpace, ::dealii::MemorySpace::Host>::value),
             ExcMessage(
               "This function is only implemented for the Host memory space"));

      AssertIndexRange(local_index,
                       partitioner->locally_owned_size() +
                         partitioner->n_ghost_indices());

      return data.values[local_index];
    }



    template <typename Number, typename MemorySpace>
    inline Number *
    Vector<Number, MemorySpace>::get_values() const
    {
      return internal::Policy<Number, MemorySpace>::get_values(data);
    }



    template <typename Number, typename MemorySpace>
    template <typename OtherNumber>
    inline void
    Vector<Number, MemorySpace>::extract_subvector_to(
      const std::vector<size_type> &indices,
      std::vector<OtherNumber> &    values) const
    {
      for (size_type i = 0; i < indices.size(); ++i)
        values[i] = operator()(indices[i]);
    }



    template <typename Number, typename MemorySpace>
    template <typename ForwardIterator, typename OutputIterator>
    inline void
    Vector<Number, MemorySpace>::extract_subvector_to(
      ForwardIterator       indices_begin,
      const ForwardIterator indices_end,
      OutputIterator        values_begin) const
    {
      while (indices_begin != indices_end)
        {
          *values_begin = operator()(*indices_begin);
          indices_begin++;
          values_begin++;
        }
    }



    template <typename Number, typename MemorySpace>
    template <typename OtherNumber>
    inline void
    Vector<Number, MemorySpace>::add(
      const std::vector<size_type> &       indices,
      const ::dealii::Vector<OtherNumber> &values)
    {
      AssertDimension(indices.size(), values.size());
      for (size_type i = 0; i < indices.size(); ++i)
        {
          Assert(
            numbers::is_finite(values[i]),
            ExcMessage(
              "The given value is not finite but either infinite or Not A Number (NaN)"));
          this->operator()(indices[i]) += values(i);
        }
    }



    template <typename Number, typename MemorySpace>
    template <typename OtherNumber>
    inline void
    Vector<Number, MemorySpace>::add(const size_type    n_elements,
                                     const size_type *  indices,
                                     const OtherNumber *values)
    {
      for (size_type i = 0; i < n_elements; ++i, ++indices, ++values)
        {
          Assert(
            numbers::is_finite(*values),
            ExcMessage(
              "The given value is not finite but either infinite or Not A Number (NaN)"));
          this->operator()(*indices) += *values;
        }
    }



    template <typename Number, typename MemorySpace>
    inline const MPI_Comm &
    Vector<Number, MemorySpace>::get_mpi_communicator() const
    {
      return partitioner->get_mpi_communicator();
    }



    template <typename Number, typename MemorySpace>
    inline const std::shared_ptr<const Utilities::MPI::Partitioner> &
    Vector<Number, MemorySpace>::get_partitioner() const
    {
      return partitioner;
    }



    template <typename Number, typename MemorySpace>
    inline void
    Vector<Number, MemorySpace>::set_ghost_state(const bool ghosted) const
    {
      vector_is_ghosted = ghosted;
    }

#endif

  } // namespace distributed
} // namespace LinearAlgebra


/**
 * 全局函数  @p swap
 * ，它重载了C++标准库的默认实现，它使用一个临时对象。该函数简单地交换了两个向量的数据。
 * @relatesalso  Vectors
 *
 *
 */
template <typename Number, typename MemorySpace>
inline void
swap(LinearAlgebra::distributed::Vector<Number, MemorySpace> &u,
     LinearAlgebra::distributed::Vector<Number, MemorySpace> &v)
{
  u.swap(v);
}


/**
 * 将 dealii::LinearAlgebra::Vector 声明为分布式向量。
 *
 *
 */
template <typename Number, typename MemorySpace>
struct is_serial_vector<LinearAlgebra::distributed::Vector<Number, MemorySpace>>
  : std::false_type
{};



namespace internal
{
  namespace LinearOperatorImplementation
  {
    template <typename>
    class ReinitHelper;

    /**
     * linear_operator.h中内部使用的一个辅助类。对
     * LinearAlgebra::distributed::Vector. 进行专业化处理。
     *
     */
    template <typename Number>
    class ReinitHelper<LinearAlgebra::distributed::Vector<Number>>
    {
    public:
      // A helper type-trait that leverage SFINAE to figure out if type T has
      // void T::get_mpi_communicator()
      template <typename T>
      struct has_get_mpi_communicator
      {
      private:
        static bool
        detect(...);

        template <typename U>
        static decltype(std::declval<U>().get_mpi_communicator())
        detect(const U &);

      public:
        static const bool value =
          !std::is_same<bool, decltype(detect(std::declval<T>()))>::value;
      };

      // A helper type-trait that leverage SFINAE to figure out if type T has
      // void T::locally_owned_domain_indices()
      template <typename T>
      struct has_locally_owned_domain_indices
      {
      private:
        static bool
        detect(...);

        template <typename U>
        static decltype(std::declval<U>().locally_owned_domain_indices())
        detect(const U &);

      public:
        static const bool value =
          !std::is_same<bool, decltype(detect(std::declval<T>()))>::value;
      };

      // A helper type-trait that leverage SFINAE to figure out if type T has
      // void T::locally_owned_range_indices()
      template <typename T>
      struct has_locally_owned_range_indices
      {
      private:
        static bool
        detect(...);

        template <typename U>
        static decltype(std::declval<U>().locally_owned_range_indices())
        detect(const U &);

      public:
        static const bool value =
          !std::is_same<bool, decltype(detect(std::declval<T>()))>::value;
      };

      // A helper type-trait that leverage SFINAE to figure out if type T has
      // void T::initialize_dof_vector(VectorType v)
      template <typename T>
      struct has_initialize_dof_vector
      {
      private:
        static bool
        detect(...);

        template <typename U>
        static decltype(std::declval<U>().initialize_dof_vector(
          std::declval<LinearAlgebra::distributed::Vector<Number> &>()))
        detect(const U &);

      public:
        static const bool value =
          !std::is_same<bool, decltype(detect(std::declval<T>()))>::value;
      };

      // Used for (Trilinos/PETSc)Wrappers::SparseMatrix
      template <typename MatrixType,
                typename std::enable_if<
                  has_get_mpi_communicator<MatrixType>::value &&
                    has_locally_owned_domain_indices<MatrixType>::value,
                  MatrixType>::type * = nullptr>
      static void
      reinit_domain_vector(MatrixType &                                mat,
                           LinearAlgebra::distributed::Vector<Number> &vec,
                           bool  /*omit_zeroing_entries*/ )
      {
        vec.reinit(mat.locally_owned_domain_indices(),
                   mat.get_mpi_communicator());
      }

      // Used for MatrixFree and DiagonalMatrix
      template <
        typename MatrixType,
        typename std::enable_if<has_initialize_dof_vector<MatrixType>::value,
                                MatrixType>::type * = nullptr>
      static void
      reinit_domain_vector(MatrixType &                                mat,
                           LinearAlgebra::distributed::Vector<Number> &vec,
                           bool omit_zeroing_entries)
      {
        mat.initialize_dof_vector(vec);
        if (!omit_zeroing_entries)
          vec = Number();
      }

      // Used for (Trilinos/PETSc)Wrappers::SparseMatrix
      template <typename MatrixType,
                typename std::enable_if<
                  has_get_mpi_communicator<MatrixType>::value &&
                    has_locally_owned_range_indices<MatrixType>::value,
                  MatrixType>::type * = nullptr>
      static void
      reinit_range_vector(MatrixType &                                mat,
                          LinearAlgebra::distributed::Vector<Number> &vec,
                          bool  /*omit_zeroing_entries*/ )
      {
        vec.reinit(mat.locally_owned_range_indices(),
                   mat.get_mpi_communicator());
      }

      // Used for MatrixFree and DiagonalMatrix
      template <
        typename MatrixType,
        typename std::enable_if<has_initialize_dof_vector<MatrixType>::value,
                                MatrixType>::type * = nullptr>
      static void
      reinit_range_vector(MatrixType &                                mat,
                          LinearAlgebra::distributed::Vector<Number> &vec,
                          bool omit_zeroing_entries)
      {
        mat.initialize_dof_vector(vec);
        if (!omit_zeroing_entries)
          vec = Number();
      }
    };

  } // namespace LinearOperatorImplementation
}  /* namespace internal */ 


DEAL_II_NAMESPACE_CLOSE

#endif


