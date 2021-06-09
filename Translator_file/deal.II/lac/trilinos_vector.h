//include/deal.II-translator/lac/trilinos_vector_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2021 by the deal.II authors
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

#ifndef dealii_trilinos_vector_h
#define dealii_trilinos_vector_h


#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_TRILINOS
#  include <deal.II/base/index_set.h>
#  include <deal.II/base/mpi.h>
#  include <deal.II/base/subscriptor.h>
#  include <deal.II/base/utilities.h>

#  include <deal.II/lac/exceptions.h>
#  include <deal.II/lac/vector.h>
#  include <deal.II/lac/vector_operation.h>
#  include <deal.II/lac/vector_type_traits.h>

#  include <Epetra_ConfigDefs.h>

#  include <memory>
#  include <utility>
#  include <vector>
#  ifdef DEAL_II_WITH_MPI // only if MPI is installed
#    include <Epetra_MpiComm.h>
#    include <mpi.h>
#  else
#    include <Epetra_SerialComm.h>
#  endif
#  include <Epetra_FEVector.h>
#  include <Epetra_LocalMap.h>
#  include <Epetra_Map.h>

DEAL_II_NAMESPACE_OPEN

// Forward declarations
#  ifndef DOXYGEN
namespace LinearAlgebra
{
  // Forward declaration
  template <typename Number>
  class ReadWriteVector;
} // namespace LinearAlgebra
#  endif

/**
 * @addtogroup  TrilinosWrappers  
     * @{ 
 *
 */

/**
 * 一个命名空间，Trilinos对象的封装类位于其中。
 *
 *
 * @ingroup TrilinosWrappers
 *
 *
 */
namespace TrilinosWrappers
{
  class SparseMatrix;

  /**
   * @cond 内部
   *
   */

  /**
   * 一个命名空间，用于TrilinosWrapper成员的内部实现细节。
   * @ingroup TrilinosWrappers
   *
   */
  namespace internal
  {
    /**
     * 声明容器大小的类型。
     *
     */
    using size_type = dealii::types::global_dof_index;

    /**
     * 这个类实现了一个用于访问Trilinos向量的包装器，就像我们访问deal.II对象一样：它被初始化为一个向量和其中的一个元素，并且有一个转换操作符来提取这个元素的标量值。它也有各种赋值运算符，用于向这一个元素写入。
     * @ingroup TrilinosWrappers
     *
     */
    class VectorReference
    {
    private:
      /**
       * 构造函数。它是私有的，以便只允许实际的向量类来创建它。
       *
       */
      VectorReference(MPI::Vector &vector, const size_type index);

    public:
      /**
       * 复制构造函数。
       *
       */
      VectorReference(const VectorReference &) = default;

      /**
       * 这看起来像一个复制操作符，但做的事情与平常不同。特别是，它并不复制这个引用的成员变量。相反，它处理的情况是，我们有两个向量
       * @p v 和 @p w,
       * ，像<tt>v(i)=w(i)</tt>中那样分配元素。这里，赋值的左手和右手都有数据类型VectorReference，但我们真正的意思是赋值这两个引用所代表的向量元素。这个操作符实现了这个操作。还要注意的是，这使得我们可以使赋值运算符成为常数。
       *
       */
      const VectorReference &
      operator=(const VectorReference &r) const;

      /**
       * 和上面一样，但对于非const的引用对象。
       *
       */
      VectorReference &
      operator=(const VectorReference &r);

      /**
       * 将向量的引用元素设置为<tt>s</tt>。
       *
       */
      const VectorReference &
      operator=(const TrilinosScalar &s) const;

      /**
       * 将<tt>s</tt>添加到矢量的引用元素中->。
       *
       */
      const VectorReference &
      operator+=(const TrilinosScalar &s) const;

      /**
       * 从向量的参考元素中减去<tt>s</tt>->。
       *
       */
      const VectorReference &
      operator-=(const TrilinosScalar &s) const;

      /**
       * 用<tt>s</tt>乘以矢量中的参考元素。
       *
       */
      const VectorReference &
      operator*=(const TrilinosScalar &s) const;

      /**
       * 将向量的参考元素除以<tt>s</tt>。
       *
       */
      const VectorReference &
      operator/=(const TrilinosScalar &s) const;

      /**
       * 将引用转换为实际值，即返回矢量中被引用元素的值。
       *
       */
      operator TrilinosScalar() const;

      /**
       * 异常情况
       *
       */
      DeclException1(ExcTrilinosError,
                     int,
                     << "An error with error number " << arg1
                     << " occurred while calling a Trilinos function");

    private:
      /**
       * 指向我们所引用的向量。
       *
       */
      MPI::Vector &vector;

      /**
       * 向量中被引用的元素的索引。
       *
       */
      const size_type index;

      // Make the vector class a friend, so that it can create objects of the
      // present type.
      friend class ::dealii::TrilinosWrappers::MPI::Vector;
    };
  } // namespace internal
  /**
   * @endcond
   *
   */

#  ifndef DEAL_II_WITH_64BIT_INDICES
    // define a helper function that queries the global ID of local ID of
  // an Epetra_BlockMap object  by calling either the 32- or 64-bit
  // function necessary.
  inline int
  gid(const Epetra_BlockMap &map, int i)
  {
    return map.GID(i);
  }
#  else
    // define a helper function that queries the global ID of local ID of
  // an Epetra_BlockMap object  by calling either the 32- or 64-bit
  // function necessary.
  inline long long int
  gid(const Epetra_BlockMap &map, int i)
  {
    return map.GID64(i);
  }
#  endif

  /**
   * 用于通过MPI并行工作的Trilinos向量类的命名空间。
   * @ingroup TrilinosWrappers
   *
   */
  namespace MPI
  {
    class BlockVector;

    /**
     * 该类实现了一个包装器，用于使用Trilinos分布式向量类Epetra_FEVector，其（并行）分区由Epetra_Map管理。
     * Epetra_FEVector正是我们一直在处理的那种矢量。
     *
     * - 我们可能从一些装配过程中得到它，在那里，不属于本地的条目也可能需要被写入，因此需要转发给所有者。        这个类的接口是以deal.II中现有的Vector类为模型的。它有几乎相同的成员函数，并且通常是可交换的。然而，由于Trilinos只支持单一的标量类型（double），所以它不是模板化的，只对这种类型起作用。        请注意，只有在向量装配后调用了函数 @p GlobalAssemble 以分配数据的情况下，Trilinos才会保证操作能达到你预期的效果。这是必要的，因为一些进程可能积累了不属于自己的元素数据，但必须被发送到拥有进程中。为了避免使用错误的数据，你需要在实际使用向量之前调用 Vector::compress() 。        <h3>Parallel communication model</h3> Trilinos的并行功能是建立在消息传递接口（MPI）之上的。MPI的通信模型建立在集体通信的基础上：如果一个进程想从另一个进程得到什么，另一个进程必须愿意接受这种通信。一个进程不能通过调用一个远程函数从另一个进程查询数据，而另一个进程不期望有这种交易。其结果是，这个类的基类中的大多数操作都必须被集体调用。例如，如果你想计算一个平行向量的l2准则， @em 所有共享这个向量的进程都必须调用 @p l2_norm 函数。如果你不这样做，而只是在一个进程上调用 @p l2_norm 函数，那么会发生以下情况。这一个进程将调用一个集体MPI函数，并等待所有其他进程加入到这个过程中。由于其他进程不调用这个函数，你将在第一个进程上得到一个超时，或者更糟糕的是，当下一个对Trilinos函数的调用在其他进程上产生MPI消息时，你将得到一个神秘的消息，即只有一个子集的进程试图进行通信。这些错误是很难搞清楚的，除非你对MPI的通信模型非常熟悉，并且知道哪些函数可能产生MPI消息。        下面将讨论一个特殊的情况，即可能会意外地产生MPI消息的情况。            <h3>Accessing individual elements of a vector</h3> Trilinos当然允许读取矢量的单个元素，但在分布式情况下，只允许读取本地存储的元素。我们通过<tt>d=vec(i)</tt>这样的调用来实现。    然而，如果你访问本地存储范围之外的元素，就会产生一个异常。        与读访问相反，Trilinos（以及相应的deal.II包装类）允许对向量的单个元素进行写入（或添加），即使它们存储在不同的进程中。你可以通过使用语法<tt>vec(i)=d</tt>或<tt>vec(i)+=d</tt>，或类似的操作写入或添加到元素中。然而，有一个问题，可能会导致非常混乱的错误信息。    Trilinos要求应用程序在从执行一组向元素添加的操作转换到执行一组向元素写入的操作时调用compress()函数。其理由是，所有进程都可能积累对元素的加法操作，即使是多个进程写到相同的元素。当我们下次调用compress()时，所有这些加法操作都被执行了。然而，如果一个进程对一个元素进行添加，而另一个进程对其进行覆盖，如果我们不确保在这中间发生与compress()的同步，那么执行的顺序将产生非确定性的行为。        为了确保这些对compress()的调用在适当的时间发生，deal.II包装器保留了一个状态变量，用于存储当前允许的操作：添加或写入。如果它遇到了相反的操作，它就会调用compress()并翻转状态。这有时会导致非常混乱的行为，例如，代码可能看起来像这样。
     * @code
     * TrilinosWrappers::MPI::Vector vector;
     * // do some write operations on the vector
     * for (size_type i=0; i<vector->size(); ++i)
     * vector(i) = i;
     *
     *                 // do some additions to vector elements, but
     *                 // only for some elements
     * for (size_type i=0; i<vector->size(); ++i)
     *   if (some_condition(i) == true)
     *     vector(i) += 1;
     *
     *                 // do another collective operation
     * const double norm = vector->l2_norm();
     * @endcode
     * 这段代码可能会遇到麻烦：当我们看到第一个加法运算时，我们需要冲刷向量的覆盖缓冲区，而deal.II库会通过调用compress()来实现。然而，它只对所有实际进行加法运算的进程进行这样的操作。
     *
     * @ingroup TrilinosWrappers
     * @ingroup Vectors  2008, 2009, 2017
     *
     */
    class Vector : public Subscriptor
    {
    public:
      /**
       * 声明一些在所有容器中使用的标准类型。这些类型与<tt>C</tt>标准库中的<tt>vector<...></tt>类中的类型平行。
       *
       */
      using value_type      = TrilinosScalar;
      using real_type       = TrilinosScalar;
      using size_type       = dealii::types::global_dof_index;
      using iterator        = value_type *;
      using const_iterator  = const value_type *;
      using reference       = internal::VectorReference;
      using const_reference = const internal::VectorReference;

      /**
       * @name  1: 基本对象处理
       *
       */
      //@{
      /**
       * 默认构造函数，生成一个空的（零大小）向量。在MPI运行的情况下，函数<tt>reinit()</tt>将必须给向量以正确的大小和在进程中的分布。
       *
       */
      Vector();

      /**
       * 使用给定向量的复制构造函数。
       *
       */
      Vector(const Vector &v);

      /**
       * 这个构造函数接受一个IndexSet，它定义了如何在MPI处理器之间分配各个组件。由于它还包括关于向量大小的信息，这就是我们生成一个%并行向量所需要的全部内容。            根据 @p parallel_partitioning 参数是否在处理器之间唯一地细分了元素，产生的向量可能有也可能没有鬼魂元素。更多信息请参见该类的一般文档。            如果提供的IndexSet形成了一个重叠的分区，就不清楚哪些元素被哪个进程拥有，local_owned_elements()将返回一个大小为0的IndexSet。              @see   @ref GlossGhostedVector  "有幽灵元素的向量"
       *
       */
      explicit Vector(const IndexSet &parallel_partitioning,
                      const MPI_Comm &communicator = MPI_COMM_WORLD);

      /**
       * 创建一个有重影的平行向量。            根据 @p ghost 参数是否在处理器之间唯一地细分元素，产生的向量可能有也可能没有鬼魂元素。更多信息请参见这个类的一般文档。              @see   @ref GlossGhostedVector  "有鬼魂元素的向量"
       *
       */
      Vector(const IndexSet &local,
             const IndexSet &ghost,
             const MPI_Comm &communicator = MPI_COMM_WORLD);

      /**
       * 从TrilinosWrappers向量类复制构造函数。由于这个类的向量不一定需要在进程间分配，用户需要向我们提供一个IndexSet和一个设置分区细节的MPI通信器。            根据 @p parallel_partitioning 参数是否在处理器之间唯一地细分元素，产生的向量可能有也可能没有鬼魂元素。更多信息请参见该类的一般文档。              @see   @ref GlossGhostedVector  "有鬼魂元素的向量"
       *
       */
      Vector(const IndexSet &parallel_partitioning,
             const Vector &  v,
             const MPI_Comm &communicator = MPI_COMM_WORLD);

      /**
       * 从deal.II向量中复制出来的构造函数。将维数设置为给定向量的维数，并复制所有的元素。            根据 @p parallel_partitioning 参数是否在处理器中唯一地划分元素，产生的向量可能有也可能没有鬼魂元素。更多信息请参见该类的一般文档。              @see   @ref GlossGhostedVector  "有鬼魂元素的向量"
       *
       */
      template <typename Number>
      Vector(const IndexSet &              parallel_partitioning,
             const dealii::Vector<Number> &v,
             const MPI_Comm &              communicator = MPI_COMM_WORLD);

      /**
       * 移动构造函数。通过窃取向量的内部数据创建一个新的向量
       * @p v.  。
       *
       */
      Vector(Vector &&v) noexcept;

      /**
       * 解构器。
       *
       */
      ~Vector() override = default;

      /**
       * 释放所有内存并返回到与调用默认构造函数后相同的状态。
       *
       */
      void
      clear();

      /**
       * 重新启动功能。这个函数将调用向量设置为输入向量的维度和平行分布，但不复制<tt>v</tt>中的元素。如果<tt>omit_zeroing_entries</tt>不是<tt>true</tt>，向量中的元素将被初始化为0。如果它被设置为<tt>true</tt>，那么向量的条目就处于一个未指定的状态，用户必须设置所有元素。在目前的实现中，如果矢量的布局与之前相比没有变化，这个方法不会触及矢量的条目，否则条目会被设置为零。
       * 请注意，这种行为在不同的版本中可能会发生变化而不被通知。
       * 这个函数有第三个参数，<tt>allow_different_maps</tt>，它允许在两个同等大小的向量之间交换数据（但在处理器之间的分布不同）。这个函数的一个微不足道的应用是在每台机器上生成整个向量的复制，当调用的向量是由每个进程上的所有索引组成的映射建立的，并且<tt>v</tt>是一个分布式向量。在这种情况下，变量<tt>omit_zeroing_entries</tt>需要被设置为<tt>false</tt>，因为在不同的并行化向量之间交换数据而不接触元素是没有意义的。
       *
       */
      void
      reinit(const Vector &v,
             const bool    omit_zeroing_entries = false,
             const bool    allow_different_maps = false);

      /**
       * Reinit功能。这个功能销毁了旧的向量内容，并根据输入的分区生成一个新的向量。 标志<tt>omit_zeroing_entries</tt>决定了向量是否应该被填充为零（false）。如果该标志被设置为<tt>true</tt>，则向量条目处于未指定的状态，用户必须设置所有元素。在目前的实现中，这个方法仍然将条目设置为零，但这可能会在没有通知的情况下在不同的版本中发生变化。            根据 @p parallel_partitioning 参数是否在处理器中唯一地细分元素，产生的向量可能有也可能没有鬼魂元素。更多信息请参见该类的一般文档。            在 @p parallel_partitioning 重叠的情况下，不清楚哪个进程应该拥有哪些元素。因此，在这种情况下，local_owned_elements()返回一个空的索引集。              @see   @ref GlossGhostedVector  "有幽灵元素的向量"
       *
       */
      void
      reinit(const IndexSet &parallel_partitioning,
             const MPI_Comm &communicator         = MPI_COMM_WORLD,
             const bool      omit_zeroing_entries = false);

      /**
       * 重新启动功能。这个功能销毁了旧的向量内容，并根据输入的分区生成一个新的向量。除了像上面所有其他方法一样只指定一个索引集外，这个方法还允许提供一个额外的鬼魂条目集。      有两个不同版本的向量可以被创建。如果标志 @p vector_writable 被设置为 @p false, ，该向量只允许读取 @p parallel_partitioning 和 @p ghost_entries. 的联合集合，那么reinit方法的效果相当于调用其他的reinit方法，其索引集包含本地拥有的条目和幽灵条目。            如果标志 @p vector_writable 被设置为true，这就为ghost元素创建了一个替代性的存储方案，允许多个线程写入向量（对于其他reinit方法，一次只允许一个线程写入ghost条目）。            根据 @p ghost_entries 参数是否在处理器之间唯一地细分元素，产生的向量可能有也可能没有鬼魂元素。更多信息请参见该类的一般文档。              @see   @ref GlossGhostedVector  "有鬼魂元素的向量"
       *
       */
      void
      reinit(const IndexSet &locally_owned_entries,
             const IndexSet &ghost_entries,
             const MPI_Comm &communicator    = MPI_COMM_WORLD,
             const bool      vector_writable = false);

      /**
       * 通过合并块状向量的成分来创建向量。
       *
       */
      void
      reinit(const BlockVector &v, const bool import_data = false);

      /**
       * 压缩Trilinos对象的底层表示，即刷新向量对象的缓冲区（如果它有的话）。在逐一写入矢量元素后，在对其进行其他操作之前，这个函数是必要的。            参数（缺省）可用于指定压缩模式（ <code>Add</code> or <code>Insert</code> ），以防在上次调用此函数后，没有对向量进行过写入。如果在上次调用compress()后，向量被添加或写入，则该参数被忽略。            更多信息见  @ref GlossCompress  "压缩分布式对象"
       * 。
       *
       */
      void
      compress(::dealii::VectorOperation::values operation);

      /**
       * 将向量的所有分量设置为给定的数字  @p s.
       * 只需将其传递给基类，但我们仍然需要声明这个函数，以使关于使构造函数显式的讨论中给出的例子发挥作用。
       * 由于将标量分配给向量的语义并不立即明确，所以只有当你想将整个向量设为0时，才能使用这个操作符。这样就可以使用直观的符号<tt>v=0</tt>。
       *
       */
      Vector &
      operator=(const TrilinosScalar s);

      /**
       * 复制给定的向量。如果有必要，可以调整现在的向量的大小。在这种情况下，设计平行分区的Epetra_Map也是从输入矢量中提取的。
       *
       */
      Vector &
      operator=(const Vector &v);

      /**
       * 移动给定的向量。该操作符通过有效地交换内部数据结构，用
       * @p v 替换当前的向量。
       *
       */
      Vector &
      operator=(Vector &&v) noexcept;

      /**
       * 另一个复制函数。这个函数接收一个deal.II向量并将其复制到一个TrilinosWrapper向量中。请注意，由于我们没有提供任何Epetra_map来告诉MPI进程之间的向量划分，TrilinosWrapper向量的大小必须与输入向量的大小相同。
       *
       */
      template <typename Number>
      Vector &
      operator=(const ::dealii::Vector<Number> &v);

      /**
       * 这个reinit函数是为了用于必须使用一些非本地数据的并行计算。需要这个函数的典型情况是并行调用
       * FEValues<dim>::get_function_values
       * 函数（或一些导数）。由于提前检索数据通常更快，这个函数可以在汇编分叉到不同处理器之前被调用。这个函数的作用如下。
       * 它获取给定矩阵的列中的信息，并寻找不同处理器之间的数据耦合。然后从输入向量中查询这些数据。注意，你不应该再向结果向量写入数据，因为有些数据可以在不同的处理器上存储多次，导致不可预测的结果。特别是，这样的向量不能用于矩阵-向量乘积，例如在线性系统的求解过程中。
       *
       */
      void
      import_nonlocal_data_for_fe(
        const dealii::TrilinosWrappers::SparseMatrix &matrix,
        const Vector &                                vector);

      /**
       * 从输入向量 @p rwv.  VectorOperation::values  @p operation
       * 中导入向量的IndexSet中存在的所有元素，用于决定 @p
       * rwv
       * 中的元素是否应该添加到当前向量中或替换当前元素。
       *
       */
      void
      import(const LinearAlgebra::ReadWriteVector<double> &rwv,
             const VectorOperation::values                 operation);


      /**
       * 检验是否相等。这个函数假定目前的向量和要比较的向量已经具有相同的大小，因为无论如何比较不同大小的向量没有什么意义。
       *
       */
      bool
      operator==(const Vector &v) const;

      /**
       * 测试不平等。这个函数假定现在的向量和要比较的向量已经有相同的大小，因为比较不同大小的向量反正没有什么意义。
       *
       */
      bool
      operator!=(const Vector &v) const;

      /**
       * 返回向量的全局尺寸。
       *
       */
      size_type
      size() const;

      /**
       * 返回向量的局部尺寸，即存储在当前MPI进程中的元素数量。对于顺序向量，这个数字与size()相同，但对于并行向量，它可能更小。
       * 要想知道到底哪些元素是存储在本地的，可以使用local_range()。
       * 如果向量包含鬼魂元素，它们会被包括在这个数字中。
       * @deprecated  这个函数已被废弃。
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
       * 返回一对索引，表明这个向量的哪些元素是本地存储的。第一个数字是存储的第一个元素的索引，第二个是本地存储的超过最后一个元素的索引。如果这是一个连续的向量，那么结果将是一对
       * <code>(0,N)</code>  ，否则将是一对 <code>(i,i+n)</code>, where
       * <code>n=local_size()</code> 和 <code>i</code>
       * 是存储在这个处理器上的向量的第一个元素，对应于半开区间
       * $[i,i+n)$
       * @note
       * 上面的描述在大多数时候是正确的，但并不总是。
       * 特别是，Trilinos向量不需要存储连续的元素范围，如
       * $[i,i+n)$
       * 。相反，它可以存储元素以任意的方式分布在所有处理器上的向量，每个处理器只是存储一个特定的子集，不一定是连续的。在这种情况下，这个函数显然是没有意义的，因为它最多只能返回一个包括本地存储的所有元素的范围。因此，这个函数只有在本地存储的范围确实是连续的情况下才会成功。如果向量的本地部分不是连续的，它将触发一个断言。
       *
       */
      std::pair<size_type, size_type>
      local_range() const;

      /**
       * 返回 @p index 是否在本地范围内，另见local_range()。
       * @note
       * 这个函数的适用性限制与local_range()的文档中所列相同。
       *
       */
      bool
      in_local_range(const size_type index) const;

      /**
       * 返回一个索引集，描述这个向量的哪些元素是由当前处理器拥有的。请注意，这个索引集不包括这个向量可能在本地存储为幽灵元素，但实际上由另一个处理器拥有的元素。因此，如果这是一个分布式向量，在不同处理器上返回的索引集将形成不相交的集合，加起来就是完整的索引集。
       * 很明显，如果一个向量只在一个处理器上创建，那么结果将满足
       * @code
       * vec.locally_owned_elements() == complete_index_set (vec.size())
       * @endcode
       *
       *
       */
      IndexSet
      locally_owned_elements() const;

      /**
       * 如果向量包含鬼魂元素，则返回。如果至少有一个进程上有鬼魂元素，这个答案就是真的。              @see   @ref GlossGhostedVector  "含有幽灵元素的向量"
       *
       */
      bool
      has_ghost_elements() const;

      /**
       * 这个函数只是为了与 @p   LinearAlgebra::distributed::Vector
       * 类兼容而存在，并不做任何事情：这个类以不同的方式实现了鬼魂值的更新，与底层的Trilinos向量对象更加匹配。
       *
       */
      void
      update_ghost_values() const;

      /**
       * 返回两个向量的标量（内）积。这些向量必须具有相同的大小。
       *
       */
      TrilinosScalar operator*(const Vector &vec) const;

      /**
       * 返回 $l_2$ -norm的平方。
       *
       */
      real_type
      norm_sqr() const;

      /**
       * 这个向量的元素的平均值。
       *
       */
      TrilinosScalar
      mean_value() const;

      /**
       * 计算这个向量的元素的最小值。
       *
       */
      TrilinosScalar
      min() const;

      /**
       * 计算这个向量的元素的最大值。
       *
       */
      TrilinosScalar
      max() const;

      /**
       * $l_1$  -矢量的规范。 绝对值的总和。
       *
       */
      real_type
      l1_norm() const;

      /**
       * $l_2$  - 矢量的规范。 各个元素的平方根之和。
       *
       */
      real_type
      l2_norm() const;

      /**
       * $l_p$  -
       * 矢量的规范。元素绝对值的<i>p</i>次幂之和的<i>p</i>根。
       *
       */
      real_type
      lp_norm(const TrilinosScalar p) const;

      /**
       * 元素的最大绝对值。
       *
       */
      real_type
      linfty_norm() const;

      /**
       * 执行矢量加法和后续内积的组合操作，返回内积的值。换句话说，这个函数的结果与用户调用的
       * @code
       * this->add(a, V);
       * return_value =this W;
       * @endcode
       * 这个函数存在的原因是为了与deal.II自己的向量类兼容，后者可以用较少的内存传输实现这个功能。然而，对于Trilinos向量来说，原生不支持这样的组合操作，因此成本完全等同于单独调用这两个方法。
       * 对于复值向量，第二步中的标量乘积被实现为
       * $\left<v,w\right>=\sum_i v_i \bar{w_i}$  .
       *
       */
      TrilinosScalar
      add_and_dot(const TrilinosScalar a, const Vector &V, const Vector &W);

      /**
       * 返回该向量是否只包含数值为0的元素。这是一个集体操作。这个函数很昂贵，因为有可能所有的元素都要被检查。
       *
       */
      bool
      all_zero() const;

      /**
       * 如果向量没有负的条目，即所有条目都是零或正，则返回
       * @p true
       * 。例如，这个函数用于检查细化指标是否真的都是正的（或零）。
       *
       */
      bool
      is_non_negative() const;
      //@}


      /**
       * @name  2：数据-获取
       *
       */
      //@{

      /**
       * 提供对一个给定元素的访问，包括读和写。
       * 当使用一个用MPI分布的向量时，这个操作只对调用处理器上实际存在的元素有意义。
       * 否则，会产生一个异常。
       *
       */
      reference
      operator()(const size_type index);

      /**
       * 提供对一个元素的只读访问。
       * 当使用MPI分布的向量时，这个操作只对调用处理器上实际存在的元素有意义。
       * 否则，会产生一个异常。
       *
       */
      TrilinosScalar
      operator()(const size_type index) const;

      /**
       * 提供对一个给定元素的访问，包括读和写。
       * 与operator()完全相同。
       *
       */
      reference operator[](const size_type index);

      /**
       * 提供对一个元素的只读访问。
       * 与operator()完全相同。
       *
       */
      TrilinosScalar operator[](const size_type index) const;

      /**
       * 与通过operator()获取向量中的单个元素不同，这个函数允许一次性获取一整组元素。要读取的元素的索引在第一个参数中说明，相应的值在第二个参数中返回。
       * 如果当前的向量被称为 @p v,
       * ，那么这个函数就等同于代码
       * @code
       * for (unsigned int i=0; i<indices.size(); ++i)
       *   values[i] = v[indices[i]];
       * @endcode
       * @pre  @p indices 和 @p values 数组的大小必须是一致的。
       *
       */
      void
      extract_subvector_to(const std::vector<size_type> &indices,
                           std::vector<TrilinosScalar> & values) const;

      /**
       * 这个函数不是通过operator()获得向量的单个元素，而是允许一次获得整个元素集。与前一个函数不同的是，这个函数通过取消引用前两个参数提供的迭代器范围内的所有元素来获得元素的索引，并将向量的值放入通过取消引用从第三个参数指向的位置开始的迭代器范围获得的内存位置。
       * 如果当前的向量被称为 @p v,
       * ，那么这个函数就等同于代码
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
       * 通过返回该向量中本地拥有的元素的开始和结束的迭代器，使Vector类有点像C++标准库中的<tt>vector<>/tt>类。如果向量是由IndexSet或deal.II中的其他方法构造的，则本地元素的排序与全局索引给出的排序一致（注意，Epetra_Map可以包含任意排序的元素）。
       * 它认为end()
       *
       * - begin() == local_size()。
       *
       */
      iterator
      begin();

      /**
       * 返回常数迭代器到向量中本地拥有的元素的起点。
       *
       */
      const_iterator
      begin() const;

      /**
       * 返回一个迭代器，指向本地拥有的条目数组结束后的元素。
       *
       */
      iterator
      end();

      /**
       * 返回一个恒定的迭代器，指向本地拥有的条目数组结束后的元素。
       *
       */
      const_iterator
      end() const;

      //@}


      /**
       * @name  3: 修改向量
       *
       */
      //@{

      /**
       * 一个集体的设置操作：这个函数允许一次性设置一整组元素，而不是设置一个向量中的单个元素。
       * 要设置的元素的索引在第一个参数中说明，相应的值在第二个参数中说明。
       *
       */
      void
      set(const std::vector<size_type> &     indices,
          const std::vector<TrilinosScalar> &values);

      /**
       * 这是第二个集体设置操作。作为区别，这个函数需要一个deal.II的数值向量。
       *
       */
      void
      set(const std::vector<size_type> &          indices,
          const ::dealii::Vector<TrilinosScalar> &values);

      /**
       * 这个集合操作级别较低，可以处理其他任何事情&mdash；你唯一需要提供的是一个存放所有索引的地址和要设置的元素数量。
       *
       */
      void
      set(const size_type       n_elements,
          const size_type *     indices,
          const TrilinosScalar *values);

      /**
       * 一个集体添加操作。这个函数将存储在 @p values
       * 中的整组值添加到 @p indices. 指定的向量成分中。
       *
       */
      void
      add(const std::vector<size_type> &     indices,
          const std::vector<TrilinosScalar> &values);

      /**
       * 这是第二次集体添加操作。作为区别，这个函数需要一个deal.II的数值向量。
       *
       */
      void
      add(const std::vector<size_type> &          indices,
          const ::dealii::Vector<TrilinosScalar> &values);

      /**
       * 取一个<tt>n_elements</tt>连续存储的地址，并将其添加到向量中。处理上述其他两个<tt>add()</tt>函数未涵盖的所有情况。
       *
       */
      void
      add(const size_type       n_elements,
          const size_type *     indices,
          const TrilinosScalar *values);

      /**
       * 将整个向量乘以一个固定系数。
       *
       */
      Vector &
      operator*=(const TrilinosScalar factor);

      /**
       * 将整个向量除以一个固定的因子。
       *
       */
      Vector &
      operator/=(const TrilinosScalar factor);

      /**
       * 将给定的向量添加到当前的向量中。
       *
       */
      Vector &
      operator+=(const Vector &V);

      /**
       * 从现在的向量中减去给定的向量。
       *
       */
      Vector &
      operator-=(const Vector &V);

      /**
       * 将 @p s 加到所有组件上。注意 @p s
       * 是一个标量而不是一个向量。
       *
       */
      void
      add(const TrilinosScalar s);

      /**
       * 简单的向量加法，等于<tt>运算符+=</tt>。
       * 不过，如果第二个参数<tt>allow_different_maps</tt>被设置，那么就有可能从一个使用不同地图的向量中添加数据，也就是说，一个向量的元素被不同的处理器分割。例如，这可能包括有幽灵元素的向量。
       * 然而，一般来说，添加具有不同元素对处理器映射的向量需要在处理器之间进行数据通信，因此，比使用相同映射的向量的操作要慢。
       *
       */
      void
      add(const Vector &V, const bool allow_different_maps = false);

      /**
       * 矢量的倍数的简单加法，即<tt>*this += a*V</tt>。
       *
       */
      void
      add(const TrilinosScalar a, const Vector &V);

      /**
       * 缩放向量的多重加法，即<tt>*this += a*V + b*W</tt>。
       *
       */
      void
      add(const TrilinosScalar a,
          const Vector &       V,
          const TrilinosScalar b,
          const Vector &       W);

      /**
       * 缩放和简单的向量相加，即<tt>*this = s*(*this) + V</tt>。
       *
       */
      void
      sadd(const TrilinosScalar s, const Vector &V);

      /**
       * 缩放和简单加法，即：<tt>*this = s*(*this) + a*V</tt>。
       *
       */
      void
      sadd(const TrilinosScalar s, const TrilinosScalar a, const Vector &V);

      /**
       * 用参数中的相应元素来缩放这个向量的每个元素。这个函数主要是为了模拟对角线缩放矩阵的乘法（和立即重新分配）。
       *
       */
      void
      scale(const Vector &scaling_factors);

      /**
       * 赋值 <tt>*this = a*V</tt>.
       *
       */
      void
      equ(const TrilinosScalar a, const Vector &V);
      //@}

      /**
       * @name  4：混合东西
       *
       */
      //@{

      /**
       * 返回一个对底层Trilinos Epetra_MultiVector类的常量引用。
       *
       */
      const Epetra_MultiVector &
      trilinos_vector() const;

      /**
       * 返回一个（可修改的）对底层Trilinos
       * Epetra_FEVector类的引用。
       *
       */
      Epetra_FEVector &
      trilinos_vector();

      /**
       * 返回一个对底层Trilinos
       * Epetra_BlockMap的常量引用，它设置了向量的平行分区。
       *
       */
      const Epetra_BlockMap &
      trilinos_partitioner() const;

      /**
       * 打印到一个流。  @p precision
       * 表示打印数值所需的精度， @p scientific
       * 是否应使用科学符号。如果 @p across 是 @p true
       * ，那么向量将被打印在一行中，而如果 @p false
       * 则元素被打印在单独的一行中。
       *
       */
      void
      print(std::ostream &     out,
            const unsigned int precision  = 3,
            const bool         scientific = true,
            const bool         across     = true) const;

      /**
       * 交换这个向量和另一个向量的内容  @p v.
       * 人们可以用一个临时变量和复制数据元素来完成这个操作，但是这个函数明显更有效率，因为它只交换了两个向量的数据指针，因此不需要分配临时存储和移动数据。注意，两个向量需要有相同的大小，并以同一地图为基础。
       * 这个函数类似于所有C++标准容器的 @p swap
       * 函数。此外，还有一个全局函数<tt>swap(u,v)</tt>，它简单地调用<tt>u.swap(v)</tt>，同样与标准函数相类似。
       *
       */
      void
      swap(Vector &v);

      /**
       * 对内存消耗的估计，以字节为单位。
       *
       */
      std::size_t
      memory_consumption() const;

      /**
       * 返回一个对与此对象一起使用的MPI通信器对象的引用。
       *
       */
      const MPI_Comm &
      get_mpi_communicator() const;
      //@}

      /**
       * 异常情况
       *
       */
      DeclException0(ExcDifferentParallelPartitioning);

      /**
       * 异常情况
       *
       */
      DeclException1(ExcTrilinosError,
                     int,
                     << "An error with error number " << arg1
                     << " occurred while calling a Trilinos function");

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
        << "You are trying to access element " << arg1
        << " of a distributed vector, but this element is not stored "
        << "on the current processor. Note: There are " << arg2
        << " elements stored "
        << "on the current processor from within the range [" << arg3 << ","
        << arg4 << "] but Trilinos vectors need not store contiguous "
        << "ranges on each processor, and not every element in "
        << "this range may in fact be stored locally."
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
       * Trilinos不允许混合添加矩阵条目和覆盖它们（以使并行计算的同步更简单）。我们的方法是，对于每个访问操作，存储它是插入还是添加。如果之前的操作是不同的类型，那么我们首先要刷新Trilinos缓冲区；否则，我们可以简单地继续下去。
       * 幸运的是，Trilinos有一个这样的对象，在这种情况下已经完成了所有的并行通信，所以我们只需使用他们的模型，存储最后一个操作是加法还是插入。
       *
       */
      Epetra_CombineMode last_action;

      /**
       * 一个布尔变量，用来保存向量是否被压缩的信息。
       *
       */
      bool compressed;

      /**
       * 这个向量是否有鬼魂元素。这在所有的处理器上都是真的，即使只有一个处理器有任何鬼魂元素。
       *
       */
      bool has_ghosts;

      /**
       * 指向实际Epetra向量对象的指针。这可能代表一个事实上分布在多个处理器上的向量。该对象在设置时需要一个现有的Epetra_Map来存储数据。
       *
       */
      std::unique_ptr<Epetra_FEVector> vector;

      /**
       * Trilinos中的一个向量对象，用于收集非本地元素，如果该向量是用描述鬼魂元素的额外IndexSet构造的。
       *
       */
      std::unique_ptr<Epetra_MultiVector> nonlocal_vector;

      /**
       * 一个IndexSet，存储这个向量专门拥有的索引。
       *
       */
      IndexSet owned_elements;

      // Make the reference class a friend.
      friend class internal::VectorReference;
    };



    // ------------------- inline and template functions --------------


    /**
     * 全局函数 @p swap
     * ，它重载了C++标准库的默认实现，它使用一个临时对象。该函数简单地交换了两个向量的数据。
     * @relatesalso   TrilinosWrappers::MPI::Vector .
     *
     */
    inline void
    swap(Vector &u, Vector &v)
    {
      u.swap(v);
    }
  } // namespace MPI

#  ifndef DOXYGEN

  namespace internal
  {
    inline VectorReference::VectorReference(MPI::Vector &   vector,
                                            const size_type index)
      : vector(vector)
      , index(index)
    {}


    inline const VectorReference &
    VectorReference::operator=(const VectorReference &r) const
    {
      // as explained in the class
      // documentation, this is not the copy
      // operator. so simply pass on to the
      // "correct" assignment operator
      *this = static_cast<TrilinosScalar>(r);

      return *this;
    }



    inline VectorReference &
    VectorReference::operator=(const VectorReference &r)
    {
      // as above
      *this = static_cast<TrilinosScalar>(r);

      return *this;
    }


    inline const VectorReference &
    VectorReference::operator=(const TrilinosScalar &value) const
    {
      vector.set(1, &index, &value);
      return *this;
    }



    inline const VectorReference &
    VectorReference::operator+=(const TrilinosScalar &value) const
    {
      vector.add(1, &index, &value);
      return *this;
    }



    inline const VectorReference &
    VectorReference::operator-=(const TrilinosScalar &value) const
    {
      TrilinosScalar new_value = -value;
      vector.add(1, &index, &new_value);
      return *this;
    }



    inline const VectorReference &
    VectorReference::operator*=(const TrilinosScalar &value) const
    {
      TrilinosScalar new_value = static_cast<TrilinosScalar>(*this) * value;
      vector.set(1, &index, &new_value);
      return *this;
    }



    inline const VectorReference &
    VectorReference::operator/=(const TrilinosScalar &value) const
    {
      TrilinosScalar new_value = static_cast<TrilinosScalar>(*this) / value;
      vector.set(1, &index, &new_value);
      return *this;
    }
  } // namespace internal

  namespace MPI
  {
    inline bool
    Vector::in_local_range(const size_type index) const
    {
      std::pair<size_type, size_type> range = local_range();

      return ((index >= range.first) && (index < range.second));
    }



    inline IndexSet
    Vector::locally_owned_elements() const
    {
      Assert(owned_elements.size() == size(),
             ExcMessage(
               "The locally owned elements have not been properly initialized!"
               " This happens for example if this object has been initialized"
               " with exactly one overlapping IndexSet."));
      return owned_elements;
    }



    inline bool
    Vector::has_ghost_elements() const
    {
      return has_ghosts;
    }



    inline void
    Vector::update_ghost_values() const
    {}



    inline internal::VectorReference
    Vector::operator()(const size_type index)
    {
      return internal::VectorReference(*this, index);
    }



    inline internal::VectorReference Vector::operator[](const size_type index)
    {
      return operator()(index);
    }



    inline TrilinosScalar Vector::operator[](const size_type index) const
    {
      return operator()(index);
    }



    inline void
    Vector::extract_subvector_to(const std::vector<size_type> &indices,
                                 std::vector<TrilinosScalar> & values) const
    {
      for (size_type i = 0; i < indices.size(); ++i)
        values[i] = operator()(indices[i]);
    }



    template <typename ForwardIterator, typename OutputIterator>
    inline void
    Vector::extract_subvector_to(ForwardIterator       indices_begin,
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



    inline Vector::iterator
    Vector::begin()
    {
      return (*vector)[0];
    }



    inline Vector::iterator
    Vector::end()
    {
      return (*vector)[0] + locally_owned_size();
    }



    inline Vector::const_iterator
    Vector::begin() const
    {
      return (*vector)[0];
    }



    inline Vector::const_iterator
    Vector::end() const
    {
      return (*vector)[0] + locally_owned_size();
    }



    inline void
    Vector::set(const std::vector<size_type> &     indices,
                const std::vector<TrilinosScalar> &values)
    {
      // if we have ghost values, do not allow
      // writing to this vector at all.
      Assert(!has_ghost_elements(), ExcGhostsPresent());

      AssertDimension(indices.size(), values.size());

      set(indices.size(), indices.data(), values.data());
    }



    inline void
    Vector::set(const std::vector<size_type> &          indices,
                const ::dealii::Vector<TrilinosScalar> &values)
    {
      // if we have ghost values, do not allow
      // writing to this vector at all.
      Assert(!has_ghost_elements(), ExcGhostsPresent());

      AssertDimension(indices.size(), values.size());

      set(indices.size(), indices.data(), values.begin());
    }



    inline void
    Vector::set(const size_type       n_elements,
                const size_type *     indices,
                const TrilinosScalar *values)
    {
      // if we have ghost values, do not allow
      // writing to this vector at all.
      Assert(!has_ghost_elements(), ExcGhostsPresent());

      if (last_action == Add)
        {
          const int ierr = vector->GlobalAssemble(Add);
          AssertThrow(ierr == 0, ExcTrilinosError(ierr));
        }

      if (last_action != Insert)
        last_action = Insert;

      for (size_type i = 0; i < n_elements; ++i)
        {
          const TrilinosWrappers::types::int_type row = indices[i];
          const TrilinosWrappers::types::int_type local_row =
            vector->Map().LID(row);
          if (local_row != -1)
            (*vector)[0][local_row] = values[i];
          else
            {
              const int ierr = vector->ReplaceGlobalValues(1, &row, &values[i]);
              AssertThrow(ierr == 0, ExcTrilinosError(ierr));
              compressed = false;
            }
          // in set operation, do not use the pre-allocated vector for nonlocal
          // entries even if it exists. This is to ensure that we really only
          // set the elements touched by the set() method and not all contained
          // in the nonlocal entries vector (there is no way to distinguish them
          // on the receiving processor)
        }
    }



    inline void
    Vector::add(const std::vector<size_type> &     indices,
                const std::vector<TrilinosScalar> &values)
    {
      // if we have ghost values, do not allow
      // writing to this vector at all.
      Assert(!has_ghost_elements(), ExcGhostsPresent());
      AssertDimension(indices.size(), values.size());

      add(indices.size(), indices.data(), values.data());
    }



    inline void
    Vector::add(const std::vector<size_type> &          indices,
                const ::dealii::Vector<TrilinosScalar> &values)
    {
      // if we have ghost values, do not allow
      // writing to this vector at all.
      Assert(!has_ghost_elements(), ExcGhostsPresent());
      AssertDimension(indices.size(), values.size());

      add(indices.size(), indices.data(), values.begin());
    }



    inline void
    Vector::add(const size_type       n_elements,
                const size_type *     indices,
                const TrilinosScalar *values)
    {
      // if we have ghost values, do not allow
      // writing to this vector at all.
      Assert(!has_ghost_elements(), ExcGhostsPresent());

      if (last_action != Add)
        {
          if (last_action == Insert)
            {
              const int ierr = vector->GlobalAssemble(Insert);
              AssertThrow(ierr == 0, ExcTrilinosError(ierr));
            }
          last_action = Add;
        }

      for (size_type i = 0; i < n_elements; ++i)
        {
          const size_type                         row       = indices[i];
          const TrilinosWrappers::types::int_type local_row = vector->Map().LID(
            static_cast<TrilinosWrappers::types::int_type>(row));
          if (local_row != -1)
            (*vector)[0][local_row] += values[i];
          else if (nonlocal_vector.get() == nullptr)
            {
              const int ierr = vector->SumIntoGlobalValues(
                1,
                reinterpret_cast<const TrilinosWrappers::types::int_type *>(
                  &row),
                &values[i]);
              AssertThrow(ierr == 0, ExcTrilinosError(ierr));
              compressed = false;
            }
          else
            {
              // use pre-allocated vector for non-local entries if it exists for
              // addition operation
              const TrilinosWrappers::types::int_type my_row =
                nonlocal_vector->Map().LID(
                  static_cast<TrilinosWrappers::types::int_type>(row));
              Assert(my_row != -1,
                     ExcMessage(
                       "Attempted to write into off-processor vector entry "
                       "that has not be specified as being writable upon "
                       "initialization"));
              (*nonlocal_vector)[0][my_row] += values[i];
              compressed = false;
            }
        }
    }



    inline Vector::size_type
    Vector::size() const
    {
#    ifndef DEAL_II_WITH_64BIT_INDICES
      return vector->Map().MaxAllGID() + 1 - vector->Map().MinAllGID();
#    else
      return vector->Map().MaxAllGID64() + 1 - vector->Map().MinAllGID64();
#    endif
    }



    inline Vector::size_type
    Vector::local_size() const
    {
      return vector->Map().NumMyElements();
    }



    inline Vector::size_type
    Vector::locally_owned_size() const
    {
      return owned_elements.n_elements();
    }



    inline std::pair<Vector::size_type, Vector::size_type>
    Vector::local_range() const
    {
#    ifndef DEAL_II_WITH_64BIT_INDICES
      const TrilinosWrappers::types::int_type begin = vector->Map().MinMyGID();
      const TrilinosWrappers::types::int_type end =
        vector->Map().MaxMyGID() + 1;
#    else
      const TrilinosWrappers::types::int_type begin =
        vector->Map().MinMyGID64();
      const TrilinosWrappers::types::int_type end =
        vector->Map().MaxMyGID64() + 1;
#    endif

      Assert(
        end - begin == vector->Map().NumMyElements(),
        ExcMessage(
          "This function only makes sense if the elements that this "
          "vector stores on the current processor form a contiguous range. "
          "This does not appear to be the case for the current vector."));

      return std::make_pair(begin, end);
    }



    inline TrilinosScalar Vector::operator*(const Vector &vec) const
    {
      Assert(vector->Map().SameAs(vec.vector->Map()),
             ExcDifferentParallelPartitioning());
      Assert(!has_ghost_elements(), ExcGhostsPresent());

      TrilinosScalar result;

      const int ierr = vector->Dot(*(vec.vector), &result);
      AssertThrow(ierr == 0, ExcTrilinosError(ierr));

      return result;
    }



    inline Vector::real_type
    Vector::norm_sqr() const
    {
      const TrilinosScalar d = l2_norm();
      return d * d;
    }



    inline TrilinosScalar
    Vector::mean_value() const
    {
      Assert(!has_ghost_elements(), ExcGhostsPresent());

      TrilinosScalar mean;
      const int      ierr = vector->MeanValue(&mean);
      AssertThrow(ierr == 0, ExcTrilinosError(ierr));

      return mean;
    }



    inline TrilinosScalar
    Vector::min() const
    {
      TrilinosScalar min_value;
      const int      ierr = vector->MinValue(&min_value);
      AssertThrow(ierr == 0, ExcTrilinosError(ierr));

      return min_value;
    }



    inline TrilinosScalar
    Vector::max() const
    {
      TrilinosScalar max_value;
      const int      ierr = vector->MaxValue(&max_value);
      AssertThrow(ierr == 0, ExcTrilinosError(ierr));

      return max_value;
    }



    inline Vector::real_type
    Vector::l1_norm() const
    {
      Assert(!has_ghost_elements(), ExcGhostsPresent());

      TrilinosScalar d;
      const int      ierr = vector->Norm1(&d);
      AssertThrow(ierr == 0, ExcTrilinosError(ierr));

      return d;
    }



    inline Vector::real_type
    Vector::l2_norm() const
    {
      Assert(!has_ghost_elements(), ExcGhostsPresent());

      TrilinosScalar d;
      const int      ierr = vector->Norm2(&d);
      AssertThrow(ierr == 0, ExcTrilinosError(ierr));

      return d;
    }



    inline Vector::real_type
    Vector::lp_norm(const TrilinosScalar p) const
    {
      Assert(!has_ghost_elements(), ExcGhostsPresent());

      TrilinosScalar  norm    = 0;
      TrilinosScalar  sum     = 0;
      const size_type n_local = locally_owned_size();

      // loop over all the elements because
      // Trilinos does not support lp norms
      for (size_type i = 0; i < n_local; ++i)
        sum += std::pow(std::fabs((*vector)[0][i]), p);

      norm = std::pow(sum, static_cast<TrilinosScalar>(1. / p));

      return norm;
    }



    inline Vector::real_type
    Vector::linfty_norm() const
    {
      // while we disallow the other
      // norm operations on ghosted
      // vectors, this particular norm
      // is safe to run even in the
      // presence of ghost elements
      TrilinosScalar d;
      const int      ierr = vector->NormInf(&d);
      AssertThrow(ierr == 0, ExcTrilinosError(ierr));

      return d;
    }



    inline TrilinosScalar
    Vector::add_and_dot(const TrilinosScalar a,
                        const Vector &       V,
                        const Vector &       W)
    {
      this->add(a, V);
      return *this * W;
    }



    // inline also scalar products, vector
    // additions etc. since they are all
    // representable by a single Trilinos
    // call. This reduces the overhead of the
    // wrapper class.
    inline Vector &
    Vector::operator*=(const TrilinosScalar a)
    {
      AssertIsFinite(a);

      const int ierr = vector->Scale(a);
      AssertThrow(ierr == 0, ExcTrilinosError(ierr));

      return *this;
    }



    inline Vector &
    Vector::operator/=(const TrilinosScalar a)
    {
      AssertIsFinite(a);

      const TrilinosScalar factor = 1. / a;

      AssertIsFinite(factor);

      const int ierr = vector->Scale(factor);
      AssertThrow(ierr == 0, ExcTrilinosError(ierr));

      return *this;
    }



    inline Vector &
    Vector::operator+=(const Vector &v)
    {
      AssertDimension(size(), v.size());
      Assert(vector->Map().SameAs(v.vector->Map()),
             ExcDifferentParallelPartitioning());

      const int ierr = vector->Update(1.0, *(v.vector), 1.0);
      AssertThrow(ierr == 0, ExcTrilinosError(ierr));

      return *this;
    }



    inline Vector &
    Vector::operator-=(const Vector &v)
    {
      AssertDimension(size(), v.size());
      Assert(vector->Map().SameAs(v.vector->Map()),
             ExcDifferentParallelPartitioning());

      const int ierr = vector->Update(-1.0, *(v.vector), 1.0);
      AssertThrow(ierr == 0, ExcTrilinosError(ierr));

      return *this;
    }



    inline void
    Vector::add(const TrilinosScalar s)
    {
      // if we have ghost values, do not allow
      // writing to this vector at all.
      Assert(!has_ghost_elements(), ExcGhostsPresent());
      AssertIsFinite(s);

      size_type n_local = locally_owned_size();
      for (size_type i = 0; i < n_local; ++i)
        (*vector)[0][i] += s;
    }



    inline void
    Vector::add(const TrilinosScalar a, const Vector &v)
    {
      // if we have ghost values, do not allow
      // writing to this vector at all.
      Assert(!has_ghost_elements(), ExcGhostsPresent());
      AssertDimension(locally_owned_size(), v.locally_owned_size());

      AssertIsFinite(a);

      const int ierr = vector->Update(a, *(v.vector), 1.);
      AssertThrow(ierr == 0, ExcTrilinosError(ierr));
    }



    inline void
    Vector::add(const TrilinosScalar a,
                const Vector &       v,
                const TrilinosScalar b,
                const Vector &       w)
    {
      // if we have ghost values, do not allow
      // writing to this vector at all.
      Assert(!has_ghost_elements(), ExcGhostsPresent());
      AssertDimension(locally_owned_size(), v.locally_owned_size());
      AssertDimension(locally_owned_size(), w.locally_owned_size());

      AssertIsFinite(a);
      AssertIsFinite(b);

      const int ierr = vector->Update(a, *(v.vector), b, *(w.vector), 1.);

      AssertThrow(ierr == 0, ExcTrilinosError(ierr));
    }



    inline void
    Vector::sadd(const TrilinosScalar s, const Vector &v)
    {
      // if we have ghost values, do not allow
      // writing to this vector at all.
      Assert(!has_ghost_elements(), ExcGhostsPresent());
      AssertDimension(size(), v.size());

      AssertIsFinite(s);

      // We assume that the vectors have the same Map
      // if the local size is the same and if the vectors are not ghosted
      if (locally_owned_size() == v.locally_owned_size() &&
          !v.has_ghost_elements())
        {
          Assert(this->vector->Map().SameAs(v.vector->Map()) == true,
                 ExcDifferentParallelPartitioning());
          const int ierr = vector->Update(1., *(v.vector), s);
          AssertThrow(ierr == 0, ExcTrilinosError(ierr));
        }
      else
        {
          (*this) *= s;
          this->add(v, true);
        }
    }



    inline void
    Vector::sadd(const TrilinosScalar s,
                 const TrilinosScalar a,
                 const Vector &       v)
    {
      // if we have ghost values, do not allow
      // writing to this vector at all.
      Assert(!has_ghost_elements(), ExcGhostsPresent());
      AssertDimension(size(), v.size());
      AssertIsFinite(s);
      AssertIsFinite(a);

      // We assume that the vectors have the same Map
      // if the local size is the same and if the vectors are not ghosted
      if (locally_owned_size() == v.locally_owned_size() &&
          !v.has_ghost_elements())
        {
          Assert(this->vector->Map().SameAs(v.vector->Map()) == true,
                 ExcDifferentParallelPartitioning());
          const int ierr = vector->Update(a, *(v.vector), s);
          AssertThrow(ierr == 0, ExcTrilinosError(ierr));
        }
      else
        {
          (*this) *= s;
          Vector tmp = v;
          tmp *= a;
          this->add(tmp, true);
        }
    }



    inline void
    Vector::scale(const Vector &factors)
    {
      // if we have ghost values, do not allow
      // writing to this vector at all.
      Assert(!has_ghost_elements(), ExcGhostsPresent());
      AssertDimension(locally_owned_size(), factors.locally_owned_size());

      const int ierr = vector->Multiply(1.0, *(factors.vector), *vector, 0.0);
      AssertThrow(ierr == 0, ExcTrilinosError(ierr));
    }



    inline void
    Vector::equ(const TrilinosScalar a, const Vector &v)
    {
      // if we have ghost values, do not allow
      // writing to this vector at all.
      Assert(!has_ghost_elements(), ExcGhostsPresent());
      AssertIsFinite(a);

      // If we don't have the same map, copy.
      if (vector->Map().SameAs(v.vector->Map()) == false)
        {
          this->sadd(0., a, v);
        }
      else
        {
          // Otherwise, just update
          int ierr = vector->Update(a, *v.vector, 0.0);
          AssertThrow(ierr == 0, ExcTrilinosError(ierr));

          last_action = Zero;
        }
    }



    inline const Epetra_MultiVector &
    Vector::trilinos_vector() const
    {
      return static_cast<const Epetra_MultiVector &>(*vector);
    }



    inline Epetra_FEVector &
    Vector::trilinos_vector()
    {
      return *vector;
    }



    inline const Epetra_BlockMap &
    Vector::trilinos_partitioner() const
    {
      return vector->Map();
    }



    inline const MPI_Comm &
    Vector::get_mpi_communicator() const
    {
      static MPI_Comm comm;

#    ifdef DEAL_II_WITH_MPI

      const Epetra_MpiComm *mpi_comm =
        dynamic_cast<const Epetra_MpiComm *>(&vector->Map().Comm());
      comm = mpi_comm->Comm();

#    else

      comm = MPI_COMM_SELF;

#    endif

      return comm;
    }

    template <typename number>
    Vector::Vector(const IndexSet &              parallel_partitioner,
                   const dealii::Vector<number> &v,
                   const MPI_Comm &              communicator)
    {
      *this =
        Vector(parallel_partitioner.make_trilinos_map(communicator, true), v);
      owned_elements = parallel_partitioner;
    }



    inline Vector &
    Vector::operator=(const TrilinosScalar s)
    {
      AssertIsFinite(s);

      int ierr = vector->PutScalar(s);
      AssertThrow(ierr == 0, ExcTrilinosError(ierr));

      if (nonlocal_vector.get() != nullptr)
        {
          ierr = nonlocal_vector->PutScalar(0.);
          AssertThrow(ierr == 0, ExcTrilinosError(ierr));
        }

      return *this;
    }
  }  /* end of namespace MPI */ 

#  endif  /* DOXYGEN */ 

}  /* end of namespace TrilinosWrappers */ 

 /*@}*/ 


namespace internal
{
  namespace LinearOperatorImplementation
  {
    template <typename>
    class ReinitHelper;

    /**
     * linear_operator.h中内部使用的一个辅助类。对
     * TrilinosWrappers::MPI::Vector. 的特殊化。
     *
     */
    template <>
    class ReinitHelper<TrilinosWrappers::MPI::Vector>
    {
    public:
      template <typename Matrix>
      static void
      reinit_range_vector(const Matrix &                 matrix,
                          TrilinosWrappers::MPI::Vector &v,
                          bool                           omit_zeroing_entries)
      {
        v.reinit(matrix.locally_owned_range_indices(),
                 matrix.get_mpi_communicator(),
                 omit_zeroing_entries);
      }

      template <typename Matrix>
      static void
      reinit_domain_vector(const Matrix &                 matrix,
                           TrilinosWrappers::MPI::Vector &v,
                           bool                           omit_zeroing_entries)
      {
        v.reinit(matrix.locally_owned_domain_indices(),
                 matrix.get_mpi_communicator(),
                 omit_zeroing_entries);
      }
    };

  } // namespace LinearOperatorImplementation
}  /* namespace internal */ 



/**
 * 将 dealii::TrilinosWrappers::MPI::Vector 声明为分布式向量。
 *
 *
 */
template <>
struct is_serial_vector<TrilinosWrappers::MPI::Vector> : std::false_type
{};


DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_TRILINOS

 /*----------------------------   trilinos_vector.h ---------------------------*/ 

#endif // dealii_trilinos_vector_h


