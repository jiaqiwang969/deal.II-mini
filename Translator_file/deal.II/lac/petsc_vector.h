//include/deal.II-translator/lac/petsc_vector_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2021 by the deal.II authors
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

#ifndef dealii_petsc_vector_h
#  define dealii_petsc_vector_h


#  include <deal.II/base/config.h>

#  ifdef DEAL_II_WITH_PETSC

#    include <deal.II/base/index_set.h>
#    include <deal.II/base/subscriptor.h>

#    include <deal.II/lac/exceptions.h>
#    include <deal.II/lac/petsc_vector_base.h>
#    include <deal.II/lac/vector.h>
#    include <deal.II/lac/vector_operation.h>
#    include <deal.II/lac/vector_type_traits.h>

DEAL_II_NAMESPACE_OPEN


/*!   @addtogroup  PETScWrappers  
     * @{ 

 
*
*/
namespace PETScWrappers
{
  /**
   * 用于在MPI上并行工作的PETSc类的命名空间，例如分布式向量和矩阵。
   * @ingroup PETScWrappers
   *
   */
  namespace MPI
  {
    /**
     * 实现一个基于PETSC的并行向量类，并使用MPI通信来同步分布式操作。所有的功能实际上都在基类中，除了生成并行向量的调用。这是可能的，因为PETSc只在一个抽象的向量类型上工作，并在内部分配给根据实际向量类型进行实际工作的函数（很像使用虚拟函数）。只有创建特定类型向量的函数不同，并在这个特定的类中实现。
     * <h3>Parallel communication model</h3>
     * PETSc的并行功能是建立在消息传递接口（MPI）之上的。MPI的通信模型建立在集体通信的基础上：如果一个进程想从另一个进程那里得到什么，另一个进程必须愿意接受这种通信。一个进程不能通过调用一个远程函数从另一个进程查询数据，而另一个进程不期望有这种交易。其结果是，这个类的基类中的大多数操作都必须被集体调用。例如，如果你想计算一个平行向量的l2准则，
     * @em 所有共享这个向量的进程都必须调用 @p l2_norm
     * 函数。如果你不这样做，而只是在一个进程上调用 @p
     * l2_norm
     * 函数，那么会发生以下情况。这一个进程将调用一个集体MPI函数，并等待所有其他进程加入到这个过程中。由于其他进程不调用这个函数，你将在第一个进程上得到一个超时，或者更糟糕的是，当下一个对PETSc函数的调用在其他进程上产生MPI消息时，你将得到一个神秘的消息，即只有一个子集的进程试图进行通信。这些错误是很难搞清楚的，除非你很熟悉MPI的通信模型，并且知道哪些函数可能产生MPI消息。
     * 下面将讨论一个特殊的情况，即可能会意外地产生MPI消息的情况。
     * <h3>Accessing individual elements of a vector</h3>
     * PETSc确实允许对向量的单个元素进行读取访问，但在分布式情况下，只允许读取本地存储的元素。我们通过<tt>d=vec(i)</tt>这样的调用来实现。然而，如果你访问本地存储范围之外的元素，就会产生一个异常。
     * 与读访问相反，PETSc（和相应的deal.II包装类）允许对向量的单个元素进行写入（或添加），即使它们存储在不同的进程中。你可以这样写，例如，<tt>vec(i)=d</tt>或<tt>vec(i)+=d</tt>，或类似的操作。但是有一个问题，可能会导致非常混乱的错误信息。PETSc要求应用程序在从对元素的添加转换到对元素的写入时调用compress()函数。其理由是，所有进程都可能积累对元素的加法操作，即使是多个进程对相同的元素进行写入。当我们下一次调用compress()时，所有这些加法操作都已执行完毕。然而，如果一个进程对一个元素进行添加，而另一个进程对其进行覆盖，如果我们不确保在两者之间发生与compress()的同步，执行的顺序将产生非确定性的行为。
     * 为了确保这些对compress()的调用在适当的时间发生，deal.II包装器保留了一个状态变量，用于存储当前允许的操作：添加或写入。如果它遇到了相反的操作，它就会调用compress()并翻转状态。这有时会导致非常混乱的行为，例如，代码可能看起来像这样。
     * @code
     * PETScWrappers::MPI::Vector vector;
     * ...
     * // do some write operations on the vector
     * for (unsigned int i=0; i<vector.size(); ++i)
     *   vector(i) = i;
     *
     * // do some additions to vector elements, but only for some elements
     * for (unsigned int i=0; i<vector.size(); ++i)
     *   if (some_condition(i) == true)
     *     vector(i) += 1;
     *
     * // do another collective operation
     * const double norm = vector.l2_norm();
     * @endcode
     * 这段代码可能会遇到麻烦：当我们看到第一个加法运算时，我们需要冲刷向量的覆盖缓冲区，而deal.II库会通过调用compress()来实现。然而，它只对所有实际进行加法运算的进程进行这样的操作。
     *
     * 如果其中一个进程的条件永远不为真，那么这个进程就不会得到实际的compress()调用，而所有其他的进程都会。    这就给我们带来了麻烦，因为所有其他的进程都会在刷新写缓冲区的调用中挂起，而另一个进程则会推进到计算l2准则的调用。这时，你会得到一个错误，即某些操作只被一个子集的进程尝试了。这种行为似乎令人惊讶，除非你知道对单个元素的写/添加操作可能会触发这种行为。        这里描述的问题可以通过对compress()放置额外的调用来避免，或者确保所有进程在同一时间做相同类型的操作，例如，如果有必要，可以放置零加法。          @see   @ref GlossGhostedVector  "有鬼魂元素的向量"
     * @ingroup PETScWrappers
     * @ingroup Vectors
     *
     */
    class Vector : public VectorBase
    {
    public:
      /**
       * 声明容器尺寸的类型。
       *
       */
      using size_type = types::global_dof_index;

      /**
       * 默认构造函数。将向量初始化为空。
       *
       */
      Vector();

      /**
       * 构造函数。设置维度为 @p n
       * 并将所有元素初始化为零。              @arg
       * locally_owned_size表示应存储在本进程中的块的大小。
       * @arg
       * communicator表示MPI通信器，向量的不同部分应通过该通信器进行通信。
       * 构造器是明确的，以避免类似的意外。
       * <tt>v=0;</tt>。据推测，用户想把向量的每个元素都设置为0，但结果却是这样的调用。
       * <tt>v=向量 @<number@>(0);</tt>,
       * 即向量被一个长度为零的向量取代。
       *
       */
      explicit Vector(const MPI_Comm &communicator,
                      const size_type n,
                      const size_type locally_owned_size);

      /**
       * 从deal.II向量中复制构造。设置维度为给定向量的维度，并复制所有元素。
       * @arg
       * local_owned_size表示应存储在当前进程中的块的大小。
       * @arg
       * communicator表示MPI通信器，向量的不同部分将通过该通信器通信。
       *
       */
      template <typename Number>
      explicit Vector(const MPI_Comm &              communicator,
                      const dealii::Vector<Number> &v,
                      const size_type               locally_owned_size);


      /**
       * Copy-constructor表示来自PETSc包装向量类的值。
       * @arg  local_size表示应存储在当前进程中的块的大小。
       * @arg
       * communicator表示MPI通信器，矢量的不同部分应通过该通信器进行通信
       * @deprecated
       * 明确使用VectorBase类型的对象已被弃用：使用
       * PETScWrappers::MPI::Vector 代替。
       *
       */
      DEAL_II_DEPRECATED
      explicit Vector(const MPI_Comm &  communicator,
                      const VectorBase &v,
                      const size_type   local_size);

      /**
       * 从IndexSets中构造一个新的并行重影PETSc向量。
       * 注意 @p local 必须是升序和1:1，见
       * IndexSet::is_ascending_and_one_to_one(). 特别是 @p local
       * 中的DoF需要是连续的，这意味着你只能从有几个有限元分量的DoFHandler中创建向量，如果它们没有按分量重新排序（否则使用
       * PETScWrappers::BlockVector ）。
       * 矢量的全局大小由local.size()决定。 @p ghost
       * 中的全局索引是作为鬼魂索引提供的，这样就可以在本地读取它们。
       * 请注意， @p ghost
       * 的IndexSet可能是空的，并且在构造过程中，任何已经包含在
       * @p local
       * 中的指数都会被忽略。这样，鬼魂参数可以等于本地相关自由度的集合，见
       * step-32  。
       * @note  这个操作总是创建一个重影向量，它被认为是只读的。              @see   @ref GlossGhostedVector  "有鬼魂元素的向量"
       *
       */
      Vector(const IndexSet &local,
             const IndexSet &ghost,
             const MPI_Comm &communicator);

      /**
       * 从IndexSet中构造一个新的没有重影元素的并行PETSc向量。
       * 注意 @p local 必须是升序和1:1，见
       * IndexSet::is_ascending_and_one_to_one().  特别是， @p local
       * 中的DoF需要是连续的，这意味着你只能从一个有几个有限元分量的DoFHandler中创建向量，如果它们没有按分量重新排序的话（否则使用
       * PETScWrappers::BlockVector ）。
       *
       */
      explicit Vector(const IndexSet &local, const MPI_Comm &communicator);

      /**
       * 复制构造函数。
       *
       */
      Vector(const Vector &v);

      /**
       * 释放所有内存并返回到与调用默认构造函数后相同的状态。
       *
       */
      virtual void
      clear() override;

      /**
       * 拷贝给定的向量。如果有必要的话，调整当前向量的大小。同时接管
       * @p v. 的MPI通信器。
       *
       */
      Vector &
      operator=(const Vector &v);

      /**
       * 将向量的所有分量设置为给定的数字  @p s.
       * 简单地将这个传给基类，但我们仍然需要声明这个函数，使讨论中给出的关于使构造函数显式的例子生效。
       *
       */
      Vector &
      operator=(const PetscScalar s);

      /**
       * 将一个deal.II向量的值（相对于PETSc向量包装类的值）复制到这个对象中。
       * 与顺序向量的情况相反，这个操作者要求现在的向量已经有正确的大小，因为我们需要有一个分区和一个通讯器存在，否则我们无法从源向量中得到。
       *
       */
      template <typename number>
      Vector &
      operator=(const dealii::Vector<number> &v);

      /**
       * 将向量的尺寸改为 @p N.
       * 没有说明调整向量的尺寸如何影响这个对象的内存分配；也就是说，不能保证将其调整到较小的尺寸实际上也减少了内存消耗，或者为了提高效率，使用相同的内存
       * @p locally_owned_size  表示在本进程中应储存多少 @p N
       * 的值。 对于较少的数据。              @p communicator
       * 表示此后将用于此向量的MPI通信器。            如果 @p
       * omit_zeroing_entries 为假，则向量由零填充。
       * 否则，这些元素将被留作未指定的状态。
       *
       */
      void
      reinit(const MPI_Comm &communicator,
             const size_type N,
             const size_type locally_owned_size,
             const bool      omit_zeroing_entries = false);

      /**
       * 将维度改为向量 @p v,
       * 的维度，同时接管对局部大小的划分，以及MPI通信器的划分。
       * 这与其他 @p reinit 函数同样适用。             @p v
       * 的元素没有被复制，即这个函数与调用<tt>reinit(v.size(),
       * v.local_owned_size(), omit_zeroing_entries)</tt>相同。
       *
       */
      void
      reinit(const Vector &v, const bool omit_zeroing_entries = false);

      /**
       * 重置为一个带有鬼魂元素的向量。更多细节请见相同签名的构造函数。              @see   @ref GlossGhostedVector  "有鬼元素的向量"
       *
       */
      void
      reinit(const IndexSet &local,
             const IndexSet &ghost,
             const MPI_Comm &communicator);

      /**
       * 作为一个没有鬼魂元素的向量重新启动。更多细节请参见具有相同签名的构造函数。              @see   @ref GlossGhostedVector  "有幽灵元素的向量"
       *
       */
      void
      reinit(const IndexSet &local, const MPI_Comm &communicator);

      /**
       * 返回一个对该向量使用的MPI通信器对象的引用。
       *
       */
      const MPI_Comm &
      get_mpi_communicator() const override;

      /**
       * 打印到一个流。  @p precision
       * 表示打印数值所需的精度， @p scientific
       * 是否应使用科学符号。如果 @p across 是 @p true
       * ，那么向量将被打印在一行中，而如果 @p false
       * 则元素被打印在单独的一行中。
       * @note
       * 这个函数重载了基类中的函数，以确保对分布在处理器上的并行向量发生正确的事情。
       *
       */
      void
      print(std::ostream &     out,
            const unsigned int precision  = 3,
            const bool         scientific = true,
            const bool         across     = true) const;

      /**
       * @copydoc   PETScWrappers::VectorBase::all_zero()
       * @note
       * 这个函数重载了基类中的函数，使之成为一个集体操作。
       *
       */
      bool
      all_zero() const;

    protected:
      /**
       * 创建一个长度为 @p n. 的向量
       * 对于这个类，我们创建一个平行向量。  @p n
       * 表示要创建的向量的总大小。  @p
       * local_owned_size表示这些元素中有多少应被存储在本地。
       *
       */
      virtual void
      create_vector(const size_type n, const size_type locally_owned_size);



      /**
       * 创建一个全局长度为 @p n,  本地大小为 @p
       * local_owned_size的向量，并具有指定的ghost索引。注意，在访问这些索引之前，你需要调用update_ghost_values()。
       *
       */
      virtual void
      create_vector(const size_type n,
                    const size_type locally_owned_size,
                    const IndexSet &ghostnodes);


    private:
      /**
       * 将用于该并行向量的通信器对象的副本。
       *
       */
      MPI_Comm communicator;
    };


    // ------------------ template and inline functions -------------


    /**
     * 全局函数 @p swap
     * ，它重载了C++标准库的默认实现，它使用一个临时对象。该函数简单地交换了两个向量的数据。
     * @relatesalso   PETScWrappers::MPI::Vector .
     *
     */
    inline void
    swap(Vector &u, Vector &v)
    {
      u.swap(v);
    }


#    ifndef DOXYGEN

    template <typename number>
    Vector::Vector(const MPI_Comm &              communicator,
                   const dealii::Vector<number> &v,
                   const size_type               locally_owned_size)
      : communicator(communicator)
    {
      Vector::create_vector(v.size(), locally_owned_size);

      *this = v;
    }



    inline Vector &
    Vector::operator=(const PetscScalar s)
    {
      VectorBase::operator=(s);

      return *this;
    }



    template <typename number>
    inline Vector &
    Vector::operator=(const dealii::Vector<number> &v)
    {
      Assert(size() == v.size(), ExcDimensionMismatch(size(), v.size()));

      // FIXME: the following isn't necessarily fast, but this is due to
      // the fact that PETSc doesn't offer an inlined access operator.
      //
      // if someone wants to contribute some code: to make this code
      // faster, one could either first convert all values to PetscScalar,
      // and then set them all at once using VecSetValues. This has the
      // drawback that it could take quite some memory, if the vector is
      // large, and it would in addition allocate memory on the heap, which
      // is expensive. an alternative would be to split the vector into
      // chunks of, say, 128 elements, convert a chunk at a time and set it
      // in the output vector using VecSetValues. since 128 elements is
      // small enough, this could easily be allocated on the stack (as a
      // local variable) which would make the whole thing much more
      // efficient.
      //
      // a second way to make things faster is for the special case that
      // number==PetscScalar. we could then declare a specialization of
      // this template, and omit the conversion. the problem with this is
      // that the best we can do is to use VecSetValues, but this isn't
      // very efficient either: it wants to see an array of indices, which
      // in this case a) again takes up a whole lot of memory on the heap,
      // and b) is totally dumb since its content would simply be the
      // sequence 0,1,2,3,...,n. the best of all worlds would probably be a
      // function in Petsc that would take a pointer to an array of
      // PetscScalar values and simply copy n elements verbatim into the
      // vector...
      for (size_type i = 0; i < v.size(); ++i)
        (*this)(i) = v(i);

      compress(::dealii::VectorOperation::insert);

      return *this;
    }



    inline const MPI_Comm &
    Vector::get_mpi_communicator() const
    {
      return communicator;
    }

#    endif // DOXYGEN
  }        // namespace MPI
} // namespace PETScWrappers

namespace internal
{
  namespace LinearOperatorImplementation
  {
    template <typename>
    class ReinitHelper;

    /**
     * linear_operator.h中内部使用的一个辅助类。对
     * PETScWrappers::MPI::Vector. 的特殊化。
     *
     */
    template <>
    class ReinitHelper<PETScWrappers::MPI::Vector>
    {
    public:
      template <typename Matrix>
      static void
      reinit_range_vector(const Matrix &              matrix,
                          PETScWrappers::MPI::Vector &v,
                          bool  /*omit_zeroing_entries*/ )
      {
        v.reinit(matrix.locally_owned_range_indices(),
                 matrix.get_mpi_communicator());
      }

      template <typename Matrix>
      static void
      reinit_domain_vector(const Matrix &              matrix,
                           PETScWrappers::MPI::Vector &v,
                           bool  /*omit_zeroing_entries*/ )
      {
        v.reinit(matrix.locally_owned_domain_indices(),
                 matrix.get_mpi_communicator());
      }
    };

  } // namespace LinearOperatorImplementation
}  /* namespace internal */ 

 /**@}*/ 


/**
 * 将 dealii::PETScWrappers::MPI::Vector 声明为分布式向量。
 *
 *
 */
template <>
struct is_serial_vector<PETScWrappers::MPI::Vector> : std::false_type
{};


DEAL_II_NAMESPACE_CLOSE

#  endif // DEAL_II_WITH_PETSC

#endif
 /*------------------------- petsc_vector.h -------------------------*/ 


