//include/deal.II-translator/lac/trilinos_parallel_block_vector_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2020 by the deal.II authors
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

#ifndef dealii_trilinos_parallel_block_vector_h
#define dealii_trilinos_parallel_block_vector_h


#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_TRILINOS

#  include <deal.II/lac/block_indices.h>
#  include <deal.II/lac/block_vector_base.h>
#  include <deal.II/lac/exceptions.h>
#  include <deal.II/lac/trilinos_vector.h>

#  include <functional>

DEAL_II_NAMESPACE_OPEN

// forward declaration
#  ifndef DOXYGEN
template <typename Number>
class BlockVectorBase;

namespace TrilinosWrappers
{
  // forward declaration
  namespace MPI
  {
    class BlockVector;
  }
  class BlockSparseMatrix;
} // namespace TrilinosWrappers
#  endif

/*!   @addtogroup  TrilinosWrappers  
     * @{ 

 
*
*/

namespace TrilinosWrappers
{
  namespace MPI
  {
    /**
     * 一个基于TrilinosWrappers中实现的向量类的块向量实现。虽然基类提供了大部分的接口，但这个类处理了向量的实际分配，并提供了底层向量类型特有的函数。
     * 数据的分配模式是这样的，每个块都分布在MPI通信器中命名的所有MPI进程中。
     * 也就是说，我们不只是分发整个向量，而是分发每个组件。因此，在构造函数和reinit()函数中，不仅要指定各个块的大小，还要指定这些块中每个元素在本地进程中的存储数量。
     * @ingroup Vectors
     * @ingroup TrilinosWrappers @see   @ref GlossBlockLA  "块（线性代数）"
     *
     */
    class BlockVector : public dealii::BlockVectorBase<MPI::Vector>
    {
    public:
      /**
       * 对基类进行类型化定义，以便更简单地访问它自己的别名。
       *
       */
      using BaseClass = dealii::BlockVectorBase<MPI::Vector>;

      /**
       * 类型化底层向量的类型。
       *
       */
      using BlockType = BaseClass::BlockType;

      /**
       * 从基类中导入别名。
       *
       */
      using value_type      = BaseClass::value_type;
      using pointer         = BaseClass::pointer;
      using const_pointer   = BaseClass::const_pointer;
      using reference       = BaseClass::reference;
      using const_reference = BaseClass::const_reference;
      using size_type       = BaseClass::size_type;
      using iterator        = BaseClass::iterator;
      using const_iterator  = BaseClass::const_iterator;

      /**
       * 默认构造函数。生成一个没有任何块的空向量。
       *
       */
      BlockVector() = default;

      /**
       * 构造函数。生成一个块向量，其块数与 @p partitioning.
       * 中的条目一样多。每个IndexSet与MPI通信器一起包含MPI进程之间的数据分配布局。
       *
       */
      explicit BlockVector(const std::vector<IndexSet> &parallel_partitioning,
                           const MPI_Comm &communicator = MPI_COMM_WORLD);

      /**
       * 创建一个带有幽灵元素的BlockVector。更多细节见各自的reinit()方法。
       * @p ghost_values 可以包含 @p parallel_partitioning,
       * 中的任何元素，它们将被忽略。
       *
       */
      BlockVector(const std::vector<IndexSet> &parallel_partitioning,
                  const std::vector<IndexSet> &ghost_values,
                  const MPI_Comm &             communicator,
                  const bool                   vector_writable = false);

      /**
       * 复制-构造器。将平行向量的所有属性设置为给定参数的属性，并复制这些元素。
       *
       */
      BlockVector(const BlockVector &v);

      /**
       * 移动构造函数。通过窃取向量的内部数据创建一个新的向量
       * @p v. 。
       *
       */
      BlockVector(BlockVector &&v) noexcept;

      /**
       * 创建一个由<tt>num_blocks</tt>组件组成的块状向量，但各个组件中没有内容，用户必须使用块状向量的重定位来填充适当的数据。
       *
       */
      explicit BlockVector(const size_type num_blocks);

      /**
       * 销毁器。清除内存
       *
       */
      ~BlockVector() override = default;

      /**
       * 复制操作符：用给定的标量值填充本地存储的向量的所有组件。
       *
       */
      BlockVector &
      operator=(const value_type s);

      /**
       * 对相同类型的参数进行复制操作。
       *
       */
      BlockVector &
      operator=(const BlockVector &v);

      /**
       * 移动给定的向量。这个操作符通过有效地交换内部数据结构，将目前的向量替换为
       * @p v 。
       *
       */
      BlockVector &
      operator=(BlockVector &&v) noexcept;

      /**
       * 另一个复制函数。这个函数接收一个deal.II块向量并将其复制到一个TrilinosWrappers块向量中。注意，块的数量必须与输入向量中的相同。
       * 使用 reinit() 命令来调整 BlockVector
       * 的大小或改变块组件的内部结构。
       * 由于Trilinos只在双数上工作，这个函数被限制为在deal.II向量中只接受一种可能的数字类型。
       *
       */
      template <typename Number>
      BlockVector &
      operator=(const ::dealii::BlockVector<Number> &v);

      /**
       * 重新初始化BlockVector，使其包含与输入参数中给出的索引集一样多的块，根据地图中描述的各个组件的平行分布。
       * 如果<tt>omit_zeroing_entries==false</tt>，该向量将被填充为零。
       *
       */
      void
      reinit(const std::vector<IndexSet> &parallel_partitioning,
             const MPI_Comm &             communicator         = MPI_COMM_WORLD,
             const bool                   omit_zeroing_entries = false);

      /**
       * Reinit功能。这个函数销毁了旧的向量内容，并根据输入的分区生成一个新的向量。除了像上面所有其他方法一样只指定一个索引集外，这个方法还允许提供一个额外的幽灵条目集。
       * 有两个不同版本的向量可以被创建。如果标志 @p
       * vector_writable 被设置为 @p false, ，该向量只允许读取 @p
       * parallel_partitioning 和 @p ghost_entries.
       * 的联合集合，那么reinit方法的效果相当于调用其他的reinit方法，其索引集包含本地拥有的条目和幽灵条目。
       * 如果标志 @p vector_writable 被设置为
       * "true"，这就为ghost元素创建了一个替代性的存储方案，允许多个线程向向量中写入数据（对于其他reinit方法，一次只允许一个线程向ghost条目写入数据）。
       *
       */
      void
      reinit(const std::vector<IndexSet> &partitioning,
             const std::vector<IndexSet> &ghost_values,
             const MPI_Comm &             communicator    = MPI_COMM_WORLD,
             const bool                   vector_writable = false);


      /**
       * 将尺寸改为向量<tt>V</tt>的尺寸。这与其他reinit()函数的情况相同。
       * <tt>V</tt>的元素不会被复制，也就是说，这个函数与调用<tt>reinit
       * (V.size(), omit_zeroing_entries)</tt>相同。
       * 注意，你必须调用这个（或其他reinit()函数）函数，而不是调用单个块的reinit()函数，以允许块向量更新它的向量大小缓存。如果你在其中一个块上调用reinit()，那么对这个对象的后续操作可能会产生不可预测的结果，因为它们可能会被路由到错误的块上。
       *
       */
      void
      reinit(const BlockVector &V, const bool omit_zeroing_entries = false);

      /**
       * 将块的数量改为<tt>num_blocks</tt>。各个区块将被初始化为零大小，所以假定用户自己以适当的方式调整各个区块的大小，并在之后调用<tt>collect_sizes</tt>。
       *
       */
      void
      reinit(const size_type num_blocks);

      /**
       * 这个reinit函数是为了用于需要使用一些非本地数据的并行计算。需要这个函数的典型情况是并行调用
       * FEValues<dim>::get_function_values
       * 函数（或一些导数）。由于提前检索数据通常更快，这个函数可以在汇编分叉到不同处理器之前被调用。这个函数的作用如下。
       * 它获取给定矩阵的列中的信息，并寻找不同处理器之间的数据耦合。然后从输入向量中查询这些数据。注意，你不应该再向结果向量写入数据，因为有些数据可以在不同的处理器上存储多次，导致不可预测的结果。特别是，这样的向量不能用于矩阵-向量乘积，例如在线性系统的求解过程中。
       *
       */
      void
      import_nonlocal_data_for_fe(const TrilinosWrappers::BlockSparseMatrix &m,
                                  const BlockVector &                        v);

      /**
       * 如果这个向量包含鬼魂元素，则返回。              @see   @ref GlossGhostedVector  "有鬼元素的向量"
       *
       */
      bool
      has_ghost_elements() const;

      /**
       * 交换这个向量和另一个向量<tt>v</tt>的内容。我们可以用一个临时变量和复制数据元素来完成这个操作，但这个函数明显更有效率，因为它只交换两个向量的数据指针，因此不需要分配临时存储空间和移动数据。
       * 限制：现在这个函数只在两个向量有相同数量的块的情况下工作。如果需要，也应该交换块的数量。
       * 这个函数类似于所有C++标准容器的swap()函数。此外，还有一个全局函数swap(u,v)，它简单地调用<tt>u.swap(v)</tt>，同样与标准函数类比。
       *
       */
      void
      swap(BlockVector &v);

      /**
       * 打印到一个流。
       *
       */
      void
      print(std::ostream &     out,
            const unsigned int precision  = 3,
            const bool         scientific = true,
            const bool         across     = true) const;

      /**
       * 异常情况
       *
       */
      DeclException0(ExcIteratorRangeDoesNotMatchVectorSize);

      /**
       * 异常情况
       *
       */
      DeclException0(ExcNonMatchingBlockVectors);
    };



     /*-------------------------- Inline functions ---------------------------*/ 
    inline BlockVector::BlockVector(
      const std::vector<IndexSet> &parallel_partitioning,
      const MPI_Comm &             communicator)
    {
      reinit(parallel_partitioning, communicator, false);
    }



    inline BlockVector::BlockVector(
      const std::vector<IndexSet> &parallel_partitioning,
      const std::vector<IndexSet> &ghost_values,
      const MPI_Comm &             communicator,
      const bool                   vector_writable)
    {
      reinit(parallel_partitioning,
             ghost_values,
             communicator,
             vector_writable);
    }



    inline BlockVector::BlockVector(const size_type num_blocks)
    {
      reinit(num_blocks);
    }



    inline BlockVector::BlockVector(const BlockVector &v)
      : dealii::BlockVectorBase<MPI::Vector>()
    {
      this->components.resize(v.n_blocks());
      this->block_indices = v.block_indices;

      for (size_type i = 0; i < this->n_blocks(); ++i)
        this->components[i] = v.components[i];
    }



    inline BlockVector::BlockVector(BlockVector &&v) noexcept
    {
      // initialize a minimal, valid object and swap
      reinit(0);
      swap(v);
    }



    template <typename Number>
    BlockVector &
    BlockVector::operator=(const ::dealii::BlockVector<Number> &v)
    {
      if (n_blocks() != v.n_blocks())
        {
          std::vector<size_type> block_sizes(v.n_blocks(), 0);
          block_indices.reinit(block_sizes);
          if (components.size() != n_blocks())
            components.resize(n_blocks());
        }

      for (size_type i = 0; i < this->n_blocks(); ++i)
        this->components[i] = v.block(i);

      collect_sizes();

      return *this;
    }



    inline bool
    BlockVector::has_ghost_elements() const
    {
      bool ghosted = block(0).has_ghost_elements();
#  ifdef DEBUG
      for (unsigned int i = 0; i < this->n_blocks(); ++i)
        Assert(block(i).has_ghost_elements() == ghosted, ExcInternalError());
#  endif
      return ghosted;
    }



    inline void
    BlockVector::swap(BlockVector &v)
    {
      std::swap(this->components, v.components);

      dealii::swap(this->block_indices, v.block_indices);
    }



    /**
     * 全局函数，它重载了C++标准库的默认实现，它使用一个临时对象。该函数简单地交换了两个向量的数据。
     * @relatesalso   TrilinosWrappers::MPI::BlockVector .
     *
     */
    inline void
    swap(BlockVector &u, BlockVector &v)
    {
      u.swap(v);
    }

  }  /* namespace MPI */ 

}  /* namespace TrilinosWrappers */ 

 /*@}*/ 


namespace internal
{
  namespace LinearOperatorImplementation
  {
    template <typename>
    class ReinitHelper;

    /**
     * linear_operator.h中内部使用的一个辅助类。对
     * TrilinosWrappers::MPI::BlockVector. 的特殊化。
     *
     */
    template <>
    class ReinitHelper<TrilinosWrappers::MPI::BlockVector>
    {
    public:
      template <typename Matrix>
      static void
      reinit_range_vector(const Matrix &                      matrix,
                          TrilinosWrappers::MPI::BlockVector &v,
                          bool omit_zeroing_entries)
      {
        v.reinit(matrix.locally_owned_range_indices(),
                 matrix.get_mpi_communicator(),
                 omit_zeroing_entries);
      }

      template <typename Matrix>
      static void
      reinit_domain_vector(const Matrix &                      matrix,
                           TrilinosWrappers::MPI::BlockVector &v,
                           bool omit_zeroing_entries)
      {
        v.reinit(matrix.locally_owned_domain_indices(),
                 matrix.get_mpi_communicator(),
                 omit_zeroing_entries);
      }
    };

  } // namespace LinearOperatorImplementation
}  /* namespace internal */ 


/**
 * 将 dealii::TrilinosWrappers::MPI::BlockVector 声明为分布式向量。
 *
 *
 */
template <>
struct is_serial_vector<TrilinosWrappers::MPI::BlockVector> : std::false_type
{};

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_TRILINOS

#endif


