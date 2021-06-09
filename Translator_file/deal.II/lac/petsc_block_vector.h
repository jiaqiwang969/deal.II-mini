//include/deal.II-translator/lac/petsc_block_vector_0.txt
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

#ifndef dealii_petsc_block_vector_h
#define dealii_petsc_block_vector_h


#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_PETSC

#  include <deal.II/lac/block_indices.h>
#  include <deal.II/lac/block_vector_base.h>
#  include <deal.II/lac/exceptions.h>
#  include <deal.II/lac/petsc_vector.h>
#  include <deal.II/lac/vector_type_traits.h>

DEAL_II_NAMESPACE_OPEN


namespace PETScWrappers
{
  // forward declaration
  class BlockVector;

  namespace MPI
  {
    /*!   @addtogroup  PETScWrappers  
     * @{     
*
*/

    /**
     * 基于PETScWrappers中实现的平行向量类的块向量的实现。虽然基类提供了大部分的接口，但这个类处理了向量的实际分配，并提供了底层向量类型的特定函数。
     * 数据的分配模式是这样的：每个块都分布在MPI通信器中命名的所有MPI进程中。
     * 也就是说，我们不只是分发整个向量，而是分发每个组件。因此，在构造函数和reinit()函数中，不仅要指定各个块的大小，还要指定这些块中每个元素在本地进程中的存储数量。
     * @ingroup Vectors @see   @ref GlossBlockLA  "块（线性代数）"
     *
     */
    class BlockVector : public BlockVectorBase<Vector>
    {
    public:
      /**
       * 对基类进行类型化定义，以便更简单地访问它自己的别名。
       *
       */
      using BaseClass = BlockVectorBase<Vector>;

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
       * 构造函数。生成一个有 @p n_blocks
       * 个块的向量，每个块都是跨越 @p communicator
       * 的平行向量，有 @p block_size 个元素，其中 @p
       * locally_owned_size 个元素被存储在当前进程中。
       *
       */
      explicit BlockVector(const unsigned int n_blocks,
                           const MPI_Comm &   communicator,
                           const size_type    block_size,
                           const size_type    locally_owned_size);

      /**
       * 复制构造函数。将平行向量的所有属性设置为给定参数的属性，并复制这些元素。
       *
       */
      BlockVector(const BlockVector &V);

      /**
       * 构造函数。将块的数量设置为<tt>block_sizes.size()</tt>，并以<tt>block_sizes[i]</tt>零元素初始化每个块。
       * 各个块分布在给定的通信器上，并且每个块都在当前进程上存储<tt>local_elements[i]</tt>元素。
       *
       */
      BlockVector(const std::vector<size_type> &block_sizes,
                  const MPI_Comm &              communicator,
                  const std::vector<size_type> &local_elements);

      /**
       * 用parallel_partitioning.size()块创建一个BlockVector，每个块用给定的IndexSet初始化。
       *
       */
      explicit BlockVector(const std::vector<IndexSet> &parallel_partitioning,
                           const MPI_Comm &communicator = MPI_COMM_WORLD);

      /**
       * 与上述相同，但包括幽灵元素
       *
       */
      BlockVector(const std::vector<IndexSet> &parallel_partitioning,
                  const std::vector<IndexSet> &ghost_indices,
                  const MPI_Comm &             communicator);



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
      operator=(const BlockVector &V);

      /**
       * 重新初始化BlockVector，使其包含大小为 @p block_size的
       * @p n_blocks ，每个都在本地存储 @p locally_owned_size
       * 元素。 @p communicator
       * 参数表示这些块中的每个块应与哪个MPI通道通信。
       * 如果<tt>omit_zeroing_entries==false</tt>，该向量将被填充为零。
       *
       */
      void
      reinit(const unsigned int n_blocks,
             const MPI_Comm &   communicator,
             const size_type    block_size,
             const size_type    locally_owned_size,
             const bool         omit_zeroing_entries = false);

      /**
       * 重新初始化BlockVector，使其包含<tt>block_sizes.size()</tt>块。每个块都被重新初始化为<tt>block_sizes[i]</tt>尺寸。它们中的每一个都在当前进程中存储<tt>locally_owned_sizes[i]</tt>元素。
       * 如果块的数量与调用此函数前相同，所有的向量都保持不变，并为每个向量调用
       * reinit()。
       * 如果<tt>omit_zeroing_entries==false</tt>，则向量被填充为零。
       * 注意，你必须调用这个（或其他reinit()函数）函数，而不是调用单个块的reinit()函数，以允许块向量更新其缓存的向量大小。如果你调用其中一个块的reinit()，那么对这个对象的后续操作可能会产生不可预测的结果，因为它们可能被路由到错误的块。
       *
       */
      void
      reinit(const std::vector<size_type> &block_sizes,
             const MPI_Comm &              communicator,
             const std::vector<size_type> &locally_owned_sizes,
             const bool                    omit_zeroing_entries = false);

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
       * 使用IndexSets重新初始化BlockVector。详细情况请参见具有相同参数的构造函数。
       *
       */
      void
      reinit(const std::vector<IndexSet> &parallel_partitioning,
             const MPI_Comm &             communicator);

      /**
       * 与上述相同，但包括幽灵条目。
       *
       */
      void
      reinit(const std::vector<IndexSet> &parallel_partitioning,
             const std::vector<IndexSet> &ghost_entries,
             const MPI_Comm &             communicator);

      /**
       * 将块的数量改为<tt>num_blocks</tt>。各个区块将被初始化为零大小，所以假定用户自己以适当的方式调整各个区块的大小，并在之后调用<tt>collect_sizes</tt>。
       *
       */
      void
      reinit(const unsigned int num_blocks);

      /**
       * 如果这个向量是一个重影向量（因此是只读的），则返回。
       *
       */
      bool
      has_ghost_elements() const;

      /**
       * 返回一个对与此向量一起使用的MPI通信器对象的引用。
       *
       */
      const MPI_Comm &
      get_mpi_communicator() const;

      /**
       * 交换这个向量和另一个向量<tt>v</tt>的内容。我们可以用一个临时变量和复制数据元素来完成这个操作，但是这个函数明显更有效率，因为它只交换两个向量的数据指针，因此不需要分配临时存储和移动数据。
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

     /*@}*/ 

     /*--------------------- Inline functions --------------------------------*/ 

    inline BlockVector::BlockVector(const unsigned int n_blocks,
                                    const MPI_Comm &   communicator,
                                    const size_type    block_size,
                                    const size_type    locally_owned_size)
    {
      reinit(n_blocks, communicator, block_size, locally_owned_size);
    }



    inline BlockVector::BlockVector(
      const std::vector<size_type> &block_sizes,
      const MPI_Comm &              communicator,
      const std::vector<size_type> &local_elements)
    {
      reinit(block_sizes, communicator, local_elements, false);
    }


    inline BlockVector::BlockVector(const BlockVector &v)
      : BlockVectorBase<Vector>()
    {
      this->components.resize(v.n_blocks());
      this->block_indices = v.block_indices;

      for (unsigned int i = 0; i < this->n_blocks(); ++i)
        this->components[i] = v.components[i];
    }

    inline BlockVector::BlockVector(
      const std::vector<IndexSet> &parallel_partitioning,
      const MPI_Comm &             communicator)
    {
      reinit(parallel_partitioning, communicator);
    }

    inline BlockVector::BlockVector(
      const std::vector<IndexSet> &parallel_partitioning,
      const std::vector<IndexSet> &ghost_indices,
      const MPI_Comm &             communicator)
    {
      reinit(parallel_partitioning, ghost_indices, communicator);
    }

    inline BlockVector &
    BlockVector::operator=(const value_type s)
    {
      BaseClass::operator=(s);
      return *this;
    }

    inline BlockVector &
    BlockVector::operator=(const BlockVector &v)
    {
      // we only allow assignment to vectors with the same number of blocks
      // or to an empty BlockVector
      Assert(n_blocks() == 0 || n_blocks() == v.n_blocks(),
             ExcDimensionMismatch(n_blocks(), v.n_blocks()));

      if (this->n_blocks() != v.n_blocks())
        reinit(v.n_blocks());

      for (size_type i = 0; i < this->n_blocks(); ++i)
        this->components[i] = v.block(i);

      collect_sizes();

      return *this;
    }



    inline void
    BlockVector::reinit(const unsigned int n_blocks,
                        const MPI_Comm &   communicator,
                        const size_type    block_size,
                        const size_type    locally_owned_size,
                        const bool         omit_zeroing_entries)
    {
      reinit(std::vector<size_type>(n_blocks, block_size),
             communicator,
             std::vector<size_type>(n_blocks, locally_owned_size),
             omit_zeroing_entries);
    }



    inline void
    BlockVector::reinit(const std::vector<size_type> &block_sizes,
                        const MPI_Comm &              communicator,
                        const std::vector<size_type> &locally_owned_sizes,
                        const bool                    omit_zeroing_entries)
    {
      this->block_indices.reinit(block_sizes);
      if (this->components.size() != this->n_blocks())
        this->components.resize(this->n_blocks());

      for (unsigned int i = 0; i < this->n_blocks(); ++i)
        this->components[i].reinit(communicator,
                                   block_sizes[i],
                                   locally_owned_sizes[i],
                                   omit_zeroing_entries);
    }


    inline void
    BlockVector::reinit(const BlockVector &v, const bool omit_zeroing_entries)
    {
      this->block_indices = v.get_block_indices();
      if (this->components.size() != this->n_blocks())
        this->components.resize(this->n_blocks());

      for (unsigned int i = 0; i < this->n_blocks(); ++i)
        block(i).reinit(v.block(i), omit_zeroing_entries);
    }

    inline void
    BlockVector::reinit(const std::vector<IndexSet> &parallel_partitioning,
                        const MPI_Comm &             communicator)
    {
      std::vector<size_type> sizes(parallel_partitioning.size());
      for (unsigned int i = 0; i < parallel_partitioning.size(); ++i)
        sizes[i] = parallel_partitioning[i].size();

      this->block_indices.reinit(sizes);
      if (this->components.size() != this->n_blocks())
        this->components.resize(this->n_blocks());

      for (unsigned int i = 0; i < this->n_blocks(); ++i)
        block(i).reinit(parallel_partitioning[i], communicator);
    }

    inline void
    BlockVector::reinit(const std::vector<IndexSet> &parallel_partitioning,
                        const std::vector<IndexSet> &ghost_entries,
                        const MPI_Comm &             communicator)
    {
      std::vector<types::global_dof_index> sizes(parallel_partitioning.size());
      for (unsigned int i = 0; i < parallel_partitioning.size(); ++i)
        sizes[i] = parallel_partitioning[i].size();

      this->block_indices.reinit(sizes);
      if (this->components.size() != this->n_blocks())
        this->components.resize(this->n_blocks());

      for (unsigned int i = 0; i < this->n_blocks(); ++i)
        block(i).reinit(parallel_partitioning[i],
                        ghost_entries[i],
                        communicator);
    }



    inline const MPI_Comm &
    BlockVector::get_mpi_communicator() const
    {
      return block(0).get_mpi_communicator();
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

      ::dealii::swap(this->block_indices, v.block_indices);
    }



    inline void
    BlockVector::print(std::ostream &     out,
                       const unsigned int precision,
                       const bool         scientific,
                       const bool         across) const
    {
      for (unsigned int i = 0; i < this->n_blocks(); ++i)
        {
          if (across)
            out << 'C' << i << ':';
          else
            out << "Component " << i << std::endl;
          this->components[i].print(out, precision, scientific, across);
        }
    }



    /**
     * 全局函数，它重载了C++标准库的默认实现，它使用一个临时对象。该函数简单地交换了两个向量的数据。
     * @relatesalso   PETScWrappers::MPI::BlockVector .
     *
     */
    inline void
    swap(BlockVector &u, BlockVector &v)
    {
      u.swap(v);
    }

  } // namespace MPI

} // namespace PETScWrappers

namespace internal
{
  namespace LinearOperatorImplementation
  {
    template <typename>
    class ReinitHelper;

    /**
     * linear_operator.h中内部使用的一个辅助类。对
     * PETScWrappers::MPI::BlockVector. 的特殊化。
     *
     */
    template <>
    class ReinitHelper<PETScWrappers::MPI::BlockVector>
    {
    public:
      template <typename Matrix>
      static void
      reinit_range_vector(const Matrix &                   matrix,
                          PETScWrappers::MPI::BlockVector &v,
                          bool  /*omit_zeroing_entries*/ )
      {
        v.reinit(matrix.locally_owned_range_indices(),
                 matrix.get_mpi_communicator());
      }

      template <typename Matrix>
      static void
      reinit_domain_vector(const Matrix &                   matrix,
                           PETScWrappers::MPI::BlockVector &v,
                           bool  /*omit_zeroing_entries*/ )
      {
        v.reinit(matrix.locally_owned_domain_indices(),
                 matrix.get_mpi_communicator());
      }
    };

  } // namespace LinearOperatorImplementation
}  /* namespace internal */ 


/**
 * 将 dealii::PETScWrappers::MPI::BlockVector 声明为分布式向量。
 *
 *
 */
template <>
struct is_serial_vector<PETScWrappers::MPI::BlockVector> : std::false_type
{};


DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_PETSC

#endif


