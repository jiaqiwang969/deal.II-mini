//include/deal.II-translator/lac/petsc_block_sparse_matrix_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2020 by the deal.II authors
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

#ifndef dealii_petsc_block_sparse_matrix_h
#define dealii_petsc_block_sparse_matrix_h


#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_PETSC

#  include <deal.II/lac/block_matrix_base.h>
#  include <deal.II/lac/block_sparsity_pattern.h>
#  include <deal.II/lac/exceptions.h>
#  include <deal.II/lac/petsc_block_vector.h>
#  include <deal.II/lac/petsc_sparse_matrix.h>

#  include <cmath>

DEAL_II_NAMESPACE_OPEN



namespace PETScWrappers
{
  namespace MPI
  {
    /*!   @addtogroup  PETScWrappers  
     * @{     
*
*/

    /**
     * 基于 PETScWrappers::MPI::SparseMatrix
     * 类的阻塞式稀疏矩阵。这个类实现了封锁稀疏矩阵的PETSc
     * SparseMatrix基对象所特有的函数，并将实际工作中对各个块的大部分调用留给基类中实现的函数。关于这个类何时有用的描述，也请参见这里。
     * 与deal.II-type
     * SparseMatrix类相比，PETSc矩阵没有外部对象来表示稀疏性模式。因此，我们不能通过附加一个块状稀疏模式来确定这种类型的块状矩阵的各个块的大小，而是通过调用
     * reinit()
     * 来设置块的数量，然后分别设置每个块的大小。为了固定块矩阵的数据结构，有必要让它知道我们已经改变了基础矩阵的大小。为此，我们必须调用collect_sizes()函数，其原因与BlockSparsityPattern类所记载的大致相同。
     * @ingroup Matrix1 @see   @ref GlossBlockLA  "块（线性代数）"
     *
     */
    class BlockSparseMatrix : public BlockMatrixBase<SparseMatrix>
    {
    public:
      /**
       * 对基类进行类型化定义，以便更简单地访问它自己的别名。
       *
       */
      using BaseClass = BlockMatrixBase<SparseMatrix>;

      /**
       * 对底层矩阵的类型进行类型化定义。
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
       * 构造函数；将矩阵初始化为空，没有任何结构，也就是说，矩阵根本无法使用。因此，这个构造函数只对作为类的成员的矩阵有用。所有其他的矩阵都应该在数据流中的一个点上创建，在那里所有必要的信息都是可用的。
       * 你必须在使用前用reinit(BlockSparsityPattern)初始化矩阵。然后每行和每列的块数由该函数决定。
       *
       */
      BlockSparseMatrix() = default;

      /**
       * 解构器。
       *
       */
      ~BlockSparseMatrix() override = default;

      /**
       * 伪拷贝操作符只拷贝空对象。块状矩阵的大小需要相同。
       *
       */
      BlockSparseMatrix &
      operator=(const BlockSparseMatrix &);

      /**
       * 这个操作符将一个标量分配给一个矩阵。因为这通常没有什么意义（我们应该把所有的矩阵条目都设置为这个值吗？
       * 仅仅是稀疏模式的非零条目？），这个操作只允许在实际要分配的值为零的情况下进行。这个操作符的存在只是为了允许明显的符号<tt>matrix=0</tt>，它将矩阵的所有元素设置为零，但保留之前使用的稀疏模式。
       *
       */
      BlockSparseMatrix &
      operator=(const double d);

      /**
       * 调整矩阵的大小，通过设置块的行数和列数。
       * 这将删除所有的块，并用未初始化的块代替，也就是那些尚未设置大小的块。你必须通过调用块本身的
       * @p reinit
       * 函数来做到这一点。不要忘了之后在这个对象上调用collect_sizes()。
       * 你必须自己设置块的大小的原因是，大小可能是变化的，每行的最大元素数可能是变化的，等等。在这里不复制SparsityPattern类的接口是比较简单的，而是让用户调用他们想要的任何函数。
       *
       */
      void
      reinit(const size_type n_block_rows, const size_type n_block_columns);


      /**
       * 有效地重新引用块状矩阵进行并行计算。只有简单类型的BlockSparsityPattern可以有效地并行存储大型稀疏模式，所以这是唯一支持的参数。IndexSets描述了每个块的本地拥有的DoF的范围。注意，IndexSets需要升序和1:1。
       * 对于一个对称的结构，前两个参数使用同一个向量。
       *
       */
      void
      reinit(const std::vector<IndexSet> &      rows,
             const std::vector<IndexSet> &      cols,
             const BlockDynamicSparsityPattern &bdsp,
             const MPI_Comm &                   com);


      /**
       * 与上述相同，但只针对对称结构。
       *
       */
      void
      reinit(const std::vector<IndexSet> &      sizes,
             const BlockDynamicSparsityPattern &bdsp,
             const MPI_Comm &                   com);



      /**
       * 矩阵-向量乘法：让 $dst = M*src$ 与 $M$ 为该矩阵。
       *
       */
      void
      vmult(BlockVector &dst, const BlockVector &src) const;

      /**
       * 矩阵-向量乘法。就像前面的函数一样，但只适用于矩阵只有一个块列的情况。
       *
       */
      void
      vmult(BlockVector &dst, const Vector &src) const;

      /**
       * 矩阵-向量乘法。就像前面的函数，但只适用于矩阵只有一个块行的情况。
       *
       */
      void
      vmult(Vector &dst, const BlockVector &src) const;

      /**
       * 矩阵-向量乘法。就像前面的函数，但只适用于矩阵只有一个块的情况。
       *
       */
      void
      vmult(Vector &dst, const Vector &src) const;

      /**
       * 矩阵-向量乘法：让 $dst = M^T*src$ 与 $M$
       * 为这个矩阵。这个函数与vmult()的作用相同，但需要转置的矩阵。
       *
       */
      void
      Tvmult(BlockVector &dst, const BlockVector &src) const;

      /**
       * 矩阵-向量乘法。就像前面的函数一样，但只适用于矩阵只有一个块行的情况。
       *
       */
      void
      Tvmult(BlockVector &dst, const Vector &src) const;

      /**
       * 矩阵-向量乘法。就像前面的函数一样，但只适用于矩阵只有一个块列的情况。
       *
       */
      void
      Tvmult(Vector &dst, const BlockVector &src) const;

      /**
       * 矩阵-向量乘法。就像前面的函数，但只适用于矩阵只有一个块的情况。
       *
       */
      void
      Tvmult(Vector &dst, const Vector &src) const;

      /**
       * 这个函数收集了子对象的大小，并将其存储在内部数组中，以便能够将矩阵的全局索引转为子对象的索引。在你改变子对象的大小后，你必须*每次都调用这个函数。
       *
       */
      void
      collect_sizes();

      /**
       * 返回该矩阵的域空间的划分，即该矩阵必须与之相乘的向量的划分。
       *
       */
      std::vector<IndexSet>
      locally_owned_domain_indices() const;

      /**
       * 返回该矩阵的范围空间的划分，即由矩阵-向量乘积产生的向量的划分。
       *
       */
      std::vector<IndexSet>
      locally_owned_range_indices() const;

      /**
       * 返回对与该矩阵一起使用的MPI通信器对象的一个引用。
       *
       */
      const MPI_Comm &
      get_mpi_communicator() const;

      /**
       * 使基类中的clear()函数可见，尽管它是受保护的。
       *
       */
      using BlockMatrixBase<SparseMatrix>::clear;
    };



     /*@}*/ 

    // ------------- inline and template functions -----------------

    inline BlockSparseMatrix &
    BlockSparseMatrix::operator=(const double d)
    {
      Assert(d == 0, ExcScalarAssignmentOnlyForZeroValue());

      for (size_type r = 0; r < this->n_block_rows(); ++r)
        for (size_type c = 0; c < this->n_block_cols(); ++c)
          this->block(r, c) = d;

      return *this;
    }



    inline void
    BlockSparseMatrix::vmult(BlockVector &dst, const BlockVector &src) const
    {
      BaseClass::vmult_block_block(dst, src);
    }



    inline void
    BlockSparseMatrix::vmult(BlockVector &dst, const Vector &src) const
    {
      BaseClass::vmult_block_nonblock(dst, src);
    }



    inline void
    BlockSparseMatrix::vmult(Vector &dst, const BlockVector &src) const
    {
      BaseClass::vmult_nonblock_block(dst, src);
    }



    inline void
    BlockSparseMatrix::vmult(Vector &dst, const Vector &src) const
    {
      BaseClass::vmult_nonblock_nonblock(dst, src);
    }


    inline void
    BlockSparseMatrix::Tvmult(BlockVector &dst, const BlockVector &src) const
    {
      BaseClass::Tvmult_block_block(dst, src);
    }



    inline void
    BlockSparseMatrix::Tvmult(BlockVector &dst, const Vector &src) const
    {
      BaseClass::Tvmult_block_nonblock(dst, src);
    }



    inline void
    BlockSparseMatrix::Tvmult(Vector &dst, const BlockVector &src) const
    {
      BaseClass::Tvmult_nonblock_block(dst, src);
    }



    inline void
    BlockSparseMatrix::Tvmult(Vector &dst, const Vector &src) const
    {
      BaseClass::Tvmult_nonblock_nonblock(dst, src);
    }

  } // namespace MPI

} // namespace PETScWrappers


DEAL_II_NAMESPACE_CLOSE


#endif // DEAL_II_WITH_PETSC

#endif // dealii_petsc_block_sparse_matrix_h


