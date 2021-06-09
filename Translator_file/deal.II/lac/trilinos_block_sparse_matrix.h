//include/deal.II-translator/lac/trilinos_block_sparse_matrix_0.txt
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

#ifndef dealii_trilinos_block_sparse_matrix_h
#define dealii_trilinos_block_sparse_matrix_h


#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_TRILINOS

#  include <deal.II/base/template_constraints.h>

#  include <deal.II/lac/block_matrix_base.h>
#  include <deal.II/lac/exceptions.h>
#  include <deal.II/lac/full_matrix.h>
#  include <deal.II/lac/trilinos_parallel_block_vector.h>
#  include <deal.II/lac/trilinos_sparse_matrix.h>

#  include <cmath>

DEAL_II_NAMESPACE_OPEN

// forward declarations
#  ifndef DOXYGEN
class BlockSparsityPattern;
template <typename number>
class BlockSparseMatrix;
#  endif

namespace TrilinosWrappers
{
  /*!   @addtogroup  TrilinosWrappers  
     * @{   
*
*/

  /**
   * 基于 TrilinosWrappers::SparseMatrix 类的阻塞式稀疏矩阵。
   * 这个类实现了Trilinos
   * SparseMatrix基对象的特定函数，用于阻塞式稀疏矩阵，并将实际工作中对各个块的大部分调用留给基类中实现的函数。关于这个类何时有用的描述，也请参见那里。
   * 与deal.II-type
   * SparseMatrix类相比，Trilinos矩阵没有外部对象来表示稀疏性模式。因此，人们并不是通过附加一个块状稀疏模式来决定这种类型的块状矩阵的各个块的大小，而是通过调用
   * reinit()
   * 来设置块的数量，然后再分别设置每个块的大小。为了固定块矩阵的数据结构，有必要让它知道我们已经改变了基础矩阵的大小。为此，我们必须调用collect_sizes()函数，其原因与BlockSparsityPattern类所记载的大致相同。
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
    ~BlockSparseMatrix() override;

    /**
     * 伪拷贝操作符只拷贝空对象。块状矩阵的大小需要相同。
     *
     */
    BlockSparseMatrix &
    operator=(const BlockSparseMatrix &) = default;

    /**
     * 这个操作符将一个标量分配给一个矩阵。因为这通常没有什么意义（我们应该把所有的矩阵条目都设置为这个值吗？仅仅是稀疏模式的非零条目？），这个操作只允许在实际要分配的值为零的情况下进行。这个操作符的存在只是为了允许明显的符号<tt>matrix=0</tt>，它将矩阵的所有元素设置为零，但保留之前使用的稀疏模式。
     *
     */
    BlockSparseMatrix &
    operator=(const double d);

    /**
     * 调整矩阵的大小，通过设置块的行数和列数。
     * 这将删除所有的块，并用未初始化的块代替，也就是那些尚未设置大小的块。你必须通过调用块本身的
     * @p reinit
     * 函数来做到这一点。不要忘了之后在这个对象上调用collect_sizes()。
     * 你必须自己设置块的大小的原因是，大小可能是变化的，每行的最大元素数可能是变化的，等等。在这里不复制
     * @p
     * SparsityPattern类的接口比较简单，而是让用户调用他们想要的任何函数。
     *
     */
    void
    reinit(const size_type n_block_rows, const size_type n_block_columns);

    /**
     * 调整矩阵的大小，通过使用一个索引集阵列来确定各个矩阵的%平行分布。这个函数假定生成了一个二次方块矩阵。
     *
     */
    template <typename BlockSparsityPatternType>
    void
    reinit(const std::vector<IndexSet> &   input_maps,
           const BlockSparsityPatternType &block_sparsity_pattern,
           const MPI_Comm &                communicator  = MPI_COMM_WORLD,
           const bool                      exchange_data = false);

    /**
     * 调整矩阵的大小，并通过给定的稀疏模式进行初始化。
     * 由于没有给出分布图，结果是一个块状矩阵，所有的元素都存储在本地。
     *
     */
    template <typename BlockSparsityPatternType>
    void
    reinit(const BlockSparsityPatternType &block_sparsity_pattern);

    /**
     * 这个函数使用deal.II稀疏矩阵和存储在其中的条目来初始化Trilinos矩阵。它使用一个阈值，只复制模数大于阈值的元素（因此deal.II矩阵中的零可以被过滤掉）。
     *
     */
    void
    reinit(
      const std::vector<IndexSet> &              parallel_partitioning,
      const ::dealii::BlockSparseMatrix<double> &dealii_block_sparse_matrix,
      const MPI_Comm &                           communicator = MPI_COMM_WORLD,
      const double                               drop_tolerance = 1e-13);

    /**
     * 这个函数使用deal.II稀疏矩阵和存储在其中的条目来初始化Trilinos矩阵。它使用一个阈值，只复制模数大于阈值的元素（所以deal.II矩阵中的零可以被过滤掉）。由于没有给出Epetra_Map，所有的元素将被本地存储。
     *
     */
    void
    reinit(const ::dealii::BlockSparseMatrix<double> &deal_ii_sparse_matrix,
           const double                               drop_tolerance = 1e-13);

    /**
     * 返回矩阵的状态，即在需要数据交换的操作之后是否需要调用compress()。只有在<tt>debug</tt>模式下使用时才返回非真值，因为跟踪所有导致需要压缩()的操作是相当昂贵的。
     *
     */
    bool
    is_compressed() const;

    /**
     * 这个函数收集了子对象的大小，并将其存储在内部数组中，以便能够将矩阵的全局索引转为子对象的索引。你必须*在你改变了子对象的大小之后，每次都调用这个函数。注意，这是一个集体操作，即需要在所有MPI进程中调用。这个命令在内部调用<tt>compress()</tt>方法，所以在使用<tt>collect_sizes()</tt>的情况下，你不需要调用这个函数。
     *
     */
    void
    collect_sizes();

    /**
     * 返回这个矩阵的非零元素的总数（所有MPI进程的总和）。
     *
     */
    size_type
    n_nonzero_elements() const;

    /**
     * 返回与该矩阵一起使用的MPI通信器对象。
     *
     */
    MPI_Comm
    get_mpi_communicator() const;

    /**
     * 返回该矩阵各个块的域空间的划分，即该矩阵要与之相乘的块向量的划分。
     *
     */
    std::vector<IndexSet>
    locally_owned_domain_indices() const;

    /**
     * 返回该矩阵各个块的范围空间的划分，即由矩阵-向量乘积产生的块向量的划分。
     *
     */
    std::vector<IndexSet>
    locally_owned_range_indices() const;

    /**
     * 矩阵-向量乘法：让 $dst = M*src$ 与 $M$
     * 是这个矩阵。矢量类型可以是块状矢量或非块状矢量（只有在矩阵只有一行或一列的情况下才可以），并且需要定义
     * TrilinosWrappers::SparseMatrix::vmult. 。
     *
     */
    template <typename VectorType1, typename VectorType2>
    void
    vmult(VectorType1 &dst, const VectorType2 &src) const;

    /**
     * 矩阵-向量乘法：让 $dst = M^T*src$ 与 $M$
     * 是这个矩阵。这个函数与vmult()的作用相同，但需要转置的矩阵。
     *
     */
    template <typename VectorType1, typename VectorType2>
    void
    Tvmult(VectorType1 &dst, const VectorType2 &src) const;

    /**
     * 计算方程<i>Mx=b</i>的残差，其中残差被定义为<i>r=b-Mx</i>。将残差写入
     * @p dst.  返回残差向量的<i>l<sub>2</sub></i>准则。
     * 源<i>x</i>和目的<i>dst</i>不能是同一个向量。
     * 注意，这两个向量必须是使用与矩阵相同的Map生成的分布式向量。
     * 这个函数只适用于矩阵只有一个块行的情况。
     *
     */
    TrilinosScalar
    residual(MPI::BlockVector &      dst,
             const MPI::BlockVector &x,
             const MPI::BlockVector &b) const;

    /**
     * 计算方程<i>Mx=b</i>的残差，其中残差被定义为<i>r=b-Mx</i>。将残差写进
     * @p dst.  返回残差向量的<i>l<sub>2</sub></i>准则。
     * 这个函数只适用于矩阵只有一个块行的情况。
     *
     */
    TrilinosScalar
    residual(MPI::BlockVector &      dst,
             const MPI::Vector &     x,
             const MPI::BlockVector &b) const;

    /**
     * 计算方程<i>Mx=b</i>的残差，其中残差被定义为<i>r=b-Mx</i>。将残差写入
     * @p dst.  返回残差向量的<i>l<sub>2</sub></i>准则。
     * 这个函数只适用于矩阵只有一个块列的情况。
     *
     */
    TrilinosScalar
    residual(MPI::Vector &           dst,
             const MPI::BlockVector &x,
             const MPI::Vector &     b) const;

    /**
     * 计算方程<i>Mx=b</i>的残差，其中残差被定义为<i>r=b-Mx</i>。将残差写入
     * @p dst.  返回残差向量的<i>l<sub>2</sub></i>准则。
     * 这个函数只适用于矩阵只有一个块的情况。
     *
     */
    TrilinosScalar
    residual(MPI::Vector &      dst,
             const MPI::Vector &x,
             const MPI::Vector &b) const;

    /**
     * 使基类中的clear()函数可见，尽管它是受保护的。
     *
     */
    using BlockMatrixBase<SparseMatrix>::clear;

    /**
     * @addtogroup  Exceptions  @{
     *
     */

    /**
     * 异常情况
     *
     */
    DeclException4(ExcIncompatibleRowNumbers,
                   int,
                   int,
                   int,
                   int,
                   << "The blocks [" << arg1 << ',' << arg2 << "] and [" << arg3
                   << ',' << arg4 << "] have differing row numbers.");

    /**
     * 异常情况
     *
     */
    DeclException4(ExcIncompatibleColNumbers,
                   int,
                   int,
                   int,
                   int,
                   << "The blocks [" << arg1 << ',' << arg2 << "] and [" << arg3
                   << ',' << arg4 << "] have differing column numbers.");
    ///@}

  private:
    /**
     * 内部版本的(T)vmult有两个块向量
     *
     */
    template <typename VectorType1, typename VectorType2>
    void
    vmult(VectorType1 &      dst,
          const VectorType2 &src,
          const bool         transpose,
          const std::integral_constant<bool, true>,
          const std::integral_constant<bool, true>) const;

    /**
     * (T)vmult的内部版本，其中源向量是一个块向量，但目的向量是一个非块向量
     *
     */
    template <typename VectorType1, typename VectorType2>
    void
    vmult(VectorType1 &      dst,
          const VectorType2 &src,
          const bool         transpose,
          const std::integral_constant<bool, false>,
          const std::integral_constant<bool, true>) const;

    /**
     * (T)vmult的内部版本，其中源向量是一个非块向量，但目的向量是一个块向量
     *
     */
    template <typename VectorType1, typename VectorType2>
    void
    vmult(VectorType1 &      dst,
          const VectorType2 &src,
          const bool         transpose,
          const std::integral_constant<bool, true>,
          const std::integral_constant<bool, false>) const;

    /**
     * (T)vmult的内部版本，其中源向量和目的向量均为非块向量（仅在矩阵仅由一个块组成时定义）。
     *
     */
    template <typename VectorType1, typename VectorType2>
    void
    vmult(VectorType1 &      dst,
          const VectorType2 &src,
          const bool         transpose,
          const std::integral_constant<bool, false>,
          const std::integral_constant<bool, false>) const;
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



  inline bool
  BlockSparseMatrix::is_compressed() const
  {
    bool compressed = true;
    for (size_type row = 0; row < n_block_rows(); ++row)
      for (size_type col = 0; col < n_block_cols(); ++col)
        if (block(row, col).is_compressed() == false)
          {
            compressed = false;
            break;
          }

    return compressed;
  }



  template <typename VectorType1, typename VectorType2>
  inline void
  BlockSparseMatrix::vmult(VectorType1 &dst, const VectorType2 &src) const
  {
    vmult(dst,
          src,
          false,
          std::integral_constant<bool, IsBlockVector<VectorType1>::value>(),
          std::integral_constant<bool, IsBlockVector<VectorType2>::value>());
  }



  template <typename VectorType1, typename VectorType2>
  inline void
  BlockSparseMatrix::Tvmult(VectorType1 &dst, const VectorType2 &src) const
  {
    vmult(dst,
          src,
          true,
          std::integral_constant<bool, IsBlockVector<VectorType1>::value>(),
          std::integral_constant<bool, IsBlockVector<VectorType2>::value>());
  }



  template <typename VectorType1, typename VectorType2>
  inline void
  BlockSparseMatrix::vmult(VectorType1 &      dst,
                           const VectorType2 &src,
                           const bool         transpose,
                           std::integral_constant<bool, true>,
                           std::integral_constant<bool, true>) const
  {
    if (transpose == true)
      BaseClass::Tvmult_block_block(dst, src);
    else
      BaseClass::vmult_block_block(dst, src);
  }



  template <typename VectorType1, typename VectorType2>
  inline void
  BlockSparseMatrix::vmult(VectorType1 &      dst,
                           const VectorType2 &src,
                           const bool         transpose,
                           std::integral_constant<bool, false>,
                           std::integral_constant<bool, true>) const
  {
    if (transpose == true)
      BaseClass::Tvmult_nonblock_block(dst, src);
    else
      BaseClass::vmult_nonblock_block(dst, src);
  }



  template <typename VectorType1, typename VectorType2>
  inline void
  BlockSparseMatrix::vmult(VectorType1 &      dst,
                           const VectorType2 &src,
                           const bool         transpose,
                           std::integral_constant<bool, true>,
                           std::integral_constant<bool, false>) const
  {
    if (transpose == true)
      BaseClass::Tvmult_block_nonblock(dst, src);
    else
      BaseClass::vmult_block_nonblock(dst, src);
  }



  template <typename VectorType1, typename VectorType2>
  inline void
  BlockSparseMatrix::vmult(VectorType1 &      dst,
                           const VectorType2 &src,
                           const bool         transpose,
                           std::integral_constant<bool, false>,
                           std::integral_constant<bool, false>) const
  {
    if (transpose == true)
      BaseClass::Tvmult_nonblock_nonblock(dst, src);
    else
      BaseClass::vmult_nonblock_nonblock(dst, src);
  }



  inline std::vector<IndexSet>
  BlockSparseMatrix::locally_owned_domain_indices() const
  {
    Assert(this->n_block_cols() != 0, ExcNotInitialized());
    Assert(this->n_block_rows() != 0, ExcNotInitialized());

    std::vector<IndexSet> domain_indices;
    for (size_type c = 0; c < this->n_block_cols(); ++c)
      domain_indices.push_back(
        this->sub_objects[0][c]->locally_owned_domain_indices());

    return domain_indices;
  }



  inline std::vector<IndexSet>
  BlockSparseMatrix::locally_owned_range_indices() const
  {
    Assert(this->n_block_cols() != 0, ExcNotInitialized());
    Assert(this->n_block_rows() != 0, ExcNotInitialized());

    std::vector<IndexSet> range_indices;
    for (size_type r = 0; r < this->n_block_rows(); ++r)
      range_indices.push_back(
        this->sub_objects[r][0]->locally_owned_range_indices());

    return range_indices;
  }



  namespace internal
  {
    namespace BlockLinearOperatorImplementation
    {
      /**
       * 这是一个BlockLinearOperators的扩展类，用于Trilinos块稀疏矩阵。
       * @note
       * 这个类目前做得很少，只是检查每个子块的正确Payload类型是否被正确选择。为了给BlockLinearOperators增加更多的功能，同时保持与Trilinos稀疏矩阵和预处理类的兼容性，未来可能需要对该类进行进一步的扩展。
       * @ingroup TrilinosWrappers
       *
       */
      template <typename PayloadBlockType>
      class TrilinosBlockPayload
      {
      public:
        /**
         * 每个子块所持有的有效载荷的类型
         *
         */
        using BlockType = PayloadBlockType;

        /**
         * 默认构造函数
         * 这只是检查每个块的有效载荷是否被正确选择（即属于TrilinosPayload类型）。除此以外，这个类不做任何特别的事情，也不需要特别的配置，我们只有一个通用的构造函数，可以在任何条件下调用。
         *
         */
        template <typename... Args>
        TrilinosBlockPayload(const Args &...)
        {
          static_assert(
            std::is_same<
              PayloadBlockType,
              internal::LinearOperatorImplementation::TrilinosPayload>::value,
            "TrilinosBlockPayload can only accept a payload of type TrilinosPayload.");
        }
      };

    } // namespace BlockLinearOperatorImplementation
  }    /* namespace internal */ 


}  /* namespace TrilinosWrappers */ 


DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_TRILINOS

#endif // dealii_trilinos_block_sparse_matrix_h


