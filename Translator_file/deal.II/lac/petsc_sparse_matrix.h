//include/deal.II-translator/lac/petsc_sparse_matrix_0.txt
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

#ifndef dealii_petsc_sparse_matrix_h
#  define dealii_petsc_sparse_matrix_h


#  include <deal.II/base/config.h>

#  ifdef DEAL_II_WITH_PETSC

#    include <deal.II/lac/exceptions.h>
#    include <deal.II/lac/petsc_matrix_base.h>
#    include <deal.II/lac/petsc_vector.h>

#    include <vector>

DEAL_II_NAMESPACE_OPEN
// forward declaration
#    ifndef DOXYGEN
template <typename MatrixType>
class BlockMatrixBase;
#    endif

namespace PETScWrappers
{
  /**
   * 实现一个基于PETSc的连续稀疏矩阵类。所有的功能实际上都在基类中，除了生成连续稀疏矩阵的调用。这是可能的，因为PETSc只在一个抽象的矩阵类型上工作，并在内部根据实际的矩阵类型分配给做实际工作的函数（很像使用虚拟函数）。只有创建特定类型矩阵的函数不同，并在这个特定的类中实现。
   * @ingroup PETScWrappers
   * @ingroup Matrix1
   *
   */
  class SparseMatrix : public MatrixBase
  {
  public:
    /**
     * 一个描述这个类在运行时行为方面的一些特征的结构。其他一些以一个或其他矩阵类作为模板参数的类（如块状矩阵类），可以根据这个类中的变量来调整其行为。
     *
     */
    struct Traits
    {
      /**
       * 对该矩阵的单个元素进行零的添加是安全的。
       *
       */
      static const bool zero_addition_can_be_elided = true;
    };

    /**
     * 默认构造函数。创建一个空矩阵。
     *
     */
    SparseMatrix();

    /**
     * 创建一个尺寸为 @p m 乘以 @p n, 的稀疏矩阵，每行有 @p
     * n_nonzero_per_row
     * 个非零元素的初始猜测。PETSc能够应对后来为某一行分配的元素超过这个数量的情况，但这涉及到复制数据，因此很昂贵。
     * @p is_symmetric
     * 标志决定了我们是否应该告诉PETSc矩阵将是对称的（如调用<tt>MatSetOption(mat,
     * MAT_SYMMETRIC)</tt>所示。请注意，PETSc文档指出，对于这个标志被设置为
     * @p true,
     * 的矩阵，我们不能形成ILU分解，只能形成ICC。这个标志的默认值是
     * @p false.  。
     *
     */
    SparseMatrix(const size_type m,
                 const size_type n,
                 const size_type n_nonzero_per_row,
                 const bool      is_symmetric = false);

    /**
     * 初始化一个具有 @p m 行和 @p n 列的矩形矩阵。
     * 每行的最大非零条目数分别由 @p row_lengths 数组给出。
     * 正如其他构造函数一样。PETSc能够应对后来为某一行分配的元素超过这个数量的情况，但这涉及到复制数据，因此很昂贵。
     * @p is_symmetric
     * 标志决定了我们是否应该告诉PETSc矩阵将是对称的（如调用<tt>MatSetOption(mat,
     * MAT_SYMMETRIC)</tt>所示。请注意，PETSc文档中指出，对于这个标志被设置为
     * @p true,
     * 的矩阵，不能形成ILU分解，只能形成ICC。这个标志的默认值是
     * @p false.  。
     *
     */
    SparseMatrix(const size_type               m,
                 const size_type               n,
                 const std::vector<size_type> &row_lengths,
                 const bool                    is_symmetric = false);

    /**
     * 使用给定的稀疏模式初始化一个稀疏矩阵。
     * 请注意，如果你不向PETSc提供对行长度的良好估计，它可能会非常慢。使用本函数是一个非常有效的方法，因为它通过使用给定的稀疏模式参数，为矩阵的每一行使用准确的非零条目数。如果
     * @p preset_nonzero_locations 标志是 @p true,
     * ，这个函数不仅预先设置了正确的行大小，而且还预先分配了矩阵中正确的非零条目。
     * PETsc允许以后向矩阵添加额外的非零条目，只需向这些元素写入即可。然而，这将导致额外的内存分配，效率非常低，将大大降低你的程序速度。因此，从一开始就做好内存分配是非常有效的。
     *
     */
    template <typename SparsityPatternType>
    explicit SparseMatrix(const SparsityPatternType &sparsity_pattern,
                          const bool preset_nonzero_locations = true);

    /**
     * 这个操作符将一个标量分配给一个矩阵。由于这通常没有什么意义（我们应该把所有的矩阵条目都设置为这个值吗？仅仅是稀疏模式的非零条目？），这个操作只允许在实际要分配的值为零时进行。这个操作符的存在只是为了允许明显的符号<tt>matrix=0</tt>，它将矩阵的所有元素设置为零，但保留之前使用的稀疏模式。
     *
     */
    SparseMatrix &
    operator=(const double d);

    /**
     * 复制构造函数被删除。
     *
     */
    SparseMatrix(const SparseMatrix &) = delete;

    /**
     * 拷贝赋值运算符被删除。
     *
     */
    SparseMatrix &
    operator=(const SparseMatrix &) = delete;

    /**
     * 扔掉现在的矩阵，并生成一个具有相同属性的矩阵，就像它是由这个类的构造函数创建的，参数列表与现在的函数相同。
     *
     */
    void
    reinit(const size_type m,
           const size_type n,
           const size_type n_nonzero_per_row,
           const bool      is_symmetric = false);

    /**
     * 扔掉当前的矩阵，并生成一个具有相同属性的矩阵，就像它是由这个类的构造函数创建的一样，参数列表与当前函数相同。
     *
     */
    void
    reinit(const size_type               m,
           const size_type               n,
           const std::vector<size_type> &row_lengths,
           const bool                    is_symmetric = false);

    /**
     * 使用给定的稀疏模式初始化一个稀疏矩阵。
     * 请注意，如果你不向PETSc提供对行长度的良好估计，它可能会非常慢。使用本函数是一个非常有效的方法，因为它通过使用给定的稀疏模式参数，对矩阵的每一行使用准确的非零条目数。如果
     * @p preset_nonzero_locations 标志是 @p true,
     * ，这个函数不仅预先设置了正确的行大小，而且还预先分配了矩阵中正确的非零条目。
     * PETsc允许以后向矩阵添加额外的非零条目，只需向这些元素写入即可。然而，这将导致额外的内存分配，效率非常低，将大大降低你的程序速度。因此，从一开始就做好内存分配是非常有效的。
     * 尽管这似乎是一个明显的胜利，但将 @p
     * preset_nonzero_locations 标志设置为 @p true
     * 似乎并没有加速程序。相反，它似乎能够在一定程度上减慢整个程序。这是令人惊讶的，因为我们可以使用高效的函数调用到PETSc中，允许一次创建多个条目；尽管如此，考虑到这是低效的，各自的标志的默认值等于
     * @p false.  。
     *
     */
    template <typename SparsityPatternType>
    void
    reinit(const SparsityPatternType &sparsity_pattern,
           const bool                 preset_nonzero_locations = true);

    /**
     * 返回一个对该矩阵使用的MPI通信器对象的引用。由于这是一个连续的矩阵，它返回MPI_COMM_SELF通信器。
     *
     */
    virtual const MPI_Comm &
    get_mpi_communicator() const override;

    /**
     * 返回这个矩阵的行数。
     *
     */
    size_t
    m() const;

    /**
     * 返回此矩阵的列数。
     *
     */
    size_t
    n() const;

    /**
     * 执行矩阵与矩阵的乘法  $C = AB$  ，或者， $C = A
     * \text{diag}(V) B$  给出一个兼容的矢量  $V$  。
     * 这个函数调用  MatrixBase::mmult()  来做实际工作。
     *
     */
    void
    mmult(SparseMatrix &      C,
          const SparseMatrix &B,
          const MPI::Vector & V = MPI::Vector()) const;

    /**
     * 与<tt>this</tt>的转置进行矩阵-矩阵乘法，即 $C = A^T B$
     * ，或者， $C = A^T \text{diag}(V) B$ 给定一个兼容的向量 $V$
     * 。        这个函数调用 MatrixBase::Tmmult()
     * 来做实际工作。
     *
     */
    void
    Tmmult(SparseMatrix &      C,
           const SparseMatrix &B,
           const MPI::Vector & V = MPI::Vector()) const;

  private:
    /**
     * 为各自的reinit()函数和匹配的构造函数做实际工作，即创建一个矩阵。摆脱之前的矩阵是留给调用者的。
     *
     */
    void
    do_reinit(const size_type m,
              const size_type n,
              const size_type n_nonzero_per_row,
              const bool      is_symmetric = false);

    /**
     * 与之前的函数相同。
     *
     */
    void
    do_reinit(const size_type               m,
              const size_type               n,
              const std::vector<size_type> &row_lengths,
              const bool                    is_symmetric = false);

    /**
     * 同前面的函数。
     *
     */
    template <typename SparsityPatternType>
    void
    do_reinit(const SparsityPatternType &sparsity_pattern,
              const bool                 preset_nonzero_locations);

    // To allow calling protected prepare_add() and prepare_set().
    friend class BlockMatrixBase<SparseMatrix>;
  };

  namespace MPI
  {
    /**
     * 实现一个基于PETSc的并行稀疏矩阵类，矩阵的行分布在一个MPI网络上。所有的功能实际上都在基类中，除了用于生成并行稀疏矩阵的调用。这是可能的，因为PETSc只在抽象矩阵类型上工作，并在内部分配给根据实际矩阵类型进行实际工作的函数（很像使用虚拟函数）。只有创建特定类型矩阵的函数不同，并在这个特定的类中实现。
     * 在平行向量类的文档中，有许多关于通信模型以及访问单个元素的评论。这些评论在这里也适用。
     * <h3>Partitioning of matrices</h3>
     * PETSc对并行矩阵进行分区，以便每个MPI进程 "拥有
     * "一定数量的行（即只有这个进程在这些行中存储各自的条目）。每个进程拥有的行数必须通过参数
     * @p  local_rows传递给构造函数和
     * reinit()函数。在所有MPI进程中作为 @p local_rows
     * 传递的单个值当然要与矩阵的全局行数相加。
     * 除此之外，PETSc还对它所拥有的矩阵的矩形块（即矩阵中
     * @p local_rows
     * 乘以n()的元素）进行分区，这样就可以有效地进行矩阵向量乘法。因此，这种列的划分必须与矩阵被乘以的向量的划分相匹配，就像行的划分必须与目的向量的划分相匹配。这个分区是通过
     * @p local_columns
     * 变量传递给构造函数和reinit()函数的，这个变量也必须与矩阵中的全局列数相加。local_columns
     * @p
     * 这个名字可能命名不当，因为它并没有反映出只有这些列被存储在本地，但它反映出这些列是传入向量的元素被存储在本地的事实。
     * 让事情变得更加复杂的是，PETSc需要对每一行要存储的元素数量有一个非常好的估计，以提高工作效率。
     * 否则，它就会把大部分时间花在分配小块内存上，这个过程如果经常发生，会使程序变得缓慢。如果对每行的条目数进行良好的估计还不够，它甚至需要将其拆分如下：对于它拥有的每一行，它需要估计这一行中属于为这个过程而设置的列的元素数（见上文），以及其他列中的元素数。
     * 因为一般来说，这个信息不是很容易得到，所以这个类的大多数初始化函数都假设你作为参数给
     * @p n_nonzero_per_row 或由 @p
     * row_lengths的所有元素数量都落入这个过程 "拥有
     * "的列中，而没有落入其他列中。这对大多数行来说是一个合理的猜测，因为在一个好的域划分中，节点只与同一子域内的节点进行交互。然而，对于子域界面上的节点来说，这并不成立，对于这些节点对应的行，PETSc将不得不分配额外的内存，这是一个昂贵的过程。
     * 避免这种情况的唯一方法是告诉PETSc矩阵的实际条目将在哪里。为此，这个类有一些构造函数和
     * reinit() 函数，它们接收一个包含所有这些信息的
     * DynamicSparsityPattern
     * 对象。虽然在一般情况下，如果构造函数和reinit()函数知道本地行和列的数量就足够了，但获得稀疏模式的函数也需要知道本地行的数量（
     * @p local_rows_per_process)  和列的数量（  @p
     * local_columns_per_process)
     * 对于所有其他过程，以便计算矩阵的哪些部分是哪些。因此，仅仅计算属于某个进程的自由度数量是不够的，你必须在所有进程中都有所有进程的数字。
     * @ingroup PETScWrappers
     * @ingroup Matrix1
     *
     */
    class SparseMatrix : public MatrixBase
    {
    public:
      /**
       * 声明容器大小的类型。
       *
       */
      using size_type = types::global_dof_index;

      /**
       * 一个描述该类在运行时行为方面的一些特征的结构。其他一些以一个或其他矩阵类作为其模板参数的类（如块状矩阵类），可以根据这个类中的变量来调整其行为。
       *
       */
      struct Traits
      {
        /**
         * 对这个矩阵的单个元素进行零的添加并不安全。原因是对矩阵的加法可能会触发集体操作，在多个进程中同步缓冲区。如果一个进程的加法被省略了，这可能会导致其他进程陷入无限的等待循环。
         *
         */
        static const bool zero_addition_can_be_elided = false;
      };

      /**
       * 默认构造函数。创建一个空矩阵。
       *
       */
      SparseMatrix();

      /**
       * 解构器释放PETSc对象。
       *
       */
      ~SparseMatrix() override;

      /**
       * 使用给定的稀疏模式进行初始化，通信发生在所提供的
       * @p communicator.   @p local_rows_per_process  和  @p
       * local_columns_per_process 参数的含义，请参见类文档。
       * 请注意，如果你不向PETSc提供对行的长度的良好估计，它的速度会非常慢。使用本函数是一个非常有效的方法，因为它通过使用给定的稀疏模式参数，为矩阵的每一行使用准确的非零条目数。如果
       * @p preset_nonzero_locations 标志是 @p true,
       * ，这个函数不仅预先设置了正确的行大小，而且还预先分配了矩阵中正确的非零条目。
       * PETsc允许以后向矩阵添加额外的非零条目，只需向这些元素写入即可。然而，这将导致额外的内存分配，效率非常低，将大大降低你的程序速度。因此，从一开始就做好内存分配是非常有效的。
       *
       */
      template <typename SparsityPatternType>
      SparseMatrix(const MPI_Comm &              communicator,
                   const SparsityPatternType &   sparsity_pattern,
                   const std::vector<size_type> &local_rows_per_process,
                   const std::vector<size_type> &local_columns_per_process,
                   const unsigned int            this_process,
                   const bool preset_nonzero_locations = true);

      /**
       * 这个操作符将一个标量分配给一个矩阵。由于这通常没有什么意义（我们应该把所有的矩阵条目都设置为这个值吗？
       * 仅仅是稀疏模式的非零条目？），这个操作只允许在实际要分配的值为零时进行。这个操作符的存在只是为了允许明显的符号<tt>matrix=0</tt>，它将矩阵的所有元素设置为零，但保留之前使用的稀疏模式。
       *
       */
      SparseMatrix &
      operator=(const value_type d);


      /**
       * 制作PETSc矩阵的副本  @p other.
       * 假设两个矩阵具有相同的SparsityPattern。
       *
       */
      void
      copy_from(const SparseMatrix &other);

      /**
       * 使用给定的稀疏模式进行初始化，通信发生在所提供的
       * @p communicator.
       * 上。注意，如果你不向PETSc提供对行长度的良好估计，它可能会非常慢。使用本函数是一个非常有效的方法，因为它通过使用给定的稀疏模式参数，使用矩阵每一行的非零条目的确切数量。如果
       * @p preset_nonzero_locations 标志是 @p true,
       * ，这个函数不仅预先设置了正确的行大小，而且还预先分配了矩阵中正确的非零条目。
       * PETsc允许以后向矩阵添加额外的非零条目，只需向这些元素写入即可。然而，这将导致额外的内存分配，效率非常低，将大大降低你的程序速度。因此，从一开始就做好内存分配是非常有效的。
       *
       */
      template <typename SparsityPatternType>
      void
      reinit(const MPI_Comm &              communicator,
             const SparsityPatternType &   sparsity_pattern,
             const std::vector<size_type> &local_rows_per_process,
             const std::vector<size_type> &local_columns_per_process,
             const unsigned int            this_process,
             const bool                    preset_nonzero_locations = true);

      /**
       * 创建一个矩阵，其中IndexSets的size()决定了全局的行和列的数量，IndexSet的条目给出了调用处理器的行和列。注意，只支持升序的1:1
       * IndexSets。
       *
       */
      template <typename SparsityPatternType>
      void
      reinit(const IndexSet &           local_rows,
             const IndexSet &           local_columns,
             const SparsityPatternType &sparsity_pattern,
             const MPI_Comm &           communicator);

      /**
       * 初始化这个矩阵，使其具有与 @p other. 相同的结构
       * 这不会复制其他矩阵的值，但你可以使用copy_from()来实现。
       *
       */
      void
      reinit(const SparseMatrix &other);

      /**
       * 返回一个对与该矩阵一起使用的MPI通信器对象的引用。
       *
       */
      virtual const MPI_Comm &
      get_mpi_communicator() const override;

      /**
       * @addtogroup  Exceptions  
     * @{ 
       *
       */
      /**
       * 异常情况
       *
       */
      DeclException2(ExcLocalRowsTooLarge,
                     int,
                     int,
                     << "The number of local rows " << arg1
                     << " must be larger than the total number of rows "
                     << arg2);
      //@}

      /**
       * 返回向量  $v$  相对于该矩阵诱导的法线的平方，即
       * $\left(v^\ast,Mv\right)$
       * 。这很有用，例如在有限元背景下，一个函数的 $L_2$
       * 规范等于相对于代表有限元函数节点值的向量的质量矩阵的矩阵规范。
       * 很明显，对于这个操作，矩阵需要是二次的。
       * 这个函数的实现没有deal.II中使用的 @p MatrixBase
       * 类（即原始的，而不是PETSc封装类）的效率高，因为PETSc不支持这个操作，需要一个临时矢量。
       *
       */
      PetscScalar
      matrix_norm_square(const Vector &v) const;

      /**
       * 计算矩阵标量乘积  $\left(u^\ast,Mv\right)$  。
       * 这个函数的实现没有deal.II中使用的 @p MatrixBase
       * 类（即原始的函数，而不是PETSc包装类）的效率高，因为PETSc不支持这个操作，需要一个临时向量。
       *
       */
      PetscScalar
      matrix_scalar_product(const Vector &u, const Vector &v) const;

      /**
       * 返回这个矩阵的域空间的划分，即这个矩阵要与之相乘的向量的划分。
       *
       */
      IndexSet
      locally_owned_domain_indices() const;

      /**
       * 返回该矩阵的范围空间的划分，即由矩阵-向量乘积产生的向量的划分。
       *
       */
      IndexSet
      locally_owned_range_indices() const;

      /**
       * 执行矩阵-矩阵乘法  $C = AB$  ，或者， $C = A
       * \text{diag}(V) B$  给出一个兼容的向量  $V$  。
       * 这个函数调用  MatrixBase::mmult()  来做实际工作。
       *
       */
      void
      mmult(SparseMatrix &      C,
            const SparseMatrix &B,
            const MPI::Vector & V = MPI::Vector()) const;

      /**
       * 与<tt>this</tt>的转置进行矩阵-矩阵乘法，即 $C = A^T B$
       * ，或者， $C = A^T \text{diag}(V) B$ 给定一个兼容的向量
       * $V$ 。            这个函数调用 MatrixBase::Tmmult()
       * 来做实际工作。
       *
       */
      void
      Tmmult(SparseMatrix &      C,
             const SparseMatrix &B,
             const MPI::Vector & V = MPI::Vector()) const;

    private:
      /**
       * 将用于该并行矢量的通信器对象的副本。
       *
       */
      MPI_Comm communicator;

      /**
       * 与之前的函数相同。
       *
       */
      template <typename SparsityPatternType>
      void
      do_reinit(const SparsityPatternType &   sparsity_pattern,
                const std::vector<size_type> &local_rows_per_process,
                const std::vector<size_type> &local_columns_per_process,
                const unsigned int            this_process,
                const bool                    preset_nonzero_locations);

      /**
       * 与以前的函数相同。
       *
       */
      template <typename SparsityPatternType>
      void
      do_reinit(const IndexSet &           local_rows,
                const IndexSet &           local_columns,
                const SparsityPatternType &sparsity_pattern);

      // To allow calling protected prepare_add() and prepare_set().
      friend class BlockMatrixBase<SparseMatrix>;
    };



    // -------- template and inline functions ----------

    inline const MPI_Comm &
    SparseMatrix::get_mpi_communicator() const
    {
      return communicator;
    }
  } // namespace MPI
} // namespace PETScWrappers

DEAL_II_NAMESPACE_CLOSE

#  endif // DEAL_II_WITH_PETSC

#endif
 /*--------------------------- petsc_sparse_matrix.h -------------------------*/ 


