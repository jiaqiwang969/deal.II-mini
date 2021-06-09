//include/deal.II-translator/lac/trilinos_sparse_matrix_0.txt
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

#ifndef dealii_trilinos_sparse_matrix_h
#  define dealii_trilinos_sparse_matrix_h


#  include <deal.II/base/config.h>

#  ifdef DEAL_II_WITH_TRILINOS

#    include <deal.II/base/index_set.h>
#    include <deal.II/base/subscriptor.h>

#    include <deal.II/lac/exceptions.h>
#    include <deal.II/lac/full_matrix.h>
#    include <deal.II/lac/trilinos_epetra_vector.h>
#    include <deal.II/lac/trilinos_index_access.h>
#    include <deal.II/lac/trilinos_tpetra_vector.h>
#    include <deal.II/lac/trilinos_vector.h>
#    include <deal.II/lac/vector_memory.h>
#    include <deal.II/lac/vector_operation.h>

#    include <Epetra_Comm.h>
#    include <Epetra_CrsGraph.h>
#    include <Epetra_Export.h>
#    include <Epetra_FECrsMatrix.h>
#    include <Epetra_Map.h>
#    include <Epetra_MultiVector.h>
#    include <Epetra_Operator.h>

#    include <cmath>
#    include <memory>
#    include <type_traits>
#    include <vector>
#    ifdef DEAL_II_WITH_MPI
#      include <Epetra_MpiComm.h>
#      include <mpi.h>
#    else
#      include <Epetra_SerialComm.h>
#    endif

DEAL_II_NAMESPACE_OPEN

// forward declarations
#    ifndef DOXYGEN
template <typename MatrixType>
class BlockMatrixBase;

template <typename number>
class SparseMatrix;
class SparsityPattern;
class DynamicSparsityPattern;

namespace TrilinosWrappers
{
  class SparseMatrix;
  class SparsityPattern;

  namespace SparseMatrixIterators
  {
    template <bool Constness>
    class Iterator;
  }
} // namespace TrilinosWrappers
#    endif

namespace TrilinosWrappers
{
  /**
   * Trilinos矩阵的迭代器
   *
   */
  namespace SparseMatrixIterators
  {
    /**
     * 异常情况
     *
     */
    DeclException0(ExcBeyondEndOfMatrix);

    /**
     * 异常情况
     *
     */
    DeclException3(ExcAccessToNonlocalRow,
                   std::size_t,
                   std::size_t,
                   std::size_t,
                   << "You tried to access row " << arg1
                   << " of a distributed sparsity pattern, "
                   << " but only rows " << arg2 << " through " << arg3
                   << " are stored locally and can be accessed.");

    /**
     * 处理常数和非常数访问器对象的指数 对于普通的
     * dealii::SparseMatrix,
     * ，我们将使用一个访问器来处理稀疏模式。对于Trilinos矩阵，这似乎并不那么简单，因此，我们在这里写一个小的基类。
     *
     */
    class AccessorBase
    {
    public:
      /**
       * 声明容器大小的类型。
       *
       */
      using size_type = dealii::types::global_dof_index;

      /**
       * 构造函数。
       *
       */
      AccessorBase(SparseMatrix *  matrix,
                   const size_type row,
                   const size_type index);

      /**
       * 这个对象所代表的元素的行数。
       *
       */
      size_type
      row() const;

      /**
       * 这个对象所代表的元素在行中的索引。
       *
       */
      size_type
      index() const;

      /**
       * 这个对象所代表的元素的列号。
       *
       */
      size_type
      column() const;

    protected:
      /**
       * 指向矩阵对象的指针。这个对象应该被适当的派生类作为常量指针或非常量来处理。为了能够同时实现这两个对象，这里不是常数，所以要小心处理
       *
       */
      mutable SparseMatrix *matrix;
      /**
       * 当前的行数。
       *
       */
      size_type a_row;

      /**
       * 当前在行中的索引。
       *
       */
      size_type a_index;

      /**
       * 丢弃旧的行缓存（它们可能仍然被其他访问器使用），并为这个访问器目前所指向的行生成新的行缓存。
       *
       */
      void
      visit_present_row();

      /**
       * 缓存，我们在这里存储当前行的列索引。这是必要的，因为Trilinos对其矩阵元素的访问相当困难，当我们进入某一行时，一次性复制该行的所有列项比反复向Trilinos索取单个列项要有效得多。这也有一定的意义，因为无论如何，我们很可能会按顺序访问它们。
       * 为了使迭代器/存取器的复制具有可接受的性能，我们为这些条目保留了一个共享指针，以便在必要时有多个存取器可以访问这些数据。
       *
       */
      std::shared_ptr<std::vector<size_type>> colnum_cache;

      /**
       * 该行的值的缓存。
       *
       */
      std::shared_ptr<std::vector<TrilinosScalar>> value_cache;
    };

    /**
     * 稀疏矩阵访问器的通用模板。第一个模板参数表示底层数字类型，第二个表示矩阵的常数。
     * 一般的模板没有被实现，只有针对第二个模板参数的两个可能值的特殊化。因此，这里列出的接口只是作为提供的模板，因为doxygen并没有链接这些特殊化。
     *
     */
    template <bool Constess>
    class Accessor : public AccessorBase
    {
      /**
       * 这个矩阵条目的值。
       *
       */
      TrilinosScalar
      value() const;

      /**
       * 此矩阵条目的值。
       *
       */
      TrilinosScalar &
      value();
    };

    /**
     * const Accessor的特殊化。
     *
     */
    template <>
    class Accessor<true> : public AccessorBase
    {
    public:
      /**
       * 这里要使用的矩阵的类型（包括constness）的类型定义。
       *
       */
      using MatrixType = const SparseMatrix;

      /**
       * 构造器。因为我们只使用访问器进行读取访问，所以一个常数矩阵指针就足够了。
       *
       */
      Accessor(MatrixType *matrix, const size_type row, const size_type index);

      /**
       * 复制构造函数，从一个常量或非常量访问器获取到一个常量访问器。
       *
       */
      template <bool Other>
      Accessor(const Accessor<Other> &a);

      /**
       * 这个矩阵条目的值。
       *
       */
      TrilinosScalar
      value() const;

    private:
      // Make iterator class a friend.
      template <bool>
      friend class Iterator;
    };

    /**
     * 可变访问器的特殊化。
     *
     */
    template <>
    class Accessor<false> : public AccessorBase
    {
      class Reference
      {
      public:
        /**
         * 构造函数。
         *
         */
        Reference(const Accessor<false> &accessor);

        /**
         * 对矩阵数据类型的转换操作。
         *
         */
        operator TrilinosScalar() const;

        /**
         * 将我们目前指向的矩阵的元素设置为 @p n. 。
         *
         */
        const Reference &
        operator=(const TrilinosScalar n) const;

        /**
         * 将 @p n 添加到我们目前指向的矩阵元素中。
         *
         */
        const Reference &
        operator+=(const TrilinosScalar n) const;

        /**
         * 从我们现在指向的矩阵元素中减去 @p n 。
         *
         */
        const Reference &
        operator-=(const TrilinosScalar n) const;

        /**
         * 将我们现在指向的矩阵元素乘以 @p n. 。
         *
         */
        const Reference &
        operator*=(const TrilinosScalar n) const;

        /**
         * 用我们现在指向的矩阵的元素除以 @p n. 。
         *
         */
        const Reference &
        operator/=(const TrilinosScalar n) const;

      private:
        /**
         * 指向访问器的指针，表示我们目前指向哪个元素。
         *
         */
        Accessor &accessor;
      };

    public:
      /**
       * 这里要使用的矩阵的类型（包括常数）的类型定义。
       *
       */
      using MatrixType = SparseMatrix;

      /**
       * 构造函数。因为我们只使用访问器进行读取访问，所以一个常数矩阵指针就足够了。
       *
       */
      Accessor(MatrixType *matrix, const size_type row, const size_type index);

      /**
       * 这个矩阵条目的值。
       *
       */
      Reference
      value() const;

    private:
      // Make iterator class a friend.
      template <bool>
      friend class Iterator;

      // Make Reference object a friend.
      friend class Reference;
    };

    /**
     * 这个类作为一个迭代器，在特里诺斯矩阵的元素上行走。这个类的实现与PETSc矩阵的实现类似。
     * 请注意，Trilinos以升序的方式存储每一行的元素。这与deal.II稀疏矩阵风格相反，在这种风格中，对角线元素（如果存在的话）被存储在所有其他数值之前，而在PETSc稀疏矩阵中，人们无法保证元素的一定顺序。
     * @ingroup TrilinosWrappers
     *
     */
    template <bool Constness>
    class Iterator
    {
    public:
      /**
       * 声明容器大小的类型。
       *
       */
      using size_type = dealii::types::global_dof_index;

      /**
       * 为我们要操作的矩阵类型（包括constness）提供类型定义。
       *
       */
      using MatrixType = typename Accessor<Constness>::MatrixType;

      /**
       * 构造函数。在矩阵 @p matrix
       * 中创建一个迭代器，用于给定行和其中的索引。
       *
       */
      Iterator(MatrixType *matrix, const size_type row, const size_type index);

      /**
       * 复制构造函数，可选择改变常态。
       *
       */
      template <bool Other>
      Iterator(const Iterator<Other> &other);

      /**
       * 前缀增量。
       *
       */
      Iterator<Constness> &
      operator++();

      /**
       * 后缀增量。
       *
       */
      Iterator<Constness>
      operator++(int);

      /**
       * 撤消运算符。
       *
       */
      const Accessor<Constness> &operator*() const;

      /**
       * 解除引用操作符。
       *
       */
      const Accessor<Constness> *operator->() const;

      /**
       * 比较。真，如果两个迭代器都指向同一个矩阵位置。
       *
       */
      bool
      operator==(const Iterator<Constness> &) const;

      /**
       * <tt>==</tt>的倒数。
       *
       */
      bool
      operator!=(const Iterator<Constness> &) const;

      /**
       * 比较运算符。如果第一行数字较小，或者行数字相等且第一个索引较小，则结果为真。
       *
       */
      bool
      operator<(const Iterator<Constness> &) const;

      /**
       * 比较运算符。与前一个运算符相反。
       *
       */
      bool
      operator>(const Iterator<Constness> &) const;

      /**
       * 异常情况
       *
       */
      DeclException2(ExcInvalidIndexWithinRow,
                     size_type,
                     size_type,
                     << "Attempt to access element " << arg2 << " of row "
                     << arg1 << " which doesn't have that many elements.");

    private:
      /**
       * 存储一个访问器类的对象。
       *
       */
      Accessor<Constness> accessor;

      template <bool Other>
      friend class Iterator;
    };

  } // namespace SparseMatrixIterators


  /**
   * 这个类实现了一个包装器，用于使用Trilinos分布式稀疏矩阵类Epetra_FECrsMatrix。这正是我们一直在处理的那种矩阵
   *
   * - 我们很可能从一些装配过程中得到它，其中不属于本地的条目也可能需要被写入，因此需要转发给所有者过程。 这个类被设计用来在分布式内存架构中使用，底部有一个MPI编译器，但也同样适用于串行进程。  该类工作的唯一要求是，Trilinos已经安装了与生成deal.II相同的编译器。    这个类的接口是以deal.II中现有的SparseMatrix类为模型的。它有几乎相同的成员函数，而且通常是可以交换的。然而，由于Trilinos只支持单一的标量类型（double），所以它没有模板化，只对double起作用。    请注意，Trilinos只保证在矩阵装配后调用了函数 @p GlobalAssemble 的情况下，操作才会达到你的期望。  因此，你需要在实际使用矩阵之前调用 SparseMatrix::compress() 。这也会调用 @p FillComplete ，通过丢弃未使用的元素来压缩稀疏矩阵的存储格式。  不过，Trilinos允许在调用这些函数后继续组装矩阵。    <h3>Thread safety of Trilinos matrices</h3> 当从共享内存中的几个线程向Trilinos矩阵写入时，必须记住几件事，因为这个类中没有内置锁来防止数据竞赛。在同一时间内同时访问同一矩阵行会导致数据竞赛，用户必须明确避免。然而，在以下三个条件下，可以从几个线程同时访问矩阵的<b>different</b>行。    <ul>   <li>  矩阵只使用一个MPI进程。    <li>  矩阵已经用reinit()方法进行了初始化，并采用了DynamicSparsityPattern（包括本地相关行的集合，即汇编程序可能写入的行）。    <li>  矩阵已经从一个 TrilinosWrappers::SparsityPattern 对象中初始化，该对象又被reinit函数初始化，指定了三个索引集，一个用于行，一个用于列，以及更大的 @p 可写行集，操作是一个加法。在未来的某个时间点，Trilinos的支持可能足够完整，以至于从一个已经被类似于 DoFTools::make_sparsity_pattern 的函数填充的 TrilinosWrappers::SparsityPattern 初始化，总是会产生一个允许多个进程写入同一矩阵行的矩阵。  然而，Trilinos至少在11.12版本之前并不正确支持这一功能。    </ul>  请注意， TrilinosWrappers::SparsityPattern 的所有其他reinit方法和构造函数都会导致矩阵需要按需分配非进程条目，这就破坏了线程安全。当然，对块状Trilinos稀疏模式和块状矩阵使用各自的reinit方法也会导致线程安全。
   * @ingroup TrilinosWrappers
   * @ingroup Matrix1
   *
   */
  class SparseMatrix : public Subscriptor
  {
  public:
    /**
     * 声明容器大小的类型。
     *
     */
    using size_type = dealii::types::global_dof_index;

    /**
     * 异常情况
     *
     */
    DeclException1(ExcAccessToNonlocalRow,
                   std::size_t,
                   << "You tried to access row " << arg1
                   << " of a non-contiguous locally owned row set."
                   << " The row " << arg1
                   << " is not stored locally and can't be accessed.");

    /**
     * 一个描述这个类在运行时行为方面的一些特征的结构。其他一些以一个或其他矩阵类作为其模板参数的类（如块状矩阵类），可以根据这个类中的变量来调整其行为。
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
     * 为迭代器类声明一个别名。
     *
     */
    using iterator = SparseMatrixIterators::Iterator<false>;

    /**
     * 为常量迭代器类声明一个别名。
     *
     */
    using const_iterator = SparseMatrixIterators::Iterator<true>;

    /**
     * 为所有其他的容器类声明一个别名，以示类比。
     *
     */
    using value_type = TrilinosScalar;

    /**
     * @name  构造函数和初始化。
     *
     */
    //@{
    /**
     * 默认构造函数。生成一个空的（零大小）矩阵。
     *
     */
    SparseMatrix();

    /**
     * 生成一个完全存储在本地的矩阵，有#m行和#n列。
     * 每行的列条目数被指定为最大条目数参数。
     *
     */
    SparseMatrix(const size_type    m,
                 const size_type    n,
                 const unsigned int n_max_entries_per_row);

    /**
     * 生成一个完全存储在本地的矩阵，有#m行和#n列。
     * 向量<tt>n_entries_per_row</tt>指定了每一行的条目数。
     *
     */
    SparseMatrix(const size_type                  m,
                 const size_type                  n,
                 const std::vector<unsigned int> &n_entries_per_row);

    /**
     * 从Trilinos稀疏模式对象生成一个矩阵。
     *
     */
    SparseMatrix(const SparsityPattern &InputSparsityPattern);

    /**
     * 移动构造函数。通过窃取内部数据创建一个新的稀疏矩阵。
     *
     */
    SparseMatrix(SparseMatrix &&other) noexcept;

    /**
     * 复制构造函数被删除。
     *
     */
    SparseMatrix(const SparseMatrix &) = delete;

    /**
     * operator=被删除。
     *
     */
    SparseMatrix &
    operator=(const SparseMatrix &) = delete;

    /**
     * 解构器。做成了虚拟，这样就可以使用指向这个类的指针。
     *
     */
    virtual ~SparseMatrix() override = default;

    /**
     * 这个函数用deal.II稀疏模式初始化Trilinos矩阵，也就是说，它使Trilinos
     * Epetra矩阵根据稀疏模式知道非零项的位置。这个函数是为了在串行程序中使用，不需要指定矩阵如何在不同处理器之间分配。这个函数也可以在%并行中使用，但建议使用Epetra_Map手动指定矩阵的%并行分区。
     * 当以%并行方式运行时，目前需要每个处理器都持有sparsity_pattern结构，因为每个处理器都要设置其行数。
     * 这是一个集体操作，需要在所有处理器上调用，以避免出现死锁。
     *
     */
    template <typename SparsityPatternType>
    void
    reinit(const SparsityPatternType &sparsity_pattern);

    /**
     * 这个函数从一个（可能是分布式的）Trilinos稀疏模式重新初始化Trilinos稀疏矩阵。
     * 这是一个集体操作，需要在所有处理器上调用，以避免出现死锁。
     * 如果你想从几个线程写到矩阵并使用MPI，你需要使用这个reinit方法，并使用已经创建的明确说明可写行的稀疏度模式。在所有其他情况下，你不能把MPI和多线程写进矩阵中。
     *
     */
    void
    reinit(const SparsityPattern &sparsity_pattern);

    /**
     * 这个函数将 @p sparse_matrix
     * 的布局复制到调用矩阵中。值不会被复制，但你可以使用copy_from()来实现。
     * 这是一个集体操作，需要在所有处理器上调用，以避免出现死锁。
     *
     */
    void
    reinit(const SparseMatrix &sparse_matrix);

    /**
     * 这个函数使用deal.II稀疏矩阵和存储在其中的条目来初始化Trilinos矩阵。它使用一个阈值，只复制模数大于阈值的元素（所以deal.II矩阵中的零可以被过滤掉）。
     * 可选参数<tt>copy_values</tt>决定是只使用输入矩阵的稀疏结构，还是也要复制矩阵的条目。
     * 这是一个集体操作，需要在所有处理器上调用，以避免出现死锁。
     * @note
     * 如果在最后一个参数中给出了不同的稀疏模式（即与第一个参数中给出的稀疏矩阵中使用的模式不同），那么生成的特里诺斯矩阵将具有如此给出的稀疏模式。当然，这也意味着给定的矩阵中所有不属于这个单独的稀疏模式的条目实际上将被删除。
     *
     */
    template <typename number>
    void
    reinit(const ::dealii::SparseMatrix<number> &dealii_sparse_matrix,
           const double                          drop_tolerance    = 1e-13,
           const bool                            copy_values       = true,
           const ::dealii::SparsityPattern *     use_this_sparsity = nullptr);

    /**
     * 这个reinit函数将一个Trilinos
     * Epetra_CrsMatrix作为输入，并复制其稀疏性模式。如果有此要求，甚至连内容（值）也将被复制。
     *
     */
    void
    reinit(const Epetra_CrsMatrix &input_matrix, const bool copy_values = true);
    //@}

    /**
     * @name  使用IndexSet描述的构造器和初始化
     *
     */
    //@{
    /**
     * 使用一个IndexSet和一个MPI通信器的构造器来描述%并行分区。参数
     * @p n_max_entries_per_row
     * 设置将被分配的每一行中的非零条目的数量。注意，这个数字不需要精确，甚至允许实际的矩阵结构有比构造函数中指定的更多的非零条目。但是在这里提供良好的估计仍然是有利的，因为这将大大增加矩阵设置的性能。然而，对矩阵-向量乘积的性能没有影响，因为Trilinos在使用前（在compress()步骤中）重新组织了矩阵内存。
     *
     */
    SparseMatrix(const IndexSet &   parallel_partitioning,
                 const MPI_Comm &   communicator          = MPI_COMM_WORLD,
                 const unsigned int n_max_entries_per_row = 0);

    /**
     * 与之前相同，但现在单独设置每个矩阵行中的非零的数量。因为在这种情况下，我们确切地知道矩阵中的元素数量，所以我们已经可以分配合适的内存数量，这使得创建过程（包括通过各自的
     * SparseMatrix::reinit 调用插入非零元素）大大加快。
     *
     */
    SparseMatrix(const IndexSet &                 parallel_partitioning,
                 const MPI_Comm &                 communicator,
                 const std::vector<unsigned int> &n_entries_per_row);

    /**
     * 这个构造函数与上面的构造函数类似，但现在它为行和列采取了两个不同的
     * IndexSet
     * 分区。这个接口是用来生成矩形矩阵的，其中第一个索引集描述了与矩阵行相关的自由度的平行分区，第二个索引集描述了矩阵列的分区。第二个索引集指定了这个矩阵要与之相乘的向量的分区，而不是实际出现在矩阵中的元素的分布。
     * 参数 @p n_max_entries_per_row
     * 定义了将为每一行分配多少内存。这个数字不需要准确，因为结构会在compress()调用中被重组。
     *
     */
    SparseMatrix(const IndexSet &row_parallel_partitioning,
                 const IndexSet &col_parallel_partitioning,
                 const MPI_Comm &communicator          = MPI_COMM_WORLD,
                 const size_type n_max_entries_per_row = 0);

    /**
     * 这个构造函数与上面的构造函数类似，但它现在需要两个不同的Epetra映射来表示行和列。这个接口是用来生成矩形矩阵的，其中一个映射指定了与矩阵行相关的自由度的平行分布%，第二个映射指定了与矩阵中的列相关的自由度的平行分布%。第二个映射也为矩阵向量乘积中的内部排列提供信息（即这个矩阵要与之相乘的向量分布），但不用于列的分布&ndash；相反，在任何情况下，一行的所有列元素都存储在同一个处理器上。向量<tt>n_entries_per_row</tt>指定了新生成的矩阵中每一行的条目数。
     *
     */
    SparseMatrix(const IndexSet &                 row_parallel_partitioning,
                 const IndexSet &                 col_parallel_partitioning,
                 const MPI_Comm &                 communicator,
                 const std::vector<unsigned int> &n_entries_per_row);

    /**
     * 这个函数是根据指定的稀疏性_模式初始化Tridinos
     * Epetra矩阵，同时根据用户提供的索引集和%并行通信器，将矩阵的行数重新分配给不同进程。在遵循教程程序风格的程序中，这个函数（以及对矩形矩阵的相应调用）是初始化矩阵大小、其在MPI进程中的分布（如果在%parallel中运行）以及非零元素位置的自然方式。Trilinos在内部存储了稀疏模式，所以在这个调用之后就不再需要它了，与deal.II自己的对象相反。可选的参数
     * @p exchange_data
     * 可用于重新初始化未完全构建的稀疏度模式。这个功能只对动态稀疏模式类型的输入稀疏模式实现。如果没有设置该标志，每个处理器只是设置疏散度模式中属于其行的元素。
     * 这是一个集体操作，需要在所有处理器上调用，以避免出现死锁。
     *
     */
    template <typename SparsityPatternType>
    void
    reinit(const IndexSet &           parallel_partitioning,
           const SparsityPatternType &sparsity_pattern,
           const MPI_Comm &           communicator  = MPI_COMM_WORLD,
           const bool                 exchange_data = false);

    /**
     * 这个函数类似于上面的另一个初始化函数，但现在也根据两个用户提供的索引集重新分配矩阵的行和列。
     * 要用于矩形矩阵。可选的参数 @p exchange_data
     * 可用于用未完全构建的稀疏模式进行重新初始化。这个功能只对动态稀疏模式类型的输入稀疏模式实现。
     * 这是一个集体操作，需要在所有处理器上调用，以避免出现死锁。
     *
     */
    template <typename SparsityPatternType>
    typename std::enable_if<
      !std::is_same<SparsityPatternType,
                    dealii::SparseMatrix<double>>::value>::type
    reinit(const IndexSet &           row_parallel_partitioning,
           const IndexSet &           col_parallel_partitioning,
           const SparsityPatternType &sparsity_pattern,
           const MPI_Comm &           communicator  = MPI_COMM_WORLD,
           const bool                 exchange_data = false);

    /**
     * 这个函数使用deal.II稀疏矩阵和存储在其中的条目来初始化Trilinos矩阵。它使用一个阈值，只复制模数大于阈值的元素（所以deal.II矩阵中的零可以被过滤掉）。与其他带有deal.II稀疏矩阵参数的reinit函数不同的是，该函数采用用户指定的%并行分区，而不是内部生成的。
     * 可选参数<tt>copy_values</tt>决定是只使用输入矩阵的稀疏结构还是也要复制矩阵条目。
     * 这是一个集体操作，需要在所有处理器上调用，以避免出现死锁。
     *
     */
    template <typename number>
    void
    reinit(const IndexSet &                      parallel_partitioning,
           const ::dealii::SparseMatrix<number> &dealii_sparse_matrix,
           const MPI_Comm &                      communicator = MPI_COMM_WORLD,
           const double                          drop_tolerance    = 1e-13,
           const bool                            copy_values       = true,
           const ::dealii::SparsityPattern *     use_this_sparsity = nullptr);

    /**
     * 这个函数类似于上面的另一个初始化函数与deal.II稀疏矩阵输入，但现在需要矩阵的行和列的索引集。选用于矩形矩阵。
     * 可选参数<tt>copy_values</tt>决定是只使用输入矩阵的稀疏结构还是也复制矩阵条目。
     * 这是一个集体操作，需要在所有处理器上调用，以避免出现死锁。
     *
     */
    template <typename number>
    void
    reinit(const IndexSet &                      row_parallel_partitioning,
           const IndexSet &                      col_parallel_partitioning,
           const ::dealii::SparseMatrix<number> &dealii_sparse_matrix,
           const MPI_Comm &                      communicator = MPI_COMM_WORLD,
           const double                          drop_tolerance    = 1e-13,
           const bool                            copy_values       = true,
           const ::dealii::SparsityPattern *     use_this_sparsity = nullptr);
    //@}
    /**
     * @name  矩阵的信息
     *
     */
    //@{

    /**
     * 返回这个矩阵的行数。
     *
     */
    size_type
    m() const;

    /**
     * 返回此矩阵中的列数。
     *
     */
    size_type
    n() const;

    /**
     * 返回矩阵的本地维度，即存储在当前MPI进程中的行数。对于顺序矩阵，这个数字与m()相同，但对于%并行矩阵，这个数字可能更小。
     * 要想知道哪些元素被存储在本地，可以使用local_range()。
     *
     */
    unsigned int
    local_size() const;

    /**
     * 返回一对指数，表明该矩阵的哪些行是本地存储的。第一个数字是存储的第一行的索引，第二个数字是本地存储的最后一行之后的那一行的索引。如果这是一个连续的矩阵，那么结果将是一对(0,m())，否则将是一对(i,i+n)，其中<tt>n=local_size()</tt>。
     *
     */
    std::pair<size_type, size_type>
    local_range() const;

    /**
     * 返回 @p index 是否在本地范围内，另见local_range()。
     *
     */
    bool
    in_local_range(const size_type index) const;

    /**
     * 返回这个矩阵的非零元素的总数（所有MPI进程的总和）。
     *
     */
    size_type
    n_nonzero_elements() const;

    /**
     * 特定行中的条目数。
     *
     */
    unsigned int
    row_length(const size_type row) const;

    /**
     * 返回矩阵的状态，即在需要数据交换的操作之后是否需要调用compress()。当方法set()被调用时，也需要调用compress()（即使在串行工作时）。
     *
     */
    bool
    is_compressed() const;

    /**
     * 确定这个对象的内存消耗（以字节为单位）的估计值。注意，如果在基于MPI的程序中调用这个方法，则只返回当前处理器上保留的内存。
     *
     */
    size_type
    memory_consumption() const;

    /**
     * 返回与该矩阵一起使用的MPI通信器对象。
     *
     */
    MPI_Comm
    get_mpi_communicator() const;

    //@}
    /**
     * @name  修改条目
     *
     */
    //@{

    /**
     * 这个操作符将一个标量分配给一个矩阵。因为这通常没有什么意义（我们应该把所有的矩阵条目都设置为这个值吗？
     * 仅仅是稀疏模式的非零条目？），只有当要分配的实际值为零时，才允许这个操作。这个操作符的存在只是为了允许明显的符号<tt>matrix=0</tt>，它将矩阵的所有元素设置为零，但保留了之前使用的稀疏模式。
     *
     */
    SparseMatrix &
    operator=(const double d);

    /**
     * 释放所有内存并返回到与调用默认构造函数后相同的状态。
     * 这是一个集体操作，需要在所有处理器上调用，以避免出现死锁。
     *
     */
    void
    clear();

    /**
     * 这个命令做两件事。      <ul>   <li>  如果矩阵在初始化时没有稀疏模式，则使用set()命令手动添加元素。当这个过程完成后，对compress()的调用重组了内部数据结构（稀疏模式），这样就可以在矩阵-向量产品中快速访问数据。      <li>  如果矩阵结构已经被固定（通过初始化稀疏模式或在设置阶段调用compress()），该命令将进行数据的%并行交换。    当我们在一个以上的（MPI）进程上进行装配时，这是必要的，因为这时一些非本地的行数据会积累在属于当前的处理器元素的节点上，但实际上是由另一个处理器持有的。这个命令通常在所有元素都被遍历后调用。      </ul>  在这两种情况下，这个函数都会压缩数据结构，并允许产生的矩阵用于所有其他操作，如矩阵-向量乘积。这是一个集体操作，也就是说，在%并行使用时，需要在所有处理器上运行。        更多信息见 @ref GlossCompress "压缩分布式对象"。
     *
     */
    void
    compress(::dealii::VectorOperation::values operation);

    /**
     * 将元素（<i>i,j</i>）设置为 @p value.
     * 只要compress()没有被调用，这个函数就能够在矩阵中插入新的元素，所以稀疏模式将被扩展。当compress()第一次被调用时（或者在矩阵被初始化为稀疏模式的情况下），不能添加新的元素，在未被初始化的位置插入元素会引发异常。
     * 如果矩阵是在没有稀疏模式的情况下构建的，并且新的矩阵条目是按需添加的，请注意底层Epetra_FECrsMatrix数据结构施加的以下行为。
     * 如果同一个矩阵条目被多次插入，那么即使
     * VectorOperation::insert
     * 被指定为compress()的参数，矩阵条目也会在调用compress()时被添加（因为Epetra不会在最终调用compress()前跟踪同一个条目的值）。如果你不能确保矩阵条目只被设置一次，那么在插入元素之前，用稀疏模式初始化矩阵，固定矩阵结构。
     *
     */
    void
    set(const size_type i, const size_type j, const TrilinosScalar value);

    /**
     * 将FullMatrix<double>中给出的所有元素设置到<tt>indices</tt>给出的稀疏矩阵位置。换句话说，这个函数将<tt>full_matrix</tt>中的元素写入调用的矩阵中，对矩阵的行和列都使用<tt>indices</tt>指定的本地到全球的索引。这个函数假设一个二次稀疏矩阵和一个二次全矩阵，这是FE计算中的通常情况。
     * 只要没有调用compress()，这个函数就能够在矩阵中插入新的元素，所以稀疏模式将被扩展。在第一次调用compress()后，或者矩阵被初始化为稀疏模式后，扩展稀疏模式就不再可能了，在未初始化的位置插入元素会产生一个异常。
     * 可选参数<tt>elide_zero_values</tt>可以用来指定零值是否应该被插入，或者应该被过滤掉。默认值是<tt>false</tt>，也就是说，即使是零值也要插入/替换。
     * 对于矩阵的构建没有稀疏模式，新的矩阵条目是按需添加的情况，请注意底层Epetra_FECrsMatrix数据结构施加的以下行为。
     * 如果同一个矩阵条目被多次插入，那么即使
     * VectorOperation::insert
     * 被指定为compress()的参数，矩阵条目也会在调用compress()时被添加（因为Epetra不会在最终调用compress()前跟踪同一个条目的值）。如果你不能确保矩阵条目只被设置一次，那么在插入元素之前，用稀疏模式初始化矩阵，固定矩阵结构。
     *
     */
    void
    set(const std::vector<size_type> &    indices,
        const FullMatrix<TrilinosScalar> &full_matrix,
        const bool                        elide_zero_values = false);

    /**
     * 与之前的功能相同，但现在包括了使用矩形full_matrices的可能性，以及在行和列上分别使用不同的局部到全局的索引。
     *
     */
    void
    set(const std::vector<size_type> &    row_indices,
        const std::vector<size_type> &    col_indices,
        const FullMatrix<TrilinosScalar> &full_matrix,
        const bool                        elide_zero_values = false);

    /**
     * 将矩阵指定行中的几个元素与<tt>col_indices</tt>给出的列索引设置为相应的值。
     * 只要没有调用compress()，这个函数就能够在矩阵中插入新的元素，所以稀疏模式将被扩展。在第一次调用compress()后，或者矩阵被初始化为稀疏模式后，扩展稀疏模式就不再可能了，在未初始化的位置插入元素将引发异常。
     * 可选参数<tt>elide_zero_values</tt>可以用来指定零值是否应该被插入，或者应该被过滤掉。默认值是<tt>false</tt>，也就是说，即使是零值也要插入/替换。
     * 对于矩阵的构建没有稀疏模式，新的矩阵条目是按需添加的情况，请注意底层Epetra_FECrsMatrix数据结构施加的以下行为。
     * 如果同一个矩阵条目被多次插入，那么即使
     * VectorOperation::insert
     * 被指定为compress()的参数，矩阵条目也会在调用compress()时被添加（因为Epetra不会在最终调用compress()前跟踪同一条目的值）。如果你不能确保矩阵条目只被设置一次，那么在插入元素之前，用稀疏模式初始化矩阵，固定矩阵结构。
     *
     */
    void
    set(const size_type                    row,
        const std::vector<size_type> &     col_indices,
        const std::vector<TrilinosScalar> &values,
        const bool                         elide_zero_values = false);

    /**
     * 将几个元素设置为<tt>values</tt>所给的值，在稀疏矩阵的col_indices所给的列中的某一行。
     * 只要没有调用compress()，这个函数就能够在矩阵中插入新的元素，所以稀疏模式将被扩展。在第一次调用compress()后，或者矩阵被初始化为稀疏模式后，扩展稀疏模式就不再可能了，在未初始化的位置插入元素将引发异常。
     * 可选参数<tt>elide_zero_values</tt>可以用来指定零值是否应该被插入，或者应该被过滤掉。默认值是<tt>false</tt>，也就是说，即使是零值也要插入/替换。
     * 对于矩阵的构建没有稀疏模式，新的矩阵条目是按需添加的情况，请注意底层Epetra_FECrsMatrix数据结构施加的以下行为。
     * 如果同一个矩阵条目被多次插入，那么即使
     * VectorOperation::insert
     * 被指定为compress()的参数，矩阵条目也会在调用compress()时被添加（因为Epetra不会在最终调用compress()前跟踪同一个条目的值）。如果你不能确保矩阵条目只被设置一次，那么在插入元素之前，用稀疏模式初始化矩阵，固定矩阵结构。
     *
     */
    template <typename Number>
    void
    set(const size_type  row,
        const size_type  n_cols,
        const size_type *col_indices,
        const Number *   values,
        const bool       elide_zero_values = false);

    /**
     * 在元素（<i>i,j</i>）上添加 @p value 。        就像deal.II
     * SparseMatrix<Number>类中的相应调用一样（但与基于PETSc的矩阵的情况不同），如果稀疏模式中不存在一个条目，这个函数会抛出一个异常。
     * 此外，如果<tt>value</tt>不是一个有限的数字，也会抛出一个异常。
     *
     */
    void
    add(const size_type i, const size_type j, const TrilinosScalar value);

    /**
     * 将FullMatrix<double>中给出的所有元素添加到由<tt>indices</tt>给出的稀疏矩阵位置。换句话说，这个函数将<tt>full_matrix</tt>中的元素添加到调用矩阵的相应条目中，使用<tt>indices</tt>为矩阵的行和列指定的本地到全球索引。这个函数假设了一个二次稀疏矩阵和一个二次全矩阵，这是FE计算中的通常情况。
     * 就像deal.II
     * SparseMatrix<Number>类中的相应调用一样（但与基于PETSc的矩阵的情况不同），如果稀疏模式中不存在一个条目，这个函数会抛出一个异常。
     * 可选参数<tt>elide_zero_values</tt>可以用来指定是无论如何都要添加零值，还是要过滤掉这些零值，只添加非零数据。默认值是<tt>true</tt>，也就是说，零值不会被添加到矩阵中。
     *
     */
    void
    add(const std::vector<size_type> &    indices,
        const FullMatrix<TrilinosScalar> &full_matrix,
        const bool                        elide_zero_values = true);

    /**
     * 与之前的函数相同，但现在包括了使用矩形full_matrices的可能性，以及在行和列上分别使用不同的本地到全球索引。
     *
     */
    void
    add(const std::vector<size_type> &    row_indices,
        const std::vector<size_type> &    col_indices,
        const FullMatrix<TrilinosScalar> &full_matrix,
        const bool                        elide_zero_values = true);

    /**
     * 将矩阵的指定行中的几个元素与<tt>col_indices</tt>给出的列索引设置为相应的值。
     * 就像deal.II
     * SparseMatrix<Number>类中的相应调用一样（但与基于PETSc的矩阵的情况不同），如果稀疏模式中不存在一个条目，这个函数会抛出一个异常。
     * 可选参数<tt>elide_zero_values</tt>可以用来指定是无论如何都要添加零值，还是要过滤掉这些零值，只添加非零数据。默认值是<tt>true</tt>，也就是说，零值不会被添加到矩阵中。
     *
     */
    void
    add(const size_type                    row,
        const std::vector<size_type> &     col_indices,
        const std::vector<TrilinosScalar> &values,
        const bool                         elide_zero_values = true);

    /**
     * 在给定的全局矩阵行中，在稀疏矩阵中由col_indices指定的列中添加一个由<tt>values</tt>给出的数值阵列。
     * 就像deal.II
     * SparseMatrix<Number>类中的相应调用一样（但与基于PETSc的矩阵的情况不同），如果在稀疏模式中不存在条目，这个函数会抛出一个异常。
     * 可选参数<tt>elide_zero_values</tt>可以用来指定是无论如何都要添加零值，还是要过滤掉这些零值，只添加非零数据。默认值是<tt>true</tt>，也就是说，零值不会被添加到矩阵中。
     *
     */
    void
    add(const size_type       row,
        const size_type       n_cols,
        const size_type *     col_indices,
        const TrilinosScalar *values,
        const bool            elide_zero_values      = true,
        const bool            col_indices_are_sorted = false);

    /**
     * 将整个矩阵乘以一个固定系数。
     *
     */
    SparseMatrix &
    operator*=(const TrilinosScalar factor);

    /**
     * 用整个矩阵除以一个固定系数。
     *
     */
    SparseMatrix &
    operator/=(const TrilinosScalar factor);

    /**
     * 复制给定的（Trilinos）矩阵（稀疏模式和条目）。
     *
     */
    void
    copy_from(const SparseMatrix &source);

    /**
     * 将<tt>matrix</tt>按<tt>factor</tt>的比例添加到该矩阵中，即<tt>factor*matrix</tt>的矩阵被添加到<tt>this</tt>。如果调用矩阵的稀疏性模式不包含输入矩阵的稀疏性模式中的所有元素，这个函数将抛出一个异常。
     *
     */
    void
    add(const TrilinosScalar factor, const SparseMatrix &matrix);

    /**
     * 将此<tt>行</tt>中的所有元素设置为零，将其删除。这个函数并不修改分配的非零条目的数量，它只是将这些条目设置为零。
     * 这个操作用于消除约束（例如由于挂起的节点），并确保我们可以将这个修改写入矩阵，而不必从矩阵中读取条目（例如非零元素的位置）&mdash；如果没有这个操作，消除%并行矩阵的约束是一个相当复杂的程序。
     * 第二个参数可以用来将该行的对角线条目设置为一个不同于零的值。默认是将其设置为零。
     * @note
     * 如果矩阵是使用MPI在多个处理器之间并行存储的，这个函数只触及本地存储的行，而简单地忽略所有其他行的索引。此外，在并行计算的背景下，如果你清除了某一行，而其他处理器对同一行仍有待处理的写入或添加，你会陷入麻烦。换句话说，如果另一个处理器仍然想向某行的某个元素添加东西，而你调用这个函数将该行清零，那么你下次调用compress()时可能会将远程值添加到你刚刚创建的零上。因此，你要在对矩阵进行最后一次修改后，在开始清空行之前调用compress()。
     *
     */
    void
    clear_row(const size_type row, const TrilinosScalar new_diag_value = 0);

    /**
     * 与clear_row()相同，不同的是，它一次对若干行进行操作。
     * 第二个参数可以用来将所有被清除的行的对角线项设置为不同于0的内容。请注意，所有这些对角线项都得到相同的值
     *
     * - 如果你想要不同的对角线条目的值，你必须手动设置它们。
     * @note
     * 如果矩阵是用MPI在多个处理器之间并行存储的，这个函数只触及本地存储的行，而简单地忽略所有其他行的索引。此外，在并行计算的背景下，如果你清除了某一行，而其他处理器对同一行仍有待处理的写入或添加，你会陷入麻烦。换句话说，如果另一个处理器仍然想向某行的某个元素添加东西，而你调用这个函数将该行清零，那么你下次调用compress()时可能会将远程值添加到你刚刚创建的零上。因此，你要在对矩阵进行最后一次修改后，在开始清空行之前调用compress()。
     *
     */
    void
    clear_rows(const std::vector<size_type> &rows,
               const TrilinosScalar          new_diag_value = 0);

    /**
     * 设置一个内部标志，使矩阵进行的所有操作，即乘法，都以转置的顺序进行。然而，这并不能直接将矩阵重塑为转置的形式，所以在使用这个标志时应该注意。
     * @note
     * 连续调用此函数的任何偶数次，都将使对象返回到其原始状态。
     *
     */
    void
    transpose();

    //@}
    /**
     * @name  入口访问
     *
     */
    //@{

    /**
     * 返回条目的值（<i>i,j</i>）。
     * 这可能是一个昂贵的操作，你应该始终注意在哪里调用这个函数。
     * 正如在deal.II稀疏矩阵类中，如果相应的条目不存在于该类的稀疏模式中，我们会抛出一个异常，这是从Trilinos中请求的。此外，如果要求的元素没有保存在调用过程中，也会抛出一个异常。
     *
     */
    TrilinosScalar
    operator()(const size_type i, const size_type j) const;

    /**
     * 返回矩阵条目的值（<i>i,j</i>）。如果这个条目在稀疏模式中不存在，那么将返回0。虽然这在某些情况下可能很方便，但请注意，由于没有使用矩阵的稀疏性，写出的算法与最优解相比很简单，很慢。
     * 另一方面，如果你想确定条目存在，你应该使用operator()来代替。
     * 如果你有一个并行的矩阵，这个函数中缺乏错误检查，也会产生令人惊讶的结果。在这种情况下，你从这个函数中得到一个零的结果，并不意味着该条目在稀疏模式中不存在，或者它存在但数值为零。相反，也可能是它根本就没有存储在当前的处理器上；在这种情况下，它可能被存储在另一个处理器上，而且可能是以非零值存储的。
     *
     */
    TrilinosScalar
    el(const size_type i, const size_type j) const;

    /**
     * 返回第<i>i</i>行的主对角线元素。如果矩阵不是二次方的，这个函数会抛出一个错误，如果<i>(i,i)</i>不是本地矩阵的元素，它也会抛出一个错误。
     * 参见trilinos_sparse_matrix.cc中的注释。
     *
     */
    TrilinosScalar
    diag_element(const size_type i) const;

    //@}
    /**
     * @name  乘法运算
     *
     */
    //@{

    /**
     * 矩阵-向量乘法：让<i>dst = M*src</i>与<i>M</i>为该矩阵。        源和目的不能是同一个向量。        这个函数可以用几种类型的向量对象来调用，即 @p VectorType 可以是 <ul>   <li>   TrilinosWrappers::MPI::Vector,   <li>   LinearAlgebra::EpetraWrappers::Vector,   <li>   LinearAlgebra::TpetraWrappers::Vector,   <li>  Vector<double>， <li>   LinearAlgebra::distributed::Vector<double>.   </ul>  当使用类型为 TrilinosWrappers::MPI::Vector,  ]的向量 @p dst 必须用用于矩阵行索引的相同IndexSet进行初始化，向量 @p src 必须用用于矩阵列索引的相同IndexSet进行初始化。        当矩阵对象和向量对象的底层数字类型相同时，这个函数将被调用。    尽管看起来很复杂，但其返回类型只是 "void"。        如果是一个串行向量，这个函数只有在一个处理器上运行时才会起作用，因为矩阵对象本身就是分布的。否则，会产生一个异常。
     *
     */
    template <typename VectorType>
    typename std::enable_if<std::is_same<typename VectorType::value_type,
                                         TrilinosScalar>::value>::type
    vmult(VectorType &dst, const VectorType &src) const;

    /**
     * 与上面的函数相同，适用于矩阵对象的底层数字类型和向量对象的数字类型不一致的情况。
     * 这种情况没有实现。调用它将导致一个运行时错误。
     * 尽管看起来很复杂，但其返回类型只是 "void"。
     *
     */
    template <typename VectorType>
    typename std::enable_if<!std::is_same<typename VectorType::value_type,
                                          TrilinosScalar>::value>::type
    vmult(VectorType &dst, const VectorType &src) const;

    /**
     * 矩阵-向量乘法：让<i>dst =
     * M<sup>T</sup>*src</i>与<i>M</i>是这个矩阵。这个函数与vmult()的作用相同，但需要转置的矩阵。
     * 源和目的不能是同一个向量。
     * 这个函数可以用几种类型的向量对象调用，见vmult()中关于
     * @p VectorType 的讨论。
     * 当矩阵对象的底层数字类型和向量对象的底层数字类型相同时，该函数将被调用。
     * 尽管看起来很复杂，但其返回类型只是 "void"。
     *
     */
    template <typename VectorType>
    typename std::enable_if<std::is_same<typename VectorType::value_type,
                                         TrilinosScalar>::value>::type
    Tvmult(VectorType &dst, const VectorType &src) const;

    /**
     * 与上面的函数相同，适用于矩阵对象的底层数字类型和向量对象的底层数字类型不一致的情况。
     * 这种情况没有实现。调用它将导致一个运行时错误。
     * 尽管看起来很复杂，但其返回类型只是 "void"。
     *
     */
    template <typename VectorType>
    typename std::enable_if<!std::is_same<typename VectorType::value_type,
                                          TrilinosScalar>::value>::type
    Tvmult(VectorType &dst, const VectorType &src) const;

    /**
     * 添加矩阵-向量的乘法。在<i>dst</i>上添加<i>M*src</i>，<i>M</i>为该矩阵。
     * 源和目的不能是同一个向量。
     * 这个函数可以用几种类型的向量对象调用，见vmult()中关于
     * @p VectorType 的讨论。
     *
     */
    template <typename VectorType>
    void
    vmult_add(VectorType &dst, const VectorType &src) const;

    /**
     * 添加矩阵-向量的乘法。将<i>M<sup>T</sup>*src</i>加到<i>dst</i>，<i>M</i>是这个矩阵。这个函数的作用与vmult_add()相同，但取的是转置的矩阵。
     * 来源和目的地不能是同一个向量。
     * 这个函数可以用几种类型的向量对象调用，见vmult()中关于
     * @p VectorType 的讨论。
     *
     */
    template <typename VectorType>
    void
    Tvmult_add(VectorType &dst, const VectorType &src) const;

    /**
     * 返回向量 $v$ 相对于该矩阵引起的规范的平方，即
     * $\left(v,Mv\right)$
     * 。这很有用，例如在有限元背景下，一个函数的 $L_2$
     * 规范等于相对于代表有限元函数节点值的向量的质量矩阵的矩阵规范。
     * 很明显，对于这个操作，矩阵需要是二次的。
     * 这个函数的实现没有deal.II中使用的 @p SparseMatrix
     * 类（即原始的，而不是Trilinos包装类）的效率高，因为Trilinos不支持这个操作，需要一个临时向量。
     * 矢量必须用矩阵初始化的相同IndexSet进行初始化。
     * 如果是一个本地化的Vector，这个函数只有在一个处理器上运行时才会起作用，因为矩阵对象本身就是分布的。否则，将抛出一个异常。
     *
     */
    TrilinosScalar
    matrix_norm_square(const MPI::Vector &v) const;

    /**
     * 计算矩阵标量乘积  $\left(u,Mv\right)$  。
     * 这个函数的实现没有deal.II中使用的 @p SparseMatrix
     * 类（即原始的，而不是Trilinos包装类）的效率高，因为Trilinos不支持这个操作，需要一个临时矢量。
     * 矢量 @p u
     * 必须用用于矩阵行索引的相同IndexSet进行初始化，矢量
     * @p v 必须用用于矩阵列索引的相同IndexSet进行初始化。
     * 如果是一个本地化的Vector，这个函数只有在一个处理器上运行时才会起作用，因为矩阵对象本身是分布式的。否则，将抛出一个异常。
     * 这个函数只对方形矩阵实现。
     *
     */
    TrilinosScalar
    matrix_scalar_product(const MPI::Vector &u, const MPI::Vector &v) const;

    /**
     * 计算方程<i>Mx=b</i>的残差，其中残差被定义为<i>r=b-Mx</i>。将残差写入
     * @p dst.  返回残差向量的<i>l<sub>2</sub></i>准则。
     * 源<i>x</i>和目的<i>dst</i>不能是同一个向量。        向量
     * @p dst 和 @p b
     * 必须用用于矩阵行索引的相同IndexSet进行初始化，向量
     * @p x 必须用用于矩阵列索引的相同IndexSet进行初始化。
     * 如果是一个本地化的Vector，这个函数只有在一个处理器上运行时才会起作用，因为矩阵对象本身是分布式的。否则，将抛出一个异常。
     *
     */
    TrilinosScalar
    residual(MPI::Vector &      dst,
             const MPI::Vector &x,
             const MPI::Vector &b) const;

    /**
     * 执行矩阵-矩阵乘法<tt>C = A
     * B</tt>，或者，如果给出一个可选的矢量参数，<tt>C = A
     * diag(V)
     * B</tt>，其中<tt>diag(V)</tt>定义了一个带有矢量条目的对角矩阵。
     * 这个函数假定调用矩阵<tt>A</tt>和<tt>B</tt>的大小是兼容的。<tt>C</tt>的大小将在本函数中设置。
     * 矩阵C的内容和稀疏模式将被这个函数改变，所以要确保稀疏模式没有在你的程序中其他地方使用。这是一个昂贵的操作，所以在使用这个函数之前请三思。
     *
     */
    void
    mmult(SparseMatrix &      C,
          const SparseMatrix &B,
          const MPI::Vector & V = MPI::Vector()) const;


    /**
     * 用<tt>this</tt>的转置进行矩阵-矩阵相乘，即<tt>C =
     * A<sup>T</sup>
     * B</tt>，或者，如果给出了可选的矢量参数，<tt>C =
     * A<sup>T</sup> diag(V)
     * B</tt>，其中<tt>diag(V)</tt>定义了一个带有矢量项的对角矩阵。
     * 这个函数假定调用矩阵<tt>A</tt>和<tt>B</tt>的大小兼容。<tt>C</tt>的大小将在本函数中设置。
     * 矩阵C的内容和稀疏模式将被这个函数改变，所以要确保稀疏模式没有在你的程序中其他地方使用。这是一个昂贵的操作，所以在使用这个函数之前请三思。
     *
     */
    void
    Tmmult(SparseMatrix &      C,
           const SparseMatrix &B,
           const MPI::Vector & V = MPI::Vector()) const;

    //@}
    /**
     * @name  矩阵规范
     *
     */
    //@{

    /**
     * 返回矩阵的<i>l</i><sub>1</sub>规范，即 $|M|_1=
     * \max_{\mathrm{all\ columns\ } j} \sum_{\mathrm{all\ rows\ } i}
     * |M_{ij}|$  ，（最大列之和）。
     * 这是自然的矩阵准则，与向量的l1准则兼容，即 $|Mv|_1
     * \leq |M|_1 |v|_1$  。 (参看Haemmerlin-Hoffmann: Numerische
     * Mathematik)
     *
     */
    TrilinosScalar
    l1_norm() const;

    /**
     * 返回矩阵的linfty-norm，即 $|M|_\infty=\max_{\mathrm{all\ rows\ }
     * i}\sum_{\mathrm{all\ columns\ } j} |M_{ij}|$  , (行的最大和)。
     * 这是一个自然的矩阵规范，与向量的linfty-norm兼容，即
     * $|Mv|_\infty \leq |M|_\infty |v|_\infty$  。
     * (参看Haemmerlin-Hoffmann: Numerische Mathematik)
     *
     */
    TrilinosScalar
    linfty_norm() const;

    /**
     * 返回矩阵的frobenius
     * norm，即矩阵中所有条目的平方之和的平方根。
     *
     */
    TrilinosScalar
    frobenius_norm() const;

    //@}
    /**
     * @name  访问底层Trilinos数据
     *
     */
    //@{

    /**
     * 返回对底层Trilinos Epetra_CrsMatrix数据的一个常量引用。
     *
     */
    const Epetra_CrsMatrix &
    trilinos_matrix() const;

    /**
     * 返回一个对底层Trilinos
     * Epetra_CrsGraph数据的常量引用，该数据存储了矩阵的稀疏性模式。
     *
     */
    const Epetra_CrsGraph &
    trilinos_sparsity_pattern() const;

    //@}

    /**
     * @name  分割器
     *
     */
    //@{

    /**
     * 返回该矩阵的域空间的分区，即该矩阵要与之相乘的向量的分区。
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

    //@}

    /**
     * @name  迭代器
     *
     */
    //@{

    /**
     * 返回一个指向矩阵的第一个元素的迭代器。
     * 迭代器在每一行中访问的元素是按照Trilinos的存储方式来排序的，尽管实现上保证了一行的所有元素在下一行的元素之前被访问。如果你的算法依赖于访问一行中的元素，你将需要咨询Trilinos文档，了解它存储数据的顺序。然而，如果你对元素进行迭代，依靠接收元素的顺序通常不是一个好的和长期稳定的想法。
     * 当你遍历一个并行矩阵的元素时，你将只能访问本地拥有的行。(你也可以访问其他行，但它们看起来是空的。)在这种情况下，你可能想调用begin()函数，该函数将行作为一个参数，以限制要循环的元素范围。
     *
     */
    const_iterator
    begin() const;

    /**
     * 像上面的函数一样，但对于非恒定矩阵。
     *
     */
    iterator
    begin();

    /**
     * 返回一个迭代器，指向这个矩阵的最后一个以上的元素。
     *
     */
    const_iterator
    end() const;

    /**
     * 像上面的函数一样，但对于非静态矩阵。
     *
     */
    iterator
    end();

    /**
     * 返回一个指向行 @p r. 第一个元素的迭代器
     * 注意，如果给定的行是空的，即不包含任何非零条目，那么这个函数返回的迭代器就等于<tt>end(r)</tt>。在这种情况下，如果行
     * @p r
     * 和下面的任何一行都不包含任何非零条目，那么返回的迭代器可能无法被解除引用。
     * 迭代器在每一行中访问的元素是按照Trilinos的存储方式排序的，尽管实现上保证了一行的所有元素在下一行的元素之前被访问。如果你的算法依赖于访问一行中的元素，你将需要咨询Trilinos文档，了解它存储数据的顺序。然而，如果你对元素进行迭代，依靠接收元素的顺序通常不是一个好的和长期稳定的想法。
     * @note
     * 当你访问一个并行矩阵的元素时，你只能访问实际存储在本地的行的元素。(你也可以访问其他行，但它们看起来是空的。)即使如此，如果另一个处理器后来写进或增加了存储在当前处理器上的矩阵元素，那么你仍然会看到这个条目的旧值，除非你在远程处理器上修改矩阵元素和在当前处理器上访问它之间调用compress()。更多信息请参见compress()函数的文档。
     *
     */
    const_iterator
    begin(const size_type r) const;

    /**
     * 像上面的函数一样，但是对于非恒定矩阵。
     *
     */
    iterator
    begin(const size_type r);

    /**
     * 返回一个指向第 @p r
     * 行最后一个元素的迭代器，如果第 @p r
     * 行之后没有任何条目，则指向整个稀疏模式的末端。
     * 请注意，结束迭代器不一定是可被解读的。特别是如果它是一个矩阵的最后一行的结束迭代器，情况更是如此。
     *
     */
    const_iterator
    end(const size_type r) const;

    /**
     * 像上面的函数一样，但是对于非恒定矩阵。
     *
     */
    iterator
    end(const size_type r);

    //@}
    /**
     * @name  输入/输出
     *
     */
    //@{

    /**
     * 抽象的Trilinos对象，帮助以ASCII格式查看其他Trilinos对象。目前这个功能还没有实现。
     * TODO：没有实现。
     *
     */
    void
    write_ascii();

    /**
     * 打印矩阵到给定的流中，使用格式<tt>(line,col)
     * value</tt>，即每行一个非零的矩阵条目。可选的标志以Trilinos风格输出稀疏模式，在打印到流中时，数据根据处理器编号进行排序，以及矩阵的摘要，如全局大小。
     *
     */
    void
    print(std::ostream &out,
          const bool    write_extended_trilinos_info = false) const;

    //@}
    /**
     * @addtogroup  Exceptions
     *
     */
    //@{
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
    DeclException2(ExcInvalidIndex,
                   size_type,
                   size_type,
                   << "The entry with index <" << arg1 << ',' << arg2
                   << "> does not exist.");

    /**
     * 异常情况
     *
     */
    DeclExceptionMsg(ExcSourceEqualsDestination,
                     "You are attempting an operation on two matrices that "
                     "are the same object, but the operation requires that the "
                     "two objects are in fact different.");

    /**
     * 异常情况
     *
     */
    DeclException0(ExcMatrixNotCompressed);

    /**
     * 异常情况
     *
     */
    DeclException4(ExcAccessToNonLocalElement,
                   size_type,
                   size_type,
                   size_type,
                   size_type,
                   << "You tried to access element (" << arg1 << "/" << arg2
                   << ")"
                   << " of a distributed matrix, but only rows in range ["
                   << arg3 << "," << arg4
                   << "] are stored locally and can be accessed.");

    /**
     * 异常情况
     *
     */
    DeclException2(ExcAccessToNonPresentElement,
                   size_type,
                   size_type,
                   << "You tried to access element (" << arg1 << "/" << arg2
                   << ")"
                   << " of a sparse matrix, but it appears to not"
                   << " exist in the Trilinos sparsity pattern.");
    //@}



  protected:
    /**
     * 对于某些矩阵存储格式，特别是PETSc分布式块矩阵，对单个元素的设置和添加操作不能自由混合。相反，当我们想从设置元素切换到添加元素时，我们必须同步操作。
     * BlockMatrixBase通过为每个块调用这个辅助函数来自动同步访问。
     * 这个函数确保矩阵处于一个允许添加元素的状态；如果它之前已经处于这个状态，那么这个函数就不会做任何事情。
     *
     */
    void
    prepare_add();

    /**
     * 与prepare_add()相同，但如果该类中的元素表示法需要这样的操作，则为设置元素准备矩阵。
     *
     */
    void
    prepare_set();



  private:
    /**
     * 指向用户提供的矩阵列的Epetra
     * Trilinos映射的指针，该映射将矩阵的一部分分配给各个进程。
     *
     */
    std::unique_ptr<Epetra_Map> column_space_map;

    /**
     * Trilinos中的一个稀疏矩阵对象，用于基于有限元的问题，允许组装成非局部元素。
     * 实际的类型是稀疏矩阵，在构造函数中设置。
     *
     */
    std::unique_ptr<Epetra_FECrsMatrix> matrix;

    /**
     * Trilinos中的一个稀疏矩阵对象，用于收集非局部元素，如果该矩阵是由Trilinos稀疏模式构建的，并有相应的选项。
     *
     */
    std::unique_ptr<Epetra_CrsMatrix> nonlocal_matrix;

    /**
     * 一个用于交流非本地矩阵的导出对象。
     *
     */
    std::unique_ptr<Epetra_Export> nonlocal_matrix_exporter;

    /**
     * Trilinos不允许混合添加矩阵条目和覆盖它们（以使%并行计算的同步更简单）。我们的方法是，对于每个访问操作，存储它是插入还是添加。如果前一个是不同的类型，那么我们首先要刷新Trilinos缓冲区；否则，我们可以简单地继续下去。幸运的是，Trilinos有一个这样的对象，在这种情况下已经完成了所有的%并行通信，所以我们只需使用他们的模型，它存储了上一个操作是加法还是插入。
     *
     */
    Epetra_CombineMode last_action;

    /**
     * 一个布尔变量，用来保存关于向量是否被压缩的信息。
     *
     */
    bool compressed;

    // To allow calling protected prepare_add() and prepare_set().
    friend class BlockMatrixBase<SparseMatrix>;
  };



  // forwards declarations
  class SolverBase;
  class PreconditionBase;

  namespace internal
  {
    inline void
    check_vector_map_equality(const Epetra_CrsMatrix &  mtrx,
                              const Epetra_MultiVector &src,
                              const Epetra_MultiVector &dst,
                              const bool                transpose)
    {
      if (transpose == false)
        {
          Assert(src.Map().SameAs(mtrx.DomainMap()) == true,
                 ExcMessage(
                   "Column map of matrix does not fit with vector map!"));
          Assert(dst.Map().SameAs(mtrx.RangeMap()) == true,
                 ExcMessage("Row map of matrix does not fit with vector map!"));
        }
      else
        {
          Assert(src.Map().SameAs(mtrx.RangeMap()) == true,
                 ExcMessage(
                   "Column map of matrix does not fit with vector map!"));
          Assert(dst.Map().SameAs(mtrx.DomainMap()) == true,
                 ExcMessage("Row map of matrix does not fit with vector map!"));
        }
      (void)mtrx; // removes -Wunused-variable in optimized mode
      (void)src;
      (void)dst;
    }

    inline void
    check_vector_map_equality(const Epetra_Operator &   op,
                              const Epetra_MultiVector &src,
                              const Epetra_MultiVector &dst,
                              const bool                transpose)
    {
      if (transpose == false)
        {
          Assert(src.Map().SameAs(op.OperatorDomainMap()) == true,
                 ExcMessage(
                   "Column map of operator does not fit with vector map!"));
          Assert(dst.Map().SameAs(op.OperatorRangeMap()) == true,
                 ExcMessage(
                   "Row map of operator does not fit with vector map!"));
        }
      else
        {
          Assert(src.Map().SameAs(op.OperatorRangeMap()) == true,
                 ExcMessage(
                   "Column map of operator does not fit with vector map!"));
          Assert(dst.Map().SameAs(op.OperatorDomainMap()) == true,
                 ExcMessage(
                   "Row map of operator does not fit with vector map!"));
        }
      (void)op; // removes -Wunused-variable in optimized mode
      (void)src;
      (void)dst;
    }

    namespace LinearOperatorImplementation
    {
      /**
       * 这是对Trilinos稀疏矩阵和预处理类型的LinearOperators的一个扩展类。它提供了对Trilinos向量类型进行基本操作（<tt>vmult</tt>和<tt>Tvmult</tt>）的接口。它满足了将调用Epetra_Operator函数的Trilinos求解器包装成LinearOperator的必要条件。
       * @note  这个有效载荷所包装的 TrilinosWrappers::SparseMatrix
       * 或 TrilinosWrappers::PreconditionBase
       * 是通过引用传递给<tt>vmult</tt>和<tt>Tvmult</tt>函数。当它或它所引用的Trilinos对象上设置了转置标志时，这个对象不是线程安全的。更多细节请参见
       * TrilinosWrappers::internal::LinearOperatorImplementation::TrilinosPayload::SetUseTranspose()
       * 函数的文档。
       * @ingroup TrilinosWrappers
       *
       */
      class TrilinosPayload : public Epetra_Operator
      {
      public:
        /**
         * 内部支持的向量类型的定义。
         *
         */
        using VectorType = Epetra_MultiVector;

        /**
         * 运算器领域空间的向量类型的定义。
         *
         */
        using Range = VectorType;

        /**
         * 运算器的范围空间的向量类型的定义。
         *
         */
        using Domain = VectorType;

        /**
         * @name  构造函数/解构函数
         *
         */
        //@{

        /**
         * 默认构造函数
         * @note
         * 根据设计，由于没有足够的信息来构建域和范围图，因此产生的对象是不可操作的。
         *
         */
        TrilinosPayload();

        /**
         * 基于典范矩阵的稀疏矩阵的构造函数
         *
         */
        TrilinosPayload(const TrilinosWrappers::SparseMatrix &matrix_exemplar,
                        const TrilinosWrappers::SparseMatrix &matrix);

        /**
         * 基于典范矩阵的预处理器的构造器
         *
         */
        TrilinosPayload(
          const TrilinosWrappers::SparseMatrix &    matrix_exemplar,
          const TrilinosWrappers::PreconditionBase &preconditioner);

        /**
         * 基于示范性预处理的预处理器的构造器
         *
         */
        TrilinosPayload(
          const TrilinosWrappers::PreconditionBase &preconditioner_exemplar,
          const TrilinosWrappers::PreconditionBase &preconditioner);

        /**
         * 默认的复制构造器
         *
         */
        TrilinosPayload(const TrilinosPayload &payload);

        /**
         * 复合复制构造函数
         * 这是PackagedOperations所需要的，因为它设置了域和范围图，以及基于两个操作的复合<tt>vmult</tt>和<tt>Tvmult</tt>操作的组合操作
         *
         */
        TrilinosPayload(const TrilinosPayload &first_op,
                        const TrilinosPayload &second_op);

        /**
         * 破坏器
         *
         */
        virtual ~TrilinosPayload() override = default;

        /**
         * 返回一个为身份操作配置的有效载荷
         *
         */
        TrilinosPayload
        identity_payload() const;

        /**
         * 返回一个为空操作而配置的有效载荷
         *
         */
        TrilinosPayload
        null_payload() const;

        /**
         * 返回一个为转置操作配置的有效载荷
         *
         */
        TrilinosPayload
        transpose_payload() const;

        /**
         * 返回一个为反转操作配置的有效载荷
         * 调用这个工厂函数将配置两个额外的函数，即<tt>inv_vmult</tt>和<tt>inv_Tvmult</tt>，这两个函数都包裹了反转操作。<tt>vmult</tt>和<tt>Tvmult</tt>操作保留了从
         * @p op. 继承的标准定义。
         * @note
         * 只有当求解器和预处理器派生自各自的TrilinosWrappers基类时，该功能才会启用。
         * 因此，C++编译器只有在满足以下标准的情况下才会考虑这个功能。
         * 1.  @p Solver 派生自 TrilinosWrappers::SolverBase, ，2.  @p
         * Preconditioner 派生自 TrilinosWrappers::PreconditionBase.  。
         *
         */
        template <typename Solver, typename Preconditioner>
        typename std::enable_if<
          std::is_base_of<TrilinosWrappers::SolverBase, Solver>::value &&
            std::is_base_of<TrilinosWrappers::PreconditionBase,
                            Preconditioner>::value,
          TrilinosPayload>::type
        inverse_payload(Solver &, const Preconditioner &) const;

        /**
         * 返回一个为逆运算配置的有效载荷
         * 调用这个工厂函数将配置两个额外的函数，即<tt>inv_vmult</tt>和<tt>inv_Tvmult</tt>，这两个函数被禁用，因为
         * @p Solver 或 @p Preconditioner 与Epetra_MultiVector不兼容。
         * <tt>vmult</tt>和<tt>Tvmult</tt>操作保留了从 @p op.
         * 继承的标准定义。
         * @note
         * 只有满足以下标准，C++编译器才会考虑这个函数。
         * 1.  @p Solver 不从 TrilinosWrappers::SolverBase, 派生，2.  @p
         * Preconditioner 不从 TrilinosWrappers::PreconditionBase. 派生。
         *
         */
        template <typename Solver, typename Preconditioner>
        typename std::enable_if<
          !(std::is_base_of<TrilinosWrappers::SolverBase, Solver>::value &&
            std::is_base_of<TrilinosWrappers::PreconditionBase,
                            Preconditioner>::value),
          TrilinosPayload>::type
        inverse_payload(Solver &, const Preconditioner &) const;

        //@}

        /**
         * @name  LinearOperator功能
         *
         */
        //@{

        /**
         * 返回一个IndexSet，它定义了这个矩阵的域空间的分区，即这个矩阵要与之相乘/操作的向量的分区。
         *
         */
        IndexSet
        locally_owned_domain_indices() const;

        /**
         * 返回一个IndexSet，它定义了这个矩阵的范围空间的划分，即由矩阵-向量积产生的向量的划分。
         *
         */
        IndexSet
        locally_owned_range_indices() const;

        /**
         * 返回与该有效载荷一起使用的MPI通信器对象。
         *
         */
        MPI_Comm
        get_mpi_communicator() const;

        /**
         * 设置一个内部标志，使矩阵进行的所有操作，即乘法，都以转置的顺序进行。
         * @note
         * 这并不直接将矩阵重塑为转置形式，所以在使用这个标志时要注意。
         *
         */
        void
        transpose();

        /**
         * 当Apply被调用时，将由有效载荷执行的标准矩阵-向量操作。
         * @note
         * 这不是由LinearOperator调用的，而是由Trilinos函数调用的，这些函数希望以此来模仿LinearOperator的动作。
         *
         */
        std::function<void(VectorType &, const VectorType &)> vmult;

        /**
         * 当Apply被调用时，标准的转置矩阵-向量操作将由有效载荷执行。
         * @note
         * 这不是由LinearOperator调用的，而是由Trilinos函数调用的，该函数希望以此来模仿LinearOperator的动作。
         *
         */
        std::function<void(VectorType &, const VectorType &)> Tvmult;

        /**
         * 当ApplyInverse被调用时，有效载荷将进行矩阵-向量的逆运算。
         * @note
         * 这不是由LinearOperator调用的，而是由Trilinos函数调用的，这些函数希望以此来模仿InverseOperator的动作。
         *
         */
        std::function<void(VectorType &, const VectorType &)> inv_vmult;

        /**
         * 当ApplyInverse被调用时，将由有效载荷执行的矩阵-向量的反转操作。
         * @note
         * 这不是由LinearOperator调用的，而是由Trilinos函数调用的，该函数期望以此来模仿InverseOperator的动作。
         *
         */
        std::function<void(VectorType &, const VectorType &)> inv_Tvmult;

        //@}

        /**
         * @name  核心Epetra_Operator功能
         *
         */
        //@{

        /**
         * 返回该运算符的转置标志状态
         * 这是对Trilinos类Epetra_Operator的同一函数的重载。
         *
         */
        virtual bool
        UseTranspose() const override;

        /**
         * 设置一个内部标志，使矩阵进行的所有操作，即乘法，都以转置的顺序进行。
         * 这是对Trilinos类Epetra_Operator的相同函数的重载。
         * @note
         * 这并不直接将矩阵重塑为转置形式，所以在使用这个标志时要注意。当该标志被设置为
         * "true
         * "时（无论是在这里还是直接在底层的Tridinos对象本身），该对象不再是线程安全的。从本质上讲，它不可能确保LinearOperator和底层Trilinos对象的转置状态在可能同时发生在不同线程上的所有操作中保持同步。
         *
         */
        virtual int
        SetUseTranspose(bool UseTranspose) override;

        /**
         * 对一个向量 @p X
         * （内部定义的VectorType类型）应用vmult操作，并将结果存储在向量中
         * @p Y.
         * 这是从Trilinos类Epetra_Operator中重载的相同函数。
         * @note
         * 预定的操作取决于内部转置标志的状态。如果这个标志被设置为
         * "真"，那么结果将相当于执行一个Tvmult操作。
         *
         */
        virtual int
        Apply(const VectorType &X, VectorType &Y) const override;

        /**
         * 对一个向量 @p X
         * （内部定义的VectorType类型）应用vmult逆运算，并将结果存储在向量
         * @p Y. 中。
         * 实际上，只有当包裹的对象作为预处理程序时，才会从特里诺斯求解器调用这个函数。
         * 这是从Tridinos类Epetra_Operator中重载的相同函数。
         * @note
         * 只有当有效载荷被InverseOperator初始化，或者是一个预处理程序的包装器时，这个函数才可操作。如果没有，那么使用这个函数将导致抛出一个错误。
         * @note
         * 预期的操作取决于内部转置标志的状态。如果这个标志被设置为
         * "true"，那么结果将等同于执行一个Tvmult操作。
         *
         */
        virtual int
        ApplyInverse(const VectorType &Y, VectorType &X) const override;
        //@}

        /**
         * @name  附加Epetra_Operator功能
         *
         */
        //@{

        /**
         * 返回一个标签来描述这个类。
         * 这重载了Trilinos类Epetra_Operator中的相同功能。
         *
         */
        virtual const char *
        Label() const override;

        /**
         * 返回这个对象的底层MPI通信器的引用。
         * 这是从Trilinos类Epetra_Operator中重载的相同函数。
         *
         */
        virtual const Epetra_Comm &
        Comm() const override;

        /**
         * 返回该矩阵的域空间的划分，即该矩阵必须与之相乘的向量的划分。
         * 这是对Trilinos类Epetra_Operator的相同函数的重载。
         *
         */
        virtual const Epetra_Map &
        OperatorDomainMap() const override;

        /**
         * 返回该矩阵的范围空间的划分，即由矩阵-向量乘积产生的向量的划分。
         * 这是对Trilinos类Epetra_Operator的相同函数的重载。
         *
         */
        virtual const Epetra_Map &
        OperatorRangeMap() const override;
        //@}

      private:
        /**
         * 一个记录运算器是执行标准的矩阵-向量乘法，还是转置运算的标志。
         *
         */
        bool use_transpose;

        /**
         * 内部通信模式，以防矩阵需要从deal.II格式中复制。
         *
         */
#    ifdef DEAL_II_WITH_MPI
        Epetra_MpiComm communicator;
#    else
        Epetra_SerialComm communicator;
#    endif

        /**
         * Epetra_Map，用于设置该运算器的域空间的划分。
         *
         */
        Epetra_Map domain_map;

        /**
         * Epetra_Map，设置此运算符的范围空间的划分。
         *
         */
        Epetra_Map range_map;

        /**
         * 返回一个标志，描述该算子是否可以返回无穷大准则的计算。因为一般情况下不是这样的，所以总是返回一个否定的结果。
         * 这是对Trilinos类Epetra_Operator的相同函数的重载。
         *
         */
        virtual bool
        HasNormInf() const override;

        /**
         * 返回该运算符的无穷大规范。
         * 抛出一个错误，因为一般来说，我们无法计算这个值。
         * 这是对Trilinos类Epetra_Operator中相同函数的重载。
         *
         */
        virtual double
        NormInf() const override;
      };

      /**
       * 返回一个运算器，该运算器返回一个配置为支持两个LinearOperator相加的有效载荷
       *
       */
      TrilinosPayload
      operator+(const TrilinosPayload &first_op,
                const TrilinosPayload &second_op);

      /**
       * 返回一个操作符，该操作符返回一个被配置为支持两个LinearOperator的乘法的有效载荷。
       *
       */
      TrilinosPayload operator*(const TrilinosPayload &first_op,
                                const TrilinosPayload &second_op);

    } // namespace LinearOperatorImplementation
  }    /* namespace internal */ 



  // ----------------------- inline and template functions --------------------

#    ifndef DOXYGEN

  namespace SparseMatrixIterators
  {
    inline AccessorBase::AccessorBase(SparseMatrix *matrix,
                                      size_type     row,
                                      size_type     index)
      : matrix(matrix)
      , a_row(row)
      , a_index(index)
    {
      visit_present_row();
    }


    inline AccessorBase::size_type
    AccessorBase::row() const
    {
      Assert(a_row < matrix->m(), ExcBeyondEndOfMatrix());
      return a_row;
    }


    inline AccessorBase::size_type
    AccessorBase::column() const
    {
      Assert(a_row < matrix->m(), ExcBeyondEndOfMatrix());
      return (*colnum_cache)[a_index];
    }


    inline AccessorBase::size_type
    AccessorBase::index() const
    {
      Assert(a_row < matrix->m(), ExcBeyondEndOfMatrix());
      return a_index;
    }


    inline Accessor<true>::Accessor(MatrixType *    matrix,
                                    const size_type row,
                                    const size_type index)
      : AccessorBase(const_cast<SparseMatrix *>(matrix), row, index)
    {}


    template <bool Other>
    inline Accessor<true>::Accessor(const Accessor<Other> &other)
      : AccessorBase(other)
    {}


    inline TrilinosScalar
    Accessor<true>::value() const
    {
      Assert(a_row < matrix->m(), ExcBeyondEndOfMatrix());
      return (*value_cache)[a_index];
    }


    inline Accessor<false>::Reference::Reference(const Accessor<false> &acc)
      : accessor(const_cast<Accessor<false> &>(acc))
    {}


    inline Accessor<false>::Reference::operator TrilinosScalar() const
    {
      return (*accessor.value_cache)[accessor.a_index];
    }

    inline const Accessor<false>::Reference &
    Accessor<false>::Reference::operator=(const TrilinosScalar n) const
    {
      (*accessor.value_cache)[accessor.a_index] = n;
      accessor.matrix->set(accessor.row(),
                           accessor.column(),
                           static_cast<TrilinosScalar>(*this));
      return *this;
    }


    inline const Accessor<false>::Reference &
    Accessor<false>::Reference::operator+=(const TrilinosScalar n) const
    {
      (*accessor.value_cache)[accessor.a_index] += n;
      accessor.matrix->set(accessor.row(),
                           accessor.column(),
                           static_cast<TrilinosScalar>(*this));
      return *this;
    }


    inline const Accessor<false>::Reference &
    Accessor<false>::Reference::operator-=(const TrilinosScalar n) const
    {
      (*accessor.value_cache)[accessor.a_index] -= n;
      accessor.matrix->set(accessor.row(),
                           accessor.column(),
                           static_cast<TrilinosScalar>(*this));
      return *this;
    }


    inline const Accessor<false>::Reference &
    Accessor<false>::Reference::operator*=(const TrilinosScalar n) const
    {
      (*accessor.value_cache)[accessor.a_index] *= n;
      accessor.matrix->set(accessor.row(),
                           accessor.column(),
                           static_cast<TrilinosScalar>(*this));
      return *this;
    }


    inline const Accessor<false>::Reference &
    Accessor<false>::Reference::operator/=(const TrilinosScalar n) const
    {
      (*accessor.value_cache)[accessor.a_index] /= n;
      accessor.matrix->set(accessor.row(),
                           accessor.column(),
                           static_cast<TrilinosScalar>(*this));
      return *this;
    }


    inline Accessor<false>::Accessor(MatrixType *    matrix,
                                     const size_type row,
                                     const size_type index)
      : AccessorBase(matrix, row, index)
    {}


    inline Accessor<false>::Reference
    Accessor<false>::value() const
    {
      Assert(a_row < matrix->m(), ExcBeyondEndOfMatrix());
      return {*this};
    }



    template <bool Constness>
    inline Iterator<Constness>::Iterator(MatrixType *    matrix,
                                         const size_type row,
                                         const size_type index)
      : accessor(matrix, row, index)
    {}


    template <bool Constness>
    template <bool Other>
    inline Iterator<Constness>::Iterator(const Iterator<Other> &other)
      : accessor(other.accessor)
    {}


    template <bool Constness>
    inline Iterator<Constness> &
    Iterator<Constness>::operator++()
    {
      Assert(accessor.a_row < accessor.matrix->m(), ExcIteratorPastEnd());

      ++accessor.a_index;

      // If at end of line: do one
      // step, then cycle until we
      // find a row with a nonzero
      // number of entries.
      if (accessor.a_index >= accessor.colnum_cache->size())
        {
          accessor.a_index = 0;
          ++accessor.a_row;

          while ((accessor.a_row < accessor.matrix->m()) &&
                 ((accessor.matrix->in_local_range(accessor.a_row) == false) ||
                  (accessor.matrix->row_length(accessor.a_row) == 0)))
            ++accessor.a_row;

          accessor.visit_present_row();
        }
      return *this;
    }


    template <bool Constness>
    inline Iterator<Constness>
    Iterator<Constness>::operator++(int)
    {
      const Iterator<Constness> old_state = *this;
      ++(*this);
      return old_state;
    }



    template <bool Constness>
    inline const Accessor<Constness> &Iterator<Constness>::operator*() const
    {
      return accessor;
    }



    template <bool Constness>
    inline const Accessor<Constness> *Iterator<Constness>::operator->() const
    {
      return &accessor;
    }



    template <bool Constness>
    inline bool
    Iterator<Constness>::operator==(const Iterator<Constness> &other) const
    {
      return (accessor.a_row == other.accessor.a_row &&
              accessor.a_index == other.accessor.a_index);
    }



    template <bool Constness>
    inline bool
    Iterator<Constness>::operator!=(const Iterator<Constness> &other) const
    {
      return !(*this == other);
    }



    template <bool Constness>
    inline bool
    Iterator<Constness>::operator<(const Iterator<Constness> &other) const
    {
      return (accessor.row() < other.accessor.row() ||
              (accessor.row() == other.accessor.row() &&
               accessor.index() < other.accessor.index()));
    }


    template <bool Constness>
    inline bool
    Iterator<Constness>::operator>(const Iterator<Constness> &other) const
    {
      return (other < *this);
    }

  } // namespace SparseMatrixIterators



  inline SparseMatrix::const_iterator
  SparseMatrix::begin() const
  {
    return begin(0);
  }



  inline SparseMatrix::const_iterator
  SparseMatrix::end() const
  {
    return const_iterator(this, m(), 0);
  }



  inline SparseMatrix::const_iterator
  SparseMatrix::begin(const size_type r) const
  {
    AssertIndexRange(r, m());
    if (in_local_range(r) && (row_length(r) > 0))
      return const_iterator(this, r, 0);
    else
      return end(r);
  }



  inline SparseMatrix::const_iterator
  SparseMatrix::end(const size_type r) const
  {
    AssertIndexRange(r, m());

    // place the iterator on the first entry
    // past this line, or at the end of the
    // matrix
    for (size_type i = r + 1; i < m(); ++i)
      if (in_local_range(i) && (row_length(i) > 0))
        return const_iterator(this, i, 0);

    // if there is no such line, then take the
    // end iterator of the matrix
    return end();
  }



  inline SparseMatrix::iterator
  SparseMatrix::begin()
  {
    return begin(0);
  }



  inline SparseMatrix::iterator
  SparseMatrix::end()
  {
    return iterator(this, m(), 0);
  }



  inline SparseMatrix::iterator
  SparseMatrix::begin(const size_type r)
  {
    AssertIndexRange(r, m());
    if (in_local_range(r) && (row_length(r) > 0))
      return iterator(this, r, 0);
    else
      return end(r);
  }



  inline SparseMatrix::iterator
  SparseMatrix::end(const size_type r)
  {
    AssertIndexRange(r, m());

    // place the iterator on the first entry
    // past this line, or at the end of the
    // matrix
    for (size_type i = r + 1; i < m(); ++i)
      if (in_local_range(i) && (row_length(i) > 0))
        return iterator(this, i, 0);

    // if there is no such line, then take the
    // end iterator of the matrix
    return end();
  }



  inline bool
  SparseMatrix::in_local_range(const size_type index) const
  {
    TrilinosWrappers::types::int_type begin, end;
#      ifndef DEAL_II_WITH_64BIT_INDICES
    begin = matrix->RowMap().MinMyGID();
    end   = matrix->RowMap().MaxMyGID() + 1;
#      else
    begin = matrix->RowMap().MinMyGID64();
    end   = matrix->RowMap().MaxMyGID64() + 1;
#      endif

    return ((index >= static_cast<size_type>(begin)) &&
            (index < static_cast<size_type>(end)));
  }



  inline bool
  SparseMatrix::is_compressed() const
  {
    return compressed;
  }



  // Inline the set() and add() functions, since they will be called
  // frequently, and the compiler can optimize away some unnecessary loops
  // when the sizes are given at compile time.
  template <>
  void
  SparseMatrix::set<TrilinosScalar>(const size_type       row,
                                    const size_type       n_cols,
                                    const size_type *     col_indices,
                                    const TrilinosScalar *values,
                                    const bool            elide_zero_values);



  template <typename Number>
  void
  SparseMatrix::set(const size_type  row,
                    const size_type  n_cols,
                    const size_type *col_indices,
                    const Number *   values,
                    const bool       elide_zero_values)
  {
    std::vector<TrilinosScalar> trilinos_values(n_cols);
    std::copy(values, values + n_cols, trilinos_values.begin());
    this->set(
      row, n_cols, col_indices, trilinos_values.data(), elide_zero_values);
  }



  inline void
  SparseMatrix::set(const size_type      i,
                    const size_type      j,
                    const TrilinosScalar value)
  {
    AssertIsFinite(value);

    set(i, 1, &j, &value, false);
  }



  inline void
  SparseMatrix::set(const std::vector<size_type> &    indices,
                    const FullMatrix<TrilinosScalar> &values,
                    const bool                        elide_zero_values)
  {
    Assert(indices.size() == values.m(),
           ExcDimensionMismatch(indices.size(), values.m()));
    Assert(values.m() == values.n(), ExcNotQuadratic());

    for (size_type i = 0; i < indices.size(); ++i)
      set(indices[i],
          indices.size(),
          indices.data(),
          &values(i, 0),
          elide_zero_values);
  }



  inline void
  SparseMatrix::add(const size_type      i,
                    const size_type      j,
                    const TrilinosScalar value)
  {
    AssertIsFinite(value);

    if (value == 0)
      {
        // we have to check after Insert/Add in any case to be consistent
        // with the MPI communication model, but we can save some
        // work if the addend is zero. However, these actions are done in case
        // we pass on to the other function.

        // TODO: fix this (do not run compress here, but fail)
        if (last_action == Insert)
          {
            int ierr;
            ierr = matrix->GlobalAssemble(*column_space_map,
                                          matrix->RowMap(),
                                          false);

            Assert(ierr == 0, ExcTrilinosError(ierr));
            (void)ierr; // removes -Wunused-but-set-variable in optimized mode
          }

        last_action = Add;

        return;
      }
    else
      add(i, 1, &j, &value, false);
  }



  // inline "simple" functions that are called frequently and do only involve
  // a call to some Trilinos function.
  inline SparseMatrix::size_type
  SparseMatrix::m() const
  {
#      ifndef DEAL_II_WITH_64BIT_INDICES
    return matrix->NumGlobalRows();
#      else
    return matrix->NumGlobalRows64();
#      endif
  }



  inline SparseMatrix::size_type
  SparseMatrix::n() const
  {
    // If the matrix structure has not been fixed (i.e., we did not have a
    // sparsity pattern), it does not know about the number of columns so we
    // must always take this from the additional column space map
    Assert(column_space_map.get() != nullptr, ExcInternalError());
    return n_global_elements(*column_space_map);
  }



  inline unsigned int
  SparseMatrix::local_size() const
  {
    return matrix->NumMyRows();
  }



  inline std::pair<SparseMatrix::size_type, SparseMatrix::size_type>
  SparseMatrix::local_range() const
  {
    size_type begin, end;
#      ifndef DEAL_II_WITH_64BIT_INDICES
    begin = matrix->RowMap().MinMyGID();
    end   = matrix->RowMap().MaxMyGID() + 1;
#      else
    begin = matrix->RowMap().MinMyGID64();
    end   = matrix->RowMap().MaxMyGID64() + 1;
#      endif

    return std::make_pair(begin, end);
  }



  inline SparseMatrix::size_type
  SparseMatrix::n_nonzero_elements() const
  {
#      ifndef DEAL_II_WITH_64BIT_INDICES
    return matrix->NumGlobalNonzeros();
#      else
    return matrix->NumGlobalNonzeros64();
#      endif
  }



  template <typename SparsityPatternType>
  inline void
  SparseMatrix::reinit(const IndexSet &           parallel_partitioning,
                       const SparsityPatternType &sparsity_pattern,
                       const MPI_Comm &           communicator,
                       const bool                 exchange_data)
  {
    reinit(parallel_partitioning,
           parallel_partitioning,
           sparsity_pattern,
           communicator,
           exchange_data);
  }



  template <typename number>
  inline void
  SparseMatrix::reinit(const IndexSet &parallel_partitioning,
                       const ::dealii::SparseMatrix<number> &sparse_matrix,
                       const MPI_Comm &                      communicator,
                       const double                          drop_tolerance,
                       const bool                            copy_values,
                       const ::dealii::SparsityPattern *     use_this_sparsity)
  {
    Epetra_Map map =
      parallel_partitioning.make_trilinos_map(communicator, false);
    reinit(parallel_partitioning,
           parallel_partitioning,
           sparse_matrix,
           drop_tolerance,
           copy_values,
           use_this_sparsity);
  }



  inline const Epetra_CrsMatrix &
  SparseMatrix::trilinos_matrix() const
  {
    return static_cast<const Epetra_CrsMatrix &>(*matrix);
  }



  inline const Epetra_CrsGraph &
  SparseMatrix::trilinos_sparsity_pattern() const
  {
    return matrix->Graph();
  }



  inline IndexSet
  SparseMatrix::locally_owned_domain_indices() const
  {
    return IndexSet(matrix->DomainMap());
  }



  inline IndexSet
  SparseMatrix::locally_owned_range_indices() const
  {
    return IndexSet(matrix->RangeMap());
  }



  inline void
  SparseMatrix::prepare_add()
  {
    // nothing to do here
  }



  inline void
  SparseMatrix::prepare_set()
  {
    // nothing to do here
  }


  namespace internal
  {
    namespace LinearOperatorImplementation
    {
      template <typename Solver, typename Preconditioner>
      typename std::enable_if<
        std::is_base_of<TrilinosWrappers::SolverBase, Solver>::value &&
          std::is_base_of<TrilinosWrappers::PreconditionBase,
                          Preconditioner>::value,
        TrilinosPayload>::type
      TrilinosPayload::inverse_payload(
        Solver &              solver,
        const Preconditioner &preconditioner) const
      {
        const auto &payload = *this;

        TrilinosPayload return_op(payload);

        // Capture by copy so the payloads are always valid

        return_op.inv_vmult = [payload, &solver, &preconditioner](
                                TrilinosPayload::Domain &     tril_dst,
                                const TrilinosPayload::Range &tril_src) {
          // Duplicated from TrilinosWrappers::PreconditionBase::vmult
          // as well as from TrilinosWrappers::SparseMatrix::Tvmult
          Assert(&tril_src != &tril_dst,
                 TrilinosWrappers::SparseMatrix::ExcSourceEqualsDestination());
          internal::check_vector_map_equality(payload,
                                              tril_src,
                                              tril_dst,
                                              !payload.UseTranspose());
          solver.solve(payload, tril_dst, tril_src, preconditioner);
        };

        return_op.inv_Tvmult = [payload, &solver, &preconditioner](
                                 TrilinosPayload::Range &       tril_dst,
                                 const TrilinosPayload::Domain &tril_src) {
          // Duplicated from TrilinosWrappers::PreconditionBase::vmult
          // as well as from TrilinosWrappers::SparseMatrix::Tvmult
          Assert(&tril_src != &tril_dst,
                 TrilinosWrappers::SparseMatrix::ExcSourceEqualsDestination());
          internal::check_vector_map_equality(payload,
                                              tril_src,
                                              tril_dst,
                                              payload.UseTranspose());

          const_cast<TrilinosPayload &>(payload).transpose();
          solver.solve(payload, tril_dst, tril_src, preconditioner);
          const_cast<TrilinosPayload &>(payload).transpose();
        };

        // If the input operator is already setup for transpose operations, then
        // we must do similar with its inverse.
        if (return_op.UseTranspose() == true)
          std::swap(return_op.inv_vmult, return_op.inv_Tvmult);

        return return_op;
      }

      template <typename Solver, typename Preconditioner>
      typename std::enable_if<
        !(std::is_base_of<TrilinosWrappers::SolverBase, Solver>::value &&
          std::is_base_of<TrilinosWrappers::PreconditionBase,
                          Preconditioner>::value),
        TrilinosPayload>::type
      TrilinosPayload::inverse_payload(Solver &, const Preconditioner &) const
      {
        TrilinosPayload return_op(*this);

        return_op.inv_vmult = [](TrilinosPayload::Domain &,
                                 const TrilinosPayload::Range &) {
          AssertThrow(false,
                      ExcMessage("Payload inv_vmult disabled because of "
                                 "incompatible solver/preconditioner choice."));
        };

        return_op.inv_Tvmult = [](TrilinosPayload::Range &,
                                  const TrilinosPayload::Domain &) {
          AssertThrow(false,
                      ExcMessage("Payload inv_vmult disabled because of "
                                 "incompatible solver/preconditioner choice."));
        };

        return return_op;
      }
    } // namespace LinearOperatorImplementation
  }   // namespace internal

  template <>
  void
  SparseMatrix::set<TrilinosScalar>(const size_type       row,
                                    const size_type       n_cols,
                                    const size_type *     col_indices,
                                    const TrilinosScalar *values,
                                    const bool            elide_zero_values);
#    endif // DOXYGEN

}  /* namespace TrilinosWrappers */ 


DEAL_II_NAMESPACE_CLOSE


#  endif // DEAL_II_WITH_TRILINOS


 /*-----------------------   trilinos_sparse_matrix.h     --------------------*/ 

#endif
 /*-----------------------   trilinos_sparse_matrix.h     --------------------*/ 


