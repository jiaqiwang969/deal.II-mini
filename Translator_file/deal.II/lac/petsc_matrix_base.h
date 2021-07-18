//include/deal.II-translator/lac/petsc_matrix_base_0.txt
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

#ifndef dealii_petsc_matrix_base_h
#  define dealii_petsc_matrix_base_h


#  include <deal.II/base/config.h>

#  ifdef DEAL_II_WITH_PETSC

#    include <deal.II/base/subscriptor.h>

#    include <deal.II/lac/exceptions.h>
#    include <deal.II/lac/full_matrix.h>
#    include <deal.II/lac/petsc_compatibility.h>
#    include <deal.II/lac/petsc_vector_base.h>
#    include <deal.II/lac/vector_operation.h>

#    include <petscmat.h>

#    include <cmath>
#    include <memory>
#    include <vector>

DEAL_II_NAMESPACE_OPEN

// Forward declarations
#    ifndef DOXYGEN
template <typename Matrix>
class BlockMatrixBase;
#    endif


namespace PETScWrappers
{
  // forward declarations
  class MatrixBase;

  namespace MatrixIterators
  {
    /**
     * 这个类作为一个迭代器，在PETSc矩阵的元素上行走。由于PETSc为所有类型的矩阵提供了一个统一的接口，这个迭代器可以用来访问稀疏和完整的矩阵。
     * 请注意，PETSc并不保证每一行中元素的顺序。还要注意的是，访问全矩阵的元素，竟然只显示矩阵的非零元素，而不是所有元素。
     * @ingroup PETScWrappers
     *
     */
    class const_iterator
    {
    private:
      /**
       * 迭代器的访问器类
       *
       */
      class Accessor
      {
      public:
        /**
         * 声明容器大小的类型。
         *
         */
        using size_type = types::global_dof_index;

        /**
         * 构造器。因为我们只使用访问器进行读取访问，所以一个常量矩阵指针就足够了。
         *
         */
        Accessor(const MatrixBase *matrix,
                 const size_type   row,
                 const size_type   index);

        /**
         * 这个对象所代表的元素的行号。
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

        /**
         * 这个矩阵条目的值。
         *
         */
        PetscScalar
        value() const;

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
                       int,
                       int,
                       int,
                       << "You tried to access row " << arg1
                       << " of a distributed matrix, but only rows " << arg2
                       << " through " << arg3
                       << " are stored locally and can be accessed.");

      private:
        /**
         * 访问的矩阵。
         *
         */
        mutable MatrixBase *matrix;

        /**
         * 当前行数。
         *
         */
        size_type a_row;

        /**
         * 当前行的索引。
         *
         */
        size_type a_index;

        /**
         * 缓存，我们在这里存储当前行的列索引。这是必要的，因为PETSc对其矩阵元素的访问是相当困难的，当我们进入某一行时，一次性复制该行的所有列项，要比反复向PETSc索取单个列项更有效率。这也有一定道理，因为我们很可能会按顺序访问它们。
         * 为了使迭代器/存取器的复制具有可接受的性能，我们为这些条目保留一个共享指针，以便在必要时有多个存取器可以访问这些数据。
         *
         */
        std::shared_ptr<const std::vector<size_type>> colnum_cache;

        /**
         * 这个行的值的类似缓存。
         *
         */
        std::shared_ptr<const std::vector<PetscScalar>> value_cache;

        /**
         * 丢弃旧的行缓存（它们可能仍然被其他访问器使用），并为这个访问器目前所指向的行生成新的行缓存。
         *
         */
        void
        visit_present_row();

        // Make enclosing class a friend.
        friend class const_iterator;
      };

    public:
      /**
       * 声明容器大小的类型。
       *
       */
      using size_type = types::global_dof_index;

      /**
       * 构造函数。在矩阵 @p matrix
       * 中为给定的行和其中的索引创建一个迭代器。
       *
       */
      const_iterator(const MatrixBase *matrix,
                     const size_type   row,
                     const size_type   index);

      /**
       * 前缀增量。
       *
       */
      const_iterator &
      operator++();

      /**
       * 后缀增量。
       *
       */
      const_iterator
      operator++(int);

      /**
       * 撤消运算符。
       *
       */
      const Accessor &operator*() const;

      /**
       * 解除引用操作符。
       *
       */
      const Accessor *operator->() const;

      /**
       * 比较。真，如果两个迭代器都指向同一个矩阵位置。
       *
       */
      bool
      operator==(const const_iterator &) const;
      /**
       * <tt>==</tt>的倒数。
       *
       */
      bool
      operator!=(const const_iterator &) const;

      /**
       * 比较运算符。如果第一行数字较小，或者行数字相等且第一个索引较小，则结果为真。
       *
       */
      bool
      operator<(const const_iterator &) const;

      /**
       * 异常情况
       *
       */
      DeclException2(ExcInvalidIndexWithinRow,
                     int,
                     int,
                     << "Attempt to access element " << arg2 << " of row "
                     << arg1 << " which doesn't have that many elements.");

    private:
      /**
       * 存储一个访问器类的对象。
       *
       */
      Accessor accessor;
    };

  } // namespace MatrixIterators


  /**
   * 所有在PETSc矩阵类型之上实现的矩阵类的基类。由于在PETSc中，所有的矩阵类型（即顺序和平行，稀疏，阻塞等）都是通过填充一个抽象对象的内容来建立的，而这个抽象对象只能通过一个独立于实际矩阵类型的指针来引用，所以我们可以在这个基类中实现几乎所有的矩阵功能。然后，派生类将只需要提供创建一种或另一种矩阵的功能。
   * 这个类的接口是以deal.II中现有的SparseMatrix类为模型的。它有几乎相同的成员函数，而且通常是可以交换的。然而，由于PETSc只支持单一的标量类型（要么是双数、浮点数，要么是复杂的数据类型），所以它没有模板化，只能与你的PETSc安装中定义的数据类型PetscScalar一起工作。
   * 请注意，只有在矩阵装配后调用了函数 @p MatAssemblyBegin
   * 和 @p MatAssemblyEnd
   * 的情况下，PETSc才能保证操作符合你的期望。因此，你需要在实际使用矩阵之前调用
   * SparseMatrix::compress()  。这也会调用 @p MatCompress
   * ，通过丢弃未使用的元素来压缩稀疏矩阵的存储格式。PETSc允许在调用这些函数后继续装配矩阵，但由于此后不再有可用的条目，所以最好在装配阶段结束后，在主动使用矩阵之前，只调用一次
   * SparseMatrix::compress() 。
   * @ingroup PETScWrappers
   * @ingroup Matrix1
   *
   */
  class MatrixBase : public Subscriptor
  {
  public:
    /**
     * 为迭代器类声明一个别名。
     *
     */
    using const_iterator = MatrixIterators::const_iterator;

    /**
     * 声明容器大小的类型。
     *
     */
    using size_type = types::global_dof_index;

    /**
     * 声明一个类似于所有其他容器类的别名。
     *
     */
    using value_type = PetscScalar;

    /**
     * 默认的构造函数。
     *
     */
    MatrixBase();

    /**
     * 复制构造函数。它被删除了，因为在不知道所存储的矩阵的具体种类的情况下复制这个基类，可能会错过重要的细节，而且如果矩阵很大的话，也会很昂贵。
     *
     */
    MatrixBase(const MatrixBase &) = delete;

    /**
     * 复制操作符。它被删除了，因为在不知道所存储的矩阵的具体种类的情况下复制这个基类可能会错过重要的细节，而且如果矩阵很大的话，费用也很高。
     *
     */
    MatrixBase &
    operator=(const MatrixBase &) = delete;

    /**
     * 销毁器。虚化，以便人们可以使用指向这个类的指针。
     *
     */
    virtual ~MatrixBase() override;

    /**
     * 这个操作符将一个标量分配给一个矩阵。因为这通常没有什么意义（我们应该把所有的矩阵条目都设置为这个值吗？仅仅是稀疏模式的非零条目？），这个操作只允许在实际要分配的值为零时进行。这个操作符的存在只是为了允许明显的符号<tt>matrix=0</tt>，它将矩阵的所有元素设置为零，但保留了之前使用的稀疏模式。
     *
     */
    MatrixBase &
    operator=(const value_type d);
    /**
     * 释放所有内存并返回到与调用默认构造函数后相同的状态。
     *
     */
    void
    clear();

    /**
     * 将元素(<i>i,j</i>)设置为 @p value.
     * 如果现在的对象(来自这个对象的派生类)恰好是一个稀疏矩阵，那么这个函数就会向矩阵添加一个新的条目，如果该条目之前不存在的话，这与SparseMatrix类形成了很大的反差，如果该条目不存在则会抛出错误。如果<tt>value</tt>不是一个有限的数字，就会抛出一个异常。
     *
     */
    void
    set(const size_type i, const size_type j, const PetscScalar value);

    /**
     * 将FullMatrix<double>中给出的所有元素设置为<tt>indices</tt>给出的稀疏矩阵位置。换句话说，这个函数将<tt>full_matrix</tt>中的元素写入调用的矩阵中，对矩阵的行和列都使用<tt>indices</tt>指定的本地到全球索引。这个函数假设一个二次稀疏矩阵和一个二次全矩阵，这是FE计算中的通常情况。
     * 如果现在的对象（来自这个对象的派生类）恰好是一个稀疏矩阵，那么这个函数就会向矩阵添加一些新的条目，如果这些条目之前不存在的话，这与SparseMatrix类非常不同，后者在条目不存在的时候会抛出一个错误。
     * 可选的参数<tt>elide_zero_values</tt>可以用来指定零值是否应该被插入，还是应该被过滤掉。默认值是<tt>false</tt>，也就是说，即使是零值也要插入/替换。
     *
     */
    void
    set(const std::vector<size_type> & indices,
        const FullMatrix<PetscScalar> &full_matrix,
        const bool                     elide_zero_values = false);

    /**
     * 与之前的功能相同，但现在包括了使用矩形full_matrices的可能性，以及在行和列上分别使用不同的局部到全局的索引。
     *
     */
    void
    set(const std::vector<size_type> & row_indices,
        const std::vector<size_type> & col_indices,
        const FullMatrix<PetscScalar> &full_matrix,
        const bool                     elide_zero_values = false);

    /**
     * 将矩阵的指定行中的几个元素与<tt>col_indices</tt>给出的列索引设置为相应的值。
     * 如果现在的对象（来自这个对象的派生类）恰好是一个稀疏矩阵，那么这个函数就会向矩阵添加一些新的条目，如果这些条目之前不存在的话，这与SparseMatrix类非常不同，后者在条目不存在的时候会抛出一个错误。
     * 可选的参数<tt>elide_zero_values</tt>可以用来指定零值是否应该被插入，还是应该被过滤掉。默认值是<tt>false</tt>，也就是说，即使是零值也要插入/替换。
     *
     */
    void
    set(const size_type                 row,
        const std::vector<size_type> &  col_indices,
        const std::vector<PetscScalar> &values,
        const bool                      elide_zero_values = false);

    /**
     * 将几个元素设置为由<tt>values</tt>给定的值，在由col_indices给定的列中设置为稀疏矩阵的行。
     * 如果现在的对象（来自这个对象的派生类）恰好是一个稀疏矩阵，那么这个函数就会向矩阵添加一些新的条目，如果这些条目之前不存在的话，这与SparseMatrix类非常不同，后者在条目不存在的时候会抛出一个错误。
     * 可选的参数<tt>elide_zero_values</tt>可以用来指定零值是否应该被插入，还是应该被过滤掉。默认值是<tt>false</tt>，也就是说，即使是零值也要插入/替换。
     *
     */
    void
    set(const size_type    row,
        const size_type    n_cols,
        const size_type *  col_indices,
        const PetscScalar *values,
        const bool         elide_zero_values = false);

    /**
     * 将 @p value 添加到元素（<i>i,j</i>）。
     * 如果现在的对象（来自这个对象的派生类）恰好是一个稀疏矩阵，那么这个函数就会向矩阵添加一个新的条目，如果该条目之前不存在的话，这与SparseMatrix类非常不同，后者在该条目不存在的情况下会抛出一个错误。如果<tt>value</tt>不是一个有限的数字，就会抛出一个异常。
     *
     */
    void
    add(const size_type i, const size_type j, const PetscScalar value);

    /**
     * 将FullMatrix<double>中给出的所有元素添加到由<tt>indices</tt>给出的稀疏矩阵位置。换句话说，这个函数将<tt>full_matrix</tt>中的元素添加到调用矩阵的相应条目中，使用<tt>indices</tt>为矩阵的行和列指定的本地到全球索引。这个函数假设一个二次稀疏矩阵和一个二次全矩阵，这是FE计算中的通常情况。
     * 如果现在的对象（来自这个对象的派生类）恰好是一个稀疏矩阵，那么这个函数就会向矩阵添加一些新的条目，如果这些条目之前不存在的话，这与SparseMatrix类非常不同，后者在条目不存在的时候会抛出一个错误。
     * 可选参数<tt>elide_zero_values</tt>可以用来指定是无论如何都要添加零值，还是要过滤掉这些零值，只添加非零数据。默认值是<tt>true</tt>，也就是说，零值不会被添加到矩阵中。
     *
     */
    void
    add(const std::vector<size_type> & indices,
        const FullMatrix<PetscScalar> &full_matrix,
        const bool                     elide_zero_values = true);

    /**
     * 与之前的函数相同，但现在包括了使用矩形full_matrices的可能性，以及在行和列上分别使用不同的本地到全球索引。
     *
     */
    void
    add(const std::vector<size_type> & row_indices,
        const std::vector<size_type> & col_indices,
        const FullMatrix<PetscScalar> &full_matrix,
        const bool                     elide_zero_values = true);

    /**
     * 将矩阵的指定行中的几个元素与<tt>col_indices</tt>给出的列索引设置为相应的值。
     * 如果现在的对象（来自这个对象的派生类）恰好是一个稀疏矩阵，那么这个函数就会向矩阵添加一些新的条目，如果这些条目之前不存在的话，这与SparseMatrix类非常不同，后者在条目不存在的时候会抛出一个错误。
     * 可选参数<tt>elide_zero_values</tt>可以用来指定是无论如何都要添加零值，还是要过滤掉这些零值，只添加非零数据。默认值是<tt>true</tt>，也就是说，零值不会被添加到矩阵中。
     *
     */
    void
    add(const size_type                 row,
        const std::vector<size_type> &  col_indices,
        const std::vector<PetscScalar> &values,
        const bool                      elide_zero_values = true);

    /**
     * 在给定的全局矩阵行中，在稀疏矩阵中由col_indices指定的列中添加一个由<tt>values</tt>给出的数值阵列。
     * 如果现在的对象（来自这个对象的派生类）恰好是一个稀疏矩阵，那么这个函数就会向矩阵添加一些新的条目，如果这些条目之前不存在的话，这与SparseMatrix类形成鲜明对比，后者在条目不存在的时候会抛出一个错误。
     * 可选参数<tt>elide_zero_values</tt>可以用来指定是无论如何都要添加零值，还是要过滤掉这些零值，只添加非零数据。默认值是<tt>true</tt>，也就是说，零值不会被添加到矩阵中。
     *
     */
    void
    add(const size_type    row,
        const size_type    n_cols,
        const size_type *  col_indices,
        const PetscScalar *values,
        const bool         elide_zero_values      = true,
        const bool         col_indices_are_sorted = false);

    /**
     * 将此<tt>行</tt>中的所有元素设置为零，将其删除。这个函数并不修改分配的非零条目的数量，它只是将一些条目设置为零。不过，它可能会将它们从稀疏模式中删除（但会保留分配的内存，以备以后再次添加新的条目）。
     * 这个操作用于消除约束（例如由于挂起的节点），并确保我们可以将这个修改写入矩阵，而不需要从矩阵中读取条目（例如非零元素的位置）。
     *
     * - 如果没有这个操作，消除平行矩阵的约束是一个相当复杂的过程。        第二个参数可以用来将该行的对角线条目设置为一个不同于零的值。默认是将其设置为零。
     *
     */
    void
    clear_row(const size_type row, const PetscScalar new_diag_value = 0);

    /**
     * 与clear_row()相同，只是它同时作用于若干行。
     * 第二个参数可以用来将所有被清除的行的对角线条目设置为不同于0的内容。请注意，所有这些对角线项都得到相同的值
     *
     * - 如果你想要不同的对角线条目的值，你必须手动设置它们。
     *
     */
    void
    clear_rows(const std::vector<size_type> &rows,
               const PetscScalar             new_diag_value = 0);

    /**
     * PETSc矩阵存储了它们自己的稀疏性模式。因此，与我们自己的SparsityPattern类相类似，这个函数压缩了稀疏模式，并允许将得到的矩阵用于所有其他操作，而以前只允许使用汇编函数。因此，一旦你组装了矩阵，就必须调用这个函数。        更多信息请参见  @ref GlossCompress  "压缩分布式对象"
     * 。
     *
     */
    void
    compress(const VectorOperation::values operation);

    /**
     * 返回条目的值（<i>i,j</i>）。
     * 这可能是一个昂贵的操作，你应该始终注意在哪里调用这个函数。
     * 与 @p MatrixBase
     * 类中的相应函数相比，如果相应的条目不存在于该类的稀疏模式中，我们不会抛出一个异常，因为PETSc并不传输这一信息。
     * 因此这个函数完全等同于<tt>el()</tt>函数。
     *
     */
    PetscScalar
    operator()(const size_type i, const size_type j) const;

    /**
     * 返回矩阵条目的值（<i>i,j</i>）。如果这个条目不存在于稀疏模式中，那么就返回0。虽然这在某些情况下可能很方便，但请注意，由于没有使用矩阵的稀疏性，写出的算法与最优解相比很简单，很慢。
     *
     */
    PetscScalar
    el(const size_type i, const size_type j) const;

    /**
     * 返回第<i>i</i>行中的主对角线元素。如果矩阵不是二次的，这个函数会抛出一个错误。
     * 由于我们不能直接访问底层数据结构，这个函数并不比使用el()函数的元素访问快。然而，我们提供这个函数是为了与SparseMatrix类兼容。
     *
     */
    PetscScalar
    diag_element(const size_type i) const;

    /**
     * 返回这个矩阵的行数。
     *
     */
    size_type
    m() const;

    /**
     * 返回这个矩阵中的列数。
     *
     */
    size_type
    n() const;

    /**
     * 返回矩阵的本地维度，即存储在当前MPI进程中的行数。对于顺序矩阵，这个数字与m()相同，但对于并行矩阵，这个数字可能更小。
     * 要想知道到底哪些元素被存储在本地，可以使用local_range()。
     *
     */
    size_type
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
     * 返回对与该矩阵一起使用的MPI通信器对象的引用。这个函数必须在派生类中实现。
     *
     */
    virtual const MPI_Comm &
    get_mpi_communicator() const = 0;

    /**
     * 返回这个矩阵的非零元素的数量。实际上，它返回的是稀疏模式中的条目数；如果任何一个条目碰巧是零，无论如何都会被计算在内。
     *
     */
    size_type
    n_nonzero_elements() const;

    /**
     * 特定行中的条目数。
     *
     */
    size_type
    row_length(const size_type row) const;

    /**
     * 返回矩阵的l1准则，即 $|M|_1=max_{all columns j}\sum_{all rows
     * i} |M_ij|$
     * ，（最大列数之和）。这是一个自然的矩阵准则，与向量的l1准则兼容，即
     * $|Mv|_1\leq |M|_1 |v|_1$  。(参看Haemmerlin-Hoffmann: Numerische
     * Mathematik)
     *
     */
    PetscReal
    l1_norm() const;

    /**
     * 返回矩阵的linfty-norm，即 $|M|_infty=max_{all rows i}\sum_{all
     * columns j} |M_ij|$
     * ，（最大行数之和）。这是一个自然的矩阵规范，与向量的linfty-norm兼容，即
     * $|Mv|_infty \leq |M|_infty |v|_infty$  。(参看Haemmerlin-Hoffmann:
     * Numerische Mathematik)
     *
     */
    PetscReal
    linfty_norm() const;

    /**
     * 返回矩阵的frobenius
     * norm，即矩阵中所有条目的平方之和的平方根。
     *
     */
    PetscReal
    frobenius_norm() const;


    /**
     * 返回向量 $v$ 相对于该矩阵诱导的准则的平方，即
     * $\left(v,Mv\right)$
     * 。这很有用，例如在有限元背景下，一个函数的 $L_2$
     * 规范等于相对于代表有限元函数节点值的向量的质量矩阵的矩阵规范。
     * 很明显，对于这个操作，矩阵需要是二次的。
     * 这个函数的实现没有deal.II中使用的 @p MatrixBase
     * 类（即原始的，而不是PETSc封装类）的效率高，因为PETSc不支持这个操作，需要一个临时向量。
     * 注意，如果当前对象代表一个并行的分布式矩阵（类型为
     * PETScWrappers::MPI::SparseMatrix),
     * ，那么给出的向量也必须是一个分布式向量。反之，如果矩阵不是分布式的，那么向量也不可能是。
     *
     */
    PetscScalar
    matrix_norm_square(const VectorBase &v) const;


    /**
     * 计算矩阵标量乘积  $\left(u,Mv\right)$  。
     * 这个函数的实现不如deal.II中使用的 @p MatrixBase
     * 类（即原始函数，而不是PETSc封装类）的效率高，因为PETSc不支持这个操作，需要一个临时矢量。
     * 注意，如果当前对象代表一个平行分布式矩阵（类型为
     * PETScWrappers::MPI::SparseMatrix),
     * ，那么两个向量也必须是分布式向量。反之，如果矩阵不是分布式的，那么两个向量都不可能是。
     *
     */
    PetscScalar
    matrix_scalar_product(const VectorBase &u, const VectorBase &v) const;

    /**
     * 返回矩阵的轨迹，即矩阵中所有对角线项的总和。
     *
     */
    PetscScalar
    trace() const;

    /**
     * 将整个矩阵乘以一个固定系数。
     *
     */
    MatrixBase &
    operator*=(const PetscScalar factor);

    /**
     * 用整个矩阵除以一个固定系数。
     *
     */
    MatrixBase &
    operator/=(const PetscScalar factor);


    /**
     * 将矩阵 @p other 按系数 @p factor
     * 的比例添加到当前矩阵中。
     *
     */
    MatrixBase &
    add(const PetscScalar factor, const MatrixBase &other);

    /**
     * 矩阵-向量乘法：让<i>dst =
     * M*src</i>与<i>M</i>是这个矩阵。
     * 源和目的不能是同一个向量。
     * 注意，如果当前对象代表一个平行分布式矩阵（类型为
     * PETScWrappers::MPI::SparseMatrix),
     * ，那么两个向量也必须是分布式向量。反之，如果矩阵不是分布式的，那么两个向量都不可以是。
     *
     */
    void
    vmult(VectorBase &dst, const VectorBase &src) const;

    /**
     * 矩阵-向量乘法：让<i>dst =
     * M<sup>T</sup>*src</i>与<i>M</i>为这个矩阵。这个函数与vmult()的作用相同，但需要转置的矩阵。
     * 源和目的不能是同一个向量。
     * 注意，如果当前对象代表一个并行分布式矩阵（类型为
     * PETScWrappers::MPI::SparseMatrix),
     * ，那么两个向量也必须是分布式向量。反之，如果矩阵不是分布式的，那么两个向量都不可以是。
     *
     */
    void
    Tvmult(VectorBase &dst, const VectorBase &src) const;

    /**
     * 加法
     * 矩阵-向量乘法。在<i>dst</i>上添加<i>M*src</i>，<i>M</i>为该矩阵。
     * 源和目的不能是同一个向量。
     * 注意，如果当前对象代表一个平行分布式矩阵（类型为
     * PETScWrappers::MPI::SparseMatrix),
     * ，那么两个向量也必须是分布式向量。反之，如果矩阵不是分布式的，那么两个向量都不可以是。
     *
     */
    void
    vmult_add(VectorBase &dst, const VectorBase &src) const;

    /**
     * 加法
     * 矩阵-向量乘法。将<i>M<sup>T</sup>*src</i>加到<i>dst</i>，<i>M</i>是这个矩阵。这个函数与vmult_add()的作用相同，但需要转置的矩阵。
     * 来源和目的地不能是同一个向量。
     * 注意，如果当前对象代表一个并行分布式矩阵（类型为
     * PETScWrappers::MPI::SparseMatrix),
     * ，那么两个向量也必须是分布式向量。反之，如果矩阵不是分布式的，那么两个向量都不可以是。
     *
     */
    void
    Tvmult_add(VectorBase &dst, const VectorBase &src) const;

    /**
     * 计算一个方程<i>Mx=b</i>的残差，其中残差被定义为<i>r=b-Mx</i>。将残差写入
     * @p dst.  返回残差向量的<i>l<sub>2</sub></i>准则。
     * 源<i>x</i>和目的<i>dst</i>不能是同一个向量。
     * 注意，如果当前对象代表一个平行分布式矩阵（类型为
     * PETScWrappers::MPI::SparseMatrix),
     * ，那么所有的向量也必须是分布式向量。反之，如果矩阵不是分布式的，那么两个向量都不能是。
     *
     */
    PetscScalar
    residual(VectorBase &dst, const VectorBase &x, const VectorBase &b) const;

    /**
     * 迭代器从第一个条目开始。这只能在拥有整个矩阵的处理器上调用。在所有其他情况下，请参考以行号为参数的begin()的版本。
     *
     */
    const_iterator
    begin() const;

    /**
     * 最后的迭代器。这只能在拥有整个矩阵的处理器上调用。在所有其他情况下，请参考end()的版本，以一个行号作为参数。
     *
     */
    const_iterator
    end() const;

    /**
     * 迭代器从第 @p r. 行的第一个条目开始
     * 注意，如果给定的行是空的，即不包含任何非零条目，那么这个函数返回的迭代器就等于<tt>end(r)</tt>。还要注意的是，在这种情况下，迭代器可能不能被解除引用。
     *
     */
    const_iterator
    begin(const size_type r) const;

    /**
     * 行<tt>r</tt>的最终迭代器。它指向超过行 @p r,
     * 末尾的第一个元素，或者超过整个稀疏模式的末尾。
     * 请注意，结束迭代器不一定是可被解除引用的。特别是如果它是一个矩阵的最后一行的结束迭代器，情况更是如此。
     *
     */
    const_iterator
    end(const size_type r) const;

    /**
     * 转换操作符，以获得对底层PETSc类型的访问。如果你这样做，你就切断了这个类可能需要的一些信息，所以这个转换操作符应该只在你知道你要做什么的时候使用。特别是，它应该只用于对矩阵的只读操作。
     *
     */
    operator Mat() const;

    /**
     * 返回一个对底层PETSc类型的引用。它可以用来修改底层数据，所以只有在你知道你在做什么的时候才使用它。
     *
     */
    Mat &
    petsc_matrix();

    /**
     * 对一个矩阵进行原地转置。
     *
     */
    void
    transpose();

    /**
     * 测试矩阵是否是对称的。 默认公差为 $1000\times32$
     * -位机器精度。
     *
     */
    PetscBool
    is_symmetric(const double tolerance = 1.e-12);

    /**
     * 测试一个矩阵是否是赫米特的，即它是其转置的复共轭。默认公差为
     * $1000\times32$  -位机器精度。
     *
     */
    PetscBool
    is_hermitian(const double tolerance = 1.e-12);

    /**
     * 使用PETSc内部矩阵查看器功能<tt>MatView</tt>打印PETSc矩阵对象的值。默认格式是打印非零矩阵元素。对于其他有效的查看格式，请参考http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Mat/MatView.html
     *
     */
    void
    write_ascii(const PetscViewerFormat format = PETSC_VIEWER_DEFAULT);

    /**
     * 打印矩阵的元素到给定的输出流。          @param[in,out]
     * out 要写入的输出流。      @param[in]  alternative_output
     * 这个参数被忽略。它的存在是为了与其他矩阵类中的类似函数兼容。
     *
     */
    void
    print(std::ostream &out, const bool alternative_output = false) const;

    /**
     * 返回该矩阵在该CPU上所消耗的字节数。
     *
     */
    std::size_t
    memory_consumption() const;

    /**
     * 异常情况
     *
     */
    DeclExceptionMsg(ExcSourceEqualsDestination,
                     "You are attempting an operation on two matrices that "
                     "are the same object, but the operation requires that the "
                     "two objects are in fact different.");

    /**
     * 异常情况。
     *
     */
    DeclException2(ExcWrongMode,
                   int,
                   int,
                   << "You tried to do a "
                   << (arg1 == 1 ? "'set'" : (arg1 == 2 ? "'add'" : "???"))
                   << " operation but the matrix is currently in "
                   << (arg2 == 1 ? "'set'" : (arg2 == 2 ? "'add'" : "???"))
                   << " mode. You first have to call 'compress()'.");

  protected:
    /**
     * 一个PETSc中的通用矩阵对象。实际的类型是稀疏矩阵，在构造函数中设置。
     *
     */
    Mat matrix;

    /**
     * 存储最后一个动作是写还是加操作。
     *
     */
    VectorOperation::values last_action;

    /**
     * 确保此调用后的动作所需的添加/设置模式与当前模式兼容。应该从所有访问矩阵元素的内部函数中调用。
     *
     */
    void
    prepare_action(const VectorOperation::values new_action);

    /**
     * 内部函数，检查是否有未决的插入/添加操作。否则会抛出一个异常。在调用任何修改矩阵的PETSc内部函数之前，都是有用的。
     *
     */
    void
    assert_is_compressed();

    /**
     * 对于某些矩阵存储格式，特别是PETSc分布式块状矩阵，单个元素的设置和添加操作不能自由混合。相反，当我们想从设置元素切换到添加元素时，我们必须同步操作。BlockMatrixBase通过为每个块调用这个辅助函数来自动同步访问。这个函数确保矩阵处于允许添加元素的状态；如果它之前已经处于这种状态，那么这个函数什么也不做。
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

    /**
     * 执行矩阵-矩阵乘法 $C = AB$
     * 的基础函数，或者，如果给出一个大小与B兼容的向量
     * $V$ ，则为 $C = A \text{diag}(V) B$ ，其中 $\text{diag}(V)$
     * 定义了一个带有向量项的对角矩阵。
     * 这个函数假定调用矩阵 $A$ 和 $B$ 的大小兼容。 $C$
     * 的大小将在本函数中设置。        矩阵 $C$
     * 的内容和稀疏模式将被这个函数重置，所以要确保稀疏模式没有在你的程序中其他地方使用。这是一个昂贵的操作，所以在使用这个函数之前请三思。
     *
     */
    void
    mmult(MatrixBase &C, const MatrixBase &B, const VectorBase &V) const;

    /**
     * 基准函数，用于执行矩阵与<tt>this</tt>的转置相乘，即
     * $C = A^T B$ ，或者，如果给出一个可选矢量 $V$
     * ，其大小与 $B$ 兼容，则为 $C = A^T \text{diag}(V) B$
     * ，其中 $\text{diag}(V)$
     * 定义了一个带有矢量项的对角矩阵。
     * 这个函数假设调用矩阵 $A$ 和 $B$ 的大小兼容。 $C$
     * 的大小将在本函数中设置。        矩阵 $C$
     * 的内容和稀疏模式将被这个函数改变，所以要确保稀疏模式没有在你的程序中其他地方使用。这是一个昂贵的操作，所以在使用这个函数之前请三思。
     *
     */
    void
    Tmmult(MatrixBase &C, const MatrixBase &B, const VectorBase &V) const;

  private:
    /**
     * 一个内部的整数数组，当向（大的）稀疏矩阵添加/插入本地数据时，用于存储列索引。
     * 这个变量并不存储矩阵对象的任何
     * "状态"。相反，它只是被这个类的一些成员函数用作临时缓冲区。与所有的
     * @p mutable
     * 成员变量一样，除非有突变器的保护，否则对这个变量的使用不是线程安全的。然而，由于PETSc矩阵操作无论如何都不是线程安全的，所以没有必要试图使事情变得线程安全，所以没有与这个变量相关的突变。
     *
     */
    mutable std::vector<PetscInt> column_indices;

    /**
     * 一个内部的双值数组，在向（大的）稀疏矩阵添加/插入本地数据时，用于存储列索引。
     * 与上面 @p column_indices 变量的注释相同。
     *
     */
    mutable std::vector<PetscScalar> column_values;


    // To allow calling protected prepare_add() and prepare_set().
    template <class>
    friend class dealii::BlockMatrixBase;
  };



#    ifndef DOXYGEN
  // ---------------------- inline and template functions ---------------------


  namespace MatrixIterators
  {
    inline const_iterator::Accessor::Accessor(const MatrixBase *matrix,
                                              const size_type   row,
                                              const size_type   index)
      : matrix(const_cast<MatrixBase *>(matrix))
      , a_row(row)
      , a_index(index)
    {
      visit_present_row();
    }



    inline const_iterator::Accessor::size_type
    const_iterator::Accessor::row() const
    {
      Assert(a_row < matrix->m(), ExcBeyondEndOfMatrix());
      return a_row;
    }


    inline const_iterator::Accessor::size_type
    const_iterator::Accessor::column() const
    {
      Assert(a_row < matrix->m(), ExcBeyondEndOfMatrix());
      return (*colnum_cache)[a_index];
    }


    inline const_iterator::Accessor::size_type
    const_iterator::Accessor::index() const
    {
      Assert(a_row < matrix->m(), ExcBeyondEndOfMatrix());
      return a_index;
    }


    inline PetscScalar
    const_iterator::Accessor::value() const
    {
      Assert(a_row < matrix->m(), ExcBeyondEndOfMatrix());
      return (*value_cache)[a_index];
    }


    inline const_iterator::const_iterator(const MatrixBase *matrix,
                                          const size_type   row,
                                          const size_type   index)
      : accessor(matrix, row, index)
    {}



    inline const_iterator &
    const_iterator::operator++()
    {
      Assert(accessor.a_row < accessor.matrix->m(), ExcIteratorPastEnd());

      ++accessor.a_index;

      // if at end of line: do one step, then cycle until we find a
      // row with a nonzero number of entries
      if (accessor.a_index >= accessor.colnum_cache->size())
        {
          accessor.a_index = 0;
          ++accessor.a_row;

          while ((accessor.a_row < accessor.matrix->m()) &&
                 (accessor.a_row < accessor.matrix->local_range().second) &&
                 (accessor.matrix->row_length(accessor.a_row) == 0))
            ++accessor.a_row;

          accessor.visit_present_row();
        }
      return *this;
    }


    inline const_iterator
    const_iterator::operator++(int)
    {
      const const_iterator old_state = *this;
      ++(*this);
      return old_state;
    }


    inline const const_iterator::Accessor &const_iterator::operator*() const
    {
      return accessor;
    }


    inline const const_iterator::Accessor *const_iterator::operator->() const
    {
      return &accessor;
    }


    inline bool
    const_iterator::operator==(const const_iterator &other) const
    {
      return (accessor.a_row == other.accessor.a_row &&
              accessor.a_index == other.accessor.a_index);
    }


    inline bool
    const_iterator::operator!=(const const_iterator &other) const
    {
      return !(*this == other);
    }


    inline bool
    const_iterator::operator<(const const_iterator &other) const
    {
      return (accessor.row() < other.accessor.row() ||
              (accessor.row() == other.accessor.row() &&
               accessor.index() < other.accessor.index()));
    }

  } // namespace MatrixIterators



  // Inline the set() and add()
  // functions, since they will be
  // called frequently, and the
  // compiler can optimize away
  // some unnecessary loops when
  // the sizes are given at
  // compile time.
  inline void
  MatrixBase::set(const size_type i, const size_type j, const PetscScalar value)
  {
    AssertIsFinite(value);

    set(i, 1, &j, &value, false);
  }



  inline void
  MatrixBase::set(const std::vector<size_type> & indices,
                  const FullMatrix<PetscScalar> &values,
                  const bool                     elide_zero_values)
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
  MatrixBase::set(const std::vector<size_type> & row_indices,
                  const std::vector<size_type> & col_indices,
                  const FullMatrix<PetscScalar> &values,
                  const bool                     elide_zero_values)
  {
    Assert(row_indices.size() == values.m(),
           ExcDimensionMismatch(row_indices.size(), values.m()));
    Assert(col_indices.size() == values.n(),
           ExcDimensionMismatch(col_indices.size(), values.n()));

    for (size_type i = 0; i < row_indices.size(); ++i)
      set(row_indices[i],
          col_indices.size(),
          col_indices.data(),
          &values(i, 0),
          elide_zero_values);
  }



  inline void
  MatrixBase::set(const size_type                 row,
                  const std::vector<size_type> &  col_indices,
                  const std::vector<PetscScalar> &values,
                  const bool                      elide_zero_values)
  {
    Assert(col_indices.size() == values.size(),
           ExcDimensionMismatch(col_indices.size(), values.size()));

    set(row,
        col_indices.size(),
        col_indices.data(),
        values.data(),
        elide_zero_values);
  }



  inline void
  MatrixBase::set(const size_type    row,
                  const size_type    n_cols,
                  const size_type *  col_indices,
                  const PetscScalar *values,
                  const bool         elide_zero_values)
  {
    prepare_action(VectorOperation::insert);

    const PetscInt  petsc_i = row;
    PetscInt const *col_index_ptr;

    PetscScalar const *col_value_ptr;
    int                n_columns;

    // If we don't elide zeros, the pointers are already available...
    if (elide_zero_values == false)
      {
        col_index_ptr = reinterpret_cast<const PetscInt *>(col_indices);
        col_value_ptr = values;
        n_columns     = n_cols;
      }
    else
      {
        // Otherwise, extract nonzero values in each row and get the
        // respective index.
        if (column_indices.size() < n_cols)
          {
            column_indices.resize(n_cols);
            column_values.resize(n_cols);
          }

        n_columns = 0;
        for (size_type j = 0; j < n_cols; ++j)
          {
            const PetscScalar value = values[j];
            AssertIsFinite(value);
            if (value != PetscScalar())
              {
                column_indices[n_columns] = col_indices[j];
                column_values[n_columns]  = value;
                n_columns++;
              }
          }
        AssertIndexRange(n_columns, n_cols + 1);

        col_index_ptr = column_indices.data();
        col_value_ptr = column_values.data();
      }

    const PetscErrorCode ierr = MatSetValues(matrix,
                                             1,
                                             &petsc_i,
                                             n_columns,
                                             col_index_ptr,
                                             col_value_ptr,
                                             INSERT_VALUES);
    AssertThrow(ierr == 0, ExcPETScError(ierr));
  }



  inline void
  MatrixBase::add(const size_type i, const size_type j, const PetscScalar value)
  {
    AssertIsFinite(value);

    if (value == PetscScalar())
      {
        // we have to check after using Insert/Add in any case to be
        // consistent with the MPI communication model, but we can save
        // some work if the addend is zero. However, these actions are done
        // in case we pass on to the other function.
        prepare_action(VectorOperation::add);

        return;
      }
    else
      add(i, 1, &j, &value, false);
  }



  inline void
  MatrixBase::add(const std::vector<size_type> & indices,
                  const FullMatrix<PetscScalar> &values,
                  const bool                     elide_zero_values)
  {
    Assert(indices.size() == values.m(),
           ExcDimensionMismatch(indices.size(), values.m()));
    Assert(values.m() == values.n(), ExcNotQuadratic());

    for (size_type i = 0; i < indices.size(); ++i)
      add(indices[i],
          indices.size(),
          indices.data(),
          &values(i, 0),
          elide_zero_values);
  }



  inline void
  MatrixBase::add(const std::vector<size_type> & row_indices,
                  const std::vector<size_type> & col_indices,
                  const FullMatrix<PetscScalar> &values,
                  const bool                     elide_zero_values)
  {
    Assert(row_indices.size() == values.m(),
           ExcDimensionMismatch(row_indices.size(), values.m()));
    Assert(col_indices.size() == values.n(),
           ExcDimensionMismatch(col_indices.size(), values.n()));

    for (size_type i = 0; i < row_indices.size(); ++i)
      add(row_indices[i],
          col_indices.size(),
          col_indices.data(),
          &values(i, 0),
          elide_zero_values);
  }



  inline void
  MatrixBase::add(const size_type                 row,
                  const std::vector<size_type> &  col_indices,
                  const std::vector<PetscScalar> &values,
                  const bool                      elide_zero_values)
  {
    Assert(col_indices.size() == values.size(),
           ExcDimensionMismatch(col_indices.size(), values.size()));

    add(row,
        col_indices.size(),
        col_indices.data(),
        values.data(),
        elide_zero_values);
  }



  inline void
  MatrixBase::add(const size_type    row,
                  const size_type    n_cols,
                  const size_type *  col_indices,
                  const PetscScalar *values,
                  const bool         elide_zero_values,
                  const bool  /*col_indices_are_sorted*/ )
  {
    (void)elide_zero_values;

    prepare_action(VectorOperation::add);

    const PetscInt  petsc_i = row;
    PetscInt const *col_index_ptr;

    PetscScalar const *col_value_ptr;
    int                n_columns;

    // If we don't elide zeros, the pointers are already available...
    if (elide_zero_values == false)
      {
        col_index_ptr = reinterpret_cast<const PetscInt *>(col_indices);
        col_value_ptr = values;
        n_columns     = n_cols;
      }
    else
      {
        // Otherwise, extract nonzero values in each row and get the
        // respective index.
        if (column_indices.size() < n_cols)
          {
            column_indices.resize(n_cols);
            column_values.resize(n_cols);
          }

        n_columns = 0;
        for (size_type j = 0; j < n_cols; ++j)
          {
            const PetscScalar value = values[j];
            AssertIsFinite(value);
            if (value != PetscScalar())
              {
                column_indices[n_columns] = col_indices[j];
                column_values[n_columns]  = value;
                n_columns++;
              }
          }
        AssertIndexRange(n_columns, n_cols + 1);

        col_index_ptr = column_indices.data();
        col_value_ptr = column_values.data();
      }

    const PetscErrorCode ierr = MatSetValues(
      matrix, 1, &petsc_i, n_columns, col_index_ptr, col_value_ptr, ADD_VALUES);
    AssertThrow(ierr == 0, ExcPETScError(ierr));
  }



  inline PetscScalar
  MatrixBase::operator()(const size_type i, const size_type j) const
  {
    return el(i, j);
  }



  inline MatrixBase::const_iterator
  MatrixBase::begin() const
  {
    Assert(
      (in_local_range(0) && in_local_range(m() - 1)),
      ExcMessage(
        "begin() and end() can only be called on a processor owning the entire matrix. If this is a distributed matrix, use begin(row) and end(row) instead."));

    // find the first non-empty row in order to make sure that the returned
    // iterator points to something useful
    size_type first_nonempty_row = 0;
    while ((first_nonempty_row < m()) && (row_length(first_nonempty_row) == 0))
      ++first_nonempty_row;

    return const_iterator(this, first_nonempty_row, 0);
  }


  inline MatrixBase::const_iterator
  MatrixBase::end() const
  {
    Assert(
      (in_local_range(0) && in_local_range(m() - 1)),
      ExcMessage(
        "begin() and end() can only be called on a processor owning the entire matrix. If this is a distributed matrix, use begin(row) and end(row) instead."));

    return const_iterator(this, m(), 0);
  }


  inline MatrixBase::const_iterator
  MatrixBase::begin(const size_type r) const
  {
    Assert(in_local_range(r),
           ExcIndexRange(r, local_range().first, local_range().second));

    if (row_length(r) > 0)
      return const_iterator(this, r, 0);
    else
      return end(r);
  }


  inline MatrixBase::const_iterator
  MatrixBase::end(const size_type r) const
  {
    Assert(in_local_range(r),
           ExcIndexRange(r, local_range().first, local_range().second));

    // place the iterator on the first entry past this line, or at the
    // end of the matrix
    //
    // in the parallel case, we need to put it on the first entry of
    // the first row after the locally owned range. this of course
    // doesn't exist, but we can nevertheless create such an
    // iterator. we need to check whether 'i' is past the locally
    // owned range of rows first, before we ask for the length of the
    // row since the latter query leads to an exception in case we ask
    // for a row that is not locally owned
    for (size_type i = r + 1; i < m(); ++i)
      if (i == local_range().second || (row_length(i) > 0))
        return const_iterator(this, i, 0);

    // if there is no such line, then take the
    // end iterator of the matrix
    // we don't allow calling end() directly for distributed matrices so we need
    // to copy the code without the assertion.
    return {this, m(), 0};
  }



  inline bool
  MatrixBase::in_local_range(const size_type index) const
  {
    PetscInt begin, end;

    const PetscErrorCode ierr =
      MatGetOwnershipRange(static_cast<const Mat &>(matrix), &begin, &end);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    return ((index >= static_cast<size_type>(begin)) &&
            (index < static_cast<size_type>(end)));
  }



  inline void
  MatrixBase::prepare_action(const VectorOperation::values new_action)
  {
    if (last_action == VectorOperation::unknown)
      last_action = new_action;

    Assert(last_action == new_action, ExcWrongMode(last_action, new_action));
  }



  inline void
  MatrixBase::assert_is_compressed()
  {
    // compress() sets the last action to none, which allows us to check if
    // there are pending add/insert operations:
    AssertThrow(last_action == VectorOperation::unknown,
                ExcMessage("Error: missing compress() call."));
  }



  inline void
  MatrixBase::prepare_add()
  {
    prepare_action(VectorOperation::add);
  }



  inline void
  MatrixBase::prepare_set()
  {
    prepare_action(VectorOperation::insert);
  }

#    endif // DOXYGEN
} // namespace PETScWrappers


DEAL_II_NAMESPACE_CLOSE


#  endif // DEAL_II_WITH_PETSC

#endif
 /*---------------------------- petsc_matrix_base.h --------------------------*/ 


