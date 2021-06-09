//include/deal.II-translator/lac/trilinos_sparsity_pattern_0.txt
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

#ifndef dealii_trilinos_sparsity_pattern_h
#  define dealii_trilinos_sparsity_pattern_h


#  include <deal.II/base/config.h>

#  ifdef DEAL_II_WITH_TRILINOS

#    include <deal.II/base/index_set.h>
#    include <deal.II/base/subscriptor.h>

#    include <deal.II/lac/exceptions.h>

#    include <Epetra_FECrsGraph.h>
#    include <Epetra_Map.h>

#    include <cmath>
#    include <memory>
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
class SparsityPattern;
class DynamicSparsityPattern;

namespace TrilinosWrappers
{
  class SparsityPattern;
  class SparseMatrix;

  namespace SparsityPatternIterators
  {
    class Iterator;
  }
} // namespace TrilinosWrappers
#    endif

namespace TrilinosWrappers
{
  namespace SparsityPatternIterators
  {
    /**
     * 迭代器进入稀疏性模式的访问器类。这个类也是进入稀疏矩阵的常量和非常量访问器类的基类。
     * 请注意，这个类只允许对元素进行读取访问，提供它们的行和列号。它不允许修改稀疏模式本身。
     * @ingroup TrilinosWrappers
     *
     */
    class Accessor
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
      Accessor(const SparsityPattern *sparsity_pattern,
               const size_type        row,
               const size_type        index);

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

      /**
       * 异常情况
       *
       */
      DeclException0(ExcBeyondEndOfSparsityPattern);

      /**
       * 异常情况
       *
       */
      DeclException3(ExcAccessToNonlocalRow,
                     size_type,
                     size_type,
                     size_type,
                     << "You tried to access row " << arg1
                     << " of a distributed sparsity pattern, "
                     << " but only rows " << arg2 << " through " << arg3
                     << " are stored locally and can be accessed.");

    private:
      /**
       * 访问的矩阵。
       *
       */
      mutable SparsityPattern *sparsity_pattern;

      /**
       * 当前行数。
       *
       */
      size_type a_row;

      /**
       * 当前行中的索引。
       *
       */
      size_type a_index;

      /**
       * 缓存，我们在这里存储当前行的列索引。这是必要的，因为Trilinos对其矩阵元素的访问相当困难，当我们进入某一行时，一次性复制该行的所有列项，比反复向Trilinos索取单个列项要高效得多。这也有一定的意义，因为无论如何，我们很可能会按顺序访问它们。
       * 为了使迭代器/存取器的复制具有可接受的性能，我们为这些条目保留了一个共享指针，以便在必要时有多个存取器可以访问这些数据。
       *
       */
      std::shared_ptr<const std::vector<size_type>> colnum_cache;

      /**
       * 丢弃旧的行缓存（它们可能仍然被其他访问器使用），并为这个访问器目前所指向的行生成新的行缓存。
       *
       */
      void
      visit_present_row();

      // Make enclosing class a friend.
      friend class Iterator;
    };

    /**
     * 类型为 TrilinosWrappers::SparsityPattern.
     * 的稀疏模式的迭代器类，对稀疏模式的单个元素的访问由该命名空间的访问器类处理。
     *
     */
    class Iterator
    {
    public:
      /**
       * 声明容器大小的类型。
       *
       */
      using size_type = dealii::types::global_dof_index;

      /**
       * 构造器。为给定的行和其中的索引创建一个进入矩阵
       * @p matrix 的迭代器。
       *
       */
      Iterator(const SparsityPattern *sparsity_pattern,
               const size_type        row,
               const size_type        index);

      /**
       * 复制构造函数。
       *
       */
      Iterator(const Iterator &i);

      /**
       * 前缀增量。
       *
       */
      Iterator &
      operator++();

      /**
       * 后缀增量。
       *
       */
      Iterator
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
      operator==(const Iterator &) const;

      /**
       * <tt>==</tt>的倒数。
       *
       */
      bool
      operator!=(const Iterator &) const;

      /**
       * 比较运算符。如果第一行数字较小，或者行数字相等且第一个索引较小，则结果为真。
       *
       */
      bool
      operator<(const Iterator &) const;

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
      Accessor accessor;

      friend class TrilinosWrappers::SparsityPattern;
    };

  } // namespace SparsityPatternIterators


  /**
   * 该类实现了一个包装类，用于使用Trilinos分布式稀疏性模式类Epetra_FECrsGraph。该类被设计用于构建%平行的Tridinos矩阵。该类的功能是以现有的稀疏模式类为模型，不同的是，该类可以根据稀疏模式行的划分，完全并行地工作。
   * 这个类与DynamicSparsityPattern有很多相似之处，因为它可以动态地添加元素到模式中，而不需要事先为其保留任何内存。然而，它也有一个方法
   * SparsityPattern::compress(),
   * 可以最终确定该模式，并使其能够与特里诺斯稀疏矩阵一起使用。
   * @ingroup TrilinosWrappers
   * @ingroup Sparsity
   *
   */
  class SparsityPattern : public Subscriptor
  {
  public:
    /**
     * 声明容器尺寸的类型。
     *
     */
    using size_type = dealii::types::global_dof_index;

    /**
     * 为迭代器类声明一个别名。
     *
     */
    using const_iterator = SparsityPatternIterators::Iterator;

    /**
     * @name  基本构造函数和初始化
     *
     */
    //@{
    /**
     * 默认构造函数。生成一个空的（零大小）稀疏模式。
     *
     */
    SparsityPattern();

    /**
     * 生成一个完全本地存储的稀疏模式，有 $m$ 行和 $n$
     * 列。产生的矩阵也将完全存储在本地。
     * 可以使用可选的 @p n_entries_per_row
     * 参数来指定每行的列条目数。然而，这个值不需要准确，甚至根本不需要给出，因为在建立稀疏模式之前，人们通常没有这种信息（通常情况下，当函数
     * DoFTools::make_sparsity_pattern()
     * 被调用）。条目是以类似于deal.II
     * DynamicSparsityPattern类的方式动态分配的。然而，一个好的估计将减少稀疏模式的设置时间。
     *
     */
    SparsityPattern(const size_type m,
                    const size_type n,
                    const size_type n_entries_per_row = 0);

    /**
     * 生成一个完全存储在本地的稀疏度模式，有 $m$ 行和
     * $n$ 列。产生的矩阵也将完全存储在本地。
     * 向量<tt>n_entries_per_row</tt>指定了每一行的条目数（不过这个信息通常是不可用的）。
     *
     */
    SparsityPattern(const size_type               m,
                    const size_type               n,
                    const std::vector<size_type> &n_entries_per_row);

    /**
     * 移动构造函数。通过窃取内部数据创建一个新的稀疏矩阵。
     *
     */
    SparsityPattern(SparsityPattern &&other) noexcept;

    /**
     * 复制构造函数。将调用的稀疏模式设置为与输入的稀疏模式相同。
     *
     */
    SparsityPattern(const SparsityPattern &input_sparsity_pattern);

    /**
     * 解除构造函数。虚化，以便人们可以使用指向该类的指针。
     *
     */
    virtual ~SparsityPattern() override = default;

    /**
     * 初始化一个完全存储在本地的稀疏度模式，有 $m$ 行和
     * $n$ 列。由此产生的矩阵将被完全存储在本地。
     * 每行的列条目数被指定为最大条目数参数。
     * 这不需要是一个准确的数字，因为条目是以类似于deal.II
     * DynamicSparsityPattern类的方式动态分配的，但是一个好的估计将减少稀疏模式的设置时间。
     *
     */
    void
    reinit(const size_type m,
           const size_type n,
           const size_type n_entries_per_row = 0);

    /**
     * 初始化一个完全存储在本地的稀疏度模式，有 $m$ 行和
     * $n$ 列。由此产生的矩阵将被完全存储在本地。
     * 向量<tt>n_entries_per_row</tt>指定每一行的条目数。
     *
     */
    void
    reinit(const size_type               m,
           const size_type               n,
           const std::vector<size_type> &n_entries_per_row);

    /**
     * 复制功能。将调用的稀疏度模式设置为与输入的稀疏度模式相同。
     *
     */
    void
    copy_from(const SparsityPattern &input_sparsity_pattern);

    /**
     * 拷贝函数，来自一个deal.II稀疏度模式。如果并行使用，该函数使用行和列的临时划分。
     *
     */
    template <typename SparsityPatternType>
    void
    copy_from(const SparsityPatternType &nontrilinos_sparsity_pattern);

    /**
     * 复制操作。这个操作只允许用于空对象，以避免由编译器自动合成的潜在的非常昂贵的操作。如果你知道你真的想复制一个具有非琐碎内容的稀疏模式，请使用copy_from()代替。
     *
     */
    SparsityPattern &
    operator=(const SparsityPattern &input_sparsity_pattern);

    /**
     * 释放所有内存并返回到与调用默认构造函数后相同的状态。
     * 这是一个集体操作，需要在所有处理器上调用，以避免出现死锁。
     *
     */
    void
    clear();

    /**
     * 与我们自己的SparsityPattern类相类似，这个函数压缩了稀疏模式，并允许产生的模式用于实际生成一个（基于Trilinos）矩阵。这个函数还交换了在添加新元素过程中可能积累的非局部数据。因此，一旦结构被固定，就必须调用这个函数。这是一个集体操作，即在并行使用时需要在所有处理器上运行。
     *
     */
    void
    compress();
    //@}

    /**
     * @name  使用IndexSet描述的构造器和初始化
     *
     */
    //@{

    /**
     * 使用一个IndexSet和一个MPI通信器的方形稀疏模式的构造器，用于描述%并行分区。
     * 此外，可以指定稀疏模式的行中非零条目的数量。请注意，这个数字不需要精确，甚至允许实际的稀疏结构有比构造函数中指定的更多非零条目。但是在这里提供好的估计值仍然是有利的，因为一个好的值可以避免重复分配内存，这大大增加了创建稀疏度模式时的性能。
     *
     */
    SparsityPattern(const IndexSet &parallel_partitioning,
                    const MPI_Comm &communicator      = MPI_COMM_WORLD,
                    const size_type n_entries_per_row = 0);

    /**
     * 与之前相同，但现在使用每m行中非零的确切数量。
     * 由于在这种情况下我们确切地知道稀疏模式中的元素数，我们已经可以分配适当数量的内存，这使得各自的
     * SparsityPattern::reinit
     * 调用的创建过程大大加快。然而，这是一个相当不寻常的情况，因为知道每一行的条目数通常与知道非零条目的指数有关，而稀疏模式就是为了描述这些非零条目。
     *
     */
    SparsityPattern(const IndexSet &              parallel_partitioning,
                    const MPI_Comm &              communicator,
                    const std::vector<size_type> &n_entries_per_row);

    /**
     * 这个构造函数与上面的构造函数类似，但是它现在需要两个不同的索引集来描述行和列的%平行分割。这个接口是为了用于生成矩形稀疏模式。请注意，沿着列没有真正的并行性；拥有某一行的处理器总是拥有所有的列元素，不管它们可能分散得多远。第二个Epetra_Map仅用于指定列数，以及在基于该列图与向量进行矩阵-向量乘积时的内部安排。
     * 每行的列条目数被指定为最大条目数的参数。
     *
     */
    SparsityPattern(const IndexSet &row_parallel_partitioning,
                    const IndexSet &col_parallel_partitioning,
                    const MPI_Comm &communicator      = MPI_COMM_WORLD,
                    const size_type n_entries_per_row = 0);

    /**
     * 这个构造函数与上面的构造函数类似，但它现在需要两个不同的行和列的索引集。这个接口是用来生成矩形矩阵的，其中一个映射指定了行的%平行分布，第二个映射指定了与矩阵列相关的自由度分布。然而，这第二个映射并不用于列本身的分布&ndash；相反，一行的所有列元素都存储在同一个处理器上。向量<tt>n_entries_per_row</tt>指定了新生成的矩阵中每一行的条目数。
     *
     */
    SparsityPattern(const IndexSet &              row_parallel_partitioning,
                    const IndexSet &              col_parallel_partitioning,
                    const MPI_Comm &              communicator,
                    const std::vector<size_type> &n_entries_per_row);

    /**
     * 这个构造函数可以构造一般的稀疏模式，可能是非方形的。通过这种方式构建稀疏模式，用户可以明确地指定我们要添加元素的行。
     * 这个集合被要求是第一个索引集 @p
     * row_parallel_partitioning的超集，其中也包括被另一个处理器拥有的行（ghost
     * rows）。注意，元素只能被添加到 @p writable_rows.
     * 指定的行中。当处理器要写入的行可以在实际插入元素到矩阵中之前被确定时，这种方法是有益的。对于deal.II中使用的典型的
     * parallel::distributed::Triangulation
     * 类，我们知道处理器只为我们所说的本地相关的道夫添加行元素（见
     * DoFTools::extract_locally_relevant_dofs). ）。
     * 其他构造函数方法使用一般的Trilinos设施，允许向任意的行添加元素（正如所有其他reinit函数所做的那样）。然而，这种灵活性是有代价的，最突出的是，只要使用MPI，从共享内存中的多个线程向同一矩阵添加元素是不安全的。对于这些设置，目前的方法是可以选择的。它将把处理器外的数据存储为一个额外的稀疏模式（然后通过reinit方法传递给Trilinos矩阵），其组织方式可以确保线程安全（当然，只要用户确保不同时写入同一矩阵行）。
     *
     */
    SparsityPattern(const IndexSet &row_parallel_partitioning,
                    const IndexSet &col_parallel_partitioning,
                    const IndexSet &writable_rows,
                    const MPI_Comm &communicator      = MPI_COMM_WORLD,
                    const size_type n_entries_per_row = 0);

    /**
     * 重新初始化函数，用于生成一个方形稀疏模式，使用IndexSet和MPI通信器来描述%并行分区和稀疏模式的行中非零项的数量。请注意，这个数字不需要精确，甚至允许实际的稀疏结构有比构造函数中指定的更多的非零条目。然而，在这里提供良好的估计仍然是有利的，因为这将大大增加创建稀疏模式时的性能。
     * 这个函数本身不创建任何条目，但提供了正确的数据结构，可以被相应的add()函数使用。
     *
     */
    void
    reinit(const IndexSet &parallel_partitioning,
           const MPI_Comm &communicator      = MPI_COMM_WORLD,
           const size_type n_entries_per_row = 0);

    /**
     * 与之前相同，但现在使用每m行中非零的确切数量。
     * 由于在这种情况下我们确切地知道稀疏模式中的元素数，我们已经可以分配适当数量的内存，这使得向稀疏模式添加条目的过程大大加快。然而，这是一个相当不寻常的情况，因为知道每一行的条目数通常与知道非零条目的指数有关，而稀疏模式就是为了描述这些非零条目。
     *
     */
    void
    reinit(const IndexSet &              parallel_partitioning,
           const MPI_Comm &              communicator,
           const std::vector<size_type> &n_entries_per_row);

    /**
     * 这个reinit函数与上面的函数类似，但它现在需要两个不同的行和列的索引集。这个接口是用来生成矩形稀疏模式的，其中一个索引集描述了与稀疏模式行相关的%平行分割，另一个是稀疏模式列的平行分割。请注意，沿着列没有真正的并行性；拥有某一行的处理器总是拥有所有的列元素，不管它们可能分散得多远。第二个IndexSet仅用于指定列数，以及在与基于该IndexSet的EpetraMap的向量进行矩阵-向量乘积时用于内部安排。
     * 每行的列条目数由参数<tt>n_entries_per_row</tt>指定。
     *
     */
    void
    reinit(const IndexSet &row_parallel_partitioning,
           const IndexSet &col_parallel_partitioning,
           const MPI_Comm &communicator      = MPI_COMM_WORLD,
           const size_type n_entries_per_row = 0);

    /**
     * 这个reinit函数用于指定一般的矩阵，可能是非方形的。除了上述其他reinit方法的参数外，它允许用户明确指定我们要添加元素的行。这个集合是第一个索引集
     * @p row_parallel_partitioning
     * 的超集，也包括被另一个处理器拥有的行（ghost
     * rows）。
     * 当一个处理器要写入的行可以在实际插入元素到矩阵之前确定时，这种方法是有益的。对于deal.II中使用的典型
     * parallel::distributed::Triangulation
     * 类，我们知道处理器只为我们所说的本地相关的道夫添加行元素（见
     * DoFTools::extract_locally_relevant_dofs).
     * ）Trilinos矩阵允许向任意的行添加元素（正如所有其他reinit函数所做的那样），这也是所有其他reinit方法所做的。
     * 然而，这种灵活性是有代价的，最突出的是在使用MPI时，从共享内存的多个线程向同一个矩阵添加元素是不安全的。对于这些设置，当前的方法是值得选择的。它将把处理器外的数据存储为一个额外的稀疏模式（然后通过reinit方法传递给Trilinos矩阵），其组织方式可以确保线程安全（当然，只要用户确保不同时写入同一矩阵行）。
     *
     */
    void
    reinit(const IndexSet &row_parallel_partitioning,
           const IndexSet &col_parallel_partitioning,
           const IndexSet &writeable_rows,
           const MPI_Comm &communicator      = MPI_COMM_WORLD,
           const size_type n_entries_per_row = 0);

    /**
     * 和以前一样，但现在使用一个向量<tt>n_entries_per_row</tt>来指定稀疏模式的每一行的条目数。
     *
     */
    void
    reinit(const IndexSet &              row_parallel_partitioning,
           const IndexSet &              col_parallel_partitioning,
           const MPI_Comm &              communicator,
           const std::vector<size_type> &n_entries_per_row);

    /**
     * Reinit函数。接受一个deal.II稀疏模式和由两个索引集指定的行和列的%平行分割，以及一个用于初始化当前Trilinos稀疏模式的%平行通信器。可选参数
     * @p exchange_data
     * 可用于用未完全构建的稀疏模式进行重新初始化。这个功能只对动态稀疏模式类型的输入稀疏模式实现。
     *
     */
    template <typename SparsityPatternType>
    void
    reinit(const IndexSet &           row_parallel_partitioning,
           const IndexSet &           col_parallel_partitioning,
           const SparsityPatternType &nontrilinos_sparsity_pattern,
           const MPI_Comm &           communicator  = MPI_COMM_WORLD,
           const bool                 exchange_data = false);

    /**
     * Reinit函数。接受一个deal.II稀疏度模式和一个%平行分割的行和列，用于初始化当前Trilinos稀疏度模式。可选参数
     * @p
     * exchange_data可用于重新初始化未完全构建的稀疏度模式。这个功能只对动态稀疏模式类型的输入稀疏模式实现。
     *
     */
    template <typename SparsityPatternType>
    void
    reinit(const IndexSet &           parallel_partitioning,
           const SparsityPatternType &nontrilinos_sparsity_pattern,
           const MPI_Comm &           communicator  = MPI_COMM_WORLD,
           const bool                 exchange_data = false);
    //@}
    /**
     * @name  稀疏度模式的信息
     *
     */
    //@{

    /**
     * 返回稀疏性模式的状态，即在需要数据交换的操作之后是否需要调用compress()。
     *
     */
    bool
    is_compressed() const;

    /**
     * 返回当前处理器上每行的最大条目数。
     *
     */
    unsigned int
    max_entries_per_row() const;

    /**
     * 返回该稀疏模式中的行数。
     *
     */
    size_type
    n_rows() const;

    /**
     * 返回这个稀疏模式中的列数。
     *
     */
    size_type
    n_cols() const;

    /**
     * 返回稀疏模式的局部维度，即存储在当前MPI进程中的行数。在顺序情况下，这个数字与n_rows()相同，但是对于并行矩阵，这个数字可能更小。
     * 要想知道到底哪些元素被存储在本地，可以使用local_range()。
     *
     */
    unsigned int
    local_size() const;

    /**
     * 返回一对指数，表明该稀疏模式的哪些行被存储在本地。第一个数字是存储的第一行的索引，第二个数字是本地存储的最后一行之后的那一行的索引。如果这是一个顺序矩阵，那么结果将是一对（0,n_rows()），否则将是一对（i,i+n），其中<tt>n=local_size()</tt>。
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
     * 返回这个稀疏模式的非零元素的数量。
     *
     */
    size_type
    n_nonzero_elements() const;

    /**
     * 返回给定行中的条目数。
     * 在一个并行的环境中，有关的行当然可能不存储在当前的处理器上，在这种情况下，就不可能查询其中的条目数。在这种情况下，返回值是`static_cast<size_type>(-1)`。
     *
     */
    size_type
    row_length(const size_type row) const;

    /**
     * 计算这个结构所代表的矩阵的带宽。该带宽是 $|i-j|$
     * 的最大值，其中索引对 $(i,j)$
     * 代表矩阵的非零条目。因此， $n\times m$
     * 矩阵的最大带宽是 $\max\{n-1,m-1\}$  。
     *
     */
    size_type
    bandwidth() const;

    /**
     * 返回该对象是否为空。如果没有分配内存，它就是空的，这与两个维度都是0时的情况相同。
     *
     */
    bool
    empty() const;

    /**
     * 返回索引（<i>i,j</i>）是否存在于稀疏模式中（即它可能是非零）。
     *
     */
    bool
    exists(const size_type i, const size_type j) const;

    /**
     * 返回给定的 @p row
     * 是否被存储在这个进程的当前对象中。
     *
     */
    bool
    row_is_stored_locally(const size_type i) const;

    /**
     * 确定这个对象的内存消耗（以字节为单位）的估计值。目前这个类没有实现。
     *
     */
    std::size_t
    memory_consumption() const;

    //@}
    /**
     * @name 添加条目
     *
     */
    //@{
    /**
     * 将元素（<i>i,j</i>）添加到稀疏模式中。
     *
     */
    void
    add(const size_type i, const size_type j);


    /**
     * 在一行中添加几个元素到稀疏模式中。
     *
     */
    template <typename ForwardIterator>
    void
    add_entries(const size_type row,
                ForwardIterator begin,
                ForwardIterator end,
                const bool      indices_are_sorted = false);
    //@}
    /**
     * @name 访问底层的Trilinos数据
     *
     */
    //@{

    /**
     * 返回一个对存储稀疏模式的底层Trilinos
     * Epetra_CrsGraph数据的常量引用。
     *
     */
    const Epetra_FECrsGraph &
    trilinos_sparsity_pattern() const;

    /**
     * 返回一个对底层Trilinos
     * Epetra_Map的常量引用，该Map设置了该稀疏模式的域空间的平行分区，即基于该稀疏模式的向量矩阵的分区与之相乘。
     *
     */
    const Epetra_Map &
    domain_partitioner() const;

    /**
     * 返回一个对底层Trilinos
     * Epetra_Map的常量引用，该Map设置了该稀疏模式的范围空间的划分，即由矩阵-向量乘积产生的向量的划分。
     *
     */
    const Epetra_Map &
    range_partitioner() const;

    /**
     * 返回与该矩阵一起使用的MPI通信器对象。
     *
     */
    MPI_Comm
    get_mpi_communicator() const;
    //@}

    /**
     * @name  分割器
     *
     */
    //@{

    /**
     * 返回该模式的域空间的分区，即基于该稀疏模式的矩阵必须与之相乘的向量的分区。
     *
     */
    IndexSet
    locally_owned_domain_indices() const;

    /**
     * 返回该模式的范围空间划分，即基于该模式的矩阵向量乘积所产生的向量的划分。
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
     * 迭代器从第一个条目开始。
     *
     */
    const_iterator
    begin() const;

    /**
     * 最后的迭代器。
     *
     */
    const_iterator
    end() const;

    /**
     * 从第 @p r. 行的第一个条目开始的迭代器
     * 注意，如果给定的行是空的，即不包含任何非零条目，那么这个函数返回的迭代器就等于<tt>end(r)</tt>。还要注意的是，在这种情况下，迭代器可能不能被解除引用。
     *
     */
    const_iterator
    begin(const size_type r) const;

    /**
     * 行<tt>r</tt>的最终迭代器。它指向过了 @p r,
     * 行末尾的第一个元素，或者过了整个稀疏模式的末尾。
     * 请注意，结束迭代器不一定是可被解除引用的。特别是如果它是一个矩阵的最后一行的结束迭代器，情况更是如此。
     *
     */
    const_iterator
    end(const size_type r) const;

    //@}
    /**
     * @name  输入/输出
     *
     */
    //@{

    /**
     * 抽象的Trilinos对象，帮助在ASCII中查看其他Trilinos对象。目前这个功能还没有实现。
     * TODO：未实现。
     *
     */
    void
    write_ascii();

    /**
     * 使用<tt>(line,col)</tt>格式，向给定的流打印（本地拥有的部分）稀疏模式。可选的标志是以Trilinos风格输出稀疏模式，在实际写入条目之前，甚至根据处理器的编号也会打印到流中，以及一个摘要。
     *
     */
    void
    print(std::ostream &out,
          const bool    write_extended_trilinos_info = false) const;

    /**
     * 以<tt>gnuplot</tt>能理解的格式打印矩阵的稀疏度，该格式可用于以图形方式绘制稀疏度模式。该格式由成对的<tt>i
     * j</tt>非零元素组成，每个元素代表该矩阵的一个条目，输出文件中每行一个。指数从零开始计算，和平常一样。由于稀疏模式的打印方式与矩阵的显示方式相同，我们打印的是列索引的负数，这意味着<tt>(0,0)</tt>元素在左上角而不是左下角。
     * 在gnuplot中通过将数据样式设置为点或点来打印稀疏模式，并使用<tt>plot</tt>命令。
     *
     */
    void
    print_gnuplot(std::ostream &out) const;

    //@}
    /**
     * @addtogroup  Exceptions  @{
     *
     */
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
    DeclExceptionMsg(
      ExcSourceEqualsDestination,
      "You are attempting an operation on two sparsity patterns that "
      "are the same object, but the operation requires that the "
      "two objects are in fact different.");

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
  private:
    /**
     * 指向用户提供的矩阵列的Epetra
     * Trilinos映射的指针，它将矩阵的一部分分配给各个进程。
     *
     */
    std::unique_ptr<Epetra_Map> column_space_map;

    /**
     * Trilinos中的稀疏模式对象，用于基于有限元的问题，允许向模式添加非局部元素。
     *
     */
    std::unique_ptr<Epetra_FECrsGraph> graph;

    /**
     * 稀疏模式的非本地部分的稀疏模式对象，将被发送到拥有的处理器。只有在设置了特定的构造函数或带有writable_rows参数的reinit方法时才使用
     *
     */
    std::unique_ptr<Epetra_CrsGraph> nonlocal_graph;

    friend class TrilinosWrappers::SparseMatrix;
    friend class SparsityPatternIterators::Accessor;
    friend class SparsityPatternIterators::Iterator;
  };



  // ----------------------- inline and template functions --------------------


#    ifndef DOXYGEN

  namespace SparsityPatternIterators
  {
    inline Accessor::Accessor(const SparsityPattern *sp,
                              const size_type        row,
                              const size_type        index)
      : sparsity_pattern(const_cast<SparsityPattern *>(sp))
      , a_row(row)
      , a_index(index)
    {
      visit_present_row();
    }



    inline Accessor::size_type
    Accessor::row() const
    {
      Assert(a_row < sparsity_pattern->n_rows(),
             ExcBeyondEndOfSparsityPattern());
      return a_row;
    }



    inline Accessor::size_type
    Accessor::column() const
    {
      Assert(a_row < sparsity_pattern->n_rows(),
             ExcBeyondEndOfSparsityPattern());
      return (*colnum_cache)[a_index];
    }



    inline Accessor::size_type
    Accessor::index() const
    {
      Assert(a_row < sparsity_pattern->n_rows(),
             ExcBeyondEndOfSparsityPattern());
      return a_index;
    }



    inline Iterator::Iterator(const SparsityPattern *sp,
                              const size_type        row,
                              const size_type        index)
      : accessor(sp, row, index)
    {}



    inline Iterator::Iterator(const Iterator &) = default;



    inline Iterator &
    Iterator::operator++()
    {
      Assert(accessor.a_row < accessor.sparsity_pattern->n_rows(),
             ExcIteratorPastEnd());

      ++accessor.a_index;

      // If at end of line: do one step, then cycle until we find a row with a
      // nonzero number of entries that is stored locally.
      if (accessor.a_index >= accessor.colnum_cache->size())
        {
          accessor.a_index = 0;
          ++accessor.a_row;

          while (accessor.a_row < accessor.sparsity_pattern->n_rows())
            {
              const auto row_length =
                accessor.sparsity_pattern->row_length(accessor.a_row);
              if (row_length == 0 ||
                  !accessor.sparsity_pattern->row_is_stored_locally(
                    accessor.a_row))
                ++accessor.a_row;
              else
                break;
            }

          accessor.visit_present_row();
        }
      return *this;
    }



    inline Iterator
    Iterator::operator++(int)
    {
      const Iterator old_state = *this;
      ++(*this);
      return old_state;
    }



    inline const Accessor &Iterator::operator*() const
    {
      return accessor;
    }



    inline const Accessor *Iterator::operator->() const
    {
      return &accessor;
    }



    inline bool
    Iterator::operator==(const Iterator &other) const
    {
      return (accessor.a_row == other.accessor.a_row &&
              accessor.a_index == other.accessor.a_index);
    }



    inline bool
    Iterator::operator!=(const Iterator &other) const
    {
      return !(*this == other);
    }



    inline bool
    Iterator::operator<(const Iterator &other) const
    {
      return (accessor.row() < other.accessor.row() ||
              (accessor.row() == other.accessor.row() &&
               accessor.index() < other.accessor.index()));
    }

  } // namespace SparsityPatternIterators



  inline SparsityPattern::const_iterator
  SparsityPattern::begin() const
  {
    const size_type first_valid_row = this->local_range().first;
    return const_iterator(this, first_valid_row, 0);
  }



  inline SparsityPattern::const_iterator
  SparsityPattern::end() const
  {
    return const_iterator(this, n_rows(), 0);
  }



  inline SparsityPattern::const_iterator
  SparsityPattern::begin(const size_type r) const
  {
    AssertIndexRange(r, n_rows());
    if (row_length(r) > 0)
      return const_iterator(this, r, 0);
    else
      return end(r);
  }



  inline SparsityPattern::const_iterator
  SparsityPattern::end(const size_type r) const
  {
    AssertIndexRange(r, n_rows());

    // place the iterator on the first entry
    // past this line, or at the end of the
    // matrix
    for (size_type i = r + 1; i < n_rows(); ++i)
      if (row_length(i) > 0)
        return const_iterator(this, i, 0);

    // if there is no such line, then take the
    // end iterator of the matrix
    return end();
  }



  inline bool
  SparsityPattern::in_local_range(const size_type index) const
  {
    TrilinosWrappers::types::int_type begin, end;
#      ifndef DEAL_II_WITH_64BIT_INDICES
    begin = graph->RowMap().MinMyGID();
    end   = graph->RowMap().MaxMyGID() + 1;
#      else
    begin = graph->RowMap().MinMyGID64();
    end   = graph->RowMap().MaxMyGID64() + 1;
#      endif

    return ((index >= static_cast<size_type>(begin)) &&
            (index < static_cast<size_type>(end)));
  }



  inline bool
  SparsityPattern::is_compressed() const
  {
    return graph->Filled();
  }



  inline bool
  SparsityPattern::empty() const
  {
    return ((n_rows() == 0) && (n_cols() == 0));
  }



  inline void
  SparsityPattern::add(const size_type i, const size_type j)
  {
    add_entries(i, &j, &j + 1);
  }



  template <typename ForwardIterator>
  inline void
  SparsityPattern::add_entries(const size_type row,
                               ForwardIterator begin,
                               ForwardIterator end,
                               const bool  /*indices_are_sorted*/ )
  {
    if (begin == end)
      return;

    // verify that the size of the data type Trilinos expects matches that the
    // iterator points to. we allow for some slippage between signed and
    // unsigned and only compare that they are both either 32 or 64 bit. to
    // write this test properly, not that we cannot compare the size of
    // '*begin' because 'begin' may be an iterator and '*begin' may be an
    // accessor class. consequently, we need to somehow get an actual value
    // from it which we can by evaluating an expression such as when
    // multiplying the value produced by 2
    Assert(sizeof(TrilinosWrappers::types::int_type) == sizeof((*begin) * 2),
           ExcNotImplemented());

    TrilinosWrappers::types::int_type *col_index_ptr =
      reinterpret_cast<TrilinosWrappers::types::int_type *>(
        const_cast<typename std::decay<decltype(*begin)>::type *>(&*begin));
    // Check at least for the first index that the conversion actually works
    AssertDimension(*col_index_ptr, *begin);
    TrilinosWrappers::types::int_type trilinos_row_index = row;
    const int                         n_cols = static_cast<int>(end - begin);

    int ierr;
    if (row_is_stored_locally(row))
      ierr =
        graph->InsertGlobalIndices(trilinos_row_index, n_cols, col_index_ptr);
    else if (nonlocal_graph.get() != nullptr)
      {
        // this is the case when we have explicitly set the off-processor rows
        // and want to create a separate matrix object for them (to retain
        // thread-safety)
        Assert(nonlocal_graph->RowMap().LID(
                 static_cast<TrilinosWrappers::types::int_type>(row)) != -1,
               ExcMessage("Attempted to write into off-processor matrix row "
                          "that has not be specified as being writable upon "
                          "initialization"));
        ierr = nonlocal_graph->InsertGlobalIndices(trilinos_row_index,
                                                   n_cols,
                                                   col_index_ptr);
      }
    else
      ierr = graph->InsertGlobalIndices(1,
                                        &trilinos_row_index,
                                        n_cols,
                                        col_index_ptr);

    AssertThrow(ierr >= 0, ExcTrilinosError(ierr));
  }



  inline const Epetra_FECrsGraph &
  SparsityPattern::trilinos_sparsity_pattern() const
  {
    return *graph;
  }



  inline IndexSet
  SparsityPattern::locally_owned_domain_indices() const
  {
    return IndexSet(graph->DomainMap());
  }



  inline IndexSet
  SparsityPattern::locally_owned_range_indices() const
  {
    return IndexSet(graph->RangeMap());
  }

#    endif // DOXYGEN
} // namespace TrilinosWrappers


DEAL_II_NAMESPACE_CLOSE


#  endif // DEAL_II_WITH_TRILINOS


 /*--------------------   trilinos_sparsity_pattern.h     --------------------*/ 

#endif
 /*--------------------   trilinos_sparsity_pattern.h     --------------------*/ 


