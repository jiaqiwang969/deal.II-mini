//include/deal.II-translator/base/table_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2002 - 2021 by the deal.II authors
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

#ifndef dealii_table_h
#define dealii_table_h

#include <deal.II/base/config.h>

#include <deal.II/base/aligned_vector.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/linear_index_iterator.h>
#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/base/table_indices.h>

#include <algorithm>
#include <cstddef>

DEAL_II_NAMESPACE_OPEN

// forward declaration
#ifndef DOXYGEN
template <int N, typename T>
class TableBase;
template <int N, typename T>
class Table;
template <typename T>
class TransposeTable;
template <typename T>
class Table<1, T>;
template <typename T>
class Table<2, T>;
template <typename T>
class Table<3, T>;
template <typename T>
class Table<4, T>;
template <typename T>
class Table<5, T>;
template <typename T>
class Table<6, T>;
#endif



namespace internal
{
  /**
   * @internal
   * 有一个命名空间，我们在其中声明一些类，用于使用<tt>operator[]<tt>访问表的元素。这些是相当有技术含量的，因为它们必须以递归的方式进行工作（由于索引的数量不知道，如果我们访问一个对象，我们必须返回一个迭代器到下一个较低维度的对象，直到我们在最低层，实际上可以返回一个对存储数据类型本身的引用）。
   * 这是非常技术性的，你通常根本不会想看这些类，除了可能出于教育的原因。
   * 这里的任何一个类都没有你应该在你的程序中明确使用的接口（当然，除了通过用<tt>operator[]</tt>访问表的元素，生成这个命名空间的类型的临时对象）。
   *
   */
  namespace TableBaseAccessors
  {
    /**
     * @internal
     * 有一个类，根据其模板参数，声明一些嵌套的别名。一般的模板什么都不声明，但有更多有用的特殊化，关于最后一个参数表示表的常数，访问器对象将在这个命名空间中生成。
     *
     */
    template <int N, typename T, bool Constness>
    struct Types
    {};

    /**
     * @internal
     * 拥有一个声明一些嵌套别名的类，这取决于它的模板参数。对常量对象的访问器进行专业化处理。
     *
     */
    template <int N, typename T>
    struct Types<N, T, true>
    {
      using value_type = const T;
      using TableType  = const TableBase<N, T>;

      using iterator       = typename AlignedVector<T>::const_iterator;
      using const_iterator = typename AlignedVector<T>::const_iterator;

      using reference       = typename AlignedVector<T>::const_reference;
      using const_reference = typename AlignedVector<T>::const_reference;
    };

    /**
     * @internal
     * 拥有一个宣布一些嵌套别名的类，这取决于它的模板参数。对非常量对象的访问器的特殊化。
     *
     */
    template <int N, typename T>
    struct Types<N, T, false>
    {
      using value_type = T;
      using TableType  = TableBase<N, T>;

      using iterator       = typename AlignedVector<T>::iterator;
      using const_iterator = typename AlignedVector<T>::const_iterator;

      using reference       = typename AlignedVector<T>::reference;
      using const_reference = typename AlignedVector<T>::const_reference;
    };


    /**
     * @internal
     * 作为<tt>Table<N,T></tt>类型表的子对象访问器的类。模板参数<tt>C</tt>可以是true或false，并表示所处理的对象是否是常量（即只有当值为false时才允许写访问）。
     * 因为对于<tt>N</tt>索引，应用<tt>operator[]</tt>的效果是获得对<tt>N-1</tt>索引的访问，我们必须递归地实现这些访问器类，当我们只剩下一个索引时就停止。对于后一种情况，下面声明了这个类的特殊化，在这里调用<tt>operator[]</tt>可以访问表实际存储的对象。在给索引运算符的值中，需要检查它是否在其范围内，为此我们需要知道我们目前实际访问的是表的哪个索引。这是通过模板参数<tt>P</tt>来实现的：它表明还有多少个索引。对于一个向量，<tt>P</tt>可能只有一个（然后使用下面的特殊化）。对于一个表来说，这个值可能是两个，当使用<tt>operator[]</tt>时，会出现一个<tt>P=1</tt>的对象。
     * <tt>P</tt>的值也被用来确定stride：这个对象存储了一个指针，表示它可能访问的对象范围的开始。当我们在这个对象上应用<tt>operator[]</tt>时，产生的新访问器只能访问这些元素的一个子集，为了知道哪个子集，我们需要知道表的尺寸和现在的索引，这个索引由<tt>P</tt>表示。
     * 正如对整个命名空间所说的那样，你通常不需要直接处理这些类，也不应该试图直接使用它们的接口，因为它可能在没有通知的情况下发生变化。事实上，由于构造函数是私有的，你甚至不能生成这个类的对象，因为它们只被认为是访问表类元素的临时性的，而不是作为函数的参数来传递的，等等。
     *
     */
    template <int N, typename T, bool C, unsigned int P>
    class Accessor
    {
    public:
      using TableType = typename Types<N, T, C>::TableType;

      using iterator       = typename Types<N, T, C>::iterator;
      using const_iterator = typename Types<N, T, C>::const_iterator;

      using size_type       = size_t;
      using difference_type = ptrdiff_t;

    private:
      /**
       * 构造函数。取一个指向表对象的指针来了解各个维度的大小，以及一个指向我们可能访问的数据子集的指针。
       *
       */
      Accessor(const TableType &table, const iterator data);

    public:
      /**
       * 复制构造函数。这个构造函数是公开的，这样就可以把子表作为参数传递给函数，就像
       * <code>f(table[i])</code>  中的那样。
       * 如果访问器的存储时间比它所指向的表长，使用这个构造函数是有风险的。不要这样做。
       *
       */
      Accessor(const Accessor &a);

      /**
       * 索引操作符。执行一个范围检查。
       *
       */
      Accessor<N, T, C, P - 1> operator[](const size_type i) const;

      /**
       * 范围检查的异常。不要使用全局异常，因为这样我们可以输出哪个索引是错误的。
       *
       */
      DeclException3(ExcIndexRange,
                     size_type,
                     size_type,
                     size_type,
                     << "Index " << N - P + 1 << "has a value of " << arg1
                     << " but needs to be in the range [" << arg2 << "," << arg3
                     << "[.");

    private:
      /**
       * 存储给构造函数的数据。这个类没有非常量的成员函数，所以没有理由不把这些元素变成常量。
       *
       */
      const TableType &table;
      const iterator   data;

      // declare some other classes
      // as friends. make sure to
      // work around bugs in some
      // compilers
      template <int N1, typename T1>
      friend class dealii::Table;
      template <int N1, typename T1, bool C1, unsigned int P1>
      friend class Accessor;
      friend class dealii::Table<N, T>;
      friend class Accessor<N, T, C, P + 1>;
    };



    /**
     * @internal
     * 表格的访问器类。这是最后一个索引的特殊化，它实际上允许对表的元素进行访问，而不是递归地返回进一步子集的访问对象。
     * 这个特殊化与一般模板的情况相同；更多信息请看那里。
     *
     */
    template <int N, typename T, bool C>
    class Accessor<N, T, C, 1>
    {
    public:
      /**
       * 对这一行的元素进行常量和非常量迭代器类型的类型化，以及标准库算法通常需要的所有其他类型。
       *
       */
      using value_type = typename Types<N, T, C>::value_type;

      using iterator       = typename Types<N, T, C>::iterator;
      using const_iterator = typename Types<N, T, C>::const_iterator;

      using reference       = typename Types<N, T, C>::reference;
      using const_reference = typename Types<N, T, C>::const_reference;

      using size_type       = size_t;
      using difference_type = ptrdiff_t;

      /**
       * 从上面的switch类中导入一个别名。
       *
       */
      using TableType = typename Types<N, T, C>::TableType;

    private:
      /**
       * 构造函数。取一个迭代器到表对象，以了解各个维度的大小，以及一个迭代器到我们可能访问的数据子集（在这个特殊的例子中，只有一行）。
       * 构造函数是私有的，以防止你身边有这样的对象。创建这种对象的唯一方法是通过<tt>Table</tt>类，它只生成临时对象。
       * 这保证了访问器对象比母对象更早退出范围，避免了数据一致性的问题。
       *
       */
      Accessor(const TableType &table, const iterator data);

    public:
      /**
       * 复制构造函数。这个构造函数是公开的，这样就可以像
       * <code>f(table[i])</code>
       * 中那样，将子表作为参数传递给函数。
       * 如果访问器的存储时间比它所指向的表长，使用这个构造函数是有风险的。不要这样做。
       *
       */
      Accessor(const Accessor &a);

      /**
       * 索引操作符。执行一个范围检查。
       *
       */
      reference operator[](const size_type) const;

      /**
       * 返回一行的长度，即表对象的最后一个索引所对应的元素数。
       *
       */
      size_type
      size() const;

      /**
       * 返回该行的第一个元素的迭代器。
       *
       */
      iterator
      begin() const;

      /**
       * 返回一个迭代器到超过此行终点的元素。
       *
       */
      iterator
      end() const;

    private:
      /**
       * 存储给构造函数的数据。这个类没有非常量的成员函数，所以没有理由不使这些元素成为常量。
       *
       */
      const TableType &table;
      const iterator   data;

      // declare some other classes
      // as friends. make sure to
      // work around bugs in some
      // compilers
      template <int N1, typename T1>
      friend class dealii::Table;
      template <int N1, typename T1, bool C1, unsigned int P1>
      friend class Accessor;
      friend class dealii::Table<2, T>;
      friend class Accessor<N, T, C, 2>;
    };
  } // namespace TableBaseAccessors

} // namespace internal



/**
 * 一个持有多维数组的模板类型的对象的类。如果表示维数的模板参数是一个，那么这个类或多或少代表了一个向量；如果是两个，那么它就是一个矩阵；以此类推。
 * 这个类特别取代了对高维数组的尝试，如
 * <tt>std::vector<std::vector<T>></tt>,
 * 或甚至更高的嵌套结构。这些结构体的缺点是很难初始化，最重要的是，如果一个矩阵或高维表的所有行都有相同的大小（这是通常的情况），那么它们的效率就非常低，因为此时每一行的内存都是独立分配的，既浪费时间又浪费内存。这可以通过为整个对象只分配一大块内存而变得更有效率，这就是当前类的做法。
 *
 *  <h3>Comparison with the Tensor class</h3>
 * 在某种程度上，这个类类似于张量类，因为它对维数进行模板化。然而，有两个主要区别。第一是Tensor类只存储数值（如<tt>double</tt>s），而Table类存储任意的对象。第二是张量类在每个维度都有固定的尺寸，也是作为模板参数给出的，而这个类可以在每个维度处理任意的、不同的尺寸。
 * 这有两个后果。首先，由于在编译时不知道大小，它必须进行显式内存分配。其次，各个元素的布局在编译时并不知道，所以访问速度比张量类慢，在张量类中，元素的数量和它们的位置在编译时就已经知道了，而且编译器可以利用这些知识进行优化（例如在展开循环的时候）。另一方面，这个类当然更灵活，例如，当你想要一个二维表，其行数等于单元格的自由度数，列数等于正交点的数量。这两个数字可能只有在运行时才知道，所以这里需要一个灵活的表格。此外，你可能想存储，例如，形状函数的梯度，所以数据类型不是单一的标量值，而是张量本身。
 *
 *  <h3>Dealing with large data sets</h3>
 * 表类（派生自该类）经常被用来存储大型数据表。在
 * step-53 中给出了一个适度的例子，我们存储了一个 $380
 * \times 220$
 * 非洲地区的地理海拔数据的表，这个数据需要大约670
 * kB，如果内存；然而，存储三维或更多维度数据的表（例如，关于地球内部的密度、压力和温度的信息，在`（纬度、经度、深度）`点的规则网格上）可以轻松达到数百兆字节或更多。这些表格通常被提供给诸如InterpolatedTensorProductGridData或InterpolatedUniformGridData等类。
 * 如果你需要在单处理器（或多线程）作业中加载这样的表，那么你对这些表的大小无能为力。该表只需适合内存即可。但是，如果你的程序是通过MPI并行化的，那么典型的第一种实现方式是在每个进程上创建一个表对象，并通过从文件中读取数据在每个MPI进程上填充它。这从两个方面来说是低效的。
 *
 *
 *
 * - 你会有很多进程同时试图从同一个文件中读取数据。
 *
 *
 * - 在大多数情况下，每个进程上存储的数据都是一样的，虽然每个进程都需要能够从表中读取数据，但没有必要每个进程都存储自己的表。碰巧位于同一台机器上的所有MPI进程不妨只存储一个副本，并通过[共享内存](https://en.wikipedia.org/wiki/Shared_memory)使其彼此可用；在这种模式下，每台机器只有一个MPI进程需要存储数据，然后所有其他进程可以访问它。
 * 这两种用例都是由内部基于
 * AlignedVector::replicate_across_communicator(). 的
 * TableBase::replicate_across_communicator()
 * 函数实现的，该函数允许像下面这样的工作流程，我们让那个等级为0的MPI进程负责读取数据（但它也可以是任何其他
 * "根等级"）。
 *
 * @code
 *  const unsigned int N=..., M=...;     // table sizes, assumed known
 *  Table<2,double>    data_table;
 *  const unsigned int root_rank = 0;
 *
 *  if (Utilities::MPI::this_mpi_process(mpi_communicator) == root_rank)
 *  {
 *    data_table.resize (N,M);
 *
 *    std::ifstream input_file ("data_file.dat");
 *    ...;                               // read the data from the file
 *  }
 *
 *  // Now distribute to all processes
 *  data_table.replicate_across_communicator (mpi_communicator, root_rank);
 * @endcode
 *
 * 这段代码中的最后一个调用确保数据在所有非根进程中可用，方法是在其他进程的内存空间中重新创建一个表的副本，或者，如果可能的话，在共享内存中为位于MPI作业使用的每台机器上的所有进程创建一次副本。
 *
 *
 * @ingroup data
 *
 */
template <int N, typename T>
class TableBase : public Subscriptor
{
public:
  using value_type = T;

  /**
   * 用来计算这个容器中的元素数量的整数类型。
   *
   */
  using size_type = typename AlignedVector<T>::size_type;


  /**
   * 默认构造函数。将所有尺寸设置为零。
   *
   */
  TableBase() = default;

  /**
   * 构造函数。用每个索引组件中的给定尺寸初始化数组。
   *
   */
  explicit TableBase(const TableIndices<N> &sizes);

  /**
   * 构造函数。在每个索引组件中用给定的尺寸初始化数组，然后通过调用fill(
   * entries,C_style_indexing)，用第二个和第三个参数初始化表格中的元素。
   *
   */
  template <typename InputIterator>
  TableBase(const TableIndices<N> &sizes,
            InputIterator          entries,
            const bool             C_style_indexing = true);

  /**
   * 拷贝构造函数。执行一个深度拷贝。
   *
   */
  TableBase(const TableBase<N, T> &src);

  /**
   * 拷贝构造函数。从存储其他数据类型的表对象中执行深度拷贝。
   *
   */
  template <typename T2>
  TableBase(const TableBase<N, T2> &src);

  /**
   * 移动构造函数。转移另一个表的内容。
   *
   */
  TableBase(TableBase<N, T> &&src) noexcept;

  /**
   * 解除构造函数。释放分配的内存。
   *
   */
  ~TableBase() override = default;

  /**
   * 赋值运算符。将<tt>src</tt>的所有元素复制到矩阵中。
   * 如果需要，大小会被调整。
   * 我们不能使用其他的、模板化的版本，因为如果我们不声明这个版本，编译器会很高兴地生成一个预定义的复制操作符，这不是我们想要的。
   *
   */
  TableBase<N, T> &
  operator=(const TableBase<N, T> &src);

  /**
   * 拷贝操作符。将<tt>src</tt>的所有元素复制到数组中。如果需要的话，大小会被调整。
   * 这个函数要求<tt>T2</tt>的类型可以转换为<tt>T</tt>。
   *
   */
  template <typename T2>
  TableBase<N, T> &
  operator=(const TableBase<N, T2> &src);

  /**
   * 移动赋值运算符。将<tt>src</tt>的所有元素转移到表中。
   *
   */
  TableBase<N, T> &
  operator=(TableBase<N, T> &&src) noexcept;

  /**
   * 测试两个表的相等。
   *
   */
  bool
  operator==(const TableBase<N, T> &T2) const;

  /**
   * 将所有条目设置为默认值（即用默认构造的对象将它们复制过来）。但不要改变表的大小。
   *
   */
  void
  reset_values();

  /**
   * 将此对象的尺寸设置为第一个参数中给出的尺寸，并为表项分配所需的内存以适应这些尺寸。如果
   * @p omit_default_initialization 被设置为 @p false,
   * ，表的所有元素都被设置为元素类型的默认构造对象。否则，内存将处于未初始化或其他未定义的状态。
   *
   */
  void
  reinit(const TableIndices<N> &new_size,
         const bool             omit_default_initialization = false);

  /**
   * 表的大小在<tt>i</tt>方向。
   *
   */
  size_type
  size(const unsigned int i) const;

  /**
   * 返回这个对象在每个方向上的大小。
   *
   */
  const TableIndices<N> &
  size() const;

  /**
   * 返回存储在此对象中的元素数量，它是每个维度上的扩展量的乘积。
   *
   */
  size_type
  n_elements() const;

  /**
   * 返回该对象是否为空，即其中一个方向为零。
   * 这等同于<tt>n_elements()==0</tt>。
   *
   */
  bool
  empty() const;

  /**
   * 通过解引用给定的前向迭代器（例如，可以是一个指向数组第一个元素的指针，或者一个插入
   * std::istream_iterator).
   * 第二个参数表示所指向的元素的排列方式是对应于最后一个索引运行最快还是最慢。默认情况下，使用C风格的索引，即最后一个索引运行最快（与Fortran风格相反，当遍历多维数组时，第一个索引运行最快。例如，如果你试图填充一个Table<2,T>类型的对象，那么用第二个参数的默认值来调用这个函数将导致等同于做
   * @code
   * Table<2,T> t;
   * for (unsigned int i=0; i<t.size(0); ++i)
   *   for (unsigned int j=0; j<t.size(1); ++j)
   *     t[i][j] =entries++;
   * @endcode
   * 另一方面，如果这个函数的第二个参数是假的，那么这将导致以下形式的代码。
   * @code
   * Table<2,T> t;
   * for (unsigned int j=0; j<t.size(1); ++j)
   *   for (unsigned int i=0; i<t.size(0); ++i)
   *     t[i][j] =entries++;
   * @endcode
   * 注意我们通过遍历给定的迭代器集合来填充表元素的转换顺序。
   * @param  entries
   * 一个迭代器，用于初始化这个表的元素集。假设迭代器可以被递增和取消引用足够多的次数来填充这个表。
   * @param  C_style_indexing
   * 如果是true，当我们解指输入范围的后续元素时，以最后一个索引最快的速度运行该表的元素。如果为假，则最快改变第一个索引。
   *
   */
  template <typename InputIterator>
  void
  fill(InputIterator entries, const bool C_style_indexing = true);

  /**
   * 用相同的值填充所有的表项。
   *
   */
  void
  fill(const T &value);

  /**
   * 返回一个对指定元素的读写引用。
   *
   */
  typename AlignedVector<T>::reference
  operator()(const TableIndices<N> &indices);

  /**
   * 返回指定元素的值作为只读引用。
   * 我们将请求的值作为一个常量引用而不是按值返回，因为这个对象可能持有的数据类型可能很大，而且我们在这里不知道复制是否昂贵。
   *
   */
  typename AlignedVector<T>::const_reference
  operator()(const TableIndices<N> &indices) const;

  /**
   * 这个函数在MPI通信器的所有进程中复制由 @p root_process
   * 指示的进程上发现的状态。在 @p root_process
   * 以外的任何进程中发现的当前状态都会在这个进程中丢失。我们可以想象这个操作就像从根进程到所有其他进程对
   * Utilities::MPI::broadcast()
   * 的调用，尽管在实践中这个函数可能试图将数据移动到每个承载MPI进程的机器上的共享内存区域，然后让这个机器上的所有MPI进程访问这个共享内存区域，而不是保留自己的副本。请看这个类的一般文档中的代码例子。
   * 这个函数的意图是将大的数组从一个进程快速交换给其他进程，而不是在所有进程上计算或创建它。这特别适用于从磁盘加载的数据
   *
   * - 说，大的数据表
   *
   * 比起让每个进程自己从磁盘上读取数据，通过读取一次，然后在MPI宇宙中的所有进程中分发，更容易处理。
   * 具体来说，共享内存区域的使用允许在MPI宇宙中每个多核机器上只复制一次数据，而不是为每个MPI进程复制一次数据。如果今天的机器上的数据很大，每个共享内存空间可以很容易地容纳几十个MPI进程，这就可以节省大量内存。
   * 这个功能并不意味着保持不同进程的数据同步，就像
   * parallel::distributed::Vector
   * 和其他矢量类所做的那样，存在一个由每个进程拥有的矢量的某些元素的概念，可能还有从其拥有的进程镜像到其他进程的幽灵元素。相反，当前对象的元素被简单地复制到其他进程中，把这个操作看作是在所有进程中创建一组`const`AlignedVector对象，在复制操作之后不应该再被改变，这是确保向量在所有进程中保持一致的唯一方法。这尤其是因为共享内存区域的使用，在一个MPI进程上对一个向量元素的任何修改也可能导致对其他进程上可见元素的修改，假设它们位于一个共享内存节点内。
   * @note
   * 在MPI进程之间使用共享内存需要检测的MPI安装支持必要的操作。
   * 这对于MPI 3.0和更高版本来说是这样的。
   * @note  这个功能并不便宜。它需要创建所提供 @p
   * communicator
   * 对象的子通信器，这通常是一个昂贵的操作。同样地，共享内存空间的生成也不是一个便宜的操作。因此，当目标是在进程之间共享大的只读数据表时，这个功能主要是有意义的；例子是在启动时加载数据表，然后在程序的运行时间内使用。
   * 在这种情况下，运行这个函数的启动成本可以随着时间的推移而摊销，而且在具有大核心数的机器上，许多MPI进程在同一台机器上运行时，不必在每个进程上存储表所带来的潜在内存节省可能是相当大的。
   * @note  这个函数只有在数据类型`T`是 "自足
   * "的情况下才有意义，也就是说，它的所有信息都存储在其成员变量中，并且没有一个成员变量是指向内存的其他部分的指针。这是因为如果一个类型`T`确实有指向内存其他部分的指针，那么将`T`移到共享内存空间不会导致其他进程访问该对象用其成员变量指针指向的数据。这些数据仍然只存在于一个进程中，并且通常在其他进程无法访问的内存区域。
   * 因此，这个函数的通常使用情况是共享简单对象的数组，如`double's或`int's。
   * @note
   * 调用该函数后，不同MPI进程的对象共享一个共同的状态。这意味着某些操作变得
   * "集体"，即必须在所有参与的处理器上同时调用。特别是，你不能再在一个MPI进程上调用resize()、reserve()或clear()。
   *
   * - 你必须在所有进程上同时这样做，因为它们必须为这些操作进行通信。如果你不这样做，你很可能会得到一个死锁，可能很难调试。推而广之，这个只集体调整大小的规则也延伸到这个函数本身。你不能连续调用它两次，因为这意味着首先除了`root_process'以外的所有进程都要扔掉他们的数据，这不是一个集体操作。一般来说，这些关于可以做什么和不可以做什么的限制，暗示了上面评论的正确性。你应该把一个当前函数被调用的AlignedVector视为`const'，在调用析构器之前，不能对其进行进一步的操作。
   *
   */
  void
  replicate_across_communicator(const MPI_Comm &   communicator,
                                const unsigned int root_process);

  /**
   * 交换这个表和另一个表的内容  @p v.
   * 人们可以用一个临时变量和复制过来的数据元素来完成这个操作，但是这个函数明显更有效率，因为它只交换了两个向量的数据指针，因此不需要分配临时存储和移动数据。
   * 这个函数类似于所有C++标准容器的 @p swap
   * 函数。此外，还有一个全局函数<tt>swap(u,v)</tt>，它简单地调用<tt>u.swap(v)</tt>，同样与标准函数相类似。
   *
   */
  void
  swap(TableBase<N, T> &v);

  /**
   * 确定这个对象的内存消耗（以字节为单位）的估计值。
   *
   */
  std::size_t
  memory_consumption() const;

  /**
   * 使用[BOOST序列化库](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html)将此对象的数据写入或读出到一个流中，以便进行序列化。
   *
   */
  template <class Archive>
  void
  serialize(Archive &ar, const unsigned int version);

protected:
  /**
   * 返回指定元素在一个接一个存储的元素阵列中的位置。这个函数不做索引检查。
   *
   */
  size_type
  position(const TableIndices<N> &indices) const;

  /**
   * 返回一个对指定元素的读写引用。
   * 这个函数不做边界检查，只在内部和已经检查过的函数中使用。
   *
   */
  typename AlignedVector<T>::reference
  el(const TableIndices<N> &indices);

  /**
   * 返回指定元素的值，作为一个只读引用。
   * 这个函数不做边界检查，只在内部和已经检查过的函数中使用。
   * 我们将请求的值作为一个常数引用而不是按值返回，因为这个对象可能持有的数据类型可能很大，而且我们在这里不知道复制是否昂贵。
   *
   */
  typename AlignedVector<T>::const_reference
  el(const TableIndices<N> &indices) const;

protected:
  /**
   * 组件-数组。
   *
   */
  AlignedVector<T> values;

  /**
   * 表的每个方向上的大小。
   *
   */
  TableIndices<N> table_size;

  // Make all other table classes friends.
  template <int, typename>
  friend class TableBase;
};


/**
 * 一个代表具有任意但固定索引数的表的类。这个一般的模板在TableBase类提供的功能之外，还实现了一些额外的功能，例如索引函数采取正确的参数数量，等等。
 * 与其说是这个通用模板，不如说是在这个类的部分特殊化中实现了这些功能，有固定数量的维数。更多信息请见那里，以及基类的文档中。
 *
 *
 * @ingroup data
 *
 *
 */
template <int N, typename T>
class Table;


/**
 * 一个代表一维表的类，也就是一个类似矢量的类。这个类的大部分接口是在TableBase基类中实现的。关于这个类的原理和接口的概要，见那里。
 *
 *
 * @ingroup data
 *
 *
 */
template <typename T>
class Table<1, T> : public TableBase<1, T>
{
public:
  /**
   * 用来计算这个容器中的元素数量的整数类型。
   *
   */
  using size_type = typename TableBase<1, T>::size_type;

  /**
   * 默认构造函数。将所有尺寸设置为零。
   *
   */
  Table() = default;

  /**
   * 构造函数。将给定的维度传给基类。
   *
   */
  explicit Table(const size_type size);

  /**
   * 构造函数。创建一个具有给定尺寸的表，并从一组迭代器中初始化它。
   * 这个函数完全等同于创建一个给定尺寸的表 <code>t</code>
   * ，然后调用
   * @code
   * t.fill (entries, C_style_indexing);
   * @endcode
   * 上，使用 TableBase::fill()
   * 函数，其中的参数会有更详细的解释。但问题是，这只有在运行构造函数后可以改变表的情况下才有可能，而调用当前的构造函数可以立即确定对象的大小和初始化，这样就可以将其标记为常量。
   * 使用这个构造函数，你可以做这样的事情。
   * @code
   * const double values[] = { 1, 2, 3 };
   * const Table<1,double> t(3, entries, true);
   * @endcode
   * 你也可以使用输入迭代器，从文件中直接初始化一个表。
   * @code
   * std::ifstream input ("myfile");
   * const Table<1,double> t(3,
   *                         std::istream_iterator<double>(input),
   *                         true);
   * @endcode
   * @param  size 这个一维表的大小。    @param  entries
   * 一组元素的迭代器，用于初始化这个表。假设迭代器可以被递增和取消引用足够多的次数来填充这个表。
   * @param  C_style_indexing
   * 如果是true，当我们解指输入范围的后续元素时，以最后一个索引变化最快的方式运行该表的元素。如果为假，则最快改变第一个索引。
   *
   */
  template <typename InputIterator>
  Table(const size_type size,
        InputIterator   entries,
        const bool      C_style_indexing = true);

  /**
   * 访问操作符。由于这是一个一维对象，这只是访问所请求的数据元素。返回一个只读的引用。
   *
   */
  typename AlignedVector<T>::const_reference
  operator[](const size_type i) const;

  /**
   * 访问操作符。由于这是一个一维的对象，这只是访问所请求的数据元素。返回一个读-写引用。
   *
   */
  typename AlignedVector<T>::reference operator[](const size_type i);

  /**
   * 访问操作符。由于这是一个一维的对象，这只是访问所请求的数据元素。返回一个只读的引用。
   *
   */
  typename AlignedVector<T>::const_reference
  operator()(const size_type i) const;

  /**
   * 访问操作符。由于这是一个一维的对象，这只是访问所请求的数据元素。返回一个读-写引用。
   *
   */
  typename AlignedVector<T>::reference
  operator()(const size_type i);

  /**
   * 使基类中的`operator()`的变化可用。
   *
   */
  using TableBase<1, T>::operator();
};



/**
 * 一个用于Table<2,
 * T>和TransposeTable的迭代器和访问器的命名空间。这些类有特殊的访问器（也就是说，与Table<3,
 * T>相比），因为它们有一个类似矩阵的结构；也就是说，访问器也提供行和列的信息，并被设计为与SparseMatrix和SparsityPattern迭代器类兼容。
 *
 *
 */
namespace MatrixTableIterators
{
  /**
   * @brief Enumeration describing the storage order (i.e., the in-memory
   * 表类的布局）。
   *
   */
  enum class Storage
  {
    /**
     * 数据以行为主（即C风格）的顺序组织。
     *
     */
    row_major,

    /**
     * 数据以列为主（即Fortran风格）的顺序组织。
     *
     */
    column_major
  };

  // Forward declaration of the iterator class.
  template <typename TableType, bool Constness, Storage storage_order>
  class Iterator;

  /**
   * @brief %Accessor class template. This class is partially specialized for
   * <code>Constness</code> 的两个值。
   *
   */
  template <typename TableType, bool Constness, Storage storage_order>
  class Accessor;

  /**
   * @brief %Accessor base class for Table<2, T> and TransposeTable.
   * 该类与LinearIndexIterator中描述的对%Accessor的要求兼容。关于迭代器和访问器之间的分割描述，请参见该类的文档。
   * @tparam  TableType %Table的类型，例如，Table<2,
   * T>或TransposeTable。      @tparam  Constness
   * 这个对象是否存储一个常量指针，并可以修改所提供的表对象。
   * @tparam  storage_order 底层表的存储方案，例如，
   * Storage::row_major  为Table<2, T>。
   *
   */
  template <typename TableType, bool Constness, Storage storage_order>
  class AccessorBase
  {
  public:
    /**
     * 存储指向表的指针的类型。
     *
     */
    using container_pointer_type = typename std::
      conditional<Constness, const TableType *, TableType *>::type;

    /**
     * 底层容器的值类型。
     *
     */
    using value_type = typename TableType::value_type;

    /**
     * 表的行和列索引的数字类型。
     *
     */
    using size_type = typename TableType::size_type;

    /**
     * 默认构造函数。
     *
     */
    AccessorBase();

    /**
     * 设置终端迭代器的构造函数。
     *
     */
    AccessorBase(const container_pointer_type table);

    /**
     * 从非const Accessor复制构造函数。
     *
     */
    AccessorBase(const AccessorBase<TableType, false, storage_order> &);

    /**
     * 接受数组索引的构造函数。
     *
     */
    AccessorBase(const container_pointer_type table,
                 const std::ptrdiff_t         linear_index);

    /**
     * 比较运算符。
     *
     */
    template <bool OtherConstness>
    friend bool
    operator==(
      const AccessorBase<TableType, Constness, storage_order> &     left,
      const AccessorBase<TableType, OtherConstness, storage_order> &right)
    {
      return left.container == right.container &&
             left.linear_index == right.linear_index;
    }

    /**
     * 获取该访问器所代表的元素的值的常数引用。
     *
     */
    const value_type &
    value() const;

    /**
     * 转换操作符，返回元素的常数引用。
     *
     */
    operator const value_type &() const;

    /**
     * 返回当前条目的行。
     *
     */
    size_type
    row() const;

    /**
     * 返回当前条目的列。
     *
     */
    size_type
    column() const;

  protected:
    /**
     * 指向表格的指针。
     *
     */
    container_pointer_type container;

    /**
     * 当前的索引。
     *
     */
    std::ptrdiff_t linear_index;

    /**
     * 检查 <code>linear_index</code>
     * 是否对应于表实际存储的条目（即断言
     * <code>linear_index</code> 为非负数且小于
     * <code>container->size()</code> ）。
     *
     */
    void
    assert_valid_linear_index() const;

    // Make the const version a friend for copying.
    friend class AccessorBase<TableType, true, storage_order>;

    // Make the underlying iterator class a friend.
    friend class LinearIndexIterator<
      Iterator<TableType, Constness, storage_order>,
      Accessor<TableType, Constness, storage_order>>;
  };

  /**
   * @brief %Accessor class offering read-only access to elements of a
   * 表。这与基类相同。
   *
   */
  template <typename TableType, Storage storage_order>
  class Accessor<TableType, true, storage_order>
    : public AccessorBase<TableType, true, storage_order>
  {
  public:
    /**
     * 使用基类的价值类型。
     *
     */
    using value_type =
      typename AccessorBase<TableType, true, storage_order>::value_type;

    /**
     * 使用基类的大小类型。
     *
     */
    using size_type =
      typename AccessorBase<TableType, true, storage_order>::size_type;

    /**
     * 继承基类的构造函数。
     *
     */
    using AccessorBase<TableType, true, storage_order>::AccessorBase;
  };

  /**
   * @brief %Accessor class offering read and write access to the elements of
   * 一个表格。
   *
   */
  template <typename TableType, Storage storage_order>
  class Accessor<TableType, false, storage_order>
    : public AccessorBase<TableType, false, storage_order>
  {
  public:
    /**
     * 使用基类的值类型。
     *
     */
    using value_type =
      typename AccessorBase<TableType, true, storage_order>::value_type;

    /**
     * 使用基类大小类型。
     *
     */
    using size_type =
      typename AccessorBase<TableType, true, storage_order>::size_type;

    /**
     * 继承基类的构造函数。
     *
     */
    using AccessorBase<TableType, false, storage_order>::AccessorBase;

    /**
     * 赋值运算符。这将给当前行和列坐标的表项分配一个新值。
     *
     */
    const Accessor<TableType, false, storage_order> &
    operator=(const value_type &) const;

    /**
     * 移动赋值运算符。在当前的行和列坐标上为表项分配一个新的值。
     *
     */
    const Accessor<TableType, false, storage_order> &
    operator=(value_type &&) const;

    /**
     * 由于我们重载了value()，我们必须明确地使用基类的版本。
     *
     */
    using AccessorBase<TableType, false, storage_order>::value;

    /**
     * 获取该访问器所代表的元素的值的引用。
     *
     */
    value_type &
    value() const;

    /**
     * 返回该元素的引用的转换操作符。
     *
     */
    operator value_type &();
  };

  /**
   * @brief %Iterator class for both matrix-like tables, i.e., Table<2, T> and
   * TransposeTable.       @tparam  TableType
   * %Table的类型，例如，Table<2, T>或TransposeTable。      @tparam
   * Constness 这是否是一个常量迭代器。      @tparam
   * storage_order 底层表的存储方案，例如， Storage::row_major
   * 为Table<2, T>。
   *
   */
  template <typename TableType, bool Constness, Storage storage_order>
  class Iterator
    : public LinearIndexIterator<Iterator<TableType, Constness, storage_order>,
                                 Accessor<TableType, Constness, storage_order>>
  {
  public:
    /**
     * 底层表使用的大小类型。
     *
     */
    using size_type = typename TableType::size_type;

    /**
     * 存储到表的指针的类型。
     *
     */
    using container_pointer_type = typename std::
      conditional<Constness, const TableType *, TableType *>::type;

    /**
     * 来自访问器的构造函数。
     *
     */
    Iterator(const Accessor<TableType, Constness, storage_order> &accessor);

    /**
     * 构造函数。创建一个表的终端迭代器。
     *
     */
    Iterator(const container_pointer_type object);

    /**
     * 为一个特定的表项的构造函数。
     *
     */
    Iterator(const container_pointer_type object,
             const size_type              row,
             const size_type              column);

    /**
     * 从一个非const迭代器复制构造器。
     *
     */
    Iterator(const Iterator<TableType, false, storage_order> &i);

    /**
     * 具有特定线性索引的条目的构造函数。
     *
     */
    Iterator(const container_pointer_type container,
             const std::ptrdiff_t         linear_index);
  };
} // namespace MatrixTableIterators



/**
 * 一个代表二维表的类，即一个对象的矩阵（不一定只有数字）。这个类的大部分接口是在TableBase基类中实现的。关于这个类的原理和接口的概要，请看那里。
 * 这个类也是FullMatrix类的基类，因此它有一些专门针对矩阵及其需求的功能。
 *
 *
 * @ingroup data
 *
 *
 */
template <typename T>
class Table<2, T> : public TableBase<2, T>
{
public:
  /**
   * 用来计算这个容器中的元素数量的整数类型。
   *
   */
  using size_type = typename TableBase<2, T>::size_type;

  /**
   * 用于表内数值的类型定义。
   *
   */
  using value_type = typename AlignedVector<T>::value_type;

  /**
   * 表中引用的类型定义。
   *
   */
  using reference = typename AlignedVector<T>::reference;

  /**
   * 为表中的常数引用提供类型定义。
   *
   */
  using const_reference = typename AlignedVector<T>::const_reference;

  /**
   * 为一个常数迭代器提供类型定义，该迭代器以列为主的顺序遍历表。
   *
   */
  using const_iterator = MatrixTableIterators::
    Iterator<Table<2, T>, true, MatrixTableIterators::Storage::row_major>;

  /**
   * 一个迭代器的类型定义，该迭代器以列的主次顺序遍历表。
   *
   */
  using iterator = MatrixTableIterators::
    Iterator<Table<2, T>, false, MatrixTableIterators::Storage::row_major>;

  /**
   * 默认构造函数。将所有的尺寸设置为零。
   *
   */
  Table() = default;

  /**
   * 构造函数。将给定的尺寸传递给基类。
   *
   */
  Table(const size_type size1, const size_type size2);

  /**
   * 构造函数。创建一个具有给定尺寸的表，并从一组迭代器中初始化它。
   * 这个函数完全等同于创建一个给定尺寸的表 <code>t</code>
   * ，然后调用
   * @code
   * t.fill (entries, C_style_indexing);
   * @endcode
   * 上，使用 TableBase::fill()
   * 函数，其中的参数会有更详细的解释。但问题是，这只有在运行构造函数后可以改变表的情况下才有可能，而调用当前的构造函数可以立即确定对象的大小和初始化，这样就可以将其标记为常量。
   * 使用这个构造函数，你可以做这样的事情。
   * @code
   * const double values[] = { 1, 2, 3, 4, 5, 6 };
   * const Table<2,double> t(2, 3, entries, true);
   * @endcode
   * 你也可以使用输入迭代器，从文件中直接初始化一个表。
   * @code
   * std::ifstream input ("myfile");
   * const Table<2,double> t(2, 3,
   *                         std::istream_iterator<double>(input),
   *                         true);
   * @endcode
   * @param  size1 这个表的第一维的大小。    @param  size2
   * 这个表在第二维中的大小。    @param  entries
   * 一个迭代器，用于初始化该表的元素集。假设迭代器可以被递增和取消引用足够多的次数来填充这个表。
   * @param  C_style_indexing
   * 如果是true，当我们解指输入范围的后续元素时，以最后一个索引变化最快的方式运行该表的元素。如果为假，则最快改变第一个索引。
   *
   */
  template <typename InputIterator>
  Table(const size_type size1,
        const size_type size2,
        InputIterator   entries,
        const bool      C_style_indexing = true);

  /**
   * 重新初始化该对象。这个函数在这里主要是为了与早期的<tt>vector2d</tt>类兼容。通过将参数转换为基类所要求的数据类型而向下传递给基类。
   *
   */
  void
  reinit(const size_type size1,
         const size_type size2,
         const bool      omit_default_initialization = false);

  using TableBase<2, T>::reinit;

  /**
   * 访问操作符。生成一个对象，访问这个二维表的请求行。会进行范围检查。
   * 这个版本的函数只允许读取访问。
   *
   */
  dealii::internal::TableBaseAccessors::Accessor<2, T, true, 1>
  operator[](const size_type i) const;

  /**
   * 访问操作符。生成一个对象，访问这个二维表的请求行。会进行范围检查。
   * 这个版本的函数允许读-写访问。
   *
   */
  dealii::internal::TableBaseAccessors::Accessor<2, T, false, 1>
  operator[](const size_type i);

  /**
   * 通过同时指定所有索引，直接访问表中的一个元素。会进行范围检查。
   * 这个版本的函数只允许读访问。
   *
   */
  typename AlignedVector<T>::const_reference
  operator()(const size_type i, const size_type j) const;

  /**
   * 通过同时指定所有索引，直接访问表中的一个元素。会进行范围检查。
   * 这个版本的函数允许读-写访问。
   *
   */
  typename AlignedVector<T>::reference
  operator()(const size_type i, const size_type j);

  /**
   * 使基类中的`operator()`的变化可用。
   *
   */
  using TableBase<2, T>::operator();

  /**
   * 行的数量。由于我们这里有一个二维的对象，所以这个函数真的只有意义。
   *
   */
  size_type
  n_rows() const;

  /**
   * 列的数量。这个函数真的很有意义，因为我们在这里有一个二维的对象。
   *
   */
  size_type
  n_cols() const;

  /**
   * 返回一个指向第一个条目的迭代器。
   *
   */
  iterator
  begin();

  /**
   * 返回一个指向第一个条目的常数迭代器。
   *
   */
  const_iterator
  begin() const;

  /**
   * 返回一个指向超过最后一个条目的迭代器。
   *
   */
  iterator
  end();

  /**
   * 返回一个常数迭代器，指向最后一个条目之后的一个条目。
   *
   */
  const_iterator
  end() const;

protected:
  /**
   * 返回一个对元素<tt>(i,j)</tt>的读写引用。
   * 这个函数不做边界检查，只在内部和已经检查过的函数中使用。
   * 这些函数在这里主要是为了与这些表类以前对2D数组的实现兼容，当时称为<tt>vector2d</tt>。
   *
   */
  typename AlignedVector<T>::reference
  el(const size_type i, const size_type j);

  /**
   * 返回元素<tt>(i,j)</tt>的值，作为一个只读的引用。
   * 这个函数不做边界检查，只在内部和已经检查过的函数中使用。
   * 我们将请求的值作为常量引用而不是按值返回，因为这个对象可能持有的数据类型可能很大，而我们在这里不知道复制是否昂贵。
   * 这些函数在这里主要是为了与这些表类以前对2D数组的实现兼容，当时称为<tt>vector2d</tt>。
   *
   */
  typename AlignedVector<T>::const_reference
  el(const size_type i, const size_type j) const;

  // Make the AccessorBase class a friend so that it may directly index into
  // the array.
  friend class MatrixTableIterators::
    AccessorBase<Table<2, T>, true, MatrixTableIterators::Storage::row_major>;
  friend class MatrixTableIterators::
    AccessorBase<Table<2, T>, false, MatrixTableIterators::Storage::row_major>;

  // Make the mutable accessor class a friend so that we can write to array
  // entries with iterators.
  friend class MatrixTableIterators::
    Accessor<Table<2, T>, false, MatrixTableIterators::Storage::row_major>;
};



/**
 * 一个代表对象（不一定只有数字）的三维表的类。这个类的大部分接口是在TableBase基类中实现的。关于这个类的原理和接口的概要，见那里。
 *
 *
 * @ingroup data
 *
 *
 */
template <typename T>
class Table<3, T> : public TableBase<3, T>
{
public:
  /**
   * 用来计算这个容器中的元素数量的整数类型。
   *
   */
  using size_type = typename TableBase<3, T>::size_type;

  /**
   * 默认构造函数。将所有尺寸设置为零。
   *
   */
  Table() = default;

  /**
   * 构造函数。将给定的尺寸传递给基类。
   *
   */
  Table(const size_type size1, const size_type size2, const size_type size3);

  /**
   * 构造函数。创建一个具有给定尺寸的表，并从一组迭代器中初始化它。
   * 这个函数完全等同于创建一个给定尺寸的表 <code>t</code>
   * ，然后调用
   * @code
   * t.fill (entries, C_style_indexing);
   * @endcode
   * 上，使用 TableBase::fill()
   * 函数，其中的参数会有更详细的解释。但问题是，这只有在运行构造函数后可以改变表的情况下才有可能，而调用当前的构造函数可以立即确定对象的大小和初始化，这样就可以将其标记为常量。
   * 使用这个构造函数，你可以做这样的事情（这里显示的是一个二维表，但同样适用于当前类）。
   * @code
   * const double values[] = { 1, 2, 3, 4, 5, 6 };
   * const Table<2,double> t(2, 3, entries, true);
   * @endcode
   * 你也可以使用输入迭代器，从文件中直接初始化一个表。
   * @code
   * std::ifstream input ("myfile");
   * const Table<2,double> t(2, 3,
   *                         std::istream_iterator<double>(input),
   *                         true);
   * @endcode
   * @param  size1 这个表在第一维的大小。    @param  size2
   * 该表在第二维中的大小。    @param  size3
   * 此表在第三维中的大小。    @param  entries
   * 一个迭代器，用于初始化这个表的元素集。假设迭代器可以被递增和取消引用足够多的次数来填充这个表。
   * @param  C_style_indexing
   * 如果是true，当我们解指输入范围的后续元素时，以最后一个索引变化最快的方式运行该表的元素。如果为假，则最快改变第一个索引。
   *
   */
  template <typename InputIterator>
  Table(const size_type size1,
        const size_type size2,
        const size_type size3,
        InputIterator   entries,
        const bool      C_style_indexing = true);

  /**
   * 访问操作符。生成一个对象，访问这个三维表的请求的二维子对象。会进行范围检查。
   * 这个版本的函数只允许读取访问。
   *
   */
  dealii::internal::TableBaseAccessors::Accessor<3, T, true, 2>
  operator[](const size_type i) const;

  /**
   * 访问操作符。生成一个对象，访问这个三维表的请求的二维子对象。会进行范围检查。
   * 这个版本的函数允许读写访问。
   *
   */
  dealii::internal::TableBaseAccessors::Accessor<3, T, false, 2>
  operator[](const size_type i);

  /**
   * 通过同时指定所有索引，直接访问表中的一个元素。会进行范围检查。
   * 这个版本的函数只允许读访问。
   *
   */
  typename AlignedVector<T>::const_reference
  operator()(const size_type i, const size_type j, const size_type k) const;


  /**
   * 通过同时指定所有索引，直接访问表中的一个元素。会进行范围检查。
   * 这个版本的函数允许读-写访问。
   *
   */
  typename AlignedVector<T>::reference
  operator()(const size_type i, const size_type j, const size_type k);

  /**
   * 使基类中的`operator()`的变化可用。
   *
   */
  using TableBase<3, T>::operator();
};



/**
 * 一个代表对象的四维表（不一定只有数字）的类。这个类的大部分接口是在TableBase基类中实现的。关于这个类的原理和接口的概要，见那里。
 *
 *
 * @ingroup data
 *
 *
 */
template <typename T>
class Table<4, T> : public TableBase<4, T>
{
public:
  /**
   * 用来计算这个容器中的元素数量的整数类型。
   *
   */
  using size_type = typename TableBase<4, T>::size_type;

  /**
   * 默认构造函数。将所有尺寸设置为零。
   *
   */
  Table() = default;

  /**
   * 构造函数。将给定的尺寸传递给基类。
   *
   */
  Table(const size_type size1,
        const size_type size2,
        const size_type size3,
        const size_type size4);

  /**
   * 访问操作符。生成一个对象，访问这个四维表的请求的三维子对象。会进行范围检查。
   * 这个版本的函数只允许读取访问。
   *
   */
  dealii::internal::TableBaseAccessors::Accessor<4, T, true, 3>
  operator[](const size_type i) const;

  /**
   * 访问操作符。生成一个对象，访问这个四维表的所要求的三维子对象。会进行范围检查。
   * 这个版本的函数允许读-写访问。
   *
   */
  dealii::internal::TableBaseAccessors::Accessor<4, T, false, 3>
  operator[](const size_type i);

  /**
   * 通过同时指定所有索引，直接访问表中的一个元素。会进行范围检查。
   * 这个版本的函数只允许读访问。
   *
   */
  typename AlignedVector<T>::const_reference
  operator()(const size_type i,
             const size_type j,
             const size_type k,
             const size_type l) const;


  /**
   * 通过同时指定所有索引，直接访问表中的一个元素。会进行范围检查。
   * 这个版本的函数允许读-写访问。
   *
   */
  typename AlignedVector<T>::reference
  operator()(const size_type i,
             const size_type j,
             const size_type k,
             const size_type l);

  /**
   * 使基类中的`operator()`的变化可用。
   *
   */
  using TableBase<4, T>::operator();
};



/**
 * 一个代表对象的五维表（不一定只有数字）的类。这个类的大部分接口是在TableBase基类中实现的。关于这个类的原理和接口的概要，见那里。
 *
 *
 * @ingroup data
 *
 *
 */
template <typename T>
class Table<5, T> : public TableBase<5, T>
{
public:
  /**
   * 用来计算这个容器中的元素数量的整数类型。
   *
   */
  using size_type = typename TableBase<5, T>::size_type;


  /**
   * 默认构造函数。将所有尺寸设置为零。
   *
   */
  Table() = default;

  /**
   * 构造函数。将给定的尺寸传递给基类。
   *
   */
  Table(const size_type size1,
        const size_type size2,
        const size_type size3,
        const size_type size4,
        const size_type size5);

  /**
   * 访问操作符。生成一个对象，访问这个五维表的请求的四维子对象。会进行范围检查。
   * 这个版本的函数只允许读取访问。
   *
   */
  dealii::internal::TableBaseAccessors::Accessor<5, T, true, 4>
  operator[](const size_type i) const;

  /**
   * 访问操作符。生成一个对象，访问这个五维表的四维子对象。会进行范围检查。
   * 这个版本的函数允许读写访问。
   *
   */
  dealii::internal::TableBaseAccessors::Accessor<5, T, false, 4>
  operator[](const size_type i);

  /**
   * 通过同时指定所有索引，直接访问表中的一个元素。会进行范围检查。
   * 这个版本的函数只允许读访问。
   *
   */
  typename AlignedVector<T>::const_reference
  operator()(const size_type i,
             const size_type j,
             const size_type k,
             const size_type l,
             const size_type m) const;

  /**
   * 通过同时指定所有索引，直接访问表中的一个元素。会进行范围检查。
   * 这个版本的函数允许读-写访问。
   *
   */
  typename AlignedVector<T>::reference
  operator()(const size_type i,
             const size_type j,
             const size_type k,
             const size_type l,
             const size_type m);

  /**
   * 使基类中的`operator()`的变化可用。
   *
   */
  using TableBase<5, T>::operator();
};



/**
 * 一个代表对象的六维表（不一定只有数字）的类。这个类的大部分接口是在TableBase基类中实现的。关于这个类的原理和接口的概要，见那里。
 *
 *
 * @ingroup data
 *
 *
 */
template <typename T>
class Table<6, T> : public TableBase<6, T>
{
public:
  /**
   * 用来计算这个容器中的元素数量的整数类型。
   *
   */
  using size_type = typename TableBase<6, T>::size_type;

  /**
   * 默认构造函数。将所有尺寸设置为零。
   *
   */
  Table() = default;

  /**
   * 构造函数。将给定的尺寸传递给基类。
   *
   */
  Table(const size_type size1,
        const size_type size2,
        const size_type size3,
        const size_type size4,
        const size_type size5,
        const size_type size6);

  /**
   * 访问操作符。生成一个对象，访问这个六维表的五维子对象的请求。会进行范围检查。
   * 这个版本的函数只允许读取访问。
   *
   */
  dealii::internal::TableBaseAccessors::Accessor<6, T, true, 5>
  operator[](const size_type i) const;

  /**
   * 访问操作符。生成一个对象，访问这个六维表的五维子对象的请求。会进行范围检查。
   * 这个版本的函数允许读写访问。
   *
   */
  dealii::internal::TableBaseAccessors::Accessor<6, T, false, 5>
  operator[](const size_type i);

  /**
   * 通过同时指定所有索引，直接访问表中的一个元素。会进行范围检查。
   * 这个版本的函数只允许读访问。
   *
   */
  typename AlignedVector<T>::const_reference
  operator()(const size_type i,
             const size_type j,
             const size_type k,
             const size_type l,
             const size_type m,
             const size_type n) const;

  /**
   * 通过同时指定所有索引，直接访问表中的一个元素。会进行范围检查。
   * 这个版本的函数允许读-写访问。
   *
   */
  typename AlignedVector<T>::reference
  operator()(const size_type i,
             const size_type j,
             const size_type k,
             const size_type l,
             const size_type m,
             const size_type n);

  /**
   * 使基类中的`operator()`的变化可用。
   *
   */
  using TableBase<6, T>::operator();
};


/**
 * 一个代表对象的七维表（不一定只有数字）的类。这个类的大部分接口是在TableBase基类中实现的。关于这个类的原理和接口的概要，见那里。
 *
 *
 * @ingroup data
 *
 *
 */
template <typename T>
class Table<7, T> : public TableBase<7, T>
{
public:
  /**
   * 用来计算这个容器中的元素数量的整数类型。
   *
   */
  using size_type = typename TableBase<7, T>::size_type;

  /**
   * 默认构造函数。将所有尺寸设置为零。
   *
   */
  Table() = default;

  /**
   * 构造函数。将给定的尺寸传递给基类。
   *
   */
  Table(const size_type size1,
        const size_type size2,
        const size_type size3,
        const size_type size4,
        const size_type size5,
        const size_type size6,
        const size_type size7);

  /**
   * 访问操作符。生成一个对象，访问这个七维表的六维子对象的请求。会进行范围检查。
   * 这个版本的函数只允许读取访问。
   *
   */
  dealii::internal::TableBaseAccessors::Accessor<7, T, true, 6>
  operator[](const size_type i) const;

  /**
   * 访问操作符。生成一个对象，访问这个七维表的六维子对象的请求。会进行范围检查。
   * 这个版本的函数允许读写访问。
   *
   */
  dealii::internal::TableBaseAccessors::Accessor<7, T, false, 6>
  operator[](const size_type i);

  /**
   * 通过同时指定所有索引，直接访问表中的一个元素。会进行范围检查。
   * 这个版本的函数只允许读访问。
   *
   */
  typename AlignedVector<T>::const_reference
  operator()(const size_type i,
             const size_type j,
             const size_type k,
             const size_type l,
             const size_type m,
             const size_type n,
             const size_type o) const;

  /**
   * 通过同时指定所有索引，直接访问表中的一个元素。会进行范围检查。
   * 这个版本的函数允许读-写访问。
   *
   */
  typename AlignedVector<T>::reference
  operator()(const size_type i,
             const size_type j,
             const size_type k,
             const size_type l,
             const size_type m,
             const size_type n,
             const size_type o);

  /**
   * 使基类中的`operator()`的变化可用。
   *
   */
  using TableBase<7, T>::operator();
};


/**
 * 一个代表转置二维表的类，即一个对象（不一定只有数字）的矩阵，采用列首编号（FORTRAN惯例）。因此，唯一真正的区别其实是在存储格式上。
 * 这个类复制了Table<2,T>的功能，但是元素访问和尺寸将用于TableBase中的数据字段的转置排序。
 *
 *
 * @ingroup data
 *
 */
template <typename T>
class TransposeTable : public TableBase<2, T>
{
public:
  /**
   * 用来计算这个容器中的元素数量的整数类型。
   *
   */
  using size_type = typename TableBase<2, T>::size_type;

  /**
   * 用于表内数值的类型化定义。
   *
   */
  using value_type = typename AlignedVector<T>::value_type;

  /**
   * 表中引用的类型定义。
   *
   */
  using reference = typename AlignedVector<T>::reference;

  /**
   * 为表中的常数引用提供类型定义。
   *
   */
  using const_reference = typename AlignedVector<T>::const_reference;

  /**
   * 为一个常数迭代器提供类型定义，该迭代器以列为主的顺序遍历表。
   *
   */
  using const_iterator =
    MatrixTableIterators::Iterator<TransposeTable<T>,
                                   true,
                                   MatrixTableIterators::Storage::column_major>;

  /**
   * 一个迭代器的类型定义，该迭代器以列的主次顺序遍历表。
   *
   */
  using iterator =
    MatrixTableIterators::Iterator<TransposeTable<T>,
                                   false,
                                   MatrixTableIterators::Storage::column_major>;

  /**
   * 默认构造函数。将所有的尺寸设置为零。
   *
   */
  TransposeTable() = default;

  /**
   * 构造函数。将给定的尺寸传递给基类。
   *
   */
  TransposeTable(const size_type size1, const size_type size2);

  /**
   * 重新初始化该对象。这个函数在这里主要是为了与早期的<tt>vector2d</tt>类兼容。通过将参数转换为基类所要求的数据类型而向下传递给基类。
   *
   */
  void
  reinit(const size_type size1,
         const size_type size2,
         const bool      omit_default_initialization = false);

  /**
   * 通过同时指定所有索引，直接访问表中的一个元素。会进行范围检查。
   * 这个版本的函数只允许读取访问。
   *
   */
  const_reference
  operator()(const size_type i, const size_type j) const;

  /**
   * 通过同时指定所有索引，直接访问表中的一个元素。会进行范围检查。
   * 这个版本的函数允许读-写访问。
   *
   */
  reference
  operator()(const size_type i, const size_type j);

  /**
   * 行的数量。由于我们这里有一个二维对象，所以这个函数真的只有意义。
   *
   */
  size_type
  n_rows() const;

  /**
   * 列的数量。这个函数真的很有意义，因为我们在这里有一个二维的对象。
   *
   */
  size_type
  n_cols() const;

  /**
   * 返回一个指向第一个条目的迭代器。
   *
   */
  iterator
  begin();

  /**
   * 返回一个指向第一个条目的常数迭代器。
   *
   */
  const_iterator
  begin() const;

  /**
   * 返回一个指向超过最后一个条目的迭代器。
   *
   */
  iterator
  end();

  /**
   * 返回一个常数迭代器，指向最后一个条目之后的一个条目。
   *
   */
  const_iterator
  end() const;

protected:
  /**
   * 返回一个对元素<tt>(i,j)</tt>的读写引用。
   * 这个函数不做边界检查，只在内部和已经检查过的函数中使用。
   * 这些函数在这里主要是为了与这些表类以前对2D数组的实现兼容，当时称为<tt>vector2d</tt>。
   *
   */
  reference
  el(const size_type i, const size_type j);

  /**
   * 返回元素<tt>(i,j)</tt>的值，作为一个只读的引用。
   * 这个函数不做边界检查，只在内部和已经检查过的函数中使用。
   * 我们将请求的值作为常量引用而不是按值返回，因为这个对象可能持有的数据类型可能很大，而我们在这里不知道复制是否昂贵。
   * 这些函数在这里主要是为了与这些表类以前对2D数组的实现兼容，当时称为<tt>vector2d</tt>。
   *
   */
  const_reference
  el(const size_type i, const size_type j) const;

  // Make the AccessorBase class a friend so that it may directly index into
  // the array.
  friend class MatrixTableIterators::AccessorBase<
    TransposeTable<T>,
    true,
    MatrixTableIterators::Storage::column_major>;
  friend class MatrixTableIterators::AccessorBase<
    TransposeTable<T>,
    false,
    MatrixTableIterators::Storage::column_major>;

  // Make the mutable accessor class a friend so that we can write to array
  // entries with iterators.
  friend class MatrixTableIterators::Accessor<
    TransposeTable<T>,
    false,
    MatrixTableIterators::Storage::column_major>;
};



 /* --------------------- Template and inline functions ---------------- */ 

#ifndef DOXYGEN

template <int N, typename T>
TableBase<N, T>::TableBase(const TableIndices<N> &sizes)
{
  reinit(sizes);
}



template <int N, typename T>
template <typename InputIterator>
TableBase<N, T>::TableBase(const TableIndices<N> &sizes,
                           InputIterator          entries,
                           const bool             C_style_indexing)
{
  reinit(sizes);
  fill(entries, C_style_indexing);
}



template <int N, typename T>
TableBase<N, T>::TableBase(const TableBase<N, T> &src)
  : Subscriptor()
{
  reinit(src.table_size, true);
  values = src.values;
}



template <int N, typename T>
template <typename T2>
TableBase<N, T>::TableBase(const TableBase<N, T2> &src)
{
  reinit(src.table_size);
  if (src.n_elements() != 0)
    std::copy(src.values.begin(), src.values.end(), values.begin());
}



template <int N, typename T>
TableBase<N, T>::TableBase(TableBase<N, T> &&src) noexcept
  : Subscriptor(std::move(src))
  , values(std::move(src.values))
  , table_size(src.table_size)
{
  src.table_size = TableIndices<N>();
}



template <int N, typename T>
template <class Archive>
inline void
TableBase<N, T>::serialize(Archive &ar, const unsigned int)
{
  ar &static_cast<Subscriptor &>(*this);

  ar &values &table_size;
}



namespace internal
{
  namespace TableBaseAccessors
  {
    template <int N, typename T, bool C, unsigned int P>
    inline Accessor<N, T, C, P>::Accessor(const TableType &table,
                                          const iterator   data)
      : table(table)
      , data(data)
    {}



    template <int N, typename T, bool C, unsigned int P>
    inline Accessor<N, T, C, P>::Accessor(const Accessor &a)
      : table(a.table)
      , data(a.data)
    {}



    template <int N, typename T, bool C, unsigned int P>
    inline Accessor<N, T, C, P - 1> Accessor<N, T, C, P>::
                                    operator[](const size_type i) const
    {
      AssertIndexRange(i, table.size()[N - P]);

      // access i-th
      // subobject. optimize on the
      // case i==0
      if (i == 0)
        return Accessor<N, T, C, P - 1>(table, data);
      else
        {
          // note: P>1, otherwise the
          // specialization would have
          // been taken!
          size_type subobject_size = table.size()[N - 1];
          for (int p = P - 1; p > 1; --p)
            subobject_size *= table.size()[N - p];
          const iterator new_data = data + i * subobject_size;
          return Accessor<N, T, C, P - 1>(table, new_data);
        }
    }



    template <int N, typename T, bool C>
    inline Accessor<N, T, C, 1>::Accessor(const TableType &table,
                                          const iterator   data)
      : table(table)
      , data(data)
    {}



    template <int N, typename T, bool C>
    inline Accessor<N, T, C, 1>::Accessor(const Accessor &a)
      : table(a.table)
      , data(a.data)
    {}



    template <int N, typename T, bool C>
    inline typename Accessor<N, T, C, 1>::reference Accessor<N, T, C, 1>::
                                                    operator[](const size_type i) const
    {
      AssertIndexRange(i, table.size()[N - 1]);
      return *(data + i);
    }



    template <int N, typename T, bool C>
    inline typename Accessor<N, T, C, 1>::size_type
    Accessor<N, T, C, 1>::size() const
    {
      return table.size()[N - 1];
    }



    template <int N, typename T, bool C>
    inline typename Accessor<N, T, C, 1>::iterator
    Accessor<N, T, C, 1>::begin() const
    {
      return data;
    }



    template <int N, typename T, bool C>
    inline typename Accessor<N, T, C, 1>::iterator
    Accessor<N, T, C, 1>::end() const
    {
      return data + table.size()[N - 1];
    }
  } // namespace TableBaseAccessors
} // namespace internal



template <int N, typename T>
inline TableBase<N, T> &
TableBase<N, T>::operator=(const TableBase<N, T> &m)
{
  if (!m.empty())
    values = m.values;
  reinit(m.size(), true);

  return *this;
}



template <int N, typename T>
template <typename T2>
inline TableBase<N, T> &
TableBase<N, T>::operator=(const TableBase<N, T2> &m)
{
  reinit(m.size(), true);
  if (!empty())
    std::copy(m.values.begin(),
              m.values.begin() + n_elements(),
              values.begin());

  return *this;
}



template <int N, typename T>
inline TableBase<N, T> &
TableBase<N, T>::operator=(TableBase<N, T> &&m) noexcept
{
  static_cast<Subscriptor &>(*this) = std::move(static_cast<Subscriptor &>(m));
  values                            = std::move(m.values);
  table_size                        = m.table_size;
  m.table_size                      = TableIndices<N>();

  return *this;
}



template <int N, typename T>
inline bool
TableBase<N, T>::operator==(const TableBase<N, T> &T2) const
{
  return (values == T2.values);
}



template <int N, typename T>
inline void
TableBase<N, T>::reset_values()
{
  // use parallel set operation
  if (n_elements() != 0)
    values.fill(T());
}



template <int N, typename T>
inline void
TableBase<N, T>::fill(const T &value)
{
  if (n_elements() != 0)
    values.fill(value);
}



template <int N, typename T>
inline void
TableBase<N, T>::replicate_across_communicator(const MPI_Comm &   communicator,
                                               const unsigned int root_process)
{
  // Replicate first the actual data, then also exchange the
  // extents of the table
  values.replicate_across_communicator(communicator, root_process);

  table_size =
    Utilities::MPI::broadcast(communicator, table_size, root_process);
}



template <int N, typename T>
inline void
TableBase<N, T>::reinit(const TableIndices<N> &new_sizes,
                        const bool             omit_default_initialization)
{
  table_size = new_sizes;

  const size_type new_size = n_elements();

  // if zero size was given: free all memory
  if (new_size == 0)
    {
      values.resize(0);
      // set all sizes to zero, even
      // if one was previously
      // nonzero. This simplifies
      // some assertions.
      table_size = TableIndices<N>();

      return;
    }

  // adjust values field. If it was empty before, we can simply call resize(),
  // which can set all the data fields. Otherwise, select the fast resize and
  // manually fill in all the elements as required by the design of this
  // class. (Selecting another code for the empty case ensures that we touch
  // the memory only once for non-trivial classes that need to initialize the
  // memory also in resize_fast.)
  if (!omit_default_initialization)
    {
      if (values.empty())
        values.resize(new_size);
      else
        {
          values.resize_fast(new_size);
          values.fill();
        }
    }
  else
    values.resize_fast(new_size);
}



template <int N, typename T>
inline const TableIndices<N> &
TableBase<N, T>::size() const
{
  return table_size;
}



template <int N, typename T>
inline typename TableBase<N, T>::size_type
TableBase<N, T>::size(const unsigned int i) const
{
  AssertIndexRange(i, N);
  return table_size[i];
}



template <int N, typename T>
inline typename TableBase<N, T>::size_type
TableBase<N, T>::n_elements() const
{
  size_type s = 1;
  for (unsigned int n = 0; n < N; ++n)
    s *= table_size[n];
  return s;
}



template <int N, typename T>
inline bool
TableBase<N, T>::empty() const
{
  return (n_elements() == 0);
}



namespace internal
{
  namespace TableImplementation
  {
    template <typename InputIterator, typename T>
    void
    fill_Fortran_style(InputIterator entries, TableBase<1, T> &table)
    {
      using size_type = typename TableBase<1, T>::size_type;
      for (size_type i = 0; i < table.size()[0]; ++i)
        table(TableIndices<1>(i)) = *entries++;
    }


    template <typename InputIterator, typename T>
    void
    fill_Fortran_style(InputIterator entries, TableBase<2, T> &table)
    {
      using size_type = typename TableBase<2, T>::size_type;
      for (size_type j = 0; j < table.size()[1]; ++j)
        for (size_type i = 0; i < table.size()[0]; ++i)
          table(TableIndices<2>(i, j)) = *entries++;
    }


    template <typename InputIterator, typename T>
    void
    fill_Fortran_style(InputIterator entries, TableBase<3, T> &table)
    {
      using size_type = typename TableBase<3, T>::size_type;
      for (size_type k = 0; k < table.size()[2]; ++k)
        for (size_type j = 0; j < table.size()[1]; ++j)
          for (size_type i = 0; i < table.size()[0]; ++i)
            table(TableIndices<3>(i, j, k)) = *entries++;
    }


    template <typename InputIterator, typename T, int N>
    void
    fill_Fortran_style(InputIterator, TableBase<N, T> &)
    {
      Assert(false, ExcNotImplemented());
    }
  } // namespace TableImplementation
} // namespace internal


template <int N, typename T>
template <typename InputIterator>
inline void
TableBase<N, T>::fill(InputIterator entries, const bool C_style_indexing)
{
  Assert(n_elements() != 0, ExcMessage("Trying to fill an empty matrix."));

  if (C_style_indexing)
    for (typename AlignedVector<T>::iterator p = values.begin();
         p != values.end();
         ++p)
      *p = *entries++;
  else
    internal::TableImplementation::fill_Fortran_style(entries, *this);
}



template <int N, typename T>
inline void
TableBase<N, T>::swap(TableBase<N, T> &v)
{
  values.swap(v.values);
  std::swap(table_size, v.table_size);
}



template <int N, typename T>
inline std::size_t
TableBase<N, T>::memory_consumption() const
{
  return sizeof(*this) + MemoryConsumption::memory_consumption(values);
}



template <int N, typename T>
inline typename TableBase<N, T>::size_type
TableBase<N, T>::position(const TableIndices<N> &indices) const
{
  // specialize this for the
  // different numbers of dimensions,
  // to make the job somewhat easier
  // for the compiler. have the
  // general formula nevertheless:
  switch (N)
    {
      case 1:
        return indices[0];
      case 2:
        return size_type(indices[0]) * table_size[1] + indices[1];
      case 3:
        return ((size_type(indices[0]) * table_size[1] + indices[1]) *
                  table_size[2] +
                indices[2]);
      default:
        {
          unsigned int s = indices[0];
          for (unsigned int n = 1; n < N; ++n)
            s = s * table_size[n] + indices[n];
          return s;
        }
    }
}



template <int N, typename T>
inline typename AlignedVector<T>::const_reference
TableBase<N, T>::operator()(const TableIndices<N> &indices) const
{
  for (unsigned int n = 0; n < N; ++n)
    AssertIndexRange(indices[n], table_size[n]);
  return el(indices);
}



template <int N, typename T>
inline typename AlignedVector<T>::reference
TableBase<N, T>::operator()(const TableIndices<N> &indices)
{
  for (unsigned int n = 0; n < N; ++n)
    AssertIndexRange(indices[n], table_size[n]);
  return el(indices);
}



template <int N, typename T>
inline typename AlignedVector<T>::const_reference
TableBase<N, T>::el(const TableIndices<N> &indices) const
{
  return values[position(indices)];
}



template <int N, typename T>
inline typename AlignedVector<T>::reference
TableBase<N, T>::el(const TableIndices<N> &indices)
{
  AssertIndexRange(position(indices), values.size());
  return values[position(indices)];
}



template <typename T>
inline Table<1, T>::Table(const size_type size)
  : TableBase<1, T>(TableIndices<1>(size))
{}



template <typename T>
template <typename InputIterator>
inline Table<1, T>::Table(const size_type size,
                          InputIterator   entries,
                          const bool      C_style_indexing)
  : TableBase<1, T>(TableIndices<1>(size), entries, C_style_indexing)
{}



template <typename T>
inline typename AlignedVector<T>::const_reference Table<1, T>::
                                                  operator[](const size_type i) const
{
  AssertIndexRange(i, this->table_size[0]);
  return this->values[i];
}



template <typename T>
inline typename AlignedVector<T>::reference Table<1, T>::
                                            operator[](const size_type i)
{
  AssertIndexRange(i, this->table_size[0]);
  return this->values[i];
}



template <typename T>
inline typename AlignedVector<T>::const_reference
Table<1, T>::operator()(const size_type i) const
{
  AssertIndexRange(i, this->table_size[0]);
  return this->values[i];
}



template <typename T>
inline typename AlignedVector<T>::reference
Table<1, T>::operator()(const size_type i)
{
  AssertIndexRange(i, this->table_size[0]);
  return this->values[i];
}


//---------------------------------------------------------------------------



template <typename T>
inline Table<2, T>::Table(const size_type size1, const size_type size2)
  : TableBase<2, T>(TableIndices<2>(size1, size2))
{}



template <typename T>
template <typename InputIterator>
inline Table<2, T>::Table(const size_type size1,
                          const size_type size2,
                          InputIterator   entries,
                          const bool      C_style_indexing)
  : TableBase<2, T>(TableIndices<2>(size1, size2), entries, C_style_indexing)
{}



template <typename T>
inline void
Table<2, T>::reinit(const size_type size1,
                    const size_type size2,
                    const bool      omit_default_initialization)
{
  this->TableBase<2, T>::reinit(TableIndices<2>(size1, size2),
                                omit_default_initialization);
}



template <typename T>
inline dealii::internal::TableBaseAccessors::Accessor<2, T, true, 1>
  Table<2, T>::operator[](const size_type i) const
{
  AssertIndexRange(i, this->table_size[0]);
  return dealii::internal::TableBaseAccessors::Accessor<2, T, true, 1>(
    *this, this->values.begin() + size_type(i) * n_cols());
}



template <typename T>
inline dealii::internal::TableBaseAccessors::Accessor<2, T, false, 1>
  Table<2, T>::operator[](const size_type i)
{
  AssertIndexRange(i, this->table_size[0]);
  return dealii::internal::TableBaseAccessors::Accessor<2, T, false, 1>(
    *this, this->values.begin() + size_type(i) * n_cols());
}



template <typename T>
inline typename AlignedVector<T>::const_reference
Table<2, T>::operator()(const size_type i, const size_type j) const
{
  AssertIndexRange(i, this->table_size[0]);
  AssertIndexRange(j, this->table_size[1]);
  return this->values[size_type(i) * this->table_size[1] + j];
}



template <typename T>
inline typename AlignedVector<T>::reference
Table<2, T>::operator()(const size_type i, const size_type j)
{
  AssertIndexRange(i, this->table_size[0]);
  AssertIndexRange(j, this->table_size[1]);
  return this->values[size_type(i) * this->table_size[1] + j];
}



template <typename T>
inline typename AlignedVector<T>::const_reference
Table<2, T>::el(const size_type i, const size_type j) const
{
  return this->values[size_type(i) * this->table_size[1] + j];
}



template <typename T>
inline typename AlignedVector<T>::reference
Table<2, T>::el(const size_type i, const size_type j)
{
  return this->values[size_type(i) * this->table_size[1] + j];
}



template <typename T>
inline typename Table<2, T>::size_type
Table<2, T>::n_rows() const
{
  return this->table_size[0];
}



template <typename T>
inline typename Table<2, T>::size_type
Table<2, T>::n_cols() const
{
  return this->table_size[1];
}



template <typename T>
inline typename Table<2, T>::iterator
Table<2, T>::begin()
{
  return typename Table<2, T>::iterator(this, 0, 0);
}



template <typename T>
inline typename Table<2, T>::const_iterator
Table<2, T>::begin() const
{
  return typename Table<2, T>::const_iterator(this, 0, 0);
}



template <typename T>
inline typename Table<2, T>::iterator
Table<2, T>::end()
{
  return typename Table<2, T>::iterator(this);
}



template <typename T>
inline typename Table<2, T>::const_iterator
Table<2, T>::end() const
{
  return typename Table<2, T>::const_iterator(this);
}



//---------------------------------------------------------------------------
namespace MatrixTableIterators
{
  namespace internal
  {
    // Internal calculation routines for AccessorBase. These do not do any
    // checking whatsoever.
    template <typename TableType, Storage storage_order>
    inline std::ptrdiff_t
    get_row_index(const std::ptrdiff_t   linear_index,
                  const TableType *const container)
    {
      switch (storage_order)
        {
          case Storage::row_major:
            return linear_index / container->n_cols();
          case Storage::column_major:
            return linear_index % container->n_rows();
          default:
            Assert(false, ExcInternalError());
        }
      return {};
    }



    template <typename TableType, Storage storage_order>
    inline std::ptrdiff_t
    get_column_index(const std::ptrdiff_t   linear_index,
                     const TableType *const container)
    {
      switch (storage_order)
        {
          case Storage::row_major:
            return linear_index % container->n_cols();
          case Storage::column_major:
            return linear_index / container->n_rows();
          default:
            Assert(false, ExcInternalError());
        }
      return {};
    }
  } // namespace internal



  template <typename TableType, bool Constness, Storage storage_order>
  inline AccessorBase<TableType, Constness, storage_order>::AccessorBase()
    : container(nullptr)
    , linear_index(std::numeric_limits<decltype(linear_index)>::max())
  {}



  template <typename TableType, bool Constness, Storage storage_order>
  inline AccessorBase<TableType, Constness, storage_order>::AccessorBase(
    const container_pointer_type table)
    : container(table)
    , linear_index(container->values.size())
  {}



  template <typename TableType, bool Constness, Storage storage_order>
  inline AccessorBase<TableType, Constness, storage_order>::AccessorBase(
    const AccessorBase<TableType, false, storage_order> &a)
    : container(a.container)
    , linear_index(a.linear_index)
  {}



  template <typename TableType, bool Constness, Storage storage_order>
  inline AccessorBase<TableType, Constness, storage_order>::AccessorBase(
    const container_pointer_type table,
    const std::ptrdiff_t         index)
    : container(table)
    , linear_index(index)
  {
    Assert(0 <= linear_index &&
             std::size_t(linear_index) < container->values.size() + 1,
           ExcMessage("The current iterator points outside of the table and is "
                      "not the end iterator."));
  }



  template <typename TableType, bool Constness, Storage storage_order>
  inline const typename AccessorBase<TableType, Constness, storage_order>::
    value_type &
    AccessorBase<TableType, Constness, storage_order>::value() const
  {
    assert_valid_linear_index();
    return this->container->values[linear_index];
  }



  template <typename TableType, bool Constness, Storage storage_order>
  inline AccessorBase<TableType, Constness, storage_order>::
  operator const value_type &() const
  {
    assert_valid_linear_index();
    return this->container->values[linear_index];
  }



  template <typename TableType, bool Constness, Storage storage_order>
  inline typename AccessorBase<TableType, Constness, storage_order>::size_type
  AccessorBase<TableType, Constness, storage_order>::row() const
  {
    assert_valid_linear_index();
    return static_cast<std::size_t>(
      internal::get_row_index<TableType, storage_order>(linear_index,
                                                        container));
  }



  template <typename TableType, bool Constness, Storage storage_order>
  inline typename AccessorBase<TableType, Constness, storage_order>::size_type
  AccessorBase<TableType, Constness, storage_order>::column() const
  {
    assert_valid_linear_index();
    return static_cast<std::size_t>(
      internal::get_column_index<TableType, storage_order>(linear_index,
                                                           container));
  }



  template <typename TableType, bool Constness, Storage storage_order>
  inline void
  AccessorBase<TableType, Constness, storage_order>::assert_valid_linear_index()
    const
  {
#  ifdef DEBUG // avoid unused variable warnings by guarding everything
    Assert(container != nullptr,
           ExcMessage("This accessor has been default-constructed and does not "
                      "have a corresponding table."));
    Assert(!container->empty(),
           ExcMessage("An empty table has no rows or columns."));
    Assert(0 <= linear_index &&
             std::size_t(linear_index) < container->values.size(),
           ExcMessage("The current iterator points outside of the table."));
    const std::ptrdiff_t row_n =
      internal::get_row_index<TableType, storage_order>(linear_index,
                                                        container);
    const std::ptrdiff_t column_n =
      internal::get_column_index<TableType, storage_order>(linear_index,
                                                           container);
    Assert(0 <= column_n && std::size_t(column_n) < container->n_cols(),
           ExcMessage("The current iterator points outside the table."));
    Assert(0 <= row_n && std::size_t(row_n) < container->n_rows(),
           ExcMessage("The current iterator points outside the table."));
#  endif
  }



  template <typename TableType, Storage storage_order>
  inline const Accessor<TableType, false, storage_order> &
  Accessor<TableType, false, storage_order>::operator=(
    const typename Accessor<TableType, false, storage_order>::value_type &t)
    const
  {
    this->assert_valid_linear_index();
    this->container->values[this->linear_index] = t;
    return *this;
  }



  template <typename TableType, Storage storage_order>
  inline const Accessor<TableType, false, storage_order> &
  Accessor<TableType, false, storage_order>::operator=(
    typename Accessor<TableType, false, storage_order>::value_type &&t) const
  {
    this->assert_valid_linear_index();
    this->container->values[this->linear_index] = t;
    return *this;
  }



  template <typename TableType, Storage storage_order>
  inline typename Accessor<TableType, false, storage_order>::value_type &
  Accessor<TableType, false, storage_order>::value() const
  {
    this->assert_valid_linear_index();
    return this->container->values[this->linear_index];
  }



  template <typename TableType, Storage storage_order>
  inline Accessor<TableType, false, storage_order>::operator value_type &()
  {
    this->assert_valid_linear_index();
    return this->container->values[this->linear_index];
  }



  template <typename TableType, bool Constness, Storage storage_order>
  inline Iterator<TableType, Constness, storage_order>::Iterator(
    const Accessor<TableType, Constness, storage_order> &a)
    : LinearIndexIterator<Iterator<TableType, Constness, storage_order>,
                          Accessor<TableType, Constness, storage_order>>(a)
  {}



  template <typename TableType, bool Constness, Storage storage_order>
  inline Iterator<TableType, Constness, storage_order>::Iterator(
    const container_pointer_type table)
    : LinearIndexIterator<Iterator<TableType, Constness, storage_order>,
                          Accessor<TableType, Constness, storage_order>>(
        Accessor<TableType, Constness, storage_order>(table))
  {}



  template <typename TableType, bool Constness, Storage storage_order>
  inline Iterator<TableType, Constness, storage_order>::Iterator(
    const Iterator<TableType, false, storage_order> &i)
    : LinearIndexIterator<Iterator<TableType, Constness, storage_order>,
                          Accessor<TableType, Constness, storage_order>>(*i)
  {}



  template <typename TableType, bool Constness, Storage storage_order>
  inline Iterator<TableType, Constness, storage_order>::Iterator(
    const container_pointer_type table,
    const size_type              row_n,
    const size_type              col_n)
    : Iterator(table,
               storage_order == Storage::row_major ?
                 table->n_cols() * row_n + col_n :
                 table->n_rows() * col_n + row_n)
  {}



  template <typename TableType, bool Constness, Storage storage_order>
  inline Iterator<TableType, Constness, storage_order>::Iterator(
    const container_pointer_type table,
    const std::ptrdiff_t         linear_index)
    : LinearIndexIterator<Iterator<TableType, Constness, storage_order>,
                          Accessor<TableType, Constness, storage_order>>(
        Accessor<TableType, Constness, storage_order>(table, linear_index))
  {}
} // namespace MatrixTableIterators



//---------------------------------------------------------------------------
template <typename T>
inline TransposeTable<T>::TransposeTable(const size_type size1,
                                         const size_type size2)
  : TableBase<2, T>(TableIndices<2>(size2, size1))
{}



template <typename T>
inline void
TransposeTable<T>::reinit(const size_type size1,
                          const size_type size2,
                          const bool      omit_default_initialization)
{
  this->TableBase<2, T>::reinit(TableIndices<2>(size2, size1),
                                omit_default_initialization);
}



template <typename T>
inline typename TransposeTable<T>::const_reference
TransposeTable<T>::operator()(const size_type i, const size_type j) const
{
  AssertIndexRange(i, this->table_size[1]);
  AssertIndexRange(j, this->table_size[0]);
  return this->values[size_type(j) * this->table_size[1] + i];
}



template <typename T>
inline typename TransposeTable<T>::reference
TransposeTable<T>::operator()(const size_type i, const size_type j)
{
  AssertIndexRange(i, this->table_size[1]);
  AssertIndexRange(j, this->table_size[0]);
  return this->values[size_type(j) * this->table_size[1] + i];
}



template <typename T>
inline typename TransposeTable<T>::const_reference
TransposeTable<T>::el(const size_type i, const size_type j) const
{
  return this->values[size_type(j) * this->table_size[1] + i];
}



template <typename T>
inline typename TransposeTable<T>::reference
TransposeTable<T>::el(const size_type i, const size_type j)
{
  return this->values[size_type(j) * this->table_size[1] + i];
}



template <typename T>
inline typename TransposeTable<T>::size_type
TransposeTable<T>::n_rows() const
{
  return this->table_size[1];
}



template <typename T>
inline typename TransposeTable<T>::size_type
TransposeTable<T>::n_cols() const
{
  return this->table_size[0];
}



template <typename T>
inline typename TransposeTable<T>::iterator
TransposeTable<T>::begin()
{
  return typename TransposeTable<T>::iterator(this, 0, 0);
}



template <typename T>
inline typename TransposeTable<T>::const_iterator
TransposeTable<T>::begin() const
{
  return typename TransposeTable<T>::const_iterator(this, 0, 0);
}



template <typename T>
inline typename TransposeTable<T>::iterator
TransposeTable<T>::end()
{
  return typename TransposeTable<T>::iterator(this);
}



template <typename T>
inline typename TransposeTable<T>::const_iterator
TransposeTable<T>::end() const
{
  return typename TransposeTable<T>::const_iterator(this);
}
//---------------------------------------------------------------------------



template <typename T>
inline Table<3, T>::Table(const size_type size1,
                          const size_type size2,
                          const size_type size3)
  : TableBase<3, T>(TableIndices<3>(size1, size2, size3))
{}



template <typename T>
template <typename InputIterator>
inline Table<3, T>::Table(const size_type size1,
                          const size_type size2,
                          const size_type size3,
                          InputIterator   entries,
                          const bool      C_style_indexing)
  : TableBase<3, T>(TableIndices<3>(size1, size2, size3),
                    entries,
                    C_style_indexing)
{}



template <typename T>
inline dealii::internal::TableBaseAccessors::Accessor<3, T, true, 2>
  Table<3, T>::operator[](const size_type i) const
{
  AssertIndexRange(i, this->table_size[0]);
  const size_type subobject_size =
    size_type(this->table_size[1]) * this->table_size[2];
  return (dealii::internal::TableBaseAccessors::Accessor<3, T, true, 2>(
    *this, this->values.begin() + i * subobject_size));
}



template <typename T>
inline dealii::internal::TableBaseAccessors::Accessor<3, T, false, 2>
  Table<3, T>::operator[](const size_type i)
{
  AssertIndexRange(i, this->table_size[0]);
  const size_type subobject_size =
    size_type(this->table_size[1]) * this->table_size[2];
  return (dealii::internal::TableBaseAccessors::Accessor<3, T, false, 2>(
    *this, this->values.begin() + i * subobject_size));
}



template <typename T>
inline typename AlignedVector<T>::const_reference
Table<3, T>::
operator()(const size_type i, const size_type j, const size_type k) const
{
  AssertIndexRange(i, this->table_size[0]);
  AssertIndexRange(j, this->table_size[1]);
  AssertIndexRange(k, this->table_size[2]);
  return this
    ->values[(size_type(i) * this->table_size[1] + j) * this->table_size[2] +
             k];
}



template <typename T>
inline typename AlignedVector<T>::reference
Table<3, T>::operator()(const size_type i, const size_type j, const size_type k)
{
  AssertIndexRange(i, this->table_size[0]);
  AssertIndexRange(j, this->table_size[1]);
  AssertIndexRange(k, this->table_size[2]);
  return this
    ->values[(size_type(i) * this->table_size[1] + j) * this->table_size[2] +
             k];
}



template <typename T>
inline Table<4, T>::Table(const size_type size1,
                          const size_type size2,
                          const size_type size3,
                          const size_type size4)
  : TableBase<4, T>(TableIndices<4>(size1, size2, size3, size4))
{}



template <typename T>
inline dealii::internal::TableBaseAccessors::Accessor<4, T, true, 3>
  Table<4, T>::operator[](const size_type i) const
{
  AssertIndexRange(i, this->table_size[0]);
  const size_type subobject_size =
    size_type(this->table_size[1]) * this->table_size[2] * this->table_size[3];
  return (dealii::internal::TableBaseAccessors::Accessor<4, T, true, 3>(
    *this, this->values.begin() + i * subobject_size));
}



template <typename T>
inline dealii::internal::TableBaseAccessors::Accessor<4, T, false, 3>
  Table<4, T>::operator[](const size_type i)
{
  AssertIndexRange(i, this->table_size[0]);
  const size_type subobject_size =
    size_type(this->table_size[1]) * this->table_size[2] * this->table_size[3];
  return (dealii::internal::TableBaseAccessors::Accessor<4, T, false, 3>(
    *this, this->values.begin() + i * subobject_size));
}



template <typename T>
inline typename AlignedVector<T>::const_reference
Table<4, T>::operator()(const size_type i,
                        const size_type j,
                        const size_type k,
                        const size_type l) const
{
  AssertIndexRange(i, this->table_size[0]);
  AssertIndexRange(j, this->table_size[1]);
  AssertIndexRange(k, this->table_size[2]);
  AssertIndexRange(l, this->table_size[3]);
  return this
    ->values[((size_type(i) * this->table_size[1] + j) * this->table_size[2] +
              k) *
               this->table_size[3] +
             l];
}



template <typename T>
inline typename AlignedVector<T>::reference
Table<4, T>::operator()(const size_type i,
                        const size_type j,
                        const size_type k,
                        const size_type l)
{
  AssertIndexRange(i, this->table_size[0]);
  AssertIndexRange(j, this->table_size[1]);
  AssertIndexRange(k, this->table_size[2]);
  AssertIndexRange(l, this->table_size[3]);
  return this
    ->values[((size_type(i) * this->table_size[1] + j) * this->table_size[2] +
              k) *
               this->table_size[3] +
             l];
}



template <typename T>
inline Table<5, T>::Table(const size_type size1,
                          const size_type size2,
                          const size_type size3,
                          const size_type size4,
                          const size_type size5)
  : TableBase<5, T>(TableIndices<5>(size1, size2, size3, size4, size5))
{}



template <typename T>
inline dealii::internal::TableBaseAccessors::Accessor<5, T, true, 4>
  Table<5, T>::operator[](const size_type i) const
{
  AssertIndexRange(i, this->table_size[0]);
  const size_type subobject_size = size_type(this->table_size[1]) *
                                   this->table_size[2] * this->table_size[3] *
                                   this->table_size[4];
  return (dealii::internal::TableBaseAccessors::Accessor<5, T, true, 4>(
    *this, this->values.begin() + i * subobject_size));
}



template <typename T>
inline dealii::internal::TableBaseAccessors::Accessor<5, T, false, 4>
  Table<5, T>::operator[](const size_type i)
{
  AssertIndexRange(i, this->table_size[0]);
  const size_type subobject_size = size_type(this->table_size[1]) *
                                   this->table_size[2] * this->table_size[3] *
                                   this->table_size[4];
  return (dealii::internal::TableBaseAccessors::Accessor<5, T, false, 4>(
    *this, this->values.begin() + i * subobject_size));
}



template <typename T>
inline typename AlignedVector<T>::const_reference
Table<5, T>::operator()(const size_type i,
                        const size_type j,
                        const size_type k,
                        const size_type l,
                        const size_type m) const
{
  AssertIndexRange(i, this->table_size[0]);
  AssertIndexRange(j, this->table_size[1]);
  AssertIndexRange(k, this->table_size[2]);
  AssertIndexRange(l, this->table_size[3]);
  AssertIndexRange(m, this->table_size[4]);
  return this
    ->values[(((size_type(i) * this->table_size[1] + j) * this->table_size[2] +
               k) *
                this->table_size[3] +
              l) *
               this->table_size[4] +
             m];
}



template <typename T>
inline typename AlignedVector<T>::reference
Table<5, T>::operator()(const size_type i,
                        const size_type j,
                        const size_type k,
                        const size_type l,
                        const size_type m)
{
  AssertIndexRange(i, this->table_size[0]);
  AssertIndexRange(j, this->table_size[1]);
  AssertIndexRange(k, this->table_size[2]);
  AssertIndexRange(l, this->table_size[3]);
  AssertIndexRange(m, this->table_size[4]);
  return this
    ->values[(((size_type(i) * this->table_size[1] + j) * this->table_size[2] +
               k) *
                this->table_size[3] +
              l) *
               this->table_size[4] +
             m];
}



template <typename T>
inline Table<6, T>::Table(const size_type size1,
                          const size_type size2,
                          const size_type size3,
                          const size_type size4,
                          const size_type size5,
                          const size_type size6)
  : TableBase<6, T>()
{
  TableIndices<6> table_indices;
  table_indices[0] = size1;
  table_indices[1] = size2;
  table_indices[2] = size3;
  table_indices[3] = size4;
  table_indices[4] = size5;
  table_indices[5] = size6;

  TableBase<6, T>::reinit(table_indices);
}



template <typename T>
inline dealii::internal::TableBaseAccessors::Accessor<6, T, true, 5>
  Table<6, T>::operator[](const size_type i) const
{
  AssertIndexRange(i, this->table_size[0]);
  const size_type subobject_size = size_type(this->table_size[1]) *
                                   this->table_size[2] * this->table_size[3] *
                                   this->table_size[4] * this->table_size[5];
  return (dealii::internal::TableBaseAccessors::Accessor<6, T, true, 5>(
    *this, this->values.begin() + i * subobject_size));
}



template <typename T>
inline dealii::internal::TableBaseAccessors::Accessor<6, T, false, 5>
  Table<6, T>::operator[](const size_type i)
{
  AssertIndexRange(i, this->table_size[0]);
  const size_type subobject_size = size_type(this->table_size[1]) *
                                   this->table_size[2] * this->table_size[3] *
                                   this->table_size[4] * this->table_size[5];
  return (dealii::internal::TableBaseAccessors::Accessor<6, T, false, 5>(
    *this, this->values.begin() + i * subobject_size));
}



template <typename T>
inline typename AlignedVector<T>::const_reference
Table<6, T>::operator()(const size_type i,
                        const size_type j,
                        const size_type k,
                        const size_type l,
                        const size_type m,
                        const size_type n) const
{
  AssertIndexRange(i, this->table_size[0]);
  AssertIndexRange(j, this->table_size[1]);
  AssertIndexRange(k, this->table_size[2]);
  AssertIndexRange(l, this->table_size[3]);
  AssertIndexRange(m, this->table_size[4]);
  AssertIndexRange(n, this->table_size[5]);
  return this
    ->values[((((size_type(i) * this->table_size[1] + j) * this->table_size[2] +
                k) *
                 this->table_size[3] +
               l) *
                this->table_size[4] +
              m) *
               this->table_size[5] +
             n];
}



template <typename T>
inline typename AlignedVector<T>::reference
Table<6, T>::operator()(const size_type i,
                        const size_type j,
                        const size_type k,
                        const size_type l,
                        const size_type m,
                        const size_type n)
{
  AssertIndexRange(i, this->table_size[0]);
  AssertIndexRange(j, this->table_size[1]);
  AssertIndexRange(k, this->table_size[2]);
  AssertIndexRange(l, this->table_size[3]);
  AssertIndexRange(m, this->table_size[4]);
  AssertIndexRange(n, this->table_size[5]);
  return this
    ->values[((((size_type(i) * this->table_size[1] + j) * this->table_size[2] +
                k) *
                 this->table_size[3] +
               l) *
                this->table_size[4] +
              m) *
               this->table_size[5] +
             n];
}



template <typename T>
inline Table<7, T>::Table(const size_type size1,
                          const size_type size2,
                          const size_type size3,
                          const size_type size4,
                          const size_type size5,
                          const size_type size6,
                          const size_type size7)
  : TableBase<7, T>()
{
  TableIndices<7> table_indices;
  table_indices[0] = size1;
  table_indices[1] = size2;
  table_indices[2] = size3;
  table_indices[3] = size4;
  table_indices[4] = size5;
  table_indices[5] = size6;
  table_indices[6] = size7;

  TableBase<7, T>::reinit(table_indices);
}



template <typename T>
inline dealii::internal::TableBaseAccessors::Accessor<7, T, true, 6>
  Table<7, T>::operator[](const size_type i) const
{
  AssertIndexRange(i, this->table_size[0]);
  const size_type subobject_size =
    size_type(this->table_size[1]) * this->table_size[2] * this->table_size[3] *
    this->table_size[4] * this->table_size[5] * this->table_size[6];
  return (dealii::internal::TableBaseAccessors::Accessor<7, T, true, 6>(
    *this, this->values.begin() + i * subobject_size));
}



template <typename T>
inline dealii::internal::TableBaseAccessors::Accessor<7, T, false, 6>
  Table<7, T>::operator[](const size_type i)
{
  AssertIndexRange(i, this->table_size[0]);
  const size_type subobject_size =
    size_type(this->table_size[1]) * this->table_size[2] * this->table_size[3] *
    this->table_size[4] * this->table_size[5] * this->table_size[6];
  return (dealii::internal::TableBaseAccessors::Accessor<7, T, false, 6>(
    *this, this->values.begin() + i * subobject_size));
}



template <typename T>
inline typename AlignedVector<T>::const_reference
Table<7, T>::operator()(const size_type i,
                        const size_type j,
                        const size_type k,
                        const size_type l,
                        const size_type m,
                        const size_type n,
                        const size_type o) const
{
  AssertIndexRange(i, this->table_size[0]);
  AssertIndexRange(j, this->table_size[1]);
  AssertIndexRange(k, this->table_size[2]);
  AssertIndexRange(l, this->table_size[3]);
  AssertIndexRange(m, this->table_size[4]);
  AssertIndexRange(n, this->table_size[5]);
  AssertIndexRange(o, this->table_size[6]);
  return this->values
    [(((((size_type(i) * this->table_size[1] + j) * this->table_size[2] + k) *
          this->table_size[3] +
        l) *
         this->table_size[4] +
       m) *
        this->table_size[5] +
      n) *
       this->table_size[6] +
     o];
}



template <typename T>
inline typename AlignedVector<T>::reference
Table<7, T>::operator()(const size_type i,
                        const size_type j,
                        const size_type k,
                        const size_type l,
                        const size_type m,
                        const size_type n,
                        const size_type o)
{
  AssertIndexRange(i, this->table_size[0]);
  AssertIndexRange(j, this->table_size[1]);
  AssertIndexRange(k, this->table_size[2]);
  AssertIndexRange(l, this->table_size[3]);
  AssertIndexRange(m, this->table_size[4]);
  AssertIndexRange(n, this->table_size[5]);
  AssertIndexRange(o, this->table_size[6]);
  return this->values
    [(((((size_type(i) * this->table_size[1] + j) * this->table_size[2] + k) *
          this->table_size[3] +
        l) *
         this->table_size[4] +
       m) *
        this->table_size[5] +
      n) *
       this->table_size[6] +
     o];
}


#endif // DOXYGEN



/**
 * 全局函数 @p swap
 * ，它重载了C++标准库的默认实现，它使用一个临时对象。该函数简单地交换了两个表的数据。
 *
 *
 */
template <int N, typename T>
inline void
swap(TableBase<N, T> &u, TableBase<N, T> &v)
{
  u.swap(v);
}

DEAL_II_NAMESPACE_CLOSE

#endif


