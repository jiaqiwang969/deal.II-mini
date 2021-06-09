//include/deal.II-translator/lac/sparse_matrix_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2021 by the deal.II authors
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

#ifndef dealii_sparse_matrix_h
#  define dealii_sparse_matrix_h


#  include <deal.II/base/config.h>

#  include <deal.II/base/smartpointer.h>
#  include <deal.II/base/subscriptor.h>

#  include <deal.II/lac/exceptions.h>
#  include <deal.II/lac/identity_matrix.h>
#  include <deal.II/lac/sparsity_pattern.h>
#  include <deal.II/lac/vector_operation.h>
#  ifdef DEAL_II_WITH_MPI
#    include <mpi.h>
#  endif

#  include <memory>


DEAL_II_NAMESPACE_OPEN

// Forward declarations
#  ifndef DOXYGEN
template <typename number>
class Vector;
template <typename number>
class FullMatrix;
template <typename Matrix>
class BlockMatrixBase;
template <typename number>
class SparseILU;
#    ifdef DEAL_II_WITH_MPI
namespace Utilities
{
  namespace MPI
  {
    template <typename Number>
    void
    sum(const SparseMatrix<Number> &, const MPI_Comm &, SparseMatrix<Number> &);
  }
} // namespace Utilities
#    endif

#    ifdef DEAL_II_WITH_TRILINOS
namespace TrilinosWrappers
{
  class SparseMatrix;
}
#    endif
#  endif

/**
 * @addtogroup  Matrix1  
     * @{ 
 *
 *
 */

/**
 * 一个命名空间，我们在其中声明对稀疏矩阵元素的迭代器。
 *
 *
 */
namespace SparseMatrixIterators
{
  /**
   * 声明容器大小的类型。
   *
   */
  using size_type = types::global_dof_index;

  // forward declaration
  template <typename number, bool Constness>
  class Iterator;

  /**
   * 稀疏矩阵访问器的通用模板。第一个模板参数表示底层数字类型，第二个表示矩阵的常数。
   * 通用模板没有被实现，只有针对第二个模板参数的两个可能值的特殊化。因此，这里列出的接口只是作为提供的模板，因为doxygen并没有链接这些特殊化。
   *
   */
  template <typename number, bool Constness>
  class Accessor : public SparsityPatternIterators::Accessor
  {
  public:
    /**
     * 这个矩阵条目的值。
     *
     */
    number
    value() const;

    /**
     * 此矩阵条目的值。
     *
     */
    number &
    value();

    /**
     * 返回该访问器所指向的矩阵的引用。注意，在本例中，这是一个常数引用。
     *
     */
    const SparseMatrix<number> &
    get_matrix() const;
  };



  /**
   * 用于常数矩阵的访问器类，在const_iterators中使用。这个类建立在用于稀疏模式的访问器类的基础上，在所有非零条目上循环，只增加了访问器函数，以获得存储在某一位置的实际值。
   *
   */
  template <typename number>
  class Accessor<number, true> : public SparsityPatternIterators::Accessor
  {
  public:
    /**
     * 这里要使用的矩阵的类型（包括常数）的类型定义。
     *
     */
    using MatrixType = const SparseMatrix<number>;

    /**
     * 构造函数。
     *
     */
    Accessor(MatrixType *matrix, const std::size_t index_within_matrix);

    /**
     * 构造器。构建给定矩阵的终端访问器。
     *
     */
    Accessor(MatrixType *matrix);

    /**
     * 复制构造器，从一个非常量访问器到一个常量访问器。
     *
     */
    Accessor(const SparseMatrixIterators::Accessor<number, false> &a);

    /**
     * 这个矩阵条目的值。
     *
     */
    number
    value() const;

    /**
     * 返回该访问器所指向的矩阵的引用。注意，在本例中，这是一个常量引用。
     *
     */
    const MatrixType &
    get_matrix() const;

  private:
    /**
     * 指向我们使用的矩阵的指针。
     *
     */
    MatrixType *matrix;

    /**
     * 使基类的前进函数可用。
     *
     */
    using SparsityPatternIterators::Accessor::advance;

    // Make iterator class a friend.
    template <typename, bool>
    friend class Iterator;
  };


  /**
   * 用于非恒定矩阵的访问器类，在迭代器中使用。这个类建立在用于稀疏模式的访问器类的基础上，在所有非零条目上循环，只增加了访问器函数，以获得存储在某一位置的实际值。
   *
   */
  template <typename number>
  class Accessor<number, false> : public SparsityPatternIterators::Accessor
  {
  private:
    /**
     * 参考类。这是访问器类在你调用value()函数时返回的东西。引用的作用就像它是对矩阵条目实际值的引用一样，也就是说，你可以读写它，可以对它进行加法和乘法，等等，但是由于矩阵并没有给出这个矩阵条目的地址，我们必须通过函数来完成这一切。
     * 构造函数需要一个指向访问器对象的指针，描述它指向矩阵的哪个元素。当人们写像iterator->value()=0（而不是iterator->value()=0.0）这样的代码时，就会产生歧义，因为右边是一个整数，既可以转换为<tt>数字</tt>（即最常见的双数），也可以转换为另一个<tt>参考</tt>类型的对象。然后编译器抱怨说不知道该采取哪种转换。
     * 由于某些原因，添加另一个重载operator=(int)似乎并不能解决这个问题。然而，我们通过向Reference构造函数添加第二个假参数来避免这个问题，该参数未被使用，但可以确保没有第二个匹配的转换序列，使用单参数的右侧。
     * 测试案例oliver_01检查了这一方法是否真的能按预期工作。
     *
     */
    class Reference
    {
    public:
      /**
       * 构造函数。关于第二个参数，请参见一般的类文档。
       *
       */
      Reference(const Accessor *accessor, const bool dummy);

      /**
       * 对矩阵的数据类型的转换操作。
       *
       */
      operator number() const;

      /**
       * 将我们目前指向的矩阵的元素设置为 @p n. 。
       *
       */
      const Reference &
      operator=(const number n) const;

      /**
       * 将 @p n 添加到我们目前指向的矩阵元素中。
       *
       */
      const Reference &
      operator+=(const number n) const;

      /**
       * 从我们现在指向的矩阵元素中减去 @p n 。
       *
       */
      const Reference &
      operator-=(const number n) const;

      /**
       * 将我们现在指向的矩阵元素乘以 @p n. 。
       *
       */
      const Reference &
      operator*=(const number n) const;

      /**
       * 用我们现在指向的矩阵的元素除以 @p n. 。
       *
       */
      const Reference &
      operator/=(const number n) const;

    private:
      /**
       * 指向访问器的指针，表示我们目前指向哪个元素。
       *
       */
      const Accessor *accessor;
    };

  public:
    /**
     * 这里要使用的矩阵的类型（包括常数）的类型定义。
     *
     */
    using MatrixType = SparseMatrix<number>;

    /**
     * 构造函数。
     *
     */
    Accessor(MatrixType *matrix, const std::size_t index);

    /**
     * 构造器。构建给定矩阵的终端访问器。
     *
     */
    Accessor(MatrixType *matrix);

    /**
     * 该矩阵条目的值，作为可读可写的引用返回。
     *
     */
    Reference
    value() const;

    /**
     * 返回该访问器所指向的矩阵的引用。注意，在本例中，这是一个非常量的引用。
     *
     */
    MatrixType &
    get_matrix() const;

  private:
    /**
     * 指向我们使用的矩阵的指针。
     *
     */
    MatrixType *matrix;

    /**
     * 使基类的前进函数可用。
     *
     */
    using SparsityPatternIterators::Accessor::advance;

    // Make iterator class a friend.
    template <typename, bool>
    friend class Iterator;
  };



  /**
   * 用于常数和非常数矩阵的迭代器。
   * 这些迭代器的典型用途是迭代稀疏矩阵的元素或单个行的元素。注意，不能保证行的元素实际上是按照列单调增加的顺序来遍历的。更多信息请参见SparsityPattern类的文档。
   * 第一个模板参数表示基础数字类型，第二个模板参数表示矩阵的常数。
   * 因为这个类有一个专门的<tt>Constness=false</tt>，所以这个类是用于常数矩阵的迭代器。
   * @note
   * 该类直接对SparsityPattern和SparseMatrix类的内部数据结构进行操作。因此，有些操作很便宜，有些则不然。特别是，访问列索引和指向的条目的值是很便宜的。另一方面，访问行索引是很昂贵的（对于一个有
   * $N$ 行的矩阵，这需要 $O(\log(N))$
   * 次操作）。因此，当你设计使用这些迭代器的算法时，通常的做法是不一次性循环稀疏矩阵的<i>all</i>个元素，而是在所有行上有一个外循环，并在这个循环中迭代这个行的元素。这样，你只需要解除对迭代器的引用来获得列索引和值，而通过使用循环索引可以避免对行索引的（昂贵）查找。
   *
   */
  template <typename number, bool Constness>
  class Iterator
  {
  public:
    /**
     * 我们要操作的矩阵类型（包括常数）的类型定义。
     *
     */
    using MatrixType = typename Accessor<number, Constness>::MatrixType;

    /**
     * 当你解除对当前类型的迭代器的定义时，你得到的类型的别名。
     *
     */
    using value_type = const Accessor<number, Constness> &;

    /**
     * 构造函数。在矩阵 @p matrix
     * 中创建一个迭代器，用于完整矩阵中的给定索引（从第2个条目开始计算）。
     *
     */
    Iterator(MatrixType *matrix, const std::size_t index_within_matrix);

    /**
     * 构造函数。为给定的矩阵创建终端迭代器。
     *
     */
    Iterator(MatrixType *matrix);

    /**
     * 转换构造函数，从非常量迭代器到常量迭代器。
     *
     */
    Iterator(const SparseMatrixIterators::Iterator<number, false> &i);

    /**
     * 从非定常迭代器到定常迭代器的复制赋值操作符。
     *
     */
    const Iterator<number, Constness> &
    operator=(const SparseMatrixIterators::Iterator<number, false> &i);

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
    const Accessor<number, Constness> &operator*() const;

    /**
     * 解除引用操作符。
     *
     */
    const Accessor<number, Constness> *operator->() const;

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
     * 比较运算符。如果第一个行号较小，或者行号相等且第一个索引较小，则结果为真。
     * 这个函数只有在两个迭代器都指向同一个矩阵时才有效。
     *
     */
    bool
    operator<(const Iterator &) const;

    /**
     * 比较运算符。与上述运算符的工作方式相同，只是反过来了。
     *
     */
    bool
    operator>(const Iterator &) const;

    /**
     * 返回当前迭代器和参数之间的距离。这个距离是通过对当前迭代器应用operator++多少次才能得到参数（对于正的返回值），或者operator--（对于负的返回值）而给出。
     *
     */
    int
    operator-(const Iterator &p) const;

    /**
     * 返回一个比当前迭代器领先 @p n 的迭代器。
     *
     */
    Iterator
    operator+(const size_type n) const;

  private:
    /**
     * 存储一个访问器类的对象。
     *
     */
    Accessor<number, Constness> accessor;
  };

} // namespace SparseMatrixIterators

/**
 * @}
 *
 *
 */


// TODO: Add multithreading to the other vmult functions.

/**
 * 稀疏矩阵。这个类实现了在SparsityPattern表示的位置存储矩阵入口值的功能。参见 @ref
 * Sparsity 中关于稀疏模式和矩阵之间分离的讨论。
 * SparseMatrix的元素是按照SparsityPattern类存储其条目的相同顺序来存储的。在每一行中，元素通常以列索引递增的顺序从左到右存储；这一规则的例外是，如果矩阵是正方形（m()
 * ==
 * n()），那么对角线条目就会被存储为每一行的第一个元素，以使应用雅可比或SSOR预处理程序等操作更快。因此，如果你用迭代器遍历稀疏矩阵的某一行的元素（使用
 * SparseMatrix::begin 和 SparseMatrix::end)
 * ，你会发现只要矩阵是平方的，每一行的元素就不会按列索引排序。
 *
 *
 * @note
 * 这个模板的实例化提供给<tt>  @<float@>  和  @<double@></tt>;
 * 其他的可以在应用程序中生成（见手册中 @ref Instantiations
 * 一节）。
 *
 *
 * @ingroup Matrix1
 *
 *
 */
template <typename number>
class SparseMatrix : public virtual Subscriptor
{
public:
  /**
   * 声明容器尺寸的类型。
   *
   */
  using size_type = types::global_dof_index;

  /**
   * 矩阵条目的类型。这个别名类似于标准库容器中的<tt>value_type</tt>。
   *
   */
  using value_type = number;

  /**
   * 声明一个类型，该类型持有与本类的模板参数相同精度的实值数。如果这个类的模板参数是一个实数数据类型，那么real_type就等于模板参数。
   * 如果模板参数是一个 std::complex
   * 类型，那么real_type等于复数的基础类型。
   * 这个别名被用来表示规范的返回类型。
   *
   */
  using real_type = typename numbers::NumberTraits<number>::real_type;

  /**
   * 一个迭代器类的类型定义，在这个矩阵的所有非零项上行走。这个迭代器不能改变矩阵的值。
   *
   */
  using const_iterator = SparseMatrixIterators::Iterator<number, true>;

  /**
   * 走过该矩阵所有非零项的迭代器类的类型定义。这个迭代器
   * @em
   * 可以改变矩阵的值，但当然不能改变稀疏模式，因为一旦稀疏矩阵被附加到它上面，这个模式就固定了。
   *
   */
  using iterator = SparseMatrixIterators::Iterator<number, false>;

  /**
   * 一个描述这个类在运行时行为方面的一些特征的结构。其他一些以一个或其他矩阵类作为模板参数的类（如块状矩阵类）可以根据这个类中的变量来调整其行为。
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
   * @name  构造函数和初始化
   *
   */
  //@{
  /**
   * 构造函数；将矩阵初始化为空，没有任何结构，也就是说，矩阵根本无法使用。因此，这个构造函数只对作为类的成员的矩阵有用。所有其他的矩阵都应该在数据流中的一个点上创建，在那里所有必要的信息都是可用的。
   * 你必须在使用前用reinit(const SparsityPattern&)初始化矩阵。
   *
   */
  SparseMatrix();

  /**
   * 复制构造函数。只有当要复制的矩阵为空时，才允许调用该构造函数。这与SparsityPattern的原因相同，详见那里。
   * 如果你真的想复制一个完整的矩阵，你可以通过使用copy_from()函数来实现。
   *
   */
  SparseMatrix(const SparseMatrix &);

  /**
   * 移动构造函数。通过将矩阵的内部数据 @p m
   * 转移到一个新的对象中，构造一个新的稀疏矩阵。
   * 移动构造允许一个对象从一个函数中返回或打包成一个元组，即使该类不能被复制构造。
   *
   */
  SparseMatrix(SparseMatrix<number> &&m) noexcept;

  /**
   * 构造器。接受给定的矩阵稀疏度结构来表示这个矩阵的稀疏度模式。你可以在以后通过调用reinit(const
   * SparsityPattern&)函数来改变稀疏性模式。
   * 你必须确保稀疏结构的寿命至少与该矩阵的寿命一样长，或者只要reinit(const
   * SparsityPattern&)没有被调用新的稀疏模式。
   * 构造函数被明确标记，以便不允许有人将稀疏模式代替稀疏矩阵传递给某个函数，这样就会生成一个空矩阵。
   *
   */
  explicit SparseMatrix(const SparsityPattern &sparsity);

  /**
   * 拷贝构造函数：用身份矩阵初始化矩阵。如果稀疏模式和身份矩阵的大小不一致，或者如果稀疏模式没有在整个对角线上提供非零条目，这个构造函数将抛出一个异常。
   *
   */
  SparseMatrix(const SparsityPattern &sparsity, const IdentityMatrix &id);

  /**
   * 解构器。释放所有内存，但不释放稀疏结构的内存。
   *
   */
  virtual ~SparseMatrix() override;

  /**
   * 复制操作符。由于复制整个稀疏矩阵是一个非常昂贵的操作，我们不允许这样做，除了大小为0的空矩阵这一特殊情况。这看起来不是特别有用，但如果想有一个
   * <code>std::vector@<SparseMatrix@<double@> @></code>
   * ，这正是人们所需要的：在这种情况下，人们可以创建一个空矩阵的向量（需要复制对象的能力），然后用有用的东西来填充。
   *
   */
  SparseMatrix<number> &
  operator=(const SparseMatrix<number> &);

  /**
   * 移动赋值运算符。这个操作符通过转移 @p m.
   * 的内部数据，将目前的矩阵替换成 @p m 。
   *
   */
  SparseMatrix<number> &
  operator=(SparseMatrix<number> &&m) noexcept;

  /**
   * 复制运算符：用身份矩阵初始化矩阵。如果稀疏模式和身份矩阵的大小不一致，或者如果稀疏模式没有规定整个对角线上的非零条目，这个操作符将抛出一个异常。
   *
   */
  SparseMatrix<number> &
  operator=(const IdentityMatrix &id);

  /**
   * 这个操作符将一个标量分配给一个矩阵。因为这通常没有什么意义（我们应该把所有的矩阵条目都设置为这个值吗？
   * 仅仅是稀疏模式的非零条目？），这个操作只允许在实际要分配的值为零时进行。这个操作符的存在只是为了允许明显的符号<tt>matrix=0</tt>，它将矩阵的所有元素设置为零，但保留之前使用的稀疏模式。
   * @dealiiOperationIsMultithreaded
   *
   */
  SparseMatrix &
  operator=(const double d);

  /**
   * 用给定的稀疏模式重新初始化稀疏矩阵。后者告诉矩阵需要保留多少个非零元素。
   * 关于内存分配，和上面说的一样。
   * 你必须确保稀疏结构的寿命至少与该矩阵的寿命一样长，或者只要reinit(const
   * SparsityPattern &)没有被调用，就不会有新的稀疏结构。
   * 矩阵的元素被这个函数设置为零。
   *
   */
  virtual void
  reinit(const SparsityPattern &sparsity);

  /**
   * 释放所有内存并返回到与调用默认构造函数后相同的状态。它也会忘记它之前绑定的稀疏模式。
   *
   */
  virtual void
  clear();
  //@}
  /**
   * @name  矩阵的信息
   *
   */
  //@{
  /**
   * 返回该对象是否为空。如果两个维度都是零或者没有关联的SparsityPattern，它就是空的。
   *
   */
  bool
  empty() const;

  /**
   * 返回共域（或范围）空间的维度。注意，矩阵的维度是
   * $m \times n$  。
   *
   */
  size_type
  m() const;

  /**
   * 返回域空间的维度。请注意，矩阵的维度是 $m \times n$  .
   *
   */
  size_type
  n() const;

  /**
   * 返回特定行中的条目数。
   *
   */
  size_type
  get_row_length(const size_type row) const;

  /**
   * 返回该矩阵的非零元素的数量。实际上，它返回的是稀疏模式中的条目数；如果任何一个条目恰好是零，无论如何都会被计算在内。
   *
   */
  std::size_t
  n_nonzero_elements() const;

  /**
   * 返回这个矩阵中实际非零元素的数量。可以指定参数<tt>threshold</tt>，以便只计算绝对值大于阈值的元素。
   * 注意，这个函数（与n_nonzero_elements()相反）不计算稀疏模式的所有条目，而只计算非零的（或绝对值大于阈值的）。
   *
   */
  std::size_t
  n_actually_nonzero_elements(const double threshold = 0.) const;

  /**
   * 返回一个对该矩阵底层稀疏性模式的（常量）引用。
   * 尽管返回值被声明为<tt>const</tt>，但你应该注意，如果你调用任何对其进行操作的对象的非常量函数，它可能会改变。
   *
   */
  const SparsityPattern &
  get_sparsity_pattern() const;

  /**
   * 确定此对象的内存消耗（以字节为单位）的估计值。参见MemoryConsumption。
   *
   */
  std::size_t
  memory_consumption() const;

  /**
   * 虚函数，与分布式并行矩阵兼容。
   *
   */
  void compress(::dealii::VectorOperation::values);

  //@}
  /**
   * @name  修改条目
   *
   */
  //@{
  /**
   * 设置元素（<i>i,j</i>）为<tt>value</tt>。如果条目不存在或者<tt>value</tt>不是一个有限的数字，则抛出一个错误。尽管如此，它仍然允许在不存在的字段中存储零值。
   *
   */
  void
  set(const size_type i, const size_type j, const number value);

  /**
   * 将FullMatrix中给出的所有元素设置到<tt>indices</tt>给出的稀疏矩阵位置。换句话说，这个函数将<tt>full_matrix</tt>中的元素写入调用的矩阵中，对矩阵的行和列都使用<tt>indices</tt>指定的本地到全球的索引。这个函数假设一个二次稀疏矩阵和一个二次全矩阵，这是FE计算中的通常情况。
   * 可选参数<tt>elide_zero_values</tt>可以用来指定是无论如何都要设置零值，还是要过滤掉零值（如果存在的话，不改变相应元素中的先前内容）。默认值是<tt>false</tt>，也就是说，即使是零值也要处理。
   *
   */
  template <typename number2>
  void
  set(const std::vector<size_type> &indices,
      const FullMatrix<number2> &   full_matrix,
      const bool                    elide_zero_values = false);

  /**
   * 与之前的函数相同，但现在包括了使用矩形full_matrices的可能性，以及在行和列上分别使用不同的本地到全球索引。
   *
   */
  template <typename number2>
  void
  set(const std::vector<size_type> &row_indices,
      const std::vector<size_type> &col_indices,
      const FullMatrix<number2> &   full_matrix,
      const bool                    elide_zero_values = false);

  /**
   * 将矩阵的指定行中的几个元素与<tt>col_indices</tt>给出的列索引设置为相应的值。
   * 可选的参数<tt>elide_zero_values</tt>可以用来指定是无论如何都要设置零值，还是要过滤掉零值（如果存在的话，不改变相应元素中的先前内容）。默认值是<tt>false</tt>，也就是说，即使是零值也要处理。
   *
   */
  template <typename number2>
  void
  set(const size_type               row,
      const std::vector<size_type> &col_indices,
      const std::vector<number2> &  values,
      const bool                    elide_zero_values = false);

  /**
   * 将几个元素设置为由<tt>values</tt>给出的值，在给定的行和col_indices给出的列中设置为稀疏矩阵。
   * 可选的参数<tt>elide_zero_values</tt>可以用来指定是无论如何都要插入零值还是要过滤掉它们。默认值是<tt>false</tt>，也就是说，即使是零值也要插入/替换。
   *
   */
  template <typename number2>
  void
  set(const size_type  row,
      const size_type  n_cols,
      const size_type *col_indices,
      const number2 *  values,
      const bool       elide_zero_values = false);

  /**
   * 向元素添加<tt>value</tt>（<i>i,j</i>）。
   * 如果该条目不存在或者<tt>value</tt>不是一个有限的数字，则抛出一个错误。尽管如此，它仍然允许在不存在的字段中存储零值。
   *
   */
  void
  add(const size_type i, const size_type j, const number value);

  /**
   * 将FullMatrix<double>中给出的所有元素添加到由<tt>indices</tt>给出的稀疏矩阵位置。换句话说，这个函数将<tt>full_matrix</tt>中的元素添加到调用矩阵的相应条目中，使用<tt>indices</tt>为矩阵的行和列指定的本地到全球索引。这个函数假定一个二次稀疏矩阵和一个二次全矩阵，这是FE计算中通常的情况。
   * 可选参数<tt>elide_zero_values</tt>可以用来指定是无论如何都要添加零值，还是要过滤掉这些零值，只添加非零数据。默认值是<tt>true</tt>，也就是说，零值不会被添加到矩阵中。
   *
   */
  template <typename number2>
  void
  add(const std::vector<size_type> &indices,
      const FullMatrix<number2> &   full_matrix,
      const bool                    elide_zero_values = true);

  /**
   * 与之前的函数相同，但现在包括了使用矩形full_matrices的可能性，以及在行和列上分别使用不同的本地到全球索引。
   *
   */
  template <typename number2>
  void
  add(const std::vector<size_type> &row_indices,
      const std::vector<size_type> &col_indices,
      const FullMatrix<number2> &   full_matrix,
      const bool                    elide_zero_values = true);

  /**
   * 将矩阵的指定行中的几个元素与<tt>col_indices</tt>给出的列索引设置为相应的值。
   * 可选参数<tt>elide_zero_values</tt>可以用来指定是无论如何都要添加零值，还是要过滤掉这些零值，只添加非零数据。默认值是<tt>true</tt>，也就是说，零值不会被添加到矩阵中。
   *
   */
  template <typename number2>
  void
  add(const size_type               row,
      const std::vector<size_type> &col_indices,
      const std::vector<number2> &  values,
      const bool                    elide_zero_values = true);

  /**
   * 在给定的全局矩阵行中，在稀疏矩阵中由col_indices指定的列中添加一个由<tt>values</tt>给出的数值阵列。
   * 可选的参数<tt>elide_zero_values</tt>可以用来指定是无论如何都要添加零值，还是要过滤掉这些数据，只添加非零值。默认值是<tt>true</tt>，也就是说，零值不会被添加到矩阵中。
   *
   */
  template <typename number2>
  void
  add(const size_type  row,
      const size_type  n_cols,
      const size_type *col_indices,
      const number2 *  values,
      const bool       elide_zero_values      = true,
      const bool       col_indices_are_sorted = false);

  /**
   * 将整个矩阵乘以一个固定系数。
   *
   */
  SparseMatrix &
  operator*=(const number factor);

  /**
   * 用整个矩阵除以一个固定系数。
   *
   */
  SparseMatrix &
  operator/=(const number factor);

  /**
   * 通过形成现有矩阵和其转置之间的平均值来对称矩阵，
   * $A = \frac 12(A+A^T)$  。
   * 这个操作假设底层的稀疏模式代表一个对称的对象。如果不是这样，那么这个操作的结果将不是一个对称矩阵，因为出于效率的考虑，它只通过在左下角的三角形部分进行循环来明确地进行对称；如果右上角的三角形有条目，那么这些元素在对称过程中会被遗漏。稀疏模式的对称化可以通过
   * SparsityPattern::symmetrize(). 得到。
   *
   */
  void
  symmetrize();

  /**
   * 将作为参数给出的矩阵复制到当前对象中。
   * 复制矩阵是一个昂贵的操作，我们不希望通过编译器生成的代码
   * <code>operator=</code>
   * 而意外发生（例如，如果不小心声明了当前类型的函数参数为<i>by
   * value</i>而不是<i>by
   * reference</i>，就会发生这种情况）。复制矩阵的功能是在这个成员函数中实现的。因此，该类型对象的所有复制操作都需要一个明确的函数调用。
   * 源矩阵可以是一个任意类型的矩阵，只要其数据类型可以转换为该矩阵的数据类型。
   * 该函数返回一个对<tt>*this</tt>的引用。
   *
   */
  template <typename somenumber>
  SparseMatrix<number> &
  copy_from(const SparseMatrix<somenumber> &source);

  /**
   * 这个函数完全类似于 SparsityPattern::copy_from()
   * 函数，它允许在一个步骤中初始化整个矩阵。关于参数类型及其含义的更多信息请见那里。你还可以在那里找到一个关于如何使用这个函数的小例子。
   * 与引用的函数唯一不同的是，内部迭代器指向的对象需要是
   * <tt>std::pair<unsigned  int,
   * value</tt>类型，其中<tt>value</tt>需要可转换为该类的元素类型，由<tt>number</tt>模板参数指定。
   * 矩阵以前的内容被覆盖。注意，由输入参数指定的条目不一定要覆盖矩阵的所有元素。未覆盖的元素保持不动。
   *
   */
  template <typename ForwardIterator>
  void
  copy_from(const ForwardIterator begin, const ForwardIterator end);

  /**
   * 将一个完整矩阵的非零条目复制到此对象中。之前的内容被删除。
   * 请注意，底层的稀疏模式必须适合容纳全矩阵的非零条目。这可以使用
   * SparsityPattern::copy_from()
   * 的那个版本来实现，该版本以FullMatrix作为参数。
   *
   */
  template <typename somenumber>
  void
  copy_from(const FullMatrix<somenumber> &matrix);

#  ifdef DEAL_II_WITH_TRILINOS
  /**
   * 将给定的特里诺斯矩阵复制到这个矩阵中。如果当前对象的稀疏模式不包含给定参数的非零条目的位置，该操作会触发一个断言。
   * 这个函数假设两个矩阵有相同的大小。
   * 该函数返回一个对<tt>*this</tt>的引用。
   *
   */
  SparseMatrix<number> &
  copy_from(const TrilinosWrappers::SparseMatrix &matrix);
#  endif

  /**
   * 将<tt>matrix</tt>按<tt>factor</tt>的比例添加到这个矩阵中，也就是说，将<tt>factor*matrix</tt>的矩阵添加到<tt>this</tt>。如果所涉及的两个矩阵的稀疏性模式不指向同一个对象，这个函数会抛出一个错误，因为在这种情况下，操作会比较便宜。
   * 源矩阵可以是一个任意底层标量类型的稀疏矩阵，只要其数据类型可以转换为这个矩阵的数据类型。
   *
   */
  template <typename somenumber>
  void
  add(const number factor, const SparseMatrix<somenumber> &matrix);

  //@}
  /**
   * @name 条目访问
   *
   */
  //@{

  /**
   * 返回条目的值（<i>i,j</i>）。
   * 这可能是一个昂贵的操作，你应该始终注意在哪里调用这个函数。为了避免滥用，如果所需元素在矩阵中不存在，该函数会抛出一个异常。
   * 如果你想要一个返回零的函数（对于不在矩阵的稀疏模式中的条目），请使用el()函数。
   * 如果你要在所有元素上循环，可以考虑使用一个迭代器类来代替，因为它们更适合稀疏的矩阵结构。
   *
   */
  const number &
  operator()(const size_type i, const size_type j) const;

  /**
   * 与上面那个相反，这个函数允许修改对象。
   *
   */
  number &
  operator()(const size_type i, const size_type j);

  /**
   * 这个函数主要像operator()()，它返回矩阵条目的值（<i>i,j</i>）。唯一的区别是，如果这个条目不存在于稀疏模式中，那么就不会引发异常，而是返回0。虽然这在某些情况下可能很方便，但请注意，由于没有使用矩阵的稀疏性，所以很容易写出与最优解相比很慢的算法。
   * 如果你要在所有元素上循环，可以考虑使用一个迭代器类来代替，因为它们更适合稀疏的矩阵结构。
   *
   */
  number
  el(const size_type i, const size_type j) const;

  /**
   * 返回第<i>i</i>行中的主对角线元素。如果矩阵不是二次方的，这个函数会抛出一个错误。
   * 这个函数比operator()()快得多，因为对于二次矩阵来说，对角线条目可能是每行中第一个被存储的，因此访问时不需要搜索正确的列号。
   *
   */
  number
  diag_element(const size_type i) const;

  /**
   * 和上面一样，但返回一个可写的引用。你确定你知道你在做什么吗？
   *
   */
  number &
  diag_element(const size_type i);

  //@}
  /**
   * @name  乘法运算
   *
   */
  //@{
  /**
   * 矩阵-向量乘法：让<i>dst =
   * M*src</i>与<i>M</i>是这个矩阵。
   * 注意，虽然这个函数可以对所有提供迭代器类的向量进行操作，但它只对类型为
   * @ref Vector  的对象真正有效。
   * 对于所有迭代元素或随机成员访问昂贵的类来说，这个函数并不高效。特别是，如果你想与BlockVector对象相乘，你应该考虑同时使用BlockSparseMatrix。
   * 源和目的不能是同一个向量。
   * @dealiiOperationIsMultithreaded
   *
   */
  template <class OutVector, class InVector>
  void
  vmult(OutVector &dst, const InVector &src) const;

  /**
   * 矩阵-向量乘法：让<i>dst =
   * M<sup>T</sup>*src</i>与<i>M</i>是这个矩阵。这个函数与vmult()的作用相同，但需要转置的矩阵。
   * 注意，虽然这个函数可以对所有提供迭代器类的向量进行操作，但它只对类型为
   * @ref Vector  的对象真正有效。
   * 对于所有迭代元素或随机成员访问昂贵的类来说，这个函数并不高效。特别是，如果你想与BlockVector对象相乘，你应该考虑同时使用BlockSparseMatrix。
   * 源和目的不能是同一个向量。
   *
   */
  template <class OutVector, class InVector>
  void
  Tvmult(OutVector &dst, const InVector &src) const;

  /**
   * 添加矩阵-向量的乘法。在<i>dst</i>上添加<i>M*src</i>，<i>M</i>是这个矩阵。
   * 注意，虽然这个函数可以对所有提供迭代器类的向量进行操作，但它只对类型为
   * @ref Vector  的对象真正有效。
   * 对于所有迭代元素或随机成员访问昂贵的类来说，这个函数并不高效。特别是，如果你想与BlockVector对象相乘，你应该考虑同时使用BlockSparseMatrix。
   * 源和目的不能是同一个向量。
   * @dealiiOperationIsMultithreaded
   *
   */
  template <class OutVector, class InVector>
  void
  vmult_add(OutVector &dst, const InVector &src) const;

  /**
   * 添加矩阵-向量的乘法。将<i>M<sup>T</sup>*src</i>加到<i>dst</i>，<i>M</i>是这个矩阵。这个函数与vmult_add()的操作相同，但取的是转置的矩阵。
   * 注意，虽然这个函数可以对所有提供迭代器类的向量进行操作，但它只对类型为
   * @ref Vector  的对象真正有效。
   * 对于所有迭代元素或随机成员访问昂贵的类来说，这个函数并不高效。特别是，如果你想与BlockVector对象相乘，你应该考虑同时使用BlockSparseMatrix。
   * 源和目的不能是同一个向量。
   *
   */
  template <class OutVector, class InVector>
  void
  Tvmult_add(OutVector &dst, const InVector &src) const;

  /**
   * 返回向量  $v$  相对于该矩阵诱导的法线的平方，即
   * $\left(v,Mv\right)$
   * 。这很有用，例如在有限元背景下，一个函数的 $L_2$
   * 规范等于相对于代表有限元函数节点值的向量的质量矩阵的矩阵规范。
   * 显然，对于这个操作来说，矩阵需要是二次的，而且为了使结果真正成为一个规范，它还需要是实数对称的或复数隐式的。
   * 该矩阵和给定向量的底层模板类型应该都是实值或复值，但不是混合的，这样这个函数才有意义。
   * @dealiiOperationIsMultithreaded
   *
   */
  template <typename somenumber>
  somenumber
  matrix_norm_square(const Vector<somenumber> &v) const;

  /**
   * 计算矩阵标量乘积  $\left(u,Mv\right)$  。
   * @dealiiOperationIsMultithreaded
   *
   */
  template <typename somenumber>
  somenumber
  matrix_scalar_product(const Vector<somenumber> &u,
                        const Vector<somenumber> &v) const;

  /**
   * 计算方程<i>Mx=b</i>的残差，其中残差被定义为<i>r=b-Mx</i>。将残差写入<tt>dst</tt>。残差向量的<i>l<sub>2</sub></i>准则被返回。
   * 源<i>x</i>和目的<i>dst</i>不能是同一个向量。
   * @dealiiOperationIsMultithreaded
   *
   */
  template <typename somenumber>
  somenumber
  residual(Vector<somenumber> &      dst,
           const Vector<somenumber> &x,
           const Vector<somenumber> &b) const;

  /**
   * 执行矩阵-矩阵乘法<tt>C = A
   * B</tt>，或者，如果给出一个可选的矢量参数，则<tt>C = A
   * diag(V)
   * B</tt>，其中<tt>diag(V)</tt>定义了一个带有矢量项的对角矩阵。
   * 这个函数假定调用矩阵 @p A 和参数 @p B
   * 的大小兼容。默认情况下，输出矩阵 @p C
   * 将被适当调整大小。    默认情况下，即如果可选的参数
   * @p rebuild_sparsity_pattern 是 @p true,
   * ，矩阵C的稀疏模式将被改变，以确保由乘积 $AB$
   * 产生的所有条目可以被存储在 $C$
   * 中。这是一个昂贵的操作，如果有办法预先预测稀疏模式，你可能应该在以
   * @p false
   * 为最后一个参数调用这个函数之前自己建立它。在这种情况下，疏密模式的重建被绕过了。
   * 当把 @p rebuild_sparsity_pattern 设置为 @p true
   * 时（即把它留在默认值），必须意识到作为第一个参数传递的矩阵
   * @p C
   * 仍然需要用稀疏模式进行初始化（可以在创建稀疏矩阵对象时，或者通过
   * SparseMatrix::reinit()
   * 函数）。这是因为我们可以在当前函数中创建一个稀疏模式，然后将
   * @p C
   * 与之关联，但一旦当前函数结束，就没有办法将这个稀疏模式的所有权转移给任何人。因此，该函数要求
   * @p C
   * 已经与一个稀疏模式对象相关联，然后该对象被重置为适合
   * @p A 和 @p B. 的乘积。 然而，作为其结果，还必须认识到
   * @p C
   * 的稀疏模式被修改，这将使碰巧<i>also</i>使用该稀疏模式对象的<i>all
   * other SparseMatrix objects</i>无效。
   *
   */
  template <typename numberB, typename numberC>
  void
  mmult(SparseMatrix<numberC> &      C,
        const SparseMatrix<numberB> &B,
        const Vector<number> &       V = Vector<number>(),
        const bool                   rebuild_sparsity_pattern = true) const;

  /**
   * 用<tt>this</tt>的转置执行矩阵-矩阵乘法，即<tt>C =
   * A<sup>T</sup>
   * B</tt>，或者，如果给出了可选的矢量参数，<tt>C =
   * A<sup>T</sup> diag(V)
   * B</tt>，其中<tt>diag(V)</tt>定义了一个带有矢量项的对角矩阵。
   * 这个函数假定调用矩阵<tt>A</tt>和<tt>B</tt>的大小兼容。<tt>C</tt>的大小将在本函数中设置。
   * 矩阵C的内容和稀疏模式将被这个函数改变，所以要确保稀疏模式没有在你的程序中其他地方使用。这是一个昂贵的操作，所以在你使用这个函数之前要三思而后行。
   * 有一个可选的标志<tt>rebuild_sparsity_pattern</tt>，可以用来绕过创建一个新的稀疏度模式，而使用存储在<tt>C</tt>中的稀疏度模式。在这种情况下，要确保它真的适合。默认情况下是重建稀疏度模式。
   * @note
   * 重建稀疏度模式需要改变它。这意味着所有与该稀疏性模式相关的其他矩阵将有无效的条目。
   *
   */
  template <typename numberB, typename numberC>
  void
  Tmmult(SparseMatrix<numberC> &      C,
         const SparseMatrix<numberB> &B,
         const Vector<number> &       V = Vector<number>(),
         const bool                   rebuild_sparsity_pattern = true) const;

  //@}
  /**
   * @name  矩阵规范
   *
   */
  //@{

  /**
   * 返回矩阵的 $l_1$ 规范，即 $|M|_1=\max_{\mathrm{all\ columns\
   * }j}\sum_{\mathrm{all\ rows\ } i} |M_{ij}|$
   * ，（最大列数之和）。 这是自然的矩阵准则，与向量的
   * $l_1$ 准则兼容，即 $|Mv|_1\leq |M|_1 |v|_1$  。(参见Haemmerlin-
   * Hoffmann: Numerische Mathematik)
   *
   */
  real_type
  l1_norm() const;

  /**
   * 返回矩阵的 $l_\infty$ 准则，即 $|M|_\infty=\max_{\mathrm{all\
   * rows\ }i}\sum_{\mathrm{all\ columns\ }j} |M_{ij}|$  ,
   * (行的最大和)。 这是自然的矩阵准则，与向量的
   * $l_\infty$ 准则兼容，即 $|Mv|_\infty \leq |M|_\infty |v|_\infty$
   * 。 (参看Haemmerlin-Hoffmann: Numerische Mathematik)
   *
   */
  real_type
  linfty_norm() const;

  /**
   * 返回矩阵的frobenius
   * norm，即矩阵中所有条目的平方之和的平方根。
   *
   */
  real_type
  frobenius_norm() const;
  //@}
  /**
   * @name  预处理方法
   *
   */
  //@{

  /**
   * 应用雅可比预处理，将<tt>src</tt>向量的每个元素乘以各自对角线元素的逆值，并将结果乘以松弛因子<tt>omega</tt>。
   *
   */
  template <typename somenumber>
  void
  precondition_Jacobi(Vector<somenumber> &      dst,
                      const Vector<somenumber> &src,
                      const number              omega = 1.) const;

  /**
   * 对<tt>src</tt>应用SSOR预处理，阻尼<tt>omega</tt>。
   * 可选的参数<tt>pos_right_of_diagonal</tt>应该提供一个数组，其中每个条目指定全局非零点阵列中对角线的右边位置。
   *
   */
  template <typename somenumber>
  void
  precondition_SSOR(Vector<somenumber> &            dst,
                    const Vector<somenumber> &      src,
                    const number                    omega = 1.,
                    const std::vector<std::size_t> &pos_right_of_diagonal =
                      std::vector<std::size_t>()) const;

  /**
   * 将SOR预处理矩阵应用于<tt>src</tt>。
   *
   */
  template <typename somenumber>
  void
  precondition_SOR(Vector<somenumber> &      dst,
                   const Vector<somenumber> &src,
                   const number              om = 1.) const;

  /**
   * 对<tt>src</tt>应用转置的SOR预处理矩阵。
   *
   */
  template <typename somenumber>
  void
  precondition_TSOR(Vector<somenumber> &      dst,
                    const Vector<somenumber> &src,
                    const number              om = 1.) const;

  /**
   * 就地执行SSOR预处理。
   * 应用预处理矩阵而不复制到第二个向量。
   * <tt>omega</tt>是放松参数。
   *
   */
  template <typename somenumber>
  void
  SSOR(Vector<somenumber> &v, const number omega = 1.) const;

  /**
   * 就地执行SOR预处理。 <tt>omega</tt>是松弛参数。
   *
   */
  template <typename somenumber>
  void
  SOR(Vector<somenumber> &v, const number om = 1.) const;

  /**
   * 就地进行转置SOR预处理。 <tt>omega</tt>是松弛参数。
   *
   */
  template <typename somenumber>
  void
  TSOR(Vector<somenumber> &v, const number om = 1.) const;

  /**
   * 就地进行移置的SOR预处理。
   * 标准的SOR方法是按照<tt>permutation</tt>规定的顺序应用的，即首先是行<tt>permutation[0]</tt>，然后是<tt>permutation[1]</tt>，依此类推。出于效率的考虑，需要排列组合以及它的逆向排列。
   * <tt>omega</tt>是放松参数。
   *
   */
  template <typename somenumber>
  void
  PSOR(Vector<somenumber> &          v,
       const std::vector<size_type> &permutation,
       const std::vector<size_type> &inverse_permutation,
       const number                  om = 1.) const;

  /**
   * 就地进行转置的包络SOR预处理。
   * 转置的SOR方法按照<tt>permutation</tt>规定的顺序应用，即首先是行<tt>permutation[m()-1]</tt>，然后是<tt>permutation[m()-2]</tt>，依此类推。出于效率的考虑，需要用到permutation以及它的逆向。
   * <tt>omega</tt>是放松参数。
   *
   */
  template <typename somenumber>
  void
  TPSOR(Vector<somenumber> &          v,
        const std::vector<size_type> &permutation,
        const std::vector<size_type> &inverse_permutation,
        const number                  om = 1.) const;

  /**
   * 对<tt>v</tt>做一个雅可比步骤。
   * 对<tt>b</tt>做一个直接的雅可比步骤，右手边<tt>b</tt>。这个函数需要一个辅助向量，它从GrowingVectorMemory中获取。
   *
   */
  template <typename somenumber>
  void
  Jacobi_step(Vector<somenumber> &      v,
              const Vector<somenumber> &b,
              const number              om = 1.) const;

  /**
   * 对<tt>v</tt>做一个SOR步骤。
   * 对右手边的<tt>b</tt>直接执行SOR步骤。
   *
   */
  template <typename somenumber>
  void
  SOR_step(Vector<somenumber> &      v,
           const Vector<somenumber> &b,
           const number              om = 1.) const;

  /**
   * 对<tt>v</tt>做一个邻接的SOR步骤。
   * 对<tt>b</tt>做一个直接的TSOR步骤，右手边<tt>b</tt>。
   *
   */
  template <typename somenumber>
  void
  TSOR_step(Vector<somenumber> &      v,
            const Vector<somenumber> &b,
            const number              om = 1.) const;

  /**
   * 对<tt>v</tt>做一个SSOR步骤。
   * 通过在SOR之后执行TSOR，对右手边的<tt>b</tt>直接执行SSOR步骤。
   *
   */
  template <typename somenumber>
  void
  SSOR_step(Vector<somenumber> &      v,
            const Vector<somenumber> &b,
            const number              om = 1.) const;
  //@}
  /**
   * @name  迭代器
   *
   */
  //@{

  /**
   * 返回一个指向矩阵的第一个元素的迭代器。
   * 注意这个类的一般文档中关于元素访问顺序的讨论。
   *
   */
  const_iterator
  begin() const;

  /**
   * 像上面的函数一样，但是对于非恒定矩阵。
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
   * 注意，如果给定的行是空的，即不包含任何非零条目，那么这个函数返回的迭代器等于<tt>end(r)</tt>。在这种情况下，如果行
   * @p r
   * 和以下任何一行都不包含任何非零条目，则返回的迭代器可能无法被解除引用。
   *
   */
  const_iterator
  begin(const size_type r) const;

  /**
   * 像上面的函数一样，但对于非恒定矩阵。
   *
   */
  iterator
  begin(const size_type r);

  /**
   * 返回一个指向第 @p r 行最后一个元素的迭代器，如果 @p
   * r
   * 之后的行根本不包含任何条目，则指向整个稀疏模式的末端。
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
   * 打印矩阵到给定的流，使用格式<tt>(row,column)
   * value</tt>，即每行打印矩阵的一个非零条目。如果<tt>across</tt>为真，则在单行上打印所有条目，使用格式row,column:value。
   * 如果参数<tt>diagonal_first</tt>为真，则二次方矩阵的对角线元素在其行中首先打印，对应于内部存储方案。如果它是假的，一行中的元素将按升列顺序写入。
   *
   */
  template <class StreamType>
  void
  print(StreamType &out,
        const bool  across         = false,
        const bool  diagonal_first = true) const;

  /**
   * 以通常的格式打印矩阵，即作为矩阵，而不是作为非零元素的列表。为了提高可读性，不在矩阵中的元素显示为空白，而明确设置为零的矩阵元素则显示为空白。
   * 参数允许对输出格式进行灵活设置。
   * <tt>precision</tt>和<tt>scientific</tt>用于确定数字格式，其中<tt>scientific
   * = false</tt>表示固定点符号。
   * <tt>width</tt>的一个零条目使函数计算出一个宽度，但如果输出粗略的话，可以将其改为一个正值。
   * 此外，还可以指定一个空值的字符。
   * 最后，整个矩阵可以与一个共同的分母相乘，以产生更可读的输出，甚至是整数。
   * @attention
   * 如果应用于一个大的矩阵，这个函数可能会产生<b>large</b>量的输出!
   *
   */
  void
  print_formatted(std::ostream &     out,
                  const unsigned int precision   = 3,
                  const bool         scientific  = true,
                  const unsigned int width       = 0,
                  const char *       zero_string = " ",
                  const double       denominator = 1.) const;

  /**
   * 打印矩阵的实际模式。对于每个绝对值大于阈值的条目，打印一个'*'，对于每个较小的数值打印一个':'，对于每个未分配的条目打印一个'.'。
   *
   */
  void
  print_pattern(std::ostream &out, const double threshold = 0.) const;

  /**
   * 将矩阵打印到输出流 @p out 中，其格式可以被
   * numpy::readtxt(). 读取 要在python中加载矩阵，只需做<code>
   * [data, row, column] = numpy.loadtxt('my_matrix.txt') sparse_matrix =
   * scipy.sparse.csr_matrix((data, (row, column))   </code>
   *
   */
  void
  print_as_numpy_arrays(std::ostream &     out,
                        const unsigned int precision = 9) const;

  /**
   * 把这个对象的数据全部写到一个文件中。这是以二进制模式进行的，所以输出的数据既不能被人类阅读，也不能（可能）被其他使用不同操作系统的数字格式的计算机阅读。
   * 这个函数的目的是，如果你的内存不足，想在不同的程序之间进行交流，或者允许对象在程序的不同运行中持续存在，你可以把矩阵和稀疏模式换出来。
   *
   */
  void
  block_write(std::ostream &out) const;

  /**
   * 从文件中读取先前由block_write()写入的数据。
   * 这是用上述函数的逆运算来完成的，所以它的速度相当快，因为除了前面的几个数字，比特流是不被解释的。
   * 在这个操作中，对象被调整了大小，所有以前的内容都会丢失。然而，请注意，没有对新数据和底层的SparsityPattern对象是否适合在一起进行检查。你有责任确保稀疏度模式和要读取的数据是匹配的。
   * 一个原始形式的错误检查会被执行，它将识别最直白的尝试，即把一些数据解释为一个矩阵，以比特方式存储到一个实际上不是这样创建的文件中，但不是更多。
   *
   */
  void
  block_read(std::istream &in);
  //@}
  /**
   * @addtogroup  Exceptions  
     * @{ 
   *
   */

  /**
   * 异常情况
   *
   */
  DeclException2(ExcInvalidIndex,
                 int,
                 int,
                 << "You are trying to access the matrix entry with index <"
                 << arg1 << ',' << arg2
                 << ">, but this entry does not exist in the sparsity pattern "
                    "of this matrix."
                    "\n\n"
                    "The most common cause for this problem is that you used "
                    "a method to build the sparsity pattern that did not "
                    "(completely) take into account all of the entries you "
                    "will later try to write into. An example would be "
                    "building a sparsity pattern that does not include "
                    "the entries you will write into due to constraints "
                    "on degrees of freedom such as hanging nodes or periodic "
                    "boundary conditions. In such cases, building the "
                    "sparsity pattern will succeed, but you will get errors "
                    "such as the current one at one point or other when "
                    "trying to write into the entries of the matrix.");
  /**
   * 异常情况
   *
   */
  DeclExceptionMsg(ExcDifferentSparsityPatterns,
                   "When copying one sparse matrix into another, "
                   "or when adding one sparse matrix to another, "
                   "both matrices need to refer to the same "
                   "sparsity pattern.");
  /**
   * 异常情况
   *
   */
  DeclException2(ExcIteratorRange,
                 int,
                 int,
                 << "The iterators denote a range of " << arg1
                 << " elements, but the given number of rows was " << arg2);
  /**
   * 异常情况
   *
   */
  DeclExceptionMsg(ExcSourceEqualsDestination,
                   "You are attempting an operation on two matrices that "
                   "are the same object, but the operation requires that the "
                   "two objects are in fact different.");
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
   * 指向该矩阵使用的稀疏模式的指针。为了保证它在使用中不被删除，我们使用SmartPointer类来订阅它。
   *
   */
  SmartPointer<const SparsityPattern, SparseMatrix<number>> cols;

  /**
   * 所有非零条目的数值数组。一个条目在矩阵中的位置，也就是这个数组中给定值的行号和列号，只能用稀疏模式来推导。同样的道理也适用于更常见的通过坐标寻找一个条目的操作。
   *
   */
  std::unique_ptr<number[]> val;

  /**
   * 拨出的#val的大小。如果在过去的某个时候，通过使用reinit()函数将具有较小尺寸的稀疏模式关联到该对象，从而减少了矩阵的尺寸，那么这个尺寸可能大于实际使用的部分。
   *
   */
  std::size_t max_len;

  // make all other sparse matrices friends
  template <typename somenumber>
  friend class SparseMatrix;
  template <typename somenumber>
  friend class SparseLUDecomposition;
  template <typename>
  friend class SparseILU;

  // To allow it calling private prepare_add() and prepare_set().
  template <typename>
  friend class BlockMatrixBase;

  // Also give access to internal details to the iterator/accessor classes.
  template <typename, bool>
  friend class SparseMatrixIterators::Iterator;
  template <typename, bool>
  friend class SparseMatrixIterators::Accessor;

#  ifdef DEAL_II_WITH_MPI
  // Give access to internal datastructures to perform MPI operations.
  template <typename Number>
  friend void
  Utilities::MPI::sum(const SparseMatrix<Number> &,
                      const MPI_Comm &,
                      SparseMatrix<Number> &);
#  endif
};

#  ifndef DOXYGEN
 /*---------------------- Inline functions -----------------------------------*/ 



template <typename number>
inline typename SparseMatrix<number>::size_type
SparseMatrix<number>::m() const
{
  Assert(cols != nullptr, ExcNeedsSparsityPattern());
  return cols->rows;
}


template <typename number>
inline typename SparseMatrix<number>::size_type
SparseMatrix<number>::n() const
{
  Assert(cols != nullptr, ExcNeedsSparsityPattern());
  return cols->cols;
}


// Inline the set() and add() functions, since they will be called frequently.
template <typename number>
inline void
SparseMatrix<number>::set(const size_type i,
                          const size_type j,
                          const number    value)
{
  AssertIsFinite(value);

  const size_type index = cols->operator()(i, j);

  // it is allowed to set elements of the matrix that are not part of the
  // sparsity pattern, if the value to which we set it is zero
  if (index == SparsityPattern::invalid_entry)
    {
      Assert((index != SparsityPattern::invalid_entry) || (value == number()),
             ExcInvalidIndex(i, j));
      return;
    }

  val[index] = value;
}



template <typename number>
template <typename number2>
inline void
SparseMatrix<number>::set(const std::vector<size_type> &indices,
                          const FullMatrix<number2> &   values,
                          const bool                    elide_zero_values)
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



template <typename number>
template <typename number2>
inline void
SparseMatrix<number>::set(const std::vector<size_type> &row_indices,
                          const std::vector<size_type> &col_indices,
                          const FullMatrix<number2> &   values,
                          const bool                    elide_zero_values)
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



template <typename number>
template <typename number2>
inline void
SparseMatrix<number>::set(const size_type               row,
                          const std::vector<size_type> &col_indices,
                          const std::vector<number2> &  values,
                          const bool                    elide_zero_values)
{
  Assert(col_indices.size() == values.size(),
         ExcDimensionMismatch(col_indices.size(), values.size()));

  set(row,
      col_indices.size(),
      col_indices.data(),
      values.data(),
      elide_zero_values);
}



template <typename number>
inline void
SparseMatrix<number>::add(const size_type i,
                          const size_type j,
                          const number    value)
{
  AssertIsFinite(value);

  if (value == number())
    return;

  const size_type index = cols->operator()(i, j);

  // it is allowed to add elements to the matrix that are not part of the
  // sparsity pattern, if the value to which we set it is zero
  if (index == SparsityPattern::invalid_entry)
    {
      Assert((index != SparsityPattern::invalid_entry) || (value == number()),
             ExcInvalidIndex(i, j));
      return;
    }

  val[index] += value;
}



template <typename number>
template <typename number2>
inline void
SparseMatrix<number>::add(const std::vector<size_type> &indices,
                          const FullMatrix<number2> &   values,
                          const bool                    elide_zero_values)
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



template <typename number>
template <typename number2>
inline void
SparseMatrix<number>::add(const std::vector<size_type> &row_indices,
                          const std::vector<size_type> &col_indices,
                          const FullMatrix<number2> &   values,
                          const bool                    elide_zero_values)
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



template <typename number>
template <typename number2>
inline void
SparseMatrix<number>::add(const size_type               row,
                          const std::vector<size_type> &col_indices,
                          const std::vector<number2> &  values,
                          const bool                    elide_zero_values)
{
  Assert(col_indices.size() == values.size(),
         ExcDimensionMismatch(col_indices.size(), values.size()));

  add(row,
      col_indices.size(),
      col_indices.data(),
      values.data(),
      elide_zero_values);
}



template <typename number>
inline SparseMatrix<number> &
SparseMatrix<number>::operator*=(const number factor)
{
  Assert(cols != nullptr, ExcNeedsSparsityPattern());
  Assert(val != nullptr, ExcNotInitialized());

  number *            val_ptr = val.get();
  const number *const end_ptr = val.get() + cols->n_nonzero_elements();

  while (val_ptr != end_ptr)
    *val_ptr++ *= factor;

  return *this;
}



template <typename number>
inline SparseMatrix<number> &
SparseMatrix<number>::operator/=(const number factor)
{
  Assert(cols != nullptr, ExcNeedsSparsityPattern());
  Assert(val != nullptr, ExcNotInitialized());
  Assert(factor != number(), ExcDivideByZero());

  const number factor_inv = number(1.) / factor;

  number *            val_ptr = val.get();
  const number *const end_ptr = val.get() + cols->n_nonzero_elements();

  while (val_ptr != end_ptr)
    *val_ptr++ *= factor_inv;

  return *this;
}



template <typename number>
inline const number &
SparseMatrix<number>::operator()(const size_type i, const size_type j) const
{
  Assert(cols != nullptr, ExcNeedsSparsityPattern());
  Assert(cols->operator()(i, j) != SparsityPattern::invalid_entry,
         ExcInvalidIndex(i, j));
  return val[cols->operator()(i, j)];
}



template <typename number>
inline number &
SparseMatrix<number>::operator()(const size_type i, const size_type j)
{
  Assert(cols != nullptr, ExcNeedsSparsityPattern());
  Assert(cols->operator()(i, j) != SparsityPattern::invalid_entry,
         ExcInvalidIndex(i, j));
  return val[cols->operator()(i, j)];
}



template <typename number>
inline number
SparseMatrix<number>::el(const size_type i, const size_type j) const
{
  Assert(cols != nullptr, ExcNeedsSparsityPattern());
  const size_type index = cols->operator()(i, j);

  if (index != SparsityPattern::invalid_entry)
    return val[index];
  else
    return 0;
}



template <typename number>
inline number
SparseMatrix<number>::diag_element(const size_type i) const
{
  Assert(cols != nullptr, ExcNeedsSparsityPattern());
  Assert(m() == n(), ExcNotQuadratic());
  AssertIndexRange(i, m());

  // Use that the first element in each row of a quadratic matrix is the main
  // diagonal
  return val[cols->rowstart[i]];
}



template <typename number>
inline number &
SparseMatrix<number>::diag_element(const size_type i)
{
  Assert(cols != nullptr, ExcNeedsSparsityPattern());
  Assert(m() == n(), ExcNotQuadratic());
  AssertIndexRange(i, m());

  // Use that the first element in each row of a quadratic matrix is the main
  // diagonal
  return val[cols->rowstart[i]];
}



template <typename number>
template <typename ForwardIterator>
void
SparseMatrix<number>::copy_from(const ForwardIterator begin,
                                const ForwardIterator end)
{
  Assert(static_cast<size_type>(std::distance(begin, end)) == m(),
         ExcIteratorRange(std::distance(begin, end), m()));

  // for use in the inner loop, we define an alias to the type of the inner
  // iterators
  using inner_iterator =
    typename std::iterator_traits<ForwardIterator>::value_type::const_iterator;
  size_type row = 0;
  for (ForwardIterator i = begin; i != end; ++i, ++row)
    {
      const inner_iterator end_of_row = i->end();
      for (inner_iterator j = i->begin(); j != end_of_row; ++j)
        // write entries
        set(row, j->first, j->second);
    };
}


//---------------------------------------------------------------------------


namespace SparseMatrixIterators
{
  template <typename number>
  inline Accessor<number, true>::Accessor(const MatrixType *matrix,
                                          const std::size_t index_within_matrix)
    : SparsityPatternIterators::Accessor(&matrix->get_sparsity_pattern(),
                                         index_within_matrix)
    , matrix(matrix)
  {}



  template <typename number>
  inline Accessor<number, true>::Accessor(const MatrixType *matrix)
    : SparsityPatternIterators::Accessor(&matrix->get_sparsity_pattern())
    , matrix(matrix)
  {}



  template <typename number>
  inline Accessor<number, true>::Accessor(
    const SparseMatrixIterators::Accessor<number, false> &a)
    : SparsityPatternIterators::Accessor(a)
    , matrix(&a.get_matrix())
  {}



  template <typename number>
  inline number
  Accessor<number, true>::value() const
  {
    AssertIndexRange(linear_index, matrix->n_nonzero_elements());
    return matrix->val[linear_index];
  }



  template <typename number>
  inline const typename Accessor<number, true>::MatrixType &
  Accessor<number, true>::get_matrix() const
  {
    return *matrix;
  }



  template <typename number>
  inline Accessor<number, false>::Reference::Reference(const Accessor *accessor,
                                                       const bool)
    : accessor(accessor)
  {}


  template <typename number>
  inline Accessor<number, false>::Reference::operator number() const
  {
    AssertIndexRange(accessor->linear_index,
                     accessor->matrix->n_nonzero_elements());
    return accessor->matrix->val[accessor->linear_index];
  }



  template <typename number>
  inline const typename Accessor<number, false>::Reference &
  Accessor<number, false>::Reference::operator=(const number n) const
  {
    AssertIndexRange(accessor->linear_index,
                     accessor->matrix->n_nonzero_elements());
    accessor->matrix->val[accessor->linear_index] = n;
    return *this;
  }



  template <typename number>
  inline const typename Accessor<number, false>::Reference &
  Accessor<number, false>::Reference::operator+=(const number n) const
  {
    AssertIndexRange(accessor->linear_index,
                     accessor->matrix->n_nonzero_elements());
    accessor->matrix->val[accessor->linear_index] += n;
    return *this;
  }



  template <typename number>
  inline const typename Accessor<number, false>::Reference &
  Accessor<number, false>::Reference::operator-=(const number n) const
  {
    AssertIndexRange(accessor->linear_index,
                     accessor->matrix->n_nonzero_elements());
    accessor->matrix->val[accessor->linear_index] -= n;
    return *this;
  }



  template <typename number>
  inline const typename Accessor<number, false>::Reference &
  Accessor<number, false>::Reference::operator*=(const number n) const
  {
    AssertIndexRange(accessor->linear_index,
                     accessor->matrix->n_nonzero_elements());
    accessor->matrix->val[accessor->linear_index] *= n;
    return *this;
  }



  template <typename number>
  inline const typename Accessor<number, false>::Reference &
  Accessor<number, false>::Reference::operator/=(const number n) const
  {
    AssertIndexRange(accessor->linear_index,
                     accessor->matrix->n_nonzero_elements());
    accessor->matrix->val[accessor->linear_index] /= n;
    return *this;
  }



  template <typename number>
  inline Accessor<number, false>::Accessor(MatrixType *      matrix,
                                           const std::size_t index)
    : SparsityPatternIterators::Accessor(&matrix->get_sparsity_pattern(), index)
    , matrix(matrix)
  {}



  template <typename number>
  inline Accessor<number, false>::Accessor(MatrixType *matrix)
    : SparsityPatternIterators::Accessor(&matrix->get_sparsity_pattern())
    , matrix(matrix)
  {}



  template <typename number>
  inline typename Accessor<number, false>::Reference
  Accessor<number, false>::value() const
  {
    return Reference(this, true);
  }



  template <typename number>
  inline typename Accessor<number, false>::MatrixType &
  Accessor<number, false>::get_matrix() const
  {
    return *matrix;
  }



  template <typename number, bool Constness>
  inline Iterator<number, Constness>::Iterator(MatrixType *      matrix,
                                               const std::size_t index)
    : accessor(matrix, index)
  {}



  template <typename number, bool Constness>
  inline Iterator<number, Constness>::Iterator(MatrixType *matrix)
    : accessor(matrix)
  {}



  template <typename number, bool Constness>
  inline Iterator<number, Constness>::Iterator(
    const SparseMatrixIterators::Iterator<number, false> &i)
    : accessor(*i)
  {}



  template <typename number, bool Constness>
  inline const Iterator<number, Constness> &
  Iterator<number, Constness>::
  operator=(const SparseMatrixIterators::Iterator<number, false> &i)
  {
    accessor = *i;
    return *this;
  }



  template <typename number, bool Constness>
  inline Iterator<number, Constness> &
  Iterator<number, Constness>::operator++()
  {
    accessor.advance();
    return *this;
  }


  template <typename number, bool Constness>
  inline Iterator<number, Constness>
  Iterator<number, Constness>::operator++(int)
  {
    const Iterator iter = *this;
    accessor.advance();
    return iter;
  }


  template <typename number, bool Constness>
  inline const Accessor<number, Constness> &Iterator<number, Constness>::
                                            operator*() const
  {
    return accessor;
  }


  template <typename number, bool Constness>
  inline const Accessor<number, Constness> *Iterator<number, Constness>::
                                            operator->() const
  {
    return &accessor;
  }


  template <typename number, bool Constness>
  inline bool
  Iterator<number, Constness>::operator==(const Iterator &other) const
  {
    return (accessor == other.accessor);
  }


  template <typename number, bool Constness>
  inline bool
  Iterator<number, Constness>::operator!=(const Iterator &other) const
  {
    return !(*this == other);
  }


  template <typename number, bool Constness>
  inline bool
  Iterator<number, Constness>::operator<(const Iterator &other) const
  {
    Assert(&accessor.get_matrix() == &other.accessor.get_matrix(),
           ExcInternalError());

    return (accessor < other.accessor);
  }


  template <typename number, bool Constness>
  inline bool
  Iterator<number, Constness>::operator>(const Iterator &other) const
  {
    return (other < *this);
  }


  template <typename number, bool Constness>
  inline int
  Iterator<number, Constness>::operator-(const Iterator &other) const
  {
    Assert(&accessor.get_matrix() == &other.accessor.get_matrix(),
           ExcInternalError());

    return (*this)->linear_index - other->linear_index;
  }



  template <typename number, bool Constness>
  inline Iterator<number, Constness>
  Iterator<number, Constness>::operator+(const size_type n) const
  {
    Iterator x = *this;
    for (size_type i = 0; i < n; ++i)
      ++x;

    return x;
  }

} // namespace SparseMatrixIterators



template <typename number>
inline typename SparseMatrix<number>::const_iterator
SparseMatrix<number>::begin() const
{
  return const_iterator(this, 0);
}


template <typename number>
inline typename SparseMatrix<number>::const_iterator
SparseMatrix<number>::end() const
{
  return const_iterator(this);
}


template <typename number>
inline typename SparseMatrix<number>::iterator
SparseMatrix<number>::begin()
{
  return iterator(this, 0);
}


template <typename number>
inline typename SparseMatrix<number>::iterator
SparseMatrix<number>::end()
{
  return iterator(this, cols->rowstart[cols->rows]);
}


template <typename number>
inline typename SparseMatrix<number>::const_iterator
SparseMatrix<number>::begin(const size_type r) const
{
  AssertIndexRange(r, m());

  return const_iterator(this, cols->rowstart[r]);
}



template <typename number>
inline typename SparseMatrix<number>::const_iterator
SparseMatrix<number>::end(const size_type r) const
{
  AssertIndexRange(r, m());

  return const_iterator(this, cols->rowstart[r + 1]);
}



template <typename number>
inline typename SparseMatrix<number>::iterator
SparseMatrix<number>::begin(const size_type r)
{
  AssertIndexRange(r, m());

  return iterator(this, cols->rowstart[r]);
}



template <typename number>
inline typename SparseMatrix<number>::iterator
SparseMatrix<number>::end(const size_type r)
{
  AssertIndexRange(r, m());

  return iterator(this, cols->rowstart[r + 1]);
}



template <typename number>
template <class StreamType>
inline void
SparseMatrix<number>::print(StreamType &out,
                            const bool  across,
                            const bool  diagonal_first) const
{
  Assert(cols != nullptr, ExcNeedsSparsityPattern());
  Assert(val != nullptr, ExcNotInitialized());

  bool   hanging_diagonal = false;
  number diagonal         = number();

  for (size_type i = 0; i < cols->rows; ++i)
    {
      for (size_type j = cols->rowstart[i]; j < cols->rowstart[i + 1]; ++j)
        {
          if (!diagonal_first && i == cols->colnums[j])
            {
              diagonal         = val[j];
              hanging_diagonal = true;
            }
          else
            {
              if (hanging_diagonal && cols->colnums[j] > i)
                {
                  if (across)
                    out << ' ' << i << ',' << i << ':' << diagonal;
                  else
                    out << '(' << i << ',' << i << ") " << diagonal
                        << std::endl;
                  hanging_diagonal = false;
                }
              if (across)
                out << ' ' << i << ',' << cols->colnums[j] << ':' << val[j];
              else
                out << "(" << i << "," << cols->colnums[j] << ") " << val[j]
                    << std::endl;
            }
        }
      if (hanging_diagonal)
        {
          if (across)
            out << ' ' << i << ',' << i << ':' << diagonal;
          else
            out << '(' << i << ',' << i << ") " << diagonal << std::endl;
          hanging_diagonal = false;
        }
    }
  if (across)
    out << std::endl;
}


template <typename number>
inline void
SparseMatrix<number>::prepare_add()
{
  // nothing to do here
}



template <typename number>
inline void
SparseMatrix<number>::prepare_set()
{
  // nothing to do here
}

#  endif // DOXYGEN


 /*----------------------------   sparse_matrix.h ---------------------------*/ 

DEAL_II_NAMESPACE_CLOSE

#endif
 /*----------------------------   sparse_matrix.h ---------------------------*/ 


