//include/deal.II-translator/lac/chunk_sparse_matrix_0.txt
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

#ifndef dealii_chunk_sparse_matrix_h
#  define dealii_chunk_sparse_matrix_h


#  include <deal.II/base/config.h>

#  include <deal.II/base/smartpointer.h>
#  include <deal.II/base/subscriptor.h>

#  include <deal.II/lac/chunk_sparsity_pattern.h>
#  include <deal.II/lac/exceptions.h>
#  include <deal.II/lac/identity_matrix.h>

#  include <memory>


DEAL_II_NAMESPACE_OPEN

// Forward declarations
#  ifndef DOXYGEN
template <typename number>
class Vector;
template <typename number>
class FullMatrix;
#  endif

/*!   @addtogroup  Matrix1  
     * @{ 

* 
*
*/

/**
 * 一个命名空间，我们在其中声明对稀疏矩阵元素的迭代器。
 *
 *
 */
namespace ChunkSparseMatrixIterators
{
  // forward declaration
  template <typename number, bool Constness>
  class Iterator;

  /**
   * 稀疏矩阵访问器的通用模板。第一个模板参数表示底层数字类型，第二个表示矩阵的常数。
   * 通用模板没有被实现，只有针对第二个模板参数的两个可能值的特殊化。因此，这里列出的接口只是作为提供的模板，因为doxygen并没有链接这些特殊化。
   *
   */
  template <typename number, bool Constness>
  class Accessor : public ChunkSparsityPatternIterators::Accessor
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
    const ChunkSparseMatrix<number> &
    get_matrix() const;
  };



  /**
   * 用于常数矩阵的访问器类，在const_iterators中使用。这个类建立在用于稀疏模式的访问器类的基础上，在所有非零条目上循环，只增加了访问器函数，以获得存储在某一位置的实际值。
   *
   */
  template <typename number>
  class Accessor<number, true> : public ChunkSparsityPatternIterators::Accessor
  {
  public:
    /**
     * 这里要使用的矩阵的类型（包括常数）的类型定义。
     *
     */
    using MatrixType = const ChunkSparseMatrix<number>;

    /**
     * 构造函数。
     *
     */
    Accessor(MatrixType *matrix, const unsigned int row);

    /**
     * 构造器。构建给定矩阵的终端访问器。
     *
     */
    Accessor(MatrixType *matrix);

    /**
     * 复制构造器，从一个非常量访问器到一个常量访问器。
     *
     */
    Accessor(const ChunkSparseMatrixIterators::Accessor<number, false> &a);

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
    using ChunkSparsityPatternIterators::Accessor::advance;

    // Make iterator class a friend.
    template <typename, bool>
    friend class Iterator;
  };


  /**
   * 用于非恒定矩阵的访问器类，在迭代器中使用。这个类建立在用于稀疏模式的访问器类的基础上，在所有非零条目上循环，只增加了访问器函数，以获得存储在某一位置的实际值。
   *
   */
  template <typename number>
  class Accessor<number, false> : public ChunkSparsityPatternIterators::Accessor
  {
  private:
    /**
     * 参考类。这是访问器类在你调用value()函数时返回的东西。引用的作用就像它是对矩阵条目实际值的引用一样，也就是说，你可以读写它，可以对它进行加法和乘法，等等，但是由于矩阵并没有给出这个矩阵条目的地址，我们必须通过函数来完成这一切。
     * 构造函数需要一个指向访问器对象的指针，描述它指向矩阵的哪个元素。当人们写像iterator->value()=0（而不是iterator->value()=0.0）这样的代码时，就会产生歧义，因为右边是一个整数，既可以转换为<tt>数字</tt>（即最常见的双数），也可以转换为另一个<tt>参考</tt>类型的对象。然后编译器抱怨说不知道该采取哪种转换。
     * 由于某些原因，添加另一个重载operator=(int)似乎并不能解决这个问题。然而，我们通过给Reference构造函数添加第二个假参数来避免这个问题，该参数未被使用，但可以确保没有第二个匹配的转换序列，并使用单参数的右侧。
     *
     */
    class Reference
    {
    public:
      /**
       * 构造函数。关于第二个参数，见一般的类文档。
       *
       */
      Reference(const Accessor *accessor, const bool dummy);

      /**
       * 对矩阵的数据类型的转换操作符。
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
    using MatrixType = ChunkSparseMatrix<number>;

    /**
     * 构造函数。
     *
     */
    Accessor(MatrixType *matrix, const unsigned int row);

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
    using ChunkSparsityPatternIterators::Accessor::advance;

    // Make iterator class a friend.
    template <typename, bool>
    friend class Iterator;
  };



  /**
   * 用于常数和非常数矩阵的迭代器。
   * 第一个模板参数表示底层数字类型，第二个表示矩阵的常数性。
   * 因为这个类有一个针对<tt>Constness=false</tt>的特殊化，这个类是用于常数矩阵的迭代器。
   *
   */
  template <typename number, bool Constness>
  class Iterator
  {
  public:
    /**
     * 我们要操作的矩阵类型（包括恒定性）的类型化定义。
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
     * 中创建一个迭代器，用于给定行和其中的索引。
     *
     */
    Iterator(MatrixType *matrix, const unsigned int row);

    /**
     * 构造函数。为给定的矩阵创建终端迭代器。
     *
     */
    Iterator(MatrixType *matrix);

    /**
     * 转换构造函数，从非常量迭代器到常量迭代器。
     *
     */
    Iterator(const ChunkSparseMatrixIterators::Iterator<number, false> &i);

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
    operator+(const unsigned int n) const;

  private:
    /**
     * 存储一个访问器类的对象。
     *
     */
    Accessor<number, Constness> accessor;
  };

} // namespace ChunkSparseMatrixIterators



/**
 * 稀疏矩阵。该类实现了在一个由SparsityPattern表示的稀疏矩阵的位置上存储数值的功能。稀疏模式和数值的分离是由于人们可以在这些位置上存储不同类型的数据元素而不需要SparsityPattern知道这一点，更重要的是人们可以将一个以上的矩阵与同一个稀疏模式联系起来。
 * 这个类的使用在 step-51 中得到了证明。
 *
 *
 * @note
 * 这个模板的实例化提供给<tt>  @<float@>  和  @<double@></tt>;
 * 其他人可以在应用程序中生成（见手册中的 @ref
 * Instantiations 部分）。
 *
 *
 */
template <typename number>
class ChunkSparseMatrix : public virtual Subscriptor
{
public:
  /**
   * 声明容器大小的类型。
   *
   */
  using size_type = types::global_dof_index;

  /**
   * 矩阵条目的类型。这个别名类似于标准库中容器的<tt>value_type</tt>。
   *
   */
  using value_type = number;

  /**
   * 声明一个类型，该类型持有与该类的模板参数相同精度的实值数。如果这个类的模板参数是一个实数数据类型，那么real_type就等于模板参数。
   * 如果模板参数是一个 std::complex
   * 类型，那么real_type等于复数的基础类型。
   * 这个别名被用来表示规范的返回类型。
   *
   */
  using real_type = typename numbers::NumberTraits<number>::real_type;

  /**
   * 一个迭代器类的类型定义，在这个矩阵的所有非零条目上行走。这个迭代器不能改变矩阵的值。
   *
   */
  using const_iterator = ChunkSparseMatrixIterators::Iterator<number, true>;

  /**
   * 遍历该矩阵所有非零项的迭代器类的类型定义。这个迭代器
   * @em
   * 可以改变矩阵的值，但当然不能改变稀疏模式，因为一旦稀疏矩阵被附加到它上面，这就固定了。
   *
   */
  using iterator = ChunkSparseMatrixIterators::Iterator<number, false>;

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
   * @name  构造器和初始化。
   *
   */
  //@{
  /**
   * 构造函数；将矩阵初始化为空，没有任何结构，也就是说，矩阵根本无法使用。因此，这个构造函数只对作为类的成员的矩阵有用。所有其他的矩阵都应该在数据流中的一个点上创建，在那里所有必要的信息都是可用的。
   * 你必须在使用前用reinit(const
   * ChunkSparsityPattern&)初始化矩阵。
   *
   */
  ChunkSparseMatrix();

  /**
   * 复制构造函数。只有当要复制的矩阵为空时，才允许调用这个构造函数。这与ChunkSparsityPattern的原因相同，详见那里。
   * 如果你真的想复制一个完整的矩阵，你可以使用copy_from()函数来实现。
   *
   */
  ChunkSparseMatrix(const ChunkSparseMatrix &);

  /**
   * 构造函数。使用给定的矩阵稀疏度结构来表示该矩阵的稀疏度模式。你可以在以后通过调用reinit(const
   * ChunkSparsityPattern&)函数来改变稀疏性模式。
   * 你必须确保稀疏度结构的寿命至少和这个矩阵的寿命一样长，或者只要reinit(const
   * ChunkSparsityPattern&)没有被调用新的稀疏度模式。
   * 构造函数被明确标记，以便不允许有人将稀疏模式代替稀疏矩阵传递给某个函数，这样就会生成一个空矩阵。
   *
   */
  explicit ChunkSparseMatrix(const ChunkSparsityPattern &sparsity);

  /**
   * 拷贝构造函数：用身份矩阵初始化矩阵。如果稀疏模式和身份矩阵的大小不一致，或者如果稀疏模式没有在整个对角线上提供非零条目，这个构造函数将抛出一个异常。
   *
   */
  ChunkSparseMatrix(const ChunkSparsityPattern &sparsity,
                    const IdentityMatrix &      id);

  /**
   * 解构器。释放所有内存，但不释放稀疏结构的内存。
   *
   */
  virtual ~ChunkSparseMatrix() override;

  /**
   * 复制操作符。由于复制整个稀疏矩阵是一个非常昂贵的操作，我们不允许这样做，除了大小为0的空矩阵这一特殊情况。这看起来不是特别有用，但如果想有一个
   * <code>std::vector@<ChunkSparseMatrix@<double@> @></code>
   * ，这正是人们所需要的：在这种情况下，人们可以创建一个空矩阵的向量（需要复制对象的能力），然后用有用的东西来填充。
   *
   */
  ChunkSparseMatrix<number> &
  operator=(const ChunkSparseMatrix<number> &);

  /**
   * 复制操作符：用身份矩阵来初始化矩阵。如果稀疏模式和身份矩阵的大小不一致，或者如果稀疏模式没有在整个对角线上提供非零条目，这个操作符将抛出一个异常。
   *
   */
  ChunkSparseMatrix<number> &
  operator=(const IdentityMatrix &id);

  /**
   * 这个操作符将一个标量分配给一个矩阵。因为这通常没有什么意义（我们应该把所有的矩阵条目都设置为这个值吗？
   * 仅仅是稀疏模式的非零条目？），这个操作只允许在实际要分配的值为零时进行。这个操作符的存在只是为了允许明显的符号<tt>matrix=0</tt>，它将矩阵的所有元素设置为零，但保留之前使用的稀疏模式。
   *
   */
  ChunkSparseMatrix &
  operator=(const double d);

  /**
   * 用给定的稀疏模式重新初始化稀疏矩阵。后者告诉矩阵需要保留多少个非零元素。
   * 关于内存分配，和上面说的一样。
   * 你必须确保稀疏结构的寿命至少与该矩阵的寿命一样长，或者只要reinit(const
   * ChunkSparsityPattern
   * &)没有被调用，就不会有新的稀疏结构。
   * 矩阵的元素被这个函数设置为零。
   *
   */
  virtual void
  reinit(const ChunkSparsityPattern &sparsity);

  /**
   * 释放所有内存并返回到与调用默认构造函数后一样的状态。它也会忘记它之前绑定的稀疏模式。
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
   * 返回该对象是否为空。如果两个维度都是零或者没有关联ChunkSparsityPattern，那么它就是空的。
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
   * 返回域空间的维度。注意，矩阵的维度是 $m \times n$  .
   *
   */
  size_type
  n() const;

  /**
   * 返回这个矩阵的非零元素的数量。实际上，它返回的是稀疏模式中的条目数；如果任何一个条目碰巧是零，无论如何都会被计算在内。
   *
   */
  size_type
  n_nonzero_elements() const;

  /**
   * 返回这个矩阵中实际非零元素的数量。
   * 注意，这个函数（与n_nonzero_elements()相反）不计算稀疏模式的所有条目，只计算非零的条目。
   *
   */
  size_type
  n_actually_nonzero_elements() const;

  /**
   * 返回一个对该矩阵底层稀疏性模式的（常数）引用。
   * 虽然返回值被声明为<tt>const</tt>，但你应该注意，如果你调用任何对其进行操作的对象的非常量函数，它可能会改变。
   *
   */
  const ChunkSparsityPattern &
  get_sparsity_pattern() const;

  /**
   * 确定此对象的内存消耗（以字节为单位）的估计值。参见MemoryConsumption。
   *
   */
  std::size_t
  memory_consumption() const;

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
   * 向元素添加<tt>value</tt>（<i>i,j</i>）。
   * 如果该条目不存在或者<tt>value</tt>不是一个有限的数字，则抛出一个错误。尽管如此，它仍然允许在不存在的字段中存储零值。
   *
   */
  void
  add(const size_type i, const size_type j, const number value);

  /**
   * 在给定的全局矩阵行中，在稀疏矩阵中由col_indices指定的列中添加一个由<tt>values</tt>给出的数值数组。
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
  ChunkSparseMatrix &
  operator*=(const number factor);

  /**
   * 用整个矩阵除以一个固定系数。
   *
   */
  ChunkSparseMatrix &
  operator/=(const number factor);

  /**
   * 通过形成现有矩阵和其转置之间的平均值来对称矩阵，
   * $A = \frac 12(A+A^T)$  。
   * 这个操作假设底层的稀疏模式代表一个对称的对象。如果不是这样，那么这个操作的结果将不是一个对称矩阵，因为出于效率的考虑，它只通过在左下角的三角形部分进行循环来明确地进行对称；如果右上角的三角形有条目，那么这些元素在对称过程中会被遗漏。稀疏模式的对称化可以通过
   * ChunkSparsityPattern::symmetrize(). 得到。
   *
   */
  void
  symmetrize();

  /**
   * 将作为参数给出的矩阵复制到当前对象中。
   * 复制矩阵是一个昂贵的操作，我们不希望通过编译器生成的代码意外发生
   * <code>operator=</code>
   * 。（例如，如果不小心声明了当前类型的函数参数<i>by
   * value</i>而不是<i>by
   * reference</i>，就会发生这种情况）。复制矩阵的功能是在这个成员函数中实现的。因此，该类型对象的所有复制操作都需要一个明确的函数调用。
   * 源矩阵可以是一个任意类型的矩阵，只要其数据类型可以转换为该矩阵的数据类型。
   * 该函数返回一个对<tt>*this</tt>的引用。
   *
   */
  template <typename somenumber>
  ChunkSparseMatrix<number> &
  copy_from(const ChunkSparseMatrix<somenumber> &source);

  /**
   * 这个函数完全类似于 ChunkSparsityPattern::copy_from()
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
   * 将一个完整矩阵的非零条目复制到此对象中。之前的内容被删除。请注意，底层的稀疏模式必须适合容纳全矩阵的非零条目。
   *
   */
  template <typename somenumber>
  void
  copy_from(const FullMatrix<somenumber> &matrix);

  /**
   * 将<tt>matrix</tt>按<tt>factor</tt>的比例添加到这个矩阵中，也就是说，矩阵<tt>factor*matrix</tt>被添加到<tt>this</tt>。如果所涉及的两个矩阵的稀疏性模式不指向同一个对象，这个函数会抛出一个错误，因为在这种情况下，操作会比较便宜。
   * 源矩阵可以是一个任意底层标量类型的稀疏矩阵，只要其数据类型可以转换为这个矩阵的数据类型。
   *
   */
  template <typename somenumber>
  void
  add(const number factor, const ChunkSparseMatrix<somenumber> &matrix);

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
  number
  operator()(const size_type i, const size_type j) const;

  /**
   * 这个函数主要像operator()()，它返回矩阵条目的值（<i>i,j</i>）。唯一的区别是，如果这个条目不存在于稀疏模式中，那么就不会引发异常，而是返回0。虽然这在某些情况下可能很方便，但请注意，由于没有使用矩阵的稀疏性，所以写出的算法与最优解相比很简单，很慢。
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
   * 提取给定矩阵行中数值和索引的副本。
   * 用户应该传递数组column_indices和values的长度，这就提供了一个检查我们是否写到未分配的内存的方法。这个方法是由Trilinos行矩阵中的一个类似方法激发的，与迭代器相比，这个方法可以更快地访问矩阵中的条目，而迭代器对这种矩阵类型来说是相当慢的。
   *
   */
  void
  extract_row_copy(const size_type row,
                   const size_type array_length,
                   size_type &     row_length,
                   size_type *     column_indices,
                   number *        values) const;

  //@}
  /**
   * @name  矩阵矢量乘法
   *
   */
  //@{
  /**
   * 矩阵向量乘法：让<i>dst =
   * M*src</i>与<i>M</i>是这个矩阵。
   * 注意，虽然这个函数可以对所有提供迭代器类的向量进行操作，但它只对类型为
   * @ref Vector  的对象真正有效。
   * 对于所有迭代元素或随机成员访问昂贵的类来说，这个函数并不高效。特别是，如果你想与BlockVector对象相乘，你应该考虑同时使用BlockChunkSparseMatrix。
   * 来源和目的地不能是同一个向量。
   *
   */
  template <class OutVector, class InVector>
  void
  vmult(OutVector &dst, const InVector &src) const;

  /**
   * 矩阵-向量乘法：让<i>dst =
   * M<sup>T</sup>*src</i>与<i>M</i>为这个矩阵。这个函数与vmult()的作用相同，但需要转置的矩阵。
   * 注意，虽然这个函数可以对所有提供迭代器类的向量进行操作，但它只对类型为
   * @ref Vector  的对象真正有效。
   * 对于所有迭代元素或随机成员访问昂贵的类来说，这个函数并不高效。特别是，如果你想与BlockVector对象相乘，你应该考虑同时使用BlockChunkSparseMatrix。
   * 来源和目的地不能是同一个向量。
   *
   */
  template <class OutVector, class InVector>
  void
  Tvmult(OutVector &dst, const InVector &src) const;

  /**
   * 添加矩阵-向量的乘法。在<i>dst</i>上添加<i>M*src</i>，<i>M</i>是这个矩阵。
   * 注意，虽然这个函数可以对所有提供迭代器类的向量进行操作，但它只对
   * @ref Vector  类型的对象真正有效。
   * 对于所有迭代元素或随机成员访问昂贵的类来说，这个函数并不高效。特别是，如果你想与BlockVector对象相乘，你应该考虑同时使用BlockChunkSparseMatrix。
   * 来源和目的地不能是同一个向量。
   *
   */
  template <class OutVector, class InVector>
  void
  vmult_add(OutVector &dst, const InVector &src) const;

  /**
   * 添加矩阵-向量的乘法。将<i>M<sup>T</sup>*src</i>加到<i>dst</i>，<i>M</i>是这个矩阵。这个函数与vmult_add()的操作相同，但取的是转置的矩阵。
   * 注意，虽然这个函数可以对所有提供迭代器类的向量进行操作，但它只对类型为
   * @ref Vector  的对象真正有效。
   * 对于所有迭代元素或随机成员访问昂贵的类来说，这个函数并不高效。特别是，如果你想与BlockVector对象相乘，你应该考虑同时使用BlockChunkSparseMatrix。
   * 来源和目的地不能是同一个向量。
   *
   */
  template <class OutVector, class InVector>
  void
  Tvmult_add(OutVector &dst, const InVector &src) const;

  /**
   * 返回向量 $v$ 相对于该矩阵诱导的法线的平方，即
   * $\left(v,Mv\right)$
   * 。这很有用，例如在有限元背景下，一个函数的 $L_2$
   * 规范等于相对于代表有限元函数节点值的向量的质量矩阵的矩阵规范。
   * 显然，对于这个操作来说，矩阵需要是二次的，而且为了使结果真正成为一个规范，它还需要是实数对称的或复数隐式的。
   * 该矩阵和给定向量的基础模板类型应该都是实值或复值，但不是混合的，这样这个函数才有意义。
   *
   */
  template <typename somenumber>
  somenumber
  matrix_norm_square(const Vector<somenumber> &v) const;

  /**
   * 计算矩阵标量乘积  $\left(u,Mv\right)$  。
   *
   */
  template <typename somenumber>
  somenumber
  matrix_scalar_product(const Vector<somenumber> &u,
                        const Vector<somenumber> &v) const;
  /**
   * 计算方程<i>Mx=b</i>的残差，其中残差被定义为<i>r=b-Mx</i>。将残差写入<tt>dst</tt>。残差向量的<i>l<sub>2</sub></i>准则被返回。
   * 源<i>x</i>和目的<i>dst</i>不能是同一个向量。
   *
   */
  template <typename somenumber>
  somenumber
  residual(Vector<somenumber> &      dst,
           const Vector<somenumber> &x,
           const Vector<somenumber> &b) const;

  //@}
  /**
   * @name  矩阵规范
   *
   */
  //@{

  /**
   * 返回矩阵的l1准则，即 $|M|_1=max_{all columns j}\sum_{all rows i}
   * |M_ij|$  ，（最大列数之和）。
   * 这是自然的矩阵规范，与向量的l1规范兼容，即
   * $|Mv|_1\leq |M|_1 |v|_1$  。 (参见Haemmerlin-Hoffmann : Numerische
   * Mathematik)
   *
   */
  real_type
  l1_norm() const;

  /**
   * 返回矩阵的linfty-norm，即 $|M|_infty=max_{all rows i}\sum_{all
   * columns j} |M_ij|$  , (行的最大和)。
   * 这是一个自然的矩阵规范，与向量的linfty-norm兼容，即
   * $|Mv|_infty \leq |M|_infty |v|_infty$  。 (参见Haemmerlin-Hoffmann :
   * Numerische Mathematik)
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
   * 应用雅可比预处理方法，将<tt>src</tt>向量的每个元素乘以各自对角线元素的逆值，并将结果与松弛因子<tt>omega</tt>相乘。
   *
   */
  template <typename somenumber>
  void
  precondition_Jacobi(Vector<somenumber> &      dst,
                      const Vector<somenumber> &src,
                      const number              omega = 1.) const;

  /**
   * 对<tt>src</tt>应用SSOR预处理。
   *
   */
  template <typename somenumber>
  void
  precondition_SSOR(Vector<somenumber> &      dst,
                    const Vector<somenumber> &src,
                    const number              om = 1.) const;

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
   * 标准的SOR方法按照<tt>permutation</tt>规定的顺序应用，即先是行<tt>permutation[0]</tt>，然后是<tt>permutation[1]</tt>等等。出于效率的考虑，需要排列组合以及它的逆向排列。
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
   * 转置的SOR方法按照<tt>permutation</tt>规定的顺序应用，即首先是行<tt>permutation[m()-1]</tt>，然后是<tt>permutation[m()-2]</tt>等等。出于效率的考虑，需要用到permutation以及它的逆向。
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
   * 对<tt>v</tt>做一个SOR步骤。
   * 对右边的<tt>b</tt>进行直接的SOR步骤。
   *
   */
  template <typename somenumber>
  void
  SOR_step(Vector<somenumber> &      v,
           const Vector<somenumber> &b,
           const number              om = 1.) const;

  /**
   * 对<tt>v</tt>做一个邻接的SOR步骤。
   * 对<tt>b</tt>的右手边做一个直接的TSOR步骤。
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
   * 迭代器从矩阵的第一个条目开始。这是对常数矩阵的版本。
   * 请注意，由于ChunkSparseMatrix中的布局，对矩阵条目的迭代要比稀疏矩阵的迭代慢得多，因为迭代器是逐行旅行的，而数据是以若干行和列的块状形式存储的。
   *
   */
  const_iterator
  begin() const;

  /**
   * 最后的迭代器。这是用于常数矩阵的版本。
   * 请注意，由于ChunkSparseMatrix的布局，对矩阵条目的迭代要比稀疏矩阵慢得多，因为迭代器是逐行进行的，而数据是以若干行和列的块来存储的。
   *
   */
  const_iterator
  end() const;

  /**
   * 迭代器从矩阵的第一个条目开始。这是对非常数矩阵的版本。
   * 请注意，由于ChunkSparseMatrix的布局，对矩阵条目的迭代比稀疏矩阵的迭代要慢得多，因为迭代器是逐行进行的，而数据是以几行几列的块状存储的。
   *
   */
  iterator
  begin();

  /**
   * 最终迭代器。这是用于非恒定矩阵的版本。
   * 请注意，由于ChunkSparseMatrix的布局，对矩阵条目的迭代要比稀疏矩阵的迭代慢得多，因为迭代器是逐行进行的，而数据是以若干行和列的块来存储的。
   *
   */
  iterator
  end();

  /**
   * 迭代器从<tt>r</tt>行的第一个条目开始。这是对常数矩阵的版本。
   * 注意，如果给定的行是空的，即不包含任何非零条目，那么这个函数返回的迭代器就等于<tt>end(r)</tt>。还要注意的是，在这种情况下，迭代器可能无法被解除引用。
   * 请注意，由于ChunkSparseMatrix的布局，对矩阵条目的迭代要比稀疏矩阵慢得多，因为迭代器是逐行进行的，而数据是以几行几列的块状存储的。
   *
   */
  const_iterator
  begin(const unsigned int r) const;

  /**
   * 行<tt>r</tt>的最终迭代器。它指向过了 @p r,
   * 行末尾的第一个元素或过了整个稀疏模式的末尾。这是对常数矩阵的版本。
   * 请注意，结束迭代器不一定是可取消引用的。特别是如果它是一个矩阵最后一行的结束迭代器，情况更是如此。
   * 请注意，由于ChunkSparseMatrix的布局，对矩阵条目的迭代要比稀疏矩阵的迭代慢得多，因为迭代器是逐行进行的，而数据是以若干行和列的块来存储的。
   *
   */
  const_iterator
  end(const unsigned int r) const;

  /**
   * 迭代器从<tt>r</tt>行的第一个条目开始。这是用于非恒定矩阵的版本。
   * 注意，如果给定的行是空的，即不包含任何非零条目，那么这个函数返回的迭代器就等于<tt>end(r)</tt>。还要注意的是，在这种情况下，迭代器可能无法被解除引用。
   * 请注意，由于ChunkSparseMatrix的布局，对矩阵条目的迭代要比稀疏矩阵慢得多，因为迭代器是逐行进行的，而数据是以几行几列的块状存储的。
   *
   */
  iterator
  begin(const unsigned int r);

  /**
   * 行<tt>r</tt>的最终迭代器。它指向超过行 @p r,
   * 末尾的第一个元素，或者超过整个稀疏模式的末尾。这是对非恒定矩阵的版本。
   * 请注意，结束迭代器不一定是可取消引用的。特别是如果它是一个矩阵最后一行的结束迭代器，情况更是如此。
   * 请注意，由于ChunkSparseMatrix的布局，对矩阵条目的迭代要比稀疏矩阵的迭代慢得多，因为迭代器是逐行进行的，而数据是以若干行和列的块来存储的。
   *
   */
  iterator
  end(const unsigned int r);
  //@}
  /**
   * @name  输入/输出
   *
   */
  //@{

  /**
   * 打印矩阵到给定的流，使用格式<tt>(line,col)
   * value</tt>，即每行一个非零的矩阵条目。
   *
   */
  void
  print(std::ostream &out) const;

  /**
   * 以通常的格式打印矩阵，即作为矩阵而不是作为非零元素的列表。为了提高可读性，不在矩阵中的元素显示为空白，而明确设置为零的矩阵元素则显示为空白。
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
   * 将此对象的数据全部写到文件中。这是以二进制模式进行的，所以输出的数据既不能被人类阅读，也不能（可能）被其他使用不同操作系统或数字格式的计算机阅读。
   * 这个函数的目的是，如果你的内存不足，想在不同的程序之间进行交流，或者允许对象在程序的不同运行中持续存在，你可以把矩阵和稀疏模式换出来。
   *
   */
  void
  block_write(std::ostream &out) const;

  /**
   * 从文件中读取先前由block_write()写入的数据。
   * 这是用上述函数的逆运算来完成的，所以它的速度相当快，因为除了前面的几个数字，比特流是不被解释的。
   * 在这个操作中，对象被调整了大小，所有以前的内容都会丢失。然而，请注意，没有检查新数据和底层的ChunkSparsityPattern对象是否适合在一起。你有责任确保稀疏度模式和要读取的数据是匹配的。
   * 一个原始形式的错误检查会被执行，它将识别最直白的尝试，即把一些数据解释为一个矩阵，以比特方式存储到一个实际上不是以这种方式创建的文件，但不会更多。
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
  DeclException0(ExcDifferentChunkSparsityPatterns);
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
private:
  /**
   * 指向用于该矩阵的稀疏模式的指针。为了保证它在使用中不被删除，我们使用SmartPointer类来订阅它。
   *
   */
  SmartPointer<const ChunkSparsityPattern, ChunkSparseMatrix<number>> cols;

  /**
   * 所有非零条目的数值数组。一个条目在矩阵中的位置，也就是这个数组中给定值的行号和列号，只能用稀疏模式来推导。同样的道理也适用于更常见的通过坐标寻找一个条目的操作。
   *
   */
  std::unique_ptr<number[]> val;

  /**
   * 拨出的#val的大小。如果在过去的某个时候，通过使用reinit()函数，将具有较小尺寸的稀疏模式与该对象相关联，从而减少了矩阵的尺寸，则该尺寸可能大于实际使用的部分。
   *
   */
  size_type max_len;

  /**
   * 返回Val数组中条目 $(i,j)$ 的位置。
   *
   */
  size_type
  compute_location(const size_type i, const size_type j) const;

  // make all other sparse matrices friends
  template <typename somenumber>
  friend class ChunkSparseMatrix;

  // Also give access to internal details to the iterator/accessor classes.
  template <typename, bool>
  friend class ChunkSparseMatrixIterators::Iterator;
  template <typename, bool>
  friend class ChunkSparseMatrixIterators::Accessor;
};

 /*@}*/ 

#  ifndef DOXYGEN
 /*---------------------- Inline functions -----------------------------------*/ 



template <typename number>
inline typename ChunkSparseMatrix<number>::size_type
ChunkSparseMatrix<number>::m() const
{
  Assert(cols != nullptr, ExcNeedsSparsityPattern());
  return cols->rows;
}


template <typename number>
inline typename ChunkSparseMatrix<number>::size_type
ChunkSparseMatrix<number>::n() const
{
  Assert(cols != nullptr, ExcNeedsSparsityPattern());
  return cols->cols;
}



template <typename number>
inline const ChunkSparsityPattern &
ChunkSparseMatrix<number>::get_sparsity_pattern() const
{
  Assert(cols != nullptr, ExcNeedsSparsityPattern());
  return *cols;
}



template <typename number>
inline typename ChunkSparseMatrix<number>::size_type
ChunkSparseMatrix<number>::compute_location(const size_type i,
                                            const size_type j) const
{
  const size_type chunk_size = cols->get_chunk_size();
  const size_type chunk_index =
    cols->sparsity_pattern(i / chunk_size, j / chunk_size);

  if (chunk_index == ChunkSparsityPattern::invalid_entry)
    return ChunkSparsityPattern::invalid_entry;
  else
    {
      return (chunk_index * chunk_size * chunk_size +
              (i % chunk_size) * chunk_size + (j % chunk_size));
    }
}


template <typename number>
inline void
ChunkSparseMatrix<number>::set(const size_type i,
                               const size_type j,
                               const number    value)
{
  AssertIsFinite(value);

  Assert(cols != nullptr, ExcNeedsSparsityPattern());
  // it is allowed to set elements of the matrix that are not part of the
  // sparsity pattern, if the value to which we set it is zero
  const size_type index = compute_location(i, j);
  Assert((index != SparsityPattern::invalid_entry) || (value == 0.),
         ExcInvalidIndex(i, j));

  if (index != SparsityPattern::invalid_entry)
    val[index] = value;
}



template <typename number>
inline void
ChunkSparseMatrix<number>::add(const size_type i,
                               const size_type j,
                               const number    value)
{
  AssertIsFinite(value);

  Assert(cols != nullptr, ExcNeedsSparsityPattern());

  if (std::abs(value) != 0.)
    {
      const size_type index = compute_location(i, j);
      Assert((index != ChunkSparsityPattern::invalid_entry),
             ExcInvalidIndex(i, j));

      val[index] += value;
    }
}



template <typename number>
template <typename number2>
inline void
ChunkSparseMatrix<number>::add(const size_type  row,
                               const size_type  n_cols,
                               const size_type *col_indices,
                               const number2 *  values,
                               const bool  /*elide_zero_values*/ ,
                               const bool  /*col_indices_are_sorted*/ )
{
  // TODO: could be done more efficiently...
  for (size_type col = 0; col < n_cols; ++col)
    add(row, col_indices[col], static_cast<number>(values[col]));
}



template <typename number>
inline ChunkSparseMatrix<number> &
ChunkSparseMatrix<number>::operator*=(const number factor)
{
  Assert(cols != nullptr, ExcNeedsSparsityPattern());
  Assert(val != nullptr, ExcNotInitialized());

  const size_type chunk_size = cols->get_chunk_size();

  // multiply all elements of the matrix with the given factor. this includes
  // the padding elements in chunks that overlap the boundaries of the actual
  // matrix -- but since multiplication with a number does not violate the
  // invariant of keeping these elements at zero nothing can happen
  number *            val_ptr = val.get();
  const number *const end_ptr =
    val.get() +
    cols->sparsity_pattern.n_nonzero_elements() * chunk_size * chunk_size;
  while (val_ptr != end_ptr)
    *val_ptr++ *= factor;

  return *this;
}



template <typename number>
inline ChunkSparseMatrix<number> &
ChunkSparseMatrix<number>::operator/=(const number factor)
{
  Assert(cols != nullptr, ExcNeedsSparsityPattern());
  Assert(val != nullptr, ExcNotInitialized());
  Assert(std::abs(factor) != 0, ExcDivideByZero());

  const number factor_inv = 1. / factor;

  const size_type chunk_size = cols->get_chunk_size();

  // multiply all elements of the matrix with the given factor. this includes
  // the padding elements in chunks that overlap the boundaries of the actual
  // matrix -- but since multiplication with a number does not violate the
  // invariant of keeping these elements at zero nothing can happen
  number *            val_ptr = val.get();
  const number *const end_ptr =
    val.get() +
    cols->sparsity_pattern.n_nonzero_elements() * chunk_size * chunk_size;

  while (val_ptr != end_ptr)
    *val_ptr++ *= factor_inv;

  return *this;
}



template <typename number>
inline number
ChunkSparseMatrix<number>::operator()(const size_type i,
                                      const size_type j) const
{
  Assert(cols != nullptr, ExcNeedsSparsityPattern());
  AssertThrow(compute_location(i, j) != SparsityPattern::invalid_entry,
              ExcInvalidIndex(i, j));
  return val[compute_location(i, j)];
}



template <typename number>
inline number
ChunkSparseMatrix<number>::el(const size_type i, const size_type j) const
{
  Assert(cols != nullptr, ExcNeedsSparsityPattern());
  const size_type index = compute_location(i, j);

  if (index != ChunkSparsityPattern::invalid_entry)
    return val[index];
  else
    return 0;
}



template <typename number>
inline number
ChunkSparseMatrix<number>::diag_element(const size_type i) const
{
  Assert(cols != nullptr, ExcNeedsSparsityPattern());
  Assert(m() == n(), ExcNotQuadratic());
  AssertIndexRange(i, m());

  // Use that the first element in each row of a quadratic matrix is the main
  // diagonal of the chunk sparsity pattern
  const size_type chunk_size = cols->get_chunk_size();
  return val[cols->sparsity_pattern.rowstart[i / chunk_size] * chunk_size *
               chunk_size +
             (i % chunk_size) * chunk_size + (i % chunk_size)];
}



template <typename number>
template <typename ForwardIterator>
inline void
ChunkSparseMatrix<number>::copy_from(const ForwardIterator begin,
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
    }
}



//---------------------------------------------------------------------------


namespace ChunkSparseMatrixIterators
{
  template <typename number>
  inline Accessor<number, true>::Accessor(const MatrixType * matrix,
                                          const unsigned int row)
    : ChunkSparsityPatternIterators::Accessor(&matrix->get_sparsity_pattern(),
                                              row)
    , matrix(matrix)
  {}



  template <typename number>
  inline Accessor<number, true>::Accessor(const MatrixType *matrix)
    : ChunkSparsityPatternIterators::Accessor(&matrix->get_sparsity_pattern())
    , matrix(matrix)
  {}



  template <typename number>
  inline Accessor<number, true>::Accessor(
    const ChunkSparseMatrixIterators::Accessor<number, false> &a)
    : ChunkSparsityPatternIterators::Accessor(a)
    , matrix(&a.get_matrix())
  {}



  template <typename number>
  inline number
  Accessor<number, true>::value() const
  {
    const unsigned int chunk_size =
      matrix->get_sparsity_pattern().get_chunk_size();
    return matrix->val[reduced_index() * chunk_size * chunk_size +
                       chunk_row * chunk_size + chunk_col];
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
    const unsigned int chunk_size =
      accessor->matrix->get_sparsity_pattern().get_chunk_size();
    return accessor->matrix
      ->val[accessor->reduced_index() * chunk_size * chunk_size +
            accessor->chunk_row * chunk_size + accessor->chunk_col];
  }



  template <typename number>
  inline const typename Accessor<number, false>::Reference &
  Accessor<number, false>::Reference::operator=(const number n) const
  {
    const unsigned int chunk_size =
      accessor->matrix->get_sparsity_pattern().get_chunk_size();
    accessor->matrix
      ->val[accessor->reduced_index() * chunk_size * chunk_size +
            accessor->chunk_row * chunk_size + accessor->chunk_col] = n;
    return *this;
  }



  template <typename number>
  inline const typename Accessor<number, false>::Reference &
  Accessor<number, false>::Reference::operator+=(const number n) const
  {
    const unsigned int chunk_size =
      accessor->matrix->get_sparsity_pattern().get_chunk_size();
    accessor->matrix
      ->val[accessor->reduced_index() * chunk_size * chunk_size +
            accessor->chunk_row * chunk_size + accessor->chunk_col] += n;
    return *this;
  }



  template <typename number>
  inline const typename Accessor<number, false>::Reference &
  Accessor<number, false>::Reference::operator-=(const number n) const
  {
    const unsigned int chunk_size =
      accessor->matrix->get_sparsity_pattern().get_chunk_size();
    accessor->matrix
      ->val[accessor->reduced_index() * chunk_size * chunk_size +
            accessor->chunk_row * chunk_size + accessor->chunk_col] -= n;
    return *this;
  }



  template <typename number>
  inline const typename Accessor<number, false>::Reference &
  Accessor<number, false>::Reference::operator*=(const number n) const
  {
    const unsigned int chunk_size =
      accessor->matrix->get_sparsity_pattern().get_chunk_size();
    accessor->matrix
      ->val[accessor->reduced_index() * chunk_size * chunk_size +
            accessor->chunk_row * chunk_size + accessor->chunk_col] *= n;
    return *this;
  }



  template <typename number>
  inline const typename Accessor<number, false>::Reference &
  Accessor<number, false>::Reference::operator/=(const number n) const
  {
    const unsigned int chunk_size =
      accessor->matrix->get_sparsity_pattern().get_chunk_size();
    accessor->matrix
      ->val[accessor->reduced_index() * chunk_size * chunk_size +
            accessor->chunk_row * chunk_size + accessor->chunk_col] /= n;
    return *this;
  }



  template <typename number>
  inline Accessor<number, false>::Accessor(MatrixType *       matrix,
                                           const unsigned int row)
    : ChunkSparsityPatternIterators::Accessor(&matrix->get_sparsity_pattern(),
                                              row)
    , matrix(matrix)
  {}



  template <typename number>
  inline Accessor<number, false>::Accessor(MatrixType *matrix)
    : ChunkSparsityPatternIterators::Accessor(&matrix->get_sparsity_pattern())
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
  inline Iterator<number, Constness>::Iterator(MatrixType *       matrix,
                                               const unsigned int row)
    : accessor(matrix, row)
  {}



  template <typename number, bool Constness>
  inline Iterator<number, Constness>::Iterator(MatrixType *matrix)
    : accessor(matrix)
  {}



  template <typename number, bool Constness>
  inline Iterator<number, Constness>::Iterator(
    const ChunkSparseMatrixIterators::Iterator<number, false> &i)
    : accessor(*i)
  {}



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

    // TODO: can be optimized
    int difference = 0;
    if (*this < other)
      {
        Iterator copy = *this;
        while (copy != other)
          {
            ++copy;
            --difference;
          }
      }
    else
      {
        Iterator copy = other;
        while (copy != *this)
          {
            ++copy;
            ++difference;
          }
      }
    return difference;
  }



  template <typename number, bool Constness>
  inline Iterator<number, Constness>
  Iterator<number, Constness>::operator+(const unsigned int n) const
  {
    Iterator x = *this;
    for (unsigned int i = 0; i < n; ++i)
      ++x;

    return x;
  }

} // namespace ChunkSparseMatrixIterators



template <typename number>
inline typename ChunkSparseMatrix<number>::const_iterator
ChunkSparseMatrix<number>::begin() const
{
  return const_iterator(this, 0);
}


template <typename number>
inline typename ChunkSparseMatrix<number>::const_iterator
ChunkSparseMatrix<number>::end() const
{
  return const_iterator(this);
}


template <typename number>
inline typename ChunkSparseMatrix<number>::iterator
ChunkSparseMatrix<number>::begin()
{
  return iterator(this, 0);
}


template <typename number>
inline typename ChunkSparseMatrix<number>::iterator
ChunkSparseMatrix<number>::end()
{
  return iterator(this);
}


template <typename number>
inline typename ChunkSparseMatrix<number>::const_iterator
ChunkSparseMatrix<number>::begin(const unsigned int r) const
{
  AssertIndexRange(r, m());
  return const_iterator(this, r);
}



template <typename number>
inline typename ChunkSparseMatrix<number>::const_iterator
ChunkSparseMatrix<number>::end(const unsigned int r) const
{
  AssertIndexRange(r, m());
  return const_iterator(this, r + 1);
}



template <typename number>
inline typename ChunkSparseMatrix<number>::iterator
ChunkSparseMatrix<number>::begin(const unsigned int r)
{
  AssertIndexRange(r, m());
  return iterator(this, r);
}



template <typename number>
inline typename ChunkSparseMatrix<number>::iterator
ChunkSparseMatrix<number>::end(const unsigned int r)
{
  AssertIndexRange(r, m());
  return iterator(this, r + 1);
}



#  endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif
 /*--------------------------- chunk_sparse_matrix.h -------------------------*/ 


