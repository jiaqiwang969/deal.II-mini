//include/deal.II-translator/lac/sparse_matrix_ez_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2002 - 2020 by the deal.II authors
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

#ifndef dealii_sparse_matrix_ez_h
#  define dealii_sparse_matrix_ez_h


#  include <deal.II/base/config.h>

#  include <deal.II/base/smartpointer.h>
#  include <deal.II/base/subscriptor.h>

#  include <deal.II/lac/exceptions.h>

#  include <vector>

DEAL_II_NAMESPACE_OPEN

// Forward declarations
#  ifndef DOXYGEN
template <typename number>
class Vector;
template <typename number>
class FullMatrix;
#  endif

/**
 * @addtogroup  Matrix1  
     * @{ 
 *
 *
 */

/**
 * 无稀疏模式的稀疏矩阵。
 * 这个矩阵没有使用预先组装好的稀疏模式，而是在飞行中建立了模式。填充矩阵可能比
 * @p SparseMatrix,
 * 消耗更多的时间，因为当新的矩阵元素被插入到矩阵中间的某个地方时，可能会涉及到大量的内存移动，而目前没有未使用的内存位置可用于插入新条目的行。为了帮助优化，可以向构造函数提供一个预期的行长度，以及一个行的增量大小。
 * 该类使用一个存储结构，与通常的稀疏矩阵格式类似，只存储非零元素。这些元素被存储在整个矩阵的一个数据数组中，并按行排序，在每行中按列号排序。一个单独的数组描述了每一行在长数据数组中的起始位置以及它的长度。
 * 由于这种结构，行与行之间可能出现空隙。每当必须创建一个新条目时，就会尝试使用其行中的间隙。如果没有空隙，该行必须被扩展，所有后续的行都必须向后移位。这是一个非常昂贵的操作，解释了这种数据结构的低效率，以及为什么像SparseMatrix类那样预先分配一个稀疏模式是有用的。
 * 这就是提供给构造函数或 reinit()函数的优化参数的作用。
 * @p default_row_length
 * 是初始化时为每行分配的条目数（行的实际长度仍为0）。这意味着，
 * @p default_row_length
 * 个条目可以被添加到这一行，而不会转移其他行。如果添加的条目较少，额外的内存当然就会被浪费掉。
 * 如果一个行的空间不够，那么它将被扩大 @p default_increment
 * 个条目。这样一来，后续的行就不会经常被单项移位了。
 * 最后， @p default_reserve
 * 在数据数组的末端分配了额外的空间。这个空间在任何必须扩大的行中都会被使用。这很重要，因为否则不仅下面的行必须被移动，而且在为整个数据阵列分配了足够多的空间后，实际上<i>all</i>行也必须被移动。
 * 建议的设置。  @p default_row_length
 * 应该是一个典型行的长度，例如网格中常规部分的网板尺寸。然后，
 * @p default_increment
 * 可以是有一个挂起的节点给该行增加的预期条目量。这样一来，就应该在内存消耗和速度之间取得一个很好的折中。
 * @p default_reserve 应该是对悬挂节点数量的估计乘以 @p
 * default_increment。 让 @p default_increment
 * 为零会导致每当有行溢出时出现异常。
 * 如果行被期望或多或少地从头到尾填满，使用 @p
 * default_row_length 为零可能不是一个坏主意。
 *
 *
 * @note
 * 这个类的名字用美国人的方式发音是有意义的，其中 "EZ
 * "的发音与 "easy "的发音相同。
 *
 *
 */
template <typename number>
class SparseMatrixEZ : public Subscriptor
{
public:
  /**
   * 声明容器尺寸的类型。
   *
   */
  using size_type = types::global_dof_index;

  /**
   * 用于存储条目的列号及其值的类。
   *
   */
  struct Entry
  {
    /**
     * 标准构造函数。设置 @p column 到 @p invalid. 。
     *
     */
    Entry();

    /**
     * 构造函数。填充列和值。
     *
     */
    Entry(const size_type column, const number &value);

    /**
     * 列的编号。
     *
     */
    size_type column;

    /**
     * 那里的值。
     *
     */
    number value;

    /**
     * 不存在的列号。
     *
     */
    static const size_type invalid = numbers::invalid_size_type;
  };

  /**
   * 用于存储矩阵行的信息的结构。每行有一个对象，将被存储在矩阵中。
   *
   */
  struct RowInfo
  {
    /**
     * 构造函数。
     *
     */
    RowInfo(const size_type start = Entry::invalid);

    /**
     * 数据域中行的第一个条目的索引。
     *
     */
    size_type start;
    /**
     * 该行的条目数。
     *
     */
    unsigned short length;
    /**
     * 对角线元素相对于起始索引的位置。
     *
     */
    unsigned short diagonal;
    /**
     * 不存在的对角线的值。
     *
     */
    static const unsigned short invalid_diagonal =
      static_cast<unsigned short>(-1);
  };

public:
  /**
   * 符合标准的迭代器。
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
       * 构造器。因为我们只使用访问器进行读取访问，一个常量矩阵指针就足够了。
       *
       */
      Accessor(const SparseMatrixEZ<number> *matrix,
               const size_type               row,
               const unsigned short          index);

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
      unsigned short
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
      number
      value() const;

    protected:
      /**
       * 访问的矩阵。
       *
       */
      const SparseMatrixEZ<number> *matrix;

      /**
       * 当前行数。
       *
       */
      size_type a_row;

      /**
       * 当前行的索引。
       *
       */
      unsigned short a_index;

      // Make enclosing class a friend.
      friend class const_iterator;
    };

  public:
    /**
     * 构造函数。
     *
     */
    const_iterator(const SparseMatrixEZ<number> *matrix,
                   const size_type               row,
                   const unsigned short          index);

    /**
     * 前缀增量。这总是返回一个有效的条目或<tt>end()</tt>。
     *
     */
    const_iterator &
    operator++();

    /**
     * 去引用操作符。
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

  private:
    /**
     * 存储一个访问器类的对象。
     *
     */
    Accessor accessor;
  };

  /**
   * 矩阵条目的类型。这个别名类似于标准库容器中的<tt>value_type</tt>。
   *
   */
  using value_type = number;

  /**
   * @name  构造函数和初始化
   *
   */
  //@{
  /**
   * 构造函数。初始化一个尺寸为0乘以0的空矩阵。
   *
   */
  SparseMatrixEZ();

  /**
   * 假的复制构造函数。这是在容器中使用的。它只能为空对象调用。
   * 如果你真的想复制一个完整的矩阵，你可以使用 @p
   * copy_from函数来实现。
   *
   */
  SparseMatrixEZ(const SparseMatrixEZ &);

  /**
   * 构造函数。生成一个给定大小的矩阵，准备被填充。
   * 可选参数 @p default_row_length 和 @p default_increment
   * 允许预先分配内存。适当地提供这些参数对于有效地组装矩阵是至关重要的。
   *
   */
  explicit SparseMatrixEZ(const size_type    n_rows,
                          const size_type    n_columns,
                          const size_type    default_row_length = 0,
                          const unsigned int default_increment  = 1);

  /**
   * 销毁器。释放所有的内存。
   *
   */
  ~SparseMatrixEZ() override = default;

  /**
   * 只复制空对象的伪操作符。
   *
   */
  SparseMatrixEZ<number> &
  operator=(const SparseMatrixEZ<number> &);

  /**
   * 这个操作符将一个标量分配给一个矩阵。因为这通常没有什么意义（我们应该把所有的矩阵条目都设置为这个值吗？仅仅是稀疏模式的非零条目？），这个操作只允许在实际要分配的值为零的情况下进行。这个操作符的存在只是为了允许明显的符号<tt>matrix=0</tt>，它将矩阵的所有元素设置为零，但保留之前使用的稀疏模式。
   *
   */
  SparseMatrixEZ<number> &
  operator=(const double d);

  /**
   * 将稀疏矩阵重新初始化为所提供的尺寸。矩阵在此时将没有任何条目。可选参数
   * @p default_row_length， @p default_increment 和 @p reserve
   * 允许预先分配内存。正确地提供这些参数对于有效地组装矩阵是至关重要的。
   *
   */
  void
  reinit(const size_type n_rows,
         const size_type n_columns,
         size_type       default_row_length = 0,
         unsigned int    default_increment  = 1,
         size_type       reserve            = 0);

  /**
   * 释放所有内存并返回到与调用默认构造函数后相同的状态。它也会忘记其稀疏模式。
   *
   */
  void
  clear();
  //@}
  /**
   * @name  矩阵的信息
   *
   */
  //@{
  /**
   * 返回该对象是否为空。如果两个维度都是零，它就是空的。
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
   * 返回该矩阵的非零元素的数量。
   *
   */
  size_type
  n_nonzero_elements() const;

  /**
   * 确定此对象的内存消耗（以字节为单位）的估计值。
   *
   */
  std::size_t
  memory_consumption() const;

  /**
   * 打印统计数据。如果 @p full 是 @p true,
   * ，则打印所有现有行长和分配行长的直方图。否则，只显示已分配和已使用条目的关系。
   *
   */
  template <class StreamType>
  void
  print_statistics(StreamType &s, bool full = false);

  /**
   * 计算条目的数量。
   * 在前三个参数中，该函数返回该矩阵使用、分配和保留的条目数。
   * 如果最后一个参数为真，每行的条目数也会被打印出来。
   *
   */
  void
  compute_statistics(size_type &             used,
                     size_type &             allocated,
                     size_type &             reserved,
                     std::vector<size_type> &used_by_line,
                     const bool              compute_by_line) const;
  //@}
  /**
   * @name  修改条目
   *
   */
  //@{
  /**
   * 将元素<tt>(i,j)</tt>设置为  @p value.
   * 如果<tt>value</tt>不是一个有限的数字，就会产生异常。
   * 可选参数<tt>elide_zero_values</tt>可以用来指定是无论如何都要添加零值，还是要过滤掉这些零值，只添加非零数据。默认值是<tt>true</tt>，也就是说，零值不会被添加到矩阵中。
   * 如果这个函数设置了一个尚不存在的元素的值，那么它将为其分配一个条目。(除非如上所述，`elide_zero_values`是`true`)。
   * @note
   * 如果你想保持矩阵的对称稀疏模式，你可能需要插入零元素。
   *
   */
  void
  set(const size_type i,
      const size_type j,
      const number    value,
      const bool      elide_zero_values = true);

  /**
   * 将 @p value 添加到元素<tt>(i,j)</tt>。
   * 如果这个函数添加到一个尚不存在的元素的值中，那么它将为其分配一个条目。
   * 该函数会自动过滤掉零，也就是说，当向一个目前不存在条目的矩阵元素添加零时，它不会创建新条目。
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
   * 将作为参数给出的矩阵复制到当前对象中。
   * 复制矩阵是一个昂贵的操作，我们不希望通过编译器生成的代码意外发生
   * <code>operator=</code>
   * 。（例如，如果不小心声明了一个当前类型为<i>by
   * value</i>而非<i>by
   * reference</i>的函数参数，就会发生这种情况）。复制矩阵的功能是在这个成员函数中实现的。因此，该类型对象的所有复制操作都需要一个明确的函数调用。
   * 源矩阵可以是一个任意类型的矩阵，只要其数据类型可以转换为该矩阵的数据类型。
   * 可选参数<tt>elide_zero_values</tt>可以用来指定是无论如何都要添加零值，还是要过滤掉这些零值，只添加非零数据。默认值是<tt>true</tt>，也就是说，零值不会被添加到矩阵中。
   * 该函数返回一个 @p this. 的引用。
   *
   */
  template <typename MatrixType>
  SparseMatrixEZ<number> &
  copy_from(const MatrixType &source, const bool elide_zero_values = true);

  /**
   * 将 @p matrix 按 @p factor 的比例添加到该矩阵中。
   * 源矩阵可以是一个任意类型的矩阵，只要其数据类型可以转换为该矩阵的数据类型，并且具有标准的
   * @p const_iterator. 。
   *
   */
  template <typename MatrixType>
  void
  add(const number factor, const MatrixType &matrix);
  //@}
  /**
   * @name 条目访问
   *
   */
  //@{
  /**
   * 返回条目(i,j)的值。
   * 这可能是一个昂贵的操作，你应该始终注意在哪里调用这个函数。
   * 为了避免滥用，如果所需元素在矩阵中不存在，该函数会抛出一个异常。
   * 如果你想要一个返回零的函数（对于不在矩阵的稀疏模式中的条目），请使用
   * @p el 函数。
   *
   */
  number
  operator()(const size_type i, const size_type j) const;

  /**
   * 返回条目(i,j)的值。对所有不存在的条目返回零。
   *
   */
  number
  el(const size_type i, const size_type j) const;
  //@}
  /**
   * @name  乘法运算
   *
   */
  //@{
  /**
   * 矩阵-向量乘法：让 $dst = M*src$ 与 $M$ 为该矩阵。
   *
   */
  template <typename somenumber>
  void
  vmult(Vector<somenumber> &dst, const Vector<somenumber> &src) const;

  /**
   * 矩阵-向量乘法：让 $dst = M^T*src$ 与 $M$
   * 为这个矩阵。这个函数与 @p vmult
   * 的作用相同，但需要转置的矩阵。
   *
   */
  template <typename somenumber>
  void
  Tvmult(Vector<somenumber> &dst, const Vector<somenumber> &src) const;

  /**
   * 加法 矩阵-向量乘法。在 $dst$ 上添加 $M*src$ ， $M$
   * 为该矩阵。
   *
   */
  template <typename somenumber>
  void
  vmult_add(Vector<somenumber> &dst, const Vector<somenumber> &src) const;

  /**
   * 添加矩阵-向量乘法。将 $M^T*src$ 加到 $dst$ ， $M$
   * 是这个矩阵。这个函数与 @p vmult_add
   * 的作用相同，但需要转置的矩阵。
   *
   */
  template <typename somenumber>
  void
  Tvmult_add(Vector<somenumber> &dst, const Vector<somenumber> &src) const;
  //@}
  /**
   * @name  矩阵的准则
   *
   */
  //@{
  /**
   * 矩阵的Frobenius-norm。
   *
   */
  number
  l2_norm() const;
  //@}
  /**
   * @name  预处理方法
   *
   */
  //@{
  /**
   * 应用雅可比预处理方法，将 @p
   * src向量的每个元素乘以各自对角线元素的逆值，并将结果与阻尼系数
   * @p omega. 相乘。
   *
   */
  template <typename somenumber>
  void
  precondition_Jacobi(Vector<somenumber> &      dst,
                      const Vector<somenumber> &src,
                      const number              omega = 1.) const;

  /**
   * 对 @p src. 应用SSOR预处理。
   *
   */
  template <typename somenumber>
  void
  precondition_SSOR(Vector<somenumber> &            dst,
                    const Vector<somenumber> &      src,
                    const number                    om = 1.,
                    const std::vector<std::size_t> &pos_right_of_diagonal =
                      std::vector<std::size_t>()) const;

  /**
   * 将SOR预处理矩阵应用于 @p src. ，该方法的结果是 $dst =
   * (om D
   *
   * - L)^{-1} src$  。
   *
   */
  template <typename somenumber>
  void
  precondition_SOR(Vector<somenumber> &      dst,
                   const Vector<somenumber> &src,
                   const number              om = 1.) const;

  /**
   * 对 @p src. 应用转置的SOR预处理矩阵，该方法的结果是
   * $dst = (om D
   *
   * - U)^{-1} src$  。
   *
   */
  template <typename somenumber>
  void
  precondition_TSOR(Vector<somenumber> &      dst,
                    const Vector<somenumber> &src,
                    const number              om = 1.) const;

  /**
   * 将由 @p B, 共轭的矩阵 @p A 即 $B A B^T$
   * 添加到该对象中。如果参数 @p transpose 为真，计算 $B^T A
   * B$  。    这个函数要求 @p B 有一个 @p const_iterator
   * 遍历所有矩阵条目，并且 @p A
   * 有一个函数<tt>el(i,j)</tt>用于访问特定条目。
   *
   */
  template <typename MatrixTypeA, typename MatrixTypeB>
  void
  conjugate_add(const MatrixTypeA &A,
                const MatrixTypeB &B,
                const bool         transpose = false);
  //@}
  /**
   * @name  迭代器
   *
   */
  //@{
  /**
   * 迭代器从第一个现有条目开始。
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
   * 迭代器从第 @p r.
   * 行的第一个条目开始，如果这一行是空的，结果是<tt>end(r)</tt>，它并不指向第
   * @p r. 行。
   *
   */
  const_iterator
  begin(const size_type r) const;

  /**
   * 行 @p r. 的最终迭代器 结果可能与<tt>end()</tt>不同!
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
   * 打印矩阵到给定的流，使用格式<tt>(line,col)
   * value</tt>，即每行有一个非零的矩阵条目。
   *
   */
  void
  print(std::ostream &out) const;

  /**
   * 以通常的格式打印矩阵，即作为矩阵而不是作为非零元素的列表。为了提高可读性，不在矩阵中的元素显示为空白，而明确设置为零的矩阵元素则显示为空白。
   * 参数允许对输出格式进行灵活设置。  @p 精度和 @p
   * scientific 用于确定数字格式，其中 @p scientific = @p false
   * 表示定点符号。  @p width
   * 的零条目使函数计算出一个宽度，但如果输出是粗略的，它可能被改变为一个正值。
   * 此外，还可以指定一个空值的字符。
   * 最后，整个矩阵可以与一个共同的分母相乘，产生更可读的输出，甚至是整数。
   * 如果应用于一个大的矩阵，这个函数可能会产生 @em
   * 大量的输出!
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
   * 以二进制模式将此对象的数据写到文件中。
   * 注意，这种二进制格式与平台有关。
   *
   */
  void
  block_write(std::ostream &out) const;

  /**
   * 读取之前由 @p block_write.
   * 写入的数据，该对象在此操作中被调整大小，所有之前的内容都会丢失。
   * 执行一种原始形式的错误检查，它将识别最直白的尝试，将一些数据解释为按位数存储到文件中的向量，但不会有更多。
   *
   */
  void
  block_read(std::istream &in);
  //@}

  /**
   * @addtogroup  Exceptions  @{
   *
   */

  /**
   * 缺少对角线条目的异常情况。
   *
   */
  DeclException0(ExcNoDiagonal);

  /**
   * 例外情况
   *
   */
  DeclException2(ExcInvalidEntry,
                 int,
                 int,
                 << "The entry with index (" << arg1 << ',' << arg2
                 << ") does not exist.");

  DeclException2(ExcEntryAllocationFailure,
                 int,
                 int,
                 << "An entry with index (" << arg1 << ',' << arg2
                 << ") cannot be allocated.");
  //@}
private:
  /**
   * 找到一个条目并返回一个常数指针。如果该条目不存在，则返回一个零指针。
   *
   */
  const Entry *
  locate(const size_type row, const size_type col) const;

  /**
   * 找到一个条目并返回一个可写指针。如果该条目不存在，则返回一个零指针。
   *
   */
  Entry *
  locate(const size_type row, const size_type col);

  /**
   * 查找一个条目或生成它。
   *
   */
  Entry *
  allocate(const size_type row, const size_type col);

  /**
   * @p vmult
   * 的版本，只对<tt>[begin_row,end_row)</tt>定义的区域进行操作。在启用多线程的情况下，这个函数被
   * @p vmult 调用。
   *
   */
  template <typename somenumber>
  void
  threaded_vmult(Vector<somenumber> &      dst,
                 const Vector<somenumber> &src,
                 const size_type           begin_row,
                 const size_type           end_row) const;

  /**
   * @p matrix_norm_square
   * 的版本，只对<tt>[begin_row,end_row)</tt>定义的区域进行操作。在启用多线程的情况下，这个函数被
   * @p matrix_norm_square 调用。
   *
   */
  template <typename somenumber>
  void
  threaded_matrix_norm_square(const Vector<somenumber> &v,
                              const size_type           begin_row,
                              const size_type           end_row,
                              somenumber *              partial_sum) const;

  /**
   * @p matrix_scalar_product
   * 的版本，只对<tt>[begin_row,end_row)</tt>定义的区域进行操作。在启用多线程的情况下，这个函数被
   * @p matrix_scalar_product 调用。
   *
   */
  template <typename somenumber>
  void
  threaded_matrix_scalar_product(const Vector<somenumber> &u,
                                 const Vector<somenumber> &v,
                                 const size_type           begin_row,
                                 const size_type           end_row,
                                 somenumber *              partial_sum) const;

  /**
   * 列的数量。这仅用于检查向量尺寸。
   *
   */
  size_type n_columns;

  /**
   * 每一行的信息结构。
   *
   */
  std::vector<RowInfo> row_info;

  /**
   * 数据存储。
   *
   */
  std::vector<Entry> data;

  /**
   * 当某行增长时进行递增。
   *
   */
  unsigned int increment;

  /**
   * 记住用户提供的默认行长度。
   *
   */
  unsigned int saved_default_row_length;
};

/**
 * @}
 *
 *
 */
 /*---------------------- Inline functions -----------------------------------*/ 

template <typename number>
inline SparseMatrixEZ<number>::Entry::Entry(const size_type column,
                                            const number &  value)
  : column(column)
  , value(value)
{}



template <typename number>
inline SparseMatrixEZ<number>::Entry::Entry()
  : column(invalid)
  , value(0)
{}


template <typename number>
inline SparseMatrixEZ<number>::RowInfo::RowInfo(const size_type start)
  : start(start)
  , length(0)
  , diagonal(invalid_diagonal)
{}


//---------------------------------------------------------------------------
template <typename number>
inline SparseMatrixEZ<number>::const_iterator::Accessor::Accessor(
  const SparseMatrixEZ<number> *matrix,
  const size_type               r,
  const unsigned short          i)
  : matrix(matrix)
  , a_row(r)
  , a_index(i)
{}


template <typename number>
inline typename SparseMatrixEZ<number>::size_type
SparseMatrixEZ<number>::const_iterator::Accessor::row() const
{
  return a_row;
}


template <typename number>
inline typename SparseMatrixEZ<number>::size_type
SparseMatrixEZ<number>::const_iterator::Accessor::column() const
{
  return matrix->data[matrix->row_info[a_row].start + a_index].column;
}


template <typename number>
inline unsigned short
SparseMatrixEZ<number>::const_iterator::Accessor::index() const
{
  return a_index;
}



template <typename number>
inline number
SparseMatrixEZ<number>::const_iterator::Accessor::value() const
{
  return matrix->data[matrix->row_info[a_row].start + a_index].value;
}


template <typename number>
inline SparseMatrixEZ<number>::const_iterator::const_iterator(
  const SparseMatrixEZ<number> *matrix,
  const size_type               r,
  const unsigned short          i)
  : accessor(matrix, r, i)
{
  // Finish if this is the end()
  if (r == accessor.matrix->m() && i == 0)
    return;

  // Make sure we never construct an
  // iterator pointing to a
  // non-existing entry

  // If the index points beyond the
  // end of the row, try the next
  // row.
  if (accessor.a_index >= accessor.matrix->row_info[accessor.a_row].length)
    {
      do
        {
          ++accessor.a_row;
        }
      // Beware! If the next row is
      // empty, iterate until a
      // non-empty row is found or we
      // hit the end of the matrix.
      while (accessor.a_row < accessor.matrix->m() &&
             accessor.matrix->row_info[accessor.a_row].length == 0);
    }
}


template <typename number>
inline typename SparseMatrixEZ<number>::const_iterator &
SparseMatrixEZ<number>::const_iterator::operator++()
{
  Assert(accessor.a_row < accessor.matrix->m(), ExcIteratorPastEnd());

  // Increment column index
  ++(accessor.a_index);
  // If index exceeds number of
  // entries in this row, proceed
  // with next row.
  if (accessor.a_index >= accessor.matrix->row_info[accessor.a_row].length)
    {
      accessor.a_index = 0;
      // Do this loop to avoid
      // elements in empty rows
      do
        {
          ++accessor.a_row;
        }
      while (accessor.a_row < accessor.matrix->m() &&
             accessor.matrix->row_info[accessor.a_row].length == 0);
    }
  return *this;
}


template <typename number>
inline const typename SparseMatrixEZ<number>::const_iterator::Accessor &
  SparseMatrixEZ<number>::const_iterator::operator*() const
{
  return accessor;
}


template <typename number>
inline const typename SparseMatrixEZ<number>::const_iterator::Accessor *
  SparseMatrixEZ<number>::const_iterator::operator->() const
{
  return &accessor;
}


template <typename number>
inline bool
SparseMatrixEZ<number>::const_iterator::
operator==(const const_iterator &other) const
{
  return (accessor.row() == other.accessor.row() &&
          accessor.index() == other.accessor.index());
}


template <typename number>
inline bool
SparseMatrixEZ<number>::const_iterator::
operator!=(const const_iterator &other) const
{
  return !(*this == other);
}


template <typename number>
inline bool
SparseMatrixEZ<number>::const_iterator::
operator<(const const_iterator &other) const
{
  return (accessor.row() < other.accessor.row() ||
          (accessor.row() == other.accessor.row() &&
           accessor.index() < other.accessor.index()));
}


//---------------------------------------------------------------------------
template <typename number>
inline typename SparseMatrixEZ<number>::size_type
SparseMatrixEZ<number>::m() const
{
  return row_info.size();
}


template <typename number>
inline typename SparseMatrixEZ<number>::size_type
SparseMatrixEZ<number>::n() const
{
  return n_columns;
}


template <typename number>
inline typename SparseMatrixEZ<number>::Entry *
SparseMatrixEZ<number>::locate(const size_type row, const size_type col)
{
  AssertIndexRange(row, m());
  AssertIndexRange(col, n());

  const RowInfo & r   = row_info[row];
  const size_type end = r.start + r.length;
  for (size_type i = r.start; i < end; ++i)
    {
      Entry *const entry = &data[i];
      if (entry->column == col)
        return entry;
      if (entry->column == Entry::invalid)
        return nullptr;
    }
  return nullptr;
}



template <typename number>
inline const typename SparseMatrixEZ<number>::Entry *
SparseMatrixEZ<number>::locate(const size_type row, const size_type col) const
{
  SparseMatrixEZ<number> *t = const_cast<SparseMatrixEZ<number> *>(this);
  return t->locate(row, col);
}


template <typename number>
inline typename SparseMatrixEZ<number>::Entry *
SparseMatrixEZ<number>::allocate(const size_type row, const size_type col)
{
  AssertIndexRange(row, m());
  AssertIndexRange(col, n());

  RowInfo &       r   = row_info[row];
  const size_type end = r.start + r.length;

  size_type i = r.start;
  // If diagonal exists and this
  // column is higher, start only
  // after diagonal.
  if (r.diagonal != RowInfo::invalid_diagonal && col >= row)
    i += r.diagonal;
  // Find position of entry
  while (i < end && data[i].column < col)
    ++i;

  // entry found
  if (i != end && data[i].column == col)
    return &data[i];

  // Now, we must insert the new
  // entry and move all successive
  // entries back.

  // If no more space is available
  // for this row, insert new
  // elements into the vector.
  // TODO:[GK] We should not extend this row if i<end
  if (row != row_info.size() - 1)
    {
      if (end >= row_info[row + 1].start)
        {
          // Failure if increment 0
          Assert(increment != 0, ExcEntryAllocationFailure(row, col));

          // Insert new entries
          data.insert(data.begin() + end, increment, Entry());
          // Update starts of
          // following rows
          for (size_type rn = row + 1; rn < row_info.size(); ++rn)
            row_info[rn].start += increment;
        }
    }
  else
    {
      if (end >= data.size())
        {
          // Here, appending a block
          // does not increase
          // performance.
          data.push_back(Entry());
        }
    }

  Entry *entry = &data[i];
  // Save original entry
  Entry temp = *entry;
  // Insert new entry here to
  // make sure all entries
  // are ordered by column
  // index
  entry->column = col;
  entry->value  = 0;
  // Update row_info
  ++r.length;
  if (col == row)
    r.diagonal = i - r.start;
  else if (col < row && r.diagonal != RowInfo::invalid_diagonal)
    ++r.diagonal;

  if (i == end)
    return entry;

  // Move all entries in this
  // row up by one
  for (size_type j = i + 1; j < end; ++j)
    {
      // There should be no invalid
      // entry below end
      Assert(data[j].column != Entry::invalid, ExcInternalError());

      // TODO[GK]: This could be done more efficiently by moving starting at the
      // top rather than swapping starting at the bottom
      std::swap(data[j], temp);
    }
  Assert(data[end].column == Entry::invalid, ExcInternalError());

  data[end] = temp;

  return entry;
}



template <typename number>
inline void
SparseMatrixEZ<number>::set(const size_type i,
                            const size_type j,
                            const number    value,
                            const bool      elide_zero_values)
{
  AssertIsFinite(value);

  AssertIndexRange(i, m());
  AssertIndexRange(j, n());

  if (elide_zero_values && value == 0.)
    {
      Entry *entry = locate(i, j);
      if (entry != nullptr)
        entry->value = 0.;
    }
  else
    {
      Entry *entry = allocate(i, j);
      entry->value = value;
    }
}



template <typename number>
inline void
SparseMatrixEZ<number>::add(const size_type i,
                            const size_type j,
                            const number    value)
{
  AssertIsFinite(value);

  AssertIndexRange(i, m());
  AssertIndexRange(j, n());

  // ignore zero additions
  if (std::abs(value) == 0.)
    return;

  Entry *entry = allocate(i, j);
  entry->value += value;
}


template <typename number>
template <typename number2>
void
SparseMatrixEZ<number>::add(const std::vector<size_type> &indices,
                            const FullMatrix<number2> &   full_matrix,
                            const bool                    elide_zero_values)
{
  // TODO: This function can surely be made more efficient
  for (size_type i = 0; i < indices.size(); ++i)
    for (size_type j = 0; j < indices.size(); ++j)
      if ((full_matrix(i, j) != 0) || (elide_zero_values == false))
        add(indices[i], indices[j], full_matrix(i, j));
}



template <typename number>
template <typename number2>
void
SparseMatrixEZ<number>::add(const std::vector<size_type> &row_indices,
                            const std::vector<size_type> &col_indices,
                            const FullMatrix<number2> &   full_matrix,
                            const bool                    elide_zero_values)
{
  // TODO: This function can surely be made more efficient
  for (size_type i = 0; i < row_indices.size(); ++i)
    for (size_type j = 0; j < col_indices.size(); ++j)
      if ((full_matrix(i, j) != 0) || (elide_zero_values == false))
        add(row_indices[i], col_indices[j], full_matrix(i, j));
}



template <typename number>
template <typename number2>
void
SparseMatrixEZ<number>::add(const size_type               row,
                            const std::vector<size_type> &col_indices,
                            const std::vector<number2> &  values,
                            const bool                    elide_zero_values)
{
  // TODO: This function can surely be made more efficient
  for (size_type j = 0; j < col_indices.size(); ++j)
    if ((values[j] != 0) || (elide_zero_values == false))
      add(row, col_indices[j], values[j]);
}



template <typename number>
template <typename number2>
void
SparseMatrixEZ<number>::add(const size_type  row,
                            const size_type  n_cols,
                            const size_type *col_indices,
                            const number2 *  values,
                            const bool       elide_zero_values,
                            const bool  /*col_indices_are_sorted*/ )
{
  // TODO: This function can surely be made more efficient
  for (size_type j = 0; j < n_cols; ++j)
    if ((std::abs(values[j]) != 0) || (elide_zero_values == false))
      add(row, col_indices[j], values[j]);
}



template <typename number>
inline number
SparseMatrixEZ<number>::el(const size_type i, const size_type j) const
{
  const Entry *entry = locate(i, j);
  if (entry)
    return entry->value;
  return 0.;
}



template <typename number>
inline number
SparseMatrixEZ<number>::operator()(const size_type i, const size_type j) const
{
  const Entry *entry = locate(i, j);
  if (entry)
    return entry->value;
  Assert(false, ExcInvalidEntry(i, j));
  return 0.;
}


template <typename number>
inline typename SparseMatrixEZ<number>::const_iterator
SparseMatrixEZ<number>::begin() const
{
  const_iterator result(this, 0, 0);
  return result;
}

template <typename number>
inline typename SparseMatrixEZ<number>::const_iterator
SparseMatrixEZ<number>::end() const
{
  return const_iterator(this, m(), 0);
}

template <typename number>
inline typename SparseMatrixEZ<number>::const_iterator
SparseMatrixEZ<number>::begin(const size_type r) const
{
  AssertIndexRange(r, m());
  const_iterator result(this, r, 0);
  return result;
}

template <typename number>
inline typename SparseMatrixEZ<number>::const_iterator
SparseMatrixEZ<number>::end(const size_type r) const
{
  AssertIndexRange(r, m());
  const_iterator result(this, r + 1, 0);
  return result;
}

template <typename number>
template <typename MatrixType>
inline SparseMatrixEZ<number> &
SparseMatrixEZ<number>::copy_from(const MatrixType &M,
                                  const bool        elide_zero_values)
{
  reinit(M.m(), M.n(), this->saved_default_row_length, this->increment);

  // loop over the elements of the argument matrix row by row, as suggested
  // in the documentation of the sparse matrix iterator class, and
  // copy them into the current object
  for (size_type row = 0; row < M.m(); ++row)
    {
      const typename MatrixType::const_iterator end_row = M.end(row);
      for (typename MatrixType::const_iterator entry = M.begin(row);
           entry != end_row;
           ++entry)
        set(row, entry->column(), entry->value(), elide_zero_values);
    }

  return *this;
}

template <typename number>
template <typename MatrixType>
inline void
SparseMatrixEZ<number>::add(const number factor, const MatrixType &M)
{
  Assert(M.m() == m(), ExcDimensionMismatch(M.m(), m()));
  Assert(M.n() == n(), ExcDimensionMismatch(M.n(), n()));

  if (factor == 0.)
    return;

  // loop over the elements of the argument matrix row by row, as suggested
  // in the documentation of the sparse matrix iterator class, and
  // add them into the current object
  for (size_type row = 0; row < M.m(); ++row)
    {
      const typename MatrixType::const_iterator end_row = M.end(row);
      for (typename MatrixType::const_iterator entry = M.begin(row);
           entry != end_row;
           ++entry)
        if (entry->value() != 0)
          add(row, entry->column(), factor * entry->value());
    }
}



template <typename number>
template <typename MatrixTypeA, typename MatrixTypeB>
inline void
SparseMatrixEZ<number>::conjugate_add(const MatrixTypeA &A,
                                      const MatrixTypeB &B,
                                      const bool         transpose)
{
  // Compute the result
  // r_ij = \sum_kl b_ik b_jl a_kl

  //    Assert (n() == B.m(), ExcDimensionMismatch(n(), B.m()));
  //    Assert (m() == B.m(), ExcDimensionMismatch(m(), B.m()));
  //    Assert (A.n() == B.n(), ExcDimensionMismatch(A.n(), B.n()));
  //    Assert (A.m() == B.n(), ExcDimensionMismatch(A.m(), B.n()));

  // Somehow, we have to avoid making
  // this an operation of complexity
  // n^2. For the transpose case, we
  // can go through the non-zero
  // elements of A^-1 and use the
  // corresponding rows of B only.
  // For the non-transpose case, we
  // must find a trick.
  typename MatrixTypeB::const_iterator       b1      = B.begin();
  const typename MatrixTypeB::const_iterator b_final = B.end();
  if (transpose)
    while (b1 != b_final)
      {
        const size_type                      i  = b1->column();
        const size_type                      k  = b1->row();
        typename MatrixTypeB::const_iterator b2 = B.begin();
        while (b2 != b_final)
          {
            const size_type j = b2->column();
            const size_type l = b2->row();

            const typename MatrixTypeA::value_type a = A.el(k, l);

            if (a != 0.)
              add(i, j, a * b1->value() * b2->value());
            ++b2;
          }
        ++b1;
      }
  else
    {
      // Determine minimal and
      // maximal row for a column in
      // advance.

      std::vector<size_type> minrow(B.n(), B.m());
      std::vector<size_type> maxrow(B.n(), 0);
      while (b1 != b_final)
        {
          const size_type r = b1->row();
          if (r < minrow[b1->column()])
            minrow[b1->column()] = r;
          if (r > maxrow[b1->column()])
            maxrow[b1->column()] = r;
          ++b1;
        }

      typename MatrixTypeA::const_iterator       ai = A.begin();
      const typename MatrixTypeA::const_iterator ae = A.end();

      while (ai != ae)
        {
          const typename MatrixTypeA::value_type a = ai->value();
          // Don't do anything if
          // this entry is zero.
          if (a == 0.)
            continue;

          // Now, loop over all rows
          // having possibly a
          // nonzero entry in column
          // ai->row()
          b1 = B.begin(minrow[ai->row()]);
          const typename MatrixTypeB::const_iterator be1 =
            B.end(maxrow[ai->row()]);
          const typename MatrixTypeB::const_iterator be2 =
            B.end(maxrow[ai->column()]);

          while (b1 != be1)
            {
              const double b1v = b1->value();
              // We need the product
              // of both. If it is
              // zero, we can save
              // the work
              if (b1->column() == ai->row() && (b1v != 0.))
                {
                  const size_type i = b1->row();

                  typename MatrixTypeB::const_iterator b2 =
                    B.begin(minrow[ai->column()]);
                  while (b2 != be2)
                    {
                      if (b2->column() == ai->column())
                        {
                          const size_type j = b2->row();
                          add(i, j, a * b1v * b2->value());
                        }
                      ++b2;
                    }
                }
              ++b1;
            }
          ++ai;
        }
    }
}


template <typename number>
template <class StreamType>
inline void
SparseMatrixEZ<number>::print_statistics(StreamType &out, bool full)
{
  size_type              used;
  size_type              allocated;
  size_type              reserved;
  std::vector<size_type> used_by_line;

  compute_statistics(used, allocated, reserved, used_by_line, full);

  out << "SparseMatrixEZ:used      entries:" << used << std::endl
      << "SparseMatrixEZ:allocated entries:" << allocated << std::endl
      << "SparseMatrixEZ:reserved  entries:" << reserved << std::endl;

  if (full)
    {
      for (size_type i = 0; i < used_by_line.size(); ++i)
        if (used_by_line[i] != 0)
          out << "SparseMatrixEZ:entries\t" << i << "\trows\t"
              << used_by_line[i] << std::endl;
    }
}


DEAL_II_NAMESPACE_CLOSE

#endif
 /*----------------------------   sparse_matrix.h ---------------------------*/ 


