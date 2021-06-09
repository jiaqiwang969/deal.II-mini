//include/deal.II-translator/lac/full_matrix_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2020 by the deal.II authors
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

#ifndef dealii_full_matrix_h
#define dealii_full_matrix_h


#include <deal.II/base/config.h>

#include <deal.II/base/numbers.h>
#include <deal.II/base/table.h>
#include <deal.II/base/tensor.h>

#include <deal.II/differentiation/ad/ad_number_traits.h>

#include <deal.II/lac/exceptions.h>
#include <deal.II/lac/identity_matrix.h>

#include <cstring>
#include <iomanip>
#include <vector>

DEAL_II_NAMESPACE_OPEN


// forward declarations
#ifndef DOXYGEN
template <typename number>
class Vector;
template <typename number>
class LAPACKFullMatrix;
#endif

/*!   @addtogroup  Matrix1  
     * @{ 

 
*
*/


/**
 * 实现一个经典的矩形数字方案。条目的数据类型由模板参数<tt>number</tt>提供。
 * 该接口相当肥大，事实上，每次需要新的功能时，它都会增长。所以，提供了很多功能。
 * 内部计算通常是以函数的向量参数的精度来完成的。如果没有数字类型的参数，则使用矩阵数字类型。
 *
 *
 * @note
 * 本模板的实例化提供给<tt>  @<float@>,   @<double@>,
 * @<std::complex@<float@>@>,   @<std::complex@<double@>@></tt>.
 * 其他可以在应用程序中生成，详见 @ref Instantiations  。
 *
 *
 */
template <typename number>
class FullMatrix : public Table<2, number>
{
public:
  // The assertion in full_matrix.templates.h for whether or not a number is
  // finite is not compatible for AD number types.
  static_assert(
    !Differentiation::AD::is_ad_number<number>::value,
    "The FullMatrix class does not support auto-differentiable numbers.");

  /**
   * 用于对该容器进行索引的一种类型。
   *
   */
  using size_type = std::size_t;

  /**
   * 矩阵条目的类型。这个别名类似于标准库容器中的<tt>value_type</tt>。
   *
   */
  using value_type = number;

  /**
   * 使用基类的可变迭代器类型。
   *
   */
  using iterator = typename Table<2, number>::iterator;

  /**
   * 使用基类的常数迭代器类型。
   *
   */
  using const_iterator = typename Table<2, number>::const_iterator;

  /**
   * 使用基类的迭代器函数。
   *
   */
  using Table<2, number>::begin;

  /**
   * 使用基类迭代器函数
   *
   */
  using Table<2, number>::end;

  /**
   * 声明一个类型，该类型持有与本类的模板参数相同精度的实值数。如果这个类的模板参数是一个实数数据类型，那么real_type就等于模板参数。
   * 如果模板参数是一个 std::complex
   * 类型，那么real_type等于复数的基础类型。
   * 这个别名被用来表示规范的返回类型。
   *
   */
  using real_type = typename numbers::NumberTraits<number>::real_type;

  /**
   * @name  构造函数和初始化。 也见基类Table。
   *
   */
  //@{

  /**
   * 构造函数。将矩阵初始化为一个维度为<tt>n</tt>的正方形矩阵。
   * 为了避免整数和其他类型的矩阵的隐式转换，这个构造函数被声明为<tt>explicit</tt>。
   * 默认情况下，不分配内存。
   *
   */
  explicit FullMatrix(const size_type n = 0);

  /**
   * 构造函数。将矩阵初始化为一个矩形矩阵。
   *
   */
  FullMatrix(const size_type rows, const size_type cols);

  /**
   * 构造函数从一个数组中初始化。该数组被逐行排列。不进行范围检查。
   *
   */
  FullMatrix(const size_type rows, const size_type cols, const number *entries);

  /**
   * 构建一个全矩阵，等于参数大小的身份矩阵。使用这个构造函数，我们可以很容易地创建一个大小为
   * <code>n</code> 的身份矩阵，方法是说
   * @code
   * FullMatrix<double> M(IdentityMatrix(n));
   * @endcode
   *
   *
   */
  FullMatrix(const IdentityMatrix &id);
  /**
   * @}
   *
   */

  /**
   * @name  复制到其他矩阵中或从其他矩阵中复制出来
   *
   */
  /**
   * @{
   *
   */

  /**
   * 变量赋值运算符。
   *
   */
  template <typename number2>
  FullMatrix<number> &
  operator=(const FullMatrix<number2> &);

  /**
   * 这个操作符将一个标量分配给一个矩阵。为了避免与这个函数的语义相混淆，0是<tt>d</tt>唯一允许的值，允许你以一种直观的方式清除一个矩阵。
   * @dealiiOperationIsMultithreaded
   *
   */
  FullMatrix<number> &
  operator=(const number d);

  /**
   * 复制操作符，创建一个完整的矩阵，等于参数大小的身份矩阵。这样一来，人们可以很容易地创建一个大小为
   * <code>n</code> 的身份矩阵，只要说
   * @code
   * M = IdentityMatrix(n);
   * @endcode
   *
   *
   */
  FullMatrix<number> &
  operator=(const IdentityMatrix &id);

  /**
   * LapackFullMatrix的赋值操作。调用的矩阵必须与LAPACK矩阵的大小相同。
   *
   */
  template <typename number2>
  FullMatrix<number> &
  operator=(const LAPACKFullMatrix<number2> &);


  /**
   * 来自不同矩阵类的赋值。这个赋值运算符使用类型名MatrixType的迭代器。因此，稀疏矩阵是可能的来源。
   *
   */
  template <typename MatrixType>
  void
  copy_from(const MatrixType &);

  /**
   * 来自不同矩阵类的转置赋值。这个赋值运算符使用了typename
   * MatrixType的迭代器。因此，稀疏矩阵是可能的来源。
   *
   */
  template <typename MatrixType>
  void
  copy_transposed(const MatrixType &);

  /**
   * 用从张量中提取的元素填充矩阵，取包括在<tt>r_i</tt>和<tt>r_j</tt>之间的行以及<tt>c_i</tt>和<tt>c_j</tt>之间的列。然后将得到的矩阵插入目标矩阵的<tt>(dst_r,
   * dst_c)</tt>位置，对索引进行检查。
   *
   */
  template <int dim>
  void
  copy_from(const Tensor<2, dim> &T,
            const unsigned int    src_r_i = 0,
            const unsigned int    src_r_j = dim - 1,
            const unsigned int    src_c_i = 0,
            const unsigned int    src_c_j = dim - 1,
            const size_type       dst_r   = 0,
            const size_type       dst_c   = 0);

  /**
   * 将一个子矩阵（也是矩形的）插入张量中，将其左上角的元素放在指定的位置<tt>(dst_r,
   * dst_c)</tt>，其他元素也随之插入。默认值的选择是，如果张量的大小和矩阵的大小一致，则不需要指定参数。
   *
   */
  template <int dim>
  void copy_to(Tensor<2, dim> &   T,
               const size_type    src_r_i = 0,
               const size_type    src_r_j = dim - 1,
               const size_type    src_c_i = 0,
               const size_type    src_c_j = dim - 1,
               const unsigned int dst_r   = 0,
               const unsigned int dst_c   = 0) const;

  /**
   * 将另一个矩阵的行和列的一个子集复制到当前对象中。
   * @param  矩阵 要从中提取子集的矩阵。    @param  row_index_set
   * 要提取的 @p matrix 的行的集合。    @param  column_index_set
   * 要提取的 @p matrix 的列的集合。  @pre   @p row_index_set  和
   * @p
   * column_index_set中的元素数量应等于当前对象的行和列的数量。换句话说，在这个操作中，当前对象的大小不会被调整。
   *
   */
  template <typename MatrixType, typename index_type>
  void
  extract_submatrix_from(const MatrixType &             matrix,
                         const std::vector<index_type> &row_index_set,
                         const std::vector<index_type> &column_index_set);

  /**
   * 将当前矩阵对象的元素复制到另一个矩阵的指定行和列集合中。因此，这是一个散点操作。
   * @param  row_index_set 要写入的 @p matrix 的行。    @param
   * column_index_set 要写入的 @p matrix 的列。    @param  matrix
   * 某些元素要被替换的矩阵。  @pre   @p row_index_set 和 @p
   * column_index_set中的元素数量应等于当前对象的行和列的数量。换句话说，在这个操作中，当前对象的大小不会被调整。
   *
   */
  template <typename MatrixType, typename index_type>
  void
  scatter_matrix_to(const std::vector<index_type> &row_index_set,
                    const std::vector<index_type> &column_index_set,
                    MatrixType &                   matrix) const;

  /**
   * 填充矩形块。
   * 矩阵<tt>src</tt>的一个矩形块被复制到<tt>this</tt>。被复制的块的左上角是<tt>(src_offset_i,src_offset_j)</tt>。
   * 被复制块的左上角是<tt>(dst_offset_i,dst_offset_j)</tt>。
   * 被复制的矩形块的尺寸是可能的最大尺寸，由<tt>this</tt>或<tt>src</tt>的尺寸决定。
   *
   */
  template <typename number2>
  void
  fill(const FullMatrix<number2> &src,
       const size_type            dst_offset_i = 0,
       const size_type            dst_offset_j = 0,
       const size_type            src_offset_i = 0,
       const size_type            src_offset_j = 0);


  /**
   * 使基类的功能可用。
   *
   */
  template <typename number2>
  void
  fill(const number2 *);

  /**
   * 用另一个矩阵的permutation来填充。
   * 矩阵<tt>src</tt>被复制到目标中。两个互换<tt>p_r</tt>和<tt>p_c</tt>的操作方式，使得<tt>result(i,j)
   * = src(p_r[i], p_c[j]) </tt>。
   * 如果矩阵<tt>src</tt>更大的话，向量也可能是从更大的整数集中选择出来的。也可以通过这种方法来复制行或列。
   *
   */
  template <typename number2>
  void
  fill_permutation(const FullMatrix<number2> &   src,
                   const std::vector<size_type> &p_rows,
                   const std::vector<size_type> &p_cols);

  /**
   * 将矩阵的某一特定条目设置为一个值。因此，调用
   * <code>A.set(1,2,3.141);</code>  完全等同于操作  <code>A(1,2) =
   * 3.141;</code>
   * 。这个函数的存在是为了与各种稀疏矩阵对象兼容。
   * @param  i 要设置的元素的行索引。    @param  j
   * 要设置的元素的列索引。    @param  value
   * 要写进元素的值。
   *
   */
  void
  set(const size_type i, const size_type j, const number value);
  /**
   * @}
   *
   */
  /**
   * @name  非修饰性操作符
   *
   */
  /**
   * @{
   *
   */

  /**
   * 比较运算符。对这个东西要小心，它可能会吞噬大量的计算时间
   * 它最常用于程序的内部一致性检查。
   *
   */
  bool
  operator==(const FullMatrix<number> &) const;

  /**
   * 这个矩阵的行数。 注意，这个矩阵的维度是<i>m x n</i>。
   *
   */
  size_type
  m() const;

  /**
   * 该矩阵的列数。 请注意，该矩阵的维度为<i>m x n</i>。
   *
   */
  size_type
  n() const;

  /**
   * 返回该矩阵是否只包含数值为0的元素。这个函数主要用于内部一致性检查，在非调试模式下应该很少使用，因为它需要花费相当多的时间。
   *
   */
  bool
  all_zero() const;

  /**
   * 返回由该矩阵引起的向量<tt>v</tt>的规范的平方，即<i>(v,Mv)</i>。这很有用，例如在有限元背景下，一个函数的<i>L<sup>2</sup></i>规范等于相对于代表有限元函数节点值的向量矩阵的矩阵规范。
   * 显然，对于这个操作，矩阵需要是二次的，而且为了使结果真正成为一个规范，它还需要是实数对称的或复数隐式的。
   * 该矩阵和给定向量的基础模板类型应该都是实值或复值，但不是混合的，这样这个函数才有意义。
   *
   */
  template <typename number2>
  number2
  matrix_norm_square(const Vector<number2> &v) const;

  /**
   * 建立矩阵标量乘积<tt>u<sup>T</sup>M
   * v</tt>。这个函数在有限元背景下建立两个函数的单元标量积时大多有用。
   * 这个矩阵和给定向量的基本模板类型应该都是实数或复数，但不是混合的，这样这个函数才有意义。
   *
   */
  template <typename number2>
  number2
  matrix_scalar_product(const Vector<number2> &u,
                        const Vector<number2> &v) const;

  /**
   * 返回矩阵的<i>l<sub>1</sub></i>-norm，其中 $||M||_1 = \max_j
   * \sum_i |M_{ij}|$ （列上之和的最大值）。
   *
   */
  real_type
  l1_norm() const;

  /**
   * 返回矩阵的 $l_\infty$ -Norm，其中 $||M||_\infty = \max_i \sum_j
   * |M_{ij}|$ （行上和的最大值）。
   *
   */
  real_type
  linfty_norm() const;

  /**
   * 计算矩阵的Frobenius规范。
   * 返回值是所有矩阵条目的平方和的根。
   * @note
   * 对于我们中的胆小鬼：这个规范不是与向量空间的<i>l<sub>2</sub></i>规范兼容的规范。
   *
   */
  real_type
  frobenius_norm() const;

  /**
   * 计算歪斜对称部分的相对规范。返回值是矩阵的倾斜对称部分的Frobenius规范除以矩阵的规范。
   * 这个函数的主要目的是检查，一个矩阵是否在一定的精度范围内是对称的。
   *
   */
  real_type
  relative_symmetry_norm2() const;

  /**
   * 计算矩阵的行列式。
   * 这个功能只在一维、二维和三维空间实现，因为对于更高的维度，数值工作就会爆发。
   * 很明显，对于这个函数来说，矩阵需要是二次的。
   *
   */
  number
  determinant() const;

  /**
   * 返回矩阵的轨迹，即对角线值的总和（恰好也等于矩阵的特征值的总和）。
   * 很明显，对于这个函数来说，矩阵需要是二次的。
   *
   */
  number
  trace() const;

  /**
   * 以用户定义的格式输出由指定的精度和宽度给出的矩阵。这个函数在为输出设置这些给定值之前保存了流的宽度和精度，并在输出后恢复之前的值。
   *
   */
  template <class StreamType>
  void
  print(StreamType &       s,
        const unsigned int width     = 5,
        const unsigned int precision = 2) const;

  /**
   * 打印矩阵并允许对条目进行格式化。
   * 参数允许对输出格式进行灵活设置。      @arg
   * <tt>precision</tt>表示尾数的数量。      @arg
   * <tt>scientific</tt>用于确定数字格式，其中<tt>scientific</tt> =
   * <tt>false</tt>表示固定点符号。      @arg
   * <tt>width</tt>表示每列的与。<tt>width</tt>的零条目使函数计算出一个宽度，但如果输出粗略，可以将其改为正值。
   * @arg  <tt>zero_string</tt>为零条目指定一个打印的字符串。
   * @arg  <tt>denominator</tt>
   * 将整个矩阵乘以这个共同的分母，得到更漂亮的数字。
   * @arg
   * <tt>阈值</tt>：所有绝对值小于此值的条目都被视为零。
   *
   */
  void
  print_formatted(std::ostream &     out,
                  const unsigned int precision   = 3,
                  const bool         scientific  = true,
                  const unsigned int width       = 0,
                  const char *       zero_string = " ",
                  const double       denominator = 1.,
                  const double       threshold   = 0.) const;

  /**
   * 确定这个对象的内存消耗（以字节为单位）的估计值。
   *
   */
  std::size_t
  memory_consumption() const;

  //@}
  ///@name Iterator functions
  //@{

  /**
   * 可变的迭代器，从<tt>r</tt>行的第一个条目开始。
   *
   */
  iterator
  begin(const size_type r);

  /**
   * 超过行<tt>r</tt>末尾的一个可变迭代器。
   *
   */
  iterator
  end(const size_type r);

  /**
   * 从<tt>r</tt>行的第一个条目开始的恒定迭代器。
   *
   */
  const_iterator
  begin(const size_type r) const;

  /**
   * 从第<tt>r</tt>行的最后一条开始的恒定迭代器。
   *
   */
  const_iterator
  end(const size_type r) const;

  //@}
  ///@name Modifying operators
  //@{

  /**
   * 用一个固定的因子来缩放整个矩阵。
   *
   */
  FullMatrix &
  operator*=(const number factor);

  /**
   * 用给定系数的逆数来缩放整个矩阵。
   *
   */
  FullMatrix &
  operator/=(const number factor);

  /**
   * 缩放后的矩阵的简单加法，即<tt>*this += a*A</tt>。
   * 矩阵<tt>A</tt>可以是一个任意底层标量类型上的全矩阵，只要其数据类型可以转换为该矩阵的数据类型。
   *
   */
  template <typename number2>
  void
  add(const number a, const FullMatrix<number2> &A);

  /**
   * 缩放矩阵的多重加法，即<tt>*this += a*A + b*B</tt>。
   * 矩阵<tt>A</tt>和<tt>B</tt>可以是一个任意底层标量类型上的全矩阵，只要其数据类型可以转换为该矩阵的数据类型。
   *
   */
  template <typename number2>
  void
  add(const number               a,
      const FullMatrix<number2> &A,
      const number               b,
      const FullMatrix<number2> &B);

  /**
   * 缩放矩阵的多重加法，即<tt>*this += a*A + b*B + c*C</tt>。
   * 矩阵<tt>A</tt>,
   * <tt>B</tt>和<tt>C</tt>可以是一个任意底层标量类型上的全矩阵，只要其数据类型可以转换为这个矩阵的数据类型。
   *
   */
  template <typename number2>
  void
  add(const number               a,
      const FullMatrix<number2> &A,
      const number               b,
      const FullMatrix<number2> &B,
      const number               c,
      const FullMatrix<number2> &C);

  /**
   * 添加矩形块。
   * 矩阵<tt>src</tt>的一个矩形块被添加到<tt>this</tt>。
   * 被复制的块的左上角是<tt>(src_offset_i,src_offset_j)/tt>。
   * 被复制块的左上角是<tt>(dst_offset_i,dst_offset_j)</tt>。
   * 被复制的矩形块的尺寸是可能的最大尺寸，由<tt>this</tt>或<tt>src</tt>的尺寸和给定的偏移量决定。
   *
   */
  template <typename number2>
  void
  add(const FullMatrix<number2> &src,
      const number               factor,
      const size_type            dst_offset_i = 0,
      const size_type            dst_offset_j = 0,
      const size_type            src_offset_i = 0,
      const size_type            src_offset_j = 0);

  /**
   * 将<tt>B</tt>的转置加到<tt>this</tt>。    <i>A += s
   * B<sup>T</sup></i>
   *
   */
  template <typename number2>
  void
  Tadd(const number s, const FullMatrix<number2> &B);

  /**
   * 增加一个矩形块的转置。
   * 矩阵<tt>src</tt>的一个矩形块被转置并添加到<tt>this</tt>。被复制的块的左上角是<tt>(src_offset_i,src_offset_j)</tt>在<b>non</b>转置的矩阵的坐标。
   * 被复制块的左上角是<tt>(dst_offset_i,dst_offset_j)</tt>。
   * 被复制的矩形块的尺寸是可能的最大尺寸，由<tt>this</tt>或<tt>src</tt>的尺寸决定。
   *
   */
  template <typename number2>
  void
  Tadd(const FullMatrix<number2> &src,
       const number               factor,
       const size_type            dst_offset_i = 0,
       const size_type            dst_offset_j = 0,
       const size_type            src_offset_i = 0,
       const size_type            src_offset_j = 0);

  /**
   * 在给定的位置添加一个单一元素。
   *
   */
  void
  add(const size_type row, const size_type column, const number value);

  /**
   * 在给定的全局矩阵行中，在全局矩阵中由col_indices指定的列上添加一个由<tt>values</tt>给出的数值数组。这个函数的出现是为了与deal.II中的各种稀疏矩阵兼容。特别是，两个布尔字段
   * @p elide_zero_values 和 @p col_indices_are_sorted
   * 并不影响这个例程的性能，与稀疏矩阵的情况不同，在实现中确实被忽略了。
   *
   */
  template <typename number2, typename index_type>
  void
  add(const size_type   row,
      const size_type   n_cols,
      const index_type *col_indices,
      const number2 *   values,
      const bool        elide_zero_values      = true,
      const bool        col_indices_are_sorted = false);

  /**
   * <i>A(i,1...n) += s*A(j,1...n)</i>.  简单地增加这一行的数量
   *
   */
  void
  add_row(const size_type i, const number s, const size_type j);

  /**
   * <i>A(i,1...n) += s*A(j,1...n) + t*A(k,1...n)</i>.
   * 多次增加这的行数。
   *
   */
  void
  add_row(const size_type i,
          const number    s,
          const size_type j,
          const number    t,
          const size_type k);

  /**
   * <i>A(1...n,i) += s*A(1...n,j)</i>.  简单添加这的列。
   *
   */
  void
  add_col(const size_type i, const number s, const size_type j);

  /**
   * <i>A(1...n,i) += s*A(1...n,j) + t*A(1...n,k)</i>.  多次增加此列。
   *
   */
  void
  add_col(const size_type i,
          const number    s,
          const size_type j,
          const number    t,
          const size_type k);

  /**
   * 交换 <i>A(i,1...n) <-> A(j,1...n)</i>.
   * 交换这个的第i行和第j行
   *
   */
  void
  swap_row(const size_type i, const size_type j);

  /**
   * 交换<i>A(1...n,i) <-> A(1...n,j)</i>。
   * 交换这个的第i列和第j列
   *
   */
  void
  swap_col(const size_type i, const size_type j);

  /**
   * 给这个的对角线元素添加常数，即添加身份矩阵的一个倍数。
   *
   */
  void
  diagadd(const number s);

  /**
   * 赋值 <tt>*this = a*A</tt>.
   *
   */
  template <typename number2>
  void
  equ(const number a, const FullMatrix<number2> &A);

  /**
   * 分配 <tt>*this = a*A + b*B</tt>.
   *
   */
  template <typename number2>
  void
  equ(const number               a,
      const FullMatrix<number2> &A,
      const number               b,
      const FullMatrix<number2> &B);

  /**
   * 分配 <tt>*this = a*A + b*B + c*C</tt>.
   *
   */
  template <typename number2>
  void
  equ(const number               a,
      const FullMatrix<number2> &A,
      const number               b,
      const FullMatrix<number2> &B,
      const number               c,
      const FullMatrix<number2> &C);

  /**
   * 通过形成现有矩阵和它的转置<i>A =
   * 1/2(A+A<sup>T</sup>)</i>之间的平均值，对矩阵进行对称。
   * 很明显，矩阵必须是二次方的，才能进行这一操作。
   *
   */
  void
  symmetrize();

  /**
   * A=Inverse(A)。A必须是一个方形矩阵。
   * 通过高斯-乔丹算法对该矩阵进行反演，并进行部分透视。
   * 这个过程对于正定矩阵来说表现良好，但要注意在不确定的情况下的舍入误差。
   * 在deal.II被配置为LAPACK的情况下，函数Xgetrf和Xgetri建立了一个LU因式分解，并在该因式分解的基础上反转矩阵，提供了最好的性能，直到有几百个行和列的矩阵。
   * 对一个<tt>n x
   * n</tt>矩阵进行反转的数值努力是<tt>n*3</tt>。
   *
   */
  void
  gauss_jordan();

  /**
   * 将给定矩阵的逆运算分配给<tt>*this</tt>。这个函数是针对一到四维的二次元矩阵的硬编码。然而，由于所需的代码量增长很快，如果维数更大，则隐含地调用gauss_jordan()方法。
   *
   */
  template <typename number2>
  void
  invert(const FullMatrix<number2> &M);

  /**
   * 将给定矩阵 $A$ 的Cholesky分解 $A=:L L^T$
   * 分配给<tt>*this</tt>，其中 $L$
   * 为下三角矩阵。给定的矩阵必须是对称正定的。
   * 如果矩阵不是正定的，就会抛出ExcMatrixNotPositiveDefinite。
   *
   */
  template <typename number2>
  void
  cholesky(const FullMatrix<number2> &A);

  /**
   * <tt>*this(i,j)</tt> =  $V(i) W(j)$  其中 $V,W$
   * 是相同长度的向量。
   *
   */
  template <typename number2>
  void
  outer_product(const Vector<number2> &V, const Vector<number2> &W);

  /**
   * 将给定矩阵的左逆分配给<tt>*this</tt>。正在进行的计算是<i>(A<sup>T</sup>*A)<sup>-1</sup>
   * A<sup>T</sup></i>。
   *
   */
  template <typename number2>
  void
  left_invert(const FullMatrix<number2> &M);

  /**
   * 将给定矩阵的右逆分配给<tt>*this</tt>。正在执行的计算是<i>A<sup>T</sup>*(A*A<sup>T</sup>)
   * <sup>-1</sup></i>。
   *
   */
  template <typename number2>
  void
  right_invert(const FullMatrix<number2> &M);

  //@}
  ///@name Multiplications
  //@{

  /**
   * 矩阵-矩阵-乘法。
   * 可选参数<tt>adding</tt>决定，结果是存储在<tt>C</tt>中，还是添加到<tt>C</tt>中。
   * if (adding) <i>C += A*B</i> if (!adding) <i>C = A*B</i>
   * 假设<tt>A</tt>和<tt>B</tt>的大小兼容，并且<tt>C</tt>已经具有正确大小。
   * 如果三个矩阵维度的乘积大于300，并且在配置过程中检测到BLAS，则该函数使用BLAS函数Xgemm。使用BLAS通常会带来相当大的性能提升。
   *
   */
  template <typename number2>
  void
  mmult(FullMatrix<number2> &      C,
        const FullMatrix<number2> &B,
        const bool                 adding = false) const;

  /**
   * 矩阵-矩阵-乘法，使用<tt>this</tt>的转置。
   * 可选的参数<tt>adding</tt>决定了结果是存储在<tt>C</tt>中还是添加到<tt>C</tt>中。
   * if (adding) <i>C += A<sup>T</sup>*B</i> if (!adding) <i>C =
   * A<sup>T</sup>*B</i>
   * 假设<tt>A</tt>和<tt>B</tt>具有兼容的大小，并且<tt>C</tt>已经具有正确大小。
   * 如果三个矩阵维度的乘积大于300，并且在配置过程中检测到BLAS，则该函数使用BLAS函数Xgemm。使用BLAS通常会带来相当大的性能提升。
   *
   */
  template <typename number2>
  void
  Tmmult(FullMatrix<number2> &      C,
         const FullMatrix<number2> &B,
         const bool                 adding = false) const;

  /**
   * 使用<tt>B</tt>的转置进行矩阵-矩阵乘法。
   * 可选的参数<tt>adding</tt>决定了结果是存储在<tt>C</tt>中还是添加到<tt>C</tt>中。
   * if (adding) <i>C += A*B<sup>T</sup></i> if (!adding) <i>C =
   * A*B<sup>T</sup></i>
   * 假设<tt>A</tt>和<tt>B</tt>具有兼容的大小，并且<tt>C</tt>已经具有正确大小。
   * 如果三个矩阵维度的乘积大于300，并且在配置过程中检测到BLAS，则该函数使用BLAS函数Xgemm。使用BLAS通常会带来相当大的性能提升。
   *
   */
  template <typename number2>
  void
  mTmult(FullMatrix<number2> &      C,
         const FullMatrix<number2> &B,
         const bool                 adding = false) const;

  /**
   * 使用<tt>this</tt>和<tt>B</tt>的转置进行矩阵-矩阵乘法。
   * 可选的参数<tt>adding</tt>决定了结果是存储在<tt>C</tt>中还是添加到<tt>C</tt>中。
   * if (adding) <i>C += A<sup>T</sup>*B<sup>T</sup></i> if (!adding) <i>C =
   * A<sup>T</sup>*B<sup>T</sup></i>
   * 假设<tt>A</tt>和<tt>B</tt>具有兼容的大小，并且<tt>C</tt>已经具有正确大小。
   * 如果三个矩阵维度的乘积大于300，并且在配置过程中检测到BLAS，则该函数使用BLAS函数Xgemm。使用BLAS通常会带来相当大的性能提升。
   *
   */
  template <typename number2>
  void
  TmTmult(FullMatrix<number2> &      C,
          const FullMatrix<number2> &B,
          const bool                 adding = false) const;

  /**
   * 在当前矩阵中加入三乘积<b>B A
   * D</b>。可以选择使用矩阵的转置<b>B</b>和<b>D</b>。缩放因子对整个乘积进行缩放，这在向矩阵添加三乘积的倍数时很有帮助。
   * 这个乘积是在考虑到舒尔补码<b>B<sup>T</sup> A<sup>-1</sup>
   * D</b>的情况下编写的。
   * 请注意，在这种情况下，<tt>A</tt>的参数必须是矩阵<b>A</b>的逆值。
   *
   */
  void
  triple_product(const FullMatrix<number> &A,
                 const FullMatrix<number> &B,
                 const FullMatrix<number> &D,
                 const bool                transpose_B = false,
                 const bool                transpose_D = false,
                 const number              scaling     = number(1.));

  /**
   * 矩阵-向量-乘法。
   * 可选参数<tt>adding</tt>决定，结果是存储在<tt>w</tt>中，还是添加到<tt>w</tt>中。
   * if (adding) <i>w += A*v</i> if (!adding) <i>w = A*v</i>
   * 源和目的必须不是同一个向量。
   *
   */
  template <typename number2>
  void
  vmult(Vector<number2> &      w,
        const Vector<number2> &v,
        const bool             adding = false) const;

  /**
   * 加法 矩阵-向量-乘法。 <i>w += A*v</i>
   * 源和目的必须不是同一个向量。
   *
   */
  template <typename number2>
  void
  vmult_add(Vector<number2> &w, const Vector<number2> &v) const;

  /**
   * 转置矩阵-向量-乘法。
   * 可选参数<tt>adding</tt>决定，结果是存储在<tt>w</tt>中，还是添加到<tt>w</tt>中。
   * if (adding) <i>w += A<sup>T</sup>*v</i> if (!adding) <i>w =
   * A<sup>T</sup>*v</i> 源和目的必须不是同一个向量。
   *
   */
  template <typename number2>
  void
  Tvmult(Vector<number2> &      w,
         const Vector<number2> &v,
         const bool             adding = false) const;

  /**
   * 添加转置矩阵-向量-乘法。 <i>w += A<sup>T</sup>*v</i>
   * 源和目的必须不是同一个向量。
   *
   */
  template <typename number2>
  void
  Tvmult_add(Vector<number2> &w, const Vector<number2> &v) const;

  /**
   * 应用雅可比预处理程序，将<tt>src</tt>向量的每个元素乘以各自对角线元素的逆值，并将结果与阻尼系数<tt>omega</tt>相乘。
   *
   */
  template <typename somenumber>
  void
  precondition_Jacobi(Vector<somenumber> &      dst,
                      const Vector<somenumber> &src,
                      const number              omega = 1.) const;

  /**
   * <i>dst=b-A*x</i>.
   * 残差计算，返回<i>l<sub>2</sub></i>-norm|<i>dst</i>|。
   * 源<i>x</i>和目的<i>dst</i>不能是同一个向量。
   *
   */
  template <typename number2, typename number3>
  number
  residual(Vector<number2> &      dst,
           const Vector<number2> &x,
           const Vector<number3> &b) const;

  /**
   * 正向消除下三角。
   * 对于一个给定的右手边，反转矩形矩阵的下三角。
   * 如果矩阵的列数多于行数，该函数只对左边的二次元子矩阵进行操作。如果行数多，则考虑矩阵的上二次方部分。
   * @note  对 @p dst 和 @p src. 使用同一个对象是安全的。
   *
   */
  template <typename number2>
  void
  forward(Vector<number2> &dst, const Vector<number2> &src) const;

  /**
   * 向后消除上三角。    参见forward()
   * @note  对 @p dst 和 @p src. 使用同一对象是安全的。
   *
   */
  template <typename number2>
  void
  backward(Vector<number2> &dst, const Vector<number2> &src) const;

  //@}

  /**
   * @addtogroup  Exceptions 
     * @{ 
   *
   */

  /**
   * 例外情况
   *
   */
  DeclException0(ExcEmptyMatrix);

  /**
   * 异常情况
   *
   */
  DeclException1(
    ExcNotRegular,
    number,
    << "The maximal pivot is " << arg1
    << ", which is below the threshold. The matrix may be singular.");
  /**
   * 异常情况
   *
   */
  DeclException3(ExcInvalidDestination,
                 size_type,
                 size_type,
                 size_type,
                 << "Target region not in matrix: size in this direction="
                 << arg1 << ", size of new matrix=" << arg2
                 << ", offset=" << arg3);
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
  DeclException0(ExcMatrixNotPositiveDefinite);
  //@}
};

 /**@}*/ 

#ifndef DOXYGEN
 /*-------------------------Inline functions -------------------------------*/ 



template <typename number>
inline typename FullMatrix<number>::size_type
FullMatrix<number>::m() const
{
  return this->n_rows();
}



template <typename number>
inline typename FullMatrix<number>::size_type
FullMatrix<number>::n() const
{
  return this->n_cols();
}



template <typename number>
FullMatrix<number> &
FullMatrix<number>::operator=(const number d)
{
  Assert(d == number(0), ExcScalarAssignmentOnlyForZeroValue());
  (void)d; // removes -Wunused-parameter warning in optimized mode

  if (this->n_elements() != 0)
    this->reset_values();

  return *this;
}



template <typename number>
template <typename number2>
inline void
FullMatrix<number>::fill(const number2 *src)
{
  Table<2, number>::fill(src);
}



template <typename number>
template <typename MatrixType>
void
FullMatrix<number>::copy_from(const MatrixType &M)
{
  this->reinit(M.m(), M.n());

  // loop over the elements of the argument matrix row by row, as suggested
  // in the documentation of the sparse matrix iterator class, and
  // copy them into the current object
  for (size_type row = 0; row < M.m(); ++row)
    {
      const typename MatrixType::const_iterator end_row = M.end(row);
      for (typename MatrixType::const_iterator entry = M.begin(row);
           entry != end_row;
           ++entry)
        this->el(row, entry->column()) = entry->value();
    }
}



template <typename number>
template <typename MatrixType>
void
FullMatrix<number>::copy_transposed(const MatrixType &M)
{
  this->reinit(M.n(), M.m());

  // loop over the elements of the argument matrix row by row, as suggested
  // in the documentation of the sparse matrix iterator class, and
  // copy them into the current object
  for (size_type row = 0; row < M.m(); ++row)
    {
      const typename MatrixType::const_iterator end_row = M.end(row);
      for (typename MatrixType::const_iterator entry = M.begin(row);
           entry != end_row;
           ++entry)
        this->el(entry->column(), row) = entry->value();
    }
}



template <typename number>
template <typename MatrixType, typename index_type>
inline void
FullMatrix<number>::extract_submatrix_from(
  const MatrixType &             matrix,
  const std::vector<index_type> &row_index_set,
  const std::vector<index_type> &column_index_set)
{
  AssertDimension(row_index_set.size(), this->n_rows());
  AssertDimension(column_index_set.size(), this->n_cols());

  const size_type n_rows_submatrix = row_index_set.size();
  const size_type n_cols_submatrix = column_index_set.size();

  for (size_type sub_row = 0; sub_row < n_rows_submatrix; ++sub_row)
    for (size_type sub_col = 0; sub_col < n_cols_submatrix; ++sub_col)
      (*this)(sub_row, sub_col) =
        matrix.el(row_index_set[sub_row], column_index_set[sub_col]);
}



template <typename number>
template <typename MatrixType, typename index_type>
inline void
FullMatrix<number>::scatter_matrix_to(
  const std::vector<index_type> &row_index_set,
  const std::vector<index_type> &column_index_set,
  MatrixType &                   matrix) const
{
  AssertDimension(row_index_set.size(), this->n_rows());
  AssertDimension(column_index_set.size(), this->n_cols());

  const size_type n_rows_submatrix = row_index_set.size();
  const size_type n_cols_submatrix = column_index_set.size();

  for (size_type sub_row = 0; sub_row < n_rows_submatrix; ++sub_row)
    for (size_type sub_col = 0; sub_col < n_cols_submatrix; ++sub_col)
      matrix.set(row_index_set[sub_row],
                 column_index_set[sub_col],
                 (*this)(sub_row, sub_col));
}


template <typename number>
inline void
FullMatrix<number>::set(const size_type i,
                        const size_type j,
                        const number    value)
{
  (*this)(i, j) = value;
}



template <typename number>
template <typename number2>
void
FullMatrix<number>::vmult_add(Vector<number2> &      w,
                              const Vector<number2> &v) const
{
  vmult(w, v, true);
}


template <typename number>
template <typename number2>
void
FullMatrix<number>::Tvmult_add(Vector<number2> &      w,
                               const Vector<number2> &v) const
{
  Tvmult(w, v, true);
}


//---------------------------------------------------------------------------
template <typename number>
inline typename FullMatrix<number>::iterator
FullMatrix<number>::begin(const size_type r)
{
  AssertIndexRange(r, m());
  return iterator(this, r, 0);
}



template <typename number>
inline typename FullMatrix<number>::iterator
FullMatrix<number>::end(const size_type r)
{
  AssertIndexRange(r, m());
  return iterator(this, r + 1, 0);
}



template <typename number>
inline typename FullMatrix<number>::const_iterator
FullMatrix<number>::begin(const size_type r) const
{
  AssertIndexRange(r, m());
  return const_iterator(this, r, 0);
}



template <typename number>
inline typename FullMatrix<number>::const_iterator
FullMatrix<number>::end(const size_type r) const
{
  AssertIndexRange(r, m());
  return const_iterator(this, r + 1, 0);
}



template <typename number>
inline void
FullMatrix<number>::add(const size_type r, const size_type c, const number v)
{
  AssertIndexRange(r, this->m());
  AssertIndexRange(c, this->n());

  this->operator()(r, c) += v;
}



template <typename number>
template <typename number2, typename index_type>
inline void
FullMatrix<number>::add(const size_type   row,
                        const size_type   n_cols,
                        const index_type *col_indices,
                        const number2 *   values,
                        const bool,
                        const bool)
{
  AssertIndexRange(row, this->m());
  for (size_type col = 0; col < n_cols; ++col)
    {
      AssertIndexRange(col_indices[col], this->n());
      this->operator()(row, col_indices[col]) += values[col];
    }
}


template <typename number>
template <class StreamType>
inline void
FullMatrix<number>::print(StreamType &       s,
                          const unsigned int w,
                          const unsigned int p) const
{
  Assert(!this->empty(), ExcEmptyMatrix());

  // save the state of out stream
  const std::streamsize old_precision = s.precision(p);
  const std::streamsize old_width     = s.width(w);

  for (size_type i = 0; i < this->m(); ++i)
    {
      for (size_type j = 0; j < this->n(); ++j)
        {
          s.width(w);
          s.precision(p);
          s << this->el(i, j);
        }
      s << std::endl;
    }

  // reset output format
  s.precision(old_precision);
  s.width(old_width);
}


#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif


