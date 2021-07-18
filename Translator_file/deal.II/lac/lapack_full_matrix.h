//include/deal.II-translator/lac/lapack_full_matrix_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2020 by the deal.II authors
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

#ifndef dealii_lapack_full_matrix_h
#define dealii_lapack_full_matrix_h


#include <deal.II/base/config.h>

#include <deal.II/base/smartpointer.h>
#include <deal.II/base/table.h>
#include <deal.II/base/thread_management.h>

#include <deal.II/lac/lapack_support.h>
#include <deal.II/lac/vector_memory.h>

#include <complex>
#include <memory>
#include <vector>

DEAL_II_NAMESPACE_OPEN

// forward declarations
#ifndef DOXYGEN
template <typename number>
class Vector;
template <typename number>
class BlockVector;
template <typename number>
class FullMatrix;
template <typename number>
class SparseMatrix;
#endif

/**
 * FullMatrix的一个变种，尽可能使用LAPACK函数。为了做到这一点，矩阵是以转置的顺序存储的。元素访问函数通过恢复转置来隐藏这一事实。
 *
 *
 * @note
 * 为了执行LAPACK函数，该类在私有部分包含了很多辅助数据。这些数据向量的名称通常是LAPACK文档中为参数选择的名称。
 *
 *
 * @ingroup Matrix1
 *
 *
 */
template <typename number>
class LAPACKFullMatrix : public TransposeTable<number>
{
public:
  /**
   * 声明容器尺寸的类型。
   *
   */
  using size_type = std::make_unsigned<types::blas_int>::type;

  /**
   * 构造函数。将矩阵初始化为一个维度为 @p size.
   * 的正方形矩阵，为了避免整数和其他类型的矩阵的隐式转换，这个构造函数被声明为<tt>explicit</tt>。
   * 默认情况下，不分配内存。
   *
   */
  explicit LAPACKFullMatrix(const size_type size = 0);


  /**
   * 构造函数。将矩阵初始化为一个矩形矩阵  $\rm{rows} \times
   * \rm{cols}$  。
   *
   */
  LAPACKFullMatrix(const size_type rows, const size_type cols);


  /**
   * 拷贝构造函数。这个构造函数对矩阵进行了深度复制。
   * 因此，它带来了一个可能的效率问题，例如，如果函数参数是通过值而不是通过引用传递的。
   * 不幸的是，我们不能将这个复制构造函数<tt>explicit</tt>，因为这样就不能在容器中使用这个类，例如
   * <code>std::vector</code>
   * 。因此，检查程序性能的责任必须由这个类的使用者来承担。
   *
   */
  LAPACKFullMatrix(const LAPACKFullMatrix &);

  /**
   * 赋值运算符。
   *
   */
  LAPACKFullMatrix<number> &
  operator=(const LAPACKFullMatrix<number> &);

  /**
   * 来自普通FullMatrix的赋值运算符。
   * @note
   * 由于LAPACK希望矩阵以转置的顺序出现，所以这里包括了这种转置。
   *
   */
  template <typename number2>
  LAPACKFullMatrix<number> &
  operator=(const FullMatrix<number2> &);

  /**
   * 来自普通稀疏矩阵的赋值操作。
   * @note
   * 由于LAPACK希望矩阵以转置的顺序排列，所以这里包括了这个转置。
   *
   */
  template <typename number2>
  LAPACKFullMatrix<number> &
  operator=(const SparseMatrix<number2> &);

  /**
   * 这个操作符将一个标量分配给一个矩阵。为了避免与构造函数混淆，0（当投到
   * @p number 类型时）是 @p d. 的唯一允许值。
   *
   */
  LAPACKFullMatrix<number> &
  operator=(const number d);

  /**
   * 这个操作符将所有条目乘以一个固定的 @p factor. 。
   *
   */
  LAPACKFullMatrix<number> &
  operator*=(const number factor);

  /**
   * 该运算符将所有条目除以一个固定的 @p factor. 。
   *
   */
  LAPACKFullMatrix<number> &
  operator/=(const number factor);

  /**
   * 将矩阵的一个特定条目设置为  @p value.  因此，调用
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
   * 缩放矩阵的简单加法，即  $\mathbf A \mathrel{+}= a \, \mathbf
   * B$  。
   *
   */
  void
  add(const number a, const LAPACKFullMatrix<number> &B);

  /**
   * 对一个对称矩阵进行秩-1更新  $ \mathbf A \leftarrow \mathbf A
   * + a \, \mathbf v \mathbf v^T $  。
   * 这个函数也适用于Cholesky分解。  在这种情况下，更新（
   * $a>0$ ）是通过Givens旋转进行的，而降维（ $a<0$
   * ）是通过双曲旋转进行的。请注意，后一种情况可能会导致一个负定矩阵，在这种情况下，错误将被抛出（因为Cholesky分解只对对称和正定矩阵有效）。
   *
   */
  void
  rank1_update(const number a, const Vector<number> &v);

  /**
   * 将<a href="https://en.wikipedia.org/wiki/Givens_rotation">Givens
   * rotation</a> @p csr
   * （余弦、正弦和半径的三要素，旋转矩阵 $\mathbf G$
   * 的定义见 Utilities::LinearAlgebra::givens_rotation() ）在 @p i'th
   * 和 @p k'th 单位向量所跨越的平面内应用于该矩阵。
   * 如果 @p left 是 <code>true</code> ，则旋转从左边 $\mathbf A
   * \leftarrow \mathbf G \mathbf A$ 开始，只有行 @p i 和 @p k
   * 受到影响。  否则，旋转矩阵的转置将从右边 $\mathbf A
   * \leftarrow \mathbf A \mathbf G^T$ 开始，只有列 @p i 和 @p k
   * 受到影响。
   *
   */
  void
  apply_givens_rotation(const std::array<number, 3> &csr,
                        const size_type              i,
                        const size_type              k,
                        const bool                   left = true);

  /**
   * 从不同的矩阵类进行赋值，执行LAPACK所期望的转置格式的常规转换。这个赋值运算符使用类型名MatrixType的迭代器。因此，稀疏矩阵是可能的来源。
   *
   */
  template <typename MatrixType>
  void
  copy_from(const MatrixType &);

  /**
   * 用一个具有相同属性的矩阵重新生成当前的矩阵，如果它是由本类的构造函数创建的，参数列表与本函数相同。
   *
   */
  void
  reinit(const size_type size);

  /**
   * 和上面一样，但在调整大小时将保留矩阵的值。
   * 矩阵的原始值在增加大小时被保留 \f[ \mathbf A \rightarrow
   * \left( \begin{array}{cc} \mathbf A & \mathbf 0 \\ \mathbf 0 & \mathbf 0
   * \end{array} \right) \f]
   * 而如果新的大小较小，矩阵将包含原始矩阵的左上块 \f[
   * \left( \begin{array}{cc} \mathbf A_{11} & \mathbf A_{12} \\ \mathbf
   * A_{21} & \mathbf A_{22} \end{array} \right) \rightarrow \mathbf A_{11}
   * \f] 。
   *
   */
  void
  grow_or_shrink(const size_type size);

  /**
   * 从矩阵中移除行 @p row 和列 @p col 。  \f[ \left(
   * \begin{array}{ccc} \mathbf A_{11}   & \mathbf a_{12} & \mathbf A_{13}
   * \\ \mathbf a_{21}^T & a_{22}         & \mathbf a_{23}^T \\ \mathbf A_{31}
   * & \mathbf a_{32} & \mathbf A_{33} \end{array} \right) \rightarrow \left(
   * \begin{array}{cc} \mathbf A_{11} & \mathbf A_{13} \\ \mathbf A_{31} &
   * \mathbf A_{33} \end{array} \right) \f]
   *
   */
  void
  remove_row_and_column(const size_type row, const size_type col);

  /**
   * 用一个具有相同属性的矩阵重新生成当前的矩阵，如果它是由本类的构造函数创建的，参数列表与本函数相同。
   *
   */
  void
  reinit(const size_type rows, const size_type cols);

  /**
   * 将 @p property 分配给这个矩阵。
   *
   */
  void
  set_property(const LAPACKSupport::Property property);

  /**
   * 返回共域（或范围）空间的维度。
   * @note  矩阵的维度为 $m \times n$  。
   *
   */
  size_type
  m() const;

  /**
   * 返回域空间的维数。
   * @note  矩阵的维数是 $m \times n$  .
   *
   */
  size_type
  n() const;

  /**
   * 填充矩形块。
   * 矩阵<tt>src</tt>的一个矩形块被复制到<tt>this</tt>。被复制的块的左上角是<tt>(src_offset_i,src_offset_j)</tt>。
   * 被复制块的左上角是<tt>(dst_offset_i,dst_offset_j)</tt>。
   * 被复制的矩形块的尺寸是可能的最大尺寸，由<tt>this</tt>或<tt>src</tt>的尺寸决定。
   * 最后两个参数允许输入源的倍数或其转置。
   *
   */
  template <typename MatrixType>
  void
  fill(const MatrixType &src,
       const size_type   dst_offset_i = 0,
       const size_type   dst_offset_j = 0,
       const size_type   src_offset_i = 0,
       const size_type   src_offset_j = 0,
       const number      factor       = 1.,
       const bool        transpose    = false);


  /**
   * 矩阵-向量-乘法。    根据#state中记录的以前的变换，这个函数的结果是 <ul>   <li>  如果#state是 LAPACKSupport::matrix 或 LAPACKSupport::inverse_matrix, ，这将是一个使用LAPACK gemv()的常规矩阵向量乘积。    <li>  如果#state是 LAPACKSupport::svd 或 LAPACKSupport::inverse_svd, ，该函数首先与右变换矩阵相乘，然后与奇异值的对角线矩阵或其倒数值相乘，最后与左变换矩阵相乘。    </ul>  可选参数 @p adding 决定了结果是存储在向量中 $\mathbf w = \mathbf A \cdot \mathbf v$ 还是添加到向量中 $\mathbf w \mathrel{+}= \mathbf A \cdot \mathbf v$  。
   * @note  来源和目的地不能是同一个向量。
   * @note  带有 @p number2
   * 的模板只存在于与FullMatrix的编译时兼容。由于底层LAPACK接口的限制，只有
   * @p number2  =  @p number
   * 的情况被实现。所有其他的变体在调用时都会出现错误。
   *
   */
  template <typename number2>
  void
  vmult(Vector<number2> &      w,
        const Vector<number2> &v,
        const bool             adding = false) const;

  /**
   * 上述函数的特殊化，用于兼容 Vector::value_type. 。
   *
   */
  void
  vmult(Vector<number> &      w,
        const Vector<number> &v,
        const bool            adding = false) const;

  /**
   * 增加了矩阵-向量-乘法  $\mathbf w \mathrel{+}= \mathbf A \cdot
   * \mathbf v$  。    有关实现的细节，请参见vmult()的文档。
   *
   */
  template <typename number2>
  void
  vmult_add(Vector<number2> &w, const Vector<number2> &v) const;

  /**
   * 上述函数的特殊化，用于兼容  Vector::value_type.  。
   *
   */
  void
  vmult_add(Vector<number> &w, const Vector<number> &v) const;

  /**
   * 转置矩阵-向量-乘法。    可选参数  @p adding
   * 决定了结果是存储在向量中  $\mathbf w = \mathbf A^T \cdot
   * \mathbf v$  还是添加到向量中  $\mathbf w \mathrel{+}= \mathbf A^T
   * \cdot \mathbf v$  。
   * 关于实现的细节，请参见vmult()的文档。
   *
   */
  template <typename number2>
  void
  Tvmult(Vector<number2> &      w,
         const Vector<number2> &v,
         const bool             adding = false) const;

  /**
   * 上述函数的特殊化，用于兼容  Vector::value_type.  。
   *
   */
  void
  Tvmult(Vector<number> &      w,
         const Vector<number> &v,
         const bool            adding = false) const;

  /**
   * 增加转置矩阵-向量-乘法  $\mathbf w \mathrel{+}= \mathbf A^T
   * \cdot \mathbf v$  。
   * 有关实现的细节，请参见vmult()的文档。
   *
   */
  template <typename number2>
  void
  Tvmult_add(Vector<number2> &w, const Vector<number2> &v) const;

  /**
   * 上述函数的特殊化，用于兼容  Vector::value_type.  。
   *
   */
  void
  Tvmult_add(Vector<number> &w, const Vector<number> &v) const;


  /**
   * 矩阵-矩阵-乘法。    可选参数  @p adding
   * 决定了结果是存储在矩阵中  $\mathbf C            = \mathbf A
   * \cdot \mathbf B$  还是添加到矩阵中  $\mathbf C \mathrel{+}=
   * \mathbf A \cdot \mathbf B$  。
   * @note  假设 @p A 和 @p B 具有兼容的大小，并且 @p C
   * 已经具有正确的大小。      @p This
   * 函数使用BLAS函数Xgemm。
   *
   */
  void
  mmult(LAPACKFullMatrix<number> &      C,
        const LAPACKFullMatrix<number> &B,
        const bool                      adding = false) const;

  /**
   * 和之前一样，但将结果存储在FullMatrix中，而不是LAPACKFullMatrix中。
   *
   */
  void
  mmult(FullMatrix<number> &            C,
        const LAPACKFullMatrix<number> &B,
        const bool                      adding = false) const;

  /**
   * 使用<tt>this</tt>的转置进行矩阵-矩阵乘法。
   * 可选的参数  @p adding  决定了结果是存储在矩阵中
   * $\mathbf C = \mathbf A^T \cdot \mathbf B$  还是添加到矩阵中
   * $\mathbf C \mathrel{+}= \mathbf A^T \cdot \mathbf B$  。
   * @note  假设 @p A 和 @p B 具有兼容的大小，并且 @p C
   * 已经具有正确的大小。
   * @note  这个函数使用BLAS函数Xgemm。
   *
   */
  void
  Tmmult(LAPACKFullMatrix<number> &      C,
         const LAPACKFullMatrix<number> &B,
         const bool                      adding = false) const;

  /**
   * 和之前一样，但将结果存储在FullMatrix中，而不是LAPACKFullMatrix中。
   *
   */
  void
  Tmmult(FullMatrix<number> &            C,
         const LAPACKFullMatrix<number> &B,
         const bool                      adding = false) const;

  /**
   * 矩阵-矩阵-乘法，使用<tt>this</tt>的转置和对角线矢量  @p
   * V.  如果  <code>adding=false</code>  则将结果存储在矩阵中
   * $\mathbf C = \mathbf A^T \cdot \rm{diag}(\mathbf V) \cdot \mathbf B$
   * 否则将增加  $\mathbf C \mathrel{+}= \mathbf A^T \cdot
   * \rm{diag}(\mathbf V) \cdot \mathbf B$  。
   * @note  假设 @p A,   @p B 和 @p V 具有兼容的大小，并且 @p C
   * 已经具有正确的大小。
   * @note  这个函数不是由LAPACK提供的。该函数首先形成
   * $\rm{diag}(\mathbf V) \cdot \mathbf B$ 积，然后使用Xgemm函数。
   *
   */
  void
  Tmmult(LAPACKFullMatrix<number> &      C,
         const LAPACKFullMatrix<number> &B,
         const Vector<number> &          V,
         const bool                      adding = false) const;

  /**
   * 矩阵-矩阵-使用 @p B. 的转置进行乘法 可选参数 @p adding
   * 决定，结果是存储在矩阵中 $\mathbf C = \mathbf A \cdot \mathbf
   * B^T$ 还是加到矩阵中 $\mathbf C \mathrel{+}= \mathbf A \cdot
   * \mathbf B^T$  。
   * @note  假设 @p A 和 @p B 具有兼容的大小，并且 @p C
   * 已经具有正确的大小。
   * @note  这个函数使用BLAS函数Xgemm。
   *
   */
  void
  mTmult(LAPACKFullMatrix<number> &      C,
         const LAPACKFullMatrix<number> &B,
         const bool                      adding = false) const;

  /**
   * 和之前一样，但将结果存储在FullMatrix中，而不是LAPACKFullMatrix中。
   *
   */
  void
  mTmult(FullMatrix<number> &            C,
         const LAPACKFullMatrix<number> &B,
         const bool                      adding = false) const;

  /**
   * 使用<tt>this</tt>和 @p B. 的转置进行矩阵乘法 可选参数 @p
   * adding 决定了结果是存储在矩阵中 $\mathbf C = \mathbf A^T
   * \cdot \mathbf B^T$ 还是添加到矩阵中 $\mathbf C \mathrel{+}=
   * \mathbf A^T \cdot \mathbf B^T$  。
   * @note  假设 @p A 和 @p B 具有兼容的尺寸，并且 @p C
   * 已经具有正确的尺寸。
   * @note  这个函数使用BLAS函数Xgemm。
   *
   */
  void
  TmTmult(LAPACKFullMatrix<number> &      C,
          const LAPACKFullMatrix<number> &B,
          const bool                      adding = false) const;

  /**
   * 和之前一样，但将结果存储在FullMatrix中，而不是LAPACKFullMatrix中。
   *
   */
  void
  TmTmult(FullMatrix<number> &            C,
          const LAPACKFullMatrix<number> &B,
          const bool                      adding = false) const;

  /**
   * 执行外置换位。  矩阵 @p B 应该有适当的大小。
   * @note  对于复数类型，将执行共轭转置。
   * @note
   * 如果deal.II配置了Intel-MKL，将使用`mkl_?omatcopy`，否则将逐元素进行转置。
   *
   */
  void
  transpose(LAPACKFullMatrix<number> &B) const;

  /**
   * 将该矩阵的行数按 @p V
   * 进行缩放。这相当于与对角线矩阵进行预乘法  $\mathbf
   * A\leftarrow {\rm diag}(\mathbf V)\mathbf A$  。
   *
   */
  void
  scale_rows(const Vector<number> &V);

  /**
   * 使用LAPACK函数Xgetrf计算矩阵的LU因子化。
   *
   */
  void
  compute_lu_factorization();

  /**
   * 使用LAPACK函数Xpotrf计算矩阵的Cholesky因式分解。
   * @note 因式分解被存储在矩阵的下三角部分。
   *
   */
  void
  compute_cholesky_factorization();

  /**
   * 使用Cholesky因式分解法估计对称正定矩阵的条件数
   * $1/k(\mathbf A)$ 在 $L_1$ 准则下的倒数（ $1/(||\mathbf A||_1 \,
   * ||\mathbf A^{-1}||_1)$
   * ）。只有在矩阵已经被分解的情况下才能调用这个函数。
   * @note  条件数  $k(\mathbf A)$
   * 可用于估计与矩阵反演或线性代数方程组的解有关的数值误差，如
   * <code>error = std::numeric_limits<Number>::epsilon k</code>  。
   * 或者可以得到准确数字的数量
   * <code>std::floor(std::log10(k))</code>  。
   * @note
   * 该函数计算条件数的倒数，以避免在矩阵接近单数时可能的溢出。
   * @param[in]  l1_norm 是调用Cholesky分解前矩阵的 $L_1$
   * 规范。它可以通过调用l1_norm()获得。
   *
   */
  number
  reciprocal_condition_number(const number l1_norm) const;

  /**
   * 估计三角矩阵的条件数 $1/k(\mathbf A)$ 在 $L_1$
   * 规范中的倒数。矩阵必须将 LAPACKSupport::Property 设置为
   * LAPACKSupport::Property::upper_triangular 或
   * LAPACKSupport::Property::lower_triangular, 见set_property()。
   *
   */
  number
  reciprocal_condition_number() const;

  /**
   * 计算一个矩阵的行列式。由于它需要对矩阵进行LU因子化，这个函数只能在调用compute_lu_factorization()之后才能被调用。
   *
   */
  number
  determinant() const;

  /**
   * 计算 $L_1$ 规范。
   *
   */
  number
  l1_norm() const;

  /**
   * 计算 $L_\infty$ 规范。
   *
   */
  number
  linfty_norm() const;

  /**
   * 计算Frobenius规范
   *
   */
  number
  frobenius_norm() const;

  /**
   * 计算矩阵的轨迹，即对角线值的总和。
   * 很明显，对于这个函数来说，矩阵需要是二次的。
   *
   */
  number
  trace() const;

  /**
   * 首先用LAPACK函数Xgetrf/Xpotrf计算LU/Cholesky因子，然后用Xgetri/Xpotri建立实际的逆矩阵。
   *
   */
  void
  invert();

  /**
   * 解出右手边 @p v 的线性系统，并将解放回 @p v.
   * ，矩阵应该是三角形的，或者之前已经计算过LU/Cholesky因式分解。
   * 转置的标志表示是否要进行转置系统的求解。
   *
   */
  void
  solve(Vector<number> &v, const bool transposed = false) const;

  /**
   * 与上述相同，但对于多个右手边（与矩阵中的列数一样多
   * @p B). ）。
   *
   */
  void
  solve(LAPACKFullMatrix<number> &B, const bool transposed = false) const;

  /**
   * 计算矩阵的特征值。调用此例程后，可以使用eigenvalue()函数检索特征值。在此操作后，矩阵本身将是
   * LAPACKSupport::unusable 。
   * 可选的参数也允许计算左和右的特征向量。
   * 注意，该函数不会立即返回计算出的特征值，因为这涉及到在LAPACK函数的输出数组和任何返回数组之间复制数据。这通常是不必要的，因为人们可能不会对所有的特征值都感兴趣，而是只对极端的特征值感兴趣。在这种情况下，只让这个函数计算特征值，并有一个单独的函数来返回所要求的任何特征值，是比较便宜的做法。
   * @note 调用LAPACK函数Xgeev。
   *
   */
  void
  compute_eigenvalues(const bool right_eigenvectors = false,
                      const bool left_eigenvectors  = false);

  /**
   * 计算一个实数对称矩阵的特征值和特征向量。只计算区间
   * $(\rm{lower\_bound}, \rm{upper\_bound}]$
   * 内的特征值，绝对公差为 $\rm abs\_accuracy$
   * 。当一个近似的特征值被确定位于宽度小于或等于
   * $\rm{abs\_accuracy} + eps \rm{max}(|a|,|b|)$ 的区间 $[a,b]$
   * 内，其中 $eps$ 为机器精度，则被接受为收敛。 如果
   * $\rm{abs\_accuracy}$ 小于或等于零，那么 $eps\,|\mathbf{T}|_1$
   * 将被用来代替，其中 $|\mathbf{T}|_1$ 是通过将 $\mathbf A$
   * 还原为三对角形式得到的三对角矩阵的1-norm。当
   * $\rm{abs\_accuracy}$
   * 被设置为下溢阈值的两倍而不是零时，特征值将被最准确地计算出来。
   * 调用此例程后， $(\rm{lower\_bound}, \rm{upper\_bound}]$
   * 中的所有特征值将被存储在特征值中，相应的特征向量将被存储在特征向量的列中，其维度被相应设置。
   * @note 调用LAPACK函数Xsyevx。
   *
   */
  void
  compute_eigenvalues_symmetric(const number        lower_bound,
                                const number        upper_bound,
                                const number        abs_accuracy,
                                Vector<number> &    eigenvalues,
                                FullMatrix<number> &eigenvectors);

  /**
   * 计算形式为实数的广义对称特征问题的广义特征值和特征向量
   *
   *
   *
   *
   * - itype = 1:  $\mathbf A \cdot \mathbf x=\lambda \mathbf B \cdot \mathbf x$  。
   *
   *
   *
   *
   *
   *
   * - itype = 2:  $\mathbf A \cdot \mathbf B \cdot \mathbf x=\lambda \mathbf x$  。
   *
   *
   *
   *
   *
   *
   * - itype = 3:  $\mathbf B \cdot \mathbf A \cdot \mathbf x=\lambda \mathbf x$  其中 $\mathbf A$ 是这个矩阵。  $\mathbf A$ 和 $\mathbf B$ 被假定为对称的，而 $\mathbf B$ 必须是正定的。只有在区间 $(\rm{lower\_bound},
   * \rm{upper\_bound}]$ 内的特征值才以绝对公差
   * $\rm{abs\_accuracy}$ 进行计算。
   * 当一个近似的特征值被确定位于宽度小于或等于
   * $\rm{abs\_accuracy} + eps \rm{max}( |a|,|b| )$ 的区间 $[a,b]$
   * 内，其中 $eps$ 是机器精度，则被接受为收敛。如果
   * $\rm{abs\_accuracy}$ 小于或等于零，那么 $eps \, |\mathbf{T}|_1$
   * 将被用来代替，其中 $|\mathbf{T}|_1$ 是通过将 $\mathbf A$
   * 还原为三对角形式得到的三对角矩阵的1-norm。当
   * $\rm{abs\_accuracy}$
   * 被设置为下溢阈值的两倍而不是零时，特征值将被最准确地计算出来。调用此例程后，
   * $(\rm{lower\_bound}, \rm{upper\_bound}]$
   * 中的所有特征值将被存储在特征值中，相应的特征向量将被存储在特征向量中，其维度被相应设置。
   * @note 调用LAPACK函数Xsygvx。
   *
   */
  void
  compute_generalized_eigenvalues_symmetric(
    LAPACKFullMatrix<number> &   B,
    const number                 lower_bound,
    const number                 upper_bound,
    const number                 abs_accuracy,
    Vector<number> &             eigenvalues,
    std::vector<Vector<number>> &eigenvectors,
    const types::blas_int        itype = 1);

  /**
   * 与其他compute_generalized_eigenvalues_symmetric函数相同，只是所有的特征值都被计算出来，并且自动设置了公差。
   * 注意这个函数不会立即返回计算的特征值，因为这涉及到在LAPACK函数的输出数组和任何返回数组之间复制数据。这通常是不必要的，因为人们可能不会对所有的特征值都感兴趣，而是只对极端的特征值感兴趣。在这种情况下，只让这个函数计算特征值，并有一个单独的函数来返回所要求的任何特征值，是比较便宜的做法。特征值可以用eigenvalue()函数来检索。
   * 计算出的特征向量的数量等于eigenvectors.size()
   * @note 调用LAPACK函数Xsygv。
   *
   */
  void
  compute_generalized_eigenvalues_symmetric(
    LAPACKFullMatrix<number> &   B,
    std::vector<Vector<number>> &eigenvectors,
    const types::blas_int        itype = 1);

  /**
   * 使用LAPACK函数Xgesdd计算矩阵的奇异值分解。
   * 要求#状态为 LAPACKSupport::matrix,
   * 填充数据成员#wr、#svd_u和#svd_vt，并使对象处于#状态
   * LAPACKSupport::svd.
   * 奇异值分解将提供的矩阵（A）分解为三部分。U、sigma和V的转置（V^T），这样A=U
   * sigma
   * V^T。Sigma是一个MxN矩阵，包含A在对角线上的奇异值，而所有其他元素都是零。U是一个MxM的正交矩阵，包含与A的奇异值相对应的左侧奇异向量。
   * V是一个NxN的正交矩阵，包含与A的奇异值相对应的右侧奇异向量。注意，变量#svd_vt包含V的转置，可以通过get_svd_vt()访问，而U可以通过get_svd_u()访问。
   *
   */
  void
  compute_svd();

  /**
   * 通过奇异值分解计算矩阵的逆值。    要求#state是
   * LAPACKSupport::matrix 或 LAPACKSupport::svd.
   * 在第一种情况下，该函数调用compute_svd()。
   * 在此函数之后，对象将具有#state  LAPACKSupport::inverse_svd.
   * 对于奇异值分解，通过将所有奇异值替换为它们的倒数来简单地计算逆值。如果矩阵没有最大等级，则不触及奇异值0，从而计算矩阵的最小规范右逆。
   * 参数 @p threshold
   * 决定了何时应将奇异值视为零。它是最小的非零奇异值与最大的非零奇异值的比率
   * $s_{max}$  。因此，所有小于 $s_{max}/\rm{threshold}$
   * 的奇异值的倒数将被设置为零。
   *
   */
  void
  compute_inverse_svd(const double threshold = 0.);

  /**
   * 与上述相同，但提供内核的大小而不是阈值，即 @p
   * kernel_size 最小的特征值。
   *
   */
  void
  compute_inverse_svd_with_kernel(const unsigned int kernel_size);

  /**
   * 检索调用compute_eigenvalues()后的特征值。
   *
   */
  std::complex<number>
  eigenvalue(const size_type i) const;

  /**
   * 在调用compute_svd()或compute_inverse_svd()后检索奇异值。
   *
   */
  number
  singular_value(const size_type i) const;

  /**
   * 在调用compute_svd()或compute_inverse_svd()后检索矩阵#svd_u。
   *
   */
  inline const LAPACKFullMatrix<number> &
  get_svd_u() const;

  /**
   * 在调用compute_svd()或compute_inverse_svd()后检索矩阵#svd_vt。
   *
   */
  inline const LAPACKFullMatrix<number> &
  get_svd_vt() const;

  /**
   * 打印矩阵并允许对条目进行格式化。
   * 参数允许对输出格式进行灵活设置。      @param  out
   * 这指定了要写入的流。      @param  precision
   * 表示尾数的数量。      @param  scientific
   * 用于确定数字格式，其中 <code>scientific = false</code>
   * 表示定点符号。      @param  宽度表示每列的宽度。 @p
   * width
   * 的零条目使函数计算出一个宽度，但如果输出是粗略的，可以将其改为正值。
   * @param  zero_string指定了一个为零条目打印的字符串。
   * @param  分母
   * 将整个矩阵乘以这个共同分母，得到更漂亮的数字。
   * @param  门槛 所有绝对值小于此值的条目都被视为零。
   * @note 只有当状态为 LAPACKSupport::matrix 或
   * LAPACK::inverse_matrix.
   * 时，存储的条目才类似于一个矩阵。否则，不允许调用此函数。
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

private:
  /**
   * 计算各种规范的内部函数。
   *
   */
  number
  norm(const char type) const;

  /**
   * 由于LAPACK的操作出了名地改变了矩阵项的含义，我们在这里记录了最后一次操作后的当前状态。
   *
   */
  LAPACKSupport::State state;

  /**
   * 矩阵的其他属性，可能有助于选择更有效的LAPACK函数。
   *
   */
  LAPACKSupport::Property property;

  /**
   * 用于某些LAPACK函数的工作阵列。
   *
   */
  mutable std::vector<number> work;

  /**
   * 用于某些LAPACK函数的整数工作数组。
   *
   */
  mutable std::vector<types::blas_int> iwork;

  /**
   * 储存用于LU-因式分解中枢轴的排列方式的向量。
   * 也用作LAPACK函数中需要的刮擦数组WORK。
   *
   */
  std::vector<types::blas_int> ipiv;

  /**
   * 用于计算LU因子化的逆矩阵的工作区。
   *
   */
  std::vector<number> inv_work;

  /**
   * 特征值的实部或奇异值。由 compute_eigenvalues() 或
   * compute_svd() 填充。
   *
   */
  std::vector<typename numbers::NumberTraits<number>::real_type> wr;

  /**
   * 特征值的虚部，或者，在复数标量的情况下，特征值本身。由
   * compute_eigenvalues 填充。
   *
   */
  std::vector<number> wi;

  /**
   * 可以存储左翼特征向量的空间。
   *
   */
  std::vector<number> vl;

  /**
   * 可以存储右特征向量的空间。
   *
   */
  std::vector<number> vr;

  /**
   * 矩阵 $\mathbf U$ 的奇异值分解 $\mathbf U \cdot \mathbf S \cdot
   * \mathbf V^T$  .
   *
   */
  std::unique_ptr<LAPACKFullMatrix<number>> svd_u;

  /**
   * 矩阵 $\mathbf V^T$ 在奇异值分解 $\mathbf U \cdot \mathbf S \cdot
   * \mathbf V^T$ 中。
   *
   */
  std::unique_ptr<LAPACKFullMatrix<number>> svd_vt;

  /**
   * 线程互斥。
   *
   */
  mutable std::mutex mutex;
};



/**
 * 一个基于LAPACKFullMatrix的LU因子化的预处理程序。
 *
 *
 * @ingroup Preconditioners
 *
 *
 */
template <typename number>
class PreconditionLU : public Subscriptor
{
public:
  void
  initialize(const LAPACKFullMatrix<number> &);
  void
  initialize(const LAPACKFullMatrix<number> &, VectorMemory<Vector<number>> &);
  void
  vmult(Vector<number> &, const Vector<number> &) const;
  void
  Tvmult(Vector<number> &, const Vector<number> &) const;
  void
  vmult(BlockVector<number> &, const BlockVector<number> &) const;
  void
  Tvmult(BlockVector<number> &, const BlockVector<number> &) const;

private:
  SmartPointer<const LAPACKFullMatrix<number>, PreconditionLU<number>> matrix;
  SmartPointer<VectorMemory<Vector<number>>, PreconditionLU<number>>   mem;
};

 /*---------------------- Inline functions -----------------------------------*/ 

template <typename number>
inline void
LAPACKFullMatrix<number>::set(const size_type i,
                              const size_type j,
                              const number    value)
{
  (*this)(i, j) = value;
}


template <typename number>
inline typename LAPACKFullMatrix<number>::size_type
LAPACKFullMatrix<number>::m() const
{
  return static_cast<size_type>(this->n_rows());
}

template <typename number>
inline typename LAPACKFullMatrix<number>::size_type
LAPACKFullMatrix<number>::n() const
{
  return static_cast<size_type>(this->n_cols());
}

template <typename number>
template <typename MatrixType>
inline void
LAPACKFullMatrix<number>::copy_from(const MatrixType &M)
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

  state = LAPACKSupport::matrix;
}



template <typename number>
template <typename MatrixType>
inline void
LAPACKFullMatrix<number>::fill(const MatrixType &M,
                               const size_type   dst_offset_i,
                               const size_type   dst_offset_j,
                               const size_type   src_offset_i,
                               const size_type   src_offset_j,
                               const number      factor,
                               const bool        transpose)
{
  // loop over the elements of the argument matrix row by row, as suggested
  // in the documentation of the sparse matrix iterator class
  for (size_type row = src_offset_i; row < M.m(); ++row)
    {
      const typename MatrixType::const_iterator end_row = M.end(row);
      for (typename MatrixType::const_iterator entry = M.begin(row);
           entry != end_row;
           ++entry)
        {
          const size_type i = transpose ? entry->column() : row;
          const size_type j = transpose ? row : entry->column();

          const size_type dst_i = dst_offset_i + i - src_offset_i;
          const size_type dst_j = dst_offset_j + j - src_offset_j;
          if (dst_i < this->n_rows() && dst_j < this->n_cols())
            (*this)(dst_i, dst_j) = factor * entry->value();
        }
    }

  state = LAPACKSupport::matrix;
}


template <typename number>
template <typename number2>
void
LAPACKFullMatrix<number>::vmult(Vector<number2> &,
                                const Vector<number2> &,
                                const bool) const
{
  Assert(false,
         ExcMessage("LAPACKFullMatrix<number>::vmult must be called with a "
                    "matching Vector<double> vector type."));
}


template <typename number>
template <typename number2>
void
LAPACKFullMatrix<number>::vmult_add(Vector<number2> &,
                                    const Vector<number2> &) const
{
  Assert(false,
         ExcMessage("LAPACKFullMatrix<number>::vmult_add must be called with a "
                    "matching Vector<double> vector type."));
}


template <typename number>
template <typename number2>
void
LAPACKFullMatrix<number>::Tvmult(Vector<number2> &,
                                 const Vector<number2> &,
                                 const bool) const
{
  Assert(false,
         ExcMessage("LAPACKFullMatrix<number>::Tvmult must be called with a "
                    "matching Vector<double> vector type."));
}


template <typename number>
template <typename number2>
void
LAPACKFullMatrix<number>::Tvmult_add(Vector<number2> &,
                                     const Vector<number2> &) const
{
  Assert(false,
         ExcMessage(
           "LAPACKFullMatrix<number>::Tvmult_add must be called with a "
           "matching Vector<double> vector type."));
}


template <typename number>
inline std::complex<number>
LAPACKFullMatrix<number>::eigenvalue(const size_type i) const
{
  Assert(state & LAPACKSupport::eigenvalues, ExcInvalidState());
  Assert(wr.size() == this->n_rows(), ExcInternalError());
  Assert(wi.size() == this->n_rows(), ExcInternalError());
  AssertIndexRange(i, this->n_rows());

  if (numbers::NumberTraits<number>::is_complex)
    return std::complex<number>(wi[i]);
  else
    return std::complex<number>(wr[i], wi[i]);
}


template <typename number>
inline number
LAPACKFullMatrix<number>::singular_value(const size_type i) const
{
  Assert(state == LAPACKSupport::svd || state == LAPACKSupport::inverse_svd,
         LAPACKSupport::ExcState(state));
  AssertIndexRange(i, wr.size());

  return wr[i];
}


template <typename number>
inline const LAPACKFullMatrix<number> &
LAPACKFullMatrix<number>::get_svd_u() const
{
  Assert(state == LAPACKSupport::svd || state == LAPACKSupport::inverse_svd,
         LAPACKSupport::ExcState(state));

  return *svd_u;
}


template <typename number>
inline const LAPACKFullMatrix<number> &
LAPACKFullMatrix<number>::get_svd_vt() const
{
  Assert(state == LAPACKSupport::svd || state == LAPACKSupport::inverse_svd,
         LAPACKSupport::ExcState(state));

  return *svd_vt;
}



DEAL_II_NAMESPACE_CLOSE

#endif


