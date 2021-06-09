//include/deal.II-translator/lac/sparse_direct_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2001 - 2020 by the deal.II authors
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

#ifndef dealii_sparse_direct_h
#define dealii_sparse_direct_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/subscriptor.h>

#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparse_matrix_ez.h>
#include <deal.II/lac/vector.h>

#ifdef DEAL_II_WITH_UMFPACK
#  include <umfpack.h>
#endif

DEAL_II_NAMESPACE_OPEN

namespace types
{
  /**
   * UMFPACK的索引类型。SuiteSparse_long必须在这里用于Windows
   * 64构建。
   *
   */
#ifdef SuiteSparse_long
  using suitesparse_index = SuiteSparse_long;
#else
  using suitesparse_index = long int;
#endif
} // namespace types

/**
 * 这个类为稀疏直接求解器UMFPACK提供了一个接口，它是SuiteSparse库的一部分（见<a
 * href="http://faculty.cse.tamu.edu/davis/suitesparse.html">this
 * link</a>）。UMFPACK是一套用于解决非对称稀疏线性系统的例程，Ax=b，使用非对称模式的多前沿方法和直接稀疏LU分解。矩阵可以有对称或不对称的稀疏模式，也可以有不对称的条目。这个类的使用在
 * step-22 和 step-29 的教程程序中解释。
 * 这个矩阵类实现了预处理程序的常规接口，即用于初始化的函数initialize(const
 * SparseMatrix<double>&matrix,const
 * AdditionalData)和所有矩阵共有的整套vmult()函数。这里只实现了vmult()和vmult_add()，它们执行与逆矩阵的乘法。此外，该类还提供了一个较早的接口，由函数
 * factorize() 和 solve() 组成。这两个接口是可以互换的。
 *
 *
 * @note  如果在配置过程中没有明确关闭<a
 * href="http://faculty.cse.tamu.edu/davis/suitesparse.html">UMFPACK</a>接口，则存在该类。
 *
 *
 * @note
 * UMFPACK有自己的许可证，独立于deal.II的许可证。如果你想使用UMFPACK，你必须接受该许可证。它可以从deal.II的ReadMe文件中链接到。UMFPACK是由其作者<a
 * href="http://faculty.cse.tamu.edu/davis/welcome.html">Timothy A.
 * Davis</a>提供的。
 *
 *  <h4>Instantiations</h4> 这个类有SparseMatrix<double>,
 * SparseMatrix<float>, SparseMatrixEZ<float>, SparseMatrixEZ<double>,
 * BlockSparseMatrix<double>, 和BlockSparseMatrix<float>的实例化。
 *
 *
 * @ingroup Solvers Preconditioners
 *
 *
 */
class SparseDirectUMFPACK : public Subscriptor
{
public:
  /**
   * 声明容器尺寸的类型。
   *
   */
  using size_type = types::global_dof_index;

  /**
   * 先决条件器的通常初始化接口所需的假类。
   *
   */
  class AdditionalData
  {};


  /**
   * 构造函数。关于这个函数的参数的含义，请参见这个类的文档。
   *
   */
  SparseDirectUMFPACK();

  /**
   * 解构器。
   *
   */
  ~SparseDirectUMFPACK() override;

  /**
   * @name  设置一个稀疏因子化
   *
   */
  /**
   * @{
   *
   */

  /**
   * 这个函数什么都不做。它在这里只是为了提供一个与其他稀疏直接求解器一致的接口。
   *
   */
  void
  initialize(const SparsityPattern &sparsity_pattern);

  /**
   * 对矩阵进行因子化。在这个类的对象被初始化为某种稀疏模式后，这个函数可以针对不同的矩阵被多次调用。因此，如果你想用相同的稀疏模式反转几个矩阵，你可以节省一些计算时间。然而，请注意，大部分的计算时间实际上是花在因式分解上的，所以这个功能可能并不总是有很大的好处。
   * 与其他直接求解器类相比，初始化方法什么都不做。因此，当初始化步骤尚未执行时，初始化不会被这个方法自动调用。
   * 这个函数将矩阵的内容复制到它自己的存储器中；因此，在这个操作之后，即使需要进行后续的求解，矩阵也可以被删除。
   *
   */
  template <class Matrix>
  void
  factorize(const Matrix &matrix);

  /**
   * 初始化内存并调用 SparseDirectUMFPACK::factorize. 。
   *
   */
  template <class Matrix>
  void
  initialize(const Matrix &       matrix,
             const AdditionalData additional_data = AdditionalData());

  /**
   * @}
   *
   */

  /**
   * @name  表示矩阵的逆的函数
   *
   */
  /**
   * @{
   *
   */

  /**
   * 预调器接口函数。通常情况下，给定源向量，该方法返回<i>Ax
   * =
   * b</i>的近似解。由于这个类提供了一个直接求解器的包装，这里实际上是精确解（当然是在数值精度范围内的精确）。
   * 换句话说，这个函数实际上是与矩阵的精确逆值相乘，
   * $A^{-1}$  .
   *
   */
  void
  vmult(Vector<double> &dst, const Vector<double> &src) const;

  /**
   * 和以前一样，但对于块状向量。
   *
   */
  void
  vmult(BlockVector<double> &dst, const BlockVector<double> &src) const;

  /**
   * 同前，但使用矩阵的转置，即此函数与  $A^{-T}$  相乘。
   *
   */
  void
  Tvmult(Vector<double> &dst, const Vector<double> &src) const;

  /**
   * 与之前相同，但对于块状向量
   *
   */
  void
  Tvmult(BlockVector<double> &dst, const BlockVector<double> &src) const;

  /**
   * 返回共域（或范围）空间的维数。注意，矩阵的维度是
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
   * @}
   *
   */

  /**
   * @name  解决线性系统的函数
   *
   */
  /**
   * @{
   *
   */

  /**
   * 求解某一个右手边的向量。在矩阵被分解后，这个函数可以针对不同的右手向量被多次调用。这可以节省大量的计算时间，因为与矩阵的因式分解相比，实际的解是很快的。
   * 解决方案将被返回，以代替右手边的向量。
   * @param[in,out]  rhs_and_solution 调用此函数时包含线性系统
   * $Ax=b$ 的右手边 $b$ ，调用此函数后包含线性系统的解 $x$
   * 的一个向量。    @param[in]  转置 如果设置为
   * "true"，该函数将解决线性  $A^T x = b$  而不是  $Ax=b$  。
   * @pre  在调用这个函数之前，你需要先调用 factorize() 。
   *
   */
  void
  solve(Vector<double> &rhs_and_solution, const bool transpose = false) const;

  /**
   * 和前面的函数一样，但是对于复值的右手边和解向量。
   * 如果之前被分解的矩阵有复值项，那么在这个函数返回时，`rhs_and_solution`向量将简单地包含线性系统的解
   * $Ax=b$
   * 。如果矩阵是实值的，那么这也是正确的，但解决方案将简单地通过将
   * $A^{-1}$
   * 的因子化应用于右侧向量的实部和虚部而计算出来。
   *
   */
  void
  solve(Vector<std::complex<double>> &rhs_and_solution,
        const bool                    transpose = false) const;

  /**
   * 与之前相同，但对于块状向量。
   *
   */
  void
  solve(BlockVector<double> &rhs_and_solution,
        const bool           transpose = false) const;

  /**
   * 同前，但适用于复值块状向量。
   *
   */
  void
  solve(BlockVector<std::complex<double>> &rhs_and_solution,
        const bool                         transpose = false) const;

  /**
   * 依次调用两个函数factize()和solve()，即对给定的右手边向量执行整个求解过程。
   * 解的结果将代替右手边的向量返回。
   *
   */
  template <class Matrix>
  void
  solve(const Matrix &  matrix,
        Vector<double> &rhs_and_solution,
        const bool      transpose = false);

  /**
   * 与之前相同，但对于复值的解向量。
   *
   */
  template <class Matrix>
  void
  solve(const Matrix &                matrix,
        Vector<std::complex<double>> &rhs_and_solution,
        const bool                    transpose = false);

  /**
   * 同前，但用于块状向量。
   *
   */
  template <class Matrix>
  void
  solve(const Matrix &       matrix,
        BlockVector<double> &rhs_and_solution,
        const bool           transpose = false);

  /**
   * 同前，但对于复值块向量。
   *
   */
  template <class Matrix>
  void
  solve(const Matrix &                     matrix,
        BlockVector<std::complex<double>> &rhs_and_solution,
        const bool                         transpose = false);

  /**
   * @}
   *
   */

  /**
   * 其中一个UMFPack例程出了一个错误。错误代码包含在输出中，可以在UMFPack用户手册中查询。例程的名称也包括在内，供参考。
   *
   */
  DeclException2(
    ExcUMFPACKError,
    std::string,
    int,
    << "UMFPACK routine " << arg1 << " returned error status " << arg2 << "."
    << "\n\n"
    << ("A complete list of error codes can be found in the file "
        "<bundled/umfpack/UMFPACK/Include/umfpack.h>."
        "\n\n"
        "That said, the two most common errors that can happen are "
        "that your matrix cannot be factorized because it is "
        "rank deficient, and that UMFPACK runs out of memory "
        "because your problem is too large."
        "\n\n"
        "The first of these cases most often happens if you "
        "forget terms in your bilinear form necessary to ensure "
        "that the matrix has full rank, or if your equation has a "
        "spatially variable coefficient (or nonlinearity) that is "
        "supposed to be strictly positive but, for whatever "
        "reasons, is negative or zero. In either case, you probably "
        "want to check your assembly procedure. Similarly, a "
        "matrix can be rank deficient if you forgot to apply the "
        "appropriate boundary conditions. For example, the "
        "Laplace equation for a problem where only Neumann boundary "
        "conditions are posed (or where you forget to apply Dirichlet "
        "boundary conditions) has exactly one eigenvalue equal to zero "
        "and its rank is therefore deficient by one. Finally, the matrix "
        "may be rank deficient because you are using a quadrature "
        "formula with too few quadrature points."
        "\n\n"
        "The other common situation is that you run out of memory. "
        "On a typical laptop or desktop, it should easily be possible "
        "to solve problems with 100,000 unknowns in 2d. If you are "
        "solving problems with many more unknowns than that, in "
        "particular if you are in 3d, then you may be running out "
        "of memory and you will need to consider iterative "
        "solvers instead of the direct solver employed by "
        "UMFPACK."));

private:
  /**
   * 范围空间的尺寸，即矩阵的行数。
   *
   */
  size_type n_rows;

  /**
   * 域空间的维度，即矩阵的列数。
   *
   */
  size_type n_cols;

  /**
   * UMFPACK例程分配了一些对象，在这些对象中存储了关于分解的符号和数字值的信息。这些对象的实际数据类型是不透明的，只作为无效指针传递。
   *
   */
  void *symbolic_decomposition;
  void *numeric_decomposition;

  /**
   * 释放所有还没有被释放的内存。
   *
   */
  void
  clear();

  /**
   * 确保数组Ai和Ap在每一行都是排序的。UMFPACK希望它是这样的。我们需要有三个版本的这个函数，一个用于通常的稀疏矩阵，一个用于SparseMatrixEZ，一个用于BlockSparseMatrix类。
   *
   */
  template <typename number>
  void
  sort_arrays(const SparseMatrixEZ<number> &);

  template <typename number>
  void
  sort_arrays(const SparseMatrix<number> &);

  template <typename number>
  void
  sort_arrays(const BlockSparseMatrix<number> &);

  /**
   * 我们为求解器存储数据的数组。这些在umfpack_*_symbolic()和umfpack_*_numeric()函数的描述中都有记载，但简单来说。
   *
   *
   *
   *
   *
   *
   * - `Ap`是数组，表示哪一行在`Ai`中开始。
   *
   *
   *
   *
   *
   *
   * - `Ai`是存储非零条目的列索引的数组
   *
   *
   *
   *
   *
   * - `Ax`是存储非零条目值的数组；如果矩阵是复值的，那么它存储的是实数部分
   *
   *
   *
   * - `Az`是存储非零项的虚部的数组，只在矩阵是复数值时使用。
   *
   */
  std::vector<types::suitesparse_index> Ap;
  std::vector<types::suitesparse_index> Ai;
  std::vector<double>                   Ax;
  std::vector<double>                   Az;

  /**
   * 解算器程序的控制和工作数组。
   *
   */
  std::vector<double> control;
};

DEAL_II_NAMESPACE_CLOSE

#endif // dealii_sparse_direct_h


