//include/deal.II-translator/lac/scalapack_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2020 by the deal.II authors
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

#ifndef dealii_scalapack_h
#define dealii_scalapack_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_SCALAPACK

#  include <deal.II/base/exceptions.h>
#  include <deal.II/base/mpi.h>
#  include <deal.II/base/process_grid.h>
#  include <deal.II/base/thread_management.h>

#  include <deal.II/lac/full_matrix.h>
#  include <deal.II/lac/lapack_full_matrix.h>
#  include <deal.II/lac/lapack_support.h>

#  include <mpi.h>

#  include <memory>

DEAL_II_NAMESPACE_OPEN

/**
 * 一个围绕ScaLAPACK并行密集线性代数的封装类。
 * ScaLAPACK假设矩阵是按照块-循环分解方案来分布的。一个
 * $M$ 乘以 $N$ 的矩阵首先被分解成 $\lceil M / MB \rceil$ 乘以
 * $\lceil N / NB \rceil$ 的块，然后均匀地分布在具有 $p q \le Np$
 * 进程的二维进程网格上，其中 $p,q$ 是网格尺寸， $Np$
 * 是进程的总数。参数MB和NB被称为行和列块大小，决定了块-循环分布的颗粒度。
 * 下面显示了 $10 \times 9$ 矩阵在 $3\times 3$
 * 笛卡尔进程网格上的块-循环分布，块大小为
 * $\text{MB}=\text{NB}=2$ 。
*\htmlonly <style>div.image img[src="scalapack_block_cycling.png"]{width:35%;}</style>\endhtmlonly  @image html scalapack_block_cycling.png "Block-Cyclic Distribution"
 * 请注意，进程P2、P5和P8所拥有的本地矩阵的列数为奇数，这说明
 * $N=9$ 不是 $\text{NB}=2$ 的整数倍。
 * 块大小的选择是一个折中，既要有足够大的大小来实现高效的本地/串行BLAS，又要有足够小的大小来实现良好的并行负载平衡。
 * 下面我们展示了 ScaLAPACKMatrix::invert()
 * 在多达5个节点上的强扩展实例，每个节点由两个英特尔至强2660v2
 * IvyBridge插座2.20GHz，10个核心/插座组成。计算在方形处理器网格1x1,
 * 2x2, 3x3, 4x4, 5x5, 6x6, 7x7, 8x8, 9x9, 10x10进行。
*  @image html scalapack_invert.png
 * @ingroup Matrix1
 *
 *
 */
template <typename NumberType>
class ScaLAPACKMatrix : protected TransposeTable<NumberType>
{
public:
  /**
   * 声明容器尺寸的类型。
   *
   */
  using size_type = unsigned int;

  /**
   * 矩形矩阵的构造函数，具有 @p n_rows 和 @p n_cols
   * 并使用网格 @p process_grid. 分布。参数 @p row_block_size 和 @p
   * column_block_size 是用于矩阵的块循环分布的块大小。
   * 一般来说，建议使用 $2$ 的幂，例如 $16,32,64, \dots$  。
   *
   */
  ScaLAPACKMatrix(
    const size_type                                           n_rows,
    const size_type                                           n_columns,
    const std::shared_ptr<const Utilities::MPI::ProcessGrid> &process_grid,
    const size_type               row_block_size    = 32,
    const size_type               column_block_size = 32,
    const LAPACKSupport::Property property = LAPACKSupport::Property::general);

  /**
   * 一个大小为 @p size, 的正方形矩阵的构造函数，使用 @p
   * process_grid. 中的过程网格进行分布。
   * 矩阵的行和列使用相同的块大小。
   * 一般来说，建议使用 $2$ 的幂，例如 $16,32,64, \dots$  。
   *
   */
  ScaLAPACKMatrix(
    const size_type                                           size,
    const std::shared_ptr<const Utilities::MPI::ProcessGrid> &process_grid,
    const size_type                                           block_size = 32,
    const LAPACKSupport::Property                             property =
      LAPACKSupport::Property::symmetric);

  /**
   * 一般矩形矩阵的构造函数，从文件 @p filename
   * 中读取并使用网格 @p process_grid. 分布。 使用HDF5从文件
   * @p filename 加载矩阵。
   * 如果在建立deal.II时没有使用HDF5，调用这个函数将导致一个异常。
   * 参数 @p row_block_size 和 @p column_block_size
   * 是用于矩阵的块循环分布的块大小。
   * 一般来说，建议使用 $2$ 的幂，例如 $16,32,64, \dots$  。
   *
   */
  ScaLAPACKMatrix(
    const std::string &                                       filename,
    const std::shared_ptr<const Utilities::MPI::ProcessGrid> &process_grid,
    const size_type row_block_size    = 32,
    const size_type column_block_size = 32);

  /**
   * 破坏器
   *
   */
  ~ScaLAPACKMatrix() override = default;

  /**
   * 用 @p n_rows 和 @p n_cols 初始化矩形矩阵，并使用网格 @p
   * process_grid. 分布。参数 @p row_block_size 和 @p column_block_size
   * 是用于矩阵的块循环分布的块大小。
   * 一般来说，建议使用 $2$ 的幂，例如 $16,32,64, \dots$  。
   *
   */
  void
  reinit(
    const size_type                                           n_rows,
    const size_type                                           n_columns,
    const std::shared_ptr<const Utilities::MPI::ProcessGrid> &process_grid,
    const size_type               row_block_size    = 32,
    const size_type               column_block_size = 32,
    const LAPACKSupport::Property property = LAPACKSupport::Property::general);

  /**
   * 初始化大小为 @p size 的正方形矩阵，并使用网格 @p
   * process_grid. 进行分布。参数 @p block_size
   * 用于矩阵的块循环分布。
   * 矩阵的行和列使用相同的块大小。
   * 一般来说，建议使用 $2$ 的幂，例如 $16,32,64, \dots$  。
   *
   */
  void
  reinit(const size_type                                           size,
         const std::shared_ptr<const Utilities::MPI::ProcessGrid> &process_grid,
         const size_type               block_size = 32,
         const LAPACKSupport::Property property =
           LAPACKSupport::Property::symmetric);

  /**
   * 将 @p property 分配给这个矩阵。
   *
   */
  void
  set_property(const LAPACKSupport::Property property);

  /**
   * 返回该矩阵的当前 @p property 。
   *
   */
  LAPACKSupport::Property
  get_property() const;

  /**
   * 返回此矩阵的当前 @p state 。
   *
   */
  LAPACKSupport::State
  get_state() const;

  /**
   * 来自普通FullMatrix的赋值操作。
   * @note
   * 这个函数应该只用于相对较小的矩阵尺寸。它主要是为了调试的目的。
   *
   */
  ScaLAPACKMatrix<NumberType> &
  operator=(const FullMatrix<NumberType> &);

  /**
   * 将本地拥有的 @p matrix 的内容复制到分布式矩阵中。
   * 分布式矩阵和进程 @p rank 上的 @p matrix
   * 必须有匹配的尺寸。
   * 对于所有进程来说，除了具有等级 @p rank 的进程，序列
   * @p matrix 不被引用。  用户必须确保所有进程都以相同的
   * @p rank. 来调用。 @p rank
   * 是指用于创建分布式矩阵的进程网格的MPI通信器的一个进程。
   *
   */
  void
  copy_from(const LAPACKFullMatrix<NumberType> &matrix,
            const unsigned int                  rank);

  /**
   * 将分布式矩阵的内容复制到 @p matrix. 中。
   * @note
   * 这个函数应该只用于相对较小的矩阵尺寸。它主要是为了调试的目的。
   *
   */
  void
  copy_to(FullMatrix<NumberType> &matrix) const;

  /**
   * 将分布式矩阵的内容复制到本地复制的 @p matrix
   * 上，该进程的等级为 @p rank.  对于除 @p rank
   * 外的所有进程， @p matrix 不被引用。  分布式矩阵和进程
   * @p rank 上的 @p matrix 必须有匹配的尺寸。
   * 用户必须确保所有进程都以相同的 @p rank. 来调用。 @p
   * rank
   * 是指用于创建分布式矩阵的进程网格的MPI通信器的一个进程。
   *
   */
  void
  copy_to(LAPACKFullMatrix<NumberType> &matrix, const unsigned int rank) const;

  /**
   * 将分布式矩阵的内容复制到不同的分布式矩阵中  @p dest.
   * 该函数也适用于具有不同进程网格或块循环分布的矩阵。
   *
   */
  void
  copy_to(ScaLAPACKMatrix<NumberType> &dest) const;

  /**
   * 将分布式矩阵A的一个子矩阵（子集）复制到分布式矩阵的一个子矩阵
   * @p B.  。
   *
   *
   *
   *
   *
   * - 子矩阵A的第一个元素的全局行和列索引由 @p offset_A 提供，行索引=  <code>offset_A.first</code>  ，列索引=  <code>offset_A.second</code>  。
   *
   *
   *
   *
   *
   * - 子矩阵B的第一个元素的全局行和列索引由 @p offset_B 提供，行索引=  <code>offset_B.first</code>  ，列索引=  <code>offset_B.second</code>  。
   *
   *
   *
   *
   *
   *
   * - 要复制的子矩阵的尺寸由 @p submatrix_size 给出，行数=  <code>submatrix_size.first</code>  ，列数=  <code>submatrix_size.second</code>  。      如果需要复制具有相同块状循环分布的完整矩阵，请使用  ScaLAPACKMatrix<NumberType>::copy_to(ScaLAPACKMatrix<NumberType>  &dest），只需一个参数，以避免通信。    矩阵 @p A 和 @p B 的底层进程网格必须是用同一个MPI通信器建立的。
   *
   */
  void
  copy_to(ScaLAPACKMatrix<NumberType> &                B,
          const std::pair<unsigned int, unsigned int> &offset_A,
          const std::pair<unsigned int, unsigned int> &offset_B,
          const std::pair<unsigned int, unsigned int> &submatrix_size) const;

  /**
   * 转置赋值。  $\mathbf{A} = \mathbf{B}^T$  矩阵 $\mathbf{A}$ 和
   * $\mathbf{B}$ 必须具有相同的进程网格。
   * 必须满足以下对齐条件。  $MB_A=NB_B$  和  $NB_A=MB_B$  .
   *
   */
  void
  copy_transposed(const ScaLAPACKMatrix<NumberType> &B);

  /**
   * 基于输入参数 @p transpose_B 和对齐条件的操作在下表中进行了总结。    | transpose_B | 块大小 | 操作 | :---------: | :--------------------------: | :-------------------------------------------: | false |  $MB_A=MB_B$   <br>   $NB_A=NB_B$  |  $\mathbf{A} = a \mathbf{A} + b \mathbf{B}$  | | true |  $MB_A=NB_B$   <br>   $NB_A=MB_B$  |  $\mathbf{A} = a \mathbf{A} + b \mathbf{B}^T$  | 矩阵 $\mathbf{A}$  和 $\mathbf{B}$  必须具有相同的过程网格。
   *
   */
  void
  add(const ScaLAPACKMatrix<NumberType> &B,
      const NumberType                   a           = 0.,
      const NumberType                   b           = 1.,
      const bool                         transpose_B = false);

  /**
   * 矩阵-加法。    $\mathbf{A} = \mathbf{A} + b\, \mathbf{B}$  矩阵
   * $\mathbf{A}$ 和 $\mathbf{B}$ 必须有相同的过程网格。
   * 必须满足以下对齐条件。  $MB_A=MB_B$  和  $NB_A=NB_B$  .
   *
   */
  void
  add(const NumberType b, const ScaLAPACKMatrix<NumberType> &B);

  /**
   * 矩阵-加法。    $\mathbf{A} = \mathbf{A} + b\, \mathbf{B}^T$  矩阵
   * $\mathbf{A}$ 和 $\mathbf{B}$ 必须具有相同的过程网格。
   * 必须满足以下对齐条件。  $MB_A=NB_B$  和  $NB_A=MB_B$  。
   *
   */
  void
  Tadd(const NumberType b, const ScaLAPACKMatrix<NumberType> &B);

  /**
   * 矩阵-矩阵-乘法。    基于输入参数和排列条件的操作总结在下表中。    | 转置_A | 转置_B | 块大小 | 操作 | :---------: | :---------: | :-------------------------------------------: | :-------------------------------------------------------------: | 虚假 | 虚假 |  $MB_A=MB_C$   <br>   $NB_A=MB_B$   <br>   $NB_B=NB_C$  |  $\mathbf{C} = b \mathbf{A} \cdot \mathbf{B} + c \mathbf{C}$  | 虚假 | 真实 |  $MB_A=MB_C$  ]  <br>   $NB_A=NB_B$   <br>   $MB_B=NB_C$  |  $\mathbf{C} = b \mathbf{A} \cdot \mathbf{B}^T + c \mathbf{C}$  | | 真实 | 错误 |  $MB_A=MB_B$   <br>   $NB_A=MB_C$   <br>  ]  $NB_B=NB_C$  |  $\mathbf{C} = b \mathbf{A}^T \cdot \mathbf{B} + c \mathbf{C}$  | | 真 | 真 |  $MB_A=NB_B$   <br>   $NB_A=MB_C$   <br>   $MB_B=NB_C$  |  $\mathbf{C} = b \mathbf{A}^T \cdot \mathbf{B}^T + c \mathbf{C}$  | 假设 $\mathbf{A}$ 和 $\mathbf{B}$ 的大小兼容，并且 $\mathbf{C}$ 已经具有正确大小。    矩阵 $\mathbf{A}$ ， $\mathbf{B}$ 和 $\mathbf{C}$ 必须具有相同的过程网格。
   *
   */
  void
  mult(const NumberType                   b,
       const ScaLAPACKMatrix<NumberType> &B,
       const NumberType                   c,
       ScaLAPACKMatrix<NumberType> &      C,
       const bool                         transpose_A = false,
       const bool                         transpose_B = false) const;

  /**
   * 矩阵-矩阵-乘法。    可选参数 @p adding
   * 决定了结果是存储在 $\mathbf{C}$ 中还是添加到 $\mathbf{C}$
   * 中。 如果（  @p adding)   $\mathbf{C} = \mathbf{C} + \mathbf{A}
   * \cdot \mathbf{B}$  否则 $\mathbf{C} = \mathbf{A} \cdot \mathbf{B}$
   * 假设 $\mathbf{A}$ 和 $\mathbf{B}$ 的大小兼容，并且
   * $\mathbf{C}$ 已经有正确大小。
   * 必须满足以下对齐条件。  $MB_A=MB_C$  ,  $NB_A=MB_B$  和
   * $NB_B=NB_C$  。
   *
   */
  void
  mmult(ScaLAPACKMatrix<NumberType> &      C,
        const ScaLAPACKMatrix<NumberType> &B,
        const bool                         adding = false) const;

  /**
   * 使用  $\mathbf{A}$  的转置进行矩阵-矩阵-乘法。
   * 可选参数 @p adding 决定了结果是存储在 $\mathbf{C}$
   * 中还是添加到 $\mathbf{C}$ 中。 如果（  @p adding)
   * $\mathbf{C} = \mathbf{C} + \mathbf{A}^T \cdot \mathbf{B}$  否则
   * $\mathbf{C} = \mathbf{A}^T \cdot \mathbf{B}$  假设 $\mathbf{A}$ 和
   * $\mathbf{B}$ 的大小兼容，并且 $\mathbf{C}$
   * 已经具有正确大小。    必须满足以下对齐条件。
   * $MB_A=MB_B$  ,  $NB_A=MB_C$  和  $NB_B=NB_C$  .
   *
   */
  void
  Tmmult(ScaLAPACKMatrix<NumberType> &      C,
         const ScaLAPACKMatrix<NumberType> &B,
         const bool                         adding = false) const;

  /**
   * 使用  $\mathbf{B}$  的转置进行矩阵-矩阵-乘法。
   * 可选参数 @p adding 决定了结果是存储在 $\mathbf{C}$
   * 中还是添加到 $\mathbf{C}$ 中。 如果（  @p adding)
   * $\mathbf{C} = \mathbf{C} + \mathbf{A} \cdot \mathbf{B}^T$  否则
   * $\mathbf{C} = \mathbf{A} \cdot \mathbf{B}^T$  假设 $\mathbf{A}$ 和
   * $\mathbf{B}$ 的大小兼容，并且 $\mathbf{C}$
   * 已经具有正确大小。    必须满足以下对齐条件。
   * $MB_A=MB_C$  ,  $NB_A=NB_B$  和  $MB_B=NB_C$  。
   *
   */
  void
  mTmult(ScaLAPACKMatrix<NumberType> &      C,
         const ScaLAPACKMatrix<NumberType> &B,
         const bool                         adding = false) const;

  /**
   * 使用  $\mathbf{A}$  和  $\mathbf{B}$
   * 的转置进行矩阵-矩阵-乘法。    可选参数 @p adding
   * 决定了结果是存储在 $\mathbf{C}$ 中还是添加到 $\mathbf{C}$
   * 中。 如果（  @p adding)   $\mathbf{C} = \mathbf{C} + \mathbf{A}^T
   * \cdot \mathbf{B}^T$  否则 $\mathbf{C} = \mathbf{A}^T \cdot
   * \mathbf{B}^T$  假设 $\mathbf{A}$ 和 $\mathbf{B}$
   * 的大小兼容，并且 $\mathbf{C}$ 已经有正确大小。
   * 必须满足以下对齐条件。  $MB_A=NB_B$  ,  $NB_A=MB_C$  和
   * $MB_B=NB_C$  。
   *
   */
  void
  TmTmult(ScaLAPACKMatrix<NumberType> &      C,
          const ScaLAPACKMatrix<NumberType> &B,
          const bool                         adding = false) const;

  /**
   * 使用HDF5将分布式矩阵存储在 @p filename 中。
   * 如果在建立deal.II时没有使用HDF5，调用这个函数将导致一个异常。
   * 如果HDF5是用MPI构建的，那么将使用并行I/O来保存矩阵。
   * 否则，只有一个进程会进行输出。这意味着在内部，分布式矩阵被复制到一个进程中，由该进程进行输出。因此，矩阵必须适合一个进程的内存。
   * 为了调整I/O性能，特别是平行I/O，用户可以定义可选的参数
   * @p chunk_size.  所有MPI进程都需要以相同的值调用该函数。
   * 矩阵是分块写入文件的，因此系统的属性决定了最佳的分块大小。在内部，HDF5将矩阵分割成<tt>chunk_size.first</tt>
   * x
   * <tt>chunk_size.second</tt>大小的块，其中<tt>chunk_size.first</tt>是一个块的行数，<tt>chunk_size.second</tt>是列的数量。
   *
   */
  void
  save(const std::string &                          filename,
       const std::pair<unsigned int, unsigned int> &chunk_size =
         std::make_pair(numbers::invalid_unsigned_int,
                        numbers::invalid_unsigned_int)) const;

  /**
   * 使用HDF5从文件 @p filename 加载分布式矩阵。
   * 如果在建立deal.II时没有使用HDF5，调用这个函数将导致一个异常。
   * 矩阵的尺寸必须与存储在文件中的矩阵相同。
   * 如果HDF5是用MPI构建的，那么将使用并行I/O来加载矩阵。
   * 否则，只有一个进程会从存储中加载矩阵，并随后将内容分发给其他进程。
   *
   */
  void
  load(const std::string &filename);

  /**
   * 使用ScaLAPACK函数计算矩阵的Cholesky因子化
   * <code>pXpotrf</code>
   * 。因式分解的结果被存储在这个对象中。
   *
   */
  void
  compute_cholesky_factorization();

  /**
   * 使用ScaLAPACK函数 <code>pXgetrf</code>
   * 计算矩阵的LU因式分解，并通过行间交换进行部分透视。
   * 因式分解的结果被存储在这个对象中。
   *
   */
  void
  compute_lu_factorization();

  /**
   * 通过首先计算对称矩阵的Cholesky或一般矩阵的LU因子化来反转矩阵，然后使用
   * <code>pXpotri</code>  或  <code>pXgetri</code>
   * 建立实际的反转。如果矩阵是三角形的，则跳过LU因子化步骤，直接使用
   * <code>pXtrtri</code> 。
   * 如果之前已经应用了Cholesky或LU因子化，则直接调用
   * <code>pXpotri</code> or <code>pXgetri</code> 。
   * 逆运算被存储在这个对象中。
   *
   */
  void
  invert();

  /**
   * 计算选定的特征值和可选的实数对称矩阵的特征向量
   * $\mathbf{A} \in \mathbb{R}^{M \times M}$  。
   * 特征值/特征向量是通过规定一个指数范围来选择的  @p
   * index_limits.  如果成功，计算出的特征值将按升序排列。
   * 特征向量被存储在矩阵的列中，从而覆盖了矩阵的原始内容。
   * 如果必须计算所有的特征值/特征向量，在 @p index_limits.
   * 中传递封闭区间 $ \left[ 0, M-1 \right] $ 如果需要 $r$
   * 最大的特征值/特征向量，则传递封闭区间 $ \left[ M-r, M-1
   * \right] $ 。
   *
   */
  std::vector<NumberType>
  eigenpairs_symmetric_by_index(
    const std::pair<unsigned int, unsigned int> &index_limits,
    const bool                                   compute_eigenvectors);

  /**
   * 计算选定的特征值，并可选择计算特征向量。
   * 通过规定特征值的取值范围 @p value_limits
   * 来选择特征值/特征向量。
   * 如果成功，计算出的特征值将按升序排列。
   * 特征向量被存储在矩阵的列中，从而覆盖了矩阵的原始内容。
   *
   */
  std::vector<NumberType>
  eigenpairs_symmetric_by_value(
    const std::pair<NumberType, NumberType> &value_limits,
    const bool                               compute_eigenvectors);

  /**
   * 使用MRRR算法计算选定的特征值和可选的实数对称矩阵
   * $\mathbf{A} \in \mathbb{R}^{M \times M}$ 的特征向量。
   * 通过规定指数范围选择特征值/特征向量  @p index_limits.
   * 如果成功，计算出的特征值将按升序排列。
   * 特征向量被存储在矩阵的列中，从而覆盖了矩阵的原始内容。
   * 如果必须计算所有的特征值/特征向量，在 @p index_limits.
   * 中传递封闭区间 $ \left[ 0, M-1 \right] $ 如果希望得到 $r$
   * 最大的特征值/特征向量，则传递封闭区间 $ \left[ M-r, M-1
   * \right] $ 。
   *
   */
  std::vector<NumberType>
  eigenpairs_symmetric_by_index_MRRR(
    const std::pair<unsigned int, unsigned int> &index_limits,
    const bool                                   compute_eigenvectors);

  /**
   * 使用MRRR算法计算选定的特征值和可选的实数对称矩阵
   * $\mathbf{A} \in \mathbb{R}^{M \times M}$ 的特征向量。
   * 特征值/特征向量是通过规定特征值的取值范围 @p
   * value_limits 来选择。
   * 如果成功，计算出的特征值将按升序排列。
   * 特征向量被存储在矩阵的列中，从而覆盖了矩阵的原始内容。
   *
   */
  std::vector<NumberType>
  eigenpairs_symmetric_by_value_MRRR(
    const std::pair<NumberType, NumberType> &value_limits,
    const bool                               compute_eigenvectors);

  /**
   * 计算矩阵 $\mathbf{A} \in \mathbb{R}^{M \times N}$
   * 的奇异值分解（SVD），可选择计算左和/或右奇异向量。SVD写成
   * $\mathbf{A} = \mathbf{U} \cdot \mathbf{\Sigma} \cdot \mathbf{V}^T$ ，
   * $\mathbf{\Sigma} \in \mathbb{R}^{M \times N}$ 为对角矩阵，
   * $\mathbf{U} \in \mathbb{R}^{M \times M}$ 和 $\mathbf{V} \in
   * \mathbb{R}^{M \times M}$ 为正交矩阵。 $\mathbf{\Sigma}$
   * 的对角线元素是 $A$ 的奇异值， $\mathbf{U}$ 和 $\mathbf{V}$
   * 的列分别是相应的左和右奇异向量。奇异值按递减顺序返回，只计算
   * $\min(M,N)$ 的第一列和 $\mathbf{V}^T$ 的行。
   * 返回时，矩阵的内容是不可用的。  矩阵 $\mathbf{A}$
   * 的行和列必须具有相同的块循环分布。
   * 如果需要左奇异向量，矩阵 $\mathbf{A}$ 和 $\mathbf{U}$
   * 必须以相同的过程网格和块循环分布构建。如果需要右奇异向量，矩阵
   * $\mathbf{A}$ 和 $\mathbf{V}^T$
   * 必须以相同的过程网格和块循环分布构建。
   * 为了避免计算左和/或右奇异向量，该函数接受
   * <code>nullptr</code> 为 @p U 和/或 @p VT. 。
   *
   */
  std::vector<NumberType>
  compute_SVD(ScaLAPACKMatrix<NumberType> *U  = nullptr,
              ScaLAPACKMatrix<NumberType> *VT = nullptr);

  /**
   * 解决涉及矩阵 $\mathbf{A} \in \mathbb{R}^{M \times N}$ 或其转置
   * $\mathbf{A}^T$ 的过定或欠定实数线性系统，使用 $\mathbf{A}$
   * 的QR或LQ因式分解，用于矩阵 $N_{\rm RHS}$ 列中的RHS向量
   * 假设 $\mathbf{A}$ 具有全等级。  $\rm{rank}(\mathbf{A}) =
   * \min(M,N)$  .     支持以下选项。
   *
   *
   *
   *
   *
   *
   * - If(!transpose) and  $M \geq N$  : 超定系统 $\min \Vert \mathbf{B}
   *
   * - \mathbf{A}\cdot \mathbf{X}\Vert$ 的最小平方解 .n 退出后， $0$
   * 至 $N-1$
   * 的行包含最小平方解向量。每一列的剩余平方和由该列中
   * $N$ 至 $M-1$ 元素的平方和给出。
   *
   *
   *
   *
   *
   *
   * - If(!transpose) and  $M < N$  : 查找欠定系统的最小规范解  $\mathbf{A} \cdot \mathbf{X} = \mathbf{B}$  .n 退出后， $\mathbf{B}$  的列包含最小规范解向量。
   *
   *
   *
   *
   *
   *
   * - 如果(转置)和 $M \geq N$ ：找到欠定系统 $ \mathbf{A}^\top \cdot \mathbf{X} = \mathbf{B}$ 的最小规范解 .n 退出后， $\mathbf{B}$ 的列包含最小规范解向量。
   *
   *
   *
   *
   *
   *
   * - 如果(转置)和  $M < N$  : 超定系统的最小平方解  $\min \Vert \mathbf{B}
   *
   * - \mathbf{A}^\top \cdot \mathbf{X}\Vert$  .n 退出后，行  $0$  到
   * $M-1$
   * 包含最小平方解向量。每一列的剩余平方和由该列中元素
   * $M$ 至 $N-1$ 的平方和给出。    If(!tranpose) then  $\mathbf{B}
   * \in \mathbb{R}^{M \times N_{\rm RHS}}$  , otherwise  $\mathbf{B} \in
   * \mathbb{R}^{N \times N_{\rm RHS}}$  。  矩阵 $\mathbf{A}$ 和
   * $\mathbf{B}$ 的行和列必须有相同的块循环分布。
   *
   */
  void
  least_squares(ScaLAPACKMatrix<NumberType> &B, const bool transpose = false);

  /**
   * 使用奇异值分解 $\mathbf{A} = \mathbf{U} \cdot \mathbf{\Sigma}
   * \cdot \mathbf{V}^T$ 计算实数矩阵 $\mathbf{A} \in \mathbb{R}^{M
   * \times N}$ 的假逆 $\mathbf{A}^+ \in \mathbb{R}^{N \times M}$
   * （Moore-Penrose逆）。    与逆不同，假逆  $\mathbf{A}^+ =
   * \mathbf{V} \cdot \mathbf{\Sigma}^+ \cdot \mathbf{U}^T$
   * 既适用于矩形矩阵，也适用于奇异矩阵  $\mathbf{A}$  。
   * 对于一个矩形 $\mathbf{\Sigma}$
   * ，伪逆的计算方法是取对角线上每个非零元素的倒数，保留零的位置，然后转置
   * $\mathbf{\Sigma}$  。在数值计算中，只有奇异值 $\sigma_i >
   * \sigma_{\text{max}} \, \text{ratio}$
   * 被考虑在内。成功退出后，该函数返回满足该条件的奇异值的数量。这个值可以解释为
   * $\mathbf{A}$  的等级。    返回时该对象包含伪逆
   * $\mathbf{A}^+ \in \mathbb{R}^{N \times M}$  。
   * 必须满足以下排列条件。  $MB_A = NB_A$  .
   *
   */
  unsigned int
  pseudoinverse(const NumberType ratio);

  /**
   * 估计一个SPD矩阵在 $l_1$ -norm中的条件数。
   * 该矩阵必须处于Cholesky状态（见compute_cholesky_factorization()）。条件数的倒数被返回，以避免在条件数非常大时出现溢出的可能性。
   * @p a_norm 在调用Cholesky因式分解之前，必须包含矩阵的
   * $l_1$  -norm（见l1_norm()）。
   * @note  另一种方法是明确计算矩阵的逆值，并手动构建
   * $k_1 = ||\mathbf{A}||_1 \, ||\mathbf{A}^{-1}||_1$  。
   *
   */
  NumberType
  reciprocal_condition_number(const NumberType a_norm) const;

  /**
   * 计算矩阵的 $l_1$ -norm。
   *
   */
  NumberType
  l1_norm() const;

  /**
   * 计算矩阵的 $l_{\infty}$ 准则。
   *
   */
  NumberType
  linfty_norm() const;

  /**
   * 计算矩阵的Frobenius准则。
   *
   */
  NumberType
  frobenius_norm() const;

  /**
   * $M \times N$ 矩阵的行数。
   *
   */
  size_type
  m() const;

  /**
   * $M \times N$ 矩阵的列数。
   *
   */
  size_type
  n() const;

  /**
   * 这个MPI进程上的本地行数。
   *
   */
  unsigned int
  local_m() const;

  /**
   * 这个MPI进程上的本地列数。
   *
   */
  unsigned int
  local_n() const;

  /**
   * 返回给定本地行的全局行号  @p loc_row  。
   *
   */
  unsigned int
  global_row(const unsigned int loc_row) const;

  /**
   * 为给定的本地列返回全局列号  @p loc_column.  。
   *
   */
  unsigned int
  global_column(const unsigned int loc_column) const;

  /**
   * 读取对本地元素的访问。
   *
   */
  NumberType
  local_el(const unsigned int loc_row, const unsigned int loc_column) const;

  /**
   * 对本地元素的写访问。
   *
   */
  NumberType &
  local_el(const unsigned int loc_row, const unsigned int loc_column);

  /**
   * 通过数组 @p factors.
   * 中提供的标量对分布式矩阵的列进行缩放。数组 @p
   * factors 必须有与矩阵列同样多的条目。     @p factors
   * 的副本必须在底层MPI通信器的所有进程中可用。
   * @note   @p InputVector
   * 的基本前提是必须能够从中创建ArrayView； @p std::vector
   * 和Vector类满足这一要求。
   *
   */
  template <class InputVector>
  void
  scale_columns(const InputVector &factors);

  /**
   * 通过数组 @p factors.
   * 中提供的标量来缩放分布式矩阵的行，数组 @p factors
   * 的条目数必须与矩阵的行一样多。     @p factors
   * 的副本必须在底层MPI通信器的所有进程中可用。
   * @note   @p InputVector
   * 的基本前提是必须能够从中创建ArrayView； @p std::vector
   * 和Vector类满足这一要求。
   *
   */
  template <class InputVector>
  void
  scale_rows(const InputVector &factors);

private:
  /**
   * 使用ScaLAPACK的内部函数计算一个分布式对称密集矩阵的规范。
   *
   */
  NumberType
  norm_symmetric(const char type) const;

  /**
   * 使用ScaLAPACK的内部函数计算一个分布式密集矩阵的规范。
   *
   */
  NumberType
  norm_general(const char type) const;

  /**
   * 计算选定的特征值和可选的特征向量。
   * 特征值/特征向量是通过为特征值规定一个指数范围 @p
   * index_limits 或一个值范围 @p value_limits
   * 来选择的。如果同时规定了两个范围（意味着两个范围都与默认值不同），该函数将抛出一个异常，因为这种模糊性是被禁止的。如果成功，计算出的特征值将按升序排列。特征向量被存储在矩阵的列中，从而覆盖了矩阵的原始内容。
   *
   */
  std::vector<NumberType>
  eigenpairs_symmetric(
    const bool                                   compute_eigenvectors,
    const std::pair<unsigned int, unsigned int> &index_limits =
      std::make_pair(numbers::invalid_unsigned_int,
                     numbers::invalid_unsigned_int),
    const std::pair<NumberType, NumberType> &value_limits =
      std::make_pair(std::numeric_limits<NumberType>::quiet_NaN(),
                     std::numeric_limits<NumberType>::quiet_NaN()));

  /**
   * 使用MRRR算法计算选定的特征值和可选的实数对称矩阵
   * $\mathbf{A} \in \mathbb{R}^{M \times M}$ 的特征向量。
   * 特征值/特征向量是通过规定指数范围 @p index_limits
   * 或特征值的取值范围 @p value_limits
   * 来选择。如果两个范围都被规定了（意味着两个范围都与默认值不同），该函数将抛出一个异常，因为这种模糊性是被禁止的。
   * 通过调用这个函数，矩阵的原始内容将被覆盖。如果要求，特征向量将被存储在矩阵的列中。同样，如果只需要特征值，矩阵的内容也将被覆盖。
   * 如果成功，计算出的特征值将按升序排列。
   * @note
   * 由于Netlib-ScaLAPACK中的一个错误，可以计算所有的或没有特征向量。
   * 因此，必须对输入 @p index_limits
   * 进行相应设置。使用Intel-MKL则不需要这种限制。
   *
   */
  std::vector<NumberType>
  eigenpairs_symmetric_MRRR(
    const bool                                   compute_eigenvectors,
    const std::pair<unsigned int, unsigned int> &index_limits =
      std::make_pair(numbers::invalid_unsigned_int,
                     numbers::invalid_unsigned_int),
    const std::pair<NumberType, NumberType> &value_limits =
      std::make_pair(std::numeric_limits<NumberType>::quiet_NaN(),
                     std::numeric_limits<NumberType>::quiet_NaN()));

  /* 使用串行例程在 @p filename 中存储分布矩阵。 
*
*/
  void
  save_serial(const std::string &                          filename,
              const std::pair<unsigned int, unsigned int> &chunk_size) const;

  /* 使用串行例程从文件 @p filename 加载分布式矩阵。 
*
*/
  void
  load_serial(const std::string &filename);

  /* 使用并行例程将分布式矩阵存储在 @p filename 中。 
*
*/
  void
  save_parallel(const std::string &                          filename,
                const std::pair<unsigned int, unsigned int> &chunk_size) const;

  /* 使用并行例程从文件 @p filename 中加载分布式矩阵。 
*
*/
  void
  load_parallel(const std::string &filename);

  /**
   * 由于ScaLAPACK操作会改变矩阵项的含义，所以我们在这里记录最后一次操作后的当前状态。
   *
   */
  LAPACKSupport::State state;

  /**
   * 矩阵的其他属性，可能有助于选择更有效的ScaLAPACK函数。
   *
   */
  LAPACKSupport::Property property;

  /**
   * 一个指向 Utilities::MPI::ProcessGrid
   * 对象的共享指针，它包含一个BLACS上下文和一个MPI通信器，以及其他必要的数据结构。
   *
   */
  std::shared_ptr<const Utilities::MPI::ProcessGrid> grid;

  /**
   * 矩阵中的行数。
   *
   */
  int n_rows;

  /**
   * 矩阵中的列数。
   *
   */
  int n_columns;

  /**
   * 行块大小。
   *
   */
  int row_block_size;

  /**
   * 列块大小。
   *
   */
  int column_block_size;

  /**
   * 当前进程拥有的矩阵中的行数。
   *
   */
  int n_local_rows;

  /**
   * 当前进程拥有的矩阵中的列数。
   *
   */
  int n_local_columns;

  /**
   * ScaLAPACK描述向量。
   *
   */
  int descriptor[9];

  /**
   * 工作区阵列。
   *
   */
  mutable std::vector<NumberType> work;

  /**
   * 整数工作区数组。
   *
   */
  mutable std::vector<int> iwork;

  /**
   * 保存ScaLAPACK的矩阵分解程序所需的枢轴信息的整数阵列。
   *
   */
  std::vector<int> ipiv;

  /**
   * 一个字符，用于定义元素的存储位置，以防ScaLAPACK操作支持这个。
   *
   */
  const char uplo;

  /**
   * 过程网格的过程行，全局矩阵的第一行分布于此。
   *
   */
  const int first_process_row;

  /**
   * 全局矩阵第一列所分布的过程网格的过程列。
   *
   */
  const int first_process_column;

  /**
   * 全局行索引，决定在哪里开始一个子矩阵。
   * 目前这等于unity，因为我们不使用子矩阵。
   *
   */
  const int submatrix_row;

  /**
   * 全局列索引，决定子矩阵的起始位置。
   * 目前这等于unity，因为我们不使用子矩阵。
   *
   */
  const int submatrix_column;

  /**
   * 线程互斥。
   *
   */
  mutable Threads::Mutex mutex;
};

// ----------------------- inline functions ----------------------------

#  ifndef DOXYGEN

template <typename NumberType>
inline NumberType
ScaLAPACKMatrix<NumberType>::local_el(const unsigned int loc_row,
                                      const unsigned int loc_column) const
{
  return (*this)(loc_row, loc_column);
}



template <typename NumberType>
inline NumberType &
ScaLAPACKMatrix<NumberType>::local_el(const unsigned int loc_row,
                                      const unsigned int loc_column)
{
  return (*this)(loc_row, loc_column);
}


template <typename NumberType>
inline unsigned int
ScaLAPACKMatrix<NumberType>::m() const
{
  return n_rows;
}



template <typename NumberType>
inline unsigned int
ScaLAPACKMatrix<NumberType>::n() const
{
  return n_columns;
}



template <typename NumberType>
unsigned int
ScaLAPACKMatrix<NumberType>::local_m() const
{
  return n_local_rows;
}



template <typename NumberType>
unsigned int
ScaLAPACKMatrix<NumberType>::local_n() const
{
  return n_local_columns;
}


#  endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_SCALAPACK

#endif


