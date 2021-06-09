//include/deal.II-translator/lac/cuda_sparse_matrix_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2020 by the deal.II authors
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

#ifndef dealii_cuda_sparse_matrix_h
#define dealii_cuda_sparse_matrix_h

#include <deal.II/base/config.h>

#include <deal.II/base/subscriptor.h>

#include <iomanip>

#ifdef DEAL_II_COMPILER_CUDA_AWARE
#  include <deal.II/base/cuda.h>

#  include <deal.II/lac/cuda_vector.h>
#  include <deal.II/lac/sparse_matrix.h>

#  include <cusparse.h>

DEAL_II_NAMESPACE_OPEN

namespace CUDAWrappers
{
  /**
   * 这个类是cuSPARSE
   * csr稀疏矩阵的一个包装器。与deal.II自己的SparseMatrix不同，每一行中的所有元素都是按照列索引递增的顺序存储。
   * @note  这个模板的实例化是为<tt>  @<float@>  和
   * @<double@></tt>.  提供的。
   * @ingroup Matrix1
   *
   */
  template <typename Number>
  class SparseMatrix : public virtual Subscriptor
  {
  public:
    /**
     * 声明容器尺寸的类型。
     *
     */
    using size_type = int;

    /**
     * 矩阵条目的类型。
     *
     */
    using value_type = Number;

    /**
     * 声明一个保存实值数的类型，其精度与该类的模板参数相同。
     *
     */
    using real_type = Number;

    /**
     * @name  构造函数和初始化
     *
     */
    //@{
    /**
     * 构造函数。将矩阵初始化为空，没有任何结构，也就是说，矩阵根本无法使用。因此，这个构造函数只对属于一个类的成员的矩阵有用。
     * 在使用之前，你必须先用reinit初始化矩阵。
     *
     */
    SparseMatrix();

    /**
     * 构造函数。接受一个 Utilities::CUDA::Handle
     * 和一个主机上的稀疏矩阵。主机上的稀疏矩阵被复制到设备上，并且根据cuSPARSE支持的格式对元素进行重新排序。
     *
     */
    SparseMatrix(Utilities::CUDA::Handle &             handle,
                 const ::dealii::SparseMatrix<Number> &sparse_matrix_host);

    /**
     * 移动构造函数。通过窃取内部数据来创建一个新的SparseMatrix。
     *
     */
    SparseMatrix(CUDAWrappers::SparseMatrix<Number> &&);

    /**
     * 复制构造函数被删除。
     *
     */
    SparseMatrix(const CUDAWrappers::SparseMatrix<Number> &) = delete;

    /**
     * 解除构造函数。释放所有的内存。
     *
     */
    ~SparseMatrix();

    /**
     * 移动赋值运算符。
     *
     */
    SparseMatrix &
    operator=(CUDAWrappers::SparseMatrix<Number> &&);

    /**
     * 复制赋值被删除。
     *
     */
    SparseMatrix &
    operator=(const CUDAWrappers::SparseMatrix<Number> &) = delete;

    /**
     * 重新初始化稀疏矩阵。主机上的稀疏矩阵被复制到设备上，并且根据cuSPARSE支持的格式对元素进行重新排序。
     *
     */
    void
    reinit(Utilities::CUDA::Handle &             handle,
           const ::dealii::SparseMatrix<Number> &sparse_matrix_host);
    //@}

    /**
     * @name  矩阵的信息
     *
     */
    //@{
    /**
     * 返回共域（或范围）空间的维数。注意，矩阵的维度是
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
    std::size_t
    n_nonzero_elements() const;

    /**
     * 打印矩阵到给定的流，使用格式<tt>(row,column)
     * value</tt>，即每行打印矩阵的一个非零条目。如果<tt>across</tt>为真，则在单行上打印所有条目，使用格式row,column:value。
     * 如果参数<tt>diagonal_first</tt>为真，二次元矩阵的对角线元素将在其行中首先打印。如果它是假的，一行中的元素将按升列顺序写入。
     *
     */
    template <class StreamType>
    void
    print(StreamType &out,
          const bool  across         = false,
          const bool  diagonal_first = true) const;

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
    //@}

    /**
     * @name  修改条目
     *
     */
    //@{
    /**
     * 将整个矩阵乘以一个固定系数。
     *
     */
    SparseMatrix &
    operator*=(const Number factor);

    /**
     * 用整个矩阵除以一个固定的因子。
     *
     */
    SparseMatrix &
    operator/=(const Number factor);
    //@}

    /**
     * @name  乘法运算
     *
     */
    //@{
    /**
     * 矩阵-向量乘法：让 $dst = M \cdot src$ 与 $M$ 为该矩阵。
     *
     */
    void
    vmult(LinearAlgebra::CUDAWrappers::Vector<Number> &      dst,
          const LinearAlgebra::CUDAWrappers::Vector<Number> &src) const;

    /**
     * 矩阵-向量乘法：让 $dst = M^T \cdot src$ 与 $M$
     * 为这个矩阵。这个函数与vmult()的作用相同，但需要这个转置的矩阵。
     *
     */
    void
    Tvmult(LinearAlgebra::CUDAWrappers::Vector<Number> &      dst,
           const LinearAlgebra::CUDAWrappers::Vector<Number> &src) const;

    /**
     * 添加矩阵-向量的乘法。在 $dst$ 上添加 $M \cdot src$ ，
     * $M$ 是这个矩阵。
     *
     */
    void
    vmult_add(LinearAlgebra::CUDAWrappers::Vector<Number> &      dst,
              const LinearAlgebra::CUDAWrappers::Vector<Number> &src) const;

    /**
     * 添加矩阵-向量的乘法。将 $M^T \cdot src$ 加到 $dst$ ， $M$
     * 是这个矩阵。这个函数的作用与vmult_add()相同，但需要转置的矩阵。
     *
     */
    void
    Tvmult_add(LinearAlgebra::CUDAWrappers::Vector<Number> &      dst,
               const LinearAlgebra::CUDAWrappers::Vector<Number> &src) const;

    /**
     * 返回向量 $v$ 相对于该矩阵诱导的准则的平方，即
     * $\left(v,Mv\right)$
     * 。这很有用，例如，在有限背景下，一个函数的 $L_2$
     * 规范等于相对于代表有限元函数节点值的向量的质量矩阵的矩阵规范。
     * 很明显，对于这个操作，矩阵需要是二次的。
     *
     */
    Number
    matrix_norm_square(
      const LinearAlgebra::CUDAWrappers::Vector<Number> &v) const;

    /**
     * 计算矩阵标量乘积  $\left(u,Mv\right)$  。
     *
     */
    Number
    matrix_scalar_product(
      const LinearAlgebra::CUDAWrappers::Vector<Number> &u,
      const LinearAlgebra::CUDAWrappers::Vector<Number> &v) const;

    /**
     * 计算方程 $M \cdot x=b$ 的残差，其中残差被定义为 $r=b-M
     * \cdot x$ 。将残差写进  $dst$  。返回残差向量的 $l_2$
     * 准则。        源 $x$ 和目的 $dst$ 不能是同一个向量。
     *
     */
    Number
    residual(LinearAlgebra::CUDAWrappers::Vector<Number> &      dst,
             const LinearAlgebra::CUDAWrappers::Vector<Number> &x,
             const LinearAlgebra::CUDAWrappers::Vector<Number> &b) const;
    //@}

    /**
     * @name  矩阵规范
     *
     */
    //@{
    /**
     * 返回矩阵的 $l_1$ -规范，即 $|M|_1=\max_{\mathrm{all\ columns\
     * }j}\sum_{\mathrm{all\ rows\ }i} |M_{ij}|$
     * ，（最大列数之和）。这是自然的矩阵准则，与向量的
     * $l_1$ 准则兼容，即 $|Mv|_1\leq |M|_1 |v|_1$  。
     *
     */
    Number
    l1_norm() const;

    /**
     * 返回矩阵的 $l_\infty$ 准则，即 $|M|_\infty=\max_{\mathrm{all\
     * rows\ }i}\sum_{\mathrm{all\ columns\ }j} |M_{ij}|$
     * ，（行的最大和）。这是与向量的 $l_\infty$
     * 规范兼容的自然规范，即 $|Mv|_\infty \leq |M|_\infty
     * |v|_\infty$  。
     *
     */
    Number
    linfty_norm() const;

    /**
     * 返回矩阵的frobenius规范，即矩阵中所有条目的平方和的平方根。
     *
     */
    Number
    frobenius_norm() const;
    //@}

    /**
     * @name  访问底层CUDA数据
     *
     */
    //@{
    /**
     * 返回一个元组，包含矩阵值的指针，列索引的指针，行指针的指针，cuSPARSE矩阵描述，以及cuSPARSE
     * SP矩阵描述。
     *
     */
    std::tuple<Number *, int *, int *, cusparseMatDescr_t, cusparseSpMatDescr_t>
    get_cusparse_matrix() const;
    //@}

  private:
    /**
     * 用于调用cuSPARSE函数的cuSPARSE句柄。
     *
     */
    cusparseHandle_t cusparse_handle;

    /**
     * 稀疏矩阵中非零元素的数量。
     *
     */
    int nnz;

    /**
     * 稀疏矩阵的行数。
     *
     */
    int n_rows;

    /**
     * 稀疏矩阵的列数。
     *
     */
    int n_cols;

    /**
     * 指向稀疏矩阵的值（在设备上）的指针。
     *
     */
    std::unique_ptr<Number[], void (*)(Number *)> val_dev;

    /**
     * 指向稀疏矩阵的列索引（在设备上）的指针。
     *
     */
    std::unique_ptr<int[], void (*)(int *)> column_index_dev;

    /**
     * 指向稀疏矩阵的行指针的指针（在设备上）。
     *
     */
    std::unique_ptr<int[], void (*)(int *)> row_ptr_dev;

    /**
     * 矩阵的cuSPARSE描述。
     *
     */
    cusparseMatDescr_t descr;

    /**
     * cuSPARSE描述的稀疏矩阵。
     *
     */
    cusparseSpMatDescr_t sp_descr;
  };



  template <typename Number>
  inline SparseMatrix<Number>::size_type
  SparseMatrix<Number>::m() const
  {
    return n_rows;
  }



  template <typename Number>
  inline SparseMatrix<Number>::size_type
  SparseMatrix<Number>::n() const
  {
    return n_cols;
  }



  template <typename Number>
  inline std::size_t
  SparseMatrix<Number>::n_nonzero_elements() const
  {
    return nnz;
  }



  template <typename Number>
  template <class StreamType>
  inline void
  SparseMatrix<Number>::print(StreamType &out,
                              const bool  across,
                              const bool  diagonal_first) const
  {
    Assert(column_index_dev != nullptr, ExcNotInitialized());
    Assert(val_dev != nullptr, ExcNotInitialized());
    Assert(row_ptr_dev != nullptr, ExcNotInitialized());

    std::vector<int>    rows(n_rows + 1);
    std::vector<int>    cols(nnz);
    std::vector<double> val(nnz);
    Utilities::CUDA::copy_to_host(row_ptr_dev.get(), rows);
    Utilities::CUDA::copy_to_host(column_index_dev.get(), cols);
    Utilities::CUDA::copy_to_host(val_dev.get(), val);

    bool   has_diagonal = false;
    Number diagonal     = Number();

    for (size_type i = 0; i < n_rows; ++i)
      {
        if (diagonal_first)
          {
            // find the diagonal and print if it exists
            for (size_type j = rows[i]; j < rows[i + 1] && cols[j] <= i; ++j)
              {
                if (i == cols[j])
                  {
                    diagonal     = val[j];
                    has_diagonal = true;
                    if (across)
                      out << ' ' << i << ',' << i << ':' << diagonal;
                    else
                      out << '(' << i << ',' << i << ") " << diagonal
                          << std::endl;
                    break;
                  }
              }
          }
        for (size_type j = rows[i]; j < rows[i + 1]; ++j)
          {
            if (has_diagonal && i == cols[j])
              continue;
            if (across)
              out << ' ' << i << ',' << cols[j] << ':' << val[j];
            else
              out << "(" << i << "," << cols[j] << ") " << val[j] << std::endl;
          }
      }
    if (across)
      out << std::endl;
  }



  template <typename Number>
  void
  SparseMatrix<Number>::print_formatted(std::ostream &     out,
                                        const unsigned int precision,
                                        const bool         scientific,
                                        const unsigned int width_,
                                        const char *       zero_string,
                                        const double       denominator) const
  {
    Assert(column_index_dev != nullptr, ExcNotInitialized());
    Assert(val_dev != nullptr, ExcNotInitialized());
    Assert(row_ptr_dev != nullptr, ExcNotInitialized());

    std::vector<int>    rows(n_rows + 1);
    std::vector<int>    cols(nnz);
    std::vector<Number> val(nnz);
    Utilities::CUDA::copy_to_host(row_ptr_dev.get(), rows);
    Utilities::CUDA::copy_to_host(column_index_dev.get(), cols);
    Utilities::CUDA::copy_to_host(val_dev.get(), val);

    unsigned int width = width_;

    std::ios::fmtflags old_flags     = out.flags();
    unsigned int       old_precision = out.precision(precision);

    if (scientific)
      {
        out.setf(std::ios::scientific, std::ios::floatfield);
        if (!width)
          width = precision + 7;
      }
    else
      {
        out.setf(std::ios::fixed, std::ios::floatfield);
        if (!width)
          width = precision + 2;
      }

    for (size_type i = 0; i < n_rows; ++i)
      {
        size_type j = rows[i];
        for (size_type k = 0; k < n_cols; ++k)
          {
            if (k == cols[j])
              {
                out << std::setw(width) << val[j] * Number(denominator) << ' ';
                ++j;
              }
            else
              out << std::setw(width) << zero_string << ' ';
          }
        out << std::endl;
      };
    AssertThrow(out, ExcIO());

    // reset output format
    out.precision(old_precision);
    out.flags(old_flags);
  }
} // namespace CUDAWrappers

DEAL_II_NAMESPACE_CLOSE

#endif
#endif


