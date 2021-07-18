//include/deal.II-translator/lac/cuda_precondition_0.txt
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

#ifndef dealii_cuda_precondition_h
#define dealii_cuda_precondition_h

#include <deal.II/base/config.h>

#include <deal.II/base/cuda.h>
#include <deal.II/base/smartpointer.h>

#include <memory>

#ifdef DEAL_II_COMPILER_CUDA_AWARE

DEAL_II_NAMESPACE_OPEN

// forward-definition
#  ifndef DOXYGEN
namespace LinearAlgebra
{
  namespace CUDAWrappers
  {
    template <typename Number>
    class Vector;
  }
} // namespace LinearAlgebra
#  endif

namespace CUDAWrappers
{
  // forward definition
  template <typename Number>
  class SparseMatrix;

  /**
   * 该类实现了 @em 对称 CUDAWrappers::SparseMatrix
   * 矩阵的不完全Cholesky因子化（IC）预处理。
   * 该实现与cuSPARSE文档（https://docs.nvidia.com/cuda/cusparse/index.html#cusparse-lt-t-gt-csric02）中记载的实现非常接近。
   * @note  本模板的实例化提供给<tt>  @<float@>  和
   * @<double@></tt>.  。
   * @ingroup Preconditioners CUDAWrappers
   *
   */
  template <typename Number>
  class PreconditionIC
  {
  public:
    /**
     * 声明容器尺寸的类型。
     *
     */
    using size_type = int;

    /**
     * 标准化的数据结构，用于向预处理程序输送附加标志。
     *
     */
    struct AdditionalData
    {
      /**
       * 构造函数。 cuSPARSE允许计算和使用级别信息。
       * 根据文档，这可能会提高性能。
       * 建议尝试这两个选项。
       *
       */
      AdditionalData(bool use_level_analysis = true);

      /**
       * 决定在创建和应用预处理程序时是否使用电平信息的标志。更多信息请参见cusparseSolvePolicy_t的文档，网址为https://docs.nvidia.com/cuda/cusparse/index.html#cusparsesolvepolicy_t。
       *
       */
      bool use_level_analysis;
    };

    /**
     * 构造函数。
     *
     */
    PreconditionIC(const Utilities::CUDA::Handle &handle);

    /**
     * 拷贝构造函数被删除。
     *
     */
    PreconditionIC(const PreconditionIC<Number> &) = delete;

    /**
     * 拷贝赋值运算符被删除。
     *
     */
    PreconditionIC &
    operator=(const PreconditionIC<Number> &) = delete;

    /**
     * 解构器。释放所有在这个类中被初始化的资源。
     *
     */
    ~PreconditionIC();

    /**
     * 初始化这个对象。特别是，给定的矩阵被复制以就地修改。对于底层的稀疏模式指针被存储。具体来说，这意味着只要
     * @p matrix
     * 是有效的，并且在调用此函数后没有被改变，当前对象就能被可靠地使用。
     * @p additional_data 决定了是否使用级别信息。
     *
     */
    void
    initialize(const SparseMatrix<Number> &matrix,
               const AdditionalData &      additional_data = AdditionalData());

    /**
     * 应用预处理程序。
     *
     */
    void
    vmult(LinearAlgebra::CUDAWrappers::Vector<Number> &      dst,
          const LinearAlgebra::CUDAWrappers::Vector<Number> &src) const;

    /**
     * 应用预处理程序。由于预处理器是对称的，这与vmult()相同。
     *
     */
    void
    Tvmult(LinearAlgebra::CUDAWrappers::Vector<Number> &      dst,
           const LinearAlgebra::CUDAWrappers::Vector<Number> &src) const;

    /**
     * 返回共域（或范围）空间的维度。注意，矩阵是正方形的，其维度为
     * $m \times m$  。
     * @note
     * 只有在预处理程序已经初始化的情况下才可以调用这个函数。
     *
     */
    size_type
    m() const;

    /**
     * 返回共域（或范围）空间的维数。注意，矩阵是正方形的，其维度为
     * $n \times n$  。
     * @note
     * 只有在预处理程序已被初始化的情况下才可调用此函数。
     *
     */
    size_type
    n() const;

  private:
    /**
     * 用于调用cuSPARSE函数的cuSPARSE句柄。
     *
     */
    cusparseHandle_t cusparse_handle;

    /**
     * cuSPARSE对稀疏矩阵的描述  $M=LL^T$  。
     *
     */
    cusparseMatDescr_t descr_M;

    /**
     * cuSPARSE描述的下三角矩阵  $L$  .
     *
     */
    cusparseMatDescr_t descr_L;

    /**
     * 解决和分析结构为  $M=LL^T$  .
     *
     */
    csric02Info_t info_M;

    /**
     * 下三角矩阵的求解和分析结构  $L$  .
     *
     */
    csrsv2Info_t info_L;

    /**
     * 上三角矩阵的求解和分析结构  $L^T$  .
     *
     */
    csrsv2Info_t info_Lt;

    /**
     * 指向此对象初始化的矩阵的指针。
     *
     */
    SmartPointer<const SparseMatrix<Number>> matrix_pointer;

    /**
     * 指向计算出的预处理矩阵的值（在设备上）的指针。
     *
     */
    std::unique_ptr<Number[], void (*)(Number *)> P_val_dev;

    /**
     * 指向此对象初始化的稀疏矩阵的行指针的指针（在设备上）。由matrix_pointer守护。
     *
     */
    const int *P_row_ptr_dev;

    /**
     * 指向此对象初始化的稀疏矩阵的列索引（在设备上）的指针。由
     * matrix_pointer 保护。
     *
     */
    const int *P_column_index_dev;

    /**
     * 指向vmult()中使用的临时(辅助)向量的值(在设备上)的指针。
     *
     */
    std::unique_ptr<Number[], void (*)(Number *)> tmp_dev;

    /**
     * 指向一个内部缓冲区的指针（在设备上），用于计算分解。
     *
     */
    std::unique_ptr<void, void (*)(void *)> buffer_dev;

    /**
     * 确定是否应该为下三角矩阵生成级别信息  $L$
     * 。这个值可以通过一个AdditionalData对象来修改。
     *
     */
    cusparseSolvePolicy_t policy_L;

    /**
     * 确定是否应该为上三角矩阵生成水平信息  $L^T$
     * 。这个值可以通过一个AdditionalData对象来修改。
     *
     */
    cusparseSolvePolicy_t policy_Lt;

    /**
     * 确定是否应该为  $M=LL^T$
     * 生成电平信息。这个值可以通过一个AdditionalData对象来修改。
     *
     */
    cusparseSolvePolicy_t policy_M;

    /**
     * 行数与此对象初始化的矩阵的行数相同。
     *
     */
    int n_rows;

    /**
     * 非零元素的数量与此对象初始化的矩阵相同。
     *
     */
    int n_nonzero_elements;
  };

  /**
   * 该类为 CUDAWrappers::SparseMatrix
   * 矩阵实现了一个不完整的LU因子化预处理程序。
   * 该实现与cuSPARSE文档（https://docs.nvidia.com/cuda/cusparse/index.html#cusparse-lt-t-gt-csrilu02）中记载的实现非常接近。
   * @note  这个模板的实例化提供给<tt>  @<float@>  和
   * @<double@></tt>.  。
   * @ingroup Preconditioners CUDAWrappers
   *
   */
  template <typename Number>
  class PreconditionILU
  {
  public:
    /**
     * 声明容器大小的类型。
     *
     */
    using size_type = int;

    /**
     * 标准化的数据结构，用于向预处理程序输送额外的标志。
     *
     */
    struct AdditionalData
    {
      /**
       * cuSPARSE允许计算和使用级别信息。根据文档，这可能会提高性能。
       * 建议尝试这两个选项。
       *
       */
      AdditionalData(bool use_level_analysis = true);

      /**
       * 决定在创建和应用预处理程序时是否使用电平信息的标志。更多信息请参见cusparseSolvePolicy_t的文档，网址为https://docs.nvidia.com/cuda/cusparse/index.html#cusparsesolvepolicy_t。
       *
       */
      bool use_level_analysis;
    };

    /**
     * 构造函数。
     *
     */
    PreconditionILU(const Utilities::CUDA::Handle &handle);

    /**
     * 拷贝构造函数被删除。
     *
     */
    PreconditionILU(const PreconditionILU<Number> &) = delete;

    /**
     * 拷贝赋值运算符被删除。
     *
     */
    PreconditionILU &
    operator=(const PreconditionILU<Number> &) = delete;

    /**
     * 解构器。释放所有在这个类中被初始化的资源。
     *
     */
    ~PreconditionILU();

    /**
     * 初始化这个对象。特别是，给定的矩阵被复制以就地修改。对于底层的稀疏模式指针被存储。具体来说，这意味着只要
     * @p matrix
     * 是有效的，并且在调用此函数后没有被改变，当前对象就能被可靠地使用。
     * @p additional_data 决定了是否使用级别信息。
     *
     */
    void
    initialize(const SparseMatrix<Number> &matrix,
               const AdditionalData &      additional_data = AdditionalData());

    /**
     * 应用预处理程序。
     *
     */
    void
    vmult(LinearAlgebra::CUDAWrappers::Vector<Number> &      dst,
          const LinearAlgebra::CUDAWrappers::Vector<Number> &src) const;

    /**
     * 应用转置的预处理程序。尚未实施。
     *
     */
    void
    Tvmult(LinearAlgebra::CUDAWrappers::Vector<Number> &      dst,
           const LinearAlgebra::CUDAWrappers::Vector<Number> &src) const;

    /**
     * 返回码域（或范围）空间的维数。注意，矩阵是正方形的，其维数为
     * $m \times m$  。
     * @note
     * 只有在预处理程序已经初始化的情况下才可以调用这个函数。
     *
     */
    size_type
    m() const;

    /**
     * 返回共域（或范围）空间的维数。注意，矩阵是正方形的，其维度为
     * $n \times n$  。
     * @note
     * 只有在预处理程序已经初始化的情况下才可以调用这个函数。
     *
     */
    size_type
    n() const;

  private:
    /**
     * 用于调用cuSPARSE函数的cuSPARSE句柄。
     *
     */
    cusparseHandle_t cusparse_handle;

    /**
     * cuSPARSE对稀疏矩阵的描述  $M=LU$  。
     *
     */
    cusparseMatDescr_t descr_M;

    /**
     * cuSPARSE描述的下三角矩阵  $L$  .
     *
     */
    cusparseMatDescr_t descr_L;

    /**
     * cuSPARSE描述的上三角矩阵  $U$  .
     *
     */
    cusparseMatDescr_t descr_U;

    /**
     * 对 $M=LU$ 的解算和分析结构 .
     *
     */
    csrilu02Info_t info_M;

    /**
     * 下三角矩阵 $L$ 的求解和分析结构。
     *
     */
    csrsv2Info_t info_L;

    /**
     * 上三角矩阵的求解和分析结构  $U$  .
     *
     */
    csrsv2Info_t info_U;

    /**
     * 指向此对象初始化的矩阵的指针。
     *
     */
    SmartPointer<const SparseMatrix<Number>> matrix_pointer;

    /**
     * 指向计算出的预处理矩阵的值（在设备上）的指针。
     *
     */
    std::unique_ptr<Number[], void (*)(Number *)> P_val_dev;

    /**
     * 指向此对象初始化的稀疏矩阵的行指针的指针（在设备上）。由matrix_pointer守护。
     *
     */
    const int *P_row_ptr_dev;

    /**
     * 指向此对象初始化的稀疏矩阵的列索引（在设备上）的指针。由
     * matrix_pointer 保护。
     *
     */
    const int *P_column_index_dev;

    /**
     * 指向vmult()中使用的临时(辅助)向量的值(在设备上)的指针。
     *
     */
    std::unique_ptr<Number[], void (*)(Number *)> tmp_dev;

    /**
     * 指向一个内部缓冲区的指针（在设备上），用于计算分解。
     *
     */
    std::unique_ptr<void, void (*)(void *)> buffer_dev;

    /**
     * 确定是否应该为下三角矩阵生成级别信息  $L$
     * 。这个值可以通过一个AdditionalData对象来修改。
     *
     */
    cusparseSolvePolicy_t policy_L;

    /**
     * 确定是否应该为上三角矩阵生成水平信息  $U$
     * 。这个值可以通过一个AdditionalData对象来修改。
     *
     */
    cusparseSolvePolicy_t policy_U;

    /**
     * 确定是否应该为  $M=LU$
     * 生成电平信息。这个值可以通过一个AdditionalData对象来修改。
     *
     */
    cusparseSolvePolicy_t policy_M;

    /**
     * 行数与此对象初始化的矩阵的行数相同。
     *
     */
    int n_rows;

    /**
     * 非零元素的数量与此对象初始化的矩阵相同。
     *
     */
    int n_nonzero_elements;
  };

   /*--------------------------- inline functions ----------------------------*/ 

#  ifndef DOXYGEN
  template <typename Number>
  inline typename PreconditionIC<Number>::size_type
  PreconditionIC<Number>::m() const
  {
    return n_rows;
  }



  template <typename Number>
  inline typename PreconditionIC<Number>::size_type
  PreconditionIC<Number>::n() const
  {
    return n_rows;
  }



  template <typename Number>
  inline typename PreconditionILU<Number>::size_type
  PreconditionILU<Number>::m() const
  {
    return n_rows;
  }



  template <typename Number>
  inline typename PreconditionILU<Number>::size_type
  PreconditionILU<Number>::n() const
  {
    return n_rows;
  }
#  endif // DOXYGEN

} // namespace CUDAWrappers

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_CUDA

#endif // dealii_cuda_precondition_h


