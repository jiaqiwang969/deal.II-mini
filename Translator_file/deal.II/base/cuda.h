//include/deal.II-translator/base/cuda_0.txt
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

#ifndef dealii_cuda_h
#define dealii_cuda_h

#include <deal.II/base/config.h>

#include <deal.II/base/array_view.h>
#include <deal.II/base/exceptions.h>

#ifdef DEAL_II_COMPILER_CUDA_AWARE
#  include <cusolverDn.h>
#  include <cusolverSp.h>
#  include <cusparse.h>

#  include <vector>

DEAL_II_NAMESPACE_OPEN
namespace Utilities
{
  /**
   * 一个用于CUDA的实用结构的命名空间。
   *
   */
  namespace CUDA
  {
    /**
     * 各种CUDA
     * API需要一个对象来存储内部数据。这个结构为deal.II内部使用的各CUDA库创建、初始化、存储和销毁这些所谓的句柄。
     *
     */
    struct Handle
    {
      /**
       * 构造函数。初始化不同库的句柄。
       *
       */
      Handle();

      /**
       * 复制构造函数被删除。
       *
       */
      Handle(Handle const &) = delete;

      /**
       * 销毁器。销毁句柄。
       *
       */
      ~Handle();

      /**
       * 指针指向一个不透明的 cuSolverDN 上下文。
       * 该句柄必须传递给每个cuSolverDN库函数。
       *
       */
      cusolverDnHandle_t cusolver_dn_handle;

      /**
       * 指向一个不透明的cuSolverSP上下文的指针。
       * 该句柄必须被传递给每个cuSolverSP库函数。
       *
       */
      cusolverSpHandle_t cusolver_sp_handle;

      /**
       * 指向一个不透明的cuSPARSE上下文的指针。
       * 该句柄必须被传递给每个cuSPARSE库函数。
       *
       */
      cusparseHandle_t cusparse_handle;
    };

    /**
     * 在设备上分配 @p n_elements 。
     *
     */
    template <typename T>
    inline void
    malloc(T *&pointer, const unsigned int n_elements)
    {
      cudaError_t cuda_error_code =
        cudaMalloc(&pointer, n_elements * sizeof(T));
      AssertCuda(cuda_error_code);
    }

    /**
     * 释放设备上的内存。
     *
     */
    template <typename T>
    inline void
    free(T *&pointer)
    {
      cudaError_t cuda_error_code = cudaFree(pointer);
      AssertCuda(cuda_error_code);
      pointer = nullptr;
    }

    /**
     * 指向设备内存的 `std::unique_ptr` 要使用的分配器。
     *
     */
    template <typename Number>
    Number *
    allocate_device_data(const std::size_t size)
    {
      Number *device_ptr;
      Utilities::CUDA::malloc(device_ptr, size);
      return device_ptr;
    }

    /**
     * 用于指向设备内存的 `std::unique_ptr` 的删除器。
     *
     */
    template <typename Number>
    void
    delete_device_data(Number *device_ptr) noexcept
    {
      const cudaError_t error_code = cudaFree(device_ptr);
      AssertNothrowCuda(error_code);
    }

    /**
     * 将设备ArrayView  @p in 复制到主机ArrayView  @p out. 。
     *
     */
    template <typename T>
    inline void
    copy_to_host(const ArrayView<const T, MemorySpace::CUDA> &in,
                 ArrayView<T, MemorySpace::Host> &            out)
    {
      AssertDimension(in.size(), out.size());
      cudaError_t cuda_error_code = cudaMemcpy(out.data(),
                                               in.data(),
                                               in.size() * sizeof(T),
                                               cudaMemcpyDeviceToHost);
      AssertCuda(cuda_error_code);
    }

    /**
     * 复制主机ArrayView @p in 到设备ArrayView @p out. 。
     *
     */
    template <typename T>
    inline void
    copy_to_dev(const ArrayView<const T, MemorySpace::Host> &in,
                ArrayView<T, MemorySpace::CUDA> &            out)
    {
      AssertDimension(in.size(), out.size());
      cudaError_t cuda_error_code = cudaMemcpy(out.data(),
                                               in.data(),
                                               in.size() * sizeof(T),
                                               cudaMemcpyHostToDevice);
      AssertCuda(cuda_error_code);
    }

    /**
     * 把 @p pointer_dev 中的元素复制到 @p vector_host.
     * 中的主机上。
     *
     */
    template <typename T>
    inline void
    copy_to_host(const T *pointer_dev, std::vector<T> &vector_host)
    {
      ArrayView<const T, MemorySpace::CUDA> in(pointer_dev, vector_host.size());
      auto                                  out = make_array_view(vector_host);
      copy_to_host(in, out);
    }

    /**
     * 将 @p vector_host 中的元素复制到 @p pointer_dev.
     * 中的设备上
     * 在调用这个函数之前，需要在设备上分配内存。
     *
     */
    template <typename T>
    inline void
    copy_to_dev(const std::vector<T> &vector_host, T *pointer_dev)
    {
      auto                            in = make_array_view(vector_host);
      ArrayView<T, MemorySpace::CUDA> out(pointer_dev, vector_host.size());
      copy_to_dev(in, out);
    }
  } // namespace CUDA
} // namespace Utilities

DEAL_II_NAMESPACE_CLOSE
#endif
#endif


