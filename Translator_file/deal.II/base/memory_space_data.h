//include/deal.II-translator/base/memory_space_data_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2020 by the deal.II authors
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


#ifndef dealii_memory_space_data_h
#define dealii_memory_space_data_h

#include <deal.II/base/config.h>

#include <deal.II/base/cuda.h>
#include <deal.II/base/exceptions.h>

#include <functional>
#include <memory>

DEAL_II_NAMESPACE_OPEN

 /**
  */
namespace MemorySpace
{
  /**
   * 数据结构
   *
   */
  template <typename Number, typename MemorySpace>
  struct MemorySpaceData
  {
    MemorySpaceData()
    {
      static_assert(std::is_same<MemorySpace, Host>::value ||
                      std::is_same<MemorySpace, CUDA>::value,
                    "MemorySpace should be Host or CUDA");
    }

    /**
     * 将活动数据（对于Host来说是values，对于CUDA来说是values_dev）复制到
     * @p begin.  如果数据在设备上，它将被移到主机上。
     *
     */
    void
    copy_to(Number *begin, std::size_t n_elements)
    {
      (void)begin;
      (void)n_elements;
    }

    /**
     * 将 @p begin
     * 中的数据复制到结构的活动数据中（对于Host来说是数值，对于CUDA来说是数值_dev）。指针
     * @p begin 必须在主机上。
     *
     */
    void
    copy_from(Number *begin, std::size_t n_elements)
    {
      (void)begin;
      (void)n_elements;
    }

    /**
     * 指向主机上的数据的指针。
     *
     */
    std::unique_ptr<Number[], std::function<void(Number *)>> values;

    /**
     * 指向设备上的数据的指针。
     *
     */
    std::unique_ptr<Number[]> values_dev;

    /**
     * 指向共享同一内存的进程的数据的指针。
     *
     */
    std::vector<ArrayView<const Number>> values_sm;
  };



  /**
   * 类似于 std::swap. 的互换功能。
   *
   */
  template <typename Number, typename MemorySpace>
  inline void
  swap(MemorySpaceData<Number, MemorySpace> &,
       MemorySpaceData<Number, MemorySpace> &)
  {
    static_assert(std::is_same<MemorySpace, Host>::value ||
                    std::is_same<MemorySpace, CUDA>::value,
                  "MemorySpace should be Host or CUDA");
  }

#ifndef DOXYGEN

  template <typename Number>
  struct MemorySpaceData<Number, Host>
  {
    MemorySpaceData()
      : values(nullptr, &std::free)
    {}

    void
    copy_to(Number *begin, std::size_t n_elements)
    {
      std::copy(values.get(), values.get() + n_elements, begin);
    }

    void
    copy_from(Number *begin, std::size_t n_elements)
    {
      std::copy(begin, begin + n_elements, values.get());
    }

    std::unique_ptr<Number[], std::function<void(Number *)>> values;

    // This is not used but it allows to simplify the code until we start using
    // CUDA-aware MPI.
    std::unique_ptr<Number[]> values_dev;

    std::vector<ArrayView<const Number>> values_sm;
  };



  template <typename Number>
  inline void
  swap(MemorySpaceData<Number, Host> &u, MemorySpaceData<Number, Host> &v)
  {
    std::swap(u.values, v.values);
  }



#  ifdef DEAL_II_COMPILER_CUDA_AWARE

  template <typename Number>
  struct MemorySpaceData<Number, CUDA>
  {
    MemorySpaceData()
      : values(nullptr, &std::free)
      , values_dev(nullptr, Utilities::CUDA::delete_device_data<Number>)
    {}

    void
    copy_to(Number *begin, std::size_t n_elements)
    {
      const cudaError_t cuda_error_code =
        cudaMemcpy(begin,
                   values_dev.get(),
                   n_elements * sizeof(Number),
                   cudaMemcpyDeviceToHost);
      AssertCuda(cuda_error_code);
    }

    void
    copy_from(Number *begin, std::size_t n_elements)
    {
      const cudaError_t cuda_error_code =
        cudaMemcpy(values_dev.get(),
                   begin,
                   n_elements * sizeof(Number),
                   cudaMemcpyHostToDevice);
      AssertCuda(cuda_error_code);
    }

    std::unique_ptr<Number[], std::function<void(Number *)>> values;
    std::unique_ptr<Number[], void (*)(Number *)>            values_dev;

    /**
     * 目前没有使用这个功能。
     *
     */
    std::vector<ArrayView<const Number>> values_sm;
  };



  template <typename Number>
  inline void
  swap(MemorySpaceData<Number, CUDA> &u, MemorySpaceData<Number, CUDA> &v)
  {
    std::swap(u.values, v.values);
    std::swap(u.values_dev, v.values_dev);
  }

#  endif

#endif

} // namespace MemorySpace

DEAL_II_NAMESPACE_CLOSE

#endif


