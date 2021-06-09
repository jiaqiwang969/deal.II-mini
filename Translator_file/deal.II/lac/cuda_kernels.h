//include/deal.II-translator/lac/cuda_kernels_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2019 by the deal.II authors
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

#ifndef dealii_cuda_kernels_h
#define dealii_cuda_kernels_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_COMPILER_CUDA_AWARE


#  include <deal.II/base/cuda_size.h>
#  include <deal.II/base/types.h>

#  include <deal.II/lac/cuda_atomic.h>

#  include <assert.h>

DEAL_II_NAMESPACE_OPEN

namespace LinearAlgebra
{
  namespace CUDAWrappers
  {
    /**
     * 包含CUDA内核的命名空间。
     *
     */
    namespace kernel
    {
      using ::dealii::CUDAWrappers::block_size;
      using ::dealii::CUDAWrappers::chunk_size;
      using ::dealii::CUDAWrappers::warp_size;
      using size_type = types::global_dof_index;

      /**
       * 将大小为 @p N 的 @p val 的每个条目乘以 @p a.  。
       * @ingroup CUDAWrappers
       *
       */
      template <typename Number>
      __global__ void
      vec_scale(Number *val, const Number a, const size_type N);



      /**
       * 定义两个数字的加法的向量。
       * @ingroup CUDAWrappers
       *
       */
      template <typename Number>
      struct Binop_Addition
      {
        __device__ static inline Number
        operation(const Number a, const Number b)
        {
          return a + b;
        }
      };

      template <typename Number>
      struct Binop_Addition<std::complex<Number>>
      {
        __device__ static inline std::complex<Number>
        operation(const std::complex<Number> a, const std::complex<Number>)
        {
          printf("This function is not implemented for std::complex<Number>!");
          assert(false);
          return a;
        }
      };



      /**
       * 定义两个数字的减法的向量。
       * @ingroup CUDAWrappers
       *
       */
      template <typename Number>
      struct Binop_Subtraction
      {
        __device__ static inline Number
        operation(const Number a, const Number b)
        {
          return a - b;
        }
      };

      template <typename Number>
      struct Binop_Subtraction<std::complex<Number>>
      {
        __device__ static inline std::complex<Number>
        operation(const std::complex<Number> a, const std::complex<Number> b)
        {
          printf("This function is not implemented for std::complex<Number>!");
          assert(false);
          return a;
        }
      };



      /**
       * 定义两个数的最大值的向量。
       * @ingroup CUDAWrappers
       *
       */
      template <typename Number>
      struct Binop_Max
      {
        __device__ static inline Number
        operation(const Number a, const Number b)
        {
          return a > b ? a : b;
        }
      };

      template <typename Number>
      struct Binop_Max<std::complex<Number>>
      {
        __device__ static inline std::complex<Number>
        operation(const std::complex<Number> a, const std::complex<Number>)
        {
          printf("This function is not implemented for std::complex<Number>!");
          assert(false);
          return a;
        }
      };



      /**
       * 定义两个数字的最大值的向量。
       * @ingroup CUDAWrappers
       *
       */
      template <typename Number>
      struct Binop_Min
      {
        __device__ static inline Number
        operation(const Number a, const Number b)
        {
          return a > b ? b : a;
        }
      };

      template <typename Number>
      struct Binop_Min<std::complex<Number>>
      {
        __device__ static inline std::complex<Number>
        operation(const std::complex<Number> a, const std::complex<Number>)
        {
          printf("This function is not implemented for std::complex<Number>!");
          assert(false);
          return a;
        }
      };



      /**
       * 对 @p v1 和 @p v2. 中的每个元素应用向量 @p Binop 。
       * @ingroup CUDAWrappers
       *
       */
      template <typename Number, template <typename> class Binop>
      __global__ void
      vector_bin_op(Number *v1, const Number *v2, const size_type N);



      /**
       * 对指数在 @p mask 和 @p v2. 中的 @p v1 的元素应用漏斗 @p
       * Binop ， @p mask 的大小应该大于 @p v1. 的大小 @p mask 和
       * @p v2 的大小应该相同 @p  N。
       * @ingroup CUDAWrappers
       *
       */
      template <typename Number, template <typename> class Binop>
      __global__ void
      masked_vector_bin_op(const unsigned int *mask,
                           Number *            v1,
                           const Number *      v2,
                           const size_type     N);



      /**
       * 实现在使用缩减时用于添加元素的函数的结构。
       * @ingroup CUDAWrappers
       *
       */
      template <typename Number>
      struct ElemSum
      {
        __device__ static Number
        reduction_op(const Number a, const Number b);

        __device__ static Number
        atomic_op(Number *dst, const Number a);

        __device__ static Number
        element_wise_op(const Number a);

        __device__ static Number
        null_value();
      };



      /**
       * 实现在使用还原时计算L1准则的函数的结构。
       * @ingroup CUDAWrappers
       *
       */
      template <typename Number>
      struct L1Norm
      {
        __device__ static Number
        reduction_op(const Number a, const Number b);

        __device__ static Number
        atomic_op(Number *dst, const Number a);

        __device__ static Number
        element_wise_op(const Number a);

        __device__ static Number
        null_value();
      };



      /**
       * 实现在使用还原法时用于计算L-无穷大规范的函数的结构。
       * @ingroup CUDAWrappers
       *
       */
      template <typename Number>
      struct LInfty
      {
        __device__ static Number
        reduction_op(const Number a, const Number b);

        __device__ static Number
        atomic_op(Number *dst, const Number a);

        __device__ static Number
        element_wise_op(const Number a);

        __device__ static Number
        null_value();
      };



      /**
       * 使用 @p Operation. 对 @p v 进行还原。
       * @ingroup CUDAWrappers
       *
       */
      template <typename Number, typename Operation>
      __global__ void
      reduction(Number *result, const Number *v, const size_type N);



      /**
       * 在使用双向量还原时，实现用于计算点乘法线的函数的结构。
       * @ingroup CUDAWrappers
       *
       */
      template <typename Number>
      struct DotProduct
      {
        __device__ static Number
        binary_op(const Number a, const Number b);

        __device__ static Number
        reduction_op(const Number a, const Number b);

        __device__ static Number
        atomic_op(Number *dst, const Number a);

        __device__ static Number
        null_value();
      };



      /**
       * 对 @p v1 和 @p v2
       * 中的每个元素进行二进制运算，然后对所得数组进行还原。
       * @ingroup CUDAWrappers
       *
       */
      template <typename Number, typename Operation>
      __global__ void
      double_vector_reduction(Number *        result,
                              const Number *  v1,
                              const Number *  v2,
                              const size_type N);



      /**
       * 将 @p a 添加到 @p val. 的每个元素中。
       * @ingroup CUDAWrappers
       *
       */
      template <typename Number>
      __global__ void
      vec_add(Number *val, const Number a, const size_type N);



      /**
       * 矢量的倍数相加，即：<tt>val += a*V_val</tt>。
       * @ingroup CUDAWrappers
       *
       */
      template <typename Number>
      __global__ void
      add_aV(Number *        val,
             const Number    a,
             const Number *  V_val,
             const size_type N);



      /**
       * 多重缩放向量的加法，即：<tt>val += a*V_val +
       * b*W_val</tt>。
       * @ingroup CUDAWrappers
       *
       */
      template <typename Number>
      __global__ void
      add_aVbW(Number *        val,
               const Number    a,
               const Number *  V_val,
               const Number    b,
               const Number *  W_val,
               const size_type N);



      /**
       * 缩放和简单添加一个向量的倍数，即<tt>val = = s*val +
       * a*V_val</tt>。
       * @ingroup CUDAWrappers
       *
       */
      template <typename Number>
      __global__ void
      sadd(const Number    s,
           Number *        val,
           const Number    a,
           const Number *  V_val,
           const size_type N);



      /**
       * 缩放和多次添加缩放的向量，即<tt>val = = s*val + a*V_val
       * + b*W_val</tt>。
       * @ingroup CUDAWrappers
       *
       */
      template <typename Number>
      __global__ void
      sadd(const Number    s,
           Number *        val,
           const Number    a,
           const Number *  V_val,
           const Number    b,
           const Number *  W_val,
           const size_type N);



      /**
       * 用参数中的相应元素来缩放这个向量的每个元素。
       * @ingroup CUDAWrappers
       *
       */
      template <typename Number>
      __global__ void
      scale(Number *val, const Number *V_val, const size_type N);



      /**
       * 赋值 <tt>val = a*V_val</tt>.
       * @ingroup CUDAWrappers
       *
       */
      template <typename Number>
      __global__ void
      equ(Number *val, const Number a, const Number *V_val, const size_type N);



      /**
       * 赋值 <tt>val = a*V_val + b*W_val</tt>.
       * @ingroup CUDAWrappers
       *
       */
      template <typename Number>
      __global__ void
      equ(Number *        val,
          const Number    a,
          const Number *  V_val,
          const Number    b,
          const Number *  W_val,
          const size_type N);



      /**
       * 执行一个向量加法和后续内积的组合操作，返回内积的值。
       * @ingroup CUDAWrappers
       *
       */
      template <typename Number>
      __global__ void
      add_and_dot(Number *        res,
                  Number *        v1,
                  const Number *  v2,
                  const Number *  v3,
                  const Number    a,
                  const size_type N);



      /**
       * 将 @p val 的每个元素设置为 @p s. 。
       * @ingroup CUDAWrappers
       *
       */
      template <typename Number>
      __global__ void
      set(Number *val, const Number s, const size_type N);


      /**
       * 将 @p val 中的每个元素设置为 @p v ，使用 @p indices
       * 作为permutation，即：<tt>val[indices[i]] = v[i]</tt>。
       * @ingroup CUDAWrappers
       *
       */
      template <typename Number, typename IndexType>
      __global__ void
      set_permutated(const IndexType *indices,
                     Number *         val,
                     const Number *   v,
                     const IndexType  N);



      /**
       * 将 @p val 中的每个元素设置为 @p v ，使用 @p indices
       * 作为permutation，即，<tt>val[i]=v[indices[i]]</tt>。
       * @ingroup CUDAWrappers
       *
       */
      template <typename Number, typename IndexType>
      __global__ void
      gather(Number *         val,
             const IndexType *indices,
             const Number *   v,
             const IndexType  N);



      /**
       * 将 @p val 中的每个元素添加到 @p v 中，使用 @p indices
       * 作为permutation，即：<tt>val[indices[i]] += v[i]</tt>。+=
       * v[i]/tt>。
       * @ingroup CUDAWrappers
       *
       */
      template <typename Number>
      __global__ void
      add_permutated(const size_type *indices,
                     Number *         val,
                     const Number *   v,
                     const size_type  N);
    } // namespace kernel
  }   // namespace CUDAWrappers
} // namespace LinearAlgebra

DEAL_II_NAMESPACE_CLOSE

#endif

#endif


