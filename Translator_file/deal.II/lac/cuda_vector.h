//include/deal.II-translator/lac/cuda_vector_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2021 by the deal.II authors
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

#ifndef dealii_cuda_vector_h
#define dealii_cuda_vector_h

#include <deal.II/base/config.h>

#include <deal.II/base/communication_pattern_base.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/index_set.h>

#include <deal.II/lac/vector_operation.h>
#include <deal.II/lac/vector_space_vector.h>

#ifdef DEAL_II_WITH_CUDA

DEAL_II_NAMESPACE_OPEN

// Forward declarations
#  ifndef DOXYGEN
template <typename Number>
class ReadWriteVector;
#  endif

namespace LinearAlgebra
{
  /**
   * 一个用于CUDA向量的命名空间。
   *
   */
  namespace CUDAWrappers
  {
    /**
     * 该类使用CUDA实现了一个在Nvidia
     * GPU上使用的向量。这个类是由
     * LinearAlgebra::VectorSpaceVector 类派生的。
     * @note  只支持float和double。          @see  CUDAWrappers
     * @ingroup Vectors
     *
     */
    template <typename Number>
    class Vector : public VectorSpaceVector<Number>
    {
    public:
      using value_type = typename VectorSpaceVector<Number>::value_type;
      using size_type  = typename VectorSpaceVector<Number>::size_type;
      using real_type  = typename VectorSpaceVector<Number>::real_type;

      /**
       * 构造函数。创建一个维度为0的向量。
       *
       */
      Vector();

      /**
       * 复制构造函数。
       *
       */
      Vector(const Vector<Number> &V);

      /**
       * 移动构造函数。
       *
       */
      Vector(Vector<Number> &&) noexcept = default;

      /**
       * 构造函数。设置维度为 @p n
       * 并将所有元素初始化为0。
       * 构造函数是明确的，以避免像这样的意外。
       * <tt>v=0;</tt>。据推测，用户希望将向量的每一个元素都设置为零，但相反，发生的是这样的调用。
       * <tt>v=向量 @<Number@>(0);</tt>,
       * 即向量被一个长度为零的向量所取代。
       *
       */
      explicit Vector(const size_type n);

      /**
       * 复制赋值运算符。
       *
       */
      Vector &
      operator=(const Vector<Number> &v);

      /**
       * 移动赋值运算符。
       *
       */
      Vector &
      operator=(Vector<Number> &&v) noexcept = default;

      /**
       * 交换这个向量和另一个向量的内容  @p v.
       * 人们可以用一个临时变量和复制数据元素来完成这个操作，但是这个函数明显更有效率，因为它只交换了两个向量的数据指针，因此不需要分配临时存储和移动数据。
       * 这个函数类似于所有C++标准容器的 @p swap
       * 函数。另外，还有一个全局函数<tt>swap(u,v)</tt>，它简单地调用<tt>u.swap(v)</tt>，同样与标准函数类似。
       * 这个函数是虚拟的，以便允许派生类单独处理内存。
       *
       */
      virtual void
      swap(Vector<Number> &v);

      /**
       * 重新启动功能。标志<tt>omit_zeroing_entries</tt>决定了向量是否应该被填充为零（false）或保持不动（true）。
       *
       */
      void
      reinit(const size_type n, const bool omit_zeroing_entries = false);

      /**
       * 将维度改为向量V的维度，V的元素不会被复制。
       *
       */
      virtual void
      reinit(const VectorSpaceVector<Number> &V,
             const bool omit_zeroing_entries = false) override;

      /**
       * 从输入向量中导入所有元素  @p V.   VectorOperation::values
       * @p operation  用于决定int  @p V
       * 的元素是否应该添加到当前向量中，或者替换当前元素。最后一个参数不使用。它只用于分布式向量。这是一个应该用来复制一个向量到GPU的函数。
       *
       */
      virtual void
      import(const ReadWriteVector<Number> &V,
             VectorOperation::values        operation,
             std::shared_ptr<const Utilities::MPI::CommunicationPatternBase>
               communication_pattern = {}) override;

      /**
       * 将向量的所有元素设置为标量 @p s.  只有在 @p s
       * 等于0的情况下才允许进行这个操作。
       *
       */
      virtual Vector<Number> &
      operator=(const Number s) override;

      /**
       * 将entive向量乘以一个固定系数。
       *
       */
      virtual Vector<Number> &
      operator*=(const Number factor) override;

      /**
       * 用整个向量除以一个固定因子。
       *
       */
      virtual Vector<Number> &
      operator/=(const Number factor) override;

      /**
       * 将向量 @p V 添加到现在的向量中。
       *
       */
      virtual Vector<Number> &
      operator+=(const VectorSpaceVector<Number> &V) override;

      /**
       * 从现在的向量中减去向量 @p V 。
       *
       */
      virtual Vector<Number> &
      operator-=(const VectorSpaceVector<Number> &V) override;

      /**
       * 返回两个向量的标量乘积。
       *
       */
      virtual Number
      operator*(const VectorSpaceVector<Number> &V) const override;

      /**
       * 将 @p to 的所有成分相加。注意 @p a
       * 是一个标量而不是一个向量。
       *
       */
      virtual void
      add(const Number a) override;

      /**
       * 矢量的倍数的简单相加，即<tt>*this += a*V</tt>。
       *
       */
      virtual void
      add(const Number a, const VectorSpaceVector<Number> &V) override;

      /**
       * 缩放向量的多次加法，即：<tt>*this += a*V+b*W</tt>。
       *
       */
      virtual void
      add(const Number                     a,
          const VectorSpaceVector<Number> &V,
          const Number                     b,
          const VectorSpaceVector<Number> &W) override;

      /**
       * 缩放和简单的向量倍数加法，即<tt>*this =
       * s*(*this)+a*V</tt>。
       *
       */
      virtual void
      sadd(const Number                     s,
           const Number                     a,
           const VectorSpaceVector<Number> &V) override;

      /**
       * 用参数中的相应元素来缩放这个向量的每个元素。这个函数主要是为了模拟对角线缩放矩阵的乘法（和立即重新分配）。
       *
       */
      virtual void
      scale(const VectorSpaceVector<Number> &scaling_factors) override;

      /**
       * 赋值 <tt>*this = a*V</tt>.
       *
       */
      virtual void
      equ(const Number a, const VectorSpaceVector<Number> &V) override;

      /**
       * 返回向量是否只包含值为0的元素。
       *
       */
      virtual bool
      all_zero() const override;

      /**
       * 返回这个向量的所有条目的平均值。
       *
       */
      virtual value_type
      mean_value() const override;

      /**
       * 返回该向量的l<sub>1</sub>准则（即所有条目在所有处理器中的绝对值之和）。
       *
       */
      virtual real_type
      l1_norm() const override;

      /**
       * 返回向量的l<sub>2</sub>准则（即所有处理器中所有条目的平方之和的平方根）。
       *
       */
      virtual real_type
      l2_norm() const override;

      /**
       * 返回 $l_2$ -norm的平方。
       *
       */
      real_type
      norm_sqr() const;

      /**
       * 返回向量的最大规范（即所有条目和所有处理器中的最大绝对值）。
       *
       */
      virtual real_type
      linfty_norm() const override;

      /**
       * 执行一个向量加法和后续内积的组合操作，返回内积的值。换句话说，这个函数的结果与用户调用的
       * @code
       * this->add(a, V);
       * return_value =this W;
       * @endcode
       * 这个函数存在的原因是这个操作比单独调用这两个函数涉及的内存转移要少。这个方法只需要加载三个向量，
       * @p this,   @p V,   @p W,
       * ，而调用单独的方法意味着要加载调用向量 @p this
       * 两次。由于大多数向量操作都有内存传输限制，这就使时间减少了25\%（如果
       * @p W 等于 @p this). ，则减少50\%）。
       * 对于复值向量，第二步的标量乘法被实现为
       * $\left<v,w\right>=\sum_i v_i \bar{w_i}$  。
       *
       */
      virtual Number
      add_and_dot(const Number                     a,
                  const VectorSpaceVector<Number> &V,
                  const VectorSpaceVector<Number> &W) override;

      /**
       * 返回指向底层数组的指针。所有权仍然在这个类中。
       *
       */
      Number *
      get_values() const;

      /**
       * 返回向量的大小。
       *
       */
      virtual size_type
      size() const override;

      /**
       * 返回一个索引集，描述这个向量的哪些元素为当前处理器所拥有，即[0,
       * size]。
       *
       */
      virtual dealii::IndexSet
      locally_owned_elements() const override;

      /**
       * 打印向量到输出流  @p out.  。
       *
       */
      virtual void
      print(std::ostream &     out,
            const unsigned int precision  = 2,
            const bool         scientific = true,
            const bool         across     = true) const override;

      /**
       * 返回这个类的内存消耗，单位是字节。
       *
       */
      virtual std::size_t
      memory_consumption() const override;

      /**
       * 试图在两个不兼容的向量类型之间进行操作。
       * @ingroup Exceptions
       *
       */
      DeclException0(ExcVectorTypeNotCompatible);

    private:
      /**
       * 指向此向量的元素数组的指针。
       *
       */
      std::unique_ptr<Number[], void (*)(Number *)> val;

      /**
       * 矢量中的元素数量。
       *
       */
      size_type n_elements;
    };
  } // namespace CUDAWrappers
} // namespace LinearAlgebra

// ---------------------------- Inline functions --------------------------

/**
 * 全局函数 @p swap
 * ，它重载了C++标准库的默认实现，它使用一个临时对象。该函数简单地交换了两个向量的数据。
 * @relatesalso  Vectors
 *
 *
 */
template <typename Number>
inline void
swap(LinearAlgebra::CUDAWrappers::Vector<Number> &u,
     LinearAlgebra::CUDAWrappers::Vector<Number> &v)
{
  u.swap(v);
}

namespace LinearAlgebra
{
  namespace CUDAWrappers
  {
    template <typename Number>
    inline Number *
    Vector<Number>::get_values() const
    {
      return val.get();
    }



    template <typename Number>
    inline typename Vector<Number>::size_type
    Vector<Number>::size() const
    {
      return n_elements;
    }


    template <typename Number>
    inline IndexSet
    Vector<Number>::locally_owned_elements() const
    {
      return complete_index_set(n_elements);
    }



    template <typename Number>
    inline void
    Vector<Number>::swap(Vector<Number> &v)
    {
      std::swap(val, v.val);
      std::swap(n_elements, v.n_elements);
    }
  } // namespace CUDAWrappers
} // namespace LinearAlgebra

DEAL_II_NAMESPACE_CLOSE

#endif

#endif


