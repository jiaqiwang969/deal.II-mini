//include/deal.II-translator/lac/vector_space_vector_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2021 by the deal.II authors
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

#ifndef dealii_vector_space_vector_h
#define dealii_vector_space_vector_h

#include <deal.II/base/config.h>

#include <deal.II/base/communication_pattern_base.h>
#include <deal.II/base/numbers.h>

#include <deal.II/lac/vector_operation.h>

#include <memory>


DEAL_II_NAMESPACE_OPEN

// Forward declarations
#ifndef DOXYGEN
class IndexSet;
namespace LinearAlgebra
{
  template <typename Number>
  class ReadWriteVector;
} // namespace LinearAlgebra
#endif

namespace LinearAlgebra
{
  /*!   @addtogroup Vectors  
     * @{   
*
*/

  /**
   * VectorSpaceVector是一个抽象类，用于定义向量类在想要实现全局操作时需要实现的接口。这个类是ReadWriteVector的补充，它允许访问单个元素，但不允许全局操作。
   *
   */
  template <typename Number>
  class VectorSpaceVector
  {
  public:
    using value_type = Number;
    using size_type  = types::global_dof_index;
    using real_type  = typename numbers::NumberTraits<Number>::real_type;

    /**
     * 将维度改为向量V的维度，V的元素不会被复制。
     *
     */
    virtual void
    reinit(const VectorSpaceVector<Number> &V,
           const bool                       omit_zeroing_entries = false) = 0;

    /**
     * 将向量的所有元素设置为标量 @p s. ，只有当 @p s
     * 等于零时，才允许该操作。
     *
     */
    virtual VectorSpaceVector<Number> &
    operator=(const Number s) = 0;

    /**
     * 将整个向量乘以一个固定系数。
     *
     */
    virtual VectorSpaceVector<Number> &
    operator*=(const Number factor) = 0;

    /**
     * 用整个向量除以一个固定的因子。
     *
     */
    virtual VectorSpaceVector<Number> &
    operator/=(const Number factor) = 0;

    /**
     * 将向量 @p V 添加到现在的向量中。
     *
     */
    virtual VectorSpaceVector<Number> &
    operator+=(const VectorSpaceVector<Number> &V) = 0;

    /**
     * 从现在的向量中减去向量 @p V 。
     *
     */
    virtual VectorSpaceVector<Number> &
    operator-=(const VectorSpaceVector<Number> &V) = 0;

    /**
     * 从输入向量 @p V.  VectorOperation::values  @p operation
     * 中导入向量的IndexSet中存在的所有元素，用于决定 @p V
     * 中的元素是否应该添加到当前向量中，或者替换当前元素。如果多次使用同一通信模式，可以使用最后一个参数。这可以用来提高性能。
     *
     */
    virtual void
    import(const ReadWriteVector<Number> &V,
           VectorOperation::values        operation,
           std::shared_ptr<const Utilities::MPI::CommunicationPatternBase>
             communication_pattern = {}) = 0;

    /**
     * 返回两个向量的标量乘积。
     *
     */
    virtual Number operator*(const VectorSpaceVector<Number> &V) const = 0;

    /**
     * 将 @p a 添加到所有组件中。注意 @p a
     * 是一个标量而不是一个向量。
     *
     */
    virtual void
    add(const Number a) = 0;

    /**
     * 矢量的倍数的简单加法，即：<tt>*this += a*V</tt>。
     *
     */
    virtual void
    add(const Number a, const VectorSpaceVector<Number> &V) = 0;

    /**
     * 缩放向量的多重加法，即<tt>*this += a*V+b*W</tt>。
     *
     */
    virtual void
    add(const Number                     a,
        const VectorSpaceVector<Number> &V,
        const Number                     b,
        const VectorSpaceVector<Number> &W) = 0;

    /**
     * 缩放和简单的向量倍数加法，即<tt>*this =
     * s*(*this)+a*V</tt>。
     *
     */
    virtual void
    sadd(const Number                     s,
         const Number                     a,
         const VectorSpaceVector<Number> &V) = 0;

    /**
     * 用参数中的相应元素来缩放这个向量的每个元素。这个函数主要是为了模拟对角线缩放矩阵的乘法（和立即重新分配）。
     *
     */
    virtual void
    scale(const VectorSpaceVector<Number> &scaling_factors) = 0;

    /**
     * 赋值 <tt>*this = a*V</tt>.
     *
     */
    virtual void
    equ(const Number a, const VectorSpaceVector<Number> &V) = 0;

    /**
     * 返回向量是否只包含值为0的元素。
     *
     */
    virtual bool
    all_zero() const = 0;

    /**
     * 返回这个向量的所有条目的平均值。
     *
     */
    virtual value_type
    mean_value() const = 0;

    /**
     * 返回该向量的l<sub>1</sub>准则（即所有条目在所有处理器中的绝对值之和）。
     *
     */
    virtual real_type
    l1_norm() const = 0;

    /**
     * 返回向量的l<sub>2</sub>准则（即所有处理器中所有条目的平方之和的平方根）。
     *
     */
    virtual real_type
    l2_norm() const = 0;

    /**
     * 返回向量的最大规范（即所有条目和所有处理器之间的最大绝对值）。
     *
     */
    virtual real_type
    linfty_norm() const = 0;

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
     * 对于复值向量，第二步的标量乘法实现为
     * $\left<v,w\right>=\sum_i v_i \bar{w_i}$  。
     *
     */
    virtual Number
    add_and_dot(const Number                     a,
                const VectorSpaceVector<Number> &V,
                const VectorSpaceVector<Number> &W) = 0;

    /**
     * 这个函数什么都不做，只是为了向后兼容而存在。
     *
     */
    virtual void compress(VectorOperation::values)
    {}

    /**
     * 返回向量的全局大小，等于所有处理器中本地拥有的索引数之和。
     *
     */
    virtual size_type
    size() const = 0;

    /**
     * 返回一个索引集，描述这个向量的哪些元素是由当前处理器拥有的。因此，如果这是一个分布式向量，在不同处理器上返回的索引集将形成不相交的集合，加起来就是完整的索引集。很明显，如果一个向量只在一个处理器上创建，那么结果将满足
     * @code
     * vec.locally_owned_elements() == complete_index_set(vec.size())
     * @endcode
     *
     *
     */
    virtual dealii::IndexSet
    locally_owned_elements() const = 0;

    /**
     * 将向量打印到输出流中  @p out.  。
     *
     */
    virtual void
    print(std::ostream &     out,
          const unsigned int precision  = 3,
          const bool         scientific = true,
          const bool         across     = true) const = 0;

    /**
     * 以字节为单位返回这个类的内存消耗。
     *
     */
    virtual std::size_t
    memory_consumption() const = 0;

    /**
     * 解构器。声明为虚拟，以便继承类（可能管理自己的内存）被正确销毁。
     *
     */
    virtual ~VectorSpaceVector() = default;
  };
   /*@}*/ 
} // namespace LinearAlgebra

// ---------------------------- Free functions --------------------------

namespace LinearAlgebra
{
  /**
   * 将向量的所有条目移动一个常数因子，使向量的平均值变为零。
   *
   */
  template <typename Number>
  void
  set_zero_mean_value(VectorSpaceVector<Number> &vector)
  {
    vector.add(-vector.mean_value());
  }
} // namespace LinearAlgebra

DEAL_II_NAMESPACE_CLOSE

#endif


