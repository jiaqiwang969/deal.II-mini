//include/deal.II-translator/lac/la_vector_0.txt
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

#ifndef dealii_la_vector_h
#define dealii_la_vector_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/logstream.h>

#include <deal.II/lac/read_write_vector.h>
#include <deal.II/lac/vector_operation.h>
#include <deal.II/lac/vector_space_vector.h>
#include <deal.II/lac/vector_type_traits.h>

// boost::serialization::make_array used to be in array.hpp, but was
// moved to a different file in BOOST 1.64
#include <boost/version.hpp>
#if BOOST_VERSION >= 106400
#  include <boost/serialization/array_wrapper.hpp>
#else
#  include <boost/serialization/array.hpp>
#endif
#include <boost/serialization/split_member.hpp>

#include <cstdio>
#include <cstring>
#include <iostream>
#include <vector>


DEAL_II_NAMESPACE_OPEN

/**
 * 一个用于向量类的命名空间。
 * 这个命名空间包含了各种类，它们为来自不同的外部库（如Trilinos
 * (EPetra)或PETSc）和本地实现（如
 * LinearAlgebra::distributed::Vector. ）的向量类提供包装。
 * 不同的向量类派生自VectorSpaceVector，为向量空间操作提供一个联合接口，或派生自ReadWriteVector（或ReadWriteVector本身），或两者皆有。通过VectorSpaceVector的向量空间操作（如规范或向量添加）和通过ReadWriteVector的元素访问的分离是设计上的，可以提高性能。
 *
 *
 */
namespace LinearAlgebra
{
  /*!   @addtogroup Vectors  
     * @{   
*
*/

  /**
   * 数值化的数据向量。该类同时派生自
   * ::dealii::LinearAlgebra::ReadWriteVector 和
   * ::dealii::LinearAlgebra::VectorSpaceVector.
   * 与C++标准库中的数组不同，该类实现了适合数值计算的向量空间的元素。
   *
   */
  template <typename Number>
  class Vector : public ReadWriteVector<Number>,
                 public VectorSpaceVector<Number>
  {
  public:
    using size_type  = types::global_dof_index;
    using value_type = typename ReadWriteVector<Number>::value_type;

    /**
     * 构造函数。创建一个维数为0的向量。
     *
     */
    Vector() = default;

    /**
     * 复制构造函数。将维数设置为给定的向量的维数，并复制所有元素。
     *
     */
    Vector(const Vector<Number> &V);

    /**
     * 构造函数。将维度设置为 @p n
     * 并将所有元素初始化为零。
     * 构造函数是明确的，以避免像这样的意外。
     * <tt>v=0;</tt>。据推测，用户希望将向量的每一个元素都设置为零，但相反，发生的是这样的调用。
     * <tt>v=向量 @<Number@>(0);</tt>,
     * 即向量被一个长度为0的向量所取代。
     *
     */
    explicit Vector(const size_type n);

    /**
     * 用迭代器指向的给定范围的值初始化向量。这个函数的存在是与
     * @p std::vector 类相类似的。
     *
     */
    template <typename InputIterator>
    Vector(const InputIterator first, const InputIterator last);

    /**
     * 将向量的全局大小设置为 @p size.
     * 存储的元素的索引在[0,size]。        如果标志 @p
     * omit_zeroing_entries
     * 被设置为false，内存将被初始化为0，否则内存将不被触动（用户在使用前必须确保用合理的数据填充它）。
     *
     */
    virtual void
    reinit(const size_type size,
           const bool      omit_zeroing_entries = false) override;

    /**
     * 使用与输入向量 @p in_vector
     * 相同的IndexSet，并为该向量分配内存。        如果标志
     * @p omit_zeroing_entries
     * 被设置为false，内存将被初始化为0，否则内存将不被触动（用户在使用前必须确保用合理的数据填充它）。
     *
     */
    template <typename Number2>
    void
    reinit(const ReadWriteVector<Number2> &in_vector,
           const bool                      omit_zeroing_entries = false);

    /**
     * 初始化向量。指数由 @p  locally_stored_indices指定。
     * 如果标志 @p omit_zeroing_entries
     * 被设置为false，内存将被初始化为0，否则内存将不被触动（用户在使用前必须确保用合理的数据填充它）。
     * local_stored_indices。
     *
     */
    virtual void
    reinit(const IndexSet &locally_stored_indices,
           const bool      omit_zeroing_entries = false) override;


    /**
     * 将维度改为向量V的维度，V的元素不会被复制。
     *
     */
    virtual void
    reinit(const VectorSpaceVector<Number> &V,
           const bool omit_zeroing_entries = false) override;

    /**
     * 返回`false`，因为这是一个串行向量。
     * 这个功能只有在使用基于MPI的向量时才需要调用，并且为了兼容而存在于其他对象中。
     *
     */
    bool
    has_ghost_elements() const;

    /**
     * 复制输入向量的数据  @p in_vector. .
     *
     */
    Vector<Number> &
    operator=(const Vector<Number> &in_vector);

    /**
     * 复制输入向量 @p in_vector. 的数据。
     *
     */
    template <typename Number2>
    Vector<Number> &
    operator=(const Vector<Number2> &in_vector);

    /**
     * 将向量的所有元素设置为标量  @p s.  只有当 @p s
     * 等于零时才允许进行此操作。
     *
     */
    virtual Vector<Number> &
    operator=(const Number s) override;

    /**
     * 将整个向量乘以一个固定系数。
     *
     */
    virtual Vector<Number> &
    operator*=(const Number factor) override;

    /**
     * 用整个向量除以一个固定的因子。
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
    virtual Number operator*(const VectorSpaceVector<Number> &V) const override;

    /**
     * 这个函数没有实现，将抛出一个异常。
     *
     */
    virtual void
    import(const ReadWriteVector<Number> &V,
           VectorOperation::values        operation,
           std::shared_ptr<const Utilities::MPI::CommunicationPatternBase>
             communication_pattern = {}) override;

    /**
     * 将 @p a 添加到所有分量中。注意， @p a
     * 是一个标量，而不是一个矢量。
     *
     */
    virtual void
    add(const Number a) override;

    /**
     * 向量的倍数的简单加法，即<tt>*this += a*V</tt>。
     *
     */
    virtual void
    add(const Number a, const VectorSpaceVector<Number> &V) override;

    /**
     * 矢量的倍数加法，即：<tt>*this += a*V+b*W</tt>。
     *
     */
    virtual void
    add(const Number                     a,
        const VectorSpaceVector<Number> &V,
        const Number                     b,
        const VectorSpaceVector<Number> &W) override;

    /**
     * 一个向量的倍数的缩放和简单加法，即<tt>*this =
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
     * 返回该向量的l<sub>1</sub>准则（即所有条目的绝对值之和）。
     *
     */
    virtual typename VectorSpaceVector<Number>::real_type
    l1_norm() const override;

    /**
     * 返回向量的l<sub>2</sub>准则（即所有处理器中所有条目的平方根之和）。
     *
     */
    virtual typename VectorSpaceVector<Number>::real_type
    l2_norm() const override;

    /**
     * 返回向量的最大规范（即所有条目和所有处理器之间的最大绝对值）。
     *
     */
    virtual typename VectorSpaceVector<Number>::real_type
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
     * @p W 等于 @p this).
     * ，则减少50\%）对于复值向量，第二步中的标量乘法被实现为
     * $\left<v,w\right>=\sum_i v_i \bar{w_i}$  。
     *
     */
    virtual Number
    add_and_dot(const Number                     a,
                const VectorSpaceVector<Number> &V,
                const VectorSpaceVector<Number> &W) override;

    /**
     * 返回向量的全局大小，等于所有处理器中本地拥有的索引数之和。
     *
     */
    virtual size_type
    size() const override;

    /**
     * 返回一个索引集，描述这个向量的哪些元素是由当前处理器拥有的。因此，如果这是一个分布式向量，在不同处理器上返回的索引集将形成不相交的集合，加起来就是完整的索引集。很明显，如果一个向量只在一个处理器上创建，那么结果将满足
     * @code
     * vec.locally_owned_elements() == complete_index_set(vec.size())
     * @endcode
     *
     *
     */
    virtual dealii::IndexSet
    locally_owned_elements() const override;

    /**
     * 将向量打印到输出流中  @p out.  。
     *
     */
    virtual void
    print(std::ostream &     out,
          const unsigned int precision  = 3,
          const bool         scientific = true,
          const bool         across     = true) const override;

    /**
     * 打印向量到输出流 @p out ，其格式可以被 numpy::readtxt().
     * 读取
     * 注意，IndexSet不被打印，而只是存储在Vector中的值。要在python中加载向量，只需执行<code>
     * vector = numpy.loadtxt('my_vector.txt') </code>。
     *
     */
    void
    print_as_numpy_array(std::ostream &     out,
                         const unsigned int precision = 9) const;

    /**
     * 将向量全部写入文件。这是在二进制模式下完成的，所以输出结果既不能被人类阅读，也不能（可能）被其他使用不同操作系统或数字格式的计算机阅读。
     *
     */
    void
    block_write(std::ostream &out) const;

    /**
     * 从一个文件中读取一个矢量en块。这是用上述函数的逆运算来完成的，所以它的速度相当快，因为位流没有被解释。
     * 如果有必要，矢量会被调整大小。
     * 一个原始形式的错误检查被执行，它将识别最直截了当的尝试，将一些数据解释为存储在文件中的位流向量，但不会更多。
     *
     */
    void
    block_read(std::istream &in);

    /**
     * 返回这个类的内存消耗，单位是字节。
     *
     */
    virtual std::size_t
    memory_consumption() const override;

    /**
     * 使用[BOOST序列化库](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html)从一个流中写入和读取此对象的数据，以便进行序列化。
     *
     */
    template <typename Archive>
    void
    serialize(Archive &ar, const unsigned int version);

    /**
     * 试图在两个不兼容的向量类型之间执行操作。
     * @ingroup Exceptions
     *
     */
    DeclException0(ExcVectorTypeNotCompatible);

  private:
    // Make all other ReadWriteVector types friends.
    template <typename Number2>
    friend class Vector;
  };

   /*@}*/ 
   /*--------------------------- Inline functions ----------------------------*/ 

  template <typename Number>
  inline Vector<Number>::Vector(const Vector<Number> &V)
    : ReadWriteVector<Number>(V)
  {}



  template <typename Number>
  inline Vector<Number>::Vector(const size_type n)
    : ReadWriteVector<Number>(n)
  {}



  template <typename Number>
  template <typename InputIterator>
  inline Vector<Number>::Vector(const InputIterator first,
                                const InputIterator last)
  {
    this->reinit(complete_index_set(std::distance(first, last)), true);
    std::copy(first, last, this->begin());
  }



  template <typename Number>
  inline typename Vector<Number>::size_type
  Vector<Number>::size() const
  {
    return ReadWriteVector<Number>::size();
  }



  template <typename Number>
  inline dealii::IndexSet
  Vector<Number>::locally_owned_elements() const
  {
    return IndexSet(ReadWriteVector<Number>::get_stored_elements());
  }



  template <typename Number>
  inline void
  Vector<Number>::print(std::ostream &     out,
                        const unsigned int precision,
                        const bool         scientific,
                        const bool) const
  {
    ReadWriteVector<Number>::print(out, precision, scientific);
  }



  template <typename Number>
  template <typename Archive>
  inline void
  Vector<Number>::serialize(Archive &ar, const unsigned int)
  {
    size_type current_size = this->size();
    ar &static_cast<Subscriptor &>(*this);
    ar & this->stored_elements;
    // If necessary, resize the vector during a read operation
    if (this->size() != current_size)
      this->reinit(this->size());
    ar &boost::serialization::make_array(this->values.get(), this->size());
  }



  template <typename Number>
  inline std::size_t
  Vector<Number>::memory_consumption() const
  {
    return ReadWriteVector<Number>::memory_consumption();
  }
} // end of namespace LinearAlgebra


/**
 * 将 dealii::LinearAlgebra::Vector 声明为串行向量。
 *
 *
 */
template <typename Number>
struct is_serial_vector<LinearAlgebra::Vector<Number>> : std::true_type
{};

#ifndef DOXYGEN
 /*----------------------- Inline functions ----------------------------------*/ 

namespace LinearAlgebra
{
  template <typename Number>
  inline bool
  Vector<Number>::has_ghost_elements() const
  {
    return false;
  }
} // namespace LinearAlgebra

#endif


DEAL_II_NAMESPACE_CLOSE

#ifdef DEAL_II_MSVC
#  include <deal.II/lac/la_vector.templates.h>
#endif

#endif


