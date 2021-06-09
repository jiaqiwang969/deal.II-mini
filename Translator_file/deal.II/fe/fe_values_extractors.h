//include/deal.II-translator/fe/fe_values_extractors_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2020 by the deal.II authors
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

#ifndef dealii_fe_values_extractors_h
#define dealii_fe_values_extractors_h


#include <deal.II/base/config.h>

DEAL_II_NAMESPACE_OPEN


/**
 * 一个命名空间，我们在其中声明
 * "提取器"，即当作为FEValues、FEFaceValues和FESubfaceValues对象上的operator[]表达式的子标时，提取矢量值元素的某些成分的类。对这些对象应用提取器的结果是一个来自命名空间FEValuesViews的具有相应类型的对象。有一些提取器用于单个标量分量、由
 * <code>dim</code> 元素组成的向量分量和由 <code>(dim*dim +
 * dim)/2</code>
 * 分量组成的二阶对称张量，以及二阶非对称张量。
 * 我们可以把提取器看作是相当于一个索引，或者一个索引范围。在标量提取器（即
 * FEValuesExtractors::Scalar 类）的情况下，创建一个像（见
 * step-20 的用途）的对象
 *
 * @code
 * const FEValuesExtractors::Scalar pressure(dim);
 * @endcode
 * 可以被认为是创建一个具有`dim`值的单一索引。就其本身而言，一个索引不知道它是什么索引，所以它需要相当于一个数组的东西来提取。因此，假设有一个至少有`dim+1`个向量分量的有限元素（在
 * step-20
 * 中确实有），以及一个对其进行操作的FEValues对象，那么就可以编写
 *
 * @code
 * fe_values[pressure]
 * @endcode
 * 的结果是一个对象，它只代表整个元素的第`dim`个分量的形状函数。在这个例子中，这些将是压力形状函数的值，或者更准确地说：所有形状函数的（标量）压力值（即使是与压力无关的形状函数，但例如速度）。在上面的例子中，如图所示，在`fe_values`对象上使用`operator[]`的结果是
 * FEValuesViews::Scalar. 类型。 同样的，当使用
 *
 * @code
 * const FEValuesExtractors::Vector velocities(0);
 * @endcode
 * 那么这样创建的对象可以被认为是一个<i>index
 * range</i>，从零开始，正好延伸到`dim`组件上。用Matlab的符号，可以写成`0:dim-1`。然后，写
 *
 * @code
 * fe_values[velocities]
 * @endcode
 * 就像在Matlab中写 "array(3:7)
 * "会返回一个长度为5的数组，这个数组是通过查看索引3到7（包括7）从原始数组中提取的。
 * 参见 @ref
 * vector_valued
 * 模块的描述，以了解如何使用该命名空间的功能的例子。
 *
 *
 * @ingroup feaccess vector_valued
 *
 *
 */
namespace FEValuesExtractors
{
  /**
   * 矢量值元素的单一标量分量的提取器。将这种类型的对象应用于FEValues、FEFaceValues或FESubfaceValues对象的结果是 FEValuesViews::Scalar.
   * 提取器的概念在命名空间FEValuesExtractors的文档和 @ref
   * vector_valued 模块中定义。
   * @ingroup feaccess vector_valued
   *
   */
  struct Scalar
  {
    /**
     * 矢量的选定标量分量。
     *
     */
    unsigned int component;

    /**
     * 默认构造函数。用一个无效的分量初始化该对象。
     * 这导致了一个不能使用的对象，但它允许将这种对象放入数组中，在调整数组大小时需要一个默认的构造函数，然后再将一个合适的对象分配给数组的每个元素。
     *
     */
    Scalar();

    /**
     * 构造函数。以选定的向量分量作为参数。
     *
     */
    Scalar(const unsigned int component);

    /**
     * 返回一个字符串，唯一标识这个有限元提取器。
     *
     */
    std::string
    get_name() const;
  };


  /**
   * 向量值元素的 <code>spacedim</code> 分量的提取器。 <code>spacedim</code> 的值由该提取器应用的FEValues对象定义。将这种类型的对象应用于FEValues、FEFaceValues或FESubfaceValues对象的结果是 FEValuesViews::Vector. 类型。 提取器的概念在命名空间FEValuesExtractors的文档和 @ref vector_valued 模块中定义。    请注意，在当前的上下文中，矢量是指物理学上使用的矢量：它具有 <code>spacedim</code> 分量，在坐标系变换下以特定的方式表现出来。例子包括速度或位移场。这与数学中使用 "向量 "
   * 一词的方式相反（以及我们在库中的其他上下文中使用这个词的方式，例如在向量类中），在那里它真正代表了一个数字的集合。后者的一个例子是火焰中化学物种浓度的集合；然而，这些实际上只是标量变量的集合，因为如果坐标系被旋转，它们不会改变，不像速度矢量的分量，因此，这个类不应该被用于这种情况。
   * @ingroup feaccess vector_valued
   *
   */
  struct Vector
  {
    /**
     * 矢量视图的第一个分量。
     *
     */
    unsigned int first_vector_component;

    /**
     * 默认构造函数。用一个无效的分量初始化对象。
     * 这导致了一个不能使用的对象，但它允许将这种对象放入数组中，在调整数组大小时需要一个默认的构造函数，然后再将一个合适的对象分配给数组的每个元素。
     *
     */
    Vector();

    /**
     * 构造函数。以FEValues对象内部所选向量的第一个分量作为参数。
     *
     */
    Vector(const unsigned int first_vector_component);

    /**
     * 返回一个唯一标识这个有限元提取器的字符串。
     *
     */
    std::string
    get_name() const;
  };


  /**
   * 对由模板参数指定的等级的对称张量进行提取。对于一个二阶对称张量，这表示一个矢量元素的 <code>(dim*dim
   * + dim)/2</code> 分量的集合。 <code>dim</code>
   * 的值由提取器所应用的FEValues对象定义。将这种类型的对象应用于FEValues、FEFaceValues或FESubfaceValues对象的结果是
   * FEValuesViews::SymmetricTensor.
   * 提取器的概念在命名空间FEValuesExtractors的文档和 @ref
   * vector_valued 模块中定义。
   * @ingroup feaccess vector_valued
   *
   */
  template <int rank>
  struct SymmetricTensor
  {
    /**
     * 张量视图的第一个组成部分。
     *
     */
    unsigned int first_tensor_component;

    /**
     * 默认构造函数。用一个无效的组件初始化对象。
     * 这导致了一个不能使用的对象，但它允许将这种对象放入数组中，在调整数组大小时需要一个默认的构造函数，然后再将一个合适的对象分配给数组的每个元素。
     *
     */
    SymmetricTensor();

    /**
     * 构造函数。以FEValues对象内部所选张量的第一个分量作为参数。
     *
     */
    SymmetricTensor(const unsigned int first_tensor_component);

    /**
     * 返回一个唯一标识这个有限元提取器的字符串。
     *
     */
    std::string
    get_name() const;
  };


  /**
   * 对由模板参数指定的给定等级的一般张量的提取器。对于一个二阶张量，这代表一个矢量值元素的 <code>(dim*dim)</code>
   * 分量的集合。 <code>dim</code>
   * 的值由提取器所应用的FEValues对象定义。将这种类型的对象应用于FEValues、FEFaceValues或FESubfaceValues对象的结果是
   * FEValuesViews::Tensor.
   * 提取器的概念在命名空间FEValuesExtractors的文档和 @ref
   * vector_valued 模块中定义。
   * @ingroup feaccess vector_valued
   *
   */
  template <int rank>
  struct Tensor
  {
    /**
     * 张量视图的第一个组成部分。
     *
     */
    unsigned int first_tensor_component;

    /**
     * 默认构造函数。用一个无效的组件初始化对象。
     * 这导致了一个不能使用的对象，但它允许将这种对象放入数组中，在调整数组大小时需要一个默认的构造函数，然后再将一个合适的对象分配给数组的每个元素。
     *
     */
    Tensor();

    /**
     * 构造函数。以FEValues对象内部所选张量的第一个分量作为参数。
     *
     */
    Tensor(const unsigned int first_tensor_component);

    /**
     * 返回一个唯一标识此有限元提取器的字符串。
     *
     */
    std::string
    get_name() const;
  };
} // namespace FEValuesExtractors


 /*-------------- Inline functions: namespace FEValuesExtractors -------------*/ 

namespace FEValuesExtractors
{
  inline Scalar::Scalar()
    : component(numbers::invalid_unsigned_int)
  {}



  inline Scalar::Scalar(const unsigned int component)
    : component(component)
  {}



  inline Vector::Vector()
    : first_vector_component(numbers::invalid_unsigned_int)
  {}


  inline Vector::Vector(const unsigned int first_vector_component)
    : first_vector_component(first_vector_component)
  {}


  template <int rank>
  inline SymmetricTensor<rank>::SymmetricTensor()
    : first_tensor_component(numbers::invalid_unsigned_int)
  {}


  template <int rank>
  inline SymmetricTensor<rank>::SymmetricTensor(
    const unsigned int first_tensor_component)
    : first_tensor_component(first_tensor_component)
  {}


  template <int rank>
  inline Tensor<rank>::Tensor()
    : first_tensor_component(numbers::invalid_unsigned_int)
  {}


  template <int rank>
  inline Tensor<rank>::Tensor(const unsigned int first_tensor_component)
    : first_tensor_component(first_tensor_component)
  {}
} // namespace FEValuesExtractors



DEAL_II_NAMESPACE_CLOSE

#endif


