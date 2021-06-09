//include/deal.II-translator/base/tensor_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2021 by the deal.II authors
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

#ifndef dealii_tensor_h
#define dealii_tensor_h

#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/numbers.h>
#include <deal.II/base/table_indices.h>
#include <deal.II/base/template_constraints.h>
#include <deal.II/base/tensor_accessors.h>
#include <deal.II/base/utilities.h>

#ifdef DEAL_II_WITH_ADOLC
#  include <adolc/adouble.h> // Taped double
#endif

#include <cmath>
#include <ostream>
#include <utility>
#include <vector>


DEAL_II_NAMESPACE_OPEN

// Forward declarations:
#ifndef DOXYGEN
template <typename ElementType, typename MemorySpace>
class ArrayView;
template <int dim, typename Number>
class Point;
template <int rank_, int dim, typename Number = double>
class Tensor;
template <typename Number>
class Vector;
template <typename number>
class FullMatrix;
namespace Differentiation
{
  namespace SD
  {
    class Expression;
  }
} // namespace Differentiation
#endif


/**
 * 这个类是<tt>Tensor<rank,dim,Number></tt>类的一个专门版本。它处理等级为零的张量，即标量。第二个模板参数
 * @p dim 被忽略了。
 * 这个类的存在是因为在某些情况下，我们想构造张量
 * @<spacedim-dim,dim,Number@>,
 * 类型的对象，它应该扩展为标量、向量、矩阵等，这取决于模板参数
 * @p dim 和 @p spacedim.
 * 的值，因此我们需要一个充当标量的类（即 @p Number)
 * 用于所有目的，但属于Tensor模板家族。
 * @tparam  dim
 * 一个整数，表示该张量所处空间的维度。当然，这等于识别一个点和秩-1张量的坐标数。由于当前对象是一个秩0张量（标量），这个模板参数对这个类没有意义。
 * @tparam  Number
 * 要存储张量元素的数据类型。几乎在所有情况下，这只是默认的
 * @p double,
 * ，但在某些情况下，人们可能希望以不同的（而且总是标量）类型来存储元素。它可以用来将张量建立在
 * @p float 或 @p complex
 * 数字或任何其他实现基本算术运算的数据类型上。另一个例子是允许自动微分的类型（例如，见
 * step-33
 * 中使用的Sacado类型），从而可以生成一个以张量为参数的函数的分析（空间）导数。
 *
 *
 * @ingroup geomprimitives
 *
 */
template <int dim, typename Number>
class Tensor<0, dim, Number>
{
public:
  static_assert(dim >= 0,
                "Tensors must have a dimension greater than or equal to one.");

  /**
   * 提供一种方法来获取一个对象的尺寸，而不需要明确知道它的数据类型。实现这种方式而不是提供一个函数<tt>dimension()</tt>，因为现在有可能在编译时获得尺寸，而不需要内联函数的扩展和预评估；因此编译器可能产生更有效的代码，你可以使用这个值来声明其他数据类型。
   *
   */
  static constexpr unsigned int dimension = dim;

  /**
   * 向外界公布这个张量的等级。
   *
   */
  static constexpr unsigned int rank = 0;

  /**
   * 秩为0的张量的独立成分的数量。
   *
   */
  static constexpr unsigned int n_independent_components = 1;

  /**
   * 声明一个类型，该类型持有与该类的模板参数相同精度的实值数。对于
   * std::complex<number>,
   * ，这对应于类型number，对于所有其他情况，它等于Number。也请参见Vector<Number>中的相应字段。
   * 这个别名是用来表示规范的返回类型的。
   *
   */
  using real_type = typename numbers::NumberTraits<Number>::real_type;

  /**
   * 由这个容器封装并由operator[]()返回的对象的类型。这是一个等级为0的张量的标量数字类型。
   *
   */
  using value_type = Number;

  /**
   * 声明一个数组类型，可以用来静态地初始化这个类型的对象。如果是一个等级为0的张量，这只是标量数字类型Number。
   *
   */
  using array_type = Number;

  /**
   * 构造函数。设置为零。
   * @note  这个函数也可以在CUDA设备代码中使用。
   *
   */
  constexpr DEAL_II_CUDA_HOST_DEV
  Tensor();

  /**
   * 来自不同底层标量类型的张量的构造器。这显然要求 @p
   * OtherNumber 类型可转换为 @p 数字。
   * @note  这个函数也可以在CUDA设备代码中使用。
   *
   */
  template <typename OtherNumber>
  constexpr DEAL_II_CUDA_HOST_DEV
  Tensor(const Tensor<0, dim, OtherNumber> &initializer);

  /**
   * 构造函数，数据从一个C风格的数组中复制出来。
   * @note  这个函数也可以在CUDA设备代码中使用。
   *
   */
  template <typename OtherNumber>
  constexpr DEAL_II_CUDA_HOST_DEV
  Tensor(const OtherNumber &initializer);

#if __GNUC__ >= 11 || defined __INTEL_COMPILER
  /**
   * 拷贝构造函数
   *
   */
  constexpr DEAL_II_CUDA_HOST_DEV
  Tensor(const Tensor<0, dim, Number> &other);

  /**
   * 拷贝赋值操作符
   *
   */
  constexpr DEAL_II_CUDA_HOST_DEV Tensor<0, dim, Number> &
                                  operator=(const Tensor<0, dim, Number> &other);

  /**
   * 移动构造函数
   *
   */
  constexpr DEAL_II_CUDA_HOST_DEV
    Tensor(Tensor<0, dim, Number> &&other) noexcept;

  /**
   * 移动赋值运算符
   *
   */
  constexpr DEAL_II_CUDA_HOST_DEV Tensor<0, dim, Number> &
                                  operator=(Tensor<0, dim, Number> &&other) noexcept;
#endif

  /**
   * 返回一个指向底层存储的第一个元素的指针。
   *
   */
  Number *
  begin_raw();

  /**
   * 返回一个指向底层存储的第一个元素的常量指针。
   *
   */
  const Number *
  begin_raw() const;

  /**
   * 返回一个指向底层存储结束后的元素的指针。
   *
   */
  Number *
  end_raw();

  /**
   * 返回一个超过底层存储结束的元素的常量指针。
   *
   */
  const Number *
  end_raw() const;

  /**
   * 返回一个对封装的Number对象的引用。由于秩0张量是标量，这是一个自然的操作。
   * 这是返回一个可写引用的非const转换操作。
   * @note  这个函数也可以在CUDA设备代码中使用。
   *
   */
  constexpr DEAL_II_CUDA_HOST_DEV operator Number &();

  /**
   * 返回一个对封装的Number对象的引用。由于秩0张量是标量，这是一个自然的操作。
   * 这是一个const转换操作，返回一个只读的引用。
   * @note  这个函数也可以在CUDA设备代码中使用。
   *
   */
  constexpr DEAL_II_CUDA_HOST_DEV operator const Number &() const;

  /**
   * 从具有不同底层标量类型的张量进行赋值。这显然要求
   * @p OtherNumber 类型可以转换为 @p 数。
   * @note 这个函数也可以在CUDA设备代码中使用。
   *
   */
  template <typename OtherNumber>
  constexpr DEAL_II_CUDA_HOST_DEV Tensor &
                                  operator=(const Tensor<0, dim, OtherNumber> &rhs);

  /**
   * 这个操作符将一个标量分配给一个张量。这显然要求 @p
   * OtherNumber 类型可以转换为 @p Number. 。
   * @note  这个函数也可以在CUDA设备代码中使用。
   *
   */
  template <typename OtherNumber>
  constexpr DEAL_II_CUDA_HOST_DEV Tensor &
                                  operator=(const OtherNumber &d);

  /**
   * 测试两个张量的相等。
   *
   */
  template <typename OtherNumber>
  constexpr bool
  operator==(const Tensor<0, dim, OtherNumber> &rhs) const;

  /**
   * 测试两个张量的不等式。
   *
   */
  template <typename OtherNumber>
  constexpr bool
  operator!=(const Tensor<0, dim, OtherNumber> &rhs) const;

  /**
   * 添加另一个标量。
   * @note 这个函数也可以在CUDA设备代码中使用。
   *
   */
  template <typename OtherNumber>
  constexpr DEAL_II_CUDA_HOST_DEV Tensor &
                                  operator+=(const Tensor<0, dim, OtherNumber> &rhs);

  /**
   * 减去另一个标量。
   * @note  这个函数也可以在CUDA设备代码中使用。
   *
   */
  template <typename OtherNumber>
  constexpr DEAL_II_CUDA_HOST_DEV Tensor &
                                  operator-=(const Tensor<0, dim, OtherNumber> &rhs);

  /**
   * 将标量乘以一个<tt>因子</tt>。
   * @note  这个函数也可以在CUDA设备代码中使用。
   *
   */
  template <typename OtherNumber>
  constexpr DEAL_II_CUDA_HOST_DEV Tensor &
                                  operator*=(const OtherNumber &factor);

  /**
   * 将标量除以<tt>因子</tt>。
   * @note  这个函数也可以在CUDA设备代码中使用。
   *
   */
  template <typename OtherNumber>
  constexpr DEAL_II_CUDA_HOST_DEV Tensor &
                                  operator/=(const OtherNumber &factor);

  /**
   * 带有反转项的张量。
   * @note  这个函数也可以在CUDA设备代码中使用。
   *
   */
  constexpr DEAL_II_CUDA_HOST_DEV Tensor
                                  operator-() const;

  /**
   * 将所有数值重置为零。    请注意，这与标准库容器的 @p
   * clear()成员函数和deal.II内的其他几个类的语义部分不一致，它们不仅将存储元素的值重置为零，而且释放所有内存并将对象返回到处女状态。然而，由于本类型的对象的大小是由其模板参数决定的，所以调整大小是不可能的，事实上，所有元素的值都为零的状态就是这样一个对象构造后的状态。
   *
   */
  constexpr void
  clear();

  /**
   * 返回张量的Frobenius-norm，即所有条目的绝对平方之和的平方根。对于目前秩-1张量的情况，这等于通常的<tt>l<sub>2</sub></tt>向量的规范。
   *
   */
  real_type
  norm() const;

  /**
   * 返回张量的Frobenius-norm的平方，即所有条目的绝对平方之和。
   * @note 这个函数也可以在CUDA设备代码中使用。
   *
   */
  constexpr DEAL_II_CUDA_HOST_DEV real_type
                                  norm_square() const;

  /**
   * 使用[BOOST序列化库](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html)将此对象的数据读入或写入一个流中，以便进行序列化。
   *
   */
  template <class Archive>
  void
  serialize(Archive &ar, const unsigned int version);

  /**
   * 内部类型声明，用于对Tensor<1,dim,Number>的operator[]()的返回类型进行专业化。
   *
   */
  using tensor_type = Number;

private:
  /**
   * 这个标量对象的值。
   *
   */
  Number value;

  /**
   * 用于unroll的内部辅助函数。
   *
   */
  template <typename OtherNumber>
  void
  unroll_recursion(Vector<OtherNumber> &result,
                   unsigned int &       start_index) const;

  // Allow an arbitrary Tensor to access the underlying values.
  template <int, int, typename>
  friend class Tensor;
};



/**
 * 一个具有任意等级的一般张量类，即具有任意数量的索引。张量类提供了一个索引操作符和一些基础结构，但大多数功能是递归到等级为1的张量中，或者放到外部模板函数中，例如<tt>contract</tt>族。
 * 张量的等级规定了它可以代表哪些类型的物理量。  <ul>   <li>  一个等级为0的张量是一个标量，可以存储诸如温度或压力等量。这些标量在本文档中显示为简单的小写拉丁字母，例如：  $a, b, c, \dots$  。    </li>   <li>  等级1张量是一个有 @p dim 分量的矢量，它可以表示矢量，如速度、位移、电场等。它们也可以描述一个标量场的梯度。    秩-1张量使用的符号是粗体小写拉丁字母，例如  $\mathbf a, \mathbf b, \mathbf c, \dots$  。    一个秩-1张量的成分，如 $\mathbf a$ ，表示为 $a_i$ ，其中 $i$ 是0和<tt>dim-1</tt>之间的一个索引。    </li>   <li>  等级2张量是一种线性算子，可以将一个向量转化为另一个向量。这些张量类似于具有 $\text{dim} \times \text{dim}$ 成分的矩阵。有一个相关的类SymmetricTensor<2,dim>，用于等级2的张量，其元素是对称的。等级2张量通常用粗体的大写拉丁字母表示，如 $\mathbf A, \mathbf B, \dots$ 或粗体的希腊字母，如 $\boldsymbol{\varepsilon}, \boldsymbol{\sigma}$  。    等级2张量的组成部分如 $\mathbf A$ 用两个指数 $(i,j)$ 表示为 $A_{ij}$  。这些张量通常描述矢量场的梯度（变形梯度、速度梯度等）或标量场的Hessians。此外，机械应力张量是秩-2张量，将内部表面的单位法向量映射为局部牵引力（单位面积的力）向量。    </li>   <li>  等级大于2的张量也是以一致的方式定义的。它们有 $\text{dim}^{\text{rank}}$ 个分量，识别一个分量所需的指数数等于<tt>秩</tt>。对于等级为4的张量，存在一个称为SymmetricTensor<4,dim>的对称变体。    </li>   </ul> 。
 * 对秩为2的对象使用这个张量类，在许多情况下比矩阵更有优势，因为编译器知道维度以及数据的位置。因此有可能产生比具有与运行时间相关的维度的矩阵更有效的代码。这也使得代码更容易阅读，因为张量（一个与坐标系相关的对象，并且在坐标旋转和变换方面具有变换属性）和矩阵（我们认为它是与线性代数事物相关的任意向量空间上的运算符）之间存在语义上的差异。
 * @tparam  rank_
 * 一个整数，表示这个张量的等级。对于秩0张量，该类存在一个特殊化。
 * @tparam  dim
 * 一个整数，表示该张量所处空间的维度。当然，这等于识别一个点和秩-1张量的坐标数。
 * @tparam  Number
 * 用来存储张量元素的数据类型。在几乎所有的情况下，这只是默认的
 * @p double,
 * ，但在有些情况下，人们可能希望以不同的（而且总是标量的）类型存储元素。它可以用来将张量建立在
 * @p float 或 @p complex
 * 数字或任何其他实现基本算术操作的数据类型上。另一个例子是允许自动微分的类型（例如，见
 * step-33
 * 中使用的Sacado类型），从而可以生成一个以张量为参数的函数的分析（空间）导数。
 *
 *
 * @ingroup geomprimitives
 *
 */
template <int rank_, int dim, typename Number>
class Tensor
{
public:
  static_assert(rank_ >= 1,
                "Tensors must have a rank greater than or equal to one.");
  static_assert(dim >= 0,
                "Tensors must have a dimension greater than or equal to one.");
  /**
   * 提供一种方法来获取一个对象的尺寸，而不需要明确知道它的数据类型。实现这种方式而不是提供一个函数<tt>dimension()</tt>，因为现在有可能在编译时获得尺寸，而不需要内联函数的扩展和预评估；因此编译器可能产生更有效的代码，你可以使用这个值来声明其他数据类型。
   *
   */
  static constexpr unsigned int dimension = dim;

  /**
   * 向外界公布这个张量的等级。
   *
   */
  static constexpr unsigned int rank = rank_;

  /**
   * 当前等级的张量的独立成分的数量。这是dim乘以每个子张量的独立分量的数量。
   *
   */
  static constexpr unsigned int n_independent_components =
    Tensor<rank_ - 1, dim>::n_independent_components * dim;

  /**
   * 由这个容器封装的、由operator[]()返回的对象的类型。对于一般的张量来说，这是一个低等级的张量，对于Tensor<1,dim,Number>来说，这是一个标量数字类型。
   *
   */
  using value_type = typename Tensor<rank_ - 1, dim, Number>::tensor_type;

  /**
   * 声明一个数组类型，可以用来静态地初始化这个类型的对象。对于`dim
   * == 0`，它的大小是1，否则，它是`dim`。
   *
   */
  using array_type =
    typename Tensor<rank_ - 1, dim, Number>::array_type[(dim != 0) ? dim : 1];

  /**
   * 构造函数。将所有条目初始化为零。
   * @note 这个函数也可以在CUDA设备代码中使用。
   *
   */
  constexpr DEAL_II_ALWAYS_INLINE DEAL_II_CUDA_HOST_DEV
                                  Tensor();

  /**
   * 一个构造函数，数据从一个C风格的数组中复制出来。
   * @note  这个函数也可以在CUDA设备代码中使用。
   *
   */
  constexpr DEAL_II_CUDA_HOST_DEV explicit Tensor(
    const array_type &initializer);

  /**
   * 一个构造函数，数据从一个ArrayView对象中复制。
   * 显然，ArrayView对象必须代表一个大小为`dim`<sup>`rank`</sup>的数据延伸。参数`initializer`的顺序排列的元素被解释为unrolled_to_component_index()所描述的那样。
   * 这个构造函数显然要求 @p ElementType 类型等于 @p Number,
   * 或可转换为 @p Number. 数。
   * @note  这个函数也可以在CUDA设备代码中使用。
   *
   */
  template <typename ElementType, typename MemorySpace>
  constexpr DEAL_II_CUDA_HOST_DEV explicit Tensor(
    const ArrayView<ElementType, MemorySpace> &initializer);

  /**
   * 来自不同底层标量类型的张量的构造函数。这显然要求
   * @p OtherNumber 类型可转换为 @p 数字。
   * @note 这个函数也可以在CUDA设备代码中使用。
   *
   */
  template <typename OtherNumber>
  constexpr DEAL_II_CUDA_HOST_DEV
  Tensor(const Tensor<rank_, dim, OtherNumber> &initializer);

  /**
   * 可以从 "张量的张量 "转换的构造函数。
   *
   */
  template <typename OtherNumber>
  constexpr Tensor(
    const Tensor<1, dim, Tensor<rank_ - 1, dim, OtherNumber>> &initializer);

  /**
   * 到张量的张量的转换操作符。
   *
   */
  template <typename OtherNumber>
  constexpr
  operator Tensor<1, dim, Tensor<rank_ - 1, dim, OtherNumber>>() const;

#if __GNUC__ >= 11 || defined __INTEL_COMPILER
  /**
   * 复制构造函数
   *
   */
  constexpr Tensor(const Tensor<rank_, dim, Number> &);

  /**
   * 拷贝赋值运算符
   *
   */
  constexpr Tensor<rank_, dim, Number> &
  operator=(const Tensor<rank_, dim, Number> &);

  /**
   * 移动构造函数
   *
   */
  constexpr Tensor(Tensor<rank_, dim, Number> &&) noexcept;

  /**
   * 移动赋值运算符
   *
   */
  constexpr Tensor<rank_, dim, Number> &
  operator=(Tensor<rank_, dim, Number> &&) noexcept;
#endif

  /**
   * 读写访问操作符。
   * @note  这个函数也可以在CUDA设备代码中使用。
   *
   */
  constexpr DEAL_II_CUDA_HOST_DEV value_type &operator[](const unsigned int i);

  /**
   * 只读访问操作符。
   * @note  这个函数也可以在CUDA设备代码中使用。
   *
   */
  constexpr DEAL_II_CUDA_HOST_DEV const value_type &
                                        operator[](const unsigned int i) const;

  /**
   * 使用TableIndices <tt>indices</tt>进行读取访问
   *
   */
  constexpr const Number &operator[](const TableIndices<rank_> &indices) const;

  /**
   * 使用TableIndices<tt>indices</tt>进行读写访问
   *
   */
  constexpr Number &operator[](const TableIndices<rank_> &indices);

  /**
   * 返回一个指向底层存储的第一个元素的指针。
   *
   */
  Number *
  begin_raw();

  /**
   * 返回一个指向底层存储的第一个元素的常量指针。
   *
   */
  const Number *
  begin_raw() const;

  /**
   * 返回一个指向底层存储结束后的元素的指针。
   *
   */
  Number *
  end_raw();

  /**
   * 返回一个指向超过底层存储结束的元素的指针。
   *
   */
  const Number *
  end_raw() const;

  /**
   * 来自具有不同底层标量类型的张量的赋值运算符。
   * 这显然要求 @p OtherNumber 类型可以转换为 @p 数。
   * @note 这个函数也可以在CUDA设备代码中使用。
   *
   */
  template <typename OtherNumber>
  constexpr DEAL_II_CUDA_HOST_DEV Tensor &
                                  operator=(const Tensor<rank_, dim, OtherNumber> &rhs);

  /**
   * 这个操作符将一个标量分配给一个张量。为了避免混淆向张量分配标量值的确切含义，零是<tt>d</tt>唯一允许的值，允许用直观的符号<tt>t=0</tt>将张量的所有元素重置为零。
   *
   */
  constexpr Tensor &
  operator=(const Number &d);

  /**
   * 测试两个张量的相等。
   *
   */
  template <typename OtherNumber>
  constexpr bool
  operator==(const Tensor<rank_, dim, OtherNumber> &) const;

  /**
   * 测试两个张量的不平等。
   *
   */
  template <typename OtherNumber>
  constexpr bool
  operator!=(const Tensor<rank_, dim, OtherNumber> &) const;

  /**
   * 添加另一个张量。
   * @note 这个函数也可以在CUDA设备代码中使用。
   *
   */
  template <typename OtherNumber>
  constexpr DEAL_II_CUDA_HOST_DEV Tensor &
                                  operator+=(const Tensor<rank_, dim, OtherNumber> &);

  /**
   * 减去另一个张量。
   * @note  这个函数也可以在CUDA设备代码中使用。
   *
   */
  template <typename OtherNumber>
  constexpr DEAL_II_CUDA_HOST_DEV Tensor &
                                  operator-=(const Tensor<rank_, dim, OtherNumber> &);

  /**
   * 用<tt>因子</tt>缩放张量，即用<tt>因子</tt>乘以所有组件。
   * @note 这个函数也可以在CUDA设备代码中使用。
   *
   */
  template <typename OtherNumber>
  constexpr DEAL_II_CUDA_HOST_DEV Tensor &
                                  operator*=(const OtherNumber &factor);

  /**
   * 用<tt>1/factor</tt>缩放向量。
   * @note  这个函数也可以在CUDA设备代码中使用。
   *
   */
  template <typename OtherNumber>
  constexpr DEAL_II_CUDA_HOST_DEV Tensor &
                                  operator/=(const OtherNumber &factor);

  /**
   * 单元减法运算符。负掉张量的所有条目。
   * @note  这个函数也可以在CUDA设备代码中使用。
   *
   */
  constexpr DEAL_II_CUDA_HOST_DEV Tensor
                                  operator-() const;

  /**
   * 将所有值重置为零。    请注意，这与标准库容器的 @p
   * clear()成员函数以及deal.II中的其他几个类的语义部分不一致，它们不仅将存储元素的值重置为零，而且释放所有内存并将对象返回到处女状态。然而，由于本类型的对象的大小是由其模板参数决定的，所以调整大小是不可能的，事实上，所有元素的值都为零的状态就是这样一个对象构建后的状态。
   *
   */
  constexpr void
  clear();

  /**
   * 返回张量的Frobenius-norm，即所有条目的绝对平方之和的平方根。对于目前秩-1张量的情况，这等于通常的<tt>l<sub>2</sub></tt>向量的规范。
   * @note  这个函数也可以在CUDA设备代码中使用。
   *
   */
  DEAL_II_CUDA_HOST_DEV
  typename numbers::NumberTraits<Number>::real_type
  norm() const;

  /**
   * 返回张量的Frobenius-norm的平方，即所有条目的绝对平方之和。
   * @note  这个函数也可以在CUDA设备代码中使用。
   *
   */
  constexpr DEAL_II_CUDA_HOST_DEV
    typename numbers::NumberTraits<Number>::real_type
    norm_square() const;

  /**
   * 用所有张量元素填充一个向量。
   * 这个函数将所有的张量条目展开为一个单一的、线性编号的向量。正如C++中的惯例，张量的最右边的索引行进得最快。
   *
   */
  template <typename OtherNumber>
  void
  unroll(Vector<OtherNumber> &result) const;

  /**
   * 为函数的参数所索引的张量元素返回一个范围为
   * $[0,\text{dim}^{\text{rank}}-1]$ 的未滚动索引。
   *
   */
  static constexpr unsigned int
  component_to_unrolled_index(const TableIndices<rank_> &indices);

  /**
   * 与 component_to_unrolled_index 相反。对于 $[0,
   * \text{dim}^{\text{rank}}-1]$
   * 范围内的一个索引，返回它所对应的索引集。
   *
   */
  static constexpr TableIndices<rank_>
  unrolled_to_component_indices(const unsigned int i);

  /**
   * 确定这个对象的内存消耗（以字节为单位）的估计值。
   *
   */
  static constexpr std::size_t
  memory_consumption();

  /**
   * 使用[BOOST序列化库](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html)将此对象的数据读入或写入一个流中，以便进行序列化。
   *
   */
  template <class Archive>
  void
  serialize(Archive &ar, const unsigned int version);

  /**
   * 内部类型声明，用于对Tensor<1,dim,Number>的operator[]()的返回类型进行专业化。
   *
   */
  using tensor_type = Tensor<rank_, dim, Number>;

private:
  /**
   * 持有子元素的张量数组。
   *
   */
  Tensor<rank_ - 1, dim, Number> values[(dim != 0) ? dim : 1];
  // ... avoid a compiler warning in case of dim == 0 and ensure that the
  // array always has positive size.

  /**
   * 用于unroll的内部辅助函数。
   *
   */
  template <typename OtherNumber>
  void
  unroll_recursion(Vector<OtherNumber> &result,
                   unsigned int &       start_index) const;

  /**
   * 这个构造函数是供内部使用的。它提供了一种方法来创建Tensor<rank,
   * dim, Number>的constexpr构造函数。
   * @note 这个函数也可以在CUDA设备代码中使用。
   *
   */
  template <typename ArrayLike, std::size_t... Indices>
  constexpr DEAL_II_CUDA_HOST_DEV
  Tensor(const ArrayLike &initializer, std::index_sequence<Indices...>);

  // Allow an arbitrary Tensor to access the underlying values.
  template <int, int, typename>
  friend class Tensor;

  // Point is allowed access to the coordinates. This is supposed to improve
  // speed.
  friend class Point<dim, Number>;
};


#ifndef DOXYGEN
namespace internal
{
  // Workaround: The following 4 overloads are necessary to be able to
  // compile the library with Apple Clang 8 and older. We should remove
  // these overloads again when we bump the minimal required version to
  // something later than clang-3.6 / Apple Clang 6.3.
  template <int rank, int dim, typename T, typename U>
  struct ProductTypeImpl<Tensor<rank, dim, T>, std::complex<U>>
  {
    using type =
      Tensor<rank, dim, std::complex<typename ProductType<T, U>::type>>;
  };

  template <int rank, int dim, typename T, typename U>
  struct ProductTypeImpl<Tensor<rank, dim, std::complex<T>>, std::complex<U>>
  {
    using type =
      Tensor<rank, dim, std::complex<typename ProductType<T, U>::type>>;
  };

  template <typename T, int rank, int dim, typename U>
  struct ProductTypeImpl<std::complex<T>, Tensor<rank, dim, U>>
  {
    using type =
      Tensor<rank, dim, std::complex<typename ProductType<T, U>::type>>;
  };

  template <int rank, int dim, typename T, typename U>
  struct ProductTypeImpl<std::complex<T>, Tensor<rank, dim, std::complex<U>>>
  {
    using type =
      Tensor<rank, dim, std::complex<typename ProductType<T, U>::type>>;
  };
  // end workaround

  /**
   * 下面的结构需要用来初始化嵌套的Tensor对象。
   * 另外请看numbers.h中的另一个特殊化。
   *
   */
  template <int rank, int dim, typename T>
  struct NumberType<Tensor<rank, dim, T>>
  {
    static constexpr DEAL_II_ALWAYS_INLINE const Tensor<rank, dim, T> &
                                                 value(const Tensor<rank, dim, T> &t)
    {
      return t;
    }

    static constexpr DEAL_II_ALWAYS_INLINE Tensor<rank, dim, T>
                                           value(const T &t)
    {
      Tensor<rank, dim, T> tmp;
      tmp = t;
      return tmp;
    }
  };
} // namespace internal


 /*---------------------- Inline functions: Tensor<0,dim> ---------------------*/ 


template <int dim, typename Number>
constexpr DEAL_II_ALWAYS_INLINE DEAL_II_CUDA_HOST_DEV
                                Tensor<0, dim, Number>::Tensor()
  // Some auto-differentiable numbers need explicit
  // zero initialization such as adtl::adouble.
  : Tensor{0.0}
{}



template <int dim, typename Number>
template <typename OtherNumber>
constexpr DEAL_II_ALWAYS_INLINE DEAL_II_CUDA_HOST_DEV
                                Tensor<0, dim, Number>::Tensor(const OtherNumber &initializer)
  : value(internal::NumberType<Number>::value(initializer))
{}



template <int dim, typename Number>
template <typename OtherNumber>
constexpr DEAL_II_ALWAYS_INLINE DEAL_II_CUDA_HOST_DEV
                                Tensor<0, dim, Number>::Tensor(const Tensor<0, dim, OtherNumber> &p)
  : Tensor{p.value}
{}



#  if __GNUC__ >= 11 || defined __INTEL_COMPILER
template <int dim, typename Number>
constexpr DEAL_II_ALWAYS_INLINE DEAL_II_CUDA_HOST_DEV
                                Tensor<0, dim, Number>::Tensor(const Tensor<0, dim, Number> &other)
  : value{other.value}
{}



template <int dim, typename Number>
constexpr DEAL_II_ALWAYS_INLINE DEAL_II_CUDA_HOST_DEV Tensor<0, dim, Number> &
Tensor<0, dim, Number>::operator=(const Tensor<0, dim, Number> &other)
{
  value = other.value;
  return *this;
}



template <int dim, typename Number>
constexpr DEAL_II_ALWAYS_INLINE DEAL_II_CUDA_HOST_DEV
                                Tensor<0, dim, Number>::Tensor(Tensor<0, dim, Number> &&other) noexcept
  : value{std::move(other.value)}
{}



template <int dim, typename Number>
constexpr DEAL_II_ALWAYS_INLINE DEAL_II_CUDA_HOST_DEV Tensor<0, dim, Number> &
  Tensor<0, dim, Number>::operator=(Tensor<0, dim, Number> &&other) noexcept
{
  value = std::move(other.value);
  return *this;
}
#  endif



template <int dim, typename Number>
inline Number *
Tensor<0, dim, Number>::begin_raw()
{
  return std::addressof(value);
}



template <int dim, typename Number>
inline const Number *
Tensor<0, dim, Number>::begin_raw() const
{
  return std::addressof(value);
}



template <int dim, typename Number>
inline Number *
Tensor<0, dim, Number>::end_raw()
{
  return begin_raw() + n_independent_components;
}



template <int dim, typename Number>
const Number *
Tensor<0, dim, Number>::end_raw() const
{
  return begin_raw() + n_independent_components;
}



template <int dim, typename Number>
constexpr inline DEAL_II_ALWAYS_INLINE
  DEAL_II_CUDA_HOST_DEV Tensor<0, dim, Number>::operator Number &()
{
  // We cannot use Assert inside a CUDA kernel
#  ifndef __CUDA_ARCH__
  Assert(dim != 0,
         ExcMessage("Cannot access an object of type Tensor<0,0,Number>"));
#  endif
  return value;
}


template <int dim, typename Number>
constexpr inline DEAL_II_ALWAYS_INLINE
  DEAL_II_CUDA_HOST_DEV Tensor<0, dim, Number>::operator const Number &() const
{
  // We cannot use Assert inside a CUDA kernel
#  ifndef __CUDA_ARCH__
  Assert(dim != 0,
         ExcMessage("Cannot access an object of type Tensor<0,0,Number>"));
#  endif
  return value;
}


template <int dim, typename Number>
template <typename OtherNumber>
constexpr inline DEAL_II_ALWAYS_INLINE
  DEAL_II_CUDA_HOST_DEV Tensor<0, dim, Number> &
  Tensor<0, dim, Number>::operator=(const Tensor<0, dim, OtherNumber> &p)
{
  value = internal::NumberType<Number>::value(p);
  return *this;
}


template <int dim, typename Number>
template <typename OtherNumber>
constexpr inline DEAL_II_ALWAYS_INLINE
  DEAL_II_CUDA_HOST_DEV Tensor<0, dim, Number> &
  Tensor<0, dim, Number>::operator=(const OtherNumber &d)
{
  value = internal::NumberType<Number>::value(d);
  return *this;
}


template <int dim, typename Number>
template <typename OtherNumber>
constexpr inline bool
Tensor<0, dim, Number>::operator==(const Tensor<0, dim, OtherNumber> &p) const
{
#  if defined(DEAL_II_ADOLC_WITH_ADVANCED_BRANCHING)
  Assert(!(std::is_same<Number, adouble>::value ||
           std::is_same<OtherNumber, adouble>::value),
         ExcMessage(
           "The Tensor equality operator for ADOL-C taped numbers has not yet "
           "been extended to support advanced branching."));
#  endif

  return numbers::values_are_equal(value, p.value);
}


template <int dim, typename Number>
template <typename OtherNumber>
constexpr bool
Tensor<0, dim, Number>::operator!=(const Tensor<0, dim, OtherNumber> &p) const
{
  return !((*this) == p);
}


template <int dim, typename Number>
template <typename OtherNumber>
constexpr inline DEAL_II_ALWAYS_INLINE
  DEAL_II_CUDA_HOST_DEV Tensor<0, dim, Number> &
  Tensor<0, dim, Number>::operator+=(const Tensor<0, dim, OtherNumber> &p)
{
  value += p.value;
  return *this;
}


template <int dim, typename Number>
template <typename OtherNumber>
constexpr inline DEAL_II_ALWAYS_INLINE
  DEAL_II_CUDA_HOST_DEV Tensor<0, dim, Number> &
  Tensor<0, dim, Number>::operator-=(const Tensor<0, dim, OtherNumber> &p)
{
  value -= p.value;
  return *this;
}



namespace internal
{
  namespace ComplexWorkaround
  {
    template <typename Number, typename OtherNumber>
    constexpr inline DEAL_II_ALWAYS_INLINE DEAL_II_CUDA_HOST_DEV void
                                           multiply_assign_scalar(Number &val, const OtherNumber &s)
    {
      val *= s;
    }

#  ifdef __CUDA_ARCH__
    template <typename Number, typename OtherNumber>
    constexpr inline DEAL_II_ALWAYS_INLINE DEAL_II_CUDA_HOST_DEV void
                                           multiply_assign_scalar(std::complex<Number> &, const OtherNumber &)
    {
      printf("This function is not implemented for std::complex<Number>!\n");
      assert(false);
    }
#  endif
  } // namespace ComplexWorkaround
} // namespace internal


template <int dim, typename Number>
template <typename OtherNumber>
constexpr inline DEAL_II_ALWAYS_INLINE
  DEAL_II_CUDA_HOST_DEV Tensor<0, dim, Number> &
  Tensor<0, dim, Number>::operator*=(const OtherNumber &s)
{
  internal::ComplexWorkaround::multiply_assign_scalar(value, s);
  return *this;
}



template <int dim, typename Number>
template <typename OtherNumber>
constexpr inline DEAL_II_CUDA_HOST_DEV Tensor<0, dim, Number> &
Tensor<0, dim, Number>::operator/=(const OtherNumber &s)
{
  value /= s;
  return *this;
}


template <int dim, typename Number>
constexpr DEAL_II_ALWAYS_INLINE DEAL_II_CUDA_HOST_DEV Tensor<0, dim, Number>
Tensor<0, dim, Number>::operator-() const
{
  return -value;
}


template <int dim, typename Number>
inline DEAL_II_ALWAYS_INLINE typename Tensor<0, dim, Number>::real_type
Tensor<0, dim, Number>::norm() const
{
  Assert(dim != 0,
         ExcMessage("Cannot access an object of type Tensor<0,0,Number>"));
  return numbers::NumberTraits<Number>::abs(value);
}


template <int dim, typename Number>
constexpr DEAL_II_CUDA_HOST_DEV inline DEAL_II_ALWAYS_INLINE
  typename Tensor<0, dim, Number>::real_type
  Tensor<0, dim, Number>::norm_square() const
{
  // We cannot use Assert inside a CUDA kernel
#  ifndef __CUDA_ARCH__
  Assert(dim != 0,
         ExcMessage("Cannot access an object of type Tensor<0,0,Number>"));
#  endif
  return numbers::NumberTraits<Number>::abs_square(value);
}


template <int dim, typename Number>
template <typename OtherNumber>
inline void
Tensor<0, dim, Number>::unroll_recursion(Vector<OtherNumber> &result,
                                         unsigned int &       index) const
{
  Assert(dim != 0,
         ExcMessage("Cannot unroll an object of type Tensor<0,0,Number>"));
  result[index] = value;
  ++index;
}


template <int dim, typename Number>
constexpr inline void
Tensor<0, dim, Number>::clear()
{
  // Some auto-differentiable numbers need explicit
  // zero initialization.
  value = internal::NumberType<Number>::value(0.0);
}


template <int dim, typename Number>
template <class Archive>
inline void
Tensor<0, dim, Number>::serialize(Archive &ar, const unsigned int)
{
  ar &value;
}


template <int dim, typename Number>
constexpr unsigned int Tensor<0, dim, Number>::n_independent_components;


 /*-------------------- Inline functions: Tensor<rank,dim> --------------------*/ 

template <int rank_, int dim, typename Number>
template <typename ArrayLike, std::size_t... indices>
constexpr DEAL_II_ALWAYS_INLINE DEAL_II_CUDA_HOST_DEV
                                Tensor<rank_, dim, Number>::Tensor(const ArrayLike &initializer,
                                   std::index_sequence<indices...>)
  : values{Tensor<rank_ - 1, dim, Number>(initializer[indices])...}
{
  static_assert(sizeof...(indices) == dim,
                "dim should match the number of indices");
}



template <int rank_, int dim, typename Number>
constexpr DEAL_II_ALWAYS_INLINE DEAL_II_CUDA_HOST_DEV
                                Tensor<rank_, dim, Number>::Tensor()
  // We would like to use =default, but this causes compile errors with some
  // MSVC versions and internal compiler errors with -O1 in gcc 5.4.
  : values{}
{}



template <int rank_, int dim, typename Number>
constexpr DEAL_II_ALWAYS_INLINE DEAL_II_CUDA_HOST_DEV
                                Tensor<rank_, dim, Number>::Tensor(const array_type &initializer)
  : Tensor(initializer, std::make_index_sequence<dim>{})
{}



template <int rank_, int dim, typename Number>
template <typename ElementType, typename MemorySpace>
constexpr DEAL_II_ALWAYS_INLINE DEAL_II_CUDA_HOST_DEV
                                Tensor<rank_, dim, Number>::Tensor(
  const ArrayView<ElementType, MemorySpace> &initializer)
{
  AssertDimension(initializer.size(), n_independent_components);

  for (unsigned int i = 0; i < n_independent_components; ++i)
    (*this)[unrolled_to_component_indices(i)] = initializer[i];
}



template <int rank_, int dim, typename Number>
template <typename OtherNumber>
constexpr DEAL_II_ALWAYS_INLINE DEAL_II_CUDA_HOST_DEV
                                Tensor<rank_, dim, Number>::Tensor(
  const Tensor<rank_, dim, OtherNumber> &initializer)
  : Tensor(initializer, std::make_index_sequence<dim>{})
{}


template <int rank_, int dim, typename Number>
template <typename OtherNumber>
constexpr DEAL_II_ALWAYS_INLINE
Tensor<rank_, dim, Number>::Tensor(
  const Tensor<1, dim, Tensor<rank_ - 1, dim, OtherNumber>> &initializer)
  : Tensor(initializer, std::make_index_sequence<dim>{})
{}


template <int rank_, int dim, typename Number>
template <typename OtherNumber>
constexpr DEAL_II_ALWAYS_INLINE Tensor<rank_, dim, Number>::
                                operator Tensor<1, dim, Tensor<rank_ - 1, dim, OtherNumber>>() const
{
  return Tensor<1, dim, Tensor<rank_ - 1, dim, Number>>(values);
}


#  if __GNUC__ >= 11 || defined __INTEL_COMPILER
template <int rank_, int dim, typename Number>
constexpr DEAL_II_ALWAYS_INLINE
Tensor<rank_, dim, Number>::Tensor(const Tensor<rank_, dim, Number> &other)
{
  for (unsigned int i = 0; i < dim; ++i)
    values[i] = other.values[i];
}


template <int rank_, int dim, typename Number>
constexpr DEAL_II_ALWAYS_INLINE Tensor<rank_, dim, Number> &
Tensor<rank_, dim, Number>::operator=(const Tensor<rank_, dim, Number> &other)
{
  for (unsigned int i = 0; i < dim; ++i)
    values[i] = other.values[i];
  return *this;
}


template <int rank_, int dim, typename Number>
constexpr DEAL_II_ALWAYS_INLINE
Tensor<rank_, dim, Number>::Tensor(Tensor<rank_, dim, Number> &&other) noexcept
{
  for (unsigned int i = 0; i < dim; ++i)
    values[i] = other.values[i];
}


template <int rank_, int dim, typename Number>
constexpr DEAL_II_ALWAYS_INLINE Tensor<rank_, dim, Number> &
                                Tensor<rank_, dim, Number>::
                                operator=(Tensor<rank_, dim, Number> &&other) noexcept
{
  for (unsigned int i = 0; i < dim; ++i)
    values[i] = other.values[i];
  return *this;
}
#  endif


namespace internal
{
  namespace TensorSubscriptor
  {
    template <typename ArrayElementType, int dim>
    constexpr inline DEAL_II_ALWAYS_INLINE
      DEAL_II_CUDA_HOST_DEV ArrayElementType &
                            subscript(ArrayElementType * values,
                                      const unsigned int i,
                                      std::integral_constant<int, dim>)
    {
      // We cannot use Assert in a CUDA kernel
#  ifndef __CUDA_ARCH__
      AssertIndexRange(i, dim);
#  endif
      return values[i];
    }

    // The variables within this struct will be referenced in the next function.
    // It is a workaround that allows returning a reference to a static variable
    // while allowing constexpr evaluation of the function.
    // It has to be defined outside the function because constexpr functions
    // cannot define static variables
    template <typename ArrayElementType>
    struct Uninitialized
    {
      static ArrayElementType value;
    };

    template <typename Type>
    Type Uninitialized<Type>::value;

    template <typename ArrayElementType>
    constexpr inline DEAL_II_ALWAYS_INLINE
      DEAL_II_CUDA_HOST_DEV ArrayElementType &
                            subscript(ArrayElementType *,
                                      const unsigned int,
                                      std::integral_constant<int, 0>)
    {
      // We cannot use Assert in a CUDA kernel
#  ifndef __CUDA_ARCH__
      Assert(
        false,
        ExcMessage(
          "Cannot access elements of an object of type Tensor<rank,0,Number>."));
#  endif
      return Uninitialized<ArrayElementType>::value;
    }
  } // namespace TensorSubscriptor
} // namespace internal


template <int rank_, int dim, typename Number>
constexpr inline DEAL_II_ALWAYS_INLINE             DEAL_II_CUDA_HOST_DEV
  typename Tensor<rank_, dim, Number>::value_type &Tensor<rank_, dim, Number>::
                                                   operator[](const unsigned int i)
{
  return dealii::internal::TensorSubscriptor::subscript(
    values, i, std::integral_constant<int, dim>());
}


template <int rank_, int dim, typename Number>
constexpr DEAL_II_ALWAYS_INLINE
    DEAL_II_CUDA_HOST_DEV const typename Tensor<rank_, dim, Number>::value_type &
    Tensor<rank_, dim, Number>::operator[](const unsigned int i) const
{
#  ifndef DEAL_II_COMPILER_CUDA_AWARE
  AssertIndexRange(i, dim);
#  endif

  return values[i];
}


template <int rank_, int dim, typename Number>
constexpr inline DEAL_II_ALWAYS_INLINE const Number &
                                             Tensor<rank_, dim, Number>::
                                             operator[](const TableIndices<rank_> &indices) const
{
#  ifndef DEAL_II_COMPILER_CUDA_AWARE
  Assert(dim != 0,
         ExcMessage("Cannot access an object of type Tensor<rank_,0,Number>"));
#  endif

  return TensorAccessors::extract<rank_>(*this, indices);
}



template <int rank_, int dim, typename Number>
constexpr inline DEAL_II_ALWAYS_INLINE Number &Tensor<rank_, dim, Number>::
                                               operator[](const TableIndices<rank_> &indices)
{
#  ifndef DEAL_II_COMPILER_CUDA_AWARE
  Assert(dim != 0,
         ExcMessage("Cannot access an object of type Tensor<rank_,0,Number>"));
#  endif

  return TensorAccessors::extract<rank_>(*this, indices);
}



template <int rank_, int dim, typename Number>
inline Number *
Tensor<rank_, dim, Number>::begin_raw()
{
  return std::addressof(
    this->operator[](this->unrolled_to_component_indices(0)));
}



template <int rank_, int dim, typename Number>
inline const Number *
Tensor<rank_, dim, Number>::begin_raw() const
{
  return std::addressof(
    this->operator[](this->unrolled_to_component_indices(0)));
}



template <int rank_, int dim, typename Number>
inline Number *
Tensor<rank_, dim, Number>::end_raw()
{
  return begin_raw() + n_independent_components;
}



template <int rank_, int dim, typename Number>
inline const Number *
Tensor<rank_, dim, Number>::end_raw() const
{
  return begin_raw() + n_independent_components;
}



template <int rank_, int dim, typename Number>
template <typename OtherNumber>
constexpr inline DEAL_II_ALWAYS_INLINE Tensor<rank_, dim, Number> &
Tensor<rank_, dim, Number>::operator=(const Tensor<rank_, dim, OtherNumber> &t)
{
  // The following loop could be written more concisely using std::copy, but
  // that function is only constexpr from C++20 on.
  for (unsigned int i = 0; i < dim; ++i)
    values[i] = t.values[i];
  return *this;
}


template <int rank_, int dim, typename Number>
constexpr inline DEAL_II_ALWAYS_INLINE Tensor<rank_, dim, Number> &
Tensor<rank_, dim, Number>::operator=(const Number &d)
{
  Assert(numbers::value_is_zero(d), ExcScalarAssignmentOnlyForZeroValue());
  (void)d;

  for (unsigned int i = 0; i < dim; ++i)
    values[i] = internal::NumberType<Number>::value(0.0);
  return *this;
}


template <int rank_, int dim, typename Number>
template <typename OtherNumber>
constexpr inline bool
Tensor<rank_, dim, Number>::
operator==(const Tensor<rank_, dim, OtherNumber> &p) const
{
  for (unsigned int i = 0; i < dim; ++i)
    if (values[i] != p.values[i])
      return false;
  return true;
}


// At some places in the library, we have Point<0> for formal reasons
// (e.g., we sometimes have Quadrature<dim-1> for faces, so we have
// Quadrature<0> for dim=1, and then we have Point<0>). To avoid warnings
// in the above function that the loop end check always fails, we
// implement this function here
template <>
template <>
constexpr inline bool
Tensor<1, 0, double>::operator==(const Tensor<1, 0, double> &) const
{
  return true;
}


template <int rank_, int dim, typename Number>
template <typename OtherNumber>
constexpr bool
Tensor<rank_, dim, Number>::
operator!=(const Tensor<rank_, dim, OtherNumber> &p) const
{
  return !((*this) == p);
}


template <int rank_, int dim, typename Number>
template <typename OtherNumber>
constexpr inline DEAL_II_ALWAYS_INLINE
  DEAL_II_CUDA_HOST_DEV Tensor<rank_, dim, Number> &
                        Tensor<rank_, dim, Number>::
                        operator+=(const Tensor<rank_, dim, OtherNumber> &p)
{
  for (unsigned int i = 0; i < dim; ++i)
    values[i] += p.values[i];
  return *this;
}


template <int rank_, int dim, typename Number>
template <typename OtherNumber>
constexpr inline DEAL_II_ALWAYS_INLINE
  DEAL_II_CUDA_HOST_DEV Tensor<rank_, dim, Number> &
                        Tensor<rank_, dim, Number>::
                        operator-=(const Tensor<rank_, dim, OtherNumber> &p)
{
  for (unsigned int i = 0; i < dim; ++i)
    values[i] -= p.values[i];
  return *this;
}


template <int rank_, int dim, typename Number>
template <typename OtherNumber>
constexpr inline DEAL_II_ALWAYS_INLINE
  DEAL_II_CUDA_HOST_DEV Tensor<rank_, dim, Number> &
  Tensor<rank_, dim, Number>::operator*=(const OtherNumber &s)
{
  for (unsigned int i = 0; i < dim; ++i)
    values[i] *= s;
  return *this;
}


namespace internal
{
  namespace TensorImplementation
  {
    template <int rank,
              int dim,
              typename Number,
              typename OtherNumber,
              typename std::enable_if<
                !std::is_integral<
                  typename ProductType<Number, OtherNumber>::type>::value &&
                  !std::is_same<Number, Differentiation::SD::Expression>::value,
                int>::type = 0>
    constexpr DEAL_II_CUDA_HOST_DEV inline DEAL_II_ALWAYS_INLINE void
    division_operator(Tensor<rank, dim, Number> (&t)[dim],
                      const OtherNumber &factor)
    {
      const Number inverse_factor = Number(1.) / factor;
      // recurse over the base objects
      for (unsigned int d = 0; d < dim; ++d)
        t[d] *= inverse_factor;
    }


    template <int rank,
              int dim,
              typename Number,
              typename OtherNumber,
              typename std::enable_if<
                std::is_integral<
                  typename ProductType<Number, OtherNumber>::type>::value ||
                  std::is_same<Number, Differentiation::SD::Expression>::value,
                int>::type = 0>
    constexpr DEAL_II_CUDA_HOST_DEV inline DEAL_II_ALWAYS_INLINE void
    division_operator(Tensor<rank, dim, Number> (&t)[dim],
                      const OtherNumber &factor)
    {
      // recurse over the base objects
      for (unsigned int d = 0; d < dim; ++d)
        t[d] /= factor;
    }
  } // namespace TensorImplementation
} // namespace internal


template <int rank_, int dim, typename Number>
template <typename OtherNumber>
constexpr inline DEAL_II_ALWAYS_INLINE
  DEAL_II_CUDA_HOST_DEV Tensor<rank_, dim, Number> &
  Tensor<rank_, dim, Number>::operator/=(const OtherNumber &s)
{
  internal::TensorImplementation::division_operator(values, s);
  return *this;
}


template <int rank_, int dim, typename Number>
constexpr inline DEAL_II_ALWAYS_INLINE
  DEAL_II_CUDA_HOST_DEV Tensor<rank_, dim, Number>
  Tensor<rank_, dim, Number>::operator-() const
{
  Tensor<rank_, dim, Number> tmp;

  for (unsigned int i = 0; i < dim; ++i)
    tmp.values[i] = -values[i];

  return tmp;
}


template <int rank_, int dim, typename Number>
inline typename numbers::NumberTraits<Number>::real_type
Tensor<rank_, dim, Number>::norm() const
{
  return std::sqrt(norm_square());
}


template <int rank_, int dim, typename Number>
constexpr inline DEAL_II_ALWAYS_INLINE DEAL_II_CUDA_HOST_DEV
  typename numbers::NumberTraits<Number>::real_type
  Tensor<rank_, dim, Number>::norm_square() const
{
  typename numbers::NumberTraits<Number>::real_type s = internal::NumberType<
    typename numbers::NumberTraits<Number>::real_type>::value(0.0);
  for (unsigned int i = 0; i < dim; ++i)
    s += values[i].norm_square();

  return s;
}


template <int rank_, int dim, typename Number>
template <typename OtherNumber>
inline void
Tensor<rank_, dim, Number>::unroll(Vector<OtherNumber> &result) const
{
  AssertDimension(result.size(),
                  (Utilities::fixed_power<rank_, unsigned int>(dim)));

  unsigned int index = 0;
  unroll_recursion(result, index);
}


template <int rank_, int dim, typename Number>
template <typename OtherNumber>
inline void
Tensor<rank_, dim, Number>::unroll_recursion(Vector<OtherNumber> &result,
                                             unsigned int &       index) const
{
  for (unsigned int i = 0; i < dim; ++i)
    values[i].unroll_recursion(result, index);
}


template <int rank_, int dim, typename Number>
constexpr inline unsigned int
Tensor<rank_, dim, Number>::component_to_unrolled_index(
  const TableIndices<rank_> &indices)
{
  unsigned int index = 0;
  for (int r = 0; r < rank_; ++r)
    index = index * dim + indices[r];

  return index;
}



namespace internal
{
  // unrolled_to_component_indices is instantiated from DataOut for dim==0
  // and rank=2. Make sure we don't have compiler warnings.

  template <int dim>
  inline constexpr unsigned int
  mod(const unsigned int x)
  {
    return x % dim;
  }

  template <>
  inline unsigned int
  mod<0>(const unsigned int x)
  {
    Assert(false, ExcInternalError());
    return x;
  }

  template <int dim>
  inline constexpr unsigned int
  div(const unsigned int x)
  {
    return x / dim;
  }

  template <>
  inline unsigned int
  div<0>(const unsigned int x)
  {
    Assert(false, ExcInternalError());
    return x;
  }

} // namespace internal



template <int rank_, int dim, typename Number>
constexpr inline TableIndices<rank_>
Tensor<rank_, dim, Number>::unrolled_to_component_indices(const unsigned int i)
{
  AssertIndexRange(i, n_independent_components);

  TableIndices<rank_> indices;

  unsigned int remainder = i;
  for (int r = rank_ - 1; r >= 0; --r)
    {
      indices[r] = internal::mod<dim>(remainder);
      remainder  = internal::div<dim>(remainder);
    }
  Assert(remainder == 0, ExcInternalError());

  return indices;
}


template <int rank_, int dim, typename Number>
constexpr inline void
Tensor<rank_, dim, Number>::clear()
{
  for (unsigned int i = 0; i < dim; ++i)
    values[i] = internal::NumberType<Number>::value(0.0);
}


template <int rank_, int dim, typename Number>
constexpr std::size_t
Tensor<rank_, dim, Number>::memory_consumption()
{
  return sizeof(Tensor<rank_, dim, Number>);
}


template <int rank_, int dim, typename Number>
template <class Archive>
inline void
Tensor<rank_, dim, Number>::serialize(Archive &ar, const unsigned int)
{
  ar &values;
}


template <int rank_, int dim, typename Number>
constexpr unsigned int Tensor<rank_, dim, Number>::n_independent_components;

#endif // DOXYGEN

 /* ----------------- Non-member functions operating on tensors. ------------ */ 

/**
 * @name  Tensor对象的输出函数
 *
 *
 */
//@{

/**
 * 张量的输出运算符。连续打印元素，中间有一个空格，等级1的子张量之间有两个空格，等级2之间有三个空格，以此类推。
 * @relatesalso  Tensor
 *
 *
 */
template <int rank_, int dim, typename Number>
inline std::ostream &
operator<<(std::ostream &out, const Tensor<rank_, dim, Number> &p)
{
  for (unsigned int i = 0; i < dim; ++i)
    {
      out << p[i];
      if (i != dim - 1)
        out << ' ';
    }

  return out;
}


/**
 * 秩为0的张量的输出算子。由于这种张量是标量，我们只需打印这一个值。
 * @relatesalso  Tensor
 *
 *
 */
template <int dim, typename Number>
inline std::ostream &
operator<<(std::ostream &out, const Tensor<0, dim, Number> &p)
{
  out << static_cast<const Number &>(p);
  return out;
}


//@}
/**
 * @name  张量对象的矢量空间操作。
 *
 *
 */
//@{


/**
 * 秩为0的张量与左边的对象进行标量乘法。
 * 这个函数将存储在张量中的底层 @p Number
 * 解开，并与之相乘 @p object 。
 *
 *
 * @note  这个函数也可以在CUDA设备代码中使用。
 * @relatesalso  Tensor
 *
 *
 */
template <int dim, typename Number, typename Other>
constexpr DEAL_II_CUDA_HOST_DEV inline DEAL_II_ALWAYS_INLINE
  typename ProductType<Other, Number>::type
  operator*(const Other &object, const Tensor<0, dim, Number> &t)
{
  return object * static_cast<const Number &>(t);
}



/**
 * 秩为0的张量与来自右边的对象进行标量乘法。
 * 这个函数将存储在张量中的底层 @p Number
 * 解开，并与之相乘 @p object 。
 *
 *
 * @note  这个函数也可以在CUDA设备代码中使用。
 * @relatesalso  Tensor
 *
 *
 */
template <int dim, typename Number, typename Other>
constexpr DEAL_II_CUDA_HOST_DEV inline DEAL_II_ALWAYS_INLINE
  typename ProductType<Number, Other>::type
  operator*(const Tensor<0, dim, Number> &t, const Other &object)
{
  return static_cast<const Number &>(t) * object;
}


/**
 * 两个等级为0的张量的标量乘法。
 * 这个函数将存储在张量中的 @p Number 和 @p
 * OtherNumber类型的底层对象解包，并将它们相乘。它返回一个解包的乘积类型的数字。
 *
 *
 * @note  这个函数也可以在CUDA设备代码中使用。
 * @relatesalso  Tensor
 *
 *
 */
template <int dim, typename Number, typename OtherNumber>
DEAL_II_CUDA_HOST_DEV constexpr DEAL_II_ALWAYS_INLINE
  typename ProductType<Number, OtherNumber>::type
  operator*(const Tensor<0, dim, Number> &     src1,
            const Tensor<0, dim, OtherNumber> &src2)
{
  return static_cast<const Number &>(src1) *
         static_cast<const OtherNumber &>(src2);
}


/**
 * 等级为0的张量除以一个标量数字。
 *
 *
 * @note  这个函数也可以在CUDA设备代码中使用。
 * @relatesalso  Tensor
 *
 *
 */
template <int dim, typename Number, typename OtherNumber>
DEAL_II_CUDA_HOST_DEV constexpr DEAL_II_ALWAYS_INLINE
  Tensor<0,
         dim,
         typename ProductType<Number,
                              typename EnableIfScalar<OtherNumber>::type>::type>
  operator/(const Tensor<0, dim, Number> &t, const OtherNumber &factor)
{
  return static_cast<const Number &>(t) / factor;
}


/**
 * 添加两个等级为0的张量。
 *
 *
 * @note  这个函数也可以在CUDA设备代码中使用。
 * @relatesalso  Tensor
 *
 *
 */
template <int dim, typename Number, typename OtherNumber>
constexpr DEAL_II_ALWAYS_INLINE DEAL_II_CUDA_HOST_DEV
                                Tensor<0, dim, typename ProductType<Number, OtherNumber>::type>
                                operator+(const Tensor<0, dim, Number> &     p,
            const Tensor<0, dim, OtherNumber> &q)
{
  return static_cast<const Number &>(p) + static_cast<const OtherNumber &>(q);
}


/**
 * 减去两个等级为0的张量。
 *
 *
 * @note  这个函数也可以在CUDA设备代码中使用。
 * @relatesalso  Tensor
 *
 *
 */
template <int dim, typename Number, typename OtherNumber>
constexpr DEAL_II_ALWAYS_INLINE DEAL_II_CUDA_HOST_DEV
                                Tensor<0, dim, typename ProductType<Number, OtherNumber>::type>
                                operator-(const Tensor<0, dim, Number> &     p,
            const Tensor<0, dim, OtherNumber> &q)
{
  return static_cast<const Number &>(p) - static_cast<const OtherNumber &>(q);
}


/**
 * 一般等级的张量与来自右边的标量的乘法。
 * 只允许与标量数类型（即浮点数、复数浮点数等）的乘法，详情请参见EnableIfScalar的文档。
 *
 *
 * @note  这个函数也可以在CUDA设备代码中使用。
 * @relatesalso  Tensor
 *
 *
 */
template <int rank, int dim, typename Number, typename OtherNumber>
constexpr DEAL_II_CUDA_HOST_DEV inline DEAL_II_ALWAYS_INLINE
  Tensor<rank,
         dim,
         typename ProductType<Number,
                              typename EnableIfScalar<OtherNumber>::type>::type>
  operator*(const Tensor<rank, dim, Number> &t, const OtherNumber &factor)
{
  // recurse over the base objects
  Tensor<rank, dim, typename ProductType<Number, OtherNumber>::type> tt;
  for (unsigned int d = 0; d < dim; ++d)
    tt[d] = t[d] * factor;
  return tt;
}


/**
 * 一般等级的张量与左边的标量数字的乘法。
 * 只允许与标量数类型（即浮点数、复数浮点数等）的乘法，详情请参见EnableIfScalar的文档。
 *
 *
 * @note  这个函数也可以在CUDA设备代码中使用。
 * @relatesalso  Tensor
 *
 *
 */
template <int rank, int dim, typename Number, typename OtherNumber>
DEAL_II_CUDA_HOST_DEV constexpr inline DEAL_II_ALWAYS_INLINE
  Tensor<rank,
         dim,
         typename ProductType<typename EnableIfScalar<Number>::type,
                              OtherNumber>::type>
  operator*(const Number &factor, const Tensor<rank, dim, OtherNumber> &t)
{
  // simply forward to the operator above
  return t * factor;
}


namespace internal
{
  namespace TensorImplementation
  {
    template <int rank,
              int dim,
              typename Number,
              typename OtherNumber,
              typename std::enable_if<
                !std::is_integral<
                  typename ProductType<Number, OtherNumber>::type>::value,
                int>::type = 0>
    constexpr DEAL_II_CUDA_HOST_DEV inline DEAL_II_ALWAYS_INLINE
      Tensor<rank, dim, typename ProductType<Number, OtherNumber>::type>
      division_operator(const Tensor<rank, dim, Number> &t,
                        const OtherNumber &              factor)
    {
      Tensor<rank, dim, typename ProductType<Number, OtherNumber>::type> tt;
      const Number inverse_factor = Number(1.) / factor;
      // recurse over the base objects
      for (unsigned int d = 0; d < dim; ++d)
        tt[d] = t[d] * inverse_factor;
      return tt;
    }


    template <int rank,
              int dim,
              typename Number,
              typename OtherNumber,
              typename std::enable_if<
                std::is_integral<
                  typename ProductType<Number, OtherNumber>::type>::value,
                int>::type = 0>
    constexpr DEAL_II_CUDA_HOST_DEV inline DEAL_II_ALWAYS_INLINE
      Tensor<rank, dim, typename ProductType<Number, OtherNumber>::type>
      division_operator(const Tensor<rank, dim, Number> &t,
                        const OtherNumber &              factor)
    {
      Tensor<rank, dim, typename ProductType<Number, OtherNumber>::type> tt;
      // recurse over the base objects
      for (unsigned int d = 0; d < dim; ++d)
        tt[d] = t[d] / factor;
      return tt;
    }
  } // namespace TensorImplementation
} // namespace internal


/**
 * 一般等级的张量与标量数字的除法。关于模板参数和返回类型的更多信息，请参见上面关于operator*()的讨论。
 *
 *
 * @note  这个函数也可以在CUDA设备代码中使用。
 * @relatesalso  Tensor
 *
 *
 */
template <int rank, int dim, typename Number, typename OtherNumber>
constexpr DEAL_II_CUDA_HOST_DEV inline DEAL_II_ALWAYS_INLINE
  Tensor<rank,
         dim,
         typename ProductType<Number,
                              typename EnableIfScalar<OtherNumber>::type>::type>
  operator/(const Tensor<rank, dim, Number> &t, const OtherNumber &factor)
{
  return internal::TensorImplementation::division_operator(t, factor);
}


/**
 * 两个一般等级的张量的相加。
 * @tparam  rank 两个张量的等级。
 *
 *
 * @note  这个函数也可以在CUDA设备代码中使用。
 * @relatesalso  Tensor
 *
 *
 */
template <int rank, int dim, typename Number, typename OtherNumber>
constexpr DEAL_II_CUDA_HOST_DEV inline DEAL_II_ALWAYS_INLINE
  Tensor<rank, dim, typename ProductType<Number, OtherNumber>::type>
  operator+(const Tensor<rank, dim, Number> &     p,
            const Tensor<rank, dim, OtherNumber> &q)
{
  Tensor<rank, dim, typename ProductType<Number, OtherNumber>::type> tmp(p);

  for (unsigned int i = 0; i < dim; ++i)
    tmp[i] += q[i];

  return tmp;
}


/**
 * 两个一般等级的张量的减法。
 * @tparam  rank 两个张量的等级。
 *
 *
 * @note  这个函数也可以在CUDA设备代码中使用。
 * @relatesalso  Tensor
 *
 *
 */
template <int rank, int dim, typename Number, typename OtherNumber>
constexpr DEAL_II_CUDA_HOST_DEV inline DEAL_II_ALWAYS_INLINE
  Tensor<rank, dim, typename ProductType<Number, OtherNumber>::type>
  operator-(const Tensor<rank, dim, Number> &     p,
            const Tensor<rank, dim, OtherNumber> &q)
{
  Tensor<rank, dim, typename ProductType<Number, OtherNumber>::type> tmp(p);

  for (unsigned int i = 0; i < dim; ++i)
    tmp[i] -= q[i];

  return tmp;
}

/**
 * 两个等级为0的张量对象的进位乘法（即两个标量值的乘法）。
 * @relatesalso Tensor
 *
 *
 */
template <int dim, typename Number, typename OtherNumber>
inline constexpr DEAL_II_ALWAYS_INLINE
  Tensor<0, dim, typename ProductType<Number, OtherNumber>::type>
  schur_product(const Tensor<0, dim, Number> &     src1,
                const Tensor<0, dim, OtherNumber> &src2)
{
  Tensor<0, dim, typename ProductType<Number, OtherNumber>::type> tmp(src1);

  tmp *= src2;

  return tmp;
}

/**
 * 两个一般等级的张量对象的入境乘法。
 * 这种乘法也被称为 "Hadamard-product"（参看https://en.wikipedia.org/wiki/Hadamard_product_(matrices)），并产生一个大小为<rank, dim>的新张量。@f[
 * \text{result}_{i, j}
 * = \text{left}_{i, j}\circ
 *   \text{right}_{i, j}
 * @f]
 * @tparam  rank 两个张量的等级。
 * @relatesalso  Tensor
 *
 *
 */
template <int rank, int dim, typename Number, typename OtherNumber>
inline constexpr DEAL_II_ALWAYS_INLINE
  Tensor<rank, dim, typename ProductType<Number, OtherNumber>::type>
  schur_product(const Tensor<rank, dim, Number> &     src1,
                const Tensor<rank, dim, OtherNumber> &src2)
{
  Tensor<rank, dim, typename ProductType<Number, OtherNumber>::type> tmp;

  for (unsigned int i = 0; i < dim; ++i)
    tmp[i] = schur_product(Tensor<rank - 1, dim, Number>(src1[i]),
                           Tensor<rank - 1, dim, OtherNumber>(src2[i]));

  return tmp;
}

//@}
/**
 * @name  张量对象的收缩操作和外积
 *
 */
//@{


/**
 * 张量的点积（单一收缩）。返回一个等级为 $(\text{rank}_1 +
 * \text{rank}_2
 *
 * - 2)$ 的张量，它是等级为 @p rank_1 的张量 @p src1 的最后一个索引与等级为 @p rank_2: 的张量@f[
 * \text{result}_{i_1,\ldots,i_{r1},j_1,\ldots,j_{r2}}
 * = \sum_{k}
 *   \text{left}_{i_1,\ldots,i_{r1}, k}
 *   \text{right}_{k, j_1,\ldots,j_{r2}}
 * @f]的第一个索引的收缩。
 *
 *
 * @note
 * 对于张量类，乘法运算符只对一对指数进行收缩。这与SymmetricTensor的乘法运算符相反，后者做的是双重收缩。
 *
 *
 * @note
 * 如果收缩产生一个等级为0的张量，那么标量数将作为一个未包装的数字类型返回。
 * @relatesalso  Tensor
 *
 *
 */
template <int rank_1,
          int rank_2,
          int dim,
          typename Number,
          typename OtherNumber,
          typename = typename std::enable_if<rank_1 >= 1 && rank_2 >= 1>::type>
constexpr inline DEAL_II_ALWAYS_INLINE
  typename Tensor<rank_1 + rank_2 - 2,
                  dim,
                  typename ProductType<Number, OtherNumber>::type>::tensor_type
  operator*(const Tensor<rank_1, dim, Number> &     src1,
            const Tensor<rank_2, dim, OtherNumber> &src2)
{
  typename Tensor<rank_1 + rank_2 - 2,
                  dim,
                  typename ProductType<Number, OtherNumber>::type>::tensor_type
    result{};

  TensorAccessors::internal::
    ReorderedIndexView<0, rank_2, const Tensor<rank_2, dim, OtherNumber>>
      reordered = TensorAccessors::reordered_index_view<0, rank_2>(src2);
  TensorAccessors::contract<1, rank_1, rank_2, dim>(result, src1, reordered);

  return result;
}


/**
 * 对任意等级的两个张量的一对索引进行通用收缩。返回一个等级为
 * $(\text{rank}_1 + \text{rank}_2
 *
 * - 2)$ 的张量，它是等级为 @p rank_1 的张量 @p src1 的索引 @p index_2 与等级为 @p rank_2: 的张量@f[
 * \text{result}_{i_1,\ldots,i_{r1},j_1,\ldots,j_{r2}}
 * = \sum_{k}
 *   \text{left}_{i_1,\ldots,k,\ldots,i_{r1}}
 *   \text{right}_{j_1,\ldots,k,\ldots,j_{r2}}
 * @f]之间的收缩。
 * 例如，如果张量 <code>index_1==0</code> 的第一个索引（
 * <code>index_1==0</code> ）应与第三个索引（
 * <code>index_2==2</code>) of a tensor <code>t2</code>
 * ）收缩，该函数应被调用为
 *
 * @code
 * contract<0, 2>(t1, t2);
 * @endcode
 *
 *
 *
 * @note  索引的位置是从0开始计算的，即
 * $0\le\text{index}_i<\text{range}_i$  。
 *
 *
 * @note
 * 如果收缩产生一个等级为0的张量，那么标量数将作为一个未包装的数字类型返回。
 * @relatesalso  Tensor
 *
 *
 */
template <int index_1,
          int index_2,
          int rank_1,
          int rank_2,
          int dim,
          typename Number,
          typename OtherNumber>
constexpr inline DEAL_II_ALWAYS_INLINE
  typename Tensor<rank_1 + rank_2 - 2,
                  dim,
                  typename ProductType<Number, OtherNumber>::type>::tensor_type
  contract(const Tensor<rank_1, dim, Number> &     src1,
           const Tensor<rank_2, dim, OtherNumber> &src2)
{
  Assert(0 <= index_1 && index_1 < rank_1,
         ExcMessage(
           "The specified index_1 must lie within the range [0,rank_1)"));
  Assert(0 <= index_2 && index_2 < rank_2,
         ExcMessage(
           "The specified index_2 must lie within the range [0,rank_2)"));

  using namespace TensorAccessors;
  using namespace TensorAccessors::internal;

  // Reorder index_1 to the end of src1:
  const ReorderedIndexView<index_1, rank_1, const Tensor<rank_1, dim, Number>>
    reord_01 = reordered_index_view<index_1, rank_1>(src1);

  // Reorder index_2 to the end of src2:
  const ReorderedIndexView<index_2,
                           rank_2,
                           const Tensor<rank_2, dim, OtherNumber>>
    reord_02 = reordered_index_view<index_2, rank_2>(src2);

  typename Tensor<rank_1 + rank_2 - 2,
                  dim,
                  typename ProductType<Number, OtherNumber>::type>::tensor_type
    result{};
  TensorAccessors::contract<1, rank_1, rank_2, dim>(result, reord_01, reord_02);
  return result;
}


/**
 * 对任意等级的两个张量的两对索引进行通用收缩。返回一个等级为
 * $(\text{rank}_1 + \text{rank}_2
 *
 * - 4)$ 的张量，它是索引 @p index_1 与索引 @p index_2, 和索引 @p index_3与索引 @p index_4 的等级为 @p rank_1 的张量 @p src2 的收缩@f[
 * \text{result}_{i_1,\ldots,i_{r1},j_1,\ldots,j_{r2}}
 * = \sum_{k, l}
 *   \text{left}_{i_1,\ldots,k,\ldots,l,\ldots,i_{r1}}
 *   \text{right}_{j_1,\ldots,k,\ldots,l\ldots,j_{r2}}
 * @f]
 * 例如，如果第一个索引（ <code>index_1==0</code>
 * ）应与第三个索引（ <code>index_2==2</code>
 * ）收缩，第二个索引（ <code>index_3==1</code>
 * ）与第一个索引（ <code>index_4==0</code>) of a tensor
 * <code>t2</code> ）收缩，这个函数应被调用为
 *
 * @code
 * double_contract<0, 2, 1, 0>(t1, t2);
 * @endcode
 *
 *
 * @note  索引的位置是从0开始计算的，即
 * $0\le\text{index}_i<\text{range}_i$  。
 *
 *
 * @note
 * 如果收缩产生一个等级为0的张量，那么标量数将作为一个未包装的数字类型返回。
 * @relatesalso  Tensor
 *
 *
 */
template <int index_1,
          int index_2,
          int index_3,
          int index_4,
          int rank_1,
          int rank_2,
          int dim,
          typename Number,
          typename OtherNumber>
constexpr inline
  typename Tensor<rank_1 + rank_2 - 4,
                  dim,
                  typename ProductType<Number, OtherNumber>::type>::tensor_type
  double_contract(const Tensor<rank_1, dim, Number> &     src1,
                  const Tensor<rank_2, dim, OtherNumber> &src2)
{
  Assert(0 <= index_1 && index_1 < rank_1,
         ExcMessage(
           "The specified index_1 must lie within the range [0,rank_1)"));
  Assert(0 <= index_3 && index_3 < rank_1,
         ExcMessage(
           "The specified index_3 must lie within the range [0,rank_1)"));
  Assert(index_1 != index_3,
         ExcMessage("index_1 and index_3 must not be the same"));
  Assert(0 <= index_2 && index_2 < rank_2,
         ExcMessage(
           "The specified index_2 must lie within the range [0,rank_2)"));
  Assert(0 <= index_4 && index_4 < rank_2,
         ExcMessage(
           "The specified index_4 must lie within the range [0,rank_2)"));
  Assert(index_2 != index_4,
         ExcMessage("index_2 and index_4 must not be the same"));

  using namespace TensorAccessors;
  using namespace TensorAccessors::internal;

  // Reorder index_1 to the end of src1:
  ReorderedIndexView<index_1, rank_1, const Tensor<rank_1, dim, Number>>
    reord_1 = TensorAccessors::reordered_index_view<index_1, rank_1>(src1);

  // Reorder index_2 to the end of src2:
  ReorderedIndexView<index_2, rank_2, const Tensor<rank_2, dim, OtherNumber>>
    reord_2 = TensorAccessors::reordered_index_view<index_2, rank_2>(src2);

  // Now, reorder index_3 to the end of src1. We have to make sure to
  // preserve the original ordering: index_1 has been removed. If
  // index_3 > index_1, we have to use (index_3 - 1) instead:
  ReorderedIndexView<
    (index_3 < index_1 ? index_3 : index_3 - 1),
    rank_1,
    ReorderedIndexView<index_1, rank_1, const Tensor<rank_1, dim, Number>>>
    reord_3 =
      TensorAccessors::reordered_index_view < index_3 < index_1 ? index_3 :
                                                                  index_3 - 1,
    rank_1 > (reord_1);

  // Now, reorder index_4 to the end of src2. We have to make sure to
  // preserve the original ordering: index_2 has been removed. If
  // index_4 > index_2, we have to use (index_4 - 1) instead:
  ReorderedIndexView<
    (index_4 < index_2 ? index_4 : index_4 - 1),
    rank_2,
    ReorderedIndexView<index_2, rank_2, const Tensor<rank_2, dim, OtherNumber>>>
    reord_4 =
      TensorAccessors::reordered_index_view < index_4 < index_2 ? index_4 :
                                                                  index_4 - 1,
    rank_2 > (reord_2);

  typename Tensor<rank_1 + rank_2 - 4,
                  dim,
                  typename ProductType<Number, OtherNumber>::type>::tensor_type
    result{};
  TensorAccessors::contract<2, rank_1, rank_2, dim>(result, reord_3, reord_4);
  return result;
}


/**
 * 标量积，或者两个等阶的张量的（广义）Frobenius内积。返回一个标量，它是张量 @p left 和 @p right:  @f[
 * \sum_{i_1,\ldots,i_r}
 * \text{left}_{i_1,\ldots,i_r}
 * \text{right}_{i_1,\ldots,i_r}
 * @f]完全收缩的结果。
 * @relatesalso  Tensor
 *
 *
 */
template <int rank, int dim, typename Number, typename OtherNumber>
constexpr inline DEAL_II_ALWAYS_INLINE
  typename ProductType<Number, OtherNumber>::type
  scalar_product(const Tensor<rank, dim, Number> &     left,
                 const Tensor<rank, dim, OtherNumber> &right)
{
  typename ProductType<Number, OtherNumber>::type result{};
  TensorAccessors::contract<rank, rank, rank, dim>(result, left, right);
  return result;
}


/**
 * 对三个张量进行完全收缩。返回一个标量，它是等级为 @p rank_1, 的张量 @p left 、等级为 $(\text{rank}_1+\text{rank}_2)$ 的张量 @p middle 和等级为 @p rank_2: 的张量 @p 完全缩减的结果 @f[
 * \sum_{i_1,\ldots,i_{r1},j_1,\ldots,j_{r2}}
 * \text{left}_{i_1,\ldots,i_{r1}}
 * \text{middle}_{i_1,\ldots,i_{r1},j_1,\ldots,j_{r2}}
 * \text{right}_{j_1,\ldots,j_{r2}}
 * @f]
 *
 * @note
 * 三个输入张量中的每一个都可以是Tensor或SymmetricTensor。
 * @relatesalso  Tensor
 *
 *
 */
template <template <int, int, typename> class TensorT1,
          template <int, int, typename> class TensorT2,
          template <int, int, typename> class TensorT3,
          int rank_1,
          int rank_2,
          int dim,
          typename T1,
          typename T2,
          typename T3>
constexpr inline DEAL_II_ALWAYS_INLINE
  typename ProductType<T1, typename ProductType<T2, T3>::type>::type
  contract3(const TensorT1<rank_1, dim, T1> &         left,
            const TensorT2<rank_1 + rank_2, dim, T2> &middle,
            const TensorT3<rank_2, dim, T3> &         right)
{
  using return_type =
    typename ProductType<T1, typename ProductType<T2, T3>::type>::type;
  return TensorAccessors::contract3<rank_1, rank_2, dim, return_type>(left,
                                                                      middle,
                                                                      right);
}


/**
 * @p rank_1 和 @p rank_2: 两个张量的外积 返回一个等级为 $(\text{rank}_1 + \text{rank}_2)$ 的张量：@f[
 * \text{result}_{i_1,\ldots,i_{r1},j_1,\ldots,j_{r2}}
 * = \text{left}_{i_1,\ldots,i_{r1}}\,\text{right}_{j_1,\ldots,j_{r2}.}
 * @f]
 * @relatesalso  Tensor
 *
 */
template <int rank_1,
          int rank_2,
          int dim,
          typename Number,
          typename OtherNumber>
constexpr inline DEAL_II_ALWAYS_INLINE
  Tensor<rank_1 + rank_2, dim, typename ProductType<Number, OtherNumber>::type>
  outer_product(const Tensor<rank_1, dim, Number> &     src1,
                const Tensor<rank_2, dim, OtherNumber> &src2)
{
  typename Tensor<rank_1 + rank_2,
                  dim,
                  typename ProductType<Number, OtherNumber>::type>::tensor_type
    result{};
  TensorAccessors::contract<0, rank_1, rank_2, dim>(result, src1, src2);
  return result;
}


//@}
/**
 * @name  对等级为1的张量的特殊操作
 *
 */
//@{


/**
 * 返回2d中的交叉积。这只是顺时针旋转90度来计算切向矢量的外法线。这个函数是为所有空间维度定义的，以允许独立于维度的编程（例如在空间维度的开关内），但只有当参数的实际维度为2时才可以调用（例如来自开关中的<tt>dim==2</tt>情况）。
 * @relatesalso  Tensor
 *
 *
 */
template <int dim, typename Number>
constexpr inline DEAL_II_ALWAYS_INLINE Tensor<1, dim, Number>
                                       cross_product_2d(const Tensor<1, dim, Number> &src)
{
  Assert(dim == 2, ExcInternalError());

  Tensor<1, dim, Number> result;

  result[0] = src[1];
  result[1] = -src[0];

  return result;
}


/**
 * 返回2个向量在3D中的交叉积。这个函数是为所有空间维度定义的，以允许独立于维度的编程（例如在空间维度的开关内），但只有当参数的实际维度是3时才可以调用（例如来自开关中的<tt>dim==3</tt>情况）。
 * @relatesalso  Tensor
 *
 *
 */
template <int dim, typename Number1, typename Number2>
constexpr inline DEAL_II_ALWAYS_INLINE
  Tensor<1, dim, typename ProductType<Number1, Number2>::type>
  cross_product_3d(const Tensor<1, dim, Number1> &src1,
                   const Tensor<1, dim, Number2> &src2)
{
  Assert(dim == 3, ExcInternalError());

  Tensor<1, dim, typename ProductType<Number1, Number2>::type> result;

  // avoid compiler warnings
  constexpr int s0 = 0 % dim;
  constexpr int s1 = 1 % dim;
  constexpr int s2 = 2 % dim;

  result[s0] = src1[s1] * src2[s2] - src1[s2] * src2[s1];
  result[s1] = src1[s2] * src2[s0] - src1[s0] * src2[s2];
  result[s2] = src1[s0] * src2[s1] - src1[s1] * src2[s0];

  return result;
}


//@}
/**
 * @name  对等级2的张量的特殊操作
 *
 *
 */
//@{


/**
 * 计算一个张量或等级2的行列式。
 * @relatesalso  Tensor
 *
 *
 */
template <int dim, typename Number>
constexpr inline DEAL_II_ALWAYS_INLINE Number
                                       determinant(const Tensor<2, dim, Number> &t)
{
  // Compute the determinant using the Laplace expansion of the
  // determinant. We expand along the last row.
  Number det = internal::NumberType<Number>::value(0.0);

  for (unsigned int k = 0; k < dim; ++k)
    {
      Tensor<2, dim - 1, Number> minor;
      for (unsigned int i = 0; i < dim - 1; ++i)
        for (unsigned int j = 0; j < dim - 1; ++j)
          minor[i][j] = t[i][j < k ? j : j + 1];

      const Number cofactor = ((k % 2 == 0) ? -1. : 1.) * determinant(minor);

      det += t[dim - 1][k] * cofactor;
    }

  return ((dim % 2 == 0) ? 1. : -1.) * det;
}

/**
 * dim==1的特殊化。
 * @relatesalso  Tensor
 *
 */
template <typename Number>
constexpr DEAL_II_ALWAYS_INLINE Number
                                determinant(const Tensor<2, 1, Number> &t)
{
  return t[0][0];
}

/**
 * 对dim==2的特化。
 * @relatesalso  Tensor
 *
 */
template <typename Number>
constexpr DEAL_II_ALWAYS_INLINE Number
                                determinant(const Tensor<2, 2, Number> &t)
{
  // hard-coded for efficiency reasons
  return t[0][0] * t[1][1] - t[1][0] * t[0][1];
}

/**
 * 对dim的特殊化==3。
 * @relatesalso  Tensor
 *
 */
template <typename Number>
constexpr DEAL_II_ALWAYS_INLINE Number
                                determinant(const Tensor<2, 3, Number> &t)
{
  // hard-coded for efficiency reasons
  const Number C0 = internal::NumberType<Number>::value(t[1][1] * t[2][2]) -
                    internal::NumberType<Number>::value(t[1][2] * t[2][1]);
  const Number C1 = internal::NumberType<Number>::value(t[1][2] * t[2][0]) -
                    internal::NumberType<Number>::value(t[1][0] * t[2][2]);
  const Number C2 = internal::NumberType<Number>::value(t[1][0] * t[2][1]) -
                    internal::NumberType<Number>::value(t[1][1] * t[2][0]);
  return t[0][0] * C0 + t[0][1] * C1 + t[0][2] * C2;
}


/**
 * 计算并返回等级为2的张量的轨迹，即其对角线项的总和。
 * @relatesalso  Tensor
 *
 *
 */
template <int dim, typename Number>
constexpr inline DEAL_II_ALWAYS_INLINE Number
                                       trace(const Tensor<2, dim, Number> &d)
{
  Number t = d[0][0];
  for (unsigned int i = 1; i < dim; ++i)
    t += d[i][i];
  return t;
}


/**
 * 计算并返回给定张量的逆值。由于编译器可以进行返回值的优化，并且由于返回对象的大小是已知的，所以可以接受通过值来返回结果，而不是通过引用作为一个参数。
 * @relatesalso  Tensor
 *
 *
 */
template <int dim, typename Number>
constexpr inline Tensor<2, dim, Number>
invert(const Tensor<2, dim, Number> &)
{
  Number return_tensor[dim][dim];

  // if desired, take over the
  // inversion of a 4x4 tensor
  // from the FullMatrix
  AssertThrow(false, ExcNotImplemented());

  return Tensor<2, dim, Number>(return_tensor);
}


#ifndef DOXYGEN

template <typename Number>
constexpr inline DEAL_II_ALWAYS_INLINE Tensor<2, 1, Number>
                                       invert(const Tensor<2, 1, Number> &t)
{
  Tensor<2, 1, Number> return_tensor;

  return_tensor[0][0] = internal::NumberType<Number>::value(1.0 / t[0][0]);

  return return_tensor;
}


template <typename Number>
constexpr inline DEAL_II_ALWAYS_INLINE Tensor<2, 2, Number>
                                       invert(const Tensor<2, 2, Number> &t)
{
  Tensor<2, 2, Number> return_tensor;

  const Number inv_det_t = internal::NumberType<Number>::value(
    1.0 / (t[0][0] * t[1][1] - t[1][0] * t[0][1]));
  return_tensor[0][0] = t[1][1];
  return_tensor[0][1] = -t[0][1];
  return_tensor[1][0] = -t[1][0];
  return_tensor[1][1] = t[0][0];
  return_tensor *= inv_det_t;

  return return_tensor;
}


template <typename Number>
constexpr inline DEAL_II_ALWAYS_INLINE Tensor<2, 3, Number>
                                       invert(const Tensor<2, 3, Number> &t)
{
  Tensor<2, 3, Number> return_tensor;

  return_tensor[0][0] = internal::NumberType<Number>::value(t[1][1] * t[2][2]) -
                        internal::NumberType<Number>::value(t[1][2] * t[2][1]);
  return_tensor[0][1] = internal::NumberType<Number>::value(t[0][2] * t[2][1]) -
                        internal::NumberType<Number>::value(t[0][1] * t[2][2]);
  return_tensor[0][2] = internal::NumberType<Number>::value(t[0][1] * t[1][2]) -
                        internal::NumberType<Number>::value(t[0][2] * t[1][1]);
  return_tensor[1][0] = internal::NumberType<Number>::value(t[1][2] * t[2][0]) -
                        internal::NumberType<Number>::value(t[1][0] * t[2][2]);
  return_tensor[1][1] = internal::NumberType<Number>::value(t[0][0] * t[2][2]) -
                        internal::NumberType<Number>::value(t[0][2] * t[2][0]);
  return_tensor[1][2] = internal::NumberType<Number>::value(t[0][2] * t[1][0]) -
                        internal::NumberType<Number>::value(t[0][0] * t[1][2]);
  return_tensor[2][0] = internal::NumberType<Number>::value(t[1][0] * t[2][1]) -
                        internal::NumberType<Number>::value(t[1][1] * t[2][0]);
  return_tensor[2][1] = internal::NumberType<Number>::value(t[0][1] * t[2][0]) -
                        internal::NumberType<Number>::value(t[0][0] * t[2][1]);
  return_tensor[2][2] = internal::NumberType<Number>::value(t[0][0] * t[1][1]) -
                        internal::NumberType<Number>::value(t[0][1] * t[1][0]);
  const Number inv_det_t = internal::NumberType<Number>::value(
    1.0 / (t[0][0] * return_tensor[0][0] + t[0][1] * return_tensor[1][0] +
           t[0][2] * return_tensor[2][0]));
  return_tensor *= inv_det_t;

  return return_tensor;
}

#endif  /* DOXYGEN */ 


/**
 * 返回给定张量的转置。
 * @relatesalso  Tensor
 *
 *
 */
template <int dim, typename Number>
constexpr inline DEAL_II_ALWAYS_INLINE Tensor<2, dim, Number>
                                       transpose(const Tensor<2, dim, Number> &t)
{
  Tensor<2, dim, Number> tt;
  for (unsigned int i = 0; i < dim; ++i)
    {
      tt[i][i] = t[i][i];
      for (unsigned int j = i + 1; j < dim; ++j)
        {
          tt[i][j] = t[j][i];
          tt[j][i] = t[i][j];
        };
    }
  return tt;
}


/**
 * 返回等级为2的给定张量的邻接值。张量的邻接 $\mathbf A$ 被定义为@f[
 * \textrm{adj}\mathbf A
 * \dealcoloneq \textrm{det}\mathbf A \; \mathbf{A}^{-1}
 * \; .
 * @f]。
 *
 * @note  这要求张量是可逆的。
 * @relatesalso  Tensor
 *
 *
 */
template <int dim, typename Number>
constexpr Tensor<2, dim, Number>
adjugate(const Tensor<2, dim, Number> &t)
{
  return determinant(t) * invert(t);
}


/**
 * 返回给定等级为2的张量的协因子。张量的协因子 $\mathbf A$ 定义为@f[
 * \textrm{cof}\mathbf A
 * \dealcoloneq \textrm{det}\mathbf A \; \mathbf{A}^{-T}
 *  = \left[ \textrm{adj}\mathbf A \right]^{T} \; .
 * @f] 。
 *
 * @note  这要求该张量是可逆的。
 * @relatesalso  Tensor
 *
 *
 */
template <int dim, typename Number>
constexpr Tensor<2, dim, Number>
cofactor(const Tensor<2, dim, Number> &t)
{
  return transpose(adjugate(t));
}


/**
 * 通过结合给定输入 ${\mathbf A}$ 的奇异值分解（SVD）
 * ${\mathbf A}=\mathbf U  \mathbf S \mathbf V^T$
 * 的产物，返回最近的正交矩阵 $\hat {\mathbf A}=\mathbf U
 * \mathbf{V}^T$ ，有效地用身份矩阵代替 $\mathbf S$ 。
 * 这是一个（非线性）[投影操作](https://en.wikipedia.org/wiki/Projection_(mathematics))，因为当应用两次时，我们有
 * $\hat{\hat{\mathbf A}}=\hat{\mathbf A}$
 * ，这很容易看到。(这是因为 $\hat {\mathbf A}$ 的SVD仅仅是
 * $\mathbf U \mathbf I \mathbf{V}^T$ 。) 此外， $\hat {\mathbf A}$
 * 实际上是一个正交矩阵，因为正交矩阵必须满足 ${\hat
 * {\mathbf A}}^T \hat {\mathbf A}={\mathbf I}$ ，这里意味着
 *
 * @f{align*}{
 * {\hat {\mathbf A}}^T \hat {\mathbf A}
 * &=
 * \left(\mathbf U \mathbf{V}^T\right)^T\left(\mathbf U \mathbf{V}^T\right)
 * \\
 * &=
 * \mathbf V \mathbf{U}^T
 * \mathbf U \mathbf{V}^T
 * \\
 * &=
 * \mathbf V \left(\mathbf{U}^T
 * \mathbf U\right) \mathbf{V}^T
 * \\
 * &=
 * \mathbf V \mathbf I \mathbf{V}^T
 * \\
 * &=
 * \mathbf V \mathbf{V}^T
 * \\
 * &=
 * \mathbf I
 * @f}
 * 由于从SVD出来的 $\mathbf U$ 和 $\mathbf V$
 * 因子本身就是正交矩阵。
 * @param  A 要为其寻找最接近的正交张量。  @tparam  Number
 * 用来存储张量的条目的类型。  必须是`float`或`double`。
 * @pre  为了使用这个函数，这个程序必须与LAPACK库相连。
 * @pre   @p A
 * 不得为单数。这是因为，从概念上讲，这里要解决的问题是试图找到一个矩阵
 * $\hat{\mathbf A}$ ，使其与 $\mathbf A$
 * 的某种距离最小，同时满足二次约束 ${\hat {\mathbf A}}^T \hat
 * {\mathbf A}={\mathbf I}$  。这与想要找到一个矢量 $\hat{\mathbf
 * x}\in{\mathbb R}^n$ 的问题没有什么不同，该矢量对于给定的
 * $\mathbf x$ 来说，在约束条件 $\|\mathbf x\|^2=1$
 * 下最小化二次目标函数 $\|\hat {\mathbf x}
 *
 * - \mathbf x\|^2$  。
 *
 * --换句话说，我们正在寻找单位球面上最接近 $\mathbf x$
 * 的点 $\hat{\mathbf x}$ 。这个问题对所有的 $\mathbf x$
 * 都有解，除了 $\mathbf x=0$
 * 。我们在这里考虑的问题的相应条件是， $\mathbf A$
 * 不能有零特征值。
 * @relatesalso  Tensor
 *
 *
 */
template <int dim, typename Number>
Tensor<2, dim, Number>
project_onto_orthogonal_tensors(const Tensor<2, dim, Number> &A);


/**
 * 返回给定等级2张量的 $l_1$ 规范，其中 $\|\mathbf T\|_1 =
 * \max_j \sum_i |T_{ij}|$ （列上之和的最大值）。
 * @relatesalso Tensor
 *
 *
 */
template <int dim, typename Number>
inline Number
l1_norm(const Tensor<2, dim, Number> &t)
{
  Number max = internal::NumberType<Number>::value(0.0);
  for (unsigned int j = 0; j < dim; ++j)
    {
      Number sum = internal::NumberType<Number>::value(0.0);
      for (unsigned int i = 0; i < dim; ++i)
        sum += std::fabs(t[i][j]);

      if (sum > max)
        max = sum;
    }

  return max;
}


/**
 * 返回给定的2级张量的 $l_\infty$ 规范，其中 $\|\mathbf
 * T\|_\infty = \max_i \sum_j |T_{ij}|$ （行上和的最大值）。
 * @relatesalso  Tensor
 *
 *
 */
template <int dim, typename Number>
inline Number
linfty_norm(const Tensor<2, dim, Number> &t)
{
  Number max = internal::NumberType<Number>::value(0.0);
  for (unsigned int i = 0; i < dim; ++i)
    {
      Number sum = internal::NumberType<Number>::value(0.0);
      for (unsigned int j = 0; j < dim; ++j)
        sum += std::fabs(t[i][j]);

      if (sum > max)
        max = sum;
    }

  return max;
}

//@}


#ifndef DOXYGEN


#  ifdef DEAL_II_ADOLC_WITH_ADVANCED_BRANCHING

// Specialization of functions for ADOL-C number types when
// the advanced branching feature is used
template <int dim>
inline adouble
l1_norm(const Tensor<2, dim, adouble> &t)
{
  adouble max = internal::NumberType<adouble>::value(0.0);
  for (unsigned int j = 0; j < dim; ++j)
    {
      adouble sum = internal::NumberType<adouble>::value(0.0);
      for (unsigned int i = 0; i < dim; ++i)
        sum += std::fabs(t[i][j]);

      condassign(max, (sum > max), sum, max);
    }

  return max;
}


template <int dim>
inline adouble
linfty_norm(const Tensor<2, dim, adouble> &t)
{
  adouble max = internal::NumberType<adouble>::value(0.0);
  for (unsigned int i = 0; i < dim; ++i)
    {
      adouble sum = internal::NumberType<adouble>::value(0.0);
      for (unsigned int j = 0; j < dim; ++j)
        sum += std::fabs(t[i][j]);

      condassign(max, (sum > max), sum, max);
    }

  return max;
}

#  endif // DEAL_II_ADOLC_WITH_ADVANCED_BRANCHING


#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif


