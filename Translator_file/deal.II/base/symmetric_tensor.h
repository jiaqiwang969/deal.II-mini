//include/deal.II-translator/base/symmetric_tensor_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2021 by the deal.II authors
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

#ifndef dealii_symmetric_tensor_h
#define dealii_symmetric_tensor_h


#include <deal.II/base/config.h>

#include <deal.II/base/numbers.h>
#include <deal.II/base/table_indices.h>
#include <deal.II/base/template_constraints.h>
#include <deal.II/base/tensor.h>

#include <algorithm>
#include <array>
#include <functional>

DEAL_II_NAMESPACE_OPEN

// Forward declaration
#ifndef DOXYGEN
template <int rank, int dim, typename Number = double>
class SymmetricTensor;
#endif

template <int dim, typename Number>
DEAL_II_CONSTEXPR inline DEAL_II_ALWAYS_INLINE SymmetricTensor<2, dim, Number>
                                               unit_symmetric_tensor();

template <int dim, typename Number>
DEAL_II_CONSTEXPR inline DEAL_II_ALWAYS_INLINE SymmetricTensor<4, dim, Number>
                                               deviator_tensor();

template <int dim, typename Number>
DEAL_II_CONSTEXPR inline DEAL_II_ALWAYS_INLINE SymmetricTensor<4, dim, Number>
                                               identity_tensor();

template <int dim, typename Number>
constexpr DEAL_II_ALWAYS_INLINE SymmetricTensor<2, dim, Number>
                                invert(const SymmetricTensor<2, dim, Number> &);

template <int dim, typename Number>
constexpr DEAL_II_ALWAYS_INLINE SymmetricTensor<4, dim, Number>
                                invert(const SymmetricTensor<4, dim, Number> &);

template <int dim2, typename Number>
constexpr inline DEAL_II_ALWAYS_INLINE Number
                                       trace(const SymmetricTensor<2, dim2, Number> &);

template <int dim, typename Number>
constexpr inline DEAL_II_ALWAYS_INLINE SymmetricTensor<2, dim, Number>
                                       deviator(const SymmetricTensor<2, dim, Number> &);

template <int dim, typename Number>
constexpr inline DEAL_II_ALWAYS_INLINE Number
                                       determinant(const SymmetricTensor<2, dim, Number> &);



namespace internal
{
  // Workaround: The following 4 overloads are necessary to be able to
  // compile the library with Apple Clang 8 and older. We should remove
  // these overloads again when we bump the minimal required version to
  // something later than clang-3.6 / Apple Clang 6.3.
  template <int rank, int dim, typename T, typename U>
  struct ProductTypeImpl<SymmetricTensor<rank, dim, T>, std::complex<U>>
  {
    using type =
      SymmetricTensor<rank,
                      dim,
                      std::complex<typename ProductType<T, U>::type>>;
  };

  template <int rank, int dim, typename T, typename U>
  struct ProductTypeImpl<SymmetricTensor<rank, dim, std::complex<T>>,
                         std::complex<U>>
  {
    using type =
      SymmetricTensor<rank,
                      dim,
                      std::complex<typename ProductType<T, U>::type>>;
  };

  template <typename T, int rank, int dim, typename U>
  struct ProductTypeImpl<std::complex<T>, SymmetricTensor<rank, dim, U>>
  {
    using type =
      SymmetricTensor<rank,
                      dim,
                      std::complex<typename ProductType<T, U>::type>>;
  };

  template <int rank, int dim, typename T, typename U>
  struct ProductTypeImpl<std::complex<T>,
                         SymmetricTensor<rank, dim, std::complex<U>>>
  {
    using type =
      SymmetricTensor<rank,
                      dim,
                      std::complex<typename ProductType<T, U>::type>>;
  };
  // end workaround

  /**
   * 一个命名空间，用于SymmetricTensor类（及其相关函数）工作方式的内部函数和类。
   *
   */
  namespace SymmetricTensorImplementation
  {
    /**
     * 计算通用 @p rank,  @p dim 和 @p Number
     * 类型的对称张量的逆运算。
     *
     */
    template <int rank, int dim, typename Number>
    struct Inverse;
  } // namespace SymmetricTensorImplementation

  /**
   * 一个命名空间，用于SymmetricTensor类工作方式的内部类。
   *
   */
  namespace SymmetricTensorAccessors
  {
    /**
     * 创建一个TableIndices<2>对象，其中到<tt>position-1</tt>的第一个条目来自previous_indices，而new_index被放在<tt>position</tt>的位置。其余的指数仍处于无效状态。
     *
     */
    DEAL_II_CONSTEXPR inline DEAL_II_ALWAYS_INLINE TableIndices<2>
                                                   merge(const TableIndices<2> &previous_indices,
                                                         const unsigned int     new_index,
                                                         const unsigned int     position)
    {
      AssertIndexRange(position, 2);

      if (position == 0)
        return {new_index, numbers::invalid_unsigned_int};
      else
        return {previous_indices[0], new_index};
    }



    /**
     * 创建一个TableIndices<4>对象，其中到<tt>position-1</tt>的第一个条目取自previous_indices，new_index被放在<tt>position</tt>的位置。其余的指数仍处于无效状态。
     *
     */
    DEAL_II_CONSTEXPR inline DEAL_II_ALWAYS_INLINE TableIndices<4>
                                                   merge(const TableIndices<4> &previous_indices,
                                                         const unsigned int     new_index,
                                                         const unsigned int     position)
    {
      AssertIndexRange(position, 4);

      switch (position)
        {
          case 0:
            return {new_index,
                    numbers::invalid_unsigned_int,
                    numbers::invalid_unsigned_int,
                    numbers::invalid_unsigned_int};
          case 1:
            return {previous_indices[0],
                    new_index,
                    numbers::invalid_unsigned_int,
                    numbers::invalid_unsigned_int};
          case 2:
            return {previous_indices[0],
                    previous_indices[1],
                    new_index,
                    numbers::invalid_unsigned_int};
          case 3:
            return {previous_indices[0],
                    previous_indices[1],
                    previous_indices[2],
                    new_index};
          default:
            Assert(false, ExcInternalError());
            return {};
        }
    }


    /**
     * Typedef模板魔法，表示两个张量或等级rank1和rank2之间的双倍收缩结果。一般来说，这是一个等级<tt>rank1+rank2-4</tt>的张量，但如果这是零的话，就是一个单一的标量Number。对于这种情况，我们有一个特殊化。
     *
     */
    template <int rank1,
              int rank2,
              int dim,
              typename Number,
              typename OtherNumber = Number>
    struct double_contraction_result
    {
      using value_type = typename ProductType<Number, OtherNumber>::type;
      using type =
        ::dealii::SymmetricTensor<rank1 + rank2 - 4, dim, value_type>;
    };


    /**
     * Typedef模板魔法表示两个张量或等级rank1和rank2之间的双重收缩的结果。一般来说，这是一个等级<tt>rank1+rank2-4</tt>的张量，但如果这是零的话，就是一个单一的标量数。对于这种情况，我们有一个特殊化。
     *
     */
    template <int dim, typename Number, typename OtherNumber>
    struct double_contraction_result<2, 2, dim, Number, OtherNumber>
    {
      using type = typename ProductType<Number, OtherNumber>::type;
    };



    /**
     * 声明用于存储对称张量的数据结构类型的别名。例如，对于等级2的对称张量，我们使用一个平面向量来存储所有的元素。另一方面，对称秩-4张量是从对称秩-2张量映射到对称秩-2张量，所以它们可以被表示为矩阵等。
     * 除了需要这些信息的访问器类之外，其他人可能对这些信息不感兴趣。特别是，你不应该在你的应用程序中对存储格式做任何假设。
     *
     */
    template <int rank, int dim, typename Number>
    struct StorageType;

    /**
     * 用于秩2张量的StorageType的特殊化。
     *
     */
    template <int dim, typename Number>
    struct StorageType<2, dim, Number>
    {
      /**
       * 等级2的对称张量的独立成分的数量。我们只存储它的右上半部分。
       *
       */
      static const unsigned int n_independent_components =
        (dim * dim + dim) / 2;

      /**
       * 声明我们实际存储数据的类型。
       *
       */
      using base_tensor_type = Tensor<1, n_independent_components, Number>;
    };



    /**
     * 用于等级4张量的StorageType的特殊化。
     *
     */
    template <int dim, typename Number>
    struct StorageType<4, dim, Number>
    {
      /**
       * 等级2的对称张量的独立分量的数量。
       * 由于秩-4张量是这类对象之间的映射，我们需要这些信息。
       *
       */
      static const unsigned int n_rank2_components = (dim * dim + dim) / 2;

      /**
       * 等级4的对称张量的独立分量的数量。
       *
       */
      static const unsigned int n_independent_components =
        (n_rank2_components *
         StorageType<2, dim, Number>::n_independent_components);

      /**
       * 声明我们实际存储数据的类型。对称秩4张量是对称秩2张量之间的映射，所以如果我们将秩2张量表示为向量，我们可以将数据表示为矩阵。
       *
       */
      using base_tensor_type = Tensor<2, n_rank2_components, Number>;
    };



    /**
     * 切换类型，选择一个秩2和维度<tt>dim</tt>的张量，切换张量是否应该是常数。
     *
     */
    template <int rank, int dim, bool constness, typename Number>
    struct AccessorTypes;

    /**
     * 切换类型选择一个等级为2，维度为<tt>dim</tt>的张量，切换张量是否应该是常数。
     * 对常数张量的特殊化。
     *
     */
    template <int rank, int dim, typename Number>
    struct AccessorTypes<rank, dim, true, Number>
    {
      using tensor_type = const ::dealii::SymmetricTensor<rank, dim, Number>;

      using reference = Number;
    };

    /**
     * 切换类型选择等级2和维度<tt>dim</tt>的张量，切换张量是否应该是常数。
     * 对非常数张量的特殊化。
     *
     */
    template <int rank, int dim, typename Number>
    struct AccessorTypes<rank, dim, false, Number>
    {
      using tensor_type = ::dealii::SymmetricTensor<rank, dim, Number>;

      using reference = Number &;
    };


    /**
     * @internal
     * 作为SymmetricTensor类型元素的访问器的类。模板参数<tt>constness</tt>可以是true或false，表示所处理的对象是否是常数（即只有在值为false时才允许写访问）。
     * 因为在<tt>N</tt>索引中，应用<tt>operator[]</tt>的效果是获得对<tt>N-1</tt>索引的访问，我们必须递归地实现这些访问器类，当我们只剩下一个索引时就停止。对于后一种情况，下面声明了这个类的特殊化，调用<tt>operator[]</tt>可以访问张量实际存储的对象；张量类还确保只访问我们实际存储的那些元素，也就是说，必要时它会重新排列索引。模板参数<tt>P</tt>表示有多少个剩余的索引。对于一个等级为2的张量，<tt>P</tt>可能是两个，而当使用<tt>operator[]</tt>时，会出现一个<tt>P=1</tt>的对象。
     * 正如对整个命名空间所说的，你通常不会直接与这些类打交道，也不应该试图直接使用它们的接口，因为它可能会在没有通知的情况下改变。事实上，由于构造函数是私有的，你甚至不能生成这个类的对象，因为它们只被认为是访问表类元素的暂时性的，而不是作为函数的参数来传递它们，等等。
     * 这个类是一个用于表类的类似类的改编。
     *
     */
    template <int rank, int dim, bool constness, int P, typename Number>
    class Accessor
    {
    public:
      /**
       * 从上面的switch类中导入两个别名。
       *
       */
      using reference =
        typename AccessorTypes<rank, dim, constness, Number>::reference;
      using tensor_type =
        typename AccessorTypes<rank, dim, constness, Number>::tensor_type;

    private:
      /**
       * 构造函数。取一个我们要访问的张量对象的引用。
       * 第二个参数表示进入张量的前几个索引的值。例如，对于一个等级为4的张量，如果P=2，那么我们将已经有两个连续的元素选择（例如通过<tt>tensor[1][2]</tt>），这两个索引值必须被储存在某个地方。因此，这个类只利用这个数组的第一个等级-P元素，但用P-1将其传递给下一级，由其填充下一个条目，以此类推。
       * 构造函数是私有的，以防止你身边有这样的对象。创建这种对象的唯一方法是通过<tt>Table</tt>类，它只作为临时对象生成它们。
       * 这保证了访问器对象比母对象更早退出范围，避免了数据一致性的问题。
       *
       */
      constexpr Accessor(tensor_type &             tensor,
                         const TableIndices<rank> &previous_indices);

      /**
       * 复制构造函数。
       *
       */
      constexpr DEAL_II_ALWAYS_INLINE
      Accessor(const Accessor &) = default;

    public:
      /**
       * 索引操作符。
       *
       */
      constexpr Accessor<rank, dim, constness, P - 1, Number>
      operator[](const unsigned int i);

      /**
       * 索引操作符。
       *
       */
      constexpr Accessor<rank, dim, constness, P - 1, Number>
      operator[](const unsigned int i) const;

    private:
      /**
       * 存储给构造函数的数据。
       *
       */
      tensor_type &            tensor;
      const TableIndices<rank> previous_indices;

      // Declare some other classes as friends. Make sure to work around bugs
      // in some compilers:
      template <int, int, typename>
      friend class dealii::SymmetricTensor;
      template <int, int, bool, int, typename>
      friend class Accessor;
      friend class ::dealii::SymmetricTensor<rank, dim, Number>;
      friend class Accessor<rank, dim, constness, P + 1, Number>;
    };



    /**
     * @internal
     * SymmetricTensor的访问器类。这是最后一个索引的特殊化，它实际上允许访问表中的元素，而不是递归地返回进一步子集的访问对象。这种特殊化与一般模板的情况相同；更多信息请看那里。
     *
     */
    template <int rank, int dim, bool constness, typename Number>
    class Accessor<rank, dim, constness, 1, Number>
    {
    public:
      /**
       * 从上面的开关类中导入两个别名。
       *
       */
      using reference =
        typename AccessorTypes<rank, dim, constness, Number>::reference;
      using tensor_type =
        typename AccessorTypes<rank, dim, constness, Number>::tensor_type;

    private:
      /**
       * 构造函数。取一个我们要访问的张量对象的引用。
       * 第二个参数表示进入张量的前几个索引的值。例如，对于一个等级为4的张量，如果P=2，那么我们将已经有两个连续的元素选择（例如通过<tt>tensor[1][2]</tt>），这两个索引值必须被储存在某个地方。因此，这个类只利用这个数组的第一个等级-P元素，但用P-1将其传递给下一级，而P-1则填充了下一个条目，以此类推。
       * 对于这个特殊的特殊化，即对于P==1，除了最后一个索引，其他的都已经被填充了。
       * 构造函数是私有的，以防止你身边有这样的对象。创建这种对象的唯一方法是通过<tt>Table</tt>类，它只生成临时对象。
       * 这保证了访问器对象比母对象更早退出范围，避免了数据一致性的问题。
       *
       */
      constexpr Accessor(tensor_type &             tensor,
                         const TableIndices<rank> &previous_indices);

      /**
       * 复制构造函数。
       *
       */
      constexpr DEAL_II_ALWAYS_INLINE
      Accessor(const Accessor &) = default;

    public:
      /**
       * 索引操作符。
       *
       */
      constexpr reference operator[](const unsigned int);

      /**
       * 索引运算器。
       *
       */
      constexpr reference operator[](const unsigned int) const;

    private:
      /**
       * 存储给构造函数的数据。
       *
       */
      tensor_type &            tensor;
      const TableIndices<rank> previous_indices;

      // Declare some other classes as friends. Make sure to work around bugs
      // in some compilers:
      template <int, int, typename>
      friend class dealii::SymmetricTensor;
      template <int, int, bool, int, typename>
      friend class SymmetricTensorAccessors::Accessor;
      friend class ::dealii::SymmetricTensor<rank, dim, Number>;
      friend class SymmetricTensorAccessors::
        Accessor<rank, dim, constness, 2, Number>;
    };
  } // namespace SymmetricTensorAccessors
} // namespace internal



/**
 * 提供一个能有效地存储等级为2,4,...的对称张量的类，即只存储全张量中那些不多余的对角线元素。例如，对于对称
 * $2\times 2$
 * 张量，这将是元素11、22和12，而元素21等于12的元素。在本文档中，二阶对称张量用粗体大写拉丁字母表示，如
 * $\mathbf A, \mathbf B, \dots$  或粗体希腊字母，如
 * $\boldsymbol{\varepsilon}$  ,  $\boldsymbol{\sigma}$  。二阶张量如
 * $\mathbf A$ 的直角坐标表示为 $A_{ij}$ ，其中 $i,j$
 * 是0到<tt>dim-1</tt>的指数。
 * 在许多情况下，对秩为2的对称张量使用这个类比矩阵有优势，因为编译器知道维度以及数据的位置。因此有可能产生比运行时间相关的维度的矩阵更有效的代码。它也比使用更通用的<tt>Tensor</tt>类更有效，因为存储的元素更少，而且该类会自动确保张量代表一个对称对象。
 * 对于更高等级的张量，存储方面的节省甚至更高。例如，对于等级为4的
 * $3 \times 3 \times 3 \times 3$
 * 张量，只需要存储36个而不是全部81个条目。这些等级4的张量用黑板风格的大写拉丁字母表示，如
 * $\mathbb A$ ，其组成部分 $\mathcal{A}_{ijkl}$  。
 * 虽然对称秩-2张量的定义很明显，但如果秩4张量是将对称秩-2张量映射到对称秩-2张量上的算子，则被视为对称的。这种所谓的秩4张量的次要对称性要求对于每一组四个指数
 * $i, j, k, l$ ，身份 $\mathcal{C}_{ijkl} = \mathcal{C}_{jikl} =
 * \mathcal{C}_{ijlk}$ 成立。然而，它并不意味着关系
 * $\mathcal{C}_{ijkl} = \mathcal{C}_{klij}$
 * 。因此，这里所理解的等级4的对称张子只是将对称张子映射到对称张子上的张子，但它们不一定能引起对称标量积
 * $\mathbf A : \mathbb C : \mathbf B = \mathbf B : \mathbb C : \mathbf A$
 * ，甚至是正（半）定式 $\mathbf A : \mathbb C : \mathbf A$
 * ，其中 $\mathbf A, \mathbf B$
 * 是对称等级2张子，冒号表示常见的双指数收缩，作为对称张子的标量积作用。
 * 对称张量最常用于结构和流体力学，其中应变和应力通常是对称张量，应力-应变关系由对称秩4张量给出。
 *
 *
 * @note
 * 对称张量只存在于偶数的指数。换句话说，你可以使用的对象只有<tt>SymmetricTensor<2,dim></tt>,
 * <tt>SymmetricTensor<4,dim></tt>等，但是<tt>SymmetricTensor<1,dim></tt>和<tt>SymmetricTensor<3,dim></tt>并不存在，使用它们很可能导致编译器错误。
 *
 *  <h3>Accessing elements</h3> 张量 $\mathbb C$
 * 的元素可以使用括号运算符来访问，即对于等级为4的张量，<tt>C[0][1][0][1]</tt>可以访问元素
 * $\mathcal{C}_{0101}$
 * 。这种访问可以用于读和写（如果张量至少是非常数）。你也可以对它进行其他操作，尽管这可能会导致混乱的情况，因为张量的几个元素被存储在同一个位置。例如，对于一个假定开始时为零的秩-2张量，写<tt>A[0][1]+=1;
 * A[1][0]+=1;</tt>将导致同一个元素被增加1 <em> 两次 </em>
 * ，因为即使访问使用不同的索引，被访问的元素是对称的，因此存储在同一位置。因此，在应用程序中，将对单个元素的操作限制为简单的读或写可能是有用的。
 *
 *
 * @ingroup geomprimitives
 *
 *
 */
template <int rank_, int dim, typename Number>
class SymmetricTensor
{
public:
  static_assert(rank_ % 2 == 0, "A SymmetricTensor must have even rank!");

  /**
   * 提供一种方法来获取一个对象的尺寸，而不需要明确知道它的数据类型。实现这种方式而不是提供一个函数<tt>dimension()</tt>，因为现在有可能在编译时获得尺寸，而不需要内联函数的扩展和预评估；因此编译器可能产生更有效的代码，你可以使用这个值来声明其他数据类型。
   *
   */
  static const unsigned int dimension = dim;

  /**
   * 向外界公布这个张量的等级。
   *
   */
  static const unsigned int rank = rank_;

  /**
   * 一个整数，表示完全描述一个对称张量的独立成分的数量。在
   * $d$ 空间维度中，这个数字等于 $\frac 12 (d^2+d)$
   * 的对称张量的等级2。
   *
   */
  static constexpr unsigned int n_independent_components =
    internal::SymmetricTensorAccessors::StorageType<rank_, dim, Number>::
      n_independent_components;

  /**
   * 默认构造函数。创建一个所有条目都等于零的张量。
   *
   */
  constexpr DEAL_II_ALWAYS_INLINE
  SymmetricTensor() = default;

  /**
   * 构造函数。从一个一般的张量生成一个对称的张量。假设
   * @p t
   * 已经是对称的，在调试模式下，这实际上是检查的。
   * 注意，没有规定保证张量只在四舍五入的情况下是对称的：如果传入的张量不是完全对称的，那么会抛出一个异常。如果你知道传入的张量只到四舍五入为止是对称的，那么你可能想先调用<tt>symmetrize()</tt>函数。如果你不确定，在调用<tt>symmetrize()</tt>之前进行检查是很好的做法。
   * 因为我们通过非constexpr函数调用来检查对称性，所以你必须在constexpr上下文中使用symmetrize()函数来代替。
   *
   */
  template <typename OtherNumber>
  explicit SymmetricTensor(const Tensor<2, dim, OtherNumber> &t);

  /**
   * 一个构造函数，用于从持有独立元素的数组中创建一个对称的张量。使用这个构造函数假定调用者知道元素在对称张量中的存储顺序；因此不鼓励使用它，但是如果你认为你想使用它，你可以使用unrolled_index()函数查询元素的顺序。
   * 这个构造函数目前只针对等级为2的对称张量实现。
   * 传递的数组的大小等于
   * SymmetricTensor<rank_,dim>::n_independent_components;
   * 使用内部命名空间的对象的原因是为了解决一些旧编译器的错误。
   *
   */
  constexpr SymmetricTensor(const Number (&array)[n_independent_components]);

  /**
   * 从具有不同底层标量类型的张量复制构造函数。这显然要求
   * @p OtherNumber 类型可以转换为 @p 数。
   *
   */
  template <typename OtherNumber>
  constexpr explicit SymmetricTensor(
    const SymmetricTensor<rank_, dim, OtherNumber> &initializer);

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
   * 来自具有不同底层标量类型的对称张量的赋值运算符。
   * 这显然要求 @p OtherNumber 类型可以转换为 @p Number. 。
   *
   */
  template <typename OtherNumber>
  constexpr SymmetricTensor &
  operator=(const SymmetricTensor<rank_, dim, OtherNumber> &rhs);

  /**
   * 这个操作符将一个标量分配给一个张量。为了避免混淆将标量值分配给张量的确切含义，零是<tt>d</tt>唯一允许的值，允许用直观的符号
   * $\mathbf A = 0$ 将张量的所有元素重置为零。
   *
   */
  constexpr SymmetricTensor &
  operator=(const Number &d);

  /**
   * 将目前的对称张量转换成具有相同元素的全张量，但使用全张量的不同存储方案。
   *
   */
  constexpr operator Tensor<rank_, dim, Number>() const;

  /**
   * 测试两个张量的相等性。
   *
   */
  constexpr bool
  operator==(const SymmetricTensor &) const;

  /**
   * 测试两个张量的不等式。
   *
   */
  constexpr bool
  operator!=(const SymmetricTensor &) const;

  /**
   * 添加另一个张量。
   *
   */
  template <typename OtherNumber>
  constexpr SymmetricTensor &
  operator+=(const SymmetricTensor<rank_, dim, OtherNumber> &);

  /**
   * 减去另一个张量。
   *
   */
  template <typename OtherNumber>
  constexpr SymmetricTensor &
  operator-=(const SymmetricTensor<rank_, dim, OtherNumber> &);

  /**
   * 用<tt>因子</tt>缩放张量，即用<tt>因子</tt>乘以所有组件。
   *
   */
  template <typename OtherNumber>
  constexpr SymmetricTensor &
  operator*=(const OtherNumber &factor);

  /**
   * 用<tt>1/factor</tt>缩放张量。
   *
   */
  template <typename OtherNumber>
  constexpr SymmetricTensor &
  operator/=(const OtherNumber &factor);

  /**
   * 单元减法运算符。负掉张量的所有条目。
   *
   */
  constexpr SymmetricTensor
  operator-() const;

  /**
   * 本对称张量与一个等级为2的张量之间的双重收缩积。例如，如果当前对象是对称秩2张量
   * $\mathbf{A}$ ，它与另一个对称秩2张量 $\mathbf{B}$
   * 相乘，那么结果是标量积双收缩  $\mathbf A : \mathbf B =
   * \sum_{i,j} A_{ij} B_{ij}$  。
   * 在这种情况下，返回值被评估为一个单一的标量。虽然有可能定义其他标量积（以及相关的诱导规范），但这个似乎是最合适的一个。
   * 如果当前对象是一个秩-4张量，比如 $\mathbb A$
   * ，那么结果是一个秩-2张量 $\mathbf C = \mathbb A : \mathbf B$
   * ，也就是说，该操作在当前对象的最后两个索引和参数的索引上收缩，结果是一个秩2的张量（
   * $C_{ij} = \sum_{k,l} \mathcal{A}_{ijkl} B_{kl}$ ）。
   * 请注意，对称张量的乘法算子被定义为对两个索引的双重收缩，而对于普通<tt>张量</tt>对象，它被定义为对一个索引的单一收缩。因此，对于对称张量，它的作用方式在数学文献中通常用
   * "冒号乘法 "来表示。
   * 有一些全局函数<tt>double_contract</tt>做了与这个操作符相同的工作，但它们不是将结果作为返回值，而是将其写入函数的第一个参数中。
   *
   */
  template <typename OtherNumber>
  DEAL_II_CONSTEXPR typename internal::SymmetricTensorAccessors::
    double_contraction_result<rank_, 2, dim, Number, OtherNumber>::type
    operator*(const SymmetricTensor<2, dim, OtherNumber> &s) const;

  /**
   * 在本对象的两个索引上收缩，并将等级4的对称张量作为参数给出。
   *
   */
  template <typename OtherNumber>
  DEAL_II_CONSTEXPR typename internal::SymmetricTensorAccessors::
    double_contraction_result<rank_, 4, dim, Number, OtherNumber>::type
    operator*(const SymmetricTensor<4, dim, OtherNumber> &s) const;

  /**
   * 返回一个对指定元素的读写引用。
   *
   */
  constexpr Number &
  operator()(const TableIndices<rank_> &indices);

  /**
   * 返回一个 @p const 对参数所指的值的引用。
   *
   */
  constexpr const Number &
  operator()(const TableIndices<rank_> &indices) const;

  /**
   * 访问这个对称张量的某一行的元素。这个函数是为常数张量调用的。
   *
   */
  constexpr internal::SymmetricTensorAccessors::
    Accessor<rank_, dim, true, rank_ - 1, Number>
    operator[](const unsigned int row) const;

  /**
   * 访问这个对称张量的一行元素。对于非常数张量，这个函数被调用。
   *
   */
  constexpr internal::SymmetricTensorAccessors::
    Accessor<rank_, dim, false, rank_ - 1, Number>
    operator[](const unsigned int row);

  /**
   * 返回一个 @p const 对参数所指数值的引用。
   * 与operator()完全相同。
   *
   */
  constexpr const Number &operator[](const TableIndices<rank_> &indices) const;

  /**
   * 返回一个对指定元素的读写引用。
   * 与operator()完全相同。
   *
   */
  constexpr Number &operator[](const TableIndices<rank_> &indices);

  /**
   * 根据unrolled
   * index访问一个元素。函数<tt>s.access_raw_entry(unrolled_index)</tt>与<tt>s[s.unrolled_to_component_indices(unrolled_index)]</tt>的作用相同，但更有效。
   *
   */
  constexpr const Number &
  access_raw_entry(const unsigned int unrolled_index) const;

  /**
   * 根据unrolled
   * index来访问一个元素。函数<tt>s.access_raw_entry(unrolled_index)</tt>与<tt>s[s.unrolled_to_component_indices(unrolled_index)]</tt>的作用相同，但效率更高。
   *
   */
  constexpr Number &
  access_raw_entry(const unsigned int unrolled_index);

  /**
   * 返回张量的Frobenius-norm，即所有条目平方之和的平方根。这个准则是由上面为两个对称张量定义的标量乘引起的。请注意，它包括张量的<i>all</i>个条目，计算对称性，而不仅仅是唯一的条目（例如，对于等级2的张量，这个准则包括将右上和左下条目的平方相加，而不仅仅是其中之一，尽管它们对于对称张量是相等的）。
   *
   */
  constexpr typename numbers::NumberTraits<Number>::real_type
  norm() const;

  /**
   * 张量对象可以通过简单地将所有元素粘贴到一个长矢量中来展开，但为此必须定义一个元素的顺序。对于对称张量，该函数返回对称张量中给定条目在
   * <code>[0,n_independent_components)</code> 范围内的哪个索引。
   *
   */
  static constexpr unsigned int
  component_to_unrolled_index(const TableIndices<rank_> &indices);

  /**
   * 与前一个函数相反：给定张量的解卷形式中的一个索引
   * $i$ ，返回与之相对应的索引集 $(k,l)$
   * （对于秩-2张量）或 $(k,l,m,n)$ （对于秩-4张量）。
   *
   */
  static constexpr TableIndices<rank_>
  unrolled_to_component_indices(const unsigned int i);

  /**
   * 将所有数值重置为零。    请注意，这与标准库容器的 @p
   * clear()成员函数和deal.II中的其他几个类的语义部分不一致，它们不仅将存储元素的值重置为零，而且释放所有内存并将对象返回到处女状态。然而，由于本类型的对象的大小是由其模板参数决定的，所以调整大小是不可能的，事实上，所有元素的值都为零的状态就是这样一个对象构造后的状态。
   *
   */
  constexpr void
  clear();

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

private:
  /**
   * 一个描述基数张量的属性的结构。
   *
   */
  using base_tensor_descriptor =
    internal::SymmetricTensorAccessors::StorageType<rank_, dim, Number>;

  /**
   * 对称张量的数据存储类型。
   *
   */
  using base_tensor_type = typename base_tensor_descriptor::base_tensor_type;

  /**
   * 我们存储张量数据的地方。
   *
   */
  base_tensor_type data;

  // Make all other symmetric tensors friends.
  template <int, int, typename>
  friend class SymmetricTensor;

  // Make a few more functions friends.
  template <int dim2, typename Number2>
  friend constexpr Number2
  trace(const SymmetricTensor<2, dim2, Number2> &d);

  template <int dim2, typename Number2>
  friend constexpr Number2
  determinant(const SymmetricTensor<2, dim2, Number2> &t);

  template <int dim2, typename Number2>
  friend constexpr SymmetricTensor<2, dim2, Number2>
  deviator(const SymmetricTensor<2, dim2, Number2> &t);

  template <int dim2, typename Number2>
  friend DEAL_II_CONSTEXPR SymmetricTensor<2, dim2, Number2>
                           unit_symmetric_tensor();

  template <int dim2, typename Number2>
  friend DEAL_II_CONSTEXPR SymmetricTensor<4, dim2, Number2>
                           deviator_tensor();

  template <int dim2, typename Number2>
  friend DEAL_II_CONSTEXPR SymmetricTensor<4, dim2, Number2>
                           identity_tensor();


  // Make a few helper classes friends as well.
  friend struct internal::SymmetricTensorImplementation::
    Inverse<2, dim, Number>;

  friend struct internal::SymmetricTensorImplementation::
    Inverse<4, dim, Number>;
};



// ------------------------- inline functions ------------------------

#ifndef DOXYGEN

// provide declarations for static members
template <int rank, int dim, typename Number>
const unsigned int SymmetricTensor<rank, dim, Number>::dimension;

template <int rank_, int dim, typename Number>
constexpr unsigned int
  SymmetricTensor<rank_, dim, Number>::n_independent_components;

namespace internal
{
  namespace SymmetricTensorAccessors
  {
    template <int rank_, int dim, bool constness, int P, typename Number>
    constexpr DEAL_II_ALWAYS_INLINE
    Accessor<rank_, dim, constness, P, Number>::Accessor(
      tensor_type &              tensor,
      const TableIndices<rank_> &previous_indices)
      : tensor(tensor)
      , previous_indices(previous_indices)
    {}



    template <int rank_, int dim, bool constness, int P, typename Number>
    constexpr inline DEAL_II_ALWAYS_INLINE
        Accessor<rank_, dim, constness, P - 1, Number>
        Accessor<rank_, dim, constness, P, Number>::
        operator[](const unsigned int i)
    {
      return Accessor<rank_, dim, constness, P - 1, Number>(
        tensor, merge(previous_indices, i, rank_ - P));
    }



    template <int rank_, int dim, bool constness, int P, typename Number>
    constexpr DEAL_II_ALWAYS_INLINE
        Accessor<rank_, dim, constness, P - 1, Number>
        Accessor<rank_, dim, constness, P, Number>::
        operator[](const unsigned int i) const
    {
      return Accessor<rank_, dim, constness, P - 1, Number>(
        tensor, merge(previous_indices, i, rank_ - P));
    }



    template <int rank_, int dim, bool constness, typename Number>
    constexpr DEAL_II_ALWAYS_INLINE
    Accessor<rank_, dim, constness, 1, Number>::Accessor(
      tensor_type &              tensor,
      const TableIndices<rank_> &previous_indices)
      : tensor(tensor)
      , previous_indices(previous_indices)
    {}



    template <int rank_, int dim, bool constness, typename Number>
    constexpr inline DEAL_II_ALWAYS_INLINE
      typename Accessor<rank_, dim, constness, 1, Number>::reference
        Accessor<rank_, dim, constness, 1, Number>::
        operator[](const unsigned int i)
    {
      return tensor(merge(previous_indices, i, rank_ - 1));
    }


    template <int rank_, int dim, bool constness, typename Number>
    constexpr DEAL_II_ALWAYS_INLINE
      typename Accessor<rank_, dim, constness, 1, Number>::reference
        Accessor<rank_, dim, constness, 1, Number>::
        operator[](const unsigned int i) const
    {
      return tensor(merge(previous_indices, i, rank_ - 1));
    }
  } // namespace SymmetricTensorAccessors
} // namespace internal



template <int rank_, int dim, typename Number>
template <typename OtherNumber>
inline DEAL_II_ALWAYS_INLINE
SymmetricTensor<rank_, dim, Number>::SymmetricTensor(
  const Tensor<2, dim, OtherNumber> &t)
{
  static_assert(rank == 2, "This function is only implemented for rank==2");
  for (unsigned int d = 0; d < dim; ++d)
    for (unsigned int e = 0; e < d; ++e)
      Assert(t[d][e] == t[e][d],
             ExcMessage("The incoming Tensor must be exactly symmetric."));

  for (unsigned int d = 0; d < dim; ++d)
    data[d] = t[d][d];

  for (unsigned int d = 0, c = 0; d < dim; ++d)
    for (unsigned int e = d + 1; e < dim; ++e, ++c)
      data[dim + c] = t[d][e];
}



template <int rank_, int dim, typename Number>
template <typename OtherNumber>
constexpr DEAL_II_ALWAYS_INLINE
SymmetricTensor<rank_, dim, Number>::SymmetricTensor(
  const SymmetricTensor<rank_, dim, OtherNumber> &initializer)
  : data(initializer.data)
{}



template <int rank_, int dim, typename Number>
constexpr inline DEAL_II_ALWAYS_INLINE
SymmetricTensor<rank_, dim, Number>::SymmetricTensor(
  const Number (&array)[n_independent_components])
  : data(
      *reinterpret_cast<const typename base_tensor_type::array_type *>(array))
{
  // ensure that the reinterpret_cast above actually works
  Assert(sizeof(typename base_tensor_type::array_type) == sizeof(array),
         ExcInternalError());
}



template <int rank_, int dim, typename Number>
template <typename OtherNumber>
constexpr inline DEAL_II_ALWAYS_INLINE SymmetricTensor<rank_, dim, Number> &
                                       SymmetricTensor<rank_, dim, Number>::
                                       operator=(const SymmetricTensor<rank_, dim, OtherNumber> &t)
{
  data = t.data;
  return *this;
}



template <int rank_, int dim, typename Number>
constexpr inline DEAL_II_ALWAYS_INLINE SymmetricTensor<rank_, dim, Number> &
SymmetricTensor<rank_, dim, Number>::operator=(const Number &d)
{
  Assert(numbers::value_is_zero(d),
         ExcMessage("Only assignment with zero is allowed"));
  (void)d;

  data = internal::NumberType<Number>::value(0.0);

  return *this;
}


namespace internal
{
  namespace SymmetricTensorImplementation
  {
    template <int dim, typename Number>
    constexpr inline DEAL_II_ALWAYS_INLINE dealii::Tensor<2, dim, Number>
                                           convert_to_tensor(const dealii::SymmetricTensor<2, dim, Number> &s)
    {
      dealii::Tensor<2, dim, Number> t;

      // diagonal entries are stored first
      for (unsigned int d = 0; d < dim; ++d)
        t[d][d] = s.access_raw_entry(d);

      // off-diagonal entries come next, row by row
      for (unsigned int d = 0, c = 0; d < dim; ++d)
        for (unsigned int e = d + 1; e < dim; ++e, ++c)
          {
            t[d][e] = s.access_raw_entry(dim + c);
            t[e][d] = s.access_raw_entry(dim + c);
          }
      return t;
    }


    template <int dim, typename Number>
    constexpr dealii::Tensor<4, dim, Number>
    convert_to_tensor(const dealii::SymmetricTensor<4, dim, Number> &st)
    {
      // utilize the symmetry properties of SymmetricTensor<4,dim>
      // discussed in the class documentation to avoid accessing all
      // independent elements of the input tensor more than once
      dealii::Tensor<4, dim, Number> t;

      for (unsigned int i = 0; i < dim; ++i)
        for (unsigned int j = i; j < dim; ++j)
          for (unsigned int k = 0; k < dim; ++k)
            for (unsigned int l = k; l < dim; ++l)
              t[TableIndices<4>(i, j, k, l)] = t[TableIndices<4>(i, j, l, k)] =
                t[TableIndices<4>(j, i, k, l)] =
                  t[TableIndices<4>(j, i, l, k)] =
                    st[TableIndices<4>(i, j, k, l)];

      return t;
    }


    template <typename Number>
    struct Inverse<2, 1, Number>
    {
      constexpr static inline DEAL_II_ALWAYS_INLINE
        dealii::SymmetricTensor<2, 1, Number>
        value(const dealii::SymmetricTensor<2, 1, Number> &t)
      {
        dealii::SymmetricTensor<2, 1, Number> tmp;

        tmp[0][0] = 1.0 / t[0][0];

        return tmp;
      }
    };


    template <typename Number>
    struct Inverse<2, 2, Number>
    {
      constexpr static inline DEAL_II_ALWAYS_INLINE
        dealii::SymmetricTensor<2, 2, Number>
        value(const dealii::SymmetricTensor<2, 2, Number> &t)
      {
        dealii::SymmetricTensor<2, 2, Number> tmp;

        // Sympy result: ([
        // [ t11/(t00*t11 - t01**2), -t01/(t00*t11 - t01**2)],
        // [-t01/(t00*t11 - t01**2),  t00/(t00*t11 - t01**2)]  ])
        const TableIndices<2> idx_00(0, 0);
        const TableIndices<2> idx_01(0, 1);
        const TableIndices<2> idx_11(1, 1);
        const Number          inv_det_t =
          1.0 / (t[idx_00] * t[idx_11] - t[idx_01] * t[idx_01]);
        tmp[idx_00] = t[idx_11];
        tmp[idx_01] = -t[idx_01];
        tmp[idx_11] = t[idx_00];
        tmp *= inv_det_t;

        return tmp;
      }
    };


    template <typename Number>
    struct Inverse<2, 3, Number>
    {
      constexpr static dealii::SymmetricTensor<2, 3, Number>
      value(const dealii::SymmetricTensor<2, 3, Number> &t)
      {
        dealii::SymmetricTensor<2, 3, Number> tmp;

        // Sympy result: ([
        // [  (t11*t22 - t12**2)/(t00*t11*t22 - t00*t12**2 - t01**2*t22 +
        //                        2*t01*t02*t12 - t02**2*t11),
        //    (-t01*t22 + t02*t12)/(t00*t11*t22 - t00*t12**2 - t01**2*t22 +
        //                          2*t01*t02*t12 - t02**2*t11),
        //    (t01*t12 - t02*t11)/(t00*t11*t22 -  t00*t12**2 - t01**2*t22 +
        //                         2*t01*t02*t12 - t02**2*t11)],
        // [  (-t01*t22 + t02*t12)/(t00*t11*t22 - t00*t12**2 - t01**2*t22 +
        //                          2*t01*t02*t12 - t02**2*t11),
        //    (t00*t22 - t02**2)/(t00*t11*t22 - t00*t12**2 - t01**2*t22 +
        //                        2*t01*t02*t12 - t02**2*t11),
        //    (t00*t12 - t01*t02)/(-t00*t11*t22 + t00*t12**2 + t01**2*t22 -
        //                         2*t01*t02*t12 + t02**2*t11)],
        // [  (t01*t12 - t02*t11)/(t00*t11*t22 - t00*t12**2 - t01**2*t22 +
        //                         2*t01*t02*t12 - t02**2*t11),
        //    (t00*t12 - t01*t02)/(-t00*t11*t22 + t00*t12**2 + t01**2*t22 -
        //                         2*t01*t02*t12 + t02**2*t11),
        //    (-t00*t11 + t01**2)/(-t00*t11*t22 + t00*t12**2 + t01**2*t22 -
        //                         2*t01*t02*t12 + t02**2*t11)]  ])
        //
        // =
        //
        // [  (t11*t22 - t12**2)/det_t,
        //    (-t01*t22 + t02*t12)/det_t,
        //    (t01*t12 - t02*t11)/det_t],
        // [  (-t01*t22 + t02*t12)/det_t,
        //    (t00*t22 - t02**2)/det_t,
        //    (-t00*t12 + t01*t02)/det_t],
        // [  (t01*t12 - t02*t11)/det_t,
        //    (-t00*t12 + t01*t02)/det_t,
        //    (t00*t11 - t01**2)/det_t]   ])
        //
        // with det_t = (t00*t11*t22 - t00*t12**2 - t01**2*t22 +
        //               2*t01*t02*t12 - t02**2*t11)
        const TableIndices<2> idx_00(0, 0);
        const TableIndices<2> idx_01(0, 1);
        const TableIndices<2> idx_02(0, 2);
        const TableIndices<2> idx_11(1, 1);
        const TableIndices<2> idx_12(1, 2);
        const TableIndices<2> idx_22(2, 2);
        const Number          inv_det_t =
          1.0 / (t[idx_00] * t[idx_11] * t[idx_22] -
                 t[idx_00] * t[idx_12] * t[idx_12] -
                 t[idx_01] * t[idx_01] * t[idx_22] +
                 2.0 * t[idx_01] * t[idx_02] * t[idx_12] -
                 t[idx_02] * t[idx_02] * t[idx_11]);
        tmp[idx_00] = t[idx_11] * t[idx_22] - t[idx_12] * t[idx_12];
        tmp[idx_01] = -t[idx_01] * t[idx_22] + t[idx_02] * t[idx_12];
        tmp[idx_02] = t[idx_01] * t[idx_12] - t[idx_02] * t[idx_11];
        tmp[idx_11] = t[idx_00] * t[idx_22] - t[idx_02] * t[idx_02];
        tmp[idx_12] = -t[idx_00] * t[idx_12] + t[idx_01] * t[idx_02];
        tmp[idx_22] = t[idx_00] * t[idx_11] - t[idx_01] * t[idx_01];
        tmp *= inv_det_t;

        return tmp;
      }
    };


    template <typename Number>
    struct Inverse<4, 1, Number>
    {
      constexpr static inline dealii::SymmetricTensor<4, 1, Number>
      value(const dealii::SymmetricTensor<4, 1, Number> &t)
      {
        dealii::SymmetricTensor<4, 1, Number> tmp;
        tmp.data[0][0] = 1.0 / t.data[0][0];
        return tmp;
      }
    };


    template <typename Number>
    struct Inverse<4, 2, Number>
    {
      constexpr static inline dealii::SymmetricTensor<4, 2, Number>
      value(const dealii::SymmetricTensor<4, 2, Number> &t)
      {
        dealii::SymmetricTensor<4, 2, Number> tmp;

        // Inverting this tensor is a little more complicated than necessary,
        // since we store the data of 't' as a 3x3 matrix t.data, but the
        // product between a rank-4 and a rank-2 tensor is really not the
        // product between this matrix and the 3-vector of a rhs, but rather
        //
        // B.vec = t.data * mult * A.vec
        //
        // where mult is a 3x3 matrix with entries [[1,0,0],[0,1,0],[0,0,2]] to
        // capture the fact that we need to add up both the c_ij12*a_12 and the
        // c_ij21*a_21 terms.
        //
        // In addition, in this scheme, the identity tensor has the matrix
        // representation mult^-1.
        //
        // The inverse of 't' therefore has the matrix representation
        //
        // inv.data = mult^-1 * t.data^-1 * mult^-1
        //
        // in order to compute it, let's first compute the inverse of t.data and
        // put it into tmp.data; at the end of the function we then scale the
        // last row and column of the inverse by 1/2, corresponding to the left
        // and right multiplication with mult^-1.
        const Number t4  = t.data[0][0] * t.data[1][1],
                     t6  = t.data[0][0] * t.data[1][2],
                     t8  = t.data[0][1] * t.data[1][0],
                     t00 = t.data[0][2] * t.data[1][0],
                     t01 = t.data[0][1] * t.data[2][0],
                     t04 = t.data[0][2] * t.data[2][0],
                     t07 = 1.0 / (t4 * t.data[2][2] - t6 * t.data[2][1] -
                                  t8 * t.data[2][2] + t00 * t.data[2][1] +
                                  t01 * t.data[1][2] - t04 * t.data[1][1]);
        tmp.data[0][0] =
          (t.data[1][1] * t.data[2][2] - t.data[1][2] * t.data[2][1]) * t07;
        tmp.data[0][1] =
          -(t.data[0][1] * t.data[2][2] - t.data[0][2] * t.data[2][1]) * t07;
        tmp.data[0][2] =
          -(-t.data[0][1] * t.data[1][2] + t.data[0][2] * t.data[1][1]) * t07;
        tmp.data[1][0] =
          -(t.data[1][0] * t.data[2][2] - t.data[1][2] * t.data[2][0]) * t07;
        tmp.data[1][1] = (t.data[0][0] * t.data[2][2] - t04) * t07;
        tmp.data[1][2] = -(t6 - t00) * t07;
        tmp.data[2][0] =
          -(-t.data[1][0] * t.data[2][1] + t.data[1][1] * t.data[2][0]) * t07;
        tmp.data[2][1] = -(t.data[0][0] * t.data[2][1] - t01) * t07;
        tmp.data[2][2] = (t4 - t8) * t07;

        // scale last row and column as mentioned
        // above
        tmp.data[2][0] /= 2;
        tmp.data[2][1] /= 2;
        tmp.data[0][2] /= 2;
        tmp.data[1][2] /= 2;
        tmp.data[2][2] /= 4;

        return tmp;
      }
    };


    template <typename Number>
    struct Inverse<4, 3, Number>
    {
      static dealii::SymmetricTensor<4, 3, Number>
      value(const dealii::SymmetricTensor<4, 3, Number> &t)
      {
        dealii::SymmetricTensor<4, 3, Number> tmp = t;

        // This function follows the exact same scheme as the 2d case, except
        // that hardcoding the inverse of a 6x6 matrix is pretty wasteful.
        // Instead, we use the Gauss-Jordan algorithm implemented for
        // FullMatrix. For historical reasons the following code is copied from
        // there, with the tangential benefit that we do not need to copy the
        // tensor entries to and from the FullMatrix.
        const unsigned int N = 6;

        // First get an estimate of the size of the elements of this matrix,
        // for later checks whether the pivot element is large enough, or
        // whether we have to fear that the matrix is not regular.
        Number diagonal_sum = internal::NumberType<Number>::value(0.0);
        for (unsigned int i = 0; i < N; ++i)
          diagonal_sum += std::fabs(tmp.data[i][i]);
        const Number typical_diagonal_element =
          diagonal_sum / static_cast<double>(N);
        (void)typical_diagonal_element;

        unsigned int p[N];
        for (unsigned int i = 0; i < N; ++i)
          p[i] = i;

        for (unsigned int j = 0; j < N; ++j)
          {
            // Pivot search: search that part of the line on and right of the
            // diagonal for the largest element.
            Number       max = std::fabs(tmp.data[j][j]);
            unsigned int r   = j;
            for (unsigned int i = j + 1; i < N; ++i)
              if (std::fabs(tmp.data[i][j]) > max)
                {
                  max = std::fabs(tmp.data[i][j]);
                  r   = i;
                }

            // Check whether the pivot is too small
            Assert(max > 1.e-16 * typical_diagonal_element,
                   ExcMessage("This tensor seems to be noninvertible"));

            // Row interchange
            if (r > j)
              {
                for (unsigned int k = 0; k < N; ++k)
                  std::swap(tmp.data[j][k], tmp.data[r][k]);

                std::swap(p[j], p[r]);
              }

            // Transformation
            const Number hr = 1. / tmp.data[j][j];
            tmp.data[j][j]  = hr;
            for (unsigned int k = 0; k < N; ++k)
              {
                if (k == j)
                  continue;
                for (unsigned int i = 0; i < N; ++i)
                  {
                    if (i == j)
                      continue;
                    tmp.data[i][k] -= tmp.data[i][j] * tmp.data[j][k] * hr;
                  }
              }
            for (unsigned int i = 0; i < N; ++i)
              {
                tmp.data[i][j] *= hr;
                tmp.data[j][i] *= -hr;
              }
            tmp.data[j][j] = hr;
          }

        // Column interchange
        Number hv[N];
        for (unsigned int i = 0; i < N; ++i)
          {
            for (unsigned int k = 0; k < N; ++k)
              hv[p[k]] = tmp.data[i][k];
            for (unsigned int k = 0; k < N; ++k)
              tmp.data[i][k] = hv[k];
          }

        // Scale rows and columns. The mult matrix
        // here is diag[1, 1, 1, 1/2, 1/2, 1/2].
        for (unsigned int i = 3; i < 6; ++i)
          for (unsigned int j = 0; j < 3; ++j)
            tmp.data[i][j] /= 2;

        for (unsigned int i = 0; i < 3; ++i)
          for (unsigned int j = 3; j < 6; ++j)
            tmp.data[i][j] /= 2;

        for (unsigned int i = 3; i < 6; ++i)
          for (unsigned int j = 3; j < 6; ++j)
            tmp.data[i][j] /= 4;

        return tmp;
      }
    };

  } // namespace SymmetricTensorImplementation
} // namespace internal



template <int rank_, int dim, typename Number>
constexpr DEAL_II_ALWAYS_INLINE SymmetricTensor<rank_, dim, Number>::
                                operator Tensor<rank_, dim, Number>() const
{
  return internal::SymmetricTensorImplementation::convert_to_tensor(*this);
}



template <int rank_, int dim, typename Number>
constexpr bool
SymmetricTensor<rank_, dim, Number>::
operator==(const SymmetricTensor<rank_, dim, Number> &t) const
{
  return data == t.data;
}



template <int rank_, int dim, typename Number>
constexpr bool
SymmetricTensor<rank_, dim, Number>::
operator!=(const SymmetricTensor<rank_, dim, Number> &t) const
{
  return data != t.data;
}



template <int rank_, int dim, typename Number>
template <typename OtherNumber>
constexpr inline DEAL_II_ALWAYS_INLINE SymmetricTensor<rank_, dim, Number> &
                                       SymmetricTensor<rank_, dim, Number>::
                                       operator+=(const SymmetricTensor<rank_, dim, OtherNumber> &t)
{
  data += t.data;
  return *this;
}



template <int rank_, int dim, typename Number>
template <typename OtherNumber>
constexpr inline DEAL_II_ALWAYS_INLINE SymmetricTensor<rank_, dim, Number> &
                                       SymmetricTensor<rank_, dim, Number>::
                                       operator-=(const SymmetricTensor<rank_, dim, OtherNumber> &t)
{
  data -= t.data;
  return *this;
}



template <int rank_, int dim, typename Number>
template <typename OtherNumber>
constexpr inline DEAL_II_ALWAYS_INLINE SymmetricTensor<rank_, dim, Number> &
SymmetricTensor<rank_, dim, Number>::operator*=(const OtherNumber &d)
{
  data *= d;
  return *this;
}



template <int rank_, int dim, typename Number>
template <typename OtherNumber>
constexpr inline DEAL_II_ALWAYS_INLINE SymmetricTensor<rank_, dim, Number> &
SymmetricTensor<rank_, dim, Number>::operator/=(const OtherNumber &d)
{
  data /= d;
  return *this;
}



template <int rank_, int dim, typename Number>
constexpr inline DEAL_II_ALWAYS_INLINE SymmetricTensor<rank_, dim, Number>
SymmetricTensor<rank_, dim, Number>::operator-() const
{
  SymmetricTensor tmp = *this;
  tmp.data            = -tmp.data;
  return tmp;
}



template <int rank_, int dim, typename Number>
constexpr inline DEAL_II_ALWAYS_INLINE void
SymmetricTensor<rank_, dim, Number>::clear()
{
  data.clear();
}



template <int rank_, int dim, typename Number>
constexpr std::size_t
SymmetricTensor<rank_, dim, Number>::memory_consumption()
{
  // all memory consists of statically allocated memory of the current
  // object, no pointers
  return sizeof(SymmetricTensor<rank_, dim, Number>);
}



namespace internal
{
  template <int dim, typename Number, typename OtherNumber = Number>
  DEAL_II_CONSTEXPR inline DEAL_II_ALWAYS_INLINE
    typename SymmetricTensorAccessors::
      double_contraction_result<2, 2, dim, Number, OtherNumber>::type
      perform_double_contraction(
        const typename SymmetricTensorAccessors::StorageType<2, dim, Number>::
          base_tensor_type &data,
        const typename SymmetricTensorAccessors::
          StorageType<2, dim, OtherNumber>::base_tensor_type &sdata)
  {
    using result_type = typename SymmetricTensorAccessors::
      double_contraction_result<2, 2, dim, Number, OtherNumber>::type;

    switch (dim)
      {
        case 1:
          return data[0] * sdata[0];
        default:
          // Start with the non-diagonal part to avoid some multiplications by
          // 2.

          result_type sum = data[dim] * sdata[dim];
          for (unsigned int d = dim + 1; d < (dim * (dim + 1) / 2); ++d)
            sum += data[d] * sdata[d];
          sum += sum; // sum = sum * 2.;

          // Now add the contributions from the diagonal
          for (unsigned int d = 0; d < dim; ++d)
            sum += data[d] * sdata[d];
          return sum;
      }
  }



  template <int dim, typename Number, typename OtherNumber = Number>
  DEAL_II_CONSTEXPR inline DEAL_II_ALWAYS_INLINE
    typename SymmetricTensorAccessors::
      double_contraction_result<4, 2, dim, Number, OtherNumber>::type
      perform_double_contraction(
        const typename SymmetricTensorAccessors::StorageType<4, dim, Number>::
          base_tensor_type &data,
        const typename SymmetricTensorAccessors::
          StorageType<2, dim, OtherNumber>::base_tensor_type &sdata)
  {
    using result_type = typename SymmetricTensorAccessors::
      double_contraction_result<4, 2, dim, Number, OtherNumber>::type;
    using value_type = typename SymmetricTensorAccessors::
      double_contraction_result<4, 2, dim, Number, OtherNumber>::value_type;

    const unsigned int data_dim = SymmetricTensorAccessors::
      StorageType<2, dim, value_type>::n_independent_components;
    value_type tmp[data_dim]{};
    for (unsigned int i = 0; i < data_dim; ++i)
      tmp[i] =
        perform_double_contraction<dim, Number, OtherNumber>(data[i], sdata);
    return result_type(tmp);
  }



  template <int dim, typename Number, typename OtherNumber = Number>
  DEAL_II_CONSTEXPR inline DEAL_II_ALWAYS_INLINE
    typename SymmetricTensorAccessors::StorageType<
      2,
      dim,
      typename SymmetricTensorAccessors::
        double_contraction_result<2, 4, dim, Number, OtherNumber>::value_type>::
      base_tensor_type
      perform_double_contraction(
        const typename SymmetricTensorAccessors::StorageType<2, dim, Number>::
          base_tensor_type &data,
        const typename SymmetricTensorAccessors::
          StorageType<4, dim, OtherNumber>::base_tensor_type &sdata)
  {
    using value_type = typename SymmetricTensorAccessors::
      double_contraction_result<2, 4, dim, Number, OtherNumber>::value_type;
    using base_tensor_type = typename SymmetricTensorAccessors::
      StorageType<2, dim, value_type>::base_tensor_type;

    base_tensor_type tmp;
    for (unsigned int i = 0; i < tmp.dimension; ++i)
      {
        // Start with the non-diagonal part
        value_type sum = data[dim] * sdata[dim][i];
        for (unsigned int d = dim + 1; d < (dim * (dim + 1) / 2); ++d)
          sum += data[d] * sdata[d][i];
        sum += sum; // sum = sum * 2.;

        // Now add the contributions from the diagonal
        for (unsigned int d = 0; d < dim; ++d)
          sum += data[d] * sdata[d][i];
        tmp[i] = sum;
      }
    return tmp;
  }



  template <int dim, typename Number, typename OtherNumber = Number>
  DEAL_II_CONSTEXPR inline DEAL_II_ALWAYS_INLINE
    typename SymmetricTensorAccessors::StorageType<
      4,
      dim,
      typename SymmetricTensorAccessors::
        double_contraction_result<4, 4, dim, Number, OtherNumber>::value_type>::
      base_tensor_type
      perform_double_contraction(
        const typename SymmetricTensorAccessors::StorageType<4, dim, Number>::
          base_tensor_type &data,
        const typename SymmetricTensorAccessors::
          StorageType<4, dim, OtherNumber>::base_tensor_type &sdata)
  {
    using value_type = typename SymmetricTensorAccessors::
      double_contraction_result<4, 4, dim, Number, OtherNumber>::value_type;
    using base_tensor_type = typename SymmetricTensorAccessors::
      StorageType<4, dim, value_type>::base_tensor_type;

    const unsigned int data_dim = SymmetricTensorAccessors::
      StorageType<2, dim, value_type>::n_independent_components;
    base_tensor_type tmp;
    for (unsigned int i = 0; i < data_dim; ++i)
      for (unsigned int j = 0; j < data_dim; ++j)
        {
          // Start with the non-diagonal part
          for (unsigned int d = dim; d < (dim * (dim + 1) / 2); ++d)
            tmp[i][j] += data[i][d] * sdata[d][j];
          tmp[i][j] += tmp[i][j]; // tmp[i][j] = tmp[i][j] * 2;

          // Now add the contributions from the diagonal
          for (unsigned int d = 0; d < dim; ++d)
            tmp[i][j] += data[i][d] * sdata[d][j];
        }
    return tmp;
  }

} // end of namespace internal



template <int rank_, int dim, typename Number>
template <typename OtherNumber>
DEAL_II_CONSTEXPR inline DEAL_II_ALWAYS_INLINE
  typename internal::SymmetricTensorAccessors::
    double_contraction_result<rank_, 2, dim, Number, OtherNumber>::type
      SymmetricTensor<rank_, dim, Number>::
      operator*(const SymmetricTensor<2, dim, OtherNumber> &s) const
{
  // need to have two different function calls
  // because a scalar and rank-2 tensor are not
  // the same data type (see internal function
  // above)
  return internal::perform_double_contraction<dim, Number, OtherNumber>(data,
                                                                        s.data);
}



template <int rank_, int dim, typename Number>
template <typename OtherNumber>
DEAL_II_CONSTEXPR inline typename internal::SymmetricTensorAccessors::
  double_contraction_result<rank_, 4, dim, Number, OtherNumber>::type
    SymmetricTensor<rank_, dim, Number>::
    operator*(const SymmetricTensor<4, dim, OtherNumber> &s) const
{
  typename internal::SymmetricTensorAccessors::
    double_contraction_result<rank_, 4, dim, Number, OtherNumber>::type tmp;
  tmp.data =
    internal::perform_double_contraction<dim, Number, OtherNumber>(data,
                                                                   s.data);
  return tmp;
}



// internal namespace to switch between the
// access of different tensors. There used to
// be explicit instantiations before for
// different ranks and dimensions, but since
// we now allow for templates on the data
// type, and since we cannot partially
// specialize the implementation, this got
// into a separate namespace
namespace internal
{
  // The variables within this struct will be referenced in the next functions.
  // It is a workaround that allows returning a reference to a static variable
  // while allowing constexpr evaluation of the function.
  // It has to be defined outside the function because constexpr functions
  // cannot define static variables.
  // A similar struct has also been defined in tensor.h
  template <typename Type>
  struct Uninitialized
  {
    static Type value;
  };

  template <typename Type>
  Type Uninitialized<Type>::value;

  template <int dim, typename Number>
  constexpr inline DEAL_II_ALWAYS_INLINE Number &
                                         symmetric_tensor_access(const TableIndices<2> &indices,
                                                                 typename SymmetricTensorAccessors::
                                                                   StorageType<2, dim, Number>::base_tensor_type &data)
  {
    // 1d is very simple and done first
    if (dim == 1)
      return data[0];

    // first treat the main diagonal elements, which are stored consecutively
    // at the beginning
    if (indices[0] == indices[1])
      return data[indices[0]];

    // the rest is messier and requires a few switches.
    switch (dim)
      {
        case 2:
          // at least for the 2x2 case it is reasonably simple
          Assert(((indices[0] == 1) && (indices[1] == 0)) ||
                   ((indices[0] == 0) && (indices[1] == 1)),
                 ExcInternalError());
          return data[2];

        default:
          // to do the rest, sort our indices before comparing
          {
            TableIndices<2> sorted_indices(std::min(indices[0], indices[1]),
                                           std::max(indices[0], indices[1]));
            for (unsigned int d = 0, c = 0; d < dim; ++d)
              for (unsigned int e = d + 1; e < dim; ++e, ++c)
                if ((sorted_indices[0] == d) && (sorted_indices[1] == e))
                  return data[dim + c];
            Assert(false, ExcInternalError());
          }
      }

    // The code should never reach there.
    // Returns a dummy reference to a dummy variable just to make the
    // compiler happy.
    return Uninitialized<Number>::value;
  }



  template <int dim, typename Number>
  constexpr inline DEAL_II_ALWAYS_INLINE const Number &
                                               symmetric_tensor_access(const TableIndices<2> &indices,
                                                                       const typename SymmetricTensorAccessors::
                                                                         StorageType<2, dim, Number>::base_tensor_type &data)
  {
    // 1d is very simple and done first
    if (dim == 1)
      return data[0];

    // first treat the main diagonal elements, which are stored consecutively
    // at the beginning
    if (indices[0] == indices[1])
      return data[indices[0]];

    // the rest is messier and requires a few switches.
    switch (dim)
      {
        case 2:
          // at least for the 2x2 case it is reasonably simple
          Assert(((indices[0] == 1) && (indices[1] == 0)) ||
                   ((indices[0] == 0) && (indices[1] == 1)),
                 ExcInternalError());
          return data[2];

        default:
          // to do the rest, sort our indices before comparing
          {
            TableIndices<2> sorted_indices(std::min(indices[0], indices[1]),
                                           std::max(indices[0], indices[1]));
            for (unsigned int d = 0, c = 0; d < dim; ++d)
              for (unsigned int e = d + 1; e < dim; ++e, ++c)
                if ((sorted_indices[0] == d) && (sorted_indices[1] == e))
                  return data[dim + c];
            Assert(false, ExcInternalError());
          }
      }

    // The code should never reach there.
    // Returns a dummy reference to a dummy variable just to make the
    // compiler happy.
    return Uninitialized<Number>::value;
  }



  template <int dim, typename Number>
  constexpr inline Number &
  symmetric_tensor_access(const TableIndices<4> &indices,
                          typename SymmetricTensorAccessors::
                            StorageType<4, dim, Number>::base_tensor_type &data)
  {
    switch (dim)
      {
        case 1:
          return data[0][0];

        case 2:
          // each entry of the tensor can be thought of as an entry in a
          // matrix that maps the rolled-out rank-2 tensors into rolled-out
          // rank-2 tensors. this is the format in which we store rank-4
          // tensors. determine which position the present entry is
          // stored in
          {
            constexpr std::size_t base_index[2][2] = {{0, 2}, {2, 1}};
            return data[base_index[indices[0]][indices[1]]]
                       [base_index[indices[2]][indices[3]]];
          }
        case 3:
          // each entry of the tensor can be thought of as an entry in a
          // matrix that maps the rolled-out rank-2 tensors into rolled-out
          // rank-2 tensors. this is the format in which we store rank-4
          // tensors. determine which position the present entry is
          // stored in
          {
            constexpr std::size_t base_index[3][3] = {{0, 3, 4},
                                                      {3, 1, 5},
                                                      {4, 5, 2}};
            return data[base_index[indices[0]][indices[1]]]
                       [base_index[indices[2]][indices[3]]];
          }

        default:
          Assert(false, ExcNotImplemented());
      }

    // The code should never reach there.
    // Returns a dummy reference to a dummy variable just to make the
    // compiler happy.
    return Uninitialized<Number>::value;
  }


  template <int dim, typename Number>
  constexpr inline DEAL_II_ALWAYS_INLINE const Number &
                                               symmetric_tensor_access(const TableIndices<4> &indices,
                                                                       const typename SymmetricTensorAccessors::
                                                                         StorageType<4, dim, Number>::base_tensor_type &data)
  {
    switch (dim)
      {
        case 1:
          return data[0][0];

        case 2:
          // each entry of the tensor can be thought of as an entry in a
          // matrix that maps the rolled-out rank-2 tensors into rolled-out
          // rank-2 tensors. this is the format in which we store rank-4
          // tensors. determine which position the present entry is
          // stored in
          {
            constexpr std::size_t base_index[2][2] = {{0, 2}, {2, 1}};
            return data[base_index[indices[0]][indices[1]]]
                       [base_index[indices[2]][indices[3]]];
          }
        case 3:
          // each entry of the tensor can be thought of as an entry in a
          // matrix that maps the rolled-out rank-2 tensors into rolled-out
          // rank-2 tensors. this is the format in which we store rank-4
          // tensors. determine which position the present entry is
          // stored in
          {
            constexpr std::size_t base_index[3][3] = {{0, 3, 4},
                                                      {3, 1, 5},
                                                      {4, 5, 2}};
            return data[base_index[indices[0]][indices[1]]]
                       [base_index[indices[2]][indices[3]]];
          }

        default:
          Assert(false, ExcNotImplemented());
      }

    // The code should never reach there.
    // Returns a dummy reference to a dummy variable just to make the
    // compiler happy.
    return Uninitialized<Number>::value;
  }

} // end of namespace internal



template <int rank_, int dim, typename Number>
constexpr inline DEAL_II_ALWAYS_INLINE Number &
                                       SymmetricTensor<rank_, dim, Number>::
                                       operator()(const TableIndices<rank_> &indices)
{
  for (unsigned int r = 0; r < rank; ++r)
    AssertIndexRange(indices[r], dimension);
  return internal::symmetric_tensor_access<dim, Number>(indices, data);
}



template <int rank_, int dim, typename Number>
constexpr inline DEAL_II_ALWAYS_INLINE const Number &
                                             SymmetricTensor<rank_, dim, Number>::
                                             operator()(const TableIndices<rank_> &indices) const
{
  for (unsigned int r = 0; r < rank; ++r)
    AssertIndexRange(indices[r], dimension);
  return internal::symmetric_tensor_access<dim, Number>(indices, data);
}



namespace internal
{
  namespace SymmetricTensorImplementation
  {
    template <int rank_>
    constexpr TableIndices<rank_>
    get_partially_filled_indices(const unsigned int row,
                                 const std::integral_constant<int, 2> &)
    {
      return TableIndices<rank_>(row, numbers::invalid_unsigned_int);
    }


    template <int rank_>
    constexpr TableIndices<rank_>
    get_partially_filled_indices(const unsigned int row,
                                 const std::integral_constant<int, 4> &)
    {
      return TableIndices<rank_>(row,
                                 numbers::invalid_unsigned_int,
                                 numbers::invalid_unsigned_int,
                                 numbers::invalid_unsigned_int);
    }
  } // namespace SymmetricTensorImplementation
} // namespace internal


template <int rank_, int dim, typename Number>
constexpr DEAL_II_ALWAYS_INLINE internal::SymmetricTensorAccessors::
  Accessor<rank_, dim, true, rank_ - 1, Number>
    SymmetricTensor<rank_, dim, Number>::
    operator[](const unsigned int row) const
{
  return internal::SymmetricTensorAccessors::
    Accessor<rank_, dim, true, rank_ - 1, Number>(
      *this,
      internal::SymmetricTensorImplementation::get_partially_filled_indices<
        rank_>(row, std::integral_constant<int, rank_>()));
}



template <int rank_, int dim, typename Number>
constexpr inline DEAL_II_ALWAYS_INLINE internal::SymmetricTensorAccessors::
  Accessor<rank_, dim, false, rank_ - 1, Number>
    SymmetricTensor<rank_, dim, Number>::operator[](const unsigned int row)
{
  return internal::SymmetricTensorAccessors::
    Accessor<rank_, dim, false, rank_ - 1, Number>(
      *this,
      internal::SymmetricTensorImplementation::get_partially_filled_indices<
        rank_>(row, std::integral_constant<int, rank_>()));
}



template <int rank_, int dim, typename Number>
constexpr DEAL_II_ALWAYS_INLINE const Number &
                                      SymmetricTensor<rank_, dim, Number>::
                                      operator[](const TableIndices<rank_> &indices) const
{
  return operator()(indices);
}



template <int rank_, int dim, typename Number>
constexpr inline DEAL_II_ALWAYS_INLINE Number &
                                       SymmetricTensor<rank_, dim, Number>::
                                       operator[](const TableIndices<rank_> &indices)
{
  return operator()(indices);
}



template <int rank_, int dim, typename Number>
inline Number *
SymmetricTensor<rank_, dim, Number>::begin_raw()
{
  return std::addressof(this->access_raw_entry(0));
}



template <int rank_, int dim, typename Number>
inline const Number *
SymmetricTensor<rank_, dim, Number>::begin_raw() const
{
  return std::addressof(this->access_raw_entry(0));
}



template <int rank_, int dim, typename Number>
inline Number *
SymmetricTensor<rank_, dim, Number>::end_raw()
{
  return begin_raw() + n_independent_components;
}



template <int rank_, int dim, typename Number>
inline const Number *
SymmetricTensor<rank_, dim, Number>::end_raw() const
{
  return begin_raw() + n_independent_components;
}



namespace internal
{
  namespace SymmetricTensorImplementation
  {
    template <int dim, typename Number>
    constexpr unsigned int
    entry_to_indices(const dealii::SymmetricTensor<2, dim, Number> &,
                     const unsigned int index)
    {
      return index;
    }


    template <int dim, typename Number>
    constexpr dealii::TableIndices<2>
    entry_to_indices(const dealii::SymmetricTensor<4, dim, Number> &,
                     const unsigned int index)
    {
      return internal::SymmetricTensorAccessors::StorageType<4, dim, Number>::
        base_tensor_type::unrolled_to_component_indices(index);
    }

  } // namespace SymmetricTensorImplementation
} // namespace internal



template <int rank_, int dim, typename Number>
constexpr inline const Number &
SymmetricTensor<rank_, dim, Number>::access_raw_entry(
  const unsigned int index) const
{
  AssertIndexRange(index, n_independent_components);
  return data[internal::SymmetricTensorImplementation::entry_to_indices(*this,
                                                                        index)];
}



template <int rank_, int dim, typename Number>
constexpr inline Number &
SymmetricTensor<rank_, dim, Number>::access_raw_entry(const unsigned int index)
{
  AssertIndexRange(index, n_independent_components);
  return data[internal::SymmetricTensorImplementation::entry_to_indices(*this,
                                                                        index)];
}



namespace internal
{
  template <int dim, typename Number>
  constexpr inline typename numbers::NumberTraits<Number>::real_type
  compute_norm(const typename SymmetricTensorAccessors::
                 StorageType<2, dim, Number>::base_tensor_type &data)
  {
    switch (dim)
      {
        case 1:
          return numbers::NumberTraits<Number>::abs(data[0]);

        case 2:
          return std::sqrt(
            numbers::NumberTraits<Number>::abs_square(data[0]) +
            numbers::NumberTraits<Number>::abs_square(data[1]) +
            2. * numbers::NumberTraits<Number>::abs_square(data[2]));

        case 3:
          return std::sqrt(
            numbers::NumberTraits<Number>::abs_square(data[0]) +
            numbers::NumberTraits<Number>::abs_square(data[1]) +
            numbers::NumberTraits<Number>::abs_square(data[2]) +
            2. * numbers::NumberTraits<Number>::abs_square(data[3]) +
            2. * numbers::NumberTraits<Number>::abs_square(data[4]) +
            2. * numbers::NumberTraits<Number>::abs_square(data[5]));

        default:
          {
            typename numbers::NumberTraits<Number>::real_type return_value =
              typename numbers::NumberTraits<Number>::real_type();

            for (unsigned int d = 0; d < dim; ++d)
              return_value +=
                numbers::NumberTraits<Number>::abs_square(data[d]);
            for (unsigned int d = dim; d < (dim * dim + dim) / 2; ++d)
              return_value +=
                2. * numbers::NumberTraits<Number>::abs_square(data[d]);

            return std::sqrt(return_value);
          }
      }
  }



  template <int dim, typename Number>
  constexpr inline typename numbers::NumberTraits<Number>::real_type
  compute_norm(const typename SymmetricTensorAccessors::
                 StorageType<4, dim, Number>::base_tensor_type &data)
  {
    switch (dim)
      {
        case 1:
          return numbers::NumberTraits<Number>::abs(data[0][0]);

        default:
          {
            typename numbers::NumberTraits<Number>::real_type return_value =
              typename numbers::NumberTraits<Number>::real_type();

            const unsigned int n_independent_components = data.dimension;

            for (unsigned int i = 0; i < dim; ++i)
              for (unsigned int j = 0; j < dim; ++j)
                return_value +=
                  numbers::NumberTraits<Number>::abs_square(data[i][j]);
            for (unsigned int i = 0; i < dim; ++i)
              for (unsigned int j = dim; j < n_independent_components; ++j)
                return_value +=
                  2. * numbers::NumberTraits<Number>::abs_square(data[i][j]);
            for (unsigned int i = dim; i < n_independent_components; ++i)
              for (unsigned int j = 0; j < dim; ++j)
                return_value +=
                  2. * numbers::NumberTraits<Number>::abs_square(data[i][j]);
            for (unsigned int i = dim; i < n_independent_components; ++i)
              for (unsigned int j = dim; j < n_independent_components; ++j)
                return_value +=
                  4. * numbers::NumberTraits<Number>::abs_square(data[i][j]);

            return std::sqrt(return_value);
          }
      }
  }

} // end of namespace internal



template <int rank_, int dim, typename Number>
constexpr typename numbers::NumberTraits<Number>::real_type
SymmetricTensor<rank_, dim, Number>::norm() const
{
  return internal::compute_norm<dim, Number>(data);
}



namespace internal
{
  namespace SymmetricTensorImplementation
  {
    // a function to do the unrolling from a set of indices to a
    // scalar index into the array in which we store the elements of
    // a symmetric tensor
    //
    // this function is for rank-2 tensors
    template <int dim>
    constexpr inline DEAL_II_ALWAYS_INLINE unsigned int
    component_to_unrolled_index(const TableIndices<2> &indices)
    {
      AssertIndexRange(indices[0], dim);
      AssertIndexRange(indices[1], dim);

      switch (dim)
        {
          case 1:
            {
              return 0;
            }

          case 2:
            {
              constexpr unsigned int table[2][2] = {{0, 2}, {2, 1}};
              return table[indices[0]][indices[1]];
            }

          case 3:
            {
              constexpr unsigned int table[3][3] = {{0, 3, 4},
                                                    {3, 1, 5},
                                                    {4, 5, 2}};
              return table[indices[0]][indices[1]];
            }

          case 4:
            {
              constexpr unsigned int table[4][4] = {{0, 4, 5, 6},
                                                    {4, 1, 7, 8},
                                                    {5, 7, 2, 9},
                                                    {6, 8, 9, 3}};
              return table[indices[0]][indices[1]];
            }

          default:
            // for the remainder, manually figure out the numbering
            {
              if (indices[0] == indices[1])
                return indices[0];

              TableIndices<2> sorted_indices(indices);
              sorted_indices.sort();

              for (unsigned int d = 0, c = 0; d < dim; ++d)
                for (unsigned int e = d + 1; e < dim; ++e, ++c)
                  if ((sorted_indices[0] == d) && (sorted_indices[1] == e))
                    return dim + c;

              // should never get here:
              Assert(false, ExcInternalError());
              return 0;
            }
        }
    }

    // a function to do the unrolling from a set of indices to a
    // scalar index into the array in which we store the elements of
    // a symmetric tensor
    //
    // this function is for tensors of ranks not already handled
    // above
    template <int dim, int rank_>
    constexpr inline unsigned int
    component_to_unrolled_index(const TableIndices<rank_> &indices)
    {
      (void)indices;
      Assert(false, ExcNotImplemented());
      return numbers::invalid_unsigned_int;
    }
  } // namespace SymmetricTensorImplementation
} // namespace internal


template <int rank_, int dim, typename Number>
constexpr unsigned int
SymmetricTensor<rank_, dim, Number>::component_to_unrolled_index(
  const TableIndices<rank_> &indices)
{
  return internal::SymmetricTensorImplementation::component_to_unrolled_index<
    dim>(indices);
}



namespace internal
{
  namespace SymmetricTensorImplementation
  {
    // a function to do the inverse of the unrolling from a set of
    // indices to a scalar index into the array in which we store
    // the elements of a symmetric tensor. in other words, it goes
    // from the scalar index into the array to a set of indices of
    // the tensor
    //
    // this function is for rank-2 tensors
    template <int dim>
    constexpr inline DEAL_II_ALWAYS_INLINE TableIndices<2>
                                           unrolled_to_component_indices(const unsigned int i,
                                                                         const std::integral_constant<int, 2> &)
    {
      Assert(
        (i < dealii::SymmetricTensor<2, dim, double>::n_independent_components),
        ExcIndexRange(
          i,
          0,
          dealii::SymmetricTensor<2, dim, double>::n_independent_components));
      switch (dim)
        {
          case 1:
            {
              return {0, 0};
            }

          case 2:
            {
              const TableIndices<2> table[3] = {TableIndices<2>(0, 0),
                                                TableIndices<2>(1, 1),
                                                TableIndices<2>(0, 1)};
              return table[i];
            }

          case 3:
            {
              const TableIndices<2> table[6] = {TableIndices<2>(0, 0),
                                                TableIndices<2>(1, 1),
                                                TableIndices<2>(2, 2),
                                                TableIndices<2>(0, 1),
                                                TableIndices<2>(0, 2),
                                                TableIndices<2>(1, 2)};
              return table[i];
            }

          default:
            if (i < dim)
              return {i, i};

            for (unsigned int d = 0, c = dim; d < dim; ++d)
              for (unsigned int e = d + 1; e < dim; ++e, ++c)
                if (c == i)
                  return {d, e};

            // should never get here:
            Assert(false, ExcInternalError());
            return {0, 0};
        }
    }

    // a function to do the inverse of the unrolling from a set of
    // indices to a scalar index into the array in which we store
    // the elements of a symmetric tensor. in other words, it goes
    // from the scalar index into the array to a set of indices of
    // the tensor
    //
    // this function is for tensors of a rank not already handled
    // above
    template <int dim, int rank_>
    constexpr inline
      typename std::enable_if<rank_ != 2, TableIndices<rank_>>::type
      unrolled_to_component_indices(const unsigned int i,
                                    const std::integral_constant<int, rank_> &)
    {
      (void)i;
      Assert(
        (i <
         dealii::SymmetricTensor<rank_, dim, double>::n_independent_components),
        ExcIndexRange(i,
                      0,
                      dealii::SymmetricTensor<rank_, dim, double>::
                        n_independent_components));
      Assert(false, ExcNotImplemented());
      return TableIndices<rank_>();
    }

  } // namespace SymmetricTensorImplementation
} // namespace internal

template <int rank_, int dim, typename Number>
constexpr DEAL_II_ALWAYS_INLINE TableIndices<rank_>
                                SymmetricTensor<rank_, dim, Number>::unrolled_to_component_indices(
  const unsigned int i)
{
  return internal::SymmetricTensorImplementation::unrolled_to_component_indices<
    dim>(i, std::integral_constant<int, rank_>());
}



template <int rank_, int dim, typename Number>
template <class Archive>
inline void
SymmetricTensor<rank_, dim, Number>::serialize(Archive &ar, const unsigned int)
{
  ar &data;
}


#endif // DOXYGEN

 /* ----------------- Non-member functions operating on tensors. ------------ */ 


/**
 * 两个等级相同的对称张量相加。结果是另一个SymmetricTensor，它的数字类型与该操作兼容。
 * 如果可能的话（例如，当 @p Number 和 @p OtherNumber
 * 是同一类型时，或者如果 <code>Number() + OtherNumber()</code>
 * 的结果是另一个 @p Number),
 * ，你应该使用<tt>操作符+=</tt>来代替，因为这不需要创建一个临时变量。
 * @relatesalso  SymmetricTensor
 *
 *
 */
template <int rank_, int dim, typename Number, typename OtherNumber>
constexpr inline DEAL_II_ALWAYS_INLINE
  SymmetricTensor<rank_, dim, typename ProductType<Number, OtherNumber>::type>
  operator+(const SymmetricTensor<rank_, dim, Number> &     left,
            const SymmetricTensor<rank_, dim, OtherNumber> &right)
{
  SymmetricTensor<rank_, dim, typename ProductType<Number, OtherNumber>::type>
    tmp = left;
  tmp += right;
  return tmp;
}


/**
 * 减去两个等级相同的对称张量。其结果是另一个SymmetricTensor，它的数字类型与该操作兼容。
 * 如果可能的话（例如，当 @p Number 和 @p OtherNumber
 * 是相同的类型，或者如果<code>Number()的结果
 *
 * - OtherNumber()</code>是另一个 @p Number), ，你应该使用<tt>operator-=</tt>代替，因为这不需要创建一个临时变量。
 * @relatesalso  SymmetricTensor
 *
 *
 */
template <int rank_, int dim, typename Number, typename OtherNumber>
constexpr inline DEAL_II_ALWAYS_INLINE
  SymmetricTensor<rank_, dim, typename ProductType<Number, OtherNumber>::type>
  operator-(const SymmetricTensor<rank_, dim, Number> &     left,
            const SymmetricTensor<rank_, dim, OtherNumber> &right)
{
  SymmetricTensor<rank_, dim, typename ProductType<Number, OtherNumber>::type>
    tmp = left;
  tmp -= right;
  return tmp;
}


/**
 * 将一个对称张量和一个相同等级的一般张量相加。其结果是一个一般的张量，其数字类型与该操作兼容。
 * @relatesalso  SymmetricTensor
 *
 *
 */
template <int rank_, int dim, typename Number, typename OtherNumber>
constexpr DEAL_II_ALWAYS_INLINE
  Tensor<rank_, dim, typename ProductType<Number, OtherNumber>::type>
  operator+(const SymmetricTensor<rank_, dim, Number> &left,
            const Tensor<rank_, dim, OtherNumber> &    right)
{
  return Tensor<rank_, dim, Number>(left) + right;
}


/**
 * 一个普通张量与一个相同等级的对称张量相加。其结果是一个一般的张量，其数字类型与该操作兼容。
 * @relatesalso  SymmetricTensor
 *
 *
 */
template <int rank_, int dim, typename Number, typename OtherNumber>
constexpr DEAL_II_ALWAYS_INLINE
  Tensor<rank_, dim, typename ProductType<Number, OtherNumber>::type>
  operator+(const Tensor<rank_, dim, Number> &              left,
            const SymmetricTensor<rank_, dim, OtherNumber> &right)
{
  return left + Tensor<rank_, dim, OtherNumber>(right);
}


/**
 * 从一个相同等级的对称张量中减去一个一般张量。其结果是一个一般的张量，其数字类型与该操作兼容。
 * @relatesalso  SymmetricTensor
 *
 *
 */
template <int rank_, int dim, typename Number, typename OtherNumber>
constexpr DEAL_II_ALWAYS_INLINE
  Tensor<rank_, dim, typename ProductType<Number, OtherNumber>::type>
  operator-(const SymmetricTensor<rank_, dim, Number> &left,
            const Tensor<rank_, dim, OtherNumber> &    right)
{
  return Tensor<rank_, dim, Number>(left) - right;
}


/**
 * 从一个相同等级的一般张量中减去一个对称张量。其结果是一个一般的张量，其数字类型与该操作兼容。
 * @relatesalso  SymmetricTensor
 *
 *
 */
template <int rank_, int dim, typename Number, typename OtherNumber>
constexpr DEAL_II_ALWAYS_INLINE
  Tensor<rank_, dim, typename ProductType<Number, OtherNumber>::type>
  operator-(const Tensor<rank_, dim, Number> &              left,
            const SymmetricTensor<rank_, dim, OtherNumber> &right)
{
  return left - Tensor<rank_, dim, OtherNumber>(right);
}



/**
 * 计算一个等级2的对称张量的行列式。行列式通常也被称为等级2张量的第三不变量。
 * 对于一个一维张量，行列式等于唯一的元素，因此等同于轨迹。
 * 为了简化符号学，还有一个<tt>third_invariant()</tt>函数可以返回一个张量的行列式。
 * @relatesalso  SymmetricTensor
 *
 *
 */
template <int dim, typename Number>
constexpr inline DEAL_II_ALWAYS_INLINE Number
                                       determinant(const SymmetricTensor<2, dim, Number> &t)
{
  switch (dim)
    {
      case 1:
        return t.data[0];
      case 2:
        return (t.data[0] * t.data[1] - t.data[2] * t.data[2]);
      case 3:
        {
          // in analogy to general tensors, but
          // there's something to be simplified for
          // the present case
          const Number tmp = t.data[3] * t.data[4] * t.data[5];
          return (tmp + tmp + t.data[0] * t.data[1] * t.data[2] -
                  t.data[0] * t.data[5] * t.data[5] -
                  t.data[1] * t.data[4] * t.data[4] -
                  t.data[2] * t.data[3] * t.data[3]);
        }
      default:
        Assert(false, ExcNotImplemented());
        return internal::NumberType<Number>::value(0.0);
    }
}



/**
 * 计算一个等级2的对称张量的行列式。因此，这个函数计算的值与<tt>determinant()</tt>函数相同，提供这个函数只是为了在符号上更加简单（因为还有first_invariant()和second_invariant()函数）。\f[
 * I_3 (\mathbf A) = III (\mathbf A) = \det (\mathbf A) \f]
 * @relatesalso  SymmetricTensor
 *
 *
 */
template <int dim, typename Number>
constexpr DEAL_II_ALWAYS_INLINE Number
                                third_invariant(const SymmetricTensor<2, dim, Number> &t)
{
  return determinant(t);
}



/**
 * 计算并返回等级为2的张量的轨迹，即其对角线项的总和。跟踪是秩2张量的第一个不变量。\f[
 * \text{tr} \mathbf A = \sum_i A_{ii} \f]
 * @relatesalso  SymmetricTensor
 *
 *
 */
template <int dim, typename Number>
constexpr inline DEAL_II_ALWAYS_INLINE Number
                                       trace(const SymmetricTensor<2, dim, Number> &d)
{
  Number t = d.data[0];
  for (unsigned int i = 1; i < dim; ++i)
    t += d.data[i];
  return t;
}


/**
 * 计算一个秩2对称张量的迹。因此，这个函数计算的值与<tt>trace()</tt>函数相同，提供这个函数只是为了在符号上更加简单（因为还有second_invariant()和third_invariant()函数）。\f[
 * I_1 (\mathbf A) = I (\mathbf A) = \text{tr} \mathbf A = \sum_i A_{ii} \f]
 * @relatesalso  SymmetricTensor
 *
 *
 */
template <int dim, typename Number>
constexpr Number
first_invariant(const SymmetricTensor<2, dim, Number> &t)
{
  return trace(t);
}


/**
 * 计算一个等级为2的张量的第二个不变式。张量 $\mathbf A$
 * 的第二不变量定义为 $I_2 (\mathbf A) = II(\mathbf A) = \frac 12
 * \left[ (\text{tr} \mathbf A)^2
 *
 * - \text{tr} (\mathbf{A}^2) \right]$  。
 * 对于这个函数的参数种类，即大小为1的秩2张量，其结果只是零。
 * @relatesalso  SymmetricTensor
 *
 *
 */
template <typename Number>
constexpr DEAL_II_ALWAYS_INLINE Number
                                second_invariant(const SymmetricTensor<2, 1, Number> &)
{
  return internal::NumberType<Number>::value(0.0);
}



/**
 * 计算一个等级为2的张量的第二不变量。张量 $\mathbf A$
 * 的第二不变量定义为 $I_2 (\mathbf A) = II(\mathbf A) = \frac 12
 * \left[ (\text{tr} \mathbf A)^2
 *
 * - \text{tr} (\mathbf{A}^2) \right]$  。
 * 对于这个函数的参数种类，即大小为2的对称秩-2张量，其结果是（从1开始计算指数）
 * $I_2(\mathbf A) = II(\mathbf A) = \frac 12 \left[ (A_{11} + A_{22})^2
 *
 * - (A_{11}^2+2 A_{12}^2+ A_{22}^2) \right] = A_{11} A_{22}
 *
 * - A_{12}^2$  。正如预期的那样，对于这个函数处理的
 * $2\times 2$ 对称张量，这等于张量的行列式。这是因为对于
 * $2\times 2$
 * 对称张量，实际上只有两个不变式，所以第二个和第三个不变式是一样的；行列式是第三个不变式）。
 * @relatesalso SymmetricTensor
 *
 *
 */
template <typename Number>
constexpr DEAL_II_ALWAYS_INLINE Number
                                second_invariant(const SymmetricTensor<2, 2, Number> &t)
{
  return t[0][0] * t[1][1] - t[0][1] * t[0][1];
}



/**
 * 计算一个等级为2的张量的第二个不变式。张量 $\mathbf A$
 * 的第二不变量定义为 $I_2 (\mathbf A) = II(\mathbf A) = \frac 12
 * \left[ (\text{tr} \mathbf A)^2
 *
 * - \text{tr} (\mathbf{A}^2) \right]$  。
 * @relatesalso SymmetricTensor
 *
 */
template <typename Number>
constexpr DEAL_II_ALWAYS_INLINE Number
                                second_invariant(const SymmetricTensor<2, 3, Number> &t)
{
  return (t[0][0] * t[1][1] + t[1][1] * t[2][2] + t[2][2] * t[0][0] -
          t[0][1] * t[0][1] - t[0][2] * t[0][2] - t[1][2] * t[1][2]);
}



/**
 * 返回一个对称的 $1 \times 1$
 * 张量的特征值。当然，张量的（单一）条目等于（单一）特征值。
 * @relatesalso  SymmetricTensor
 *
 *
 */
template <typename Number>
std::array<Number, 1>
eigenvalues(const SymmetricTensor<2, 1, Number> &T);



/**
 * 返回一个对称的 $2\times 2$
 * 张量的特征值。特征值的数组按降序排序。 对于 $2\times
 * 2$ 张量，张量 $\mathbf T$ 的特征值是<a
 * href="https://en.wikipedia.org/wiki/Eigenvalue_algorithm#2.C3.972_matrices">the
 * characteristic polynomial</a>  $0 = \lambda^2
 *
 *
 *
 * - \lambda\;\text{tr}\mathbf{T} + \det \mathbf{T}$ 的根，由 $\lambda_1,
 * \lambda_2 = \frac{1}{2} \left[ \text{tr} \mathbf{T} \pm \sqrt{(\text{tr}
 * \mathbf{T})^2
 *
 * - 4 \det \mathbf{T}} \right]$ 给出。
 * @warning
 * 这里采用的算法通过计算特征多项式的根来确定特征值。在存在共根的情况下（特征值相等），计算结果为<a
 * href="https://scicomp.stackexchange.com/q/23686">subject to round-off
 * errors</a>阶的 $\sqrt{\epsilon}$
 * 。作为一种替代方法，eigenvectors()函数提供了一种更稳健但昂贵的方法来计算对称张量的特征值。
 * @relatesalso  SymmetricTensor
 *
 *
 */
template <typename Number>
std::array<Number, 2>
eigenvalues(const SymmetricTensor<2, 2, Number> &T);



/**
 * 返回一个对称 $3\times 3$
 * 张量的特征值。特征值的数组按降序排序。 对于 $3\times
 * 3$ 张量，张量 $\mathbf T$ 的特征值是<a
 * href="https://en.wikipedia.org/wiki/Eigenvalue_algorithm#3.C3.973_matrices">the
 * characteristic polynomial</a>  $0 = \lambda^3
 *
 * - \lambda^2\;\text{tr}\mathbf T
 *
 *
 *
 * - \frac{1}{2} \lambda \left[\text{tr}(\mathbf{T}^2)
 *
 * - (\text{tr}\mathbf T)^2\right]
 *
 * - \det \mathbf T$  的根。
 * @warning
 * 这里采用的算法通过计算特征多项式的根来确定特征值。在存在共根的情况下（特征值相等），计算是<a
 * href="https://scicomp.stackexchange.com/q/23686">subject to round-off
 * errors</a>的顺序 $\sqrt{\epsilon}$
 * 。作为一种替代方法，eigenvectors()函数提供了一种更稳健但昂贵的方法来计算一个对称张量的特征值。
 * @relatesalso  SymmetricTensor
 *
 *
 */
template <typename Number>
std::array<Number, 3>
eigenvalues(const SymmetricTensor<2, 3, Number> &T);



namespace internal
{
  namespace SymmetricTensorImplementation
  {
    /**
     * 使用Householder方法对一个等级2的对称张量进行三对角化。
     * 这里实现的专门算法是在
     * @code{.bib}
     * @article{Kopp2008,
     * title       = {Efficient numerical diagonalization of hermitian 3x3
     *                matrices},
     * author      = {Kopp, J.},
     * journal     = {International Journal of Modern Physics C},
     * year        = {2008},
     * volume      = {19},
     * number      = {3},
     * pages       = {523--548},
     * doi         = {10.1142/S0129183108012303},
     * eprinttype  = {arXiv},
     * eprint      = {physics/0610206v3},
     * eprintclass = {physics.comp-ph},
     * url         =
     * {https://www.mpi-hd.mpg.de/personalhomes/globes/3x3/index.html}
     * }
     * @endcode
     * 中给出的，并且是基于11.3.2节中提出的通用算法。
     * @code{.bib}
     * @book{Press2007,
     * title   = {Numerical recipes 3rd edition: The art of scientific
     *            computing},
     * author  = {Press, W. H.},
     * journal = {Cambridge university press},
     * year    = {2007}
     * }
     * @endcode
     * @param[in]  A 要被三对角化的张量  @param[out]  Q
     * 实现转换的正交矩阵  @param[out]  d
     * 三对角矩阵的对角线元素  @param[out]  e
     * 三对角矩阵的非对角线元素
     *
     */
    template <int dim, typename Number>
    void
    tridiagonalize(const dealii::SymmetricTensor<2, dim, Number> &A,
                   dealii::Tensor<2, dim, Number> &               Q,
                   std::array<Number, dim> &                      d,
                   std::array<Number, dim - 1> &                  e);



    /**
     * 使用QL算法计算实值等级为2的对称张量的特征值和特征向量，并采用隐式移位的方法。
     * 这里实现的专门算法是在
     * @code{.bib}
     * @article{Kopp2008,
     * title       = {Efficient numerical diagonalization of hermitian 3x3
     *                matrices},
     * author      = {Kopp, J.},
     * journal     = {International Journal of Modern Physics C},
     * year        = {2008},
     * volume      = {19},
     * number      = {3},
     * pages       = {523--548},
     * doi         = {10.1142/S0129183108012303},
     * eprinttype  = {arXiv},
     * eprint      = {physics/0610206v3},
     * eprintclass = {physics.comp-ph},
     * url         =
     * {https://www.mpi-hd.mpg.de/personalhomes/globes/3x3/index.html}
     * }
     * @endcode
     * 的第11.4.3节中提出的通用算法的基础上。
     * @code{.bib}
     * @book{Press2007,
     * title   = {Numerical recipes 3rd edition: The art of scientific
     *            computing},
     * author  = {Press, W. H.},
     * journal = {Cambridge university press},
     * year    = {2007}
     * }
     * @endcode
     * @param[in]  A 要计算特征向量和特征值的张量。
     * @return  一个包含特征向量和相关特征值的数组。
     * 该数组不按任何特定顺序排序。
     *
     */
    template <int dim, typename Number>
    std::array<std::pair<Number, Tensor<1, dim, Number>>, dim>
    ql_implicit_shifts(const dealii::SymmetricTensor<2, dim, Number> &A);



    /**
     * 使用雅可比算法计算实值秩2对称张量的特征值和特征向量。
     * 这里实现的专门算法是在
     * @code{.bib}
     * @article{Kopp2008,
     * title       = {Efficient numerical diagonalization of hermitian 3x3
     *                matrices},
     * author      = {Kopp, J.},
     * journal     = {International Journal of Modern Physics C},
     * year        = {2008},
     * volume      = {19},
     * number      = {3},
     * pages       = {523--548},
     * doi         = {10.1142/S0129183108012303},
     * eprinttype  = {arXiv},
     * eprint      = {physics/0610206v3},
     * eprintclass = {physics.comp-ph},
     * url         =
     * {https://www.mpi-hd.mpg.de/personalhomes/globes/3x3/index.html}
     * }
     * @endcode
     * 的第11.4.3节中提出的通用算法的基础上。
     * @code{.bib}
     * @book{Press2007,
     * title   = {Numerical recipes 3rd edition: The art of scientific
     *            computing},
     * author  = {Press, W. H.},
     * journal = {Cambridge university press},
     * year    = {2007}
     * }
     * @endcode
     * @param[in]  A 要计算的特征向量和特征值的张量。
     * @return  一个包含特征向量和相关特征值的数组。
     * 该数组不按任何特定顺序排序。
     *
     */
    template <int dim, typename Number>
    std::array<std::pair<Number, Tensor<1, dim, Number>>, dim>
      jacobi(dealii::SymmetricTensor<2, dim, Number> A);



    /**
     * 计算实值等级为2的对称2x2张量的特征值和特征向量，使用特征方程计算特征值，并使用基于交叉积的分析方法计算特征向量。如果计算结果被认为太不准确，那么该方法就会退回到ql_implicit_shifts。
     * @param[in]  A 将要计算特征向量和特征值的张量。
     * @return  一个包含特征向量和相关特征值的数组。
     * 该数组不按任何特定顺序排序。
     *
     */
    template <typename Number>
    std::array<std::pair<Number, Tensor<1, 2, Number>>, 2>
    hybrid(const dealii::SymmetricTensor<2, 2, Number> &A);



    /**
     * 计算实值等级为2的对称3x3张量的特征值和特征向量，使用特征方程计算特征值，并使用基于交叉积的分析方法计算特征向量。如果计算结果被认为太不准确，那么该方法就会退回到ql_implicit_shifts。这里实现的专门算法是在
     * @code{.bib}
     * @article{Kopp2008,
     * title       = {Efficient numerical diagonalization of hermitian 3x3
     *                matrices},
     * author      = {Kopp, J.},
     * journal     = {International Journal of Modern Physics C},
     * year        = {2008},
     * volume      = {19},
     * number      = {3},
     * pages       = {523--548},
     * doi         = {10.1142/S0129183108012303},
     * eprinttype  = {arXiv},
     * eprint      = {physics/0610206v3},
     * eprintclass = {physics.comp-ph},
     * url         =
     * {https://www.mpi-hd.mpg.de/personalhomes/globes/3x3/index.html}
     * }
     * @endcode
     * @param[in]  A 要计算特征向量和特征值的张量。
     * @return  一个包含特征向量和相关特征值的数组。
     * 该数组不按任何特定顺序排序。
     *
     */
    template <typename Number>
    std::array<std::pair<Number, Tensor<1, 3, Number>>, 3>
    hybrid(const dealii::SymmetricTensor<2, 3, Number> &A);

    /**
     * 一个结构，用于对eign=envalues和特征向量的数组进行排序。排序是按照特征值的降序进行的。
     *
     */
    template <int dim, typename Number>
    struct SortEigenValuesVectors
    {
      using EigValsVecs = std::pair<Number, Tensor<1, dim, Number>>;
      bool
      operator()(const EigValsVecs &lhs, const EigValsVecs &rhs)
      {
        return lhs.first > rhs.first;
      }
    };

  } // namespace SymmetricTensorImplementation

} // namespace internal



// The line below is to ensure that doxygen puts the full description
// of this global enumeration into the documentation
// See https://stackoverflow.com/a/1717984
 /** @file */ 
/**
 * 枚举算法，用于在对SymmetricTensor对象使用eigenvalues()和eigenvectors()方法进行归一化特征向量及其相应特征值的计算时使用。
 * 在计算特征向量时使用的专门算法在下文中介绍。
 *
 * @code{.bib}
 * @article{Kopp2008,
 * title       = {Efficient numerical diagonalization of hermitian 3x3
 *                matrices},
 * author      = {Kopp, J.},
 * journal     = {International Journal of Modern Physics C},
 * year        = {2008},
 * volume      = {19},
 * number      = {3},
 * pages       = {523--548},
 * doi         = {10.1142/S0129183108012303},
 * eprinttype  = {arXiv},
 * eprint      = {physics/0610206v3},
 * eprintclass = {physics.comp-ph},
 * url         =
 * {https://www.mpi-hd.mpg.de/personalhomes/globes/3x3/index.html}
 * }
 * @endcode
 *
 *
 *
 */
enum struct SymmetricTensorEigenvectorMethod
{
  /**
   * 一种混合方法，优先使用特征方程来计算特征值，并使用基于交叉积的分析方法来计算特征向量。如果计算结果被认为太不准确，那么该方法就会退回到ql_implicit_shifts。
   * 如果没有遇到病态的情况，这种方法可能会提供最快速的计算。
   *
   */
  hybrid,
  /**
   * 迭代QL算法，在使用Householder方法对张量进行三对角化后应用隐式移位。
   * 这种方法在计算速度和它的稳健性之间提供了一个折中。当
   * $T$
   * 的元素有很大的变化幅度时，这种方法特别有用，这通常会导致在计算较小的特征值时失去准确性。
   *
   */
  ql_implicit_shifts,
  /**
   * 迭代雅可比算法。
   * 这种方法提供的是现有选项中最稳健的，即使是最病态的情况也能得到可靠的结果。然而，它是所有实施的算法中最慢的算法。
   *
   */
  jacobi
};



/**
 * 返回实值等级2的对称张量的特征值和特征向量  $\mathbf T$
 * 。匹配的特征值和特征向量对阵列按降序排序（由特征值决定）。
 * 在计算特征向量时利用的专门算法见于
 *
 * @code{.bib}
 * @article{Kopp2008,
 * title       = {Efficient numerical diagonalization of hermitian 3x3
 *                matrices},
 * author      = {Kopp, J.},
 * journal     = {International Journal of Modern Physics C},
 * year        = {2008},
 * volume      = {19},
 * number      = {3},
 * pages       = {523--548},
 * doi         = {10.1142/S0129183108012303},
 * eprinttype  = {arXiv},
 * eprint      = {physics/0610206v3},
 * eprintclass = {physics.comp-ph},
 * url         =
 * {https://www.mpi-hd.mpg.de/personalhomes/globes/3x3/index.html}
 * }
 * @endcode
 *
 * @relatesalso  SymmetricTensor
 *
 *
 */
template <int dim, typename Number>
std::array<std::pair<Number, Tensor<1, dim, Number>>,
           std::integral_constant<int, dim>::value>
eigenvectors(const SymmetricTensor<2, dim, Number> &T,
             const SymmetricTensorEigenvectorMethod method =
               SymmetricTensorEigenvectorMethod::ql_implicit_shifts);



/**
 * 返回给定的对称张量的转置。由于我们正在处理对称对象，转置当然与原始张量相同。这个函数的存在主要是为了与张量类兼容。
 * @relatesalso  SymmetricTensor
 *
 *
 */
template <int rank_, int dim, typename Number>
constexpr DEAL_II_ALWAYS_INLINE SymmetricTensor<rank_, dim, Number>
                                transpose(const SymmetricTensor<rank_, dim, Number> &t)
{
  return t;
}



/**
 * 计算对称张量的偏差，其定义为  $\text{dev} \mathbf T = \mathbf
 * T
 *
 * - \frac{1}{\text{dim}} \text{tr}\mathbf T \; \mathbf I$  ，其中
 * $\mathbf I$
 * 是身份算子。这个量等于原始张量减去其收缩或扩张成分，指的是弹性等方面的剪切力。
 * @relatesalso SymmetricTensor
 *
 *
 */
template <int dim, typename Number>
constexpr inline DEAL_II_ALWAYS_INLINE SymmetricTensor<2, dim, Number>
                                       deviator(const SymmetricTensor<2, dim, Number> &t)
{
  SymmetricTensor<2, dim, Number> tmp = t;

  // subtract scaled trace from the diagonal
  const Number tr = trace(t) / dim;
  for (unsigned int i = 0; i < dim; ++i)
    tmp.data[i] -= tr;

  return tmp;
}



/**
 * 返回一个等级为2的单位对称张量，即
 * $\text{dim}\times\text{dim}$ 身份矩阵 $\mathbf I$  。
 * @relatesalso  SymmetricTensor
 *
 *
 */
template <int dim, typename Number>
DEAL_II_CONSTEXPR inline DEAL_II_ALWAYS_INLINE SymmetricTensor<2, dim, Number>
                                               unit_symmetric_tensor()
{
  // create a default constructed matrix filled with
  // zeros, then set the diagonal elements to one
  SymmetricTensor<2, dim, Number> tmp;
  switch (dim)
    {
      case 1:
        tmp.data[0] = internal::NumberType<Number>::value(1.);
        break;
      case 2:
        tmp.data[0] = tmp.data[1] = internal::NumberType<Number>::value(1.);
        break;
      case 3:
        tmp.data[0] = tmp.data[1] = tmp.data[2] =
          internal::NumberType<Number>::value(1.);
        break;
      default:
        for (unsigned int d = 0; d < dim; ++d)
          tmp.data[d] = internal::NumberType<Number>::value(1.);
    }
  return tmp;
}



/**
 * unit_symmetric_tensor<dim>()是unit_symmetric_tensor<dim,Number>()函数的特殊化，它使用
 * <code>double</code>  作为元素的数据类型。
 * @relatesalso  SymmetricTensor
 *
 *
 */
template <int dim>
DEAL_II_CONSTEXPR inline DEAL_II_ALWAYS_INLINE SymmetricTensor<2, dim>
                                               unit_symmetric_tensor()
{
  return unit_symmetric_tensor<dim, double>();
}



/**
 * 返回等级4的张量，当与对称等级2的张量 $\mathbf T$
 * 相乘时，返回偏差 $\text{dev}\ \mathbf T$
 * 。它是线性偏差算子 $\mathbb P$
 * 的算子表示，也被称为体积投影张量，计算公式为：。\f{align*}{
 * \mathbb{P} &=\mathbb{I}
 *
 * -\frac{1}{\text{dim}} \mathbf I \otimes \mathbf I \\ \mathcal{P}_{ijkl} &=
 * \frac 12 \left(\delta_{ik} \delta_{jl} + \delta_{il} \delta_{jk} \right)
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * - \frac{1}{\text{dim}} \delta_{ij} \delta_{kl} \f}
 * 对于每一个张量<tt>T</tt>，都有一个特征<tt>deviator<dim,Number>(T)
 * == deviator_tensor<dim,Number>() T</tt>，直到数字上的舍入。\f[
 * \text{dev}\mathbf T = \mathbb P : \mathbf T \f]
 *
 *
 * @note
 * 提供这种运算符表示的原因是为了简化对张量的偏离部分的求导。\f[
 * \frac{\partial \text{dev}\mathbf{T}}{\partial \mathbf T} = \mathbb P. \f]
 * @relatesalso  SymmetricTensor
 *
 *
 */
template <int dim, typename Number>
DEAL_II_CONSTEXPR inline SymmetricTensor<4, dim, Number>
deviator_tensor()
{
  SymmetricTensor<4, dim, Number> tmp;

  // fill the elements treating the diagonal
  for (unsigned int i = 0; i < dim; ++i)
    for (unsigned int j = 0; j < dim; ++j)
      tmp.data[i][j] =
        internal::NumberType<Number>::value((i == j ? 1. : 0.) - 1. / dim);

  // then fill the ones that copy over the
  // non-diagonal elements. note that during
  // the double-contraction, we handle the
  // off-diagonal elements twice, so simply
  // copying requires a weight of 1/2
  for (unsigned int i = dim;
       i < internal::SymmetricTensorAccessors::StorageType<4, dim, Number>::
             n_rank2_components;
       ++i)
    tmp.data[i][i] = internal::NumberType<Number>::value(0.5);

  return tmp;
}



/**
 * 这个版本的deviator_tensor<dim>()函数是deviator_tensor<dim,Number>()的特殊化，使用<tt>double</tt>作为张量元素的数据类型。
 * @relatesalso  SymmetricTensor
 *
 *
 */
template <int dim>
DEAL_II_CONSTEXPR inline DEAL_II_ALWAYS_INLINE SymmetricTensor<4, dim>
                                               deviator_tensor()
{
  return deviator_tensor<dim, double>();
}



/**
 * 返回四阶对称身份张量 $\mathbb I$
 * ，它将对称的二阶张量，如 $\mathbf A$
 * ，映射到它们自己。\f[ \mathbb I : \mathbf A = \mathbf A \f]
 * 请注意，这个张量，尽管它是同一性的，但它的形式有些滑稽，特别是不只由零和一组成。例如，对于<tt>dim=2</tt>，身份张量除了\f[
 * \mathcal{I}_{0000} = \mathcal{I}_{1111} = 1 \f] \f[ \mathcal{I}_{0101} =
 * \mathcal{I}_{0110} = \mathcal{I}_{1001} = \mathcal{I}_{1010} = \frac 12.
 * \f] 以指数符号表示，我们可以写出一般的形式\f[
 * \mathcal{I}_{ijkl} = \frac 12 \left( \delta_{ik} \delta_{jl} + \delta_{il}
 * \delta_{jl} \right). \f] 为了了解为什么 $1 / 2$
 * 的这个因子是必要的，考虑计算 $\mathbf A= \mathbb I : \mathbf
 * B$  。对于元素  $A_{01}$  我们有  $A_{01} = \mathcal{I}_{0100}
 * B_{00} + \mathcal{I}_{0111} B_{11} + \mathcal{I}_{0101} B_{01} +
 * \mathcal{I}_{0110} B_{10}$  。另一方面，我们需要有 $A_{01} =
 * B_{01}$ ，对称性意味着 $B_{01}=B_{10}$ ，导致 $A_{01} =
 * (\mathcal{I}_{0101} + \mathcal{I}_{0110}) B_{01}$
 * ，或者，同样通过对称性， $\mathcal{I}_{0101} =
 * \mathcal{I}_{0110} = \frac 12$
 * 。类似的考虑也适用于三维的情况。 这个问题在 step-44
 * 的介绍中也有解释。
 * @relatesalso  SymmetricTensor
 *
 *
 */
template <int dim, typename Number>
DEAL_II_CONSTEXPR inline DEAL_II_ALWAYS_INLINE SymmetricTensor<4, dim, Number>
                                               identity_tensor()
{
  SymmetricTensor<4, dim, Number> tmp;

  // fill the elements treating the diagonal
  for (unsigned int i = 0; i < dim; ++i)
    tmp.data[i][i] = internal::NumberType<Number>::value(1.);

  // then fill the ones that copy over the
  // non-diagonal elements. note that during
  // the double-contraction, we handle the
  // off-diagonal elements twice, so simply
  // copying requires a weight of 1/2
  for (unsigned int i = dim;
       i < internal::SymmetricTensorAccessors::StorageType<4, dim, Number>::
             n_rank2_components;
       ++i)
    tmp.data[i][i] = internal::NumberType<Number>::value(0.5);

  return tmp;
}



/**
 * 这个版本的identity_tensor<dim>()函数是identity_tensor<dim,Number>()的特殊化，它使用<tt>double</tt>作为张量的元素的数据类型。
 * @relatesalso  SymmetricTensor
 *
 *
 */
template <int dim>
DEAL_II_CONSTEXPR inline DEAL_II_ALWAYS_INLINE SymmetricTensor<4, dim>
                                               identity_tensor()
{
  return identity_tensor<dim, double>();
}



/**
 * 反转一个对称的秩2张量。
 *
 *
 * @note
 * 如果一个张量不是可反转的，那么其结果是不明确的，但很可能包含除以0的结果或至少是一个非常小的数字。
 * @relatesalso  SymmetricTensor
 *
 *
 */
template <int dim, typename Number>
constexpr DEAL_II_ALWAYS_INLINE SymmetricTensor<2, dim, Number>
                                invert(const SymmetricTensor<2, dim, Number> &t)
{
  return internal::SymmetricTensorImplementation::Inverse<2, dim, Number>::
    value(t);
}



/**
 * 反转一个对称的秩-4张量。由于对称秩-4张量是对称秩-2张量的映射，它们可以有一个反转。
 * 如果一个张量不是可逆的，那么其结果是不明确的，但可能包含除以0的结果或至少是一个非常小的数字。
 * @relatesalso  SymmetricTensor
 *
 *
 */
template <int dim, typename Number>
constexpr SymmetricTensor<4, dim, Number>
invert(const SymmetricTensor<4, dim, Number> &t)
{
  return internal::SymmetricTensorImplementation::Inverse<4, dim, Number>::
    value(t);
}



/**
 * 返回等级为4的张量，它是作为参数给出的两个张量的外积，即结果
 * $\mathbb A = \mathbf{T}_1 \otimes \mathbf{T}_2$ 满足 $\mathbb A : \mathbf
 * B = (\mathbf{T}_2 : \mathbf B) \mathbf{T}_1$ 的所有对称张量 $\mathbf
 * B$  。在索引符号中\f[ \mathcal{A}_{ijkl} = (T_1)_{ij} (T_2)_{kl} \f]
 * 。 例如，偏差张量 $\mathbb P = \mathbb I
 *
 * - \frac{1}{\text{dim}} \mathbf I \otimes \mathbf I$
 * 可以被计算为<tt>identity_tensor<dim>()
 *
 * - 1/d outer_product (unit_symmetric_tensor<dim>(), unit_symmetric_tensor<dim>())</tt>，因为与单位张量的（双）收缩会产生对称张量的迹线（ $\mathbf I : \mathbf B = \text{tr} \mathbf B$  ）。
 * @relatesalso  SymmetricTensor
 *
 *
 */
template <int dim, typename Number>
constexpr inline SymmetricTensor<4, dim, Number>
outer_product(const SymmetricTensor<2, dim, Number> &t1,
              const SymmetricTensor<2, dim, Number> &t2)
{
  SymmetricTensor<4, dim, Number> tmp;

  // fill only the elements really needed
  for (unsigned int i = 0; i < dim; ++i)
    for (unsigned int j = i; j < dim; ++j)
      for (unsigned int k = 0; k < dim; ++k)
        for (unsigned int l = k; l < dim; ++l)
          tmp[i][j][k][l] = t1[i][j] * t2[k][l];

  return tmp;
}



/**
 * 返回全秩2张量的对称版本，即 $\text{sym}\mathbf A = \frac 12
 * \left(\mathbf A + \mathbf{A}^T\right)$
 * ，作为对称秩2张量。这是对一般维度的版本。
 * @relatesalso  SymmetricTensor
 *
 *
 */
template <int dim, typename Number>
constexpr inline DEAL_II_ALWAYS_INLINE SymmetricTensor<2, dim, Number>
                                       symmetrize(const Tensor<2, dim, Number> &t)
{
  SymmetricTensor<2, dim, Number> result;
  for (unsigned int d = 0; d < dim; ++d)
    result[d][d] = t[d][d];

  const Number half = internal::NumberType<Number>::value(0.5);
  for (unsigned int d = 0; d < dim; ++d)
    for (unsigned int e = d + 1; e < dim; ++e)
      result[d][e] = (t[d][e] + t[e][d]) * half;
  return result;
}



/**
 * 一般等级的对称张量与一个来自右边的标量相乘。如果标量的数据类型与用于存储对称张量元素的数据类型相同，则使用该运算符的这个版本。
 * @relatesalso  SymmetricTensor
 *
 *
 */
template <int rank_, int dim, typename Number>
constexpr inline DEAL_II_ALWAYS_INLINE SymmetricTensor<rank_, dim, Number>
                                       operator*(const SymmetricTensor<rank_, dim, Number> &t, const Number &factor)
{
  SymmetricTensor<rank_, dim, Number> tt = t;
  tt *= factor;
  return tt;
}



/**
 * 一般等级的对称张量与来自左边的标量相乘。如果标量的数据类型与用于存储对称张量元素的数据类型相同，则使用该运算符的这个版本。
 * @relatesalso  SymmetricTensor
 *
 *
 */
template <int rank_, int dim, typename Number>
constexpr DEAL_II_ALWAYS_INLINE SymmetricTensor<rank_, dim, Number>
                                operator*(const Number &factor, const SymmetricTensor<rank_, dim, Number> &t)
{
  // simply forward to the other operator
  return t * factor;
}



/**
 * 对称张量与来自右边的标量数字的乘法。
 * 这个运算符的目的是只实现张量与标量数（即浮点数、复数浮点数等）的乘法。该函数的写法是，只有在第二个参数确实是标量数的情况下，编译器才会考虑这个函数
 *
 * - 换句话说， @p OtherNumber 将不匹配，例如 <code>std::vector@<double@></code> ，因为张量和矢量的乘积显然没有意义。在EnableIfScalar类的文档中解释了编译器禁止考虑这个操作符与非标量类型的乘法的机制。
 * 函数的返回类型被选择为与张量和标量参数的类型一致。例如，如果你用
 * <code>SymmetricTensor@<2,dim,double@></code>  乘以
 * <code>std::complex@<double@></code>  ，那么结果将是
 * <code>SymmetricTensor@<2,dim,std::complex@<double@>@></code>
 * 。换句话说，返回的张量存储其组件的类型等于你用标量因子乘以输入张量的单个组件后得到的类型。
 * @relatesalso  SymmetricTensor  @relatesalso  EnableIfScalar
 *
 *
 */
template <int rank_, int dim, typename Number, typename OtherNumber>
constexpr inline DEAL_II_ALWAYS_INLINE SymmetricTensor<
  rank_,
  dim,
  typename ProductType<Number,
                       typename EnableIfScalar<OtherNumber>::type>::type>
operator*(const SymmetricTensor<rank_, dim, Number> &t,
          const OtherNumber &                        factor)
{
  // form the product. we have to convert the two factors into the final
  // type via explicit casts because, for awkward reasons, the C++
  // standard committee saw it fit to not define an
  //   operator*(float,std::complex<double>)
  // (as well as with switched arguments and double<->float).
  using product_type = typename ProductType<Number, OtherNumber>::type;
  SymmetricTensor<rank_, dim, product_type> tt(t);
  tt *= internal::NumberType<product_type>::value(factor);
  return tt;
}



/**
 * 对称张量与左边的标量数字相乘。关于模板参数和返回类型的更多信息，请参见与带交换参数的运算符的讨论。
 * @relatesalso  SymmetricTensor  @relatesalso  EnableIfScalar
 *
 *
 */
template <int rank_, int dim, typename Number, typename OtherNumber>
constexpr inline DEAL_II_ALWAYS_INLINE SymmetricTensor<
  rank_,
  dim,
  typename ProductType<OtherNumber,
                       typename EnableIfScalar<Number>::type>::type>
operator*(const Number &                                  factor,
          const SymmetricTensor<rank_, dim, OtherNumber> &t)
{
  // simply forward to the other operator with switched arguments
  return (t * factor);
}



/**
 * 一般等级的对称张量除以一个标量。
 * @relatesalso  SymmetricTensor
 *
 *
 */
template <int rank_, int dim, typename Number, typename OtherNumber>
constexpr inline SymmetricTensor<
  rank_,
  dim,
  typename ProductType<Number,
                       typename EnableIfScalar<OtherNumber>::type>::type>
operator/(const SymmetricTensor<rank_, dim, Number> &t,
          const OtherNumber &                        factor)
{
  using product_type = typename ProductType<Number, OtherNumber>::type;
  SymmetricTensor<rank_, dim, product_type> tt(t);
  tt /= internal::NumberType<product_type>::value(factor);
  return tt;
}



/**
 * 一般等级的对称张量与一个标量从右边相乘。
 * @relatesalso  SymmetricTensor
 *
 *
 */
template <int rank_, int dim>
constexpr inline DEAL_II_ALWAYS_INLINE SymmetricTensor<rank_, dim>
                                       operator*(const SymmetricTensor<rank_, dim> &t, const double factor)
{
  SymmetricTensor<rank_, dim> tt(t);
  tt *= factor;
  return tt;
}



/**
 * 一般等级的对称张量与来自左边的标量相乘。
 * @relatesalso  SymmetricTensor
 *
 *
 */
template <int rank_, int dim>
constexpr inline DEAL_II_ALWAYS_INLINE SymmetricTensor<rank_, dim>
                                       operator*(const double factor, const SymmetricTensor<rank_, dim> &t)
{
  SymmetricTensor<rank_, dim> tt(t);
  tt *= factor;
  return tt;
}



/**
 * 一般等级的对称张量除以一个标量。
 * @relatesalso  SymmetricTensor
 *
 *
 */
template <int rank_, int dim>
constexpr inline SymmetricTensor<rank_, dim>
operator/(const SymmetricTensor<rank_, dim> &t, const double factor)
{
  SymmetricTensor<rank_, dim> tt(t);
  tt /= factor;
  return tt;
}

/**
 * 计算两个秩为2的张量 $\mathbf A, \mathbf B$ 之间的标量乘积
 * $\mathbf A: \mathbf B=\sum_{i,j} A_{ij}B_{ij}$
 * 。在当前两个参数都是对称张量的情况下，这相当于调用表达式
 * <code>A*B</code> ，它使用 <code>SymmetricTensor::operator*()</code>
 * 。
 * @relatesalso  SymmetricTensor
 *
 *
 */
template <int dim, typename Number, typename OtherNumber>
constexpr DEAL_II_ALWAYS_INLINE typename ProductType<Number, OtherNumber>::type
scalar_product(const SymmetricTensor<2, dim, Number> &     t1,
               const SymmetricTensor<2, dim, OtherNumber> &t2)
{
  return (t1 * t2);
}


/**
 * 计算两个等级为2的张量 $\mathbf A, \mathbf B$
 * 之间的标量乘积 $\mathbf A: \mathbf B=\sum_{i,j} A_{ij}B_{ij}$
 * 。我们不使用 <code>operator*</code>
 * 进行这个操作，因为两个张量之间的乘积通常被认为是对第一个张量的最后一个索引和第二个张量的第一个索引的收缩。例如，如果<tt>B</tt>是一个张量，调用<tt>A*B</tt>（而不是<tt>scalar_product(A,B)</tt>）提供
 * $(\mathbf A \cdot\mathbf B)_{ij}=\sum_k A_{ik}B_{kj}$  。
 * @relatesalso  Tensor  @relatesalso  SymmetricTensor
 *
 *
 */
template <int dim, typename Number, typename OtherNumber>
constexpr inline DEAL_II_ALWAYS_INLINE
  typename ProductType<Number, OtherNumber>::type
  scalar_product(const SymmetricTensor<2, dim, Number> &t1,
                 const Tensor<2, dim, OtherNumber> &    t2)
{
  typename ProductType<Number, OtherNumber>::type s = internal::NumberType<
    typename ProductType<Number, OtherNumber>::type>::value(0.0);
  for (unsigned int i = 0; i < dim; ++i)
    for (unsigned int j = 0; j < dim; ++j)
      s += t1[i][j] * t2[i][j];
  return s;
}


/**
 * 计算两个等级为2的张量 $\mathbf A, \mathbf B$
 * 之间的标量乘积 $\mathbf A:\mathbf B=\sum_{i,j} A_{ij}B_{ij}$
 * 。我们不使用 <code>operator*</code>
 * 进行这个操作，因为两个张量之间的乘积通常被认为是对第一个张量的最后一个索引和第二个张量的第一个索引的收缩。例如，如果<tt>A</tt>是一个张量，调用<tt>A*B</tt>（而不是<tt>scalar_product(A,B)</tt>）提供
 * $(\mathbf A \cdot\mathbf B)_{ij}=\sum_k A_{ik}B_{kj}$  。
 * @relatesalso  Tensor  @relatesalso  SymmetricTensor
 *
 *
 */
template <int dim, typename Number, typename OtherNumber>
constexpr DEAL_II_ALWAYS_INLINE typename ProductType<Number, OtherNumber>::type
scalar_product(const Tensor<2, dim, Number> &              t1,
               const SymmetricTensor<2, dim, OtherNumber> &t2)
{
  return scalar_product(t2, t1);
}


/**
 * 在等级4和等级2的对称张量之间进行双重收缩，产生等级2的对称张量，作为该函数的第一个参数给出。这个操作是矩阵-向量乘法的对称张量类似物。
 * 这个函数的作用与 SymmetricTensor::operator*().
 * 相同，但是不应该使用它，因为成员运算器知道实际的数据存储格式，并且至少快两个数量级。这个函数的存在主要是为了与一般的张量类兼容。
 * @relatesalso  SymmetricTensor
 *
 *
 */
template <typename Number, typename OtherNumber>
constexpr inline DEAL_II_ALWAYS_INLINE void double_contract(
  SymmetricTensor<2, 1, typename ProductType<Number, OtherNumber>::type> &tmp,
  const SymmetricTensor<4, 1, Number> &                                   t,
  const SymmetricTensor<2, 1, OtherNumber> &                              s)
{
  tmp[0][0] = t[0][0][0][0] * s[0][0];
}



/**
 * 在等级4和等级2的对称张量之间进行双重收缩，产生等级2的对称张量，作为该函数的第一个参数给出。这个操作是矩阵-向量乘法的对称张量类似物。
 * 这个函数的作用与 SymmetricTensor::operator*().
 * 相同，但是不应该使用它，因为成员运算器知道实际的数据存储格式，并且至少快两个数量级。这个函数的存在主要是为了与一般的张量类兼容。
 * @relatesalso  SymmetricTensor
 *
 *
 */
template <typename Number, typename OtherNumber>
constexpr inline void double_contract(
  SymmetricTensor<2, 1, typename ProductType<Number, OtherNumber>::type> &tmp,
  const SymmetricTensor<2, 1, Number> &                                   s,
  const SymmetricTensor<4, 1, OtherNumber> &                              t)
{
  tmp[0][0] = t[0][0][0][0] * s[0][0];
}



/**
 * 在等级4和等级2的对称张量之间进行双重收缩，产生等级2的对称张量，作为这个函数的第一个参数给出。这个操作是矩阵-向量乘法的对称张量类似物。
 * 这个函数的作用与 SymmetricTensor::operator*().
 * 相同，但是不应该使用它，因为成员运算器知道实际的数据存储格式，并且至少快两个数量级。这个函数的存在主要是为了与一般的张量类兼容。
 * @relatesalso  SymmetricTensor
 *
 *
 */
template <typename Number, typename OtherNumber>
constexpr inline void double_contract(
  SymmetricTensor<2, 2, typename ProductType<Number, OtherNumber>::type> &tmp,
  const SymmetricTensor<4, 2, Number> &                                   t,
  const SymmetricTensor<2, 2, OtherNumber> &                              s)
{
  const unsigned int dim = 2;

  for (unsigned int i = 0; i < dim; ++i)
    for (unsigned int j = i; j < dim; ++j)
      tmp[i][j] = t[i][j][0][0] * s[0][0] + t[i][j][1][1] * s[1][1] +
                  2 * t[i][j][0][1] * s[0][1];
}



/**
 * 在等级4和等级2的对称张量之间进行双重收缩，产生等级2的对称张量，作为这个函数的第一个参数给出。这个操作是矩阵-向量乘法的对称张量类似物。
 * 这个函数的作用与 SymmetricTensor::operator*().
 * 相同，但是不应该使用它，因为成员运算器知道实际的数据存储格式，并且至少快两个数量级。这个函数的存在主要是为了与一般的张量类兼容。
 * @relatesalso  SymmetricTensor
 *
 *
 */
template <typename Number, typename OtherNumber>
constexpr inline void double_contract(
  SymmetricTensor<2, 2, typename ProductType<Number, OtherNumber>::type> &tmp,
  const SymmetricTensor<2, 2, Number> &                                   s,
  const SymmetricTensor<4, 2, OtherNumber> &                              t)
{
  const unsigned int dim = 2;

  for (unsigned int i = 0; i < dim; ++i)
    for (unsigned int j = i; j < dim; ++j)
      tmp[i][j] = s[0][0] * t[0][0][i][j] * +s[1][1] * t[1][1][i][j] +
                  2 * s[0][1] * t[0][1][i][j];
}



/**
 * 在等级4和等级2的对称张量之间进行双重收缩，产生等级2的对称张量，作为这个函数的第一个参数给出。这个操作是矩阵-向量乘法的对称张量类似物。
 * 这个函数的作用与 SymmetricTensor::operator*().
 * 相同，但是不应该使用它，因为成员运算器知道实际的数据存储格式，并且至少要快两个数量级。这个函数的存在主要是为了与一般的张量类兼容。
 * @relatesalso  SymmetricTensor
 *
 *
 */
template <typename Number, typename OtherNumber>
constexpr inline void double_contract(
  SymmetricTensor<2, 3, typename ProductType<Number, OtherNumber>::type> &tmp,
  const SymmetricTensor<4, 3, Number> &                                   t,
  const SymmetricTensor<2, 3, OtherNumber> &                              s)
{
  const unsigned int dim = 3;

  for (unsigned int i = 0; i < dim; ++i)
    for (unsigned int j = i; j < dim; ++j)
      tmp[i][j] = t[i][j][0][0] * s[0][0] + t[i][j][1][1] * s[1][1] +
                  t[i][j][2][2] * s[2][2] + 2 * t[i][j][0][1] * s[0][1] +
                  2 * t[i][j][0][2] * s[0][2] + 2 * t[i][j][1][2] * s[1][2];
}



/**
 * 在等级4和等级2的对称张量之间进行双重收缩，产生等级2的对称张量，作为这个函数的第一个参数给出。这个操作是矩阵-向量乘法的对称张量类似物。
 * 这个函数的作用与 SymmetricTensor::operator*().
 * 相同，但是不应该使用它，因为成员运算器知道实际的数据存储格式，并且至少要快两个数量级。这个函数的存在主要是为了与一般的张量类兼容。
 * @relatesalso  SymmetricTensor
 *
 *
 */
template <typename Number, typename OtherNumber>
constexpr inline void double_contract(
  SymmetricTensor<2, 3, typename ProductType<Number, OtherNumber>::type> &tmp,
  const SymmetricTensor<2, 3, Number> &                                   s,
  const SymmetricTensor<4, 3, OtherNumber> &                              t)
{
  const unsigned int dim = 3;

  for (unsigned int i = 0; i < dim; ++i)
    for (unsigned int j = i; j < dim; ++j)
      tmp[i][j] = s[0][0] * t[0][0][i][j] + s[1][1] * t[1][1][i][j] +
                  s[2][2] * t[2][2][i][j] + 2 * s[0][1] * t[0][1][i][j] +
                  2 * s[0][2] * t[0][2][i][j] + 2 * s[1][2] * t[1][2][i][j];
}



/**
 * 将一个对称的秩2张量（即一个矩阵）乘以一个秩1张量（即一个向量）。其结果是一个秩-1张量（即一个矢量）。
 * @relatesalso  SymmetricTensor
 *
 *
 */
template <int dim, typename Number, typename OtherNumber>
constexpr Tensor<1, dim, typename ProductType<Number, OtherNumber>::type>
operator*(const SymmetricTensor<2, dim, Number> &src1,
          const Tensor<1, dim, OtherNumber> &    src2)
{
  Tensor<1, dim, typename ProductType<Number, OtherNumber>::type> dest;
  for (unsigned int i = 0; i < dim; ++i)
    for (unsigned int j = 0; j < dim; ++j)
      dest[i] += src1[i][j] * src2[j];
  return dest;
}


/**
 * 将一个秩-1张量（即一个向量）与一个对称的秩-2张量（即一个矩阵）相乘。其结果是一个秩-1张量（即一个向量）。
 * @relatesalso  SymmetricTensor
 *
 *
 */
template <int dim, typename Number, typename OtherNumber>
constexpr Tensor<1, dim, typename ProductType<Number, OtherNumber>::type>
operator*(const Tensor<1, dim, Number> &              src1,
          const SymmetricTensor<2, dim, OtherNumber> &src2)
{
  // this is easy for symmetric tensors:
  return src2 * src1;
}



/**
 * 张量的点乘（单一收缩）。返回一个等级为 $(\text{rank}_1 +
 * \text{rank}_2
 *
 * - 2)$ 的张量，它是等级为 @p src1 的张量的最后一个索引与等级为 @p rank_2: 的张量的第一个索引的收缩@f[
 * \text{result}_{i_1,\ldots,i_{r1},j_1,\ldots,j_{r2}}
 * = \sum_{k}
 *   \text{left}_{i_1,\ldots,i_{r1}, k}
 *   \text{right}_{k, j_1,\ldots,j_{r2}}
 * @f] 。
 *
 * @note
 * 由于一个操作数是张量，乘法运算符只对一对索引执行收缩。这与SymmetricTensor的乘法运算符不同，后者做的是双倍收缩。
 * @relatesalso  SymmetricTensor
 *
 *
 */
template <int rank_1,
          int rank_2,
          int dim,
          typename Number,
          typename OtherNumber>
constexpr DEAL_II_ALWAYS_INLINE
  typename Tensor<rank_1 + rank_2 - 2,
                  dim,
                  typename ProductType<Number, OtherNumber>::type>::tensor_type
  operator*(const Tensor<rank_1, dim, Number> &              src1,
            const SymmetricTensor<rank_2, dim, OtherNumber> &src2)
{
  return src1 * Tensor<rank_2, dim, OtherNumber>(src2);
}



/**
 * 张量的点积（单收缩）。返回一个等级为 $(\text{rank}_1 +
 * \text{rank}_2
 *
 * - 2)$ 的张量，它是等级为 @p rank_1 的张量 @p src1 的最后一个索引与等级为 @p rank_2: 的张量@f[
 * \text{result}_{i_1,\ldots,i_{r1},j_1,\ldots,j_{r2}}
 * = \sum_{k}
 *   \text{left}_{i_1,\ldots,i_{r1}, k}
 *   \text{right}_{k, j_1,\ldots,j_{r2}}
 * @f]的第一个索引的收缩。
 *
 * @note
 * 由于一个操作数是张量，乘法运算符只对一对索引执行收缩。这与SymmetricTensor的乘法运算符不同，后者做的是双倍收缩。
 * @relatesalso  SymmetricTensor
 *
 *
 */
template <int rank_1,
          int rank_2,
          int dim,
          typename Number,
          typename OtherNumber>
constexpr DEAL_II_ALWAYS_INLINE
  typename Tensor<rank_1 + rank_2 - 2,
                  dim,
                  typename ProductType<Number, OtherNumber>::type>::tensor_type
  operator*(const SymmetricTensor<rank_1, dim, Number> &src1,
            const Tensor<rank_2, dim, OtherNumber> &    src2)
{
  return Tensor<rank_1, dim, Number>(src1) * src2;
}



/**
 * 秩为2的对称张量的输出运算符。连续打印元素，中间有一个空格，等级1的子张量之间有两个空格，等级2之间有三个空格，以此类推。在输出中不做特殊的修改来表示对称性，例如只输出唯一的条目。
 * @relatesalso  SymmetricTensor
 *
 *
 */
template <int dim, typename Number>
inline std::ostream &
operator<<(std::ostream &out, const SymmetricTensor<2, dim, Number> &t)
{
  // make our lives a bit simpler by outputting
  // the tensor through the operator for the
  // general Tensor class
  Tensor<2, dim, Number> tt;

  for (unsigned int i = 0; i < dim; ++i)
    for (unsigned int j = 0; j < dim; ++j)
      tt[i][j] = t[i][j];

  return out << tt;
}



/**
 * 秩为4的对称张量的输出运算符。连续打印元素，中间有一个空格，等级1的子张量之间有两个空格，等级2之间有三个空格，以此类推。在输出中不做特殊的修改以表示对称性，例如只输出唯一的条目。
 * @relatesalso  SymmetricTensor
 *
 *
 */
template <int dim, typename Number>
inline std::ostream &
operator<<(std::ostream &out, const SymmetricTensor<4, dim, Number> &t)
{
  // make our lives a bit simpler by outputting
  // the tensor through the operator for the
  // general Tensor class
  Tensor<4, dim, Number> tt;

  for (unsigned int i = 0; i < dim; ++i)
    for (unsigned int j = 0; j < dim; ++j)
      for (unsigned int k = 0; k < dim; ++k)
        for (unsigned int l = 0; l < dim; ++l)
          tt[i][j][k][l] = t[i][j][k][l];

  return out << tt;
}


DEAL_II_NAMESPACE_CLOSE

#endif


