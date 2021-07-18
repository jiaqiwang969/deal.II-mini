//include/deal.II-translator/base/template_constraints_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2021 by the deal.II authors
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

#ifndef dealii_template_constraints_h
#define dealii_template_constraints_h


#include <deal.II/base/config.h>

#include <deal.II/base/complex_overloads.h>

#include <complex>
#include <iterator>
#include <utility>

DEAL_II_NAMESPACE_OPEN

namespace internal
{
  namespace TemplateConstraints
  {
    // helper struct for is_base_of_all and all_same_as
    template <bool... Values>
    struct BoolStorage;


    /**
     * 一个辅助类，其`value`成员是真还是假，取决于给定的布尔模板参数是否全部为真。
     *
     */
    template <bool... Values>
    struct all_true
    {
      static constexpr bool value =
        std::is_same<BoolStorage<Values..., true>,
                     BoolStorage<true, Values...>>::value;
    };
  } // namespace TemplateConstraints
} // namespace internal

/**
 * 这个结构是 std::is_base_of<Base,
 * Derived>对模板参数包的泛化，测试所有Derived...类是否以Base为基类，或者本身就是Base。结果存储在成员变量值中。
 *
 *
 */
template <class Base, class... Derived>
struct is_base_of_all
{
  static constexpr bool value = internal::TemplateConstraints::all_true<
    std::is_base_of<Base, Derived>::value...>::value;
};



/**
 * 这个结构是  std::is_same  对模板参数包的泛化，测试
 * `Types...`
 * 参数包中的所有类型是否等于作为第一个模板参数的
 * `Type`。结果存储在成员变量值中。
 *
 *
 */
template <class Type, class... Types>
struct all_same_as
{
  static constexpr bool value = internal::TemplateConstraints::all_true<
    std::is_same<Type, Types>::value...>::value;
};



/*  `std::enable_if` 的概括，只有在给定的布尔模板参数中的<i>all</i>为真时才起作用。

* 
*
*/
template <bool... Values>
struct enable_if_all
  : std::enable_if<internal::TemplateConstraints::all_true<Values...>::value>
{};



/**
 * 一个类型特征，检查一个类是否表现为一个有开始和结束的可迭代容器。这意味着该类要么定义了`begin()`和`end()`函数，要么是一个C风格的数组。
 *
 *
 */
template <typename T>
class has_begin_and_end
{
  template <typename C>
  static std::false_type
  test(...);

  template <typename C>
  static auto
  test(int) -> decltype(std::begin(std::declval<C>()),
                        std::end(std::declval<C>()),
                        std::true_type());

public:
  using type = decltype(test<T>(0));

  static const bool value = type::value;
};



/**
 * 一个模板类，简单地将其模板参数导出为一个本地别名。这个类，虽然一开始看起来毫无用处，但在以下情况下是有意义的：如果你有一个如下的函数模板。
 *
 * @code
 * template <typename T>
 * void f(T, T);
 * @endcode
 * 那么它就不能在 <code>f(1, 3.141)</code>
 * 这样的表达式中被调用，因为模板的类型 <code>T</code>
 * 不能以唯一的方式从参数的类型中推导出来。然而，如果模板被写成
 *
 * @code
 * template <typename T>
 * void f(T, typename identity<T>::type);
 * @endcode
 * 那么这个调用就变得有效了： <code>T</code>
 * 的类型不能从函数的第二个参数中推导出来，所以只有第一个参数参与了模板类型解析。
 * 这个特征的背景如下：考虑
 *
 * @code
 * template <typename RT, typename A>
 * void forward_call(RT (*p) (A), A a)
 * {
 * p(a);
 * }
 *
 * void h (double);
 *
 * void g()
 * {
 * forward_call(&h, 1);
 * }
 * @endcode
 * 这段代码不能编译，因为编译器不能决定模板类型
 * <code>A</code> should be <code>double</code>
 * （来自作为第一参数给定的函数的签名，因为表达式
 * <code>1</code>
 * 有该类型。当然，我们希望编译器做的是简单地将
 * <code>1</code> to <code>double</code>
 * 投入。我们可以通过编写以下代码来实现这一点。
 *
 * @code
 * template <typename RT, typename A>
 * void forward_call(RT (*p) (A), typename identity<A>::type a)
 * {
 * p(a);
 * }
 *
 * void h (double);
 *
 * void g()
 * {
 * forward_call(&h, 1);
 * }
 * @endcode
 *
 *
 *
 */
template <typename T>
struct identity
{
  using type = T;
};



/**
 * 一个总是返回一个给定值的类。这需要作为一些编译器难以处理的用作默认参数的lambdas的变通方法。
 *
 *
 */
template <typename ArgType, typename ValueType>
struct always_return
{
  ValueType value;
  ValueType
  operator()(const ArgType &)
  {
    return value;
  }
};



/**
 * 一个用于对任意指针进行平等比较的类。在某些情况下，人们希望确保一个函数的两个参数不是同一个对象。在这种情况下，我们会确保它们的地址不一样。然而，有时这两个参数的类型可能是模板类型，它们可能是同一类型，也可能不是。在这种情况下，像<tt>&object1
 * !=
 * &object2</tt>这样的简单比较只有在两个对象的类型相同的情况下才起作用，但如果不相同，编译器会barf。然而，在后一种情况下，由于两个对象的类型不同，我们可以确定这两个对象不可能是相同的。
 * 这个类实现了一个比较函数，如果它的两个参数的类型不同，总是返回
 * @p false ，否则返回<tt>p1 == p2</tt>。
 *
 *
 */
struct PointerComparison
{
  /**
   * 相同类型的指针的比较函数。如果两个指针相等，则返回
   * @p true 。
   *
   */
  template <typename T>
  static bool
  equal(const T *p1, const T *p2)
  {
    return (p1 == p2);
  }


  /**
   * 用于不同类型的指针的比较函数。C++语言不允许使用<tt>operator==</tt>来比较这些指针。
   * 然而，由于这两个指针的类型不同，我们知道它们不可能相同，所以我们总是返回
   * @p false.  。
   *
   */
  template <typename T, typename U>
  static bool
  equal(const T *, const U *)
  {
    return false;
  }
};



namespace internal
{
  /**
   * 一个结构，实现了两种类型相乘产生的默认乘积类型。
   * @note 当 @p T 或 @p U 有限定词（ @p const 或 @p volatile)
   * 或者是 @p lvalue 或 @p rvalue 的引用时，应该注意!
   * 建议只对未限定的（完全剥离的）类型进行该类的专业化处理，并使用ProductType类来确定对（潜在的）限定类型进行操作的结果。
   *
   */
  template <typename T, typename U>
  struct ProductTypeImpl
  {
    using type = decltype(std::declval<T>() * std::declval<U>());
  };

} // namespace internal



/**
 * 一个具有本地别名的类，它代表了由类型为 @p T 和 @p U.
 * 的两个变量的乘积所产生的类型。换句话说，我们想在这样的代码中推断
 * <code>product</code> 变量的类型。
 *
 * @code
 * T t;
 * U u;
 * auto product = t*u;
 * @endcode
 * 这个结构的本地别名代表了变量 <code>product</code>
 * 会有的类型。
 *
 *  <h3>Where is this useful</h3>
 * 这个类的目的主要是表示人们需要用来表示正交点的有限元场的值或梯度的类型。例如，假设你在Vector<float>中存储未知数的值
 * $U_j$ ，那么在正交点评估 $u_h(x_q) = \sum_j U_j \varphi_j(x_q)$
 * 的结果是需要存储为 $u_h(x_q)$ 的变量，因为 $U_j$ 是 @p
 * float 值， $\varphi_j(x_q)$ 被计算为 @p double 值，然后积为 @p
 * double 值。另一方面，如果你将未知数 $U_j$ 存储为
 * <code>std::complex@<double@></code> 值，并试图在正交点评估
 * $\nabla u_h(x_q) = \sum_j U_j \nabla\varphi_j(x_q)$ ，那么梯度 $\nabla
 * u_h(x_q)$ 需要存储为
 * <code>Tensor@<1,dim,std::complex@<double@>@></code>
 * 类型的对象，因为当你用复数乘以 <code>Tensor@<1,dim@></code>
 * （用于表示标量有限元的形状函数的梯度的类型）时，你会得到这样的对象。
 * 同样，如果你使用的是矢量值元素（有dim成分），并且
 * $U_j$ 被存储为 @p double 变量，那么 $u_h(x_q) = \sum_j U_j
 * \varphi_j(x_q)$ 需要有 <code>Tensor@<1,dim@></code>
 * 类型（因为形状函数有 <code>Tensor@<1,dim@></code>
 * 类型）。最后，如果你将 $U_j$ 存储为类型为
 * <code>std::complex@<double@></code>
 * 的对象，并且你有一个矢量值的元素，那么梯度 $\nabla
 * u_h(x_q) = \sum_j U_j \nabla\varphi_j(x_q)$ 将导致类型为
 * <code>Tensor@<2,dim,std::complex@<double@>  @></code>. 的对象。
 * 在所有这些情况下，这个类型是用来识别哪种类型需要用于计算未知数与值、梯度或形状函数的其他属性的乘积的结果。
 *
 *
 */
template <typename T, typename U>
struct ProductType
{
  using type =
    typename internal::ProductTypeImpl<typename std::decay<T>::type,
                                       typename std::decay<U>::type>::type;
};

namespace internal
{
  // Annoyingly, there is no std::complex<T>::operator*(U) for scalars U
  // other than T (not even in C++11, or C++14). We provide our own overloads
  // in base/complex_overloads.h, but in order for them to work, we have to
  // manually specify all products we want to allow:

  template <typename T>
  struct ProductTypeImpl<std::complex<T>, std::complex<T>>
  {
    using type = std::complex<T>;
  };

  template <typename T, typename U>
  struct ProductTypeImpl<std::complex<T>, std::complex<U>>
  {
    using type = std::complex<typename ProductType<T, U>::type>;
  };

  template <typename U>
  struct ProductTypeImpl<double, std::complex<U>>
  {
    using type = std::complex<typename ProductType<double, U>::type>;
  };

  template <typename T>
  struct ProductTypeImpl<std::complex<T>, double>
  {
    using type = std::complex<typename ProductType<T, double>::type>;
  };

  template <typename U>
  struct ProductTypeImpl<float, std::complex<U>>
  {
    using type = std::complex<typename ProductType<float, U>::type>;
  };

  template <typename T>
  struct ProductTypeImpl<std::complex<T>, float>
  {
    using type = std::complex<typename ProductType<T, float>::type>;
  };

} // namespace internal



/**
 * 该类提供了一个本地别名 @p type
 * ，该别名等于模板参数，但只有在模板参数对应于标量类型（即浮点类型、有符号或无符号整数或复数中的一种）的情况下才是如此。如果模板类型
 * @p T 不是标量，那么就没有声明类
 * <code>EnableIfScalar@<T@></code> ，因此也就没有本地别名可用。
 * 该类的目的是如果其中一个参数不是标量数字，则禁用某些模板函数。通过（无意义的）例子，考虑以下函数。
 * @code
 * template <typename T>
 * T multiply (const T t1, const T t2)
 * {
 *   return t1*t2;
 * }
 * @endcode
 * 这个函数可以用任何两个相同类型的参数来调用  @p T.
 * 这包括一些参数，这显然没有意义。因此，人们可能想把这个函数限制在只有标量，这可以写成
 *
 * @code
 * template <typename T>
 * typename EnableIfScalar<T>::type
 * multiply (const T t1, const T t2)
 * {
 *   return t1*t2;
 * }
 * @endcode
 * 在你调用函数的地方，编译器会从参数中推断出类型 @p T
 * 。例如，在
 *
 * @code
 * multiply(1.234, 2.345);
 * @endcode
 * 它将推断出 @p T 是 @p double, ，由于
 * <code>EnableIfScalar@<double@>::%type</code> 等于 @p double,
 * ，编译器将从上面的模板实例化一个函数<code>double
 * multiply(const double, const
 * double)</code>。另一方面，在这样的背景下
 *
 * @code
 * std::vector<char> v1, v2;
 * multiply(v1, v2);
 * @endcode
 * 编译器会推断出 @p T 是 <code>std::vector@<char@></code>
 * ，但由于 <code>EnableIfScalar@<std::vector@<char@>@>::%type</code>
 * 不存在，编译器不会考虑模板的实例化。这种技术被称为
 * "替换失败不是错误（SFINAE）"。它确保了模板函数甚至不能被调用，而不是导致后来的错误，即操作
 * <code>t1*t2</code>
 * 没有被定义（或可能导致一些无意义的结果）。它还允许声明一个函数的重载，如
 * @p
 * 乘以不同类型的参数，而不会导致编译器产生模糊的调用错误。
 *
 *
 */
template <typename T>
struct EnableIfScalar;


template <>
struct EnableIfScalar<double>
{
  using type = double;
};

template <>
struct EnableIfScalar<float>
{
  using type = float;
};

template <>
struct EnableIfScalar<long double>
{
  using type = long double;
};

template <>
struct EnableIfScalar<int>
{
  using type = int;
};

template <>
struct EnableIfScalar<unsigned int>
{
  using type = unsigned int;
};

template <typename T>
struct EnableIfScalar<std::complex<T>>
{
  using type = std::complex<T>;
};


DEAL_II_NAMESPACE_CLOSE

#endif


