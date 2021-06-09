//include/deal.II-translator/base/numbers_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2020 by the deal.II authors
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

#ifndef dealii_numbers_h
#define dealii_numbers_h


#include <deal.II/base/config.h>

#include <deal.II/base/types.h>

#ifdef DEAL_II_COMPILER_CUDA_AWARE
#  include <cuComplex.h>
#endif

#include <cmath>
#include <complex>
#include <cstddef>
#include <type_traits>

#ifdef DEAL_II_COMPILER_CUDA_AWARE
#  define DEAL_II_CUDA_HOST_DEV __host__ __device__
#else
#  define DEAL_II_CUDA_HOST_DEV
#endif

DEAL_II_NAMESPACE_OPEN

namespace internal
{
  /**
   * 一个辅助类，指定指定数据类型Number的VectorizedArray的最大向量长度，用于给定的处理器架构和优化级别。
   * 最大向量长度的值被用作VectorizedArray的默认模板参数，这样VectorizedArray<Number>就相当于VectorizedArray<Number,
   * VectorizedArrayWidthSpecifier<Number>::max_width>.  。
   * @note  该类是不支持矢量化的数据类型的默认实现。
   * @tparam  Number
   * 想找出硬件支持的向量的最大长度的基础数据类型。
   *
   */
  template <typename Number>
  struct VectorizedArrayWidthSpecifier
  {
    /**
     * 任意类型的VectorizedArray的最大向量长度。
     *
     */
    constexpr static unsigned int max_width = 1;
  };

  /**
   * 一个辅助类，指定数据类型`double`的VectorizedArray的最大向量长度，用于给定的处理器架构和优化级别。关于支持的最大向量长度的详细描述，请参见VectorizedArray的文档。
   *
   */
  template <>
  struct VectorizedArrayWidthSpecifier<double>
  {
    /**
     * Double的VectorizedArray的最大向量长度。
     *
     */
    constexpr static unsigned int max_width =
#if DEAL_II_VECTORIZATION_WIDTH_IN_BITS >= 512
      8;
#elif DEAL_II_VECTORIZATION_WIDTH_IN_BITS >= 256
      4;
#elif DEAL_II_VECTORIZATION_WIDTH_IN_BITS >= 128
      2;
#else
      1;
#endif
  };

  /**
   * 一个辅助类，指定数据类型 "float
   * "的VectorizedArray的最大向量长度，适用于给定的处理器架构和优化级别。关于支持的最大向量长度的详细描述，请参见VectorizedArray的文档。
   *
   */
  template <>
  struct VectorizedArrayWidthSpecifier<float>
  {
    /**
     * VectorizedArray的最大向量长度为float。
     *
     */
    constexpr static unsigned int max_width =
#if DEAL_II_VECTORIZATION_WIDTH_IN_BITS >= 128 && defined(__ALTIVEC__)
      4;
#elif DEAL_II_VECTORIZATION_WIDTH_IN_BITS >= 512 && defined(__AVX512F__)
      16;
#elif DEAL_II_VECTORIZATION_WIDTH_IN_BITS >= 256 && defined(__AVX__)
      8;
#elif DEAL_II_VECTORIZATION_WIDTH_IN_BITS >= 128 && defined(__SSE2__)
      4;
#else
      1;
#endif
  };


} // namespace internal

// forward declarations to support abs or sqrt operations on VectorizedArray
#ifndef DOXYGEN
template <typename Number,
          std::size_t width =
            internal::VectorizedArrayWidthSpecifier<Number>::max_width>
class VectorizedArray;
template <typename T>
struct EnableIfScalar;
#endif

DEAL_II_NAMESPACE_CLOSE

// Declare / Import auto-differentiable math functions in(to) standard
// namespace before numbers::NumberTraits is defined
#ifdef DEAL_II_WITH_ADOLC
#  include <deal.II/differentiation/ad/adolc_math.h>

#  include <adolc/adouble.h> // Taped double
#endif
// Ideally we'd like to #include <deal.II/differentiation/ad/sacado_math.h>
// but header indirectly references numbers.h. We therefore simply
// import the whole Sacado header at this point to get the math
// functions imported into the standard namespace.
#ifdef DEAL_II_TRILINOS_WITH_SACADO
#  include <Sacado.hpp>
#endif

namespace std
{
  template <typename Number, std::size_t width>
  DEAL_II_ALWAYS_INLINE ::dealii::VectorizedArray<Number, width>
  sqrt(const ::dealii::VectorizedArray<Number, width> &);
  template <typename Number, std::size_t width>
  DEAL_II_ALWAYS_INLINE ::dealii::VectorizedArray<Number, width>
  abs(const ::dealii::VectorizedArray<Number, width> &);
  template <typename Number, std::size_t width>
  DEAL_II_ALWAYS_INLINE ::dealii::VectorizedArray<Number, width>
  max(const ::dealii::VectorizedArray<Number, width> &,
      const ::dealii::VectorizedArray<Number, width> &);
  template <typename Number, std::size_t width>
  DEAL_II_ALWAYS_INLINE ::dealii::VectorizedArray<Number, width>
  min(const ::dealii::VectorizedArray<Number, width> &,
      const ::dealii::VectorizedArray<Number, width> &);
  template <typename Number, size_t width>
  ::dealii::VectorizedArray<Number, width>
  pow(const ::dealii::VectorizedArray<Number, width> &, const Number p);
  template <typename Number, size_t width>
  ::dealii::VectorizedArray<Number, width>
  sin(const ::dealii::VectorizedArray<Number, width> &);
  template <typename Number, size_t width>
  ::dealii::VectorizedArray<Number, width>
  cos(const ::dealii::VectorizedArray<Number, width> &);
  template <typename Number, size_t width>
  ::dealii::VectorizedArray<Number, width>
  tan(const ::dealii::VectorizedArray<Number, width> &);
  template <typename Number, size_t width>
  ::dealii::VectorizedArray<Number, width>
  exp(const ::dealii::VectorizedArray<Number, width> &);
  template <typename Number, size_t width>
  ::dealii::VectorizedArray<Number, width>
  log(const ::dealii::VectorizedArray<Number, width> &);
} // namespace std

DEAL_II_NAMESPACE_OPEN

/**
 * 用于声明通用常数的命名空间。由于<tt>math.h</tt>中的可用性并不总是被保证，我们把它们放在这里。由于这个文件被<tt>base/config.h</tt>所包含，它们对整个库是可用的。
 * 这里定义的常量是有时在系统包含文件<tt>math.h</tt>中声明的<tt>M_XXX</tt>常量的一个子集，但没有前缀<tt>M_</tt>。
 * 除此之外，我们声明<tt>invalid_unsigned_int</tt>为可表示的最大无符号整数；这个值在库中被广泛用作无效索引、无效数组大小的标记，以及类似的目的。
 *
 *
 */
namespace numbers
{
  /**
   * e
   *
   */
  static constexpr double E = 2.7182818284590452354;

  /**
   * log_2 e
   *
   */
  static constexpr double LOG2E = 1.4426950408889634074;

  /**
   * log_10 e
   *
   */
  static constexpr double LOG10E = 0.43429448190325182765;

  /**
   * log_e 2
   *
   */
  static constexpr double LN2 = 0.69314718055994530942;

  /**
   * 日志_e 10
   *
   */
  static constexpr double LN10 = 2.30258509299404568402;

  /**
   * 圆周率
   *
   */
  static constexpr double PI = 3.14159265358979323846;

  /**
   * pi/2
   *
   */
  static constexpr double PI_2 = 1.57079632679489661923;

  /**
   * 圆周率/4
   *
   */
  static constexpr double PI_4 = 0.78539816339744830962;

  /**
   * sqrt(2)
   *
   */
  static constexpr double SQRT2 = 1.41421356237309504880;

  /**
   * 1/sqrt(2)
   *
   */
  static constexpr double SQRT1_2 = 0.70710678118654752440;

  /**
   * 检查给定的类型是否可以在CUDA设备代码中使用。
   * 如果不能，DEAL_II_CUDA_HOST_DEV需要在使用该类型的函数中被禁用。
   *
   */
  template <typename Number, typename = void>
  struct is_cuda_compatible : std::true_type
  {};

  /**
   * std::complex  不能在CUDA设备代码中使用。
   *
   */
  template <typename Number>
  struct is_cuda_compatible<std::complex<Number>, void> : std::false_type
  {};

  /**
   * 如果给定值是一个有限的浮点数，即既不是正负无穷大，也不是NaN（不是一个数字），则返回
   * @p true 。    注意，这个函数的参数类型是
   * <code>double</code>
   * 。换句话说，如果你给出一个非常大的<code>long
   * double</code>类型的数字，这个函数可能会返回
   * <code>false</code>  ，即使这个数字相对于  <code>long
   * double</code>  类型来说是有限的。
   *
   */
  bool
  is_finite(const double x);

  /**
   * 如果给定复数的实部和虚部是有限的，则返回  @p true 。
   *
   */
  bool
  is_finite(const std::complex<double> &x);

  /**
   * 如果给定复数的实部和虚部是有限的，返回 @p true 。
   *
   */
  bool
  is_finite(const std::complex<float> &x);

  /**
   * 如果给定复数的实部和虚部是有限的，则返回  @p true 。
   * 如果实部或虚部是非常大的数字，就 <code>double</code>
   * 而言是无限的，但就 <code>long double</code>
   * 而言是有限的，也可能无法正常工作。
   *
   */
  bool
  is_finite(const std::complex<long double> &x);

  /**
   * 返回两个数字是否互相相等。
   * 对于复杂的数据类型（例如一些自动可分的数字），该函数返回输入参数存储的标量值的比较结果。
   * @note 这个函数期望 @p value_2 可以投射到 @p value_1.
   * 的类型。
   *
   */
  template <typename Number1, typename Number2>
  constexpr bool
  values_are_equal(const Number1 &value_1, const Number2 &value_2);

  /**
   * 返回两个数字是否相互不相等。
   * 对于复杂的数据类型（例如一些自动可分的数字），该函数返回输入参数存储的标量值的比较结果。
   * @note 这个函数期望 @p value_2 可以投射到 @p value_1.
   * 的类型。
   *
   */
  template <typename Number1, typename Number2>
  bool
  values_are_not_equal(const Number1 &value_1, const Number2 &value_2);

  /**
   * 返回一个值是否等于零。
   * 对于复杂的数据类型（例如一些自动可分的数字），该函数返回输入参数存储的标量值的比较结果。
   *
   */
  template <typename Number>
  constexpr bool
  value_is_zero(const Number &value);

  /**
   * 返回  @p value_1  是否小于  @p value_2.
   * 对于复杂的数据类型（例如一些自动可分的数字），此函数返回由输入参数存储的标量值的比较结果。
   * @note 这个函数期望 @p value_2 可以投射到 @p value_1.
   * 的类型。
   *
   */
  template <typename Number1, typename Number2>
  bool
  value_is_less_than(const Number1 &value_1, const Number2 &value_2);

  /**
   * 返回 @p value_1 是否小于或等于 @p value_2.
   * 对于复杂的数据类型（例如一些自动可分的数字），该函数返回输入参数存储的标量值的比较结果。
   * @note 这个函数期望 @p value_2 可投到 @p value_1. 的类型。
   *
   */
  template <typename Number1, typename Number2>
  bool
  value_is_less_than_or_equal_to(const Number1 &value_1,
                                 const Number2 &value_2);



  /**
   * 返回 @p value_1 是否大于 @p value_2.
   * 对于复杂的数据类型（例如一些自动可分的数字），该函数返回输入参数存储的标量值的比较结果。
   * @note 这个函数期望 @p value_2 可投到 @p value_1. 的类型。
   *
   */
  template <typename Number1, typename Number2>
  bool
  value_is_greater_than(const Number1 &value_1, const Number2 &value_2);

  /**
   * 返回 @p value_1 是否大于或等于 @p value_2.
   * 对于复杂的数据类型（例如一些自动可分的数字），该函数返回输入参数存储的标量值的比较结果。
   * @note 这个函数期望 @p value_2 可以投射到 @p value_1.
   * 的类型。
   *
   */
  template <typename Number1, typename Number2>
  bool
  value_is_greater_than_or_equal_to(const Number1 &value_1,
                                    const Number2 &value_2);

  /**
   * 一个结构，连同其部分特殊化 NumberTraits<std::complex<number>
   * >，提供了特质和成员函数，使得编写同时适用于实数类型和复数类型的模板成为可能。这个模板主要用于实现线性代数类，如向量和矩阵，对实数和复数都有效。
   *
   */
  template <typename number>
  struct NumberTraits
  {
    /**
     * 一个标志，指定给这个类的模板类型是复数还是实数。因为一般的模板是为非复数类型选择的，所以答案是
     * <code>false</code>  。
     *
     */
    static constexpr bool is_complex = false;

    /**
     * 对于这个数据类型，别名为相应的实数类型。由于一般模板是为所有不是
     * std::complex<T>,
     * 的特殊化的数据类型选择的，底层类型必须是实值，所以real_type等于底层类型。
     *
     */
    using real_type = number;

    /**
     * 返回给定数字的复数共轭值。因为如果number不是复数数据类型，就会选择一般的模板，这个函数只是返回给定的数字。
     * @note 这个函数也可以在CUDA设备代码中使用。
     *
     */
    static constexpr DEAL_II_CUDA_HOST_DEV const number &
                                                 conjugate(const number &x);

    /**
     * 返回给定数字的绝对值的平方。由于一般模板选择的是不等于
     * std::complex,
     * 的类型，这个函数只是返回给定数字的平方。
     * @note
     * 如果模板类型可以在CUDA设备代码中使用，这个函数也同样适用。
     *
     */
    template <typename Dummy = number>
    static constexpr DEAL_II_CUDA_HOST_DEV
      typename std::enable_if<std::is_same<Dummy, number>::value &&
                                is_cuda_compatible<Dummy>::value,
                              real_type>::type
      abs_square(const number &x);

    template <typename Dummy = number>
    static constexpr
      typename std::enable_if<std::is_same<Dummy, number>::value &&
                                !is_cuda_compatible<Dummy>::value,
                              real_type>::type
      abs_square(const number &x);

    /**
     * 返回一个数字的绝对值。
     *
     */
    static real_type
    abs(const number &x);
  };


  /**
   * 一般NumberTraits类的特殊化，如果底层数据类型为
   * std::complex<T>. ，则提供相关信息。
   *
   */
  template <typename number>
  struct NumberTraits<std::complex<number>>
  {
    /**
     * 一个标志，指定给这个类的模板类型是复数还是实数。由于一般模板的这种特殊化是为复杂类型选择的，所以答案是
     * <code>true</code>  。
     *
     */
    static constexpr bool is_complex = true;

    /**
     * 对于这个数据类型，别名为相应的实数类型。由于模板的这种特殊化是为数字类型选择的
     * std::complex<T>,
     * ，所以实数类型等于用于存储复数的两个分量的类型。
     *
     */
    using real_type = number;

    /**
     * 返回给定数的复数共轭值。
     *
     */
    static constexpr std::complex<number>
    conjugate(const std::complex<number> &x);

    /**
     * 返回给定数的绝对值的平方。由于一般模板的这种特殊化是为等于
     * std::complex,
     * 的类型选择的，这个函数返回一个数与它的复数共轭的乘积。
     *
     */
    static constexpr real_type
    abs_square(const std::complex<number> &x);


    /**
     * 返回一个复数的绝对值。
     *
     */
    static real_type
    abs(const std::complex<number> &x);
  };

  // --------------- inline and template functions ---------------- //

  inline bool
  is_nan(const double x)
  {
    return std::isnan(x);
  }



  inline bool
  is_finite(const double x)
  {
    return std::isfinite(x);
  }



  inline bool
  is_finite(const std::complex<double> &x)
  {
    // Check complex numbers for infinity
    // by testing real and imaginary part
    return (is_finite(x.real()) && is_finite(x.imag()));
  }



  inline bool
  is_finite(const std::complex<float> &x)
  {
    // Check complex numbers for infinity
    // by testing real and imaginary part
    return (is_finite(x.real()) && is_finite(x.imag()));
  }



  inline bool
  is_finite(const std::complex<long double> &x)
  {
    // Same for std::complex<long double>
    return (is_finite(x.real()) && is_finite(x.imag()));
  }


  template <typename number>
  constexpr DEAL_II_CUDA_HOST_DEV const number &
                                        NumberTraits<number>::conjugate(const number &x)
  {
    return x;
  }



  template <typename number>
  template <typename Dummy>
  constexpr DEAL_II_CUDA_HOST_DEV
    typename std::enable_if<std::is_same<Dummy, number>::value &&
                              is_cuda_compatible<Dummy>::value,
                            typename NumberTraits<number>::real_type>::type
    NumberTraits<number>::abs_square(const number &x)
  {
    return x * x;
  }



  template <typename number>
  template <typename Dummy>
  constexpr
    typename std::enable_if<std::is_same<Dummy, number>::value &&
                              !is_cuda_compatible<Dummy>::value,
                            typename NumberTraits<number>::real_type>::type
    NumberTraits<number>::abs_square(const number &x)
  {
    return x * x;
  }



  template <typename number>
  typename NumberTraits<number>::real_type
  NumberTraits<number>::abs(const number &x)
  {
    return std::abs(x);
  }



  template <typename number>
  constexpr std::complex<number>
  NumberTraits<std::complex<number>>::conjugate(const std::complex<number> &x)
  {
    return std::conj(x);
  }



  template <typename number>
  typename NumberTraits<std::complex<number>>::real_type
  NumberTraits<std::complex<number>>::abs(const std::complex<number> &x)
  {
    return std::abs(x);
  }



  template <typename number>
  constexpr typename NumberTraits<std::complex<number>>::real_type
  NumberTraits<std::complex<number>>::abs_square(const std::complex<number> &x)
  {
    return std::norm(x);
  }

} // namespace numbers


// Forward declarations
namespace Differentiation
{
  namespace AD
  {
    namespace internal
    {
      // Defined in differentiation/ad/ad_number_traits.h
      template <typename T>
      struct NumberType;
    } // namespace internal

    // Defined in differentiation/ad/ad_number_traits.h
    template <typename NumberType>
    struct is_ad_number;
  } // namespace AD
} // namespace Differentiation


namespace internal
{
  /**
   * 测试是否有可能将一种数字类型转换为另一种。
   *
   */
  template <typename From, typename To>
  struct is_explicitly_convertible
  {
    // Source: https://stackoverflow.com/a/16944130
  private:
    template <typename T>
    static void f(T);

    template <typename F, typename T>
    static constexpr auto
    test(int) -> decltype(f(static_cast<T>(std::declval<F>())), true)
    {
      return true;
    }

    template <typename F, typename T>
    static constexpr auto
    test(...) -> bool
    {
      return false;
    }

  public:
    static bool const value = test<From, To>(0);
  };

  /*下面的结构需要在一些特殊的数字类型之间进行转换。  也可以参见tensor.h了解另一种特殊化。 
*
*/
  template <typename T>
  struct NumberType
  {
    static constexpr DEAL_II_ALWAYS_INLINE DEAL_II_CUDA_HOST_DEV const T &
                                                                       value(const T &t)
    {
      return t;
    }

    // Below are generic functions that allows an overload for any
    // type U that is transformable to type T. This is particularly
    // useful when needing to cast exotic number types
    // (e.g. auto-differentiable or symbolic numbers) to a floating
    // point one, such as might  happen when converting between tensor
    // types.

    // Type T is constructible from F.
    template <typename F>
    static constexpr DEAL_II_ALWAYS_INLINE DEAL_II_CUDA_HOST_DEV T
                                                                 value(const F &f,
                                                                       typename std::enable_if<
            !std::is_same<typename std::decay<T>::type,
                          typename std::decay<F>::type>::value &&
            std::is_constructible<T, F>::value>::type * = nullptr)
    {
      return T(f);
    }

    // Type T is explicitly convertible (but not constructible) from F.
    template <typename F>
    static constexpr DEAL_II_ALWAYS_INLINE T
                                           value(const F &f,
                                                 typename std::enable_if<
            !std::is_same<typename std::decay<T>::type,
                          typename std::decay<F>::type>::value &&
            !std::is_constructible<T, F>::value &&
            is_explicitly_convertible<const F, T>::value>::type * = nullptr)
    {
      return static_cast<T>(f);
    }

    // Sacado doesn't provide any conversion operators, so we have
    // to extract the value and perform further conversions from there.
    // To be safe, we extend this to other possible AD numbers that
    // might fall into the same category.
    template <typename F>
    static T
    value(const F &f,
          typename std::enable_if<
            !std::is_same<typename std::decay<T>::type,
                          typename std::decay<F>::type>::value &&
            !std::is_constructible<T, F>::value &&
            !is_explicitly_convertible<const F, T>::value &&
            Differentiation::AD::is_ad_number<F>::value>::type * = nullptr)
    {
      return Differentiation::AD::internal::NumberType<T>::value(f);
    }
  };

  template <typename T>
  struct NumberType<std::complex<T>>
  {
    static constexpr const std::complex<T> &
    value(const std::complex<T> &t)
    {
      return t;
    }

    static constexpr std::complex<T>
    value(const T &t)
    {
      return std::complex<T>(t);
    }

    // Facilitate cast from complex<double> to complex<float>
    template <typename U>
    static constexpr std::complex<T>
    value(const std::complex<U> &t)
    {
      return std::complex<T>(NumberType<T>::value(t.real()),
                             NumberType<T>::value(t.imag()));
    }
  };

#ifdef DEAL_II_COMPILER_CUDA_AWARE
  template <>
  struct NumberType<cuComplex>
  {
    static cuComplex
    value(const float t)
    {
      return make_cuComplex(t, 0.f);
    }
  };

  template <>
  struct NumberType<cuDoubleComplex>
  {
    static cuDoubleComplex
    value(const double t)
    {
      return make_cuDoubleComplex(t, 0.);
    }
  };
#endif
} // namespace internal

namespace numbers
{
#ifdef DEAL_II_ADOLC_WITH_ADVANCED_BRANCHING

  /**
   * 返回两个数字是否相互相等。对于复杂的数据类型（例如一些自动可微分的数字），这个函数只返回输入值所存储的标量值是否相等。
   * @note 当ADOL-C被编译为 "高级分支
   * "功能时，那么这个特殊化只用于断言和其他不影响计算最终结果的代码路径中。
   *
   */
  // Defined in differentiation/ad/adolc_number_types.cc
  bool
  values_are_equal(const adouble &value_1, const adouble &value_2);


  /**
   * 返回两个数字是否彼此相等。对于复杂的数据类型（例如一些自动可分的数字），这个函数只返回输入值所存储的标量值是否相等。
   * @note 当ADOL-C被编译为 "高级分支
   * "功能时，那么这个特殊化只用于断言和其他不影响计算最终结果的代码路径中。
   *
   */
  template <typename Number>
  bool
  values_are_equal(const adouble &value_1, const Number &value_2)
  {
    // Use the specialized definition for two ADOL-C taped types
    return values_are_equal(value_1,
                            internal::NumberType<adouble>::value(value_2));
  }


  /**
   * 返回两个数字是否彼此相等。对于复杂的数据类型（例如一些自动可微分的数字），这个函数只返回输入值存储的标量值是否相等。
   * @note 当ADOL-C被编译为 "高级分支
   * "功能时，那么这个特殊化只用于断言和其他不影响计算最终结果的代码路径中。
   *
   */
  template <typename Number>
  bool
  values_are_equal(const Number &value_1, const adouble &value_2)
  {
    // Use the above definition
    return values_are_equal(value_2, value_1);
  }

  /**
   * 返回 @p value_1 是否小于 @p value_2.
   * 对于复杂的数据类型（例如一些自动可分的数字），这个函数返回输入参数存储的标量值的比较结果。
   * @note 当ADOL-C被编译为 "高级分支
   * "功能时，那么这个特殊化只用于断言和其他不影响计算最终结果的代码路径中。
   *
   */
  // Defined in differentiation/ad/adolc_number_types.cc
  bool
  value_is_less_than(const adouble &value_1, const adouble &value_2);


  /**
   * 返回 @p value_1 是否小于 @p value_2.
   * 对于复杂的数据类型（例如一些自动可分的数字），这个函数返回输入参数存储的标量值的比较结果。
   * @note 当ADOL-C被编译为 "高级分支
   * "功能时，那么这个特殊化只用于断言和其他不影响计算最终结果的代码路径中。
   *
   */
  template <typename Number>
  bool
  value_is_less_than(const adouble &value_1, const Number &value_2)
  {
    // Use the specialized definition for two ADOL-C taped types
    return value_is_less_than(value_1,
                              internal::NumberType<adouble>::value(value_2));
  }


  /**
   * 返回 @p value_1 是否小于 @p value_2.
   * 对于复杂的数据类型（例如一些自动可分的数字），该函数返回输入参数存储的标量值的比较结果。
   * @note  当ADOL-C被编译为 "高级分支
   * "功能时，那么这个特殊化只用于断言和其他不影响计算最终结果的代码路径中。
   *
   */
  template <typename Number>
  bool
  value_is_less_than(const Number &value_1, const adouble &value_2)
  {
    // Use the specialized definition for two ADOL-C taped types
    return value_is_less_than(internal::NumberType<adouble>::value(value_1),
                              value_2);
  }

#endif


  template <typename Number1, typename Number2>
  constexpr bool
  values_are_equal(const Number1 &value_1, const Number2 &value_2)
  {
    return (value_1 == internal::NumberType<Number1>::value(value_2));
  }


  template <typename Number1, typename Number2>
  inline bool
  values_are_not_equal(const Number1 &value_1, const Number2 &value_2)
  {
    return !(values_are_equal(value_1, value_2));
  }


  template <typename Number>
  constexpr bool
  value_is_zero(const Number &value)
  {
    return values_are_equal(value, 0.0);
  }


  template <typename Number1, typename Number2>
  inline bool
  value_is_less_than(const Number1 &value_1, const Number2 &value_2)
  {
    return (value_1 < internal::NumberType<Number1>::value(value_2));
  }


  template <typename Number1, typename Number2>
  inline bool
  value_is_less_than_or_equal_to(const Number1 &value_1, const Number2 &value_2)
  {
    return (value_is_less_than(value_1, value_2) ||
            values_are_equal(value_1, value_2));
  }


  template <typename Number1, typename Number2>
  bool
  value_is_greater_than(const Number1 &value_1, const Number2 &value_2)
  {
    return !(value_is_less_than_or_equal_to(value_1, value_2));
  }


  template <typename Number1, typename Number2>
  inline bool
  value_is_greater_than_or_equal_to(const Number1 &value_1,
                                    const Number2 &value_2)
  {
    return !(value_is_less_than(value_1, value_2));
  }
} // namespace numbers

DEAL_II_NAMESPACE_CLOSE

#endif


