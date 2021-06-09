//include/deal.II-translator/base/vectorization_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2020 by the deal.II authors
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


#ifndef dealii_vectorization_h
#define dealii_vectorization_h

#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/template_constraints.h>

#include <array>
#include <cmath>

// Note:
// The flag DEAL_II_VECTORIZATION_WIDTH_IN_BITS is essentially constructed
// according to the following scheme (on x86-based architectures)
// #ifdef __AVX512F__
// #define DEAL_II_VECTORIZATION_WIDTH_IN_BITS 512
// #elif defined (__AVX__)
// #define DEAL_II_VECTORIZATION_WIDTH_IN_BITS 256
// #elif defined (__SSE2__)
// #define DEAL_II_VECTORIZATION_WIDTH_IN_BITS 128
// #else
// #define DEAL_II_VECTORIZATION_WIDTH_IN_BITS 0
// #endif
// In addition to checking the flags __AVX512F__, __AVX__ and __SSE2__, a CMake
// test, 'check_01_cpu_features.cmake', ensures that these feature are not only
// present in the compilation unit but also working properly.

#if DEAL_II_VECTORIZATION_WIDTH_IN_BITS > 0

// These error messages try to detect the case that deal.II was compiled with
// a wider instruction set extension as the current compilation unit, for
// example because deal.II was compiled with AVX, but a user project does not
// add -march=native or similar flags, making it fall to SSE2. This leads to
// very strange errors as the size of data structures differs between the
// compiled deal.II code sitting in libdeal_II.so and the user code if not
// detected.
#  if DEAL_II_VECTORIZATION_WIDTH_IN_BITS >= 256 && !defined(__AVX__)
#    error \
      "Mismatch in vectorization capabilities: AVX was detected during configuration of deal.II and switched on, but it is apparently not available for the file you are trying to compile at the moment. Check compilation flags controlling the instruction set, such as -march=native."
#  endif
#  if DEAL_II_VECTORIZATION_WIDTH_IN_BITS >= 512 && !defined(__AVX512F__)
#    error \
      "Mismatch in vectorization capabilities: AVX-512F was detected during configuration of deal.II and switched on, but it is apparently not available for the file you are trying to compile at the moment. Check compilation flags controlling the instruction set, such as -march=native."
#  endif

#  if defined(_MSC_VER)
#    include <intrin.h>
#  elif defined(__ALTIVEC__)
#    include <altivec.h>

// altivec.h defines vector, pixel, bool, but we do not use them, so undefine
// them before they make trouble
#    undef vector
#    undef pixel
#    undef bool
#  else
#    include <x86intrin.h>
#  endif

#endif


DEAL_II_NAMESPACE_OPEN


// Enable the EnableIfScalar type trait for VectorizedArray<Number> such
// that it can be used as a Number type in Tensor<rank,dim,Number>, etc.

template <typename Number, std::size_t width>
struct EnableIfScalar<VectorizedArray<Number, width>>
{
  using type = VectorizedArray<typename EnableIfScalar<Number>::type, width>;
};



/**
 * VectorizedArray的一个迭代器。
 *
 *
 */
template <typename T>
class VectorizedArrayIterator
{
public:
  /**
   * 构造函数。      @param  data 实际的VectorizedArray。    @param
   * lane 一个指向当前车道的指针。
   *
   */
  VectorizedArrayIterator(T &data, const std::size_t lane)
    : data(&data)
    , lane(lane)
  {}

  /**
   * 比较是否相等。
   *
   */
  bool
  operator==(const VectorizedArrayIterator<T> &other) const
  {
    Assert(this->data == other.data,
           ExcMessage(
             "You are trying to compare iterators into different arrays."));
    return this->lane == other.lane;
  }

  /**
   * 比较不等式。
   *
   */
  bool
  operator!=(const VectorizedArrayIterator<T> &other) const
  {
    Assert(this->data == other.data,
           ExcMessage(
             "You are trying to compare iterators into different arrays."));
    return this->lane != other.lane;
  }

  /**
   * 复制赋值。
   *
   */
  VectorizedArrayIterator<T> &
  operator=(const VectorizedArrayIterator<T> &other) = default;

  /**
   * 解除引用操作符（const版本）：返回当前车道的值。
   *
   */
  const typename T::value_type &operator*() const
  {
    AssertIndexRange(lane, T::size());
    return (*data)[lane];
  }


  /**
   * 去引用操作符（非 @p const 版本）：返回当前车道的值。
   *
   */
  template <typename U = T>
  typename std::enable_if<!std::is_same<U, const U>::value,
                          typename T::value_type>::type &
  operator*()
  {
    AssertIndexRange(lane, T::size());
    return (*data)[lane];
  }

  /**
   * 前缀<tt>++</tt>运算符。<tt>++iterator</tt>。这个操作符将迭代器推进到下一个车道，并返回对<tt>*this</tt>的引用。
   *
   */
  VectorizedArrayIterator<T> &
  operator++()
  {
    AssertIndexRange(lane + 1, T::size() + 1);
    lane++;
    return *this;
  }

  /**
   * 这个操作符使迭代器前进了 @p offset
   * 个通道，并返回对<tt>*this</tt>的引用。
   *
   */
  VectorizedArrayIterator<T> &
  operator+=(const std::size_t offset)
  {
    AssertIndexRange(lane + offset, T::size() + 1);
    lane += offset;
    return *this;
  }

  /**
   * 前缀<tt>--</tt>操作符。<tt>--iterator</tt>。这个操作符使迭代器前进到前一个通道，并返回一个对<tt>*this</tt>的引用。
   *
   */
  VectorizedArrayIterator<T> &
  operator--()
  {
    Assert(
      lane > 0,
      ExcMessage(
        "You can't decrement an iterator that is already at the beginning of the range."));
    --lane;
    return *this;
  }

  /**
   * 创建新的迭代器，该迭代器被移位 @p offset. 。
   *
   */
  VectorizedArrayIterator<T>
  operator+(const std::size_t &offset) const
  {
    AssertIndexRange(lane + offset, T::size() + 1);
    return VectorizedArrayIterator<T>(*data, lane + offset);
  }

  /**
   * 计算这个迭代器和迭代器之间的距离  @p other.  。
   *
   */
  std::ptrdiff_t
  operator-(const VectorizedArrayIterator<T> &other) const
  {
    return static_cast<std::ptrdiff_t>(lane) -
           static_cast<ptrdiff_t>(other.lane);
  }

private:
  /**
   * 指向实际的VectorizedArray的指针。
   *
   */
  T *data;

  /**
   * 指向当前车道的指针。
   *
   */
  std::size_t lane;
};



/**
 * 各种VectorizedArray模板特化的基类，包含共同的功能。
 * @tparam  T
 * 实际矢量数组的类型。我们在这个类中使用Couriously
 * Recurring Template
 * Pattern（见https://en.wikipedia.org/wiki/Curiously_recurring_template_pattern），以避免求助于`虚拟`成员函数。
 *
 *
 */
template <typename T, std::size_t width>
class VectorizedArrayBase
{
public:
  /**
   * 返回数组中元素的数量。
   *
   */
  static constexpr std::size_t
  size()
  {
    return width;
  }

  /**
   * @return  指向底层数据开始的一个迭代器。
   *
   */
  VectorizedArrayIterator<T>
  begin()
  {
    return VectorizedArrayIterator<T>(static_cast<T &>(*this), 0);
  }

  /**
   * @return  一个指向基础数据开始的迭代器（`const`版本）。
   *
   */
  VectorizedArrayIterator<const T>
  begin() const
  {
    return VectorizedArrayIterator<const T>(static_cast<const T &>(*this), 0);
  }

  /**
   * @return  指向基础数据结束的迭代器。
   *
   */
  VectorizedArrayIterator<T>
  end()
  {
    return VectorizedArrayIterator<T>(static_cast<T &>(*this), width);
  }

  /**
   * @return  一个指向基础数据末端的迭代器（`const`版本）。
   *
   */
  VectorizedArrayIterator<const T>
  end() const
  {
    return VectorizedArrayIterator<const T>(static_cast<const T &>(*this),
                                            width);
  }
};



/**
 * 这个通用类定义了一个矢量数据类型的统一接口。对于一般的模板参数，这个类简单地对应于模板参数。例如，VectorizedArray<long
 * double>只不过是<tt>long
 * double</tt>的一个封装器，它的数据字段正好是<tt>long
 * double</tt>类型，并且有重载的算术操作。这意味着<tt>VectorizedArray<ComplicatedType></tt>具有与ComplicatedType类似的布局，只要ComplicatedType定义了基本的算术操作。对于浮点数和双数来说，数组被打包在一起，目的是为了以单指令/多数据（SIMD）的方式进行处理。在SIMD背景下，这种短矢量的元素通常被称为通道。包装在一起的元素的数量，即车道的数量，取决于计算机系统和用于编译deal.II的编译器标志。这些打包数据类型的基本思想是使用一条CPU指令，利用处理器的向量（SIMD）单元对整个阵列进行算术运算。按照2010年的标准，大多数计算机系统在64位操作系统上编译deal.II时，将分别使用两个双数或四个浮点数的数组（这对应于SSE/SSE2数据集）。在英特尔Sandy
 * Bridge处理器和更新的处理器或AMD
 * Bulldozer处理器和更新的处理器上，当deal.II使用gcc配置为
 * \--with-cpu=native 或
 * \--with-cpu=corei7-avx时，会使用四个双数或八个浮点。在支持AVX-512的编译器上（例如2017年开始的英特尔Skylake服务器），会使用8个双数或16个浮点。
 * 这个类的行为被做成与基本数据类型double和float类似。向量数组的定义并不初始化数据字段，而是让它未被定义，就像double和float的情况一样。然而，当调用诸如`矢量化数组<double>
 * a =矢量化数组<double>()`或`矢量化数组<double> a =
 * 0.`时，它将这个字段中的所有数字设置为零。根据C++11标准，该类属于标准布局类型，这意味着有一个等效的C表示法，例如，该类可以用
 * std::memcpy.
 * 安全地复制（也见https://en.cppreference.com/w/cpp/named_req/StandardLayoutType）。标准布局对于确保在向量中收集的数据与地址边界正确对齐也是必要的（即，当向量中的第一个元素正确对齐时，所有后续元素也将正确对齐）。
 * 请注意，为了使这个类的正常运行，必须遵守某些数据对齐规则。这是因为计算机期望VectorizedArray<double>字段的起始地址在内存中的特定地址（通常，矢量数组的地址应该是以字节为单位的数组长度的倍数）。否则，可能会出现分段故障或严重的性能损失。当在堆栈上创建一个单一的数据字段时，如
 * "VectorizedArray<double> a = 5.
 * ;"，编译器会自动处理数据对齐。然而，当分配一个VectorizedArray<double>数据的长向量时，需要尊重这些规则。为此要使用AlignedVector类或基于AlignedVector的数据容器（如Table）。这是一个与
 * std::vector
 * 非常相似的类，否则，它总是能确保数据正确对齐。
 * 用户可以通过该封装类的第二个模板参数指定通道数，明确控制特定指令集架构（ISA）扩展的宽度。例如，在英特尔Skylake服务器上，你对数据类型double有以下选择。
 *
 *
 *
 *
 *
 * - VectorizedArray<double, 1> // 无矢量化（自动优化）。
 *
 *
 *
 *
 *
 * - VectorizedArray<double, 2> // SSE2
 *
 *
 *
 *
 *
 * - VectorizedArray<double, 4> // AVX
 *
 *
 *
 *
 * - VectorizedArray<double, 8> // AVX-512 (默认)
 * 并用于英特尔Sandy Bridge、Haswell、Broadwell、AMD
 * Bulldozer和Zen/Ryzen。
 *
 *
 *
 *
 *
 * - VectorizedArray<double, 1> // 无矢量化（自动优化）。
 *
 *
 *
 *
 *
 * - VectorizedArray<double, 2> // SSE2
 *
 *
 *
 *
 *
 * - VectorizedArray<double, 4> // AVX (默认)
 * 以及对于支持AltiVec的处理器。
 *
 *
 *
 *
 *
 * - 矢量Array<double, 1> * - 矢量Array<double, 1>
 *
 *
 * 适用于旧的x86处理器或在没有添加特定处理器编译标志的情况下（即没有`-D
 * CMAKE_CXX_FLAGS=-march=native`或类似标志）。
 *
 *
 *
 *
 *
 * - VectorizedArray<double, 1> // 无矢量化（自动优化）。
 *
 *
 *
 *
 *
 * - VectorizedArray<double, 2> // SSE2
 * 类似的考虑也适用于数据类型`float`。
 * 错误地选择宽度，例如，在不支持AVX-512的处理器上选择width=3或width=8，会导致静态断言。
 * @tparam  数字基础数据类型  @tparam
 * 宽度向量长度（可选；如果不设置，则使用架构的最大宽度）。
 *
 *
 */
template <typename Number, std::size_t width>
class VectorizedArray
  : public VectorizedArrayBase<VectorizedArray<Number, width>, 1>
{
public:
  /**
   * 这给出了数组元素的类型。
   *
   */
  using value_type = Number;

  static_assert(width == 1,
                "You specified an illegal width that is not supported.");

  /**
   * 默认的空构造函数，使数据处于未初始化的状态，类似于float/double。
   *
   */
  VectorizedArray() = default;

  /**
   * 用给定的标量构建一个数组广播到所有的通道。
   *
   */
  VectorizedArray(const Number scalar)
  {
    this->operator=(scalar);
  }

  /**
   * 这个函数将一个标量分配给这个类。
   *
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray &
  operator=(const Number scalar)
  {
    data = scalar;
    return *this;
  }

  /**
   * 访问操作符（只对基类中的0号组件有效，没有特殊化）。
   *
   */
  DEAL_II_ALWAYS_INLINE
  Number &operator[](const unsigned int comp)
  {
    (void)comp;
    AssertIndexRange(comp, 1);
    return data;
  }

  /**
   * 常数访问操作符（只对基类中的0号组件有效，没有特殊化）。
   *
   */
  DEAL_II_ALWAYS_INLINE
  const Number &operator[](const unsigned int comp) const
  {
    (void)comp;
    AssertIndexRange(comp, 1);
    return data;
  }

  /**
   * 加法
   *
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray &
  operator+=(const VectorizedArray &vec)
  {
    data += vec.data;
    return *this;
  }

  /**
   * 减法
   *
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray &
  operator-=(const VectorizedArray &vec)
  {
    data -= vec.data;
    return *this;
  }

  /**
   * 乘法
   *
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray &
  operator*=(const VectorizedArray &vec)
  {
    data *= vec.data;
    return *this;
  }

  /**
   * 除法
   *
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray &
  operator/=(const VectorizedArray &vec)
  {
    data /= vec.data;
    return *this;
  }

  /**
   * 从内存中加载size()数据项到调用类中，从给定的地址开始。指针`ptr`不需要按矢量数组中的字节数对齐，与之相反的是，将双倍地址投给VectorizedArray<double>*。
   *
   */
  DEAL_II_ALWAYS_INLINE
  void
  load(const Number *ptr)
  {
    data = *ptr;
  }

  /**
   * 将调用类的内容以size()数据项的形式写入内存中，到给定的地址。指针`ptr`不需要按矢量数组中的字节数对齐，与之相反的是，将一个双倍地址投给VectorizedArray<double>*。
   *
   */
  DEAL_II_ALWAYS_INLINE
  void
  store(Number *ptr) const
  {
    *ptr = data;
  }

  /**
   * 在支持的CPU上使用 @p _mm_stream_pd
   * 存储本征，将调用类的内容以size()数据项的形式写入内存，并使用绕过处理器缓存的非时间性存储来给定地址。存储
   * @p ptr 的目标必须按矢量数组中的字节数对齐。
   * 在存储是流式的情况下，这种存储操作可以比通常的存储操作更快，因为它避免了通常在标准存储中调用的读取所有权转移。这个大约的工作原理如下（详见计算机结构方面的文献）。当一个算法将一些结果存储到一个内存地址时，一个处理器通常希望将它移到它的一些缓存中，因为它期望这些数据在某个时候会被再次使用。由于缓存是以64字节或128字节的行来组织的，但写入的数据通常较小，所以处理器在写入时必须首先加载目标缓存行，因为最初只有部分缓存行被覆盖。如果一系列的存储写入的数据比任何一个缓冲区都要大，那么这些数据最后就必须从缓冲区转移到主内存。但是，由于所有的地址都是先被读取的，这就使主内存的负载增加了一倍，这就会产生性能上的损失。此外，在多核环境下的缓存组织也需要在将某一地址写入缓存之前读取该地址，详情请参见<a
   * href="https://en.wikipedia.org/wiki/MESI_protocol">Wikipedia article on
   * the MESI
   * protocol</a>。这个函数调用的基础指令向处理器发出信号，在存储上的这两个先决条件被放宽了：首先，人们期望整个高速缓存行被覆盖（意味着内存子系统确保一起跨越高速缓存行的连续存储被合并，并适当处理只有部分高速缓存行被写入的情况），所以没有必要首先读取高速缓存行的
   * "剩余部分"。其次，该特定内存后面的数据将不受缓存一致性协议的约束，因为当同一处理器想要再次访问它时，以及在多核芯片中的任何其他处理器，它都将在主内存中。由于这种特殊的设置，随后对该函数写入的数据的任何访问都需要查询主存，这在延迟和吞吐量上都比从高速缓存中访问要慢。因此，这个命令应该只用于存储大的数组，这些数组总体上不适合放在缓存中，否则性能就会下降。关于一个典型的用例，请参见<a
   * href="https://blogs.fau.de/hager/archives/2103">this blog article</a>。
   * 注意，流式存储只在类型为 @p double 或 @p float,
   * 的VectorizedArray的专业SSE/AVX类中可用，而不是在通用基类中。
   *
   */
  DEAL_II_ALWAYS_INLINE
  void
  streaming_store(Number *ptr) const
  {
    *ptr = data;
  }

  /**
   * 从内存中加载size()数据项到调用的类中，从给定的地址和给定的偏移量开始，从偏移量开始的每个条目提供一个矢量化数组的元素。
   * 这个操作对应于以下代码（但在硬件允许的情况下，使用了更有效的实现）。
   * @code
   * for (unsigned int v=0; v<VectorizedArray<Number>::size(); ++v)
   * this->operator[](v) = base_ptr[offsets[v]];
   * @endcode
   *
   *
   */
  DEAL_II_ALWAYS_INLINE
  void
  gather(const Number *base_ptr, const unsigned int *offsets)
  {
    data = base_ptr[offsets[0]];
  }

  /**
   * 将调用类的内容以size()数据项的形式写入内存，到给定的地址和给定的偏移量，将矢量数组的元素填入每个偏移量。
   * 这个操作对应于下面的代码（但在硬件允许的情况下，使用更有效的实现）。
   * @code
   * for (unsigned int v=0; v<VectorizedArray<Number>::size(); ++v)
   * base_ptr[offsets[v]] = this->operator[](v);
   * @endcode
   *
   */
  DEAL_II_ALWAYS_INLINE
  void
  scatter(const unsigned int *offsets, Number *base_ptr) const
  {
    base_ptr[offsets[0]] = data;
  }

  /**
   * 实际的数据字段。为了与标准布局类型保持一致，并且能够与外部SIMD功能进行交互，这个成员被声明为公共的。
   *
   */
  Number data;

private:
  /**
   * 返回这个字段的平方根。不适合在用户代码中使用。使用sqrt(x)代替。
   *
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray
  get_sqrt() const
  {
    VectorizedArray res;
    res.data = std::sqrt(data);
    return res;
  }

  /**
   * 返回这个字段的绝对值。不适合在用户代码中使用。请使用abs(x)代替。
   *
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray
  get_abs() const
  {
    VectorizedArray res;
    res.data = std::fabs(data);
    return res;
  }

  /**
   * 返回这个字段和另一个字段的分量上的最大值。不适合在用户代码中使用。使用max(x,y)代替。
   *
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray
  get_max(const VectorizedArray &other) const
  {
    VectorizedArray res;
    res.data = std::max(data, other.data);
    return res;
  }

  /**
   * 返回这个字段和另一个字段的最小分量。不适合在用户代码中使用。使用min(x,y)代替。
   *
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray
  get_min(const VectorizedArray &other) const
  {
    VectorizedArray res;
    res.data = std::min(data, other.data);
    return res;
  }

  // Make a few functions friends.
  template <typename Number2, std::size_t width2>
  friend VectorizedArray<Number2, width2>
  std::sqrt(const VectorizedArray<Number2, width2> &);
  template <typename Number2, std::size_t width2>
  friend VectorizedArray<Number2, width2>
  std::abs(const VectorizedArray<Number2, width2> &);
  template <typename Number2, std::size_t width2>
  friend VectorizedArray<Number2, width2>
  std::max(const VectorizedArray<Number2, width2> &,
           const VectorizedArray<Number2, width2> &);
  template <typename Number2, std::size_t width2>
  friend VectorizedArray<Number2, width2>
  std::min(const VectorizedArray<Number2, width2> &,
           const VectorizedArray<Number2, width2> &);
};



/**
 * @name  矢量数组的打包和拆包
 *
 *
 */
//@{


/**
 * 创建一个矢量数组，将数组中的所有条目设置为给定的标量，也就是说，将标量广播到所有数组元素。
 * @relatesalso VectorizedArray
 *
 *
 */
template <typename Number,
          std::size_t width =
            internal::VectorizedArrayWidthSpecifier<Number>::max_width>
inline DEAL_II_ALWAYS_INLINE VectorizedArray<Number, width>
                             make_vectorized_array(const Number &u)
{
  VectorizedArray<Number, width> result = u;
  return result;
}



/**
 * 创建一个给定类型的矢量数组，并将标量值广播给所有数组元素。
 * @relatesalso  VectorizedArray
 *
 *
 */
template <typename VectorizedArrayType>
inline DEAL_II_ALWAYS_INLINE VectorizedArrayType
                             make_vectorized_array(const typename VectorizedArrayType::value_type &u)
{
  static_assert(
    std::is_same<VectorizedArrayType,
                 VectorizedArray<typename VectorizedArrayType::value_type,
                                 VectorizedArrayType::size()>>::value,
    "VectorizedArrayType is not a VectorizedArray.");

  VectorizedArrayType result = u;
  return result;
}



/**
 * 从内存中加载size()数据项到VectorizedArray  @p out,
 * ，从给定的地址和给定的偏移量开始，从偏移量开始的每个条目提供一个矢量化数组的元素。
 * 这个操作对应于以下代码。
 *
 * @code
 * for (unsigned int v=0; v<VectorizedArray<Number>::size(); ++v)
 * out.data[v] = ptrs[v][offset];
 * @endcode
 *
 *
 *
 */
template <typename Number, std::size_t width>
inline DEAL_II_ALWAYS_INLINE void
gather(VectorizedArray<Number, width> &   out,
       const std::array<Number *, width> &ptrs,
       const unsigned int                 offset)
{
  for (unsigned int v = 0; v < width; v++)
    out.data[v] = ptrs[v][offset];
}



/**
 * 这个方法从给定的数组 @p in. 中加载 VectorizedArray::size()
 * 数据流，输入数组的偏移量由数组 @p
 * 偏移量给出。从每个数据流中读取n_entries。然后数据被转置并存储到一个VectorizedArray类型的数组中。输出数组
 * @p  out预计是一个大小为 @p n_entries.
 * 的数组。这个方法在普通数组上操作，所以没有检查有效的数据访问。用户有责任确保给定的数组根据下面的访问布局是有效的。
 * 该操作对应于根据以下公式将一个结构数组（输入）转换为一个数组结构（输出）。
 *
 *
 * @code
 * for (unsigned int i=0; i<n_entries; ++i)
 * for (unsigned int v=0; v<VectorizedArray<Number>::size(); ++v)
 *   out[i][v] = in[offsets[v]+i];
 * @endcode
 *
 * 该代码的一个更优化的版本将用于支持的类型。
 * 这是对vectorized_transpose_and_store()的逆向操作。
 * @relatesalso  VectorizedArray
 *
 *
 */
template <typename Number, std::size_t width>
inline DEAL_II_ALWAYS_INLINE void
vectorized_load_and_transpose(const unsigned int              n_entries,
                              const Number *                  in,
                              const unsigned int *            offsets,
                              VectorizedArray<Number, width> *out)
{
  for (unsigned int i = 0; i < n_entries; ++i)
    for (unsigned int v = 0; v < VectorizedArray<Number, width>::size(); ++v)
      out[i][v] = in[offsets[v] + i];
}


/**
 * 与上面的函数相同，不同的是将一个指针数组作为输入参数传入
 * @p in.  。 与上面的函数类比，可以认为
 * "in+offset[v]"是预先计算并作为输入参数传递的。
 * 然而，如果某些函数返回一个指针数组，并且不能假设它们属于同一个数组，也就是说，它们可以在不同的内存分配中获得原点，那么也可以使用这个函数。
 *
 *
 */
template <typename Number, std::size_t width>
inline DEAL_II_ALWAYS_INLINE void
vectorized_load_and_transpose(const unsigned int                 n_entries,
                              const std::array<Number *, width> &in,
                              VectorizedArray<Number, width> *   out)
{
  for (unsigned int i = 0; i < n_entries; ++i)
    for (unsigned int v = 0; v < VectorizedArray<Number, width>::size(); ++v)
      out[i][v] = in[v][i];
}



/**
 * 该方法以转置的形式将矢量数组存储到给定的输出数组 @p
 * out 中，并给定偏移量 @p offsets.
 * ，该操作相当于将一个数组结构（输入）转换为一个数组结构（输出）。该方法对纯数组进行操作，所以没有对有效的数据访问进行检查。用户有责任确保给定的数组根据下面的访问布局是有效的。
 * 该方法假设指定的偏移量不重叠。否则，在矢量化的情况下，该行为是未定义的。用户有责任确保访问不重叠，避免未定义的行为。
 * 参数 @p add_into
 * 选择哪里的条目应该只被写入输出数组，或者结果应该被添加到输出的现有条目中。对于
 * <code>add_into == false</code>  ，假设以下代码。
 *
 *
 * @code
 * for (unsigned int i=0; i<n_entries; ++i)
 * for (unsigned int v=0; v<VectorizedArray<Number>::size(); ++v)
 *   out[offsets[v]+i] = in[i][v];
 * @endcode
 *
 * 对于  <code>add_into == true</code>  ，代码实现了以下动作。
 * @code
 * for (unsigned int i=0; i<n_entries; ++i)
 * for (unsigned int v=0; v<VectorizedArray<Number>::size(); ++v)
 *   out[offsets[v]+i] += in[i][v];
 * @endcode
 *
 * 对于支持的类型，将使用该代码的一个更优化的版本。
 * 这是对vectorized_load_and_transpose()的逆向操作。
 * @relatesalso  VectorizedArray
 *
 *
 */
template <typename Number, std::size_t width>
inline DEAL_II_ALWAYS_INLINE void
vectorized_transpose_and_store(const bool                            add_into,
                               const unsigned int                    n_entries,
                               const VectorizedArray<Number, width> *in,
                               const unsigned int *                  offsets,
                               Number *                              out)
{
  if (add_into)
    for (unsigned int i = 0; i < n_entries; ++i)
      for (unsigned int v = 0; v < VectorizedArray<Number, width>::size(); ++v)
        out[offsets[v] + i] += in[i][v];
  else
    for (unsigned int i = 0; i < n_entries; ++i)
      for (unsigned int v = 0; v < VectorizedArray<Number, width>::size(); ++v)
        out[offsets[v] + i] = in[i][v];
}


/**
 * 和上面一样，不同的是，一个指针数组被作为输入参数传入
 * @p out.  。
 * 与上面的函数相类似，可以认为`out+offset[v]`是预先计算好的，并作为输入参数传入。
 * 然而，如果某些函数返回一个指针数组，并且不能假设它们属于同一个数组，也就是说，它们可以在不同的内存分配中拥有自己的原点，那么也可以使用这个函数。
 *
 *
 */
template <typename Number, std::size_t width>
inline DEAL_II_ALWAYS_INLINE void
vectorized_transpose_and_store(const bool                            add_into,
                               const unsigned int                    n_entries,
                               const VectorizedArray<Number, width> *in,
                               std::array<Number *, width> &         out)
{
  if (add_into)
    for (unsigned int i = 0; i < n_entries; ++i)
      for (unsigned int v = 0; v < VectorizedArray<Number, width>::size(); ++v)
        out[v][i] += in[i][v];
  else
    for (unsigned int i = 0; i < n_entries; ++i)
      for (unsigned int v = 0; v < VectorizedArray<Number, width>::size(); ++v)
        out[v][i] = in[i][v];
}


//@}

#ifndef DOXYGEN

// for safety, also check that __AVX512F__ is defined in case the user manually
// set some conflicting compile flags which prevent compilation

#  if DEAL_II_VECTORIZATION_WIDTH_IN_BITS >= 512 && defined(__AVX512F__)

/**
 * 针对double和AVX-512的VectorizedArray类的专业化。
 *
 *
 */
template <>
class VectorizedArray<double, 8>
  : public VectorizedArrayBase<VectorizedArray<double, 8>, 8>
{
public:
  /**
   * 这给出了数组元素的类型。
   *
   */
  using value_type = double;

  /**
   * 默认的空构造函数，让数据处于未初始化的状态，类似于float/double。
   *
   */
  VectorizedArray() = default;

  /**
   * 用给定的标量构建一个数组，广播给所有通道。
   *
   */
  VectorizedArray(const double scalar)
  {
    this->operator=(scalar);
  }

  /**
   * 这个函数可以用来将所有的数据字段设置为一个给定的标量。
   *
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray &
  operator=(const double x)
  {
    data = _mm512_set1_pd(x);
    return *this;
  }

  /**
   * 访问操作符。
   *
   */
  DEAL_II_ALWAYS_INLINE
  double &operator[](const unsigned int comp)
  {
    AssertIndexRange(comp, 8);
    return *(reinterpret_cast<double *>(&data) + comp);
  }

  /**
   * 常数访问运算符。
   *
   */
  DEAL_II_ALWAYS_INLINE
  const double &operator[](const unsigned int comp) const
  {
    AssertIndexRange(comp, 8);
    return *(reinterpret_cast<const double *>(&data) + comp);
  }

  /**
   * 加法。
   *
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray &
  operator+=(const VectorizedArray &vec)
  {
    // if the compiler supports vector arithmetic, we can simply use +=
    // operator on the given data type. this allows the compiler to combine
    // additions with multiplication (fused multiply-add) if those
    // instructions are available. Otherwise, we need to use the built-in
    // intrinsic command for __m512d
#    ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data += vec.data;
#    else
    data = _mm512_add_pd(data, vec.data);
#    endif
    return *this;
  }

  /**
   * 减法。
   *
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray &
  operator-=(const VectorizedArray &vec)
  {
#    ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data -= vec.data;
#    else
    data = _mm512_sub_pd(data, vec.data);
#    endif
    return *this;
  }
  /**
   * 乘法。
   *
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray &
  operator*=(const VectorizedArray &vec)
  {
#    ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data *= vec.data;
#    else
    data = _mm512_mul_pd(data, vec.data);
#    endif
    return *this;
  }

  /**
   * 除法。
   *
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray &
  operator/=(const VectorizedArray &vec)
  {
#    ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data /= vec.data;
#    else
    data = _mm512_div_pd(data, vec.data);
#    endif
    return *this;
  }

  /**
   * 从内存中加载size()数据项到调用类中，从给定的地址开始。内存不需要按64字节对齐，相对于将双倍地址投给VectorizedArray<double>*来说。
   *
   */
  DEAL_II_ALWAYS_INLINE
  void
  load(const double *ptr)
  {
    data = _mm512_loadu_pd(ptr);
  }

  /**
   * 将调用类的内容以 @p
   * size()的形式写入内存到给定地址。内存不需要按64字节对齐，相对于将双倍地址投给VectorizedArray<double>*来说。
   *
   */
  DEAL_II_ALWAYS_INLINE
  void
  store(double *ptr) const
  {
    _mm512_storeu_pd(ptr, data);
  }

  /**
   * @copydoc   VectorizedArray<Number>::streaming_store() 。
   * @note  内存必须以64字节对齐。
   *
   */
  DEAL_II_ALWAYS_INLINE
  void
  streaming_store(double *ptr) const
  {
    Assert(reinterpret_cast<std::size_t>(ptr) % 64 == 0,
           ExcMessage("Memory not aligned"));
    _mm512_stream_pd(ptr, data);
  }

  /**
   * 将 @p size()
   * 从内存中加载到调用类中，从给定的地址开始，并有给定的偏移量，从偏移量开始的每个条目提供一个矢量数组的元素。
   * 这个操作对应于以下代码（但在硬件允许的情况下，使用了更有效的实现）。
   * @code
   * for (unsigned int v=0; v<VectorizedArray<Number>::size(); ++v)
   * this->operator[](v) = base_ptr[offsets[v]];
   * @endcode
   *
   *
   */
  DEAL_II_ALWAYS_INLINE
  void
  gather(const double *base_ptr, const unsigned int *offsets)
  {
    // unfortunately, there does not appear to be a 256 bit integer load, so
    // do it by some reinterpret casts here. this is allowed because the Intel
    // API allows aliasing between different vector types.
    const __m256 index_val =
      _mm256_loadu_ps(reinterpret_cast<const float *>(offsets));
    const __m256i index = *reinterpret_cast<const __m256i *>(&index_val);
    data                = _mm512_i32gather_pd(index, base_ptr, 8);
  }

  /**
   * 将调用类的内容以 @p
   * size()的形式写入内存，到给定的地址和给定的偏移量，将矢量数组的元素填入每个偏移量。
   * 这个操作对应于下面的代码（但在硬件允许的情况下使用了更有效的实现）。
   * @code
   * for (unsigned int v=0; v<VectorizedArray<Number>::size(); ++v)
   * base_ptr[offsets[v]] = this->operator[](v);
   * @endcode
   *
   */
  DEAL_II_ALWAYS_INLINE
  void
  scatter(const unsigned int *offsets, double *base_ptr) const
  {
    for (unsigned int i = 0; i < 8; ++i)
      for (unsigned int j = i + 1; j < 8; ++j)
        Assert(offsets[i] != offsets[j],
               ExcMessage("Result of scatter undefined if two offset elements"
                          " point to the same position"));

    // unfortunately, there does not appear to be a 256 bit integer load, so
    // do it by some reinterpret casts here. this is allowed because the Intel
    // API allows aliasing between different vector types.
    const __m256 index_val =
      _mm256_loadu_ps(reinterpret_cast<const float *>(offsets));
    const __m256i index = *reinterpret_cast<const __m256i *>(&index_val);
    _mm512_i32scatter_pd(base_ptr, index, data, 8);
  }

  /**
   * 实际的数据字段。为了与标准布局类型保持一致，并且能够与外部SIMD功能进行交互，这个成员被声明为公共的。
   *
   */
  __m512d data;

private:
  /**
   * 返回这个字段的平方根。不适合在用户代码中使用。使用sqrt(x)代替。
   *
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray
  get_sqrt() const
  {
    VectorizedArray res;
    res.data = _mm512_sqrt_pd(data);
    return res;
  }

  /**
   * 返回该字段的绝对值。不适合在用户代码中使用。请使用abs(x)代替。
   *
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray
  get_abs() const
  {
    // to compute the absolute value, perform bitwise andnot with -0. This
    // will leave all value and exponent bits unchanged but force the sign
    // value to +. Since there is no andnot for AVX512, we interpret the data
    // as 64 bit integers and do the andnot on those types (note that andnot
    // is a bitwise operation so the data type does not matter)
    __m512d         mask = _mm512_set1_pd(-0.);
    VectorizedArray res;
    res.data = reinterpret_cast<__m512d>(
      _mm512_andnot_epi64(reinterpret_cast<__m512i>(mask),
                          reinterpret_cast<__m512i>(data)));
    return res;
  }

  /**
   * 返回这个字段和另一个字段的分量上的最大值。不适合在用户代码中使用。使用max(x,y)代替。
   *
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray
  get_max(const VectorizedArray &other) const
  {
    VectorizedArray res;
    res.data = _mm512_max_pd(data, other.data);
    return res;
  }

  /**
   * 返回这个字段和另一个字段的最小分量。不适合在用户代码中使用。使用min(x,y)代替。
   *
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray
  get_min(const VectorizedArray &other) const
  {
    VectorizedArray res;
    res.data = _mm512_min_pd(data, other.data);
    return res;
  }

  // Make a few functions friends.
  template <typename Number2, std::size_t width2>
  friend VectorizedArray<Number2, width2>
  std::sqrt(const VectorizedArray<Number2, width2> &);
  template <typename Number2, std::size_t width2>
  friend VectorizedArray<Number2, width2>
  std::abs(const VectorizedArray<Number2, width2> &);
  template <typename Number2, std::size_t width2>
  friend VectorizedArray<Number2, width2>
  std::max(const VectorizedArray<Number2, width2> &,
           const VectorizedArray<Number2, width2> &);
  template <typename Number2, std::size_t width2>
  friend VectorizedArray<Number2, width2>
  std::min(const VectorizedArray<Number2, width2> &,
           const VectorizedArray<Number2, width2> &);
};



/**
 * 对double和AVX-512的特殊化。
 *
 *
 */
template <>
inline DEAL_II_ALWAYS_INLINE void
vectorized_load_and_transpose(const unsigned int          n_entries,
                              const double *              in,
                              const unsigned int *        offsets,
                              VectorizedArray<double, 8> *out)
{
  // do not do full transpose because the code is long and will most
  // likely not pay off because many processors have two load units
  // (for the top 8 instructions) but only 1 permute unit (for the 8
  // shuffle/unpack instructions). rather start the transposition on the
  // vectorized array of half the size with 256 bits
  const unsigned int n_chunks = n_entries / 4;
  for (unsigned int i = 0; i < n_chunks; ++i)
    {
      __m512d t0, t1, t2, t3 = {};

      t0 = _mm512_insertf64x4(t3, _mm256_loadu_pd(in + offsets[0] + 4 * i), 0);
      t0 = _mm512_insertf64x4(t0, _mm256_loadu_pd(in + offsets[2] + 4 * i), 1);
      t1 = _mm512_insertf64x4(t3, _mm256_loadu_pd(in + offsets[1] + 4 * i), 0);
      t1 = _mm512_insertf64x4(t1, _mm256_loadu_pd(in + offsets[3] + 4 * i), 1);
      t2 = _mm512_insertf64x4(t3, _mm256_loadu_pd(in + offsets[4] + 4 * i), 0);
      t2 = _mm512_insertf64x4(t2, _mm256_loadu_pd(in + offsets[6] + 4 * i), 1);
      t3 = _mm512_insertf64x4(t3, _mm256_loadu_pd(in + offsets[5] + 4 * i), 0);
      t3 = _mm512_insertf64x4(t3, _mm256_loadu_pd(in + offsets[7] + 4 * i), 1);

      __m512d v0          = _mm512_shuffle_f64x2(t0, t2, 0x88);
      __m512d v1          = _mm512_shuffle_f64x2(t0, t2, 0xdd);
      __m512d v2          = _mm512_shuffle_f64x2(t1, t3, 0x88);
      __m512d v3          = _mm512_shuffle_f64x2(t1, t3, 0xdd);
      out[4 * i + 0].data = _mm512_unpacklo_pd(v0, v2);
      out[4 * i + 1].data = _mm512_unpackhi_pd(v0, v2);
      out[4 * i + 2].data = _mm512_unpacklo_pd(v1, v3);
      out[4 * i + 3].data = _mm512_unpackhi_pd(v1, v3);
    }
  // remainder loop of work that does not divide by 4
  for (unsigned int i = 4 * n_chunks; i < n_entries; ++i)
    out[i].gather(in + i, offsets);
}



/**
 * 对double和AVX-512的特殊化。
 *
 *
 */
template <>
inline DEAL_II_ALWAYS_INLINE void
vectorized_load_and_transpose(const unsigned int             n_entries,
                              const std::array<double *, 8> &in,
                              VectorizedArray<double, 8> *   out)
{
  const unsigned int n_chunks = n_entries / 4;
  for (unsigned int i = 0; i < n_chunks; ++i)
    {
      __m512d t0, t1, t2, t3 = {};

      t0 = _mm512_insertf64x4(t3, _mm256_loadu_pd(in[0] + 4 * i), 0);
      t0 = _mm512_insertf64x4(t0, _mm256_loadu_pd(in[2] + 4 * i), 1);
      t1 = _mm512_insertf64x4(t3, _mm256_loadu_pd(in[1] + 4 * i), 0);
      t1 = _mm512_insertf64x4(t1, _mm256_loadu_pd(in[3] + 4 * i), 1);
      t2 = _mm512_insertf64x4(t3, _mm256_loadu_pd(in[4] + 4 * i), 0);
      t2 = _mm512_insertf64x4(t2, _mm256_loadu_pd(in[6] + 4 * i), 1);
      t3 = _mm512_insertf64x4(t3, _mm256_loadu_pd(in[5] + 4 * i), 0);
      t3 = _mm512_insertf64x4(t3, _mm256_loadu_pd(in[7] + 4 * i), 1);

      __m512d v0          = _mm512_shuffle_f64x2(t0, t2, 0x88);
      __m512d v1          = _mm512_shuffle_f64x2(t0, t2, 0xdd);
      __m512d v2          = _mm512_shuffle_f64x2(t1, t3, 0x88);
      __m512d v3          = _mm512_shuffle_f64x2(t1, t3, 0xdd);
      out[4 * i + 0].data = _mm512_unpacklo_pd(v0, v2);
      out[4 * i + 1].data = _mm512_unpackhi_pd(v0, v2);
      out[4 * i + 2].data = _mm512_unpacklo_pd(v1, v3);
      out[4 * i + 3].data = _mm512_unpackhi_pd(v1, v3);
    }

  for (unsigned int i = 4 * n_chunks; i < n_entries; ++i)
    gather(out[i], in, i);
}



/**
 * 用于双倍和AVX-512的特殊化。
 *
 *
 */
template <>
inline DEAL_II_ALWAYS_INLINE void
vectorized_transpose_and_store(const bool                        add_into,
                               const unsigned int                n_entries,
                               const VectorizedArray<double, 8> *in,
                               const unsigned int *              offsets,
                               double *                          out)
{
  // as for the load, we split the store operations into 256 bit units to
  // better balance between code size, shuffle instructions, and stores
  const unsigned int n_chunks = n_entries / 4;
  __m512i mask1 = _mm512_set_epi64(0xd, 0xc, 0x5, 0x4, 0x9, 0x8, 0x1, 0x0);
  __m512i mask2 = _mm512_set_epi64(0xf, 0xe, 0x7, 0x6, 0xb, 0xa, 0x3, 0x2);
  for (unsigned int i = 0; i < n_chunks; ++i)
    {
      __m512d t0   = _mm512_unpacklo_pd(in[i * 4].data, in[i * 4 + 1].data);
      __m512d t1   = _mm512_unpackhi_pd(in[i * 4].data, in[i * 4 + 1].data);
      __m512d t2   = _mm512_unpacklo_pd(in[i * 4 + 2].data, in[i * 4 + 3].data);
      __m512d t3   = _mm512_unpackhi_pd(in[i * 4 + 2].data, in[i * 4 + 3].data);
      __m512d v0   = _mm512_permutex2var_pd(t0, mask1, t2);
      __m512d v1   = _mm512_permutex2var_pd(t0, mask2, t2);
      __m512d v2   = _mm512_permutex2var_pd(t1, mask1, t3);
      __m512d v3   = _mm512_permutex2var_pd(t1, mask2, t3);
      __m256d res0 = _mm512_extractf64x4_pd(v0, 0);
      __m256d res4 = _mm512_extractf64x4_pd(v0, 1);
      __m256d res1 = _mm512_extractf64x4_pd(v2, 0);
      __m256d res5 = _mm512_extractf64x4_pd(v2, 1);
      __m256d res2 = _mm512_extractf64x4_pd(v1, 0);
      __m256d res6 = _mm512_extractf64x4_pd(v1, 1);
      __m256d res3 = _mm512_extractf64x4_pd(v3, 0);
      __m256d res7 = _mm512_extractf64x4_pd(v3, 1);

      // Cannot use the same store instructions in both paths of the 'if'
      // because the compiler cannot know that there is no aliasing
      // between pointers
      if (add_into)
        {
          res0 = _mm256_add_pd(_mm256_loadu_pd(out + 4 * i + offsets[0]), res0);
          _mm256_storeu_pd(out + 4 * i + offsets[0], res0);
          res1 = _mm256_add_pd(_mm256_loadu_pd(out + 4 * i + offsets[1]), res1);
          _mm256_storeu_pd(out + 4 * i + offsets[1], res1);
          res2 = _mm256_add_pd(_mm256_loadu_pd(out + 4 * i + offsets[2]), res2);
          _mm256_storeu_pd(out + 4 * i + offsets[2], res2);
          res3 = _mm256_add_pd(_mm256_loadu_pd(out + 4 * i + offsets[3]), res3);
          _mm256_storeu_pd(out + 4 * i + offsets[3], res3);
          res4 = _mm256_add_pd(_mm256_loadu_pd(out + 4 * i + offsets[4]), res4);
          _mm256_storeu_pd(out + 4 * i + offsets[4], res4);
          res5 = _mm256_add_pd(_mm256_loadu_pd(out + 4 * i + offsets[5]), res5);
          _mm256_storeu_pd(out + 4 * i + offsets[5], res5);
          res6 = _mm256_add_pd(_mm256_loadu_pd(out + 4 * i + offsets[6]), res6);
          _mm256_storeu_pd(out + 4 * i + offsets[6], res6);
          res7 = _mm256_add_pd(_mm256_loadu_pd(out + 4 * i + offsets[7]), res7);
          _mm256_storeu_pd(out + 4 * i + offsets[7], res7);
        }
      else
        {
          _mm256_storeu_pd(out + 4 * i + offsets[0], res0);
          _mm256_storeu_pd(out + 4 * i + offsets[1], res1);
          _mm256_storeu_pd(out + 4 * i + offsets[2], res2);
          _mm256_storeu_pd(out + 4 * i + offsets[3], res3);
          _mm256_storeu_pd(out + 4 * i + offsets[4], res4);
          _mm256_storeu_pd(out + 4 * i + offsets[5], res5);
          _mm256_storeu_pd(out + 4 * i + offsets[6], res6);
          _mm256_storeu_pd(out + 4 * i + offsets[7], res7);
        }
    }

  // remainder loop of work that does not divide by 4
  if (add_into)
    for (unsigned int i = 4 * n_chunks; i < n_entries; ++i)
      for (unsigned int v = 0; v < 8; ++v)
        out[offsets[v] + i] += in[i][v];
  else
    for (unsigned int i = 4 * n_chunks; i < n_entries; ++i)
      for (unsigned int v = 0; v < 8; ++v)
        out[offsets[v] + i] = in[i][v];
}



/**
 * 用于双倍和AVX-512的特殊化。
 *
 *
 */
template <>
inline DEAL_II_ALWAYS_INLINE void
vectorized_transpose_and_store(const bool                        add_into,
                               const unsigned int                n_entries,
                               const VectorizedArray<double, 8> *in,
                               std::array<double *, 8> &         out)
{
  // see the comments in the vectorized_transpose_and_store above

  const unsigned int n_chunks = n_entries / 4;
  __m512i mask1 = _mm512_set_epi64(0xd, 0xc, 0x5, 0x4, 0x9, 0x8, 0x1, 0x0);
  __m512i mask2 = _mm512_set_epi64(0xf, 0xe, 0x7, 0x6, 0xb, 0xa, 0x3, 0x2);
  for (unsigned int i = 0; i < n_chunks; ++i)
    {
      __m512d t0   = _mm512_unpacklo_pd(in[i * 4].data, in[i * 4 + 1].data);
      __m512d t1   = _mm512_unpackhi_pd(in[i * 4].data, in[i * 4 + 1].data);
      __m512d t2   = _mm512_unpacklo_pd(in[i * 4 + 2].data, in[i * 4 + 3].data);
      __m512d t3   = _mm512_unpackhi_pd(in[i * 4 + 2].data, in[i * 4 + 3].data);
      __m512d v0   = _mm512_permutex2var_pd(t0, mask1, t2);
      __m512d v1   = _mm512_permutex2var_pd(t0, mask2, t2);
      __m512d v2   = _mm512_permutex2var_pd(t1, mask1, t3);
      __m512d v3   = _mm512_permutex2var_pd(t1, mask2, t3);
      __m256d res0 = _mm512_extractf64x4_pd(v0, 0);
      __m256d res4 = _mm512_extractf64x4_pd(v0, 1);
      __m256d res1 = _mm512_extractf64x4_pd(v2, 0);
      __m256d res5 = _mm512_extractf64x4_pd(v2, 1);
      __m256d res2 = _mm512_extractf64x4_pd(v1, 0);
      __m256d res6 = _mm512_extractf64x4_pd(v1, 1);
      __m256d res3 = _mm512_extractf64x4_pd(v3, 0);
      __m256d res7 = _mm512_extractf64x4_pd(v3, 1);

      if (add_into)
        {
          res0 = _mm256_add_pd(_mm256_loadu_pd(out[0] + 4 * i), res0);
          _mm256_storeu_pd(out[0] + 4 * i, res0);
          res1 = _mm256_add_pd(_mm256_loadu_pd(out[1] + 4 * i), res1);
          _mm256_storeu_pd(out[1] + 4 * i, res1);
          res2 = _mm256_add_pd(_mm256_loadu_pd(out[2] + 4 * i), res2);
          _mm256_storeu_pd(out[2] + 4 * i, res2);
          res3 = _mm256_add_pd(_mm256_loadu_pd(out[3] + 4 * i), res3);
          _mm256_storeu_pd(out[3] + 4 * i, res3);
          res4 = _mm256_add_pd(_mm256_loadu_pd(out[4] + 4 * i), res4);
          _mm256_storeu_pd(out[4] + 4 * i, res4);
          res5 = _mm256_add_pd(_mm256_loadu_pd(out[5] + 4 * i), res5);
          _mm256_storeu_pd(out[5] + 4 * i, res5);
          res6 = _mm256_add_pd(_mm256_loadu_pd(out[6] + 4 * i), res6);
          _mm256_storeu_pd(out[6] + 4 * i, res6);
          res7 = _mm256_add_pd(_mm256_loadu_pd(out[7] + 4 * i), res7);
          _mm256_storeu_pd(out[7] + 4 * i, res7);
        }
      else
        {
          _mm256_storeu_pd(out[0] + 4 * i, res0);
          _mm256_storeu_pd(out[1] + 4 * i, res1);
          _mm256_storeu_pd(out[2] + 4 * i, res2);
          _mm256_storeu_pd(out[3] + 4 * i, res3);
          _mm256_storeu_pd(out[4] + 4 * i, res4);
          _mm256_storeu_pd(out[5] + 4 * i, res5);
          _mm256_storeu_pd(out[6] + 4 * i, res6);
          _mm256_storeu_pd(out[7] + 4 * i, res7);
        }
    }

  if (add_into)
    for (unsigned int i = 4 * n_chunks; i < n_entries; ++i)
      for (unsigned int v = 0; v < 8; ++v)
        out[v][i] += in[i][v];
  else
    for (unsigned int i = 4 * n_chunks; i < n_entries; ++i)
      for (unsigned int v = 0; v < 8; ++v)
        out[v][i] = in[i][v];
}



/**
 * 浮动和AVX512的特殊化。
 *
 *
 */
template <>
class VectorizedArray<float, 16>
  : public VectorizedArrayBase<VectorizedArray<float, 16>, 16>
{
public:
  /**
   * 这给出了数组元素的类型。
   *
   */
  using value_type = float;

  /**
   * 默认的空构造函数，让数据处于未初始化的状态，类似于float/double。
   *
   */
  VectorizedArray() = default;

  /**
   * 用给定的标量构建一个数组广播到所有的通道。
   *
   */
  VectorizedArray(const float scalar)
  {
    this->operator=(scalar);
  }

  /**
   * 这个函数可以用来将所有的数据字段设置为一个给定的标量。
   *
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray &
  operator=(const float x)
  {
    data = _mm512_set1_ps(x);
    return *this;
  }

  /**
   * 访问操作符。
   *
   */
  DEAL_II_ALWAYS_INLINE
  float &operator[](const unsigned int comp)
  {
    AssertIndexRange(comp, 16);
    return *(reinterpret_cast<float *>(&data) + comp);
  }

  /**
   * 常数访问运算符。
   *
   */
  DEAL_II_ALWAYS_INLINE
  const float &operator[](const unsigned int comp) const
  {
    AssertIndexRange(comp, 16);
    return *(reinterpret_cast<const float *>(&data) + comp);
  }

  /**
   * 加法。
   *
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray &
  operator+=(const VectorizedArray &vec)
  {
    // if the compiler supports vector arithmetic, we can simply use +=
    // operator on the given data type. this allows the compiler to combine
    // additions with multiplication (fused multiply-add) if those
    // instructions are available. Otherwise, we need to use the built-in
    // intrinsic command for __m512d
#    ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data += vec.data;
#    else
    data = _mm512_add_ps(data, vec.data);
#    endif
    return *this;
  }

  /**
   * 减法。
   *
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray &
  operator-=(const VectorizedArray &vec)
  {
#    ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data -= vec.data;
#    else
    data = _mm512_sub_ps(data, vec.data);
#    endif
    return *this;
  }
  /**
   * 乘法。
   *
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray &
  operator*=(const VectorizedArray &vec)
  {
#    ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data *= vec.data;
#    else
    data = _mm512_mul_ps(data, vec.data);
#    endif
    return *this;
  }

  /**
   * 除法。
   *
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray &
  operator/=(const VectorizedArray &vec)
  {
#    ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data /= vec.data;
#    else
    data = _mm512_div_ps(data, vec.data);
#    endif
    return *this;
  }

  /**
   * 从内存中加载 @p size()
   * 到调用类中，从给定的地址开始。内存不需要对齐64字节，相对于将浮点地址投给VectorizedArray<float>*。
   *
   */
  DEAL_II_ALWAYS_INLINE
  void
  load(const float *ptr)
  {
    data = _mm512_loadu_ps(ptr);
  }

  /**
   * 将调用类的内容以 @p
   * size()的形式写到内存中的给定地址。内存不需要按64字节对齐，相对于将浮点地址投给VectorizedArray<float>*来说。
   *
   */
  DEAL_II_ALWAYS_INLINE
  void
  store(float *ptr) const
  {
    _mm512_storeu_ps(ptr, data);
  }

  /**
   * @copydoc   VectorizedArray<Number>::streaming_store() 。
   * @note  内存必须以64字节对齐。
   *
   */
  DEAL_II_ALWAYS_INLINE
  void
  streaming_store(float *ptr) const
  {
    Assert(reinterpret_cast<std::size_t>(ptr) % 64 == 0,
           ExcMessage("Memory not aligned"));
    _mm512_stream_ps(ptr, data);
  }

  /**
   * 将 @p size()
   * 从内存中加载到调用类中，从给定的地址开始，并有给定的偏移量，从偏移量开始的每个条目提供一个矢量数组的元素。
   * 这个操作对应于以下代码（但在硬件允许的情况下，使用了更有效的实现）。
   * @code
   * for (unsigned int v=0; v<VectorizedArray<Number>::size(); ++v)
   * this->operator[](v) = base_ptr[offsets[v]];
   * @endcode
   *
   *
   */
  DEAL_II_ALWAYS_INLINE
  void
  gather(const float *base_ptr, const unsigned int *offsets)
  {
    // unfortunately, there does not appear to be a 512 bit integer load, so
    // do it by some reinterpret casts here. this is allowed because the Intel
    // API allows aliasing between different vector types.
    const __m512 index_val =
      _mm512_loadu_ps(reinterpret_cast<const float *>(offsets));
    const __m512i index = *reinterpret_cast<const __m512i *>(&index_val);
    data                = _mm512_i32gather_ps(index, base_ptr, 4);
  }

  /**
   * 将调用类的内容以 @p
   * size()的形式写入内存，到给定的地址和给定的偏移量，将矢量数组的元素填入每个偏移量。
   * 这个操作对应于下面的代码（但在硬件允许的情况下使用了更有效的实现）。
   * @code
   * for (unsigned int v=0; v<VectorizedArray<Number>::size(); ++v)
   * base_ptr[offsets[v]] = this->operator[](v);
   * @endcode
   *
   *
   */
  DEAL_II_ALWAYS_INLINE
  void
  scatter(const unsigned int *offsets, float *base_ptr) const
  {
    for (unsigned int i = 0; i < 16; ++i)
      for (unsigned int j = i + 1; j < 16; ++j)
        Assert(offsets[i] != offsets[j],
               ExcMessage("Result of scatter undefined if two offset elements"
                          " point to the same position"));

    // unfortunately, there does not appear to be a 512 bit integer load, so
    // do it by some reinterpret casts here. this is allowed because the Intel
    // API allows aliasing between different vector types.
    const __m512 index_val =
      _mm512_loadu_ps(reinterpret_cast<const float *>(offsets));
    const __m512i index = *reinterpret_cast<const __m512i *>(&index_val);
    _mm512_i32scatter_ps(base_ptr, index, data, 4);
  }

  /**
   * 实际的数据字段。为了与标准布局类型保持一致，并且能够与外部SIMD功能进行交互，这个成员被声明为公共的。
   *
   */
  __m512 data;

private:
  /**
   * 返回这个字段的平方根。不适合在用户代码中使用。使用sqrt(x)代替。
   *
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray
  get_sqrt() const
  {
    VectorizedArray res;
    res.data = _mm512_sqrt_ps(data);
    return res;
  }

  /**
   * 返回该字段的绝对值。不适合在用户代码中使用。使用abs(x)代替。
   *
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray
  get_abs() const
  {
    // to compute the absolute value, perform bitwise andnot with -0. This
    // will leave all value and exponent bits unchanged but force the sign
    // value to +. Since there is no andnot for AVX512, we interpret the data
    // as 32 bit integers and do the andnot on those types (note that andnot
    // is a bitwise operation so the data type does not matter)
    __m512          mask = _mm512_set1_ps(-0.f);
    VectorizedArray res;
    res.data = reinterpret_cast<__m512>(
      _mm512_andnot_epi32(reinterpret_cast<__m512i>(mask),
                          reinterpret_cast<__m512i>(data)));
    return res;
  }

  /**
   * 返回这个字段和另一个字段的分量上的最大值。不适合在用户代码中使用。使用max(x,y)代替。
   *
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray
  get_max(const VectorizedArray &other) const
  {
    VectorizedArray res;
    res.data = _mm512_max_ps(data, other.data);
    return res;
  }

  /**
   * 返回这个字段和另一个字段的最小分量。不适合在用户代码中使用。使用min(x,y)代替。
   *
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray
  get_min(const VectorizedArray &other) const
  {
    VectorizedArray res;
    res.data = _mm512_min_ps(data, other.data);
    return res;
  }

  // Make a few functions friends.
  template <typename Number2, std::size_t width2>
  friend VectorizedArray<Number2, width2>
  std::sqrt(const VectorizedArray<Number2, width2> &);
  template <typename Number2, std::size_t width2>
  friend VectorizedArray<Number2, width2>
  std::abs(const VectorizedArray<Number2, width2> &);
  template <typename Number2, std::size_t width2>
  friend VectorizedArray<Number2, width2>
  std::max(const VectorizedArray<Number2, width2> &,
           const VectorizedArray<Number2, width2> &);
  template <typename Number2, std::size_t width2>
  friend VectorizedArray<Number2, width2>
  std::min(const VectorizedArray<Number2, width2> &,
           const VectorizedArray<Number2, width2> &);
};



/**
 * 浮点数和AVX-512的特殊化。
 *
 *
 */
template <>
inline DEAL_II_ALWAYS_INLINE void
vectorized_load_and_transpose(const unsigned int          n_entries,
                              const float *               in,
                              const unsigned int *        offsets,
                              VectorizedArray<float, 16> *out)
{
  // Similar to the double case, we perform the work on smaller entities. In
  // this case, we start from 128 bit arrays and insert them into a full 512
  // bit index. This reduces the code size and register pressure because we do
  // shuffles on 4 numbers rather than 16.
  const unsigned int n_chunks = n_entries / 4;

  // To avoid warnings about uninitialized variables, need to initialize one
  // variable to a pre-exisiting value in out, which will never get used in
  // the end. Keep the initialization outside the loop because of a bug in
  // gcc-9.1 which generates a "vmovapd" instruction instead of "vmovupd" in
  // case t3 is initialized to zero (inside/outside of loop), see
  // https://gcc.gnu.org/bugzilla/show_bug.cgi?id=90991
  __m512 t0, t1, t2, t3;
  if (n_chunks > 0)
    t3 = out[0].data;
  for (unsigned int i = 0; i < n_chunks; ++i)
    {
      t0 = _mm512_insertf32x4(t3, _mm_loadu_ps(in + offsets[0] + 4 * i), 0);
      t0 = _mm512_insertf32x4(t0, _mm_loadu_ps(in + offsets[4] + 4 * i), 1);
      t0 = _mm512_insertf32x4(t0, _mm_loadu_ps(in + offsets[8] + 4 * i), 2);
      t0 = _mm512_insertf32x4(t0, _mm_loadu_ps(in + offsets[12] + 4 * i), 3);
      t1 = _mm512_insertf32x4(t3, _mm_loadu_ps(in + offsets[1] + 4 * i), 0);
      t1 = _mm512_insertf32x4(t1, _mm_loadu_ps(in + offsets[5] + 4 * i), 1);
      t1 = _mm512_insertf32x4(t1, _mm_loadu_ps(in + offsets[9] + 4 * i), 2);
      t1 = _mm512_insertf32x4(t1, _mm_loadu_ps(in + offsets[13] + 4 * i), 3);
      t2 = _mm512_insertf32x4(t3, _mm_loadu_ps(in + offsets[2] + 4 * i), 0);
      t2 = _mm512_insertf32x4(t2, _mm_loadu_ps(in + offsets[6] + 4 * i), 1);
      t2 = _mm512_insertf32x4(t2, _mm_loadu_ps(in + offsets[10] + 4 * i), 2);
      t2 = _mm512_insertf32x4(t2, _mm_loadu_ps(in + offsets[14] + 4 * i), 3);
      t3 = _mm512_insertf32x4(t3, _mm_loadu_ps(in + offsets[3] + 4 * i), 0);
      t3 = _mm512_insertf32x4(t3, _mm_loadu_ps(in + offsets[7] + 4 * i), 1);
      t3 = _mm512_insertf32x4(t3, _mm_loadu_ps(in + offsets[11] + 4 * i), 2);
      t3 = _mm512_insertf32x4(t3, _mm_loadu_ps(in + offsets[15] + 4 * i), 3);

      __m512 v0 = _mm512_shuffle_ps(t0, t1, 0x44);
      __m512 v1 = _mm512_shuffle_ps(t0, t1, 0xee);
      __m512 v2 = _mm512_shuffle_ps(t2, t3, 0x44);
      __m512 v3 = _mm512_shuffle_ps(t2, t3, 0xee);

      out[4 * i + 0].data = _mm512_shuffle_ps(v0, v2, 0x88);
      out[4 * i + 1].data = _mm512_shuffle_ps(v0, v2, 0xdd);
      out[4 * i + 2].data = _mm512_shuffle_ps(v1, v3, 0x88);
      out[4 * i + 3].data = _mm512_shuffle_ps(v1, v3, 0xdd);
    }

  // remainder loop of work that does not divide by 4
  for (unsigned int i = 4 * n_chunks; i < n_entries; ++i)
    out[i].gather(in + i, offsets);
}



/**
 * 浮点数和AVX-512的特殊化。
 *
 *
 */
template <>
inline DEAL_II_ALWAYS_INLINE void
vectorized_load_and_transpose(const unsigned int             n_entries,
                              const std::array<float *, 16> &in,
                              VectorizedArray<float, 16> *   out)
{
  // see the comments in the vectorized_load_and_transpose above

  const unsigned int n_chunks = n_entries / 4;

  __m512 t0, t1, t2, t3;
  if (n_chunks > 0)
    t3 = out[0].data;
  for (unsigned int i = 0; i < n_chunks; ++i)
    {
      t0 = _mm512_insertf32x4(t3, _mm_loadu_ps(in[0] + 4 * i), 0);
      t0 = _mm512_insertf32x4(t0, _mm_loadu_ps(in[4] + 4 * i), 1);
      t0 = _mm512_insertf32x4(t0, _mm_loadu_ps(in[8] + 4 * i), 2);
      t0 = _mm512_insertf32x4(t0, _mm_loadu_ps(in[12] + 4 * i), 3);
      t1 = _mm512_insertf32x4(t3, _mm_loadu_ps(in[1] + 4 * i), 0);
      t1 = _mm512_insertf32x4(t1, _mm_loadu_ps(in[5] + 4 * i), 1);
      t1 = _mm512_insertf32x4(t1, _mm_loadu_ps(in[9] + 4 * i), 2);
      t1 = _mm512_insertf32x4(t1, _mm_loadu_ps(in[13] + 4 * i), 3);
      t2 = _mm512_insertf32x4(t3, _mm_loadu_ps(in[2] + 4 * i), 0);
      t2 = _mm512_insertf32x4(t2, _mm_loadu_ps(in[6] + 4 * i), 1);
      t2 = _mm512_insertf32x4(t2, _mm_loadu_ps(in[10] + 4 * i), 2);
      t2 = _mm512_insertf32x4(t2, _mm_loadu_ps(in[14] + 4 * i), 3);
      t3 = _mm512_insertf32x4(t3, _mm_loadu_ps(in[3] + 4 * i), 0);
      t3 = _mm512_insertf32x4(t3, _mm_loadu_ps(in[7] + 4 * i), 1);
      t3 = _mm512_insertf32x4(t3, _mm_loadu_ps(in[11] + 4 * i), 2);
      t3 = _mm512_insertf32x4(t3, _mm_loadu_ps(in[15] + 4 * i), 3);

      __m512 v0 = _mm512_shuffle_ps(t0, t1, 0x44);
      __m512 v1 = _mm512_shuffle_ps(t0, t1, 0xee);
      __m512 v2 = _mm512_shuffle_ps(t2, t3, 0x44);
      __m512 v3 = _mm512_shuffle_ps(t2, t3, 0xee);

      out[4 * i + 0].data = _mm512_shuffle_ps(v0, v2, 0x88);
      out[4 * i + 1].data = _mm512_shuffle_ps(v0, v2, 0xdd);
      out[4 * i + 2].data = _mm512_shuffle_ps(v1, v3, 0x88);
      out[4 * i + 3].data = _mm512_shuffle_ps(v1, v3, 0xdd);
    }

  for (unsigned int i = 4 * n_chunks; i < n_entries; ++i)
    gather(out[i], in, i);
}



/**
 * 浮动和AVX-512的特殊化。
 *
 *
 */
template <>
inline DEAL_II_ALWAYS_INLINE void
vectorized_transpose_and_store(const bool                        add_into,
                               const unsigned int                n_entries,
                               const VectorizedArray<float, 16> *in,
                               const unsigned int *              offsets,
                               float *                           out)
{
  const unsigned int n_chunks = n_entries / 4;
  for (unsigned int i = 0; i < n_chunks; ++i)
    {
      __m512 t0 = _mm512_shuffle_ps(in[4 * i].data, in[1 + 4 * i].data, 0x44);
      __m512 t1 = _mm512_shuffle_ps(in[4 * i].data, in[1 + 4 * i].data, 0xee);
      __m512 t2 =
        _mm512_shuffle_ps(in[2 + 4 * i].data, in[3 + 4 * i].data, 0x44);
      __m512 t3 =
        _mm512_shuffle_ps(in[2 + 4 * i].data, in[3 + 4 * i].data, 0xee);
      __m512 u0 = _mm512_shuffle_ps(t0, t2, 0x88);
      __m512 u1 = _mm512_shuffle_ps(t0, t2, 0xdd);
      __m512 u2 = _mm512_shuffle_ps(t1, t3, 0x88);
      __m512 u3 = _mm512_shuffle_ps(t1, t3, 0xdd);

      __m128 res0  = _mm512_extractf32x4_ps(u0, 0);
      __m128 res4  = _mm512_extractf32x4_ps(u0, 1);
      __m128 res8  = _mm512_extractf32x4_ps(u0, 2);
      __m128 res12 = _mm512_extractf32x4_ps(u0, 3);
      __m128 res1  = _mm512_extractf32x4_ps(u1, 0);
      __m128 res5  = _mm512_extractf32x4_ps(u1, 1);
      __m128 res9  = _mm512_extractf32x4_ps(u1, 2);
      __m128 res13 = _mm512_extractf32x4_ps(u1, 3);
      __m128 res2  = _mm512_extractf32x4_ps(u2, 0);
      __m128 res6  = _mm512_extractf32x4_ps(u2, 1);
      __m128 res10 = _mm512_extractf32x4_ps(u2, 2);
      __m128 res14 = _mm512_extractf32x4_ps(u2, 3);
      __m128 res3  = _mm512_extractf32x4_ps(u3, 0);
      __m128 res7  = _mm512_extractf32x4_ps(u3, 1);
      __m128 res11 = _mm512_extractf32x4_ps(u3, 2);
      __m128 res15 = _mm512_extractf32x4_ps(u3, 3);

      // Cannot use the same store instructions in both paths of the 'if'
      // because the compiler cannot know that there is no aliasing between
      // pointers
      if (add_into)
        {
          res0 = _mm_add_ps(_mm_loadu_ps(out + 4 * i + offsets[0]), res0);
          _mm_storeu_ps(out + 4 * i + offsets[0], res0);
          res1 = _mm_add_ps(_mm_loadu_ps(out + 4 * i + offsets[1]), res1);
          _mm_storeu_ps(out + 4 * i + offsets[1], res1);
          res2 = _mm_add_ps(_mm_loadu_ps(out + 4 * i + offsets[2]), res2);
          _mm_storeu_ps(out + 4 * i + offsets[2], res2);
          res3 = _mm_add_ps(_mm_loadu_ps(out + 4 * i + offsets[3]), res3);
          _mm_storeu_ps(out + 4 * i + offsets[3], res3);
          res4 = _mm_add_ps(_mm_loadu_ps(out + 4 * i + offsets[4]), res4);
          _mm_storeu_ps(out + 4 * i + offsets[4], res4);
          res5 = _mm_add_ps(_mm_loadu_ps(out + 4 * i + offsets[5]), res5);
          _mm_storeu_ps(out + 4 * i + offsets[5], res5);
          res6 = _mm_add_ps(_mm_loadu_ps(out + 4 * i + offsets[6]), res6);
          _mm_storeu_ps(out + 4 * i + offsets[6], res6);
          res7 = _mm_add_ps(_mm_loadu_ps(out + 4 * i + offsets[7]), res7);
          _mm_storeu_ps(out + 4 * i + offsets[7], res7);
          res8 = _mm_add_ps(_mm_loadu_ps(out + 4 * i + offsets[8]), res8);
          _mm_storeu_ps(out + 4 * i + offsets[8], res8);
          res9 = _mm_add_ps(_mm_loadu_ps(out + 4 * i + offsets[9]), res9);
          _mm_storeu_ps(out + 4 * i + offsets[9], res9);
          res10 = _mm_add_ps(_mm_loadu_ps(out + 4 * i + offsets[10]), res10);
          _mm_storeu_ps(out + 4 * i + offsets[10], res10);
          res11 = _mm_add_ps(_mm_loadu_ps(out + 4 * i + offsets[11]), res11);
          _mm_storeu_ps(out + 4 * i + offsets[11], res11);
          res12 = _mm_add_ps(_mm_loadu_ps(out + 4 * i + offsets[12]), res12);
          _mm_storeu_ps(out + 4 * i + offsets[12], res12);
          res13 = _mm_add_ps(_mm_loadu_ps(out + 4 * i + offsets[13]), res13);
          _mm_storeu_ps(out + 4 * i + offsets[13], res13);
          res14 = _mm_add_ps(_mm_loadu_ps(out + 4 * i + offsets[14]), res14);
          _mm_storeu_ps(out + 4 * i + offsets[14], res14);
          res15 = _mm_add_ps(_mm_loadu_ps(out + 4 * i + offsets[15]), res15);
          _mm_storeu_ps(out + 4 * i + offsets[15], res15);
        }
      else
        {
          _mm_storeu_ps(out + 4 * i + offsets[0], res0);
          _mm_storeu_ps(out + 4 * i + offsets[1], res1);
          _mm_storeu_ps(out + 4 * i + offsets[2], res2);
          _mm_storeu_ps(out + 4 * i + offsets[3], res3);
          _mm_storeu_ps(out + 4 * i + offsets[4], res4);
          _mm_storeu_ps(out + 4 * i + offsets[5], res5);
          _mm_storeu_ps(out + 4 * i + offsets[6], res6);
          _mm_storeu_ps(out + 4 * i + offsets[7], res7);
          _mm_storeu_ps(out + 4 * i + offsets[8], res8);
          _mm_storeu_ps(out + 4 * i + offsets[9], res9);
          _mm_storeu_ps(out + 4 * i + offsets[10], res10);
          _mm_storeu_ps(out + 4 * i + offsets[11], res11);
          _mm_storeu_ps(out + 4 * i + offsets[12], res12);
          _mm_storeu_ps(out + 4 * i + offsets[13], res13);
          _mm_storeu_ps(out + 4 * i + offsets[14], res14);
          _mm_storeu_ps(out + 4 * i + offsets[15], res15);
        }
    }

  // remainder loop of work that does not divide by 4
  if (add_into)
    for (unsigned int i = 4 * n_chunks; i < n_entries; ++i)
      for (unsigned int v = 0; v < 16; ++v)
        out[offsets[v] + i] += in[i][v];
  else
    for (unsigned int i = 4 * n_chunks; i < n_entries; ++i)
      for (unsigned int v = 0; v < 16; ++v)
        out[offsets[v] + i] = in[i][v];
}



/**
 * 浮动和AVX-512的特殊化。
 *
 *
 */
template <>
inline DEAL_II_ALWAYS_INLINE void
vectorized_transpose_and_store(const bool                        add_into,
                               const unsigned int                n_entries,
                               const VectorizedArray<float, 16> *in,
                               std::array<float *, 16> &         out)
{
  // see the comments in the vectorized_transpose_and_store above

  const unsigned int n_chunks = n_entries / 4;
  for (unsigned int i = 0; i < n_chunks; ++i)
    {
      __m512 t0 = _mm512_shuffle_ps(in[4 * i].data, in[1 + 4 * i].data, 0x44);
      __m512 t1 = _mm512_shuffle_ps(in[4 * i].data, in[1 + 4 * i].data, 0xee);
      __m512 t2 =
        _mm512_shuffle_ps(in[2 + 4 * i].data, in[3 + 4 * i].data, 0x44);
      __m512 t3 =
        _mm512_shuffle_ps(in[2 + 4 * i].data, in[3 + 4 * i].data, 0xee);
      __m512 u0 = _mm512_shuffle_ps(t0, t2, 0x88);
      __m512 u1 = _mm512_shuffle_ps(t0, t2, 0xdd);
      __m512 u2 = _mm512_shuffle_ps(t1, t3, 0x88);
      __m512 u3 = _mm512_shuffle_ps(t1, t3, 0xdd);

      __m128 res0  = _mm512_extractf32x4_ps(u0, 0);
      __m128 res4  = _mm512_extractf32x4_ps(u0, 1);
      __m128 res8  = _mm512_extractf32x4_ps(u0, 2);
      __m128 res12 = _mm512_extractf32x4_ps(u0, 3);
      __m128 res1  = _mm512_extractf32x4_ps(u1, 0);
      __m128 res5  = _mm512_extractf32x4_ps(u1, 1);
      __m128 res9  = _mm512_extractf32x4_ps(u1, 2);
      __m128 res13 = _mm512_extractf32x4_ps(u1, 3);
      __m128 res2  = _mm512_extractf32x4_ps(u2, 0);
      __m128 res6  = _mm512_extractf32x4_ps(u2, 1);
      __m128 res10 = _mm512_extractf32x4_ps(u2, 2);
      __m128 res14 = _mm512_extractf32x4_ps(u2, 3);
      __m128 res3  = _mm512_extractf32x4_ps(u3, 0);
      __m128 res7  = _mm512_extractf32x4_ps(u3, 1);
      __m128 res11 = _mm512_extractf32x4_ps(u3, 2);
      __m128 res15 = _mm512_extractf32x4_ps(u3, 3);

      if (add_into)
        {
          res0 = _mm_add_ps(_mm_loadu_ps(out[0] + 4 * i), res0);
          _mm_storeu_ps(out[0] + 4 * i, res0);
          res1 = _mm_add_ps(_mm_loadu_ps(out[1] + 4 * i), res1);
          _mm_storeu_ps(out[1] + 4 * i, res1);
          res2 = _mm_add_ps(_mm_loadu_ps(out[2] + 4 * i), res2);
          _mm_storeu_ps(out[2] + 4 * i, res2);
          res3 = _mm_add_ps(_mm_loadu_ps(out[3] + 4 * i), res3);
          _mm_storeu_ps(out[3] + 4 * i, res3);
          res4 = _mm_add_ps(_mm_loadu_ps(out[4] + 4 * i), res4);
          _mm_storeu_ps(out[4] + 4 * i, res4);
          res5 = _mm_add_ps(_mm_loadu_ps(out[5] + 4 * i), res5);
          _mm_storeu_ps(out[5] + 4 * i, res5);
          res6 = _mm_add_ps(_mm_loadu_ps(out[6] + 4 * i), res6);
          _mm_storeu_ps(out[6] + 4 * i, res6);
          res7 = _mm_add_ps(_mm_loadu_ps(out[7] + 4 * i), res7);
          _mm_storeu_ps(out[7] + 4 * i, res7);
          res8 = _mm_add_ps(_mm_loadu_ps(out[8] + 4 * i), res8);
          _mm_storeu_ps(out[8] + 4 * i, res8);
          res9 = _mm_add_ps(_mm_loadu_ps(out[9] + 4 * i), res9);
          _mm_storeu_ps(out[9] + 4 * i, res9);
          res10 = _mm_add_ps(_mm_loadu_ps(out[10] + 4 * i), res10);
          _mm_storeu_ps(out[10] + 4 * i, res10);
          res11 = _mm_add_ps(_mm_loadu_ps(out[11] + 4 * i), res11);
          _mm_storeu_ps(out[11] + 4 * i, res11);
          res12 = _mm_add_ps(_mm_loadu_ps(out[12] + 4 * i), res12);
          _mm_storeu_ps(out[12] + 4 * i, res12);
          res13 = _mm_add_ps(_mm_loadu_ps(out[13] + 4 * i), res13);
          _mm_storeu_ps(out[13] + 4 * i, res13);
          res14 = _mm_add_ps(_mm_loadu_ps(out[14] + 4 * i), res14);
          _mm_storeu_ps(out[14] + 4 * i, res14);
          res15 = _mm_add_ps(_mm_loadu_ps(out[15] + 4 * i), res15);
          _mm_storeu_ps(out[15] + 4 * i, res15);
        }
      else
        {
          _mm_storeu_ps(out[0] + 4 * i, res0);
          _mm_storeu_ps(out[1] + 4 * i, res1);
          _mm_storeu_ps(out[2] + 4 * i, res2);
          _mm_storeu_ps(out[3] + 4 * i, res3);
          _mm_storeu_ps(out[4] + 4 * i, res4);
          _mm_storeu_ps(out[5] + 4 * i, res5);
          _mm_storeu_ps(out[6] + 4 * i, res6);
          _mm_storeu_ps(out[7] + 4 * i, res7);
          _mm_storeu_ps(out[8] + 4 * i, res8);
          _mm_storeu_ps(out[9] + 4 * i, res9);
          _mm_storeu_ps(out[10] + 4 * i, res10);
          _mm_storeu_ps(out[11] + 4 * i, res11);
          _mm_storeu_ps(out[12] + 4 * i, res12);
          _mm_storeu_ps(out[13] + 4 * i, res13);
          _mm_storeu_ps(out[14] + 4 * i, res14);
          _mm_storeu_ps(out[15] + 4 * i, res15);
        }
    }

  if (add_into)
    for (unsigned int i = 4 * n_chunks; i < n_entries; ++i)
      for (unsigned int v = 0; v < 16; ++v)
        out[v][i] += in[i][v];
  else
    for (unsigned int i = 4 * n_chunks; i < n_entries; ++i)
      for (unsigned int v = 0; v < 16; ++v)
        out[v][i] = in[i][v];
}

#  endif

#  if DEAL_II_VECTORIZATION_WIDTH_IN_BITS >= 256 && defined(__AVX__)

/**
 * 适用于双数和AVX的VectorizedArray类的特化。
 *
 *
 */
template <>
class VectorizedArray<double, 4>
  : public VectorizedArrayBase<VectorizedArray<double, 4>, 4>
{
public:
  /**
   * 这给出了数组元素的类型。
   *
   */
  using value_type = double;

  /**
   * 默认的空构造函数，让数据处于未初始化的状态，类似于float/double。
   *
   */
  VectorizedArray() = default;

  /**
   * 用给定的标量构建一个数组，广播给所有通道。
   *
   */
  VectorizedArray(const double scalar)
  {
    this->operator=(scalar);
  }

  /**
   * 这个函数可以用来将所有的数据字段设置为一个给定的标量。
   *
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray &
  operator=(const double x)
  {
    data = _mm256_set1_pd(x);
    return *this;
  }

  /**
   * 访问操作符。
   *
   */
  DEAL_II_ALWAYS_INLINE
  double &operator[](const unsigned int comp)
  {
    AssertIndexRange(comp, 4);
    return *(reinterpret_cast<double *>(&data) + comp);
  }

  /**
   * 常数访问运算符。
   *
   */
  DEAL_II_ALWAYS_INLINE
  const double &operator[](const unsigned int comp) const
  {
    AssertIndexRange(comp, 4);
    return *(reinterpret_cast<const double *>(&data) + comp);
  }

  /**
   * 加法。
   *
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray &
  operator+=(const VectorizedArray &vec)
  {
    // if the compiler supports vector arithmetic, we can simply use +=
    // operator on the given data type. this allows the compiler to combine
    // additions with multiplication (fused multiply-add) if those
    // instructions are available. Otherwise, we need to use the built-in
    // intrinsic command for __m256d
#    ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data += vec.data;
#    else
    data = _mm256_add_pd(data, vec.data);
#    endif
    return *this;
  }

  /**
   * 减法。
   *
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray &
  operator-=(const VectorizedArray &vec)
  {
#    ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data -= vec.data;
#    else
    data = _mm256_sub_pd(data, vec.data);
#    endif
    return *this;
  }
  /**
   * 乘法。
   *
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray &
  operator*=(const VectorizedArray &vec)
  {
#    ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data *= vec.data;
#    else
    data = _mm256_mul_pd(data, vec.data);
#    endif
    return *this;
  }

  /**
   * 除法。
   *
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray &
  operator/=(const VectorizedArray &vec)
  {
#    ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data /= vec.data;
#    else
    data = _mm256_div_pd(data, vec.data);
#    endif
    return *this;
  }

  /**
   * 从内存中加载 @p size()
   * 到调用类中，从给定的地址开始。内存不需要按32字节对齐，与之相对应的是将一个双倍地址投给VectorizedArray<double>*。
   *
   */
  DEAL_II_ALWAYS_INLINE
  void
  load(const double *ptr)
  {
    data = _mm256_loadu_pd(ptr);
  }

  /**
   * 将调用类的内容以 @p
   * size()的形式写入内存到给定地址。内存不需要按32字节对齐，相对于将双倍地址投给VectorizedArray<double>*来说。
   *
   */
  DEAL_II_ALWAYS_INLINE
  void
  store(double *ptr) const
  {
    _mm256_storeu_pd(ptr, data);
  }

  /**
   * @copydoc   VectorizedArray<Number>::streaming_store() 。
   * @note  内存必须以32字节对齐。
   *
   */
  DEAL_II_ALWAYS_INLINE
  void
  streaming_store(double *ptr) const
  {
    Assert(reinterpret_cast<std::size_t>(ptr) % 32 == 0,
           ExcMessage("Memory not aligned"));
    _mm256_stream_pd(ptr, data);
  }

  /**
   * 将 @p size()
   * 从内存中加载到调用的类中，从给定的地址和给定的偏移量开始，从偏移量开始的每个条目提供一个矢量数组的元素。
   * 这个操作对应于以下代码（但在硬件允许的情况下，使用了更有效的实现）。
   * @code
   * for (unsigned int v=0; v<VectorizedArray<Number>::size(); ++v)
   * this->operator[](v) = base_ptr[offsets[v]];
   * @endcode
   *
   */
  DEAL_II_ALWAYS_INLINE
  void
  gather(const double *base_ptr, const unsigned int *offsets)
  {
#    ifdef __AVX2__
    // unfortunately, there does not appear to be a 128 bit integer load, so
    // do it by some reinterpret casts here. this is allowed because the Intel
    // API allows aliasing between different vector types.
    const __m128 index_val =
      _mm_loadu_ps(reinterpret_cast<const float *>(offsets));
    const __m128i index = *reinterpret_cast<const __m128i *>(&index_val);
    data                = _mm256_i32gather_pd(base_ptr, index, 8);
#    else
    for (unsigned int i = 0; i < 4; ++i)
      *(reinterpret_cast<double *>(&data) + i) = base_ptr[offsets[i]];
#    endif
  }

  /**
   * 将调用类的内容以 @p
   * size()的形式写入内存，到给定的地址和给定的偏移量，将矢量数组的元素填入每个偏移量。
   * 这个操作对应于下面的代码（但在硬件允许的情况下使用了更有效的实现）。
   * @code
   * for (unsigned int v=0; v<VectorizedArray<Number>::size(); ++v)
   * base_ptr[offsets[v]] = this->operator[](v);
   * @endcode
   *
   *
   */
  DEAL_II_ALWAYS_INLINE
  void
  scatter(const unsigned int *offsets, double *base_ptr) const
  {
    // no scatter operation in AVX/AVX2
    for (unsigned int i = 0; i < 4; ++i)
      base_ptr[offsets[i]] = *(reinterpret_cast<const double *>(&data) + i);
  }

  /**
   * 实际的数据字段。为了与标准布局类型保持一致，并且能够与外部SIMD功能进行交互，这个成员被声明为公共的。
   *
   */
  __m256d data;

private:
  /**
   * 返回这个字段的平方根。不适合在用户代码中使用。使用sqrt(x)代替。
   *
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray
  get_sqrt() const
  {
    VectorizedArray res;
    res.data = _mm256_sqrt_pd(data);
    return res;
  }

  /**
   * 返回这个字段的绝对值。不适合在用户代码中使用。请使用abs(x)代替。
   *
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray
  get_abs() const
  {
    // to compute the absolute value, perform bitwise andnot with -0. This
    // will leave all value and exponent bits unchanged but force the sign
    // value to +.
    __m256d         mask = _mm256_set1_pd(-0.);
    VectorizedArray res;
    res.data = _mm256_andnot_pd(mask, data);
    return res;
  }

  /**
   * 返回这个字段和另一个字段的分量上的最大值。不适合在用户代码中使用。使用max(x,y)代替。
   *
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray
  get_max(const VectorizedArray &other) const
  {
    VectorizedArray res;
    res.data = _mm256_max_pd(data, other.data);
    return res;
  }

  /**
   * 返回这个字段和另一个字段的最小分量。不适合在用户代码中使用。使用min(x,y)代替。
   *
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray
  get_min(const VectorizedArray &other) const
  {
    VectorizedArray res;
    res.data = _mm256_min_pd(data, other.data);
    return res;
  }

  // Make a few functions friends.
  template <typename Number2, std::size_t width2>
  friend VectorizedArray<Number2, width2>
  std::sqrt(const VectorizedArray<Number2, width2> &);
  template <typename Number2, std::size_t width2>
  friend VectorizedArray<Number2, width2>
  std::abs(const VectorizedArray<Number2, width2> &);
  template <typename Number2, std::size_t width2>
  friend VectorizedArray<Number2, width2>
  std::max(const VectorizedArray<Number2, width2> &,
           const VectorizedArray<Number2, width2> &);
  template <typename Number2, std::size_t width2>
  friend VectorizedArray<Number2, width2>
  std::min(const VectorizedArray<Number2, width2> &,
           const VectorizedArray<Number2, width2> &);
};



/**
 * double和AVX的特殊化。
 *
 *
 */
template <>
inline DEAL_II_ALWAYS_INLINE void
vectorized_load_and_transpose(const unsigned int          n_entries,
                              const double *              in,
                              const unsigned int *        offsets,
                              VectorizedArray<double, 4> *out)
{
  const unsigned int n_chunks = n_entries / 4;
  const double *     in0      = in + offsets[0];
  const double *     in1      = in + offsets[1];
  const double *     in2      = in + offsets[2];
  const double *     in3      = in + offsets[3];

  for (unsigned int i = 0; i < n_chunks; ++i)
    {
      __m256d u0          = _mm256_loadu_pd(in0 + 4 * i);
      __m256d u1          = _mm256_loadu_pd(in1 + 4 * i);
      __m256d u2          = _mm256_loadu_pd(in2 + 4 * i);
      __m256d u3          = _mm256_loadu_pd(in3 + 4 * i);
      __m256d t0          = _mm256_permute2f128_pd(u0, u2, 0x20);
      __m256d t1          = _mm256_permute2f128_pd(u1, u3, 0x20);
      __m256d t2          = _mm256_permute2f128_pd(u0, u2, 0x31);
      __m256d t3          = _mm256_permute2f128_pd(u1, u3, 0x31);
      out[4 * i + 0].data = _mm256_unpacklo_pd(t0, t1);
      out[4 * i + 1].data = _mm256_unpackhi_pd(t0, t1);
      out[4 * i + 2].data = _mm256_unpacklo_pd(t2, t3);
      out[4 * i + 3].data = _mm256_unpackhi_pd(t2, t3);
    }

  // remainder loop of work that does not divide by 4
  for (unsigned int i = 4 * n_chunks; i < n_entries; ++i)
    out[i].gather(in + i, offsets);
}



/**
 * 双重和AVX的特殊化。
 *
 *
 */
template <>
inline DEAL_II_ALWAYS_INLINE void
vectorized_load_and_transpose(const unsigned int             n_entries,
                              const std::array<double *, 4> &in,
                              VectorizedArray<double, 4> *   out)
{
  // see the comments in the vectorized_load_and_transpose above

  const unsigned int n_chunks = n_entries / 4;
  const double *     in0      = in[0];
  const double *     in1      = in[1];
  const double *     in2      = in[2];
  const double *     in3      = in[3];

  for (unsigned int i = 0; i < n_chunks; ++i)
    {
      __m256d u0          = _mm256_loadu_pd(in0 + 4 * i);
      __m256d u1          = _mm256_loadu_pd(in1 + 4 * i);
      __m256d u2          = _mm256_loadu_pd(in2 + 4 * i);
      __m256d u3          = _mm256_loadu_pd(in3 + 4 * i);
      __m256d t0          = _mm256_permute2f128_pd(u0, u2, 0x20);
      __m256d t1          = _mm256_permute2f128_pd(u1, u3, 0x20);
      __m256d t2          = _mm256_permute2f128_pd(u0, u2, 0x31);
      __m256d t3          = _mm256_permute2f128_pd(u1, u3, 0x31);
      out[4 * i + 0].data = _mm256_unpacklo_pd(t0, t1);
      out[4 * i + 1].data = _mm256_unpackhi_pd(t0, t1);
      out[4 * i + 2].data = _mm256_unpacklo_pd(t2, t3);
      out[4 * i + 3].data = _mm256_unpackhi_pd(t2, t3);
    }

  for (unsigned int i = 4 * n_chunks; i < n_entries; ++i)
    gather(out[i], in, i);
}



/**
 * 为双数和AVX的特殊化。
 *
 *
 */
template <>
inline DEAL_II_ALWAYS_INLINE void
vectorized_transpose_and_store(const bool                        add_into,
                               const unsigned int                n_entries,
                               const VectorizedArray<double, 4> *in,
                               const unsigned int *              offsets,
                               double *                          out)
{
  const unsigned int n_chunks = n_entries / 4;
  double *           out0     = out + offsets[0];
  double *           out1     = out + offsets[1];
  double *           out2     = out + offsets[2];
  double *           out3     = out + offsets[3];
  for (unsigned int i = 0; i < n_chunks; ++i)
    {
      __m256d u0   = in[4 * i + 0].data;
      __m256d u1   = in[4 * i + 1].data;
      __m256d u2   = in[4 * i + 2].data;
      __m256d u3   = in[4 * i + 3].data;
      __m256d t0   = _mm256_permute2f128_pd(u0, u2, 0x20);
      __m256d t1   = _mm256_permute2f128_pd(u1, u3, 0x20);
      __m256d t2   = _mm256_permute2f128_pd(u0, u2, 0x31);
      __m256d t3   = _mm256_permute2f128_pd(u1, u3, 0x31);
      __m256d res0 = _mm256_unpacklo_pd(t0, t1);
      __m256d res1 = _mm256_unpackhi_pd(t0, t1);
      __m256d res2 = _mm256_unpacklo_pd(t2, t3);
      __m256d res3 = _mm256_unpackhi_pd(t2, t3);

      // Cannot use the same store instructions in both paths of the 'if'
      // because the compiler cannot know that there is no aliasing between
      // pointers
      if (add_into)
        {
          res0 = _mm256_add_pd(_mm256_loadu_pd(out0 + 4 * i), res0);
          _mm256_storeu_pd(out0 + 4 * i, res0);
          res1 = _mm256_add_pd(_mm256_loadu_pd(out1 + 4 * i), res1);
          _mm256_storeu_pd(out1 + 4 * i, res1);
          res2 = _mm256_add_pd(_mm256_loadu_pd(out2 + 4 * i), res2);
          _mm256_storeu_pd(out2 + 4 * i, res2);
          res3 = _mm256_add_pd(_mm256_loadu_pd(out3 + 4 * i), res3);
          _mm256_storeu_pd(out3 + 4 * i, res3);
        }
      else
        {
          _mm256_storeu_pd(out0 + 4 * i, res0);
          _mm256_storeu_pd(out1 + 4 * i, res1);
          _mm256_storeu_pd(out2 + 4 * i, res2);
          _mm256_storeu_pd(out3 + 4 * i, res3);
        }
    }

  // remainder loop of work that does not divide by 4
  if (add_into)
    for (unsigned int i = 4 * n_chunks; i < n_entries; ++i)
      for (unsigned int v = 0; v < 4; ++v)
        out[offsets[v] + i] += in[i][v];
  else
    for (unsigned int i = 4 * n_chunks; i < n_entries; ++i)
      for (unsigned int v = 0; v < 4; ++v)
        out[offsets[v] + i] = in[i][v];
}



/**
 * 双重和AVX的特殊化。
 *
 *
 */
template <>
inline DEAL_II_ALWAYS_INLINE void
vectorized_transpose_and_store(const bool                        add_into,
                               const unsigned int                n_entries,
                               const VectorizedArray<double, 4> *in,
                               std::array<double *, 4> &         out)
{
  // see the comments in the vectorized_transpose_and_store above

  const unsigned int n_chunks = n_entries / 4;
  double *           out0     = out[0];
  double *           out1     = out[1];
  double *           out2     = out[2];
  double *           out3     = out[3];
  for (unsigned int i = 0; i < n_chunks; ++i)
    {
      __m256d u0   = in[4 * i + 0].data;
      __m256d u1   = in[4 * i + 1].data;
      __m256d u2   = in[4 * i + 2].data;
      __m256d u3   = in[4 * i + 3].data;
      __m256d t0   = _mm256_permute2f128_pd(u0, u2, 0x20);
      __m256d t1   = _mm256_permute2f128_pd(u1, u3, 0x20);
      __m256d t2   = _mm256_permute2f128_pd(u0, u2, 0x31);
      __m256d t3   = _mm256_permute2f128_pd(u1, u3, 0x31);
      __m256d res0 = _mm256_unpacklo_pd(t0, t1);
      __m256d res1 = _mm256_unpackhi_pd(t0, t1);
      __m256d res2 = _mm256_unpacklo_pd(t2, t3);
      __m256d res3 = _mm256_unpackhi_pd(t2, t3);

      // Cannot use the same store instructions in both paths of the 'if'
      // because the compiler cannot know that there is no aliasing between
      // pointers
      if (add_into)
        {
          res0 = _mm256_add_pd(_mm256_loadu_pd(out0 + 4 * i), res0);
          _mm256_storeu_pd(out0 + 4 * i, res0);
          res1 = _mm256_add_pd(_mm256_loadu_pd(out1 + 4 * i), res1);
          _mm256_storeu_pd(out1 + 4 * i, res1);
          res2 = _mm256_add_pd(_mm256_loadu_pd(out2 + 4 * i), res2);
          _mm256_storeu_pd(out2 + 4 * i, res2);
          res3 = _mm256_add_pd(_mm256_loadu_pd(out3 + 4 * i), res3);
          _mm256_storeu_pd(out3 + 4 * i, res3);
        }
      else
        {
          _mm256_storeu_pd(out0 + 4 * i, res0);
          _mm256_storeu_pd(out1 + 4 * i, res1);
          _mm256_storeu_pd(out2 + 4 * i, res2);
          _mm256_storeu_pd(out3 + 4 * i, res3);
        }
    }

  // remainder loop of work that does not divide by 4
  if (add_into)
    for (unsigned int i = 4 * n_chunks; i < n_entries; ++i)
      for (unsigned int v = 0; v < 4; ++v)
        out[v][i] += in[i][v];
  else
    for (unsigned int i = 4 * n_chunks; i < n_entries; ++i)
      for (unsigned int v = 0; v < 4; ++v)
        out[v][i] = in[i][v];
}



/**
 * 浮动和AVX的特殊化。
 *
 *
 */
template <>
class VectorizedArray<float, 8>
  : public VectorizedArrayBase<VectorizedArray<float, 8>, 8>
{
public:
  /**
   * 这给出了数组元素的类型。
   *
   */
  using value_type = float;

  /**
   * 默认的空构造函数，使数据处于未初始化的状态，类似于float/double。
   *
   */
  VectorizedArray() = default;

  /**
   * 用给定的标量构建一个数组广播到所有的通道。
   *
   */
  VectorizedArray(const float scalar)
  {
    this->operator=(scalar);
  }

  /**
   * 这个函数可以用来将所有的数据字段设置为一个给定的标量。
   *
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray &
  operator=(const float x)
  {
    data = _mm256_set1_ps(x);
    return *this;
  }

  /**
   * 访问操作符。
   *
   */
  DEAL_II_ALWAYS_INLINE
  float &operator[](const unsigned int comp)
  {
    AssertIndexRange(comp, 8);
    return *(reinterpret_cast<float *>(&data) + comp);
  }

  /**
   * 常数访问运算符。
   *
   */
  DEAL_II_ALWAYS_INLINE
  const float &operator[](const unsigned int comp) const
  {
    AssertIndexRange(comp, 8);
    return *(reinterpret_cast<const float *>(&data) + comp);
  }

  /**
   * 加法。
   *
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray &
  operator+=(const VectorizedArray &vec)
  {
    // if the compiler supports vector arithmetic, we can simply use +=
    // operator on the given data type. this allows the compiler to combine
    // additions with multiplication (fused multiply-add) if those
    // instructions are available. Otherwise, we need to use the built-in
    // intrinsic command for __m256d
#    ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data += vec.data;
#    else
    data = _mm256_add_ps(data, vec.data);
#    endif
    return *this;
  }

  /**
   * 减法。
   *
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray &
  operator-=(const VectorizedArray &vec)
  {
#    ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data -= vec.data;
#    else
    data = _mm256_sub_ps(data, vec.data);
#    endif
    return *this;
  }
  /**
   * 乘法。
   *
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray &
  operator*=(const VectorizedArray &vec)
  {
#    ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data *= vec.data;
#    else
    data = _mm256_mul_ps(data, vec.data);
#    endif
    return *this;
  }

  /**
   * 除法。
   *
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray &
  operator/=(const VectorizedArray &vec)
  {
#    ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data /= vec.data;
#    else
    data = _mm256_div_ps(data, vec.data);
#    endif
    return *this;
  }

  /**
   * 从内存中加载 @p size()
   * 到调用类中，从给定的地址开始。内存不需要按32字节对齐，与之相对应的是将一个浮点地址投给VectorizedArray<float>*。
   *
   */
  DEAL_II_ALWAYS_INLINE
  void
  load(const float *ptr)
  {
    data = _mm256_loadu_ps(ptr);
  }

  /**
   * 将调用类的内容以 @p
   * size()的形式写入内存到给定地址。内存不需要按32字节对齐，与之相对应的是将浮点地址投给VectorizedArray<float>*。
   *
   */
  DEAL_II_ALWAYS_INLINE
  void
  store(float *ptr) const
  {
    _mm256_storeu_ps(ptr, data);
  }

  /**
   * @copydoc   VectorizedArray<Number>::streaming_store() 。
   * @note  内存必须以32字节对齐。
   *
   */
  DEAL_II_ALWAYS_INLINE
  void
  streaming_store(float *ptr) const
  {
    Assert(reinterpret_cast<std::size_t>(ptr) % 32 == 0,
           ExcMessage("Memory not aligned"));
    _mm256_stream_ps(ptr, data);
  }

  /**
   * 将 @p size()
   * 从内存中加载到调用的类中，从给定的地址和给定的偏移量开始，从偏移量开始的每个条目提供一个矢量数组的元素。
   * 这个操作对应于以下代码（但在硬件允许的情况下，使用了更有效的实现）。
   * @code
   * for (unsigned int v=0; v<VectorizedArray<Number>::size(); ++v)
   * this->operator[](v) = base_ptr[offsets[v]];
   * @endcode
   *
   */
  DEAL_II_ALWAYS_INLINE
  void
  gather(const float *base_ptr, const unsigned int *offsets)
  {
#    ifdef __AVX2__
    // unfortunately, there does not appear to be a 256 bit integer load, so
    // do it by some reinterpret casts here. this is allowed because the Intel
    // API allows aliasing between different vector types.
    const __m256 index_val =
      _mm256_loadu_ps(reinterpret_cast<const float *>(offsets));
    const __m256i index = *reinterpret_cast<const __m256i *>(&index_val);
    data                = _mm256_i32gather_ps(base_ptr, index, 4);
#    else
    for (unsigned int i = 0; i < 8; ++i)
      *(reinterpret_cast<float *>(&data) + i) = base_ptr[offsets[i]];
#    endif
  }

  /**
   * 将调用类的内容以 @p
   * size()的形式写入内存，到给定的地址和给定的偏移量，将矢量数组的元素填入每个偏移量。
   * 这个操作对应于下面的代码（但在硬件允许的情况下使用了更有效的实现）。
   * @code
   * for (unsigned int v=0; v<VectorizedArray<Number>::size(); ++v)
   * base_ptr[offsets[v]] = this->operator[](v);
   * @endcode
   *
   *
   */
  DEAL_II_ALWAYS_INLINE
  void
  scatter(const unsigned int *offsets, float *base_ptr) const
  {
    // no scatter operation in AVX/AVX2
    for (unsigned int i = 0; i < 8; ++i)
      base_ptr[offsets[i]] = *(reinterpret_cast<const float *>(&data) + i);
  }

  /**
   * 实际的数据字段。为了与标准布局类型保持一致，并能与外部SIMD功能进行交互，这个成员被声明为公共的。
   *
   */
  __m256 data;

private:
  /**
   * 返回这个字段的平方根。不适合在用户代码中使用。使用sqrt(x)代替。
   *
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray
  get_sqrt() const
  {
    VectorizedArray res;
    res.data = _mm256_sqrt_ps(data);
    return res;
  }

  /**
   * 返回这个字段的绝对值。不适合在用户代码中使用。请使用abs(x)代替。
   *
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray
  get_abs() const
  {
    // to compute the absolute value, perform bitwise andnot with -0. This
    // will leave all value and exponent bits unchanged but force the sign
    // value to +.
    __m256          mask = _mm256_set1_ps(-0.f);
    VectorizedArray res;
    res.data = _mm256_andnot_ps(mask, data);
    return res;
  }

  /**
   * 返回这个字段和另一个字段的分量上的最大值。不适合在用户代码中使用。使用max(x,y)代替。
   *
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray
  get_max(const VectorizedArray &other) const
  {
    VectorizedArray res;
    res.data = _mm256_max_ps(data, other.data);
    return res;
  }

  /**
   * 返回这个字段和另一个字段的最小分量。不适合在用户代码中使用。使用min(x,y)代替。
   *
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray
  get_min(const VectorizedArray &other) const
  {
    VectorizedArray res;
    res.data = _mm256_min_ps(data, other.data);
    return res;
  }

  // Make a few functions friends.
  template <typename Number2, std::size_t width2>
  friend VectorizedArray<Number2, width2>
  std::sqrt(const VectorizedArray<Number2, width2> &);
  template <typename Number2, std::size_t width2>
  friend VectorizedArray<Number2, width2>
  std::abs(const VectorizedArray<Number2, width2> &);
  template <typename Number2, std::size_t width2>
  friend VectorizedArray<Number2, width2>
  std::max(const VectorizedArray<Number2, width2> &,
           const VectorizedArray<Number2, width2> &);
  template <typename Number2, std::size_t width2>
  friend VectorizedArray<Number2, width2>
  std::min(const VectorizedArray<Number2, width2> &,
           const VectorizedArray<Number2, width2> &);
};



/**
 * 浮点数和AVX的特殊化。
 *
 *
 */
template <>
inline DEAL_II_ALWAYS_INLINE void
vectorized_load_and_transpose(const unsigned int         n_entries,
                              const float *              in,
                              const unsigned int *       offsets,
                              VectorizedArray<float, 8> *out)
{
  const unsigned int n_chunks = n_entries / 4;
  for (unsigned int i = 0; i < n_chunks; ++i)
    {
      // To avoid warnings about uninitialized variables, need to initialize
      // one variable with zero before using it.
      __m256 t0, t1, t2, t3 = {};
      t0 = _mm256_insertf128_ps(t3, _mm_loadu_ps(in + 4 * i + offsets[0]), 0);
      t0 = _mm256_insertf128_ps(t0, _mm_loadu_ps(in + 4 * i + offsets[4]), 1);
      t1 = _mm256_insertf128_ps(t3, _mm_loadu_ps(in + 4 * i + offsets[1]), 0);
      t1 = _mm256_insertf128_ps(t1, _mm_loadu_ps(in + 4 * i + offsets[5]), 1);
      t2 = _mm256_insertf128_ps(t3, _mm_loadu_ps(in + 4 * i + offsets[2]), 0);
      t2 = _mm256_insertf128_ps(t2, _mm_loadu_ps(in + 4 * i + offsets[6]), 1);
      t3 = _mm256_insertf128_ps(t3, _mm_loadu_ps(in + 4 * i + offsets[3]), 0);
      t3 = _mm256_insertf128_ps(t3, _mm_loadu_ps(in + 4 * i + offsets[7]), 1);

      __m256 v0           = _mm256_shuffle_ps(t0, t1, 0x44);
      __m256 v1           = _mm256_shuffle_ps(t0, t1, 0xee);
      __m256 v2           = _mm256_shuffle_ps(t2, t3, 0x44);
      __m256 v3           = _mm256_shuffle_ps(t2, t3, 0xee);
      out[4 * i + 0].data = _mm256_shuffle_ps(v0, v2, 0x88);
      out[4 * i + 1].data = _mm256_shuffle_ps(v0, v2, 0xdd);
      out[4 * i + 2].data = _mm256_shuffle_ps(v1, v3, 0x88);
      out[4 * i + 3].data = _mm256_shuffle_ps(v1, v3, 0xdd);
    }

  // remainder loop of work that does not divide by 4
  for (unsigned int i = 4 * n_chunks; i < n_entries; ++i)
    out[i].gather(in + i, offsets);
}



/**
 * 浮点和AVX的特殊化。
 *
 *
 */
template <>
inline DEAL_II_ALWAYS_INLINE void
vectorized_load_and_transpose(const unsigned int            n_entries,
                              const std::array<float *, 8> &in,
                              VectorizedArray<float, 8> *   out)
{
  // see the comments in the vectorized_load_and_transpose above

  const unsigned int n_chunks = n_entries / 4;
  for (unsigned int i = 0; i < n_chunks; ++i)
    {
      __m256 t0, t1, t2, t3 = {};
      t0 = _mm256_insertf128_ps(t3, _mm_loadu_ps(in[0] + 4 * i), 0);
      t0 = _mm256_insertf128_ps(t0, _mm_loadu_ps(in[4] + 4 * i), 1);
      t1 = _mm256_insertf128_ps(t3, _mm_loadu_ps(in[1] + 4 * i), 0);
      t1 = _mm256_insertf128_ps(t1, _mm_loadu_ps(in[5] + 4 * i), 1);
      t2 = _mm256_insertf128_ps(t3, _mm_loadu_ps(in[2] + 4 * i), 0);
      t2 = _mm256_insertf128_ps(t2, _mm_loadu_ps(in[6] + 4 * i), 1);
      t3 = _mm256_insertf128_ps(t3, _mm_loadu_ps(in[3] + 4 * i), 0);
      t3 = _mm256_insertf128_ps(t3, _mm_loadu_ps(in[7] + 4 * i), 1);

      __m256 v0           = _mm256_shuffle_ps(t0, t1, 0x44);
      __m256 v1           = _mm256_shuffle_ps(t0, t1, 0xee);
      __m256 v2           = _mm256_shuffle_ps(t2, t3, 0x44);
      __m256 v3           = _mm256_shuffle_ps(t2, t3, 0xee);
      out[4 * i + 0].data = _mm256_shuffle_ps(v0, v2, 0x88);
      out[4 * i + 1].data = _mm256_shuffle_ps(v0, v2, 0xdd);
      out[4 * i + 2].data = _mm256_shuffle_ps(v1, v3, 0x88);
      out[4 * i + 3].data = _mm256_shuffle_ps(v1, v3, 0xdd);
    }

  for (unsigned int i = 4 * n_chunks; i < n_entries; ++i)
    gather(out[i], in, i);
}



/**
 * 浮动和AVX的特殊化。
 *
 *
 */
template <>
inline DEAL_II_ALWAYS_INLINE void
vectorized_transpose_and_store(const bool                       add_into,
                               const unsigned int               n_entries,
                               const VectorizedArray<float, 8> *in,
                               const unsigned int *             offsets,
                               float *                          out)
{
  const unsigned int n_chunks = n_entries / 4;
  for (unsigned int i = 0; i < n_chunks; ++i)
    {
      __m256 u0   = in[4 * i + 0].data;
      __m256 u1   = in[4 * i + 1].data;
      __m256 u2   = in[4 * i + 2].data;
      __m256 u3   = in[4 * i + 3].data;
      __m256 t0   = _mm256_shuffle_ps(u0, u1, 0x44);
      __m256 t1   = _mm256_shuffle_ps(u0, u1, 0xee);
      __m256 t2   = _mm256_shuffle_ps(u2, u3, 0x44);
      __m256 t3   = _mm256_shuffle_ps(u2, u3, 0xee);
      u0          = _mm256_shuffle_ps(t0, t2, 0x88);
      u1          = _mm256_shuffle_ps(t0, t2, 0xdd);
      u2          = _mm256_shuffle_ps(t1, t3, 0x88);
      u3          = _mm256_shuffle_ps(t1, t3, 0xdd);
      __m128 res0 = _mm256_extractf128_ps(u0, 0);
      __m128 res4 = _mm256_extractf128_ps(u0, 1);
      __m128 res1 = _mm256_extractf128_ps(u1, 0);
      __m128 res5 = _mm256_extractf128_ps(u1, 1);
      __m128 res2 = _mm256_extractf128_ps(u2, 0);
      __m128 res6 = _mm256_extractf128_ps(u2, 1);
      __m128 res3 = _mm256_extractf128_ps(u3, 0);
      __m128 res7 = _mm256_extractf128_ps(u3, 1);

      // Cannot use the same store instructions in both paths of the 'if'
      // because the compiler cannot know that there is no aliasing between
      // pointers
      if (add_into)
        {
          res0 = _mm_add_ps(_mm_loadu_ps(out + 4 * i + offsets[0]), res0);
          _mm_storeu_ps(out + 4 * i + offsets[0], res0);
          res1 = _mm_add_ps(_mm_loadu_ps(out + 4 * i + offsets[1]), res1);
          _mm_storeu_ps(out + 4 * i + offsets[1], res1);
          res2 = _mm_add_ps(_mm_loadu_ps(out + 4 * i + offsets[2]), res2);
          _mm_storeu_ps(out + 4 * i + offsets[2], res2);
          res3 = _mm_add_ps(_mm_loadu_ps(out + 4 * i + offsets[3]), res3);
          _mm_storeu_ps(out + 4 * i + offsets[3], res3);
          res4 = _mm_add_ps(_mm_loadu_ps(out + 4 * i + offsets[4]), res4);
          _mm_storeu_ps(out + 4 * i + offsets[4], res4);
          res5 = _mm_add_ps(_mm_loadu_ps(out + 4 * i + offsets[5]), res5);
          _mm_storeu_ps(out + 4 * i + offsets[5], res5);
          res6 = _mm_add_ps(_mm_loadu_ps(out + 4 * i + offsets[6]), res6);
          _mm_storeu_ps(out + 4 * i + offsets[6], res6);
          res7 = _mm_add_ps(_mm_loadu_ps(out + 4 * i + offsets[7]), res7);
          _mm_storeu_ps(out + 4 * i + offsets[7], res7);
        }
      else
        {
          _mm_storeu_ps(out + 4 * i + offsets[0], res0);
          _mm_storeu_ps(out + 4 * i + offsets[1], res1);
          _mm_storeu_ps(out + 4 * i + offsets[2], res2);
          _mm_storeu_ps(out + 4 * i + offsets[3], res3);
          _mm_storeu_ps(out + 4 * i + offsets[4], res4);
          _mm_storeu_ps(out + 4 * i + offsets[5], res5);
          _mm_storeu_ps(out + 4 * i + offsets[6], res6);
          _mm_storeu_ps(out + 4 * i + offsets[7], res7);
        }
    }

  // remainder loop of work that does not divide by 4
  if (add_into)
    for (unsigned int i = 4 * n_chunks; i < n_entries; ++i)
      for (unsigned int v = 0; v < 8; ++v)
        out[offsets[v] + i] += in[i][v];
  else
    for (unsigned int i = 4 * n_chunks; i < n_entries; ++i)
      for (unsigned int v = 0; v < 8; ++v)
        out[offsets[v] + i] = in[i][v];
}



/**
 * 浮动和AVX的特殊化。
 *
 *
 */
template <>
inline DEAL_II_ALWAYS_INLINE void
vectorized_transpose_and_store(const bool                       add_into,
                               const unsigned int               n_entries,
                               const VectorizedArray<float, 8> *in,
                               std::array<float *, 8> &         out)
{
  // see the comments in the vectorized_transpose_and_store above

  const unsigned int n_chunks = n_entries / 4;
  for (unsigned int i = 0; i < n_chunks; ++i)
    {
      __m256 u0   = in[4 * i + 0].data;
      __m256 u1   = in[4 * i + 1].data;
      __m256 u2   = in[4 * i + 2].data;
      __m256 u3   = in[4 * i + 3].data;
      __m256 t0   = _mm256_shuffle_ps(u0, u1, 0x44);
      __m256 t1   = _mm256_shuffle_ps(u0, u1, 0xee);
      __m256 t2   = _mm256_shuffle_ps(u2, u3, 0x44);
      __m256 t3   = _mm256_shuffle_ps(u2, u3, 0xee);
      u0          = _mm256_shuffle_ps(t0, t2, 0x88);
      u1          = _mm256_shuffle_ps(t0, t2, 0xdd);
      u2          = _mm256_shuffle_ps(t1, t3, 0x88);
      u3          = _mm256_shuffle_ps(t1, t3, 0xdd);
      __m128 res0 = _mm256_extractf128_ps(u0, 0);
      __m128 res4 = _mm256_extractf128_ps(u0, 1);
      __m128 res1 = _mm256_extractf128_ps(u1, 0);
      __m128 res5 = _mm256_extractf128_ps(u1, 1);
      __m128 res2 = _mm256_extractf128_ps(u2, 0);
      __m128 res6 = _mm256_extractf128_ps(u2, 1);
      __m128 res3 = _mm256_extractf128_ps(u3, 0);
      __m128 res7 = _mm256_extractf128_ps(u3, 1);

      if (add_into)
        {
          res0 = _mm_add_ps(_mm_loadu_ps(out[0] + 4 * i), res0);
          _mm_storeu_ps(out[0] + 4 * i, res0);
          res1 = _mm_add_ps(_mm_loadu_ps(out[1] + 4 * i), res1);
          _mm_storeu_ps(out[1] + 4 * i, res1);
          res2 = _mm_add_ps(_mm_loadu_ps(out[2] + 4 * i), res2);
          _mm_storeu_ps(out[2] + 4 * i, res2);
          res3 = _mm_add_ps(_mm_loadu_ps(out[3] + 4 * i), res3);
          _mm_storeu_ps(out[3] + 4 * i, res3);
          res4 = _mm_add_ps(_mm_loadu_ps(out[4] + 4 * i), res4);
          _mm_storeu_ps(out[4] + 4 * i, res4);
          res5 = _mm_add_ps(_mm_loadu_ps(out[5] + 4 * i), res5);
          _mm_storeu_ps(out[5] + 4 * i, res5);
          res6 = _mm_add_ps(_mm_loadu_ps(out[6] + 4 * i), res6);
          _mm_storeu_ps(out[6] + 4 * i, res6);
          res7 = _mm_add_ps(_mm_loadu_ps(out[7] + 4 * i), res7);
          _mm_storeu_ps(out[7] + 4 * i, res7);
        }
      else
        {
          _mm_storeu_ps(out[0] + 4 * i, res0);
          _mm_storeu_ps(out[1] + 4 * i, res1);
          _mm_storeu_ps(out[2] + 4 * i, res2);
          _mm_storeu_ps(out[3] + 4 * i, res3);
          _mm_storeu_ps(out[4] + 4 * i, res4);
          _mm_storeu_ps(out[5] + 4 * i, res5);
          _mm_storeu_ps(out[6] + 4 * i, res6);
          _mm_storeu_ps(out[7] + 4 * i, res7);
        }
    }

  if (add_into)
    for (unsigned int i = 4 * n_chunks; i < n_entries; ++i)
      for (unsigned int v = 0; v < 8; ++v)
        out[v][i] += in[i][v];
  else
    for (unsigned int i = 4 * n_chunks; i < n_entries; ++i)
      for (unsigned int v = 0; v < 8; ++v)
        out[v][i] = in[i][v];
}

#  endif

#  if DEAL_II_VECTORIZATION_WIDTH_IN_BITS >= 128 && defined(__SSE2__)

/**
 * 对双倍数和SSE2的特殊化。
 *
 *
 */
template <>
class VectorizedArray<double, 2>
  : public VectorizedArrayBase<VectorizedArray<double, 2>, 2>
{
public:
  /**
   * 这给出了数组元素的类型。
   *
   */
  using value_type = double;

  /**
   * 默认的空构造函数，让数据处于未初始化的状态，类似于float/double。
   *
   */
  VectorizedArray() = default;

  /**
   * 用给定的标量构建一个数组广播到所有的通道。
   *
   */
  VectorizedArray(const double scalar)
  {
    this->operator=(scalar);
  }

  /**
   * 这个函数可以用来将所有的数据字段设置为一个给定的标量。
   *
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray &
  operator=(const double x)
  {
    data = _mm_set1_pd(x);
    return *this;
  }

  /**
   * 访问操作符。
   *
   */
  DEAL_II_ALWAYS_INLINE
  double &operator[](const unsigned int comp)
  {
    AssertIndexRange(comp, 2);
    return *(reinterpret_cast<double *>(&data) + comp);
  }

  /**
   * 常数访问运算符。
   *
   */
  DEAL_II_ALWAYS_INLINE
  const double &operator[](const unsigned int comp) const
  {
    AssertIndexRange(comp, 2);
    return *(reinterpret_cast<const double *>(&data) + comp);
  }

  /**
   * 加法。
   *
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray &
  operator+=(const VectorizedArray &vec)
  {
#    ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data += vec.data;
#    else
    data = _mm_add_pd(data, vec.data);
#    endif
    return *this;
  }

  /**
   * 减法。
   *
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray &
  operator-=(const VectorizedArray &vec)
  {
#    ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data -= vec.data;
#    else
    data = _mm_sub_pd(data, vec.data);
#    endif
    return *this;
  }

  /**
   * 乘法。
   *
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray &
  operator*=(const VectorizedArray &vec)
  {
#    ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data *= vec.data;
#    else
    data = _mm_mul_pd(data, vec.data);
#    endif
    return *this;
  }

  /**
   * 除法。
   *
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray &
  operator/=(const VectorizedArray &vec)
  {
#    ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data /= vec.data;
#    else
    data = _mm_div_pd(data, vec.data);
#    endif
    return *this;
  }

  /**
   * 从内存中加载 @p size()
   * 到调用类中，从给定的地址开始。内存不需要对齐16个字节，相对于将双倍地址投给VectorizedArray<double>*来说。
   *
   */
  DEAL_II_ALWAYS_INLINE
  void
  load(const double *ptr)
  {
    data = _mm_loadu_pd(ptr);
  }

  /**
   * 将调用类的内容以 @p
   * size()的形式写到内存中的给定地址。内存不需要对齐16个字节，相对于将双倍地址投给VectorizedArray<double>*来说。
   *
   */
  DEAL_II_ALWAYS_INLINE
  void
  store(double *ptr) const
  {
    _mm_storeu_pd(ptr, data);
  }

  /**
   * @copydoc   VectorizedArray<Number>::streaming_store() 。
   * @note  内存必须以16字节对齐。
   *
   */
  DEAL_II_ALWAYS_INLINE
  void
  streaming_store(double *ptr) const
  {
    Assert(reinterpret_cast<std::size_t>(ptr) % 16 == 0,
           ExcMessage("Memory not aligned"));
    _mm_stream_pd(ptr, data);
  }

  /**
   * 将 @p size()
   * 从内存中加载到调用的类中，从给定的地址和给定的偏移量开始，从偏移量开始的每个条目提供一个矢量数组的元素。
   * 这个操作对应于以下代码（但在硬件允许的情况下，使用了更有效的实现）。
   * @code
   * for (unsigned int v=0; v<VectorizedArray<Number>::size(); ++v)
   * this->operator[](v) = base_ptr[offsets[v]];
   * @endcode
   *
   */
  DEAL_II_ALWAYS_INLINE
  void
  gather(const double *base_ptr, const unsigned int *offsets)
  {
    for (unsigned int i = 0; i < 2; ++i)
      *(reinterpret_cast<double *>(&data) + i) = base_ptr[offsets[i]];
  }

  /**
   * 将调用类的内容以 @p
   * size()的形式写入内存，到给定的地址和给定的偏移量，将矢量数组的元素填入每个偏移量。
   * 这个操作对应于下面的代码（但在硬件允许的情况下使用了更有效的实现）。
   * @code
   * for (unsigned int v=0; v<VectorizedArray<Number>::size(); ++v)
   * base_ptr[offsets[v]] = this->operator[](v);
   * @endcode
   *
   */
  DEAL_II_ALWAYS_INLINE
  void
  scatter(const unsigned int *offsets, double *base_ptr) const
  {
    for (unsigned int i = 0; i < 2; ++i)
      base_ptr[offsets[i]] = *(reinterpret_cast<const double *>(&data) + i);
  }

  /**
   * 实际的数据字段。为了与标准布局类型保持一致，并且能够与外部SIMD功能进行交互，这个成员被声明为公共的。
   *
   */
  __m128d data;

private:
  /**
   * 返回这个字段的平方根。不适合在用户代码中使用。使用sqrt(x)代替。
   *
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray
  get_sqrt() const
  {
    VectorizedArray res;
    res.data = _mm_sqrt_pd(data);
    return res;
  }

  /**
   * 返回这个字段的绝对值。不适合在用户代码中使用。请使用abs(x)代替。
   *
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray
  get_abs() const
  {
    // to compute the absolute value, perform
    // bitwise andnot with -0. This will leave all
    // value and exponent bits unchanged but force
    // the sign value to +.
    __m128d         mask = _mm_set1_pd(-0.);
    VectorizedArray res;
    res.data = _mm_andnot_pd(mask, data);
    return res;
  }

  /**
   * 返回这个字段和另一个字段的分量上的最大值。不适合在用户代码中使用。使用max(x,y)代替。
   *
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray
  get_max(const VectorizedArray &other) const
  {
    VectorizedArray res;
    res.data = _mm_max_pd(data, other.data);
    return res;
  }

  /**
   * 返回这个字段和另一个字段的最小分量。不适合在用户代码中使用。使用min(x,y)代替。
   *
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray
  get_min(const VectorizedArray &other) const
  {
    VectorizedArray res;
    res.data = _mm_min_pd(data, other.data);
    return res;
  }

  // Make a few functions friends.
  template <typename Number2, std::size_t width2>
  friend VectorizedArray<Number2, width2>
  std::sqrt(const VectorizedArray<Number2, width2> &);
  template <typename Number2, std::size_t width2>
  friend VectorizedArray<Number2, width2>
  std::abs(const VectorizedArray<Number2, width2> &);
  template <typename Number2, std::size_t width2>
  friend VectorizedArray<Number2, width2>
  std::max(const VectorizedArray<Number2, width2> &,
           const VectorizedArray<Number2, width2> &);
  template <typename Number2, std::size_t width2>
  friend VectorizedArray<Number2, width2>
  std::min(const VectorizedArray<Number2, width2> &,
           const VectorizedArray<Number2, width2> &);
};



/**
 * double和SSE2的特殊化。
 *
 *
 */
template <>
inline DEAL_II_ALWAYS_INLINE void
vectorized_load_and_transpose(const unsigned int          n_entries,
                              const double *              in,
                              const unsigned int *        offsets,
                              VectorizedArray<double, 2> *out)
{
  const unsigned int n_chunks = n_entries / 2;
  for (unsigned int i = 0; i < n_chunks; ++i)
    {
      __m128d u0          = _mm_loadu_pd(in + 2 * i + offsets[0]);
      __m128d u1          = _mm_loadu_pd(in + 2 * i + offsets[1]);
      out[2 * i + 0].data = _mm_unpacklo_pd(u0, u1);
      out[2 * i + 1].data = _mm_unpackhi_pd(u0, u1);
    }

  // remainder loop of work that does not divide by 2
  for (unsigned int i = 2 * n_chunks; i < n_entries; ++i)
    for (unsigned int v = 0; v < 2; ++v)
      out[i][v] = in[offsets[v] + i];
}



/**
 * 对double和SSE2的特殊化。
 *
 *
 */
template <>
inline DEAL_II_ALWAYS_INLINE void
vectorized_load_and_transpose(const unsigned int             n_entries,
                              const std::array<double *, 2> &in,
                              VectorizedArray<double, 2> *   out)
{
  // see the comments in the vectorized_load_and_transpose above

  const unsigned int n_chunks = n_entries / 2;
  for (unsigned int i = 0; i < n_chunks; ++i)
    {
      __m128d u0          = _mm_loadu_pd(in[0] + 2 * i);
      __m128d u1          = _mm_loadu_pd(in[1] + 2 * i);
      out[2 * i + 0].data = _mm_unpacklo_pd(u0, u1);
      out[2 * i + 1].data = _mm_unpackhi_pd(u0, u1);
    }

  for (unsigned int i = 2 * n_chunks; i < n_entries; ++i)
    for (unsigned int v = 0; v < 2; ++v)
      out[i][v] = in[v][i];
}



/**
 * 双重和SSE2的特殊化。
 *
 *
 */
template <>
inline DEAL_II_ALWAYS_INLINE void
vectorized_transpose_and_store(const bool                        add_into,
                               const unsigned int                n_entries,
                               const VectorizedArray<double, 2> *in,
                               const unsigned int *              offsets,
                               double *                          out)
{
  const unsigned int n_chunks = n_entries / 2;
  if (add_into)
    {
      for (unsigned int i = 0; i < n_chunks; ++i)
        {
          __m128d u0   = in[2 * i + 0].data;
          __m128d u1   = in[2 * i + 1].data;
          __m128d res0 = _mm_unpacklo_pd(u0, u1);
          __m128d res1 = _mm_unpackhi_pd(u0, u1);
          _mm_storeu_pd(out + 2 * i + offsets[0],
                        _mm_add_pd(_mm_loadu_pd(out + 2 * i + offsets[0]),
                                   res0));
          _mm_storeu_pd(out + 2 * i + offsets[1],
                        _mm_add_pd(_mm_loadu_pd(out + 2 * i + offsets[1]),
                                   res1));
        }
      // remainder loop of work that does not divide by 2
      for (unsigned int i = 2 * n_chunks; i < n_entries; ++i)
        for (unsigned int v = 0; v < 2; ++v)
          out[offsets[v] + i] += in[i][v];
    }
  else
    {
      for (unsigned int i = 0; i < n_chunks; ++i)
        {
          __m128d u0   = in[2 * i + 0].data;
          __m128d u1   = in[2 * i + 1].data;
          __m128d res0 = _mm_unpacklo_pd(u0, u1);
          __m128d res1 = _mm_unpackhi_pd(u0, u1);
          _mm_storeu_pd(out + 2 * i + offsets[0], res0);
          _mm_storeu_pd(out + 2 * i + offsets[1], res1);
        }
      // remainder loop of work that does not divide by 2
      for (unsigned int i = 2 * n_chunks; i < n_entries; ++i)
        for (unsigned int v = 0; v < 2; ++v)
          out[offsets[v] + i] = in[i][v];
    }
}



/**
 * 双重和SSE2的特殊化。
 *
 *
 */
template <>
inline DEAL_II_ALWAYS_INLINE void
vectorized_transpose_and_store(const bool                        add_into,
                               const unsigned int                n_entries,
                               const VectorizedArray<double, 2> *in,
                               std::array<double *, 2> &         out)
{
  // see the comments in the vectorized_transpose_and_store above

  const unsigned int n_chunks = n_entries / 2;
  if (add_into)
    {
      for (unsigned int i = 0; i < n_chunks; ++i)
        {
          __m128d u0   = in[2 * i + 0].data;
          __m128d u1   = in[2 * i + 1].data;
          __m128d res0 = _mm_unpacklo_pd(u0, u1);
          __m128d res1 = _mm_unpackhi_pd(u0, u1);
          _mm_storeu_pd(out[0] + 2 * i,
                        _mm_add_pd(_mm_loadu_pd(out[0] + 2 * i), res0));
          _mm_storeu_pd(out[1] + 2 * i,
                        _mm_add_pd(_mm_loadu_pd(out[1] + 2 * i), res1));
        }

      for (unsigned int i = 2 * n_chunks; i < n_entries; ++i)
        for (unsigned int v = 0; v < 2; ++v)
          out[v][i] += in[i][v];
    }
  else
    {
      for (unsigned int i = 0; i < n_chunks; ++i)
        {
          __m128d u0   = in[2 * i + 0].data;
          __m128d u1   = in[2 * i + 1].data;
          __m128d res0 = _mm_unpacklo_pd(u0, u1);
          __m128d res1 = _mm_unpackhi_pd(u0, u1);
          _mm_storeu_pd(out[0] + 2 * i, res0);
          _mm_storeu_pd(out[1] + 2 * i, res1);
        }

      for (unsigned int i = 2 * n_chunks; i < n_entries; ++i)
        for (unsigned int v = 0; v < 2; ++v)
          out[v][i] = in[i][v];
    }
}



/**
 * 浮动和SSE2的特殊化。
 *
 *
 */
template <>
class VectorizedArray<float, 4>
  : public VectorizedArrayBase<VectorizedArray<float, 4>, 4>
{
public:
  /**
   * 这给出了数组元素的类型。
   *
   */
  using value_type = float;

  /**
   * 这个函数可以用来将所有的数据字段设置为一个给定的标量。
   *
   */

  /**
   * 默认的空构造函数，让数据处于未初始化的状态，类似于float/double。
   *
   */
  VectorizedArray() = default;

  /**
   * 构造一个数组，将给定的标量广播给所有通道。
   *
   */
  VectorizedArray(const float scalar)
  {
    this->operator=(scalar);
  }

  DEAL_II_ALWAYS_INLINE
  VectorizedArray &
  operator=(const float x)
  {
    data = _mm_set1_ps(x);
    return *this;
  }

  /**
   * 访问操作符。
   *
   */
  DEAL_II_ALWAYS_INLINE
  float &operator[](const unsigned int comp)
  {
    AssertIndexRange(comp, 4);
    return *(reinterpret_cast<float *>(&data) + comp);
  }

  /**
   * 常数访问操作符。
   *
   */
  DEAL_II_ALWAYS_INLINE
  const float &operator[](const unsigned int comp) const
  {
    AssertIndexRange(comp, 4);
    return *(reinterpret_cast<const float *>(&data) + comp);
  }

  /**
   * 加法。
   *
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray &
  operator+=(const VectorizedArray &vec)
  {
#    ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data += vec.data;
#    else
    data = _mm_add_ps(data, vec.data);
#    endif
    return *this;
  }

  /**
   * 减法。
   *
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray &
  operator-=(const VectorizedArray &vec)
  {
#    ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data -= vec.data;
#    else
    data = _mm_sub_ps(data, vec.data);
#    endif
    return *this;
  }

  /**
   * 乘法。
   *
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray &
  operator*=(const VectorizedArray &vec)
  {
#    ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data *= vec.data;
#    else
    data = _mm_mul_ps(data, vec.data);
#    endif
    return *this;
  }

  /**
   * 除法。
   *
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray &
  operator/=(const VectorizedArray &vec)
  {
#    ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data /= vec.data;
#    else
    data = _mm_div_ps(data, vec.data);
#    endif
    return *this;
  }

  /**
   * 从内存中加载 @p size()
   * 到调用类中，从给定的地址开始。内存不需要对齐16个字节，相对于将浮点地址投给VectorizedArray<float>*来说。
   *
   */
  DEAL_II_ALWAYS_INLINE
  void
  load(const float *ptr)
  {
    data = _mm_loadu_ps(ptr);
  }

  /**
   * 将调用类的内容以 @p
   * size()的形式写入内存到给定地址。内存不需要对齐16个字节，与之相对应的是将浮点地址投给VectorizedArray<float>*。
   *
   */
  DEAL_II_ALWAYS_INLINE
  void
  store(float *ptr) const
  {
    _mm_storeu_ps(ptr, data);
  }

  /**
   * @copydoc   VectorizedArray<Number>::streaming_store() 。
   * @note  内存必须以16字节对齐。
   *
   */
  DEAL_II_ALWAYS_INLINE
  void
  streaming_store(float *ptr) const
  {
    Assert(reinterpret_cast<std::size_t>(ptr) % 16 == 0,
           ExcMessage("Memory not aligned"));
    _mm_stream_ps(ptr, data);
  }

  /**
   * 将 @p size()
   * 从内存中加载到调用的类中，从给定的地址和给定的偏移量开始，从偏移量开始的每个条目提供一个矢量数组的元素。
   * 这个操作对应于以下代码（但在硬件允许的情况下，使用了更有效的实现）。
   * @code
   * for (unsigned int v=0; v<VectorizedArray<Number>::size(); ++v)
   * this->operator[](v) = base_ptr[offsets[v]];
   * @endcode
   *
   */
  DEAL_II_ALWAYS_INLINE
  void
  gather(const float *base_ptr, const unsigned int *offsets)
  {
    for (unsigned int i = 0; i < 4; ++i)
      *(reinterpret_cast<float *>(&data) + i) = base_ptr[offsets[i]];
  }

  /**
   * 将调用类的内容以 @p
   * size()的形式写入内存，到给定的地址和给定的偏移量，将矢量数组的元素填入每个偏移量。
   * 这个操作对应于下面的代码（但在硬件允许的情况下使用了更有效的实现）。
   * @code
   * for (unsigned int v=0; v<VectorizedArray<Number>::size(); ++v)
   * base_ptr[offsets[v]] = this->operator[](v);
   * @endcode
   *
   *
   */
  DEAL_II_ALWAYS_INLINE
  void
  scatter(const unsigned int *offsets, float *base_ptr) const
  {
    for (unsigned int i = 0; i < 4; ++i)
      base_ptr[offsets[i]] = *(reinterpret_cast<const float *>(&data) + i);
  }

  /**
   * 实际的数据字段。为了与标准布局类型保持一致，并且能够与外部SIMD功能进行交互，这个成员被声明为公共的。
   *
   */
  __m128 data;

private:
  /**
   * 返回这个字段的平方根。不适合在用户代码中使用。使用sqrt(x)代替。
   *
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray
  get_sqrt() const
  {
    VectorizedArray res;
    res.data = _mm_sqrt_ps(data);
    return res;
  }

  /**
   * 返回该字段的绝对值。不适合在用户代码中使用。请使用abs(x)代替。
   *
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray
  get_abs() const
  {
    // to compute the absolute value, perform bitwise andnot with -0. This
    // will leave all value and exponent bits unchanged but force the sign
    // value to +.
    __m128          mask = _mm_set1_ps(-0.f);
    VectorizedArray res;
    res.data = _mm_andnot_ps(mask, data);
    return res;
  }

  /**
   * 返回这个字段和另一个字段的分量上的最大值。不适合在用户代码中使用。使用max(x,y)代替。
   *
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray
  get_max(const VectorizedArray &other) const
  {
    VectorizedArray res;
    res.data = _mm_max_ps(data, other.data);
    return res;
  }

  /**
   * 返回这个字段和另一个字段的最小分量。不适合在用户代码中使用。使用min(x,y)代替。
   *
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray
  get_min(const VectorizedArray &other) const
  {
    VectorizedArray res;
    res.data = _mm_min_ps(data, other.data);
    return res;
  }

  // Make a few functions friends.
  template <typename Number2, std::size_t width2>
  friend VectorizedArray<Number2, width2>
  std::sqrt(const VectorizedArray<Number2, width2> &);
  template <typename Number2, std::size_t width2>
  friend VectorizedArray<Number2, width2>
  std::abs(const VectorizedArray<Number2, width2> &);
  template <typename Number2, std::size_t width2>
  friend VectorizedArray<Number2, width2>
  std::max(const VectorizedArray<Number2, width2> &,
           const VectorizedArray<Number2, width2> &);
  template <typename Number2, std::size_t width2>
  friend VectorizedArray<Number2, width2>
  std::min(const VectorizedArray<Number2, width2> &,
           const VectorizedArray<Number2, width2> &);
};



/**
 * 对float和SSE2的特殊化。
 *
 *
 */
template <>
inline DEAL_II_ALWAYS_INLINE void
vectorized_load_and_transpose(const unsigned int         n_entries,
                              const float *              in,
                              const unsigned int *       offsets,
                              VectorizedArray<float, 4> *out)
{
  const unsigned int n_chunks = n_entries / 4;
  for (unsigned int i = 0; i < n_chunks; ++i)
    {
      __m128 u0           = _mm_loadu_ps(in + 4 * i + offsets[0]);
      __m128 u1           = _mm_loadu_ps(in + 4 * i + offsets[1]);
      __m128 u2           = _mm_loadu_ps(in + 4 * i + offsets[2]);
      __m128 u3           = _mm_loadu_ps(in + 4 * i + offsets[3]);
      __m128 v0           = _mm_shuffle_ps(u0, u1, 0x44);
      __m128 v1           = _mm_shuffle_ps(u0, u1, 0xee);
      __m128 v2           = _mm_shuffle_ps(u2, u3, 0x44);
      __m128 v3           = _mm_shuffle_ps(u2, u3, 0xee);
      out[4 * i + 0].data = _mm_shuffle_ps(v0, v2, 0x88);
      out[4 * i + 1].data = _mm_shuffle_ps(v0, v2, 0xdd);
      out[4 * i + 2].data = _mm_shuffle_ps(v1, v3, 0x88);
      out[4 * i + 3].data = _mm_shuffle_ps(v1, v3, 0xdd);
    }

  // remainder loop of work that does not divide by 4
  for (unsigned int i = 4 * n_chunks; i < n_entries; ++i)
    for (unsigned int v = 0; v < 4; ++v)
      out[i][v] = in[offsets[v] + i];
}



/**
 * 对浮点和SSE2的特殊化。
 *
 *
 */
template <>
inline DEAL_II_ALWAYS_INLINE void
vectorized_load_and_transpose(const unsigned int            n_entries,
                              const std::array<float *, 4> &in,
                              VectorizedArray<float, 4> *   out)
{
  // see the comments in the vectorized_load_and_transpose above

  const unsigned int n_chunks = n_entries / 4;
  for (unsigned int i = 0; i < n_chunks; ++i)
    {
      __m128 u0           = _mm_loadu_ps(in[0] + 4 * i);
      __m128 u1           = _mm_loadu_ps(in[1] + 4 * i);
      __m128 u2           = _mm_loadu_ps(in[2] + 4 * i);
      __m128 u3           = _mm_loadu_ps(in[3] + 4 * i);
      __m128 v0           = _mm_shuffle_ps(u0, u1, 0x44);
      __m128 v1           = _mm_shuffle_ps(u0, u1, 0xee);
      __m128 v2           = _mm_shuffle_ps(u2, u3, 0x44);
      __m128 v3           = _mm_shuffle_ps(u2, u3, 0xee);
      out[4 * i + 0].data = _mm_shuffle_ps(v0, v2, 0x88);
      out[4 * i + 1].data = _mm_shuffle_ps(v0, v2, 0xdd);
      out[4 * i + 2].data = _mm_shuffle_ps(v1, v3, 0x88);
      out[4 * i + 3].data = _mm_shuffle_ps(v1, v3, 0xdd);
    }

  for (unsigned int i = 4 * n_chunks; i < n_entries; ++i)
    for (unsigned int v = 0; v < 4; ++v)
      out[i][v] = in[v][i];
}



/**
 * 浮动和SSE2的特殊化。
 *
 *
 */
template <>
inline DEAL_II_ALWAYS_INLINE void
vectorized_transpose_and_store(const bool                       add_into,
                               const unsigned int               n_entries,
                               const VectorizedArray<float, 4> *in,
                               const unsigned int *             offsets,
                               float *                          out)
{
  const unsigned int n_chunks = n_entries / 4;
  for (unsigned int i = 0; i < n_chunks; ++i)
    {
      __m128 u0 = in[4 * i + 0].data;
      __m128 u1 = in[4 * i + 1].data;
      __m128 u2 = in[4 * i + 2].data;
      __m128 u3 = in[4 * i + 3].data;
      __m128 t0 = _mm_shuffle_ps(u0, u1, 0x44);
      __m128 t1 = _mm_shuffle_ps(u0, u1, 0xee);
      __m128 t2 = _mm_shuffle_ps(u2, u3, 0x44);
      __m128 t3 = _mm_shuffle_ps(u2, u3, 0xee);
      u0        = _mm_shuffle_ps(t0, t2, 0x88);
      u1        = _mm_shuffle_ps(t0, t2, 0xdd);
      u2        = _mm_shuffle_ps(t1, t3, 0x88);
      u3        = _mm_shuffle_ps(t1, t3, 0xdd);

      // Cannot use the same store instructions in both paths of the 'if'
      // because the compiler cannot know that there is no aliasing between
      // pointers
      if (add_into)
        {
          u0 = _mm_add_ps(_mm_loadu_ps(out + 4 * i + offsets[0]), u0);
          _mm_storeu_ps(out + 4 * i + offsets[0], u0);
          u1 = _mm_add_ps(_mm_loadu_ps(out + 4 * i + offsets[1]), u1);
          _mm_storeu_ps(out + 4 * i + offsets[1], u1);
          u2 = _mm_add_ps(_mm_loadu_ps(out + 4 * i + offsets[2]), u2);
          _mm_storeu_ps(out + 4 * i + offsets[2], u2);
          u3 = _mm_add_ps(_mm_loadu_ps(out + 4 * i + offsets[3]), u3);
          _mm_storeu_ps(out + 4 * i + offsets[3], u3);
        }
      else
        {
          _mm_storeu_ps(out + 4 * i + offsets[0], u0);
          _mm_storeu_ps(out + 4 * i + offsets[1], u1);
          _mm_storeu_ps(out + 4 * i + offsets[2], u2);
          _mm_storeu_ps(out + 4 * i + offsets[3], u3);
        }
    }

  // remainder loop of work that does not divide by 4
  if (add_into)
    for (unsigned int i = 4 * n_chunks; i < n_entries; ++i)
      for (unsigned int v = 0; v < 4; ++v)
        out[offsets[v] + i] += in[i][v];
  else
    for (unsigned int i = 4 * n_chunks; i < n_entries; ++i)
      for (unsigned int v = 0; v < 4; ++v)
        out[offsets[v] + i] = in[i][v];
}



/**
 * 浮动和SSE2的特殊化。
 *
 *
 */
template <>
inline DEAL_II_ALWAYS_INLINE void
vectorized_transpose_and_store(const bool                       add_into,
                               const unsigned int               n_entries,
                               const VectorizedArray<float, 4> *in,
                               std::array<float *, 4> &         out)
{
  // see the comments in the vectorized_transpose_and_store above

  const unsigned int n_chunks = n_entries / 4;
  for (unsigned int i = 0; i < n_chunks; ++i)
    {
      __m128 u0 = in[4 * i + 0].data;
      __m128 u1 = in[4 * i + 1].data;
      __m128 u2 = in[4 * i + 2].data;
      __m128 u3 = in[4 * i + 3].data;
      __m128 t0 = _mm_shuffle_ps(u0, u1, 0x44);
      __m128 t1 = _mm_shuffle_ps(u0, u1, 0xee);
      __m128 t2 = _mm_shuffle_ps(u2, u3, 0x44);
      __m128 t3 = _mm_shuffle_ps(u2, u3, 0xee);
      u0        = _mm_shuffle_ps(t0, t2, 0x88);
      u1        = _mm_shuffle_ps(t0, t2, 0xdd);
      u2        = _mm_shuffle_ps(t1, t3, 0x88);
      u3        = _mm_shuffle_ps(t1, t3, 0xdd);

      if (add_into)
        {
          u0 = _mm_add_ps(_mm_loadu_ps(out[0] + 4 * i), u0);
          _mm_storeu_ps(out[0] + 4 * i, u0);
          u1 = _mm_add_ps(_mm_loadu_ps(out[1] + 4 * i), u1);
          _mm_storeu_ps(out[1] + 4 * i, u1);
          u2 = _mm_add_ps(_mm_loadu_ps(out[2] + 4 * i), u2);
          _mm_storeu_ps(out[2] + 4 * i, u2);
          u3 = _mm_add_ps(_mm_loadu_ps(out[3] + 4 * i), u3);
          _mm_storeu_ps(out[3] + 4 * i, u3);
        }
      else
        {
          _mm_storeu_ps(out[0] + 4 * i, u0);
          _mm_storeu_ps(out[1] + 4 * i, u1);
          _mm_storeu_ps(out[2] + 4 * i, u2);
          _mm_storeu_ps(out[3] + 4 * i, u3);
        }
    }

  if (add_into)
    for (unsigned int i = 4 * n_chunks; i < n_entries; ++i)
      for (unsigned int v = 0; v < 4; ++v)
        out[v][i] += in[i][v];
  else
    for (unsigned int i = 4 * n_chunks; i < n_entries; ++i)
      for (unsigned int v = 0; v < 4; ++v)
        out[v][i] = in[i][v];
}



#  endif // if DEAL_II_VECTORIZATION_WIDTH_IN_BITS > 0 && defined(__SSE2__)

#  if DEAL_II_VECTORIZATION_WIDTH_IN_BITS >= 128 && defined(__ALTIVEC__) && \
    defined(__VSX__)

template <>
class VectorizedArray<double, 2>
  : public VectorizedArrayBase<VectorizedArray<double, 2>, 2>
{
public:
  /**
   * 这给出了数组元素的类型。
   *
   */
  using value_type = double;

  /**
   * 默认的空构造函数，让数据处于未初始化的状态，类似于float/double。
   *
   */
  VectorizedArray() = default;

  /**
   * 用给定的标量构建一个数组广播到所有的通道。
   *
   */
  VectorizedArray(const double scalar)
  {
    this->operator=(scalar);
  }

  /**
   * 这个函数将一个标量分配给这个类。
   *
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray &
  operator=(const double x)
  {
    data = vec_splats(x);

    // Some compilers believe that vec_splats sets 'x', but that's not true.
    // They then warn about setting a variable and not using it. Suppress the
    // warning by "using" the variable:
    (void)x;
    return *this;
  }

  /**
   * 访问操作符。该组件必须是0或1。
   *
   */
  DEAL_II_ALWAYS_INLINE
  double &operator[](const unsigned int comp)
  {
    AssertIndexRange(comp, 2);
    return *(reinterpret_cast<double *>(&data) + comp);
  }

  /**
   * 常数访问运算符。
   *
   */
  DEAL_II_ALWAYS_INLINE
  const double &operator[](const unsigned int comp) const
  {
    AssertIndexRange(comp, 2);
    return *(reinterpret_cast<const double *>(&data) + comp);
  }

  /**
   * 加法。
   *
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray &
  operator+=(const VectorizedArray &vec)
  {
    data = vec_add(data, vec.data);
    return *this;
  }

  /**
   * 减法。
   *
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray &
  operator-=(const VectorizedArray &vec)
  {
    data = vec_sub(data, vec.data);
    return *this;
  }

  /**
   * 乘法。
   *
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray &
  operator*=(const VectorizedArray &vec)
  {
    data = vec_mul(data, vec.data);
    return *this;
  }

  /**
   * 除法。
   *
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray &
  operator/=(const VectorizedArray &vec)
  {
    data = vec_div(data, vec.data);
    return *this;
  }

  /**
   * 从内存中加载 @p size()
   * 到调用类中，从给定的地址开始。
   *
   */
  DEAL_II_ALWAYS_INLINE
  void
  load(const double *ptr)
  {
    data = vec_vsx_ld(0, ptr);
  }

  /**
   * 将调用类的内容以 @p
   * size()的形式写进内存，到给定的地址。
   *
   */
  DEAL_II_ALWAYS_INLINE
  void
  store(double *ptr) const
  {
    vec_vsx_st(data, 0, ptr);
  }

   /**
    * @copydoc VectorizedArray<Number>::streaming_store()
    */
  DEAL_II_ALWAYS_INLINE
  void
  streaming_store(double *ptr) const
  {
    store(ptr);
  }

   /**
    * @copydoc VectorizedArray<Number>::gather()
    */
  DEAL_II_ALWAYS_INLINE
  void
  gather(const double *base_ptr, const unsigned int *offsets)
  {
    for (unsigned int i = 0; i < 2; ++i)
      *(reinterpret_cast<double *>(&data) + i) = base_ptr[offsets[i]];
  }

   /**
    * @copydoc VectorizedArray<Number>::scatter
    */
  DEAL_II_ALWAYS_INLINE
  void
  scatter(const unsigned int *offsets, double *base_ptr) const
  {
    for (unsigned int i = 0; i < 2; ++i)
      base_ptr[offsets[i]] = *(reinterpret_cast<const double *>(&data) + i);
  }

  /**
   * 实际的数据字段。为了与标准布局类型保持一致，并且能够与外部SIMD功能进行交互，这个成员被声明为公共的。
   *
   */
  __vector double data;

private:
  /**
   * 返回这个字段的平方根。不适合在用户代码中使用。使用sqrt(x)代替。
   *
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray
  get_sqrt() const
  {
    VectorizedArray res;
    res.data = vec_sqrt(data);
    return res;
  }

  /**
   * 返回该字段的绝对值。不适合在用户代码中使用。请使用abs(x)代替。
   *
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray
  get_abs() const
  {
    VectorizedArray res;
    res.data = vec_abs(data);
    return res;
  }

  /**
   * 返回这个字段和另一个字段的分量上的最大值。不适合在用户代码中使用。使用max(x,y)代替。
   *
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray
  get_max(const VectorizedArray &other) const
  {
    VectorizedArray res;
    res.data = vec_max(data, other.data);
    return res;
  }

  /**
   * 返回这个字段和另一个字段的最小分量。不适合在用户代码中使用。使用min(x,y)代替。
   *
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray
  get_min(const VectorizedArray &other) const
  {
    VectorizedArray res;
    res.data = vec_min(data, other.data);
    return res;
  }

  // Make a few functions friends.
  template <typename Number2, std::size_t width2>
  friend VectorizedArray<Number2, width2>
  std::sqrt(const VectorizedArray<Number2, width2> &);
  template <typename Number2, std::size_t width2>
  friend VectorizedArray<Number2, width2>
  std::abs(const VectorizedArray<Number2, width2> &);
  template <typename Number2, std::size_t width2>
  friend VectorizedArray<Number2, width2>
  std::max(const VectorizedArray<Number2, width2> &,
           const VectorizedArray<Number2, width2> &);
  template <typename Number2, std::size_t width2>
  friend VectorizedArray<Number2, width2>
  std::min(const VectorizedArray<Number2, width2> &,
           const VectorizedArray<Number2, width2> &);
};



template <>
class VectorizedArray<float, 4>
  : public VectorizedArrayBase<VectorizedArray<float, 4>, 4>
{
public:
  /**
   * 这给出了数组元素的类型。
   *
   */
  using value_type = float;

  /**
   * 默认的空构造函数，让数据处于未初始化的状态，类似于float/double。
   *
   */
  VectorizedArray() = default;

  /**
   * 用给定的标量构建一个数组，广播给所有通道。
   *
   */
  VectorizedArray(const float scalar)
  {
    this->operator=(scalar);
  }

  /**
   * 这个函数将一个标量分配给这个类。
   *
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray &
  operator=(const float x)
  {
    data = vec_splats(x);

    // Some compilers believe that vec_splats sets 'x', but that's not true.
    // They then warn about setting a variable and not using it. Suppress the
    // warning by "using" the variable:
    (void)x;
    return *this;
  }

  /**
   * 访问操作符。该分量必须在0和3之间。
   *
   */
  DEAL_II_ALWAYS_INLINE
  float &operator[](const unsigned int comp)
  {
    AssertIndexRange(comp, 4);
    return *(reinterpret_cast<float *>(&data) + comp);
  }

  /**
   * 常数访问运算符。
   *
   */
  DEAL_II_ALWAYS_INLINE
  const float &operator[](const unsigned int comp) const
  {
    AssertIndexRange(comp, 4);
    return *(reinterpret_cast<const float *>(&data) + comp);
  }

  /**
   * 加法。
   *
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray &
  operator+=(const VectorizedArray &vec)
  {
    data = vec_add(data, vec.data);
    return *this;
  }

  /**
   * 减法。
   *
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray &
  operator-=(const VectorizedArray &vec)
  {
    data = vec_sub(data, vec.data);
    return *this;
  }

  /**
   * 乘法。
   *
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray &
  operator*=(const VectorizedArray &vec)
  {
    data = vec_mul(data, vec.data);
    return *this;
  }

  /**
   * 除法。
   *
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray &
  operator/=(const VectorizedArray &vec)
  {
    data = vec_div(data, vec.data);
    return *this;
  }

  /**
   * 从内存中加载 @p size()
   * 到调用类中，从给定的地址开始。
   *
   */
  DEAL_II_ALWAYS_INLINE
  void
  load(const float *ptr)
  {
    data = vec_vsx_ld(0, ptr);
  }

  /**
   * 将调用类的内容以 @p
   * size()的形式写进内存，到给定的地址。
   *
   */
  DEAL_II_ALWAYS_INLINE
  void
  store(float *ptr) const
  {
    vec_vsx_st(data, 0, ptr);
  }

   /**
    * @copydoc VectorizedArray<Number>::streaming_store()
    */
  DEAL_II_ALWAYS_INLINE
  void
  streaming_store(float *ptr) const
  {
    store(ptr);
  }

   /**
    * @copydoc VectorizedArray<Number>::gather()
    */
  DEAL_II_ALWAYS_INLINE
  void
  gather(const float *base_ptr, const unsigned int *offsets)
  {
    for (unsigned int i = 0; i < 4; ++i)
      *(reinterpret_cast<float *>(&data) + i) = base_ptr[offsets[i]];
  }

   /**
    * @copydoc VectorizedArray<Number>::scatter
    */
  DEAL_II_ALWAYS_INLINE
  void
  scatter(const unsigned int *offsets, float *base_ptr) const
  {
    for (unsigned int i = 0; i < 4; ++i)
      base_ptr[offsets[i]] = *(reinterpret_cast<const float *>(&data) + i);
  }

  /**
   * 实际的数据字段。为了与标准布局类型保持一致，并且能够与外部SIMD功能进行交互，这个成员被声明为公共的。
   *
   */
  __vector float data;

private:
  /**
   * 返回这个字段的平方根。不适合在用户代码中使用。使用sqrt(x)代替。
   *
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray
  get_sqrt() const
  {
    VectorizedArray res;
    res.data = vec_sqrt(data);
    return res;
  }

  /**
   * 返回这个字段的绝对值。不适合在用户代码中使用。使用abs(x)代替。
   *
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray
  get_abs() const
  {
    VectorizedArray res;
    res.data = vec_abs(data);
    return res;
  }

  /**
   * 返回这个字段和另一个字段的分量上的最大值。不适合在用户代码中使用。使用max(x,y)代替。
   *
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray
  get_max(const VectorizedArray &other) const
  {
    VectorizedArray res;
    res.data = vec_max(data, other.data);
    return res;
  }

  /**
   * 返回这个字段和另一个字段的最小分量。不适合在用户代码中使用。使用min(x,y)代替。
   *
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray
  get_min(const VectorizedArray &other) const
  {
    VectorizedArray res;
    res.data = vec_min(data, other.data);
    return res;
  }

  // Make a few functions friends.
  template <typename Number2, std::size_t width2>
  friend VectorizedArray<Number2, width2>
  std::sqrt(const VectorizedArray<Number2, width2> &);
  template <typename Number2, std::size_t width2>
  friend VectorizedArray<Number2, width2>
  std::abs(const VectorizedArray<Number2, width2> &);
  template <typename Number2, std::size_t width2>
  friend VectorizedArray<Number2, width2>
  std::max(const VectorizedArray<Number2, width2> &,
           const VectorizedArray<Number2, width2> &);
  template <typename Number2, std::size_t width2>
  friend VectorizedArray<Number2, width2>
  std::min(const VectorizedArray<Number2, width2> &,
           const VectorizedArray<Number2, width2> &);
};

#  endif // if DEAL_II_VECTORIZATION_LEVEL >=1 && defined(__ALTIVEC__) &&
         // defined(__VSX__)


#endif // DOXYGEN

/**
 * @name  用VectorizedArray进行算术运算
 *
 *
 */
//@{

/**
 * 对VectorizedArray的关系运算符==
 * 
 * @relatesalso  VectorizedArray 
 *
 *
 */
template <typename Number, std::size_t width>
inline DEAL_II_ALWAYS_INLINE bool
operator==(const VectorizedArray<Number, width> &lhs,
           const VectorizedArray<Number, width> &rhs)
{
  for (unsigned int i = 0; i < VectorizedArray<Number, width>::size(); ++i)
    if (lhs[i] != rhs[i])
      return false;

  return true;
}


/**
 * 用操作符+对两个矢量数组进行加法。
 * @relatesalso  VectorizedArray
 *
 *
 */
template <typename Number, std::size_t width>
inline DEAL_II_ALWAYS_INLINE VectorizedArray<Number, width>
                             operator+(const VectorizedArray<Number, width> &u,
          const VectorizedArray<Number, width> &v)
{
  VectorizedArray<Number, width> tmp = u;
  return tmp += v;
}

/**
 * 两个矢量数组的减法，用操作符
 *
 * -
 * @relatesalso VectorizedArray
 *
 *
 */
template <typename Number, std::size_t width>
inline DEAL_II_ALWAYS_INLINE VectorizedArray<Number, width>
                             operator-(const VectorizedArray<Number, width> &u,
          const VectorizedArray<Number, width> &v)
{
  VectorizedArray<Number, width> tmp = u;
  return tmp -= v;
}

/**
 * 用运算符对两个矢量数组进行乘法。
 * @relatesalso VectorizedArray
 *
 *
 */
template <typename Number, std::size_t width>
inline DEAL_II_ALWAYS_INLINE VectorizedArray<Number, width>
                             operator*(const VectorizedArray<Number, width> &u,
          const VectorizedArray<Number, width> &v)
{
  VectorizedArray<Number, width> tmp = u;
  return tmp *= v;
}

/**
 * 用运算符/对两个矢量数组进行除法。
 * @relatesalso  VectorizedArray
 *
 *
 */
template <typename Number, std::size_t width>
inline DEAL_II_ALWAYS_INLINE VectorizedArray<Number, width>
                             operator/(const VectorizedArray<Number, width> &u,
          const VectorizedArray<Number, width> &v)
{
  VectorizedArray<Number, width> tmp = u;
  return tmp /= v;
}

/**
 * 一个标量（用 @p
 * size()相等的条目扩展为一个向量数组）和一个向量数组的加法。
 * @relatesalso  VectorizedArray
 *
 *
 */
template <typename Number, std::size_t width>
inline DEAL_II_ALWAYS_INLINE VectorizedArray<Number, width>
                             operator+(const Number &u, const VectorizedArray<Number, width> &v)
{
  VectorizedArray<Number, width> tmp = u;
  return tmp += v;
}

/**
 * 在标量是双数的情况下，标量的加法（扩展为矢量数组，
 * @p
 * size()相等的条目）和矢量数组的加法（为了能够用通常是双数的常量编写简单的代码而需要）。
 * @relatesalso VectorizedArray
 *
 *
 */
template <std::size_t width>
inline DEAL_II_ALWAYS_INLINE VectorizedArray<float, width>
                             operator+(const double u, const VectorizedArray<float, width> &v)
{
  VectorizedArray<float, width> tmp = u;
  return tmp += v;
}

/**
 * 矢量数组和标量的相加（扩展为一个具有 @p size()
 * 等量项的矢量数组）。
 * @relatesalso  VectorizedArray
 *
 *
 */
template <typename Number, std::size_t width>
inline DEAL_II_ALWAYS_INLINE VectorizedArray<Number, width>
                             operator+(const VectorizedArray<Number, width> &v, const Number &u)
{
  return u + v;
}

/**
 * 在标量为双数的情况下，将一个矢量数组和一个标量相加（扩展为一个具有
 * @p size()
 * 相等条目的矢量数组）（为了能够用通常为双数的常量编写简单的代码而需要）。
 * @relatesalso VectorizedArray
 *
 *
 */
template <std::size_t width>
inline DEAL_II_ALWAYS_INLINE VectorizedArray<float, width>
                             operator+(const VectorizedArray<float, width> &v, const double u)
{
  return u + v;
}

/**
 * 从标量中减去一个矢量数组（扩展为一个具有 @p size()
 * 相等条目的矢量数组）。
 * @relatesalso  VectorizedArray
 *
 *
 */
template <typename Number, std::size_t width>
inline DEAL_II_ALWAYS_INLINE VectorizedArray<Number, width>
                             operator-(const Number &u, const VectorizedArray<Number, width> &v)
{
  VectorizedArray<Number, width> tmp = u;
  return tmp -= v;
}

/**
 * 在标量为双数的情况下，从标量中减去一个矢量数组（扩展为一个具有
 * @p size()
 * 相等条目的矢量数组）（为了能够用通常为双数的常数编写简单的代码而需要）。
 * @relatesalso VectorizedArray
 *
 *
 */
template <std::size_t width>
inline DEAL_II_ALWAYS_INLINE VectorizedArray<float, width>
                             operator-(const double u, const VectorizedArray<float, width> &v)
{
  VectorizedArray<float, width> tmp = static_cast<float>(u);
  return tmp -= v;
}

/**
 * 从一个矢量数组中减去一个标量（扩展为一个矢量数组，
 * @p  size()等于条目）。
 * @relatesalso  VectorizedArray
 *
 *
 */
template <typename Number, std::size_t width>
inline DEAL_II_ALWAYS_INLINE VectorizedArray<Number, width>
                             operator-(const VectorizedArray<Number, width> &v, const Number &u)
{
  VectorizedArray<Number, width> tmp = u;
  return v - tmp;
}

/**
 * 从一个矢量数组中减去一个标量（扩展为一个矢量数组，
 * @p
 * size()相等的条目），如果该标量是一个双数（为了能够用通常是双数的常量编写简单的代码而需要）。
 * @relatesalso VectorizedArray
 *
 *
 */
template <std::size_t width>
inline DEAL_II_ALWAYS_INLINE VectorizedArray<float, width>
                             operator-(const VectorizedArray<float, width> &v, const double u)
{
  VectorizedArray<float, width> tmp = static_cast<float>(u);
  return v - tmp;
}

/**
 * 一个标量（扩展为一个具有 @p
 * size()等分项的向量数组）和一个向量数组的乘法。
 * @relatesalso VectorizedArray
 *
 *
 */
template <typename Number, std::size_t width>
inline DEAL_II_ALWAYS_INLINE VectorizedArray<Number, width>
                             operator*(const Number &u, const VectorizedArray<Number, width> &v)
{
  VectorizedArray<Number, width> tmp = u;
  return tmp *= v;
}

/**
 * 在标量是双数的情况下，标量（扩展为一个矢量数组， @p
 * size()相等的条目）和矢量数组的乘法（为了能够用通常是双数的常量编写简单的代码而需要）。
 * @relatesalso VectorizedArray
 *
 *
 */
template <std::size_t width>
inline DEAL_II_ALWAYS_INLINE VectorizedArray<float, width>
                             operator*(const double u, const VectorizedArray<float, width> &v)
{
  VectorizedArray<float, width> tmp = static_cast<float>(u);
  return tmp *= v;
}

/**
 * 矢量数组和标量的乘法（扩展为具有 @p size()
 * 相等条目的矢量数组）。
 * @relatesalso  VectorizedArray
 *
 *
 */
template <typename Number, std::size_t width>
inline DEAL_II_ALWAYS_INLINE VectorizedArray<Number, width>
                             operator*(const VectorizedArray<Number, width> &v, const Number &u)
{
  return u * v;
}

/**
 * 在标量为双数的情况下，矢量数组和标量的乘法（扩展为具有
 * @p size()
 * 相等条目的矢量数组）（为了能够用通常为双数的常量编写简单的代码而需要）。
 * @relatesalso VectorizedArray
 *
 *
 */
template <std::size_t width>
inline DEAL_II_ALWAYS_INLINE VectorizedArray<float, width>
                             operator*(const VectorizedArray<float, width> &v, const double u)
{
  return u * v;
}

/**
 * 标量（扩展为矢量数组， @p
 * size()等于条目）与矢量数组之间的商。
 * @relatesalso VectorizedArray
 *
 *
 */
template <typename Number, std::size_t width>
inline DEAL_II_ALWAYS_INLINE VectorizedArray<Number, width>
                             operator/(const Number &u, const VectorizedArray<Number, width> &v)
{
  VectorizedArray<Number, width> tmp = u;
  return tmp /= v;
}

/**
 * 在标量是双数的情况下，标量（扩展为具有 @p
 * size()相等条目的向量数组）和向量数组之间的商（为了能够用通常为双数的常量编写简单的代码，需要）。
 * @relatesalso VectorizedArray
 *
 *
 */
template <std::size_t width>
inline DEAL_II_ALWAYS_INLINE VectorizedArray<float, width>
                             operator/(const double u, const VectorizedArray<float, width> &v)
{
  VectorizedArray<float, width> tmp = static_cast<float>(u);
  return tmp /= v;
}

/**
 * 矢量数组和标量之间的商（扩展为具有 @p size()
 * 相等条目的矢量数组）。
 * @relatesalso VectorizedArray
 *
 *
 */
template <typename Number, std::size_t width>
inline DEAL_II_ALWAYS_INLINE VectorizedArray<Number, width>
                             operator/(const VectorizedArray<Number, width> &v, const Number &u)
{
  VectorizedArray<Number, width> tmp = u;
  return v / tmp;
}

/**
 * 在标量为双数的情况下，矢量数组和标量之间的商（扩展为具有
 * @p size()
 * 相等条目的矢量数组）（为了能够用通常为双数的常量编写简单的代码而需要）。
 * @relatesalso VectorizedArray
 *
 *
 */
template <std::size_t width>
inline DEAL_II_ALWAYS_INLINE VectorizedArray<float, width>
                             operator/(const VectorizedArray<float, width> &v, const double u)
{
  VectorizedArray<float, width> tmp = static_cast<float>(u);
  return v / tmp;
}

/**
 * 矢量化数组上的单项运算符+。
 * @relatesalso VectorizedArray
 *
 *
 */
template <typename Number, std::size_t width>
inline DEAL_II_ALWAYS_INLINE VectorizedArray<Number, width>
                             operator+(const VectorizedArray<Number, width> &u)
{
  return u;
}

/**
 * 单元运算符
 *
 * - 在一个矢量化数组上。
 * @relatesalso VectorizedArray
 *
 *
 */
template <typename Number, std::size_t width>
inline DEAL_II_ALWAYS_INLINE VectorizedArray<Number, width>
                             operator-(const VectorizedArray<Number, width> &u)
{
  // to get a negative sign, subtract the input from zero (could also
  // multiply by -1, but this one is slightly simpler)
  return VectorizedArray<Number, width>() - u;
}

/**
 * 矢量化数组的输出运算符。
 * @relatesalso VectorizedArray
 *
 *
 */
template <typename Number, std::size_t width>
inline std::ostream &
operator<<(std::ostream &out, const VectorizedArray<Number, width> &p)
{
  constexpr unsigned int n = VectorizedArray<Number, width>::size();
  for (unsigned int i = 0; i < n - 1; ++i)
    out << p[i] << ' ';
  out << p[n - 1];

  return out;
}

//@}

/**
 * @name  对VectorizedArray的三元操作
 *
 *
 */
//@{


/**
 * 编码二进制操作的枚举类，用于对VectorizedArray数据类型进行组件式比较。
 *
 *
 * @note  在SIMD矢量化（sse, avx,
 * av512）的情况下，我们选择相应的有序、非信号（
 * <code>OQ</code>  ）变体。
 *
 *
 */
enum class SIMDComparison : int
{
#if DEAL_II_VECTORIZATION_WIDTH_IN_BITS >= 256 && defined(__AVX__)
  equal                 = _CMP_EQ_OQ,
  not_equal             = _CMP_NEQ_OQ,
  less_than             = _CMP_LT_OQ,
  less_than_or_equal    = _CMP_LE_OQ,
  greater_than          = _CMP_GT_OQ,
  greater_than_or_equal = _CMP_GE_OQ
#else
  equal,
  not_equal,
  less_than,
  less_than_or_equal,
  greater_than,
  greater_than_or_equal
#endif
};


/**
 * 计算以下三元操作的矢量等价物。
 *
 * @code
 * (left OP right) ? true_value : false_value
 * @endcode
 * 其中 <code>OP</code> is a binary operator (such as <code>=</code>  ,
 * <code>!=</code>, <code><</code>, <code><=</code>, <code>></code>  , 和
 * <code>>=</code>  )。
 * 当控制流本身取决于（计算的）数据时，这样的计算成语作为分支的替代是很有用的。例如，在标量数据类型的情况下，语句
 * <code>(left < right) ? true_value : false_value</code> 也可以使用
 * <code>if</code> 语句来实现。
 *
 * @code
 * if (left < right)
 *   result = true_value;
 * else
 *   result = false_value;
 * @endcode
 * 然而，这在向量化的情况下是根本不可能的，因为不同的向量条目（通道）需要不同的决定，所以必须使用第一个变体（基于三元运算符）来代替。
 *
 * @code
 * result = compare_and_apply_mask<SIMDComparison::less_than>
 *   (left, right, true_value, false_value);
 * @endcode
 * 一些更具说明性的例子（比专用的 <code>std::max</code> and
 * <code>std::abs</code> 重载效率低）。
 *
 * @code
 * VectorizedArray<double> left;
 * VectorizedArray<double> right;
 *
 * // std::max
 * const auto maximum = compare_and_apply_mask<SIMDComparison::greater_than>
 *   (left, right, left, right);
 *
 * // std::abs
 * const auto absolute = compare_and_apply_mask<SIMDComparison::less_than>
 *   (left, VectorizedArray<double>(0.),
 *
 * -left, left);
 * @endcode
 *
 * 更确切地说，这个函数首先计算一个（布尔）掩码，它是二进制运算符
 * <code>OP</code> 应用于VectorizedArray参数 @p left 和 @p right.
 * 的所有元素的结果，然后掩码被用来选择 @p true_value
 * 的相应组件（如果二进制运算相当于真），或者 @p
 * false_value. 二进制运算符通过SIMDComparison模板参数 @p
 * predicate. 编码。
 * 为了方便通用编程方法，该函数为所有VectorizedArray<Number>变体以及通用POD类型（如double和float）提供重载。
 *
 *
 * @note
 * 为了使这个函数工作，二进制操作必须通过SIMDComparison模板参数进行编码。这也解释了为什么
 * @p predicate
 * 是一个编译时常量模板参数，而不是一个常量函数参数。为了能够发出正确的低级指令，编译器必须在编译时知道比较的情况。
 *
 *
 */
template <SIMDComparison predicate, typename Number>
DEAL_II_ALWAYS_INLINE inline Number
compare_and_apply_mask(const Number &left,
                       const Number &right,
                       const Number &true_value,
                       const Number &false_value)
{
  bool mask;
  switch (predicate)
    {
      case SIMDComparison::equal:
        mask = (left == right);
        break;
      case SIMDComparison::not_equal:
        mask = (left != right);
        break;
      case SIMDComparison::less_than:
        mask = (left < right);
        break;
      case SIMDComparison::less_than_or_equal:
        mask = (left <= right);
        break;
      case SIMDComparison::greater_than:
        mask = (left > right);
        break;
      case SIMDComparison::greater_than_or_equal:
        mask = (left >= right);
        break;
    }

  return mask ? true_value : false_value;
}


/**
 * 上述函数对非矢量化的VectorizedArray<Number,
 * 1>变体的特殊化。
 *
 *
 */
template <SIMDComparison predicate, typename Number>
DEAL_II_ALWAYS_INLINE inline VectorizedArray<Number, 1>
compare_and_apply_mask(const VectorizedArray<Number, 1> &left,
                       const VectorizedArray<Number, 1> &right,
                       const VectorizedArray<Number, 1> &true_value,
                       const VectorizedArray<Number, 1> &false_value)
{
  VectorizedArray<Number, 1> result;
  result.data = compare_and_apply_mask<predicate, Number>(left.data,
                                                          right.data,
                                                          true_value.data,
                                                          false_value.data);
  return result;
}

//@}

#ifndef DOXYGEN
#  if DEAL_II_VECTORIZATION_WIDTH_IN_BITS >= 512 && defined(__AVX512F__)

template <SIMDComparison predicate>
DEAL_II_ALWAYS_INLINE inline VectorizedArray<float, 16>
compare_and_apply_mask(const VectorizedArray<float, 16> &left,
                       const VectorizedArray<float, 16> &right,
                       const VectorizedArray<float, 16> &true_values,
                       const VectorizedArray<float, 16> &false_values)
{
  const __mmask16 mask =
    _mm512_cmp_ps_mask(left.data, right.data, static_cast<int>(predicate));
  VectorizedArray<float, 16> result;
  result.data = _mm512_mask_mov_ps(false_values.data, mask, true_values.data);
  return result;
}



template <SIMDComparison predicate>
DEAL_II_ALWAYS_INLINE inline VectorizedArray<double, 8>
compare_and_apply_mask(const VectorizedArray<double, 8> &left,
                       const VectorizedArray<double, 8> &right,
                       const VectorizedArray<double, 8> &true_values,
                       const VectorizedArray<double, 8> &false_values)
{
  const __mmask16 mask =
    _mm512_cmp_pd_mask(left.data, right.data, static_cast<int>(predicate));
  VectorizedArray<double, 8> result;
  result.data = _mm512_mask_mov_pd(false_values.data, mask, true_values.data);
  return result;
}

#  endif

#  if DEAL_II_VECTORIZATION_WIDTH_IN_BITS >= 256 && defined(__AVX__)

template <SIMDComparison predicate>
DEAL_II_ALWAYS_INLINE inline VectorizedArray<float, 8>
compare_and_apply_mask(const VectorizedArray<float, 8> &left,
                       const VectorizedArray<float, 8> &right,
                       const VectorizedArray<float, 8> &true_values,
                       const VectorizedArray<float, 8> &false_values)
{
  const auto mask =
    _mm256_cmp_ps(left.data, right.data, static_cast<int>(predicate));

  VectorizedArray<float, 8> result;
  result.data = _mm256_or_ps(_mm256_and_ps(mask, true_values.data),
                             _mm256_andnot_ps(mask, false_values.data));
  return result;
}


template <SIMDComparison predicate>
DEAL_II_ALWAYS_INLINE inline VectorizedArray<double, 4>
compare_and_apply_mask(const VectorizedArray<double, 4> &left,
                       const VectorizedArray<double, 4> &right,
                       const VectorizedArray<double, 4> &true_values,
                       const VectorizedArray<double, 4> &false_values)
{
  const auto mask =
    _mm256_cmp_pd(left.data, right.data, static_cast<int>(predicate));

  VectorizedArray<double, 4> result;
  result.data = _mm256_or_pd(_mm256_and_pd(mask, true_values.data),
                             _mm256_andnot_pd(mask, false_values.data));
  return result;
}

#  endif

#  if DEAL_II_VECTORIZATION_WIDTH_IN_BITS >= 128 && defined(__SSE2__)

template <SIMDComparison predicate>
DEAL_II_ALWAYS_INLINE inline VectorizedArray<float, 4>
compare_and_apply_mask(const VectorizedArray<float, 4> &left,
                       const VectorizedArray<float, 4> &right,
                       const VectorizedArray<float, 4> &true_values,
                       const VectorizedArray<float, 4> &false_values)
{
  __m128 mask;
  switch (predicate)
    {
      case SIMDComparison::equal:
        mask = _mm_cmpeq_ps(left.data, right.data);
        break;
      case SIMDComparison::not_equal:
        mask = _mm_cmpneq_ps(left.data, right.data);
        break;
      case SIMDComparison::less_than:
        mask = _mm_cmplt_ps(left.data, right.data);
        break;
      case SIMDComparison::less_than_or_equal:
        mask = _mm_cmple_ps(left.data, right.data);
        break;
      case SIMDComparison::greater_than:
        mask = _mm_cmpgt_ps(left.data, right.data);
        break;
      case SIMDComparison::greater_than_or_equal:
        mask = _mm_cmpge_ps(left.data, right.data);
        break;
    }

  VectorizedArray<float, 4> result;
  result.data = _mm_or_ps(_mm_and_ps(mask, true_values.data),
                          _mm_andnot_ps(mask, false_values.data));

  return result;
}


template <SIMDComparison predicate>
DEAL_II_ALWAYS_INLINE inline VectorizedArray<double, 2>
compare_and_apply_mask(const VectorizedArray<double, 2> &left,
                       const VectorizedArray<double, 2> &right,
                       const VectorizedArray<double, 2> &true_values,
                       const VectorizedArray<double, 2> &false_values)
{
  __m128d mask;
  switch (predicate)
    {
      case SIMDComparison::equal:
        mask = _mm_cmpeq_pd(left.data, right.data);
        break;
      case SIMDComparison::not_equal:
        mask = _mm_cmpneq_pd(left.data, right.data);
        break;
      case SIMDComparison::less_than:
        mask = _mm_cmplt_pd(left.data, right.data);
        break;
      case SIMDComparison::less_than_or_equal:
        mask = _mm_cmple_pd(left.data, right.data);
        break;
      case SIMDComparison::greater_than:
        mask = _mm_cmpgt_pd(left.data, right.data);
        break;
      case SIMDComparison::greater_than_or_equal:
        mask = _mm_cmpge_pd(left.data, right.data);
        break;
    }

  VectorizedArray<double, 2> result;
  result.data = _mm_or_pd(_mm_and_pd(mask, true_values.data),
                          _mm_andnot_pd(mask, false_values.data));

  return result;
}

#  endif
#endif // DOXYGEN


DEAL_II_NAMESPACE_CLOSE

/**
 * 在VectorizedArray上实现来自cmath的函数。这些函数不在dealii命名空间中，以确保与cmath中的相关函数有类似的接口。相反，使用
 * std::sin. 调用它们。
 *
 *
 */
namespace std
{
  /**
   * 计算一个矢量数据域的正弦。结果将以<tt>{sin(x[0]),
   * sin(x[1]), ...,  sin(x[VectorizedArray::size()-1])}</tt>.   @relatesalso
   * VectorizedArray的形式返回。
   *
   */
  template <typename Number, std::size_t width>
  inline ::dealii::VectorizedArray<Number, width>
  sin(const ::dealii::VectorizedArray<Number, width> &x)
  {
    // put values in an array and later read in that array with an unaligned
    // read. This should save some instructions as compared to directly
    // setting the individual elements and also circumvents a compiler
    // optimization bug in gcc-4.6 with SSE2 (see also deal.II developers list
    // from April 2014, topic "matrix_free/step-48 Test").
    Number values[::dealii::VectorizedArray<Number, width>::size()];
    for (unsigned int i = 0; i < dealii::VectorizedArray<Number, width>::size();
         ++i)
      values[i] = std::sin(x[i]);
    ::dealii::VectorizedArray<Number, width> out;
    out.load(&values[0]);
    return out;
  }



  /**
   * 计算一个矢量数据域的余弦。结果以矢量数组的形式返回<tt>{cos(x[0]),
   * cos(x[1]), ..., cos(x[size()-1]) }</tt>。      @relatesalso VectorizedArray
   *
   */
  template <typename Number, std::size_t width>
  inline ::dealii::VectorizedArray<Number, width>
  cos(const ::dealii::VectorizedArray<Number, width> &x)
  {
    Number values[::dealii::VectorizedArray<Number, width>::size()];
    for (unsigned int i = 0; i < dealii::VectorizedArray<Number, width>::size();
         ++i)
      values[i] = std::cos(x[i]);
    ::dealii::VectorizedArray<Number, width> out;
    out.load(&values[0]);
    return out;
  }



  /**
   * 计算一个矢量数据域的正切。结果以矢量数组的形式返回<tt>{tan(x[0]),
   * tan(x[1]), ..., tan(x[size()-1]) }</tt>。      @relatesalso VectorizedArray
   *
   */
  template <typename Number, std::size_t width>
  inline ::dealii::VectorizedArray<Number, width>
  tan(const ::dealii::VectorizedArray<Number, width> &x)
  {
    Number values[::dealii::VectorizedArray<Number, width>::size()];
    for (unsigned int i = 0; i < dealii::VectorizedArray<Number, width>::size();
         ++i)
      values[i] = std::tan(x[i]);
    ::dealii::VectorizedArray<Number, width> out;
    out.load(&values[0]);
    return out;
  }



  /**
   * 计算一个矢量数据域的指数。结果以矢量数组的形式返回<tt>{exp(x[0]),
   * exp(x[1]), ..., exp(x[size()-1]) }</tt>。      @relatesalso
   * VectorizedArray
   *
   */
  template <typename Number, std::size_t width>
  inline ::dealii::VectorizedArray<Number, width>
  exp(const ::dealii::VectorizedArray<Number, width> &x)
  {
    Number values[::dealii::VectorizedArray<Number, width>::size()];
    for (unsigned int i = 0; i < dealii::VectorizedArray<Number, width>::size();
         ++i)
      values[i] = std::exp(x[i]);
    ::dealii::VectorizedArray<Number, width> out;
    out.load(&values[0]);
    return out;
  }



  /**
   * 计算一个矢量数据域的自然对数。结果以矢量数组的形式返回<tt>{log(x[0]),
   * log(x[1]), ..., log(x[size()-1])}</tt>。      @relatesalso
   * VectorizedArray
   *
   */
  template <typename Number, std::size_t width>
  inline ::dealii::VectorizedArray<Number, width>
  log(const ::dealii::VectorizedArray<Number, width> &x)
  {
    Number values[::dealii::VectorizedArray<Number, width>::size()];
    for (unsigned int i = 0; i < dealii::VectorizedArray<Number, width>::size();
         ++i)
      values[i] = std::log(x[i]);
    ::dealii::VectorizedArray<Number, width> out;
    out.load(&values[0]);
    return out;
  }



  /**
   * 计算一个矢量数据域的平方根。结果以矢量数组的形式返回<tt>{sqrt(x[0]),
   * sqrt(x[1]), ..., sqrt(x[size()-1]) }</tt>。      @relatesalso
   * 矢量Array
   *
   */
  template <typename Number, std::size_t width>
  inline ::dealii::VectorizedArray<Number, width>
  sqrt(const ::dealii::VectorizedArray<Number, width> &x)
  {
    return x.get_sqrt();
  }



  /**
   * 将给定的数字 @p x 提高到幂数 @p p
   * ，用于一个矢量数据域。结果以矢量数组的形式返回<tt>{pow(x[0],p),
   * pow(x[1],p), ..., pow(x[size()-1], p)}</tt>。      @relatesalso
   * 矢量Array
   *
   */
  template <typename Number, std::size_t width>
  inline ::dealii::VectorizedArray<Number, width>
  pow(const ::dealii::VectorizedArray<Number, width> &x, const Number p)
  {
    Number values[::dealii::VectorizedArray<Number, width>::size()];
    for (unsigned int i = 0; i < dealii::VectorizedArray<Number, width>::size();
         ++i)
      values[i] = std::pow(x[i], p);
    ::dealii::VectorizedArray<Number, width> out;
    out.load(&values[0]);
    return out;
  }



  /**
   * 将给定的数字 @p x 提高到幂数 @p p
   * ，用于一个矢量数据域。结果以矢量数组的形式返回<tt>{pow(x[0],p[0]),
   * pow(x[1],p[1]), ..., pow(x[size()-1],p[size()-1]) }</tt>。
   * @relatesalso VectorizedArray
   *
   */
  template <typename Number, std::size_t width>
  inline ::dealii::VectorizedArray<Number, width>
  pow(const ::dealii::VectorizedArray<Number, width> &x,
      const ::dealii::VectorizedArray<Number, width> &p)
  {
    Number values[::dealii::VectorizedArray<Number, width>::size()];
    for (unsigned int i = 0; i < dealii::VectorizedArray<Number, width>::size();
         ++i)
      values[i] = std::pow(x[i], p[i]);
    ::dealii::VectorizedArray<Number, width> out;
    out.load(&values[0]);
    return out;
  }



  /**
   * 计算一个矢量数据字段的绝对值（模数）。结果以矢量数组的形式返回<tt>{abs(x[0]),
   * abs(x[1]), ..., abs(x[size()-1]) }</tt>。      @relatesalso
   * VectorizedArray
   *
   */
  template <typename Number, std::size_t width>
  inline ::dealii::VectorizedArray<Number, width>
  abs(const ::dealii::VectorizedArray<Number, width> &x)
  {
    return x.get_abs();
  }



  /**
   * 计算两个矢量数据域的分量最大。结果以矢量数组的形式返回<tt>{max(x[0],y[0]),
   * max(x[1],y[1), ...}</tt>。      @relatesalso  VectorizedArray
   *
   */
  template <typename Number, std::size_t width>
  inline ::dealii::VectorizedArray<Number, width>
  max(const ::dealii::VectorizedArray<Number, width> &x,
      const ::dealii::VectorizedArray<Number, width> &y)
  {
    return x.get_max(y);
  }



  /**
   * 计算两个矢量数据域的分量最小值。结果以矢量数组的形式返回<tt>{min(x[0],y[0]),
   * min(x[1],y[1), ...}</tt>。      @relatesalso VectorizedArray
   *
   */
  template <typename Number, std::size_t width>
  inline ::dealii::VectorizedArray<Number, width>
  min(const ::dealii::VectorizedArray<Number, width> &x,
      const ::dealii::VectorizedArray<Number, width> &y)
  {
    return x.get_min(y);
  }



  /**
   * VectorizedArrayIterator的迭代器特质。
   *
   */
  template <class T>
  struct iterator_traits<dealii::VectorizedArrayIterator<T>>
  {
    using iterator_category = random_access_iterator_tag;
    using value_type        = typename T::value_type;
    using difference_type   = std::ptrdiff_t;
  };

} // namespace std

#endif


