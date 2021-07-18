//include/deal.II-translator/base/aligned_vector_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2021 by the deal.II authors
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


#ifndef dealii_aligned_vector_h
#define dealii_aligned_vector_h

#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/parallel.h>
#include <deal.II/base/utilities.h>

// boost::serialization::make_array used to be in array.hpp, but was
// moved to a different file in BOOST 1.64
#include <boost/version.hpp>
#if BOOST_VERSION >= 106400
#  include <boost/serialization/array_wrapper.hpp>
#else
#  include <boost/serialization/array.hpp>
#endif
#include <boost/serialization/split_member.hpp>

#include <cstring>
#include <memory>
#include <type_traits>



DEAL_II_NAMESPACE_OPEN


/**
 * 这是一个替代 std::vector
 * 的类，可与VectorizedArray和派生数据类型结合使用。它分配的内存与矢量数据类型的地址对齐（以避免分段故障，当编译器假设一个矢量数组类型的变量与某些内存地址对齐时，实际上并没有遵循这些规则）。这也可以通过证明
 * std::vector
 * 与一个用户定义的分配器来实现。另一方面，编写一个自己的小向量类可以让我们用TBB实现并行的复制和移动操作，插入deal.II风格的断言，并削减一些不必要的功能。注意，由于对齐的原因，这个向量比
 * std::vector
 * 更耗费内存，所以建议只在长向量上使用这个向量。
 *
 *
 */
template <class T>
class AlignedVector
{
public:
  /**
   * 声明所有容器中使用的标准类型。这些类型与<tt>C++</tt>标准库中的<tt>vector<...></tt>类中的类型平行。
   *
   */
  using value_type      = T;
  using pointer         = value_type *;
  using const_pointer   = const value_type *;
  using iterator        = value_type *;
  using const_iterator  = const value_type *;
  using reference       = value_type &;
  using const_reference = const value_type &;
  using size_type       = std::size_t;

  /**
   * 空构造函数。将向量大小设置为零。
   *
   */
  AlignedVector();

  /**
   * 将向量大小设置为给定的大小，并用T()初始化所有元素。
   * @dealiiOperationIsMultithreaded
   *
   */
  explicit AlignedVector(const size_type size, const T &init = T());

  /**
   * 解构器。
   *
   */
  ~AlignedVector() = default;

  /**
   * 复制构造器。      @dealiiOperationIsMultithreaded
   *
   */
  AlignedVector(const AlignedVector<T> &vec);

  /**
   * 移动构造器。通过窃取 @p vec.
   * 的内容创建一个新的对齐向量。
   *
   */
  AlignedVector(AlignedVector<T> &&vec) noexcept;

  /**
   * 赋值给输入向量 @p vec.   @dealiiOperationIsMultithreaded 。
   *
   */
  AlignedVector &
  operator=(const AlignedVector<T> &vec);

  /**
   * 移动赋值运算符。
   *
   */
  AlignedVector &
  operator=(AlignedVector<T> &&vec) noexcept;

  /**
   * 改变向量的大小。如果新的大小大于先前的大小，那么新的元素将被添加到向量的末端；如果
   * `std::is_trivial<T>`
   * 是`true`，这些元素将保持未初始化（即留在未定义的状态），如果
   * `std::is_trivial<T>` 是`false`，将被默认初始化。  关于
   * `std::is_trivial`
   * 的定义，见[here](https://en.cppreference.com/w/cpp/types/is_trivial)。
   * 如果新的大小小于先前的大小，那么如果
   * `std::is_trivial<T>`
   * 为`false`，最后的几个元素将被销毁，或者如果
   * `std::is_trivial<T>` 为`true`，将来将被简单地忽略。
   * 作为上述大纲的结果，这个函数的后缀"_fast "是指对于
   * "琐碎的
   * "类`T`，这个函数省略了构造函数/析构函数的调用，特别是新元素的初始化。
   * @note 这个方法只能对定义了默认构造函数的类 @p T
   * 进行调用， @p T(). 否则，编译会失败。
   *
   */
  void
  resize_fast(const size_type new_size);

  /**
   * 改变向量的大小。它保留以前可用的旧元素，并将每个新添加的元素初始化为一个默认构造的对象类型
   * @p T.
   * 如果新的向量大小比旧的短，除非新的大小为零，否则不会立即释放内存；但是，当前对象的大小当然会被设置为要求的值。被释放元素的析构器也被调用。
   * @dealiiOperationIsMultithreaded
   *
   */
  void
  resize(const size_type new_size);

  /**
   * 改变向量的大小。它保留以前可用的旧元素，并用提供的初始化器初始化每个新添加的元素。
   * 如果新的向量大小比旧的短，除非新的大小为零，否则不会立即释放内存；但是，当前对象的大小当然被设置为请求的值。
   * @note
   * 这个方法只能对定义了复制赋值操作符的类进行调用。否则，编译会失败。
   * @dealiiOperationIsMultithreaded
   *
   */
  void
  resize(const size_type new_size, const T &init);

  /**
   * 为 @p new_allocated_size 元素保留内存空间。    如果参数 @p
   * new_allocated_size
   * 小于当前存储元素的数量（通过调用size()表示），那么这个函数根本不做任何事情。但如果参数
   * @p new_allocated_size
   * 被设置为零，那么所有先前分配的内存就会被释放（这个操作就相当于直接调用clear()函数）。
   * 为了避免过于频繁的重新分配（这涉及到数据的复制），当给定的大小大于先前分配的大小时，这个函数将占用的内存量加倍。
   * 请注意，这个函数只改变对象可以*存储的元素数量，而不是它实际*存储的元素数量。因此，不会运行新创建的对象的构造函数或析构函数，尽管现有的元素可以被移动到一个新的位置（这涉及到在新的位置运行移动构造函数和在旧的位置运行析构函数）。
   *
   */
  void
  reserve(const size_type new_allocated_size);

  /**
   * 释放所有先前分配的内存，使向量处于相当于调用默认构造函数后的状态。
   *
   */
  void
  clear();

  /**
   * 在向量的末端插入一个元素，使向量的大小增加一个。注意，只要之前的空间不足以容纳新的元素，分配的大小就会翻倍。
   *
   */
  void
  push_back(const T in_data);

  /**
   * 返回向量的最后一个元素（读和写访问）。
   *
   */
  reference
  back();

  /**
   * 返回向量的最后一个元素（只读访问）。
   *
   */
  const_reference
  back() const;

  /**
   * 在由一个元素范围给定的向量的末尾插入几个元素。
   *
   */
  template <typename ForwardIterator>
  void
  insert_back(ForwardIterator begin, ForwardIterator end);

  /**
   * 用默认构造对象的size()拷贝来填充向量。
   * @note
   * 与其他fill()函数不同，这个方法也可以为没有定义复制赋值操作符的类调用。
   * @dealiiOperationIsMultithreaded
   *
   */
  void
  fill();

  /**
   * 用给定输入的size()拷贝来填充向量。
   * @note
   * 这个方法只能对定义了复制赋值操作符的类调用。否则，编译会失败。
   * @dealiiOperationIsMultithreaded
   *
   */
  void
  fill(const T &element);

  /**
   * 这个函数在MPI通信器的所有进程中复制由 @p root_process
   * 指示的进程上发现的状态。在 @p root_process
   * 以外的任何进程中发现的当前状态都会在这个进程中丢失。我们可以想象这个操作就像从根进程到所有其他进程对
   * Utilities::MPI::broadcast()
   * 的调用，尽管在实践中这个函数可能试图将数据移到每个承载MPI进程的机器上的共享内存区域，然后让这个机器上的所有MPI进程访问这个共享内存区域，而不是保留自己的副本。
   * 这个函数的意图是将大的数组从一个进程快速交换给其他进程，而不是在所有进程上计算或创建它。这特别适用于从磁盘加载的数据
   *
   * --比如说，大的数据表格
   *
   * 比起让每个进程自己从磁盘上读取数据，通过读取一次，然后在MPI宇宙中的所有进程中分发，更容易处理。
   * 具体来说，共享内存区域的使用允许在MPI宇宙中每个多核机器上只复制一次数据，而不是为每个MPI进程复制一次数据。如果今天的机器上的数据很大，每个共享内存空间可以很容易地容纳几十个MPI进程，这就可以节省大量内存。这个用例在TableBase类的文档中有所概述，因为当前的函数是从
   * TableBase::replicate_across_communicator().
   * 中调用的。的确，这个函数的主要原理是使基于TableBase的数据表在MPI进程中共享。
   * 这个函数并不意味着保持不同进程上的数据同步的模型，就像
   * parallel::distributed::Vector
   * 和其他向量类所做的那样，其中存在一个由每个进程拥有的向量的某些元素的概念，可能还有从其拥有的进程镜像到其他进程的幽灵元素。相反，当前对象的元素被简单地复制到其他进程中，把这个操作看作是在所有进程中创建一组`const`AlignedVector对象，在复制操作之后不应该再被改变，这是确保向量在所有进程中保持一致的唯一方法。这尤其是因为共享内存区域的使用，在一个MPI进程上对一个向量元素的任何修改也可能导致对其他进程上可见元素的修改，假设它们位于一个共享内存节点内。
   * @note
   * 在MPI进程之间使用共享内存需要检测的MPI安装支持必要的操作。
   * 这对于MPI 3.0和更高版本来说是这样的。
   * @note  这个功能并不便宜。它需要创建所提供 @p
   * communicator
   * 对象的子通信器，这通常是一个昂贵的操作。同样地，共享内存空间的生成也不是一个便宜的操作。因此，当目标是在进程之间共享大的只读数据表时，这个功能主要是有意义的；例子是在启动时加载数据表，然后在程序的运行时间内使用。
   * 在这种情况下，运行这个函数的启动成本可以随着时间的推移而摊销，而且在具有大核心数的机器上，许多MPI进程在同一台机器上运行时，不必在每个进程上存储表所带来的潜在内存节省可能是相当大的。
   * @note  这个函数只有在数据类型`T`是 "自足
   * "的情况下才有意义，也就是说，如果它的所有信息都存储在其成员变量中，并且没有一个成员变量是指向内存的其他部分。这是因为如果一个类型`T`确实有指向内存其他部分的指针，那么将`T`移到一个共享内存空间不会导致其他进程访问该对象用其成员变量指针指向的数据。这些数据仍然只存在于一个进程中，并且通常在其他进程无法访问的内存区域。
   * 因此，这个函数的通常使用情况是共享简单对象的数组，如`double's或`int's。
   * @note
   * 调用该函数后，不同MPI进程的对象共享一个共同的状态。这意味着某些操作变得
   * "集体"，即必须在所有参与的处理器上同时调用。特别是，你不能再在一个MPI进程上调用resize()、reserve()或clear()。
   *
   * - 你必须在所有进程上同时这样做，因为它们必须为这些操作进行通信。如果你不这样做，你很可能会得到一个死锁，可能很难调试。推而广之，这个只集体调整大小的规则也延伸到这个函数本身。你不能连续调用它两次，因为这意味着首先除了`root_process'以外的所有进程都要扔掉他们的数据，这不是一个集体操作。一般来说，这些关于可以做什么和不可以做什么的限制，暗示了上面评论的正确性。你应该把一个当前函数被调用的AlignedVector视为`const'，在调用析构器之前，不能对其进行进一步的操作。
   *
   */
  void
  replicate_across_communicator(const MPI_Comm &   communicator,
                                const unsigned int root_process);

  /**
   * 将给定的向量与调用的向量交换。
   *
   */
  void
  swap(AlignedVector<T> &vec);

  /**
   * 返回该向量是否为空，即其大小为零。
   *
   */
  bool
  empty() const;

  /**
   * 返回该向量的大小。
   *
   */
  size_type
  size() const;

  /**
   * 返回向量的容量，即这个向量在不重新分配的情况下所能容纳的大小。注意，容量（）>=大小（）。
   *
   */
  size_type
  capacity() const;

  /**
   * 读写访问向量中的条目 @p index 。
   *
   */
  reference operator[](const size_type index);

  /**
   * 只读访问向量中的条目 @p index 。
   *
   */
  const_reference operator[](const size_type index) const;

  /**
   * 返回一个指向底层数据缓冲区的指针。
   *
   */
  pointer
  data();

  /**
   * 返回一个指向底层数据缓冲区的常量指针。
   *
   */
  const_pointer
  data() const;

  /**
   * 返回一个指向数据阵列开头的读写指针。
   *
   */
  iterator
  begin();

  /**
   * 返回一个指向数据阵列末端的读写指针。
   *
   */
  iterator
  end();

  /**
   * 返回一个指向数据数组开始的只读指针。
   *
   */
  const_iterator
  begin() const;

  /**
   * 返回一个指向数据数组末端的只读指针。
   *
   */
  const_iterator
  end() const;

  /**
   * 返回该类中分配的内存的消耗量。如果底层类型 @p T
   * 自己分配了内存，这个内存就不计算在内。
   *
   */
  size_type
  memory_consumption() const;

  /**
   * 使用[BOOST序列化库](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html)将此对象的数据写到一个流中，以便进行序列化。
   *
   */
  template <class Archive>
  void
  save(Archive &ar, const unsigned int version) const;

  /**
   * 使用[BOOST序列化库](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html)从流中读取此对象的数据，以便进行序列化。
   *
   */
  template <class Archive>
  void
  load(Archive &ar, const unsigned int version);

#ifdef DOXYGEN
  /**
   * 使用[BOOST序列化库](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html)从流中写入和读取此对象的数据，以达到序列化的目的。
   *
   */
  template <class Archive>
  void
  serialize(Archive &archive, const unsigned int version);
#else
  // This macro defines the serialize() method that is compatible with
  // the templated save() and load() method that have been implemented.
  BOOST_SERIALIZATION_SPLIT_MEMBER()
#endif

private:
  /**
   * 指向实际数据阵列的指针。
   *
   */
  std::unique_ptr<T[], std::function<void(T *)>> elements;

  /**
   * 指向超过最后一个有效值的指针。
   *
   */
  T *used_elements_end;

  /**
   * 指向已分配内存的终点的指针。
   *
   */
  T *allocated_elements_end;
};


// ------------------------------- inline functions --------------------------

/**
 * 这个命名空间定义了AlignedVector中使用的复制和设置函数。当向量中有足够的元素时，这些函数是并行操作的。
 *
 *
 */
namespace internal
{
  /**
   * 一个给定内存位置范围的类，在这些内存位置上调用放置-新建操作符，并在那里复制构造`T`类型的对象。
   * 这个类是基于parallel.h中专门的for循环基类ParallelForLoop，其目的如下。当在AlignedVector上用apply_to_subranges调用一个并行for循环时，它为我们可能选择的每个不同的参数生成不同的代码（因为它是模板化的）。这就产生了大量的代码（例如，它使编译matrix_free.cc文件所需的内存增加了三倍，而且最终的对象大小也大了好几倍），这些代码是完全无用的。因此，这个类通过调用apply_to_subrange为所有可能的类型引导所有的复制命令，这使得复制操作更加简洁（感谢一个虚拟函数，其成本在这种情况下可以忽略不计）。
   * @relatesalso  AlignedVector
   *
   */
  template <typename T>
  class AlignedVectorCopyConstruct
    : private dealii::parallel::ParallelForInteger
  {
    static const std::size_t minimum_parallel_grain_size =
      160000 / sizeof(T) + 1;

  public:
    /**
     * 构造函数。如果有足够多的元素，就发出一个并行调用，否则就以串行方式工作。将
     * @p source_begin 和 @p source_end
     * 之间的半开区间的数据复制到 @p destination
     * 开始的数组中（通过调用带有放置new的复制构造函数）。
     * 源数组中的元素通过放置新的复制构造函数被简单地复制出来。
     *
     */
    AlignedVectorCopyConstruct(const T *const source_begin,
                               const T *const source_end,
                               T *const       destination)
      : source_(source_begin)
      , destination_(destination)
    {
      Assert(source_end >= source_begin, ExcInternalError());
      Assert(source_end == source_begin || destination != nullptr,
             ExcInternalError());
      const std::size_t size = source_end - source_begin;
      if (size < minimum_parallel_grain_size)
        AlignedVectorCopyConstruct::apply_to_subrange(0, size);
      else
        apply_parallel(0, size, minimum_parallel_grain_size);
    }

    /**
     * 这个方法在两个整数给定的子范围内将元素从源数组移动到构造函数中给定的目的地。
     *
     */
    virtual void
    apply_to_subrange(const std::size_t begin,
                      const std::size_t end) const override
    {
      if (end == begin)
        return;

      // for classes trivial assignment can use memcpy. cast element to
      // (void*) to silence compiler warning for virtual classes (they will
      // never arrive here because they are non-trivial).

      if (std::is_trivial<T>::value == true)
        std::memcpy(static_cast<void *>(destination_ + begin),
                    static_cast<const void *>(source_ + begin),
                    (end - begin) * sizeof(T));
      else
        for (std::size_t i = begin; i < end; ++i)
          new (&destination_[i]) T(source_[i]);
    }

  private:
    const T *const source_;
    T *const       destination_;
  };


  /**
   * 像AlignedVectorCopyConstruct一样，但是使用`T`的移动构造函数来创建新的对象。
   * @relatesalso  AlignedVector
   *
   */
  template <typename T>
  class AlignedVectorMoveConstruct
    : private dealii::parallel::ParallelForInteger
  {
    static const std::size_t minimum_parallel_grain_size =
      160000 / sizeof(T) + 1;

  public:
    /**
     * 构造器。如果有足够多的元素，发出一个并行调用，否则以串行方式工作。将数据从
     * @p source_begin 和 @p source_end 之间的半开区间移动到从 @p
     * destination
     * 开始的数组中（通过调用带有位置new的移动构造函数）。
     * 通过调用源区间的析构器（为后续调用free做准备），数据在两个数组之间移动。
     *
     */
    AlignedVectorMoveConstruct(T *const source_begin,
                               T *const source_end,
                               T *const destination)
      : source_(source_begin)
      , destination_(destination)
    {
      Assert(source_end >= source_begin, ExcInternalError());
      Assert(source_end == source_begin || destination != nullptr,
             ExcInternalError());
      const std::size_t size = source_end - source_begin;
      if (size < minimum_parallel_grain_size)
        AlignedVectorMoveConstruct::apply_to_subrange(0, size);
      else
        apply_parallel(0, size, minimum_parallel_grain_size);
    }

    /**
     * 这个方法在由两个整数给定的子范围上将元素从源数移到构造函数中给定的目的地。
     *
     */
    virtual void
    apply_to_subrange(const std::size_t begin,
                      const std::size_t end) const override
    {
      if (end == begin)
        return;

      // Classes with trivial assignment can use memcpy. cast element to
      // (void*) to silence compiler warning for virtual classes (they will
      // never arrive here because they are non-trivial).
      if (std::is_trivial<T>::value == true)
        std::memcpy(static_cast<void *>(destination_ + begin),
                    static_cast<void *>(source_ + begin),
                    (end - begin) * sizeof(T));
      else
        // For everything else just use the move constructor. The original
        // object remains alive and will be destroyed elsewhere.
        for (std::size_t i = begin; i < end; ++i)
          new (&destination_[i]) T(std::move(source_[i]));
    }

  private:
    T *const source_;
    T *const destination_;
  };


  /**
   * 一个给定内存位置范围的类，调用这些内存位置的place-new操作符（如果`initialize_memory==true`）或者只是将给定的初始化器复制到这个内存位置（如果`initialize_memory==false`）。后者适合于只有琐碎构造器的类，如内置类型`double`，`int`等，以及由这类类型组成的结构。
   * @tparam  initialize_memory
   * 决定set命令是否应该初始化内存（通过调用拷贝构造函数）或者使用拷贝赋值操作符。有必要使用模板来选择合适的操作，因为有些类可能只定义这两种操作中的一种。
   * @relatesalso  AlignedVector
   *
   */
  template <typename T, bool initialize_memory>
  class AlignedVectorInitialize : private dealii::parallel::ParallelForInteger
  {
    static const std::size_t minimum_parallel_grain_size =
      160000 / sizeof(T) + 1;

  public:
    /**
     * 构造函数。如果有足够多的元素，发出一个并行的调用，否则以串行方式工作。
     *
     */
    AlignedVectorInitialize(const std::size_t size,
                            const T &         element,
                            T *const          destination)
      : element_(element)
      , destination_(destination)
      , trivial_element(false)
    {
      if (size == 0)
        return;
      Assert(destination != nullptr, ExcInternalError());

      // do not use memcmp for long double because on some systems it does not
      // completely fill its memory and may lead to false positives in
      // e.g. valgrind
      if (std::is_trivial<T>::value == true &&
          std::is_same<T, long double>::value == false)
        {
          const unsigned char zero[sizeof(T)] = {};
          // cast element to (void*) to silence compiler warning for virtual
          // classes (they will never arrive here because they are
          // non-trivial).
          if (std::memcmp(zero,
                          static_cast<const void *>(&element),
                          sizeof(T)) == 0)
            trivial_element = true;
        }
      if (size < minimum_parallel_grain_size)
        AlignedVectorInitialize::apply_to_subrange(0, size);
      else
        apply_parallel(0, size, minimum_parallel_grain_size);
    }

    /**
     * 这是在一个由两个整数给定的子范围内设置元素。
     *
     */
    virtual void
    apply_to_subrange(const std::size_t begin,
                      const std::size_t end) const override
    {
      // for classes with trivial assignment of zero can use memset. cast
      // element to (void*) to silence compiler warning for virtual
      // classes (they will never arrive here because they are
      // non-trivial).
      if (std::is_trivial<T>::value == true && trivial_element)
        std::memset(static_cast<void *>(destination_ + begin),
                    0,
                    (end - begin) * sizeof(T));
      else
        copy_construct_or_assign(
          begin, end, std::integral_constant<bool, initialize_memory>());
    }

  private:
    const T &  element_;
    mutable T *destination_;
    bool       trivial_element;

    // copy assignment operation
    void
    copy_construct_or_assign(const std::size_t begin,
                             const std::size_t end,
                             std::integral_constant<bool, false>) const
    {
      for (std::size_t i = begin; i < end; ++i)
        destination_[i] = element_;
    }

    // copy constructor (memory initialization)
    void
    copy_construct_or_assign(const std::size_t begin,
                             const std::size_t end,
                             std::integral_constant<bool, true>) const
    {
      for (std::size_t i = begin; i < end; ++i)
        new (&destination_[i]) T(element_);
    }
  };



  /**
   * 像AlignedVectorInitialize一样，但使用默认构造的对象作为初始化器。
   * @tparam  initialize_memory
   * 设置set命令是否应该初始化内存（通过调用拷贝构造器），或者宁可使用拷贝赋值操作符。模板对于选择适当的操作是必要的，因为有些类可能只定义这两种操作中的一种。
   * @relatesalso  AlignedVector
   *
   */
  template <typename T, bool initialize_memory>
  class AlignedVectorDefaultInitialize
    : private dealii::parallel::ParallelForInteger
  {
    static const std::size_t minimum_parallel_grain_size =
      160000 / sizeof(T) + 1;

  public:
    /**
     * 构造函数。如果有足够多的元素，发出一个并行的调用，否则以串行方式工作。
     *
     */
    AlignedVectorDefaultInitialize(const std::size_t size, T *const destination)
      : destination_(destination)
    {
      if (size == 0)
        return;
      Assert(destination != nullptr, ExcInternalError());

      if (size < minimum_parallel_grain_size)
        AlignedVectorDefaultInitialize::apply_to_subrange(0, size);
      else
        apply_parallel(0, size, minimum_parallel_grain_size);
    }

    /**
     * 这个初始化元素的子范围由两个整数给出。
     *
     */
    virtual void
    apply_to_subrange(const std::size_t begin,
                      const std::size_t end) const override
    {
      // for classes with trivial assignment of zero can use memset. cast
      // element to (void*) to silence compiler warning for virtual
      // classes (they will never arrive here because they are
      // non-trivial).
      if (std::is_trivial<T>::value == true)
        std::memset(static_cast<void *>(destination_ + begin),
                    0,
                    (end - begin) * sizeof(T));
      else
        default_construct_or_assign(
          begin, end, std::integral_constant<bool, initialize_memory>());
    }

  private:
    mutable T *destination_;

    // copy assignment operation
    void
    default_construct_or_assign(const std::size_t begin,
                                const std::size_t end,
                                std::integral_constant<bool, false>) const
    {
      for (std::size_t i = begin; i < end; ++i)
        destination_[i] = std::move(T());
    }

    // copy constructor (memory initialization)
    void
    default_construct_or_assign(const std::size_t begin,
                                const std::size_t end,
                                std::integral_constant<bool, true>) const
    {
      for (std::size_t i = begin; i < end; ++i)
        new (&destination_[i]) T;
    }
  };

} // end of namespace internal


#ifndef DOXYGEN


template <class T>
inline AlignedVector<T>::AlignedVector()
  : elements(nullptr, [](T *) { Assert(false, ExcInternalError()); })
  , used_elements_end(nullptr)
  , allocated_elements_end(nullptr)
{}



template <class T>
inline AlignedVector<T>::AlignedVector(const size_type size, const T &init)
  : elements(nullptr, [](T *) { Assert(false, ExcInternalError()); })
  , used_elements_end(nullptr)
  , allocated_elements_end(nullptr)
{
  if (size > 0)
    resize(size, init);
}



template <class T>
inline AlignedVector<T>::AlignedVector(const AlignedVector<T> &vec)
  : elements(nullptr, [](T *) { Assert(false, ExcInternalError()); })
  , used_elements_end(nullptr)
  , allocated_elements_end(nullptr)
{
  // copy the data from vec
  reserve(vec.size());
  used_elements_end = allocated_elements_end;
  internal::AlignedVectorCopyConstruct<T>(vec.elements.get(),
                                          vec.used_elements_end,
                                          elements.get());
}



template <class T>
inline AlignedVector<T>::AlignedVector(AlignedVector<T> &&vec) noexcept
  : AlignedVector<T>()
{
  // forward to the move operator
  *this = std::move(vec);
}



template <class T>
inline AlignedVector<T> &
AlignedVector<T>::operator=(const AlignedVector<T> &vec)
{
  resize(0);
  resize_fast(vec.used_elements_end - vec.elements.get());
  internal::AlignedVectorCopyConstruct<T>(vec.elements.get(),
                                          vec.used_elements_end,
                                          elements.get());
  return *this;
}



template <class T>
inline AlignedVector<T> &
AlignedVector<T>::operator=(AlignedVector<T> &&vec) noexcept
{
  clear();

  // Move the actual data in the 'elements' object. One problem is that this
  // also moves the deleter object, but the deleter object is a lambda function
  // that references 'this' (i.e., the 'this' pointer of the *moved-from*
  // object). So what we actually do is steal the pointer via
  // std::unique_ptr::release() and then install our own deleter object that
  // mirrors the one used in reserve() below.
  elements = decltype(elements)(vec.elements.release(), [this](T *ptr) {
    if (ptr != nullptr)
      {
        Assert(this->used_elements_end != nullptr, ExcInternalError());

        if (std::is_trivial<T>::value == false)
          for (T *p = this->used_elements_end - 1; p >= ptr; --p)
            p->~T();
      }

    std::free(ptr);
  });

  // Then also steal the other pointers and clear them in the original object:
  used_elements_end      = vec.used_elements_end;
  allocated_elements_end = vec.allocated_elements_end;

  vec.used_elements_end      = nullptr;
  vec.allocated_elements_end = nullptr;

  return *this;
}



template <class T>
inline void
AlignedVector<T>::resize_fast(const size_type new_size)
{
  const size_type old_size = size();

  if (new_size == 0)
    clear();
  else if (new_size == old_size)
    {} // nothing to do here
  else if (new_size < old_size)
    {
      // call destructor on fields that are released, if the type requires it.
      // doing it backward releases the elements in reverse order as compared to
      // how they were created
      if (std::is_trivial<T>::value == false)
        for (T *p = used_elements_end - 1; p >= elements.get() + new_size; --p)
          p->~T();
      used_elements_end = elements.get() + new_size;
    }
  else // new_size > old_size
    {
      // Allocate more space, and claim that space as used
      reserve(new_size);
      used_elements_end = elements.get() + new_size;

      // need to still set the values in case the class is non-trivial because
      // virtual classes etc. need to run their (default) constructor
      if (std::is_trivial<T>::value == false)
        dealii::internal::AlignedVectorDefaultInitialize<T, true>(
          new_size - old_size, elements.get() + old_size);
    }
}



template <class T>
inline void
AlignedVector<T>::resize(const size_type new_size)
{
  const size_type old_size = size();

  if (new_size == 0)
    clear();
  else if (new_size == old_size)
    {} // nothing to do here
  else if (new_size < old_size)
    {
      // call destructor on fields that are released, if the type requires it.
      // doing it backward releases the elements in reverse order as compared to
      // how they were created
      if (std::is_trivial<T>::value == false)
        for (T *p = used_elements_end - 1; p >= elements.get() + new_size; --p)
          p->~T();
      used_elements_end = elements.get() + new_size;
    }
  else // new_size > old_size
    {
      // Allocate more space, and claim that space as used
      reserve(new_size);
      used_elements_end = elements.get() + new_size;

      // finally set the values to the default initializer
      dealii::internal::AlignedVectorDefaultInitialize<T, true>(
        new_size - old_size, elements.get() + old_size);
    }
}



template <class T>
inline void
AlignedVector<T>::resize(const size_type new_size, const T &init)
{
  const size_type old_size = size();

  if (new_size == 0)
    clear();
  else if (new_size == old_size)
    {} // nothing to do here
  else if (new_size < old_size)
    {
      // call destructor on fields that are released, if the type requires it.
      // doing it backward releases the elements in reverse order as compared to
      // how they were created
      if (std::is_trivial<T>::value == false)
        for (T *p = used_elements_end - 1; p >= elements.get() + new_size; --p)
          p->~T();
      used_elements_end = elements.get() + new_size;
    }
  else // new_size > old_size
    {
      // Allocate more space, and claim that space as used
      reserve(new_size);
      used_elements_end = elements.get() + new_size;

      // finally set the desired init values
      dealii::internal::AlignedVectorInitialize<T, true>(
        new_size - old_size, init, elements.get() + old_size);
    }
}



template <class T>
inline void
AlignedVector<T>::reserve(const size_type new_allocated_size)
{
  const size_type old_size           = used_elements_end - elements.get();
  const size_type old_allocated_size = allocated_elements_end - elements.get();
  if (new_allocated_size > old_allocated_size)
    {
      // if we continuously increase the size of the vector, we might be
      // reallocating a lot of times. therefore, try to increase the size more
      // aggressively
      const size_type new_size =
        std::max(new_allocated_size, 2 * old_allocated_size);

      // allocate and align along 64-byte boundaries (this is enough for all
      // levels of vectorization currently supported by deal.II)
      T *new_data_ptr;
      Utilities::System::posix_memalign(
        reinterpret_cast<void **>(&new_data_ptr), 64, new_size * sizeof(T));

      // Now create a deleter that encodes what should happen when the object is
      // released: We need to destroy the objects that are currently alive (in
      // reverse order, and then release the memory. Note that we catch the
      // 'this' pointer because the number of elements currently alive might
      // change over time.
      auto deleter = [this](T *ptr) {
        if (ptr != nullptr)
          {
            Assert(this->used_elements_end != nullptr, ExcInternalError());

            if (std::is_trivial<T>::value == false)
              for (T *p = this->used_elements_end - 1; p >= ptr; --p)
                p->~T();
          }

        std::free(ptr);
      };

      // copy whatever elements we need to retain
      if (new_allocated_size > 0)
        dealii::internal::AlignedVectorMoveConstruct<T>(
          elements.get(), elements.get() + old_size, new_data_ptr);

      // Now reset all of the member variables of the current object
      // based on the allocation above. Assigning to a std::unique_ptr
      // object also releases the previously pointed to memory.
      //
      // Note that at the time of releasing the old memory, 'used_elements_end'
      // still points to its previous value, and this is important for the
      // deleter object of the previously allocated array (see how it loops over
      // the to-be-destroyed elements a few lines above).
      elements               = decltype(elements)(new_data_ptr, deleter);
      used_elements_end      = elements.get() + old_size;
      allocated_elements_end = elements.get() + new_size;
    }
  else if (new_allocated_size == 0)
    clear();
  else // size_alloc < allocated_size
    {} // nothing to do here
}



template <class T>
inline void
AlignedVector<T>::clear()
{
  // Just release the memory (which also calls the destructor of the elements),
  // and then set the auxiliary pointers to invalid values.
  //
  // Note that at the time of releasing the old memory, 'used_elements_end'
  // still points to its previous value, and this is important for the
  // deleter object of the previously allocated array (see how it loops over
  // the to-be-destroyed elements a few lines above).
  elements.reset();
  used_elements_end      = nullptr;
  allocated_elements_end = nullptr;
}



template <class T>
inline void
AlignedVector<T>::push_back(const T in_data)
{
  Assert(used_elements_end <= allocated_elements_end, ExcInternalError());
  if (used_elements_end == allocated_elements_end)
    reserve(std::max(2 * capacity(), static_cast<size_type>(16)));
  if (std::is_trivial<T>::value == false)
    new (used_elements_end++) T(in_data);
  else
    *used_elements_end++ = in_data;
}



template <class T>
inline typename AlignedVector<T>::reference
AlignedVector<T>::back()
{
  AssertIndexRange(0, size());
  T *field = used_elements_end - 1;
  return *field;
}



template <class T>
inline typename AlignedVector<T>::const_reference
AlignedVector<T>::back() const
{
  AssertIndexRange(0, size());
  const T *field = used_elements_end - 1;
  return *field;
}



template <class T>
template <typename ForwardIterator>
inline void
AlignedVector<T>::insert_back(ForwardIterator begin, ForwardIterator end)
{
  const unsigned int old_size = size();
  reserve(old_size + (end - begin));
  for (; begin != end; ++begin, ++used_elements_end)
    {
      if (std::is_trivial<T>::value == false)
        new (used_elements_end) T;
      *used_elements_end = *begin;
    }
}



template <class T>
inline void
AlignedVector<T>::fill()
{
  dealii::internal::AlignedVectorDefaultInitialize<T, false>(size(),
                                                             elements.get());
}



template <class T>
inline void
AlignedVector<T>::fill(const T &value)
{
  dealii::internal::AlignedVectorInitialize<T, false>(size(),
                                                      value,
                                                      elements.get());
}



template <class T>
inline void
AlignedVector<T>::replicate_across_communicator(const MPI_Comm &   communicator,
                                                const unsigned int root_process)
{
#  ifdef DEAL_II_WITH_MPI
#    if DEAL_II_MPI_VERSION_GTE(3, 0)

  // **** Step 0 ****
  // All but the root process no longer need their data, so release the memory
  // used to store the previous elements.
  if (Utilities::MPI::this_mpi_process(communicator) != root_process)
    {
      elements.reset();
      used_elements_end      = nullptr;
      allocated_elements_end = nullptr;
    }

  // **** Step 1 ****
  // Create communicators for each group of processes that can use
  // shared memory areas. Within each of these groups, we don't care about
  // which rank each of the old processes gets except that we would like to
  // make sure that the (global) root process will have rank=0 within
  // its own sub-communicator. We can do that through the third argument of
  // MPI_Comm_split_type (the "key") which is an integer meant to indicate the
  // order of processes within the split communicators, and we should set it to
  // zero for the root processes and one for all others -- which means that
  // for all of these other processes, MPI can choose whatever order it
  // wants because they have the same key (MPI then documents that these ties
  // will be broken according to these processes' rank in the old group).
  //
  // At least that's the theory. In practice, the MPI implementation where
  // this function was developed on does not seem to do that. (Bug report
  // is here: https://github.com/open-mpi/ompi/issues/8854)
  // We work around this by letting MPI_Comm_split_type choose whatever
  // rank it wants, and then reshuffle with MPI_Comm_split in a second
  // step -- not elegant, nor efficient, but seems to work:
  MPI_Comm shmem_group_communicator;
  {
    MPI_Comm shmem_group_communicator_temp;
    int      ierr = MPI_Comm_split_type(communicator,
                                   MPI_COMM_TYPE_SHARED,
                                    /* key */  0,
                                   MPI_INFO_NULL,
                                   &shmem_group_communicator_temp);
    AssertThrowMPI(ierr);

    const int key =
      (Utilities::MPI::this_mpi_process(communicator) == root_process ? 0 : 1);
    ierr = MPI_Comm_split(shmem_group_communicator_temp,
                           /* color */  0,
                          key,
                          &shmem_group_communicator);
    AssertThrowMPI(ierr);

    // Verify the explanation from above
    if (Utilities::MPI::this_mpi_process(communicator) == root_process)
      Assert(Utilities::MPI::this_mpi_process(shmem_group_communicator) == 0,
             ExcInternalError());

    // And get rid of the temporary communicator
    ierr = MPI_Comm_free(&shmem_group_communicator_temp);
    AssertThrowMPI(ierr);
  }
  const bool is_shmem_root =
    Utilities::MPI::this_mpi_process(shmem_group_communicator) == 0;

  // **** Step 2 ****
  // We then have to send the state of the current object from the
  // root process to one exemplar in each shmem group. To this end,
  // we create another subcommunicator that includes the ranks zero
  // of all shmem groups, and because of the trick above, we know
  // that this also includes the original root process.
  //
  // There are different ways of creating a "shmem_roots_communicator".
  // The conceptually easiest way is to create an MPI_Group that only
  // includes the shmem roots and then create a communicator from this
  // via MPI_Comm_create or MPI_Comm_create_group. The problem
  // with this is that we would have to exchange among all processes
  // which ones are shmem roots and which are not. This is awkward.
  //
  // A simpler way is to use MPI_Comm_split that uses "colors" to
  // indicate which sub-communicator each process wants to be in.
  // We use color=0 to indicate the group of shmem roots, and color=1
  // for all other processes -- the latter will simply not ever do
  // anything among themselves with the communicator so created.
  //
  // Using MPI_Comm_split has the additional benefit that, just as above,
  // we can choose where each rank will end up in shmem_roots_communicator.
  // We again set key=0 for the original root_process, and key=1 for all other
  // ranks; then, the global root becomes rank=0 on the
  // shmem_roots_communicator. We don't care how the other processes are
  // ordered.
  MPI_Comm shmem_roots_communicator;
  {
    const int key =
      (Utilities::MPI::this_mpi_process(communicator) == root_process ? 0 : 1);

    const int ierr = MPI_Comm_split(communicator,
                                     /*color=*/ 
                                    (is_shmem_root ? 0 : 1),
                                    key,
                                    &shmem_roots_communicator);
    AssertThrowMPI(ierr);

    // Again verify the explanation from above
    if (Utilities::MPI::this_mpi_process(communicator) == root_process)
      Assert(Utilities::MPI::this_mpi_process(shmem_roots_communicator) == 0,
             ExcInternalError());
  }

  const unsigned int shmem_roots_root_rank = 0;
  const bool         is_shmem_roots_root =
    (is_shmem_root && (Utilities::MPI::this_mpi_process(
                         shmem_roots_communicator) == shmem_roots_root_rank));

  // Now let the original root_process broadcast the current object to all
  // shmem roots. We know that the last rank is the original root process that
  // has all of the data.
  if (is_shmem_root)
    {
      if (std::is_trivial<T>::value)
        {
          // The data is "trivial", i.e., we can copy things directly without
          // having to go through the serialization/deserialization machinery of
          // Utilities::MPI::broadcast.
          //
          // In that case, first tell all of the other shmem roots how many
          // elements we will have to deal with, and let them resize their
          // (non-shared) arrays.
          const size_type new_size =
            Utilities::MPI::broadcast(shmem_roots_communicator,
                                      size(),
                                      shmem_roots_root_rank);
          if (is_shmem_roots_root == false)
            resize(new_size);

          // Then directly copy from the root process into these buffers
          int ierr = MPI_Bcast(elements.get(),
                               sizeof(T) * new_size,
                               MPI_CHAR,
                               shmem_roots_root_rank,
                               shmem_roots_communicator);
          AssertThrowMPI(ierr);
        }
      else
        {
          // The objects to be sent around are not "trivial", and so we have
          // to go through the serialization/deserialization machinery. On all
          // but the sending process, overwrite the current state with the
          // vector just broadcast.
          //
          // On the root rank, this would lead to resetting the 'entries'
          // pointer, which would trigger the deleter which would lead to a
          // deadlock. So we just send the result of the broadcast() call to
          // nirvana on the root process and keep our current state.
          if (Utilities::MPI::this_mpi_process(shmem_roots_communicator) == 0)
            Utilities::MPI::broadcast(shmem_roots_communicator,
                                      *this,
                                      shmem_roots_root_rank);
          else
            *this = Utilities::MPI::broadcast(shmem_roots_communicator,
                                              *this,
                                              shmem_roots_root_rank);
        }
    }

  // We no longer need the shmem roots communicator, so get rid of it
  {
    const int ierr = MPI_Comm_free(&shmem_roots_communicator);
    AssertThrowMPI(ierr);
  }


  // **** Step 3 ****
  // At this point, all shmem groups have one shmem root process that has
  // a copy of the data. This is the point where each shmem group should
  // establish a shmem area to put the data into. As mentioned above,
  // we know that the shmem roots are the last rank in their respective
  // shmem_group_communicator.
  //
  // The process for all of this works as follows: While all processes in
  // the shmem group participate in the generation of the shmem memory window,
  // only the shmem root actually allocates any memory -- the rest just
  // allocate zero bytes of their own. We allocate space for exactly
  // size() elements (computed on the shmem_root that already has the data)
  // and add however many bytes are necessary so that we know that we can align
  // things to 64-byte boundaries. The worst case happens if the memory system
  // gives us a pointer to an address one byte past a desired alignment
  // boundary, and in that case aligning the memory will require us to waste the
  // first (align_by-1) bytes. So we have to ask for
  //   size() * sizeof(T) + (align_by - 1)
  // bytes.
  //
  // Before MPI 4.0, there was no way to specify that we want memory aligned to
  // a certain number of bytes. This is going to come back to bite us further
  // down below when we try to get a properly aligned pointer to our memory
  // region, see the commentary there. Starting with MPI 4.0, one can set a
  // flag in an MPI_Info structure that requests a desired alignment, so we do
  // this for forward compatibility; MPI implementations ignore flags they don't
  // know anything about, and so setting this flag is backward compatible also
  // to older MPI versions.
  //
  // There is one final piece we can already take care of here. At the beginning
  // of all of this, only the shmem_root knows how many elements there are in
  // the array. But at the end of it, all processes of course need to know. We
  // could put this information somewhere into the shmem area, along with the
  // other data, but that seems clumsy. It turns out that when calling
  // MPI_Win_allocate_shared, we are asked for the value of a parameter called
  // 'disp_unit' whose meaning is difficult to determine from the MPI
  // documentation, and that we do not actually need. So we "abuse" it a bit: On
  // the shmem root, we put the array size into it. Later on, the remaining
  // processes can query the shmem root's value of 'disp_unit', and so will be
  // able to learn about the array size that way.
  MPI_Win        shmem_window;
  void *         base_ptr;
  const MPI_Aint align_by = 64;
  const MPI_Aint alloc_size =
    Utilities::MPI::broadcast(shmem_group_communicator,
                              (size() * sizeof(T) + (align_by - 1)),
                              0);

  {
    const int disp_unit = (is_shmem_root ? size() : 1);

    int ierr;

    MPI_Info mpi_info;
    ierr = MPI_Info_create(&mpi_info);
    AssertThrowMPI(ierr);
    ierr = MPI_Info_set(mpi_info,
                        "mpi_minimum_memory_alignment",
                        std::to_string(align_by).c_str());
    AssertThrowMPI(ierr);
    ierr = MPI_Win_allocate_shared((is_shmem_root ? alloc_size : 0),
                                   disp_unit,
                                   mpi_info,
                                   shmem_group_communicator,
                                   &base_ptr,
                                   &shmem_window);
    AssertThrowMPI(ierr);

    ierr = MPI_Info_free(&mpi_info);
    AssertThrowMPI(ierr);
  }


  // **** Step 4 ****
  // The next step is to teach all non-shmem root processes what the pointer to
  // the array is that the shmem-root created. MPI has a nifty way for this
  // given that only a single process actually allocated memory in the window:
  // When calling MPI_Win_shared_query, the MPI documentation says that
  // "When rank is MPI_PROC_NULL, the pointer, disp_unit, and size returned are
  // the pointer, disp_unit, and size of the memory segment belonging the lowest
  // rank that specified size > 0. If all processes in the group attached to the
  // window specified size = 0, then the call returns size = 0 and a baseptr as
  // if MPI_ALLOC_MEM was called with size = 0."
  //
  // This will allow us to obtain the pointer to the shmem root's memory area,
  // which is the only one we care about. (None of the other processes have
  // even allocated any memory.) But this will also retrieve the shmem root's
  // disp_unit, which in step 3 above we have abused to pass along the number of
  // elements in the array.
  //
  // We don't need to do this on the shmem root process: This process has
  // already gotten its base_ptr correctly set above, and we can determine the
  // array size by just calling size().
  unsigned int array_size =
    (is_shmem_root ? size() : numbers::invalid_unsigned_int);
  if (is_shmem_root == false)
    {
      int       disp_unit;
      MPI_Aint  alloc_size; // not actually used
      const int ierr = MPI_Win_shared_query(
        shmem_window, MPI_PROC_NULL, &alloc_size, &disp_unit, &base_ptr);
      AssertThrowMPI(ierr);

      // Make sure we actually got a pointer, and also unpack the array size as
      // discussed above.
      Assert(base_ptr != nullptr, ExcInternalError());

      array_size = disp_unit;
    }


  // **** Step 5 ****
  // Now that all processes know the address of the space that is visible to
  // everyone, we need to figure out whether it is properly aligned and if not,
  // find the next aligned address.
  //
  // std::align does that, but it also modifies its last two arguments. The
  // documentation of that function at
  // https://en.cppreference.com/w/cpp/memory/align is not entirely clear, but I
  // *think* that the following should do given that we do not use base_ptr and
  // available_space any further after the call to std::align.
  std::size_t available_space       = alloc_size;
  void *      base_ptr_backup       = base_ptr;
  T *         aligned_shmem_pointer = static_cast<T *>(
    std::align(align_by, array_size * sizeof(T), base_ptr, available_space));
  Assert(aligned_shmem_pointer != nullptr, ExcInternalError());

  // There is one step to guard against. It is *conceivable* that the base_ptr
  // we have previously obtained from MPI_Win_shared_query is mapped so
  // awkwardly into the different MPI processes' memory spaces that it is
  // aligned in one memory space, but not another. In that case, different
  // processes would align base_ptr differently, and adjust available_space
  // differently. We can check that by making sure that the max (or min) over
  // all processes is equal to every process's value. If that's not the case,
  // then the whole idea of aligning above is wrong and we need to rethink what
  // it means to align data in a shared memory space.
  //
  // One might be tempted to think that this is not how MPI implementations
  // actually arrange things. Alas, when developing this functionality in 2021,
  // this is really how at least OpenMPI ends up doing things. (This is with an
  // OpenMPI implementation of MPI 3.1, so it does not support the flag we set
  // in the MPI_Info structure above when allocating the memory window.) Indeed,
  // when running this code on three processes, one ends up with base_ptr values
  // of
  //     base_ptr=0x7f0842f02108
  //     base_ptr=0x7fc0a47881d0
  //     base_ptr=0x7f64872db108
  // which, most annoyingly, are aligned to 8 and 16 byte boundaries -- so there
  // is no common offset std::align could find that leads to a 64-byte
  // aligned memory address in all three memory spaces. That's a tremendous
  // nuisance and there is really nothing we can do about this other than just
  // fall back on the (unaligned) base_ptr in that case.
  if (Utilities::MPI::min(available_space, shmem_group_communicator) !=
      Utilities::MPI::max(available_space, shmem_group_communicator))
    aligned_shmem_pointer = static_cast<T *>(base_ptr_backup);


  // **** Step 6 ****
  // If this is the shmem root process, we need to copy the data into the
  // shared memory space.
  if (is_shmem_root)
    {
      if (std::is_trivial<T>::value == true)
        std::memcpy(aligned_shmem_pointer, elements.get(), sizeof(T) * size());
      else
        for (std::size_t i = 0; i < size(); ++i)
          new (&aligned_shmem_pointer[i]) T(std::move(elements[i]));
    }

  // Make sure that the shared memory host has copied the data before we try to
  // access it.
  const int ierr = MPI_Barrier(shmem_group_communicator);
  AssertThrowMPI(ierr);

  // **** Step 7 ****
  // Finally, we need to set the pointers of this object to what we just
  // learned. This also releases all memory that may have been in use
  // previously.
  //
  // The part that is a bit tricky is how to write the deleter of this
  // shared memory object. When we want to get rid of it, we need to
  // also release the MPI_Win object along with the shmem_group_communicator
  // object. That's because as long as we use the shared memory, we still need
  // to hold on to the MPI_Win object, and the MPI_Win object is based on the
  // communicator. (The former is definitely true, the latter is not quite clear
  // from the MPI documentation, but seems reasonable.) So we need to have a
  // deleter for the pointer that ensures that upon release of the memory, we
  // not only call the destructor of these memory elements (but only once, on
  // the shmem root!) but also destroy the MPI_Win and the communicator. All of
  // that is encapsulated in the following call where the deleter makes copies
  // of the arguments in the lambda capture.
  elements =
    decltype(elements)(aligned_shmem_pointer,
                       [is_shmem_root,
                        array_size,
                        aligned_shmem_pointer,
                        shmem_group_communicator,
                        shmem_window](T *) mutable {
                         if (is_shmem_root)
                           for (unsigned int i = 0; i < array_size; ++i)
                             aligned_shmem_pointer[i].~T();

                         int ierr;
                         ierr = MPI_Win_free(&shmem_window);
                         AssertThrowMPI(ierr);

                         ierr = MPI_Comm_free(&shmem_group_communicator);
                         AssertThrowMPI(ierr);
                       });

  // We then also have to set the other two pointers that define the state of
  // the current object. Note that the new buffer size is exactly as large as
  // necessary, i.e., can store size() elements, regardless of the number of
  // allocated elements in the original objects.
  used_elements_end      = elements.get() + array_size;
  allocated_elements_end = used_elements_end;

  // **** Consistency check ****
  // At this point, each process should have a copy of the data.
  // Verify this in some sort of round-about way
#      ifdef DEBUG
  const std::vector<char> packed_data = Utilities::pack(*this);
  const int               hash =
    std::accumulate(packed_data.begin(), packed_data.end(), int(0));
  Assert(Utilities::MPI::max(hash, communicator) == hash, ExcInternalError());
#      endif



#    else
  // If we only have MPI 2.x, then simply broadcast the current object to all
  // other processes and forego the idea of using shmem
  *this = Utilities::MPI::broadcast(communicator, *this, root_process);
#    endif
#  else
  // No MPI -> nothing to replicate
  (void)communicator;
  (void)root_process;
#  endif
}



template <class T>
inline void
AlignedVector<T>::swap(AlignedVector<T> &vec)
{
  // Swap the data in the 'elements' objects. One problem is that this
  // also moves the deleter object, but the deleter object is a lambda function
  // that references 'this' (i.e., the 'this' pointer of the *moved-from*
  // object). So what we actually do is steal the pointer via
  // std::unique_ptr::release() and then install our own deleter object that
  // mirrors the one used in reserve() below.
  //
  // We have to do the same for the other object
  T *this_element_pointer = elements.release();

  elements = decltype(elements)(vec.elements.release(), [this](T *ptr) {
    if (ptr != nullptr)
      {
        Assert(this->used_elements_end != nullptr, ExcInternalError());

        if (std::is_trivial<T>::value == false)
          for (T *p = this->used_elements_end - 1; p >= ptr; --p)
            p->~T();
      }

    std::free(ptr);
  });

  vec.elements = decltype(vec.elements)(this_element_pointer, [&vec](T *ptr) {
    if (ptr != nullptr)
      {
        Assert(vec.used_elements_end != nullptr, ExcInternalError());

        if (std::is_trivial<T>::value == false)
          for (T *p = vec.used_elements_end - 1; p >= ptr; --p)
            p->~T();
      }

    std::free(ptr);
  });

  std::swap(used_elements_end, vec.used_elements_end);
  std::swap(allocated_elements_end, vec.allocated_elements_end);
}



template <class T>
inline bool
AlignedVector<T>::empty() const
{
  return used_elements_end == elements.get();
}



template <class T>
inline typename AlignedVector<T>::size_type
AlignedVector<T>::size() const
{
  return used_elements_end - elements.get();
}



template <class T>
inline typename AlignedVector<T>::size_type
AlignedVector<T>::capacity() const
{
  return allocated_elements_end - elements.get();
}



template <class T>
inline typename AlignedVector<T>::reference AlignedVector<T>::
                                            operator[](const size_type index)
{
  AssertIndexRange(index, size());
  return elements[index];
}



template <class T>
inline typename AlignedVector<T>::const_reference AlignedVector<T>::
                                                  operator[](const size_type index) const
{
  AssertIndexRange(index, size());
  return elements[index];
}



template <typename T>
inline typename AlignedVector<T>::pointer
AlignedVector<T>::data()
{
  return elements.get();
}



template <typename T>
inline typename AlignedVector<T>::const_pointer
AlignedVector<T>::data() const
{
  return elements.get();
}



template <class T>
inline typename AlignedVector<T>::iterator
AlignedVector<T>::begin()
{
  return elements.get();
}



template <class T>
inline typename AlignedVector<T>::iterator
AlignedVector<T>::end()
{
  return used_elements_end;
}



template <class T>
inline typename AlignedVector<T>::const_iterator
AlignedVector<T>::begin() const
{
  return elements.get();
}



template <class T>
inline typename AlignedVector<T>::const_iterator
AlignedVector<T>::end() const
{
  return used_elements_end;
}



template <class T>
template <class Archive>
inline void
AlignedVector<T>::save(Archive &ar, const unsigned int) const
{
  size_type vec_size = size();
  ar &      vec_size;
  if (vec_size > 0)
    ar &boost::serialization::make_array(elements.get(), vec_size);
}



template <class T>
template <class Archive>
inline void
AlignedVector<T>::load(Archive &ar, const unsigned int)
{
  size_type vec_size = 0;
  ar &      vec_size;

  if (vec_size > 0)
    {
      reserve(vec_size);
      ar &boost::serialization::make_array(elements.get(), vec_size);
      used_elements_end = elements.get() + vec_size;
    }
}



template <class T>
inline typename AlignedVector<T>::size_type
AlignedVector<T>::memory_consumption() const
{
  size_type memory = sizeof(*this);
  for (const T *t = elements.get(); t != used_elements_end; ++t)
    memory += dealii::MemoryConsumption::memory_consumption(*t);
  memory += sizeof(T) * (allocated_elements_end - used_elements_end);
  return memory;
}


#endif // ifndef DOXYGEN


/**
 * 关系运算符 == 用于AlignedVector
 * @relatesalso  AlignedVector
 *
 *
 */
template <class T>
bool
operator==(const AlignedVector<T> &lhs, const AlignedVector<T> &rhs)
{
  if (lhs.size() != rhs.size())
    return false;
  for (typename AlignedVector<T>::const_iterator lit = lhs.begin(),
                                                 rit = rhs.begin();
       lit != lhs.end();
       ++lit, ++rit)
    if (*lit != *rit)
      return false;
  return true;
}



/**
 * 关系运算符！= for AlignedVector
 * @relatesalso  AlignedVector
 *
 *
 */
template <class T>
bool
operator!=(const AlignedVector<T> &lhs, const AlignedVector<T> &rhs)
{
  return !(operator==(lhs, rhs));
}


DEAL_II_NAMESPACE_CLOSE

#endif


