//include/deal.II-translator/lac/vector_memory_0.txt
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

#ifndef dealii_vector_memory_h
#define dealii_vector_memory_h


#include <deal.II/base/config.h>

#include <deal.II/base/logstream.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/base/thread_management.h>

#include <deal.II/lac/vector.h>

#include <iostream>
#include <memory>
#include <vector>

DEAL_II_NAMESPACE_OPEN


 /*!@addtogroup VMemory */ 
 /*@{*/ 

/**
 * 向量的内存管理基类。这是一个抽象的基类，除其他地方外，还被所有迭代方法用来为辅助向量分配空间。
 * 这个类的目的如下：在迭代求解器和其他地方，人们需要为向量分配临时存储空间，例如辅助向量。我们可以每次重新分配和释放它们，但是在某些情况下，如果必须非常频繁地这样做，可能会很昂贵。一个常见的情况是，当一个迭代方法被用来在外部求解器的每次迭代中反转一个矩阵时，比如为Schur补码求解器反转一个矩阵块。(例如，
 * step-20
 * 就是这样做的，但它只是永久地保留一个向量作为临时存储。)
 * 在这种情况下，在每次调用内部求解器时重新分配和删除向量是很昂贵的，并且会导致内存碎片化。本类通过提供一个其他类可以用来分配和删除向量的接口来避免这种情况。然后不同的派生类实现不同的策略，为使用类提供临时存储向量。
 * 例如，PrimitiveVectorMemory类简单地通过操作系统设施分配和取消分配向量（即，每次要求向量时使用
 * @p new 和 @p delete)
 * 。对于只被调用一次或很少被调用的迭代求解器来说，这是一个合适的实现。
 * 另一方面，GrowingVectorMemory类在其生命周期内从不向操作系统内存管理子系统返回内存空间；它只将它们标记为未使用，并允许在下次请求矢量时重新使用它们。
 *
 *  <h3> Practical use </h3> 从这个基类派生的类通过
 * VectorMemory::alloc() 函数返回指向新向量的指针，并在通过
 * VectorMemory::free().
 * 返回向量时重新获得向量。因此，这两个函数起到了与 @p
 * new 和 @p delete.
 * 类似的作用，这包括通常的缺点。很容易忘记在使用这一设施的函数的结尾处调用
 * VectorMemory::free() ，或者在函数的 @p if
 * 分支中忘记它，而在这个分支中，人们有一个早期的 @p
 * return
 * 来自该函数。在这两种情况下，这都会导致内存泄漏：一段正确的代码必须在<i>all</i>可能的退出点为所有分配的向量调用
 * VectorMemory::free()
 * 。这包括一个函数因为在调用堆栈的更远处抛出异常而没有在这里明确处理而被留下的地方。
 * 换句话说，通过 VectorMemory::alloc() 分配的向量与通过 @p
 * new:
 * 分配的原始指针有同样的问题，很容易写出有内存泄露的代码。在原始指针的情况下，常见的解决方案是使用
 * std::unique_ptr
 * 类来代替（见http://en.cppreference.com/w/cpp/memory/unique_ptr）。在当前类的情况下，
 * VectorMemory::Pointer
 * 类是解决方案：它是一个在所有实际意义上看起来像指针的类，但在销毁时也将矢量返回到它所获得的VectorMemory对象。由于
 * VectorMemory::Pointer
 * 类的销毁发生在它超出范围的时候（无论是因为函数明确返回，还是因为控制流因异常而离开），内存泄漏不会发生：
 * VectroMemory::Pointer 对象指向的向量被<i>always</i>返回。
 *
 *
 */
template <typename VectorType = dealii::Vector<double>>
class VectorMemory : public Subscriptor
{
public:
  /**
   * 虚拟解构器。这个析构器被声明为 @p virtual
   * ，以允许通过指向该基类的指针销毁派生类型的对象。
   *
   */
  virtual ~VectorMemory() override = default;

  /**
   * 返回一个指向新向量的指针。元素的数量或其细分为块（如果适用）是未指定的，该函数的用户应该将向量重置为其适当的大小。向量的内容也是如此：它们是未指定的。换句话说，调用这个函数的地方将需要适当地调整大小或重新初始化它。
   * @warning  就像在代码中明确使用 <code>new</code> and
   * <code>delete</code>
   * 会招致内存泄露的错误（要么是因为相应的
   * <code>delete</code>
   * 被完全遗忘，要么是因为异常安全问题），明确使用alloc()和free()函数会招致编写代码时意外泄露内存。你应该考虑使用
   * VectorMemory::Pointer 类来代替，它提供了与
   * <code>std::unique</code> 相同的在堆上分配任意内存的服务。
   *
   */
  virtual VectorType *
  alloc() = 0;

  /**
   * 返回一个向量，并表明它不会被调用alloc()获得指针的地方继续使用。
   * @warning  就像在代码中明确使用 <code>new</code> and
   * <code>delete</code>
   * 会招致内存泄漏的错误（要么是因为相应的
   * <code>delete</code>
   * 被完全遗忘，要么是因为异常安全问题），明确使用alloc()和free()函数会招致编写意外泄漏内存的代码。你应该考虑使用
   * VectorMemory::Pointer 类来代替，它提供了与
   * <code>std::unique</code>
   * 相同的服务，用于在堆上分配任意内存。
   *
   */
  virtual void
  free(const VectorType *const) = 0;

  /**
   * @addtogroup  Exceptions 
     * @{ 
   *
   */

  /**
   * 向量不是从这个内存池分配的。
   *
   */
  DeclExceptionMsg(
    ExcNotAllocatedHere,
    "You are trying to deallocate a vector from a memory pool, but this "
    "vector has not actually been allocated by the same pool before.");

  //@}

  /**
   * 一个看起来像指针的类，在构造时从VectorMemory对象（或从VectorMemory派生的类的对象）中分配一个向量，并传递给这个类的构造器。然后，解构器会自动将向量的所有权返回给同一个VectorMemory对象。
   * 因此，这种类型的指针是安全的，因为当它们被销毁时，它们会自动调用
   * VectorMemory::free()
   * ，不管是发生在代码块的末尾还是因为局部变量在异常解除过程中被销毁。因此，这些类型的对象使用户不必明确地使用向量管理函数。
   * 在很多意义上，这个类的行为就像
   * <code>std::unique_ptr</code>
   * 一样，它是一个内存块的唯一所有者，它在销毁时释放了这个内存块。
   * 与 <code>std::unique_ptr</code>
   * 的主要区别是：(i)它在构造时从内存池中分配内存，以及(ii)内存不是用
   * "operator delete "销毁，而是返回到VectorMemory池中。
   *
   */
  class Pointer
    : public std::unique_ptr<VectorType, std::function<void(VectorType *)>>
  {
  public:
    /**
     * 默认构造函数。这个构造函数对应于一个 @p nullptr
     * 对象，它并不拥有一个向量。然而，它以后可以通过移动赋值分配给另一个指针对象，在这种情况下，它将偷取另一个对象所拥有的向量（就像
     * @p std::unique_ptr 那样）。
     *
     */
    Pointer() = default;

    /**
     * 移动构造函数：这通过窃取 @p p.
     * 所拥有的内部数据来创建一个新的Pointer。
     *
     */
    Pointer(Pointer &&p) noexcept = default;

    /**
     * 移动操作符：它释放当前指针所拥有的向量，然后窃取
     * @p p. 所拥有的内部数据。
     *
     */
    Pointer &
    operator=(Pointer &&p) noexcept = default;

    /**
     * 构造函数。这个构造函数自动从给定的向量内存对象中分配一个向量
     * @p mem.  。
     *
     */
    Pointer(VectorMemory<VectorType> &mem);

    /**
     * 解构器，自动从内存池中释放向量。
     *
     */
    ~Pointer() = default;
  };
};



/**
 * 简单的内存管理。参见基类的文档以了解其目的。
 * 这个类根据需要从全局堆中分配和删除向量，即不执行任何特别适应的内存管理的动作。
 *
 *
 */
template <typename VectorType = dealii::Vector<double>>
class PrimitiveVectorMemory : public VectorMemory<VectorType>
{
public:
  /**
   * 返回一个指向新向量的指针。元素的数量或其细分为块（如果适用）是未指定的，该函数的用户应该将向量重置为其适当的大小。向量的内容也是如此：它们是未指定的。换句话说，调用这个函数的地方需要适当地调整其大小或重新初始化。
   * 对于本类，调用这个函数将在堆上分配一个新的向量，并返回一个指向它的指针。之后调用free()，然后将内存返回到由操作系统管理的全局堆中。
   * @warning  就像在代码中明确使用 <code>new</code> and
   * <code>delete</code>
   * 会招致内存泄露的错误（要么是因为相应的
   * <code>delete</code>
   * 被完全遗忘，要么是因为异常安全问题），明确使用alloc()和free()函数会招致编写代码时意外泄露内存。你应该考虑使用
   * VectorMemory::Pointer 类来代替，它提供了与
   * <code>std::unique</code> 相同的在堆上分配任意内存的服务。
   *
   */
  virtual VectorType *
  alloc() override;

  /**
   * 返回一个向量，并表明它不会被调用alloc()获得指针的实例继续使用。
   * 对于本类来说，这意味着向量被返回到全局堆中。
   * @warning  就像在代码中明确使用 <code>new</code> and
   * <code>delete</code>
   * 会招致内存泄露的错误（要么是因为相应的
   * <code>delete</code>
   * 被完全遗忘，要么是因为异常安全问题），明确使用
   * alloc() 和 free()
   * 函数会招致编写意外泄露内存的代码。你应该考虑使用
   * VectorMemory::Pointer 类来代替，它提供了与
   * <code>std::unique</code> 相同的在堆上分配任意内存的服务。
   *
   */
  virtual void
  free(const VectorType *const v) override;
};



/**
 * 一个基于池的内存管理类。参见基类的文档以了解其用途。
 * 每次从这个类中请求一个向量，它都会检查它是否有一个可用的向量并返回其地址，或者在堆上分配一个新的向量。如果一个向量通过
 * GrowingVectorMemory::free()
 * 成员函数从其用户处返回，它不会将分配的内存返回到操作系统的内存子系统中，而是将其保留在周围未使用的地方，以便在以后再次调用
 * GrowingVectorMemory::alloc()
 * 时使用。因此，如果经常需要和释放临时向量，该类可以避免在堆上重复分配内存的开销；另一方面，它不会在尽可能早的时间释放曾经分配过的内存，因此可能导致整体内存消耗增加。
 * 同一向量类型的所有GrowingVectorMemory对象都使用同一个内存池。
 * (换句话说。这个类所抽取的向量池是<i>global</i>，而不是当前类的一个普通成员变量，该成员变量在周围的GrowingVectorMemory对象被销毁时被销毁）。)
 * 因此，函数可以在任何需要的时候创建这样一个GrowingVectorMemory对象，而不需要每次都创建一个新的内存池的性能惩罚。这个策略的一个缺点是，一旦分配了向量，只有在程序运行结束时才会释放。
 *
 *
 */
template <typename VectorType = dealii::Vector<double>>
class GrowingVectorMemory : public VectorMemory<VectorType>
{
public:
  /**
   * 声明容器大小的类型。
   *
   */
  using size_type = types::global_dof_index;

  /**
   * 构造器。
   * 该参数允许预先分配一定数量的向量。默认情况下是不这样做的。
   *
   */
  GrowingVectorMemory(const size_type initial_size   = 0,
                      const bool      log_statistics = false);

  /**
   * 解构器。解构器也会检查所有通过当前对象分配的向量是否都已被再次释放。
   * 然而，正如在类的文档中所讨论的那样，这并不意味着它们的内存被返回到操作系统中。
   *
   */
  virtual ~GrowingVectorMemory() override;

  /**
   * 返回一个指向新向量的指针。元素的数量或它们被细分为块（如果适用）是没有指定的，这个函数的用户应该将向量重置为适当的大小。向量的内容也是如此：它们是未指定的。换句话说，调用这个函数的地方将需要适当地调整大小或重新初始化它。
   * @warning  就像在代码中明确使用 <code>new</code> and
   * <code>delete</code>
   * 会招致内存泄露的错误（要么是因为相应的
   * <code>delete</code>
   * 被完全遗忘，要么是因为异常安全问题），明确使用alloc()和free()函数会招致编写代码时意外泄露内存。你应该考虑使用
   * VectorMemory::Pointer 类来代替，它提供了与
   * <code>std::unique</code>
   * 相同的服务，用于在堆上分配任意内存。
   *
   */
  virtual VectorType *
  alloc() override;

  /**
   * 返回一个向量，并表明它不会被调用alloc()来获取指针的实例继续使用。
   * 对于本类来说，这意味着保留该向量以便以后由alloc()方法重新使用。
   * @warning  就像在代码中明确使用 <code>new</code> and
   * <code>delete</code>
   * 会招致内存泄露的错误（要么是因为相应的
   * <code>delete</code>
   * 被完全遗忘，要么是因为异常安全问题），明确使用alloc()和free()函数会招致编写意外泄露内存的代码。你应该考虑使用
   * VectorMemory::Pointer 类来代替，它提供了与
   * <code>std::unique</code>
   * 相同的服务，用于在堆上分配任意内存。
   *
   */
  virtual void
  free(const VectorType *const) override;

  /**
   * 释放所有当前不使用的向量。
   *
   */
  static void
  release_unused_memory();

  /**
   * 本类和所有当前分配的向量所消耗的内存。
   *
   */
  virtual std::size_t
  memory_consumption() const;

private:
  /**
   * 一个描述代表此对象存储的向量的数组的这个条目的类型。这一对的第一个分量是一个标志，告诉人们是否使用该向量，第二个分量是指向该向量本身的指针。
   *
   */
  using entry_type = std::pair<bool, std::unique_ptr<VectorType>>;

  /**
   * 为内存池提供实际存储的类。
   * 这是提供GrowingVectorMemory的实际存储的地方。
   * 每个向量类型只使用一个内存池，因此从同一个存储空间分配所有的向量。
   *
   */
  struct Pool
  {
    /**
     * 标准构造函数创建一个空池
     *
     */
    Pool();

    /**
     * 解构器。
     *
     */
    ~Pool();

    /**
     * 创建数据向量；第一次初始化后不做任何事情
     *
     */
    void
    initialize(const size_type size);

    /**
     * 指向存储对象的指针
     *
     */
    std::vector<entry_type> *data;
  };

  /**
   * 返回一个分配好的向量数组。
   *
   */
  static Pool &
  get_pool();

  /**
   * 分配的总数量。只用于记账和在对象的生命周期结束时生成输出。
   *
   */
  size_type total_alloc;

  /**
   * 当前在此对象中分配的向量的数量；用于检测内存泄漏。
   *
   */
  size_type current_alloc;

  /**
   * 一个控制由析构器记录统计数据的标志。
   *
   */
  bool log_statistics;

  /**
   * Mutex，用于从多个线程同步访问此对象的内部数据。
   *
   */
  static Threads::Mutex mutex;
};



namespace internal
{
  namespace GrowingVectorMemoryImplementation
  {
    void
    release_all_unused_memory();
  }
} // namespace internal

 /*@}*/ 

#ifndef DOXYGEN
 /* --------------------- inline functions ---------------------- */ 


template <typename VectorType>
inline VectorMemory<VectorType>::Pointer::Pointer(VectorMemory<VectorType> &mem)
  : std::unique_ptr<VectorType, std::function<void(VectorType *)>>(
      mem.alloc(),
      [&mem](VectorType *v) { mem.free(v); })
{}



template <typename VectorType>
VectorType *
PrimitiveVectorMemory<VectorType>::alloc()
{
  return new VectorType();
}



template <typename VectorType>
void
PrimitiveVectorMemory<VectorType>::free(const VectorType *const v)
{
  delete v;
}



#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif


