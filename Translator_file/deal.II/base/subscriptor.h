//include/deal.II-translator/base/subscriptor_0.txt
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

#ifndef dealii_subscriptor_h
#define dealii_subscriptor_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>

#include <atomic>
#include <cstring>
#include <map>
#include <mutex>
#include <string>
#include <typeinfo>
#include <vector>

DEAL_II_NAMESPACE_OPEN

/**
 * 对订阅的处理。
 * 这个类，作为一个基类，允许使用一个特定的对象来跟踪其他对象。它被用来避免指向从Subscriptor派生的类的对象的指针在该对象被无效后被引用。在这里，无效化被认为是在对象被移出或销毁时发生的。该机制的工作方式如下。成员函数subscribe()接受一个指向布尔值的指针，该指针在失效时被修改。拥有这个指针的对象（通常是SmartPointer类型的对象）然后被期望在试图访问这个类之前检查布尔值的状态。
 * 通过为函数subscribe()和unsubscribe()提供识别字符串，这个类的效用甚至得到了加强。这些字符串被表示为
 * <code>const char</code>
 * 指针，因为底层的缓冲区来自（并由）运行时类型信息系统管理：更确切地说，这些指针是函数调用
 * <code>typeid(x).name()</code> 的结果，其中 <code>x</code>
 * 是一些对象。因此，提供给 subscribe() 和 unsubscribe()
 * 的指针必须相同。内容相同的字符串将不会被识别为相同。SmartPointer中的处理方法将照顾到这一点。该类的当前订阅者可以通过调用list_subscribers()获得。
 *
 *
 * @ingroup memory
 *
 *
 */
class Subscriptor
{
public:
  /**
   * 构造函数将计数器设置为零。
   *
   */
  Subscriptor();

  /**
   * 复制构造器。
   * 复制的计数器是零，因为引用指向原始对象。
   *
   */
  Subscriptor(const Subscriptor &);

  /**
   * 移动构造函数。
   * 继承自Subscriptor的对象只有在没有其他对象订阅它的情况下才能被移动。
   *
   */
  Subscriptor(Subscriptor &&) noexcept;

  /**
   * 解构器，断言计数器为零。
   *
   */
  virtual ~Subscriptor();

  /**
   * 赋值操作符。
   * 这也要小心处理，因为计数器必须保持不变。因此，它只是返回<tt>*this</tt>而已。
   *
   */
  Subscriptor &
  operator=(const Subscriptor &);

  /**
   * 移动赋值运算符。只会使所移动的对象无效。
   *
   */
  Subscriptor &
  operator=(Subscriptor &&) noexcept;

  /**
   * @name Subscriptor功能 从Subscriptor派生的类提供了一个订阅此对象的设施。这主要是由SmartPointer类使用。
   *
   */
  // @{

  /**
   * 通过存储指针来订阅该对象的用户  @p validity.
   * 订阅者可以通过提供的文本来识别  @p identifier.  。
   *
   */
  void
  subscribe(std::atomic<bool> *const validity,
            const std::string &      identifier = "") const;

  /**
   * 从对象中取消用户的订阅。
   * @note   @p identifier 和 @p validity
   * 的指针必须与提供给subscribe()的指针相同。
   *
   */
  void
  unsubscribe(std::atomic<bool> *const validity,
              const std::string &      identifier = "") const;

  /**
   * 返回目前对这个对象的订阅数量。这允许使用这个类来确定引用计数的寿命，其中最后一个取消订阅的人也会删除该对象。
   *
   */
  unsigned int
  n_subscriptions() const;

  /**
   * 列出输入的订阅者  @p stream.  。
   *
   */
  template <typename StreamType>
  void
  list_subscribers(StreamType &stream) const;

  /**
   * 列出输入 @p deallog. 的订阅者。
   *
   */
  void
  list_subscribers() const;

  // @}

  /**
   * @addtogroup  Exceptions  
     * @{ 
   *
   */

  /**
   * 异常情况。对象不能被删除，因为它被使用了。
   *
   */
  DeclException3(ExcInUse,
                 int,
                 std::string,
                 std::string,
                 << "Object of class " << arg2 << " is still used by " << arg1
                 << " other objects."
                 << "\n\n"
                 << "(Additional information: " << arg3 << ")\n\n"
                 << "See the entry in the Frequently Asked Questions of "
                 << "deal.II (linked to from http://www.dealii.org/) for "
                 << "a lot more information on what this error means and "
                 << "how to fix programs in which it happens.");

  /**
   * 具有 Subscriptor::unsubscribe()
   * 给出的识别字符串的订阅者没有订阅该对象。
   *
   */
  DeclException2(ExcNoSubscriber,
                 std::string,
                 std::string,
                 << "No subscriber with identifier <" << arg2
                 << "> subscribes to this object of class " << arg1
                 << ". Consequently, it cannot be unsubscribed.");
  //@}

  /**
   * 为了使用[BOOST序列化库](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html)进行序列化，将此对象的数据读入或写入一个流中。
   * 这个函数实际上并没有对这个类的任何成员变量进行序列化。原因是这个类所存储的只是谁订阅了这个对象，但在存储这个对象的内容时，谁订阅了这个对象，与恢复时谁订阅了这个对象并不一定有关系。因此，我们不希望在恢复时覆盖订阅者，那么就没有理由在一开始就把订阅者写出来。
   *
   */
  template <class Archive>
  void
  serialize(Archive &ar, const unsigned int version);

private:
  /**
   * 存储订阅此对象的对象的数量。最初，这个数字是零，在销毁时，它应再次为零（即所有订阅过的对象应再次取消订阅）。
   * 如果一个对象的创建者（和所有者）能够提供身份证明，那么他将被计算在下面的地图中。
   * 我们使用<tt>mutable</tt>关键字，以便也允许订阅常量对象。
   * 这个计数器可以在多线程代码中被同时读出和写入：因此我们使用了
   * <code>std::atomic</code> 类模板。
   *
   */
  mutable std::atomic<unsigned int> counter;

  /**
   * 在这个地图中，我们为提供给subscribe()的每个不同的标识字符串计算订阅数。
   *
   */
  mutable std::map<std::string, unsigned int> counter_map;

  /**
   * 在#counter_map中使用的数据类型。
   *
   */
  using map_value_type = decltype(counter_map)::value_type;

  /**
   * 在#counter_map中使用的迭代器类型。
   *
   */
  using map_iterator = decltype(counter_map)::iterator;

  /**
   * 在这个向量中，我们在订阅这个类的SmartPointer对象中存储指向有效性bool的指针。
   *
   */
  mutable std::vector<std::atomic<bool> *> validity_pointers;

  /**
   * 指向这个对象的typeinfo对象的指针，以后我们可以从中推导出类的名称。因为这个派生类的信息既不能在析构器中获得，也不能在构造器中获得，所以我们在两者之间获得它，并把它存储在这里。
   *
   */
  mutable const std::type_info *object_info;

  /**
   * 检查是否有任何对象订阅了这个对象。如果这个检查通过，那么销毁当前对象是安全的。如果这个检查失败，那么这个函数将中止或向deallog打印错误信息（通过使用AssertNothrow机制），但不会抛出一个异常。
   * @note
   * 因为这个函数只是一个一致性检查，所以它在释放模式下没有任何作用。
   * @note
   * 如果这个函数在有一个未捕获的异常时被调用，那么，这个函数不是中止，而是向标准错误流打印一个错误信息并返回。
   *
   */
  void
  check_no_subscribers() const noexcept;

  /**
   * 一个突变器，用于在打印出订阅者列表时确保数据的一致性。
   *
   */
  static std::mutex mutex;
};

//---------------------------------------------------------------------------

inline Subscriptor::Subscriptor()
  : counter(0)
  , object_info(nullptr)
{}



inline Subscriptor::Subscriptor(const Subscriptor &)
  : counter(0)
  , object_info(nullptr)
{}



inline Subscriptor &
Subscriptor::operator=(const Subscriptor &s)
{
  object_info = s.object_info;
  return *this;
}



inline unsigned int
Subscriptor::n_subscriptions() const
{
  return counter;
}



template <class Archive>
inline void
Subscriptor::serialize(Archive &, const unsigned int)
{
  // do nothing, as explained in the
  // documentation of this function
}

template <typename StreamType>
inline void
Subscriptor::list_subscribers(StreamType &stream) const
{
  std::lock_guard<std::mutex> lock(mutex);

  for (const auto &it : counter_map)
    stream << it.second << '/' << counter << " subscriptions from \""
           << it.first << '\"' << std::endl;
}

DEAL_II_NAMESPACE_CLOSE

#endif


