//include/deal.II-translator/base/smartpointer_0.txt
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

#ifndef dealii_smartpointer_h
#define dealii_smartpointer_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>

#include <atomic>
#include <typeinfo>

DEAL_II_NAMESPACE_OPEN

/**
 * 智能指针避免使用悬空的指针。它们可以像指针一样被使用（即使用<tt>*</tt>和<tt>-></tt>操作符和通过铸造），但要确保在使用指针的过程中，被指向的对象不会被删除或移动，要给被指向者发出使用信号。
 * 指向的对象，即类T，应该继承Subscriptor或者必须实现相同的功能。空指针是这个规则的一个例外，也是允许的。
 * 第二个模板参数P只有一个作用：如果使用了一个没有调试字符串的构造函数，那么P的名字就被用作调试字符串。
 * SmartPointer没有实现任何内存处理！尤其是删除一个SmartPointer。特别是，删除一个SmartPointer并不会删除该对象。写作
 *
 * @code
 * SmartPointer<T,P> dont_do_this = new T;
 * @endcode
 * 是一个肯定的方法来编写一个内存泄漏的程序!
 * 安全的版本是
 *
 * @code
 * T* p = new T;
 * {
 * SmartPointer<T,P> t(p);
 * ...
 * }
 * delete p;
 * @endcode
 *
 * 注意智能指针可以处理对象的<tt>const</tt>性，即一个<tt>SmartPointer<const
 * ABC></tt>真的表现得像一个指向常量对象的指针（取消引用时不允许写访问），而<tt>SmartPointer<ABC></tt>是一个可变的指针。
 *
 *
 * @ingroup memory
 *
 *
 */
template <typename T, typename P = void>
class SmartPointer
{
public:
  /**
   * 空指针的标准构造函数。这个指针的id被设置为类P的名称。
   *
   */
  SmartPointer();

  /**
   * SmartPointer的复制构造函数。我们不复制从<tt>tt</tt>订阅的对象，而是再次订阅自己的对象。
   *
   */
  template <class Q>
  SmartPointer(const SmartPointer<T, Q> &tt);

  /**
   * SmartPointer的复制构造函数。我们不复制从<tt>tt</tt>订阅的对象，而是自己再次订阅它。
   *
   */
  SmartPointer(const SmartPointer<T, P> &tt);

  /**
   * 构造函数接受一个普通的指针。如果可能的话，也就是说，如果该指针不是空指针，构造函数会订阅给定的对象以锁定它，也就是说，防止它在使用结束前被破坏。
   * <tt>id</tt>在调用 Subscriptor::subscribe(id) 时被使用，在调用
   * Subscriptor::unsubscribe(). 时被~SmartPointer()使用。
   *
   */
  SmartPointer(T *t, const std::string &id);

  /**
   * 构造函数取一个正常的指针。如果可能的话，即如果该指针不是空指针，构造函数订阅给定的对象以锁定它，即防止它在使用结束前被破坏。这个指针的id被设置为类P的名称。
   *
   */
  SmartPointer(T *t);

  /**
   * 销毁器，删除订阅。
   *
   */
  ~SmartPointer();

  /**
   * 普通指针的赋值操作符。指针会自动订阅新的对象，如果旧的对象存在，则取消订阅。它不会尝试订阅一个空指针，但仍会删除旧的订阅。
   *
   */
  SmartPointer<T, P> &
  operator=(T *tt);

  /**
   * SmartPointer的赋值操作符。指针会自动订阅新的对象，如果存在旧的对象，则会取消订阅。
   *
   */
  template <class Q>
  SmartPointer<T, P> &
  operator=(const SmartPointer<T, Q> &tt);

  /**
   * 用于SmartPointer的赋值操作符。指针会自动订阅新的对象，如果旧的对象存在，则会取消订阅。
   *
   */
  SmartPointer<T, P> &
  operator=(const SmartPointer<T, P> &tt);

  /**
   * 删除所指向的对象并将指针设为零。
   *
   */
  void
  clear();

  /**
   * 转换为普通指针。
   *
   */
  operator T *() const;

  /**
   * 解除引用操作符。如果指针是一个空指针，这个操作符会抛出一个ExcNotInitialized()。
   *
   */
  T &operator*() const;

  /**
   * 解除引用操作符。如果指针是一个空指针，该操作符会抛出一个ExcNotInitializedi()。
   *
   */
  T *operator->() const;

  /**
   * 交换这个对象和参数的指针。由于被指向的两个对象在之前和之后都被订阅了，我们不必改变它们的订阅计数器。
   * 请注意，这个函数（有两个参数）以及其中一个参数是指针，另一个参数是C型指针的相应函数是在全局命名空间实现的。
   *
   */
  template <class Q>
  void
  swap(SmartPointer<T, Q> &tt);

  /**
   * 在这个对象和给出的指针之间交换指针。由于这释放了目前所指向的对象，我们把它的订阅数减少了一个，并在我们将来要指向的对象处增加订阅数。
   * 注意，我们确实需要一个指针的引用，因为我们想改变我们所给的指针变量。
   *
   */
  void
  swap(T *&tt);

  /**
   * 返回这个类所使用的内存量的估计值（以字节为单位）。
   * 特别注意，这只包括<b>this</b>对象所使用的内存量，而不是所指向的对象。
   *
   */
  std::size_t
  memory_consumption() const;

private:
  /**
   * 指向我们要订阅的对象的指针。由于在调试时经常需要跟踪这个指针，我们特意选择了一个简短的名字。
   *
   */
  T *t;

  /**
   * 用于下标的标识。
   *
   */
  const std::string id;

  /**
   * 当所指向的对象被销毁或被移出时，Smartpointer就会被废止。
   *
   */
  std::atomic<bool> pointed_to_object_is_alive;
};


 /* --------------------- inline Template functions ------------------------- */ 


template <typename T, typename P>
inline SmartPointer<T, P>::SmartPointer()
  : t(nullptr)
  , id(typeid(P).name())
  , pointed_to_object_is_alive(false)
{}



template <typename T, typename P>
inline SmartPointer<T, P>::SmartPointer(T *t)
  : t(t)
  , id(typeid(P).name())
  , pointed_to_object_is_alive(false)
{
  if (t != nullptr)
    t->subscribe(&pointed_to_object_is_alive, id);
}



template <typename T, typename P>
inline SmartPointer<T, P>::SmartPointer(T *t, const std::string &id)
  : t(t)
  , id(id)
  , pointed_to_object_is_alive(false)
{
  if (t != nullptr)
    t->subscribe(&pointed_to_object_is_alive, id);
}



template <typename T, typename P>
template <class Q>
inline SmartPointer<T, P>::SmartPointer(const SmartPointer<T, Q> &tt)
  : t(tt.t)
  , id(tt.id)
  , pointed_to_object_is_alive(false)
{
  if (tt.pointed_to_object_is_alive && t != nullptr)
    t->subscribe(&pointed_to_object_is_alive, id);
}



template <typename T, typename P>
inline SmartPointer<T, P>::SmartPointer(const SmartPointer<T, P> &tt)
  : t(tt.t)
  , id(tt.id)
  , pointed_to_object_is_alive(false)
{
  if (tt.pointed_to_object_is_alive && t != nullptr)
    t->subscribe(&pointed_to_object_is_alive, id);
}



template <typename T, typename P>
inline SmartPointer<T, P>::~SmartPointer()
{
  if (pointed_to_object_is_alive && t != nullptr)
    t->unsubscribe(&pointed_to_object_is_alive, id);
}



template <typename T, typename P>
inline void
SmartPointer<T, P>::clear()
{
  if (pointed_to_object_is_alive && t != nullptr)
    {
      t->unsubscribe(&pointed_to_object_is_alive, id);
      delete t;
      Assert(pointed_to_object_is_alive == false, ExcInternalError());
    }
  t = nullptr;
}



template <typename T, typename P>
inline SmartPointer<T, P> &
SmartPointer<T, P>::operator=(T *tt)
{
  // optimize if no real action is
  // requested
  if (t == tt)
    return *this;

  if (pointed_to_object_is_alive && t != nullptr)
    t->unsubscribe(&pointed_to_object_is_alive, id);
  t = tt;
  if (tt != nullptr)
    tt->subscribe(&pointed_to_object_is_alive, id);
  return *this;
}



template <typename T, typename P>
template <class Q>
inline SmartPointer<T, P> &
SmartPointer<T, P>::operator=(const SmartPointer<T, Q> &tt)
{
  // if objects on the left and right
  // hand side of the operator= are
  // the same, then this is a no-op
  if (&tt == this)
    return *this;

  if (pointed_to_object_is_alive && t != nullptr)
    t->unsubscribe(&pointed_to_object_is_alive, id);
  t = static_cast<T *>(tt);
  if (tt.pointed_to_object_is_alive && tt != nullptr)
    tt->subscribe(&pointed_to_object_is_alive, id);
  return *this;
}



template <typename T, typename P>
inline SmartPointer<T, P> &
SmartPointer<T, P>::operator=(const SmartPointer<T, P> &tt)
{
  // if objects on the left and right
  // hand side of the operator= are
  // the same, then this is a no-op
  if (&tt == this)
    return *this;

  if (pointed_to_object_is_alive && t != nullptr)
    t->unsubscribe(&pointed_to_object_is_alive, id);
  t = static_cast<T *>(tt);
  if (tt.pointed_to_object_is_alive && tt != nullptr)
    tt->subscribe(&pointed_to_object_is_alive, id);
  return *this;
}



template <typename T, typename P>
inline SmartPointer<T, P>::operator T *() const
{
  return t;
}



template <typename T, typename P>
inline T &SmartPointer<T, P>::operator*() const
{
  Assert(t != nullptr, ExcNotInitialized());
  Assert(pointed_to_object_is_alive,
         ExcMessage("The object pointed to is not valid anymore."));
  return *t;
}



template <typename T, typename P>
inline T *SmartPointer<T, P>::operator->() const
{
  Assert(t != nullptr, ExcNotInitialized());
  Assert(pointed_to_object_is_alive,
         ExcMessage("The object pointed to is not valid anymore."));
  return t;
}



template <typename T, typename P>
template <class Q>
inline void
SmartPointer<T, P>::swap(SmartPointer<T, Q> &tt)
{
#ifdef DEBUG
  SmartPointer<T, P> aux(t, id);
  *this = tt;
  tt    = aux;
#else
  std::swap(t, tt.t);
#endif
}



template <typename T, typename P>
inline void
SmartPointer<T, P>::swap(T *&tt)
{
  if (pointed_to_object_is_alive && t != nullptr)
    t->unsubscribe(pointed_to_object_is_alive, id);

  std::swap(t, tt);

  if (t != nullptr)
    t->subscribe(pointed_to_object_is_alive, id);
}



template <typename T, typename P>
inline std::size_t
SmartPointer<T, P>::memory_consumption() const
{
  return sizeof(SmartPointer<T, P>);
}



// The following function is not strictly necessary but is an optimization
// for places where you call swap(p1,p2) with SmartPointer objects p1, p2.
// Unfortunately, MS Visual Studio (at least up to the 2013 edition) trips
// over it when calling std::swap(v1,v2) where v1,v2 are std::vectors of
// SmartPointer objects: it can't determine whether it should call std::swap
// or dealii::swap on the individual elements (see bug #184 on our Google Code
// site. Consequently, just take this function out of the competition for this
// compiler.
#ifndef _MSC_VER
/**
 * 全局性的函数来交换两个智能指针的内容。由于指针所指向的两个对象都保留了订阅，我们不必改变它们的订阅数。
 *
 *
 */
template <typename T, typename P, class Q>
inline void
swap(SmartPointer<T, P> &t1, SmartPointer<T, Q> &t2)
{
  t1.swap(t2);
}
#endif


/**
 * 全局函数来交换一个智能指针和一个C型指针的内容。
 * 注意，我们确实需要一个指针的引用，因为我们要改变我们所给的指针变量。
 *
 *
 */
template <typename T, typename P>
inline void
swap(SmartPointer<T, P> &t1, T *&t2)
{
  t1.swap(t2);
}



/**
 * 全局函数，用于交换一个C-style指针和一个智能指针的内容。
 * 注意，我们确实需要一个指针的引用，因为我们想改变我们所给的指针变量。
 *
 *
 */
template <typename T, typename P>
inline void
swap(T *&t1, SmartPointer<T, P> &t2)
{
  t2.swap(t1);
}

DEAL_II_NAMESPACE_CLOSE

#endif


