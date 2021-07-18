//include/deal.II-translator/base/mg_level_object_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2020 by the deal.II authors
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

#ifndef dealii_mg_level_object_h
#define dealii_mg_level_object_h

#include <deal.II/base/config.h>

#include <deal.II/base/subscriptor.h>

#include <memory>
#include <vector>

DEAL_II_NAMESPACE_OPEN


/**
 * 这个类代表了一个数组，在多级层次结构中，每个使用的层次都有一个对象，例如用于多栅格算法。与一般的
 * <code>std::vector</code>
 * 相比，这个类只允许存储一些最小和最大索引（=level）之间的对象，因为人们经常希望只在网格的一个子集上运行多层次算法（例如，因为第二或第三层最粗的层次已经足够小，在那里运行直接求解器比递归到更粗的层次要便宜）。尽管只为这些
 * "有趣的
 * "层次存储对象，该类允许简单地按层次进行索引。在内部，这当然是通过简单地将给定的索引移到我们所存储的最小级别来实现的。
 * 在这个类的一个典型的使用案例中，每个层次上存储的对象都是矩阵或向量。
 *
 *
 * @ingroup mg
 *
 * @ingroup data
 *
 */
template <class Object>
class MGLevelObject : public Subscriptor
{
public:
  /**
   * 构造函数。创建一个具有给定最小和最大级别的多级对象，并为
   * <code>maxlevel-minlevel+1</code> 级别的对象分配存储空间。
   * @note
   * 与库中的许多其他地方不同，这里的两个参数并不表示第一级和最后加一级的级别，而是实际上是一个<i>inclusive</i>的级别范围，为级别对象分配存储器。因此，这两个参数的默认值将创建一个有一个级别对象的数组，而不是一个空数组。
   * @param[in]  minlevel 为级别对象提供内存的最低级别。
   * @param[in]  maxlevel 为级别对象提供内存的最高级别。
   * @param[in]  args 传递给底层对象构造器的可选参数。
   * @pre  minlevel <= maxlevel
   *
   */
  template <class... Args>
  MGLevelObject(const unsigned int minlevel,
                const unsigned int maxlevel,
                Args &&... args);

  /**
   * 构造函数。与上述相同，但没有转给底层对象构造器的参数。
   *
   */
  MGLevelObject(const unsigned int minlevel = 0,
                const unsigned int maxlevel = 0);

  /**
   * 访问级别为 @p level. 的对象。
   *
   */
  Object &operator[](const unsigned int level);

  /**
   * 访问 @p level. 级别的对象 这个函数可以在 @p const
   * 对象上调用，并因此返回一个 @p const 引用。
   *
   */
  const Object &operator[](const unsigned int level) const;

  /**
   * 删除此对象以前的所有内容，并根据 @p new_minlevel 和 @p
   * new_maxlevel. 的值重置其大小  @param[in]  new_minlevel
   * 为级别对象提供内存的最低级别。    @param[in]  new_maxlevel
   * 为级别对象提供内存的最高级别。    @param[in]  args
   * 传递给底层对象构造器的可选参数。      @pre  minlevel <=
   * maxlevel
   *
   */
  template <class... Args>
  void
  resize(const unsigned int new_minlevel,
         const unsigned int new_maxlevel,
         Args &&... args);

  /**
   * 对这个对象存储的所有对象调用<tt>operator = (s)</tt>。
   * 这显然要求存储在每个级别上的对象允许这个操作。特别是对于向量和矩阵来说，如果
   * @p d
   * 为零，从而将所有向量或矩阵条目清零，这一点是正确的。
   *
   */
  MGLevelObject<Object> &
  operator=(const double d);

  /**
   * 对这个对象存储的所有对象调用 @p clear
   * 。这个函数只对一些 @p Object
   * 类实现，例如，矩阵类型或PreconditionBlockSOR和类似的类。如果这个类的
   * @p Object 模板类型没有提供 <code>clear()</code>
   * 的成员函数，使用这个函数会出现编译器错误。
   *
   */
  void
  clear_elements();

  /**
   * 该类存储水平对象的最粗略的水平。
   *
   */
  unsigned int
  min_level() const;

  /**
   * 该类存储水平对象的最高水平。
   *
   */
  unsigned int
  max_level() const;

  /**
   * 对存储在这里的每个对象应用 @p action 的动作。参数 @p
   * action 应该是一个接受<code> action(const unsigned int level,
   * Object &object); </code>
   * 语法的函数对象，这意味着这个函数可以接受一个lambda，一个
   * std::function, 或一个普通函数指针。
   *
   */
  template <typename ActionFunctionObjectType>
  void
  apply(ActionFunctionObjectType action);

  /**
   * 这个对象所使用的内存。
   *
   */
  std::size_t
  memory_consumption() const;

private:
  /**
   * 第一个组件的级别。
   *
   */
  unsigned int minlevel;

  /**
   * 要保存的对象的数组。
   *
   */
  std::vector<std::shared_ptr<Object>> objects;
};


 /* ------------------------------------------------------------------- */ 


template <class Object>
template <class... Args>
MGLevelObject<Object>::MGLevelObject(const unsigned int min,
                                     const unsigned int max,
                                     Args &&... args)
  : minlevel(0)
{
  resize(min, max, std::forward<Args>(args)...);
}


template <class Object>
MGLevelObject<Object>::MGLevelObject(const unsigned int min,
                                     const unsigned int max)
  : minlevel(0)
{
  resize(min, max);
}


template <class Object>
Object &MGLevelObject<Object>::operator[](const unsigned int i)
{
  Assert((i >= minlevel) && (i < minlevel + objects.size()),
         ExcIndexRange(i, minlevel, minlevel + objects.size()));
  return *objects[i - minlevel];
}


template <class Object>
const Object &MGLevelObject<Object>::operator[](const unsigned int i) const
{
  Assert((i >= minlevel) && (i < minlevel + objects.size()),
         ExcIndexRange(i, minlevel, minlevel + objects.size()));
  return *objects[i - minlevel];
}


template <class Object>
template <class... Args>
void
MGLevelObject<Object>::resize(const unsigned int new_minlevel,
                              const unsigned int new_maxlevel,
                              Args &&... args)
{
  Assert(new_minlevel <= new_maxlevel, ExcInternalError());
  // note that on clear(), the
  // shared_ptr class takes care of
  // deleting the object it points to
  // by itself
  objects.clear();

  minlevel = new_minlevel;
  for (unsigned int i = 0; i < new_maxlevel - new_minlevel + 1; ++i)
    objects.push_back(std::make_shared<Object>(std::forward<Args>(args)...));
}


template <class Object>
MGLevelObject<Object> &
MGLevelObject<Object>::operator=(const double d)
{
  typename std::vector<std::shared_ptr<Object>>::iterator v;
  for (v = objects.begin(); v != objects.end(); ++v)
    **v = d;
  return *this;
}


template <class Object>
void
MGLevelObject<Object>::clear_elements()
{
  typename std::vector<std::shared_ptr<Object>>::iterator v;
  for (v = objects.begin(); v != objects.end(); ++v)
    (*v)->clear();
}


template <class Object>
unsigned int
MGLevelObject<Object>::min_level() const
{
  return minlevel;
}


template <class Object>
unsigned int
MGLevelObject<Object>::max_level() const
{
  return minlevel + objects.size() - 1;
}

template <class Object>
template <typename ActionFunctionObjectType>
void
MGLevelObject<Object>::apply(ActionFunctionObjectType action)
{
  for (unsigned int lvl = min_level(); lvl <= max_level(); ++lvl)
    {
      action(lvl, (*this)[lvl]);
    }
}


template <class Object>
std::size_t
MGLevelObject<Object>::memory_consumption() const
{
  std::size_t result = sizeof(*this);
  using Iter = typename std::vector<std::shared_ptr<Object>>::const_iterator;
  const Iter end = objects.end();
  for (Iter o = objects.begin(); o != end; ++o)
    result += (*o)->memory_consumption();

  return result;
}

DEAL_II_NAMESPACE_CLOSE

#endif


