//include/deal.II-translator/particles/particle_iterator_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2020 by the deal.II authors
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

#ifndef dealii_particles_particle_iterator_h
#define dealii_particles_particle_iterator_h

#include <deal.II/base/config.h>

#include <deal.II/particles/particle_accessor.h>

DEAL_II_NAMESPACE_OPEN

namespace Particles
{
  // Forward declaration
#ifndef DOXYGEN
  template <int, int>
  class ParticleHandler;
#endif

  /**
   * 一个用于迭代粒子的类。与ParticleAccessor类一起，它被用来隐藏粒子类和粒子容器的内部实现。
   *
   */
  template <int dim, int spacedim = dim>
  class ParticleIterator
  {
  public:
    /**
     * 空的构造函数。这样的对象是不能使用的!
     *
     */
    ParticleIterator() = default;

    /**
     * 迭代器的构造函数。接受一个对粒子容器的引用，以及一个对细胞-粒子对的迭代器。
     *
     */
    ParticleIterator(
      const std::multimap<internal::LevelInd, Particle<dim, spacedim>> &map,
      const typename std::multimap<internal::LevelInd,
                                   Particle<dim, spacedim>>::iterator
        &particle);

    /**
     * 解除引用操作符，返回一个访问器的引用。因此，用法类似于<tt>(*i).get_id
     * ();</tt>。
     *
     */
    const ParticleAccessor<dim, spacedim> &operator*() const;

    /**
     * 去引用操作符，非 @p const 版本。
     *
     */
    ParticleAccessor<dim, spacedim> &operator*();

    /**
     * 解除引用操作符，返回所指向的粒子的指针。
     * 因此，用法类似于 <tt>i->get_id ();</tt> 有一个  @p const
     * 和一个非  @p const  版本。
     *
     */
    const ParticleAccessor<dim, spacedim> *operator->() const;

    /**
     * 去引用操作符，非 @p const 版本。
     *
     */
    ParticleAccessor<dim, spacedim> *operator->();

    /**
     * 等价比较。
     *
     */
    bool
    operator==(const ParticleIterator<dim, spacedim> &) const;

    /**
     * 不等式比较。
     *
     */
    bool
    operator!=(const ParticleIterator<dim, spacedim> &) const;

    /**
     * 前缀<tt>++</tt>运算符。<tt>++iterator</tt>。该操作符将迭代器推进到下一个元素，并返回一个对<tt>*this</tt>的引用。
     *
     */
    ParticleIterator &
    operator++();

    /**
     * 后缀<tt>++</tt>操作符。<tt>iterator++</tt>。该操作符将迭代器推进到下一个元素，但返回之前指向的元素的迭代器。
     *
     */
    ParticleIterator
    operator++(int);

    /**
     * 前缀<tt>\--</tt>操作符。<tt>\--iterator</tt>。这个操作符将迭代器移到前一个元素，并返回一个对<tt>*this</tt>的引用。
     *
     */
    ParticleIterator &
    operator--();

    /**
     * 后缀<tt>\--</tt>操作符。<tt>iterator\--</tt>。这个操作符将迭代器移动到前一个元素，但返回一个迭代器到之前指向的元素。
     *
     */
    ParticleIterator
    operator--(int);

    /**
     * 将该类标记为双向迭代器，并声明一些别名，这些别名是迭代器的标准，被算法用来查询它们所工作的迭代器的具体内容。
     *
     */
    using iterator_category = std::bidirectional_iterator_tag;
    using value_type        = ParticleAccessor<dim, spacedim>;
    using difference_type   = std::ptrdiff_t;
    using pointer           = ParticleAccessor<dim, spacedim> *;
    using reference         = ParticleAccessor<dim, spacedim> &;

  private:
    /**
     * 对实际粒子的访问器。
     *
     */
    ParticleAccessor<dim, spacedim> accessor;
  };



  // ------------------------------ inline functions -------------------------

  template <int dim, int spacedim>
  inline ParticleIterator<dim, spacedim>::ParticleIterator(
    const std::multimap<internal::LevelInd, Particle<dim, spacedim>> &map,
    const typename std::multimap<internal::LevelInd,
                                 Particle<dim, spacedim>>::iterator & particle)
    : accessor(map, particle)
  {}



  template <int dim, int spacedim>
  inline ParticleAccessor<dim, spacedim> &ParticleIterator<dim, spacedim>::
                                          operator*()
  {
    return accessor;
  }



  template <int dim, int spacedim>
  inline ParticleAccessor<dim, spacedim> *ParticleIterator<dim, spacedim>::
                                          operator->()
  {
    return &(this->operator*());
  }



  template <int dim, int spacedim>
  inline const ParticleAccessor<dim, spacedim> &
    ParticleIterator<dim, spacedim>::operator*() const
  {
    return accessor;
  }



  template <int dim, int spacedim>
  inline const ParticleAccessor<dim, spacedim> *
    ParticleIterator<dim, spacedim>::operator->() const
  {
    return &(this->operator*());
  }



  template <int dim, int spacedim>
  inline bool
  ParticleIterator<dim, spacedim>::
  operator!=(const ParticleIterator<dim, spacedim> &other) const
  {
    return accessor != other.accessor;
  }



  template <int dim, int spacedim>
  inline bool
  ParticleIterator<dim, spacedim>::
  operator==(const ParticleIterator<dim, spacedim> &other) const
  {
    return accessor == other.accessor;
  }



  template <int dim, int spacedim>
  inline ParticleIterator<dim, spacedim> &
  ParticleIterator<dim, spacedim>::operator++()
  {
    accessor.next();
    return *this;
  }



  template <int dim, int spacedim>
  inline ParticleIterator<dim, spacedim>
  ParticleIterator<dim, spacedim>::operator++(int)
  {
    ParticleIterator tmp(*this);
                     operator++();

    return tmp;
  }



  template <int dim, int spacedim>
  inline ParticleIterator<dim, spacedim> &
  ParticleIterator<dim, spacedim>::operator--()
  {
    accessor.prev();
    return *this;
  }



  template <int dim, int spacedim>
  inline ParticleIterator<dim, spacedim>
  ParticleIterator<dim, spacedim>::operator--(int)
  {
    ParticleIterator tmp(*this);
                     operator--();

    return tmp;
  }


} // namespace Particles

DEAL_II_NAMESPACE_CLOSE

#endif


