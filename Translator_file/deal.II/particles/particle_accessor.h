//include/deal.II-translator/particles/particle_accessor_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2021 by the deal.II authors
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

#ifndef dealii_particles_particle_accessor_h
#define dealii_particles_particle_accessor_h

#include <deal.II/base/config.h>

#include <deal.II/base/array_view.h>

#include <deal.II/grid/tria.h>

#include <deal.II/particles/particle.h>

DEAL_II_NAMESPACE_OPEN

namespace Particles
{
  // Forward declarations
#ifndef DOXYGEN
  template <int, int>
  class ParticleIterator;
  template <int, int>
  class ParticleHandler;
#endif

  /**
   * ParticleIterator使用的访问器类，用于访问粒子数据。
   *
   */
  template <int dim, int spacedim = dim>
  class ParticleAccessor
  {
  public:
    /**
     * @copydoc   Particle::write_particle_data_to_memory .
     *
     */
    void *
    write_particle_data_to_memory(void *data) const;


    /**
     * @copydoc   Particle::read_particle_data_from_memory .
     *
     */
    const void *
    read_particle_data_from_memory(const void *data);

    /**
     * 设置这个粒子的位置。注意，这并不检查这是否是模拟域中的有效位置。
     * @param  [in] new_location 这个粒子的新位置。
     *
     */
    void
    set_location(const Point<spacedim> &new_location);

    /**
     * 获取这个粒子的位置。          @return
     * 这个粒子的位置。
     *
     */
    const Point<spacedim> &
    get_location() const;

    /**
     * 设置此粒子的参考位置。          @param  [in]
     * new_reference_location 这个粒子的新参考位置。
     *
     */
    void
    set_reference_location(const Point<dim> &new_reference_location);

    /**
     * 返回这个粒子在其当前单元中的参考位置。
     *
     */
    const Point<dim> &
    get_reference_location() const;

    /**
     * 返回这个粒子的ID号。
     *
     */
    types::particle_index
    get_id() const;

    /**
     * 告诉粒子在哪里存储它的属性（即使它不拥有属性）。通常情况下，每个粒子只做一次，但是由于粒子生成器不知道属性，我们希望在构造时不做。这个函数的另一个用途是在粒子转移到一个新进程之后。
     *
     */
    void
    set_property_pool(PropertyPool<dim, spacedim> &property_pool);

    /**
     * 返回这个粒子是否有一个有效的属性库和一个有效的属性句柄。
     *
     */
    bool
    has_properties() const;

    /**
     * 设置此粒子的属性。          @param  [in] new_properties
     * 一个包含此粒子的新属性的向量。
     *
     */
    void
    set_properties(const std::vector<double> &new_properties);

    /**
     * 设置这个粒子的属性。          @param  [in] new_properties
     * 一个ArrayView，指向包含此粒子的新属性的内存位置。
     *
     */
    void
    set_properties(const ArrayView<const double> &new_properties);

    /**
     * 获得对这个粒子的属性的写访问权。          @return
     * 这个粒子的属性的ArrayView。
     *
     */
    const ArrayView<double>
    get_properties();

    /**
     * 获取对该粒子的属性的读取权限。          @return
     * 此粒子的属性的ArrayView。
     *
     */
    const ArrayView<const double>
    get_properties() const;

    /**
     * 返回这个粒子在所有数据被序列化的情况下所占的字节数（即这个类的write_data函数所写入的字节数）。
     *
     */
    std::size_t
    serialized_size_in_bytes() const;

    /**
     * 获取当前粒子周围的单元格的迭代器。
     * 由于粒子是以三角形的结构组织的，但是三角形本身并不存储在粒子中，所以这个操作需要对三角形的引用。
     *
     */
    typename Triangulation<dim, spacedim>::cell_iterator
    get_surrounding_cell(
      const Triangulation<dim, spacedim> &triangulation) const;

    /**
     * 使用[BOOST序列化库](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html)对该类的内容进行序列化。
     *
     */
    template <class Archive>
    void
    serialize(Archive &ar, const unsigned int version);

    /**
     * 将ParticleAccessor推进到下一个粒子。
     *
     */
    void
    next();

    /**
     * 将ParticleAccessor移动到前一个粒子。
     *
     */
    void
    prev();

    /**
     * 不等式运算符。
     *
     */
    bool
    operator!=(const ParticleAccessor<dim, spacedim> &other) const;

    /**
     * 等价运算符。
     *
     */
    bool
    operator==(const ParticleAccessor<dim, spacedim> &other) const;

  private:
    /**
     * 构建一个无效的访问器。这样的对象是不能使用的。
     *
     */
    ParticleAccessor();

    /**
     * 从对地图的引用和对地图的迭代器构造一个访问器。这个构造函数是`private`的，所以它只能被朋友类访问。
     *
     */
    ParticleAccessor(
      const std::multimap<internal::LevelInd, Particle<dim, spacedim>> &map,
      const typename std::multimap<internal::LevelInd,
                                   Particle<dim, spacedim>>::iterator
        &particle);

  private:
    /**
     * 一个指向存储粒子的容器的指针。很明显，如果容器发生变化，这个访问器就会失效。
     *
     */
    std::multimap<internal::LevelInd, Particle<dim, spacedim>> *map;

    /**
     * 一个进入粒子容器的迭代器。很明显，如果容器发生变化，这个访问器就会失效。
     *
     */
    typename std::multimap<internal::LevelInd,
                           Particle<dim, spacedim>>::iterator particle;

    // Make ParticleIterator a friend to allow it constructing
    // ParticleAccessors.
    template <int, int>
    friend class ParticleIterator;
    template <int, int>
    friend class ParticleHandler;
  };



  template <int dim, int spacedim>
  template <class Archive>
  void
  ParticleAccessor<dim, spacedim>::serialize(Archive &          ar,
                                             const unsigned int version)
  {
    return particle->second.serialize(ar, version);
  }


  // ------------------------- inline functions ------------------------------

  template <int dim, int spacedim>
  inline ParticleAccessor<dim, spacedim>::ParticleAccessor()
    : map(nullptr)
    , particle()
  {}



  template <int dim, int spacedim>
  inline ParticleAccessor<dim, spacedim>::ParticleAccessor(
    const std::multimap<internal::LevelInd, Particle<dim, spacedim>> &map,
    const typename std::multimap<internal::LevelInd,
                                 Particle<dim, spacedim>>::iterator & particle)
    : map(const_cast<
          std::multimap<internal::LevelInd, Particle<dim, spacedim>> *>(&map))
    , particle(particle)
  {}



  template <int dim, int spacedim>
  inline const void *
  ParticleAccessor<dim, spacedim>::read_particle_data_from_memory(
    const void *data)
  {
    Assert(particle != map->end(), ExcInternalError());

    return particle->second.read_particle_data_from_memory(data);
  }



  template <int dim, int spacedim>
  inline void *
  ParticleAccessor<dim, spacedim>::write_particle_data_to_memory(
    void *data) const
  {
    Assert(particle != map->end(), ExcInternalError());

    return particle->second.write_particle_data_to_memory(data);
  }



  template <int dim, int spacedim>
  inline void
  ParticleAccessor<dim, spacedim>::set_location(const Point<spacedim> &new_loc)
  {
    Assert(particle != map->end(), ExcInternalError());

    particle->second.set_location(new_loc);
  }



  template <int dim, int spacedim>
  inline const Point<spacedim> &
  ParticleAccessor<dim, spacedim>::get_location() const
  {
    Assert(particle != map->end(), ExcInternalError());

    return particle->second.get_location();
  }



  template <int dim, int spacedim>
  inline void
  ParticleAccessor<dim, spacedim>::set_reference_location(
    const Point<dim> &new_loc)
  {
    Assert(particle != map->end(), ExcInternalError());

    particle->second.set_reference_location(new_loc);
  }



  template <int dim, int spacedim>
  inline const Point<dim> &
  ParticleAccessor<dim, spacedim>::get_reference_location() const
  {
    Assert(particle != map->end(), ExcInternalError());

    return particle->second.get_reference_location();
  }



  template <int dim, int spacedim>
  inline types::particle_index
  ParticleAccessor<dim, spacedim>::get_id() const
  {
    Assert(particle != map->end(), ExcInternalError());

    return particle->second.get_id();
  }



  template <int dim, int spacedim>
  inline void
  ParticleAccessor<dim, spacedim>::set_property_pool(
    PropertyPool<dim, spacedim> &new_property_pool)
  {
    Assert(particle != map->end(), ExcInternalError());

    particle->second.set_property_pool(new_property_pool);
  }



  template <int dim, int spacedim>
  inline bool
  ParticleAccessor<dim, spacedim>::has_properties() const
  {
    Assert(particle != map->end(), ExcInternalError());

    return particle->second.has_properties();
  }



  template <int dim, int spacedim>
  inline void
  ParticleAccessor<dim, spacedim>::set_properties(
    const std::vector<double> &new_properties)
  {
    Assert(particle != map->end(), ExcInternalError());

    particle->second.set_properties(new_properties);
  }



  template <int dim, int spacedim>
  inline void
  ParticleAccessor<dim, spacedim>::set_properties(
    const ArrayView<const double> &new_properties)
  {
    Assert(particle != map->end(), ExcInternalError());

    particle->second.set_properties(new_properties);
  }



  template <int dim, int spacedim>
  inline const ArrayView<const double>
  ParticleAccessor<dim, spacedim>::get_properties() const
  {
    Assert(particle != map->end(), ExcInternalError());

    return particle->second.get_properties();
  }



  template <int dim, int spacedim>
  inline typename Triangulation<dim, spacedim>::cell_iterator
  ParticleAccessor<dim, spacedim>::get_surrounding_cell(
    const Triangulation<dim, spacedim> &triangulation) const
  {
    Assert(particle != map->end(), ExcInternalError());

    const typename Triangulation<dim, spacedim>::cell_iterator cell(
      &triangulation, particle->first.first, particle->first.second);
    return cell;
  }



  template <int dim, int spacedim>
  inline const ArrayView<double>
  ParticleAccessor<dim, spacedim>::get_properties()
  {
    Assert(particle != map->end(), ExcInternalError());

    return particle->second.get_properties();
  }



  template <int dim, int spacedim>
  inline std::size_t
  ParticleAccessor<dim, spacedim>::serialized_size_in_bytes() const
  {
    Assert(particle != map->end(), ExcInternalError());

    return particle->second.serialized_size_in_bytes();
  }



  template <int dim, int spacedim>
  inline void
  ParticleAccessor<dim, spacedim>::next()
  {
    Assert(particle != map->end(), ExcInternalError());
    ++particle;
  }



  template <int dim, int spacedim>
  inline void
  ParticleAccessor<dim, spacedim>::prev()
  {
    Assert(particle != map->begin(), ExcInternalError());
    --particle;
  }



  template <int dim, int spacedim>
  inline bool
  ParticleAccessor<dim, spacedim>::
  operator!=(const ParticleAccessor<dim, spacedim> &other) const
  {
    return (map != other.map) || (particle != other.particle);
  }



  template <int dim, int spacedim>
  inline bool
  ParticleAccessor<dim, spacedim>::
  operator==(const ParticleAccessor<dim, spacedim> &other) const
  {
    return (map == other.map) && (particle == other.particle);
  }


} // namespace Particles

DEAL_II_NAMESPACE_CLOSE

namespace boost
{
  namespace geometry
  {
    namespace index
    {
      // Forward declaration of bgi::indexable
      template <class T>
      struct indexable;

      /**
       * 确保我们可以从 Particles::ParticleAccessor
       * 对象中构建一个RTree。
       *
       */
      template <int dim, int spacedim>
      struct indexable<dealii::Particles::ParticleAccessor<dim, spacedim>>
      {
        /**
         * boost::rtree
         * 期望一个可索引对象的常量引用。对于一个
         * Particles::Particle 对象，这是它的引用位置。
         *
         */
        using result_type = const dealii::Point<spacedim> &;

        result_type
        operator()(const dealii::Particles::ParticleAccessor<dim, spacedim>
                     &accessor) const
        {
          return accessor.get_location();
        }
      };
    } // namespace index
  }   // namespace geometry
} // namespace boost

#endif


