//include/deal.II-translator/particles/particle_0.txt
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

#ifndef dealii_particles_particle_h
#define dealii_particles_particle_h

#include <deal.II/base/config.h>

#include <deal.II/base/array_view.h>
#include <deal.II/base/point.h>
#include <deal.II/base/types.h>

#include <deal.II/particles/property_pool.h>

#include <cstdint>

DEAL_II_NAMESPACE_OPEN

/**
 * 一个命名空间，包含所有与粒子实现相关的类，特别是基本的粒子类。
 *
 *
 */
namespace Particles
{
  namespace internal
  {
    /**
     * 电池级别/索引对的内部别名。
     *
     */
    using LevelInd = std::pair<int, int>;
  } // namespace internal

  /**
   * 一个代表领域中的粒子的类，该领域由某种三角形网格组成。这个类所存储的数据是粒子在整体空间中的位置，粒子在它当前所在的单元的参考坐标系中的位置，一个在所有粒子中唯一的ID号，以及一个数量可变的
   * "属性"。    附在该类每个对象上的 "属性
   * "由一个PropertyPool对象存储。这些属性被存储为一个 "双
   * "变量数组，可以通过ArrayView对象访问。例如，如果人们想给每个粒子配备一个
   * "温度 "和 "化学成分
   * "属性，这些属性会随着粒子的移动而移动（例如，可能会根据一些微分方程在不同的时间步长上发生变化），那么人们会在PropertyPool对象中给每个粒子分配两个属性。
   * 然而，在实践中，人们经常希望将属性与粒子联系起来，而不是像上面的情况那样只是独立的数字。一个例子是，如果人们想要跟踪一个粒子所受到的应力或应变
   *
   * --一个张量值的量。在这些情况下，人们会将<i>interpret</i>这些标量属性作为<i>components
   * of the stress or
   * strain</i>。换句话说，我们首先会告诉PropertyPool为每个粒子分配尽可能多的属性，就像我们想要追踪的张量中的成分一样，然后编写小的转换函数，将get_properties()函数返回的标量属性的ArrayView转换为适当类型的张量。然后，这可以在每个时间步骤中被评估和演化。第二个转换函数将从张量转换回ArrayView对象，通过set_properties()函数将更新的数据存储回粒子中。
   * 当然，在有些情况下，人们关心的属性不是实数（或者在计算机中是浮点），而是分类的。例如，人们可能想把一些粒子标记为
   * "红色"、"蓝色 "或
   * "绿色"。然后，该属性可能被表示为一个整数，或作为
   * "enum
   * "的一个元素。在这些情况下，人们需要想出一种方法来<i>represent</i>这些分类字段的浮点数字。例如，我们可以将
   * "红色 "映射到浮点数1.0，"蓝色 "映射到2.0，"绿色
   * "映射到3.0。那么，在这两种表示法之间进行转换的转换函数应该也不是很难写。
   * @ingroup Particle
   *
   */
  template <int dim, int spacedim = dim>
  class Particle
  {
  public:
    /**
     * 粒子的空构造函数，在原点创建一个粒子。
     *
     */
    Particle();

    /**
     * 粒子的构造函数。这个函数在指定的位置创建一个具有指定ID的粒子。注意，没有对重复的粒子ID进行检查，所以用户必须确保在所有进程中ID是唯一的。数据被存储在一个全局的PropertyPool对象中（对应于全局的
     * "堆"），但以后可以通过调用set_property_pool()转移到另一个属性池。
     * @param[in]  location 粒子的初始位置。      @param[in]
     * reference_location
     * 粒子在参考单元的坐标系中的初始位置。      @param[in]
     * id 粒子的全球唯一ID号。
     *
     */
    Particle(const Point<spacedim> &     location,
             const Point<dim> &          reference_location,
             const types::particle_index id);

    /**
     * 粒子的复制构造器。这个函数创建一个与输入参数的状态完全相同的粒子。复制的数据存储在一个全局的PropertyPool对象中（对应于全局的
     * "堆"），但以后可以通过调用set_property_pool()转移到另一个属性池。
     *
     */
    Particle(const Particle<dim, spacedim> &particle);

    /**
     * 粒子的构造函数。这个函数从一个数据矢量创建一个粒子。数据被存储在一个全局PropertyPool对象中（对应于全局
     * "堆"），但以后可以通过调用set_property_pool()转移到另一个属性池中。这个构造函数通常在通过调用write_data()函数对粒子进行序列化之后被调用。
     * @param[in,out]  begin_data
     * 一个指向内存位置的指针，可以从中读取完全描述一个粒子的信息。然后这个类从这个内存位置反序列化它的数据，并将指针推进到已经读过的数据之外，以初始化粒子信息。
     * @param[in,out]  property_pool
     * 一个可选的指向属性池的指针，用于管理这个粒子所使用的属性数据。如果没有提供这个参数，那么就会使用一个全局属性池；另一方面，如果提供了一个非空的指针，那么这个构造函数就会假设
     * @p begin_data 包含了由 @p property_pool.
     * 分配的相同长度和类型的序列化数据。
     * 如果这里提供的数据指针对应于有属性的粒子的数据，那么这个函数只有在提供一个属性池作为第二个参数，能够存储每个粒子正确数量的属性时才会成功。
     *
     */
    Particle(const void *&                      begin_data,
             PropertyPool<dim, spacedim> *const property_pool = nullptr);

    /**
     * 粒子的移动构造函数，通过窃取现有粒子的状态来创建一个粒子。
     *
     */
    Particle(Particle<dim, spacedim> &&particle) noexcept;

    /**
     * 复制赋值运算符。
     *
     */
    Particle<dim, spacedim> &
    operator=(const Particle<dim, spacedim> &particle);

    /**
     * 移动赋值运算符。
     *
     */
    Particle<dim, spacedim> &
    operator=(Particle<dim, spacedim> &&particle) noexcept;

    /**
     * 销毁器。如果属性句柄是有效的，则将其释放，从而为其他粒子释放内存空间。注意：内存是由属性池管理的，属性池负责对内存的处理。
     *
     */
    ~Particle();

    /**
     * 将粒子数据写入一个数据阵列。这个数组应该足够大，可以容纳这些数据，无效指针应该指向数组的第一个条目，这些数据应该被写入该数组。这个函数是用来序列化所有的粒子属性，然后通过调用适当的构造函数来反序列化这些属性
     * Particle(void&data, PropertyPoolproperty_pool = nullptr);  @param  [in]
     * data 要写进粒子数据的内存位置。          @return
     * 指向已写入数据的阵列后的下一个字节的指针。
     *
     */
    void *
    write_particle_data_to_memory(void *data) const;


    /**
     * 通过使用一个数据数组来更新与粒子相关的所有数据：ID、位置、参考位置，如果有的话，还有属性。这个数组应该足够大，以容纳这些数据，而且空指针应该指向数组的第一个条目，这些数据应该被写入这个数组中。这个函数是用来反序列化粒子数据的，不需要建立一个新的粒子类。这在ParticleHandler中被用来更新幽灵粒子，而不需要取消分配和重新分配内存。
     * @param[in]  data
     * 一个指向内存位置的指针，可以从中读取完全描述一个粒子的信息。然后这个类从这个内存位置去序列化它的数据。
     * @return
     * 一个指向从该阵列中读取数据后的下一个字节的指针。
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
     * 返回这个粒子的ID号。粒子的ID是一个即使在并行计算中也是全局唯一的属性，如果它从当前处理器拥有的单元移动到另一个处理器拥有的单元，或者它所在的单元的所有权被转移到另一个处理器，那么它的ID将与粒子的其他属性一起被转移。
     *
     */
    types::particle_index
    get_id() const;

    /**
     * 设置这个粒子的ID号。粒子的ID旨在成为一个全局唯一的属性，即使在并行计算中也是如此，如果它从当前处理器拥有的单元移动到另一个处理器拥有的单元，或者它所在的单元的所有权被转移到另一个处理器，那么它的ID将与粒子的其他属性一起转移。因此，在设置粒子的ID时，需要注意确保粒子具有全球唯一的ID。(ParticleHandler本身并不检查如此设置的粒子ID在并行设置中是否是全局唯一的，因为这将是一个非常昂贵的操作。)
     * @param[in]  new_id 这个粒子的新ID号。
     *
     */
    void
    set_id(const types::particle_index &new_id);

    /**
     * 告诉粒子在哪里存储它的属性（即使它并不拥有属性）。通常情况下，每个粒子只做一次，但是由于粒子不知道属性，我们希望在构造时不做。这个函数的另一个用途是在粒子转移到一个新进程之后。
     * 如果一个粒子已经在一个属性池中存储了属性，那么它们的值就会被保存下来，之前的属性池中的内存会被释放，而这个粒子的属性副本将被分配到新的属性池中。
     *
     */
    void
    set_property_pool(PropertyPool<dim, spacedim> &property_pool);

    /**
     * 返回这个粒子是否有一个有效的属性池和一个有效的属性句柄。
     *
     */
    bool
    has_properties() const;

    /**
     * 设置此粒子的属性。          @param  [in] new_properties
     * 一个包含此粒子的新属性的ArrayView。
     *
     */
    void
    set_properties(const ArrayView<const double> &new_properties);

    /**
     * 获得对这个粒子的属性的写入权限。如果粒子还没有属性，但可以访问PropertyPool对象，它将分配属性以允许写入它们。如果它没有属性，也不能访问一个PropertyPool，这个函数将抛出一个异常。
     * @return  这个粒子的属性的ArrayView。
     *
     */
    const ArrayView<double>
    get_properties();

    /**
     * 获取对该粒子的属性的读取权限。如果粒子没有属性，这个函数会抛出一个异常。
     * @return  此粒子的属性的ArrayView。
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
     * 使用[BOOST序列化库](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html)将此对象的数据写到一个流中，以便进行序列化。
     *
     */
    template <class Archive>
    void
    save(Archive &ar, const unsigned int version) const;

    /**
     * 为了序列化的目的，使用[BOOST序列化库](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html)从流中读取此对象的数据。
     * 请注意，为了正确地存储属性，在读取时必须知道这个粒子的属性库，即在调用这个函数之前，必须调用set_property_pool()。
     *
     */
    template <class Archive>
    void
    load(Archive &ar, const unsigned int version);

    /**
     * 释放属性池的内存
     *
     */
    void
    free_properties();

#ifdef DOXYGEN
    /**
     * 使用[BOOST序列化库](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html)从一个流中写入和读取此对象的数据，以便进行序列化。
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
     * 一个全局属性池，当一个粒子没有与属于例如ParticleHandler的属性池相关联时使用。
     *
     */
    static PropertyPool<dim, spacedim> global_property_pool;

    /**
     * 一个指向属性池的指针。从句柄到实际内存位置的转换是必要的。
     *
     */
    PropertyPool<dim, spacedim> *property_pool;

    /**
     * 一个指向所有粒子属性的句柄
     *
     */
    typename PropertyPool<dim, spacedim>::Handle property_pool_handle;
  };



   /* ---------------------- inline and template functions ------------------ */ 

  template <int dim, int spacedim>
  template <class Archive>
  inline void
  Particle<dim, spacedim>::load(Archive &ar, const unsigned int)
  {
    unsigned int n_properties = 0;

    Point<spacedim>       location;
    Point<dim>            reference_location;
    types::particle_index id;
    ar &location &reference_location &id &n_properties;

    set_location(location);
    set_reference_location(reference_location);
    set_id(id);

    if (n_properties > 0)
      {
        ArrayView<double> properties(get_properties());
        Assert(
          properties.size() == n_properties,
          ExcMessage(
            "This particle was serialized with " +
            std::to_string(n_properties) +
            " properties, but the new property handler provides space for " +
            std::to_string(properties.size()) +
            " properties. Deserializing a particle only works for matching property sizes."));

        ar &boost::serialization::make_array(properties.data(), n_properties);
      }
  }



  template <int dim, int spacedim>
  template <class Archive>
  inline void
  Particle<dim, spacedim>::save(Archive &ar, const unsigned int) const
  {
    unsigned int n_properties = 0;
    if ((property_pool != nullptr) &&
        (property_pool_handle != PropertyPool<dim, spacedim>::invalid_handle))
      n_properties = get_properties().size();

    Point<spacedim>       location           = get_location();
    Point<dim>            reference_location = get_reference_location();
    types::particle_index id                 = get_id();

    ar &location &reference_location &id &n_properties;

    if (n_properties > 0)
      ar &boost::serialization::make_array(get_properties().data(),
                                           n_properties);
  }



  template <int dim, int spacedim>
  inline void
  Particle<dim, spacedim>::set_location(const Point<spacedim> &new_loc)
  {
    property_pool->set_location(property_pool_handle, new_loc);
  }



  template <int dim, int spacedim>
  inline const Point<spacedim> &
  Particle<dim, spacedim>::get_location() const
  {
    return property_pool->get_location(property_pool_handle);
  }



  template <int dim, int spacedim>
  inline void
  Particle<dim, spacedim>::set_reference_location(const Point<dim> &new_loc)
  {
    property_pool->set_reference_location(property_pool_handle, new_loc);
  }



  template <int dim, int spacedim>
  inline const Point<dim> &
  Particle<dim, spacedim>::get_reference_location() const
  {
    return property_pool->get_reference_location(property_pool_handle);
  }



  template <int dim, int spacedim>
  inline types::particle_index
  Particle<dim, spacedim>::get_id() const
  {
    return property_pool->get_id(property_pool_handle);
  }



  template <int dim, int spacedim>
  inline void
  Particle<dim, spacedim>::set_id(const types::particle_index &new_id)
  {
    property_pool->set_id(property_pool_handle, new_id);
  }



  template <int dim, int spacedim>
  inline void
  Particle<dim, spacedim>::set_property_pool(
    PropertyPool<dim, spacedim> &new_property_pool)
  {
    // First, we do want to save any properties that may
    // have previously been set, and copy them over to the memory allocated
    // on the new pool.
    //
    // It is possible that a particle currently has no properties -- for
    // example if it has been created without an associated property
    // pool (i.e., uses the default global pool which does not store any
    // properties) but that the new pool has properties. In that case,
    // there is simply nothing to transfer -- but the register_particle()
    // call here will make sure that the newly allocated properties are
    // zero-initialized.
    const typename PropertyPool<dim, spacedim>::Handle new_handle =
      new_property_pool.register_particle();

    const Point<spacedim>       location           = get_location();
    const Point<dim>            reference_location = get_reference_location();
    const types::particle_index id                 = get_id();

    if ( /* old pool */  has_properties())
      {
        ArrayView<const double> old_properties = this->get_properties();
        ArrayView<double>       new_properties =
          new_property_pool.get_properties(new_handle);
        std::copy(old_properties.cbegin(),
                  old_properties.cend(),
                  new_properties.begin());
      }

    // Now release the old memory handle
    property_pool->deregister_particle(property_pool_handle);


    // Then set the pointer to the property pool we want to use. Also set the
    // handle to any properties.
    property_pool        = &new_property_pool;
    property_pool_handle = new_handle;

    // Now also store the saved locations
    set_location(location);
    set_reference_location(reference_location);
    set_id(id);
  }



  template <int dim, int spacedim>
  inline const ArrayView<const double>
  Particle<dim, spacedim>::get_properties() const
  {
    if (has_properties() == false)
      return {};
    else
      return property_pool->get_properties(property_pool_handle);
  }



  template <int dim, int spacedim>
  inline bool
  Particle<dim, spacedim>::has_properties() const
  {
    // Particles always have a property pool associated with them,
    // but we can access properties only if there is a valid handle.
    // The only way a particle can have no valid handle if it has
    // been moved-from -- but that leaves an object in an invalid
    // state, and so we can just assert that that can't be the case.
    Assert((property_pool_handle !=
            PropertyPool<dim, spacedim>::invalid_handle),
           ExcInternalError());
    return (property_pool->n_properties_per_slot() > 0);
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
       * 确保我们可以构建一个 Particles::Particle 对象的RTree。
       *
       */
      template <int dim, int spacedim>
      struct indexable<dealii::Particles::Particle<dim, spacedim>>
      {
        /**
         * boost::rtree
         * 期望一个可索引对象的常量引用。对于一个
         * Particles::Particle 对象，这是它的引用位置。
         *
         */
        using result_type = const dealii::Point<spacedim> &;

        result_type
        operator()(
          const dealii::Particles::Particle<dim, spacedim> &particle) const
        {
          return particle.get_location();
        }
      };

    } // namespace index
  }   // namespace geometry
} // namespace boost

#endif


