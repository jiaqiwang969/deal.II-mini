//include/deal.II-translator/particles/property_pool_0.txt
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

#ifndef dealii_particles_property_pool_h
#define dealii_particles_property_pool_h

#include <deal.II/base/config.h>

#include <deal.II/base/array_view.h>
#include <deal.II/base/point.h>


DEAL_II_NAMESPACE_OPEN

namespace types
{
   /* Type definitions */ 

#ifdef DEAL_II_WITH_64BIT_INDICES
  /**
   * 用于粒子的索引的类型。虽然在顺序计算中，40亿个32位无符号整数的索引已经足够了，但使用数百个进程的并行计算会溢出这个数字，我们需要一个更大的索引空间。我们在这里利用了控制自由度指数的同一个构建变量，因为自由度的数量和粒子的数量通常在同一个数量级上。
   * 数据类型总是表示无符号整数类型。
   *
   */
  using particle_index = uint64_t;

#  ifdef DEAL_II_WITH_MPI
  /**
   * 一个标识符，表示与 types::global_dof_index.
   * 相关的MPI类型。
   *
   */
#    define DEAL_II_PARTICLE_INDEX_MPI_TYPE MPI_UINT64_T
#  endif

#else
  /**
   * 用于粒子的索引的类型。虽然在顺序计算中，40亿个32位无符号整数的索引已经足够了，但使用数百个进程的并行计算会溢出这个数字，我们需要一个更大的索引空间。我们在这里利用了控制自由度指数的同一个构建变量，因为自由度的数量和粒子的数量通常在同一个数量级上。
   * 数据类型总是表示无符号整数类型。
   *
   */
  using particle_index = unsigned int;

#  ifdef DEAL_II_WITH_MPI
  /**
   * 一个标识符，表示与 types::global_dof_index.
   * 相关的MPI类型。
   *
   */
#    define DEAL_II_PARTICLE_INDEX_MPI_TYPE MPI_UNSIGNED
#  endif
#endif
} // namespace types

namespace Particles
{
  /**
   * 这个类管理着一个内存空间，与ParticleHandler相关的所有粒子都在其中存储它们的属性。它还存储粒子的位置和参考位置。
   * 这个类的基本原理是，因为通常每个粒子都存储相同数量的属性，并且因为算法通常会遍历所有的粒子，对所有粒子的属性进行相同的操作，所以让用于属性的内存由一个中央管理器来处理会更有效率。
   * 这样，粒子就不会存储一个指向它们存储属性的内存区域的指针，而是一个
   * "句柄"，然后PropertyPool类将其转化为指向具体内存的指针。
   * 综上所述，目前的实现只提供了这种接口，但仍然对粒子所要求的每一组属性使用简单的新建/删除分配。此外，目前的实现假设每个粒子的属性数量相同，但当然PropertyType可以包含一个指向动态分配的内存的指针，每个粒子的大小不同（这个内存不会被这个类所管理）。
   *
   */
  template <int dim, int spacedim = dim>
  class PropertyPool
  {
  public:
    /**
     * 返回给粒子的句柄的类型定义，该句柄唯一地标识了为该粒子保留的内存槽。
     *
     */
    using Handle = unsigned int;

    /**
     * 为句柄定义一个默认（无效）值。
     *
     */
    static const Handle invalid_handle;

    /**
     * 构造器。存储每个保留槽的属性数量。
     *
     */
    PropertyPool(const unsigned int n_properties_per_slot);

    /**
     * 解构器。这个函数确保之前使用allocate_properties_array()分配的所有内存也已经通过deallocate_properties_array()返回。
     *
     */
    ~PropertyPool();

    /**
     * 清除这个类所分配的动态内存。这个函数确保之前使用allocate_properties_array()分配的所有内存也已经通过deallocate_properties_array()返回。
     *
     */
    void
    clear();

    /**
     * 返回一个新的句柄，允许粒子存储诸如属性和位置等信息。这也是在这个PropertyPool变量中分配内存。
     *
     */
    Handle
    register_particle();

    /**
     * 返回一个通过register_particle()获得的句柄，并将分配给存储粒子数据的内存标记为空闲，以便重新使用。
     *
     */
    void
    deregister_particle(Handle &handle);

    /**
     * 返回由给定的`handle'标识的粒子的位置。
     *
     */
    const Point<spacedim> &
    get_location(const Handle handle) const;

    /**
     * 设置由给定的`handle'标识的粒子的位置。
     *
     */
    void
    set_location(const Handle handle, const Point<spacedim> &new_location);

    /**
     * 返回由给定的 "handle "标识的粒子的参考位置。
     *
     */
    const Point<dim> &
    get_reference_location(const Handle handle) const;

    /**
     * 设置由给定的 "handle "标识的粒子的参考位置。
     *
     */
    void
    set_reference_location(const Handle      handle,
                           const Point<dim> &new_reference_location);

    /**
     * 返回由给定的`handle'标识的这个粒子的ID号。
     *
     */
    types::particle_index
    get_id(const Handle handle) const;

    /**
     * 设置由给定的 "handle "标识的这个粒子的ID号。
     *
     */
    void
    set_id(const Handle handle, const types::particle_index &new_id);

    /**
     * 返回一个ArrayView给定句柄所对应的属性  @p handle.
     *
     */
    ArrayView<double>
    get_properties(const Handle handle);

    /**
     * 预留储存 @p size 粒子的属性所需的动态内存。
     *
     */
    void
    reserve(const std::size_t size);

    /**
     * 返回池中每个槽位存储的属性数量。
     *
     */
    unsigned int
    n_properties_per_slot() const;

  private:
    /**
     * 每个粒子保留的属性数量。
     *
     */
    const unsigned int n_properties;

    /**
     * 一个存储粒子位置的向量。它的索引方式与`reference_locations`和`properties`数组相同，也就是通过句柄。
     *
     */
    std::vector<Point<spacedim>> locations;

    /**
     * 一个存储粒子参考位置的向量。它的索引方式与`locations`和`properties`数组相同，即通过句柄。
     *
     */
    std::vector<Point<dim>> reference_locations;

    /**
     * 一个存储粒子的唯一标识符的向量。它的索引方式与`位置`和`属性`数组的索引方式相同，即通过句柄。
     *
     */
    std::vector<types::particle_index> ids;

    /**
     * 当前分配的属性（无论是分配给粒子的还是可分配的）。它的索引方式与`locations`和`reference_locations`数组相同，通过句柄。
     *
     */
    std::vector<double> properties;

    /**
     * 由allocate_properties_array()创建并由deallocate_properties_array()销毁的句柄的集合。由于内存仍然被分配，这些句柄可以被重新用于新的粒子以避免内存分配。
     *
     */
    std::vector<Handle> currently_available_handles;
  };



   /* ---------------------- inline and template functions ------------------ */ 

  template <int dim, int spacedim>
  inline const Point<spacedim> &
  PropertyPool<dim, spacedim>::get_location(const Handle handle) const
  {
    const std::vector<double>::size_type data_index =
      (handle != invalid_handle) ? handle : 0;

    // Ideally we would need to assert that 'handle' has not been deallocated
    // by searching through 'currently_available_handles'. However, this
    // is expensive and this function is performance critical, so instead
    // just check against the array range, and rely on the fact
    // that handles are invalidated when handed over to
    // deallocate_properties_array().
    Assert(data_index <= locations.size() - 1,
           ExcMessage("Invalid location handle. This can happen if the "
                      "handle was duplicated and then one copy was deallocated "
                      "before trying to access the properties."));

    return locations[data_index];
  }



  template <int dim, int spacedim>
  inline void
  PropertyPool<dim, spacedim>::set_location(const Handle           handle,
                                            const Point<spacedim> &new_location)
  {
    const std::vector<double>::size_type data_index =
      (handle != invalid_handle) ? handle : 0;

    // Ideally we would need to assert that 'handle' has not been deallocated
    // by searching through 'currently_available_handles'. However, this
    // is expensive and this function is performance critical, so instead
    // just check against the array range, and rely on the fact
    // that handles are invalidated when handed over to
    // deallocate_properties_array().
    Assert(data_index <= locations.size() - 1,
           ExcMessage("Invalid location handle. This can happen if the "
                      "handle was duplicated and then one copy was deallocated "
                      "before trying to access the properties."));

    locations[data_index] = new_location;
  }



  template <int dim, int spacedim>
  inline const Point<dim> &
  PropertyPool<dim, spacedim>::get_reference_location(const Handle handle) const
  {
    const std::vector<double>::size_type data_index =
      (handle != invalid_handle) ? handle : 0;

    // Ideally we would need to assert that 'handle' has not been deallocated
    // by searching through 'currently_available_handles'. However, this
    // is expensive and this function is performance critical, so instead
    // just check against the array range, and rely on the fact
    // that handles are invalidated when handed over to
    // deallocate_properties_array().
    Assert(data_index <= reference_locations.size() - 1,
           ExcMessage("Invalid location handle. This can happen if the "
                      "handle was duplicated and then one copy was deallocated "
                      "before trying to access the properties."));

    return reference_locations[data_index];
  }



  template <int dim, int spacedim>
  inline void
  PropertyPool<dim, spacedim>::set_reference_location(
    const Handle      handle,
    const Point<dim> &new_reference_location)
  {
    const std::vector<double>::size_type data_index =
      (handle != invalid_handle) ? handle : 0;

    // Ideally we would need to assert that 'handle' has not been deallocated
    // by searching through 'currently_available_handles'. However, this
    // is expensive and this function is performance critical, so instead
    // just check against the array range, and rely on the fact
    // that handles are invalidated when handed over to
    // deallocate_properties_array().
    Assert(data_index <= reference_locations.size() - 1,
           ExcMessage("Invalid location handle. This can happen if the "
                      "handle was duplicated and then one copy was deallocated "
                      "before trying to access the properties."));

    reference_locations[data_index] = new_reference_location;
  }



  template <int dim, int spacedim>
  inline types::particle_index
  PropertyPool<dim, spacedim>::get_id(const Handle handle) const
  {
    const std::vector<double>::size_type data_index =
      (handle != invalid_handle) ? handle : 0;

    // Ideally we would need to assert that 'handle' has not been deallocated
    // by searching through 'currently_available_handles'. However, this
    // is expensive and this function is performance critical, so instead
    // just check against the array range, and rely on the fact
    // that handles are invalidated when handed over to
    // deallocate_properties_array().
    Assert(data_index <= ids.size() - 1,
           ExcMessage("Invalid location handle. This can happen if the "
                      "handle was duplicated and then one copy was deallocated "
                      "before trying to access the properties."));

    return ids[data_index];
  }



  template <int dim, int spacedim>
  inline void
  PropertyPool<dim, spacedim>::set_id(const Handle                 handle,
                                      const types::particle_index &new_id)
  {
    const std::vector<double>::size_type data_index =
      (handle != invalid_handle) ? handle : 0;

    // Ideally we would need to assert that 'handle' has not been deallocated
    // by searching through 'currently_available_handles'. However, this
    // is expensive and this function is performance critical, so instead
    // just check against the array range, and rely on the fact
    // that handles are invalidated when handed over to
    // deallocate_properties_array().
    Assert(data_index <= ids.size() - 1,
           ExcMessage("Invalid location handle. This can happen if the "
                      "handle was duplicated and then one copy was deallocated "
                      "before trying to access the properties."));

    ids[data_index] = new_id;
  }



  template <int dim, int spacedim>
  inline ArrayView<double>
  PropertyPool<dim, spacedim>::get_properties(const Handle handle)
  {
    const std::vector<double>::size_type data_index =
      (handle != invalid_handle) ? handle * n_properties : 0;

    // Ideally we would need to assert that 'handle' has not been deallocated
    // by searching through 'currently_available_handles'. However, this
    // is expensive and this function is performance critical, so instead
    // just check against the array range, and rely on the fact
    // that handles are invalidated when handed over to
    // deallocate_properties_array().
    Assert(data_index <= properties.size() - n_properties,
           ExcMessage("Invalid property handle. This can happen if the "
                      "handle was duplicated and then one copy was deallocated "
                      "before trying to access the properties."));

    return ArrayView<double>(properties.data() + data_index, n_properties);
  }


} // namespace Particles

DEAL_II_NAMESPACE_CLOSE

#endif


