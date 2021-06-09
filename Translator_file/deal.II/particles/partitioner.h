//include/deal.II-translator/particles/partitioner_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2020 - 2021 by the deal.II authors
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

#ifndef dealii_particles_partitioner_h
#define dealii_particles_partitioner_h

#include <deal.II/base/config.h>

#include <deal.II/particles/particle_iterator.h>

DEAL_II_NAMESPACE_OPEN

namespace Particles
{
  namespace internal
  {
    /**
     * 缓存结构，用于存储跨处理器交换粒子信息（位置和属性）所需的元素，以便更新鬼魂粒子。
     * 这个结构应该只在希望使用粒子进行工作而不在每次迭代时调用sort_particles_into_subdomain_and_cells时使用。当粒子与粒子之间的相互作用发生在不同的时间尺度上时，这很有用。
     *
     */
    template <int dim, int spacedim>
    struct GhostParticlePartitioner
    {
      /**
       * 一个可以用来迭代域中所有粒子的类型。
       *
       */
      using particle_iterator = ParticleIterator<dim, spacedim>;

      /**
       * 表示是否已经建立了缓存，以防止用无效的缓存更新粒子。
       *
       */
      bool valid = false;

      /**
       * 当前子域的所有可能的邻居的子域ID的矢量。
       *
       */
      std::vector<types::subdomain_id> neighbors;

      /**
       * 大小（neighbors.size()+1）的向量，用于存储必须从当前子域到邻居的数据的起点和终点。对于邻居i，send_pointers[i]表示开始，send_pointers[i+1]表示必须发送的数据的结束。
       *
       */
      std::vector<unsigned int> send_pointers;

      /**
       * 目前居住在本地域的幽灵单元中的粒子集合，由subdomain_id组织。这些粒子相当于分布式向量中的幽灵条目。
       *
       */
      std::map<types::subdomain_id, std::vector<particle_iterator>>
        ghost_particles_by_domain;

      /**
       * 大小为(neighbors.size()+1)的向量，用于存储必须从当前子域上的邻居[i]接收的数据的起点和终点。对于邻居i，recv_pointers[i]表示开始，recv_pointers[i+1]表示必须接收的数据的结束。
       * 这个结构与邻居结合时类似于
       * Utilities::MPI::Partitioner::import_targets 。
       *
       */
      std::vector<unsigned int> recv_pointers;

      /**
       * 幽灵粒子的矢量，其顺序是它们被插入用于存储三角形上的粒子的多图中。该信息用于更新鬼魂粒子信息，而无需清除多图中的鬼魂粒子，从而大大降低了交换鬼魂粒子信息的成本。
       *
       */
      std::vector<typename std::multimap<internal::LevelInd,
                                         Particle<dim, spacedim>>::iterator>
        ghost_particles_iterators;

      /**
       * 临时存储，用于保存要发送给其他处理器的粒子数据，以便在update_ghost_particles()
       * send_recv_particles_properties_and_location()中更新幽灵粒子的信息。
       *
       */
      std::vector<char> send_data;

      /**
       * 保存粒子数据的临时存储器，以便在update_ghost_particles()
       * send_recv_particles_properties_and_location()中接收来自其他处理器的幽灵粒子信息。
       *
       */
      std::vector<char> recv_data;
    };
  } // namespace internal

} // namespace Particles

DEAL_II_NAMESPACE_CLOSE

#endif


