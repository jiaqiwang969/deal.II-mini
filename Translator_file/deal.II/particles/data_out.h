//include/deal.II-translator/particles/data_out_0.txt
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
#ifndef dealii_particles_data_out_h
#define dealii_particles_data_out_h

#include <deal.II/base/config.h>

#include <deal.II/base/data_out_base.h>

#include <deal.II/numerics/data_component_interpretation.h>

#include <string>
#include <vector>

DEAL_II_NAMESPACE_OPEN

namespace Particles
{
  template <int dim, int spacedim>
  class ParticleHandler;

  /**
   * 该类为ParticleHandler对象所存储的粒子生成图形输出。从一个粒子处理程序中，它生成的补丁可以用来编写传统的输出文件。这个类目前只支持写出粒子的位置和它们的ID，不允许写出粒子的附加属性。
   * @ingroup Particle
   *
   */
  template <int dim, int spacedim = dim>
  class DataOut : public dealii::DataOutInterface<0, spacedim>
  {
  public:
    /**
     * Particles::DataOut 类的默认构造函数。
     *
     */
    DataOut() = default;

    /**
     * Particles::DataOut 类的解构器。
     *
     */
    ~DataOut() = default;


    /**
     * 为一个给定的粒子处理程序建立补丁。          @param
     * [in] particles 一个粒子处理程序，将为其建立补丁。
     * 为每个粒子建立一个dim=0的补丁。粒子的位置被用来建立节点位置，粒子的ID被添加为一个数据元素。
     * @param  [in] data_component_names
     * 一个可选的字符串向量，描述每个粒子的属性。只有提供了这个向量，粒子的属性才会被写入。
     * @param  [in] data_component_interpretations
     * 一个可选的向量，控制粒子属性是否被解释为标量、向量或张量。必须与
     * @p data_component_names. 的长度相同。
     *
     */
    void
    build_patches(const Particles::ParticleHandler<dim, spacedim> &particles,
                  const std::vector<std::string> &data_component_names = {},
                  const std::vector<
                    DataComponentInterpretation::DataComponentInterpretation>
                    &data_component_interpretations = {});

  protected:
    /**
     * 返回由data_out类构建的补丁，该补丁是之前使用粒子处理程序构建的。
     *
     */
    virtual const std::vector<DataOutBase::Patch<0, spacedim>> &
    get_patches() const override;

    /**
     * 虚拟函数，通过该函数从该类中获得数据集的名称
     *
     */
    virtual std::vector<std::string>
    get_dataset_names() const override;


    /**
     * 各自 DataOutInterface::get_nonscalar_data_ranges()
     * 函数的重载。请看那里有更多的文档。
     * 这个函数是对函数 DataOut_DoFData::get_nonscalar_data_ranges().
     * 的重新实现。
     *
     */
    virtual std::vector<
      std::tuple<unsigned int,
                 unsigned int,
                 std::string,
                 DataComponentInterpretation::DataComponentInterpretation>>
    get_nonscalar_data_ranges() const override;

  private:
    /**
     * 这是一个补丁的向量，在每次调用 build_patches()
     * 时都会创建。这些补丁在基类的输出例程中使用。
     *
     */
    std::vector<DataOutBase::Patch<0, spacedim>> patches;

    /**
     * 储存在补丁中的所有数据组件的字段名的向量。
     *
     */
    std::vector<std::string> dataset_names;

    /**
     * 一个向量，对于当前数据集的每个数据成分，表明它们是标量字段、向量字段的一部分，还是其他任何支持的数据种类。
     *
     */
    std::vector<DataComponentInterpretation::DataComponentInterpretation>
      data_component_interpretations;
  };

} // namespace Particles

DEAL_II_NAMESPACE_CLOSE

#endif


