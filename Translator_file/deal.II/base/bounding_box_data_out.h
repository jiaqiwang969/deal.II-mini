//include/deal.II-translator/base/bounding_box_data_out_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2020 by the deal.II authors
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

#ifndef dealii_bounding_box_data_out_h
#define dealii_bounding_box_data_out_h

#include <deal.II/base/config.h>

#include <deal.II/base/data_out_base.h>
#include <deal.II/base/geometry_info.h>
#include <deal.II/base/iterator_range.h>
#include <deal.II/base/point.h>

#include <deal.II/boost_adaptors/bounding_box.h>
#include <deal.II/boost_adaptors/point.h>
#include <deal.II/boost_adaptors/segment.h>

DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
#include <boost/geometry/index/rtree.hpp>
#include <boost/geometry/strategies/strategies.hpp>
DEAL_II_ENABLE_EXTRA_DIAGNOSTICS


DEAL_II_NAMESPACE_OPEN

/**
 * 这个类为BoundingBox对象生成图形输出，从任何可以通过提升转换为BoundingBox的对象开始。
 *
 *
 */
template <int dim>
class BoundingBoxDataOut : public DataOutInterface<dim, dim>
{
public:
  /**
   * 构造函数。
   *
   */
  BoundingBoxDataOut() = default;

  /**
   * 解构器。
   *
   */
  ~BoundingBoxDataOut() = default;

  /**
   * 从一系列对象中生成补丁，可以通过boost转换为BoundingBox对象的集合。
   * 你可以把BoundingBox对象的迭代器传给这个函数，或者传给第一元素是BoundingBox对象的对或图元。
   *
   */
  template <class ConvertibleToBoundingBoxIterator>
  void
  build_patches(const ConvertibleToBoundingBoxIterator &begin,
                const ConvertibleToBoundingBoxIterator &end);

  /**
   * 从一个对象的容器中生成补丁，可以通过boost转换为BoundingBox对象的集合。
   * 你可以向这个函数传递一个BoundingBox对象的容器，或者一个第一元素是BoundingBox对象的对或图元的容器。
   *
   */
  template <class Container>
  void
  build_patches(const Container &boxes);

  /**
   * 为 build_patches() 生成的每个输出补丁附加数据。     @p
   * datasets
   * 参数应与你在调用build_patches()时使用的容器大小相同，而每个条目应与
   * @p dataset_names 参数的大小相同。      @param[in]  datasets
   * 要附加到每个补丁的实际数据  @param[in]  dataset_names
   * 数据集的每个组件的名称
   *
   */
  void
  add_datasets(const std::vector<std::vector<double>> &datasets,
               const std::vector<std::string> &        dataset_names);

protected:
  // Copy doc
  virtual const std::vector<dealii::DataOutBase::Patch<dim, dim>> &
  get_patches() const override;

  // Copy doc
  virtual std::vector<std::string>
  get_dataset_names() const override;

private:
  /**
   * 实际的方框。
   *
   */
  std::vector<DataOutBase::Patch<dim, dim>> patches;

  /**
   * 数据集的名称。
   *
   */
  std::vector<std::string> dataset_names;
};


// Template and inline functions
#ifndef DOXYGEN
template <int dim>
template <class ConvertibleToBoundingBoxIterator>
void
BoundingBoxDataOut<dim>::build_patches(
  const ConvertibleToBoundingBoxIterator &begin,
  const ConvertibleToBoundingBoxIterator &end)
{
  using Getter = boost::geometry::index::indexable<
    typename ConvertibleToBoundingBoxIterator::value_type>;
  Getter                 getter;
  constexpr unsigned int boxdim =
    boost::geometry::dimension<typename Getter::result_type>::value;
  const unsigned int N = std::distance(begin, end);
  static_assert(boxdim == dim, "Bounding boxes are of the wrong dimension!");

  dataset_names.clear();
  patches.resize(N);

  unsigned int i = 0;
  for (const auto &value :
       IteratorRange<ConvertibleToBoundingBoxIterator>(begin, end))
    {
      BoundingBox<dim> box;
      boost::geometry::convert(getter(*value), box);
      for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_cell; ++v)
        {
          patches[i].vertices[v]          = box.vertex(v);
          patches[i].patch_index          = i;
          patches[i].n_subdivisions       = 1;
          patches[i].points_are_available = false;
        }
      ++i;
    }
}



template <int dim>
template <class Container>
void
BoundingBoxDataOut<dim>::build_patches(const Container &boxes)
{
  build_patches(boxes.begin(), boxes.end());
}



template <int dim>
void
BoundingBoxDataOut<dim>::add_datasets(
  const std::vector<std::vector<double>> &datasets,
  const std::vector<std::string> &        names)
{
  AssertDimension(datasets.size(), patches.size());
  dataset_names = names;
  for (unsigned int i = 0; i < datasets.size(); ++i)
    {
      AssertDimension(datasets[i].size(), names.size());
      patches[i].data.reinit(names.size(),
                             GeometryInfo<dim>::vertices_per_cell);
      for (unsigned int j = 0; j < names.size(); ++j)
        for (unsigned int k = 0; k < GeometryInfo<dim>::vertices_per_cell; ++k)
          patches[i].data(j, k) = datasets[i][j];
    }
}



template <int dim>
const std::vector<DataOutBase::Patch<dim, dim>> &
BoundingBoxDataOut<dim>::get_patches() const
{
  return patches;
}



template <int dim>
std::vector<std::string>
BoundingBoxDataOut<dim>::get_dataset_names() const
{
  return dataset_names;
}

#endif

DEAL_II_NAMESPACE_CLOSE

#endif

