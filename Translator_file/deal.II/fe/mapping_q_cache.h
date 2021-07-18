//include/deal.II-translator/fe/mapping_q_cache_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2019 - 2021 by the deal.II authors
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

#ifndef dealii_mapping_q_cache_h
#define dealii_mapping_q_cache_h


#include <deal.II/base/config.h>

#include <deal.II/base/function.h>
#include <deal.II/base/mg_level_object.h>

#include <deal.II/fe/mapping_q_generic.h>

#include <deal.II/grid/tria.h>


DEAL_II_NAMESPACE_OPEN

// Forward declarations
#ifndef DOXYGEN
template <int, int>
class DoFHandler;
#endif


 /*!@addtogroup mapping */ 
 /*@{*/ 


/**
 * 这个类以 MappingQGeneric::compute_mapping_support_points()
 * 函数的形式为MappingQ家族的对象实现了一个缓存策略，它被用于MappingQGeneric的所有操作中。映射的信息是由
 * MappingQCache::initialize() 函数预先计算出来的。
 * 这个类的使用在  step-65  中有广泛的讨论。
 *
 *
 */
template <int dim, int spacedim = dim>
class MappingQCache : public MappingQGeneric<dim, spacedim>
{
public:
  /**
   * 构造函数。  @p polynomial_degree
   * 表示用于从参考单元映射到实数单元的多项式的度数。
   *
   */
  explicit MappingQCache(const unsigned int polynomial_degree);

  /**
   * 复制构造器。
   *
   */
  explicit MappingQCache(const MappingQCache<dim, spacedim> &mapping);

  /**
   * 解构器。
   *
   */
  ~MappingQCache();

  /**
   * clone()功能。有关文档，请参见 Mapping::clone(). 。
   *
   */
  virtual std::unique_ptr<Mapping<dim, spacedim>>
  clone() const override;

  /**
   * 返回 @p false ，因为顶点位置的保存取决于交给 reinit()
   * 函数的映射。
   *
   */
  virtual bool
  preserves_vertex_locations() const override;

  /**
   * 通过计算给定三角形的所有单元（在所有层次上）的映射支持点来初始化数据缓存。
   * @note 在底层三角图的信号 Triangulation::Signals::any_change
   * 时，缓存失效。
   *
   */
  void
  initialize(const Mapping<dim, spacedim> &      mapping,
             const Triangulation<dim, spacedim> &triangulation);

  /**
   * 通过计算给定三角形的所有单元（在所有层面）的映射支持点来初始化数据缓存。
   * @note  缓存在底层三角图的信号
   * Triangulation::Signals::any_change 时失效。      @deprecated
   * 使用上面的initialize()版本代替。
   *
   */
  DEAL_II_DEPRECATED void
  initialize(const Triangulation<dim, spacedim> &  triangulation,
             const MappingQGeneric<dim, spacedim> &mapping);

  /**
   * 通过让作为参数的函数提供给定三角形的所有单元（在所有层面）的映射支持点来初始化数据缓存。该函数必须返回一个`Point<spacedim>`的向量，其长度与多项式空间的大小相同，
   * $(p+1)^\text{dim}$  ，其中 $p$
   * 是映射的多项式程度，而且必须按照映射或FE_Q对其点的排序，即首先是所有
   * $2^\text{dim}$
   * 的顶点，然后是线、四边形和六边形的点，按照通常的分层编号。除了给定的点的数量之外，没有尝试在内部验证这些点。
   * @note
   * 如果启用了多线程，这个函数将并行运行，多次调用传入的函数。因此，在
   * MultithreadInfo::n_threads()>1,
   * 的情况下，用户代码必须确保该函数（通常是一个lambda）不会写进与其他线程共享的数据。
   * @note 缓存在底层三角化的信号
   * Triangulation::Signals::any_change 时被无效。
   *
   */
  void
  initialize(const Triangulation<dim, spacedim> &triangulation,
             const std::function<std::vector<Point<spacedim>>(
               const typename Triangulation<dim, spacedim>::cell_iterator &)>
               &compute_points_on_cell);

  /**
   * 通过计算给定三角形的所有单元（在所有层次上）和给定
   * @p mapping 的映射支持点并通过函数 @p transformation_function.
   * 转换这些点来初始化数据缓存。 bool  @p
   * function_describes_relative_displacement 表示函数 @p
   * transformation_function 映射到绝对坐标。
   * 如果参数设置为真，函数的返回值将被解释为相对变形，其结果最终将被添加到原始点上，用于本类最终使用的支持点。
   * 这个函数调用前一个函数，所以上面列出的关于线程的评论也适用于此。
   * @note 缓存在底层三角函数的信号
   * Triangulation::Signals::any_change 时被废止。
   *
   */
  void
  initialize(const Mapping<dim, spacedim> &      mapping,
             const Triangulation<dim, spacedim> &tria,
             const std::function<Point<spacedim>(
               const typename Triangulation<dim, spacedim>::cell_iterator &,
               const Point<spacedim> &)> &       transformation_function,
             const bool function_describes_relative_displacement);

  /**
   * 与上述相同，但取一个 dealii::Function 对象。
   *
   */
  void
  initialize(const Mapping<dim, spacedim> &      mapping,
             const Triangulation<dim, spacedim> &tria,
             const Function<spacedim> &          transformation_function,
             const bool function_describes_relative_displacement);

  /**
   * 通过一个离散字段（由 @p dof_handler 和 @p vector)
   * 指定，描述每个支持点的绝对或相对位置，初始化活动单元的数据缓存。
   * @note
   * 通过使用这个函数进行重新初始化，这个类的行为就像MappingFEField（vector_describes_relative_displacement
   * == false）或MappingQEulerian（vector_describes_relative_displacement
   * == true），但内部的操作要高效得多。
   *
   */
  template <typename VectorType>
  void
  initialize(const Mapping<dim, spacedim> &   mapping,
             const DoFHandler<dim, spacedim> &dof_handler,
             const VectorType &               vector,
             const bool vector_describes_relative_displacement);

  /**
   * 通过一个描述每个支持点的绝对或相对位置的解决方案（由
   * @p dof_handler 和一组 @p vectors
   * 在三角测量的所有层面上指定）来初始化所有非人工单元的数据缓存。
   * @note
   * 通过使用这个函数进行重新初始化，这个类的行为就像MappingFEField（vector_describes_relative_displacement
   * == false）或MappingQEulerian（vector_describes_relative_displacement
   * == true），但内部的操作要高效得多。
   *
   */
  template <typename VectorType>
  void
  initialize(const Mapping<dim, spacedim> &   mapping,
             const DoFHandler<dim, spacedim> &dof_handler,
             const MGLevelObject<VectorType> &vectors,
             const bool vector_describes_relative_displacement);

  /**
   * 返回缓存的内存消耗（以字节为单位）。
   *
   */
  std::size_t
  memory_consumption() const;

protected:
  /**
   * 这是从基类MappingQGeneric中重写的主要函数。
   *
   */
  virtual std::vector<Point<spacedim>>
  compute_mapping_support_points(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell)
    const override;

private:
  /**
   * 调用initialize()时填充的点缓存。它被做成一个共享指针，以允许多个实例（通过clone()创建）共享这个缓存。
   *
   */
  std::shared_ptr<std::vector<std::vector<std::vector<Point<spacedim>>>>>
    support_point_cache;

  /**
   * 与 Triangulation::signals::any
   * 的连接，一旦这个类超出范围就必须被重置。
   *
   */
  boost::signals2::connection clear_signal;

  /**
   * 指定是否已经为层次上的单元设置了support_point_cache。
   *
   */
  bool uses_level_info;
};

 /*@}*/ 

DEAL_II_NAMESPACE_CLOSE

#endif


