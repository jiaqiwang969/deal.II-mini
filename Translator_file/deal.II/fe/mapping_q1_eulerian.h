//include/deal.II-translator/fe/mapping_q1_eulerian_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2001 - 2020 by the deal.II authors
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

#ifndef dealii_mapping_q1_eulerian_h
#define dealii_mapping_q1_eulerian_h

#include <deal.II/base/config.h>

#include <deal.II/base/smartpointer.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/mapping_q1.h>

#include <array>

DEAL_II_NAMESPACE_OPEN

template <typename>
class Vector;


 /*!@addtogroup mapping */ 
 /*@{*/ 

/**
 * 该类提供了一个映射，在每个单元的位置上添加一个 $d$
 * -线性位移场。对高阶多项式的概括在MappingQEulerian类中提供）。因此，每个单元在空间中的位移是通过有限元场给映射的值来实现的。
 * <h3>Usage</h3>
 * 这个类的构造函数需要两个参数：一个是对定义从参考配置到当前配置的映射的向量的引用，一个是对DoFHandler的引用。然后，该向量应该代表一个在DoFHandler定义的节点上定义的（扁平化版本的）向量值场，其中向量场的分量数量等于空间维度的数量。因此，DoFHandler应该在一个有限元上操作，该有限元的分量与空间维数一样多。作为一个额外的要求，我们规定它的每个顶点的自由度与空间维数一样多；因为这个对象只在顶点评估有限元场，所有其他自由度的值（不与顶点相关）都被忽略了。如果给定的DoFHandler操作的有限元是由
 * @p dim
 * 连续的FE_Q()对象构造成系统元(FESystem)，则满足这些要求。
 * 在许多情况下，移位矢量也将是所研究问题的解决矢量。如果不是这种情况（即解变量的分量数不等于空间维度，例如对于<tt>dim>1</tt>中的标量问题，欧拉坐标只给出一个背景场），或者对于需要计算更多变量而不仅仅是流场的耦合问题），那么必须在给定的三角形上设置一个不同的DoFHandler，然后将位移向量与之关联。
 * 下面是一个例子。
 *
 * @code
 *  FESystem<dim> fe(FE_Q<dim>(1), dim);
 *  DoFHandler<dim> flowfield_dof_handler(triangulation);
 *  flowfield_dof_handler.distribute_dofs(fe);
 *  Vector<double> displacement_field(flowfield_dof_handler.n_dofs());
 *  MappingQ1Eulerian<dim> mymapping(flowfield_dof_handler,
 *                                   displacement_field);
 * @endcode
 *
 * 请注意，由于移位值向量和dof处理程序只在构造时与这个对象相关联，你必须确保无论何时使用这个对象，给定的对象仍然代表有效数据。
 * 为了使MappingQ1Eulerian类也能在使用PETSc或Trilinos包装类的并行代码中使用，矢量的类型可以被指定为模板参数<tt>VectorType</tt>。
 * 关于<tt>spacedim</tt>模板参数的更多信息，请查阅FiniteElement或Triangulation的文档。
 *
 *
 */
template <int dim, typename VectorType = Vector<double>, int spacedim = dim>
class MappingQ1Eulerian : public MappingQGeneric<dim, spacedim>
{
public:
  /**
   * 构造函数。      @param[in]  euler_dof_handler
   * 一个DoFHandler对象，定义了一个有限元空间。这个空间需要有精确的dim分量，这些分量将被视为相对于三角形的单元的原始位置的位移。
   * 这个DoFHandler必须基于一个 <code>FESystem(FE_Q(1),dim)</code>
   * 有限元。    @param[in]  euler_vector
   * 在第一个参数定义的空间中的一个有限元函数。这个函数的dim分量将被解释为我们在定义映射时使用的位移，相对于底层三角结构的单元位置。
   *
   */
  MappingQ1Eulerian(const DoFHandler<dim, spacedim> &euler_dof_handler,
                    const VectorType &               euler_vector);

  /**
   * 返回单元格的映射顶点。对于当前类，这个函数不使用当前单元格的几何形状中的支持点，而是在单元格的几何形状之外，评估一个外部给定的位移场。
   *
   */
  virtual boost::container::small_vector<Point<spacedim>,
                                         GeometryInfo<dim>::vertices_per_cell>
  get_vertices(const typename Triangulation<dim, spacedim>::cell_iterator &cell)
    const override;

  /**
   * 返回一个指向当前对象副本的指针。这个副本的调用者就拥有了它的所有权。
   *
   */
  virtual std::unique_ptr<Mapping<dim, spacedim>>
  clone() const override;

  /**
   * 总是返回 @p false
   * ，因为MappingQ1Eulerian一般不保留顶点位置（除非翻译矢量恰好规定顶点位置的位移为零）。
   *
   */
  virtual bool
  preserves_vertex_locations() const override;

  /**
   * 异常情况。
   *
   */
  DeclException0(ExcInactiveCell);



protected:
  /**
   * 计算一个单元的映射相关信息。关于这个函数的目的、参数和返回值的讨论，请参见
   * Mapping::fill_fe_values() 的文档。
   * 这个函数覆盖了基类中的函数，因为我们不能为这个类使用任何单元格的相似性。
   *
   */
  virtual CellSimilarity::Similarity
  fill_fe_values(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const CellSimilarity::Similarity                            cell_similarity,
    const Quadrature<dim> &                                     quadrature,
    const typename Mapping<dim, spacedim>::InternalDataBase &   internal_data,
    internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
      &output_data) const override;

  /**
   * 计算映射的支持点。对于当前类，这些是顶点，通过调用
   * Mapping::get_vertices(). 得到的，更多信息请参见
   * MappingQGeneric::compute_mapping_support_points() 的文档。
   *
   */
  virtual std::vector<Point<spacedim>>
  compute_mapping_support_points(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell)
    const override;

  /**
   * 对移位矢量的引用。
   *
   */
  SmartPointer<const VectorType, MappingQ1Eulerian<dim, VectorType, spacedim>>
    euler_transform_vectors;

  /**
   * 指向映射向量所关联的DoFHandler的指针。
   *
   */
  SmartPointer<const DoFHandler<dim, spacedim>,
               MappingQ1Eulerian<dim, VectorType, spacedim>>
    shiftmap_dof_handler;
};

 /*@}*/ 

 /*----------------------------------------------------------------------*/ 

#ifndef DOXYGEN

template <int dim, typename VectorType, int spacedim>
inline bool
MappingQ1Eulerian<dim, VectorType, spacedim>::preserves_vertex_locations() const
{
  return false;
}

#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif


