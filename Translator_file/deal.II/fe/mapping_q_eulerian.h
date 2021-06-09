//include/deal.II-translator/fe/mapping_q_eulerian_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2001 - 2021 by the deal.II authors
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


#ifndef dealii_mapping_q_eulerian_h
#define dealii_mapping_q_eulerian_h

#include <deal.II/base/config.h>

#include <deal.II/base/smartpointer.h>
#include <deal.II/base/thread_management.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/tria_iterator.h>


DEAL_II_NAMESPACE_OPEN

template <typename>
class Vector;


 /*!@addtogroup mapping */ 
 /*@{*/ 

/**
 * 这个类是MappingQ1Eulerian类对高阶 $Q_p$ 映射的扩展。
 * 当人们想在一个随着计算进行而变形的域上计算形状函数信息时，它很有用。
 * <h3>Usage</h3>
 * 这个类的构造函数需要三个参数：所需Qp映射的多项式程度，对定义从初始配置到当前配置的映射的向量的引用，以及对DoFHandler的引用。最常见的情况是使用所考虑的问题的解决方案向量作为转移向量。关键的要求是，给定矢量场的分量数量必须等于（或可能大于）空间维度的数量。如果分量多于空间维数（例如，如果我们正在处理一个耦合问题，其中有额外的求解变量），则假定第一个<tt>dim</tt>分量代表位移场，其余分量被忽略。
 * 如果这个假设不成立，我们可能需要在三角结构上设置一个单独的DoFHandler，并将所需的位移矢量与之关联。
 * 通常情况下，DoFHandler操作的有限元是由连续的FE_Q对象构成的系统元（FESystem）。下面是一个例子。
 *
 * @code
 *  FESystem<dim> fe(FE_Q<dim>(2), dim, FE_Q<dim>(1), 1);
 *  DoFHandler<dim> dof_handler(triangulation);
 *  dof_handler.distribute_dofs(fe);
 *  Vector<double> displacement_field(dof_handler.n_dofs());
 *  // ... compute displacement field somehow...
 *  MappingQEulerian<dim> q2_mapping(2, dof_handler, displacement_field);
 * @endcode
 *
 * 在这个例子中，我们的元素由<tt>(dim+1)</tt>组件组成。然而，只有第一个<tt>dim</tt>组件将被用于定义Q2映射。
 * 其余的组件将被忽略。
 * 注意，在构建映射对象之前，必须调用distribution_dofs(...)函数。
 * 还要注意的是，由于移位值向量和dof处理程序只在构造时与此对象相关联，所以你必须确保无论何时使用此对象，给定的对象仍然代表有效的数据。
 * 为了使MappingQEulerian类也能在使用PETSc或Trilinos包装类的并行代码中使用，矢量的类型可以被指定为模板参数<tt>VectorType</tt>。
 *
 *
 */
template <int dim, typename VectorType = Vector<double>, int spacedim = dim>
class MappingQEulerian : public MappingQ<dim, spacedim>
{
public:
  /**
   * 构造函数。      @param[in]  degree 所需 $Q_p$
   * 映射的多项式程度。    @param[in]  euler_dof_handler
   * 一个DoFHandler对象，它定义了一个有限元空间。这个空间需要至少有dim分量，空间的第一个dim分量将被视为相对于三角形单元的原始位置的位移。
   * @param[in]  euler_vector
   * 第二个参数定义的空间中的有限元函数。这个函数的第一个dim分量将被解释为我们在定义映射时使用的位移，相对于底层三角板的单元位置。
   * @param[in]  级别
   * 将使用映射的多网格级别。它主要用于检查 @p euler_vector
   * 的大小是否与 @p euler_dof_handler 一致。
   *
   */
  MappingQEulerian(const unsigned int               degree,
                   const DoFHandler<dim, spacedim> &euler_dof_handler,
                   const VectorType &               euler_vector,
                   const unsigned int level = numbers::invalid_unsigned_int);

  /**
   * 返回单元格的映射顶点。对于当前类，这个函数不使用当前单元格几何形状中的支持点，而是在单元格的几何形状之外，评估一个外部给定的位移场。
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
   * ，因为MappingQEulerian一般不保留顶点位置（除非翻译矢量碰巧在顶点位置提供零位移）。
   *
   */
  virtual bool
  preserves_vertex_locations() const override;

  // for documentation, see the Mapping base class
  virtual BoundingBox<spacedim>
  get_bounding_box(const typename Triangulation<dim, spacedim>::cell_iterator
                     &cell) const override;

  /**
   * 当映射在非活动单元被评估时，会抛出异常。
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
   * 对移位矢量的引用。
   *
   */
  SmartPointer<const VectorType, MappingQEulerian<dim, VectorType, spacedim>>
    euler_vector;

  /**
   * 指向映射向量所关联的DoFHandler的指针。
   *
   */
  SmartPointer<const DoFHandler<dim, spacedim>,
               MappingQEulerian<dim, VectorType, spacedim>>
    euler_dof_handler;


private:
  /**
   * 要使用映射的多网格级别。
   *
   */
  const unsigned int level;

  /**
   * 一个派生自MappingQGeneric的类，它为通用映射提供边界对象上的支持点，从而使相应的Q3映射最终成为C1。
   *
   */
  class MappingQEulerianGeneric : public MappingQGeneric<dim, spacedim>
  {
  public:
    /**
     * 构造函数。
     *
     */
    MappingQEulerianGeneric(
      const unsigned int                                 degree,
      const MappingQEulerian<dim, VectorType, spacedim> &mapping_q_eulerian);

    /**
     * 返回单元格的映射顶点。对于当前的类，这个函数不使用当前单元格的几何形状中的支持点，而是在单元格的几何形状之外，评估一个外部给定的位移场。
     *
     */
    virtual boost::container::small_vector<Point<spacedim>,
                                           GeometryInfo<dim>::vertices_per_cell>
    get_vertices(const typename Triangulation<dim, spacedim>::cell_iterator
                   &cell) const override;

    /**
     * 计算当前配置中的支撑点的位置。更多信息请参见
     * MappingQGeneric::compute_mapping_support_points() 的文档。
     *
     */
    virtual std::vector<Point<spacedim>>
    compute_mapping_support_points(
      const typename Triangulation<dim, spacedim>::cell_iterator &cell)
      const override;

    /**
     * 总是返回 @p false
     * ，因为MappingQEulerianGeneric一般不保留顶点位置（除非平移矢量刚好规定了顶点位置的零位移）。
     *
     */
    virtual bool
    preserves_vertex_locations() const override;

  private:
    /**
     * 对我们所处的周围物体的参考。
     *
     */
    const MappingQEulerian<dim, VectorType, spacedim> &mapping_q_eulerian;


    /**
     * 用于定义参考配置中的支持点的特殊正交规则。
     *
     */
    class SupportQuadrature : public Quadrature<dim>
    {
    public:
      /**
       * 构造函数，有一个参数定义了所需的多项式程度。
       *
       */
      SupportQuadrature(const unsigned int map_degree);
    };

    /**
     * 一个成员变量，按正确顺序保存正交点。
     *
     */
    const SupportQuadrature support_quadrature;

    /**
     * FEValues对象，用于查询参考配置中支持点的给定有限元场。
     * 该变量被标记为可变的，因为我们必须从compute_mapping_support_points中调用
     * FEValues::reinit ，这个函数是'常量'。
     *
     */
    mutable FEValues<dim, spacedim> fe_values;

    /**
     * 一个变量来保护对fe_values变量的访问。
     *
     */
    mutable Threads::Mutex fe_values_mutex;
  };
};

 /*@}*/ 


 /*----------------------------------------------------------------------*/ 

#ifndef DOXYGEN

template <int dim, typename VectorType, int spacedim>
inline bool
MappingQEulerian<dim, VectorType, spacedim>::preserves_vertex_locations() const
{
  return false;
}

#endif // DOXYGEN


DEAL_II_NAMESPACE_CLOSE


#endif // dealii_mapping_q_eulerian_h


