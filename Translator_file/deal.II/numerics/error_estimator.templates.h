//include/deal.II-translator/numerics/error_estimator.templates_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2021 by the deal.II authors
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

#ifndef dealii_error_estimator_templates_h
#define dealii_error_estimator_templates_h

#include <deal.II/base/config.h>

#include <deal.II/base/geometry_info.h>
#include <deal.II/base/numbers.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/work_stream.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_update_flags.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/tria_iterator.h>

#include <deal.II/hp/fe_values.h>
#include <deal.II/hp/mapping_collection.h>
#include <deal.II/hp/q_collection.h>

#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/la_parallel_block_vector.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/la_vector.h>
#include <deal.II/lac/petsc_block_vector.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/error_estimator.h>

#include <algorithm>
#include <cmath>
#include <functional>
#include <numeric>
#include <vector>

DEAL_II_NAMESPACE_OPEN

namespace internal
{
  /**
   * 误差估计器的几个函数在每个线程中需要一次的所有小型临时数据对象都聚集在这个结构中。
   * 建立这个结构的原因主要是我们有一些对单元或面进行操作的函数，需要一些小的临时数据对象。由于这些函数可能会并行运行，我们不能让这些对象成为包围类的成员变量。另一方面，在这些函数中的每一个局部声明它们将需要在我们每次访问下一个单元格或面时重新分配，我们发现如果经常发生这种情况，即使是在单线程的情况下（在我们的测量中为10-20%）也会花费大量的时间；然而，最重要的是，内存分配需要在多线程模式下进行同步。虽然这是由C++库完成的，不需要手工编码，但它还是严重损害了有效地并行运行该类函数的能力，因为它们经常被这些同步点阻断，使一切都慢了两到三倍。
   * 因此，每个线程都有一个这个类的实例来工作，不需要自己分配内存，也不需要与其他线程同步。
   * 数组的大小被初始化为hp情况下所需的最大条目数。在各个单元的循环中，我们会根据需要调整数组的大小。因为对于
   * std::vector
   * 来说，调整大小到一个较小的尺寸并不意味着内存分配，所以这很快速。
   *
   */
  template <int dim, int spacedim, typename number>
  struct ParallelData
  {
    /**
     * 要使用的有限元。
     *
     */
    const dealii::hp::FECollection<dim, spacedim> finite_element;

    /**
     * 用于面的正交公式。
     *
     */
    const dealii::hp::QCollection<dim - 1> face_quadratures;

    /**
     * FEFaceValues对象，用于对当前和可能的邻近单元的面进行积分。
     *
     */
    dealii::hp::FEFaceValues<dim, spacedim>    fe_face_values_cell;
    dealii::hp::FEFaceValues<dim, spacedim>    fe_face_values_neighbor;
    dealii::hp::FESubfaceValues<dim, spacedim> fe_subface_values;

    /**
     * 一个向量用于存储每个求解向量的正交点的法向量跳动（即一个临时值）。
     * 这个向量不是在使用它的函数内部分配的，而是全局分配的，因为内存分配很慢，特别是在有多个线程的情况下，同步化会使事情变得更慢。
     *
     */
    std::vector<std::vector<std::vector<number>>> phi;

    /**
     * 一个单元上的有限元函数梯度的向量 让psi成为<tt>a grad
     * u_h</tt>的简称，其中第三个索引是有限元的分量，第二个索引是正交点的编号。第一个索引表示解向量的索引。
     *
     */
    std::vector<std::vector<std::vector<Tensor<1, spacedim, number>>>> psi;

    /**
     * 邻近单元的相同向量
     *
     */
    std::vector<std::vector<std::vector<Tensor<1, spacedim, number>>>>
      neighbor_psi;

    /**
     * 一个面上的有限元函数的法向量
     *
     */
    std::vector<Tensor<1, spacedim>> normal_vectors;

    /**
     * 对立面的法向量。
     *
     */
    std::vector<Tensor<1, spacedim>> neighbor_normal_vectors;

    /**
     * 跳跃中的系数值所需的两个数组，如果它们被给出的话。
     *
     */
    std::vector<double>                 coefficient_values1;
    std::vector<dealii::Vector<double>> coefficient_values;

    /**
     * 用于雅各布行列式和四分仪点的权重的乘积的数组。
     *
     */
    std::vector<double> JxW_values;

    /**
     * 我们要关注的子域ID。
     *
     */
    const types::subdomain_id subdomain_id;
    /**
     * 我们要关注的材料ID。
     *
     */
    const types::material_id material_id;

    /**
     * 还有一些对 KellyErrorEstimator::estimate()
     * 函数的输入数据的引用。
     *
     */
    const std::map<types::boundary_id, const Function<spacedim, number> *>
      *                       neumann_bc;
    const ComponentMask       component_mask;
    const Function<spacedim> *coefficients;

    /**
     * 构造函数。
     *
     */
    template <class FE>
    ParallelData(const FE &                              fe,
                 const dealii::hp::QCollection<dim - 1> &face_quadratures,
                 const dealii::hp::MappingCollection<dim, spacedim> &mapping,
                 const bool                need_quadrature_points,
                 const unsigned int        n_solution_vectors,
                 const types::subdomain_id subdomain_id,
                 const types::material_id  material_id,
                 const std::map<types::boundary_id,
                                const Function<spacedim, number> *> *neumann_bc,
                 const ComponentMask &     component_mask,
                 const Function<spacedim> *coefficients);

    /**
     * 调整数组的大小，使其适合与给定的有限元索引相关的正交点数量的hp-集合。
     *
     */
    void
    resize(const unsigned int active_fe_index);
  };


  template <int dim, int spacedim, typename number>
  template <class FE>
  ParallelData<dim, spacedim, number>::ParallelData(
    const FE &                                          fe,
    const dealii::hp::QCollection<dim - 1> &            face_quadratures,
    const dealii::hp::MappingCollection<dim, spacedim> &mapping,
    const bool                                          need_quadrature_points,
    const unsigned int                                  n_solution_vectors,
    const types::subdomain_id                           subdomain_id,
    const types::material_id                            material_id,
    const std::map<types::boundary_id, const Function<spacedim, number> *>
      *                       neumann_bc,
    const ComponentMask &     component_mask,
    const Function<spacedim> *coefficients)
    : finite_element(fe)
    , face_quadratures(face_quadratures)
    , fe_face_values_cell(mapping,
                          finite_element,
                          face_quadratures,
                          update_gradients | update_JxW_values |
                            (need_quadrature_points ? update_quadrature_points :
                                                      UpdateFlags()) |
                            update_normal_vectors)
    , fe_face_values_neighbor(mapping,
                              finite_element,
                              face_quadratures,
                              update_gradients | update_normal_vectors)
    , fe_subface_values(mapping,
                        finite_element,
                        face_quadratures,
                        update_gradients | update_normal_vectors)
    , phi(n_solution_vectors,
          std::vector<std::vector<number>>(
            face_quadratures.max_n_quadrature_points(),
            std::vector<number>(fe.n_components())))
    , psi(n_solution_vectors,
          std::vector<std::vector<Tensor<1, spacedim, number>>>(
            face_quadratures.max_n_quadrature_points(),
            std::vector<Tensor<1, spacedim, number>>(fe.n_components())))
    , neighbor_psi(n_solution_vectors,
                   std::vector<std::vector<Tensor<1, spacedim, number>>>(
                     face_quadratures.max_n_quadrature_points(),
                     std::vector<Tensor<1, spacedim, number>>(
                       fe.n_components())))
    , normal_vectors(face_quadratures.max_n_quadrature_points())
    , neighbor_normal_vectors(face_quadratures.max_n_quadrature_points())
    , coefficient_values1(face_quadratures.max_n_quadrature_points())
    , coefficient_values(face_quadratures.max_n_quadrature_points(),
                         dealii::Vector<double>(fe.n_components()))
    , JxW_values(face_quadratures.max_n_quadrature_points())
    , subdomain_id(subdomain_id)
    , material_id(material_id)
    , neumann_bc(neumann_bc)
    , component_mask(component_mask)
    , coefficients(coefficients)
  {}



  template <int dim, int spacedim, typename number>
  void
  ParallelData<dim, spacedim, number>::resize(
    const unsigned int active_fe_index)
  {
    const unsigned int n_q_points   = face_quadratures[active_fe_index].size();
    const unsigned int n_components = finite_element.n_components();

    normal_vectors.resize(n_q_points);
    neighbor_normal_vectors.resize(n_q_points);
    coefficient_values1.resize(n_q_points);
    coefficient_values.resize(n_q_points);
    JxW_values.resize(n_q_points);

    for (unsigned int i = 0; i < phi.size(); ++i)
      {
        phi[i].resize(n_q_points);
        psi[i].resize(n_q_points);
        neighbor_psi[i].resize(n_q_points);

        for (unsigned int qp = 0; qp < n_q_points; ++qp)
          {
            phi[i][qp].resize(n_components);
            psi[i][qp].resize(n_components);
            neighbor_psi[i][qp].resize(n_components);
          }
      }

    for (unsigned int qp = 0; qp < n_q_points; ++qp)
      coefficient_values[qp].reinit(n_components);
  }



  /**
   * 将单个ParallelData对象的local_face_integrals地图中的数据复制到一个全局的此类地图中。这是一个WorkStream流水线的复制器阶段。
   *
   */
  template <int dim, int spacedim>
  void
  copy_local_to_global(
    const std::map<typename DoFHandler<dim, spacedim>::face_iterator,
                   std::vector<double>> &local_face_integrals,
    std::map<typename DoFHandler<dim, spacedim>::face_iterator,
             std::vector<double>> &      face_integrals)
  {
    // now copy locally computed elements into the global map
    for (typename std::map<typename DoFHandler<dim, spacedim>::face_iterator,
                           std::vector<double>>::const_iterator p =
           local_face_integrals.begin();
         p != local_face_integrals.end();
         ++p)
      {
        // double check that the element does not already exists in the
        // global map
        Assert(face_integrals.find(p->first) == face_integrals.end(),
               ExcInternalError());

        for (unsigned int i = 0; i < p->second.size(); ++i)
          {
            Assert(numbers::is_finite(p->second[i]), ExcInternalError());
            Assert(p->second[i] >= 0, ExcInternalError());
          }

        face_integrals[p->first] = p->second;
      }
  }


  /**
   * 实际上是根据ParallelData中评估的梯度来进行计算。
   *
   */
  template <int dim, int spacedim, typename number>
  std::vector<double>
  integrate_over_face(
    ParallelData<dim, spacedim, number> &                    parallel_data,
    const typename DoFHandler<dim, spacedim>::face_iterator &face,
    dealii::hp::FEFaceValues<dim, spacedim> &fe_face_values_cell)
  {
    const unsigned int n_q_points = parallel_data.psi[0].size(),
                       n_components =
                         parallel_data.finite_element.n_components(),
                       n_solution_vectors = parallel_data.psi.size();

    // now psi contains the following:
    // - for an internal face, psi=[grad u]
    // - for a neumann boundary face, psi=grad u
    // each component being the mentioned value at one of the quadrature
    // points

    // next we have to multiply this with the normal vector. Since we have
    // taken the difference of gradients for internal faces, we may chose
    // the normal vector of one cell, taking that of the neighbor would only
    // change the sign. We take the outward normal.

    parallel_data.normal_vectors =
      fe_face_values_cell.get_present_fe_values().get_normal_vectors();

    for (unsigned int n = 0; n < n_solution_vectors; ++n)
      for (unsigned int component = 0; component < n_components; ++component)
        for (unsigned int point = 0; point < n_q_points; ++point)
          parallel_data.phi[n][point][component] =
            (parallel_data.psi[n][point][component] *
             parallel_data.normal_vectors[point]);

    if (face->at_boundary() == false)
      {
        // compute the jump in the gradients

        for (unsigned int n = 0; n < n_solution_vectors; ++n)
          for (unsigned int component = 0; component < n_components;
               ++component)
            for (unsigned int p = 0; p < n_q_points; ++p)
              parallel_data.phi[n][p][component] +=
                (parallel_data.neighbor_psi[n][p][component] *
                 parallel_data.neighbor_normal_vectors[p]);
      }

    // if a coefficient was given: use that to scale the jump in the
    // gradient
    if (parallel_data.coefficients != nullptr)
      {
        // scalar coefficient
        if (parallel_data.coefficients->n_components == 1)
          {
            parallel_data.coefficients->value_list(
              fe_face_values_cell.get_present_fe_values()
                .get_quadrature_points(),
              parallel_data.coefficient_values1);
            for (unsigned int n = 0; n < n_solution_vectors; ++n)
              for (unsigned int component = 0; component < n_components;
                   ++component)
                for (unsigned int point = 0; point < n_q_points; ++point)
                  parallel_data.phi[n][point][component] *=
                    parallel_data.coefficient_values1[point];
          }
        else
          // vector-valued coefficient
          {
            parallel_data.coefficients->vector_value_list(
              fe_face_values_cell.get_present_fe_values()
                .get_quadrature_points(),
              parallel_data.coefficient_values);
            for (unsigned int n = 0; n < n_solution_vectors; ++n)
              for (unsigned int component = 0; component < n_components;
                   ++component)
                for (unsigned int point = 0; point < n_q_points; ++point)
                  parallel_data.phi[n][point][component] *=
                    parallel_data.coefficient_values[point](component);
          }
      }


    if (face->at_boundary() == true)
      // neumann boundary face. compute difference between normal derivative
      // and boundary function
      {
        const types::boundary_id boundary_id = face->boundary_id();

        Assert(parallel_data.neumann_bc->find(boundary_id) !=
                 parallel_data.neumann_bc->end(),
               ExcInternalError());
        // get the values of the boundary function at the quadrature points
        if (n_components == 1)
          {
            std::vector<number> g(n_q_points);
            parallel_data.neumann_bc->find(boundary_id)
              ->second->value_list(fe_face_values_cell.get_present_fe_values()
                                     .get_quadrature_points(),
                                   g);

            for (unsigned int n = 0; n < n_solution_vectors; ++n)
              for (unsigned int point = 0; point < n_q_points; ++point)
                parallel_data.phi[n][point][0] -= g[point];
          }
        else
          {
            std::vector<dealii::Vector<number>> g(
              n_q_points, dealii::Vector<number>(n_components));
            parallel_data.neumann_bc->find(boundary_id)
              ->second->vector_value_list(fe_face_values_cell
                                            .get_present_fe_values()
                                            .get_quadrature_points(),
                                          g);

            for (unsigned int n = 0; n < n_solution_vectors; ++n)
              for (unsigned int component = 0; component < n_components;
                   ++component)
                for (unsigned int point = 0; point < n_q_points; ++point)
                  parallel_data.phi[n][point][component] -= g[point](component);
          }
      }



    // now phi contains the following:
    // - for an internal face, phi=[a du/dn]
    // - for a neumann boundary face, phi=a du/dn-g
    // each component being the mentioned value at one of the quadrature
    // points

    parallel_data.JxW_values =
      fe_face_values_cell.get_present_fe_values().get_JxW_values();

    // take the square of the phi[i] for integration, and sum up
    std::vector<double> face_integral(n_solution_vectors, 0);
    for (unsigned int n = 0; n < n_solution_vectors; ++n)
      for (unsigned int component = 0; component < n_components; ++component)
        if (parallel_data.component_mask[component] == true)
          for (unsigned int p = 0; p < n_q_points; ++p)
            face_integral[n] += numbers::NumberTraits<number>::abs_square(
                                  parallel_data.phi[n][p][component]) *
                                parallel_data.JxW_values[p];

    return face_integral;
  }

  /**
   * 一个用于缩放边界处面的积分的因子。用于Neumann BC。
   *
   */
  template <int dim, int spacedim>
  double
  boundary_face_factor(
    const typename DoFHandler<dim, spacedim>::active_cell_iterator &cell,
    const unsigned int                                              face_no,
    const dealii::hp::FEFaceValues<dim, spacedim> &fe_face_values_cell,
    const typename KellyErrorEstimator<dim, spacedim>::Strategy strategy)
  {
    switch (strategy)
      {
        case KellyErrorEstimator<dim, spacedim>::cell_diameter_over_24:
          {
            return 1.0;
          }
        case KellyErrorEstimator<dim, spacedim>::cell_diameter:
          {
            return 1.0;
          }
        case KellyErrorEstimator<dim,
                                 spacedim>::face_diameter_over_twice_max_degree:
          {
            const double cell_degree =
              fe_face_values_cell.get_fe_collection()[cell->active_fe_index()]
                .degree;
            return cell->face(face_no)->diameter() / cell_degree;
          }
        default:
          {
            Assert(false, ExcNotImplemented());
            return -std::numeric_limits<double>::max();
          }
      }
  }


  /**
   * 一个用于缩放常规面的积分的因子。
   *
   */
  template <int dim, int spacedim>
  double
  regular_face_factor(
    const typename DoFHandler<dim, spacedim>::active_cell_iterator &cell,
    const unsigned int                                              face_no,
    const dealii::hp::FEFaceValues<dim, spacedim> &fe_face_values_cell,
    const dealii::hp::FEFaceValues<dim, spacedim> &fe_face_values_neighbor,
    const typename KellyErrorEstimator<dim, spacedim>::Strategy strategy)
  {
    switch (strategy)
      {
        case KellyErrorEstimator<dim, spacedim>::cell_diameter_over_24:
          {
            return 1.0;
          }
        case KellyErrorEstimator<dim, spacedim>::cell_diameter:
          {
            return 1.0;
          }
        case KellyErrorEstimator<dim,
                                 spacedim>::face_diameter_over_twice_max_degree:
          {
            const double cell_degree =
              fe_face_values_cell.get_fe_collection()[cell->active_fe_index()]
                .degree;
            const double neighbor_degree =
              fe_face_values_neighbor
                .get_fe_collection()[cell->neighbor(face_no)->active_fe_index()]
                .degree;
            return cell->face(face_no)->diameter() /
                   std::max(cell_degree, neighbor_degree) / 2.0;
          }
        default:
          {
            Assert(false, ExcNotImplemented());
            return -std::numeric_limits<double>::max();
          }
      }
  }

  /**
   * 对不规则面的积分进行缩放的系数。
   *
   */
  template <int dim, int spacedim>
  double
  irregular_face_factor(
    const typename DoFHandler<dim, spacedim>::active_cell_iterator &cell,
    const typename DoFHandler<dim, spacedim>::active_cell_iterator
      &                                            neighbor_child,
    const unsigned int                             face_no,
    const unsigned int                             subface_no,
    const dealii::hp::FEFaceValues<dim, spacedim> &fe_face_values,
    dealii::hp::FESubfaceValues<dim, spacedim> &   fe_subface_values,
    const typename KellyErrorEstimator<dim, spacedim>::Strategy strategy)
  {
    switch (strategy)
      {
        case KellyErrorEstimator<dim, spacedim>::cell_diameter_over_24:
          {
            return 1.0;
          }
        case KellyErrorEstimator<dim, spacedim>::cell_diameter:
          {
            return 1.0;
          }
        case KellyErrorEstimator<dim,
                                 spacedim>::face_diameter_over_twice_max_degree:
          {
            const double cell_degree =
              fe_face_values.get_fe_collection()[cell->active_fe_index()]
                .degree;
            const double neighbor_child_degree =
              fe_subface_values
                .get_fe_collection()[neighbor_child->active_fe_index()]
                .degree;
            return cell->face(face_no)->child(subface_no)->diameter() /
                   std::max(neighbor_child_degree, cell_degree) / 2.0;
          }
        default:
          {
            Assert(false, ExcNotImplemented());
            return -std::numeric_limits<double>::max();
          }
      }
  }

  /**
   * 一个用于汇总每个单元中不同面的所有贡献的系数。
   *
   */
  template <int dim, int spacedim>
  double
  cell_factor(
    const typename DoFHandler<dim, spacedim>::active_cell_iterator &cell,
    const unsigned int  /*face_no*/ ,
    const DoFHandler<dim, spacedim> &  /*dof_handler*/ ,
    const typename KellyErrorEstimator<dim, spacedim>::Strategy strategy)
  {
    switch (strategy)
      {
        case KellyErrorEstimator<dim, spacedim>::cell_diameter_over_24:
          {
            return cell->diameter() / 24;
          }
        case KellyErrorEstimator<dim, spacedim>::cell_diameter:
          {
            return cell->diameter();
          }
        case KellyErrorEstimator<dim,
                                 spacedim>::face_diameter_over_twice_max_degree:
          {
            return 1.0;
          }
        default:
          {
            Assert(false, ExcNotImplemented());
            return -std::numeric_limits<double>::max();
          }
      }
  }



  /**
   * 实际上是在一个没有悬挂节点的面（它是规则的）上进行计算，也就是说，要么在另一边有涅槃（面在边界），要么另一边的细化水平与这一边的细化水平相同，然后把这两种情况的积分一起处理。
   *
   */
  template <typename InputVector, int dim, int spacedim>
  void
  integrate_over_regular_face(
    const std::vector<const InputVector *> &solutions,
    ParallelData<dim, spacedim, typename InputVector::value_type>
      &                            parallel_data,
    std::map<typename DoFHandler<dim, spacedim>::face_iterator,
             std::vector<double>> &local_face_integrals,
    const typename DoFHandler<dim, spacedim>::active_cell_iterator &cell,
    const unsigned int                                              face_no,
    dealii::hp::FEFaceValues<dim, spacedim> &fe_face_values_cell,
    dealii::hp::FEFaceValues<dim, spacedim> &fe_face_values_neighbor,
    const typename KellyErrorEstimator<dim, spacedim>::Strategy strategy)
  {
    const typename DoFHandler<dim, spacedim>::face_iterator face =
      cell->face(face_no);
    const unsigned int n_solution_vectors = solutions.size();


    // initialize data of the restriction
    // of this cell to the present face
    fe_face_values_cell.reinit(cell, face_no, cell->active_fe_index());

    // get gradients of the finite element
    // function on this cell
    for (unsigned int n = 0; n < n_solution_vectors; ++n)
      fe_face_values_cell.get_present_fe_values().get_function_gradients(
        *solutions[n], parallel_data.psi[n]);

    double factor;
    // now compute over the other side of the face
    if (face->at_boundary() == false)
      // internal face; integrate jump of gradient across this face
      {
        Assert(cell->neighbor(face_no).state() == IteratorState::valid,
               ExcInternalError());

        const typename DoFHandler<dim, spacedim>::active_cell_iterator
          neighbor = cell->neighbor(face_no);

        // find which number the current face has relative to the
        // neighboring cell
        const unsigned int neighbor_neighbor =
          cell->neighbor_of_neighbor(face_no);
        Assert(neighbor_neighbor < GeometryInfo<dim>::faces_per_cell,
               ExcInternalError());

        // get restriction of finite element function of @p{neighbor} to the
        // common face. in the hp-case, use the quadrature formula that
        // matches the one we would use for the present cell
        fe_face_values_neighbor.reinit(neighbor,
                                       neighbor_neighbor,
                                       cell->active_fe_index());

        factor = regular_face_factor<dim, spacedim>(cell,
                                                    face_no,
                                                    fe_face_values_cell,
                                                    fe_face_values_neighbor,
                                                    strategy);

        // get gradients on neighbor cell
        for (unsigned int n = 0; n < n_solution_vectors; ++n)
          {
            fe_face_values_neighbor.get_present_fe_values()
              .get_function_gradients(*solutions[n],
                                      parallel_data.neighbor_psi[n]);
          }

        parallel_data.neighbor_normal_vectors =
          fe_face_values_neighbor.get_present_fe_values().get_normal_vectors();
      }
    else
      {
        factor = boundary_face_factor<dim, spacedim>(cell,
                                                     face_no,
                                                     fe_face_values_cell,
                                                     strategy);
      }

    // now go to the generic function that does all the other things
    local_face_integrals[face] =
      integrate_over_face(parallel_data, face, fe_face_values_cell);

    for (unsigned int i = 0; i < local_face_integrals[face].size(); i++)
      local_face_integrals[face][i] *= factor;
  }



  /**
   * 和上面的函数一样，只是积分是在 @p cell, 的 @p face_no
   * 面上，其中各自的邻居被细化，所以积分会更复杂一些。
   *
   */
  template <typename InputVector, int dim, int spacedim>
  void
  integrate_over_irregular_face(
    const std::vector<const InputVector *> &solutions,
    ParallelData<dim, spacedim, typename InputVector::value_type>
      &                            parallel_data,
    std::map<typename DoFHandler<dim, spacedim>::face_iterator,
             std::vector<double>> &local_face_integrals,
    const typename DoFHandler<dim, spacedim>::active_cell_iterator &cell,
    const unsigned int                                              face_no,
    dealii::hp::FEFaceValues<dim, spacedim> &   fe_face_values,
    dealii::hp::FESubfaceValues<dim, spacedim> &fe_subface_values,
    const typename KellyErrorEstimator<dim, spacedim>::Strategy strategy)
  {
    const auto neighbor = cell->neighbor(face_no);
    (void)neighbor;
    const unsigned int n_solution_vectors = solutions.size();
    const auto         face               = cell->face(face_no);

    Assert(neighbor.state() == IteratorState::valid, ExcInternalError());
    Assert(face->has_children(), ExcInternalError());

    // set up a vector of the gradients of the finite element function on
    // this cell at the quadrature points
    //
    // let psi be a short name for [a grad u_h], where the second index be
    // the component of the finite element, and the first index the number
    // of the quadrature point

    // store which number @p{cell} has in the list of neighbors of
    // @p{neighbor}
    const unsigned int neighbor_neighbor = cell->neighbor_of_neighbor(face_no);
    Assert(neighbor_neighbor < GeometryInfo<dim>::faces_per_cell,
           ExcInternalError());

    // loop over all subfaces
    for (unsigned int subface_no = 0; subface_no < face->n_children();
         ++subface_no)
      {
        // get an iterator pointing to the cell behind the present subface
        const typename DoFHandler<dim, spacedim>::active_cell_iterator
          neighbor_child = cell->neighbor_child_on_subface(face_no, subface_no);
        Assert(neighbor_child->is_active(), ExcInternalError());

        // restrict the finite element on the present cell to the subface
        fe_subface_values.reinit(cell,
                                 face_no,
                                 subface_no,
                                 cell->active_fe_index());

        // restrict the finite element on the neighbor cell to the common
        // @p{subface}.
        fe_face_values.reinit(neighbor_child,
                              neighbor_neighbor,
                              cell->active_fe_index());

        const double factor =
          irregular_face_factor<dim, spacedim>(cell,
                                               neighbor_child,
                                               face_no,
                                               subface_no,
                                               fe_face_values,
                                               fe_subface_values,
                                               strategy);

        // store the gradient of the solution in psi
        for (unsigned int n = 0; n < n_solution_vectors; ++n)
          fe_subface_values.get_present_fe_values().get_function_gradients(
            *solutions[n], parallel_data.psi[n]);

        // store the gradient from the neighbor's side in @p{neighbor_psi}
        for (unsigned int n = 0; n < n_solution_vectors; ++n)
          fe_face_values.get_present_fe_values().get_function_gradients(
            *solutions[n], parallel_data.neighbor_psi[n]);

        // call generic evaluate function
        parallel_data.neighbor_normal_vectors =
          fe_subface_values.get_present_fe_values().get_normal_vectors();

        local_face_integrals[neighbor_child->face(neighbor_neighbor)] =
          integrate_over_face(parallel_data, face, fe_face_values);
        for (unsigned int i = 0;
             i < local_face_integrals[neighbor_child->face(neighbor_neighbor)]
                   .size();
             i++)
          local_face_integrals[neighbor_child->face(neighbor_neighbor)][i] *=
            factor;
      }

    // finally loop over all subfaces to collect the contributions of the
    // subfaces and store them with the mother face
    std::vector<double> sum(n_solution_vectors, 0);
    for (unsigned int subface_no = 0; subface_no < face->n_children();
         ++subface_no)
      {
        Assert(local_face_integrals.find(face->child(subface_no)) !=
                 local_face_integrals.end(),
               ExcInternalError());
        Assert(local_face_integrals[face->child(subface_no)][0] >= 0,
               ExcInternalError());

        for (unsigned int n = 0; n < n_solution_vectors; ++n)
          sum[n] += local_face_integrals[face->child(subface_no)][n];
      }

    local_face_integrals[face] = sum;
  }


  /**
   * 计算单个单元的面的误差。
   * 这个函数只在二维或三维中需要。
   * 一维的误差估计器是单独实现的。
   *
   */
  template <typename InputVector, int dim, int spacedim>
  void
  estimate_one_cell(
    const typename DoFHandler<dim, spacedim>::active_cell_iterator &cell,
    ParallelData<dim, spacedim, typename InputVector::value_type>
      &                                     parallel_data,
    std::map<typename DoFHandler<dim, spacedim>::face_iterator,
             std::vector<double>> &         local_face_integrals,
    const std::vector<const InputVector *> &solutions,
    const typename KellyErrorEstimator<dim, spacedim>::Strategy strategy)
  {
    const unsigned int n_solution_vectors = solutions.size();

    const types::subdomain_id subdomain_id = parallel_data.subdomain_id;
    const unsigned int        material_id  = parallel_data.material_id;

    // empty our own copy of the local face integrals
    local_face_integrals.clear();

    // loop over all faces of this cell
    for (const unsigned int face_no : cell->face_indices())
      {
        const typename DoFHandler<dim, spacedim>::face_iterator face =
          cell->face(face_no);

        // make sure we do work only once: this face may either be regular
        // or irregular. if it is regular and has a neighbor, then we visit
        // the face twice, once from every side. let the one with the lower
        // index do the work. if it is at the boundary, or if the face is
        // irregular, then do the work below
        if ((face->has_children() == false) && !cell->at_boundary(face_no) &&
            (!cell->neighbor_is_coarser(face_no) &&
             (cell->neighbor(face_no)->index() < cell->index() ||
              (cell->neighbor(face_no)->index() == cell->index() &&
               cell->neighbor(face_no)->level() < cell->level()))))
          continue;

        // if the neighboring cell is less refined than the present one,
        // then do nothing since we integrate over the subfaces when we
        // visit the coarse cells.
        if (face->at_boundary() == false)
          if (cell->neighbor_is_coarser(face_no))
            continue;

        // if this face is part of the boundary but not of the neumann
        // boundary -> nothing to do. However, to make things easier when
        // summing up the contributions of the faces of cells, we enter this
        // face into the list of faces with contribution zero.
        if (face->at_boundary() &&
            (parallel_data.neumann_bc->find(face->boundary_id()) ==
             parallel_data.neumann_bc->end()))
          {
            local_face_integrals[face] =
              std::vector<double>(n_solution_vectors, 0.);
            continue;
          }

        // finally: note that we only have to do something if either the
        // present cell is on the subdomain we care for (and the same for
        // material_id), or if one of the neighbors behind the face is on
        // the subdomain we care for
        if (!(((subdomain_id == numbers::invalid_subdomain_id) ||
               (cell->subdomain_id() == subdomain_id)) &&
              ((material_id == numbers::invalid_material_id) ||
               (cell->material_id() == material_id))))
          {
            // ok, cell is unwanted, but maybe its neighbor behind the face
            // we presently work on? oh is there a face at all?
            if (face->at_boundary())
              continue;

            bool care_for_cell = false;
            if (face->has_children() == false)
              care_for_cell |=
                ((cell->neighbor(face_no)->subdomain_id() == subdomain_id) ||
                 (subdomain_id == numbers::invalid_subdomain_id)) &&
                ((cell->neighbor(face_no)->material_id() == material_id) ||
                 (material_id == numbers::invalid_material_id));
            else
              {
                for (unsigned int sf = 0; sf < face->n_children(); ++sf)
                  if (((cell->neighbor_child_on_subface(face_no, sf)
                          ->subdomain_id() == subdomain_id) &&
                       (material_id == numbers::invalid_material_id)) ||
                      ((cell->neighbor_child_on_subface(face_no, sf)
                          ->material_id() == material_id) &&
                       (subdomain_id == numbers::invalid_subdomain_id)))
                    {
                      care_for_cell = true;
                      break;
                    }
              }

            // so if none of the neighbors cares for this subdomain or
            // material either, then try next face
            if (care_for_cell == false)
              continue;
          }

        // so now we know that we care for this face, let's do something
        // about it. first re-size the arrays we may use to the correct
        // size:
        parallel_data.resize(cell->active_fe_index());


        // then do the actual integration
        if (face->has_children() == false)
          // if the face is a regular one, i.e.  either on the other side
          // there is nirvana (face is at boundary), or the other side's
          // refinement level is the same as that of this side, then handle
          // the integration of these both cases together
          integrate_over_regular_face(solutions,
                                      parallel_data,
                                      local_face_integrals,
                                      cell,
                                      face_no,
                                      parallel_data.fe_face_values_cell,
                                      parallel_data.fe_face_values_neighbor,
                                      strategy);

        else
          // otherwise we need to do some special computations which do not
          // fit into the framework of the above function
          integrate_over_irregular_face(solutions,
                                        parallel_data,
                                        local_face_integrals,
                                        cell,
                                        face_no,
                                        parallel_data.fe_face_values_cell,
                                        parallel_data.fe_subface_values,
                                        strategy);
      }
  }
} // namespace internal



// the following function is still independent of dimension, but it
// calls dimension dependent functions
template <int dim, int spacedim>
template <typename InputVector>
void
KellyErrorEstimator<dim, spacedim>::estimate(
  const Mapping<dim, spacedim> &   mapping,
  const DoFHandler<dim, spacedim> &dof_handler,
  const Quadrature<dim - 1> &      quadrature,
  const std::map<types::boundary_id,
                 const Function<spacedim, typename InputVector::value_type> *>
    &                       neumann_bc,
  const InputVector &       solution,
  Vector<float> &           error,
  const ComponentMask &     component_mask,
  const Function<spacedim> *coefficients,
  const unsigned int        n_threads,
  const types::subdomain_id subdomain_id,
  const types::material_id  material_id,
  const Strategy            strategy)
{
  // just pass on to the other function
  const std::vector<const InputVector *> solutions(1, &solution);
  std::vector<Vector<float> *>           errors(1, &error);
  estimate(mapping,
           dof_handler,
           quadrature,
           neumann_bc,
           solutions,
           errors,
           component_mask,
           coefficients,
           n_threads,
           subdomain_id,
           material_id,
           strategy);
}


template <int dim, int spacedim>
template <typename InputVector>
void
KellyErrorEstimator<dim, spacedim>::estimate(
  const DoFHandler<dim, spacedim> &dof_handler,
  const Quadrature<dim - 1> &      quadrature,
  const std::map<types::boundary_id,
                 const Function<spacedim, typename InputVector::value_type> *>
    &                       neumann_bc,
  const InputVector &       solution,
  Vector<float> &           error,
  const ComponentMask &     component_mask,
  const Function<spacedim> *coefficients,
  const unsigned int        n_threads,
  const types::subdomain_id subdomain_id,
  const types::material_id  material_id,
  const Strategy            strategy)
{
  estimate(get_default_linear_mapping(dof_handler.get_triangulation()),
           dof_handler,
           quadrature,
           neumann_bc,
           solution,
           error,
           component_mask,
           coefficients,
           n_threads,
           subdomain_id,
           material_id,
           strategy);
}


template <int dim, int spacedim>
template <typename InputVector>
void
KellyErrorEstimator<dim, spacedim>::estimate(
  const Mapping<dim, spacedim> &   mapping,
  const DoFHandler<dim, spacedim> &dof_handler,
  const hp::QCollection<dim - 1> & quadrature,
  const std::map<types::boundary_id,
                 const Function<spacedim, typename InputVector::value_type> *>
    &                       neumann_bc,
  const InputVector &       solution,
  Vector<float> &           error,
  const ComponentMask &     component_mask,
  const Function<spacedim> *coefficients,
  const unsigned int        n_threads,
  const types::subdomain_id subdomain_id,
  const types::material_id  material_id,
  const Strategy            strategy)
{
  // just pass on to the other function
  const std::vector<const InputVector *> solutions(1, &solution);
  std::vector<Vector<float> *>           errors(1, &error);
  estimate(mapping,
           dof_handler,
           quadrature,
           neumann_bc,
           solutions,
           errors,
           component_mask,
           coefficients,
           n_threads,
           subdomain_id,
           material_id,
           strategy);
}


template <int dim, int spacedim>
template <typename InputVector>
void
KellyErrorEstimator<dim, spacedim>::estimate(
  const DoFHandler<dim, spacedim> &dof_handler,
  const hp::QCollection<dim - 1> & quadrature,
  const std::map<types::boundary_id,
                 const Function<spacedim, typename InputVector::value_type> *>
    &                       neumann_bc,
  const InputVector &       solution,
  Vector<float> &           error,
  const ComponentMask &     component_mask,
  const Function<spacedim> *coefficients,
  const unsigned int        n_threads,
  const types::subdomain_id subdomain_id,
  const types::material_id  material_id,
  const Strategy            strategy)
{
  estimate(get_default_linear_mapping(dof_handler.get_triangulation()),
           dof_handler,
           quadrature,
           neumann_bc,
           solution,
           error,
           component_mask,
           coefficients,
           n_threads,
           subdomain_id,
           material_id,
           strategy);
}



template <int dim, int spacedim>
template <typename InputVector>
void
KellyErrorEstimator<dim, spacedim>::estimate(
  const Mapping<dim, spacedim> &   mapping,
  const DoFHandler<dim, spacedim> &dof_handler,
  const hp::QCollection<dim - 1> & face_quadratures,
  const std::map<types::boundary_id,
                 const Function<spacedim, typename InputVector::value_type> *>
    &                                     neumann_bc,
  const std::vector<const InputVector *> &solutions,
  std::vector<Vector<float> *> &          errors,
  const ComponentMask &                   component_mask,
  const Function<spacedim> *              coefficients,
  const unsigned int,
  const types::subdomain_id subdomain_id_,
  const types::material_id  material_id,
  const Strategy            strategy)
{
#ifdef DEAL_II_WITH_P4EST
  if (dynamic_cast<const parallel::distributed::Triangulation<dim, spacedim> *>(
        &dof_handler.get_triangulation()) != nullptr)
    Assert((subdomain_id_ == numbers::invalid_subdomain_id) ||
             (subdomain_id_ ==
              dynamic_cast<
                const parallel::distributed::Triangulation<dim, spacedim> &>(
                dof_handler.get_triangulation())
                .locally_owned_subdomain()),
           ExcMessage(
             "For parallel distributed triangulations, the only "
             "valid subdomain_id that can be passed here is the "
             "one that corresponds to the locally owned subdomain id."));

  const types::subdomain_id subdomain_id =
    ((dynamic_cast<const parallel::distributed::Triangulation<dim, spacedim> *>(
        &dof_handler.get_triangulation()) != nullptr) ?
       dynamic_cast<const parallel::distributed::Triangulation<dim, spacedim>
                      &>(dof_handler.get_triangulation())
         .locally_owned_subdomain() :
       subdomain_id_);
#else
  const types::subdomain_id subdomain_id = subdomain_id_;
#endif

  const unsigned int n_components = dof_handler.get_fe(0).n_components();
  (void)n_components;

  // sanity checks
  Assert(solutions.size() > 0, ExcNoSolutions());
  Assert(solutions.size() == errors.size(),
         ExcIncompatibleNumberOfElements(solutions.size(), errors.size()));

  for (const auto &boundary_function : neumann_bc)
    {
      (void)boundary_function;
      Assert(boundary_function.second->n_components == n_components,
             ExcInvalidBoundaryFunction(boundary_function.first,
                                        boundary_function.second->n_components,
                                        n_components));
    }

  Assert(component_mask.represents_n_components(n_components),
         ExcInvalidComponentMask());
  Assert(component_mask.n_selected_components(n_components) > 0,
         ExcInvalidComponentMask());

  Assert((coefficients == nullptr) ||
           (coefficients->n_components == n_components) ||
           (coefficients->n_components == 1),
         ExcInvalidCoefficient());

  for (unsigned int n = 0; n < solutions.size(); ++n)
    Assert(solutions[n]->size() == dof_handler.n_dofs(),
           ExcDimensionMismatch(solutions[n]->size(), dof_handler.n_dofs()));

  const unsigned int n_solution_vectors = solutions.size();

  // Map of integrals indexed by the corresponding face. In this map we store
  // the integrated jump of the gradient for each face.  At the end of the
  // function, we again loop over the cells and collect the contributions of
  // the different faces of the cell.
  std::map<typename DoFHandler<dim, spacedim>::face_iterator,
           std::vector<double>>
    face_integrals;

  // all the data needed in the error estimator by each of the threads is
  // gathered in the following structures
  const hp::MappingCollection<dim, spacedim> mapping_collection(mapping);
  const internal::ParallelData<dim, spacedim, typename InputVector::value_type>
    parallel_data(dof_handler.get_fe_collection(),
                  face_quadratures,
                  mapping_collection,
                  (!neumann_bc.empty() || (coefficients != nullptr)),
                  solutions.size(),
                  subdomain_id,
                  material_id,
                  &neumann_bc,
                  component_mask,
                  coefficients);
  std::map<typename DoFHandler<dim, spacedim>::face_iterator,
           std::vector<double>>
    sample_local_face_integrals;

  // now let's work on all those cells:
  WorkStream::run(
    dof_handler.begin_active(),
    static_cast<typename DoFHandler<dim, spacedim>::active_cell_iterator>(
      dof_handler.end()),
    [&solutions, strategy](
      const typename DoFHandler<dim, spacedim>::active_cell_iterator &cell,
      internal::ParallelData<dim, spacedim, typename InputVector::value_type>
        &                            parallel_data,
      std::map<typename DoFHandler<dim, spacedim>::face_iterator,
               std::vector<double>> &local_face_integrals) {
      internal::estimate_one_cell(
        cell, parallel_data, local_face_integrals, solutions, strategy);
    },
    [&face_integrals](
      const std::map<typename DoFHandler<dim, spacedim>::face_iterator,
                     std::vector<double>> &local_face_integrals) {
      internal::copy_local_to_global<dim, spacedim>(local_face_integrals,
                                                    face_integrals);
    },
    parallel_data,
    sample_local_face_integrals);

  // finally add up the contributions of the faces for each cell

  // reserve one slot for each cell and set it to zero
  for (unsigned int n = 0; n < n_solution_vectors; ++n)
    {
      (*errors[n]).reinit(dof_handler.get_triangulation().n_active_cells());
      for (unsigned int i = 0;
           i < dof_handler.get_triangulation().n_active_cells();
           ++i)
        (*errors[n])(i) = 0;
    }

  // now walk over all cells and collect information from the faces. only do
  // something if this is a cell we care for based on the subdomain id
  for (const auto &cell : dof_handler.active_cell_iterators())
    if (((subdomain_id == numbers::invalid_subdomain_id) ||
         (cell->subdomain_id() == subdomain_id)) &&
        ((material_id == numbers::invalid_material_id) ||
         (cell->material_id() == material_id)))
      {
        const unsigned int present_cell = cell->active_cell_index();

        // loop over all faces of this cell
        for (const unsigned int face_no : cell->face_indices())
          {
            Assert(face_integrals.find(cell->face(face_no)) !=
                     face_integrals.end(),
                   ExcInternalError());
            const double factor = internal::cell_factor<dim, spacedim>(
              cell, face_no, dof_handler, strategy);

            for (unsigned int n = 0; n < n_solution_vectors; ++n)
              {
                // make sure that we have written a meaningful value into this
                // slot
                Assert(face_integrals[cell->face(face_no)][n] >= 0,
                       ExcInternalError());

                (*errors[n])(present_cell) +=
                  (face_integrals[cell->face(face_no)][n] * factor);
              }
          }

        for (unsigned int n = 0; n < n_solution_vectors; ++n)
          (*errors[n])(present_cell) = std::sqrt((*errors[n])(present_cell));
      }
}



template <int dim, int spacedim>
template <typename InputVector>
void
KellyErrorEstimator<dim, spacedim>::estimate(
  const Mapping<dim, spacedim> &   mapping,
  const DoFHandler<dim, spacedim> &dof_handler,
  const Quadrature<dim - 1> &      quadrature,
  const std::map<types::boundary_id,
                 const Function<spacedim, typename InputVector::value_type> *>
    &                                     neumann_bc,
  const std::vector<const InputVector *> &solutions,
  std::vector<Vector<float> *> &          errors,
  const ComponentMask &                   component_mask,
  const Function<spacedim> *              coefficients,
  const unsigned int                      n_threads,
  const types::subdomain_id               subdomain_id,
  const types::material_id                material_id,
  const Strategy                          strategy)
{
  // forward to the function with the QCollection
  estimate(mapping,
           dof_handler,
           hp::QCollection<dim - 1>(quadrature),
           neumann_bc,
           solutions,
           errors,
           component_mask,
           coefficients,
           n_threads,
           subdomain_id,
           material_id,
           strategy);
}


template <int dim, int spacedim>
template <typename InputVector>
void
KellyErrorEstimator<dim, spacedim>::estimate(
  const DoFHandler<dim, spacedim> &dof_handler,
  const Quadrature<dim - 1> &      quadrature,
  const std::map<types::boundary_id,
                 const Function<spacedim, typename InputVector::value_type> *>
    &                                     neumann_bc,
  const std::vector<const InputVector *> &solutions,
  std::vector<Vector<float> *> &          errors,
  const ComponentMask &                   component_mask,
  const Function<spacedim> *              coefficients,
  const unsigned int                      n_threads,
  const types::subdomain_id               subdomain_id,
  const types::material_id                material_id,
  const Strategy                          strategy)
{
  estimate(get_default_linear_mapping(dof_handler.get_triangulation()),
           dof_handler,
           quadrature,
           neumann_bc,
           solutions,
           errors,
           component_mask,
           coefficients,
           n_threads,
           subdomain_id,
           material_id,
           strategy);
}



template <int dim, int spacedim>
template <typename InputVector>
void
KellyErrorEstimator<dim, spacedim>::estimate(
  const DoFHandler<dim, spacedim> &dof_handler,
  const hp::QCollection<dim - 1> & quadrature,
  const std::map<types::boundary_id,
                 const Function<spacedim, typename InputVector::value_type> *>
    &                                     neumann_bc,
  const std::vector<const InputVector *> &solutions,
  std::vector<Vector<float> *> &          errors,
  const ComponentMask &                   component_mask,
  const Function<spacedim> *              coefficients,
  const unsigned int                      n_threads,
  const types::subdomain_id               subdomain_id,
  const types::material_id                material_id,
  const Strategy                          strategy)
{
  estimate(get_default_linear_mapping(dof_handler.get_triangulation()),
           dof_handler,
           quadrature,
           neumann_bc,
           solutions,
           errors,
           component_mask,
           coefficients,
           n_threads,
           subdomain_id,
           material_id,
           strategy);
}

DEAL_II_NAMESPACE_CLOSE

#endif


