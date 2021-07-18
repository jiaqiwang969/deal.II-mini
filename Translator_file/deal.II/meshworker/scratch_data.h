//include/deal.II-translator/meshworker/scratch_data_0.txt
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

#ifndef dealii_meshworker_scratch_data_h
#define dealii_meshworker_scratch_data_h

#include <deal.II/base/config.h>

#include <deal.II/algorithms/general_data_storage.h>

#include <deal.II/differentiation/ad.h>

#include <deal.II/fe/fe_interface_values.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>

#include <boost/any.hpp>

#include <algorithm>
#include <map>
#include <vector>

DEAL_II_NAMESPACE_OPEN

namespace MeshWorker
{
  /**
   * 一个帮助类，用于简化线性和非线性问题的并行装配，以及有限元场的评估。
   * 该类是ScratchData中的一个落款，可与 WorkStream::run()
   * 函数和 MeshWorker::mesh_loop 函数（）一起使用。
   * ScratchData类有三个主要目标。
   *
   *
   *
   *
   *
   * - 为当前单元及其相邻单元按要求创建FEValues、FEFaceValues、FESSubfaceValues和FEInterfaceValues*（仅当用户提供的算法需要它们时），并在组装单元、面或子面贡献时提供统一接口访问FEValues对象。
   *
   *
   *
   *
   *
   * - 存储任意数据类型（或对任意数据类型的引用），可以通过名称检索（例如，当装配需要引用其他对象，如以前的时间步骤的解决方案向量、以前的非线性迭代向量、几何描述等），在一个类型的对象中。
   *
   *
   *
   *
   *
   *
   * - 为那些用户可能需要在正交点上的临时数据向量的用例提供一个合理的接口，允许按需构建这些临时向量*，并可方便地访问已经计算出的解向量的值、梯度等。    在 "在当前单元上工作的方法 "一节中的方法按要求初始化当前单元上的内部FEValues、FEFaceValues、FESSubfaceValues和FEInterfaceValues对象，允许使用这个类作为四个不同对象的单一替代品，用于整合和查询单元、面和子面的有限元值。    同样地，"在邻接单元上工作的方法 "一节中的方法也会根据需要初始化邻接单元上的（不同的）内部FEValues、FEFaceValues和FESSubfaceValues，允许使用这个类来替代通常需要在邻接单元、其面和子面上集成的另外三个对象（例如，在非连续Galerkin方法中）。    如果你需要在刚刚用 "在当前单元上工作的方法 "一节中的某个函数初始化的单元、面或子面上检索有限元解向量的值或梯度，你可以使用 "在当前单元上评估有限元场及其导数 "一节中的方法。    下面的代码片断给出了这个类的一个使用实例。
   * @code
   * ScratchData<dim,spacedim> scratch(fe,
   *                                 cell_quadrature,
   *                                 update_values | update_gradients);
   *
   * FEValuesExtractors::Vector velocity(0);
   * FEValuesExtractors::Scalar pressure(dim);
   *
   * ...
   *
   * for(const auto &cell : dof_handler.active_cell_iterators())
   * {
   *  const auto &fe_values = scratch.reinit(cell);
   *  scratch.extract_local_dof_values("current_solution",
   *                                   current_solution);
   *
   *  scratch.extract_local_dof_values("previous_solution",
   *                                   previous_solution);
   *
   *  const auto &local_solution_values =
   *        scratch.get_local_dof_values("current_solution",
   *                                     current_solution);
   *
   *  const auto &local_previous_solution_values =
   *        scratch.get_local_dof_values("previous_solution",
   *                                     previous_solution);
   *
   *  const auto &previous_solution_values_velocity =
   *        scratch.get_values("previous_solution", velocity);
   *  const auto &previous_solution_values_pressure =
   *        scratch.get_values("previous_solution", pressure);
   *
   *  const auto &solution_values_velocity =
   *        scratch.get_values("solution", velocity);
   *  const auto &solution_values_pressure =
   *        scratch.get_values("solution", pressure);
   *
   *  const auto &solution_symmetric_gradients_velocity =
   *        scratch.get_symmetric_gradients("solution", velocity);
   *
   *  // Do something with the above.
   * }
   * @endcode
   * 你调用这个类的函数的顺序很重要：如果你调用
   * ScratchData::reinit()
   * 函数，该函数需要一个活动单元的迭代器，那么随后调用内部需要FEValuesBase对象的方法将使用以给定单元初始化的内部FEValues对象来执行其计算。另一方面，如果你调用了
   * ScratchData::reinit()
   * 方法，该方法也需要一个面的索引，那么随后所有对需要FEValuesBase对象的方法的调用，将使用内部存储的FEFaceValues对象，用传递给
   * ScratchData::reinit()
   * 函数的单元格和面的索引初始化。这同样适用于
   * ScratchData::reinit()
   * 方法，它需要三个参数：单元格、面的索引和子面的索引。
   * 用户代码的结构应该是不交错进行单元格的工作和面的工作。
   * 例如，考虑下面的代码片段。
   * @code
   * ScratchData<dim,spacedim> scratch(...);
   * FEValuesExtractors::Scalar temperature(0);
   *
   * for(const auto &cell : dof_handler.active_cell_iterators())
   * {
   * const auto &fe_values = scratch.reinit(cell);
   * const auto &local_dof_values =
   *       scratch.extract_local_dof_values("solution", solution);
   *
   * // This will contain all values of the temperature on the cell
   * // quadrature points
   * const auto &solution_values_cell =
   *       scratch.get_values("solution", temperature);
   *
   * // Do something with values on the cell
   * ...
   *
   * // Now start working on the faces
   * for(unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
   * {
   *    const auto &fe_face_values = scratch.reinit(cell, f);
   *
   *    // Now we'll have access to the values of the temperature on the faces
   *    const auto &solution_values_face =
   *          scratch.get_values("solution", temperature);
   *
   *    // Notice how the function call is the same, but the result depends
   *    // on what was the last reinit() function that was called. In this
   *    // case, we called reinit(cell, f), triggering internal work on
   *    // faces.
   * }
   *
   * // At this point, we would like to go back and work on cells,
   * // for example querying the values of the gradients:
   * const auto &solution_gradients =
   *       scratch.get_gradients("solution", temperature);
   *
   * // This assertion would be triggered in debug mode
   * AssertDimension(solution_gradients.size(), quadrature_cell.size());
   *
   * // However, with the above call the content of solution_gradients is the
   * // gradient of the temperature on the quadrature points of the last face
   * // visited in the previous loop.
   * // If you really want to have the values of the gradients on the cell,
   * // at this point your only option is to call
   * scratch.reinit(cell);
   *
   * // (again) before querying for the gradients:
   * const auto &solution_gradients_cell =
   *       scratch.get_gradients("solution", temperature);
   *
   * // The call to reinit(cell), however, is an expensive one. You
   * // should make sure you only call it once per cell by grouping together
   * // all queries that concern cells in the same place (right after you
   * // call the reinit(cell) method).
   * // A similar argument holds for all calls on each face and subface.
   * }
   * @endcode
   * 当使用这个类时，请引用
   * @code{.bib}
   * @article{SartoriGiulianiBardelloni-2018-a,
   * Author  = {Sartori, Alberto and Giuliani, Nicola and
   *          Bardelloni, Mauro and Heltai, Luca},
   * Journal = {SoftwareX},
   * Pages   = {318--327},
   * Title   = {{deal2lkit: A toolkit library for high performance
   *            programming in deal.II}},
   * Doi     = {10.1016/j.softx.2018.09.004},
   * Volume  = {7},
   * Year    = {2018}}
   * @endcode
   *
   *
   */
  template <int dim, int spacedim = dim>
  class ScratchData
  {
  public:
    /**
     * 创建一个空的ScratchData对象。一个指向 @p mapping 和 @p fe
     * 的SmartPointer被内部存储。请确保它们的寿命比这个类实例长。
     * 构造函数不初始化任何内部的FEValues对象。
     * 这些对象会在第一次调用 reinit()
     * 函数时被初始化，使用这里传递的参数。          @param
     * mapping 在内部FEValues对象中使用的映射  @param  fe 有限元
     * @param  quadrature 单元正交  @param  update_flags
     * 当前单元FEValues和邻近单元FEValues的更新标志  @param  ]
     * face_quadrature
     * 面部正交，用于当前单元和邻居单元的FEFaceValues和FESubfaceValues
     * @param  face_update_flags
     * 用于当前单元和邻居单元的FEFaceValues和FESubfaceValues的更新标志
     *
     */
    ScratchData(
      const Mapping<dim, spacedim> &      mapping,
      const FiniteElement<dim, spacedim> &fe,
      const Quadrature<dim> &             quadrature,
      const UpdateFlags &                 update_flags,
      const Quadrature<dim - 1> &face_quadrature   = Quadrature<dim - 1>(),
      const UpdateFlags &        face_update_flags = update_default);

    /**
     * 与其他构造函数类似，但这个构造函数允许为相邻单元和面指定不同的标志。
     * @param  mapping 在内部FEValues对象中使用的映射  @param  fe
     * 有限元素  @param  quadrature 单元正交  @param  update_flags
     * 当前单元FEValues的更新标志  @param  neighbor_update_flags
     * 邻居单元FEValues的更新标志  @param  ] face_quadrature
     * 面部正交，用于当前单元和邻居单元的FEFaceValues和FESubfaceValues
     * @param  face_update_flags
     * 用于当前单元的FEFaceValues和FESubfaceValues的更新标记
     * @param  neighbor_face_update_flags
     * 用于邻居单元的FEFaceValues和FESubfaceValues的更新标记
     *
     */
    ScratchData(
      const Mapping<dim, spacedim> &      mapping,
      const FiniteElement<dim, spacedim> &fe,
      const Quadrature<dim> &             quadrature,
      const UpdateFlags &                 update_flags,
      const UpdateFlags &                 neighbor_update_flags,
      const Quadrature<dim - 1> &face_quadrature   = Quadrature<dim - 1>(),
      const UpdateFlags &        face_update_flags = update_default,
      const UpdateFlags &        neighbor_face_update_flags = update_default);

    /**
     * 与其他构造函数相同，使用默认的MappingQ1。
     * @param  fe 有限元  @param  quadrature 单元正交  @param
     * update_flags 当前单元FEValues和邻近单元FEValues的更新标志
     * @param  ] face_quadrature
     * 面部正交，用于当前单元和邻居单元的FEFaceValues和FESubfaceValues
     * @param  face_update_flags
     * 用于当前单元和邻居单元的FEFaceValues和FESubfaceValues的更新标志
     *
     */
    ScratchData(
      const FiniteElement<dim, spacedim> &fe,
      const Quadrature<dim> &             quadrature,
      const UpdateFlags &                 update_flags,
      const Quadrature<dim - 1> &face_quadrature   = Quadrature<dim - 1>(),
      const UpdateFlags &        face_update_flags = update_default);

    /**
     * 与其他构造函数相同，使用默认的MappingQ1。
     * @param  fe 有限元  @param  quadrature 单元正交  @param
     * update_flags 当前单元FEValues的更新标志  @param
     * neighbor_update_flags 邻近单元FEValues的更新标志  @param  ]
     * face_quadrature
     * 面部正交，用于当前单元和邻居单元的FEFaceValues和FESubfaceValues
     * @param  face_update_flags
     * 用于当前单元的FEFaceValues和FESubfaceValues的更新标记
     * @param  neighbor_face_update_flags
     * 用于邻居单元的FEFaceValues和FESubfaceValues的更新标记
     *
     */
    ScratchData(
      const FiniteElement<dim, spacedim> &fe,
      const Quadrature<dim> &             quadrature,
      const UpdateFlags &                 update_flags,
      const UpdateFlags &                 neighbor_update_flags,
      const Quadrature<dim - 1> &face_quadrature   = Quadrature<dim - 1>(),
      const UpdateFlags &        face_update_flags = update_default,
      const UpdateFlags &        neighbor_face_update_flags = update_default);

    /**
     * 深度复制构造函数。FEValues对象不会被复制。
     *
     */
    ScratchData(const ScratchData<dim, spacedim> &scratch);

    /**
     * @name  对当前单元工作的方法
     *
     */
     /**@{*/  // CurrentCellMethods

    /**
     * 用给定的 @p cell,
     * 初始化内部的FEValues，并返回它的引用。
     * 调用此函数后，get_current_fe_values()将返回此方法的同一对象，作为FEValuesBase引用。
     *
     */
    const FEValues<dim, spacedim> &
    reinit(
      const typename DoFHandler<dim, spacedim>::active_cell_iterator &cell);

    /**
     * 初始化内部的FEFaceValues，在给定的 @p face_no
     * 上使用，并返回一个引用。
     * 调用此函数后，get_current_fe_values()将返回此方法的同一对象，作为FEValuesBase引用。
     *
     */
    const FEFaceValues<dim, spacedim> &
    reinit(const typename DoFHandler<dim, spacedim>::active_cell_iterator &cell,
           const unsigned int face_no);

    /**
     * 初始化内部的FESubfaceValues，使其在给定的 @p subface_no,
     * 上使用 @p face_no, ，并返回一个引用。
     * 调用此函数后，get_current_fe_values()将返回此方法的同一对象，作为FEValuesBase的引用。
     * 如果 @p subface_no 是 numbers::invalid_unsigned_int,
     * ，则会调用只接收 @p cell 和 @p face_no 的reinit()函数。
     *
     */
    const FEFaceValuesBase<dim, spacedim> &
    reinit(const typename DoFHandler<dim, spacedim>::active_cell_iterator &cell,
           const unsigned int face_no,
           const unsigned int subface_no);

    /**
     * 用给定的参数初始化内部的FEInterfaceValues，并返回它的一个引用。
     * 调用此函数后，get_local_dof_indices(), get_quadrature_points(),
     * get_normal_vectors(), and
     * get_JxW_values()将被转发给本地FEInterfaceValues对象。方法get_current_fe_values()将返回与当前单元相关的FEValuesBase，而get_neighbor_fe_values()将与邻居单元相关。get_local_dof_indices()方法将返回与
     * FEInterfaceValues::get_interface_dof_to_dof_indices(),
     * 相同的结果，而get_neighbor_dof_indices()将返回邻居单元的本地dof指数。
     *
     */
    const FEInterfaceValues<dim, spacedim> &
    reinit(const typename DoFHandler<dim, spacedim>::active_cell_iterator &cell,
           const unsigned int face_no,
           const unsigned int sub_face_no,
           const typename DoFHandler<dim, spacedim>::active_cell_iterator
             &                cell_neighbor,
           const unsigned int face_no_neighbor,
           const unsigned int sub_face_no_neighbor);

    /**
     * 获取当前初始化的FEValues。
     * 如果reinit(cell)函数最后被调用，该函数将返回内部FEValues。如果调用了reinit(cell,
     * face_no)函数，则此函数返回内部的FEFaceValues，如果调用了reinit(cell,
     * face_no, subface_no)函数（有一个有效的 @p subface_no
     * 参数），则返回内部的FESubfaceValues对象。
     *
     */
    const FEValuesBase<dim, spacedim> &
    get_current_fe_values() const;

    /**
     * 返回内部FEValues对象的正交点。
     *
     */
    const std::vector<Point<spacedim>> &
    get_quadrature_points() const;

    /**
     * 返回内部FEValues对象的JxW值。
     *
     */
    const std::vector<double> &
    get_JxW_values() const;

    /**
     * 返回最后计算的法向量。
     *
     */
    const std::vector<Tensor<1, spacedim>> &
    get_normal_vectors() const;

    /**
     * 返回上次调用reinit()函数时传递的单元格的本地dof指数。
     *
     */
    const std::vector<types::global_dof_index> &
    get_local_dof_indices() const;

     /** @} */  // CurrentCellMethods

    /**
     * @name  对邻居单元工作的方法
     *
     */
     /** 
     * @{ */  // NeighborCellMethods

    /**
     * 初始化内部邻居FEValues，使用给定的 @p cell,
     * ，并返回对它的引用。
     * 调用此函数后，get_current_neighbor_fe_values()将返回此方法的同一对象，作为FEValuesBase引用。
     *
     */
    const FEValues<dim, spacedim> &
    reinit_neighbor(
      const typename DoFHandler<dim, spacedim>::active_cell_iterator &cell);

    /**
     * 初始化内部的FEFaceValues，在给定的 @p face_no
     * 上使用，并返回一个引用。
     * 调用此函数后，get_current_neighbor_fe_values()将返回此方法的同一对象，作为FEValuesBase引用。
     *
     */
    const FEFaceValues<dim, spacedim> &
    reinit_neighbor(
      const typename DoFHandler<dim, spacedim>::active_cell_iterator &cell,
      const unsigned int                                              face_no);

    /**
     * 初始化内部的FESubfaceValues，使其在给定的 @p subface_no,
     * 上使用 @p face_no, ，并返回一个引用。
     * 调用此函数后，get_current_neighbor_fe_values()将返回此方法的同一对象，作为FEValuesBase引用。
     * 如果 @p subface_no 是 numbers::invalid_unsigned_int,
     * ，则会调用只接收 @p cell 和 @p face_no 的reinit()函数。
     *
     */
    const FEFaceValuesBase<dim, spacedim> &
    reinit_neighbor(
      const typename DoFHandler<dim, spacedim>::active_cell_iterator &cell,
      const unsigned int                                              face_no,
      const unsigned int subface_no);

    /**
     * 获取当前初始化的邻居FEValues。
     * 如果最后一次调用 reinit_neighbor(cell)
     * 函数，该函数将返回邻居FEValues。如果调用了reinit_neighbor(cell,
     * face_no)函数，则此函数返回内部邻居FEFaceValues，如果调用了reinit_neighbor(cell,
     * face_no, subface_no)函数（有一个有效的 @p subface_no
     * 参数），则返回内部邻居FESubfaceValues对象。
     *
     */
    const FEValuesBase<dim, spacedim> &
    get_current_neighbor_fe_values() const;

    /**
     * 返回邻居FEValues对象的JxW值。
     *
     */
    const std::vector<double> &
    get_neighbor_JxW_values() const;

    /**
     * 返回邻居上最后计算的法向量。
     *
     */
    const std::vector<Tensor<1, spacedim>> &
    get_neighbor_normal_vectors();

    /**
     * 返回上次调用reinit_neighbor()函数时传递的邻居的局部道夫指数。
     *
     */
    const std::vector<types::global_dof_index> &
    get_neighbor_dof_indices() const;

     /** @} */  // NeighborCellMethods

    /**
     * 返回一个GeneralDataStorage对象，该对象可用于存储任何数量、任何类型的数据，然后通过一个标识符字符串进行访问。
     *
     */
    GeneralDataStorage &
    get_general_data_storage();

    /**
     * 返回一个GeneralDataStorage对象，该对象可用于存储任何类型的任何数量的数据，然后通过一个标识符字符串进行访问。
     *
     */
    const GeneralDataStorage &
    get_general_data_storage() const;

    /**
     * @name  评估当前单元上的有限元场及其导数
     *
     */
     /** 
     * @{ */  // CurrentCellEvaluation

    /**
     * 提取与内部初始化单元相关的局部道夫值。
     * 在调用这个函数之前，你必须确保你之前已经调用了reinit()函数中的一个。
     * 在每次调用这个函数时，都会产生一个新的dof值向量并在内部存储，除非找到一个先前的同名向量。如果是这样的话，那个向量的内容就会被覆盖掉。
     * 如果你给出一个唯一的 @p global_vector_name,
     * ，那么对于每个单元格，你可以保证得到一个与虚拟变量相同类型的独立的道夫向量。如果你使用一个自动微分的数字类型（如
     * Sacado::Fad::DFad<double>,
     * Sacado::Fad::DFad<Sacado::Fad::DFad<double>>,
     * 等），这个方法也会在内部初始化独立变量，允许你进行自动微分。
     * 你可以通过调用get_local_dof_values()方法来访问提取的局部道夫值，参数与你在这里传递的
     * @p global_vector_name 相同。
     * 注意，使用这种初始化策略会使这个ScratchData对象的使用与AD辅助类不兼容（因为它们会做自己的数据管理）。特别是，用户有必要自己管理所有的AD数据（包括在这个调用之前和之后）。
     *
     */
    template <typename VectorType, typename Number = double>
    void
    extract_local_dof_values(const std::string &global_vector_name,
                             const VectorType & input_vector,
                             const Number       dummy = Number(0));

    /**
     * 调用extract_local_dof_values()后，你可以通过这个方法检索存储的信息。
     * 参数 @p global_vector_name 和 @p dummy
     * 变量的类型都应该与你传递给extract_local_dof_values()函数的一致。
     *
     */
    template <typename Number = double>
    const std::vector<Number> &
    get_local_dof_values(const std::string &global_vector_name,
                         Number             dummy = Number(0)) const;

    /**
     * 对于由 @p global_vector_name,
     * 确定的解向量，计算正交点的函数值，并返回一个由你作为
     * @p variable 参数传递的 Extractor
     * 推断出的正确类型的向量。
     * 在你调用这个方法之前，你需要至少调用一次
     * extract_local_dof_values() 方法，传递相同的  @p
     * global_vector_name  字符串，以及相同类型的变量  @p dummy.
     * 如果你之前没有调用 extract_local_dof_values()
     * 方法，这个函数会抛出一个异常。
     * 为了使这个函数正常工作，你所调用的reinit()函数中的底层FEValues、FEFaceValues或FESubfaceValues对象必须已经计算了你所请求的信息。要做到这一点，update_values标志必须是你传递给这个对象的构造函数的UpdateFlags列表中的一个元素。
     * 参见FEValues类文档中的
     * "UpdateFlags、Mapping和FiniteElement的相互作用"，以了解更多信息。
     *
     */
    template <typename Extractor, typename Number = double>
    const std::vector<typename FEValuesViews::View<dim, spacedim, Extractor>::
                        template solution_value_type<Number>> &
    get_values(const std::string &global_vector_name,
               const Extractor &  variable,
               const Number       dummy = Number(0));

    /**
     * 对于 @p global_vector_name,
     * 确定的解向量，计算正交点的函数梯度，并返回一个由你作为
     * @p variable
     * 参数传递的提取器推导出的正确类型的向量。
     * 在你调用这个方法之前，你需要至少调用一次
     * extract_local_dof_values() 方法，传递相同的  @p
     * global_vector_name  字符串，以及相同类型的变量  @p dummy.
     * 如果你之前没有调用 extract_local_dof_values()
     * 方法，这个函数会抛出一个异常。
     * 为了使这个函数正常工作，你所调用的reinit()函数中的底层FEValues、FEFaceValues或FESubfaceValues对象必须已经计算了你所请求的信息。要做到这一点，update_gradients标志必须是你传递给这个对象的构造函数的UpdateFlags列表中的一个元素。
     * 参见FEValues类文档中的
     * "UpdateFlags、Mapping和FiniteElement的相互作用"，以了解更多信息。
     *
     */
    template <typename Extractor, typename Number = double>
    const std::vector<typename FEValuesViews::View<dim, spacedim, Extractor>::
                        template solution_gradient_type<Number>> &
    get_gradients(const std::string &global_vector_name,
                  const Extractor &  variable,
                  const Number       dummy = Number(0));

    /**
     * 对于 @p global_vector_name,
     * 确定的解向量，计算函数在正交点的对称梯度，并返回一个由你作为
     * @p variable
     * 参数传递的提取器推导出的正确类型的向量。
     * 在你调用这个方法之前，你需要至少调用一次
     * extract_local_dof_values() 方法，传递相同的  @p
     * global_vector_name  字符串，以及相同类型的变量  @p dummy.
     * 如果你之前没有调用 extract_local_dof_values()
     * 方法，这个函数会抛出一个异常。
     * 为了使这个函数正常工作，你所调用的reinit()函数中的底层FEValues、FEFaceValues或FESubfaceValues对象必须已经计算了你所请求的信息。要做到这一点，update_gradients标志必须是你传递给这个对象的构造函数的UpdateFlags列表中的一个元素。
     * 参见FEValues类文档中的
     * "UpdateFlags、Mapping和FiniteElement的相互作用
     * "以获得更多信息。
     *
     */
    template <typename Extractor, typename Number = double>
    const std::vector<typename FEValuesViews::View<dim, spacedim, Extractor>::
                        template solution_symmetric_gradient_type<Number>> &
    get_symmetric_gradients(const std::string &global_vector_name,
                            const Extractor &  variable,
                            const Number       dummy = Number(0));

    /**
     * 对于由 @p global_vector_name,
     * 确定的解向量，计算函数在正交点的发散，并返回一个由你作为
     * @p variable
     * 参数传递的提取器推导出的正确类型的向量。
     * 在你调用这个方法之前，你需要至少调用一次
     * extract_local_dof_values() 方法，传递相同的  @p
     * global_vector_name  字符串，以及相同类型的变量  @p dummy.
     * 如果你之前没有调用 extract_local_dof_values()
     * 方法，这个函数会抛出一个异常。
     * 为了使这个函数正常工作，你所调用的reinit()函数中的底层FEValues、FEFaceValues或FESubfaceValues对象必须已经计算了你所请求的信息。要做到这一点，update_gradients标志必须是你传递给这个对象的构造函数的UpdateFlags列表中的一个元素。
     * 参见FEValues类文档中的
     * "UpdateFlags、Mapping和FiniteElement的相互作用"，以了解更多信息。
     *
     */
    template <typename Extractor, typename Number = double>
    const std::vector<typename FEValuesViews::View<dim, spacedim, Extractor>::
                        template solution_divergence_type<Number>> &
    get_divergences(const std::string &global_vector_name,
                    const Extractor &  variable,
                    const Number       dummy = Number(0));

    /**
     * 对于由 @p global_vector_name,
     * 确定的解向量，计算正交点的函数卷曲，并返回一个由你作为
     * @p variable
     * 参数传递的提取器推导出的正确类型的向量。
     * 在你调用这个方法之前，你需要至少调用一次
     * extract_local_dof_values() 方法，传递相同的  @p
     * global_vector_name  字符串，以及相同类型的变量  @p dummy.
     * 如果你之前没有调用 extract_local_dof_values()
     * 方法，这个函数将抛出一个异常。
     * 为了使这个函数正常工作，你所调用的reinit()函数中的底层FEValues、FEFaceValues或FESubfaceValues对象必须已经计算了你所请求的信息。要做到这一点，update_gradients标志必须是你传递给这个对象的构造函数的UpdateFlags列表中的一个元素。
     * 参见FEValues类文档中的
     * "UpdateFlags、Mapping和FiniteElement的相互作用"，以了解更多信息。
     *
     */
    template <typename Extractor, typename Number = double>
    const std::vector<typename FEValuesViews::View<dim, spacedim, Extractor>::
                        template solution_curl_type<Number>> &
    get_curls(const std::string &global_vector_name,
              const Extractor &  variable,
              const Number       dummy = Number(0));

    /**
     * 对于由 @p global_vector_name,
     * 确定的解向量，计算正交点上的函数的
     * hessians，并返回一个由你作为 @p variable 参数传递的
     * Extractor 推断的正确类型的向量。
     * 在你调用这个方法之前，你需要至少调用一次
     * extract_local_dof_values() 方法，传递相同的  @p
     * global_vector_name  字符串，以及相同类型的变量  @p dummy.
     * 如果你之前没有调用 extract_local_dof_values()
     * 方法，这个函数会抛出一个异常。
     * 为了使这个函数正常工作，你所调用的reinit()函数中的底层FEValues、FEFaceValues或FESubfaceValues对象必须已经计算了你所请求的信息。要做到这一点，update_hessians标志必须是你传递给这个对象的构造函数的UpdateFlags列表中的一个元素。
     * 参见FEValues类文档中的
     * "UpdateFlags、Mapping和FiniteElement的相互作用
     * "以获得更多信息。
     *
     */
    template <typename Extractor, typename Number = double>
    const std::vector<typename FEValuesViews::View<dim, spacedim, Extractor>::
                        template solution_hessian_type<Number>> &
    get_hessians(const std::string &global_vector_name,
                 const Extractor &  variable,
                 const Number       dummy = Number(0));

    /**
     * 对于由 @p global_vector_name,
     * 确定的解向量，计算正交点上的函数的拉普拉斯，并返回一个由你作为
     * @p variable
     * 参数传递的提取器推导出的正确类型的向量。
     * 在你调用这个方法之前，你需要至少调用一次
     * extract_local_dof_values() 方法，传递相同的  @p
     * global_vector_name  字符串，以及相同类型的变量  @p dummy.
     * 如果你之前没有调用 extract_local_dof_values()
     * 方法，这个函数将抛出一个异常。
     * 为了使这个函数正常工作，你所调用的reinit()函数中的底层FEValues、FEFaceValues或FESubfaceValues对象必须已经计算了你所请求的信息。要做到这一点，update_hessians标志必须是你传递给这个对象的构造函数的UpdateFlags列表中的一个元素。
     * 参见FEValues类文档中的
     * "UpdateFlags、Mapping和FiniteElement的相互作用
     * "以获得更多信息。
     *
     */
    template <typename Extractor, typename Number = double>
    const std::vector<typename FEValuesViews::View<dim, spacedim, Extractor>::
                        template solution_laplacian_type<Number>> &
    get_laplacians(const std::string &global_vector_name,
                   const Extractor &  variable,
                   const Number       dummy = Number(0));

    /**
     * 对于由 @p global_vector_name,
     * 确定的解向量，计算正交点的函数的三次方，并返回一个由你作为
     * @p variable
     * 参数传递的提取器推导出的正确类型的向量。
     * 在你调用这个方法之前，你需要至少调用一次
     * extract_local_dof_values() 方法，传递相同的  @p
     * global_vector_name  字符串，以及相同类型的变量  @p dummy.
     * 如果你之前没有调用 extract_local_dof_values()
     * 方法，这个函数会抛出一个异常。
     * 为了使这个函数正常工作，你所调用的reinit()函数中的底层FEValues、FEFaceValues或FESubfaceValues对象必须已经计算了你所请求的信息。要做到这一点，update_3rd_derivatives标志必须是你传递给这个对象的构造函数的UpdateFlags列表中的一个元素。参见FEValues文档中的
     * "UpdateFlags、Mapping和FiniteElement的相互作用
     * "以获得更多信息。
     *
     */
    template <typename Extractor, typename Number = double>
    const std::vector<typename FEValuesViews::View<dim, spacedim, Extractor>::
                        template solution_third_derivative_type<Number>> &
    get_third_derivatives(const std::string &global_vector_name,
                          const Extractor &  variable,
                          const Number       dummy = Number(0));

     /** @} */  // CurrentCellEvaluation

    /**
     * 返回一个对所用映射的引用。
     *
     */
    const Mapping<dim, spacedim> &
    get_mapping() const;

  private:
    /**
     * 构建一个唯一的名称，以便在内部GeneralDataStorage对象中存储数值、梯度、发散等的向量。
     *
     */
    template <typename Extractor, typename Number = double>
    std::string
    get_unique_name(const std::string &global_vector_name,
                    const Extractor &  variable,
                    const std::string &object_type,
                    const unsigned int size,
                    const Number &     exemplar_number) const;

    /**
     * 构建一个唯一的名称来存储本地dof值。
     *
     */
    template <typename Number = double>
    std::string
    get_unique_dofs_name(const std::string &global_vector_name,
                         const unsigned int size,
                         const Number &     exemplar_number) const;

    /**
     * 内部FEValues所使用的映射。确保它的寿命比这个类长。
     *
     */
    SmartPointer<const Mapping<dim, spacedim>> mapping;

    /**
     * 内部FEValues使用的有限元。确保它的寿命比这个类长。
     *
     */
    SmartPointer<const FiniteElement<dim, spacedim>> fe;

    /**
     * 用于对当前单元和其相邻单元进行积分的正交公式。
     *
     */
    Quadrature<dim> cell_quadrature;

    /**
     * 用于在面、子面、相邻面和子面进行积分的正交公式。
     *
     */
    Quadrature<dim - 1> face_quadrature;

    /**
     * 在初始化单元格FEValues对象时使用的UpdateFlags。
     *
     */
    UpdateFlags cell_update_flags;

    /**
     * 在初始化相邻的单元格FEValues对象时使用的UpdateFlags。
     *
     */
    UpdateFlags neighbor_cell_update_flags;

    /**
     * 在初始化FEFaceValues和FESubfaceValues对象时使用的UpdateFlags。
     *
     */
    UpdateFlags face_update_flags;

    /**
     * 在初始化邻居FEFaceValues和FESubfaceValues对象时使用的UpdateFlags。
     *
     */
    UpdateFlags neighbor_face_update_flags;

    /**
     * 当前单元上的有限元值。
     *
     */
    std::unique_ptr<FEValues<dim, spacedim>> fe_values;

    /**
     * 当前面的有限元值。
     *
     */
    std::unique_ptr<FEFaceValues<dim, spacedim>> fe_face_values;

    /**
     * 当前子面的有限元值。
     *
     */
    std::unique_ptr<FESubfaceValues<dim, spacedim>> fe_subface_values;

    /**
     * 邻近单元上的有限元值。
     *
     */
    std::unique_ptr<FEValues<dim, spacedim>> neighbor_fe_values;

    /**
     * 邻近面的有限元值。
     *
     */
    std::unique_ptr<FEFaceValues<dim, spacedim>> neighbor_fe_face_values;

    /**
     * 邻近子面的有限元值。
     *
     */
    std::unique_ptr<FESubfaceValues<dim, spacedim>> neighbor_fe_subface_values;

    /**
     * 面上的界面值。
     *
     */
    std::unique_ptr<FEInterfaceValues<dim, spacedim>> interface_fe_values;

    /**
     * 当前单元上的Dof指数。
     *
     */
    std::vector<types::global_dof_index> local_dof_indices;

    /**
     * 邻近单元上的阻值指数。
     *
     */
    std::vector<types::global_dof_index> neighbor_dof_indices;

    /**
     * 用户数据存储。
     *
     */
    GeneralDataStorage user_data_storage;

    /**
     * 内部数据存储。
     *
     */
    GeneralDataStorage internal_data_storage;

    /**
     * 指向该单元上最后使用的FEValues/FEFaceValues或FESubfaceValues对象的指针。
     *
     */
    SmartPointer<const FEValuesBase<dim, spacedim>> current_fe_values;

    /**
     * 指向邻近单元上最后使用的FEValues/FEFaceValues或FESubfaceValues对象的指针。
     *
     */
    SmartPointer<const FEValuesBase<dim, spacedim>> current_neighbor_fe_values;
  };

#ifndef DOXYGEN
  template <int dim, int spacedim>
  template <typename Extractor, typename Number>
  std::string
  ScratchData<dim, spacedim>::get_unique_name(
    const std::string &global_vector_name,
    const Extractor &  variable,
    const std::string &object_type,
    const unsigned int size,
    const Number &     exemplar_number) const
  {
    return global_vector_name + "_" + variable.get_name() + "_" + object_type +
           "_" + Utilities::int_to_string(size) + "_" +
           Utilities::type_to_string(exemplar_number);
  }



  template <int dim, int spacedim>
  template <typename Number>
  std::string
  ScratchData<dim, spacedim>::get_unique_dofs_name(
    const std::string &global_vector_name,
    const unsigned int size,
    const Number &     exemplar_number) const
  {
    return global_vector_name + "_independent_local_dofs_" +
           Utilities::int_to_string(size) + "_" +
           Utilities::type_to_string(exemplar_number);
  }



  template <int dim, int spacedim>
  template <typename VectorType, typename Number>
  void
  ScratchData<dim, spacedim>::extract_local_dof_values(
    const std::string &global_vector_name,
    const VectorType & input_vector,
    const Number       dummy)
  {
    const unsigned int n_dofs = local_dof_indices.size();

    const std::string name =
      get_unique_dofs_name(global_vector_name, n_dofs, dummy);

    auto &independent_local_dofs =
      internal_data_storage
        .template get_or_add_object_with_name<std::vector<Number>>(name,
                                                                   n_dofs);

    AssertDimension(independent_local_dofs.size(), n_dofs);

    if (Differentiation::AD::is_tapeless_ad_number<Number>::value == true)
      for (unsigned int i = 0; i < n_dofs; ++i)
        Differentiation::AD::internal::Marking<Number>::independent_variable(
          input_vector(local_dof_indices[i]),
          i,
          n_dofs,
          independent_local_dofs[i]);
    else
      for (unsigned int i = 0; i < n_dofs; ++i)
        independent_local_dofs[i] = input_vector(local_dof_indices[i]);
  }



  template <int dim, int spacedim>
  template <typename Number>
  const std::vector<Number> &
  ScratchData<dim, spacedim>::get_local_dof_values(
    const std::string &global_vector_name,
    Number             dummy) const
  {
    const unsigned int n_dofs =
      get_current_fe_values().get_fe().n_dofs_per_cell();

    const std::string dofs_name =
      get_unique_dofs_name(global_vector_name, n_dofs, dummy);

    Assert(
      internal_data_storage.stores_object_with_name(dofs_name),
      ExcMessage(
        "You did not call yet extract_local_dof_values with the right types!"));

    return internal_data_storage
      .template get_object_with_name<std::vector<Number>>(dofs_name);
  }



  template <int dim, int spacedim>
  template <typename Extractor, typename Number>
  const std::vector<typename FEValuesViews::View<dim, spacedim, Extractor>::
                      template solution_value_type<Number>> &
  ScratchData<dim, spacedim>::get_values(const std::string &global_vector_name,
                                         const Extractor &  variable,
                                         const Number       dummy)
  {
    const std::vector<Number> &independent_local_dofs =
      get_local_dof_values(global_vector_name, dummy);

    const FEValuesBase<dim, spacedim> &fev = get_current_fe_values();

    const unsigned int n_q_points = fev.n_quadrature_points;

    const std::string name = get_unique_name(
      global_vector_name, variable, "_values_q", n_q_points, dummy);

    // Now build the return type
    using RetType =
      std::vector<typename FEValuesViews::View<dim, spacedim, Extractor>::
                    template solution_value_type<Number>>;

    RetType &ret =
      internal_data_storage.template get_or_add_object_with_name<RetType>(
        name, n_q_points);

    AssertDimension(ret.size(), n_q_points);

    fev[variable].get_function_values_from_local_dof_values(
      independent_local_dofs, ret);
    return ret;
  }



  template <int dim, int spacedim>
  template <typename Extractor, typename Number>
  const std::vector<typename FEValuesViews::View<dim, spacedim, Extractor>::
                      template solution_gradient_type<Number>> &
  ScratchData<dim, spacedim>::get_gradients(
    const std::string &global_vector_name,
    const Extractor &  variable,
    const Number       dummy)
  {
    const std::vector<Number> &independent_local_dofs =
      get_local_dof_values(global_vector_name, dummy);

    const FEValuesBase<dim, spacedim> &fev = get_current_fe_values();

    const unsigned int n_q_points = fev.n_quadrature_points;

    const std::string name = get_unique_name(
      global_vector_name, variable, "_gradients_q", n_q_points, dummy);

    // Now build the return type
    using RetType =
      std::vector<typename FEValuesViews::View<dim, spacedim, Extractor>::
                    template solution_gradient_type<Number>>;

    RetType &ret =
      internal_data_storage.template get_or_add_object_with_name<RetType>(
        name, n_q_points);

    AssertDimension(ret.size(), n_q_points);

    fev[variable].get_function_gradients_from_local_dof_values(
      independent_local_dofs, ret);
    return ret;
  }



  template <int dim, int spacedim>
  template <typename Extractor, typename Number>
  const std::vector<typename FEValuesViews::View<dim, spacedim, Extractor>::
                      template solution_hessian_type<Number>> &
  ScratchData<dim, spacedim>::get_hessians(
    const std::string &global_vector_name,
    const Extractor &  variable,
    const Number       dummy)
  {
    const std::vector<Number> &independent_local_dofs =
      get_local_dof_values(global_vector_name, dummy);

    const FEValuesBase<dim, spacedim> &fev = get_current_fe_values();

    const unsigned int n_q_points = fev.n_quadrature_points;

    const std::string name = get_unique_name(
      global_vector_name, variable, "_hessians_q", n_q_points, dummy);

    // Now build the return type
    using RetType =
      std::vector<typename FEValuesViews::View<dim, spacedim, Extractor>::
                    template solution_hessian_type<Number>>;

    RetType &ret =
      internal_data_storage.template get_or_add_object_with_name<RetType>(
        name, n_q_points);


    AssertDimension(ret.size(), n_q_points);

    fev[variable].get_function_hessians_from_local_dof_values(
      independent_local_dofs, ret);
    return ret;
  }



  template <int dim, int spacedim>
  template <typename Extractor, typename Number>
  const std::vector<typename FEValuesViews::View<dim, spacedim, Extractor>::
                      template solution_laplacian_type<Number>> &
  ScratchData<dim, spacedim>::get_laplacians(
    const std::string &global_vector_name,
    const Extractor &  variable,
    const Number       dummy)
  {
    const std::vector<Number> &independent_local_dofs =
      get_local_dof_values(global_vector_name, dummy);

    const FEValuesBase<dim, spacedim> &fev = get_current_fe_values();

    const unsigned int n_q_points = fev.n_quadrature_points;

    const std::string name = get_unique_name(
      global_vector_name, variable, "_laplacians_q", n_q_points, dummy);

    // Now build the return type
    using RetType =
      std::vector<typename FEValuesViews::View<dim, spacedim, Extractor>::
                    template solution_laplacian_type<Number>>;

    RetType &ret =
      internal_data_storage.template get_or_add_object_with_name<RetType>(
        name, n_q_points);


    AssertDimension(ret.size(), n_q_points);

    fev[variable].get_function_laplacians_from_local_dof_values(
      independent_local_dofs, ret);
    return ret;
  }



  template <int dim, int spacedim>
  template <typename Extractor, typename Number>
  const std::vector<typename FEValuesViews::View<dim, spacedim, Extractor>::
                      template solution_third_derivative_type<Number>> &
  ScratchData<dim, spacedim>::get_third_derivatives(
    const std::string &global_vector_name,
    const Extractor &  variable,
    const Number       dummy)
  {
    const std::vector<Number> &independent_local_dofs =
      get_local_dof_values(global_vector_name, dummy);

    const FEValuesBase<dim, spacedim> &fev = get_current_fe_values();

    const unsigned int n_q_points = fev.n_quadrature_points;

    const std::string name = get_unique_name(
      global_vector_name, variable, "_third_derivatives_q", n_q_points, dummy);

    // Now build the return type
    using RetType =
      std::vector<typename FEValuesViews::View<dim, spacedim, Extractor>::
                    template solution_third_derivative_type<Number>>;

    RetType &ret =
      internal_data_storage.template get_or_add_object_with_name<RetType>(
        name, n_q_points);


    AssertDimension(ret.size(), n_q_points);

    fev[variable].get_function_third_derivatives_from_local_dof_values(
      independent_local_dofs, ret);
    return ret;
  }



  template <int dim, int spacedim>
  template <typename Extractor, typename Number>
  const std::vector<typename FEValuesViews::View<dim, spacedim, Extractor>::
                      template solution_symmetric_gradient_type<Number>> &
  ScratchData<dim, spacedim>::get_symmetric_gradients(
    const std::string &global_vector_name,
    const Extractor &  variable,
    const Number       dummy)
  {
    const std::vector<Number> &independent_local_dofs =
      get_local_dof_values(global_vector_name, dummy);

    const FEValuesBase<dim, spacedim> &fev = get_current_fe_values();

    const unsigned int n_q_points = fev.n_quadrature_points;

    const std::string name = get_unique_name(
      global_vector_name, variable, "_symmetric_gradient_q", n_q_points, dummy);


    // Now build the return type
    using RetType =
      std::vector<typename FEValuesViews::View<dim, spacedim, Extractor>::
                    template solution_symmetric_gradient_type<Number>>;

    RetType &ret =
      internal_data_storage.template get_or_add_object_with_name<RetType>(
        name, n_q_points);


    AssertDimension(ret.size(), n_q_points);

    fev[variable].get_function_symmetric_gradients_from_local_dof_values(
      independent_local_dofs, ret);
    return ret;
  }


  template <int dim, int spacedim>
  template <typename Extractor, typename Number>
  const std::vector<typename FEValuesViews::View<dim, spacedim, Extractor>::
                      template solution_divergence_type<Number>> &
  ScratchData<dim, spacedim>::get_divergences(
    const std::string &global_vector_name,
    const Extractor &  variable,
    const Number       dummy)
  {
    const std::vector<Number> &independent_local_dofs =
      get_local_dof_values(global_vector_name, dummy);

    const FEValuesBase<dim, spacedim> &fev = get_current_fe_values();

    const unsigned int n_q_points = fev.n_quadrature_points;

    const std::string name = get_unique_name(
      global_vector_name, variable, "_divergence_q", n_q_points, dummy);

    // Now build the return type
    using RetType =
      std::vector<typename FEValuesViews::View<dim, spacedim, Extractor>::
                    template solution_divergence_type<Number>>;

    RetType &ret =
      internal_data_storage.template get_or_add_object_with_name<RetType>(
        name, n_q_points);


    AssertDimension(ret.size(), n_q_points);

    fev[variable].get_function_divergences_from_local_dof_values(
      independent_local_dofs, ret);
    return ret;
  }



  template <int dim, int spacedim>
  template <typename Extractor, typename Number>
  const std::vector<typename FEValuesViews::View<dim, spacedim, Extractor>::
                      template solution_curl_type<Number>> &
  ScratchData<dim, spacedim>::get_curls(const std::string &global_vector_name,
                                        const Extractor &  variable,
                                        const Number       dummy)
  {
    const std::vector<Number> &independent_local_dofs =
      get_local_dof_values(global_vector_name, dummy);

    const FEValuesBase<dim, spacedim> &fev = get_current_fe_values();

    const unsigned int n_q_points = fev.n_quadrature_points;

    const std::string name = get_unique_name(
      global_vector_name, variable, "_curl_q", n_q_points, dummy);

    // Now build the return type
    using RetType =
      std::vector<typename FEValuesViews::View<dim, spacedim, Extractor>::
                    template solution_curl_type<Number>>;

    RetType &ret =
      internal_data_storage.template get_or_add_object_with_name<RetType>(
        name, n_q_points);

    AssertDimension(ret.size(), n_q_points);

    fev[variable].get_function_curls_from_local_dof_values(
      independent_local_dofs, ret);
    return ret;
  }

#endif

} // namespace MeshWorker

DEAL_II_NAMESPACE_CLOSE

#endif


