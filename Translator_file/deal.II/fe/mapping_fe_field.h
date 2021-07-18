//include/deal.II-translator/fe/mapping_fe_field_0.txt
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

#ifndef dealii_mapping_fe_field_h
#define dealii_mapping_fe_field_h


#include <deal.II/base/config.h>

#include <deal.II/base/mg_level_object.h>
#include <deal.II/base/thread_management.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping.h>

#include <deal.II/lac/vector.h>

#include <array>


DEAL_II_NAMESPACE_OPEN


 /*!@addtogroup mapping */ 
 /*@{*/ 

/**
 * @deprecated  使用MappingFEField<dim, spacedim, VectorType>代替。
 *
 *
 */
template <int dim,
          int spacedim            = dim,
          typename VectorType     = Vector<double>,
          typename DoFHandlerType = void>
class MappingFEField;

#ifndef DOXYGEN
// prevent doxygen from complaining about potential recursive class relations
template <int dim, int spacedim, typename VectorType, typename DoFHandlerType>
class MappingFEField : public MappingFEField<dim, spacedim, VectorType, void>
{
public:
  DEAL_II_DEPRECATED
  MappingFEField(const DoFHandlerType &euler_dof_handler,
                 const VectorType &    euler_vector,
                 const ComponentMask & mask = ComponentMask())
    : MappingFEField<dim, spacedim, VectorType, void>(euler_dof_handler,
                                                      euler_vector,
                                                      mask)
  {}

  DEAL_II_DEPRECATED
  MappingFEField(const DoFHandlerType &         euler_dof_handler,
                 const std::vector<VectorType> &euler_vector,
                 const ComponentMask &          mask = ComponentMask())
    : MappingFEField<dim, spacedim, VectorType, void>(euler_dof_handler,
                                                      euler_vector,
                                                      mask)
  {}

  DEAL_II_DEPRECATED
  MappingFEField(const DoFHandlerType &           euler_dof_handler,
                 const MGLevelObject<VectorType> &euler_vector,
                 const ComponentMask &            mask = ComponentMask())
    : MappingFEField<dim, spacedim, VectorType, void>(euler_dof_handler,
                                                      euler_vector,
                                                      mask)
  {}

  DEAL_II_DEPRECATED
  MappingFEField(
    const MappingFEField<dim, spacedim, VectorType, DoFHandlerType> &mapping)
    : MappingFEField<dim, spacedim, VectorType, void>(mapping)
  {}
};
#endif // DOXYGEN

/**
 * MappingFEField是MappingQEulerian类的一个泛化，用于任意矢量有限元。两个主要的区别是，这个类使用绝对位置的矢量，而不是位移的矢量，而且它允许任意的有限元类型，而不是只有FE_Q。
 * 该类有效地将拓扑结构与几何结构相分离，将所有的几何信息归结为FiniteElement矢量场的某些组件。用于几何学的分量可以在构造时任意选择。
 * 我们的想法是将三角剖分视为一个参数配置空间，在此基础上构建一个任意的几何映射，使用deal.II库的工具：一个自由度矢量，一个与问题几何相关的DoFHandler，以及一个ComponentMask，告诉我们FiniteElement的哪些组件要用于映射。
 * 通常，DoFHandler操作的有限元是由连续的FE_Q()（用于等参数离散）或FE_Bernstein()（用于等几何离散）对象构建的系统元（FESystem()
 * ）。下面是一个例子。
 *
 *
 * @code
 *  const FE_Q<dim,spacedim> feq(1);
 *  const FESystem<dim,spacedim> fesystem(feq, spacedim);
 *  DoFHandler<dim,spacedim> dhq(triangulation);
 *  dhq.distribute_dofs(fesystem);
 *  const ComponentMask mask(spacedim, true);
 *  Vector<double> eulerq(dhq.n_dofs());
 *  // Fills the euler vector with information from the Triangulation
 *  VectorTools::get_position_vector(dhq, eulerq, mask);
 *  MappingFEField<dim,spacedim> map(dhq, eulerq, mask);
 * @endcode
 *
 *
 *
 */
template <int dim, int spacedim, typename VectorType>
class MappingFEField<dim, spacedim, VectorType, void>
  : public Mapping<dim, spacedim>
{
public:
  /**
   * 构造函数。第一个参数是一个VectorType，指定域从参考到当前配置的转换。
   * 一般来说，这个类将几何学与拓扑学解耦，允许用户定义仅在拓扑学上与底层三角结构等价的几何学，但其他方面可能是任意的。
   * 与MappingQEulerian不同的是，传递给构造函数的FiniteElement字段被解释为一个绝对的几何配置，因此我们必须确保euler_vector实际代表一个有效的几何体（即没有倒置的单元或没有零体积的单元）。
   * 如果底层的FiniteElement是FE_Q()的一个系统，并且euler_vector是用
   * VectorTools::get_position_vector(),
   * 初始化的，那么这个类在所有方面都与MappingQ()相同。
   * 可选的ComponentMask参数可以用来指定FiniteElement的哪些组件要用于几何变换。如果在构造时没有指定掩码，那么将使用默认的掩码，这使得该类的工作方式与MappingQEulerian()相同，即假定FiniteElement的第一个间隔分量代表问题的几何形状。
   * 请注意，如果指定了一个掩码，它的大小必须与底层的FiniteElement相匹配，并且它必须有完全间隔的非零元素，表明FiniteElement的组件（按顺序）将被用于几何。
   * 如果传递了一个不兼容的掩码，就会抛出一个异常。
   *
   */
  MappingFEField(const DoFHandler<dim, spacedim> &euler_dof_handler,
                 const VectorType &               euler_vector,
                 const ComponentMask &            mask = ComponentMask());

  /**
   * 构造函数在多网格层上获取向量，而不是只获取活动单元。向量的向量应该有与三角形中的全局层次一样多的条目，并在每个层次上提供有效的数据，即长度兼容
   * DoFHandler::n_dofs(level).  这个构造函数的前提条件是
   * DoFHandler::distribute_mg_dofs()
   * 已经被调用。除了级别向量外，还需要提供与其他构造函数相同的参数。
   *
   */
  MappingFEField(const DoFHandler<dim, spacedim> &euler_dof_handler,
                 const std::vector<VectorType> &  euler_vector,
                 const ComponentMask &            mask = ComponentMask());

  /**
   * 用MGLevelObject代替 std::vector,
   * 的构造函数，否则与上述相同。要求`euler_vector.max_level()+1`等于三角形中的全局级数。最小级别可以是0或更多&mdash;
   * 它只需要在这里设置的内容和后来用于评估映射的内容之间保持一致。
   *
   */
  MappingFEField(const DoFHandler<dim, spacedim> &euler_dof_handler,
                 const MGLevelObject<VectorType> &euler_vector,
                 const ComponentMask &            mask = ComponentMask());

  /**
   * 复制构造函数。
   *
   */
  MappingFEField(
    const MappingFEField<dim, spacedim, VectorType, void> &mapping);

  /**
   * 返回一个指向当前对象副本的指针。然后，这个副本的调用者就拥有了它的所有权。
   *
   */
  virtual std::unique_ptr<Mapping<dim, spacedim>>
  clone() const override;

  /**
   * 关于这个函数的目的，请参见
   * Mapping::preserves_vertex_locations()
   * 的文档。这个类中的实现总是返回 @p false. 。
   *
   */
  virtual bool
  preserves_vertex_locations() const override;

  virtual bool
  is_compatible_with(const ReferenceCell &reference_cell) const override;

  /**
   * 返回一个单元格的映射顶点。
   * 这个映射忽略了它所关联的三角形的顶点，并根据构造时传递的
   * @p euler_vector 来构造顶点的位置。
   *
   */
  virtual boost::container::small_vector<Point<spacedim>,
                                         GeometryInfo<dim>::vertices_per_cell>
  get_vertices(const typename Triangulation<dim, spacedim>::cell_iterator &cell)
    const override;

  /**
   * @name  参考单元和实数单元之间的映射点 
     * @{ 
   *
   */

  // for documentation, see the Mapping base class
  virtual Point<spacedim>
  transform_unit_to_real_cell(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const Point<dim> &p) const override;

  // for documentation, see the Mapping base class
  virtual Point<dim>
  transform_real_to_unit_cell(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const Point<spacedim> &p) const override;

  /**
   * @}
   *
   */

  /**
   * @name  将张量从参考坐标转换为实坐标的函数  @{
   *
   */

  // for documentation, see the Mapping base class
  virtual void
  transform(const ArrayView<const Tensor<1, dim>> &                  input,
            const MappingKind                                        kind,
            const typename Mapping<dim, spacedim>::InternalDataBase &internal,
            const ArrayView<Tensor<1, spacedim>> &output) const override;

  // for documentation, see the Mapping base class
  virtual void
  transform(const ArrayView<const DerivativeForm<1, dim, spacedim>> &input,
            const MappingKind                                        kind,
            const typename Mapping<dim, spacedim>::InternalDataBase &internal,
            const ArrayView<Tensor<2, spacedim>> &output) const override;

  // for documentation, see the Mapping base class
  virtual void
  transform(const ArrayView<const Tensor<2, dim>> &                  input,
            const MappingKind                                        kind,
            const typename Mapping<dim, spacedim>::InternalDataBase &internal,
            const ArrayView<Tensor<2, spacedim>> &output) const override;

  // for documentation, see the Mapping base class
  virtual void
  transform(const ArrayView<const DerivativeForm<2, dim, spacedim>> &input,
            const MappingKind                                        kind,
            const typename Mapping<dim, spacedim>::InternalDataBase &internal,
            const ArrayView<Tensor<3, spacedim>> &output) const override;

  // for documentation, see the Mapping base class
  virtual void
  transform(const ArrayView<const Tensor<3, dim>> &                  input,
            const MappingKind                                        kind,
            const typename Mapping<dim, spacedim>::InternalDataBase &internal,
            const ArrayView<Tensor<3, spacedim>> &output) const override;

  /**
   *
   */

  /**
   * 返回映射的程度，即传递给构造函数的值。
   *
   */
  unsigned int
  get_degree() const;

  /**
   * 返回映射的ComponentMask，即哪些组件要用于映射。
   *
   */
  ComponentMask
  get_component_mask() const;

  /**
   * 异常情况
   *
   */
  DeclException0(ExcInactiveCell);

private:
  /**
   * @name  与FEValues的接口 
     * @{ 
   *
   */

  // documentation can be found in Mapping::requires_update_flags()
  virtual UpdateFlags
  requires_update_flags(const UpdateFlags update_flags) const override;

public:
  /**
   * 此映射的内部数据的存储。见 Mapping::InternalDataBase
   * 的广泛描述。
   * 这包括在创建对象时（在get_data()中）计算一次的数据，以及该类希望从调用fill_fe_values()、fill_fe_face_values()或fill_fe_subface_values()之间存储的数据，直到以后可能从有限元调用转化()等函数。后一类成员变量与从头开始的数组一起被标记为
   * "可变"。
   *
   */
  class InternalData : public Mapping<dim, spacedim>::InternalDataBase
  {
  public:
    /**
     * 构造函数。
     *
     */
    InternalData(const FiniteElement<dim, spacedim> &fe,
                 const ComponentMask &               mask);

    /**
     * 正交点的形状函数。形状函数是按张量积顺序排列的，所以顶点必须重新排序以获得变换。
     *
     */
    const double &
    shape(const unsigned int qpoint, const unsigned int shape_nr) const;

    /**
     * 正交点的形状函数。见上文。
     *
     */
    double &
    shape(const unsigned int qpoint, const unsigned int shape_nr);

    /**
     * 形状函数在正交点的梯度。见上文。
     *
     */
    const Tensor<1, dim> &
    derivative(const unsigned int qpoint, const unsigned int shape_nr) const;

    /**
     * 形状函数在正交点的梯度。见上文。
     *
     */
    Tensor<1, dim> &
    derivative(const unsigned int qpoint, const unsigned int shape_nr);

    /**
     * 形状函数在正交点的二阶导数。见上文。
     *
     */
    const Tensor<2, dim> &
    second_derivative(const unsigned int qpoint,
                      const unsigned int shape_nr) const;

    /**
     * 形状函数在正交点的二阶导数。见上文。
     *
     */
    Tensor<2, dim> &
    second_derivative(const unsigned int qpoint, const unsigned int shape_nr);

    /**
     * 形状函数在正交点的三次导数。见上文。
     *
     */
    const Tensor<3, dim> &
    third_derivative(const unsigned int qpoint,
                     const unsigned int shape_nr) const;

    /**
     * 在正交点的形状函数的第四导数。见上文。
     *
     */
    Tensor<3, dim> &
    third_derivative(const unsigned int qpoint, const unsigned int shape_nr);

    /**
     * 形状函数在正交点的四次导数。见上文。
     *
     */
    const Tensor<4, dim> &
    fourth_derivative(const unsigned int qpoint,
                      const unsigned int shape_nr) const;

    /**
     * 形状函数在正交点的三次导数。见上文。
     *
     */
    Tensor<4, dim> &
    fourth_derivative(const unsigned int qpoint, const unsigned int shape_nr);

    /**
     * 返回这个对象的内存消耗估计值（以字节为单位）。
     *
     */
    virtual std::size_t
    memory_consumption() const override;

    /**
     * 形状函数的值。通过函数访问  @p shape.  计算一次。
     *
     */
    std::vector<double> shape_values;

    /**
     * 形状函数导数的值。通过函数访问  @p derivative.
     * 计算一次。
     *
     */
    std::vector<Tensor<1, dim>> shape_derivatives;

    /**
     * 形状函数二次导数的值。通过函数 @p
     * second_derivative访问。        计算一次。
     *
     */
    std::vector<Tensor<2, dim>> shape_second_derivatives;

    /**
     * 形状函数第三导数的值。通过函数 @p
     * third_derivative访问。        计算一次。
     *
     */
    std::vector<Tensor<3, dim>> shape_third_derivatives;

    /**
     * 形状函数第四导数的值。通过函数 @p
     * fourth_derivative访问。        计算一次。
     *
     */
    std::vector<Tensor<4, dim>> shape_fourth_derivatives;

    /**
     * 单位切向量。用于计算边界形式和法向量。
     * 这个数组有 `(dim-1)*GeometryInfo<dim>::%faces_per_cell`
     * 个条目。第一组 GeometryInfo::faces_per_cell
     * 包含每个面的第一个切向的向量；第二组
     * GeometryInfo::faces_per_cell
     * 条目包含第二个切向的向量（仅在三维中，因为每个面有两个切向），等等。
     * 填充一次。
     *
     */
    std::array<std::vector<Tensor<1, dim>>,
               GeometryInfo<dim>::faces_per_cell *(dim - 1)>
      unit_tangentials;

    /**
     * 形状函数的数量。如果这是一个Q1的映射，那么它只是每个单元格的顶点数量。然而，因为还有派生类使用这个类（例如Mapping_Q()类），所以形状函数的数量也可能不同。
     *
     */
    unsigned int n_shape_functions;

    /**
     * 存储构造时给出的掩码。如果在构造时没有指定掩码，那么将使用默认的掩码，这使得该类与MappingQEulerian()的工作方式相同，即FiniteElement的第一个间隔分量被用于euler_vector和euler_dh。
     * 如果指定了一个掩码，那么它必须与底层的FiniteElement相匹配，并且它必须有完全间隔的非零元素，表明FiniteElement的组件（按顺序）将被用于欧拉向量和欧拉道夫处理。
     *
     */
    ComponentMask mask;

    /**
     * 每个正交点的协方差变换的张量。
     * 存储的矩阵是Jacobian G^{-1}，其中G = Jacobian^{t}
     * Jacobian，是地图的第一基本形式；如果dim=spacedim，那么它还原为Jacobian矩阵的转置，它本身存储在这个结构的
     * @p contravariant 域中。        在每个单元格上计算。
     *
     */
    mutable std::vector<DerivativeForm<1, dim, spacedim>> covariant;

    /**
     * 每个正交点上的禁忌变换的张量。不变矩阵是变换的雅各布系数，即
     * $J_{ij}=dx_i/d\hat x_j$  。        在每个单元上计算。
     *
     */
    mutable std::vector<DerivativeForm<1, dim, spacedim>> contravariant;

    /**
     * 在每个正交点的雅各布系数的行列式。如果#update_volume_elements就会被填满。
     *
     */
    mutable std::vector<double> volume_elements;

    /**
     * 供内部使用的辅助向量。
     *
     */
    mutable std::vector<std::vector<Tensor<1, spacedim>>> aux;

    /**
     * 储存局部自由度的指数。
     *
     */
    mutable std::vector<types::global_dof_index> local_dof_indices;

    /**
     * 本地自由度的存储。
     *
     */
    mutable std::vector<double> local_dof_values;
  };

private:
  // documentation can be found in Mapping::get_data()
  virtual std::unique_ptr<typename Mapping<dim, spacedim>::InternalDataBase>
  get_data(const UpdateFlags, const Quadrature<dim> &quadrature) const override;

  using Mapping<dim, spacedim>::get_face_data;

  // documentation can be found in Mapping::get_face_data()
  virtual std::unique_ptr<typename Mapping<dim, spacedim>::InternalDataBase>
  get_face_data(const UpdateFlags               flags,
                const hp::QCollection<dim - 1> &quadrature) const override;

  // documentation can be found in Mapping::get_subface_data()
  virtual std::unique_ptr<typename Mapping<dim, spacedim>::InternalDataBase>
  get_subface_data(const UpdateFlags          flags,
                   const Quadrature<dim - 1> &quadrature) const override;

  // documentation can be found in Mapping::fill_fe_values()
  virtual CellSimilarity::Similarity
  fill_fe_values(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const CellSimilarity::Similarity                            cell_similarity,
    const Quadrature<dim> &                                     quadrature,
    const typename Mapping<dim, spacedim>::InternalDataBase &   internal_data,
    internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
      &output_data) const override;

  using Mapping<dim, spacedim>::fill_fe_face_values;

  // documentation can be found in Mapping::fill_fe_face_values()
  virtual void
  fill_fe_face_values(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const unsigned int                                          face_no,
    const hp::QCollection<dim - 1> &                            quadrature,
    const typename Mapping<dim, spacedim>::InternalDataBase &   internal_data,
    internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
      &output_data) const override;

  // documentation can be found in Mapping::fill_fe_subface_values()
  virtual void
  fill_fe_subface_values(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const unsigned int                                          face_no,
    const unsigned int                                          subface_no,
    const Quadrature<dim - 1> &                                 quadrature,
    const typename Mapping<dim, spacedim>::InternalDataBase &   internal_data,
    internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
      &output_data) const override;

  /**
   * @}
   *
   */

  /**
   * 指定我们是在活动自由度上（用一个欧拉向量）还是在水平自由度上（通过一个欧拉向量的向量）访问未知数。
   *
   */
  const bool uses_level_dofs;

  /**
   * 对移位矢量的参考。
   *
   */
  std::vector<SmartPointer<const VectorType,
                           MappingFEField<dim, spacedim, VectorType, void>>>
    euler_vector;

  /**
   * 指向映射向量所关联的DoFHandler的指针。
   *
   */
  SmartPointer<const DoFHandler<dim, spacedim>,
               MappingFEField<dim, spacedim, VectorType, void>>
    euler_dof_handler;

private:
  /**
   * 将单元格上的点 @p p 转换为实数单元格 @p cell 上的点 @p
   * p_real ，并返回 @p p_real.  这个函数被 @p
   * transform_unit_to_real_cell 调用，并被 @p
   * transform_real_to_unit_cell_internal(牛顿反复)多次调用。
   * 接收一个对 @p InternalData 的引用，该引用必须已经包括
   * @p p 点的形状值和单元格的映射支持点。    这个 @p
   * InternalData 参数可以避免对 @p p
   * 点的形状值进行多次计算，特别是对映射支持点进行多次计算。
   *
   */
  Point<spacedim>
  do_transform_unit_to_real_cell(const InternalData &mdata) const;


  /**
   * 通过牛顿迭代将实细胞上的点 @p p 转换为单元格 @p cell
   * 上的相应点。    取一个对 @p InternalData
   * 的引用，该引用假定是由 @p get_data 函数与 @p UpdateFlags
   * 包括 @p update_transformation_values和 @p
   * update_transformation_gradients
   * 和一个点的正交所创建的，包括对变换的初始猜测 @p
   * initial_p_unit.  ] 因此，这个函数假定 @p
   * mdata已经包括了变换的形状值和在 @p initial_p_unit.  @p mdata
   * 计算的梯度将被这个函数所改变。
   *
   */
  Point<dim>
  do_transform_real_to_unit_cell(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const Point<spacedim> &                                     p,
    const Point<dim> &                                          initial_p_unit,
    InternalData &                                              mdata) const;

  /**
   * 更新内部自由度。
   *
   */
  void
  update_internal_dofs(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const typename MappingFEField<dim, spacedim, VectorType, void>::InternalData
      &data) const;

  /**
   * 详细情况见基类的文档。
   *
   */
  virtual void
  compute_shapes_virtual(
    const std::vector<Point<dim>> &unit_points,
    typename MappingFEField<dim, spacedim, VectorType, void>::InternalData
      &data) const;

  /*在映射中使用哪些组件。 
*
*/
  const ComponentMask fe_mask;

  /**
   * 在FE空间和实空间的指数之间的映射。这个向量包含有限元空间的每个分量的一个索引。如果该索引是一个用于构建该元素的ComponentMask为false，那么将返回
   * numbers::invalid_unsigned_int
   * ，否则将返回实空间的分量。例如，如果我们使用ComponentMask(spacedim,
   * true)构建映射，那么这个向量包含{0,1,2}在spacedim=3。
   *
   */
  std::vector<unsigned int> fe_to_real;

  /**
   * FEValues对象，用于查询参考配置中支持点的给定有限元场。
   *
   */
  mutable FEValues<dim, spacedim> fe_values;

  /**
   * 一个变量用于保护对fe_values变量的访问。
   *
   */
  mutable std::mutex fe_values_mutex;

  void
  compute_data(const UpdateFlags      update_flags,
               const Quadrature<dim> &q,
               const unsigned int     n_original_q_points,
               InternalData &         data) const;

  void
  compute_face_data(const UpdateFlags      update_flags,
                    const Quadrature<dim> &q,
                    const unsigned int     n_original_q_points,
                    InternalData &         data) const;


  // Declare other MappingFEField classes friends.
  template <int, int, class, class>
  friend class MappingFEField;
};

 /*@}*/ 

 /* -------------- declaration of explicit specializations ------------- */ 

#ifndef DOXYGEN

#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif


