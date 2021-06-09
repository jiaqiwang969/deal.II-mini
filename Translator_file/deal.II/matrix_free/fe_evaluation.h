//include/deal.II-translator/matrix_free/fe_evaluation_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2021 by the deal.II authors
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


#ifndef dealii_matrix_free_fe_evaluation_h
#define dealii_matrix_free_fe_evaluation_h


#include <deal.II/base/config.h>

#include <deal.II/base/array_view.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/template_constraints.h>
#include <deal.II/base/vectorization.h>

#include <deal.II/lac/vector_operation.h>

#include <deal.II/matrix_free/evaluation_flags.h>
#include <deal.II/matrix_free/evaluation_kernels.h>
#include <deal.II/matrix_free/evaluation_selector.h>
#include <deal.II/matrix_free/evaluation_template_factory.h>
#include <deal.II/matrix_free/mapping_data_on_the_fly.h>
#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/shape_info.h>
#include <deal.II/matrix_free/tensor_product_kernels.h>
#include <deal.II/matrix_free/type_traits.h>
#include <deal.II/matrix_free/vector_access_internal.h>


DEAL_II_NAMESPACE_OPEN



namespace internal
{
  DeclException0(ExcAccessToUninitializedField);

  DeclException1(
    ExcMatrixFreeAccessToUninitializedMappingField,
    std::string,
    << "You are requesting information from an FEEvaluation/FEFaceEvaluation "
    << "object for which this kind of information has not been computed. What "
    << "information these objects compute is determined by the update_* flags you "
    << "pass to MatrixFree::reinit() via MatrixFree::AdditionalData. Here, "
    << "the operation you are attempting requires the <" << arg1
    << "> flag to be set, but it was apparently not specified "
    << "upon initialization.");
} // namespace internal

template <int dim,
          int fe_degree,
          int n_q_points_1d            = fe_degree + 1,
          int n_components_            = 1,
          typename Number              = double,
          typename VectorizedArrayType = VectorizedArray<Number>>
class FEEvaluation;



/**
 * 这个FEEvaluation和FEFaceEvaluation类的基类处理与映射相关的信息，与使用中的自由度和有限元无关。该类为用户代码提供了访问功能，但除此之外没有任何公共构造函数，是不可见的。使用方法是通过FEEvaluation类来代替。
 * 这个类有四个模板参数。
 * @tparam  该类要使用的dim尺寸
 * @tparam  Number 数字格式，通常是 @p double 或 @p float  。
 * @tparam  is_face
 * 该类是用于单元积分器（正交尺寸与空间尺寸相同）还是用于面积分器（正交尺寸少一个）。
 * @tparam  VectorizedArrayType
 * 以矢量方式处理的数组类型，默认为 VectorizedArray<Number>。
 *
 *
 * @note  目前只有VectorizedArray<Number,
 * width>被支持作为VectorizedArrayType。
 *
 *
 *
 * @ingroup matrixfree
 *
 */
template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
class FEEvaluationBaseData
{
  static_assert(
    std::is_same<Number, typename VectorizedArrayType::value_type>::value,
    "Type of Number and of VectorizedArrayType do not match.");

public:
  static constexpr unsigned int dimension = dim;

  /**
   * 解构器。
   *
   */
  ~FEEvaluationBaseData();

  /**
   * 返回 @p
   * reinit()函数所调用的单元格在几何字段中的索引偏移。这个索引可以用来访问一个字段的索引，这个字段的压缩行为与几何体的Jacobian相同，例如，存储一个有效的系数tensors，将系数与几何体结合起来，以降低内存传输，作为可用的数据字段。
   *
   */
  unsigned int
  get_mapping_data_index_offset() const;

  /**
   * 返回 @p reinit() 函数被调用的单元格的类型。  有效值是
   * @p cartesian 用于笛卡尔单元（允许相当大的数据压缩），
   * @p affine 用于具有仿射映射的单元， @p general
   * 用于没有应用任何压缩存储的一般单元。
   *
   */
  internal::MatrixFreeFunctions::GeometryType
  get_cell_type() const;

  /**
   * 返回当前正在使用的ShapeInfo对象的引用。
   *
   */
  const internal::MatrixFreeFunctions::ShapeInfo<VectorizedArrayType> &
  get_shape_info() const;

  /**
   * 返回当前正在使用的DoFInfo对象的引用。
   *
   */
  const internal::MatrixFreeFunctions::DoFInfo &
  get_dof_info() const;

  /**
   * 返回从单位到实数单元的雅各布系数乘以正交权重的行列式。
   *
   */
  VectorizedArrayType
  JxW(const unsigned int q_point) const;

  /**
   * 返回定义为  $J_{ij} = d x_i / d\hat x_j$
   * 的单位到实数单元之间映射的雅各布系数的反转和转置版本
   * $J^{-\mathrm T}$  。返回张量的 $(i,j)$ 项包含 $d\hat x_j/dx_i$
   * ，即列是指参考空间坐标，行是指实单元坐标。因此，返回的张量代表一个协变变换，在
   * FEEvaluationBase::get_gradient() 函数中用于通过乘法 $J^{-\mathrm
   * T} \hat{\nabla} u_h$
   * 将单元格梯度转换为实单元格上的梯度。
   *
   */
  Tensor<2, dim, VectorizedArrayType>
  inverse_jacobian(const unsigned int q_point) const;

  /**
   * 返回一个面的单位法向量。注意，一个面的两边都使用相同方向的法向量。对于在FaceToCellTopology中被列举为
   * "内部 "并在构造函数中被选择为 "is_interior_face=true
   * "的面，这对应于外部法向量，而对于在FaceToCellTopology中被列举为
   * "外部 "并在构造函数中被选择为 "is_interior_face=false
   * "的面，作为单一法向量的结果，法向量指向该元素中。
   * @note 只在`is_face == true`的情况下实现。
   *
   */
  Tensor<1, dim, VectorizedArrayType>
  get_normal_vector(const unsigned int q_point) const;

  /**
   * 提供一个统一的接口来访问长度为
   * MatrixFree::n_cell_batches() + MatrixFree::n_ghost_cell_batches()
   * 的向量Array字段中的数据，用于单元（普通读取）和面（间接寻址）。
   *
   */
  VectorizedArrayType
  read_cell_data(const AlignedVector<VectorizedArrayType> &array) const;

  /**
   * 为单元格（明读）和面（间接寻址）提供一个统一的接口来设置长度为
   * MatrixFree::n_cell_batches()  +  MatrixFree::n_ghost_cell_batches()
   * 的矢量Array字段中的数据。
   *
   */
  void
  set_cell_data(AlignedVector<VectorizedArrayType> &array,
                const VectorizedArrayType &         value) const;

  /**
   * 与上述相同，只是对于任意数据类型的VectorizedArrayType的长度
   * std::array 。
   *
   */
  template <typename T>
  std::array<T, VectorizedArrayType::size()>
  read_cell_data(const AlignedVector<std::array<T, VectorizedArrayType::size()>>
                   &array) const;

  /**
   * 和上面一样，只是对于任意数据类型的VectorizedArrayType的长度为
   * std::array 。
   *
   */
  template <typename T>
  void
  set_cell_data(
    AlignedVector<std::array<T, VectorizedArrayType::size()>> &array,
    const std::array<T, VectorizedArrayType::size()> &         value) const;

  /**
   * 返回这个FEEvaluation或FEFaceEvaluation所关联的单元格的id。
   *
   */
  std::array<unsigned int, VectorizedArrayType::size()>
  get_cell_ids() const;

  /**
   * 返回与此FEEvaluation/FEFaceEvaluation相关的单元格/面的id。
   *
   */
  std::array<unsigned int, VectorizedArrayType::size()>
  get_cell_or_face_ids() const;


  /**
   * 返回FEEvaluation的评估程序中局部自由度的编号，以有限元的标准编号为准。
   *
   */
  const std::vector<unsigned int> &
  get_internal_dof_numbering() const;

  /**
   * 返回一个ArrayView到内部内存，供临时使用。注意，在evaluation()和integration()调用过程中，这部分内存会被覆盖，所以不要认为它在这些调用中是稳定的。你可以写入的最大容量是3*dofs_per_cell+2*n_q_points。
   *
   */
  ArrayView<VectorizedArrayType>
  get_scratch_data() const;

  /**
   * 返回当前单元格的正交公式的编号。
   *
   */
  unsigned int
  get_quadrature_index() const;

  /**
   * 返回当前单元格或面的索引。
   *
   */
  unsigned int
  get_current_cell_index() const;

  /**
   * 返回该类的活动FE索引，以便在hp-情况下有效地进行索引。
   *
   */
  unsigned int
  get_active_fe_index() const;

  /**
   * 返回该类的活动正交索引，以便在hp情况下有效地进行索引。
   *
   */
  unsigned int
  get_active_quadrature_index() const;

  /**
   * 返回底层的MatrixFree对象。
   *
   */
  const MatrixFree<dim, Number, VectorizedArrayType> &
  get_matrix_free() const;

protected:
  /**
   * 构造函数。为了防止用户直接使用这个类，做了保护。采取所有存储在MatrixFree中的数据。如果应用于在构造`matrix_free`时选择了多个正交公式的问题，`quad_no`允许选择适当的公式。
   *
   */
  FEEvaluationBaseData(
    const MatrixFree<dim, Number, VectorizedArrayType> &matrix_free,
    const unsigned int                                  dof_no,
    const unsigned int first_selected_component,
    const unsigned int quad_no,
    const unsigned int fe_degree,
    const unsigned int n_q_points,
    const bool         is_interior_face,
    const unsigned int active_fe_index,
    const unsigned int active_quad_index,
    const unsigned int face_type);

  /**
   * 构造函数，功能减少，工作方式与FEValues类似。
   *
   */
  FEEvaluationBaseData(
    const Mapping<dim> &      mapping,
    const FiniteElement<dim> &fe,
    const Quadrature<1> &     quadrature,
    const UpdateFlags         update_flags,
    const unsigned int        first_selected_component,
    const FEEvaluationBaseData<dim, Number, is_face, VectorizedArrayType>
      *other);

  /**
   * 复制构造函数。如果FEEvaluationBase是由映射、fe、正交和更新标志构建的，基于FEValues的底层几何评估将被深度复制，以便于与线程并行使用。
   *
   */
  FEEvaluationBaseData(const FEEvaluationBaseData &other);

  /**
   * 复制赋值运算符。如果FEEvaluationBase是由映射、fe、正交和更新标志构建的，基于FEValues的底层几何评估将被深度复制，以允许与线程并行使用。
   *
   */
  FEEvaluationBaseData &
  operator=(const FEEvaluationBaseData &other);

  /**
   * 这是所有数据字段的一般阵列。
   *
   */
  AlignedVector<VectorizedArrayType> *scratch_data_array;

  /**
   * 这是 scratch_data_array 中用户可见的部分，只显示
   * scratch_data_array 的最后一部分。第一部分被 values_dofs,
   * values_quad 等消耗了。
   *
   */
  VectorizedArrayType *scratch_data;

  /**
   * 本单元格的正交公式的编号。
   *
   */
  const unsigned int quad_no;

  /**
   * 一个指向基础数据的指针。
   *
   */
  const MatrixFree<dim, Number, VectorizedArrayType> *matrix_info;

  /**
   * 一个指向底层DoF指数和约束描述的指针，用于构造时指定的组件。也包含在matrix_info中，但是如果我们存储一个对它的引用，可以简化代码。
   *
   */
  const internal::MatrixFreeFunctions::DoFInfo *dof_info;

  /**
   * 指向结构中指定的正交公式从单位到实数单元的基础转换数据的指针。也包含在matrix_info中，但如果我们存储对它的引用，可以简化代码。
   *
   */
  const internal::MatrixFreeFunctions::MappingInfoStorage<
    (is_face ? dim - 1 : dim),
    dim,
    Number,
    VectorizedArrayType> *mapping_data;

  /**
   * 该类的活动FE索引，用于在hp情况下的有效索引。
   *
   */
  const unsigned int active_fe_index;

  /**
   * 该类的活动正交索引，用于在hp情况下的有效索引。
   *
   */
  const unsigned int active_quad_index;

  /**
   * 一个指向构造时指定的底层正交公式的指针。
   * 也包含在matrix_info中，但是如果我们存储一个对它的引用，可以简化代码。
   *
   */
  const typename internal::MatrixFreeFunctions::MappingInfoStorage<
    (is_face ? dim - 1 : dim),
    dim,
    Number,
    VectorizedArrayType>::QuadratureDescriptor *descriptor;

  /**
   * 当前评估环境中正交点的数量。
   *
   */
  const unsigned int n_quadrature_points;

  /**
   * 一个指向单元格形状数据的指针，即构成张量积的正交点的值、梯度和一维的Hessians。也包含在matrix_info中，但是如果我们存储一个对它的引用，可以简化代码。
   *
   */
  const internal::MatrixFreeFunctions::ShapeInfo<VectorizedArrayType> *data;

  /**
   * 一个指向当前单元格的雅各布信息的指针。只有在非卡尔蒂斯单元上才设置为有用的值。
   *
   */
  const Tensor<2, dim, VectorizedArrayType> *jacobian;

  /**
   * 指向当前单元格的雅各布行列式的指针。如果在笛卡尔单元或具有恒定雅各布系数的单元上，这只是雅各布行列式，否则就是雅各布行列式乘以正交权。
   *
   */
  const VectorizedArrayType *J_value;

  /**
   * 一个指向面的法向量的指针。
   *
   */
  const Tensor<1, dim, VectorizedArrayType> *normal_vectors;

  /**
   * 一个指向面的法向量乘以雅各布式的指针。
   *
   */
  const Tensor<1, dim, VectorizedArrayType> *normal_x_jacobian;

  /**
   * 一个指向基础正交公式的正交权重的指针。
   *
   */
  const Number *quadrature_weights;

  /**
   * 在调用reinit()后，存储我们当前正在处理的单元格的编号。
   *
   */
  unsigned int cell;

  /**
   * 根据定义的法线方向，保存一个面是内部还是外部面的信息的标志。
   * 不用于单元格。
   *
   */
  bool is_interior_face;

  /**
   * 存储FEFaceEvaluation对象当前所指向的索引（内部面、外部面、与单元格相关的数据）。
   *
   */
  internal::MatrixFreeFunctions::DoFInfo::DoFAccessIndex dof_access_index;

  /**
   * 在`is_face==true'的情况下，存储给定单元格内一个面的当前编号，使用`0'到`2*dim'之间的值。
   *
   */
  unsigned int face_no;

  /**
   * 存储给定的面相对于标准方向的方向，如果在标准方向，则为0。
   *
   */
  unsigned int face_orientation;

  /**
   * 存储给定面的子面索引。通常情况下，这个变量的值为
   * numbers::invalid_unsigned_int
   * ，以表示对整个面的整合，但如果当前的物理面有一个更精细的邻居，它就是一个子面，必须适当地缩放ShapeInfo中的条目。
   *
   */
  unsigned int subface_index;

  /**
   * 在调用reinit()后，存储我们当前正在处理的单元格的类型。有效值是
   * @p cartesian,  @p affine 和 @p general,
   * ，它们对MappingInfo中的雅各布变换的内部存储方式有不同的影响。
   *
   */
  internal::MatrixFreeFunctions::GeometryType cell_type;

  /**
   * 可以用各自的构造函数在空中生成FEValues的几何数据。
   *
   */
  std::shared_ptr<internal::MatrixFreeFunctions::
                    MappingDataOnTheFly<dim, Number, VectorizedArrayType>>
    mapped_geometry;

  // Make FEEvaluation objects friends for access to protected member
  // mapped_geometry.
  template <int, int, int, int, typename, typename>
  friend class FEEvaluation;
};



/**
 * 这是FEEvaluation类的基类。
 * 该类通常不需要在用户代码中调用，也没有任何公共构造函数。使用方法是通过FEEvaluation类来代替。它实现了一个reinit方法，用于设置指针，以便快速执行对正交点的操作，访问
 * FEEvaluationBase::read_dof_values(),  FEEvaluationBase::set_dof_values(),
 * 和 FEEvaluationBase::distribute_local_to_global()
 * 函数的向量函数，以及访问有限元函数的值和梯度的方法。它还继承了FEEvaluationBaseData类所提供的几何体访问功能。
 * 这个类有五个模板参数。
 * @tparam  该类所使用的dim尺寸
 * @tparam  n_components
 * 解决PDEs系统时的矢量分量的数量。如果同一个操作被应用于一个PDE的几个分量（例如，一个矢量拉普拉斯方程），它们可以通过一个调用同时应用（而且通常更有效）。
 * @tparam  数字 Number 格式，通常为  @p double  或  @p float  。
 * @tparam  is_face
 * 该类是用于单元积分器（正交维度与空间维度相同）还是用于面积分器（正交维度少一个）？
 * @tparam  VectorizedArrayType
 * 以矢量方式处理的数组类型，默认为 VectorizedArray<Number>。
 *
 *
 * @note  目前只有VectorizedArray<Number,
 * width>被支持作为VectorizedArrayType。
 *
 *
 *
 * @ingroup matrixfree
 *
 */
template <int dim,
          int n_components_,
          typename Number,
          bool is_face                 = false,
          typename VectorizedArrayType = VectorizedArray<Number>>
class FEEvaluationBase
  : public FEEvaluationBaseData<dim, Number, is_face, VectorizedArrayType>
{
public:
  using number_type = Number;
  using value_type  = Tensor<1, n_components_, VectorizedArrayType>;
  using gradient_type =
    Tensor<1, n_components_, Tensor<1, dim, VectorizedArrayType>>;
  static constexpr unsigned int dimension    = dim;
  static constexpr unsigned int n_components = n_components_;

  /**
   * @name  1：从向量中读取和写入向量
   *
   */
  //@{
  /**
   * 对于向量 @p src,
   * 读出当前单元格自由度上的数值，并在内部存储。与没有约束条件时的功能
   * DoFAccessor::get_interpolated_dof_values
   * 类似，但它也包括来自悬挂节点的约束条件，所以可以把它也看作是与
   * AffineConstraints::read_dof_values
   * 类似的功能。注意，如果启用了矢量化，几个单元的DoF值会被设置。
   * 如果向量上的某些约束条件是不均匀的，则使用函数read_dof_values_plain代替，并通过调用
   * AffineConstraints::distribute.
   * 为向量提供有用的数据，也是在受约束的位置，在线性系统的求解过程中访问向量条目时，临时解应该总是有均匀的约束，这种方法是正确的。
   * 如果给定的向量模板类是块向量（通过模板函数
   * 'IsBlockVector<VectorType>::value', 确定，该函数检查从
   * dealii::BlockVectorBase) 或 std::vector<VectorType> 或
   * std::vector<VectorType>, 派生的向量，该函数从块向量的索引
   * @p first_index. 开始读取 @p n_components 块 对于非块向量， @p
   * first_index 被忽略了。
   * @note
   * 如果这个类是在没有MatrixFree对象的情况下构建的，并且信息是通过
   * DoFHandler<dim>::cell_iterator,
   * 来获取的，那么这个类只使用一个单元，这个函数会提取给定单元上的基础分量的值。这个调用比通过MatrixFree对象完成的调用要慢，并且导致一个结构在基于这些值的评估例程中不能有效地使用矢量化（相反，
   * VectorizedArray::size() 相同的副本被工作）。
   *
   */
  template <typename VectorType>
  void
  read_dof_values(const VectorType &src, const unsigned int first_index = 0);

  /**
   * 对于矢量 @p src,
   * 读出当前单元格自由度上的值，并在内部存储。与函数
   * DoFAccessor::get_interpolated_dof_values. 的功能相似
   * 相对于read_dof_values函数，这个函数从向量中读出普通条目，而不考虑存储的约束。当约束条件已经通过先前调用
   * AffineConstraints::distribute
   * 分布在向量上时，这种访问方式是合适的。当要使用不均匀的约束时，这个函数也是必要的，因为MatrixFree只能处理均匀的约束。注意，如果启用了矢量化，几个单元的DoF值会被设置。
   * 如果给定的向量模板类是块向量（通过模板函数
   * 'IsBlockVector<VectorType>::value', 确定，该函数检查从
   * dealii::BlockVectorBase) 或 std::vector<VectorType> 或
   * std::vector<VectorType>, 派生的向量，该函数从块向量的索引
   * @p first_index. 开始读取 @p n_components 块 对于非块向量， @p
   * first_index 被忽略了。
   * @note
   * 如果这个类是在没有MatrixFree对象的情况下构建的，并且信息是通过
   * DoFHandler<dim>::cell_iterator,
   * 来获取的，那么这个类只使用一个单元，这个函数提取给定单元上的基础分量的值。这个调用比通过MatrixFree对象完成的调用要慢，并且导致一个结构在基于这些值的评估例程中不能有效地使用矢量化（相反，
   * VectorizedArray::size() 相同的副本被工作）。
   *
   */
  template <typename VectorType>
  void
  read_dof_values_plain(const VectorType & src,
                        const unsigned int first_index = 0);

  /**
   * 获取内部存储在当前单元格的dof值的值，并将它们加到向量中
   * @p dst.
   * 该函数在写操作过程中也应用约束。因此，其功能与函数
   * AffineConstraints::distribute_local_to_global.
   * 相似。如果启用了矢量化，则会使用几个单元的DoF值。
   * 如果给定的向量模板类是块向量（通过模板函数
   * 'IsBlockVector<VectorType>::value', 确定，该函数检查从
   * dealii::BlockVectorBase) 或 std::vector<VectorType> 或
   * std::vector<VectorType>, 派生的向量，该函数写入块向量的 @p
   * n_components 块，从索引 @p first_index. 开始 对于非块向量，
   * @p first_index  被忽略了。     @p mask
   * 可以用来抑制对当前单元向量批中包含的一些单元的写入访问，例如，在本地时间步进的情况下，一些单元被排除在调用之外。比特集中的
   * "true "值意味着将处理相应的车道索引，而 "false
   * "值则跳过该索引。默认设置是一个包含所有1的比特集，这将把累积的积分写到批次中的所有单元。
   * @note
   * 如果这个类是在没有MatrixFree对象的情况下构建的，并且信息是通过
   * DoFHandler<dim>::cell_iterator,
   * 来获取的，那么这个类只使用一个单元格，这个函数会提取给定单元格上的底层组件的值。这个调用比通过MatrixFree对象完成的调用要慢，并且导致一个结构在基于这些值的评估例程中不能有效地使用矢量化（相反，
   * VectorizedArray::size() 相同的副本被工作）。
   *
   */
  template <typename VectorType>
  void
  distribute_local_to_global(
    VectorType &                                    dst,
    const unsigned int                              first_index = 0,
    const std::bitset<VectorizedArrayType::size()> &mask =
      std::bitset<VectorizedArrayType::size()>().flip()) const;

  /**
   * 获取内部存储在当前单元格的自由度值，并将其写入向量
   * @p dst.
   * ，该函数跳过了被约束的自由度。与distribution_local_to_global方法相反，当前单元格给出的位置上的旧值被覆盖。因此，如果一个自由度与一个以上的单元相关联（在连续有限元中很常见），这些值将被覆盖，只有最后写入的值被保留。请注意，在并行环境下，这个函数也可能触及其他MPI进程所拥有的自由度，因此，随后的更新或积累幽灵值（如
   * MatrixFree::loop()
   * 所做的）可能会使这个函数设置的自由度失效。
   * 如果给定的向量模板类是块向量（通过模板函数
   * 'IsBlockVector<VectorType>::value', 确定，该函数检查从
   * dealii::BlockVectorBase) 或 std::vector<VectorType> 或
   * std::vector<VectorType>, 派生的向量，该函数从索引 @p
   * first_index. 开始向块向量的 @p n_components 块写入。
   * 对于非块向量， @p first_index 被忽略了。     @p mask
   * 可以用来抑制对当前单元向量批中包含的一些单元的写入访问，例如，在本地时间步进的情况下，一些单元被排除在调用之外。
   * 比特集中的 "true "值意味着将处理相应的车道索引，而
   * "false
   * "值则跳过该索引。默认设置是一个包含所有1的bitset，它将把累积的积分写到批次中的所有单元。
   * @note
   * 如果这个类是在没有MatrixFree对象的情况下构建的，并且信息是通过
   * DoFHandler<dim>::cell_iterator,
   * 来获取的，那么这个类只使用一个单一的单元格，这个函数会提取给定单元格上的底层组件的值。这个调用比通过MatrixFree对象完成的调用要慢，并且导致一个结构在基于这些值的评估例程中不能有效地使用矢量化（相反，
   * VectorizedArray::size() 相同的副本被工作）。
   *
   */
  template <typename VectorType>
  void
  set_dof_values(VectorType &       dst,
                 const unsigned int first_index = 0,
                 const std::bitset<VectorizedArrayType::size()> &mask =
                   std::bitset<VectorizedArrayType::size()>().flip()) const;

  /**
   * 与set_dof_values()相同，但不解决约束。
   *
   */
  template <typename VectorType>
  void
  set_dof_values_plain(
    VectorType &                                    dst,
    const unsigned int                              first_index = 0,
    const std::bitset<VectorizedArrayType::size()> &mask =
      std::bitset<VectorizedArrayType::size()>().flip()) const;

  //@}

  /**
   * @name  2：访问正交点的数据或聚集矢量数据
   *
   */
  //@{
  /**
   * 返回索引为 @p
   * dof的局部自由度的存储值。如果对象是矢量值的，就会给出一个矢量值的返回参数。因此，参数
   * @p dof 最多可以运行到 @p  dofs_per_component，而不是 @p
   * dofs_per_cell
   * ，因为矢量值FE的不同成分会一起返回。请注意，当矢量化被启用时，来自几个单元的值被分组在一起。如果
   * @p set_dof_values
   * 是最后被调用的，那么该值就对应于那里设置的值。如果
   * @p integrate
   * 是最后被调用的，那么它就对应于具有给定索引的测试函数的综合函数的值。
   * 请注意，派生类FEEvaluationAccess为标量情况（n_components ==
   * 1）和向量情况（n_components == dim）重载了这个操作。
   *
   */
  value_type
  get_dof_value(const unsigned int dof) const;

  /**
   * 向包含自由度分量的字段写入一个值  @p dof.  向通过 @p
   * get_dof_value访问的同一字段写入一个值。因此，一旦提交了一个值，从向量中读取的原始数据就会被覆盖。
   * 请注意，派生类FEEvaluationAccess为标量情况（n_components ==
   * 1）和向量值情况（n_components == dim）重载了这个操作。
   *
   */
  void
  submit_dof_value(const value_type val_in, const unsigned int dof);

  /**
   * 在调用 FEEvaluation::evaluate() 并设置了 EvaluationFlags::values
   * 后，返回正交点号 @p q_point
   * 处的有限元函数值，或者调用 FEEvaluationBase::submit_value().
   * 时已经存储在那里的值。
   * 如果对象是矢量值的，将给出一个矢量值的返回参数。注意，当矢量化被启用时，来自几个单元的值被分组在一起。
   * 请注意，派生类FEEvaluationAccess为标量情况（n_components ==
   * 1）和向量值情况（n_components ==
   * dim）重载了这一操作的特殊性。
   *
   */
  value_type
  get_value(const unsigned int q_point) const;

  /**
   * 向包含正交点上的值的字段写一个值，其成分为  @p
   * q_point.
   * 通过get_value()访问同一个字段。如果在调用设置了
   * EvaluationFlags::values 的函数 FEEvaluation::integrate()
   * 之前应用，这将指定由当前单元上的所有基函数测试并整合的值。
   * 请注意，派生类FEEvaluationAccess为标量情况（n_components ==
   * 1）和向量情况（n_components == dim）重载了这个操作。
   *
   */
  void
  submit_value(const value_type val_in, const unsigned int q_point);

  /**
   * 在用 EvaluationFlags::gradients, 调用 FEEvaluation::evaluate()
   * 后，返回正交点 @p q_point 的有限元函数梯度，或用
   * FEEvaluationBase::submit_gradient().
   * 调用后，返回存储在那里的值
   * 注意，派生类FEEvaluationAccess为标量情况（n_components ==
   * 1）和矢量情况（n_components ==
   * dim）重载了该操作，具有特殊性。
   *
   */
  gradient_type
  get_gradient(const unsigned int q_point) const;

  /**
   * 在调用 FEEvaluation::evaluate(EvaluationFlags::gradients)
   * 面的法线方向后，返回正交点号 @p q_point
   * 的有限元函数的导数。  $\boldsymbol \nabla u(\mathbf x_q) \cdot
   * \mathbf n(\mathbf x_q)$  这个调用等同于调用get_gradient()
   * get_normal_vector()，但将使用更有效的内部数据表示。
   * 请注意，派生类FEEvaluationAccess为标量情况（n_components ==
   * 1）和向量情况（n_components == dim）重载了这个操作。
   *
   */
  value_type
  get_normal_derivative(const unsigned int q_point) const;

  /**
   * 写一个贡献，这个贡献被梯度测试到包含分量为 @p
   * q_point.
   * 的正交点上的值的字段，通过get_gradient()访问同一个字段。如果在函数
   * FEEvaluation::integrate(EvaluationFlags::gradients)
   * 被调用之前应用，这指定了当前单元上所有基函数梯度测试的内容，并对其进行积分。
   * 请注意，派生类FEEvaluationAccess为标量情况（n_components ==
   * 1）和向量情况（n_components ==
   * dim）重载了这个操作，并进行了特殊化处理。
   *
   */
  void
  submit_gradient(const gradient_type grad_in, const unsigned int q_point);

  /**
   * 编写一个贡献，该贡献由梯度测试到包含分量为 @p
   * q_point的正交点上的值的域。与通过get_gradient()或get_normal_derivative()访问相同的字段。如果在函数
   * FEEvaluation::integrate(EvaluationFlags::gradients)
   * 被调用之前应用，这指定了当前单元上所有基函数梯度的测试内容，并在此基础上进行积分。
   * @note
   * 这个操作将数据写到与submit_gradient()相同的字段。因此，只能使用这两者中的一个。通常情况下，对这个函数的潜在调用的贡献必须加到submit_gradient()的贡献中。
   * @note 派生类FEEvaluationAccess为标量情况（n_components ==
   * 1）和矢量情况（n_components == dim）重载了这个操作。
   *
   */
  void
  submit_normal_derivative(const value_type   grad_in,
                           const unsigned int q_point);

  /**
   * 在调用 FEEvaluation::evaluate(EvaluationFlags::hessians).
   * 后，返回正交点号 @p q_point
   * 处的有限元函数的Hessian。如果只需要Hessian的对角线甚至跟踪，即拉普拉斯，请使用以下其他函数。
   * 请注意，派生类FEEvaluationAccess为标量情况（n_components ==
   * 1）和向量情况（n_components == dim）重载了这个操作。
   *
   */
  Tensor<1, n_components_, Tensor<2, dim, VectorizedArrayType>>
  get_hessian(const unsigned int q_point) const;

  /**
   * 在调用 FEEvaluation::evaluate(EvaluationFlags::hessians).
   * 后，返回正交点编号为 @p q_point
   * 的有限元函数的对角线。注意，派生类FEEvaluationAccess对该操作进行了重载，并对标量情况（n_components
   * == 1）和矢量情况（n_components ==
   * dim）进行了专业化处理。
   *
   */
  gradient_type
  get_hessian_diagonal(const unsigned int q_point) const;

  /**
   * 在调用 FEEvaluation::evaluate(EvaluationFlags::hessians).
   * 后，返回正交点号 @p q_point
   * 处的有限元函数的拉普拉斯（即Hessian的踪迹）。
   * 与计算全部Hessian的情况相比，当只要求拉普拉斯时，可以节省一些操作。
   * 请注意，派生类FEEvaluationAccess为标量情况（n_components ==
   * 1）和向量情况（n_components == dim）重载了这个操作。
   *
   */
  value_type
  get_laplacian(const unsigned int q_point) const;

#ifdef DOXYGEN
  // doxygen does not anyhow mention functions coming from partial template
  // specialization of the base class, in this case FEEvaluationAccess<dim,dim>.
  // For now, hack in those functions manually only to fix documentation:

  /**
   * 在调用 @p evaluate(...,true,...). 后，返回正交点号 @p q_point
   * 处的矢量值有限元的发散。
   * @note 仅对n_components_==dim有效。
   *
   */
  VectorizedArrayType
  get_divergence(const unsigned int q_point) const;

  /**
   * 在调用 @p evaluation(...,true,...)后，返回正交点号 @p q_point
   * 处的矢量值有限元的对称梯度。它对应于<tt>0.5
   * (grad+grad<sup>T</sup>)</tt>。
   * @note 只对n_components_==dim有效。
   *
   */
  SymmetricTensor<2, dim, VectorizedArrayType>
  get_symmetric_gradient(const unsigned int q_point) const;

  /**
   * 在调用 @p evaluate(...,true,...)后，返回矢量场的卷积，
   * $\nabla \times v$ 。
   * @note  只对n_components_==dim有效。
   *
   */
  Tensor<1, (dim == 2 ? 1 : dim), VectorizedArrayType>
  get_curl(const unsigned int q_point) const;

  /**
   * 写入一个贡献，这个贡献被发散测试到包含组件 @p
   * q_point. 的正交点上的值的字段，访问与通过 @p
   * get_gradient. 相同的字段，如果在调用函数 @p
   * integrate(...,true)
   * 之前应用，这指定了由当前单元上的所有基函数梯度测试并整合的内容。
   * @note  只对n_components_==dim有效。
   * @note
   * 该操作将数据写到与submit_gradient()相同的字段。因此，只能使用这两者中的一个。通常情况下，对这个函数的潜在调用的贡献必须加到submit_gradient()的贡献的对角线中。
   *
   */
  void
  submit_divergence(const VectorizedArrayType div_in,
                    const unsigned int        q_point);

  /**
   * 写入一个贡献，该贡献被对称梯度测试到包含正交点上的值的字段，其分量为
   * @p q_point.  通过 @p get_symmetric_gradient. 访问同一字段
   * 如果在函数 @p integrate(...,true)
   * 被调用之前应用，这指定了对称梯度，它被当前单元上的所有基函数对称梯度测试并整合到。
   * @note  只对n_components_==dim有效。
   * @note
   * 该操作将数据写到与submit_gradient()相同的字段。因此，只能使用这两者中的一个。通常情况下，对这个函数的潜在调用的贡献必须加到submit_gradient()的rank-2张量的相应条目中。
   *
   */
  void
  submit_symmetric_gradient(
    const SymmetricTensor<2, dim, VectorizedArrayType> grad_in,
    const unsigned int                                 q_point);

  /**
   * 写下包含正交点上的值的curl的分量  @p q_point.  通过 @p
   * get_gradient. 访问同一数据域
   * @note  只对n_components_==dim有效。
   * @note
   * 该操作将数据写到与submit_gradient()相同的字段。因此，只能使用这两者中的一个。通常情况下，对这个函数的潜在调用的贡献必须加到submit_gradient()的rank-2张量的相应条目中。
   *
   */
  void
  submit_curl(const Tensor<1, dim == 2 ? 1 : dim, VectorizedArrayType> curl_in,
              const unsigned int                                       q_point);

#endif

  /**
   * 取正交点的值，乘以雅各布行列式和正交权重（JxW），并对单元上所有正交点的值进行求和。其结果是一个标量，代表函数在单元上的积分。如果使用了一个矢量元素，结果的分量仍然是分开的。此外，如果启用了矢量化，几个单元的积分值将包含在返回的VectorizedArray字段的槽中。
   * @note
   * 如果FEEvaluation对象被初始化为一批单元，在SIMD向量VectorizedArray中并非所有的通道都代表实际数据，这个方法在虚拟数据（从最后一个有效通道复制的）上执行计算，将没有意义。因此，用户需要确保在任何计算中不明确地使用它，比如在对几个单元格的结果求和时。
   *
   */
  value_type
  integrate_value() const;

  //@}

  /**
   * @name  3：对内部数据的访问
   *
   */
  //@{
  /**
   * 返回一个指向dof值的第一个字段的只读指针。这是read_dof_values()函数写进的数据字段。首先是第一个组件的dof值，然后是第二个组件的所有值，以此类推。这与这个类中使用的内部数据结构有关。一般来说，使用get_dof_value()函数来代替比较安全。
   *
   */
  const VectorizedArrayType *
  begin_dof_values() const;

  /**
   * 返回一个指向dof值的第一个字段的读写指针。
   * 这是read_dof_values()函数写进的数据字段。首先是第一个组件的dof值，然后是第二个组件的所有值，以此类推。这与这个类中使用的内部数据结构有关。一般来说，使用get_dof_value()函数来代替比较安全。
   *
   */
  VectorizedArrayType *
  begin_dof_values();

  /**
   * 返回一个指向正交点上函数值第一域的只读指针。首先是第一个分量的所有正交点上的函数值，然后是第二个分量的所有数值，以此类推。这与本类中使用的内部数据结构有关。调用
   * @p evaluate
   * 后的原始数据只包含单元格操作，所以可能的变换、正交权重等必须手动应用。一般来说，使用get_value()函数反而更安全，它在内部完成所有的转换。
   *
   */
  const VectorizedArrayType *
  begin_values() const;

  /**
   * 返回一个指向正交点上的函数值的第一个字段的读和写指针。首先是第一个分量的所有正交点上的函数值，然后是第二个分量的所有数值，以此类推。这与本类中使用的内部数据结构有关。调用
   * @p evaluate
   * 后的原始数据只包含单元格操作，所以可能的变换、正交权重等必须手动应用。一般来说，使用get_value()函数反而更安全，它在内部完成所有的转换。
   *
   */
  VectorizedArrayType *
  begin_values();

  /**
   * 返回一个指向正交点上的函数梯度的第一个字段的只读指针。首先是所有正交点上第一个分量的梯度的x分量，然后是y分量，以此类推。接下来是第二个分量的x分量，以此类推。这与本类中使用的内部数据结构有关。调用
   * @p evaluate
   * 后的原始数据只包含单元格操作，所以可能的变换、正交权重等必须手动应用。一般来说，使用get_gradient()函数反而更安全，它在内部完成所有的变换。
   *
   */
  const VectorizedArrayType *
  begin_gradients() const;

  /**
   * 返回一个指向正交点上函数梯度第一域的读写指针。首先是所有正交点上第一个分量的梯度的x分量，然后是y分量，以此类推。接下来是第二个分量的x分量，以此类推。这与本类中使用的内部数据结构有关。调用
   * @p evaluate
   * 后的原始数据只包含单元格操作，所以可能的变换、正交权重等必须手动应用。一般来说，使用get_gradient()函数反而更安全，它在内部完成所有的变换。
   *
   */
  VectorizedArrayType *
  begin_gradients();

  /**
   * 返回一个只读指针，指向正交点上的函数 hessians
   * 的第一个字段。首先是所有正交点上第一个分量的
   * hessians 的 xx-分量，然后是
   * yy-分量，zz-分量（3D），然后是
   * xy-分量，以此类推。接下来是第二个分量的xx-分量，以此类推。这与本类中使用的内部数据结构有关。调用
   * @p evaluate
   * 后的原始数据只包含单元格操作，所以可能的变换、正交权重等必须手动应用。一般来说，使用get_laplacian()或get_hessian()函数来代替比较安全，它在内部完成所有的变换。
   *
   */
  const VectorizedArrayType *
  begin_hessians() const;

  /**
   * 返回一个读写指针，指向正交点上的函数hesians的第一个字段。首先是所有正交点上第一个分量的
   * hessians 的 xx-分量，然后是
   * yy-分量，zz-分量（3D），然后是
   * xy-分量，以此类推。接下来是第二个分量的xx-分量，以此类推。这与本类中使用的内部数据结构有关。调用
   * @p evaluate
   * 后的原始数据只包含单元格操作，所以可能的变换、正交权重等必须手动应用。一般来说，使用get_laplacian()或get_hessian()函数来代替比较安全，它在内部完成所有的变换。
   *
   */
  VectorizedArrayType *
  begin_hessians();

  //@}

  /**
   * 返回第一个选定的分量。
   *
   */
  unsigned int
  get_first_selected_component() const;

protected:
  /**
   * 构造函数。为了防止用户直接使用这个类，做了保护。采取所有存储在MatrixFree中的数据。如果应用于有多个有限元或多个正交公式的问题，在构造
   * @p matrix_free,   @p dof_no,   @p  first_selected_component和 @p quad_no
   * 时，允许选择适当的组件。
   *
   */
  FEEvaluationBase(
    const MatrixFree<dim, Number, VectorizedArrayType> &matrix_free,
    const unsigned int                                  dof_no,
    const unsigned int first_selected_component,
    const unsigned int quad_no,
    const unsigned int fe_degree,
    const unsigned int n_q_points,
    const bool         is_interior_face,
    const unsigned int active_fe_index,
    const unsigned int active_quad_index,
    const unsigned int face_type);

  /**
   * 构造函数，带有减少的功能，工作方式与FEValues类似。参数与传递给FEValues构造函数的参数类似，但明显的区别是，FEEvaluation期望一个一维的正交公式，Quadrature<1>，而不是
   * @p dim
   * 维的。有限元既可以是标量值，也可以是矢量值，但是这个方法每次总是只选择一个标量基元（按照类模板参数指定的
   * @p n_components 副本）。对于矢量值元素，可选的参数 @p
   * first_selected_component
   * 允许指定用于评估的基元素的索引。注意，内部数据结构总是假定基元是原始的，目前不支持非原始的。
   * 正如从FEValues中得知的那样，调用带有
   * Triangulation::cell_iterator
   * 的reinit方法是必要的，以使当前类的几何和自由度为人所知。如果迭代器包括DoFHandler信息（即它是
   * DoFHandler::cell_iterator 或类似的），初始化也允许以
   * DoFHandler::active_cell_iterator
   * 类型的标准方式一次从向量中读取或写入一个单元。然而，这种方法比使用MPI的MatrixFree的路径要慢得多，因为必须进行索引转换。由于每次只使用一个单元，这种方法不会在几个元素上进行矢量化（这对矢量操作来说是最有效的），而只可能在元素内进行，如果评估/整合例程在用户代码内被结合起来（例如计算单元矩阵）。
   * 可选的FEEvaluationBaseData对象允许几个FEEvaluation对象共享几何评估，也就是说，底层映射和正交点只需要评估一次。这只有在正交公式相同的情况下才有效。否则，将创建一个新的评估对象。当你打算将FEEvaluation对象与另一个对象并行使用时，请确保不要传递一个可选的对象，因为否则打算共享的对象可能会产生竞赛条件。
   *
   */
  FEEvaluationBase(
    const Mapping<dim> &      mapping,
    const FiniteElement<dim> &fe,
    const Quadrature<1> &     quadrature,
    const UpdateFlags         update_flags,
    const unsigned int        first_selected_component,
    const FEEvaluationBaseData<dim, Number, is_face, VectorizedArrayType>
      *other);

  /**
   * 复制构造函数。如果FEEvaluationBase是由映射、fe、正交和更新标志构建的，基于FEValues的底层几何评估将被深度复制，以允许与线程并行使用。
   *
   */
  FEEvaluationBase(const FEEvaluationBase &other);

  /**
   * 复制赋值运算符。如果FEEvaluationBase是由映射、fe、正交和更新标志构建的，基于FEValues的底层几何评估将被深度复制，以便允许与线程并行使用。
   *
   */
  FEEvaluationBase &
  operator=(const FEEvaluationBase &other);

  /**
   * 一个统一的函数，根据给定的模板操作从向量中读出和写入向量。它可以对
   * @p read_dof_values,  @p distribute_local_to_global, 和 @p
   * set_dof_values. 进行操作，一次对几个向量进行操作。
   *
   */
  template <typename VectorType, typename VectorOperation>
  void
  read_write_operation(
    const VectorOperation &                        operation,
    const std::array<VectorType *, n_components_> &vectors,
    const std::array<
      const std::vector<ArrayView<const typename VectorType::value_type>> *,
      n_components_> &                              vectors_sm,
    const std::bitset<VectorizedArrayType::size()> &mask,
    const bool apply_constraints = true) const;

  /**
   * 一个统一的函数，基于给定的模板操作从向量中读出和写入向量，用于DG型方案，其中单元格上的所有自由度是连续的。它可以一次对多个向量进行read_dof_values()、distribut_local_to_global()和set_dof_values()的操作，具体取决于n_components。
   *
   */
  template <typename VectorType, typename VectorOperation>
  void
  read_write_operation_contiguous(
    const VectorOperation &                        operation,
    const std::array<VectorType *, n_components_> &vectors,
    const std::array<
      const std::vector<ArrayView<const typename VectorType::value_type>> *,
      n_components_> &                              vectors_sm,
    const std::bitset<VectorizedArrayType::size()> &mask) const;

  /**
   * 一个统一的函数，在我们没有底层MatrixFree对象的情况下，根据给定的模板操作从向量中读取和写入向量。它可以对
   * @p read_dof_values,  @p distribute_local_to_global, 和 @p
   * set_dof_values. 进行操作。
   * 它一次对几个向量进行操作，取决于n_components。
   *
   */
  template <typename VectorType, typename VectorOperation>
  void
  read_write_operation_global(
    const VectorOperation &                        operation,
    const std::array<VectorType *, n_components_> &vectors) const;

  /**
   * 这个字段存储了局部自由度的值（例如，从矢量中读出后，但在应用单元格变换前或在将它们分配到结果矢量中之前）。get_dof_value()和submit_dof_value()方法从这个字段读取或写入。
   * 这个数组的值存储在 @p scratch_data_array.
   * 的起始部分。由于其作为线程本地内存的访问，该内存可以在不同的调用之间得到重复使用。相对于在堆栈上请求内存，这种方法允许非常大的多项式程度。
   *
   */
  VectorizedArrayType *values_dofs[n_components];

  /**
   * 这个字段存储了应用单元格变换后或积分前正交点上的有限元函数的值。
   * get_value()和submit_value()方法访问这个字段。
   * 这个数组的值存储在 @p scratch_data_array.
   * 的起始部分。由于其作为线程本地内存的访问，内存可以在不同的调用之间得到重复使用。相对于在堆栈上请求内存，这种方法允许非常大的多项式程度。
   *
   */
  VectorizedArrayType *values_quad;

  /**
   * 这个字段存储了应用单元格变换后或积分前正交点上的有限元函数的梯度。get_gradient()和submit_gradient()方法（以及一些特殊的方法如get_symmetric_gradient()或get_divergence()）访问这个字段。
   * 这个数组的值存储在 @p scratch_data_array.
   * 的起始部分。由于它作为线程本地内存的访问，内存可以在不同的调用之间得到重复使用。相对于在堆栈上请求内存，这种方法允许非常大的多项式程度。
   *
   */
  VectorizedArrayType *gradients_quad;

  /**
   * 这个字段存储了应用单元格变换后正交点上的有限元函数的Hessians。get_hessian(),
   * get_laplacian(), get_hessian_diagonal()方法访问这个字段。
   * 这个数组的值存储在 @p scratch_data_array.
   * 的起始部分。由于其作为线程本地内存的访问，该内存可以在不同的调用之间得到重复使用。相对于在堆栈上请求内存，这种方法允许非常大的多项式程度。
   *
   */
  VectorizedArrayType *hessians_quad;

  /**
   * 存储在MatrixFree存储类中检测到的有限元中的组件数量，以便与模板参数进行比较。
   *
   */
  const unsigned int n_fe_components;

  /**
   * 调试信息，跟踪dof值在访问前是否已被初始化。当使用未初始化的数据时，用于控制异常。
   *
   */
  bool dof_values_initialized;

  /**
   * 调试信息，跟踪正交点的值在访问前是否已经被初始化。用于控制使用未初始化数据时的异常情况。
   *
   */
  bool values_quad_initialized;

  /**
   * 调试信息，跟踪正交点上的梯度是否在访问前被初始化。用于控制使用未初始化数据时的异常情况。
   *
   */
  bool gradients_quad_initialized;

  /**
   * 调试信息，跟踪正交点上的Hessians在访问前是否已经被初始化。用于控制使用未初始化数据时的异常情况。
   *
   */
  bool hessians_quad_initialized;

  /**
   * 调试信息跟踪正交点上的值是否在实际盯住积分之前被提交给了积分。用于控制使用未初始化数据时的异常情况。
   *
   */
  bool values_quad_submitted;

  /**
   * 调试信息，跟踪正交点的梯度在积分实际被盯住之前是否已被提交用于积分。
   * 用于控制使用未初始化数据时的异常情况。
   *
   */
  bool gradients_quad_submitted;

  /**
   * 对于一个有多个基元的FiniteElement，选择这个数据结构应该从哪个分量开始。
   *
   */
  const unsigned int first_selected_component;

  /**
   * 当初始化时没有给出MatrixFree对象时，需要一个临时数据结构来读取自由度。
   *
   */
  mutable std::vector<types::global_dof_index> local_dof_indices;

private:
  /**
   * 将数值、梯度、赫西恩的指针设置为基类的中央scratch_data_array。
   *
   */
  void
  set_data_pointers();
};



/**
 * 这个类提供对FEEvaluation类的数据字段的访问。通用访问是通过基类实现的，标量和矢量值元素的特殊化是单独定义的。
 *
 *
 * @ingroup matrixfree
 *
 *
 */
template <int dim,
          int n_components_,
          typename Number,
          bool is_face,
          typename VectorizedArrayType = VectorizedArray<Number>>
class FEEvaluationAccess : public FEEvaluationBase<dim,
                                                   n_components_,
                                                   Number,
                                                   is_face,
                                                   VectorizedArrayType>
{
  static_assert(
    std::is_same<Number, typename VectorizedArrayType::value_type>::value,
    "Type of Number and of VectorizedArrayType do not match.");

public:
  using number_type = Number;
  using value_type  = Tensor<1, n_components_, VectorizedArrayType>;
  using gradient_type =
    Tensor<1, n_components_, Tensor<1, dim, VectorizedArrayType>>;
  static constexpr unsigned int dimension    = dim;
  static constexpr unsigned int n_components = n_components_;
  using BaseClass =
    FEEvaluationBase<dim, n_components_, Number, is_face, VectorizedArrayType>;

protected:
  /**
   * 构造函数。为了防止在用户代码中进行初始化，该构造函数被保护起来。接受存储在MatrixFree中的所有数据。如果应用于有多个有限元的问题，或者在构造
   * @p matrix_free,  @p first_selected_component 和 @p
   * quad_no时选择了多个正交公式，允许选择适当的组件。
   *
   */
  FEEvaluationAccess(
    const MatrixFree<dim, Number, VectorizedArrayType> &matrix_free,
    const unsigned int                                  dof_no,
    const unsigned int first_selected_component,
    const unsigned int quad_no,
    const unsigned int fe_degree,
    const unsigned int n_q_points,
    const bool         is_interior_face  = true,
    const unsigned int active_fe_index   = numbers::invalid_unsigned_int,
    const unsigned int active_quad_index = numbers::invalid_unsigned_int,
    const unsigned int face_type         = numbers::invalid_unsigned_int);

  /**
   * 构造函数的功能减少，用于类似于FEValues的FEEvaluation的用法，包括矩阵组装。
   *
   */
  FEEvaluationAccess(
    const Mapping<dim> &      mapping,
    const FiniteElement<dim> &fe,
    const Quadrature<1> &     quadrature,
    const UpdateFlags         update_flags,
    const unsigned int        first_selected_component,
    const FEEvaluationBaseData<dim, Number, is_face, VectorizedArrayType>
      *other);

  /**
   * 复制构造函数
   *
   */
  FEEvaluationAccess(const FEEvaluationAccess &other);

  /**
   * 拷贝赋值操作符
   *
   */
  FEEvaluationAccess &
  operator=(const FEEvaluationAccess &other);
};



/**
 * 该类提供对FEEvaluation类的数据字段的访问。标量字段的部分特殊化，定义了对简单数据字段的访问，即标量的值和Tensor<1,dim>的梯度。
 *
 *
 * @ingroup matrixfree
 *
 */
template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
class FEEvaluationAccess<dim, 1, Number, is_face, VectorizedArrayType>
  : public FEEvaluationBase<dim, 1, Number, is_face, VectorizedArrayType>
{
  static_assert(
    std::is_same<Number, typename VectorizedArrayType::value_type>::value,
    "Type of Number and of VectorizedArrayType do not match.");

public:
  using number_type                       = Number;
  using value_type                        = VectorizedArrayType;
  using gradient_type                     = Tensor<1, dim, VectorizedArrayType>;
  static constexpr unsigned int dimension = dim;
  using BaseClass =
    FEEvaluationBase<dim, 1, Number, is_face, VectorizedArrayType>;

  /**
   * @copydoc   FEEvaluationBase<dim,1,Number,is_face>::get_dof_value()
   *
   */
  value_type
  get_dof_value(const unsigned int dof) const;

  /**
   * @copydoc   FEEvaluationBase<dim,1,Number,is_face>::submit_dof_value()
   * FEEvaluationBase<dim,1,Number,is_face>::submit_dof_value() .
   *
   */
  void
  submit_dof_value(const value_type val_in, const unsigned int dof);

  /**
   * @copydoc   FEEvaluationBase<dim,1,Number,is_face>::get_value()
   * FEEvaluationBase<dim,1,Number,is_face>::get_value() 。
   *
   */
  value_type
  get_value(const unsigned int q_point) const;

  /**
   * @copydoc   FEEvaluationBase<dim,1,Number,is_face>::submit_value()
   * FEEvaluationBase<dim,1,Number,is_face>::submit_value() 。
   *
   */
  void
  submit_value(const value_type val_in, const unsigned int q_point);

  /**
   * @copydoc   FEEvaluationBase<dim,1,Number,is_face>::submit_value()
   * FEEvaluationBase<dim,1,Number,is_face>::submit_value() 。
   *
   */
  void
  submit_value(const Tensor<1, 1, VectorizedArrayType> val_in,
               const unsigned int                      q_point);

  /**
   * @copydoc   FEEvaluationBase<dim,1,Number,is_face>::get_gradient()
   * FEEvaluationBase<dim,1,Number,is_face>::get_gradient() .
   *
   */
  gradient_type
  get_gradient(const unsigned int q_point) const;

  /**
   * @copydoc
   * FEEvaluationBase<dim,1,Number,is_face>::get_normal_derivative()
   * FEEvaluationBase<dim,1,Number,is_face>::get_normal_derivative() 。
   *
   */
  value_type
  get_normal_derivative(const unsigned int q_point) const;

  /**
   * @copydoc   FEEvaluationBase<dim,1,Number,is_face>::submit_gradient()
   *
   */
  void
  submit_gradient(const gradient_type grad_in, const unsigned int q_point);

  /**
   * @copydoc
   * FEEvaluationBase<dim,1,Number,is_face>::submit_normal_derivative()
   * FEEvaluationBase<dim,1,Number,is_face>::submit_normal_derivative() 。
   *
   */
  void
  submit_normal_derivative(const value_type   grad_in,
                           const unsigned int q_point);

  /**
   * @copydoc   FEEvaluationBase<dim,1,Number,is_face>::get_hessian()
   * FEEvaluationBase<dim,1,Number,is_face>::get_hessian() .
   *
   */
  Tensor<2, dim, VectorizedArrayType>
  get_hessian(unsigned int q_point) const;

  /**
   * @copydoc   FEEvaluationBase<dim,1,Number,is_face>::get_hessian_diagonal()
   *
   */
  gradient_type
  get_hessian_diagonal(const unsigned int q_point) const;

  /**
   * @copydoc   FEEvaluationBase<dim,1,Number,is_face>::get_laplacian()
   * FEEvaluationBase<dim,1,Number,is_face>::get_laplacian() .
   *
   */
  value_type
  get_laplacian(const unsigned int q_point) const;

  /**
   * @copydoc   FEEvaluationBase<dim,1,Number,is_face>::integrate_value() .
   *
   */
  value_type
  integrate_value() const;

protected:
  /**
   * 构造函数。为了避免在用户代码中进行初始化，做了保护。接受所有存储在MatrixFree中的数据。如果应用于有多个有限元的问题，或在构造
   * @p matrix_free,   @p first_selected_component  和  @p
   * quad_no时选择了多个正交公式，则允许选择适当的组件。
   *
   */
  FEEvaluationAccess(
    const MatrixFree<dim, Number, VectorizedArrayType> &matrix_free,
    const unsigned int                                  dof_no,
    const unsigned int first_selected_component,
    const unsigned int quad_no,
    const unsigned int fe_degree,
    const unsigned int n_q_points,
    const bool         is_interior_face  = true,
    const unsigned int active_fe_index   = numbers::invalid_unsigned_int,
    const unsigned int active_quad_index = numbers::invalid_unsigned_int,
    const unsigned int face_type         = numbers::invalid_unsigned_int);

  /**
   * 构造函数的功能减少，用于类似于FEValues的FEEvaluation的用法，包括矩阵组装。
   *
   */
  FEEvaluationAccess(
    const Mapping<dim> &      mapping,
    const FiniteElement<dim> &fe,
    const Quadrature<1> &     quadrature,
    const UpdateFlags         update_flags,
    const unsigned int        first_selected_component,
    const FEEvaluationBaseData<dim, Number, is_face, VectorizedArrayType>
      *other);

  /**
   * 复制构造函数
   *
   */
  FEEvaluationAccess(const FEEvaluationAccess &other);

  /**
   * 拷贝赋值操作符
   *
   */
  FEEvaluationAccess &
  operator=(const FEEvaluationAccess &other);
};



/**
 * 该类提供对FEEvaluation类的数据字段的访问。对具有与基础空间维度一样多的分量的字段进行了部分专业化处理，即值是Tensor<1,dim>类型，梯度是Tensor<2,dim>类型。提供一些额外的访问函数，如对称梯度和发散。
 *
 *
 * @ingroup matrixfree
 *
 *
 */
template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
class FEEvaluationAccess<dim, dim, Number, is_face, VectorizedArrayType>
  : public FEEvaluationBase<dim, dim, Number, is_face, VectorizedArrayType>
{
  static_assert(
    std::is_same<Number, typename VectorizedArrayType::value_type>::value,
    "Type of Number and of VectorizedArrayType do not match.");

public:
  using number_type                       = Number;
  using value_type                        = Tensor<1, dim, VectorizedArrayType>;
  using gradient_type                     = Tensor<2, dim, VectorizedArrayType>;
  static constexpr unsigned int dimension = dim;
  static constexpr unsigned int n_components = dim;
  using BaseClass =
    FEEvaluationBase<dim, dim, Number, is_face, VectorizedArrayType>;

  /**
   * @copydoc   FEEvaluationBase<dim,dim,Number,is_face>::get_gradient() .
   *
   */
  gradient_type
  get_gradient(const unsigned int q_point) const;

  /**
   * 在调用 @p evaluate(...,true,...). 后，返回正交点号 @p q_point
   * 的矢量值有限元的发散。
   *
   */
  VectorizedArrayType
  get_divergence(const unsigned int q_point) const;

  /**
   * 在调用 @p evaluation(...,true,...)后，返回在正交点号 @p
   * q_point 处的矢量值有限元的对称梯度。它对应于<tt>0.5
   * (grad+grad<sup>T</sup>)</tt>。
   *
   */
  SymmetricTensor<2, dim, VectorizedArrayType>
  get_symmetric_gradient(const unsigned int q_point) const;

  /**
   * 在调用 @p evaluate(...,true,...)之后，返回矢量场的卷积，
   * $\nabla \times v$ 。
   *
   */
  Tensor<1, (dim == 2 ? 1 : dim), VectorizedArrayType>
  get_curl(const unsigned int q_point) const;

  /**
   * @copydoc   FEEvaluationBase<dim,dim,Number,is_face>::get_hessian() 。
   *
   */
  Tensor<3, dim, VectorizedArrayType>
  get_hessian(const unsigned int q_point) const;

  /**
   * @copydoc
   * FEEvaluationBase<dim,dim,Number,is_face>::get_hessian_diagonal()
   * FEEvaluationBase<dim,dim,Number,is_face>::get_hessian_diagonal() .
   *
   */
  gradient_type
  get_hessian_diagonal(const unsigned int q_point) const;

  /**
   * @copydoc   FEEvaluationBase<dim,dim,Number,is_face>::submit_gradient() .
   *
   */
  void
  submit_gradient(const gradient_type grad_in, const unsigned int q_point);

  /**
   * 编写一个贡献，通过梯度对包含分量为 @p q_point.
   * 的正交点上的值的域进行测试，这个函数是其他submit_gradient函数的替代方案，当使用固定数量的方程组时，恰好与某些维度的维度重合，但不是全部。为了实现与维度无关的编程，可以使用这个函数来代替。
   *
   */
  void
  submit_gradient(
    const Tensor<1, dim, Tensor<1, dim, VectorizedArrayType>> grad_in,
    const unsigned int                                        q_point);

  /**
   * 编写一个贡献，通过发散测试到包含正交点上的值的字段，组件
   * @p q_point.  通过 @p get_gradient.
   * 访问相同的字段，如果在调用函数 @p integrate(...,true)
   * 之前应用，这指定了由当前单元上的所有基函数梯度测试并整合的内容。
   *
   */
  void
  submit_divergence(const VectorizedArrayType div_in,
                    const unsigned int        q_point);

  /**
   * 写入一个贡献，这个贡献被对称梯度测试到包含正交点上的值的字段，其分量为
   * @p q_point.  通过  @p get_symmetric_gradient.  访问同一个字段
   * 如果在调用函数  @p integrate(...,true)
   * 之前应用，这指定了对称梯度，它被当前单元上的所有基函数对称梯度测试并整合在一起。
   *
   */
  void
  submit_symmetric_gradient(
    const SymmetricTensor<2, dim, VectorizedArrayType> grad_in,
    const unsigned int                                 q_point);

  /**
   * 写下包含正交点 @p q_point. 上的值的curl的分量 通过 @p
   * get_gradient. 访问相同的数据域。
   *
   */
  void
  submit_curl(const Tensor<1, dim == 2 ? 1 : dim, VectorizedArrayType> curl_in,
              const unsigned int                                       q_point);

protected:
  /**
   * 构造函数。为了避免在用户代码中进行初始化，做了保护。取用存储在MatrixFree中的所有数据。如果应用于有多个有限元的问题，或在构建
   * @p matrix_free,  @p first_selected_component 和 @p
   * quad_no时选择了多个正交公式，允许选择适当的组件。
   *
   */
  FEEvaluationAccess(
    const MatrixFree<dim, Number, VectorizedArrayType> &matrix_free,
    const unsigned int                                  dof_no,
    const unsigned int first_selected_component,
    const unsigned int quad_no,
    const unsigned int dofs_per_cell,
    const unsigned int n_q_points,
    const bool         is_interior_face  = true,
    const unsigned int active_fe_index   = numbers::invalid_unsigned_int,
    const unsigned int active_quad_index = numbers::invalid_unsigned_int,
    const unsigned int face_type         = numbers::invalid_unsigned_int);

  /**
   * 构造函数的功能减少，用于类似于FEValues的FEEvaluation的用法，包括矩阵组装。
   *
   */
  FEEvaluationAccess(
    const Mapping<dim> &      mapping,
    const FiniteElement<dim> &fe,
    const Quadrature<1> &     quadrature,
    const UpdateFlags         update_flags,
    const unsigned int        first_selected_component,
    const FEEvaluationBaseData<dim, Number, is_face, VectorizedArrayType>
      *other);

  /**
   * 复制构造函数
   *
   */
  FEEvaluationAccess(const FEEvaluationAccess &other);

  /**
   * 拷贝赋值操作符
   *
   */
  FEEvaluationAccess &
  operator=(const FEEvaluationAccess &other);
};


/**
 * 该类提供对FEEvaluation类的数据字段的访问。1d中标量字段的部分特殊化，定义了对简单数据字段的访问，即标量为值，Tensor<1,1>为梯度。
 *
 *
 * @ingroup matrixfree
 *
 *
 */
template <typename Number, bool is_face, typename VectorizedArrayType>
class FEEvaluationAccess<1, 1, Number, is_face, VectorizedArrayType>
  : public FEEvaluationBase<1, 1, Number, is_face, VectorizedArrayType>
{
  static_assert(
    std::is_same<Number, typename VectorizedArrayType::value_type>::value,
    "Type of Number and of VectorizedArrayType do not match.");

public:
  using number_type                       = Number;
  using value_type                        = VectorizedArrayType;
  using gradient_type                     = Tensor<1, 1, VectorizedArrayType>;
  static constexpr unsigned int dimension = 1;
  using BaseClass =
    FEEvaluationBase<1, 1, Number, is_face, VectorizedArrayType>;

  /**
   * @copydoc   FEEvaluationBase<1,1,Number,is_face>::get_dof_value()
   * FEEvaluationBase<1,1,Number,is_face>::get_dof_value() 。
   *
   */
  value_type
  get_dof_value(const unsigned int dof) const;

  /**
   * @copydoc   FEEvaluationBase<1,1,Number,is_face>::submit_dof_value()
   * FEEvaluationBase<1,1,Number,is_face>::submit_dof_value() .
   *
   */
  void
  submit_dof_value(const value_type val_in, const unsigned int dof);

  /**
   * @copydoc   FEEvaluationBase<1,1,Number,is_face>::get_value() .
   *
   */
  value_type
  get_value(const unsigned int q_point) const;

  /**
   * @copydoc   FEEvaluationBase<1,1,Number,is_face>::submit_value() .
   *
   */
  void
  submit_value(const value_type val_in, const unsigned int q_point);

  /**
   * @copydoc   FEEvaluationBase<1,1,Number,is_face>::submit_value() .
   *
   */
  void
  submit_value(const gradient_type val_in, const unsigned int q_point);

  /**
   * @copydoc   FEEvaluationBase<1,1,Number,is_face>::get_gradient() .
   *
   */
  gradient_type
  get_gradient(const unsigned int q_point) const;

  /**
   * @copydoc   FEEvaluationBase<1,1,Number,is_face>::get_divergence() .
   *
   */
  value_type
  get_divergence(const unsigned int q_point) const;

  /**
   * @copydoc
   * FEEvaluationBase<dim,1,Number,is_face>::get_normal_derivative() .
   *
   */
  value_type
  get_normal_derivative(const unsigned int q_point) const;

  /**
   * @copydoc   FEEvaluationBase<1,1,Number,is_face>::submit_gradient()
   * FEEvaluationBase<1,1,Number,is_face>::submit_gradient() 。
   *
   */
  void
  submit_gradient(const gradient_type grad_in, const unsigned int q_point);

  /**
   * @copydoc   FEEvaluationBase<1,1,Number,is_face>::submit_gradient() .
   *
   */
  void
  submit_gradient(const value_type grad_in, const unsigned int q_point);

  /**
   * @copydoc   FEEvaluationBase<1,1,Number,is_face>::submit_gradient()
   *
   */
  void
  submit_gradient(const Tensor<2, 1, VectorizedArrayType> grad_in,
                  const unsigned int                      q_point);

  /**
   * @copydoc
   * FEEvaluationBase<1,1,Number,is_face>::submit_normal_derivative()
   * FEEvaluationBase<1,1,Number,is_face>::submit_normal_derivative() 。
   *
   */
  void
  submit_normal_derivative(const value_type   grad_in,
                           const unsigned int q_point);

  /**
   * @copydoc
   * FEEvaluationBase<1,1,Number,is_face>::submit_normal_derivative()
   * FEEvaluationBase<1,1,Number,is_face>::submit_normal_derivative() .
   *
   */
  void
  submit_normal_derivative(const gradient_type grad_in,
                           const unsigned int  q_point);

  /**
   * @copydoc   FEEvaluationBase<1,1,Number,is_face>::get_hessian()
   * FEEvaluationBase<1,1,Number,is_face>::get_hessian() 。
   *
   */
  Tensor<2, 1, VectorizedArrayType>
  get_hessian(unsigned int q_point) const;

  /**
   * @copydoc   FEEvaluationBase<1,1,Number,is_face>::get_hessian_diagonal()
   * FEEvaluationBase<1,1,Number,is_face>::get_hessian_diagonal() 。
   *
   */
  gradient_type
  get_hessian_diagonal(const unsigned int q_point) const;

  /**
   * @copydoc   FEEvaluationBase<1,1,Number,is_face>::get_laplacian()
   * FEEvaluationBase<1,1,Number,is_face>::get_laplacian() .
   *
   */
  value_type
  get_laplacian(const unsigned int q_point) const;

  /**
   * @copydoc   FEEvaluationBase<1,1,Number,is_face>::integrate_value()
   * FEEvaluationBase<1,1,Number,is_face>::integrate_value() 。
   *
   */
  value_type
  integrate_value() const;

protected:
  /**
   * 构造函数。为了避免在用户代码中进行初始化，做了保护。接受所有存储在MatrixFree中的数据。如果应用于有多个有限元的问题，或在构造
   * @p matrix_free,  @p first_selected_component 和 @p
   * quad_no时选择了多个正交公式，允许选择适当的组件。
   *
   */
  FEEvaluationAccess(
    const MatrixFree<1, Number, VectorizedArrayType> &matrix_free,
    const unsigned int                                dof_no,
    const unsigned int                                first_selected_component,
    const unsigned int                                quad_no,
    const unsigned int                                fe_degree,
    const unsigned int                                n_q_points,
    const bool                                        is_interior_face = true,
    const unsigned int active_fe_index   = numbers::invalid_unsigned_int,
    const unsigned int active_quad_index = numbers::invalid_unsigned_int,
    const unsigned int face_type         = numbers::invalid_unsigned_int);

  /**
   * 构造器的功能减少，用于类似FEValues的FEEvaluation的使用，包括矩阵组装。
   *
   */
  FEEvaluationAccess(
    const Mapping<1> &      mapping,
    const FiniteElement<1> &fe,
    const Quadrature<1> &   quadrature,
    const UpdateFlags       update_flags,
    const unsigned int      first_selected_component,
    const FEEvaluationBaseData<1, Number, is_face, VectorizedArrayType> *other);

  /**
   * 复制构造函数
   *
   */
  FEEvaluationAccess(const FEEvaluationAccess &other);

  /**
   * 拷贝赋值操作符
   *
   */
  FEEvaluationAccess &
  operator=(const FEEvaluationAccess &other);
};



/**
 * 该类提供了在正交点和单元格积分上评估函数所需的所有功能。在功能上，这个类与FEValues类似，但是，它包括很多专门的函数，使其速度更快（5到500之间，取决于多项式的程度）。对于DG中人脸项的评估，请参见类FEFaceEvaluation。
 * <h3>Usage and initialization</h3> <h4>Fast usage in combination with
 * MatrixFree</h4>
 * 首要的使用方法是通过MatrixFree对象初始化这个类，该对象缓存了所有与自由度和映射信息相关的内容。这样，就有可能使用矢量化的方式，一次为几个单元应用微分算子。
 * FEEvaluation的能力涵盖了大量的弱形式的积分任务。一般来说，有两类任务是可以完成的。一个是
 * @p evaluate 路径，从解向量插值到正交点。
 *
 *
 * @code
 * FEEvaluation<dim,fe_degree> phi(matrix_free);
 * for (unsigned int cell_index = cell_range.first;
 *    cell_index < cell_range.second; ++cell_index)
 * {
 *   phi.reinit(cell_index);
 *   phi.read_dof_values(vector);
 *   phi.evaluate(EvaluationFlags::values);   // interpolate values only
 *   for (unsigned int q=0; q<phi.n_q_points; ++q)
 *     {
 *       VectorizedArray<double> val = phi.get_value(q);
 *       // do something with val
 *     }
 * }
 * @endcode
 *
 * 同样，由 @p 矢量代表的有限元解的梯度可以通过 @p
 * phi.get_gradient(q)插值到正交点。read_dof_values()、evaluate()和get_value()的组合与
 * FEValues::get_function_values 或 FEValues::get_function_gradients
 * 所做的类似，但一般来说要快得多，因为它利用了张量积，见下面对评估例程的描述，并且可以通过矢量化一次为几个单元做这个操作。
 * FEEvaluation完成的第二类任务是右手边的积分任务。在有限元计算中，这些任务通常包括将正交点上的一个量（一个函数值，或者一个由有限元空间本身插值的场）与一组测试函数相乘，通过对每个正交点的值进行求和，再乘以正交权重和变换的雅各布行列式，对单元进行积分。如果给定一个通用的Function对象，我们想计算
 * $v_i = \int_\Omega \varphi_i f dx$
 * ，这可以通过以下单元的积分来完成。
 *
 *
 * @code
 * FEEvaluation<dim,fe_degree> phi(matrix_free);
 * Function<dim> &function = ...;
 * for (unsigned int cell_index = cell_range.first;
 *    cell_index < cell_range.second; ++cell_index)
 * {
 *   phi.reinit(cell_index);
 *   for (unsigned int q=0; q<phi.n_q_points; ++q)
 *     {
 *       Point<dim,VectorizedArray<double> > p_vect =
 *         phi.quadrature_point(q);
 *       // Need to evaluate function for each component in VectorizedArray
 *       VectorizedArray<double> f_value;
 *       for (unsigned int v=0; v<VectorizedArray<double>::size(); ++v)
 *         {
 *           Point<dim> p;
 *           for (unsigned int d=0; d<dim; ++d)
 *             p[d] = p_vect[d][v];
 *           f_value[v] = function.value(p);
 *         }
 *       phi.submit_value(f_value, q);
 *     }
 *   phi.integrate(EvaluationFlags::values);
 *   phi.distribute_local_to_global(dst);
 * }
 * @endcode
 *
 * 在这段代码中，对 @p phi.submit_value()
 * 的调用在实际积分之前为测试函数的乘法做了准备（在提交调用中，要测试的值也被乘以雅各布的行列式和正交的权重）。在
 * @p integrate()
 * 调用中，由FEEvaluation对象的每个基础函数（例如二维的FE_Q
 * @<2@>(1)
 * 的四个线性形状函数）测试的积分贡献被计算出来，这就给出了要加到
 * @p dst
 * 向量中的向量条目。需要注意的是，上面的代码需要明确地在向量数组中循环计算函数的分量，这对于与具有双参数的通用Function对象对接是必要的。简单的函数也可以直接用VectorizedArray形式实现，因为VectorizedArray提供了基本的数学运算。
 * 对于评估一个双线性形式，在源向量上的评估与涉及测试函数的积分相结合，被写入结果向量中。这种设置是无矩阵算子求值的背景，在
 * step-37 和 step-48 的教程程序中解释过。 请注意，通过
 * FEEvaluation::read_dof_values 和 FEEvaluation::distribute_local_to_global
 * 的两个向量访问，基于 MatrixFree::reinit()
 * 调用时指定的AffineConstraints对象，在飞行中解决约束。如果对自由度的值感兴趣（通常只需要正交点的值），可以通过
 * FEEvaluation::get_dof_value(i),
 * 访问这些值，其中i是基函数的索引。请注意，FEEvaluation中连续元素自由度的编号与FE_Q（或FEValues）中的排序不同，因为FEEvaluation需要以lexicographic顺序访问它们，例如FE_DGQ中使用的就是这种排序。重新索引的成本太高，因为evaluate()和integration()里面的访问是在张量评估部分的关键路径上。在evaluate()调用之前，通过read_dof_values()填充DoF值的一个替代方法是通过set_dof_value()调用手动赋值。同样，如果积分的局部结果应该被进一步处理，而不是通过distribut_local_to_global()分散到一个向量中，我们可以在调用
 * integrate()后通过get_dof_value()来访问它。在不同背景下使用积分值的一个例子是快速装配矩阵，如下一小节所示。
 * 对于大多数反复穿过网格的算子评估任务，MatrixFree的实现方式是将预先计算的映射数据（几何描述的雅各布变换）与基函数的即时评估相结合，是最有效的方式。换句话说，该框架在内存使用和对象的初始化之间选择了一种权衡，适合用无矩阵的方式替代矩阵-向量乘积或显式时间积分。
 * <h4>Usage without pre-initialized MatrixFree object</h4>
 * 第二种使用形式是通过FEValues生成的几何信息来初始化FEEvaluation。这允许在不事先初始化MatrixFree对象的情况下，即时应用积分循环。当MatrixFree的内存和初始化成本不可接受时，这可能很有用，例如，在误差计算中，不同数量的正交点应该被用于一次评估。另外，当使用这个类的例程来组装矩阵时，MatrixFree类所暗示的权衡可能是不可取的。在这种情况下，即时初始化必要的几何数据的成本是相当低的，因此避免全局对象MatrixFree是有用的。当以这种方式使用时，会使用让人想起带有单元格迭代器的FEValues的reinit方法。然而，请注意，这种模式的结果是一次只处理一个单元，几何数据在矢量化数组的所有组件中都是重复的。因此，只有在可以对不同的数据进行相同的操作时，矢量化才是有用的，例如在进行矩阵装配时。
 * 作为一个例子，考虑下面的代码来组装对拉普拉斯矩阵的贡献。
 *
 *
 * @code
 * FEEvaluation<dim,fe_degree> fe_eval (mapping, finite_element,
 *                                    QGauss<1>(fe_degree+1), flags);
 * for (const auto &cell : dof_handler.active_cell_iterators())
 * {
 *   fe_eval.reinit(cell);
 *   for (unsigned int i=0; i<dofs_per_cell;
 *        i += VectorizedArray<double>::size())
 *     {
 *       const unsigned int n_items =
 *         i+VectorizedArray<double>::size() > dofs_per_cell ?
 *         (dofs_per_cell
 *
 * - i) :
 *         VectorizedArray<double>::size();
 *
 *       // Set n_items unit vectors
 *       for (unsigned int j=0; j<dofs_per_cell; ++j)
 *         fe_eval.set_dof_value(VectorizedArray<double>(), j);
 *       for (unsigned int v=0; v<n_items; ++v)
 *         {
 *           VectorizedArray<double> one_value = VectorizedArray<double>();
 *           one_value[v] = 1.;
 *           fe_eval.set_dof_value(one_value, i+v);
 *         }
 *
 *       // Apply operator on unit vector to generate the next few matrix
 *       // columns
 *       fe_eval.evaluate(EvaluationFlags::values|EvaluationFlags::gradients);
 *       for (unsigned int q=0; q<n_q_points; ++q)
 *         {
 *           fe_eval.submit_value(10.*fe_eval.get_value(q), q);
 *           fe_eval.submit_gradient(fe_eval.get_gradient(q), q);
 *         }
 *       fe_eval.integrate(EvaluationFlags::values|EvaluationFlags::gradients);
 *
 *       // Insert computed entries in matrix
 *       for (unsigned int v=0; v<n_items; ++v)
 *         for (unsigned int j=0; j<dofs_per_cell; ++j)
 *           cell_matrix(fe_eval.get_internal_dof_numbering()[j],
 *                       fe_eval.get_internal_dof_numbering()[i+v])
 *             = fe_eval.get_dof_value(j)[v];
 *     }
 *   ...
 * }
 * @endcode
 *
 * 这段代码用上面 @p i
 * 的循环生成了单元格矩阵的列。这样做的方式如下。FEEvaluation的例程专注于有限元算子的评估，所以对于计算单元矩阵的算子评估，它被应用于单元上的所有单位向量。在单位向量上应用算子可能看起来效率不高，但这里使用的评估例程非常快，以至于它们的工作速度仍然比FEValues可能的要快得多。特别是，复杂度是
 * <code>(fe_degree+1)<sup>2*dim+1</sup> </code>  而不是
 * <code>(fe_degree+1)<sup>3*dim</sup> </code>  。
 * 由于矢量化，我们可以一次生成几个单位向量的矩阵列（例如4）。变量
 * @p n_items
 * 确保我们正确地做最后一次迭代，其中单元格的数量不被矢量化长度所除。还要注意的是，我们需要得到fe_eval应用的内部自由度编号，因为FEEvaluation内部使用的是自由度的词法编号，如上所述。
 * <h4>Internal data organization</h4>
 * 用于保存局部自由度的求解值以及正交点的内插值、梯度和Hessians的临时数据是由
 * MatrixFree::acquire_scratch_data()
 * 提供的Scratch数组，在对FEEvaluation的不同调用之间被重复使用。因此，构造一个FEEvaluation对象通常很便宜，不涉及任何昂贵的操作。在构建过程中，只有几十个指向实际数据字段的指针被设置。因此，在每个循环中多次创建一个FEEvaluation时，不会产生负面的性能影响，例如在一个`local_cell_operation`操作的顶部，该操作被分割成小块用于并行for循环，避免了在
 * @p WorkStream. 的循环中需要单独的抓取数据字段。
 * 当在多线程模式下使用FEEvaluation类时，MatrixFree中抓取数据的线程本地存储会自动确保每个线程得到它的私有数据阵列。然而，请注意，当所有的线程并行化是由外部提供的，而不是通过deal.II的例程完成的，如OpenMP，deal.II也必须被编译为支持线程。这是因为deal.II需要知道线程本地存储的符号。FEEvaluation内核已经被验证可以在OpenMP循环中工作。
 * <h4>Vectorization scheme through VectorizedArray</h4>
 * 该类旨在通过显式矢量化来执行现代CPU上存在的单指令多数据（SIMD）指令的所有算术，这些指令在deal.II中通过类VectorizedArray提供，使用配置/编译时可用的最宽矢量宽度。为了保持程序的灵活性，FEEvaluation总是在几个元素上应用矢量化。这通常是最好的妥协，因为在有限元方法中，不同元素上的计算通常是独立的（当然，除了向全局残差向量添加积分贡献的过程），在更复杂的情况下也是如此。例如，稳定参数可以定义为一个单元的所有正交点上的某些量的最大值除以该单元的体积，但不需要将结果与邻接点进行局部混合。使用计算机结构的术语，FEEvaluation的设计依赖于在典型的集成场景下对单元进行操作时不做任何跨线数据交换。
 * 当问题中的单元数不是SIMD向量中数组元素的倍数时，FEEvaluation的实现会在未使用的SIMD通道中填入一些假条目，并将它们带入周围，由于VectorizedArray的长度在编译时是固定的，所以这一选择是必要的。然而，与自动矢量化设置相比，这种方法通常会产生更好的代码，因为在自动矢量化设置中，除了在完全填充的通道上使用的矢量化版本外，还需要另一个未矢量化的代码路径，同时还有一个调度机制。在
 * @p read_dof_values,
 * 中，对不完整的一批单元的reinit()调用所产生的空道被设置为零，而
 * @p distribute_local_to_global 或 @p set_dof_values
 * 则简单地忽略了空道中的内容。实际填充的SIMD通道的数量可以通过
 * MatrixFree::n_components_filled(). 来查询。
 * 很明显，在人工通道上进行的计算（没有真实的数据）不应该与有效的结果混合。使用这个类的契约是，用户要确保通道在用户代码中不被交叉，特别是由于事先不清楚哪些单元会在矢量化中被放在一起。例如，除了通过全局向量访问方法或通过被屏蔽的访问
 * MatrixFree::n_components_filled().
 * ，一个元素上的结果不应该被添加到其他元素上的结果中。
 * 不能保证人工车道上的结果永远是可以安全地添加到其他结果中的零。JxW或Jacobian上的数据是从最后一个有效车道复制的，以避免除以零，这可能引发浮点异常或其他情况下的麻烦。
 * <h3>Description of evaluation routines</h3>
 * 这个类包含了基于张量积正交公式和类似张量积的形状函数的元素的专门评估例程，包括标准的FE_Q或FE_DGQ元素和围绕0.5对称的正交点（像高斯正交），基于截断张量积的FE_DGP元素以及高斯-洛巴托正交元素的更快情况，它给出对角线质量矩阵和内部更快评估。该类的主要优点是在所有正交中评估所有形状函数，或者在
 * <code>dim (fe_degree+1)<sup>dim+1</sup> </code>
 * 操作中对所有形状函数进行积分，而不是在FEValues的评估例程中的较慢的
 * <code> (fe_degree+1)<sup>2*dim</sup></code>
 * 复杂性。这是由一种叫做和因子化的算法完成的，该算法在沿坐标方向的评估过程中剔除了常数因子。这个算法是许多谱元算法的基础。
 * 请注意，通过这个类可以进行的许多操作都是从基类FEEvaluationBase继承的，特别是对向量的读写。另外，该类继承了FEEvaluationAccess，实现了对正交点上有限元函数的值、梯度和Hessians的访问。
 * 该类假定所考虑的有限元的形状函数 <em> 不 </em>
 * 依赖于实空间中单元的几何形状。目前，其他有限元不能用无矩阵概念处理。
 * <h4>Degree of finite element as a compile-time parameter</h4>
 * 该类FEEvaluation为两个使用模式。第一个使用模式是将多项式程度作为模板参数来指定。这保证了最大的效率。用和因数法进行评估时，会执行一些嵌套的短1D循环，其长度等于多项式度数加1。如果在编译时知道循环的边界，编译器可以根据其启发式方法认为最有效的方式展开循环。至少最里面的循环几乎总是被完全解开，避免了循环的开销。
 * 然而，将多项式度数（以及正交点的数量）作为模板参数携带，在需要考虑不同的多项式度数的代码中，例如在通过输入文件给出多项式度数的应用代码中，事情就变得更加复杂。第二种使用模式是依靠预编译的代码来处理多项式度数。虽然用户代码可以为单元格使用不同的函数（这些函数会被一些动态调度机制调用，用于不同的度数模板），但deal.II也支持基于传递给初始化的元素中的信息来使用这个类。对于这种使用模式，将多项式程度的模板参数设置为
 *
 * - 并为正交点的数量选择一个任意的数字。该代码部分包含预编译的模板代码，用于1到6之间的多项式度数和常见的正交公式，其运行速度几乎与模板版本相同。如果所选的度数没有被预编译，一个具有模板专业性的评估器对象将用于
 *
 * -被调用，根据运行时间的界限运行。
 * 下图给出了FEEvaluation的性能概览。它考虑了使用类似于
 * step-37
 * 单精度算术教程程序的代码，用连续有限元评估拉普拉斯的每个自由度所花费的时间。该时间是基于英特尔至强E5-2687W
 * v4单核的实验，运行频率为3.4GHz，在问题大小为1000万左右时测量的。该图列出了计算时间（约0.1秒）除以自由度数。
*  @image html fe_evaluation_laplacian_time_per_dof.png
 * 该图显示，模板化的计算核比非模板化的计算核快2.5到3倍。在这个设置上，最快的周转是对多项式5度的计算，每自由度7.4e-9秒，或每秒1.34亿自由度。
 *
 * - 在一个单核上。非模板化版本在多项式5度时也是最快的，每个自由度2.1e-9秒，或每秒4800万自由度。注意，使用模板`degree=-1`的FEEvaluation会选择1到6度之间的快速路径，而其他度数则选择慢速路径。
 * <h4>Pre-compiling code for more polynomial degrees</h4>
 * 也可以为不同的最大多项式度数预先编译FEEvaluation中的代码。这由类
 * internal::FEEvaluationFactory
 * 和`nclude/deal.II/matrix_free/evaluation_template_factory.templates.h`中的实现控制。通过设置宏`FE_EVAL_FACTORY_DEGREE_MAX`为所需的整数，并实例化类FEEvaluationFactory和FEFaceEvaluationFactory（后者用于FEFaceEvaluation），为可能更大的度数集创建模板函数的路径。你可以通过调用
 * FEEvaluation::fast_evaluation_supported() 或
 * FEFaceEvaluation::fast_evaluation_supported().
 * 来检查是否对一个给定的度数/n_正交点进行了快速评估/积分。
 * <h3>Handling multi-component systems</h3>
 * FEEvaluation还允许通过一个关于分量数量的模板参数来处理矢量值问题。
 *
 *
 * @code
 * FEEvaluation<dim,fe_degree,n_q_points_1d,n_components> phi(matrix_free);
 * @endcode
 *
 * 如果这样使用，可以通过调用从一个 @p std::vector<VectorType>
 * 的几个组件中收集组件
 *
 *
 * @code
 * phi.read_dof_values(src, 0);
 * @endcode
 *
 * 其中的0表示应该使用从 @p std::vector
 * 中的第2个向量开始的向量，<code>src[0], src[1], ...,
 * src[n_components-1]</code>。
 * 如果MatrixFree数据底层的DoFHandler是基于 @p
 * n_components项的FESystem，那么读取多成分系统的另一种方式是可能的。在这种情况下，为read_dof_values()和distribut_local_to_global()调用提供一个单一的向量。
 * FEEvaluation在多组件系统中的一个重要属性是在get_value()、get_gradient()或get_dof_value()调用中的多组件布局。在这种情况下，返回的不是标量字段
 * VectorizedArray  @<double@>  而是张量。
 *
 *
 * @code
 * get_value
 *
 * -> Tensor<1,n_components,VectorizedArray<double>>
 * get_gradient
 *
 * -> Tensor<1,n_components,Tensor<1,dim,VectorizedArray<double>>
 * @endcode
 *
 * 与此类似，submit_value()和submit_gradient()调用的是张量的值。请注意，存在对
 * @p n_components=1和 @p n_components=dim,
 * 的特殊化，它们是通过基类FEEvaluationAccess提供的。在标量情况下，这些提供了上述的标量返回类型。在矢量值的情况下，梯度从<code>Tensor
 * @<1,dim,Tensor@<1,dim,VectorizedArray@<double@>   @>   @></code>
 * 转换为<code>Tensor  @<2,dim,VectorizedArray@<double@>   @></code>.
 * 此外，额外的操作，如diveregence或curl也是可用的。
 * 如果不同的形状函数被结合起来，例如Stokes流中的混合有限元公式，将创建两个FEEvaluation对象，一个用于速度，一个用于压力。然后在正交点上进行组合。
 *
 *
 * @code
 * FEEvaluation<dim,degree_p+1,degree_p+2,dim> velocity (data, 0);
 * FEEvaluation<dim,degree_p,  degree_p+2,1, > pressure (data, 1);
 *
 * for (unsigned int cell=cell_range.first; cell<cell_range.second; ++cell)
 * {
 *   velocity.reinit (cell);
 *   velocity.read_dof_values (src.block(0));
 *   velocity.evaluate (EvaluationFlags::gradients);
 *   pressure.reinit (cell);
 *   pressure.read_dof_values (src.block(1));
 *   pressure.evaluate (EvaluationFlags::values);
 *
 *   for (unsigned int q=0; q<velocity.n_q_points; ++q)
 *     {
 *       SymmetricTensor<2,dim,VectorizedArray<double> > sym_grad_u =
 *         velocity.get_symmetric_gradient (q);
 *       VectorizedArray<double> pres = pressure.get_value(q);
 *       VectorizedArray<double> div =
 *
 * -trace(sym_grad_u);
 *       pressure.submit_value (div, q);
 *
 *       // subtract p I
 *       for (unsigned int d=0; d<dim; ++d)
 *         sym_grad_u[d][d]
 *
 * -= pres;
 *
 *       velocity.submit_symmetric_gradient(sym_grad_u, q);
 *    }
 *
 *   velocity.integrate (EvaluationFlags::gradients);
 *   velocity.distribute_local_to_global (dst.block(0));
 *   pressure.integrate (EvaluationFlags::values);
 *   pressure.distribute_local_to_global (dst.block(1));
 * }
 * @endcode
 *
 * 这段代码假设一个由两个分量组成的BlockVector分别描述了速度和压力分量。为了识别速度和压力的不同DoFHandler对象，FEEvaluation对象的第二个参数为速度指定各自的分量0，为压力指定1。对于矢量值问题的进一步例子，deal.II测试套件也包括一些额外的例子，例如，上面描述的斯托克斯算子可以在https://github.com/dealii/dealii/blob/master/tests/matrix_free/matrix_vector_stokes_noflux.cc
 * <h3>Handling several integration tasks and data storage in quadrature
 * points</h3>
 * FEEvaluation和MatrixFree的设计将几何学与基函数分开。因此，几个DoFHandler对象（或者同一个DoFHandler配备了不同的约束对象）可以共享相同的几何信息，就像上面的Stokes例子一样。所有的几何体在MatrixFree中被缓存一次，所以FEEvaluation不需要做昂贵的初始化调用，而是设置一些指针。这种实现是基于这样的想法：当几个字段被评估时，也只需要一次几何信息，这与FEValues不同，后者为每个字段设置了内部映射数据。例如，如果一个多分量的PDE涉及到一个分量的形状值和另一个分量的形状梯度，如果两者都基于同一个MatrixFree对象，并且更新标志指定
 * @p  update_values ,  @p update_gradients  , 和  @p update_jxw_values
 * 都给出，就不会失去效率。形状值所需数量的选择是通过evaluation()或integration调用中的标志和正交点的访问。
 *
 *
 * @code
 * phi1.evaluate(EvaluationFlags::values);
 * phi2.evaluate(EvaluationFlags::gradients);
 * for (unsigned int q=0; q<phi1.n_q_points; ++q)
 * {
 *   VectorizedArray<double> val1 = phi1.get_value(q);
 *   Tensor<1,dim,VectorizedArray<double> > grad2 = phi2.get_gradient(q);
 *   Point<dim,VectorizedArray<double> > point = phi1.quadrature_point(q);
 *   // ... some complicated formula combining those three...
 * }
 * @endcode
 *
 * 在正交点的循环中，我们可以要求两个FEEvaluation对象中的任何一个提供正交点的位置，这并不重要，因为它们只保留正交点数据的指针。
 * 这一观察结果也转化为在程序中实现不同的微分算子的情况，例如，在算法的一个阶段，质量矩阵的作用，而在另一个阶段，刚度矩阵的作用。只需要一个MatrixFree对象，通过在不同的FEEvaluation对象中使用不同的局部函数和各自的实现来保持充分的效率。换句话说，用户在为MatrixFree的初始化提供update_flags时，不需要为了效率的原因而费心保守
 *
 * - 除了在 FEEvaluation::reinit() 调用里面最多增加一两个 @p if 语句外，FEEvaluation内部不会产生任何开销。相反，从效率的角度来看，所有调用中必要的最大标志集是完全没有问题的。
 * 对于不同字段的组合，包括来自不同时间步骤的不同解向量，强制要求所有的FEEvaluation对象共享同一个MatrixFree对象。这是因为单元格被
 * MatrixFree::cell_loop()
 * 循环的方式对于不同的DoFHandler或AffineConstraints参数可能不同。更确切地说，即使在串行中布局是相同的，但在MPI情况下，对于不同的DoFHandler/AffineConstraints的排序没有任何保证。原因是该算法检测需要与MPI进行数据交换的单元，而这些单元对于不同的元素可能会发生变化，例如，具有悬挂节点约束的FE_Q比FE_DGQ元素连接到更多的邻居，而需要数据交换的单元被放在单元循环中的不同位置。当然，如果设置了完全相同的DoFHandler、AffineConstraints和选项（如线程并行的设置），那么顺序就会相同，因为算法是确定性的。
 * @tparam  该类所使用的dim维度
 * @tparam  fe_degree
 * 每个坐标方向具有fe_degree+1个自由度的张量积有限元的程度。可以设置为
 *
 * - 如果在编译时不知道度数，但性能通常会差2-3倍。
 * @tparam  n_q_points_1d
 * 一维正交公式中的点数，默认为fe_degree+1。
 * @tparam  n_components
 * 在求解PDEs系统时，向量分量的数量。如果同一个操作被应用于一个PDE的几个分量（例如，一个矢量拉普拉斯方程），它们可以通过一次调用同时应用（而且通常更有效）。默认为1。
 * @tparam  数字 数字格式，通常为 @p double 或 @p float.  默认为
 * @p 双数
 *
 *
 * @ingroup matrixfree
 *
 *
 */
template <int dim,
          int fe_degree,
          int n_q_points_1d,
          int n_components_,
          typename Number,
          typename VectorizedArrayType>
class FEEvaluation : public FEEvaluationAccess<dim,
                                               n_components_,
                                               Number,
                                               false,
                                               VectorizedArrayType>
{
  static_assert(
    std::is_same<Number, typename VectorizedArrayType::value_type>::value,
    "Type of Number and of VectorizedArrayType do not match.");

public:
  /**
   * 基类的别名。
   *
   */
  using BaseClass =
    FEEvaluationAccess<dim, n_components_, Number, false, VectorizedArrayType>;

  /**
   * 一个作为模板参数指定的底层数字类型。
   *
   */
  using number_type = Number;

  /**
   * 函数值的类型，例如 "n_components=1 "的 "VectorizedArrayType
   * "或 "n_components=dim "的 "Tensor<1,dim,VectorizedArrayType>"。
   *
   */
  using value_type = typename BaseClass::value_type;

  /**
   * 梯度的类型，例如，`Tensor<1,dim,VectorizedArrayType>`代表`n_components=1`或者`Tensor<2,dim,VectorizedArrayType
   * >`代表`n_components=dim`。
   *
   */
  using gradient_type = typename BaseClass::gradient_type;

  /**
   * 作为模板参数给出的尺寸。
   *
   */
  static constexpr unsigned int dimension = dim;

  /**
   * 作为模板参数给定的评估器的解决方案组件的数量。
   *
   */
  static constexpr unsigned int n_components = n_components_;

  /**
   * 从给定的模板参数`n_q_points_1d`确定的正交点的静态数量。请注意，如果给定了`fe_degree=-1`，并且使用了运行时的循环长度而不是编译时的循环长度，那么实际的正交点数量`n_q_points`可能不同。
   *
   */
  static constexpr unsigned int static_n_q_points =
    Utilities::pow(n_q_points_1d, dim);

  /**
   * 根据给定的模板参数`fe_degree`确定的标量分量的静态自由度数。请注意，如果给定了`fe_degree=-1`，或者如果底层的类型比通常的FE_Q或FE_DGQ更复杂，如FE_DGP，那么实际的自由度`dofs_per_component`可能会不同。
   *
   */
  static constexpr unsigned int static_dofs_per_component =
    Utilities::pow(fe_degree + 1, dim);

  /**
   * 根据给定的模板参数`fe_degree`确定的所有组件的静态自由度数。请注意，如果给定了`fe_degree=-1`，或者如果底层的类型比通常的FE_Q或FE_DGQ更复杂，比如FE_DGP，那么实际的自由度`dofs_per_cell`可能会不同。
   *
   */
  static constexpr unsigned int tensor_dofs_per_cell =
    static_dofs_per_component * n_components;

  /**
   * 根据给定的模板参数`fe_degree`确定的所有组件的静态自由度数。请注意，如果给定了`fe_degree=-1`，或者如果底层的类型比通常的FE_Q或FE_DGQ更复杂，比如FE_DGP，那么实际的自由度`dofs_per_cell`可能会不同。
   *
   */
  static constexpr unsigned int static_dofs_per_cell =
    static_dofs_per_component * n_components;

  /**
   * 构造函数。接受存储在MatrixFree中的所有数据。如果应用于在构建
   * @p matrix_free,
   * 过程中选择了一个以上的有限元或一个以上的正交公式的问题，可以通过可选的参数选择合适的组件。
   * @param  matrix_free 包含所有数据的数据对象  @param  dof_no
   * 如果matrix_free是用多个DoFHandler对象设置的，这个参数可以选择给定的评估器应该连接到哪个DoFHandler/AffineConstraints对。
   * @param  quad_no
   * 如果matrix_free被设置为多个正交对象，该参数将选择适当的正交公式编号。
   * @param  first_selected_component
   * 如果由dof_no选择的dof_handler使用一个由多个组件组成的FESystem，这个参数允许选择当前评估程序开始的组件。注意，一个评估器不支持在不同的组件中结合不同的形状函数。换句话说，FESystem的同一个基本元素需要为
   * @p first_selected_component 和
   * <code>first_selected_component+n_components_</code>
   * 之间的组件设置。      @param  active_fe_index
   * 如果matrix_free是用 hp::FECollections,
   * 的DoFHandler对象设置的，这个参数可以选择给定的评估器应该连接到哪个DoFHandler/AffineConstraints对。
   * @param  active_quad_index 如果matrix_free是用 hp::Collection
   * 对象设置的，该参数选择正交公式的适当编号。
   *
   */
  FEEvaluation(
    const MatrixFree<dim, Number, VectorizedArrayType> &matrix_free,
    const unsigned int                                  dof_no  = 0,
    const unsigned int                                  quad_no = 0,
    const unsigned int first_selected_component                 = 0,
    const unsigned int active_fe_index   = numbers::invalid_unsigned_int,
    const unsigned int active_quad_index = numbers::invalid_unsigned_int);

  /**
   * 构造函数。接受存储在MatrixFree中的所有数据，用于给定的单元格范围，这可以在p-adaptive策略的情况下自动识别active_fe_index和active_quad_index。
   * 其余的参数与上面的构造函数相同。
   *
   */
  FEEvaluation(const MatrixFree<dim, Number, VectorizedArrayType> &matrix_free,
               const std::pair<unsigned int, unsigned int> &       range,
               const unsigned int                                  dof_no  = 0,
               const unsigned int                                  quad_no = 0,
               const unsigned int first_selected_component                 = 0);

  /**
   * 构造函数，功能减少，工作方式与FEValues类似。参数与传递给FEValues的构造函数的参数类似，值得注意的是，FEEvaluation期望一个一维的正交公式Quadrature<1>，而不是一个
   * @p dim
   * 维的公式。有限元既可以是标量也可以是矢量值，但是这个方法每次总是只选择一个标量基元（按照类模板指定的
   * @p n_components 副本）。对于向量值元素，可选的参数 @p
   * first_selected_component
   * 允许指定用于评估的基础元素的索引。注意，内部数据结构总是假定基元是原始的，目前不支持非原始的。
   * 正如从FEValues中得知的那样，调用带有
   * Triangulation<dim>::cell_iterator
   * 的reinit方法是必要的，以使当前类的几何和自由度为人所知。如果迭代器包括DoFHandler信息（即它是一个
   * DoFHandler<dim>::cell_iterator 或类似的），初始化也允许以
   * DoFHandler<dim>::active_cell_iterator
   * 类型的标准方式一次从向量中读取或写入一个单元。然而，这种方法比使用MPI的MatrixFree的路径要慢得多，因为必须进行索引转换。由于每次只使用一个单元，这种方法不会在几个元素上进行矢量化（这对矢量操作来说是最有效的），而只可能在元素内进行矢量化，如果评估/整合例程在用户代码内组合的话（例如，计算单元矩阵）。
   *
   */
  FEEvaluation(const Mapping<dim> &      mapping,
               const FiniteElement<dim> &fe,
               const Quadrature<1> &     quadrature,
               const UpdateFlags         update_flags,
               const unsigned int        first_selected_component = 0);

  /**
   * 缩小功能的构造函数。这个构造函数等同于另一个构造函数，除了它使对象隐含地使用
   * $Q_1$ 映射（即，MappingQGeneric(1)类型的对象）。
   *
   */
  FEEvaluation(const FiniteElement<dim> &fe,
               const Quadrature<1> &     quadrature,
               const UpdateFlags         update_flags,
               const unsigned int        first_selected_component = 0);

  /**
   * 缩小功能的构造函数。类似于其他带有FiniteElement参数的构造函数，但使用另一个FEEvaluationBase对象来提供关于几何的信息。这允许几个FEEvaluation对象共享几何评估，即底层映射和正交点只需要评估一次。当你打算使用与给定对象平行的FEEvaluation对象时，请确保不要传递一个可选的对象，因为否则打算共享的对象可能会产生竞赛条件。
   *
   */
  FEEvaluation(
    const FiniteElement<dim> &                                           fe,
    const FEEvaluationBaseData<dim, Number, false, VectorizedArrayType> &other,
    const unsigned int first_selected_component = 0);

  /**
   * 复制构造函数。如果FEEvaluationBase是由映射、fe、正交和更新标志构建的，基于FEValues的底层几何评估将被深度复制，以便允许与线程并行使用。
   *
   */
  FEEvaluation(const FEEvaluation &other);

  /**
   * 复制赋值运算符。如果FEEvaluationBase是由映射、fe、正交和更新标志构建的，基于FEValues的底层几何评估将被深度复制，以允许与线程并行使用。
   *
   */
  FEEvaluation &
  operator=(const FEEvaluation &other);

  /**
   * 将操作指针初始化为当前的单元格批索引。与下面以单元格迭代器为参数的reinit函数和
   * FEValues::reinit()
   * 方法不同，与特定单元格相关的信息是在reinit调用中生成的，这个函数非常便宜，因为所有数据都是在
   * @p matrix_free,
   * 中预先计算的，只有少数指数需要适当设置。
   *
   */
  void
  reinit(const unsigned int cell_batch_index);

  /**
   * 如同FEValues中的惯例，使用TriaIterator对象将数据初始化到当前单元。参数类型为
   * DoFHandler::active_cell_iterator 或 DoFHandler::level_cell_iterator.
   * 只有在创建FEEvaluation对象时带有有限元、正交公式和正确的更新标志以及<b>without</b>MatrixFree对象时，该选项才能使用。这种初始化方法失去了使用矢量化的能力，另见FEEvaluation类的描述。当使用这种重启方法时，FEEvaluation也可以从向量中读取数据（但效率比来自MatrixFree的数据低）。
   *
   */
  template <bool level_dof_access>
  void
  reinit(const TriaIterator<DoFCellAccessor<dim, dim, level_dof_access>> &cell);

  /**
   * 像在FEValues中一样使用TriaIterator对象将数据初始化到当前单元。这个选项只有在FEEvaluation对象是用有限元、正交公式和正确的更新标志以及<b>without</b>MatrixFree对象创建的情况下才可用。这种初始化方法失去了使用矢量化的能力，另见FEEvaluation类的描述。当使用这种重启方法时，FEEvaluation可以<b>not</b>从矢量中读取，因为没有DoFHandler信息可用。
   *
   */
  void
  reinit(const typename Triangulation<dim>::cell_iterator &cell);

  /**
   * 检查是否支持面部评估/整合。
   *
   */
  static bool
  fast_evaluation_supported(const unsigned int given_degree,
                            const unsigned int give_n_q_points_1d);

  /**
   * 评估函数值、梯度和从输入矢量中的DoF值到单元格上的正交点的多项式内插的Hessians。
   * 函数参数指定哪些部分应被实际计算。这个函数必须首先被调用，以便访问函数
   * @p get_value(),   @p get_gradient()  或  @p
   * get_laplacian提供有用的信息（除非这些值已经被手动设置）。
   *
   */
  void
  evaluate(const EvaluationFlags::EvaluationFlags evaluation_flag);

  /**
   * 和上面一样，但有单独的bool标志。    @deprecated
   * 使用evaluate()与EvaluationFlags参数。
   *
   */
  DEAL_II_DEPRECATED_EARLY void
  evaluate(const bool evaluate_values,
           const bool evaluate_gradients,
           const bool evaluate_hessians = false);

  /**
   * 评估从输入数组 @p
   * values_array中的DoF值到单元格上的正交点的函数值、梯度和多项式内插的Hessians。如果当前的FEEvaluation对象涉及多个部件，
   * @p values_array
   * 中的排序是第一个部件的所有自由度排在前面，然后是第二个部件的所有自由度，以此类推。函数参数指定哪些部分应被实际计算。这个函数必须先被调用，这样访问函数
   * @p get_value(),  @p get_gradient() 或 @p get_laplacian
   * 才能提供有用的信息（除非这些值是手动设置的）。
   *
   */
  void
  evaluate(const VectorizedArrayType *            values_array,
           const EvaluationFlags::EvaluationFlags evaluation_flag);

  /**
   * 像上面一样，但使用单独的bool标志。    @deprecated
   * 使用evaluate()与EvaluationFlags参数。
   *
   */
  DEAL_II_DEPRECATED_EARLY void
  evaluate(const VectorizedArrayType *values_array,
           const bool                 evaluate_values,
           const bool                 evaluate_gradients,
           const bool                 evaluate_hessians = false);

  /**
   * 从输入向量中读取并评估函数值、梯度以及从与当前单元格相关的
   * @p input_vector
   * 向量条目到单元格上的正交点的多项式内插的Hessians。函数参数指定哪些部分应被实际计算。这个函数必须首先被调用，以便访问函数
   * @p get_value(),   @p get_gradient() 或 @p
   * get_laplacian提供有用的信息（除非这些值已经被手动设置）。
   * 这个调用等同于调用read_dof_values()，然后再调用evaluate()，但内部可能会使用一些额外的优化。
   *
   */
  template <typename VectorType>
  void
  gather_evaluate(const VectorType &                     input_vector,
                  const EvaluationFlags::EvaluationFlags evaluation_flag);

  /**
   * @deprecated  请使用带有EvaluationFlags参数的
   * gather_evaluate()函数。
   *
   */
  template <typename VectorType>
  DEAL_II_DEPRECATED_EARLY void
  gather_evaluate(const VectorType &input_vector,
                  const bool        evaluate_values,
                  const bool        evaluate_gradients,
                  const bool        evaluate_hessians = false);

  /**
   * 这个函数获取存储在正交点上的值和/或梯度，通过单元上的所有基函数/梯度进行测试，并执行单元积分。两个函数参数
   * @p integrate_values 和 @p integrate_gradients
   * 分别用于启用/禁用提交给数值或梯度槽的贡献求和。结果被写入内部数据域
   * @p dof_values
   * （通常由distribut_local_to_global()或set_dof_values()方法写入结果向量中）。
   *
   */
  void
  integrate(const EvaluationFlags::EvaluationFlags integration_flag);


  /**
   * @deprecated  请使用带EvaluationFlags参数的 integrate() 函数。
   *
   */
  DEAL_II_DEPRECATED_EARLY void
  integrate(const bool integrate_values, const bool integrate_gradients);

  /**
   * 该函数获取存储在正交点上的值和/或梯度，通过单元上的所有基函数/梯度进行测试，并执行单元积分。两个函数参数
   * @p  integrate_values 和  @p integrate_gradients
   * 分别用于启用/禁用提交给数值或梯度槽的贡献的求和。与其他integration()方法相反，这个调用将测试结果存储在给定的数组
   * @p values_array,
   * 中，其先前的结果被覆盖，而不是写在begin_dof_values()后面的内部数据结构中。
   *
   */
  void
  integrate(const EvaluationFlags::EvaluationFlags integration_flag,
            VectorizedArrayType *                  values_array);

  /**
   * @deprecated  请使用带EvaluationFlags参数的 integrate() 函数。
   *
   */
  DEAL_II_DEPRECATED_EARLY void
  integrate(const bool           integrate_values,
            const bool           integrate_gradients,
            VectorizedArrayType *values_array);

  /**
   * 该函数获取存储在正交点上的值和/或梯度，通过单元上的所有基函数/梯度进行测试，执行单元积分，并将结果加入与当前单元索引相关的自由度上的全局向量
   * @p output_vector 。两个函数参数 @p integrate_values 和 @p
   * integrate_gradients
   * 分别用于启用/禁用提交给数值或梯度槽的贡献求和。
   * 这个调用等同于调用 integrate() 后面的
   * distribute_local_to_global()
   * ，但可能在内部使用一些额外的优化。
   *
   */
  template <typename VectorType>
  void
  integrate_scatter(const EvaluationFlags::EvaluationFlags evaluation_flag,
                    VectorType &                           output_vector);

  /**
   * @deprecated  请使用带有EvaluationFlags参数的 integrate_scatter()
   * 函数。
   *
   */
  template <typename VectorType>
  DEAL_II_DEPRECATED_EARLY void
  integrate_scatter(const bool  integrate_values,
                    const bool  integrate_gradients,
                    VectorType &output_vector);

  /**
   * 返回存储在MappingInfo中的实坐标的第q个正交点。
   *
   */
  Point<dim, VectorizedArrayType>
  quadrature_point(const unsigned int q_point) const;

  /**
   * 底层评价对象的单元上的单个组件的自由度数。通常接近static_dofs_per_component，但这个数字取决于实际选择的元素，因此不是静态的。
   *
   */
  const unsigned int dofs_per_component;

  /**
   * 单元上的自由度数量在当前评估对象的所有元件上累积。通常接近static_dofs_per_cell
   * =
   * static_dofs_per_component*n_components，但这个数字取决于所选择的实际元素，因此不是静态的。
   *
   */
  const unsigned int dofs_per_cell;

  /**
   * 使用中的正交点的数量。如果1d中的正交点数量是作为模板给出的，这个数字只是该值的<tt>dim</tt>-次幂。如果元素度数被设置为
   *
   * -（元素度的动态选择），正交点的静态值就不准确了，必须用这个值来代替。
   *
   */
  const unsigned int n_q_points;

private:
  /**
   * 检查关于元素度的模板参数是否与初始化时使用的实际元素相一致。
   *
   */
  void
  check_template_arguments(const unsigned int fe_no,
                           const unsigned int first_selected_component);
};



/**
 * 该类提供了在正交点评估函数和面积分所需的所有函数。该类的设计与FEEvaluation相似，大部分接口与该类共享，特别是大部分访问函数来自共同的基类FEEvaluationAccess和FEEvaluationBase。此外，该类与FEEvaluation的关系类似于FEValues与FEFaceValues的关系。
 * @tparam  该类将被用于的dim尺寸
 * @tparam  fe_degree
 * 每个坐标方向具有fe_degree+1个自由度的张量积有限元的程度。如果设置为
 *
 * - ，将使用底层元素的度数，它作为一个运行时常数，而不是减缓执行速度的编译时常数。
 * @tparam  n_q_points_1d
 * 一维正交公式中的点数，通常选择为fe_degree+1
 * @tparam  n_components
 * 在求解PDEs系统时，向量分量的数量。如果相同的操作应用于一个PDE的几个分量（例如一个矢量拉普拉斯方程），它们可以通过一次调用同时应用（而且通常更有效）。
 * @tparam  数字 Number 格式，通常为  @p double  或  @p float  。
 * @tparam  VectorizedArrayType
 * 以矢量方式处理的数组类型，默认为 VectorizedArray<Number>。
 *
 *
 * @note  目前只有VectorizedArray<Number,
 * width>被支持作为VectorizedArrayType。
 *
 *
 */
template <int dim,
          int fe_degree,
          int n_q_points_1d            = fe_degree + 1,
          int n_components_            = 1,
          typename Number              = double,
          typename VectorizedArrayType = VectorizedArray<Number>>
class FEFaceEvaluation : public FEEvaluationAccess<dim,
                                                   n_components_,
                                                   Number,
                                                   true,
                                                   VectorizedArrayType>
{
  static_assert(
    std::is_same<Number, typename VectorizedArrayType::value_type>::value,
    "Type of Number and of VectorizedArrayType do not match.");

public:
  /**
   * 基类的别名。
   *
   */
  using BaseClass =
    FEEvaluationAccess<dim, n_components_, Number, true, VectorizedArrayType>;

  /**
   * 一个作为模板参数指定的底层数字类型。
   *
   */
  using number_type = Number;

  /**
   * 函数值的类型，例如 "n_components=1 "的 "VectorizedArrayType
   * "或 "n_components=dim "的 "Tensor<1,dim,VectorizedArrayType>"。
   *
   */
  using value_type = typename BaseClass::value_type;

  /**
   * 梯度的类型，例如，`Tensor<1,dim,VectorizedArrayType>`代表`n_components=1`或`Tensor<2,dim,VectorizedArrayType
   * >`代表`n_components=dim`。
   *
   */
  using gradient_type = typename BaseClass::gradient_type;

  /**
   * 作为模板参数给出的尺寸。
   *
   */
  static constexpr unsigned int dimension = dim;

  /**
   * 作为模板参数给定的评估器的解决方案组件的数量。
   *
   */
  static constexpr unsigned int n_components = n_components_;

  /**
   * 根据给定的模板参数`n_q_points_1d`确定的正交点的静态数量，取为dim-1的幂。请注意，如果给定了`fe_degree=-1`，并且使用了运行时的循环长度而不是编译时的长度，那么实际的正交点数量`n_q_points`可能不同。
   *
   */
  static constexpr unsigned int static_n_q_points =
    Utilities::pow(n_q_points_1d, dim - 1);

  /**
   * 一个单元格上具有相同正交公式的正交点的静态数量。请注意，这个值只是为了与单元格正交进行更简单的比较，因为实际的点数是由`static_n_q_points`变量赋予一个面的。
   *
   */
  static constexpr unsigned int static_n_q_points_cell =
    Utilities::pow(n_q_points_1d, dim);

  /**
   * 从给定的模板参数`fe_degree'确定的标量分量的静态自由度数。注意，如果给定`fe_degree=-1`，实际的自由度`dofs_per_component`可能不同。
   *
   */
  static constexpr unsigned int static_dofs_per_component =
    Utilities::pow(fe_degree + 1, dim);

  /**
   * 根据给定的模板参数`fe_degree`确定的所有组件的静态自由度数。注意，如果给定了`fe_degree=-1`，实际的自由度`dofs_per_cell`可能不同。
   *
   */
  static constexpr unsigned int tensor_dofs_per_cell =
    static_dofs_per_component * n_components;

  /**
   * 根据给定的模板参数`fe_degree`确定的所有组件的静态自由度数。注意，如果给定了`fe_degree=-1`，实际的自由度`dofs_per_cell`可能不同。
   *
   */
  static constexpr unsigned int static_dofs_per_cell =
    static_dofs_per_component * n_components;

  /**
   * 构造函数。接受存储在MatrixFree中的所有数据。如果应用于在构造
   * @p matrix_free,
   * 过程中选择了一个以上的有限元或一个以上的正交公式的问题，可以通过可选的参数选择适当的组件。
   * @param  matrix_free 包含所有数据的数据对象  @param
   * is_interior_face
   * 这选择了当前评估器将基于内部面的两个单元中的哪个。内部面是主要的面，法线向量是沿着它的方向的。来自另一侧的外部面提供了与内部面相同的法向量，所以如果想要该面的外部法向量，必须乘以
   *
   * - .       @param  dof_no 如果matrix_free是用多个DoFHandler对象设置的，这个参数可以选择给定的评估器应该连接到哪个DoFHandler/AffineConstraints对。      @param  quad_no 如果matrix_free被设置为多个正交对象，该参数将选择适当的正交公式编号。      @param  first_selected_component 如果由dof_no选择的dof_handler使用一个由多个基本元素组成的FESystem，这个参数选择FESystem中的基本元素的编号。请注意，由于元素中可能存在多重性，这与各元素的分量没有直接关系。      @param  active_fe_index 如果matrix_free是用 hp::FECollections, 的DoFHandler对象设置的，这个参数可以选择给定的评估器应该连接到哪个DoFHandler/AffineConstraints对。      @param  face_type 在面的情况下，指示其参考单元的类型（0代表线或四边形，1代表三角形）。      @param  active_quad_index 如果matrix_free是用 hp::Collection 对象设置的，这个参数会选择正交公式的适当编号。
   *
   */
  FEFaceEvaluation(
    const MatrixFree<dim, Number, VectorizedArrayType> &matrix_free,
    const bool                                          is_interior_face = true,
    const unsigned int                                  dof_no           = 0,
    const unsigned int                                  quad_no          = 0,
    const unsigned int first_selected_component                          = 0,
    const unsigned int active_fe_index   = numbers::invalid_unsigned_int,
    const unsigned int active_quad_index = numbers::invalid_unsigned_int,
    const unsigned int face_type         = numbers::invalid_unsigned_int);

  /**
   * 构造函数。接受存储在MatrixFree中的所有数据，用于给定的面的范围，这可以在p自适应策略的情况下自动识别active_fe_index和active_quad_index。
   * 其余的参数与上面的构造函数相同。
   *
   */
  FEFaceEvaluation(
    const MatrixFree<dim, Number, VectorizedArrayType> &matrix_free,
    const std::pair<unsigned int, unsigned int> &       range,
    const bool                                          is_interior_face = true,
    const unsigned int                                  dof_no           = 0,
    const unsigned int                                  quad_no          = 0,
    const unsigned int first_selected_component                          = 0);

  /**
   * 将操作指针初始化为当前的面。这个方法是面集成的默认选择，因为MappingInfo中存储的数据是按照这个编号存储的。与下面以单元格迭代器为参数的reinit函数和
   * FEValues::reinit()
   * 方法不同，与特定单元格相关的信息是在reinit调用中生成的，这个函数非常便宜，因为所有数据都是在
   * @p matrix_free,
   * 中预先计算的，只有少数索引和指针需要适当设置。
   *
   */
  void
  reinit(const unsigned int face_batch_number);

  /**
   * 相对于基类中的reinit()方法，这个reinit()方法是针对给定的单元数和面数进行初始化。这个方法的效率低于其他采取面数的reinit()方法，因为它需要在这个调用中把与面数相关的数据复制到单元中。
   *
   */
  void
  reinit(const unsigned int cell_batch_number, const unsigned int face_number);

  /**
   * 检查是否支持面的评估/整合。
   *
   */
  static bool
  fast_evaluation_supported(const unsigned int given_degree,
                            const unsigned int give_n_q_points_1d);

  /**
   * 评估FE函数的函数值、梯度和Laplacians，这些都是存储在内部数据域`dof_values`中的DoF值（通常由read_dof_values()方法填充），在单元格上的正交点。
   * 函数参数指定哪些部分应被实际计算。需要在函数get_value()、get_gradient()或get_normal_derivative()提供有用信息之前调用（除非这些值是通过访问内部数据指针手动设置的）。
   *
   */
  void
  evaluate(const EvaluationFlags::EvaluationFlags evaluation_flag);

  /**
   * @deprecated  请使用带EvaluationFlags参数的evaluation()函数。
   *
   */
  DEAL_II_DEPRECATED_EARLY void
  evaluate(const bool evaluate_values, const bool evaluate_gradients);

  /**
   * 评估FE函数的函数值、梯度以及在输入数组`values_array`中的DoF值在单元格上的正交点处给出的拉普拉斯。如果当前的FEEvaluation对象涉及多个组件，values_array中的排序是：第一个组件的所有自由度排在前面，然后是第二个组件的所有自由度，以此类推。函数参数指定哪些部分应被实际计算。需要在函数get_value()、get_gradient()或get_normal_derivative()提供有用信息之前调用（除非这些值已经被手动设置）。
   *
   */
  void
  evaluate(const VectorizedArrayType *            values_array,
           const EvaluationFlags::EvaluationFlags evaluation_flag);

  /**
   * @deprecated  请使用带EvaluationFlags参数的evaluation()函数。
   *
   */
  DEAL_II_DEPRECATED_EARLY void
  evaluate(const VectorizedArrayType *values_array,
           const bool                 evaluate_values,
           const bool                 evaluate_gradients);

  /**
   * 从输入向量中读取并评估FE函数在单元格上的正交点的函数值、梯度和拉普拉斯。函数参数指定哪些部分应被实际计算。需要在函数get_value(),
   * get_gradient(),
   * 或get_normal_derivative()提供有用信息之前调用。
   * 这个调用等同于调用read_dof_values()，然后再调用evaluate()，但内部可能会使用一些额外的优化措施。
   *
   */
  template <typename VectorType>
  void
  gather_evaluate(const VectorType &                     input_vector,
                  const EvaluationFlags::EvaluationFlags evaluation_flag);

  /**
   * @deprecated
   * 请使用带有EvaluationFlags参数的gather_evaluate（）函数。
   *
   */
  template <typename VectorType>
  DEAL_II_DEPRECATED_EARLY void
  gather_evaluate(const VectorType &input_vector,
                  const bool        evaluate_values,
                  const bool        evaluate_gradients);

  /**
   * 这个函数获取存储在正交点上的值和/或梯度，通过单元上的所有基函数/梯度来测试它们，并执行单元积分。两个函数参数`integrate_val`和`integrate_grad`用于启用/禁用某些值或梯度。结果被写入内部数据域`dof_values`（通常由distribut_local_to_global()或set_dof_values()方法写入结果向量中）。
   *
   */
  void
  integrate(const EvaluationFlags::EvaluationFlags evaluation_flag);

  /**
   * @deprecated  请使用带EvaluationFlags参数的 integrate() 函数。
   *
   */
  DEAL_II_DEPRECATED_EARLY void
  integrate(const bool integrate_values, const bool integrate_gradients);

  /**
   * 这个函数获取存储在正交点上的值和/或梯度，通过单元上的所有基函数/梯度进行测试，并执行单元积分。两个函数参数`integrate_val`和`integrate_grad`用于启用/禁用某些值或梯度。与其他
   * integrate()方法相反，这个调用将测试结果存储在给定的数组`values_array`中。
   *
   */
  void
  integrate(const EvaluationFlags::EvaluationFlags evaluation_flag,
            VectorizedArrayType *                  values_array);

  /**
   * @deprecated  请使用带EvaluationFlags参数的 integrate() 函数。
   *
   */
  DEAL_II_DEPRECATED_EARLY void
  integrate(const bool           integrate_values,
            const bool           integrate_gradients,
            VectorizedArrayType *values_array);

  /**
   * 这个函数获取存储在正交点上的值和/或梯度，通过单元上的所有基函数/梯度进行测试，并执行单元积分。两个函数参数`integrate_val`和`integrate_grad`用于启用/禁用某些值或梯度。
   * 这个调用等同于调用 integrate() 后面的
   * distribute_local_to_global()
   * ，但内部可能使用一些额外的优化。
   *
   */
  template <typename VectorType>
  void
  integrate_scatter(const EvaluationFlags::EvaluationFlags evaluation_flag,
                    VectorType &                           output_vector);

  /**
   * @deprecated  请使用带有EvaluationFlags参数的 integrate_scatter()
   * 函数。
   *
   */
  template <typename VectorType>
  void
  integrate_scatter(const bool  integrate_values,
                    const bool  integrate_gradients,
                    VectorType &output_vector);

  /**
   * 返回存储在MappingInfo中的实坐标的面的第q个正交点。
   *
   */
  Point<dim, VectorizedArrayType>
  quadrature_point(const unsigned int q_point) const;

  /**
   * 底层评估对象的单元格上单个组件的自由度数。通常接近于static_dofs_per_component，但是这个数字取决于实际选择的元素，因此不是静态的。
   *
   */
  const unsigned int dofs_per_component;

  /**
   * 单元上的自由度数量在当前评估对象的所有元件上累积。通常接近static_dofs_per_cell
   * =
   * static_dofs_per_component*n_components，但这个数字取决于所选择的实际元素，因此不是静态的。
   *
   */
  const unsigned int dofs_per_cell;

  /**
   * 使用中的正交点的数量。如果1d中的正交点的数量是作为模板给出的，这个数字只是该值的<tt>dim-1</tt>-次幂。如果元素度数被设置为
   *
   * -（元素度的动态选择），正交点的静态值就不准确了，必须用这个值代替。
   *
   */
  const unsigned int n_q_points;


private:
  /**
   * 返回当前面批中每个面的面号。
   *
   */
  std::array<unsigned int, VectorizedArrayType::size()>
  compute_face_no_data();

  /**
   * 确定当前面批中每个面的方向。
   *
   */
  std::array<unsigned int, VectorizedArrayType::size()>
  compute_face_orientations();
};



namespace internal
{
  namespace MatrixFreeFunctions
  {
    // a helper function to compute the number of DoFs of a DGP element at
    // compile time, depending on the degree
    template <int dim, int degree>
    struct DGP_dofs_per_component
    {
      // this division is always without remainder
      static constexpr unsigned int value =
        (DGP_dofs_per_component<dim - 1, degree>::value * (degree + dim)) / dim;
    };

    // base specialization: 1d elements have 'degree+1' degrees of freedom
    template <int degree>
    struct DGP_dofs_per_component<1, degree>
    {
      static constexpr unsigned int value = degree + 1;
    };
  } // namespace MatrixFreeFunctions
} // namespace internal


 /*----------------------- Inline functions ----------------------------------*/ 

#ifndef DOXYGEN


 /*----------------------- FEEvaluationBaseData ------------------------*/ 

template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline FEEvaluationBaseData<dim, Number, is_face, VectorizedArrayType>::
  FEEvaluationBaseData(
    const MatrixFree<dim, Number, VectorizedArrayType> &data_in,
    const unsigned int                                  dof_no,
    const unsigned int first_selected_component,
    const unsigned int quad_no_in,
    const unsigned int fe_degree,
    const unsigned int n_q_points,
    const bool         is_interior_face,
    const unsigned int active_fe_index_in,
    const unsigned int active_quad_index_in,
    const unsigned int face_type)
  : scratch_data_array(data_in.acquire_scratch_data())
  , quad_no(quad_no_in)
  , matrix_info(&data_in)
  , dof_info(&data_in.get_dof_info(dof_no))
  , mapping_data(
      internal::MatrixFreeFunctions::
        MappingInfoCellsOrFaces<dim, Number, is_face, VectorizedArrayType>::get(
          data_in.get_mapping_info(),
          quad_no))
  , active_fe_index(fe_degree != numbers::invalid_unsigned_int ?
                      data_in.get_dof_info(dof_no).fe_index_from_degree(
                        first_selected_component,
                        fe_degree) :
                      (active_fe_index_in != numbers::invalid_unsigned_int ?
                         active_fe_index_in :
                         0))
  , active_quad_index(
      fe_degree != numbers::invalid_unsigned_int ?
        (mapping_data->quad_index_from_n_q_points(n_q_points)) :
        (active_quad_index_in != numbers::invalid_unsigned_int ?
           active_quad_index_in :
           std::min<unsigned int>(active_fe_index,
                                  mapping_data->descriptor.size() - 1)))
  , descriptor(
      &mapping_data->descriptor
         [is_face ?
            (active_quad_index * std::max<unsigned int>(1, dim - 1) +
             (face_type == numbers::invalid_unsigned_int ? 0 : face_type)) :
            active_quad_index])
  , n_quadrature_points(descriptor->n_q_points)
  , data(&data_in.get_shape_info(
      dof_no,
      quad_no_in,
      dof_info->component_to_base_index[first_selected_component],
      active_fe_index,
      active_quad_index))
  , jacobian(nullptr)
  , J_value(nullptr)
  , normal_vectors(nullptr)
  , normal_x_jacobian(nullptr)
  , quadrature_weights(descriptor->quadrature_weights.begin())
  , cell(numbers::invalid_unsigned_int)
  , is_interior_face(is_interior_face)
  , dof_access_index(
      is_face ?
        (is_interior_face ?
           internal::MatrixFreeFunctions::DoFInfo::dof_access_face_interior :
           internal::MatrixFreeFunctions::DoFInfo::dof_access_face_exterior) :
        internal::MatrixFreeFunctions::DoFInfo::dof_access_cell)
  , cell_type(internal::MatrixFreeFunctions::general)
{
  Assert(matrix_info->mapping_initialized() == true, ExcNotInitialized());
  AssertDimension(matrix_info->get_task_info().vectorization_length,
                  VectorizedArrayType::size());
  AssertDimension(n_quadrature_points, descriptor->n_q_points);
}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline FEEvaluationBaseData<dim, Number, is_face, VectorizedArrayType>::
  FEEvaluationBaseData(
    const Mapping<dim> &      mapping,
    const FiniteElement<dim> &fe,
    const Quadrature<1> &     quadrature,
    const UpdateFlags         update_flags,
    const unsigned int        first_selected_component,
    const FEEvaluationBaseData<dim, Number, is_face, VectorizedArrayType>
      *other)
  : scratch_data_array(new AlignedVector<VectorizedArrayType>())
  , quad_no(numbers::invalid_unsigned_int)
  , active_fe_index(numbers::invalid_unsigned_int)
  , active_quad_index(numbers::invalid_unsigned_int)
  , descriptor(nullptr)
  , n_quadrature_points(
      Utilities::fixed_power < is_face ? dim - 1 : dim > (quadrature.size()))
  , matrix_info(nullptr)
  , dof_info(nullptr)
  , mapping_data(nullptr)
  ,
  // select the correct base element from the given FE component
  data(new internal::MatrixFreeFunctions::ShapeInfo<VectorizedArrayType>(
    Quadrature<dim - is_face>(quadrature),
    fe,
    fe.component_to_base_index(first_selected_component).first))
  , jacobian(nullptr)
  , J_value(nullptr)
  , normal_vectors(nullptr)
  , normal_x_jacobian(nullptr)
  , quadrature_weights(nullptr)
  , cell(0)
  , cell_type(internal::MatrixFreeFunctions::general)
  , is_interior_face(true)
  , dof_access_index(internal::MatrixFreeFunctions::DoFInfo::dof_access_cell)
{
  Assert(other == nullptr || other->mapped_geometry.get() != nullptr,
         ExcInternalError());
  if (other != nullptr &&
      other->mapped_geometry->get_quadrature() == quadrature)
    mapped_geometry = other->mapped_geometry;
  else
    mapped_geometry =
      std::make_shared<internal::MatrixFreeFunctions::
                         MappingDataOnTheFly<dim, Number, VectorizedArrayType>>(
        mapping, quadrature, update_flags);
  cell = 0;

  mapping_data = &mapped_geometry->get_data_storage();
  jacobian     = mapped_geometry->get_data_storage().jacobians[0].begin();
  J_value      = mapped_geometry->get_data_storage().JxW_values.begin();
}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline FEEvaluationBaseData<dim, Number, is_face, VectorizedArrayType>::
  FEEvaluationBaseData(
    const FEEvaluationBaseData<dim, Number, is_face, VectorizedArrayType>
      &other)
  : scratch_data_array(other.matrix_info == nullptr ?
                         new AlignedVector<VectorizedArrayType>() :
                         other.matrix_info->acquire_scratch_data())
  , quad_no(other.quad_no)
  , active_fe_index(other.active_fe_index)
  , active_quad_index(other.active_quad_index)
  , descriptor(other.descriptor == nullptr ? nullptr : other.descriptor)
  , n_quadrature_points(other.n_quadrature_points)
  , matrix_info(other.matrix_info)
  , dof_info(other.dof_info)
  , mapping_data(other.mapping_data)
  , data(other.matrix_info == nullptr ?
           new internal::MatrixFreeFunctions::ShapeInfo<VectorizedArrayType>(
             *other.data) :
           other.data)
  , jacobian(nullptr)
  , J_value(nullptr)
  , normal_vectors(nullptr)
  , normal_x_jacobian(nullptr)
  , quadrature_weights(other.descriptor == nullptr ?
                         nullptr :
                         descriptor->quadrature_weights.begin())
  , cell(numbers::invalid_unsigned_int)
  , cell_type(internal::MatrixFreeFunctions::general)
  , is_interior_face(other.is_interior_face)
  , dof_access_index(other.dof_access_index)
{
  // Create deep copy of mapped geometry for use in parallel...
  if (other.mapped_geometry.get() != nullptr)
    {
      mapped_geometry = std::make_shared<
        internal::MatrixFreeFunctions::
          MappingDataOnTheFly<dim, Number, VectorizedArrayType>>(
        other.mapped_geometry->get_fe_values().get_mapping(),
        other.mapped_geometry->get_quadrature(),
        other.mapped_geometry->get_fe_values().get_update_flags());
      mapping_data = &mapped_geometry->get_data_storage();
      cell         = 0;

      jacobian = mapped_geometry->get_data_storage().jacobians[0].begin();
      J_value  = mapped_geometry->get_data_storage().JxW_values.begin();
    }
}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline FEEvaluationBaseData<dim, Number, is_face, VectorizedArrayType> &
FEEvaluationBaseData<dim, Number, is_face, VectorizedArrayType>::operator=(
  const FEEvaluationBaseData<dim, Number, is_face, VectorizedArrayType> &other)
{
  AssertDimension(quad_no, other.quad_no);
  AssertDimension(active_fe_index, other.active_fe_index);
  AssertDimension(active_quad_index, other.active_quad_index);

  // release old memory
  if (matrix_info == nullptr)
    {
      delete data;
      delete scratch_data_array;
    }
  else
    {
      matrix_info->release_scratch_data(scratch_data_array);
    }

  matrix_info  = other.matrix_info;
  dof_info     = other.dof_info;
  descriptor   = other.descriptor;
  mapping_data = other.mapping_data;
  if (other.matrix_info == nullptr)
    {
      data = new internal::MatrixFreeFunctions::ShapeInfo<VectorizedArrayType>(
        *other.data);
      scratch_data_array = new AlignedVector<VectorizedArrayType>();
    }
  else
    {
      data               = other.data;
      scratch_data_array = matrix_info->acquire_scratch_data();
    }

  quadrature_weights =
    (descriptor != nullptr ? descriptor->quadrature_weights.begin() : nullptr);
  cell             = numbers::invalid_unsigned_int;
  cell_type        = internal::MatrixFreeFunctions::general;
  is_interior_face = other.is_interior_face;
  dof_access_index = other.dof_access_index;

  // Create deep copy of mapped geometry for use in parallel...
  if (other.mapped_geometry.get() != nullptr)
    {
      mapped_geometry = std::make_shared<
        internal::MatrixFreeFunctions::
          MappingDataOnTheFly<dim, Number, VectorizedArrayType>>(
        other.mapped_geometry->get_fe_values().get_mapping(),
        other.mapped_geometry->get_quadrature(),
        other.mapped_geometry->get_fe_values().get_update_flags());
      cell         = 0;
      mapping_data = &mapped_geometry->get_data_storage();
      jacobian     = mapped_geometry->get_data_storage().jacobians[0].begin();
      J_value      = mapped_geometry->get_data_storage().JxW_values.begin();
    }

  return *this;
}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline FEEvaluationBaseData<dim, Number, is_face, VectorizedArrayType>::
  ~FEEvaluationBaseData()
{
  if (matrix_info != nullptr)
    {
      try
        {
          matrix_info->release_scratch_data(scratch_data_array);
        }
      catch (...)
        {}
    }
  else
    {
      delete scratch_data_array;
      delete data;
      data = nullptr;
    }
  scratch_data_array = nullptr;
}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline unsigned int
FEEvaluationBaseData<dim, Number, is_face, VectorizedArrayType>::
  get_mapping_data_index_offset() const
{
  if (matrix_info == nullptr)
    return 0;
  else
    {
      AssertIndexRange(cell, this->mapping_data->data_index_offsets.size());
      return this->mapping_data->data_index_offsets[cell];
    }
}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline internal::MatrixFreeFunctions::GeometryType
FEEvaluationBaseData<dim, Number, is_face, VectorizedArrayType>::get_cell_type()
  const
{
  Assert(cell != numbers::invalid_unsigned_int, ExcNotInitialized());
  return cell_type;
}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline const internal::MatrixFreeFunctions::ShapeInfo<VectorizedArrayType> &
FEEvaluationBaseData<dim, Number, is_face, VectorizedArrayType>::
  get_shape_info() const
{
  Assert(data != nullptr, ExcInternalError());
  return *data;
}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline const internal::MatrixFreeFunctions::DoFInfo &
FEEvaluationBaseData<dim, Number, is_face, VectorizedArrayType>::get_dof_info()
  const
{
  Assert(dof_info != nullptr,
         ExcMessage(
           "FEEvaluation was not initialized with a MatrixFree object!"));
  return *dof_info;
}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline DEAL_II_ALWAYS_INLINE Tensor<1, dim, VectorizedArrayType>
                             FEEvaluationBaseData<dim, Number, is_face, VectorizedArrayType>::
  get_normal_vector(const unsigned int q_point) const
{
  AssertIndexRange(q_point, n_quadrature_points);
  Assert(normal_vectors != nullptr,
         internal::ExcMatrixFreeAccessToUninitializedMappingField(
           "update_normal_vectors"));
  if (this->cell_type <= internal::MatrixFreeFunctions::flat_faces)
    return normal_vectors[0];
  else
    return normal_vectors[q_point];
}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline DEAL_II_ALWAYS_INLINE VectorizedArrayType
                             FEEvaluationBaseData<dim, Number, is_face, VectorizedArrayType>::JxW(
  const unsigned int q_point) const
{
  AssertIndexRange(q_point, n_quadrature_points);
  Assert(J_value != nullptr,
         internal::ExcMatrixFreeAccessToUninitializedMappingField(
           "update_values|update_gradients"));
  if (this->cell_type <= internal::MatrixFreeFunctions::affine)
    {
      Assert(this->quadrature_weights != nullptr, ExcInternalError());
      return J_value[0] * this->quadrature_weights[q_point];
    }
  else
    return J_value[q_point];
}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline DEAL_II_ALWAYS_INLINE Tensor<2, dim, VectorizedArrayType>
                             FEEvaluationBaseData<dim, Number, is_face, VectorizedArrayType>::
  inverse_jacobian(const unsigned int q_point) const
{
  AssertIndexRange(q_point, n_quadrature_points);
  Assert(this->jacobian != nullptr,
         internal::ExcMatrixFreeAccessToUninitializedMappingField(
           "update_gradients"));
  if (this->cell_type <= internal::MatrixFreeFunctions::affine)
    return jacobian[0];
  else
    return jacobian[q_point];
}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline std::array<unsigned int, VectorizedArrayType::size()>
FEEvaluationBaseData<dim, Number, is_face, VectorizedArrayType>::get_cell_ids()
  const
{
  Assert(this->matrix_info != nullptr, ExcNotInitialized());

  const unsigned int                n_lanes = VectorizedArrayType::size();
  std::array<unsigned int, n_lanes> cells;

  // initialize array
  for (unsigned int i = 0; i < n_lanes; ++i)
    cells[i] = numbers::invalid_unsigned_int;

  if ((is_face == false) ||
      (is_face &&
       this->dof_access_index ==
         internal::MatrixFreeFunctions::DoFInfo::dof_access_cell &&
       this->is_interior_face))
    {
      // cell or interior face face (element-centric loop)
      for (unsigned int i = 0; i < n_lanes; ++i)
        cells[i] = cell * n_lanes + i;
    }
  else if (is_face &&
           this->dof_access_index ==
             internal::MatrixFreeFunctions::DoFInfo::dof_access_cell &&
           this->is_interior_face == false)
    {
      // exterior face (element-centric loop): for this case, we need to
      // look into the FaceInfo field that collects information from both
      // sides of a face once for the global mesh, and pick the face id that
      // is not the local one (cell_this).
      for (unsigned int i = 0; i < n_lanes; i++)
        {
          // compute actual (non vectorized) cell ID
          const unsigned int cell_this = this->cell * n_lanes + i;
          // compute face ID
          unsigned int face_index =
            this->matrix_info->get_cell_and_face_to_plain_faces()(this->cell,
                                                                  this->face_no,
                                                                  i);

          if (face_index == numbers::invalid_unsigned_int)
            continue; // invalid face ID: no neighbor on boundary

          // get cell ID on both sides of face
          auto cell_m = this->matrix_info->get_face_info(face_index / n_lanes)
                          .cells_interior[face_index % n_lanes];
          auto cell_p = this->matrix_info->get_face_info(face_index / n_lanes)
                          .cells_exterior[face_index % n_lanes];

          // compare the IDs with the given cell ID
          if (cell_m == cell_this)
            cells[i] = cell_p; // neighbor has the other ID
          else if (cell_p == cell_this)
            cells[i] = cell_m;
        }
    }
  else if (is_face)
    {
      // face-centric faces
      const unsigned int *cells_ =
        is_interior_face ?
          &this->matrix_info->get_face_info(cell).cells_interior[0] :
          &this->matrix_info->get_face_info(cell).cells_exterior[0];
      for (unsigned int i = 0; i < VectorizedArrayType::size(); ++i)
        if (cells_[i] != numbers::invalid_unsigned_int)
          cells[i] = cells_[i];
    }

  return cells;
}


namespace internal
{
  template <int dim,
            typename Number,
            bool is_face,
            typename VectorizedArrayType,
            typename VectorizedArrayType2,
            typename GlobalVectorType,
            typename FU>
  inline void
  process_cell_data(
    const FEEvaluationBaseData<dim, Number, is_face, VectorizedArrayType> &phi,
    const MatrixFree<dim, Number, VectorizedArrayType> *matrix_info,
    GlobalVectorType &                                  array,
    VectorizedArrayType2 &                              out,
    const FU &                                          fu)
  {
    (void)matrix_info;
    Assert(matrix_info != nullptr, ExcNotImplemented());
    AssertDimension(array.size(),
                    matrix_info->get_task_info().cell_partition_data.back());

    // 1) collect ids of cell
    const auto cells = phi.get_cell_ids();

    // 2) actually gather values
    for (unsigned int i = 0; i < VectorizedArrayType::size(); ++i)
      if (cells[i] != numbers::invalid_unsigned_int)
        fu(out[i],
           array[cells[i] / VectorizedArrayType::size()]
                [cells[i] % VectorizedArrayType::size()]);
  }
} // namespace internal



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
std::array<unsigned int, VectorizedArrayType::size()>
FEEvaluationBaseData<dim, Number, is_face, VectorizedArrayType>::
  get_cell_or_face_ids() const
{
  const unsigned int v_len = VectorizedArrayType::size();
  std::array<unsigned int, VectorizedArrayType::size()> cells;

  // initialize array
  for (unsigned int i = 0; i < v_len; ++i)
    cells[i] = numbers::invalid_unsigned_int;

  if (is_face &&
      this->dof_access_index ==
        internal::MatrixFreeFunctions::DoFInfo::dof_access_cell &&
      this->is_interior_face == false)
    {
      // cell-based face-loop: plus face
      for (unsigned int i = 0; i < v_len; i++)
        {
          // compute actual (non vectorized) cell ID
          const unsigned int cell_this = this->cell * v_len + i;
          // compute face ID
          unsigned int fn =
            this->matrix_info->get_cell_and_face_to_plain_faces()(this->cell,
                                                                  this->face_no,
                                                                  i);

          if (fn == numbers::invalid_unsigned_int)
            continue; // invalid face ID: no neighbor on boundary

          // get cell ID on both sides of face
          auto cell_m = this->matrix_info->get_face_info(fn / v_len)
                          .cells_interior[fn % v_len];
          auto cell_p = this->matrix_info->get_face_info(fn / v_len)
                          .cells_exterior[fn % v_len];

          // compare the IDs with the given cell ID
          if (cell_m == cell_this)
            cells[i] = cell_p; // neighbor has the other ID
          else if (cell_p == cell_this)
            cells[i] = cell_m;
        }
    }
  else
    {
      for (unsigned int i = 0; i < v_len; ++i)
        cells[i] = cell * v_len + i;
    }

  return cells;
}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline VectorizedArrayType
FEEvaluationBaseData<dim, Number, is_face, VectorizedArrayType>::read_cell_data(
  const AlignedVector<VectorizedArrayType> &array) const
{
  VectorizedArrayType out = Number(1.);
  internal::process_cell_data(
    *this, this->matrix_info, array, out, [](auto &local, const auto &global) {
      local = global;
    });
  return out;
}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline void
FEEvaluationBaseData<dim, Number, is_face, VectorizedArrayType>::set_cell_data(
  AlignedVector<VectorizedArrayType> &array,
  const VectorizedArrayType &         in) const
{
  internal::process_cell_data(
    *this, this->matrix_info, array, in, [](const auto &local, auto &global) {
      global = local;
    });
}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
template <typename T>
inline std::array<T, VectorizedArrayType::size()>
FEEvaluationBaseData<dim, Number, is_face, VectorizedArrayType>::read_cell_data(
  const AlignedVector<std::array<T, VectorizedArrayType::size()>> &array) const
{
  std::array<T, VectorizedArrayType::size()> out;
  internal::process_cell_data(
    *this, this->matrix_info, array, out, [](auto &local, const auto &global) {
      local = global;
    });
  return out;
}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
template <typename T>
inline void
FEEvaluationBaseData<dim, Number, is_face, VectorizedArrayType>::set_cell_data(
  AlignedVector<std::array<T, VectorizedArrayType::size()>> &array,
  const std::array<T, VectorizedArrayType::size()> &         in) const
{
  internal::process_cell_data(
    *this, this->matrix_info, array, in, [](const auto &local, auto &global) {
      global = local;
    });
}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline const std::vector<unsigned int> &
FEEvaluationBaseData<dim, Number, is_face, VectorizedArrayType>::
  get_internal_dof_numbering() const
{
  return data->lexicographic_numbering;
}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline ArrayView<VectorizedArrayType>
FEEvaluationBaseData<dim, Number, is_face, VectorizedArrayType>::
  get_scratch_data() const
{
  return ArrayView<VectorizedArrayType>(
    const_cast<VectorizedArrayType *>(scratch_data),
    scratch_data_array->end() - scratch_data);
}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline unsigned int
FEEvaluationBaseData<dim, Number, is_face, VectorizedArrayType>::
  get_quadrature_index() const
{
  return this->quad_no;
}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline unsigned int
FEEvaluationBaseData<dim, Number, is_face, VectorizedArrayType>::
  get_current_cell_index() const
{
  if (is_face && this->dof_access_index ==
                   internal::MatrixFreeFunctions::DoFInfo::dof_access_cell)
    return this->cell * GeometryInfo<dim>::faces_per_cell + this->face_no;
  else
    return this->cell;
}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline unsigned int
FEEvaluationBaseData<dim, Number, is_face, VectorizedArrayType>::
  get_active_fe_index() const
{
  return active_fe_index;
}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline unsigned int
FEEvaluationBaseData<dim, Number, is_face, VectorizedArrayType>::
  get_active_quadrature_index() const
{
  return active_quad_index;
}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline const MatrixFree<dim, Number, VectorizedArrayType> &
FEEvaluationBaseData<dim, Number, is_face, VectorizedArrayType>::
  get_matrix_free() const
{
  Assert(matrix_info != nullptr,
         ExcMessage(
           "FEEvaluation was not initialized with a MatrixFree object!"));
  return *matrix_info;
}


 /*----------------------- FEEvaluationBase ----------------------------------*/ 

template <int dim,
          int n_components_,
          typename Number,
          bool is_face,
          typename VectorizedArrayType>
inline FEEvaluationBase<dim,
                        n_components_,
                        Number,
                        is_face,
                        VectorizedArrayType>::
  FEEvaluationBase(const MatrixFree<dim, Number, VectorizedArrayType> &data_in,
                   const unsigned int                                  dof_no,
                   const unsigned int first_selected_component,
                   const unsigned int quad_no_in,
                   const unsigned int fe_degree,
                   const unsigned int n_q_points,
                   const bool         is_interior_face,
                   const unsigned int active_fe_index,
                   const unsigned int active_quad_index,
                   const unsigned int face_type)
  : FEEvaluationBaseData<dim, Number, is_face, VectorizedArrayType>(
      data_in,
      dof_no,
      first_selected_component,
      quad_no_in,
      fe_degree,
      n_q_points,
      is_interior_face,
      active_fe_index,
      active_quad_index,
      face_type)
  , n_fe_components(data_in.get_dof_info(dof_no).start_components.back())
  , dof_values_initialized(false)
  , values_quad_initialized(false)
  , gradients_quad_initialized(false)
  , hessians_quad_initialized(false)
  , values_quad_submitted(false)
  , gradients_quad_submitted(false)
  , first_selected_component(first_selected_component)
{
  set_data_pointers();
  Assert(
    this->dof_info->start_components.back() == 1 ||
      static_cast<int>(n_components_) <=
        static_cast<int>(
          this->dof_info->start_components
            [this->dof_info->component_to_base_index[first_selected_component] +
             1]) -
          first_selected_component,
    ExcMessage(
      "You tried to construct a vector-valued evaluator with " +
      std::to_string(n_components) +
      " components. However, "
      "the current base element has only " +
      std::to_string(
        this->dof_info->start_components
          [this->dof_info->component_to_base_index[first_selected_component] +
           1] -
        first_selected_component) +
      " components left when starting from local element index " +
      std::to_string(
        first_selected_component -
        this->dof_info->start_components
          [this->dof_info->component_to_base_index[first_selected_component]]) +
      " (global index " + std::to_string(first_selected_component) + ")"));

  // do not check for correct dimensions of data fields here, should be done
  // in derived classes
}



template <int dim,
          int n_components_,
          typename Number,
          bool is_face,
          typename VectorizedArrayType>
inline FEEvaluationBase<dim,
                        n_components_,
                        Number,
                        is_face,
                        VectorizedArrayType>::
  FEEvaluationBase(
    const Mapping<dim> &      mapping,
    const FiniteElement<dim> &fe,
    const Quadrature<1> &     quadrature,
    const UpdateFlags         update_flags,
    const unsigned int        first_selected_component,
    const FEEvaluationBaseData<dim, Number, is_face, VectorizedArrayType>
      *other)
  : FEEvaluationBaseData<dim, Number, is_face, VectorizedArrayType>(
      mapping,
      fe,
      quadrature,
      update_flags,
      first_selected_component,
      other)
  , n_fe_components(n_components_)
  , dof_values_initialized(false)
  , values_quad_initialized(false)
  , gradients_quad_initialized(false)
  , hessians_quad_initialized(false)
  , values_quad_submitted(false)
  , gradients_quad_submitted(false)
  // keep the number of the selected component within the current base element
  // for reading dof values
  , first_selected_component(first_selected_component)
{
  set_data_pointers();

  const unsigned int base_element_number =
    fe.component_to_base_index(first_selected_component).first;
  Assert(fe.element_multiplicity(base_element_number) == 1 ||
           fe.element_multiplicity(base_element_number) -
               first_selected_component >=
             n_components_,
         ExcMessage("The underlying element must at least contain as many "
                    "components as requested by this class"));
  (void)base_element_number;
}



template <int dim,
          int n_components_,
          typename Number,
          bool is_face,
          typename VectorizedArrayType>
inline FEEvaluationBase<dim,
                        n_components_,
                        Number,
                        is_face,
                        VectorizedArrayType>::
  FEEvaluationBase(const FEEvaluationBase<dim,
                                          n_components_,
                                          Number,
                                          is_face,
                                          VectorizedArrayType> &other)
  : FEEvaluationBaseData<dim, Number, is_face, VectorizedArrayType>(other)
  , n_fe_components(other.n_fe_components)
  , dof_values_initialized(false)
  , values_quad_initialized(false)
  , gradients_quad_initialized(false)
  , hessians_quad_initialized(false)
  , values_quad_submitted(false)
  , gradients_quad_submitted(false)
  , first_selected_component(other.first_selected_component)
{
  set_data_pointers();
}



template <int dim,
          int n_components_,
          typename Number,
          bool is_face,
          typename VectorizedArrayType>
inline FEEvaluationBase<dim,
                        n_components_,
                        Number,
                        is_face,
                        VectorizedArrayType> &
FEEvaluationBase<dim, n_components_, Number, is_face, VectorizedArrayType>::
operator=(const FEEvaluationBase<dim,
                                 n_components_,
                                 Number,
                                 is_face,
                                 VectorizedArrayType> &other)
{
  this->FEEvaluationBaseData<dim, Number, is_face, VectorizedArrayType>::
  operator=(other);
  AssertDimension(n_fe_components, other.n_fe_components);
  AssertDimension(first_selected_component, other.first_selected_component);

  return *this;
}



template <int dim,
          int n_components_,
          typename Number,
          bool is_face,
          typename VectorizedArrayType>
inline void
FEEvaluationBase<dim, n_components_, Number, is_face, VectorizedArrayType>::
  set_data_pointers()
{
  Assert(this->scratch_data_array != nullptr, ExcInternalError());

  const unsigned int tensor_dofs_per_component =
    Utilities::fixed_power<dim>(this->data->data.front().fe_degree + 1);
  const unsigned int dofs_per_component =
    this->data->dofs_per_component_on_cell;
  const unsigned int n_quadrature_points = this->n_quadrature_points;

  const unsigned int shift =
    std::max(tensor_dofs_per_component + 1, dofs_per_component) *
      n_components_ * 3 +
    2 * n_quadrature_points;
  const unsigned int allocated_size =
    shift + n_components_ * dofs_per_component +
    (n_components_ * ((dim * (dim + 1)) / 2 + dim + 1) * n_quadrature_points);
  this->scratch_data_array->resize_fast(allocated_size);

  // set the pointers to the correct position in the data array
  for (unsigned int c = 0; c < n_components_; ++c)
    {
      values_dofs[c] =
        this->scratch_data_array->begin() + c * dofs_per_component;
    }
  values_quad =
    this->scratch_data_array->begin() + n_components * dofs_per_component;
  gradients_quad = this->scratch_data_array->begin() +
                   n_components * (dofs_per_component + n_quadrature_points);
  hessians_quad =
    this->scratch_data_array->begin() +
    n_components * (dofs_per_component + (dim + 1) * n_quadrature_points);
  this->scratch_data =
    this->scratch_data_array->begin() + n_components_ * dofs_per_component +
    (n_components_ * ((dim * (dim + 1)) / 2 + dim + 1) * n_quadrature_points);
}



namespace internal
{
  // allows to select between block vectors and non-block vectors, which
  // allows to use a unified interface for extracting blocks on block vectors
  // and doing nothing on usual vectors
  template <typename VectorType, bool>
  struct BlockVectorSelector
  {};

  template <typename VectorType>
  struct BlockVectorSelector<VectorType, true>
  {
    using BaseVectorType = typename VectorType::BlockType;

    static BaseVectorType *
    get_vector_component(VectorType &vec, const unsigned int component)
    {
      AssertIndexRange(component, vec.n_blocks());
      return &vec.block(component);
    }
  };

  template <typename VectorType>
  struct BlockVectorSelector<VectorType, false>
  {
    using BaseVectorType = VectorType;

    static BaseVectorType *
    get_vector_component(VectorType &vec, const unsigned int component)
    {
      // FEEvaluation allows to combine several vectors from a scalar
      // FiniteElement into a "vector-valued" FEEvaluation object with
      // multiple components. These components can be extracted with the other
      // get_vector_component functions. If we do not get a vector of vectors
      // (std::vector<VectorType>, std::vector<VectorType*>, BlockVector), we
      // must make sure that we do not duplicate the components in input
      // and/or duplicate the resulting integrals. In such a case, we should
      // only get the zeroth component in the vector contained set nullptr for
      // the others which allows us to catch unintended use in
      // read_write_operation.
      if (component == 0)
        return &vec;
      else
        return nullptr;
    }
  };

  template <typename VectorType>
  struct BlockVectorSelector<std::vector<VectorType>, false>
  {
    using BaseVectorType = VectorType;

    static BaseVectorType *
    get_vector_component(std::vector<VectorType> &vec,
                         const unsigned int       component)
    {
      AssertIndexRange(component, vec.size());
      return &vec[component];
    }
  };

  template <typename VectorType>
  struct BlockVectorSelector<std::vector<VectorType *>, false>
  {
    using BaseVectorType = VectorType;

    static BaseVectorType *
    get_vector_component(std::vector<VectorType *> &vec,
                         const unsigned int         component)
    {
      AssertIndexRange(component, vec.size());
      return vec[component];
    }
  };
} // namespace internal



template <int dim,
          int n_components_,
          typename Number,
          bool is_face,
          typename VectorizedArrayType>
template <typename VectorType, typename VectorOperation>
inline void
FEEvaluationBase<dim, n_components_, Number, is_face, VectorizedArrayType>::
  read_write_operation(
    const VectorOperation &                        operation,
    const std::array<VectorType *, n_components_> &src,
    const std::array<
      const std::vector<ArrayView<const typename VectorType::value_type>> *,
      n_components_> &                              src_sm,
    const std::bitset<VectorizedArrayType::size()> &mask,
    const bool                                      apply_constraints) const
{
  // Case 1: No MatrixFree object given, simple case because we do not need to
  // process constraints and need not care about vectorization -> go to
  // separate function
  if (this->matrix_info == nullptr)
    {
      read_write_operation_global(operation, src);
      return;
    }

  Assert(this->dof_info != nullptr, ExcNotInitialized());
  Assert(this->matrix_info->indices_initialized() == true, ExcNotInitialized());
  if (n_fe_components == 1)
    for (unsigned int comp = 0; comp < n_components; ++comp)
      {
        Assert(src[comp] != nullptr,
               ExcMessage("The finite element underlying this FEEvaluation "
                          "object is scalar, but you requested " +
                          std::to_string(n_components) +
                          " components via the template argument in "
                          "FEEvaluation. In that case, you must pass an "
                          "std::vector<VectorType> or a BlockVector to " +
                          "read_dof_values and distribute_local_to_global."));
        internal::check_vector_compatibility(*src[comp], *this->dof_info);
      }
  else
    {
      internal::check_vector_compatibility(*src[0], *this->dof_info);
    }

  // Case 2: contiguous indices which use reduced storage of indices and can
  // use vectorized load/store operations -> go to separate function
  AssertIndexRange(
    this->cell,
    this->dof_info->index_storage_variants[this->dof_access_index].size());
  if (this->dof_info->index_storage_variants
        [is_face ? this->dof_access_index :
                   internal::MatrixFreeFunctions::DoFInfo::dof_access_cell]
        [this->cell] >=
      internal::MatrixFreeFunctions::DoFInfo::IndexStorageVariants::contiguous)
    {
      read_write_operation_contiguous(operation, src, src_sm, mask);
      return;
    }

  // Case 3: standard operation with one index per degree of freedom -> go on
  // here
  constexpr unsigned int n_lanes = VectorizedArrayType::size();
  Assert(mask.count() == n_lanes,
         ExcNotImplemented("Masking currently not implemented for "
                           "non-contiguous DoF storage"));

  std::integral_constant<bool,
                         internal::is_vectorizable<VectorType, Number>::value>
    vector_selector;

  const unsigned int dofs_per_component =
    this->data->dofs_per_component_on_cell;
  if (this->dof_info->index_storage_variants
        [is_face ? this->dof_access_index :
                   internal::MatrixFreeFunctions::DoFInfo::dof_access_cell]
        [this->cell] ==
      internal::MatrixFreeFunctions::DoFInfo::IndexStorageVariants::interleaved)
    {
      const unsigned int *dof_indices =
        this->dof_info->dof_indices_interleaved.data() +
        this->dof_info->row_starts[this->cell * n_fe_components * n_lanes]
          .first +
        this->dof_info
            ->component_dof_indices_offset[this->active_fe_index]
                                          [this->first_selected_component] *
          n_lanes;
      if (n_components == 1 || n_fe_components == 1)
        for (unsigned int i = 0; i < dofs_per_component;
             ++i, dof_indices += n_lanes)
          for (unsigned int comp = 0; comp < n_components; ++comp)
            operation.process_dof_gather(dof_indices,
                                         *src[comp],
                                         0,
                                         values_dofs[comp][i],
                                         vector_selector);
      else
        for (unsigned int comp = 0; comp < n_components; ++comp)
          for (unsigned int i = 0; i < dofs_per_component;
               ++i, dof_indices += n_lanes)
            operation.process_dof_gather(
              dof_indices, *src[0], 0, values_dofs[comp][i], vector_selector);
      return;
    }

  const unsigned int *  dof_indices[n_lanes];
  VectorizedArrayType **values_dofs =
    const_cast<VectorizedArrayType **>(&this->values_dofs[0]);

  // Assign the appropriate cell ids for face/cell case and get the pointers
  // to the dof indices of the cells on all lanes
  unsigned int        cells_copied[n_lanes];
  const unsigned int *cells;
  unsigned int        n_vectorization_actual =
    this->dof_info
      ->n_vectorization_lanes_filled[this->dof_access_index][this->cell];
  bool               has_constraints   = false;
  const unsigned int n_components_read = n_fe_components > 1 ? n_components : 1;
  if (is_face)
    {
      if (this->dof_access_index ==
          internal::MatrixFreeFunctions::DoFInfo::dof_access_cell)
        for (unsigned int v = 0; v < n_vectorization_actual; ++v)
          cells_copied[v] = this->cell * VectorizedArrayType::size() + v;
      cells =
        this->dof_access_index ==
            internal::MatrixFreeFunctions::DoFInfo::dof_access_cell ?
          &cells_copied[0] :
          (this->is_interior_face ?
             &this->matrix_info->get_face_info(this->cell).cells_interior[0] :
             &this->matrix_info->get_face_info(this->cell).cells_exterior[0]);
      for (unsigned int v = 0; v < n_vectorization_actual; ++v)
        {
          Assert(cells[v] < this->dof_info->row_starts.size() - 1,
                 ExcInternalError());
          const std::pair<unsigned int, unsigned int> *my_index_start =
            &this->dof_info->row_starts[cells[v] * n_fe_components +
                                        first_selected_component];

          // check whether any of the SIMD lanes has constraints, i.e., the
          // constraint indicator which is the second entry of row_starts
          // increments on this cell
          if (my_index_start[n_components_read].second !=
              my_index_start[0].second)
            has_constraints = true;

          dof_indices[v] =
            this->dof_info->dof_indices.data() + my_index_start[0].first;
        }
      for (unsigned int v = n_vectorization_actual; v < n_lanes; ++v)
        dof_indices[v] = nullptr;
    }
  else
    {
      AssertIndexRange((this->cell + 1) * n_lanes * n_fe_components,
                       this->dof_info->row_starts.size());
      for (unsigned int v = 0; v < n_vectorization_actual; ++v)
        {
          const std::pair<unsigned int, unsigned int> *my_index_start =
            &this->dof_info
               ->row_starts[(this->cell * n_lanes + v) * n_fe_components +
                            first_selected_component];
          if (my_index_start[n_components_read].second !=
              my_index_start[0].second)
            has_constraints = true;
          Assert(my_index_start[n_components_read].first ==
                     my_index_start[0].first ||
                   my_index_start[0].first < this->dof_info->dof_indices.size(),
                 ExcIndexRange(0,
                               my_index_start[0].first,
                               this->dof_info->dof_indices.size()));
          dof_indices[v] =
            this->dof_info->dof_indices.data() + my_index_start[0].first;
        }
      for (unsigned int v = n_vectorization_actual; v < n_lanes; ++v)
        dof_indices[v] = nullptr;
    }

  // Case where we have no constraints throughout the whole cell: Can go
  // through the list of DoFs directly
  if (!has_constraints)
    {
      if (n_vectorization_actual < n_lanes)
        for (unsigned int comp = 0; comp < n_components; ++comp)
          for (unsigned int i = 0; i < dofs_per_component; ++i)
            operation.process_empty(values_dofs[comp][i]);
      if (n_components == 1 || n_fe_components == 1)
        {
          for (unsigned int v = 0; v < n_vectorization_actual; ++v)
            for (unsigned int i = 0; i < dofs_per_component; ++i)
              for (unsigned int comp = 0; comp < n_components; ++comp)
                operation.process_dof(dof_indices[v][i],
                                      *src[comp],
                                      values_dofs[comp][i][v]);
        }
      else
        {
          for (unsigned int comp = 0; comp < n_components; ++comp)
            for (unsigned int v = 0; v < n_vectorization_actual; ++v)
              for (unsigned int i = 0; i < dofs_per_component; ++i)
                operation.process_dof(
                  dof_indices[v][comp * dofs_per_component + i],
                  *src[0],
                  values_dofs[comp][i][v]);
        }
      return;
    }

  // In the case where there are some constraints to be resolved, loop over
  // all vector components that are filled and then over local dofs. ind_local
  // holds local number on cell, index iterates over the elements of
  // index_local_to_global and dof_indices points to the global indices stored
  // in index_local_to_global
  if (n_vectorization_actual < n_lanes)
    for (unsigned int comp = 0; comp < n_components; ++comp)
      for (unsigned int i = 0; i < dofs_per_component; ++i)
        operation.process_empty(values_dofs[comp][i]);
  for (unsigned int v = 0; v < n_vectorization_actual; ++v)
    {
      const unsigned int cell_index =
        is_face ? cells[v] : this->cell * n_lanes + v;
      const unsigned int cell_dof_index =
        cell_index * n_fe_components + first_selected_component;
      const unsigned int n_components_read =
        n_fe_components > 1 ? n_components : 1;
      unsigned int index_indicators =
        this->dof_info->row_starts[cell_dof_index].second;
      unsigned int next_index_indicators =
        this->dof_info->row_starts[cell_dof_index + 1].second;

      // For read_dof_values_plain, redirect the dof_indices field to the
      // unconstrained indices
      if (apply_constraints == false &&
          this->dof_info->row_starts[cell_dof_index].second !=
            this->dof_info->row_starts[cell_dof_index + n_components_read]
              .second)
        {
          Assert(this->dof_info->row_starts_plain_indices[cell_index] !=
                   numbers::invalid_unsigned_int,
                 ExcNotInitialized());
          dof_indices[v] =
            this->dof_info->plain_dof_indices.data() +
            this->dof_info
              ->component_dof_indices_offset[this->active_fe_index]
                                            [this->first_selected_component] +
            this->dof_info->row_starts_plain_indices[cell_index];
          next_index_indicators = index_indicators;
        }

      if (n_components == 1 || n_fe_components == 1)
        {
          unsigned int ind_local = 0;
          for (; index_indicators != next_index_indicators; ++index_indicators)
            {
              const std::pair<unsigned short, unsigned short> indicator =
                this->dof_info->constraint_indicator[index_indicators];
              // run through values up to next constraint
              for (unsigned int j = 0; j < indicator.first; ++j)
                for (unsigned int comp = 0; comp < n_components; ++comp)
                  operation.process_dof(dof_indices[v][j],
                                        *src[comp],
                                        values_dofs[comp][ind_local + j][v]);

              ind_local += indicator.first;
              dof_indices[v] += indicator.first;

              // constrained case: build the local value as a linear
              // combination of the global value according to constraints
              Number value[n_components];
              for (unsigned int comp = 0; comp < n_components; ++comp)
                operation.pre_constraints(values_dofs[comp][ind_local][v],
                                          value[comp]);

              const Number *data_val =
                this->matrix_info->constraint_pool_begin(indicator.second);
              const Number *end_pool =
                this->matrix_info->constraint_pool_end(indicator.second);
              for (; data_val != end_pool; ++data_val, ++dof_indices[v])
                for (unsigned int comp = 0; comp < n_components; ++comp)
                  operation.process_constraint(*dof_indices[v],
                                               *data_val,
                                               *src[comp],
                                               value[comp]);

              for (unsigned int comp = 0; comp < n_components; ++comp)
                operation.post_constraints(value[comp],
                                           values_dofs[comp][ind_local][v]);
              ind_local++;
            }

          AssertIndexRange(ind_local, dofs_per_component + 1);

          for (; ind_local < dofs_per_component; ++dof_indices[v], ++ind_local)
            for (unsigned int comp = 0; comp < n_components; ++comp)
              operation.process_dof(*dof_indices[v],
                                    *src[comp],
                                    values_dofs[comp][ind_local][v]);
        }
      else
        {
          // case with vector-valued finite elements where all components are
          // included in one single vector. Assumption: first come all entries
          // to the first component, then all entries to the second one, and
          // so on. This is ensured by the way MatrixFree reads out the
          // indices.
          for (unsigned int comp = 0; comp < n_components; ++comp)
            {
              unsigned int ind_local = 0;

              // check whether there is any constraint on the current cell
              for (; index_indicators != next_index_indicators;
                   ++index_indicators)
                {
                  const std::pair<unsigned short, unsigned short> indicator =
                    this->dof_info->constraint_indicator[index_indicators];

                  // run through values up to next constraint
                  for (unsigned int j = 0; j < indicator.first; ++j)
                    operation.process_dof(dof_indices[v][j],
                                          *src[0],
                                          values_dofs[comp][ind_local + j][v]);
                  ind_local += indicator.first;
                  dof_indices[v] += indicator.first;

                  // constrained case: build the local value as a linear
                  // combination of the global value according to constraints
                  Number value;
                  operation.pre_constraints(values_dofs[comp][ind_local][v],
                                            value);

                  const Number *data_val =
                    this->matrix_info->constraint_pool_begin(indicator.second);
                  const Number *end_pool =
                    this->matrix_info->constraint_pool_end(indicator.second);

                  for (; data_val != end_pool; ++data_val, ++dof_indices[v])
                    operation.process_constraint(*dof_indices[v],
                                                 *data_val,
                                                 *src[0],
                                                 value);

                  operation.post_constraints(value,
                                             values_dofs[comp][ind_local][v]);
                  ind_local++;
                }

              AssertIndexRange(ind_local, dofs_per_component + 1);

              // get the dof values past the last constraint
              for (; ind_local < dofs_per_component;
                   ++dof_indices[v], ++ind_local)
                {
                  AssertIndexRange(*dof_indices[v], src[0]->size());
                  operation.process_dof(*dof_indices[v],
                                        *src[0],
                                        values_dofs[comp][ind_local][v]);
                }

              if (apply_constraints == true && comp + 1 < n_components)
                next_index_indicators =
                  this->dof_info->row_starts[cell_dof_index + comp + 2].second;
            }
        }
    }
}



template <int dim,
          int n_components_,
          typename Number,
          bool is_face,
          typename VectorizedArrayType>
template <typename VectorType, typename VectorOperation>
inline void
FEEvaluationBase<dim, n_components_, Number, is_face, VectorizedArrayType>::
  read_write_operation_global(
    const VectorOperation &                        operation,
    const std::array<VectorType *, n_components_> &src) const
{
  Assert(!local_dof_indices.empty(), ExcNotInitialized());

  unsigned int index =
    first_selected_component * this->data->dofs_per_component_on_cell;
  for (unsigned int comp = 0; comp < n_components; ++comp)
    {
      for (unsigned int i = 0; i < this->data->dofs_per_component_on_cell;
           ++i, ++index)
        {
          operation.process_empty(values_dofs[comp][i]);
          operation.process_dof_global(
            local_dof_indices[this->data->lexicographic_numbering[index]],
            *src[0],
            values_dofs[comp][i][0]);
        }
    }
}



template <int dim,
          int n_components_,
          typename Number,
          bool is_face,
          typename VectorizedArrayType>
template <typename VectorType, typename VectorOperation>
inline void
FEEvaluationBase<dim, n_components_, Number, is_face, VectorizedArrayType>::
  read_write_operation_contiguous(
    const VectorOperation &                        operation,
    const std::array<VectorType *, n_components_> &src,
    const std::array<
      const std::vector<ArrayView<const typename VectorType::value_type>> *,
      n_components_> &                              vectors_sm,
    const std::bitset<VectorizedArrayType::size()> &mask) const
{
  // This functions processes the functions read_dof_values,
  // distribute_local_to_global, and set_dof_values with the same code for
  // contiguous cell indices (DG case). The distinction between these three
  // cases is made by the input VectorOperation that either reads values from
  // a vector and puts the data into the local data field or write local data
  // into the vector. Certain operations are no-ops for the given use case.

  std::integral_constant<bool,
                         internal::is_vectorizable<VectorType, Number>::value>
                                                               vector_selector;
  const internal::MatrixFreeFunctions::DoFInfo::DoFAccessIndex ind =
    is_face ? this->dof_access_index :
              internal::MatrixFreeFunctions::DoFInfo::dof_access_cell;
  const unsigned int n_lanes = mask.count();

  const std::vector<unsigned int> &dof_indices_cont =
    this->dof_info->dof_indices_contiguous[ind];

  // Simple case: We have contiguous storage, so we can simply copy out the
  // data
  if ((this->dof_info->index_storage_variants[ind][this->cell] ==
         internal::MatrixFreeFunctions::DoFInfo::IndexStorageVariants::
           interleaved_contiguous &&
       n_lanes == VectorizedArrayType::size()) &&
      !(is_face &&
        this->dof_access_index ==
          internal::MatrixFreeFunctions::DoFInfo::dof_access_cell &&
        this->is_interior_face == false))
    {
      const unsigned int dof_index =
        dof_indices_cont[this->cell * VectorizedArrayType::size()] +
        this->dof_info->component_dof_indices_offset[this->active_fe_index]
                                                    [first_selected_component] *
          VectorizedArrayType::size();
      if (n_components == 1 || n_fe_components == 1)
        for (unsigned int comp = 0; comp < n_components; ++comp)
          operation.process_dofs_vectorized(
            this->data->dofs_per_component_on_cell,
            dof_index,
            *src[comp],
            values_dofs[comp],
            vector_selector);
      else
        operation.process_dofs_vectorized(
          this->data->dofs_per_component_on_cell * n_components,
          dof_index,
          *src[0],
          values_dofs[0],
          vector_selector);
      return;
    }

  std::array<unsigned int, VectorizedArrayType::size()> cells =
    this->get_cell_or_face_ids();

  // More general case: Must go through the components one by one and apply
  // some transformations
  const unsigned int n_filled_lanes =
    this->dof_info->n_vectorization_lanes_filled[ind][this->cell];

  const bool is_ecl =
    this->dof_access_index ==
      internal::MatrixFreeFunctions::DoFInfo::dof_access_cell &&
    this->is_interior_face == false;

  if (vectors_sm[0] != nullptr)
    {
      const auto compute_vector_ptrs = [&](const unsigned int comp) {
        std::array<typename VectorType::value_type *,
                   VectorizedArrayType::size()>
          vector_ptrs = {};

        for (unsigned int v = 0; v < n_filled_lanes; ++v)
          {
            Assert(cells[v] != numbers::invalid_unsigned_int,
                   ExcNotImplemented());
            Assert(ind < this->dof_info->dof_indices_contiguous_sm.size(),
                   ExcIndexRange(
                     ind, 0, this->dof_info->dof_indices_contiguous_sm.size()));
            Assert(cells[v] <
                     this->dof_info->dof_indices_contiguous_sm[ind].size(),
                   ExcIndexRange(
                     cells[v],
                     0,
                     this->dof_info->dof_indices_contiguous_sm[ind].size()));

            const auto &temp =
              this->dof_info->dof_indices_contiguous_sm[ind][cells[v]];

            if (temp.first != numbers::invalid_unsigned_int)
              vector_ptrs[v] = const_cast<typename VectorType::value_type *>(
                vectors_sm[comp]->operator[](temp.first).data() + temp.second +
                this->dof_info->component_dof_indices_offset
                  [this->active_fe_index][this->first_selected_component]);
            else
              vector_ptrs[v] = nullptr;
          }
        for (unsigned int v = n_filled_lanes; v < VectorizedArrayType::size();
             ++v)
          vector_ptrs[v] = nullptr;

        return vector_ptrs;
      };

      if (n_filled_lanes == VectorizedArrayType::size() &&
          n_lanes == VectorizedArrayType::size() && !is_ecl)
        {
          if (n_components == 1 || n_fe_components == 1)
            {
              for (unsigned int comp = 0; comp < n_components; ++comp)
                {
                  auto vector_ptrs = compute_vector_ptrs(comp);
                  operation.process_dofs_vectorized_transpose(
                    this->data->dofs_per_component_on_cell,
                    vector_ptrs,
                    values_dofs[comp],
                    vector_selector);
                }
            }
          else
            {
              auto vector_ptrs = compute_vector_ptrs(0);
              operation.process_dofs_vectorized_transpose(
                this->data->dofs_per_component_on_cell * n_components,
                vector_ptrs,
                &values_dofs[0][0],
                vector_selector);
            }
        }
      else
        for (unsigned int comp = 0; comp < n_components; ++comp)
          {
            auto vector_ptrs = compute_vector_ptrs(
              (n_components == 1 || n_fe_components == 1) ? comp : 0);

            for (unsigned int i = 0; i < this->data->dofs_per_component_on_cell;
                 ++i)
              operation.process_empty(values_dofs[comp][i]);

            if (n_components == 1 || n_fe_components == 1)
              {
                for (unsigned int v = 0; v < n_filled_lanes; ++v)
                  if (mask[v] == true)
                    for (unsigned int i = 0;
                         i < this->data->dofs_per_component_on_cell;
                         ++i)
                      operation.process_dof(vector_ptrs[v][i],
                                            values_dofs[comp][i][v]);
              }
            else
              {
                for (unsigned int v = 0; v < n_filled_lanes; ++v)
                  if (mask[v] == true)
                    for (unsigned int i = 0;
                         i < this->data->dofs_per_component_on_cell;
                         ++i)
                      operation.process_dof(
                        vector_ptrs[v]
                                   [i + comp * this->data
                                                 ->dofs_per_component_on_cell],
                        values_dofs[comp][i][v]);
              }
          }
      return;
    }

  unsigned int dof_indices[VectorizedArrayType::size()];

  for (unsigned int v = 0; v < n_filled_lanes; ++v)
    {
      Assert(cells[v] != numbers::invalid_unsigned_int, ExcNotImplemented());
      dof_indices[v] =
        dof_indices_cont[cells[v]] +
        this->dof_info
            ->component_dof_indices_offset[this->active_fe_index]
                                          [this->first_selected_component] *
          this->dof_info->dof_indices_interleave_strides[ind][cells[v]];
    }

  for (unsigned int v = n_filled_lanes; v < VectorizedArrayType::size(); ++v)
    dof_indices[v] = numbers::invalid_unsigned_int;

  // In the case with contiguous cell indices, we know that there are no
  // constraints and that the indices within each element are contiguous
  if (n_filled_lanes == VectorizedArrayType::size() &&
      n_lanes == VectorizedArrayType::size() && !is_ecl)
    {
      if (this->dof_info->index_storage_variants[ind][this->cell] ==
          internal::MatrixFreeFunctions::DoFInfo::IndexStorageVariants::
            contiguous)
        {
          if (n_components == 1 || n_fe_components == 1)
            for (unsigned int comp = 0; comp < n_components; ++comp)
              operation.process_dofs_vectorized_transpose(
                this->data->dofs_per_component_on_cell,
                dof_indices,
                *src[comp],
                values_dofs[comp],
                vector_selector);
          else
            operation.process_dofs_vectorized_transpose(
              this->data->dofs_per_component_on_cell * n_components,
              dof_indices,
              *src[0],
              &values_dofs[0][0],
              vector_selector);
        }
      else if (this->dof_info->index_storage_variants[ind][this->cell] ==
               internal::MatrixFreeFunctions::DoFInfo::IndexStorageVariants::
                 interleaved_contiguous_strided)
        {
          if (n_components == 1 || n_fe_components == 1)
            for (unsigned int i = 0; i < this->data->dofs_per_component_on_cell;
                 ++i)
              {
                for (unsigned int comp = 0; comp < n_components; ++comp)
                  operation.process_dof_gather(dof_indices,
                                               *src[comp],
                                               i * VectorizedArrayType::size(),
                                               values_dofs[comp][i],
                                               vector_selector);
              }
          else
            for (unsigned int comp = 0; comp < n_components; ++comp)
              for (unsigned int i = 0;
                   i < this->data->dofs_per_component_on_cell;
                   ++i)
                {
                  operation.process_dof_gather(
                    dof_indices,
                    *src[0],
                    (comp * this->data->dofs_per_component_on_cell + i) *
                      VectorizedArrayType::size(),
                    values_dofs[comp][i],
                    vector_selector);
                }
        }
      else
        {
          Assert(this->dof_info->index_storage_variants[ind][this->cell] ==
                   internal::MatrixFreeFunctions::DoFInfo::
                     IndexStorageVariants::interleaved_contiguous_mixed_strides,
                 ExcNotImplemented());
          const unsigned int *offsets =
            &this->dof_info->dof_indices_interleave_strides
               [ind][VectorizedArrayType::size() * this->cell];
          if (n_components == 1 || n_fe_components == 1)
            for (unsigned int i = 0; i < this->data->dofs_per_component_on_cell;
                 ++i)
              {
                for (unsigned int comp = 0; comp < n_components; ++comp)
                  operation.process_dof_gather(dof_indices,
                                               *src[comp],
                                               0,
                                               values_dofs[comp][i],
                                               vector_selector);
                DEAL_II_OPENMP_SIMD_PRAGMA
                for (unsigned int v = 0; v < VectorizedArrayType::size(); ++v)
                  dof_indices[v] += offsets[v];
              }
          else
            for (unsigned int comp = 0; comp < n_components; ++comp)
              for (unsigned int i = 0;
                   i < this->data->dofs_per_component_on_cell;
                   ++i)
                {
                  operation.process_dof_gather(dof_indices,
                                               *src[0],
                                               0,
                                               values_dofs[comp][i],
                                               vector_selector);
                  DEAL_II_OPENMP_SIMD_PRAGMA
                  for (unsigned int v = 0; v < VectorizedArrayType::size(); ++v)
                    dof_indices[v] += offsets[v];
                }
        }
    }
  else
    for (unsigned int comp = 0; comp < n_components; ++comp)
      {
        for (unsigned int i = 0; i < this->data->dofs_per_component_on_cell;
             ++i)
          operation.process_empty(values_dofs[comp][i]);
        if (this->dof_info->index_storage_variants[ind][this->cell] ==
            internal::MatrixFreeFunctions::DoFInfo::IndexStorageVariants::
              contiguous)
          {
            if (n_components == 1 || n_fe_components == 1)
              {
                for (unsigned int v = 0; v < n_filled_lanes; ++v)
                  if (mask[v] == true)
                    for (unsigned int i = 0;
                         i < this->data->dofs_per_component_on_cell;
                         ++i)
                      operation.process_dof(dof_indices[v] + i,
                                            *src[comp],
                                            values_dofs[comp][i][v]);
              }
            else
              {
                for (unsigned int v = 0; v < n_filled_lanes; ++v)
                  if (mask[v] == true)
                    for (unsigned int i = 0;
                         i < this->data->dofs_per_component_on_cell;
                         ++i)
                      operation.process_dof(
                        dof_indices[v] + i +
                          comp * this->data->dofs_per_component_on_cell,
                        *src[0],
                        values_dofs[comp][i][v]);
              }
          }
        else
          {
            const unsigned int *offsets =
              &this->dof_info->dof_indices_interleave_strides
                 [ind][VectorizedArrayType::size() * this->cell];
            for (unsigned int v = 0; v < n_filled_lanes; ++v)
              AssertIndexRange(offsets[v], VectorizedArrayType::size() + 1);
            if (n_components == 1 || n_fe_components == 1)
              for (unsigned int v = 0; v < n_filled_lanes; ++v)
                {
                  if (mask[v] == true)
                    for (unsigned int i = 0;
                         i < this->data->dofs_per_component_on_cell;
                         ++i)
                      operation.process_dof(dof_indices[v] + i * offsets[v],
                                            *src[comp],
                                            values_dofs[comp][i][v]);
                }
            else
              {
                for (unsigned int v = 0; v < n_filled_lanes; ++v)
                  if (mask[v] == true)
                    for (unsigned int i = 0;
                         i < this->data->dofs_per_component_on_cell;
                         ++i)
                      operation.process_dof(
                        dof_indices[v] +
                          (i + comp * this->data->dofs_per_component_on_cell) *
                            offsets[v],
                        *src[0],
                        values_dofs[comp][i][v]);
              }
          }
      }
}

namespace internal
{
  template <typename Number,
            typename VectorType,
            typename std::enable_if<!IsBlockVector<VectorType>::value,
                                    VectorType>::type * = nullptr>
  decltype(std::declval<VectorType>().begin())
  get_beginning(VectorType &vec)
  {
    return vec.begin();
  }

  template <typename Number,
            typename VectorType,
            typename std::enable_if<IsBlockVector<VectorType>::value,
                                    VectorType>::type * = nullptr>
  typename VectorType::value_type *
  get_beginning(VectorType &)
  {
    return nullptr;
  }

  template <typename VectorType,
            typename std::enable_if<has_shared_vector_data<VectorType>::value,
                                    VectorType>::type * = nullptr>
  const std::vector<ArrayView<const typename VectorType::value_type>> *
  get_shared_vector_data(VectorType &       vec,
                         const bool         is_valid_mode_for_sm,
                         const unsigned int active_fe_index,
                         const internal::MatrixFreeFunctions::DoFInfo *dof_info)
  {
    // note: no hp is supported
    if (is_valid_mode_for_sm &&
        dof_info->dof_indices_contiguous_sm[0  /*any index (<3) should work*/ ]
            .size() > 0 &&
        active_fe_index == 0)
      return &vec.shared_vector_data();
    else
      return nullptr;
  }

  template <typename VectorType,
            typename std::enable_if<!has_shared_vector_data<VectorType>::value,
                                    VectorType>::type * = nullptr>
  const std::vector<ArrayView<const typename VectorType::value_type>> *
  get_shared_vector_data(VectorType &,
                         const bool,
                         const unsigned int,
                         const internal::MatrixFreeFunctions::DoFInfo *)
  {
    return nullptr;
  }

  template <int n_components, typename VectorType>
  std::pair<
    std::array<typename internal::BlockVectorSelector<
                 typename std::remove_const<VectorType>::type,
                 IsBlockVector<typename std::remove_const<VectorType>::type>::
                   value>::BaseVectorType *,
               n_components>,
    std::array<
      const std::vector<ArrayView<const typename internal::BlockVectorSelector<
        typename std::remove_const<VectorType>::type,
        IsBlockVector<typename std::remove_const<VectorType>::type>::value>::
                                    BaseVectorType::value_type>> *,
      n_components>>
  get_vector_data(VectorType &       src,
                  const unsigned int first_index,
                  const bool         is_valid_mode_for_sm,
                  const unsigned int active_fe_index,
                  const internal::MatrixFreeFunctions::DoFInfo *dof_info)
  {
    // select between block vectors and non-block vectors. Note that the number
    // of components is checked in the internal data
    std::pair<
      std::array<typename internal::BlockVectorSelector<
                   typename std::remove_const<VectorType>::type,
                   IsBlockVector<typename std::remove_const<VectorType>::type>::
                     value>::BaseVectorType *,
                 n_components>,
      std::array<
        const std::vector<
          ArrayView<const typename internal::BlockVectorSelector<
            typename std::remove_const<VectorType>::type,
            IsBlockVector<typename std::remove_const<VectorType>::type>::
              value>::BaseVectorType::value_type>> *,
        n_components>>
      src_data;

    for (unsigned int d = 0; d < n_components; ++d)
      src_data.first[d] = internal::BlockVectorSelector<
        typename std::remove_const<VectorType>::type,
        IsBlockVector<typename std::remove_const<VectorType>::type>::value>::
        get_vector_component(
          const_cast<typename std::remove_const<VectorType>::type &>(src),
          d + first_index);

    for (unsigned int d = 0; d < n_components; ++d)
      src_data.second[d] = get_shared_vector_data(*src_data.first[d],
                                                  is_valid_mode_for_sm,
                                                  active_fe_index,
                                                  dof_info);

    return src_data;
  }
} // namespace internal



template <int dim,
          int n_components_,
          typename Number,
          bool is_face,
          typename VectorizedArrayType>
template <typename VectorType>
inline void
FEEvaluationBase<dim, n_components_, Number, is_face, VectorizedArrayType>::
  read_dof_values(const VectorType &src, const unsigned int first_index)
{
  const auto src_data = internal::get_vector_data<n_components_>(
    src,
    first_index,
    this->dof_access_index ==
      internal::MatrixFreeFunctions::DoFInfo::dof_access_cell,
    this->active_fe_index,
    this->dof_info);

  internal::VectorReader<Number, VectorizedArrayType> reader;
  read_write_operation(reader,
                       src_data.first,
                       src_data.second,
                       std::bitset<VectorizedArrayType::size()>().flip(),
                       true);

#  ifdef DEBUG
  dof_values_initialized = true;
#  endif
}



template <int dim,
          int n_components_,
          typename Number,
          bool is_face,
          typename VectorizedArrayType>
template <typename VectorType>
inline void
FEEvaluationBase<dim, n_components_, Number, is_face, VectorizedArrayType>::
  read_dof_values_plain(const VectorType &src, const unsigned int first_index)
{
  const auto src_data = internal::get_vector_data<n_components_>(
    src,
    first_index,
    this->dof_access_index ==
      internal::MatrixFreeFunctions::DoFInfo::dof_access_cell,
    this->active_fe_index,
    this->dof_info);

  internal::VectorReader<Number, VectorizedArrayType> reader;
  read_write_operation(reader,
                       src_data.first,
                       src_data.second,
                       std::bitset<VectorizedArrayType::size()>().flip(),
                       false);

#  ifdef DEBUG
  dof_values_initialized = true;
#  endif
}



template <int dim,
          int n_components_,
          typename Number,
          bool is_face,
          typename VectorizedArrayType>
template <typename VectorType>
inline void
FEEvaluationBase<dim, n_components_, Number, is_face, VectorizedArrayType>::
  distribute_local_to_global(
    VectorType &                                    dst,
    const unsigned int                              first_index,
    const std::bitset<VectorizedArrayType::size()> &mask) const
{
#  ifdef DEBUG
  Assert(dof_values_initialized == true,
         internal::ExcAccessToUninitializedField());
#  endif

  const auto dst_data = internal::get_vector_data<n_components_>(
    dst,
    first_index,
    this->dof_access_index ==
      internal::MatrixFreeFunctions::DoFInfo::dof_access_cell,
    this->active_fe_index,
    this->dof_info);

  internal::VectorDistributorLocalToGlobal<Number, VectorizedArrayType>
    distributor;
  read_write_operation(distributor, dst_data.first, dst_data.second, mask);
}



template <int dim,
          int n_components_,
          typename Number,
          bool is_face,
          typename VectorizedArrayType>
template <typename VectorType>
inline void
FEEvaluationBase<dim, n_components_, Number, is_face, VectorizedArrayType>::
  set_dof_values(VectorType &                                    dst,
                 const unsigned int                              first_index,
                 const std::bitset<VectorizedArrayType::size()> &mask) const
{
#  ifdef DEBUG
  Assert(dof_values_initialized == true,
         internal::ExcAccessToUninitializedField());
#  endif

  const auto dst_data = internal::get_vector_data<n_components_>(
    dst,
    first_index,
    this->dof_access_index ==
      internal::MatrixFreeFunctions::DoFInfo::dof_access_cell,
    this->active_fe_index,
    this->dof_info);

  internal::VectorSetter<Number, VectorizedArrayType> setter;
  read_write_operation(setter, dst_data.first, dst_data.second, mask);
}



template <int dim,
          int n_components_,
          typename Number,
          bool is_face,
          typename VectorizedArrayType>
template <typename VectorType>
inline void
FEEvaluationBase<dim, n_components_, Number, is_face, VectorizedArrayType>::
  set_dof_values_plain(
    VectorType &                                    dst,
    const unsigned int                              first_index,
    const std::bitset<VectorizedArrayType::size()> &mask) const
{
#  ifdef DEBUG
  Assert(dof_values_initialized == true,
         internal::ExcAccessToUninitializedField());
#  endif

  const auto dst_data = internal::get_vector_data<n_components_>(
    dst,
    first_index,
    this->dof_access_index ==
      internal::MatrixFreeFunctions::DoFInfo::dof_access_cell,
    this->active_fe_index,
    this->dof_info);

  internal::VectorSetter<Number, VectorizedArrayType> setter;
  read_write_operation(setter, dst_data.first, dst_data.second, mask, false);
}



 /*------------------------------ access to data fields ----------------------*/ 



template <int dim,
          int n_components,
          typename Number,
          bool is_face,
          typename VectorizedArrayType>
inline const VectorizedArrayType *
FEEvaluationBase<dim, n_components, Number, is_face, VectorizedArrayType>::
  begin_dof_values() const
{
  return &values_dofs[0][0];
}



template <int dim,
          int n_components,
          typename Number,
          bool is_face,
          typename VectorizedArrayType>
inline VectorizedArrayType *
FEEvaluationBase<dim, n_components, Number, is_face, VectorizedArrayType>::
  begin_dof_values()
{
#  ifdef DEBUG
  dof_values_initialized = true;
#  endif
  return &values_dofs[0][0];
}



template <int dim,
          int n_components,
          typename Number,
          bool is_face,
          typename VectorizedArrayType>
inline const VectorizedArrayType *
FEEvaluationBase<dim, n_components, Number, is_face, VectorizedArrayType>::
  begin_values() const
{
#  ifdef DEBUG
  Assert(values_quad_initialized || values_quad_submitted, ExcNotInitialized());
#  endif
  return values_quad;
}



template <int dim,
          int n_components,
          typename Number,
          bool is_face,
          typename VectorizedArrayType>
inline VectorizedArrayType *
FEEvaluationBase<dim, n_components, Number, is_face, VectorizedArrayType>::
  begin_values()
{
#  ifdef DEBUG
  values_quad_initialized = true;
  values_quad_submitted   = true;
#  endif
  return values_quad;
}



template <int dim,
          int n_components,
          typename Number,
          bool is_face,
          typename VectorizedArrayType>
inline const VectorizedArrayType *
FEEvaluationBase<dim, n_components, Number, is_face, VectorizedArrayType>::
  begin_gradients() const
{
#  ifdef DEBUG
  Assert(gradients_quad_initialized || gradients_quad_submitted,
         ExcNotInitialized());
#  endif
  return gradients_quad;
}



template <int dim,
          int n_components,
          typename Number,
          bool is_face,
          typename VectorizedArrayType>
inline VectorizedArrayType *
FEEvaluationBase<dim, n_components, Number, is_face, VectorizedArrayType>::
  begin_gradients()
{
#  ifdef DEBUG
  gradients_quad_submitted   = true;
  gradients_quad_initialized = true;
#  endif
  return gradients_quad;
}



template <int dim,
          int n_components,
          typename Number,
          bool is_face,
          typename VectorizedArrayType>
inline const VectorizedArrayType *
FEEvaluationBase<dim, n_components, Number, is_face, VectorizedArrayType>::
  begin_hessians() const
{
#  ifdef DEBUG
  Assert(hessians_quad_initialized, ExcNotInitialized());
#  endif
  return hessians_quad;
}



template <int dim,
          int n_components,
          typename Number,
          bool is_face,
          typename VectorizedArrayType>
inline VectorizedArrayType *
FEEvaluationBase<dim, n_components, Number, is_face, VectorizedArrayType>::
  begin_hessians()
{
#  ifdef DEBUG
  hessians_quad_initialized = true;
#  endif
  return hessians_quad;
}



template <int dim,
          int n_components,
          typename Number,
          bool is_face,
          typename VectorizedArrayType>
inline unsigned int
FEEvaluationBase<dim, n_components, Number, is_face, VectorizedArrayType>::
  get_first_selected_component() const
{
  return first_selected_component;
}



template <int dim,
          int n_components_,
          typename Number,
          bool is_face,
          typename VectorizedArrayType>
inline DEAL_II_ALWAYS_INLINE Tensor<1, n_components_, VectorizedArrayType>
                             FEEvaluationBase<dim, n_components_, Number, is_face, VectorizedArrayType>::
  get_dof_value(const unsigned int dof) const
{
  AssertIndexRange(dof, this->data->dofs_per_component_on_cell);
  Tensor<1, n_components_, VectorizedArrayType> return_value;
  for (unsigned int comp = 0; comp < n_components; comp++)
    return_value[comp] = this->values_dofs[comp][dof];
  return return_value;
}



template <int dim,
          int n_components_,
          typename Number,
          bool is_face,
          typename VectorizedArrayType>
inline DEAL_II_ALWAYS_INLINE Tensor<1, n_components_, VectorizedArrayType>
                             FEEvaluationBase<dim, n_components_, Number, is_face, VectorizedArrayType>::
  get_value(const unsigned int q_point) const
{
#  ifdef DEBUG
  Assert(this->values_quad_initialized == true,
         internal::ExcAccessToUninitializedField());
#  endif

  AssertIndexRange(q_point, this->n_quadrature_points);
  const std::size_t                             nqp = this->n_quadrature_points;
  Tensor<1, n_components_, VectorizedArrayType> return_value;
  for (unsigned int comp = 0; comp < n_components; comp++)
    return_value[comp] = values_quad[comp * nqp + q_point];
  return return_value;
}



template <int dim,
          int n_components_,
          typename Number,
          bool is_face,
          typename VectorizedArrayType>
inline DEAL_II_ALWAYS_INLINE
  Tensor<1, n_components_, Tensor<1, dim, VectorizedArrayType>>
  FEEvaluationBase<dim, n_components_, Number, is_face, VectorizedArrayType>::
    get_gradient(const unsigned int q_point) const
{
#  ifdef DEBUG
  Assert(this->gradients_quad_initialized == true,
         internal::ExcAccessToUninitializedField());
#  endif

  AssertIndexRange(q_point, this->n_quadrature_points);
  Assert(this->jacobian != nullptr,
         internal::ExcMatrixFreeAccessToUninitializedMappingField(
           "update_gradients"));
  const std::size_t nqp = this->n_quadrature_points;
  Tensor<1, n_components_, Tensor<1, dim, VectorizedArrayType>> grad_out;

  // Cartesian cell
  if (!is_face && this->cell_type == internal::MatrixFreeFunctions::cartesian)
    {
      for (unsigned int d = 0; d < dim; ++d)
        for (unsigned int comp = 0; comp < n_components; comp++)
          grad_out[comp][d] = gradients_quad[(comp * dim + d) * nqp + q_point] *
                              this->jacobian[0][d][d];
    }
  // cell with general/affine Jacobian
  else
    {
      const Tensor<2, dim, VectorizedArrayType> &jac =
        this->jacobian[this->cell_type > internal::MatrixFreeFunctions::affine ?
                         q_point :
                         0];
      for (unsigned int comp = 0; comp < n_components; comp++)
        for (unsigned int d = 0; d < dim; ++d)
          {
            grad_out[comp][d] =
              jac[d][0] * gradients_quad[(comp * dim) * nqp + q_point];
            for (unsigned int e = 1; e < dim; ++e)
              grad_out[comp][d] +=
                jac[d][e] * gradients_quad[(comp * dim + e) * nqp + q_point];
          }
    }
  return grad_out;
}



template <int dim,
          int n_components_,
          typename Number,
          bool is_face,
          typename VectorizedArrayType>
inline DEAL_II_ALWAYS_INLINE Tensor<1, n_components_, VectorizedArrayType>
                             FEEvaluationBase<dim, n_components_, Number, is_face, VectorizedArrayType>::
  get_normal_derivative(const unsigned int q_point) const
{
  AssertIndexRange(q_point, this->n_quadrature_points);
#  ifdef DEBUG
  Assert(this->gradients_quad_initialized == true,
         internal::ExcAccessToUninitializedField());
#  endif

  Assert(this->normal_x_jacobian != nullptr,
         internal::ExcMatrixFreeAccessToUninitializedMappingField(
           "update_gradients"));

  const std::size_t                            nqp = this->n_quadrature_points;
  Tensor<1, n_components, VectorizedArrayType> grad_out;

  if (this->cell_type == internal::MatrixFreeFunctions::cartesian)
    for (unsigned int comp = 0; comp < n_components; comp++)
      grad_out[comp] = gradients_quad[(comp * dim + dim - 1) * nqp + q_point] *
                       (this->normal_x_jacobian[0][dim - 1]);
  else
    {
      const std::size_t index =
        this->cell_type <= internal::MatrixFreeFunctions::affine ? 0 : q_point;
      for (unsigned int comp = 0; comp < n_components; comp++)
        {
          grad_out[comp] = gradients_quad[comp * dim * nqp + q_point] *
                           this->normal_x_jacobian[index][0];
          for (unsigned int d = 1; d < dim; ++d)
            grad_out[comp] += gradients_quad[(comp * dim + d) * nqp + q_point] *
                              this->normal_x_jacobian[index][d];
        }
    }
  return grad_out;
}



namespace internal
{
  // compute tmp = hess_unit(u) * J^T. do this manually because we do not
  // store the lower diagonal because of symmetry
  template <typename VectorizedArrayType>
  inline void
  hessian_unit_times_jac(const Tensor<2, 1, VectorizedArrayType> &jac,
                         const VectorizedArrayType *const         hessians,
                         const unsigned int,
                         VectorizedArrayType (&tmp)[1][1])
  {
    tmp[0][0] = jac[0][0] * hessians[0];
  }

  template <typename VectorizedArrayType>
  inline void
  hessian_unit_times_jac(const Tensor<2, 2, VectorizedArrayType> &jac,
                         const VectorizedArrayType *const         hessians,
                         const unsigned int                       nqp,
                         VectorizedArrayType (&tmp)[2][2])
  {
    for (unsigned int d = 0; d < 2; ++d)
      {
        tmp[0][d] = (jac[d][0] * hessians[0] + jac[d][1] * hessians[2 * nqp]);
        tmp[1][d] =
          (jac[d][0] * hessians[2 * nqp] + jac[d][1] * hessians[1 * nqp]);
      }
  }

  template <typename VectorizedArrayType>
  inline void
  hessian_unit_times_jac(const Tensor<2, 3, VectorizedArrayType> &jac,
                         const VectorizedArrayType *const         hessians,
                         const unsigned int                       nqp,
                         VectorizedArrayType (&tmp)[3][3])
  {
    for (unsigned int d = 0; d < 3; ++d)
      {
        tmp[0][d] =
          (jac[d][0] * hessians[0 * nqp] + jac[d][1] * hessians[3 * nqp] +
           jac[d][2] * hessians[4 * nqp]);
        tmp[1][d] =
          (jac[d][0] * hessians[3 * nqp] + jac[d][1] * hessians[1 * nqp] +
           jac[d][2] * hessians[5 * nqp]);
        tmp[2][d] =
          (jac[d][0] * hessians[4 * nqp] + jac[d][1] * hessians[5 * nqp] +
           jac[d][2] * hessians[2 * nqp]);
      }
  }
} // namespace internal



template <int dim,
          int n_components_,
          typename Number,
          bool is_face,
          typename VectorizedArrayType>
inline Tensor<1, n_components_, Tensor<2, dim, VectorizedArrayType>>
FEEvaluationBase<dim, n_components_, Number, is_face, VectorizedArrayType>::
  get_hessian(const unsigned int q_point) const
{
  Assert(!is_face, ExcNotImplemented());
#  ifdef DEBUG
  Assert(this->hessians_quad_initialized == true,
         internal::ExcAccessToUninitializedField());
#  endif
  AssertIndexRange(q_point, this->n_quadrature_points);

  Assert(this->jacobian != nullptr,
         internal::ExcMatrixFreeAccessToUninitializedMappingField(
           "update_hessian"));
  const Tensor<2, dim, VectorizedArrayType> &jac =
    this->jacobian[this->cell_type <= internal::MatrixFreeFunctions::affine ?
                     0 :
                     q_point];

  Tensor<1, n_components, Tensor<2, dim, VectorizedArrayType>> hessian_out;

  const std::size_t      nqp  = this->n_quadrature_points;
  constexpr unsigned int hdim = (dim * (dim + 1)) / 2;

  // Cartesian cell
  if (this->cell_type == internal::MatrixFreeFunctions::cartesian)
    {
      for (unsigned int comp = 0; comp < n_components; comp++)
        {
          for (unsigned int d = 0; d < dim; ++d)
            hessian_out[comp][d][d] =
              hessians_quad[(comp * hdim + d) * nqp + q_point] *
              (jac[d][d] * jac[d][d]);
          switch (dim)
            {
              case 1:
                break;
              case 2:
                hessian_out[comp][0][1] =
                  hessians_quad[(comp * hdim + 2) * nqp + q_point] *
                  (jac[0][0] * jac[1][1]);
                break;
              case 3:
                hessian_out[comp][0][1] =
                  hessians_quad[(comp * hdim + 3) * nqp + q_point] *
                  (jac[0][0] * jac[1][1]);
                hessian_out[comp][0][2] =
                  hessians_quad[(comp * hdim + 4) * nqp + q_point] *
                  (jac[0][0] * jac[2][2]);
                hessian_out[comp][1][2] =
                  hessians_quad[(comp * hdim + 5) * nqp + q_point] *
                  (jac[1][1] * jac[2][2]);
                break;
              default:
                Assert(false, ExcNotImplemented());
            }
          for (unsigned int d = 0; d < dim; ++d)
            for (unsigned int e = d + 1; e < dim; ++e)
              hessian_out[comp][e][d] = hessian_out[comp][d][e];
        }
    }
  // cell with general Jacobian, but constant within the cell
  else if (this->cell_type == internal::MatrixFreeFunctions::affine)
    {
      for (unsigned int comp = 0; comp < n_components; comp++)
        {
          VectorizedArrayType tmp[dim][dim];
          internal::hessian_unit_times_jac(
            jac, hessians_quad + comp * hdim * nqp + q_point, nqp, tmp);

          // compute first part of hessian, J * tmp = J * hess_unit(u) * J^T
          for (unsigned int d = 0; d < dim; ++d)
            for (unsigned int e = d; e < dim; ++e)
              {
                hessian_out[comp][d][e] = jac[d][0] * tmp[0][e];
                for (unsigned int f = 1; f < dim; ++f)
                  hessian_out[comp][d][e] += jac[d][f] * tmp[f][e];
              }

          // no J' * grad(u) part here because the Jacobian is constant
          // throughout the cell and hence, its derivative is zero

          // take symmetric part
          for (unsigned int d = 0; d < dim; ++d)
            for (unsigned int e = d + 1; e < dim; ++e)
              hessian_out[comp][e][d] = hessian_out[comp][d][e];
        }
    }
  // cell with general Jacobian
  else
    {
      const auto &jac_grad =
        this->mapping_data->jacobian_gradients
          [1 - this->is_interior_face]
          [this->mapping_data->data_index_offsets[this->cell] + q_point];
      for (unsigned int comp = 0; comp < n_components; comp++)
        {
          // compute laplacian before the gradient because it needs to access
          // unscaled gradient data
          VectorizedArrayType tmp[dim][dim];
          internal::hessian_unit_times_jac(
            jac, hessians_quad + comp * hdim * nqp + q_point, nqp, tmp);

          // compute first part of hessian, J * tmp = J * hess_unit(u) * J^T
          for (unsigned int d = 0; d < dim; ++d)
            for (unsigned int e = d; e < dim; ++e)
              {
                hessian_out[comp][d][e] = jac[d][0] * tmp[0][e];
                for (unsigned int f = 1; f < dim; ++f)
                  hessian_out[comp][d][e] += jac[d][f] * tmp[f][e];
              }

          // add diagonal part of J' * grad(u)
          for (unsigned int d = 0; d < dim; ++d)
            for (unsigned int e = 0; e < dim; ++e)
              hessian_out[comp][d][d] +=
                jac_grad[d][e] *
                gradients_quad[(comp * dim + e) * nqp + q_point];

          // add off-diagonal part of J' * grad(u)
          for (unsigned int d = 0, count = dim; d < dim; ++d)
            for (unsigned int e = d + 1; e < dim; ++e, ++count)
              for (unsigned int f = 0; f < dim; ++f)
                hessian_out[comp][d][e] +=
                  jac_grad[count][f] *
                  gradients_quad[(comp * dim + f) * nqp + q_point];

          // take symmetric part
          for (unsigned int d = 0; d < dim; ++d)
            for (unsigned int e = d + 1; e < dim; ++e)
              hessian_out[comp][e][d] = hessian_out[comp][d][e];
        }
    }
  return hessian_out;
}



template <int dim,
          int n_components_,
          typename Number,
          bool is_face,
          typename VectorizedArrayType>
inline Tensor<1, n_components_, Tensor<1, dim, VectorizedArrayType>>
FEEvaluationBase<dim, n_components_, Number, is_face, VectorizedArrayType>::
  get_hessian_diagonal(const unsigned int q_point) const
{
  Assert(!is_face, ExcNotImplemented());
#  ifdef DEBUG
  Assert(this->hessians_quad_initialized == true,
         internal::ExcAccessToUninitializedField());
#  endif
  AssertIndexRange(q_point, this->n_quadrature_points);

  Assert(this->jacobian != nullptr, ExcNotImplemented());
  const Tensor<2, dim, VectorizedArrayType> &jac =
    this->jacobian[this->cell_type <= internal::MatrixFreeFunctions::affine ?
                     0 :
                     q_point];

  const std::size_t      nqp  = this->n_quadrature_points;
  constexpr unsigned int hdim = (dim * (dim + 1)) / 2;
  Tensor<1, n_components_, Tensor<1, dim, VectorizedArrayType>> hessian_out;

  // Cartesian cell
  if (this->cell_type == internal::MatrixFreeFunctions::cartesian)
    {
      for (unsigned int comp = 0; comp < n_components; comp++)
        for (unsigned int d = 0; d < dim; ++d)
          hessian_out[comp][d] =
            hessians_quad[(comp * hdim + d) * nqp + q_point] *
            (jac[d][d] * jac[d][d]);
    }
  // cell with general Jacobian, but constant within the cell
  else if (this->cell_type == internal::MatrixFreeFunctions::affine)
    {
      for (unsigned int comp = 0; comp < n_components; comp++)
        {
          // compute laplacian before the gradient because it needs to access
          // unscaled gradient data
          VectorizedArrayType tmp[dim][dim];
          internal::hessian_unit_times_jac(
            jac, hessians_quad + comp * hdim * nqp + q_point, nqp, tmp);

          // compute only the trace part of hessian, J * tmp = J *
          // hess_unit(u) * J^T
          for (unsigned int d = 0; d < dim; ++d)
            {
              hessian_out[comp][d] = jac[d][0] * tmp[0][d];
              for (unsigned int f = 1; f < dim; ++f)
                hessian_out[comp][d] += jac[d][f] * tmp[f][d];
            }
        }
    }
  // cell with general Jacobian
  else
    {
      const Tensor<1, dim *(dim + 1) / 2, Tensor<1, dim, VectorizedArrayType>>
        &jac_grad =
          this->mapping_data->jacobian_gradients
            [0][this->mapping_data->data_index_offsets[this->cell] + q_point];
      for (unsigned int comp = 0; comp < n_components; comp++)
        {
          // compute laplacian before the gradient because it needs to access
          // unscaled gradient data
          VectorizedArrayType tmp[dim][dim];
          internal::hessian_unit_times_jac(
            jac, hessians_quad + comp * hdim * nqp + q_point, nqp, tmp);

          // compute only the trace part of hessian, J * tmp = J *
          // hess_unit(u) * J^T
          for (unsigned int d = 0; d < dim; ++d)
            {
              hessian_out[comp][d] = jac[d][0] * tmp[0][d];
              for (unsigned int f = 1; f < dim; ++f)
                hessian_out[comp][d] += jac[d][f] * tmp[f][d];
            }

          for (unsigned int d = 0; d < dim; ++d)
            for (unsigned int e = 0; e < dim; ++e)
              hessian_out[comp][d] +=
                jac_grad[d][e] *
                gradients_quad[(comp * dim + e) * nqp + q_point];
        }
    }
  return hessian_out;
}



template <int dim,
          int n_components_,
          typename Number,
          bool is_face,
          typename VectorizedArrayType>
inline Tensor<1, n_components_, VectorizedArrayType>
FEEvaluationBase<dim, n_components_, Number, is_face, VectorizedArrayType>::
  get_laplacian(const unsigned int q_point) const
{
  Assert(is_face == false, ExcNotImplemented());
#  ifdef DEBUG
  Assert(this->hessians_quad_initialized == true,
         internal::ExcAccessToUninitializedField());
#  endif
  AssertIndexRange(q_point, this->n_quadrature_points);

  Tensor<1, n_components_, VectorizedArrayType> laplacian_out;
  const auto hess_diag = get_hessian_diagonal(q_point);
  for (unsigned int comp = 0; comp < n_components; ++comp)
    {
      laplacian_out[comp] = hess_diag[comp][0];
      for (unsigned int d = 1; d < dim; ++d)
        laplacian_out[comp] += hess_diag[comp][d];
    }
  return laplacian_out;
}



template <int dim,
          int n_components_,
          typename Number,
          bool is_face,
          typename VectorizedArrayType>
inline DEAL_II_ALWAYS_INLINE void
FEEvaluationBase<dim, n_components_, Number, is_face, VectorizedArrayType>::
  submit_dof_value(const Tensor<1, n_components_, VectorizedArrayType> val_in,
                   const unsigned int                                  dof)
{
#  ifdef DEBUG
  this->dof_values_initialized = true;
#  endif
  AssertIndexRange(dof, this->data->dofs_per_component_on_cell);
  for (unsigned int comp = 0; comp < n_components; comp++)
    this->values_dofs[comp][dof] = val_in[comp];
}



template <int dim,
          int n_components_,
          typename Number,
          bool is_face,
          typename VectorizedArrayType>
inline DEAL_II_ALWAYS_INLINE void
FEEvaluationBase<dim, n_components_, Number, is_face, VectorizedArrayType>::
  submit_value(const Tensor<1, n_components_, VectorizedArrayType> val_in,
               const unsigned int                                  q_point)
{
  Assert(this->cell != numbers::invalid_unsigned_int, ExcNotInitialized());
  AssertIndexRange(q_point, this->n_quadrature_points);
  Assert(this->J_value != nullptr,
         internal::ExcMatrixFreeAccessToUninitializedMappingField(
           "update_values"));
#  ifdef DEBUG
  this->values_quad_submitted = true;
#  endif

  const std::size_t nqp = this->n_quadrature_points;
  if (this->cell_type <= internal::MatrixFreeFunctions::affine)
    {
      const VectorizedArrayType JxW =
        this->J_value[0] * this->quadrature_weights[q_point];
      for (unsigned int comp = 0; comp < n_components; ++comp)
        values_quad[comp * nqp + q_point] = val_in[comp] * JxW;
    }
  else
    {
      const VectorizedArrayType JxW = this->J_value[q_point];
      for (unsigned int comp = 0; comp < n_components; ++comp)
        values_quad[comp * nqp + q_point] = val_in[comp] * JxW;
    }
}



template <int dim,
          int n_components_,
          typename Number,
          bool is_face,
          typename VectorizedArrayType>
inline DEAL_II_ALWAYS_INLINE void
FEEvaluationBase<dim, n_components_, Number, is_face, VectorizedArrayType>::
  submit_gradient(
    const Tensor<1, n_components_, Tensor<1, dim, VectorizedArrayType>> grad_in,
    const unsigned int                                                  q_point)
{
  Assert(this->cell != numbers::invalid_unsigned_int, ExcNotInitialized());
  AssertIndexRange(q_point, this->n_quadrature_points);
  Assert(this->J_value != nullptr,
         internal::ExcMatrixFreeAccessToUninitializedMappingField(
           "update_gradients"));
  Assert(this->jacobian != nullptr,
         internal::ExcMatrixFreeAccessToUninitializedMappingField(
           "update_gradients"));
#  ifdef DEBUG
  this->gradients_quad_submitted = true;
#  endif

  const std::size_t nqp = this->n_quadrature_points;
  if (!is_face && this->cell_type == internal::MatrixFreeFunctions::cartesian)
    {
      const VectorizedArrayType JxW =
        this->J_value[0] * this->quadrature_weights[q_point];
      for (unsigned int d = 0; d < dim; ++d)
        {
          const VectorizedArrayType factor = this->jacobian[0][d][d] * JxW;
          for (unsigned int comp = 0; comp < n_components; comp++)
            gradients_quad[(comp * dim + d) * nqp + q_point] =
              grad_in[comp][d] * factor;
        }
    }
  else
    {
      const Tensor<2, dim, VectorizedArrayType> jac =
        this->cell_type > internal::MatrixFreeFunctions::affine ?
          this->jacobian[q_point] :
          this->jacobian[0];
      const VectorizedArrayType JxW =
        this->cell_type > internal::MatrixFreeFunctions::affine ?
          this->J_value[q_point] :
          this->J_value[0] * this->quadrature_weights[q_point];
      for (unsigned int comp = 0; comp < n_components; ++comp)
        for (unsigned int d = 0; d < dim; ++d)
          {
            VectorizedArrayType new_val = jac[0][d] * grad_in[comp][0];
            for (unsigned int e = 1; e < dim; ++e)
              new_val += (jac[e][d] * grad_in[comp][e]);
            gradients_quad[(comp * dim + d) * nqp + q_point] = new_val * JxW;
          }
    }
}



template <int dim,
          int n_components_,
          typename Number,
          bool is_face,
          typename VectorizedArrayType>
inline DEAL_II_ALWAYS_INLINE void
FEEvaluationBase<dim, n_components_, Number, is_face, VectorizedArrayType>::
  submit_normal_derivative(
    const Tensor<1, n_components_, VectorizedArrayType> grad_in,
    const unsigned int                                  q_point)
{
  AssertIndexRange(q_point, this->n_quadrature_points);
  Assert(this->normal_x_jacobian != nullptr,
         internal::ExcMatrixFreeAccessToUninitializedMappingField(
           "update_gradients"));
#  ifdef DEBUG
  this->gradients_quad_submitted = true;
#  endif

  const std::size_t nqp = this->n_quadrature_points;
  if (this->cell_type == internal::MatrixFreeFunctions::cartesian)
    for (unsigned int comp = 0; comp < n_components; comp++)
      {
        for (unsigned int d = 0; d < dim - 1; ++d)
          gradients_quad[(comp * dim + d) * nqp + q_point] =
            VectorizedArrayType();
        gradients_quad[(comp * dim + dim - 1) * nqp + q_point] =
          grad_in[comp] *
          (this->normal_x_jacobian[0][dim - 1] * this->J_value[0] *
           this->quadrature_weights[q_point]);
      }
  else
    {
      const unsigned int index =
        this->cell_type <= internal::MatrixFreeFunctions::affine ? 0 : q_point;
      const Tensor<1, dim, VectorizedArrayType> jac =
        this->normal_x_jacobian[index];
      for (unsigned int comp = 0; comp < n_components; comp++)
        {
          VectorizedArrayType factor = grad_in[comp] * this->J_value[index];
          if (this->cell_type <= internal::MatrixFreeFunctions::affine)
            factor = factor * this->quadrature_weights[q_point];
          for (unsigned int d = 0; d < dim; ++d)
            gradients_quad[(comp * dim + d) * nqp + q_point] = factor * jac[d];
        }
    }
}



template <int dim,
          int n_components_,
          typename Number,
          bool is_face,
          typename VectorizedArrayType>
inline Tensor<1, n_components_, VectorizedArrayType>
FEEvaluationBase<dim, n_components_, Number, is_face, VectorizedArrayType>::
  integrate_value() const
{
  Assert(this->cell != numbers::invalid_unsigned_int, ExcNotInitialized());
#  ifdef DEBUG
  Assert(this->values_quad_submitted == true,
         internal::ExcAccessToUninitializedField());
#  endif

  Tensor<1, n_components_, VectorizedArrayType> return_value;
  const std::size_t                             nqp = this->n_quadrature_points;
  for (unsigned int q = 0; q < nqp; ++q)
    for (unsigned int comp = 0; comp < n_components; ++comp)
      return_value[comp] += this->values_quad[comp * nqp + q];
  return (return_value);
}



 /*----------------------- FEEvaluationAccess --------------------------------*/ 


template <int dim,
          int n_components_,
          typename Number,
          bool is_face,
          typename VectorizedArrayType>
inline FEEvaluationAccess<dim,
                          n_components_,
                          Number,
                          is_face,
                          VectorizedArrayType>::
  FEEvaluationAccess(
    const MatrixFree<dim, Number, VectorizedArrayType> &data_in,
    const unsigned int                                  dof_no,
    const unsigned int first_selected_component,
    const unsigned int quad_no_in,
    const unsigned int fe_degree,
    const unsigned int n_q_points,
    const bool         is_interior_face,
    const unsigned int active_fe_index,
    const unsigned int active_quad_index,
    const unsigned int face_type)
  : FEEvaluationBase<dim, n_components_, Number, is_face, VectorizedArrayType>(
      data_in,
      dof_no,
      first_selected_component,
      quad_no_in,
      fe_degree,
      n_q_points,
      is_interior_face,
      active_fe_index,
      active_quad_index,
      face_type)
{}



template <int dim,
          int n_components_,
          typename Number,
          bool is_face,
          typename VectorizedArrayType>
inline FEEvaluationAccess<dim,
                          n_components_,
                          Number,
                          is_face,
                          VectorizedArrayType>::
  FEEvaluationAccess(
    const Mapping<dim> &      mapping,
    const FiniteElement<dim> &fe,
    const Quadrature<1> &     quadrature,
    const UpdateFlags         update_flags,
    const unsigned int        first_selected_component,
    const FEEvaluationBaseData<dim, Number, is_face, VectorizedArrayType>
      *other)
  : FEEvaluationBase<dim, n_components_, Number, is_face, VectorizedArrayType>(
      mapping,
      fe,
      quadrature,
      update_flags,
      first_selected_component,
      other)
{}



template <int dim,
          int n_components_,
          typename Number,
          bool is_face,
          typename VectorizedArrayType>
inline FEEvaluationAccess<dim,
                          n_components_,
                          Number,
                          is_face,
                          VectorizedArrayType>::
  FEEvaluationAccess(const FEEvaluationAccess<dim,
                                              n_components_,
                                              Number,
                                              is_face,
                                              VectorizedArrayType> &other)
  : FEEvaluationBase<dim, n_components_, Number, is_face, VectorizedArrayType>(
      other)
{}



template <int dim,
          int n_components_,
          typename Number,
          bool is_face,
          typename VectorizedArrayType>
inline FEEvaluationAccess<dim,
                          n_components_,
                          Number,
                          is_face,
                          VectorizedArrayType> &
FEEvaluationAccess<dim, n_components_, Number, is_face, VectorizedArrayType>::
operator=(const FEEvaluationAccess<dim,
                                   n_components_,
                                   Number,
                                   is_face,
                                   VectorizedArrayType> &other)
{
  this->FEEvaluationBase<dim,
                         n_components_,
                         Number,
                         is_face,
                         VectorizedArrayType>::operator=(other);
  return *this;
}



 /*-------------------- FEEvaluationAccess scalar ----------------------------*/ 


template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline FEEvaluationAccess<dim, 1, Number, is_face, VectorizedArrayType>::
  FEEvaluationAccess(
    const MatrixFree<dim, Number, VectorizedArrayType> &data_in,
    const unsigned int                                  dof_no,
    const unsigned int first_selected_component,
    const unsigned int quad_no_in,
    const unsigned int fe_degree,
    const unsigned int n_q_points,
    const bool         is_interior_face,
    const unsigned int active_fe_index,
    const unsigned int active_quad_index,
    const unsigned int face_type)
  : FEEvaluationBase<dim, 1, Number, is_face, VectorizedArrayType>(
      data_in,
      dof_no,
      first_selected_component,
      quad_no_in,
      fe_degree,
      n_q_points,
      is_interior_face,
      active_fe_index,
      active_quad_index,
      face_type)
{}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline FEEvaluationAccess<dim, 1, Number, is_face, VectorizedArrayType>::
  FEEvaluationAccess(
    const Mapping<dim> &      mapping,
    const FiniteElement<dim> &fe,
    const Quadrature<1> &     quadrature,
    const UpdateFlags         update_flags,
    const unsigned int        first_selected_component,
    const FEEvaluationBaseData<dim, Number, is_face, VectorizedArrayType>
      *other)
  : FEEvaluationBase<dim, 1, Number, is_face, VectorizedArrayType>(
      mapping,
      fe,
      quadrature,
      update_flags,
      first_selected_component,
      other)
{}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline FEEvaluationAccess<dim, 1, Number, is_face, VectorizedArrayType>::
  FEEvaluationAccess(
    const FEEvaluationAccess<dim, 1, Number, is_face, VectorizedArrayType>
      &other)
  : FEEvaluationBase<dim, 1, Number, is_face, VectorizedArrayType>(other)
{}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline FEEvaluationAccess<dim, 1, Number, is_face, VectorizedArrayType> &
FEEvaluationAccess<dim, 1, Number, is_face, VectorizedArrayType>::operator=(
  const FEEvaluationAccess<dim, 1, Number, is_face, VectorizedArrayType> &other)
{
  this->FEEvaluationBase<dim, 1, Number, is_face, VectorizedArrayType>::
  operator=(other);
  return *this;
}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline DEAL_II_ALWAYS_INLINE VectorizedArrayType
                             FEEvaluationAccess<dim, 1, Number, is_face, VectorizedArrayType>::get_dof_value(
  const unsigned int dof) const
{
  AssertIndexRange(dof, this->data->dofs_per_component_on_cell);
  return this->values_dofs[0][dof];
}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline DEAL_II_ALWAYS_INLINE VectorizedArrayType
                             FEEvaluationAccess<dim, 1, Number, is_face, VectorizedArrayType>::get_value(
  const unsigned int q_point) const
{
#  ifdef DEBUG
  Assert(this->values_quad_initialized == true,
         internal::ExcAccessToUninitializedField());
#  endif
  AssertIndexRange(q_point, this->n_quadrature_points);
  return this->values_quad[q_point];
}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline DEAL_II_ALWAYS_INLINE VectorizedArrayType
                             FEEvaluationAccess<dim, 1, Number, is_face, VectorizedArrayType>::
  get_normal_derivative(const unsigned int q_point) const
{
  return BaseClass::get_normal_derivative(q_point)[0];
}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline DEAL_II_ALWAYS_INLINE Tensor<1, dim, VectorizedArrayType>
                             FEEvaluationAccess<dim, 1, Number, is_face, VectorizedArrayType>::get_gradient(
  const unsigned int q_point) const
{
  // could use the base class gradient, but that involves too many expensive
  // initialization operations on tensors

#  ifdef DEBUG
  Assert(this->gradients_quad_initialized == true,
         internal::ExcAccessToUninitializedField());
#  endif
  AssertIndexRange(q_point, this->n_quadrature_points);

  Assert(this->jacobian != nullptr,
         internal::ExcMatrixFreeAccessToUninitializedMappingField(
           "update_gradients"));

  Tensor<1, dim, VectorizedArrayType> grad_out;

  const std::size_t nqp = this->n_quadrature_points;
  if (!is_face && this->cell_type == internal::MatrixFreeFunctions::cartesian)
    {
      for (unsigned int d = 0; d < dim; ++d)
        grad_out[d] =
          this->gradients_quad[d * nqp + q_point] * this->jacobian[0][d][d];
    }
  // cell with general/affine Jacobian
  else
    {
      const Tensor<2, dim, VectorizedArrayType> &jac =
        this->jacobian[this->cell_type > internal::MatrixFreeFunctions::affine ?
                         q_point :
                         0];
      for (unsigned int d = 0; d < dim; ++d)
        {
          grad_out[d] = jac[d][0] * this->gradients_quad[q_point];
          for (unsigned int e = 1; e < dim; ++e)
            grad_out[d] += jac[d][e] * this->gradients_quad[e * nqp + q_point];
        }
    }
  return grad_out;
}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline Tensor<2, dim, VectorizedArrayType>
FEEvaluationAccess<dim, 1, Number, is_face, VectorizedArrayType>::get_hessian(
  const unsigned int q_point) const
{
  return BaseClass::get_hessian(q_point)[0];
}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline Tensor<1, dim, VectorizedArrayType>
FEEvaluationAccess<dim, 1, Number, is_face, VectorizedArrayType>::
  get_hessian_diagonal(const unsigned int q_point) const
{
  return BaseClass::get_hessian_diagonal(q_point)[0];
}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline VectorizedArrayType
FEEvaluationAccess<dim, 1, Number, is_face, VectorizedArrayType>::get_laplacian(
  const unsigned int q_point) const
{
  return BaseClass::get_laplacian(q_point)[0];
}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline void DEAL_II_ALWAYS_INLINE
            FEEvaluationAccess<dim, 1, Number, is_face, VectorizedArrayType>::
  submit_dof_value(const VectorizedArrayType val_in, const unsigned int dof)
{
#  ifdef DEBUG
  this->dof_values_initialized = true;
  AssertIndexRange(dof, this->data->dofs_per_component_on_cell);
#  endif
  this->values_dofs[0][dof] = val_in;
}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline void DEAL_II_ALWAYS_INLINE
            FEEvaluationAccess<dim, 1, Number, is_face, VectorizedArrayType>::submit_value(
  const VectorizedArrayType val_in,
  const unsigned int        q_point)
{
  Assert(this->cell != numbers::invalid_unsigned_int, ExcNotInitialized());
  AssertIndexRange(q_point, this->n_quadrature_points);
  Assert(this->J_value != nullptr,
         internal::ExcMatrixFreeAccessToUninitializedMappingField(
           "update_value"));
#  ifdef DEBUG
  this->values_quad_submitted = true;
#  endif

  if (this->cell_type <= internal::MatrixFreeFunctions::affine)
    {
      const VectorizedArrayType JxW =
        this->J_value[0] * this->quadrature_weights[q_point];
      this->values_quad[q_point] = val_in * JxW;
    }
  else // if (this->cell_type < internal::MatrixFreeFunctions::general)
    {
      this->values_quad[q_point] = val_in * this->J_value[q_point];
    }
}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline DEAL_II_ALWAYS_INLINE void
FEEvaluationAccess<dim, 1, Number, is_face, VectorizedArrayType>::submit_value(
  const Tensor<1, 1, VectorizedArrayType> val_in,
  const unsigned int                      q_point)
{
  submit_value(val_in[0], q_point);
}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline DEAL_II_ALWAYS_INLINE void
FEEvaluationAccess<dim, 1, Number, is_face, VectorizedArrayType>::
  submit_normal_derivative(const VectorizedArrayType grad_in,
                           const unsigned int        q_point)
{
  Tensor<1, 1, VectorizedArrayType> grad;
  grad[0] = grad_in;
  BaseClass::submit_normal_derivative(grad, q_point);
}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline DEAL_II_ALWAYS_INLINE void
FEEvaluationAccess<dim, 1, Number, is_face, VectorizedArrayType>::
  submit_gradient(const Tensor<1, dim, VectorizedArrayType> grad_in,
                  const unsigned int                        q_point)
{
  Assert(this->cell != numbers::invalid_unsigned_int, ExcNotInitialized());
  AssertIndexRange(q_point, this->n_quadrature_points);
  Assert(this->J_value != nullptr,
         internal::ExcMatrixFreeAccessToUninitializedMappingField(
           "update_gradients"));
  Assert(this->jacobian != nullptr,
         internal::ExcMatrixFreeAccessToUninitializedMappingField(
           "update_gradients"));
#  ifdef DEBUG
  this->gradients_quad_submitted = true;
#  endif

  const std::size_t nqp = this->n_quadrature_points;
  if (!is_face && this->cell_type == internal::MatrixFreeFunctions::cartesian)
    {
      const VectorizedArrayType JxW =
        this->J_value[0] * this->quadrature_weights[q_point];
      for (unsigned int d = 0; d < dim; ++d)
        this->gradients_quad[d * nqp + q_point] =
          (grad_in[d] * this->jacobian[0][d][d] * JxW);
    }
  // general/affine cell type
  else
    {
      const Tensor<2, dim, VectorizedArrayType> &jac =
        this->cell_type > internal::MatrixFreeFunctions::affine ?
          this->jacobian[q_point] :
          this->jacobian[0];
      const VectorizedArrayType JxW =
        this->cell_type > internal::MatrixFreeFunctions::affine ?
          this->J_value[q_point] :
          this->J_value[0] * this->quadrature_weights[q_point];
      for (unsigned int d = 0; d < dim; ++d)
        {
          VectorizedArrayType new_val = jac[0][d] * grad_in[0];
          for (unsigned int e = 1; e < dim; ++e)
            new_val += jac[e][d] * grad_in[e];
          this->gradients_quad[d * nqp + q_point] = new_val * JxW;
        }
    }
}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline VectorizedArrayType
FEEvaluationAccess<dim, 1, Number, is_face, VectorizedArrayType>::
  integrate_value() const
{
  return BaseClass::integrate_value()[0];
}



 /*----------------- FEEvaluationAccess vector-valued ------------------------*/ 


template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline FEEvaluationAccess<dim, dim, Number, is_face, VectorizedArrayType>::
  FEEvaluationAccess(
    const MatrixFree<dim, Number, VectorizedArrayType> &data_in,
    const unsigned int                                  dof_no,
    const unsigned int first_selected_component,
    const unsigned int quad_no_in,
    const unsigned int fe_degree,
    const unsigned int n_q_points,
    const bool         is_interior_face,
    const unsigned int active_fe_index,
    const unsigned int active_quad_index,
    const unsigned int face_type)
  : FEEvaluationBase<dim, dim, Number, is_face, VectorizedArrayType>(
      data_in,
      dof_no,
      first_selected_component,
      quad_no_in,
      fe_degree,
      n_q_points,
      is_interior_face,
      active_fe_index,
      active_quad_index,
      face_type)
{}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline FEEvaluationAccess<dim, dim, Number, is_face, VectorizedArrayType>::
  FEEvaluationAccess(
    const Mapping<dim> &      mapping,
    const FiniteElement<dim> &fe,
    const Quadrature<1> &     quadrature,
    const UpdateFlags         update_flags,
    const unsigned int        first_selected_component,
    const FEEvaluationBaseData<dim, Number, is_face, VectorizedArrayType>
      *other)
  : FEEvaluationBase<dim, dim, Number, is_face, VectorizedArrayType>(
      mapping,
      fe,
      quadrature,
      update_flags,
      first_selected_component,
      other)
{}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline FEEvaluationAccess<dim, dim, Number, is_face, VectorizedArrayType>::
  FEEvaluationAccess(
    const FEEvaluationAccess<dim, dim, Number, is_face, VectorizedArrayType>
      &other)
  : FEEvaluationBase<dim, dim, Number, is_face, VectorizedArrayType>(other)
{}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline FEEvaluationAccess<dim, dim, Number, is_face, VectorizedArrayType> &
FEEvaluationAccess<dim, dim, Number, is_face, VectorizedArrayType>::operator=(
  const FEEvaluationAccess<dim, dim, Number, is_face, VectorizedArrayType>
    &other)
{
  this->FEEvaluationBase<dim, dim, Number, is_face, VectorizedArrayType>::
  operator=(other);
  return *this;
}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline DEAL_II_ALWAYS_INLINE Tensor<2, dim, VectorizedArrayType>
                             FEEvaluationAccess<dim, dim, Number, is_face, VectorizedArrayType>::
  get_gradient(const unsigned int q_point) const
{
  return BaseClass::get_gradient(q_point);
}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline DEAL_II_ALWAYS_INLINE VectorizedArrayType
                             FEEvaluationAccess<dim, dim, Number, is_face, VectorizedArrayType>::
  get_divergence(const unsigned int q_point) const
{
#  ifdef DEBUG
  Assert(this->gradients_quad_initialized == true,
         internal::ExcAccessToUninitializedField());
#  endif
  AssertIndexRange(q_point, this->n_quadrature_points);
  Assert(this->jacobian != nullptr,
         internal::ExcMatrixFreeAccessToUninitializedMappingField(
           "update_gradients"));

  VectorizedArrayType divergence;
  const std::size_t   nqp = this->n_quadrature_points;

  // Cartesian cell
  if (!is_face && this->cell_type == internal::MatrixFreeFunctions::cartesian)
    {
      divergence = this->gradients_quad[q_point] * this->jacobian[0][0][0];
      for (unsigned int d = 1; d < dim; ++d)
        divergence += this->gradients_quad[(dim * d + d) * nqp + q_point] *
                      this->jacobian[0][d][d];
    }
  // cell with general/constant Jacobian
  else
    {
      const Tensor<2, dim, VectorizedArrayType> &jac =
        this->cell_type == internal::MatrixFreeFunctions::general ?
          this->jacobian[q_point] :
          this->jacobian[0];
      divergence = jac[0][0] * this->gradients_quad[q_point];
      for (unsigned int e = 1; e < dim; ++e)
        divergence += jac[0][e] * this->gradients_quad[e * nqp + q_point];
      for (unsigned int d = 1; d < dim; ++d)
        for (unsigned int e = 0; e < dim; ++e)
          divergence +=
            jac[d][e] * this->gradients_quad[(d * dim + e) * nqp + q_point];
    }
  return divergence;
}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline DEAL_II_ALWAYS_INLINE SymmetricTensor<2, dim, VectorizedArrayType>
                             FEEvaluationAccess<dim, dim, Number, is_face, VectorizedArrayType>::
  get_symmetric_gradient(const unsigned int q_point) const
{
  // copy from generic function into dim-specialization function
  const auto          grad = get_gradient(q_point);
  VectorizedArrayType symmetrized[(dim * dim + dim) / 2];
  VectorizedArrayType half = Number(0.5);
  for (unsigned int d = 0; d < dim; ++d)
    symmetrized[d] = grad[d][d];
  switch (dim)
    {
      case 1:
        break;
      case 2:
        symmetrized[2] = grad[0][1] + grad[1][0];
        symmetrized[2] *= half;
        break;
      case 3:
        symmetrized[3] = grad[0][1] + grad[1][0];
        symmetrized[3] *= half;
        symmetrized[4] = grad[0][2] + grad[2][0];
        symmetrized[4] *= half;
        symmetrized[5] = grad[1][2] + grad[2][1];
        symmetrized[5] *= half;
        break;
      default:
        Assert(false, ExcNotImplemented());
    }
  return SymmetricTensor<2, dim, VectorizedArrayType>(symmetrized);
}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline DEAL_II_ALWAYS_INLINE
  Tensor<1, (dim == 2 ? 1 : dim), VectorizedArrayType>
  FEEvaluationAccess<dim, dim, Number, is_face, VectorizedArrayType>::get_curl(
    const unsigned int q_point) const
{
  // copy from generic function into dim-specialization function
  const Tensor<2, dim, VectorizedArrayType> grad = get_gradient(q_point);
  Tensor<1, (dim == 2 ? 1 : dim), VectorizedArrayType> curl;
  switch (dim)
    {
      case 1:
        Assert(false,
               ExcMessage(
                 "Computing the curl in 1d is not a useful operation"));
        break;
      case 2:
        curl[0] = grad[1][0] - grad[0][1];
        break;
      case 3:
        curl[0] = grad[2][1] - grad[1][2];
        curl[1] = grad[0][2] - grad[2][0];
        curl[2] = grad[1][0] - grad[0][1];
        break;
      default:
        Assert(false, ExcNotImplemented());
    }
  return curl;
}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline DEAL_II_ALWAYS_INLINE Tensor<2, dim, VectorizedArrayType>
                             FEEvaluationAccess<dim, dim, Number, is_face, VectorizedArrayType>::
  get_hessian_diagonal(const unsigned int q_point) const
{
  return BaseClass::get_hessian_diagonal(q_point);
}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline DEAL_II_ALWAYS_INLINE Tensor<3, dim, VectorizedArrayType>
                             FEEvaluationAccess<dim, dim, Number, is_face, VectorizedArrayType>::get_hessian(
  const unsigned int q_point) const
{
#  ifdef DEBUG
  Assert(this->hessians_quad_initialized == true,
         internal::ExcAccessToUninitializedField());
#  endif
  AssertIndexRange(q_point, this->n_quadrature_points);
  return BaseClass::get_hessian(q_point);
}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline DEAL_II_ALWAYS_INLINE void
FEEvaluationAccess<dim, dim, Number, is_face, VectorizedArrayType>::
  submit_gradient(const Tensor<2, dim, VectorizedArrayType> grad_in,
                  const unsigned int                        q_point)
{
  BaseClass::submit_gradient(grad_in, q_point);
}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline DEAL_II_ALWAYS_INLINE void
FEEvaluationAccess<dim, dim, Number, is_face, VectorizedArrayType>::
  submit_gradient(
    const Tensor<1, dim, Tensor<1, dim, VectorizedArrayType>> grad_in,
    const unsigned int                                        q_point)
{
  BaseClass::submit_gradient(grad_in, q_point);
}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline DEAL_II_ALWAYS_INLINE void
FEEvaluationAccess<dim, dim, Number, is_face, VectorizedArrayType>::
  submit_divergence(const VectorizedArrayType div_in,
                    const unsigned int        q_point)
{
  Assert(this->cell != numbers::invalid_unsigned_int, ExcNotInitialized());
  AssertIndexRange(q_point, this->n_quadrature_points);
  Assert(this->J_value != nullptr,
         internal::ExcMatrixFreeAccessToUninitializedMappingField(
           "update_gradients"));
  Assert(this->jacobian != nullptr,
         internal::ExcMatrixFreeAccessToUninitializedMappingField(
           "update_gradients"));
#  ifdef DEBUG
  this->gradients_quad_submitted = true;
#  endif

  const std::size_t nqp = this->n_quadrature_points;
  if (!is_face && this->cell_type == internal::MatrixFreeFunctions::cartesian)
    {
      const VectorizedArrayType fac =
        this->J_value[0] * this->quadrature_weights[q_point] * div_in;
      for (unsigned int d = 0; d < dim; ++d)
        {
          this->gradients_quad[(d * dim + d) * nqp + q_point] =
            (fac * this->jacobian[0][d][d]);
          for (unsigned int e = d + 1; e < dim; ++e)
            {
              this->gradients_quad[(d * dim + e) * nqp + q_point] =
                VectorizedArrayType();
              this->gradients_quad[(e * dim + d) * nqp + q_point] =
                VectorizedArrayType();
            }
        }
    }
  else
    {
      const Tensor<2, dim, VectorizedArrayType> &jac =
        this->cell_type == internal::MatrixFreeFunctions::general ?
          this->jacobian[q_point] :
          this->jacobian[0];
      const VectorizedArrayType fac =
        (this->cell_type == internal::MatrixFreeFunctions::general ?
           this->J_value[q_point] :
           this->J_value[0] * this->quadrature_weights[q_point]) *
        div_in;
      for (unsigned int d = 0; d < dim; ++d)
        {
          for (unsigned int e = 0; e < dim; ++e)
            this->gradients_quad[(d * dim + e) * nqp + q_point] =
              jac[d][e] * fac;
        }
    }
}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline DEAL_II_ALWAYS_INLINE void
FEEvaluationAccess<dim, dim, Number, is_face, VectorizedArrayType>::
  submit_symmetric_gradient(
    const SymmetricTensor<2, dim, VectorizedArrayType> sym_grad,
    const unsigned int                                 q_point)
{
  // could have used base class operator, but that involves some overhead
  // which is inefficient. it is nice to have the symmetric tensor because
  // that saves some operations
  Assert(this->cell != numbers::invalid_unsigned_int, ExcNotInitialized());
  AssertIndexRange(q_point, this->n_quadrature_points);
  Assert(this->J_value != nullptr,
         internal::ExcMatrixFreeAccessToUninitializedMappingField(
           "update_gradients"));
  Assert(this->jacobian != nullptr,
         internal::ExcMatrixFreeAccessToUninitializedMappingField(
           "update_gradients"));
#  ifdef DEBUG
  this->gradients_quad_submitted = true;
#  endif

  const std::size_t nqp = this->n_quadrature_points;
  if (!is_face && this->cell_type == internal::MatrixFreeFunctions::cartesian)
    {
      const VectorizedArrayType JxW =
        this->J_value[0] * this->quadrature_weights[q_point];
      for (unsigned int d = 0; d < dim; ++d)
        this->gradients_quad[(d * dim + d) * nqp + q_point] =
          (sym_grad.access_raw_entry(d) * JxW * this->jacobian[0][d][d]);
      for (unsigned int e = 0, counter = dim; e < dim; ++e)
        for (unsigned int d = e + 1; d < dim; ++d, ++counter)
          {
            const VectorizedArrayType value =
              sym_grad.access_raw_entry(counter) * JxW;
            this->gradients_quad[(e * dim + d) * nqp + q_point] =
              value * this->jacobian[0][d][d];
            this->gradients_quad[(d * dim + e) * nqp + q_point] =
              value * this->jacobian[0][e][e];
          }
    }
  // general/affine cell type
  else
    {
      const VectorizedArrayType JxW =
        this->cell_type == internal::MatrixFreeFunctions::general ?
          this->J_value[q_point] :
          this->J_value[0] * this->quadrature_weights[q_point];
      const Tensor<2, dim, VectorizedArrayType> &jac =
        this->cell_type == internal::MatrixFreeFunctions::general ?
          this->jacobian[q_point] :
          this->jacobian[0];
      VectorizedArrayType weighted[dim][dim];
      for (unsigned int i = 0; i < dim; ++i)
        weighted[i][i] = sym_grad.access_raw_entry(i) * JxW;
      for (unsigned int i = 0, counter = dim; i < dim; ++i)
        for (unsigned int j = i + 1; j < dim; ++j, ++counter)
          {
            const VectorizedArrayType value =
              sym_grad.access_raw_entry(counter) * JxW;
            weighted[i][j] = value;
            weighted[j][i] = value;
          }
      for (unsigned int comp = 0; comp < dim; ++comp)
        for (unsigned int d = 0; d < dim; ++d)
          {
            VectorizedArrayType new_val = jac[0][d] * weighted[comp][0];
            for (unsigned int e = 1; e < dim; ++e)
              new_val += jac[e][d] * weighted[comp][e];
            this->gradients_quad[(comp * dim + d) * nqp + q_point] = new_val;
          }
    }
}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline DEAL_II_ALWAYS_INLINE void
FEEvaluationAccess<dim, dim, Number, is_face, VectorizedArrayType>::submit_curl(
  const Tensor<1, dim == 2 ? 1 : dim, VectorizedArrayType> curl,
  const unsigned int                                       q_point)
{
  Tensor<2, dim, VectorizedArrayType> grad;
  switch (dim)
    {
      case 1:
        Assert(false,
               ExcMessage(
                 "Testing by the curl in 1d is not a useful operation"));
        break;
      case 2:
        grad[1][0] = curl[0];
        grad[0][1] = -curl[0];
        break;
      case 3:
        grad[2][1] = curl[0];
        grad[1][2] = -curl[0];
        grad[0][2] = curl[1];
        grad[2][0] = -curl[1];
        grad[1][0] = curl[2];
        grad[0][1] = -curl[2];
        break;
      default:
        Assert(false, ExcNotImplemented());
    }
  submit_gradient(grad, q_point);
}


 /*-------------------- FEEvaluationAccess scalar for 1d ---------------------*/ 


template <typename Number, bool is_face, typename VectorizedArrayType>
inline FEEvaluationAccess<1, 1, Number, is_face, VectorizedArrayType>::
  FEEvaluationAccess(const MatrixFree<1, Number, VectorizedArrayType> &data_in,
                     const unsigned int                                dof_no,
                     const unsigned int first_selected_component,
                     const unsigned int quad_no_in,
                     const unsigned int fe_degree,
                     const unsigned int n_q_points,
                     const bool         is_interior_face,
                     const unsigned int active_fe_index,
                     const unsigned int active_quad_index,
                     const unsigned int face_type)
  : FEEvaluationBase<1, 1, Number, is_face, VectorizedArrayType>(
      data_in,
      dof_no,
      first_selected_component,
      quad_no_in,
      fe_degree,
      n_q_points,
      is_interior_face,
      active_fe_index,
      active_quad_index,
      face_type)
{}



template <typename Number, bool is_face, typename VectorizedArrayType>
inline FEEvaluationAccess<1, 1, Number, is_face, VectorizedArrayType>::
  FEEvaluationAccess(
    const Mapping<1> &      mapping,
    const FiniteElement<1> &fe,
    const Quadrature<1> &   quadrature,
    const UpdateFlags       update_flags,
    const unsigned int      first_selected_component,
    const FEEvaluationBaseData<1, Number, is_face, VectorizedArrayType> *other)
  : FEEvaluationBase<1, 1, Number, is_face, VectorizedArrayType>(
      mapping,
      fe,
      quadrature,
      update_flags,
      first_selected_component,
      other)
{}



template <typename Number, bool is_face, typename VectorizedArrayType>
inline FEEvaluationAccess<1, 1, Number, is_face, VectorizedArrayType>::
  FEEvaluationAccess(
    const FEEvaluationAccess<1, 1, Number, is_face, VectorizedArrayType> &other)
  : FEEvaluationBase<1, 1, Number, is_face, VectorizedArrayType>(other)
{}



template <typename Number, bool is_face, typename VectorizedArrayType>
inline FEEvaluationAccess<1, 1, Number, is_face, VectorizedArrayType> &
FEEvaluationAccess<1, 1, Number, is_face, VectorizedArrayType>::operator=(
  const FEEvaluationAccess<1, 1, Number, is_face, VectorizedArrayType> &other)
{
  this->FEEvaluationBase<1, 1, Number, is_face, VectorizedArrayType>::operator=(
    other);
  return *this;
}



template <typename Number, bool is_face, typename VectorizedArrayType>
inline DEAL_II_ALWAYS_INLINE VectorizedArrayType
                             FEEvaluationAccess<1, 1, Number, is_face, VectorizedArrayType>::get_dof_value(
  const unsigned int dof) const
{
  AssertIndexRange(dof, this->data->dofs_per_component_on_cell);
  return this->values_dofs[0][dof];
}



template <typename Number, bool is_face, typename VectorizedArrayType>
inline DEAL_II_ALWAYS_INLINE VectorizedArrayType
                             FEEvaluationAccess<1, 1, Number, is_face, VectorizedArrayType>::get_value(
  const unsigned int q_point) const
{
#  ifdef DEBUG
  Assert(this->values_quad_initialized == true,
         internal::ExcAccessToUninitializedField());
#  endif
  AssertIndexRange(q_point, this->n_quadrature_points);
  return this->values_quad[q_point];
}



template <typename Number, bool is_face, typename VectorizedArrayType>
inline DEAL_II_ALWAYS_INLINE Tensor<1, 1, VectorizedArrayType>
                             FEEvaluationAccess<1, 1, Number, is_face, VectorizedArrayType>::get_gradient(
  const unsigned int q_point) const
{
  // could use the base class gradient, but that involves too many inefficient
  // initialization operations on tensors

#  ifdef DEBUG
  Assert(this->gradients_quad_initialized == true,
         internal::ExcAccessToUninitializedField());
#  endif
  AssertIndexRange(q_point, this->n_quadrature_points);

  const Tensor<2, 1, VectorizedArrayType> &jac =
    this->cell_type == internal::MatrixFreeFunctions::general ?
      this->jacobian[q_point] :
      this->jacobian[0];

  Tensor<1, 1, VectorizedArrayType> grad_out;
  grad_out[0] = jac[0][0] * this->gradients_quad[q_point];

  return grad_out;
}



template <typename Number, bool is_face, typename VectorizedArrayType>
inline DEAL_II_ALWAYS_INLINE VectorizedArrayType
                             FEEvaluationAccess<1, 1, Number, is_face, VectorizedArrayType>::get_divergence(
  const unsigned int q_point) const
{
  return get_gradient(q_point)[0];
}



template <typename Number, bool is_face, typename VectorizedArrayType>
inline DEAL_II_ALWAYS_INLINE VectorizedArrayType
                             FEEvaluationAccess<1, 1, Number, is_face, VectorizedArrayType>::
  get_normal_derivative(const unsigned int q_point) const
{
  return BaseClass::get_normal_derivative(q_point)[0];
}



template <typename Number, bool is_face, typename VectorizedArrayType>
inline DEAL_II_ALWAYS_INLINE Tensor<2, 1, VectorizedArrayType>
                             FEEvaluationAccess<1, 1, Number, is_face, VectorizedArrayType>::get_hessian(
  const unsigned int q_point) const
{
  return BaseClass::get_hessian(q_point)[0];
}



template <typename Number, bool is_face, typename VectorizedArrayType>
inline DEAL_II_ALWAYS_INLINE Tensor<1, 1, VectorizedArrayType>
                             FEEvaluationAccess<1, 1, Number, is_face, VectorizedArrayType>::
  get_hessian_diagonal(const unsigned int q_point) const
{
  return BaseClass::get_hessian_diagonal(q_point)[0];
}



template <typename Number, bool is_face, typename VectorizedArrayType>
inline DEAL_II_ALWAYS_INLINE VectorizedArrayType
                             FEEvaluationAccess<1, 1, Number, is_face, VectorizedArrayType>::get_laplacian(
  const unsigned int q_point) const
{
  return BaseClass::get_laplacian(q_point)[0];
}



template <typename Number, bool is_face, typename VectorizedArrayType>
inline DEAL_II_ALWAYS_INLINE void DEAL_II_ALWAYS_INLINE
                                  FEEvaluationAccess<1, 1, Number, is_face, VectorizedArrayType>::
  submit_dof_value(const VectorizedArrayType val_in, const unsigned int dof)
{
#  ifdef DEBUG
  this->dof_values_initialized = true;
  AssertIndexRange(dof, this->data->dofs_per_component_on_cell);
#  endif
  this->values_dofs[0][dof] = val_in;
}



template <typename Number, bool is_face, typename VectorizedArrayType>
inline DEAL_II_ALWAYS_INLINE void
FEEvaluationAccess<1, 1, Number, is_face, VectorizedArrayType>::submit_value(
  const VectorizedArrayType val_in,
  const unsigned int        q_point)
{
  Assert(this->cell != numbers::invalid_unsigned_int, ExcNotInitialized());
  AssertIndexRange(q_point, this->n_quadrature_points);
#  ifdef DEBUG
  this->values_quad_submitted = true;
#  endif

  if (this->cell_type == internal::MatrixFreeFunctions::general)
    {
      const VectorizedArrayType JxW = this->J_value[q_point];
      this->values_quad[q_point]    = val_in * JxW;
    }
  else // if (this->cell_type == internal::MatrixFreeFunctions::general)
    {
      const VectorizedArrayType JxW =
        this->J_value[0] * this->quadrature_weights[q_point];
      this->values_quad[q_point] = val_in * JxW;
    }
}



template <typename Number, bool is_face, typename VectorizedArrayType>
inline DEAL_II_ALWAYS_INLINE void
FEEvaluationAccess<1, 1, Number, is_face, VectorizedArrayType>::submit_value(
  const Tensor<1, 1, VectorizedArrayType> val_in,
  const unsigned int                      q_point)
{
  submit_value(val_in[0], q_point);
}



template <typename Number, bool is_face, typename VectorizedArrayType>
inline DEAL_II_ALWAYS_INLINE void
FEEvaluationAccess<1, 1, Number, is_face, VectorizedArrayType>::submit_gradient(
  const Tensor<1, 1, VectorizedArrayType> grad_in,
  const unsigned int                      q_point)
{
  submit_gradient(grad_in[0], q_point);
}



template <typename Number, bool is_face, typename VectorizedArrayType>
inline DEAL_II_ALWAYS_INLINE void
FEEvaluationAccess<1, 1, Number, is_face, VectorizedArrayType>::submit_gradient(
  const VectorizedArrayType grad_in,
  const unsigned int        q_point)
{
  Assert(this->cell != numbers::invalid_unsigned_int, ExcNotInitialized());
  AssertIndexRange(q_point, this->n_quadrature_points);
#  ifdef DEBUG
  this->gradients_quad_submitted = true;
#  endif

  const Tensor<2, 1, VectorizedArrayType> &jac =
    this->cell_type == internal::MatrixFreeFunctions::general ?
      this->jacobian[q_point] :
      this->jacobian[0];
  const VectorizedArrayType JxW =
    this->cell_type == internal::MatrixFreeFunctions::general ?
      this->J_value[q_point] :
      this->J_value[0] * this->quadrature_weights[q_point];

  this->gradients_quad[q_point] = jac[0][0] * grad_in * JxW;
}



template <typename Number, bool is_face, typename VectorizedArrayType>
inline DEAL_II_ALWAYS_INLINE void
FEEvaluationAccess<1, 1, Number, is_face, VectorizedArrayType>::submit_gradient(
  const Tensor<2, 1, VectorizedArrayType> grad_in,
  const unsigned int                      q_point)
{
  submit_gradient(grad_in[0][0], q_point);
}



template <typename Number, bool is_face, typename VectorizedArrayType>
inline DEAL_II_ALWAYS_INLINE void
FEEvaluationAccess<1, 1, Number, is_face, VectorizedArrayType>::
  submit_normal_derivative(const VectorizedArrayType grad_in,
                           const unsigned int        q_point)
{
  Tensor<1, 1, VectorizedArrayType> grad;
  grad[0] = grad_in;
  BaseClass::submit_normal_derivative(grad, q_point);
}



template <typename Number, bool is_face, typename VectorizedArrayType>
inline DEAL_II_ALWAYS_INLINE void
FEEvaluationAccess<1, 1, Number, is_face, VectorizedArrayType>::
  submit_normal_derivative(const Tensor<1, 1, VectorizedArrayType> grad_in,
                           const unsigned int                      q_point)
{
  BaseClass::submit_normal_derivative(grad_in, q_point);
}



template <typename Number, bool is_face, typename VectorizedArrayType>
inline VectorizedArrayType
FEEvaluationAccess<1, 1, Number, is_face, VectorizedArrayType>::
  integrate_value() const
{
  return BaseClass::integrate_value()[0];
}



 /*-------------------------- FEEvaluation -----------------------------------*/ 


template <int dim,
          int fe_degree,
          int n_q_points_1d,
          int n_components_,
          typename Number,
          typename VectorizedArrayType>
inline FEEvaluation<dim,
                    fe_degree,
                    n_q_points_1d,
                    n_components_,
                    Number,
                    VectorizedArrayType>::
  FEEvaluation(const MatrixFree<dim, Number, VectorizedArrayType> &data_in,
               const unsigned int                                  fe_no,
               const unsigned int                                  quad_no,
               const unsigned int first_selected_component,
               const unsigned int active_fe_index,
               const unsigned int active_quad_index)
  : BaseClass(data_in,
              fe_no,
              first_selected_component,
              quad_no,
              fe_degree,
              static_n_q_points,
              true  /*note: this is not a face*/ ,
              active_fe_index,
              active_quad_index)
  , dofs_per_component(this->data->dofs_per_component_on_cell)
  , dofs_per_cell(this->data->dofs_per_component_on_cell * n_components_)
  , n_q_points(this->data->n_q_points)
{
  check_template_arguments(fe_no, 0);
}



template <int dim,
          int fe_degree,
          int n_q_points_1d,
          int n_components_,
          typename Number,
          typename VectorizedArrayType>
inline FEEvaluation<dim,
                    fe_degree,
                    n_q_points_1d,
                    n_components_,
                    Number,
                    VectorizedArrayType>::
  FEEvaluation(const MatrixFree<dim, Number, VectorizedArrayType> &matrix_free,
               const std::pair<unsigned int, unsigned int> &       range,
               const unsigned int                                  dof_no,
               const unsigned int                                  quad_no,
               const unsigned int first_selected_component)
  : FEEvaluation(matrix_free,
                 dof_no,
                 quad_no,
                 first_selected_component,
                 matrix_free.get_cell_active_fe_index(range))
{}



template <int dim,
          int fe_degree,
          int n_q_points_1d,
          int n_components_,
          typename Number,
          typename VectorizedArrayType>
inline FEEvaluation<dim,
                    fe_degree,
                    n_q_points_1d,
                    n_components_,
                    Number,
                    VectorizedArrayType>::
  FEEvaluation(const Mapping<dim> &      mapping,
               const FiniteElement<dim> &fe,
               const Quadrature<1> &     quadrature,
               const UpdateFlags         update_flags,
               const unsigned int        first_selected_component)
  : BaseClass(mapping,
              fe,
              quadrature,
              update_flags,
              first_selected_component,
              nullptr)
  , dofs_per_component(this->data->dofs_per_component_on_cell)
  , dofs_per_cell(this->data->dofs_per_component_on_cell * n_components_)
  , n_q_points(this->data->n_q_points)
{
  check_template_arguments(numbers::invalid_unsigned_int, 0);
}



template <int dim,
          int fe_degree,
          int n_q_points_1d,
          int n_components_,
          typename Number,
          typename VectorizedArrayType>
inline FEEvaluation<dim,
                    fe_degree,
                    n_q_points_1d,
                    n_components_,
                    Number,
                    VectorizedArrayType>::
  FEEvaluation(const FiniteElement<dim> &fe,
               const Quadrature<1> &     quadrature,
               const UpdateFlags         update_flags,
               const unsigned int        first_selected_component)
  : BaseClass(StaticMappingQ1<dim>::mapping,
              fe,
              quadrature,
              update_flags,
              first_selected_component,
              nullptr)
  , dofs_per_component(this->data->dofs_per_component_on_cell)
  , dofs_per_cell(this->data->dofs_per_component_on_cell * n_components_)
  , n_q_points(this->data->n_q_points)
{
  check_template_arguments(numbers::invalid_unsigned_int, 0);
}



template <int dim,
          int fe_degree,
          int n_q_points_1d,
          int n_components_,
          typename Number,
          typename VectorizedArrayType>
inline FEEvaluation<dim,
                    fe_degree,
                    n_q_points_1d,
                    n_components_,
                    Number,
                    VectorizedArrayType>::
  FEEvaluation(
    const FiniteElement<dim> &                                           fe,
    const FEEvaluationBaseData<dim, Number, false, VectorizedArrayType> &other,
    const unsigned int first_selected_component)
  : BaseClass(other.mapped_geometry->get_fe_values().get_mapping(),
              fe,
              other.mapped_geometry->get_quadrature(),
              other.mapped_geometry->get_fe_values().get_update_flags(),
              first_selected_component,
              &other)
  , dofs_per_component(this->data->dofs_per_component_on_cell)
  , dofs_per_cell(this->data->dofs_per_component_on_cell * n_components_)
  , n_q_points(this->data->n_q_points)
{
  check_template_arguments(numbers::invalid_unsigned_int, 0);
}



template <int dim,
          int fe_degree,
          int n_q_points_1d,
          int n_components_,
          typename Number,
          typename VectorizedArrayType>
inline FEEvaluation<dim,
                    fe_degree,
                    n_q_points_1d,
                    n_components_,
                    Number,
                    VectorizedArrayType>::FEEvaluation(const FEEvaluation
                                                         &other)
  : BaseClass(other)
  , dofs_per_component(this->data->dofs_per_component_on_cell)
  , dofs_per_cell(this->data->dofs_per_component_on_cell * n_components_)
  , n_q_points(this->data->n_q_points)
{
  check_template_arguments(numbers::invalid_unsigned_int, 0);
}



template <int dim,
          int fe_degree,
          int n_q_points_1d,
          int n_components_,
          typename Number,
          typename VectorizedArrayType>
inline FEEvaluation<dim,
                    fe_degree,
                    n_q_points_1d,
                    n_components_,
                    Number,
                    VectorizedArrayType> &
FEEvaluation<dim,
             fe_degree,
             n_q_points_1d,
             n_components_,
             Number,
             VectorizedArrayType>::operator=(const FEEvaluation &other)
{
  BaseClass::operator=(other);
  check_template_arguments(numbers::invalid_unsigned_int, 0);
  return *this;
}



template <int dim,
          int fe_degree,
          int n_q_points_1d,
          int n_components_,
          typename Number,
          typename VectorizedArrayType>
inline void
FEEvaluation<dim,
             fe_degree,
             n_q_points_1d,
             n_components_,
             Number,
             VectorizedArrayType>::
  check_template_arguments(const unsigned int dof_no,
                           const unsigned int first_selected_component)
{
  (void)dof_no;
  (void)first_selected_component;

#  ifdef DEBUG
  // print error message when the dimensions do not match. Propose a possible
  // fix
  if ((static_cast<unsigned int>(fe_degree) != numbers::invalid_unsigned_int &&
       static_cast<unsigned int>(fe_degree) !=
         this->data->data.front().fe_degree) ||
      n_q_points != this->n_quadrature_points)
    {
      std::string message =
        "-------------------------------------------------------\n";
      message += "Illegal arguments in constructor/wrong template arguments!\n";
      message += "    Called -->   FEEvaluation<dim,";
      message += Utilities::int_to_string(fe_degree) + ",";
      message += Utilities::int_to_string(n_q_points_1d);
      message += "," + Utilities::int_to_string(n_components);
      message += ",Number>(data";
      if (first_selected_component != numbers::invalid_unsigned_int)
        {
          message += ", " + Utilities::int_to_string(dof_no) + ", ";
          message += Utilities::int_to_string(this->quad_no) + ", ";
          message += Utilities::int_to_string(first_selected_component);
        }
      message += ")\n";

      // check whether some other vector component has the correct number of
      // points
      unsigned int proposed_dof_comp  = numbers::invalid_unsigned_int,
                   proposed_fe_comp   = numbers::invalid_unsigned_int,
                   proposed_quad_comp = numbers::invalid_unsigned_int;
      if (dof_no != numbers::invalid_unsigned_int)
        {
          if (static_cast<unsigned int>(fe_degree) ==
              this->data->data.front().fe_degree)
            {
              proposed_dof_comp = dof_no;
              proposed_fe_comp  = first_selected_component;
            }
          else
            for (unsigned int no = 0; no < this->matrix_info->n_components();
                 ++no)
              for (unsigned int nf = 0;
                   nf < this->matrix_info->n_base_elements(no);
                   ++nf)
                if (this->matrix_info
                      ->get_shape_info(no, 0, nf, this->active_fe_index, 0)
                      .data.front()
                      .fe_degree == static_cast<unsigned int>(fe_degree))
                  {
                    proposed_dof_comp = no;
                    proposed_fe_comp  = nf;
                    break;
                  }
          if (n_q_points ==
              this->mapping_data->descriptor[this->active_quad_index]
                .n_q_points)
            proposed_quad_comp = this->quad_no;
          else
            for (unsigned int no = 0;
                 no < this->matrix_info->get_mapping_info().cell_data.size();
                 ++no)
              if (this->matrix_info->get_mapping_info()
                    .cell_data[no]
                    .descriptor[this->active_quad_index]
                    .n_q_points == n_q_points)
                {
                  proposed_quad_comp = no;
                  break;
                }
        }
      if (proposed_dof_comp != numbers::invalid_unsigned_int &&
          proposed_quad_comp != numbers::invalid_unsigned_int)
        {
          if (proposed_dof_comp != first_selected_component)
            message += "Wrong vector component selection:\n";
          else
            message += "Wrong quadrature formula selection:\n";
          message += "    Did you mean FEEvaluation<dim,";
          message += Utilities::int_to_string(fe_degree) + ",";
          message += Utilities::int_to_string(n_q_points_1d);
          message += "," + Utilities::int_to_string(n_components);
          message += ",Number>(data";
          if (dof_no != numbers::invalid_unsigned_int)
            {
              message +=
                ", " + Utilities::int_to_string(proposed_dof_comp) + ", ";
              message += Utilities::int_to_string(proposed_quad_comp) + ", ";
              message += Utilities::int_to_string(proposed_fe_comp);
            }
          message += ")?\n";
          std::string correct_pos;
          if (proposed_dof_comp != dof_no)
            correct_pos = " ^ ";
          else
            correct_pos = "   ";
          if (proposed_quad_comp != this->quad_no)
            correct_pos += " ^ ";
          else
            correct_pos += "   ";
          if (proposed_fe_comp != first_selected_component)
            correct_pos += " ^\n";
          else
            correct_pos += "  \n";
          message += "                                                     " +
                     correct_pos;
        }
      // ok, did not find the numbers specified by the template arguments in
      // the given list. Suggest correct template arguments
      const unsigned int proposed_n_q_points_1d = static_cast<unsigned int>(
        std::pow(1.001 * this->n_quadrature_points, 1. / dim));
      message += "Wrong template arguments:\n";
      message += "    Did you mean FEEvaluation<dim,";
      message +=
        Utilities::int_to_string(this->data->data.front().fe_degree) + ",";
      message += Utilities::int_to_string(proposed_n_q_points_1d);
      message += "," + Utilities::int_to_string(n_components);
      message += ",Number>(data";
      if (dof_no != numbers::invalid_unsigned_int)
        {
          message += ", " + Utilities::int_to_string(dof_no) + ", ";
          message += Utilities::int_to_string(this->quad_no);
          message += ", " + Utilities::int_to_string(first_selected_component);
        }
      message += ")?\n";
      std::string correct_pos;
      if (this->data->data.front().fe_degree !=
          static_cast<unsigned int>(fe_degree))
        correct_pos = " ^";
      else
        correct_pos = "  ";
      if (proposed_n_q_points_1d != n_q_points_1d)
        correct_pos += " ^\n";
      else
        correct_pos += "  \n";
      message += "                                 " + correct_pos;

      Assert(static_cast<unsigned int>(fe_degree) ==
                 this->data->data.front().fe_degree &&
               n_q_points == this->n_quadrature_points,
             ExcMessage(message));
    }
  if (dof_no != numbers::invalid_unsigned_int)
    AssertDimension(
      n_q_points,
      this->mapping_data->descriptor[this->active_quad_index].n_q_points);
#  endif
}



template <int dim,
          int fe_degree,
          int n_q_points_1d,
          int n_components_,
          typename Number,
          typename VectorizedArrayType>
inline void
FEEvaluation<dim,
             fe_degree,
             n_q_points_1d,
             n_components_,
             Number,
             VectorizedArrayType>::reinit(const unsigned int cell_index)
{
  Assert(this->mapped_geometry == nullptr,
         ExcMessage("FEEvaluation was initialized without a matrix-free object."
                    " Integer indexing is not possible"));
  if (this->mapped_geometry != nullptr)
    return;

  Assert(this->dof_info != nullptr, ExcNotInitialized());
  Assert(this->mapping_data != nullptr, ExcNotInitialized());
  this->cell = cell_index;
  this->cell_type =
    this->matrix_info->get_mapping_info().get_cell_type(cell_index);

  const unsigned int offsets =
    this->mapping_data->data_index_offsets[cell_index];
  this->jacobian = &this->mapping_data->jacobians[0][offsets];
  this->J_value  = &this->mapping_data->JxW_values[offsets];

#  ifdef DEBUG
  this->dof_values_initialized     = false;
  this->values_quad_initialized    = false;
  this->gradients_quad_initialized = false;
  this->hessians_quad_initialized  = false;
#  endif
}



template <int dim,
          int fe_degree,
          int n_q_points_1d,
          int n_components_,
          typename Number,
          typename VectorizedArrayType>
template <bool level_dof_access>
inline void
FEEvaluation<dim,
             fe_degree,
             n_q_points_1d,
             n_components_,
             Number,
             VectorizedArrayType>::
  reinit(const TriaIterator<DoFCellAccessor<dim, dim, level_dof_access>> &cell)
{
  Assert(this->matrix_info == nullptr,
         ExcMessage("Cannot use initialization from cell iterator if "
                    "initialized from MatrixFree object. Use variant for "
                    "on the fly computation with arguments as for FEValues "
                    "instead"));
  Assert(this->mapped_geometry.get() != nullptr, ExcNotInitialized());
  this->mapped_geometry->reinit(
    static_cast<typename Triangulation<dim>::cell_iterator>(cell));
  this->local_dof_indices.resize(cell->get_fe().n_dofs_per_cell());
  if (level_dof_access)
    cell->get_mg_dof_indices(this->local_dof_indices);
  else
    cell->get_dof_indices(this->local_dof_indices);
}



template <int dim,
          int fe_degree,
          int n_q_points_1d,
          int n_components_,
          typename Number,
          typename VectorizedArrayType>
inline void
FEEvaluation<dim,
             fe_degree,
             n_q_points_1d,
             n_components_,
             Number,
             VectorizedArrayType>::
  reinit(const typename Triangulation<dim>::cell_iterator &cell)
{
  Assert(this->matrix_info == 0,
         ExcMessage("Cannot use initialization from cell iterator if "
                    "initialized from MatrixFree object. Use variant for "
                    "on the fly computation with arguments as for FEValues "
                    "instead"));
  Assert(this->mapped_geometry.get() != 0, ExcNotInitialized());
  this->mapped_geometry->reinit(cell);
}



template <int dim,
          int fe_degree,
          int n_q_points_1d,
          int n_components_,
          typename Number,
          typename VectorizedArrayType>
inline Point<dim, VectorizedArrayType>
FEEvaluation<dim,
             fe_degree,
             n_q_points_1d,
             n_components_,
             Number,
             VectorizedArrayType>::quadrature_point(const unsigned int q) const
{
  if (this->matrix_info == nullptr)
    {
      Assert((this->mapped_geometry->get_fe_values().get_update_flags() |
              update_quadrature_points),
             internal::ExcMatrixFreeAccessToUninitializedMappingField(
               "update_quadrature_points"));
    }
  else
    {
      Assert(this->mapping_data->quadrature_point_offsets.empty() == false,
             internal::ExcMatrixFreeAccessToUninitializedMappingField(
               "update_quadrature_points"));
    }

  AssertIndexRange(q, n_q_points);

  const Point<dim, VectorizedArrayType> *quadrature_points =
    &this->mapping_data->quadrature_points
       [this->mapping_data->quadrature_point_offsets[this->cell]];

  // Cartesian/affine mesh: only first vertex of cell is stored, we must
  // compute it through the Jacobian (which is stored in non-inverted and
  // non-transposed form as index '1' in the jacobian field)
  if (this->cell_type <= internal::MatrixFreeFunctions::affine)
    {
      Assert(this->jacobian != nullptr, ExcNotInitialized());
      Point<dim, VectorizedArrayType> point = quadrature_points[0];

      const Tensor<2, dim, VectorizedArrayType> &jac = this->jacobian[1];
      if (this->cell_type == internal::MatrixFreeFunctions::cartesian)
        for (unsigned int d = 0; d < dim; ++d)
          point[d] += jac[d][d] * static_cast<Number>(
                                    this->descriptor->quadrature.point(q)[d]);
      else
        for (unsigned int d = 0; d < dim; ++d)
          for (unsigned int e = 0; e < dim; ++e)
            point[d] += jac[d][e] * static_cast<Number>(
                                      this->descriptor->quadrature.point(q)[e]);
      return point;
    }
  else
    return quadrature_points[q];
}



template <int dim,
          int fe_degree,
          int n_q_points_1d,
          int n_components_,
          typename Number,
          typename VectorizedArrayType>
inline void
FEEvaluation<dim,
             fe_degree,
             n_q_points_1d,
             n_components_,
             Number,
             VectorizedArrayType>::evaluate(const bool evaluate_values,
                                            const bool evaluate_gradients,
                                            const bool evaluate_hessians)
{
#  ifdef DEBUG
  Assert(this->dof_values_initialized == true,
         internal::ExcAccessToUninitializedField());
#  endif
  evaluate(this->values_dofs[0],
           evaluate_values,
           evaluate_gradients,
           evaluate_hessians);
}


template <int dim,
          int fe_degree,
          int n_q_points_1d,
          int n_components_,
          typename Number,
          typename VectorizedArrayType>
inline void
FEEvaluation<dim,
             fe_degree,
             n_q_points_1d,
             n_components_,
             Number,
             VectorizedArrayType>::
  evaluate(const EvaluationFlags::EvaluationFlags evaluation_flags)
{
#  ifdef DEBUG
  Assert(this->dof_values_initialized == true,
         internal::ExcAccessToUninitializedField());
#  endif
  evaluate(this->values_dofs[0], evaluation_flags);
}



template <int dim,
          int fe_degree,
          int n_q_points_1d,
          int n_components_,
          typename Number,
          typename VectorizedArrayType>
inline void
FEEvaluation<dim,
             fe_degree,
             n_q_points_1d,
             n_components_,
             Number,
             VectorizedArrayType>::evaluate(const VectorizedArrayType
                                              *        values_array,
                                            const bool evaluate_values,
                                            const bool evaluate_gradients,
                                            const bool evaluate_hessians)
{
  const EvaluationFlags::EvaluationFlags flag =
    ((evaluate_values) ? EvaluationFlags::values : EvaluationFlags::nothing) |
    ((evaluate_gradients) ? EvaluationFlags::gradients :
                            EvaluationFlags::nothing) |
    ((evaluate_hessians) ? EvaluationFlags::hessians :
                           EvaluationFlags::nothing);

  evaluate(values_array, flag);
}



template <int dim,
          int fe_degree,
          int n_q_points_1d,
          int n_components_,
          typename Number,
          typename VectorizedArrayType>
inline void
FEEvaluation<dim,
             fe_degree,
             n_q_points_1d,
             n_components_,
             Number,
             VectorizedArrayType>::
  evaluate(const VectorizedArrayType *            values_array,
           const EvaluationFlags::EvaluationFlags evaluation_flags)
{
  if (fe_degree > -1)
    SelectEvaluator<dim, fe_degree, n_q_points_1d, VectorizedArrayType>::
      evaluate(n_components,
               evaluation_flags,
               *this->data,
               const_cast<VectorizedArrayType *>(values_array),
               this->values_quad,
               this->gradients_quad,
               this->hessians_quad,
               this->scratch_data);
  else
    internal::FEEvaluationFactory<dim, Number, VectorizedArrayType>::evaluate(
      n_components,
      evaluation_flags,
      *this->data,
      const_cast<VectorizedArrayType *>(values_array),
      this->values_quad,
      this->gradients_quad,
      this->hessians_quad,
      this->scratch_data);

#  ifdef DEBUG
  if (evaluation_flags & EvaluationFlags::values)
    this->values_quad_initialized = true;
  if (evaluation_flags & EvaluationFlags::gradients)
    this->gradients_quad_initialized = true;
  if (evaluation_flags & EvaluationFlags::hessians)
    this->hessians_quad_initialized = true;
#  endif
}



template <int dim,
          int fe_degree,
          int n_q_points_1d,
          int n_components_,
          typename Number,
          typename VectorizedArrayType>
template <typename VectorType>
inline void
FEEvaluation<
  dim,
  fe_degree,
  n_q_points_1d,
  n_components_,
  Number,
  VectorizedArrayType>::gather_evaluate(const VectorType &input_vector,
                                        const bool        evaluate_values,
                                        const bool        evaluate_gradients,
                                        const bool        evaluate_hessians)
{
  const EvaluationFlags::EvaluationFlags flag =
    ((evaluate_values) ? EvaluationFlags::values : EvaluationFlags::nothing) |
    ((evaluate_gradients) ? EvaluationFlags::gradients :
                            EvaluationFlags::nothing) |
    ((evaluate_hessians) ? EvaluationFlags::hessians :
                           EvaluationFlags::nothing);

  gather_evaluate(input_vector, flag);
}


namespace internal
{
  /**
   * 对标准向量（有begin()方法的）的实现。
   *
   */
  template <typename Number,
            typename VectorizedArrayType,
            typename VectorType,
            typename T,
            typename std::enable_if<
              internal::has_begin<VectorType>::value &&
                std::is_same<decltype(std::declval<VectorType>().begin()),
                             Number *>::value,
              VectorType>::type * = nullptr>
  bool
  try_gather_evaluate_inplace(
    T                                             phi,
    const VectorType &                            input_vector,
    const unsigned int                            cell,
    const unsigned int                            active_fe_index,
    const unsigned int                            first_selected_component,
    const internal::MatrixFreeFunctions::DoFInfo *dof_info,
    const EvaluationFlags::EvaluationFlags        evaluation_flag)
  {
    // If the index storage is interleaved and contiguous and the vector storage
    // has the correct alignment, we can directly pass the pointer into the
    // vector to the evaluate() call, without reading the vector entries into a
    // separate data field. This saves some operations.
    if (std::is_same<typename VectorType::value_type, Number>::value &&
        dof_info->index_storage_variants
            [internal::MatrixFreeFunctions::DoFInfo::dof_access_cell][cell] ==
          internal::MatrixFreeFunctions::DoFInfo::IndexStorageVariants::
            interleaved_contiguous &&
        reinterpret_cast<std::size_t>(
          input_vector.begin() +
          dof_info->dof_indices_contiguous
            [internal::MatrixFreeFunctions::DoFInfo::dof_access_cell]
            [cell * VectorizedArrayType::size()]) %
            sizeof(VectorizedArrayType) ==
          0)
      {
        const VectorizedArrayType *vec_values =
          reinterpret_cast<const VectorizedArrayType *>(
            input_vector.begin() +
            dof_info->dof_indices_contiguous
              [internal::MatrixFreeFunctions::DoFInfo::dof_access_cell]
              [cell * VectorizedArrayType::size()] +
            dof_info->component_dof_indices_offset[active_fe_index]
                                                  [first_selected_component] *
              VectorizedArrayType::size());

        phi->evaluate(vec_values, evaluation_flag);

        return true;
      }

    return false;
  }

  /**
   * 对块状向量的实现。
   *
   */
  template <typename Number,
            typename VectorizedArrayType,
            typename VectorType,
            typename T,
            typename std::enable_if<
              !internal::has_begin<VectorType>::value ||
                !std::is_same<decltype(std::declval<VectorType>().begin()),
                              Number *>::value,
              VectorType>::type * = nullptr>
  bool
  try_gather_evaluate_inplace(T,
                              const VectorType &,
                              const unsigned int,
                              const unsigned int,
                              const unsigned int,
                              const internal::MatrixFreeFunctions::DoFInfo *,
                              const EvaluationFlags::EvaluationFlags)
  {
    return false;
  }

  /**
   * 对有begin()方法的向量的实现。
   *
   */
  template <int dim,
            int fe_degree,
            int n_q_points_1d,
            typename Number,
            typename VectorizedArrayType,
            typename VectorType,
            typename std::enable_if<
              internal::has_begin<VectorType>::value &&
                std::is_same<decltype(std::declval<VectorType>().begin()),
                             Number *>::value,
              VectorType>::type * = nullptr>
  bool
  try_integrate_scatter_inplace(
    VectorType &                                  destination,
    const unsigned int                            cell,
    const unsigned int                            n_components,
    const unsigned int                            active_fe_index,
    const unsigned int                            first_selected_component,
    const internal::MatrixFreeFunctions::DoFInfo *dof_info,
    VectorizedArrayType *                         values_quad,
    VectorizedArrayType *                         gradients_quad,
    VectorizedArrayType *                         scratch_data,
    const internal::MatrixFreeFunctions::ShapeInfo<VectorizedArrayType> *data,
    const EvaluationFlags::EvaluationFlags integration_flag)
  {
    // If the index storage is interleaved and contiguous and the vector storage
    // has the correct alignment, we can directly pass the pointer into the
    // vector to the integrate() call, without writing temporary results into a
    // separate data field that will later be added into the vector. This saves
    // some operations.
    if (std::is_same<typename VectorType::value_type, Number>::value &&
        dof_info->index_storage_variants
            [internal::MatrixFreeFunctions::DoFInfo::dof_access_cell][cell] ==
          internal::MatrixFreeFunctions::DoFInfo::IndexStorageVariants::
            interleaved_contiguous &&
        reinterpret_cast<std::size_t>(
          destination.begin() +
          dof_info->dof_indices_contiguous
            [internal::MatrixFreeFunctions::DoFInfo::dof_access_cell]
            [cell * VectorizedArrayType::size()]) %
            sizeof(VectorizedArrayType) ==
          0)
      {
        VectorizedArrayType *vec_values =
          reinterpret_cast<VectorizedArrayType *>(
            destination.begin() +
            dof_info->dof_indices_contiguous
              [internal::MatrixFreeFunctions::DoFInfo::dof_access_cell]
              [cell * VectorizedArrayType::size()] +
            dof_info->component_dof_indices_offset[active_fe_index]
                                                  [first_selected_component] *
              VectorizedArrayType::size());
        if (fe_degree > -1)
          SelectEvaluator<dim, fe_degree, n_q_points_1d, VectorizedArrayType>::
            integrate(n_components,
                      integration_flag,
                      *data,
                      vec_values,
                      values_quad,
                      gradients_quad,
                      scratch_data,
                      true);
        else
          FEEvaluationFactory<dim, Number, VectorizedArrayType>::integrate(
            n_components,
            integration_flag,
            *data,
            vec_values,
            values_quad,
            gradients_quad,
            scratch_data,
            true);

        return true;
      }

    return false;
  }

  /**
   * 对所有其他向量（如块向量）的实现。
   *
   */
  template <int dim,
            int fe_degree,
            int n_q_points_1d,
            typename Number,
            typename VectorizedArrayType,
            typename VectorType,
            typename std::enable_if<
              !internal::has_begin<VectorType>::value ||
                !std::is_same<decltype(std::declval<VectorType>().begin()),
                              Number *>::value,
              VectorType>::type * = nullptr>
  bool
  try_integrate_scatter_inplace(
    VectorType &,
    const unsigned int,
    const unsigned int,
    const unsigned int,
    const unsigned int,
    const internal::MatrixFreeFunctions::DoFInfo *,
    const VectorizedArrayType *,
    const VectorizedArrayType *,
    const VectorizedArrayType *,
    const internal::MatrixFreeFunctions::ShapeInfo<VectorizedArrayType> *,
    const EvaluationFlags::EvaluationFlags)
  {
    return false;
  }
} // namespace internal



template <int dim,
          int fe_degree,
          int n_q_points_1d,
          int n_components_,
          typename Number,
          typename VectorizedArrayType>
template <typename VectorType>
inline void
FEEvaluation<dim,
             fe_degree,
             n_q_points_1d,
             n_components_,
             Number,
             VectorizedArrayType>::
  gather_evaluate(const VectorType &                     input_vector,
                  const EvaluationFlags::EvaluationFlags evaluation_flag)
{
  if (internal::try_gather_evaluate_inplace<Number, VectorizedArrayType>(
        this,
        input_vector,
        this->cell,
        this->active_fe_index,
        this->first_selected_component,
        this->dof_info,
        evaluation_flag) == false)
    {
      this->read_dof_values(input_vector);
      evaluate(this->begin_dof_values(), evaluation_flag);
    }
}



template <int dim,
          int fe_degree,
          int n_q_points_1d,
          int n_components_,
          typename Number,
          typename VectorizedArrayType>
inline void
FEEvaluation<dim,
             fe_degree,
             n_q_points_1d,
             n_components_,
             Number,
             VectorizedArrayType>::integrate(const bool integrate_values,
                                             const bool integrate_gradients)
{
  integrate(integrate_values, integrate_gradients, this->values_dofs[0]);

#  ifdef DEBUG
  this->dof_values_initialized = true;
#  endif
}



template <int dim,
          int fe_degree,
          int n_q_points_1d,
          int n_components_,
          typename Number,
          typename VectorizedArrayType>
inline void
FEEvaluation<dim,
             fe_degree,
             n_q_points_1d,
             n_components_,
             Number,
             VectorizedArrayType>::
  integrate(const EvaluationFlags::EvaluationFlags integration_flag)
{
  integrate(integration_flag, this->values_dofs[0]);

#  ifdef DEBUG
  this->dof_values_initialized = true;
#  endif
}



template <int dim,
          int fe_degree,
          int n_q_points_1d,
          int n_components_,
          typename Number,
          typename VectorizedArrayType>
inline void
FEEvaluation<dim,
             fe_degree,
             n_q_points_1d,
             n_components_,
             Number,
             VectorizedArrayType>::integrate(const bool integrate_values,
                                             const bool integrate_gradients,
                                             VectorizedArrayType *values_array)
{
  EvaluationFlags::EvaluationFlags flag =
    (integrate_values ? EvaluationFlags::values : EvaluationFlags::nothing) |
    (integrate_gradients ? EvaluationFlags::gradients :
                           EvaluationFlags::nothing);
  integrate(flag, values_array);
}



template <int dim,
          int fe_degree,
          int n_q_points_1d,
          int n_components_,
          typename Number,
          typename VectorizedArrayType>
inline void
FEEvaluation<dim,
             fe_degree,
             n_q_points_1d,
             n_components_,
             Number,
             VectorizedArrayType>::
  integrate(const EvaluationFlags::EvaluationFlags integration_flag,
            VectorizedArrayType *                  values_array)
{
#  ifdef DEBUG
  if (integration_flag & EvaluationFlags::values)
    Assert(this->values_quad_submitted == true,
           internal::ExcAccessToUninitializedField());
  if (integration_flag & EvaluationFlags::gradients)
    Assert(this->gradients_quad_submitted == true,
           internal::ExcAccessToUninitializedField());
#  endif
  Assert(this->matrix_info != nullptr ||
           this->mapped_geometry->is_initialized(),
         ExcNotInitialized());

  Assert(
    (integration_flag &
     ~(EvaluationFlags::values | EvaluationFlags::gradients)) == 0,
    ExcMessage(
      "Only EvaluationFlags::values and EvaluationFlags::gradients are supported."));

  if (fe_degree > -1)
    SelectEvaluator<dim, fe_degree, n_q_points_1d, VectorizedArrayType>::
      integrate(n_components,
                integration_flag,
                *this->data,
                values_array,
                this->values_quad,
                this->gradients_quad,
                this->scratch_data,
                false);
  else
    internal::FEEvaluationFactory<dim, Number, VectorizedArrayType>::integrate(
      n_components,
      integration_flag,
      *this->data,
      values_array,
      this->values_quad,
      this->gradients_quad,
      this->scratch_data,
      false);

#  ifdef DEBUG
  this->dof_values_initialized = true;
#  endif
}



template <int dim,
          int fe_degree,
          int n_q_points_1d,
          int n_components_,
          typename Number,
          typename VectorizedArrayType>
template <typename VectorType>
inline void
FEEvaluation<
  dim,
  fe_degree,
  n_q_points_1d,
  n_components_,
  Number,
  VectorizedArrayType>::integrate_scatter(const bool  integrate_values,
                                          const bool  integrate_gradients,
                                          VectorType &destination)
{
  const EvaluationFlags::EvaluationFlags flag =
    ((integrate_values) ? EvaluationFlags::values : EvaluationFlags::nothing) |
    ((integrate_gradients) ? EvaluationFlags::gradients :
                             EvaluationFlags::nothing);

  integrate_scatter(flag, destination);
}



template <int dim,
          int fe_degree,
          int n_q_points_1d,
          int n_components_,
          typename Number,
          typename VectorizedArrayType>
template <typename VectorType>
inline void
FEEvaluation<dim,
             fe_degree,
             n_q_points_1d,
             n_components_,
             Number,
             VectorizedArrayType>::
  integrate_scatter(const EvaluationFlags::EvaluationFlags integration_flag,
                    VectorType &                           destination)
{
  if (internal::try_integrate_scatter_inplace<dim,
                                              fe_degree,
                                              n_q_points_1d,
                                              Number,
                                              VectorizedArrayType>(
        destination,
        this->cell,
        n_components,
        this->active_fe_index,
        this->first_selected_component,
        this->dof_info,
        this->values_quad,
        this->gradients_quad,
        this->scratch_data,
        this->data,
        integration_flag) == false)
    {
      integrate(integration_flag, this->begin_dof_values());
      this->distribute_local_to_global(destination);
    }
}



 /*-------------------------- FEFaceEvaluation ---------------------------*/ 



template <int dim,
          int fe_degree,
          int n_q_points_1d,
          int n_components_,
          typename Number,
          typename VectorizedArrayType>
inline FEFaceEvaluation<dim,
                        fe_degree,
                        n_q_points_1d,
                        n_components_,
                        Number,
                        VectorizedArrayType>::
  FEFaceEvaluation(
    const MatrixFree<dim, Number, VectorizedArrayType> &matrix_free,
    const bool                                          is_interior_face,
    const unsigned int                                  dof_no,
    const unsigned int                                  quad_no,
    const unsigned int first_selected_component,
    const unsigned int active_fe_index,
    const unsigned int active_quad_index,
    const unsigned int face_type)
  : BaseClass(matrix_free,
              dof_no,
              first_selected_component,
              quad_no,
              fe_degree,
              static_n_q_points,
              is_interior_face,
              active_fe_index,
              active_quad_index,
              face_type)
  , dofs_per_component(this->data->dofs_per_component_on_cell)
  , dofs_per_cell(this->data->dofs_per_component_on_cell * n_components_)
  , n_q_points(this->n_quadrature_points)
{}



template <int dim,
          int fe_degree,
          int n_q_points_1d,
          int n_components_,
          typename Number,
          typename VectorizedArrayType>
inline FEFaceEvaluation<dim,
                        fe_degree,
                        n_q_points_1d,
                        n_components_,
                        Number,
                        VectorizedArrayType>::
  FEFaceEvaluation(
    const MatrixFree<dim, Number, VectorizedArrayType> &matrix_free,
    const std::pair<unsigned int, unsigned int> &       range,
    const bool                                          is_interior_face,
    const unsigned int                                  dof_no,
    const unsigned int                                  quad_no,
    const unsigned int first_selected_component)
  : FEFaceEvaluation(matrix_free,
                     is_interior_face,
                     dof_no,
                     quad_no,
                     first_selected_component,
                     matrix_free.get_face_active_fe_index(range,
                                                          is_interior_face),
                     numbers::invalid_unsigned_int,
                     matrix_free.get_face_info(range.first).face_type)
{}



template <int dim,
          int fe_degree,
          int n_q_points_1d,
          int n_components_,
          typename Number,
          typename VectorizedArrayType>
inline void
FEFaceEvaluation<dim,
                 fe_degree,
                 n_q_points_1d,
                 n_components_,
                 Number,
                 VectorizedArrayType>::reinit(const unsigned int face_index)
{
  Assert(this->mapped_geometry == nullptr,
         ExcMessage("FEEvaluation was initialized without a matrix-free object."
                    " Integer indexing is not possible"));
  if (this->mapped_geometry != nullptr)
    return;

  this->cell = face_index;
  this->dof_access_index =
    this->is_interior_face ?
      internal::MatrixFreeFunctions::DoFInfo::dof_access_face_interior :
      internal::MatrixFreeFunctions::DoFInfo::dof_access_face_exterior;
  Assert(this->mapping_data != nullptr, ExcNotInitialized());
  const internal::MatrixFreeFunctions::FaceToCellTopology<
    VectorizedArrayType::size()> &faces =
    this->matrix_info->get_face_info(face_index);
  if (face_index >=
        this->matrix_info->get_task_info().face_partition_data.back() &&
      face_index <
        this->matrix_info->get_task_info().boundary_partition_data.back())
    Assert(this->is_interior_face,
           ExcMessage(
             "Boundary faces do not have a neighbor. When looping over "
             "boundary faces use FEFaceEvaluation with the parameter "
             "is_interior_face set to true. "));

  this->face_no =
    (this->is_interior_face ? faces.interior_face_no : faces.exterior_face_no);
  this->subface_index = this->is_interior_face == true ?
                          GeometryInfo<dim>::max_children_per_cell :
                          faces.subface_index;

  // First check if interior or exterior cell has non-standard orientation
  // (i.e. the third bit is one or not). Then set zero if this cell has
  // standard-orientation else copy the first three bits
  // (which is equivalent to modulo 8). See also the documentation of
  // internal::MatrixFreeFunctions::FaceToCellTopology::face_orientation.
  this->face_orientation =
    (this->is_interior_face == (faces.face_orientation >= 8)) ?
      (faces.face_orientation % 8) :
      0;

  this->cell_type = this->matrix_info->get_mapping_info().face_type[face_index];
  const unsigned int offsets =
    this->mapping_data->data_index_offsets[face_index];
  this->J_value        = &this->mapping_data->JxW_values[offsets];
  this->normal_vectors = &this->mapping_data->normal_vectors[offsets];
  this->jacobian =
    &this->mapping_data->jacobians[!this->is_interior_face][offsets];
  this->normal_x_jacobian =
    &this->mapping_data
       ->normals_times_jacobians[!this->is_interior_face][offsets];

#  ifdef DEBUG
  this->dof_values_initialized     = false;
  this->values_quad_initialized    = false;
  this->gradients_quad_initialized = false;
  this->hessians_quad_initialized  = false;
#  endif
}



template <int dim,
          int fe_degree,
          int n_q_points_1d,
          int n_components_,
          typename Number,
          typename VectorizedArrayType>
inline void
FEFaceEvaluation<dim,
                 fe_degree,
                 n_q_points_1d,
                 n_components_,
                 Number,
                 VectorizedArrayType>::reinit(const unsigned int cell_index,
                                              const unsigned int face_number)
{
  Assert(
    this->quad_no <
      this->matrix_info->get_mapping_info().face_data_by_cells.size(),
    ExcMessage(
      "You must set MatrixFree::AdditionalData::mapping_update_flags_faces_by_cells to use the present reinit method."));
  AssertIndexRange(face_number, GeometryInfo<dim>::faces_per_cell);
  AssertIndexRange(cell_index,
                   this->matrix_info->get_mapping_info().cell_type.size());
  Assert(this->mapped_geometry == nullptr,
         ExcMessage("FEEvaluation was initialized without a matrix-free object."
                    " Integer indexing is not possible"));
  if (this->mapped_geometry != nullptr)
    return;
  Assert(this->matrix_info != nullptr, ExcNotInitialized());

  this->cell_type = this->matrix_info->get_mapping_info().cell_type[cell_index];
  this->cell      = cell_index;
  this->face_orientation = 0;
  this->subface_index    = GeometryInfo<dim>::max_children_per_cell;
  this->face_no          = face_number;
  this->dof_access_index =
    internal::MatrixFreeFunctions::DoFInfo::dof_access_cell;

  const unsigned int offsets =
    this->matrix_info->get_mapping_info()
      .face_data_by_cells[this->quad_no]
      .data_index_offsets[cell_index * GeometryInfo<dim>::faces_per_cell +
                          face_number];
  AssertIndexRange(offsets,
                   this->matrix_info->get_mapping_info()
                     .face_data_by_cells[this->quad_no]
                     .JxW_values.size());
  this->J_value = &this->matrix_info->get_mapping_info()
                     .face_data_by_cells[this->quad_no]
                     .JxW_values[offsets];
  this->normal_vectors = &this->matrix_info->get_mapping_info()
                            .face_data_by_cells[this->quad_no]
                            .normal_vectors[offsets];
  this->jacobian = &this->matrix_info->get_mapping_info()
                      .face_data_by_cells[this->quad_no]
                      .jacobians[!this->is_interior_face][offsets];
  this->normal_x_jacobian =
    &this->matrix_info->get_mapping_info()
       .face_data_by_cells[this->quad_no]
       .normals_times_jacobians[!this->is_interior_face][offsets];

#  ifdef DEBUG
  this->dof_values_initialized     = false;
  this->values_quad_initialized    = false;
  this->gradients_quad_initialized = false;
  this->hessians_quad_initialized  = false;
#  endif
}



template <int dim,
          int fe_degree,
          int n_q_points_1d,
          int n_components,
          typename Number,
          typename VectorizedArrayType>
inline void
FEFaceEvaluation<dim,
                 fe_degree,
                 n_q_points_1d,
                 n_components,
                 Number,
                 VectorizedArrayType>::evaluate(const bool evaluate_values,
                                                const bool evaluate_gradients)
{
#  ifdef DEBUG
  Assert(this->dof_values_initialized, ExcNotInitialized());
#  endif

  evaluate(this->values_dofs[0], evaluate_values, evaluate_gradients);
}



template <int dim,
          int fe_degree,
          int n_q_points_1d,
          int n_components,
          typename Number,
          typename VectorizedArrayType>
inline void
FEFaceEvaluation<dim,
                 fe_degree,
                 n_q_points_1d,
                 n_components,
                 Number,
                 VectorizedArrayType>::
  evaluate(const EvaluationFlags::EvaluationFlags evaluation_flag)
{
#  ifdef DEBUG
  Assert(this->dof_values_initialized, ExcNotInitialized());
#  endif

  evaluate(this->values_dofs[0], evaluation_flag);
}



template <int dim,
          int fe_degree,
          int n_q_points_1d,
          int n_components,
          typename Number,
          typename VectorizedArrayType>
inline void
FEFaceEvaluation<dim,
                 fe_degree,
                 n_q_points_1d,
                 n_components,
                 Number,
                 VectorizedArrayType>::evaluate(const VectorizedArrayType
                                                  *        values_array,
                                                const bool evaluate_values,
                                                const bool evaluate_gradients)
{
  const EvaluationFlags::EvaluationFlags flag =
    ((evaluate_values) ? EvaluationFlags::values : EvaluationFlags::nothing) |
    ((evaluate_gradients) ? EvaluationFlags::gradients :
                            EvaluationFlags::nothing);

  evaluate(values_array, flag);
}



template <int dim,
          int fe_degree,
          int n_q_points_1d,
          int n_components,
          typename Number,
          typename VectorizedArrayType>
inline void
FEFaceEvaluation<dim,
                 fe_degree,
                 n_q_points_1d,
                 n_components,
                 Number,
                 VectorizedArrayType>::
  evaluate(const VectorizedArrayType *            values_array,
           const EvaluationFlags::EvaluationFlags evaluation_flag)
{
  Assert(
    (evaluation_flag &
     ~(EvaluationFlags::values | EvaluationFlags::gradients)) == 0,
    ExcMessage(
      "Only EvaluationFlags::values and EvaluationFlags::gradients are supported."));

  if (!(evaluation_flag & EvaluationFlags::values) &&
      !(evaluation_flag & EvaluationFlags::gradients))
    return;

  if (this->dof_access_index ==
        internal::MatrixFreeFunctions::DoFInfo::dof_access_cell &&
      this->is_interior_face == false)
    {
      const auto face_nos          = this->compute_face_no_data();
      const auto face_orientations = this->compute_face_orientations();

#  ifdef DEBUG
      // currently on structured meshes are supported -> face numbers and
      // orientations have to be the same for all filled lanes
      for (unsigned int v = 1; v < VectorizedArrayType::size(); ++v)
        {
          if (face_nos[v] != numbers::invalid_unsigned_int)
            AssertDimension(face_nos[0], face_nos[v]);
          if (face_orientations[v] != numbers::invalid_unsigned_int)
            AssertDimension(face_orientations[0], face_orientations[v]);
        }
#  endif

      internal::FEFaceEvaluationImplEvaluateSelector<dim, VectorizedArrayType>::
        template run<fe_degree, n_q_points_1d>(
          n_components,
          *this->data,
          values_array,
          this->begin_values(),
          this->begin_gradients(),
          this->scratch_data,
          evaluation_flag & EvaluationFlags::values,
          evaluation_flag & EvaluationFlags::gradients,
          face_nos[0],
          this->subface_index,
          face_orientations[0],
          this->descriptor->face_orientations);
    }
  else
    {
      if (fe_degree > -1)
        internal::FEFaceEvaluationImplEvaluateSelector<dim,
                                                       VectorizedArrayType>::
          template run<fe_degree, n_q_points_1d>(
            n_components,
            *this->data,
            values_array,
            this->begin_values(),
            this->begin_gradients(),
            this->scratch_data,
            evaluation_flag & EvaluationFlags::values,
            evaluation_flag & EvaluationFlags::gradients,
            this->face_no,
            this->subface_index,
            this->face_orientation,
            this->descriptor->face_orientations);
      else
        internal::FEFaceEvaluationFactory<dim, Number, VectorizedArrayType>::
          evaluate(n_components,
                   *this->data,
                   values_array,
                   this->begin_values(),
                   this->begin_gradients(),
                   this->scratch_data,
                   evaluation_flag & EvaluationFlags::values,
                   evaluation_flag & EvaluationFlags::gradients,
                   this->face_no,
                   this->subface_index,
                   this->face_orientation,
                   this->descriptor->face_orientations);
    }

#  ifdef DEBUG
  if (evaluation_flag & EvaluationFlags::values)
    this->values_quad_initialized = true;
  if (evaluation_flag & EvaluationFlags::gradients)
    this->gradients_quad_initialized = true;
#  endif
}



template <int dim,
          int fe_degree,
          int n_q_points_1d,
          int n_components,
          typename Number,
          typename VectorizedArrayType>
inline void
FEFaceEvaluation<dim,
                 fe_degree,
                 n_q_points_1d,
                 n_components,
                 Number,
                 VectorizedArrayType>::
  integrate(const EvaluationFlags::EvaluationFlags evaluation_flag)
{
  integrate(evaluation_flag, this->values_dofs[0]);

#  ifdef DEBUG
  this->dof_values_initialized = true;
#  endif
}



template <int dim,
          int fe_degree,
          int n_q_points_1d,
          int n_components,
          typename Number,
          typename VectorizedArrayType>
inline void
FEFaceEvaluation<dim,
                 fe_degree,
                 n_q_points_1d,
                 n_components,
                 Number,
                 VectorizedArrayType>::integrate(const bool integrate_values,
                                                 const bool integrate_gradients)
{
  integrate(integrate_values, integrate_gradients, this->values_dofs[0]);

#  ifdef DEBUG
  this->dof_values_initialized = true;
#  endif
}



template <int dim,
          int fe_degree,
          int n_q_points_1d,
          int n_components,
          typename Number,
          typename VectorizedArrayType>
inline void
FEFaceEvaluation<dim,
                 fe_degree,
                 n_q_points_1d,
                 n_components,
                 Number,
                 VectorizedArrayType>::integrate(const bool integrate_values,
                                                 const bool integrate_gradients,
                                                 VectorizedArrayType
                                                   *values_array)
{
  const EvaluationFlags::EvaluationFlags flag =
    ((integrate_values) ? EvaluationFlags::values : EvaluationFlags::nothing) |
    ((integrate_gradients) ? EvaluationFlags::gradients :
                             EvaluationFlags::nothing);

  integrate(flag, values_array);
}



template <int dim,
          int fe_degree,
          int n_q_points_1d,
          int n_components,
          typename Number,
          typename VectorizedArrayType>
inline void
FEFaceEvaluation<dim,
                 fe_degree,
                 n_q_points_1d,
                 n_components,
                 Number,
                 VectorizedArrayType>::
  integrate(const EvaluationFlags::EvaluationFlags evaluation_flag,
            VectorizedArrayType *                  values_array)
{
  Assert(
    (evaluation_flag &
     ~(EvaluationFlags::values | EvaluationFlags::gradients)) == 0,
    ExcMessage(
      "Only EvaluationFlags::values and EvaluationFlags::gradients are supported."));

  if (!(evaluation_flag & EvaluationFlags::values) &&
      !(evaluation_flag & EvaluationFlags::gradients))
    return;

  if (fe_degree > -1)
    internal::FEFaceEvaluationImplIntegrateSelector<dim, VectorizedArrayType>::
      template run<fe_degree, n_q_points_1d>(
        n_components,
        *this->data,
        values_array,
        this->begin_values(),
        this->begin_gradients(),
        this->scratch_data,
        evaluation_flag & EvaluationFlags::values,
        evaluation_flag & EvaluationFlags::gradients,
        this->face_no,
        this->subface_index,
        this->face_orientation,
        this->descriptor->face_orientations);
  else
    internal::FEFaceEvaluationFactory<dim, Number, VectorizedArrayType>::
      integrate(n_components,
                *this->data,
                values_array,
                this->begin_values(),
                this->begin_gradients(),
                this->scratch_data,
                evaluation_flag & EvaluationFlags::values,
                evaluation_flag & EvaluationFlags::gradients,
                this->face_no,
                this->subface_index,
                this->face_orientation,
                this->descriptor->face_orientations);
}



template <int dim,
          int fe_degree,
          int n_q_points_1d,
          int n_components_,
          typename Number,
          typename VectorizedArrayType>
template <typename VectorType>
inline void
FEFaceEvaluation<
  dim,
  fe_degree,
  n_q_points_1d,
  n_components_,
  Number,
  VectorizedArrayType>::gather_evaluate(const VectorType &input_vector,
                                        const bool        evaluate_values,
                                        const bool        evaluate_gradients)
{
  const EvaluationFlags::EvaluationFlags flag =
    ((evaluate_values) ? EvaluationFlags::values : EvaluationFlags::nothing) |
    ((evaluate_gradients) ? EvaluationFlags::gradients :
                            EvaluationFlags::nothing);

  gather_evaluate(input_vector, flag);
}



template <int dim,
          int fe_degree,
          int n_q_points_1d,
          int n_components_,
          typename Number,
          typename VectorizedArrayType>
template <typename VectorType>
inline void
FEFaceEvaluation<dim,
                 fe_degree,
                 n_q_points_1d,
                 n_components_,
                 Number,
                 VectorizedArrayType>::
  gather_evaluate(const VectorType &                     input_vector,
                  const EvaluationFlags::EvaluationFlags evaluation_flag)
{
  Assert(
    (evaluation_flag &
     ~(EvaluationFlags::values | EvaluationFlags::gradients)) == 0,
    ExcMessage(
      "Only EvaluationFlags::values and EvaluationFlags::gradients are supported."));

  const auto fu = [&]() {
    const auto shared_vector_data = internal::get_shared_vector_data(
      input_vector,
      this->dof_access_index ==
        internal::MatrixFreeFunctions::DoFInfo::dof_access_cell,
      this->active_fe_index,
      this->dof_info);

    if (this->dof_access_index ==
          internal::MatrixFreeFunctions::DoFInfo::dof_access_cell &&
        this->is_interior_face == false)
      {
        const auto cells             = this->get_cell_or_face_ids();
        const auto face_nos          = this->compute_face_no_data();
        const auto face_orientations = this->compute_face_orientations();

        return internal::FEFaceEvaluationImplGatherEvaluateSelector<
          dim,
          Number,
          VectorizedArrayType>::template run<fe_degree,
                                             n_q_points_1d>(
          n_components,
          VectorizedArrayType::size(),
          internal::get_beginning<Number>(input_vector),
          shared_vector_data,
          *this->data,
          *this->dof_info,
          this->begin_values(),
          this->begin_gradients(),
          this->scratch_data,
          evaluation_flag & EvaluationFlags::values,
          evaluation_flag & EvaluationFlags::gradients,
          this->active_fe_index,
          this->first_selected_component,
          cells,
          face_nos,
          this->subface_index,
          this->dof_access_index,
          face_orientations,
          this->descriptor->face_orientations);
      }
    else
      {
        // TODO: this copying should not be necessary once we have introduced
        // an internal-data structure
        std::array<unsigned int, VectorizedArrayType::size()> cells_   = {};
        std::array<unsigned int, VectorizedArrayType::size()> face_no_ = {};
        std::array<unsigned int, VectorizedArrayType::size()>
          face_orientation_ = {};

        cells_[0]            = this->cell;
        face_no_[0]          = this->face_no;
        face_orientation_[0] = this->face_orientation;

        if (fe_degree > -1)
          {
            return internal::FEFaceEvaluationImplGatherEvaluateSelector<
              dim,
              Number,
              VectorizedArrayType>::template run<fe_degree,
                                                 n_q_points_1d>(
              n_components,
              1,
              internal::get_beginning<Number>(input_vector),
              shared_vector_data,
              *this->data,
              *this->dof_info,
              this->begin_values(),
              this->begin_gradients(),
              this->scratch_data,
              evaluation_flag & EvaluationFlags::values,
              evaluation_flag & EvaluationFlags::gradients,
              this->active_fe_index,
              this->first_selected_component,
              cells_,
              face_no_,
              this->subface_index,
              this->dof_access_index,
              face_orientation_,
              this->descriptor->face_orientations);
          }
        else
          return internal::
            FEFaceEvaluationFactory<dim, Number, VectorizedArrayType>::
              gather_evaluate(n_components,
                              1,
                              internal::get_beginning<Number>(input_vector),
                              shared_vector_data,
                              *this->data,
                              *this->dof_info,
                              this->begin_values(),
                              this->begin_gradients(),
                              this->scratch_data,
                              evaluation_flag & EvaluationFlags::values,
                              evaluation_flag & EvaluationFlags::gradients,
                              this->active_fe_index,
                              this->first_selected_component,
                              cells_,
                              face_no_,
                              this->subface_index,
                              this->dof_access_index,
                              face_orientation_,
                              this->descriptor->face_orientations);
      }
  };

  if (!fu())
    {
      this->read_dof_values(input_vector);
      this->evaluate(evaluation_flag);
    }

#  ifdef DEBUG
  if (evaluation_flag & EvaluationFlags::values)
    this->values_quad_initialized = true;
  if (evaluation_flag & EvaluationFlags::gradients)
    this->gradients_quad_initialized = true;
#  endif
}



template <int dim,
          int fe_degree,
          int n_q_points_1d,
          int n_components_,
          typename Number,
          typename VectorizedArrayType>
template <typename VectorType>
inline void
FEFaceEvaluation<
  dim,
  fe_degree,
  n_q_points_1d,
  n_components_,
  Number,
  VectorizedArrayType>::integrate_scatter(const bool  integrate_values,
                                          const bool  integrate_gradients,
                                          VectorType &destination)
{
  const EvaluationFlags::EvaluationFlags flag =
    ((integrate_values) ? EvaluationFlags::values : EvaluationFlags::nothing) |
    ((integrate_gradients) ? EvaluationFlags::gradients :
                             EvaluationFlags::nothing);

  integrate_scatter(flag, destination);
}



template <int dim,
          int fe_degree,
          int n_q_points_1d,
          int n_components_,
          typename Number,
          typename VectorizedArrayType>
template <typename VectorType>
inline void
FEFaceEvaluation<dim,
                 fe_degree,
                 n_q_points_1d,
                 n_components_,
                 Number,
                 VectorizedArrayType>::
  integrate_scatter(const EvaluationFlags::EvaluationFlags evaluation_flag,
                    VectorType &                           destination)
{
  Assert((this->dof_access_index ==
            internal::MatrixFreeFunctions::DoFInfo::dof_access_cell &&
          this->is_interior_face == false) == false,
         ExcNotImplemented());

  const auto shared_vector_data = internal::get_shared_vector_data(
    destination,
    this->dof_access_index ==
      internal::MatrixFreeFunctions::DoFInfo::dof_access_cell,
    this->active_fe_index,
    this->dof_info);

  // TODO: this copying should not be necessary once we have introduced
  // an internal-data structure
  std::array<unsigned int, VectorizedArrayType::size()> cells_            = {};
  std::array<unsigned int, VectorizedArrayType::size()> face_no_          = {};
  std::array<unsigned int, VectorizedArrayType::size()> face_orientation_ = {};

  cells_[0]            = this->cell;
  face_no_[0]          = this->face_no;
  face_orientation_[0] = this->face_orientation;

  if (fe_degree > -1)
    {
      if (!internal::FEFaceEvaluationImplIntegrateScatterSelector<
            dim,
            Number,
            VectorizedArrayType>::template run<fe_degree,
                                               n_q_points_1d>(
            n_components,
            1,
            internal::get_beginning<Number>(destination),
            shared_vector_data,
            *this->data,
            *this->dof_info,
            this->begin_dof_values(),
            this->begin_values(),
            this->begin_gradients(),
            this->scratch_data,
            evaluation_flag & EvaluationFlags::values,
            evaluation_flag & EvaluationFlags::gradients,
            this->active_fe_index,
            this->first_selected_component,
            cells_,
            face_no_,
            this->subface_index,
            this->dof_access_index,
            face_orientation_,
            this->descriptor->face_orientations))
        {
          // if we arrive here, writing into the destination vector did not
          // succeed because some of the assumptions in integrate_scatter were
          // not fulfilled (e.g. an element or degree that does not support
          // direct writing), so we must do it here
          this->distribute_local_to_global(destination);
        }
    }
  else
    {
      if (!internal::FEFaceEvaluationFactory<dim, Number, VectorizedArrayType>::
            integrate_scatter(n_components,
                              1,
                              internal::get_beginning<Number>(destination),
                              shared_vector_data,
                              *this->data,
                              *this->dof_info,
                              this->begin_dof_values(),
                              this->begin_values(),
                              this->begin_gradients(),
                              this->scratch_data,
                              evaluation_flag & EvaluationFlags::values,
                              evaluation_flag & EvaluationFlags::gradients,
                              this->active_fe_index,
                              this->first_selected_component,
                              cells_,
                              face_no_,
                              this->subface_index,
                              this->dof_access_index,
                              face_orientation_,
                              this->descriptor->face_orientations))
        {
          this->distribute_local_to_global(destination);
        }
    }
}



template <int dim,
          int fe_degree,
          int n_q_points_1d,
          int n_components_,
          typename Number,
          typename VectorizedArrayType>
inline Point<dim, VectorizedArrayType>
FEFaceEvaluation<dim,
                 fe_degree,
                 n_q_points_1d,
                 n_components_,
                 Number,
                 VectorizedArrayType>::quadrature_point(const unsigned int q)
  const
{
  AssertIndexRange(q, n_q_points);
  if (this->dof_access_index < 2)
    {
      Assert(this->mapping_data->quadrature_point_offsets.empty() == false,
             internal::ExcMatrixFreeAccessToUninitializedMappingField(
               "update_quadrature_points"));
      AssertIndexRange(this->cell,
                       this->mapping_data->quadrature_point_offsets.size());
      return this->mapping_data->quadrature_points
        [this->mapping_data->quadrature_point_offsets[this->cell] + q];
    }
  else
    {
      Assert(this->matrix_info->get_mapping_info()
                 .face_data_by_cells[this->quad_no]
                 .quadrature_point_offsets.empty() == false,
             internal::ExcMatrixFreeAccessToUninitializedMappingField(
               "update_quadrature_points"));
      const unsigned int index =
        this->cell * GeometryInfo<dim>::faces_per_cell + this->face_no;
      AssertIndexRange(index,
                       this->matrix_info->get_mapping_info()
                         .face_data_by_cells[this->quad_no]
                         .quadrature_point_offsets.size());
      return this->matrix_info->get_mapping_info()
        .face_data_by_cells[this->quad_no]
        .quadrature_points[this->matrix_info->get_mapping_info()
                             .face_data_by_cells[this->quad_no]
                             .quadrature_point_offsets[index] +
                           q];
    }
}



template <int dim,
          int fe_degree,
          int n_q_points_1d,
          int n_components_,
          typename Number,
          typename VectorizedArrayType>
std::array<unsigned int, VectorizedArrayType::size()>
FEFaceEvaluation<dim,
                 fe_degree,
                 n_q_points_1d,
                 n_components_,
                 Number,
                 VectorizedArrayType>::compute_face_no_data()
{
  std::array<unsigned int, VectorizedArrayType::size()> face_no_data;

  if (this->dof_access_index !=
        internal::MatrixFreeFunctions::DoFInfo::dof_access_cell ||
      this->is_interior_face == true)
    {
      std::fill(face_no_data.begin(),
                face_no_data.begin() +
                  this->dof_info->n_vectorization_lanes_filled
                    [internal::MatrixFreeFunctions::DoFInfo::dof_access_cell]
                    [this->cell],
                this->face_no);
    }
  else
    {
      std::fill(face_no_data.begin(),
                face_no_data.end(),
                numbers::invalid_unsigned_int);

      for (unsigned int i = 0;
           i < this->dof_info->n_vectorization_lanes_filled
                 [internal::MatrixFreeFunctions::DoFInfo::dof_access_cell]
                 [this->cell];
           i++)
        {
          // compute actual (non vectorized) cell ID
          const unsigned int cell_this =
            this->cell * VectorizedArrayType::size() + i;
          // compute face ID
          const unsigned int face_index =
            this->matrix_info->get_cell_and_face_to_plain_faces()(this->cell,
                                                                  this->face_no,
                                                                  i);

          Assert(face_index != numbers::invalid_unsigned_int,
                 ExcNotInitialized());

          // get cell ID on both sides of face
          auto cell_m =
            this->matrix_info
              ->get_face_info(face_index / VectorizedArrayType::size())
              .cells_interior[face_index % VectorizedArrayType::size()];

          // compare the IDs with the given cell ID
          face_no_data[i] =
            (cell_m == cell_this) ?
              this->matrix_info
                ->get_face_info(face_index / VectorizedArrayType::size())
                .exterior_face_no :
              this->matrix_info
                ->get_face_info(face_index / VectorizedArrayType::size())
                .interior_face_no;
        }
    }

  return face_no_data;
}



template <int dim,
          int fe_degree,
          int n_q_points_1d,
          int n_components_,
          typename Number,
          typename VectorizedArrayType>
std::array<unsigned int, VectorizedArrayType::size()>
FEFaceEvaluation<dim,
                 fe_degree,
                 n_q_points_1d,
                 n_components_,
                 Number,
                 VectorizedArrayType>::compute_face_orientations()
{
  std::array<unsigned int, VectorizedArrayType::size()> face_no_data;

  if (this->dof_access_index !=
        internal::MatrixFreeFunctions::DoFInfo::dof_access_cell ||
      this->is_interior_face == true)
    {
      Assert(false, ExcNotImplemented());
    }
  else
    {
      std::fill(face_no_data.begin(),
                face_no_data.end(),
                numbers::invalid_unsigned_int);

      if (dim == 3)
        {
          for (unsigned int i = 0;
               i < this->dof_info->n_vectorization_lanes_filled
                     [internal::MatrixFreeFunctions::DoFInfo::dof_access_cell]
                     [this->cell];
               i++)
            {
              // compute actual (non vectorized) cell ID
              const unsigned int cell_this =
                this->cell * VectorizedArrayType::size() + i;
              // compute face ID
              const unsigned int face_index =
                this->matrix_info->get_cell_and_face_to_plain_faces()(
                  this->cell, this->face_no, i);

              Assert(face_index != numbers::invalid_unsigned_int,
                     ExcNotInitialized());

              const unsigned int macro =
                face_index / VectorizedArrayType::size();
              const unsigned int lane =
                face_index % VectorizedArrayType::size();

              const auto &faces = this->matrix_info->get_face_info(macro);

              // get cell ID on both sides of face
              auto cell_m = faces.cells_interior[lane];

              const bool is_interior_face = cell_m != cell_this;
              const bool fo_interior_face = faces.face_orientation >= 8;

              unsigned int face_orientation = faces.face_orientation % 8;

              if (is_interior_face != fo_interior_face)
                {
                  // invert (see also:
                  // Triangulation::update_periodic_face_map())
                  static const std::array<unsigned int, 8> table{
                    {0, 1, 0, 3, 6, 5, 4, 7}};

                  face_orientation = table[face_orientation];
                }

              // compare the IDs with the given cell ID
              face_no_data[i] = face_orientation;
            }
        }
      else
        {
          std::fill(
            face_no_data.begin(),
            face_no_data.begin() +
              this->dof_info->n_vectorization_lanes_filled
                [internal::MatrixFreeFunctions::DoFInfo::dof_access_cell]
                [this->cell],
            0);
        }
    }

  return face_no_data;
}



template <int dim,
          int fe_degree,
          int n_q_points_1d,
          int n_components_,
          typename Number,
          typename VectorizedArrayType>
bool
FEEvaluation<dim,
             fe_degree,
             n_q_points_1d,
             n_components_,
             Number,
             VectorizedArrayType>::
  fast_evaluation_supported(const unsigned int given_degree,
                            const unsigned int give_n_q_points_1d)
{
  return fe_degree == -1 ?
           internal::FEEvaluationFactory<dim, Number, VectorizedArrayType>::
             fast_evaluation_supported(given_degree, give_n_q_points_1d) :
           true;
}



template <int dim,
          int fe_degree,
          int n_q_points_1d,
          int n_components_,
          typename Number,
          typename VectorizedArrayType>
bool
FEFaceEvaluation<dim,
                 fe_degree,
                 n_q_points_1d,
                 n_components_,
                 Number,
                 VectorizedArrayType>::
  fast_evaluation_supported(const unsigned int given_degree,
                            const unsigned int give_n_q_points_1d)
{
  return fe_degree == -1 ?
           internal::FEFaceEvaluationFactory<dim, Number, VectorizedArrayType>::
             fast_evaluation_supported(given_degree, give_n_q_points_1d) :
           true;
}



 /*------------------------- end FEFaceEvaluation ------------------------- */ 


#endif // ifndef DOXYGEN


DEAL_II_NAMESPACE_CLOSE

#endif


