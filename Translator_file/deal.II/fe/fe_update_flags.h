//include/deal.II-translator/fe/fe_update_flags_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2020 by the deal.II authors
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

#ifndef dealii_fe_update_flags_h
#define dealii_fe_update_flags_h


#include <deal.II/base/config.h>

#include <deal.II/base/derivative_form.h>
#include <deal.II/base/point.h>
#include <deal.II/base/table.h>
#include <deal.II/base/tensor.h>

#include <vector>


DEAL_II_NAMESPACE_OPEN

// Forward declaration
#ifndef DOXYGEN
template <int, int>
class FiniteElement;
#endif

 /*!@addtogroup feaccess */ 
 /*@{*/ 

/**
 * 给予FEValues、FEFaceValues和FESubfaceValues的构造函数的枚举类型，告诉这些对象在每个网格单元上需要哪些数据。
 * 以限制性的方式选择这些标志对于 FEValues::reinit(),
 * FEFaceValues::reinit() 和 FESubfaceValues::reinit().
 * 的效率至关重要。因此，应该只选择实际需要的标志。相关的Mapping和FiniteElement有责任根据自己的要求添加额外的标志。例如，如果选择了#update_gradients，大多数有限元会添加#update_covariant_transformation。
 * 默认情况下，所有标志都是关闭的，也就是说，不会进行重新初始化。
 * 你可以用位数或运算符|(UpdateFlags,UpdateFlags)通过串联的方式选择多个标志。
 * <h3>Use of these flags flags</h3>
 *
 *
 */
enum UpdateFlags
{
  //! No update
  update_default = 0,
  //! Shape function values
  /**
   * 计算实空间单元上正交点的形状函数值。对于通常的拉格朗日元素，这些值等于单元格上正交点的形状函数值，但对于更复杂的元素，如FE_RaviartThomas元素，它们是不同的。
   *
   */
  update_values = 0x0001,
  //! Shape function gradients
  /**
   * 计算实单元坐标中形状函数的梯度。
   *
   */
  update_gradients = 0x0002,
  //! Second derivatives of shape functions
  /**
   * 计算实心单元坐标下的形状函数的二阶导数。
   *
   */
  update_hessians = 0x0004,
  //! Third derivatives of shape functions
  /**
   * 计算实心单元坐标中的形状函数的第三导数
   *
   */
  update_3rd_derivatives = 0x0008,
  //! Outer normal vector, not normalized
  /**
   * 切向矢量的矢量乘积，产生一个法向矢量，其长度对应于表面元素；可能比计算两者更有效。
   *
   */
  update_boundary_forms = 0x0010,
  //! Transformed quadrature points
  /**
   * 计算实际单元坐标中的正交点位置。
   * FEValues对象将参考单元上的正交点位置作为构造函数的一个参数（通过Quadrature对象）。对于大多数有限元来说，知道参考单元上正交点的位置是评估形状函数、评估映射和其他事情所必需的。另一方面，如果你想在真实单元上的正交点位置
   * $\mathbf x_q$ 评估一个右手函数 $f(\mathbf x_q)$
   * ，你需要将这个标志传递给FEValues构造函数，以确保你以后可以访问它们。
   * 在DataPostprocessor的背景下，
   * DataPostprocessorInputs::CommonInputs::evaluation_points 将被更新。
   *
   */
  update_quadrature_points = 0x0020,
  //! Transformed quadrature weights
  /**
   * 计算实数单元上的正交权重，即正交规则的权重乘以从参考单元到实数单元的转换的雅各布定理。
   *
   */
  update_JxW_values = 0x0040,
  //! Normal vectors
  /**
   * 计算法向量，可以是一个面，也可以是一个一维的单元。对任何其他对象设置这个标志都会引起一个错误。
   *
   */
  update_normal_vectors = 0x0080,
  //! Volume element
  /**
   * 计算从参考单元到实际单元的转换的雅各布系数。
   *
   */
  update_jacobians = 0x0100,
  //! Gradient of volume element
  /**
   * 计算变换的Jacobian的导数。
   *
   */
  update_jacobian_grads = 0x0200,
  //! Volume element
  /**
   * 计算从参考单元到实数单元的变换的逆雅各布系数。
   *
   */
  update_inverse_jacobians = 0x0400,
  //! Covariant transformation
  /**
   * 计算Mapping对向量进行反演变换所需的所有数值。对于像MappingCartesian这样的特殊映射，这可能比#update_inverse_jacobians更简单。
   *
   */
  update_covariant_transformation = 0x0800,
  //! Contravariant transformation
  /**
   * 计算Mapping需要的所有值，以便对向量进行逆向变换。对于像MappingCartesian这样的特殊映射，这可能比#update_jacobians更简单。
   *
   */
  update_contravariant_transformation = 0x1000,
  //! Shape function values of transformation
  /**
   * 计算由Mapping定义的变换的形状函数值。
   *
   */
  update_transformation_values = 0x2000,
  //! Shape function gradients of transformation
  /**
   * 计算由Mapping定义的变换的形状函数梯度。
   *
   */
  update_transformation_gradients = 0x4000,
  //! Determinant of the Jacobian
  /**
   * 计算每个正交点的体积元素。
   *
   */
  update_volume_elements = 0x10000,
  /**
   * 计算向前推到实际单元坐标的变换的雅各布系数的导数。
   *
   */
  update_jacobian_pushed_forward_grads = 0x100000,
  /**
   * 计算变换的Jacobian的二次导数。
   *
   */
  update_jacobian_2nd_derivatives = 0x200000,
  /**
   * 计算向前推至真实单元坐标的变换的雅各布系数的二阶导数。
   *
   */
  update_jacobian_pushed_forward_2nd_derivatives = 0x400000,
  /**
   * 计算变换的雅各布系数的第三导数。
   *
   */
  update_jacobian_3rd_derivatives = 0x800000,
  /**
   * 计算向前推至实数单元坐标的变换的雅各布系数的第三导数。
   *
   */
  update_jacobian_pushed_forward_3rd_derivatives = 0x1000000,
  //! Values needed for Piola transform
  /**
   * 结合Hdiv元素的Piola变换所需的标志。
   *
   */
  update_piola = update_volume_elements | update_contravariant_transformation,
  /**
   * 需要进行映射计算的标志的组合
   *
   */
  update_mapping =
    // Direct data
  update_quadrature_points | update_JxW_values | update_jacobians |
  update_jacobian_grads | update_jacobian_pushed_forward_grads |
  update_jacobian_2nd_derivatives |
  update_jacobian_pushed_forward_2nd_derivatives |
  update_jacobian_3rd_derivatives |
  update_jacobian_pushed_forward_3rd_derivatives | update_inverse_jacobians |
  update_boundary_forms | update_normal_vectors |
  // Transformation dependence
  update_covariant_transformation | update_contravariant_transformation |
  update_transformation_values | update_transformation_gradients |
  // Volume data
  update_volume_elements
};


/**
 * 输出操作符，它将更新标志作为一组or'd文本值输出。
 * @ref UpdateFlags
 *
 *
 */
template <class StreamType>
inline StreamType &
operator<<(StreamType &s, const UpdateFlags u)
{
  s << " UpdateFlags|";
  if (u & update_values)
    s << "values|";
  if (u & update_gradients)
    s << "gradients|";
  if (u & update_hessians)
    s << "hessians|";
  if (u & update_3rd_derivatives)
    s << "3rd_derivatives|";
  if (u & update_quadrature_points)
    s << "quadrature_points|";
  if (u & update_JxW_values)
    s << "JxW_values|";
  if (u & update_normal_vectors)
    s << "normal_vectors|";
  if (u & update_jacobians)
    s << "jacobians|";
  if (u & update_inverse_jacobians)
    s << "inverse_jacobians|";
  if (u & update_jacobian_grads)
    s << "jacobian_grads|";
  if (u & update_covariant_transformation)
    s << "covariant_transformation|";
  if (u & update_contravariant_transformation)
    s << "contravariant_transformation|";
  if (u & update_transformation_values)
    s << "transformation_values|";
  if (u & update_transformation_gradients)
    s << "transformation_gradients|";
  if (u & update_jacobian_pushed_forward_grads)
    s << "jacobian_pushed_forward_grads|";
  if (u & update_jacobian_2nd_derivatives)
    s << "jacobian_2nd_derivatives|";
  if (u & update_jacobian_pushed_forward_2nd_derivatives)
    s << "jacobian_pushed_forward_2nd_derivatives|";
  if (u & update_jacobian_3rd_derivatives)
    s << "jacobian_3rd_derivatives|";
  if (u & update_jacobian_pushed_forward_3rd_derivatives)
    s << "jacobian_pushed_forward_3rd_derivatives|";

  // TODO: check that 'u' really only has the flags set that are handled above
  return s;
}


/**
 * 全局操作符，它返回一个对象，其中所有的位都被设置为第一或第二个参数中的设置。这个操作符的存在是因为如果它不存在，那么bit-or
 * <tt>操作符|</tt>的结果将是一个整数，当我们试图将其分配给UpdateFlags类型的对象时，又会引发编译器警告。
 * @ref UpdateFlags
 *
 *
 */
inline UpdateFlags
operator|(const UpdateFlags f1, const UpdateFlags f2)
{
  return static_cast<UpdateFlags>(static_cast<unsigned int>(f1) |
                                  static_cast<unsigned int>(f2));
}



/**
 * 全局操作符，它将第二个参数的位也设置在第一个参数中。
 * @ref UpdateFlags
 *
 *
 */
inline UpdateFlags &
operator|=(UpdateFlags &f1, const UpdateFlags f2)
{
  f1 = f1 | f2;
  return f1;
}


/**
 * 全局操作符，它返回一个对象，其中所有位都被设置在第一个和第二个参数中。这个操作符的存在是因为如果它不存在，那么位和<tt>操作符&</tt>的结果将是一个整数，当我们试图将其分配给UpdateFlags类型的对象时，会引发编译器警告。
 * @ref UpdateFlags
 *
 *
 */
inline UpdateFlags operator&(const UpdateFlags f1, const UpdateFlags f2)
{
  return static_cast<UpdateFlags>(static_cast<unsigned int>(f1) &
                                  static_cast<unsigned int>(f2));
}


/**
 * 全局操作符，如果第一个参数中的所有位没有在第二个参数中设置，则将其清除。
 * @ref UpdateFlags
 *
 *
 */
inline UpdateFlags &
operator&=(UpdateFlags &f1, const UpdateFlags f2)
{
  f1 = f1 & f2;
  return f1;
}



/**
 * 这个枚举定义用于存储当前单元与先前访问的单元的相似性。这些信息用于在调用方法
 * FEValues::reinit()
 * 时重复使用数据（比如导数，如果一个单元格只是前一个单元格的翻译，则导数不会改变）。目前，这个变量只能识别一个平移和一个倒置的平移（如果dim<spacedim）。然而，这个概念使得在FEValues/FEFaceValues中添加额外的状态来检测这些相似性也变得容易。
 *
 *
 */
namespace CellSimilarity
{
  enum Similarity
  {
    /**
     * 除了平移或倒置的平移之外，单元格还存在一些差异。
     *
     */
    none,
    /**
     * 这些单元格因翻译而不同。
     *
     */
    translation,
    /**
     * 这些单元格因倒置的平移而不同。
     *
     */
    inverted_translation,
    /**
     * 下一个单元格是无效的。
     *
     */
    invalid_next_cell
  };
}


namespace internal
{
  namespace FEValuesImplementation
  {
    /**
     * 一个存储所有用于 dealii::FEValues,  dealii::FEFaceValues, 和
     * dealii::FESubfaceValues 对象中的映射相关数据的类。当
     * dealii::FEValues::reinit() 调用 Mapping::fill_fe_values()
     * 时，这种对象将作为<i>output</i>的参数提供给某个单元、面或子面。
     * 然后这里的数据将作为<i>input</i>参数在下面调用
     * FiniteElement::fill_fe_values(). 时提供。
     * @ingroup feaccess
     *
     */
    template <int dim, int spacedim = dim>
    class MappingRelatedData
    {
    public:
      /**
       * 将所有向量初始化为正确的大小。
       *
       */
      void
      initialize(const unsigned int n_quadrature_points,
                 const UpdateFlags  flags);

      /**
       * 计算并返回这个对象的内存消耗（以字节为单位）的估计值。
       *
       */
      std::size_t
      memory_consumption() const;

      /**
       * 在正交点存储一个权重乘以雅可比行列式的数组。每次调用reinit()时，这个函数都会被重置。雅可比行列式实际上是存储在这个类中的雅可比矩阵的倒数值，更多信息见这个类的一般文档。
       * 然而，如果这个对象指的是FEFaceValues或FESubfaceValues对象，那么JxW_values对应的是面的变换的雅可比，而不是单元，也就是说，维度是面的量度，而不是体积的量度。在这种情况下，它是由边界形式，而不是由雅各布矩阵计算出来的。
       *
       */
      std::vector<double> JxW_values;

      /**
       * 正交点的雅各布矩阵的阵列。
       *
       */
      std::vector<DerivativeForm<1, dim, spacedim>> jacobians;

      /**
       * 雅各布矩阵在正交点的导数数组。
       *
       */
      std::vector<DerivativeForm<2, dim, spacedim>> jacobian_grads;

      /**
       * 正交点上的反雅各布矩阵阵列。
       *
       */
      std::vector<DerivativeForm<1, spacedim, dim>> inverse_jacobians;

      /**
       * 正交点的雅各布矩阵的导数数组，向前推至实际单元坐标。
       *
       */
      std::vector<Tensor<3, spacedim>> jacobian_pushed_forward_grads;

      /**
       * 正交点的雅各布矩阵的二阶导数数组。
       *
       */
      std::vector<DerivativeForm<3, dim, spacedim>> jacobian_2nd_derivatives;

      /**
       * 在正交点的雅各布矩阵的二阶导数数组，向前推至真实单元坐标。
       *
       */
      std::vector<Tensor<4, spacedim>> jacobian_pushed_forward_2nd_derivatives;

      /**
       * 正交点上的雅各布矩阵的第三导数数组。
       *
       */
      std::vector<DerivativeForm<4, dim, spacedim>> jacobian_3rd_derivatives;

      /**
       * 在正交点的雅各布矩阵的三次导数数组，向前推至真实单元坐标。
       *
       */
      std::vector<Tensor<5, spacedim>> jacobian_pushed_forward_3rd_derivatives;

      /**
       * 正交点的数组。这个数组是在调用reinit()时设置的，包含实数元素上的正交点，而不是参考元素上的。
       *
       */
      std::vector<Point<spacedim>> quadrature_points;

      /**
       * 正交点的外向法向量列表。
       *
       */
      std::vector<Tensor<1, spacedim>> normal_vectors;

      /**
       * 正交点上的边界形式列表。
       *
       */
      std::vector<Tensor<1, spacedim>> boundary_forms;
    };


    /**
     * 一个存储所有用于 dealii::FEValues,  dealii::FEFaceValues, 和
     * dealii::FESubfaceValues
     * 对象中的形状函数相关数据的类。当
     * dealii::FEValues::reinit() 调用 FiniteElement::fill_fe_values().
     * 时，这种对象将作为<i>output</i>参数给出。
     * @ingroup feaccess
     *
     */
    template <int dim, int spacedim = dim>
    class FiniteElementRelatedData
    {
    public:
      /**
       * 将所有向量初始化为正确的大小。
       *
       */
      void
      initialize(const unsigned int                  n_quadrature_points,
                 const FiniteElement<dim, spacedim> &fe,
                 const UpdateFlags                   flags);

      /**
       * 计算并返回这个对象的内存消耗（以字节为单位）的估计值。
       *
       */
      std::size_t
      memory_consumption() const;

      /**
       * 形状值的存储类型。矩阵中的每一行表示不同点上的单个形状函数的值，列是针对单个点的不同形状函数。
       * 如果一个形状函数有一个以上的非零分量（在deal.II的字典中：它是非正数），那么我们为每个非零分量分配一行，并将随后的行向后移。
       * 因此，如果整个有限元是原始的（即所有形状函数都是原始的），为一个形状函数查询正确的行是很简单的，因为那时形状函数的编号等于行的编号。否则，使用#shape_function_to_row_table数组来获得属于这个特定形状函数的第一行，并使用
       * FiniteElement::get_nonzero_components()
       * 函数在这个形状函数的所有行中导航，该函数告诉我们哪些组件是非零的，因此在目前讨论的数组中有一行。
       *
       */
      using ShapeVector = dealii::Table<2, double>;

      /**
       * 梯度的存储类型。数据的布局与#ShapeVector数据类型相同。
       *
       */
      using GradientVector = dealii::Table<2, Tensor<1, spacedim>>;

      /**
       * 同样，对于二阶导数也是如此。
       *
       */
      using HessianVector = dealii::Table<2, Tensor<2, spacedim>>;

      /**
       * 而同样的情况也适用于三阶导数。
       *
       */
      using ThirdDerivativeVector = dealii::Table<2, Tensor<3, spacedim>>;

      /**
       * 存储正交点的形状函数的值。关于这个字段中的数据布局，请参见数据类型的描述。
       *
       */
      ShapeVector shape_values;

      /**
       * 存储形状函数在正交点的梯度。
       * 参见数据类型的描述，了解该字段的数据布局。
       *
       */
      GradientVector shape_gradients;

      /**
       * 存储形状函数在正交点的二阶导数。
       * 参见数据类型的描述，了解该字段的数据布局。
       *
       */
      HessianVector shape_hessians;

      /**
       * 存储正交点上的形状函数的3次导数。
       * 参见数据类型的描述，了解该字段的数据布局。
       *
       */
      ThirdDerivativeVector shape_3rd_derivatives;

      /**
       * 当被问及形状函数i的第c个向量分量的值（或梯度，或Hessian）时，我们需要在#shape_values、#shape_gradients和#shape_hessians数组中查找它。
       * 问题是，在这个数组中，形状函数i的c分量的数据在哪里？这就是这个表的答案。
       * 该表的格式如下。
       *
       * - 它有dofs_per_cell乘以n_components条目。
       *
       * - 对应于形状函数i，组件c的条目是 <code>i n_components + c</code>  。
       *
       * - 存储在这个位置的值表示#shape_values和其他表格中所有正交点的对应基准点的行。            在一般情况下，矢量值的情况下，分量的数量大于1，但是对于一个给定的形状函数，并非所有的矢量分量都是非零的（例如，如果一个形状函数是原始的，那么正好有一个矢量分量是非零的，而其他的都是零）。对于这种零分量，#shape_values和friends没有一行。因此，对于形状函数i为零的向量分量，在当前表中的条目是 numbers::invalid_unsigned_int.  另一方面，该表保证每个形状函数至少有一个有效的索引。特别是，对于原始有限元，每个形状函数正好有一个非零分量，因此对于每个i，在 <code>[i*n_components, (i+1)*n_components)</code> 范围内正好有一个有效索引。
       *
       */
      std::vector<unsigned int> shape_function_to_row_table;
    };
  } // namespace FEValuesImplementation
} // namespace internal


 /*@}*/ 



DEAL_II_NAMESPACE_CLOSE

#endif


