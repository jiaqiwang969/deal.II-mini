//include/deal.II-translator/fe/fe_values_0.txt
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

#ifndef dealii_fe_values_h
#define dealii_fe_values_h


#include <deal.II/base/config.h>

#include <deal.II/base/derivative_form.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/point.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/std_cxx20/iota_view.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/base/symmetric_tensor.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_update_flags.h>
#include <deal.II/fe/fe_values_extractors.h>
#include <deal.II/fe/mapping.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/hp/q_collection.h>

#include <algorithm>
#include <memory>
#include <type_traits>


// dummy include in order to have the
// definition of PetscScalar available
// without including other PETSc stuff
#ifdef DEAL_II_WITH_PETSC
#  include <petsc.h>
#endif

DEAL_II_NAMESPACE_OPEN

// Forward declaration
#ifndef DOXYGEN
template <int dim, int spacedim = dim>
class FEValuesBase;
#endif

namespace internal
{
  /**
   * 一个类，其特化用于定义一个向量值函数的卷曲对应于什么类型。
   *
   */
  template <int dim, class NumberType = double>
  struct CurlType;

  /**
   * 一个专门用于定义向量值函数的curl所对应的类型的类。
   * 在1d中，curl是一个标量。
   *
   */
  template <class NumberType>
  struct CurlType<1, NumberType>
  {
    using type = Tensor<1, 1, NumberType>;
  };

  /**
   * 一个专门用于定义向量值函数的curl所对应的类型的类。
   * 在2d中，curl是一个标量。
   *
   */
  template <class NumberType>
  struct CurlType<2, NumberType>
  {
    using type = Tensor<1, 1, NumberType>;
  };

  /**
   * 一个专门用于定义向量值函数的curl所对应的类型的类。
   * 在3D中，curl是一个矢量。
   *
   */
  template <class NumberType>
  struct CurlType<3, NumberType>
  {
    using type = Tensor<1, 3, NumberType>;
  };
} // namespace internal



/**
 * 一个用于FEValues、FEFaceValues或FESSubfaceValues对象的 "视图
 * "的命名空间。一个视图只代表整体的某一部分：而FEValues对象代表<i>all</i>向量值元素的所有分量的值、梯度或二阶导数，视图将注意力限制在单个分量或分量的一个子集上。你通常通过使用方括号操作符将FEValuesExtractors对象应用于FEValues、FEFaceValues或FESubfaceValues对象，来获得这个命名空间中定义的类的对象。
 * 有一些类为单个标量组件、由 <code>dim</code>
 * 元素组成的向量组件和由 <code>(dim*dim + dim)/2</code>
 * 元素组成的对称二阶张量组件提供视图。
 * 参见 @ref
 * vector_valued
 * 模块的描述，了解如何使用该命名空间的特征的例子。
 *
 *
 * @ingroup feaccess vector_valued
 *
 *
 */
namespace FEValuesViews
{
  /**
   * 一个代表对一个可能是矢量值的有限元的单一标量分量的视图的类。视图将在 @ref
   * vector_valued 模块中讨论。    如果你将
   * FEValuesExtractors::Scalar
   * 应用于FEValues、FEFaceValues或FESSubfaceValues对象，你会得到这种类型的对象。
   * @ingroup feaccess vector_valued
   *
   */
  template <int dim, int spacedim = dim>
  class Scalar
  {
  public:
    /**
     * 这个类所代表的视图值的数据类型的别名。由于我们处理的是一个单一的组件，所以值的类型是一个标量的双数。
     *
     */
    using value_type = double;

    /**
     * 该类代表的视图的梯度类型的别名。
     * 这里，对于有限元的标量分量，梯度是一个
     * <code>Tensor@<1,dim@></code>  。
     *
     */
    using gradient_type = dealii::Tensor<1, spacedim>;

    /**
     * 该类代表的视图的二阶导数类型的别名。这里，对于有限元的标量分量，Hessian是一个
     * <code>Tensor@<2,dim@></code>  .
     *
     */
    using hessian_type = dealii::Tensor<2, spacedim>;

    /**
     * 该类所代表的视图的第三导数的类型的别名。这里，对于有限元的标量分量，第三导数是一个
     * <code>Tensor@<3,dim@></code>  。
     *
     */
    using third_derivative_type = dealii::Tensor<3, spacedim>;

    /**
     * a  @p Number
     * 与本类提供的视图值的乘积的数据类型的别名。这是一个有限元场的标量分量的数据类型，其自由度由元素类型为
     * @p Number. 的向量描述。
     *
     */
    template <typename Number>
    using solution_value_type = typename ProductType<Number, value_type>::type;

    /**
     * @p Number
     * 和该类提供的视图梯度的乘积的数据类型的别名。这是一个有限元场的标量分量的数据类型，其自由度由元素类型为
     * @p Number. 的向量描述。
     *
     */
    template <typename Number>
    using solution_gradient_type =
      typename ProductType<Number, gradient_type>::type;

    /**
     * 一个 @p Number
     * 和该类提供的视图的拉普拉斯的乘积的数据类型的别名。这是一个有限元场的标量分量的数据类型，其自由度由元素类型为
     * @p Number. 的向量描述。
     *
     */
    template <typename Number>
    using solution_laplacian_type =
      typename ProductType<Number, value_type>::type;

    /**
     * 一个 @p Number
     * 和该类提供的视图的赫西的乘积的数据类型的别名。这是一个有限元场的标量分量的数据类型，其自由度由元素类型为
     * @p Number. 的向量描述。
     *
     */
    template <typename Number>
    using solution_hessian_type =
      typename ProductType<Number, hessian_type>::type;

    /**
     * 一个 @p Number
     * 和该类提供的视图的第三导数的乘积的数据类型的别名。这是一个有限元场的标量分量的数据类型，其自由度由元素类型为
     * @p Number. 的向量描述。
     *
     */
    template <typename Number>
    using solution_third_derivative_type =
      typename ProductType<Number, third_derivative_type>::type;

    /**
     * 一个结构，为标量视图和任何 @p Number
     * 类型的基函数的值和导数的乘积提供输出类型。
     * @deprecated  使用周围类中定义的类型来代替。
     *
     */
    template <typename Number>
    struct DEAL_II_DEPRECATED OutputType
    {
      /**
       * 一个 @p Number
       * 和Scalar类的视图值的乘积的数据类型的别名。
       *
       */
      using value_type =
        typename ProductType<Number,
                             typename Scalar<dim, spacedim>::value_type>::type;

      /**
       * @p Number
       * 和Scalar类视图的梯度的乘积的数据类型的别名。
       *
       */
      using gradient_type = typename ProductType<
        Number,
        typename Scalar<dim, spacedim>::gradient_type>::type;

      /**
       * 一个 @p Number
       * 和拉普拉斯类视图的乘积的数据类型的别名。
       *
       */
      using laplacian_type =
        typename ProductType<Number,
                             typename Scalar<dim, spacedim>::value_type>::type;

      /**
       * @p Number
       * 和Scalar类视图的Hessians之积的数据类型的别名。
       *
       */
      using hessian_type = typename ProductType<
        Number,
        typename Scalar<dim, spacedim>::hessian_type>::type;

      /**
       * 一个 @p Number
       * 的乘积和Scalar类视图的第三导数的数据类型的别名。
       *
       */
      using third_derivative_type = typename ProductType<
        Number,
        typename Scalar<dim, spacedim>::third_derivative_type>::type;
    };

    /**
     * 一个结构，对于每个形状函数，我们预先计算出一堆数据，这将使以后的访问变得更加便宜。
     *
     */
    struct ShapeFunctionData
    {
      /**
       * 对于每个形状函数，存储所选向量分量是否可能为非零。对于原始形状函数，我们可以肯定地知道某个给定形状函数的某个标量分量是否为非零，而对于非原始形状函数，这可能并不完全清楚（例如，对于RT元素，它取决于单元格的形状）。
       *
       */
      bool is_nonzero_shape_function_component;

      /**
       * 对于每个形状函数，在shape_values、shape_gradients和shape_hessians表中存储行索引（列索引是正交点的索引）。如果形状函数是原始的，那么我们可以从FEValues对象的shape_function_to_row_table中获得这些信息；否则，我们必须花点功夫来计算这些信息。
       *
       */
      unsigned int row_index;
    };

    /**
     * 默认构造函数。创建一个无效的对象。
     *
     */
    Scalar();

    /**
     * 表示FEValuesBase对象（或从FEValuesBase派生的类之一）的单个标量分量的对象的构造函数。
     *
     */
    Scalar(const FEValuesBase<dim, spacedim> &fe_values_base,
           const unsigned int                 component);

    /**
     * 复制构造函数。这不是一个轻量级的对象，所以我们不允许复制，如果调用这个函数会产生一个编译时错误。
     *
     */
    Scalar(const Scalar<dim, spacedim> &) = delete;

    /**
     * 移动构造函数。
     *
     */
    // NOLINTNEXTLINE OSX does not compile with noexcept
    Scalar(Scalar<dim, spacedim> &&) = default;

    /**
     * 解构器。
     *
     */
    ~Scalar() = default;

    /**
     * 拷贝操作符。这不是一个轻量级的对象，所以我们不允许复制，如果调用这个函数会产生一个编译时错误。
     *
     */
    Scalar &
    operator=(const Scalar<dim, spacedim> &) = delete;

    /**
     * 移动赋值运算符。
     *
     */
    Scalar &
    operator=(Scalar<dim, spacedim> &&) noexcept = default;

    /**
     * 返回该视图所选择的矢量分量的值，用于参数所选择的形状函数和正交点。
     * @param  shape_function 要评估的形状函数的编号。
     * 请注意，这个数字从零到dofs_per_cell，即使是在FEFaceValues或FESubfaceValues对象的情况下。
     * @param  q_point 要评估函数的正交点的编号。
     * @dealiiRequiresUpdateFlags{update_values}
     *
     */
    value_type
    value(const unsigned int shape_function, const unsigned int q_point) const;

    /**
     * 返回该视图所选择的矢量分量的梯度（等级为1的张量），用于形状函数和参数选择的正交点。
     * @note 参数的含义与value()函数的记录相同。
     * @dealiiRequiresUpdateFlags{update_gradients}
     *
     */
    gradient_type
    gradient(const unsigned int shape_function,
             const unsigned int q_point) const;

    /**
     * 对于参数选择的形状函数和正交点，返回该视图选择的向量分量的Hessian（所有二次导数的等级2的张量）。
     * @note  参数的含义与value()函数的记载相同。
     * @dealiiRequiresUpdateFlags{update_hessians}
     *
     */
    hessian_type
    hessian(const unsigned int shape_function,
            const unsigned int q_point) const;

    /**
     * 返回该视图选择的向量分量的所有三次导数的等级3的张量，用于参数选择的形状函数和正交点。
     * @note 参数的含义与value()函数的记录相同。
     * @dealiiRequiresUpdateFlags{update_third_derivatives}
     *
     */
    third_derivative_type
    third_derivative(const unsigned int shape_function,
                     const unsigned int q_point) const;

    /**
     * 返回上次调用FEValues对象的<tt>reinit</tt>函数时选择的单元、面或子面的正交点上由<tt>fe_function</tt>表征的有限元函数的标量分量的值。
     * 这个函数等同于 FEValuesBase::get_function_values
     * 函数，但它只对选定的标量分量起作用。
     * 输出向量所存储的数据类型必须是你将形状函数的值相乘后得到的（即
     * @p value_type) 乘以用于存储你的有限元向量 $U$ 的未知数
     * $U_j$ 值的类型（由 @p fe_function 参数代表）。
     * @dealiiRequiresUpdateFlags{update_values}
     *
     */
    template <class InputVector>
    void
    get_function_values(
      const InputVector &fe_function,
      std::vector<solution_value_type<typename InputVector::value_type>>
        &values) const;

    /**
     * 与上述相同，但使用局部自由度值的向量。换句话说，这个函数不是从与DoFHandler对象相关的全局向量中提取位于当前单元上的自由度的节点值（如上面的函数），而是通过其第一个参数获取这些局部节点值。获得这样一个向量的典型方法是通过调用如下代码
     * @code
     * cell->get_dof_values (dof_values, local_dof_values);
     * @endcode
     * 当前函数的意义在于，人们可以先修改这些局部值，例如应用限制器或确保所有节点值为正，然后再评估与当前单元上这些局部值相对应的有限元场。另一种应用是，人们希望将一个单元上的解后处理为每个单元上的不同的有限元空间，而不需要实际创建一个相应的DoFHandler
     *
     * 在这种情况下，我们所要计算的是该后处理函数的局部表示，其特征是节点值；然后该函数允许在正交点评估该表示。
     * @param[in]  dof_values
     * 一个本地节点值的向量。该向量的长度必须等于当前单元上的DoF数量，并且必须按照参考单元上自由度的编号顺序排列。
     * @param[out]  值
     * 给定的有限元场的值的向量，在当前对象上的正交点。
     * @tparam  InputVector  @p InputVector
     * 类型必须允许从中创建ArrayView对象；这一点由
     * `std::vector` 类和其他类满足。
     *
     */
    template <class InputVector>
    void
    get_function_values_from_local_dof_values(
      const InputVector &dof_values,
      std::vector<solution_value_type<typename InputVector::value_type>>
        &values) const;

    /**
     * 返回上次调用FEValues对象的<tt>reinit</tt>函数时选择的单元、面或子面的正交点上由<tt>fe_function</tt>表征的有限元函数的标量分量的梯度。
     * 这个函数等同于 FEValuesBase::get_function_gradients
     * 函数，但它只对选定的标量分量起作用。
     * 输出向量存储的数据类型必须是你将形状函数的梯度相乘后得到的数据类型（即
     * @p gradient_type) 乘以用于存储你的有限元向量 $U$
     * 的未知数 $U_j$ 值的类型（由 @p fe_function 参数代表）。
     * @dealiiRequiresUpdateFlags{update_gradients}
     *
     */
    template <class InputVector>
    void
    get_function_gradients(
      const InputVector &fe_function,
      std::vector<solution_gradient_type<typename InputVector::value_type>>
        &gradients) const;

    /**
     * 这个函数与get_function_gradients()的关系，与get_function_values_from_local_dof_values()与get_function_values()的关系相同。参见get_function_values_from_local_dof_values()的文档以了解更多信息。
     *
     */
    template <class InputVector>
    void
    get_function_gradients_from_local_dof_values(
      const InputVector &dof_values,
      std::vector<solution_gradient_type<typename InputVector::value_type>>
        &gradients) const;

    /**
     * 返回上次调用FEValues对象的<tt>reinit</tt>函数时选择的单元、面或子面的正交点上由<tt>fe_function</tt>表征的有限元函数的Hessians。
     * 这个函数等同于 FEValuesBase::get_function_hessians
     * 函数，但它只对选定的标量分量起作用。
     * 输出向量所存储的数据类型必须是你将形状函数的Hessians相乘后得到的（即
     * @p hessian_type) 乘以用于存储你的有限元向量 $U$
     * 的未知数 $U_j$ 值的类型（由 @p fe_function 参数代表）。
     * @dealiiRequiresUpdateFlags{update_hessians}
     *
     */
    template <class InputVector>
    void
    get_function_hessians(
      const InputVector &fe_function,
      std::vector<solution_hessian_type<typename InputVector::value_type>>
        &hessians) const;

    /**
     * 这个函数与get_function_hessians()的关系，与get_function_values_from_local_dof_values()与get_function_values()的关系相同。参见get_function_values_from_local_dof_values()的文档以了解更多信息。
     *
     */
    template <class InputVector>
    void
    get_function_hessians_from_local_dof_values(
      const InputVector &dof_values,
      std::vector<solution_hessian_type<typename InputVector::value_type>>
        &hessians) const;


    /**
     * 返回上次调用FEValues对象的<tt>reinit</tt>函数时选择的单元、面或子面的正交点上由<tt>fe_function</tt>表征的有限元函数标量分量的拉普拉斯。Laplacians是Hessians的轨迹。
     * 这个函数相当于 FEValuesBase::get_function_laplacians
     * 函数，但它只对选定的标量分量起作用。
     * 输出向量所存储的数据类型必须是你将形状函数的拉普拉斯系数相乘后得到的（即
     * @p value_type) 乘以用于存储你的有限元向量 $U_j$
     * 的未知数值的类型 $U$ （由 @p fe_function 参数代表）。
     * @dealiiRequiresUpdateFlags{update_hessians}
     *
     */
    template <class InputVector>
    void
    get_function_laplacians(
      const InputVector &fe_function,
      std::vector<solution_laplacian_type<typename InputVector::value_type>>
        &laplacians) const;

    /**
     * 这个函数与get_function_laplacians()的关系，与get_function_values_from_local_dof_values()与get_function_values()的关系相同。参见get_function_values_from_local_dof_values()的文档以了解更多信息。
     *
     */
    template <class InputVector>
    void
    get_function_laplacians_from_local_dof_values(
      const InputVector &dof_values,
      std::vector<solution_laplacian_type<typename InputVector::value_type>>
        &laplacians) const;


    /**
     * 返回上次调用FEValues对象的<tt>reinit</tt>函数时选择的单元、面或子面的正交点上由<tt>fe_function</tt>表征的有限元函数的标量分量的三次导数。
     * 这个函数相当于 FEValuesBase::get_function_third_derivatives
     * 函数，但它只对选定的标量分量起作用。
     * 输出向量所存储的数据类型必须是你将形状函数的三次导数（即
     * @p  third_derivative_type）乘以你的有限元向量 $U$ （由 @p
     * fe_function参数表示）用于存储未知数 $U_j$
     * 的值的类型。
     * @dealiiRequiresUpdateFlags{update_third_derivatives}
     *
     */
    template <class InputVector>
    void
    get_function_third_derivatives(
      const InputVector &fe_function,
      std::vector<
        solution_third_derivative_type<typename InputVector::value_type>>
        &third_derivatives) const;

    /**
     * 这个函数与get_function_third_derivatives()的关系，与get_function_values_from_local_dof_values()与get_function_values()的关系相同。参见get_function_values_from_local_dof_values()的文档以了解更多信息。
     *
     */
    template <class InputVector>
    void
    get_function_third_derivatives_from_local_dof_values(
      const InputVector &dof_values,
      std::vector<
        solution_third_derivative_type<typename InputVector::value_type>>
        &third_derivatives) const;


  private:
    /**
     * 一个指向我们操作的FEValuesBase对象的指针。
     *
     */
    const SmartPointer<const FEValuesBase<dim, spacedim>> fe_values;

    /**
     * 这个视图代表FEValuesBase对象的单一标量组件。
     *
     */
    const unsigned int component;

    /**
     * 存储有关形状函数的数据。
     *
     */
    std::vector<ShapeFunctionData> shape_function_data;
  };



  /**
   * 一个代表对一组 <code>spacedim</code> 分量的视图的类，这些分量构成了一个矢量值的有限元的矢量部分。视图将在 @ref vector_valued 模块中讨论。    请注意，在目前的上下文中，矢量是指物理学上使用的意义：它有 <code>spacedim</code> 个分量，在坐标系变换下以特定的方式表现出来。例子包括速度或位移场。这与数学中使用 "向量 "
   * 一词的方式相反（以及我们在库中的其他上下文中使用这个词的方式，例如在向量类中），在那里它真正代表了一个数字的集合。后者的一个例子是火焰中化学物种浓度的集合；然而，这些实际上只是标量变量的集合，因为如果坐标系被旋转，它们不会改变，不像速度矢量的分量，因此，这个类不应该被用于这种情况。
   * 该类允许查询代表矢量的形状函数和解决方案的（分量）的值、梯度和发散。一个向量的梯度
   * $d_{k}, 0\le k<\text{dim}$  被定义为  $S_{ij} = \frac{\partial
   * d_{i}}{\partial x_j}, 0\le i,j<\text{dim}$  。    如果你将
   * FEValuesExtractors::Vector
   * 应用于FEValues、FEFaceValues或FESubfaceValues对象，你会得到这种类型的对象。
   * @ingroup feaccess vector_valued
   *
   */
  template <int dim, int spacedim = dim>
  class Vector
  {
  public:
    /**
     * 这个类所代表的视图值的数据类型的别名。因为我们处理的是一组
     * <code>dim</code> 的组件，所以值的类型是Tensor<1,spacedim>。
     *
     */
    using value_type = dealii::Tensor<1, spacedim>;

    /**
     * 这个类所代表的视图的梯度类型的别名。
     * 这里，对于一组 <code>dim</code>
     * 分量的有限元，梯度是一个 <code>Tensor@<2,spacedim@></code>
     * 。
     * 关于向量的梯度到底是如何定义的，请看这个类的一般文档。
     *
     */
    using gradient_type = dealii::Tensor<2, spacedim>;

    /**
     * 这个类所代表的视图的对称梯度的类型的别名。这里，对于一组
     * <code>dim</code> 分量的有限元，对称梯度是一个
     * <code>SymmetricTensor@<2,spacedim@></code>  。        一个矢量场
     * $\mathbf v$ 的对称梯度定义为 $\varepsilon(\mathbf v)=\frac 12
     * (\nabla \mathbf v + \nabla \mathbf v^T)$  .
     *
     */
    using symmetric_gradient_type = dealii::SymmetricTensor<2, spacedim>;

    /**
     * 该类代表的视图发散的类型的别名。这里，对于一组
     * <code>dim</code> 分量的有限元，发散当然是一个标量。
     *
     */
    using divergence_type = double;

    /**
     * 这个类所代表的视图的卷曲类型的别名。
     * 这里，对于一组 <code>spacedim=2</code>
     * 的有限元分量，curl是一个 <code>Tensor@<1, 1@></code>
     * 。对于 <code>spacedim=3</code> it is a <code>Tensor@<1, dim@></code>
     * .
     *
     */
    using curl_type = typename dealii::internal::CurlType<spacedim>::type;

    /**
     * 这个类所代表的视图的二阶导数的类型的别名。这里，对于一组
     * <code>dim</code> 分量的有限元，Hessian是一个
     * <code>Tensor@<3,dim@></code>  。
     *
     */
    using hessian_type = dealii::Tensor<3, spacedim>;

    /**
     * 该类代表的视图的第三导数类型的别名。这里，对于有限元的一组
     * <code>dim</code> 分量，第三导数是一个
     * <code>Tensor@<4,dim@></code>  。
     *
     */
    using third_derivative_type = dealii::Tensor<4, spacedim>;

    /**
     * a  @p Number
     * 和该类提供的视图值的乘积的数据类型的别名。这是一个有限元场的矢量分量的数据类型，其自由度由元素类型为
     * @p Number. 的矢量描述。
     *
     */
    template <typename Number>
    using solution_value_type = typename ProductType<Number, value_type>::type;

    /**
     * @p Number
     * 和该类提供的视图梯度的乘积的数据类型的别名。这是一个有限元场的矢量分量的数据类型，其自由度由元素类型为
     * @p Number. 的矢量描述。
     *
     */
    template <typename Number>
    using solution_gradient_type =
      typename ProductType<Number, gradient_type>::type;

    /**
     * 一个 @p Number
     * 和该类提供的视图的对称梯度的乘积的数据类型的别名。这是一个有限元场的矢量分量的数据类型，其自由度由元素类型为
     * @p Number. 的矢量描述。
     *
     */
    template <typename Number>
    using solution_symmetric_gradient_type =
      typename ProductType<Number, symmetric_gradient_type>::type;

    /**
     * 一个 @p Number
     * 和该类提供的视图的分歧的乘积的数据类型的别名。这是一个有限元场的矢量分量的数据类型，其自由度由元素类型为
     * @p Number. 的矢量描述。
     *
     */
    template <typename Number>
    using solution_divergence_type =
      typename ProductType<Number, divergence_type>::type;

    /**
     * 一个 @p Number
     * 和该类提供的视图的拉普拉斯的乘积的数据类型的别名。这是一个有限元场的矢量分量的数据类型，其自由度由元素类型为
     * @p Number. 的矢量描述。
     *
     */
    template <typename Number>
    using solution_laplacian_type =
      typename ProductType<Number, value_type>::type;

    /**
     * 一个 @p Number
     * 和该类提供的视图的卷积的数据类型的别名。这是一个有限元场的矢量分量的数据类型，其自由度由元素类型为
     * @p Number. 的矢量描述。
     *
     */
    template <typename Number>
    using solution_curl_type = typename ProductType<Number, curl_type>::type;

    /**
     * 一个 @p Number
     * 和该类提供的视图的赫西的乘积的数据类型的别名。这是一个有限元场的矢量分量的数据类型，其自由度由元素类型为
     * @p Number. 的矢量描述。
     *
     */
    template <typename Number>
    using solution_hessian_type =
      typename ProductType<Number, hessian_type>::type;

    /**
     * 一个 @p Number
     * 和该类提供的视图的第三导数的乘积的数据类型的别名。这是一个有限元场的矢量分量的数据类型，其自由度由元素类型为
     * @p Number. 的矢量描述。
     *
     */
    template <typename Number>
    using solution_third_derivative_type =
      typename ProductType<Number, third_derivative_type>::type;

    /**
     * 一个结构，为矢量视图和任何 @p Number
     * 类型的基函数的值和导数的乘积提供输出类型。
     * @deprecated  使用周围类中定义的类型来代替。
     *
     */
    template <typename Number>
    struct DEAL_II_DEPRECATED OutputType
    {
      /**
       * 一个 @p Number
       * 和Vector类的视图值的乘积的数据类型的别名。
       *
       */
      using value_type =
        typename ProductType<Number,
                             typename Vector<dim, spacedim>::value_type>::type;

      /**
       * 一个 @p Number
       * 和视图矢量类的梯度的乘积的数据类型的别名。
       *
       */
      using gradient_type = typename ProductType<
        Number,
        typename Vector<dim, spacedim>::gradient_type>::type;

      /**
       * @p Number
       * 与向量类的对称梯度的乘积的数据类型的别名。
       *
       */
      using symmetric_gradient_type = typename ProductType<
        Number,
        typename Vector<dim, spacedim>::symmetric_gradient_type>::type;

      /**
       * 一个 @p Number
       * 和向量类视图的发散的数据类型的别名。
       *
       */
      using divergence_type = typename ProductType<
        Number,
        typename Vector<dim, spacedim>::divergence_type>::type;

      /**
       * 一个 @p Number
       * 和矢量类视图的拉普拉斯的乘积的数据类型的别名。
       *
       */
      using laplacian_type =
        typename ProductType<Number,
                             typename Vector<dim, spacedim>::value_type>::type;

      /**
       * 一个 @p Number
       * 和向量类视图的卷积的数据类型的别名。
       *
       */
      using curl_type =
        typename ProductType<Number,
                             typename Vector<dim, spacedim>::curl_type>::type;

      /**
       * 一个 @p Number
       * 和视图向量类的Hessians的乘积的数据类型的别名。
       *
       */
      using hessian_type = typename ProductType<
        Number,
        typename Vector<dim, spacedim>::hessian_type>::type;

      /**
       * 一个 @p Number
       * 和向量类视图的第三导数的乘积的数据类型的别名。
       *
       */
      using third_derivative_type = typename ProductType<
        Number,
        typename Vector<dim, spacedim>::third_derivative_type>::type;
    };

    /**
     * 一个结构，对于每个形状函数，我们预先计算出一堆数据，这将使以后的访问变得更加便宜。
     *
     */
    struct ShapeFunctionData
    {
      /**
       * 对于每一对(形状函数,向量内的分量)，存储所选向量分量是否可能为非零。对于原始形状函数，我们可以肯定地知道某个给定形状函数的某个标量分量是否为非零，而对于非原始形状函数，这可能并不完全清楚（例如，对于RT元素，它取决于单元格的形状）。
       *
       */
      bool is_nonzero_shape_function_component[spacedim];

      /**
       * 对于每一对（形状函数，矢量内的组件），在shape_values、shape_gradients和shape_hessians表中存储行索引（列索引是正交点索引）。如果形状函数是原始的，那么我们可以从FEValues对象的shape_function_to_row_table中获得这些信息；否则，我们必须花点功夫来计算这些信息。
       *
       */
      unsigned int row_index[spacedim];

      /**
       * 对于每个形状函数说如下：如果这个形状函数的is_nonzero_shape_function_component中只有一个条目是非零的，那么就存储row_index的相应值，single_nonzero_component_index代表在0和dim之间的索引，对它来说，是达到了非零的。如果多个分量不为零，那么就存储
       *
       * 如果没有组件是非零的，则存储 * - 。
       *
       * - .
       *
       */
      int          single_nonzero_component;
      unsigned int single_nonzero_component_index;
    };

    /**
     * 默认构造函数。创建一个无效的对象。
     *
     */
    Vector();

    /**
     * 表示FEValuesBase对象（或从FEValuesBase派生的类之一）的dim分量的对象的构造函数，代表一个向量值变量。
     * 第二个参数表示所选向量的第一个分量的索引。
     *
     */
    Vector(const FEValuesBase<dim, spacedim> &fe_values_base,
           const unsigned int                 first_vector_component);

    /**
     * 复制构造函数。这不是一个轻量级的对象，所以我们不允许复制，如果调用这个函数会产生一个编译时错误。
     *
     */
    Vector(const Vector<dim, spacedim> &) = delete;

    /**
     * 移动构造器。
     *
     */
    // NOLINTNEXTLINE OSX does not compile with noexcept
    Vector(Vector<dim, spacedim> &&) = default;

    /**
     * 解构器。
     *
     */
    ~Vector() = default;

    /**
     * 拷贝操作符。这不是一个轻量级的对象，所以我们不允许复制，如果调用这个函数会产生一个编译时错误。
     *
     */
    Vector &
    operator=(const Vector<dim, spacedim> &) = delete;

    /**
     * 移动赋值运算符。
     *
     */
    // NOLINTNEXTLINE OSX does not compile with noexcept
    Vector &
    operator=(Vector<dim, spacedim> &&) = default; // NOLINT

    /**
     * 返回该视图所选择的向量分量的值，用于参数所选择的形状函数和正交点。
     * 这里，由于视图代表FEValues对象中具有 <code>dim</code>
     * 分量的向量值部分，所以返回类型是具有 <code>dim</code>
     * 分量的秩1张量。          @param  shape_function
     * 要评估的形状函数的编号。
     * 注意，这个数字从零到dofs_per_cell，即使是在FEFaceValues或FESubfaceValues对象的情况下。
     * @param  q_point 要评估函数的正交点的编号。
     * @dealiiRequiresUpdateFlags{update_values}
     *
     */
    value_type
    value(const unsigned int shape_function, const unsigned int q_point) const;

    /**
     * 返回该视图所选择的矢量分量的梯度（等级为2的张量），用于参数所选择的形状函数和正交点。
     * 关于向量的梯度到底是如何定义的，请看这个类的一般文档。
     * @note  参数的含义与value()函数的记录相同。
     * @dealiiRequiresUpdateFlags{update_gradients}
     *
     */
    gradient_type
    gradient(const unsigned int shape_function,
             const unsigned int q_point) const;

    /**
     * 返回该视图选择的向量分量的对称梯度（等级为2的对称张量），用于参数选择的形状函数和正交点。
     * 对称梯度定义为 $\frac 12 [(\nabla \phi_i(x_q)) + (\nabla
     * \phi_i(x_q))^T]$  ，其中 $\phi_i$
     * 代表从FEValuesBase对象中选择的 <code>dim</code> 分量，
     * $x_q$ 是第 $q$ 个正交点的位置。
     * @note  参数的含义与value()函数的记载相同。
     * @dealiiRequiresUpdateFlags{update_gradients}
     *
     */
    symmetric_gradient_type
    symmetric_gradient(const unsigned int shape_function,
                       const unsigned int q_point) const;

    /**
     * 返回该视图所选择的向量分量的标量发散，用于参数所选择的形状函数和正交点。
     * @note  参数的含义与value()函数的记录相同。
     * @dealiiRequiresUpdateFlags{update_gradients}
     *
     */
    divergence_type
    divergence(const unsigned int shape_function,
               const unsigned int q_point) const;

    /**
     * 对于参数选择的形状函数和正交点，返回该视图选择的矢量分量的向量卷曲。    对于1d来说，这个函数没有任何意义。因此，它没有在  <code>spacedim=1</code>  中实现。 在2D中，卷曲被定义为@f{equation*}{
     * \operatorname{curl}(u) \dealcoloneq \frac{du_2}{dx}
     *
     * -\frac{du_1}{dy},
     * @f} 。
     * 而在三维中，它是由@f{equation*}{
     * \operatorname{curl}(u) \dealcoloneq \left( \begin{array}{c}
     * \frac{du_3}{dy}-\frac{du_2}{dz}\\ \frac{du_1}{dz}-\frac{du_3}{dx}\\
     * \frac{du_2}{dx}-\frac{du_1}{dy} \end{array} \right).
     * @f}给出的。
     *
     * @note  参数的含义与value()函数的记录相同。
     * @dealiiRequiresUpdateFlags{update_gradients}
     *
     */
    curl_type
    curl(const unsigned int shape_function, const unsigned int q_point) const;

    /**
     * 对于参数选择的形状函数和正交点，返回该视图选择的向量分量的Hessian（所有二次导数的等级2的张量）。
     * @note  参数的含义与value()函数的记录相同。
     * @dealiiRequiresUpdateFlags{update_hessians}
     *
     */
    hessian_type
    hessian(const unsigned int shape_function,
            const unsigned int q_point) const;

    /**
     * 返回该视图选择的向量分量的所有三次导数的等级3的张量，用于参数选择的形状函数和正交点。
     * @note 参数的含义与value()函数的记录相同。
     * @dealiiRequiresUpdateFlags{update_3rd_derivatives}
     *
     */
    third_derivative_type
    third_derivative(const unsigned int shape_function,
                     const unsigned int q_point) const;

    /**
     * 返回上次调用FEValues对象的<tt>reinit</tt>函数时选择的单元、面或子面的正交点上以<tt>fe_function</tt>为特征的有限元函数的选定向量分量的值。
     * 这个函数相当于 FEValuesBase::get_function_values
     * 函数，但它只对选定的向量分量起作用。
     * 输出向量所存储的数据类型必须是你将形状函数的值相乘后得到的（即
     * @p value_type) 乘以用于存储你的有限元向量 $U$ 的未知数
     * $U_j$ 值的类型（由 @p fe_function 参数代表）。
     * @dealiiRequiresUpdateFlags{update_values}
     *
     */
    template <class InputVector>
    void
    get_function_values(
      const InputVector &fe_function,
      std::vector<solution_value_type<typename InputVector::value_type>>
        &values) const;

    /**
     * 与上述相同，但使用局部自由度值的向量。换句话说，这个函数不是从与DoFHandler对象相关的全局向量中提取位于当前单元上的自由度的节点值（如上面的函数），而是通过其第一个参数获取这些局部节点值。获得这样一个向量的典型方法是调用如下代码
     * @code
     * cell->get_dof_values (dof_values, local_dof_values);
     * @endcode
     * 当前函数的意义在于，人们可以先修改这些局部值，例如应用限制器或确保所有节点值为正，然后再评估当前单元上与这些局部值相对应的有限元场。另一种应用是，人们希望将一个单元上的解后处理为每个单元上的不同的有限元空间，而不需要实际创建一个相应的DoFHandler
     *
     * 在这种情况下，我们所要计算的是该后处理函数的局部表示，其特征是节点值；然后该函数允许在正交点评估该表示。
     * @param[in]  dof_values
     * 一个本地节点值的向量。该向量的长度必须等于当前单元上的DoF数量，并且必须按照参考单元上自由度的编号顺序排列。
     * @param[out]  值
     * 给定的有限元场的值的向量，在当前对象上的正交点。
     * @tparam  InputVector  @p InputVector
     * 类型必须允许从中创建ArrayView对象；这一点由
     * `std::vector` 类等满足。
     *
     */
    template <class InputVector>
    void
    get_function_values_from_local_dof_values(
      const InputVector &dof_values,
      std::vector<solution_value_type<typename InputVector::value_type>>
        &values) const;

    /**
     * 返回上次调用FEValues对象的<tt>reinit</tt>函数时选择的单元、面或子面的正交点上由<tt>fe_function</tt>表征的有限元函数的选定矢量成分的梯度。
     * 这个函数相当于 FEValuesBase::get_function_gradients
     * 函数，但它只对选定的向量分量起作用。
     * 输出向量所存储的数据类型必须是你将形状函数的梯度相乘后得到的数据（即
     * @p gradient_type) 乘以用于存储你的有限元向量 $U_j$
     * 的未知数值的类型 $U$ （由 @p fe_function 参数代表）。
     * @dealiiRequiresUpdateFlags{update_gradients}
     *
     */
    template <class InputVector>
    void
    get_function_gradients(
      const InputVector &fe_function,
      std::vector<solution_gradient_type<typename InputVector::value_type>>
        &gradients) const;

    /**
     * 这个函数与get_function_gradients()的关系，与get_function_values_from_local_dof_values()与get_function_values()的关系相同。参见get_function_values_from_local_dof_values()的文档以了解更多信息。
     *
     */
    template <class InputVector>
    void
    get_function_gradients_from_local_dof_values(
      const InputVector &dof_values,
      std::vector<solution_gradient_type<typename InputVector::value_type>>
        &gradients) const;

    /**
     * 返回上次调用FEValues对象的<tt>reinit</tt>函数时选择的单元、面或子面的正交点上由<tt>fe_function</tt>描述的有限元函数的对称梯度。
     * 矢量场的对称梯度 $\mathbf v$ 定义为 $\varepsilon(\mathbf
     * v)=\frac 12 (\nabla \mathbf v + \nabla \mathbf v^T)$  。
     * @note  在FEValues类中没有
     * FEValuesBase::get_function_symmetric_gradients
     * 这样的等效函数，但当然可以从
     * FEValuesBase::get_function_gradients, 中获得相关信息。
     * 输出向量存储的数据类型必须是你将形状函数的对称梯度（即
     * @p  symmetric_gradient_type）乘以用于存储你的有限元向量
     * $U_j$ 的未知数值的类型 $U$ （由 @p
     * fe_function参数表示）时得到的。
     * @dealiiRequiresUpdateFlags{update_gradients}
     *
     */
    template <class InputVector>
    void
    get_function_symmetric_gradients(
      const InputVector &fe_function,
      std::vector<
        solution_symmetric_gradient_type<typename InputVector::value_type>>
        &symmetric_gradients) const;

    /**
     * 这个函数与get_function_symmetric_gradients()的关系，与get_function_values_from_local_dof_values()与get_function_values()的关系相同。参见get_function_values_from_local_dof_values()的文档以了解更多信息。
     *
     */
    template <class InputVector>
    void
    get_function_symmetric_gradients_from_local_dof_values(
      const InputVector &dof_values,
      std::vector<
        solution_symmetric_gradient_type<typename InputVector::value_type>>
        &symmetric_gradients) const;

    /**
     * 返回上次调用FEValues对象的<tt>reinit</tt>函数时选择的单元、面或子面的正交点上由<tt>fe_function</tt>表征的有限元函数的选定矢量分量的发散。
     * 在FEValues类中没有诸如 FEValuesBase::get_function_divergences
     * 这样的等效函数，但当然可以从
     * FEValuesBase::get_function_gradients, 中获得相关信息。
     * 输出向量所存储的数据类型必须是你将形状函数的发散量相乘后得到的（即
     * @p divergence_type) 乘以用于存储你的有限元向量 $U$ （由
     * @p fe_function 参数表示）的未知数 $U_j$ 的值的类型）。
     * @dealiiRequiresUpdateFlags{update_gradients}
     *
     */
    template <class InputVector>
    void
    get_function_divergences(
      const InputVector &fe_function,
      std::vector<solution_divergence_type<typename InputVector::value_type>>
        &divergences) const;

    /**
     * 这个函数与get_function_divergences()的关系，与get_function_values_from_local_dof_values()与get_function_values()的关系相同。更多信息请参见get_function_values_from_local_dof_values()的文档。
     *
     */
    template <class InputVector>
    void
    get_function_divergences_from_local_dof_values(
      const InputVector &dof_values,
      std::vector<solution_divergence_type<typename InputVector::value_type>>
        &divergences) const;

    /**
     * 返回上次调用FEValues对象的<tt>reinit</tt>函数时选择的单元、面或子面的正交点上由<tt>fe_function</tt>表征的有限元函数的所选向量分量的卷曲。
     * 在FEValues类中没有诸如 FEValuesBase::get_function_curls
     * 这样的等效函数，但当然可以从
     * FEValuesBase::get_function_gradients, 中获得相关信息。
     * 输出向量所存储的数据类型必须是你将形状函数的卷曲相乘后得到的（即
     * @p curl_type) 乘以用于存储你的有限元向量 $U$ （由 @p
     * fe_function 参数表示）的未知数 $U_j$ 的值的类型）。
     * @dealiiRequiresUpdateFlags{update_gradients}
     *
     */
    template <class InputVector>
    void
    get_function_curls(
      const InputVector &fe_function,
      std::vector<solution_curl_type<typename InputVector::value_type>> &curls)
      const;

    /**
     * 这个函数与get_function_curls()的关系，与get_function_values_from_local_dof_values()与get_function_values()的关系相同。参见get_function_values_from_local_dof_values()的文档以了解更多信息。
     *
     */
    template <class InputVector>
    void
    get_function_curls_from_local_dof_values(
      const InputVector &dof_values,
      std::vector<solution_curl_type<typename InputVector::value_type>> &curls)
      const;

    /**
     * 返回上次调用FEValues对象的<tt>reinit</tt>函数时选择的单元、面或子面的正交点上由<tt>fe_function</tt>表征的有限元函数的选定矢量分量的Hessians。
     * 这个函数相当于 FEValuesBase::get_function_hessians
     * 函数，但它只对选定的向量分量起作用。
     * 输出向量所存储的数据类型必须是你将形状函数的Hessians相乘后得到的（即
     * @p hessian_type) 乘以用于存储你的有限元向量 $U_j$
     * 的未知数值的类型 $U$ （由 @p fe_function 参数代表）。
     * @dealiiRequiresUpdateFlags{update_hessians}
     *
     */
    template <class InputVector>
    void
    get_function_hessians(
      const InputVector &fe_function,
      std::vector<solution_hessian_type<typename InputVector::value_type>>
        &hessians) const;

    /**
     * 这个函数与get_function_hessians()的关系，与get_function_values_from_local_dof_values()与get_function_values()的关系相同。参见get_function_values_from_local_dof_values()的文档以了解更多信息。
     *
     */
    template <class InputVector>
    void
    get_function_hessians_from_local_dof_values(
      const InputVector &dof_values,
      std::vector<solution_hessian_type<typename InputVector::value_type>>
        &hessians) const;

    /**
     * 返回上次调用FEValues对象的<tt>reinit</tt>函数时选择的单元、面或子面的正交点上由<tt>fe_function</tt>表征的有限元函数的选定向量分量的拉普拉斯。Laplacians是Hessians的轨迹。
     * 这个函数相当于 FEValuesBase::get_function_laplacians
     * 函数，但它只对选定的向量分量起作用。
     * 输出向量所存储的数据类型必须是你将形状函数的拉普拉斯系数相乘后得到的数据（即
     * @p laplacian_type) 乘以用于存储你的有限元向量 $U$ （由
     * @p fe_function 参数表示）的未知数 $U_j$ 的值的类型。
     * @dealiiRequiresUpdateFlags{update_hessians}
     *
     */
    template <class InputVector>
    void
    get_function_laplacians(
      const InputVector &fe_function,
      std::vector<solution_laplacian_type<typename InputVector::value_type>>
        &laplacians) const;

    /**
     * 这个函数与get_function_laplacians()的关系，与get_function_values_from_local_dof_values()与get_function_values()的关系相同。参见get_function_values_from_local_dof_values()的文档以了解更多信息。
     *
     */
    template <class InputVector>
    void
    get_function_laplacians_from_local_dof_values(
      const InputVector &dof_values,
      std::vector<solution_laplacian_type<typename InputVector::value_type>>
        &laplacians) const;

    /**
     * 返回上次调用FEValues对象的<tt>reinit</tt>函数时选择的单元、面或子面的正交点上由<tt>fe_function</tt>表征的有限元函数的标量分量的三次导数。
     * 这个函数相当于 FEValuesBase::get_function_third_derivatives
     * 函数，但它只对选定的标量分量起作用。
     * 输出向量所存储的数据类型必须是你将形状函数的三次导数（即
     * @p  third_derivative_type）乘以用于存储你的有限元向量 $U$
     * （由 @p  fe_function参数表示）的未知数 $U_j$
     * 的值的类型。
     * @dealiiRequiresUpdateFlags{update_third_derivatives}
     *
     */
    template <class InputVector>
    void
    get_function_third_derivatives(
      const InputVector &fe_function,
      std::vector<
        solution_third_derivative_type<typename InputVector::value_type>>
        &third_derivatives) const;

    /**
     * 这个函数与get_function_third_derivatives()的关系，与get_function_values_from_local_dof_values()与get_function_values()的关系相同。参见get_function_values_from_local_dof_values()的文档以了解更多信息。
     *
     */
    template <class InputVector>
    void
    get_function_third_derivatives_from_local_dof_values(
      const InputVector &dof_values,
      std::vector<
        solution_third_derivative_type<typename InputVector::value_type>>
        &third_derivatives) const;

  private:
    /**
     * 一个指向我们操作的FEValuesBase对象的指针。
     *
     */
    const SmartPointer<const FEValuesBase<dim, spacedim>> fe_values;

    /**
     * 这个视图代表FEValuesBase对象的向量的第一个分量。
     *
     */
    const unsigned int first_vector_component;

    /**
     * 存储有关形状函数的数据。
     *
     */
    std::vector<ShapeFunctionData> shape_function_data;
  };


  template <int rank, int dim, int spacedim = dim>
  class SymmetricTensor;

  /**
   * 一个代表对一组 <code>(dim*dim
   * + dim)/2</code>
   * 分量的视图的类，该分量形成了一个来自矢量值有限元的对称二阶张量。视图将在
   * @ref vector_valued 模块中讨论。
   * 该类允许查询代表对称张量的形状函数和解决方案的（组件）的值和发散。对称张量的发散
   * $S_{ij}, 0\le i,j<\text{dim}$  被定义为  $d_i = \sum_j \frac{\partial
   * S_{ij}}{\partial x_j}, 0\le i<\text{dim}$
   * ，由于张量的对称性，它也是  $d_i = \sum_j \frac{\partial
   * S_{ji}}{\partial x_j}$  。 换句话说，由于 $S$
   * 的对称性，我们是按行还是按列应用纳布拉算子来得到发散并不重要。
   * 如果你将 FEValuesExtractors::SymmetricTensor
   * 应用于FEValues、FEFaceValues或FESubfaceValues对象，你会得到一个这种类型的对象。
   * @ingroup feaccess vector_valued
   *
   */
  template <int dim, int spacedim>
  class SymmetricTensor<2, dim, spacedim>
  {
  public:
    /**
     * 这个类所代表的视图值的数据类型的别名。由于我们处理的是一组
     * <code>(dim*dim + dim)/2</code>
     * 成分（即对称二阶张量的唯一成分），所以数值类型是SymmetricTensor<2,spacedim>。
     *
     */
    using value_type = dealii::SymmetricTensor<2, spacedim>;

    /**
     * 这个类所代表的视图的发散类型的别名。这里，对于一组代表对称二阶张量的有限元的
     * <code>(dim*dim + dim)/2</code> 唯一分量，发散当然是一个
     * <code>Tensor@<1,dim@></code>  。
     * 关于发散的定义，请参见该类的一般讨论。
     *
     */
    using divergence_type = dealii::Tensor<1, spacedim>;

    /**
     * a  @p Number
     * 与本类提供的视图值之积的数据类型的别名。这是一个有限元场的矢量分量的数据类型，其自由度由元素类型为
     * @p Number. 的矢量描述。
     *
     */
    template <typename Number>
    using solution_value_type = typename ProductType<Number, value_type>::type;

    /**
     * 一个 @p Number
     * 和该类提供的视图的分歧的乘积的数据类型的别名。这是一个有限元场的矢量分量的数据类型，其自由度由元素类型为
     * @p Number. 的矢量描述。
     *
     */
    template <typename Number>
    using solution_divergence_type =
      typename ProductType<Number, divergence_type>::type;


    /**
     * 一个结构，为SymmetricTensor视图和任何 @p Number
     * 类型的基函数的值和导数的乘积提供输出类型。
     * @deprecated  使用周围类中定义的类型来代替。
     *
     */
    template <typename Number>
    struct DEAL_II_DEPRECATED OutputType
    {
      /**
       * @p Number
       * 与SymmetricTensor类的视图值之积的数据类型的别名。
       *
       */
      using value_type = typename ProductType<
        Number,
        typename SymmetricTensor<2, dim, spacedim>::value_type>::type;

      /**
       * 一个 @p Number
       * 的乘积的数据类型和SymmetricTensor类视图的分歧的别名。
       *
       */
      using divergence_type = typename ProductType<
        Number,
        typename SymmetricTensor<2, dim, spacedim>::divergence_type>::type;
    };

    /**
     * 一个结构，对于每个形状函数，我们预先计算出一堆数据，这将使以后的访问变得更加便宜。
     *
     */
    struct ShapeFunctionData
    {
      /**
       * 对于每一对(形状函数,向量内的分量)，存储所选向量分量是否可能为非零。对于原始形状函数，我们肯定知道某个给定形状函数的某个标量分量是否为非零，而对于非原始形状函数，这可能并不完全清楚（例如，对于RT元素，它取决于单元格的形状）。
       *
       */
      bool is_nonzero_shape_function_component
        [value_type::n_independent_components];

      /**
       * 对于每一对（形状函数，矢量内的组件），在shape_values、shape_gradients和shape_hessians表中存储行索引（列索引是正交点索引）。如果形状函数是原始的，那么我们可以从FEValues对象的shape_function_to_row_table中获得这些信息；否则，我们必须花点功夫来计算这些信息。
       *
       */
      unsigned int row_index[value_type::n_independent_components];

      /**
       * 对于每个形状函数说如下：如果这个形状函数的is_nonzero_shape_function_component中只有一个条目是非零的，那么就存储row_index的相应值，single_nonzero_component_index代表在0和（dim^2+dim）/2之间的索引，对它来说是达到的。如果多个分量为非零，那么就存储
       *
       * 如果没有分量是非零的，则存储 * - 。
       *
       * - .
       *
       */
      int single_nonzero_component;

      /**
       * @p single_nonzero_component 的索引。
       *
       */
      unsigned int single_nonzero_component_index;
    };

    /**
     * 默认构造函数。创建一个无效的对象。
     *
     */
    SymmetricTensor();

    /**
     * 表示<code>(dim*dim + dim)/2</code>
     * FEValuesBase对象（或从FEValuesBase派生的类之一）的组件的构造函数，代表构成对称二阶张量值变量的唯一组件。
     * 第二个参数表示所选对称二阶张量的第一个分量的索引。
     *
     */
    SymmetricTensor(const FEValuesBase<dim, spacedim> &fe_values_base,
                    const unsigned int                 first_tensor_component);

    /**
     * 复制构造函数。这不是一个轻量级的对象，所以我们不允许复制，如果调用这个函数会产生一个编译时错误。
     *
     */
    SymmetricTensor(const SymmetricTensor<2, dim, spacedim> &) = delete;

    /**
     * 移动构造函数。
     *
     */
    // NOLINTNEXTLINE OSX does not compile with noexcept
    SymmetricTensor(SymmetricTensor<2, dim, spacedim> &&) = default;

    /**
     * 复制操作符。这不是一个轻量级的对象，所以我们不允许复制，如果调用这个函数会产生一个编译时错误。
     *
     */
    SymmetricTensor &
    operator=(const SymmetricTensor<2, dim, spacedim> &) = delete;

    /**
     * 移动赋值运算符。
     *
     */
    SymmetricTensor &
    operator=(SymmetricTensor<2, dim, spacedim> &&) noexcept = default;

    /**
     * 返回该视图所选择的向量分量的值，用于参数所选择的形状函数和正交点。
     * 这里，由于视图代表FEValues对象的一个矢量值部分，具有
     * <code>(dim*dim + dim)/2</code>
     * 分量（对称二阶张量的唯一分量），所以返回类型是等级2的对称张量。
     * @param  shape_function 要评估的形状函数的编号。
     * 注意，这个数字从零到dofs_per_cell，即使是在FEFaceValues或FESubfaceValues对象的情况下。
     * @param  q_point 要评估函数的正交点的编号。
     * @dealiiRequiresUpdateFlags{update_values}
     *
     */
    value_type
    value(const unsigned int shape_function, const unsigned int q_point) const;

    /**
     * 返回该视图所选择的矢量分量的矢量发散，对于参数所选择的形状函数和正交点。
     * 关于发散的定义，请参见本类的一般讨论。
     * @note  参数的含义与value()函数的记录相同。
     * @dealiiRequiresUpdateFlags{update_gradients}
     *
     */
    divergence_type
    divergence(const unsigned int shape_function,
               const unsigned int q_point) const;

    /**
     * 返回上次调用FEValues对象的<tt>reinit</tt>函数时选择的单元、面或子面的正交点上以<tt>fe_function</tt>为特征的有限元函数的选定向量分量的值。
     * 这个函数相当于 FEValuesBase::get_function_values
     * 函数，但它只对选定的向量分量起作用。
     * 输出向量所存储的数据类型必须是你将形状函数的值相乘后得到的（即
     * @p value_type) 乘以用于存储你的有限元向量 $U$ 的未知数
     * $U_j$ 值的类型（由 @p fe_function 参数代表）。
     * @dealiiRequiresUpdateFlags{update_values}
     *
     */
    template <class InputVector>
    void
    get_function_values(
      const InputVector &fe_function,
      std::vector<solution_value_type<typename InputVector::value_type>>
        &values) const;

    /**
     * 与上述相同，但使用局部自由度值的向量。换句话说，这个函数不是从与DoFHandler对象相关的全局向量中提取位于当前单元上的自由度的节点值（如上面的函数），而是通过其第一个参数获取这些局部节点值。获得这样一个向量的典型方法是调用如下代码
     * @code
     * cell->get_dof_values (dof_values, local_dof_values);
     * @endcode
     * 当前函数的意义在于，人们可以先修改这些局部值，例如应用限制器或确保所有节点值为正，然后再评估与当前单元上这些局部值相对应的有限元场。另一种应用是，人们希望将一个单元上的解后处理为每个单元上的不同的有限元空间，而不需要实际创建一个相应的DoFHandler
     *
     * 在这种情况下，我们所要计算的是该后处理函数的局部表示，其特征是节点值；然后该函数允许在正交点评估该表示。
     * @param[in]  dof_values
     * 一个本地节点值的向量。该向量的长度必须等于当前单元上的DoF数量，并且必须按照参考单元上自由度的编号顺序排列。
     * @param[out]  值
     * 给定的有限元场的值的向量，在当前对象的正交点上。
     * @tparam  InputVector  @p InputVector
     * 类型必须允许从中创建ArrayView对象；这一点由
     * `std::vector` 类等满足。
     *
     */
    template <class InputVector>
    void
    get_function_values_from_local_dof_values(
      const InputVector &dof_values,
      std::vector<solution_value_type<typename InputVector::value_type>>
        &values) const;

    /**
     * 返回上次调用FEValues对象的<tt>reinit</tt>函数时选择的单元、面或子面的正交点上由<tt>fe_function</tt>表征的有限元函数的选定向量分量的发散。
     * 在FEValues类中没有诸如 FEValuesBase::get_function_divergences
     * 这样的等效函数，但当然可以从
     * FEValuesBase::get_function_gradients, 中获得信息。
     * 关于发散的定义，请参见该类的一般讨论。
     * 输出向量存储的数据类型必须是你将形状函数的发散量相乘后得到的（即
     * @p divergence_type) 乘以用于存储你的有限元向量 $U_j$
     * 的未知数值的类型 $U$ （由 @p fe_function 参数代表）。
     * @dealiiRequiresUpdateFlags{update_gradients}
     *
     */
    template <class InputVector>
    void
    get_function_divergences(
      const InputVector &fe_function,
      std::vector<solution_divergence_type<typename InputVector::value_type>>
        &divergences) const;

    /**
     * 这个函数与get_function_divergences()的关系，与get_function_values_from_local_dof_values()与get_function_values()的关系相同。参见get_function_values_from_local_dof_values()的文档以了解更多信息。
     *
     */
    template <class InputVector>
    void
    get_function_divergences_from_local_dof_values(
      const InputVector &dof_values,
      std::vector<solution_divergence_type<typename InputVector::value_type>>
        &divergences) const;

  private:
    /**
     * 一个指向我们操作的FEValuesBase对象的指针。
     *
     */
    const SmartPointer<const FEValuesBase<dim, spacedim>> fe_values;

    /**
     * 这个视图代表FEValuesBase对象的向量的第一个分量。
     *
     */
    const unsigned int first_tensor_component;

    /**
     * 存储有关形状函数的数据。
     *
     */
    std::vector<ShapeFunctionData> shape_function_data;
  };


  template <int rank, int dim, int spacedim = dim>
  class Tensor;

  /**
   * 一个代表对一组 <code>dim*dim</code>
   * 分量的视图的类，这些分量构成了一个来自矢量值有限元的二阶张量。视图将在
   * @ref vector_valued 模块中讨论。
   * 该类允许查询代表张量的形状函数和解决方案的（组件）的值、梯度和发散。张量
   * $T_{ij},\, 0\le i,j<\text{dim}$ 的发散被定义为 $d_i = \sum_j
   * \frac{\partial T_{ij}}{\partial x_j}, \, 0\le i<\text{dim}$
   * ，而它的梯度是 $G_{ijk} = \frac{\partial T_{ij}}{\partial x_k}$
   * 。    如果你将 FEValuesExtractors::Tensor
   * 应用于FEValues、FEFaceValues或FESubfaceValues对象，你会得到这种类型的对象。
   * @ingroup feaccess vector_valued
   *
   */
  template <int dim, int spacedim>
  class Tensor<2, dim, spacedim>
  {
  public:
    /**
     * 当你将这种提取器应用于一个矢量值的有限元时，你得到的数据类型。
     *
     */
    using value_type = dealii::Tensor<2, spacedim>;

    /**
     * 用于获取张量的发散的数据类型：一个矢量。
     *
     */
    using divergence_type = dealii::Tensor<1, spacedim>;

    /**
     * 用于获取二阶张量的梯度的数据类型：三阶张量。
     *
     */
    using gradient_type = dealii::Tensor<3, spacedim>;

    /**
     * @p Number
     * 与该类提供的视图值的乘积的数据类型的别名。这是一个有限元场的矢量分量的数据类型，其自由度由一个元素类型为
     * @p Number. 的矢量描述。
     *
     */
    template <typename Number>
    using solution_value_type = typename ProductType<Number, value_type>::type;

    /**
     * 一个 @p Number
     * 和该类提供的视图的分歧的乘积的数据类型的别名。这是一个有限元场的矢量分量的数据类型，其自由度由元素类型为
     * @p Number. 的矢量描述。
     *
     */
    template <typename Number>
    using solution_divergence_type =
      typename ProductType<Number, divergence_type>::type;

    /**
     * 该类提供的 @p Number
     * 与视图梯度的乘积的数据类型的别名。这是一个有限元场的矢量分量的数据类型，其自由度由元素类型为
     * @p Number. 的矢量描述。
     *
     */
    template <typename Number>
    using solution_gradient_type =
      typename ProductType<Number, gradient_type>::type;


    /**
     * 一个结构，为张量视图和任何 @p Number
     * 类型的基函数的值和导数的乘积提供输出类型。
     * @deprecated  使用周围类中定义的类型来代替。
     *
     */
    template <typename Number>
    struct DEAL_II_DEPRECATED OutputType
    {
      /**
       * 一个 @p Number
       * 和张量类的视图值的乘积的数据类型的别名。
       *
       */
      using value_type = typename ProductType<
        Number,
        typename Tensor<2, dim, spacedim>::value_type>::type;

      /**
       * 一个 @p Number
       * 和张量类视图的分歧的乘积的数据类型的别名。
       *
       */
      using divergence_type = typename ProductType<
        Number,
        typename Tensor<2, dim, spacedim>::divergence_type>::type;

      /**
       * 张量类中 @p Number 与梯度的乘积的数据类型的别名。
       *
       */
      using gradient_type = typename ProductType<
        Number,
        typename Tensor<2, dim, spacedim>::gradient_type>::type;
    };

    /**
     * 一个结构，对于每个形状函数，我们预先计算出一堆数据，这将使以后的访问变得更便宜。
     *
     */
    struct ShapeFunctionData
    {
      /**
       * 对于每一对(形状函数,向量内的分量)，存储所选向量分量是否可能为非零。对于原始形状函数，我们肯定知道某个给定形状函数的某个标量分量是否为非零，而对于非原始形状函数，这可能并不完全清楚（例如，对于RT元素，它取决于单元格的形状）。
       *
       */
      bool is_nonzero_shape_function_component
        [value_type::n_independent_components];

      /**
       * 对于每一对（形状函数，矢量内的成分），在shape_values、shape_gradients和shape_hessians表中存储行索引（列索引是正交点索引）。如果形状函数是原始的，那么我们可以从FEValues对象的shape_function_to_row_table中获得这些信息；否则，我们必须花点功夫来计算这些信息。
       *
       */
      unsigned int row_index[value_type::n_independent_components];

      /**
       * 对于每个形状函数说如下：如果这个形状函数的is_nonzero_shape_function_component中只有一个条目是非零的，那么就存储row_index的相应值，single_nonzero_component_index代表在0和（dim^2）之间的索引，对它来说是达到了。如果多个分量为非零，则存储
       *
       * 如果没有组件是非零的，则存储* - 。
       *
       * - .
       *
       */
      int single_nonzero_component;

      /**
       * @p single_nonzero_component 的索引。
       *
       */
      unsigned int single_nonzero_component_index;
    };

    /**
     * 默认构造函数。创建一个无效的对象。
     *
     */
    Tensor();

    /**
     * 复制构造函数。这不是一个轻量级的对象，所以我们不允许复制，如果调用这个函数会产生一个编译时错误。
     *
     */
    Tensor(const Tensor<2, dim, spacedim> &) = delete;

    /**
     * 移动构造函数。
     *
     */
    // NOLINTNEXTLINE OSX does not compile with noexcept
    Tensor(Tensor<2, dim, spacedim> &&) = default;

    /**
     * 解构器。
     *
     */
    ~Tensor() = default;

    /**
     * 表示FEValuesBase对象（或从FEValuesBase派生的类之一）的
     * <code>(dim*dim)</code>
     * 分量的对象的构造函数，表示构成二阶张量值变量的唯一分量。
     * 第二个参数表示所选对称二阶张量的第一个分量的索引。
     *
     */
    Tensor(const FEValuesBase<dim, spacedim> &fe_values_base,
           const unsigned int                 first_tensor_component);


    /**
     * 复制操作符。这不是一个轻量级的对象，所以我们不允许复制，如果调用这个函数会产生一个编译时错误。
     *
     */
    Tensor &
    operator=(const Tensor<2, dim, spacedim> &) = delete;

    /**
     * 移动赋值运算符。
     *
     */
    // NOLINTNEXTLINE
    Tensor &operator=(Tensor<2, dim, spacedim> &&) = default;

    /**
     * 返回该视图所选择的向量分量的值，用于参数所选择的形状函数和正交点。
     * 这里，由于视图代表了FEValues对象中具有
     * <code>(dim*dim)</code>
     * 分量（二阶张量的唯一分量）的矢量值部分，所以返回类型是等级2的张量。
     * @param  shape_function 要评估的形状函数的编号。
     * 注意，这个数字从零到dofs_per_cell，即使是在FEFaceValues或FESubfaceValues对象的情况下。
     * @param  q_point 要评估函数的正交点的编号。
     * @dealiiRequiresUpdateFlags{update_values}
     *
     */
    value_type
    value(const unsigned int shape_function, const unsigned int q_point) const;

    /**
     * 返回该视图所选择的矢量分量的矢量发散，对于参数所选择的形状函数和正交点。
     * 关于发散的定义，请参见本类的一般讨论。
     * @note  参数的含义与value()函数的记录相同。
     * @dealiiRequiresUpdateFlags{update_gradients}
     *
     */
    divergence_type
    divergence(const unsigned int shape_function,
               const unsigned int q_point) const;

    /**
     * 返回该视图选择的向量分量的梯度（3阶张量），用于参数选择的形状函数和正交点。
     * 关于梯度的定义，请参见本类的一般讨论。
     * @note  参数的含义与value()函数的记录相同。
     * @dealiiRequiresUpdateFlags{update_gradients}
     *
     */
    gradient_type
    gradient(const unsigned int shape_function,
             const unsigned int q_point) const;

    /**
     * 返回上次调用FEValues对象的<tt>reinit</tt>函数时选择的单元、面或子面的正交点上以<tt>fe_function</tt>为特征的有限元函数的选定向量分量的值。
     * 这个函数相当于 FEValuesBase::get_function_values
     * 函数，但它只对选定的向量分量起作用。
     * 输出向量所存储的数据类型必须是你将形状函数的值相乘后得到的（即
     * @p value_type) 乘以用于存储你的有限元向量 $U$ 的未知数
     * $U_j$ 的值的类型（由 @p fe_function 参数代表）。
     * @dealiiRequiresUpdateFlags{update_values}
     *
     */
    template <class InputVector>
    void
    get_function_values(
      const InputVector &fe_function,
      std::vector<solution_value_type<typename InputVector::value_type>>
        &values) const;

    /**
     * 与上述相同，但使用局部自由度值的向量。换句话说，这个函数不是从与DoFHandler对象相关的全局向量中提取位于当前单元上的自由度的节点值（如上面的函数），而是通过其第一个参数获取这些局部节点值。获得这样一个向量的典型方法是调用如下代码
     * @code
     * cell->get_dof_values (dof_values, local_dof_values);
     * @endcode
     * 当前函数的意义在于，人们可以先修改这些局部值，例如应用限制器或确保所有节点值为正，然后再评估当前单元上与这些局部值相对应的有限元场。另一种应用是，人们希望将一个单元上的解后处理为每个单元上的不同的有限元空间，而不需要实际创建一个相应的DoFHandler
     *
     * 在这种情况下，我们所要计算的是该后处理函数的局部表示，其特征是节点值；然后该函数允许在正交点评估该表示。
     * @param[in]  dof_values
     * 一个本地节点值的向量。该向量的长度必须等于当前单元上的DoF数量，并且必须按照参考单元上自由度编号的相同顺序排序。
     * @param[out]  值
     * 给定的有限元场的值的向量，在当前对象上的正交点。
     * @tparam  InputVector  @p InputVector
     * 类型必须允许从中创建ArrayView对象；这一点由
     * `std::vector` 类和其他类满足。
     *
     */
    template <class InputVector>
    void
    get_function_values_from_local_dof_values(
      const InputVector &dof_values,
      std::vector<solution_value_type<typename InputVector::value_type>>
        &values) const;

    /**
     * 返回上次调用FEValues对象的<tt>reinit</tt>函数时选择的单元、面或子面的正交点上由<tt>fe_function</tt>表征的有限元函数的选定向量分量的发散。
     * 在FEValues类中没有诸如 FEValuesBase::get_function_divergences
     * 这样的等效函数，但当然可以从
     * FEValuesBase::get_function_gradients, 中获取信息。
     * 关于发散的定义，请参见该类的一般讨论。
     * 输出向量存储的数据类型必须是你将形状函数的发散量相乘后得到的数据（即
     * @p divergence_type) 乘以用于存储你的有限元向量 $U$
     * 的未知数 $U_j$ 值的类型（由 @p fe_function 参数代表）。
     * @dealiiRequiresUpdateFlags{update_gradients}
     *
     */
    template <class InputVector>
    void
    get_function_divergences(
      const InputVector &fe_function,
      std::vector<solution_divergence_type<typename InputVector::value_type>>
        &divergences) const;

    /**
     * 这个函数与get_function_divergences()的关系，与get_function_values_from_local_dof_values()与get_function_values()的关系相同。参见get_function_values_from_local_dof_values()的文档以了解更多信息。
     *
     */
    template <class InputVector>
    void
    get_function_divergences_from_local_dof_values(
      const InputVector &dof_values,
      std::vector<solution_divergence_type<typename InputVector::value_type>>
        &divergences) const;

    /**
     * 返回上次调用FEValues对象的<tt>reinit</tt>函数时选择的单元、面或子面的正交点上由<tt>fe_function</tt>表征的有限元函数的梯度。
     * 关于梯度的定义，请参见本类的一般讨论。
     * 输出向量所存储的数据类型必须是你将形状函数的梯度相乘后得到的（即
     * @p gradient_type) 乘以用于存储你的有限元向量 $U_j$
     * 的未知数值的类型 $U$ （由 @p fe_function 参数代表）。
     * @dealiiRequiresUpdateFlags{update_gradients}
     *
     */
    template <class InputVector>
    void
    get_function_gradients(
      const InputVector &fe_function,
      std::vector<solution_gradient_type<typename InputVector::value_type>>
        &gradients) const;

    /**
     * 这个函数与get_function_gradients()的关系，与get_function_values_from_local_dof_values()与get_function_values()的关系相同。参见get_function_values_from_local_dof_values()的文档以了解更多信息。
     *
     */
    template <class InputVector>
    void
    get_function_gradients_from_local_dof_values(
      const InputVector &dof_values,
      std::vector<solution_gradient_type<typename InputVector::value_type>>
        &gradients) const;

  private:
    /**
     * 一个指向我们操作的FEValuesBase对象的指针。
     *
     */
    const SmartPointer<const FEValuesBase<dim, spacedim>> fe_values;

    /**
     * 这个视图代表FEValuesBase对象的向量的第一个分量。
     *
     */
    const unsigned int first_tensor_component;

    /**
     * 存储有关形状函数的数据。
     *
     */
    std::vector<ShapeFunctionData> shape_function_data;
  };

} // namespace FEValuesViews


namespace internal
{
  namespace FEValuesViews
  {
    /**
     * 一个类，其特化用于定义什么FEValuesViews对象与给定的FEValuesExtractors对象相对应。
     *
     */
    template <int dim, int spacedim, typename Extractor>
    struct ViewType
    {};

    /**
     * 一个类，它的特化是用来定义什么FEValuesViews对象对应于给定的FEValuesExtractors对象。
     * 当使用 FEValuesExtractors::Scalar, 时，对应的视图是一个
     * FEValuesViews::Scalar<dim, spacedim>。
     *
     */
    template <int dim, int spacedim>
    struct ViewType<dim, spacedim, FEValuesExtractors::Scalar>
    {
      using type = typename dealii::FEValuesViews::Scalar<dim, spacedim>;
    };

    /**
     * 一个类，其专门化用于定义什么FEValuesViews对象与给定的FEValuesExtractors对象相对应。
     * 当使用 FEValuesExtractors::Vector, 时，对应的视图是一个
     * FEValuesViews::Vector<dim, spacedim>。
     *
     */
    template <int dim, int spacedim>
    struct ViewType<dim, spacedim, FEValuesExtractors::Vector>
    {
      using type = typename dealii::FEValuesViews::Vector<dim, spacedim>;
    };

    /**
     * 一个类，其专门化用于定义什么FEValuesViews对象与给定的FEValuesExtractors对象相对应。
     * 当使用 FEValuesExtractors::Tensor<rank>,
     * 时，对应的视图是一个 FEValuesViews::Tensor<rank, dim,
     * spacedim>。
     *
     */
    template <int dim, int spacedim, int rank>
    struct ViewType<dim, spacedim, FEValuesExtractors::Tensor<rank>>
    {
      using type = typename dealii::FEValuesViews::Tensor<rank, dim, spacedim>;
    };

    /**
     * 一个类，其专门化用于定义什么FEValuesViews对象与给定的FEValuesExtractors对象相对应。
     * 当使用 FEValuesExtractors::SymmetricTensor<rank>,
     * 时，对应的视图是一个 FEValuesViews::SymmetricTensor<rank,
     * dim, spacedim>。
     *
     */
    template <int dim, int spacedim, int rank>
    struct ViewType<dim, spacedim, FEValuesExtractors::SymmetricTensor<rank>>
    {
      using type =
        typename dealii::FEValuesViews::SymmetricTensor<rank, dim, spacedim>;
    };

    /**
     * 一个类的对象，其中存储有 FEValuesViews::Scalar,
     * FEValuesViews::Vector,
     * 等对象的集合。FEValuesBase类在构建时使用它来生成所有可能的Views类；我们在构建时这样做，因为Views类会缓存一些信息，因此创建起来相对昂贵。
     *
     */
    template <int dim, int spacedim>
    struct Cache
    {
      /**
       * 缓存标量和矢量以及对称二阶张量值的视图。
       *
       */
      std::vector<dealii::FEValuesViews::Scalar<dim, spacedim>> scalars;
      std::vector<dealii::FEValuesViews::Vector<dim, spacedim>> vectors;
      std::vector<dealii::FEValuesViews::SymmetricTensor<2, dim, spacedim>>
        symmetric_second_order_tensors;
      std::vector<dealii::FEValuesViews::Tensor<2, dim, spacedim>>
        second_order_tensors;

      /**
       * 构造函数。
       *
       */
      Cache(const FEValuesBase<dim, spacedim> &fe_values);
    };
  } // namespace FEValuesViews
} // namespace internal

namespace FEValuesViews
{
  /**
   * 一个模板化的别名，将FEValuesViews中的相应视图关联到给定的Extractor类。
   *
   */
  template <int dim, int spacedim, typename Extractor>
  using View = typename dealii::internal::FEValuesViews::
    ViewType<dim, spacedim, Extractor>::type;
} // namespace FEValuesViews


/**
 * FEValues、FEFaceValues和FESSubfaceValues对象一方面是有限元和映射类的接口，另一方面是单元和正交规则的接口。它们允许在正交公式的正交点评估形状函数的值或导数，当通过映射从单元格投射到实空间的单元格时。这种抽象的原因是可能的优化。根据有限元和映射的类型，有些值可以在单元格上计算一次。其他的必须在每个单元上计算，但也许同时计算几个值会提供优化的方法。由于这种相互作用可能很复杂，并且取决于实际的有限元，所以不能留给应用程序的程序员。
 * FEValues、FEFaceValues和FESSubfaceValues只提供数据处理：计算留给Mapping和FiniteElement类型的对象。这些对象提供了<tt>get_*_data</tt>和<tt>fill_*_values</tt>函数，分别由<tt>FEValues*</tt>的构造器和<tt>reinit</tt>函数调用。
 * <h3>General usage</h3>
 * 通常，<tt>FEValues*</tt>的对象被用于一个三角形的所有单元（或单元的面）的积分循环中。为了充分利用优化功能，应该在循环之前构建该对象，这样就可以一劳永逸地计算不依赖于单元格位置和形状的信息（例如，这包括最常见元素的正交点的形状函数值：我们可以在单元格上评估它们，当映射到真实单元格时它们将是一样的）。然后，在所有单元的循环中，必须为每个网格单元重新初始化，以计算根据实际单元而变化的那部分信息（例如，形状函数的梯度等于单元格上的梯度
 *
 * - 可以一次性计算出来
 *
 * - 乘以单位和实际单元之间映射的雅各布矩阵，这需要为每个单元重新计算）。)
 * 一段典型的代码，把对拉普拉斯矩阵的局部贡献加起来，看起来像这样。
 *
 *
 * @code
 * FEValues values (mapping, finite_element, quadrature, flags);
 * for (const auto &cell : dof_handler.active_cell_iterators())
 * {
 *   values.reinit(cell);
 *   for (unsigned int q=0; q<quadrature.size(); ++q)
 *     for (unsigned int i=0; i<finite_element.n_dofs_per_cell(); ++i)
 *       for (unsigned int j=0; j<finite_element.n_dofs_per_cell(); ++j)
 *         A(i,j) += fe_values.shape_value(i,q)
 *                   fe_values.shape_value(j,q)
 *                   fe_values.JxW(q);
 *   ...
 * }
 * @endcode
 *
 * 这里使用的各个函数描述如下。请注意，根据设计，在FEValues对象内部使用的正交点的顺序与上述传递给FEValues对象构造函数的正交公式所定义的相同。
 * <h3>Member functions</h3>
 * 本类的函数分为不同的类别。  <ul>   <li>  shape_value(), shape_grad()等：每次返回此对象的一个值。这些函数都是内联的，所以这就是对所有有限元值的建议访问。使用优化的编译器应该不会有性能上的损失。如果有限元是矢量值的，那么这些函数返回所请求的形状函数的唯一非零分量。然而，一些有限元的形状函数有一个以上的非零分量（我们称之为非 "原始"），在这种情况下，这组函数将抛出一个异常，因为它们不能产生一个有用的结果。相反，请使用下一组函数。
 * <li>  shape_value_component(), shape_grad_component(),
 * 等等。这是与上面相同的一组函数，只是对于矢量值的有限元，它们只返回一个矢量分量。这对于那些形状函数有多个非零分量的元素来说是很有用的，因为这时不能使用上述函数，你必须使用这组函数走遍形状函数的所有（或只有非零）分量。
 * <li>  get_function_values(), get_function_gradients(), etc.:
 * 计算一个有限元函数或其在正交点的导数。
 * <li>  reinit：初始化某个单元的FEValues对象。这个函数不在本类中，只在派生类中，并且有一个变量调用语法。更多信息请参见派生类的文档。  </ul>
 *
 *  <h3>Internals about the implementation</h3>
 *
 *
 *
 * @ingroup feaccess
 *
 */
template <int dim, int spacedim>
class FEValuesBase : public Subscriptor
{
public:
  /**
   * 该对象所处的尺寸。
   *
   */
  static const unsigned int dimension = dim;

  /**
   * 该对象所处空间的尺寸。
   *
   */
  static const unsigned int space_dimension = spacedim;

  /**
   * 当前对象的正交点的数量。它的值由max_n_quadrature_points的值初始化并被更新，例如，如果为一个新的单元/面调用
   * FEFaceValues::reinit() 。
   * @note 默认值等于max_n_quadrature_points的值。
   *
   */
  const unsigned int n_quadrature_points;

  /**
   * 正交点的最大数量。这个值可能与n_quadrature_points不同，例如，如果一个具有不同面孔正交规则的QCollection被传递给初始化FEFaceValues。
   * 这对于初始化数组，分配最大的内存量是非常有用的，在以后重新调整大小的时候，可能会用到n_quadrature_points给出的当前正交点的数量。
   *
   */
  const unsigned int max_n_quadrature_points;

  /**
   * 每个单元的形状函数的数量。如果我们使用这个基类来评估单元面的有限元，这仍然是每个单元的自由度数量，而不是每个面的自由度。
   *
   */
  const unsigned int dofs_per_cell;


  /**
   * 构造函数。用<tt>n_q_points</tt>正交点，<tt>dofs_per_cell</tt>每个单元的试验函数来设置数组大小，用给定的模式在调用派生类的<tt>reinit</tt>函数时更新字段。字段本身没有被设置，这必须发生在派生类的构造函数中。
   *
   */
  FEValuesBase(const unsigned int                  n_q_points,
               const unsigned int                  dofs_per_cell,
               const UpdateFlags                   update_flags,
               const Mapping<dim, spacedim> &      mapping,
               const FiniteElement<dim, spacedim> &fe);

  /**
   * 复制赋值被删除，因为这个类的对象是不可复制的。
   *
   */
  FEValuesBase &
  operator=(const FEValuesBase &) = delete;

  /**
   * 复制构造函数被删除，因为这个类的对象是不可复制的。
   *
   */
  FEValuesBase(const FEValuesBase &) = delete;

  /**
   * 解构器。
   *
   */
  virtual ~FEValuesBase() override;


  /// @name Access to shape function values
  ///
  /// These fields are filled by the finite element.
  //@{

  /**
   * 派生类的<tt>reinit</tt>函数最后一次被调用时选择的单元格、面或子面的正交点上的形状函数的值。
   * 如果形状函数是矢量值，那么这将返回唯一的非零分量。如果形状函数有一个以上的非零分量（也就是说，它不是原始的），那么抛出一个ExcShapeFunctionNotPrimitive类型的异常。在这种情况下，使用
   * shape_value_component() 函数。      @param  function_no
   * 要被评估的形状函数的数字。请注意，这个数字从零到dofs_per_cell，即使是在FEFaceValues或FESubfaceValues对象的情况下。
   * @param  point_no 要评估函数的正交点的数目
   * @dealiiRequiresUpdateFlags{update_values}  。
   *
   */
  const double &
  shape_value(const unsigned int function_no,
              const unsigned int point_no) const;

  /**
   * 计算一个正交点的形状函数值的一个向量分量。如果有限元是标量的，那么只允许零分量，返回值等于shape_value()函数的值。如果有限元是矢量值的，但所有形状函数都是原始的（即它们只有一个分量是非零的），那么shape_value()返回的值正好等于这个函数的一个分量。因此，只有在形状函数不是基元的情况下，这个函数才更有意义，但此时它是必要的，因为其他函数不能被使用。
   * @param  function_no 要评估的形状函数的编号。      @param
   * point_no 要评估函数的正交点的编号。      @param  component
   * 要评估的向量分量。
   * @dealiiRequiresUpdateFlags{update_values}
   *
   */
  double
  shape_value_component(const unsigned int function_no,
                        const unsigned int point_no,
                        const unsigned int component) const;

  /**
   * 计算<tt>function_no</tt>第1个形状函数在<tt>quadrature_point</tt>第1个正交点的梯度，相对于实际单元坐标。
   * 如果你想得到其中一个坐标方向的导数，使用张量类的适当函数来提取这个函数返回的张量的一个分量。因为只返回对梯度值的引用，所以应该不会有大的性能缺陷。
   * 如果形状函数是矢量值的，那么它返回唯一的非零分量。如果形状函数有一个以上的非零分量（即它不是原始的），那么它将抛出一个ExcShapeFunctionNotPrimitive类型的异常。在这种情况下，请使用
   * shape_grad_component() 函数。
   * 这个函数的参数与shape_value()函数的参数相同。      @param
   * function_no 要被评估的形状函数的编号。      @param
   * quadrature_point 要评估函数的正交点的数量。
   * @dealiiRequiresUpdateFlags{update_gradients}
   *
   */
  const Tensor<1, spacedim> &
  shape_grad(const unsigned int function_no,
             const unsigned int quadrature_point) const;

  /**
   * 返回形状函数在正交点的梯度的一个向量分量。如果有限元是标量的，那么只允许零分量，返回值等于shape_grad()函数的值。如果有限元是矢量值的，但所有形状函数都是原始的（即它们只有一个分量是非零的），那么shape_grad()返回的值就等于这个函数的一个分量。因此，只有在形状函数不是基元的情况下，这个函数才更有意义，但此时它是必要的，因为其他函数不能使用。
   * 这个函数的参数与shape_value_component()函数的参数同样成立。
   * @dealiiRequiresUpdateFlags{update_gradients}
   *
   */
  Tensor<1, spacedim>
  shape_grad_component(const unsigned int function_no,
                       const unsigned int point_no,
                       const unsigned int component) const;

  /**
   * <tt>function_no</tt>第1个形状函数在<tt>point_no</tt>第1个正交点相对于实际单元坐标的二次导数。如果你想得到其中一个坐标方向的导数，使用张量类的适当函数来提取一个分量。由于只返回对Hessian值的引用，应该不会有大的性能缺陷。
   * 如果形状函数是矢量值的，那么这将返回唯一的非零分量。如果形状函数有一个以上的非零分量（即它不是原始的），那么抛出一个ExcShapeFunctionNotPrimitive类型的异常。在这种情况下，使用shape_hessian_component()函数。
   * 这个函数的参数与shape_value()函数的参数相同。
   * @dealiiRequiresUpdateFlags{update_hessians}
   *
   */
  const Tensor<2, spacedim> &
  shape_hessian(const unsigned int function_no,
                const unsigned int point_no) const;

  /**
   * 返回一个正交点上的形状函数的 hessian
   * 的一个向量分量。如果有限元是标量的，那么只允许零分量，返回值等于shape_hessian()函数的值。如果有限元是矢量的，但是所有的形状函数都是原始的（即它们只有一个分量是不为零的），那么shape_hessian()返回的值就等于这个函数的一个分量。因此，只有在形状函数不是原始函数的情况下，这个函数才更有意义，但此时它是必要的，因为其他函数不能使用。
   * 这个函数的参数与shape_value_component()函数的参数同样成立。
   * @dealiiRequiresUpdateFlags{update_hessians}
   *
   */
  Tensor<2, spacedim>
  shape_hessian_component(const unsigned int function_no,
                          const unsigned int point_no,
                          const unsigned int component) const;

  /**
   * <tt>function_no</tt>第1个形状函数在<tt>point_no</tt>第1个正交点相对于实际单元坐标的3次导数。如果你想得到其中一个坐标方向的3阶导数，请使用张量类的适当函数来提取一个分量。由于只返回对三阶导数值的引用，应该不会有大的性能缺陷。
   * 如果形状函数是矢量值的，那么这将返回唯一的非零分量。如果形状函数有一个以上的非零分量（即它不是原始的），那么抛出一个ExcShapeFunctionNotPrimitive类型的异常。在这种情况下，使用
   * shape_3rdderivative_component() 函数。
   * 这个函数的参数与shape_value()函数的参数相同。
   * @dealiiRequiresUpdateFlags{update_3rd_derivatives}
   *
   */
  const Tensor<3, spacedim> &
  shape_3rd_derivative(const unsigned int function_no,
                       const unsigned int point_no) const;

  /**
   * 返回一个正交点的形状函数的三次导数的一个向量分量。如果有限元是标量的，那么只允许零分量，返回值等于shape_3rdderivative()函数的值。如果有限元是矢量值的，但所有形状函数都是原始的（即它们只有一个分量是非零的），那么shape_3rdderivative()返回的值就等于这个函数的一个分量。因此，只有在形状函数不是原始函数的情况下，这个函数才更有意义，但此时它是必要的，因为其他函数不能使用。
   * 这个函数的参数与 shape_value_component() 函数的参数相同。
   * @dealiiRequiresUpdateFlags{update_3rd_derivatives}
   *
   */
  Tensor<3, spacedim>
  shape_3rd_derivative_component(const unsigned int function_no,
                                 const unsigned int point_no,
                                 const unsigned int component) const;

  //@}
  /// @name Access to values of global finite element fields
  //@{

  /**
   * 返回限制在当前单元、面或子面的有限元函数的值，这些单元、面或子面是在最后一次调用派生类的<tt>reinit</tt>函数时选择的正交点。
   * 如果当前的单元没有被激活，那么数值将被内插到当前的单元，并从中计算出点值。
   * 这个函数只能在使用的有限元是标量的情况下使用，即只有一个矢量分量。
   * 为了获得多分量元素的值，下面有另一个get_function_values()，返回结果的向量的向量。
   * @param[in]  fe_function
   * 一个值向量，描述（全局）该函数应在当前单元的正交点评估的有限元函数。
   * @param[out]  values
   * 由fe_function指定的函数在当前单元的正交点的值。
   * 假设该对象已经具有正确的大小。这个输出向量所存储的数据类型必须是当你将形状函数的值乘以用于存储你的有限元向量
   * $U$ （由 @p fe_function
   * 参数表示）的未知数的值的类型时得到的。这恰好等于解向量元素的类型。
   * @post   <code>values[q]</code>  将包含fe_function描述的场在 $q$
   * 第三个正交点的值。
   * @note
   * 输入矢量的实际数据类型可以是Vector&lt;T&gt;、BlockVector&lt;T&gt;或PETSc或Trilinos矢量包装类之一。它代表了与DoFHandler对象相关的全局DoF值的向量，这个FEValues对象最后被初始化。
   * @dealiiRequiresUpdateFlags{update_values}
   *
   */
  template <class InputVector>
  void
  get_function_values(
    const InputVector &                            fe_function,
    std::vector<typename InputVector::value_type> &values) const;

  /**
   * 这个函数与其他get_function_values()的作用相同，但应用于多分量（矢量值）元素。参数的含义如那里所解释。
   * @post   <code>values[q]</code> 是fe_function描述的场在 $q$
   * 个正交点的值的向量。由 <code>values[q]</code>
   * 访问的矢量的大小等于有限元的分量数，即
   * <code>values[q](c)</code> 返回 $q$ 个正交点的 $c$
   * 个矢量的值。
   *
   */
  template <class InputVector>
  void
  get_function_values(
    const InputVector &                                    fe_function,
    std::vector<Vector<typename InputVector::value_type>> &values) const;

  /**
   * 从一个任意的矢量生成函数值。这个函数与上面这个名字的第一个函数本质上是一样的，只是它没有假设输入矢量对应于描述有限元场的未知数的DoFHandler（然后我们会假设`fe_function.size()
   * ==
   * dof_handler.n_dofs()`）。相反，对应于当前单元的节点值是一个任意矢量的元素，这些元素由这个函数的第二个参数来索引。`fe_function`输入参数的其余部分对应于什么，对这个函数没有影响。
   * 鉴于此，上面的函数相当于将`fe_function`作为第一个参数传给当前函数，并将以下调用产生的`local_dof_indices`数组作为当前函数的第二个参数。
   * @code
   * cell->get_dof_indices (local_dof_indices);
   * @endcode
   * (更多信息见 DoFCellAccessor::get_dof_indices()
   * 。)同样地，上面的函数也相当于调用
   * @code
   * cell->get_dof_values (fe_function, local_dof_values);
   * @endcode
   * 然后调用当前函数，将`local_dof_values`作为第一个参数，并将一个索引为`{0,...,fe.dofs_per_cell-1}的数组作为第二个参数。
   * 当前函数的意义在于，人们有时希望在正交点评估有限元函数，其节点值没有存储在全局矢量中
   *
   * 例如，可以先修改这些局部值，例如应用限制器或确保所有节点值为正值，然后再评估当前单元上与这些局部值对应的有限元场。另一种应用是，人们希望将一个单元上的解后处理为每个单元上的不同有限元空间，而不需要实际创建一个相应的DoFHandler
   *
   * 在这种情况下，我们所要计算的是该后处理函数的局部表示，其特征是节点值；然后该函数允许在正交点评估该表示。
   * @param[in]  fe_function
   * 一个结点值的向量。这个向量可以有一个任意的大小，只要所有由
   * "indices "索引的元素可以被实际访问。      @param[in]
   * indices
   * 进入`fe_function`的索引的一个向量。这个向量的长度必须等于当前单元上的自由度数量，并且必须按照参考单元上自由度的索引顺序识别`fe_function`中的元素。
   * @param[out]  values
   * 给定的有限元场的数值向量，在当前对象上的正交点。
   * @dealiiRequiresUpdateFlags{update_values}
   *
   */
  template <class InputVector>
  void
  get_function_values(
    const InputVector &                             fe_function,
    const ArrayView<const types::global_dof_index> &indices,
    std::vector<typename InputVector::value_type> & values) const;

  /**
   * 从一个任意的矢量生成矢量函数值。
   * 这个函数与前一个函数相对应，只是针对矢量值的情况。
   * @dealiiRequiresUpdateFlags{update_values}
   *
   */
  template <class InputVector>
  void
  get_function_values(
    const InputVector &                                    fe_function,
    const ArrayView<const types::global_dof_index> &       indices,
    std::vector<Vector<typename InputVector::value_type>> &values) const;


  /**
   * 从一个任意的向量生成向量函数值。这个函数与前一个函数类似，但`indices`向量也可以是每个单元格的dofs数量的倍数。然后，<tt>value</tt>中的向量应该允许有限元的分量的相同倍数。
   * 根据最后一个参数的值，<tt>values</tt>的外向量要么有正交规则的长度（<tt>quadrature_points_fastest
   * ==
   * false</tt>），要么有要填充的元件的长度<tt>quadrature_points_fastest
   * ==
   * true</tt>。如果<tt>p</tt>是当前的正交点编号，<tt>i</tt>是所需解决方案的矢量分量，如果<tt>quadrature_points_fastest
   * ==
   * false</tt>，对<tt>values[p][i]</tt>的访问是<tt>values[i][p]</tt>，否则是<tt>values</tt>。
   * 由于这个函数允许相当普遍的参数大小组合，所以要注意对参数的检查可能不会发现错误。
   * @dealiiRequiresUpdateFlags{update_values}
   *
   */
  template <class InputVector>
  void
  get_function_values(
    const InputVector &                                      fe_function,
    const ArrayView<const types::global_dof_index> &         indices,
    ArrayView<std::vector<typename InputVector::value_type>> values,
    const bool quadrature_points_fastest) const;

  //@}
  /// @name Access to derivatives of global finite element fields
  //@{

  /**
   * 计算一个单元的正交点上的有限元梯度。这个函数等同于相应的get_function_values()函数（更多信息见那里），但评估的是有限元场的梯度而不是它的值。
   * 这个函数只能在使用的有限元是标量的情况下使用，即只有一个矢量成分。对于矢量值的有限元，有一个相同名称的相应函数。
   * @param[in]  fe_function
   * 一个值的向量，描述（全局）该函数应在当前单元的正交点评估的有限元函数。
   * @param[out]  gradients
   * 由fe_function指定的函数在当前单元的正交点的梯度。
   * 梯度是在实空间计算的（而不是在单元格上）。
   * 假设该对象已经有了正确的尺寸。这个输出向量所存储的数据类型必须是当你将形状函数的梯度乘以用于存储有限元向量
   * $U_j$ 的未知数 $U$ （由 @p
   * fe_function参数表示）的类型时所得到的。      @post
   * <code>gradients[q]</code>  将包含fe_function描述的场在 $q$
   * 第三个正交点的梯度。    <code>gradients[q][d]</code>
   * 代表坐标方向 $d$ 在正交点 $q$ 的导数。
   * @note
   * 输入矢量的实际数据类型可以是Vector&lt;T&gt;、BlockVector&lt;T&gt;或PETSc或Trilinos矢量包装类之一。它代表了与DoFHandler对象相关的全局DoF值的向量，而这个FEValues对象最后被初始化。
   * @dealiiRequiresUpdateFlags{update_gradients}
   *
   */
  template <class InputVector>
  void
  get_function_gradients(
    const InputVector &fe_function,
    std::vector<Tensor<1, spacedim, typename InputVector::value_type>>
      &gradients) const;

  /**
   * 这个函数与其他get_function_gradients()的作用相同，但应用于多分量（矢量值）元素。参数的含义与那里解释的一样。
   * @post   <code>gradients[q]</code>  是fe_function描述的场在 $q$
   * 第1个正交点的梯度的矢量。由 <code>gradients[q]</code>
   * 访问的矢量的大小等于有限元的分量数，即
   * <code>gradients[q][c]</code> 返回 $c$ 个矢量分量在 $q$
   * 个正交点的梯度。因此， <code>gradients[q][c][d]</code>
   * 是当前单元的 $c$ 矢量场的第1个矢量分量在坐标方向 $d$
   * 上的导数。      @dealiiRequiresUpdateFlags{update_gradients}
   *
   */
  template <class InputVector>
  void
  get_function_gradients(
    const InputVector &fe_function,
    std::vector<
      std::vector<Tensor<1, spacedim, typename InputVector::value_type>>>
      &gradients) const;

  /**
   * 这个函数与上述get_function_gradients()函数中的第一个函数的关系，与带有类似参数的get_function_values()函数与get_function_values()函数中的第一个函数的关系相同。更多信息见那里。
   * @dealiiRequiresUpdateFlags{update_gradients}
   *
   */
  template <class InputVector>
  void
  get_function_gradients(
    const InputVector &                             fe_function,
    const ArrayView<const types::global_dof_index> &indices,
    std::vector<Tensor<1, spacedim, typename InputVector::value_type>>
      &gradients) const;

  /**
   * 这个函数与上述get_function_gradients()函数中的第一个函数的关系，与带有类似参数的get_function_values()函数与get_function_values()函数中的第一个函数的关系相同。更多信息见那里。
   * @dealiiRequiresUpdateFlags{update_gradients}
   *
   */
  template <class InputVector>
  void
  get_function_gradients(
    const InputVector &                             fe_function,
    const ArrayView<const types::global_dof_index> &indices,
    ArrayView<
      std::vector<Tensor<1, spacedim, typename InputVector::value_type>>>
               gradients,
    const bool quadrature_points_fastest = false) const;

  //@}
  /// @name Access to second derivatives
  ///
  /// Hessian matrices and Laplacians of global finite element fields
  //@{

  /**
   * 计算单元格正交点上的有限元二阶导数的张量。这个函数等同于相应的get_function_values()函数（更多信息见那里），但评估有限元场的二阶导数而不是其值。
   * 这个函数只能在使用的有限元是标量的情况下使用，即只有一个矢量成分。对于矢量值的有限元有一个相同名称的相应函数。
   * @param[in]  fe_function
   * 一个值的向量，描述（全局）该函数应在当前单元的正交点评估的有限元函数。
   * @param[out]  hessians
   * 由fe_function指定的函数在当前单元的正交点上的Hessians。
   * Hessians是在实空间计算的（而不是在单元格上）。
   * 假设该对象已经有了正确的尺寸。这个输出向量所存储的数据类型必须是你将形状函数的Hessians乘以用于存储你的有限元向量
   * $U_j$ （由 @p
   * fe_function参数表示）的未知数值的类型时得到的。
   * @post   <code>hessians[q]</code>  将包含fe_function描述的场在 $q$
   * 第三个正交点的Hessian。    <code>hessians[q][i][j]</code>
   * 代表正交点 $q$ 的第二导数矩阵的 $(i,j)$ 分量。
   * @note
   * 输入矢量的实际数据类型可以是Vector&lt;T&gt;、BlockVector&lt;T&gt;或PETSc或Trilinos矢量封装类之一。它代表了与DoFHandler对象相关的全局DoF值的向量，这个FEValues对象最后被初始化。
   * @dealiiRequiresUpdateFlags{update_hessians}
   *
   */
  template <class InputVector>
  void
  get_function_hessians(
    const InputVector &fe_function,
    std::vector<Tensor<2, spacedim, typename InputVector::value_type>>
      &hessians) const;

  /**
   * 这个函数与其他get_function_hessians()的作用相同，但应用于多分量（矢量值）元素。参数的含义与那里解释的一样。
   * @post   <code>hessians[q]</code> 是fe_function描述的场在 $q$
   * 第1个正交点的Hessians的向量。由 <code>hessians[q]</code>
   * 访问的矢量的大小等于有限元的分量数，即
   * <code>hessians[q][c]</code> 返回 $c$ 第1个正交点的 $q$
   * 矢量分量的Hessian。因此， <code>hessians[q][c][i][j]</code>
   * 是当前单元格的正交点 $q$ 处向量场的 $c$
   * 个向量分量的第二导数矩阵的 $(i,j)$ 个分量。
   * @dealiiRequiresUpdateFlags{update_hessians}
   *
   */
  template <class InputVector>
  void
  get_function_hessians(
    const InputVector &fe_function,
    std::vector<
      std::vector<Tensor<2, spacedim, typename InputVector::value_type>>>
      &        hessians,
    const bool quadrature_points_fastest = false) const;

  /**
   * 这个函数与上述get_function_hessians()函数中的第一个函数的关系，与带有类似参数的get_function_values()函数与get_function_values()函数中的第一个函数的关系相同。更多信息见那里。
   * @dealiiRequiresUpdateFlags{update_hessians}
   *
   */
  template <class InputVector>
  void
  get_function_hessians(
    const InputVector &                             fe_function,
    const ArrayView<const types::global_dof_index> &indices,
    std::vector<Tensor<2, spacedim, typename InputVector::value_type>>
      &hessians) const;

  /**
   * 这个函数与上述get_function_hessians()函数中的第一个函数的关系，与带类似参数的get_function_values()函数与get_function_values()函数中的第一个函数的关系相同。更多信息见那里。
   * @dealiiRequiresUpdateFlags{update_hessians}
   *
   */
  template <class InputVector>
  void
  get_function_hessians(
    const InputVector &                             fe_function,
    const ArrayView<const types::global_dof_index> &indices,
    ArrayView<
      std::vector<Tensor<2, spacedim, typename InputVector::value_type>>>
               hessians,
    const bool quadrature_points_fastest = false) const;

  /**
   * 计算单元格正交点上的有限元的（标量）拉普拉斯（即二阶导数张量的迹）。这个函数等同于相应的get_function_values()函数（更多信息见那里），但评估的是有限元场的二阶导数而不是其值。
   * 这个函数只能在使用的有限元是标量的情况下使用，即只有一个矢量成分。对于矢量值的有限元有一个相同名称的相应函数。
   * @param[in]  fe_function
   * 一个值的向量，描述（全局）该函数应在当前单元的正交点评估的有限元函数。
   * @param[out]  laplacians
   * 由fe_function指定的函数在当前单元的正交点上的拉普拉斯方程。
   * 拉普拉斯是在实空间计算的（而不是在单元格上）。
   * 假设该对象已经有了正确的尺寸。这个输出向量所存储的数据类型必须是当你将形状函数的拉普拉斯系数乘以用于存储有限元向量
   * $U_j$ 的未知数 $U$ （由 @p
   * fe_function参数表示）的类型时所得到的。这恰好等于输入向量元素的类型。
   * @post   <code>laplacians[q]</code>  将包含fe_function描述的场在
   * $q$  第三个正交点的拉普拉斯。      @post
   * 对于输出向量的每个分量，都有
   * <code>laplacians[q]=trace(hessians[q])</code>
   * ，其中<tt>hessians</tt>将是get_function_hessians()函数的输出。
   * @note
   * 输入向量的实际数据类型可以是Vector&lt;T&gt;、BlockVector&lt;T&gt;或PETSc或Trilinos向量封装类之一。它代表了与DoFHandler对象相关的全局DoF值的向量，而这个FEValues对象最后被初始化。
   * @dealiiRequiresUpdateFlags{update_hessians}
   *
   */
  template <class InputVector>
  void
  get_function_laplacians(
    const InputVector &                            fe_function,
    std::vector<typename InputVector::value_type> &laplacians) const;

  /**
   * 这个函数与其他get_function_laplacians()的作用相同，但应用于多分量（矢量值）元素。参数的含义与那里解释的一样。
   * @post   <code>laplacians[q]</code>  是fe_function描述的场在 $q$
   * 个正交点的拉普拉斯矢量。由 <code>laplacians[q]</code>
   * 访问的矢量的大小等于有限元的分量数，即
   * <code>laplacians[q][c]</code> 返回 $c$ 第1个正交点的 $q$
   * 矢量的拉普拉斯。      @post
   * 对于输出向量的每个分量，持有
   * <code>laplacians[q][c]=trace(hessians[q][c])</code>
   * ，其中<tt>hessians</tt>将是get_function_hessians()
   * 函数的输出。      @dealiiRequiresUpdateFlags{update_hessians}
   *
   */
  template <class InputVector>
  void
  get_function_laplacians(
    const InputVector &                                    fe_function,
    std::vector<Vector<typename InputVector::value_type>> &laplacians) const;

  /**
   * 这个函数与上述get_function_laplacians()函数中的第一个函数的关系，与带有类似参数的get_function_values()函数与get_function_values()函数中的第一个函数的关系相同。更多信息见那里。
   * @dealiiRequiresUpdateFlags{update_hessians}
   *
   */
  template <class InputVector>
  void
  get_function_laplacians(
    const InputVector &                             fe_function,
    const ArrayView<const types::global_dof_index> &indices,
    std::vector<typename InputVector::value_type> & laplacians) const;

  /**
   * 这个函数与上述get_function_laplacians()函数中的第一个函数的关系，与带类似参数的get_function_values()函数与get_function_values()函数中的第一个函数的关系相同。更多信息见那里。
   * @dealiiRequiresUpdateFlags{update_hessians}
   *
   */
  template <class InputVector>
  void
  get_function_laplacians(
    const InputVector &                                    fe_function,
    const ArrayView<const types::global_dof_index> &       indices,
    std::vector<Vector<typename InputVector::value_type>> &laplacians) const;

  /**
   * 这个函数与上述get_function_laplacians()函数中的第一个函数的关系，与带类似参数的get_function_values()函数与get_function_values()函数中的第一个函数的关系相同。更多信息见那里。
   * @dealiiRequiresUpdateFlags{update_hessians}
   *
   */
  template <class InputVector>
  void
  get_function_laplacians(
    const InputVector &                                         fe_function,
    const ArrayView<const types::global_dof_index> &            indices,
    std::vector<std::vector<typename InputVector::value_type>> &laplacians,
    const bool quadrature_points_fastest = false) const;

  //@}
  /// @name Access to third derivatives of global finite element fields
  //@{

  /**
   * 计算一个单元的正交点上的有限元的三阶导数张量。这个函数等同于相应的get_function_values()函数（更多信息见那里），但评估的是有限元场的第三导数而不是它的值。
   * 这个函数只能在使用的有限元是标量的情况下使用，即只有一个矢量成分。对于矢量值的有限元有一个相同名称的相应函数。
   * @param[in]  fe_function
   * 一个值的向量，描述（全局）该函数应在当前单元的正交点评估的有限元函数。
   * @param[out]  third_derivatives
   * 由fe_function指定的函数在当前单元的正交点上的第三导数。
   * 三次导数是在实空间计算的（而不是在单元格上）。
   * 假设该对象已经有了正确的尺寸。这个输出向量所存储的数据类型必须是当你将形状函数的三次导数乘以用于存储你的有限元向量
   * $U_j$ 的未知数值的类型 $U$ （由 @p fe_function
   * 参数表示）时所得到的。      @post
   * <code>third_derivatives[q]</code> 将包含fe_function描述的场在 $q$
   * 第1个正交点的第三导数。
   * <code>third_derivatives[q][i][j][k]</code> 代表在正交点 $q$
   * 的三阶张量的 $(i,j,k)$ 的第三导数的第三部分。
   * @note
   * 输入矢量的实际数据类型可以是Vector&lt;T&gt;、BlockVector&lt;T&gt;或PETSc或Trilinos矢量封装类之一。它代表了与DoFHandler对象相关的全局DoF值的向量，而这个FEValues对象最后被初始化。
   * @dealiiRequiresUpdateFlags{update_3rd_derivatives}
   *
   */
  template <class InputVector>
  void
  get_function_third_derivatives(
    const InputVector &fe_function,
    std::vector<Tensor<3, spacedim, typename InputVector::value_type>>
      &third_derivatives) const;

  /**
   * 这个函数的作用与其他get_function_third_derivatives()相同，但应用于多分量（矢量值）元素。参数的含义与那里的解释相同。
   * @post   <code>third_derivatives[q]</code>  是fe_function描述的场在
   * $q$ 第1个正交点上的第三导数的向量。由
   * <code>third_derivatives[q]</code>
   * 访问的矢量的大小等于有限元的分量数，即
   * <code>third_derivatives[q][c]</code> 返回 $c$ 在 $q$
   * 第1个正交点的第1个矢量分量的三次导数。因此，
   * <code>third_derivatives[q][c][i][j][k]</code> 是当前单元的正交点
   * $q$ 处向量场的 $c$ 第3个向量分量的张量的 $(i,j,k)$ 。
   * @dealiiRequiresUpdateFlags{update_3rd_derivatives}
   *
   */
  template <class InputVector>
  void
  get_function_third_derivatives(
    const InputVector &fe_function,
    std::vector<
      std::vector<Tensor<3, spacedim, typename InputVector::value_type>>>
      &        third_derivatives,
    const bool quadrature_points_fastest = false) const;

  /**
   * 这个函数与上述get_function_third_derivatives()函数中的第一个函数的关系，与带有类似参数的get_function_values()函数与get_function_values()函数中的第一个函数的关系相同。更多信息见那里。
   * @dealiiRequiresUpdateFlags{update_3rd_derivatives}
   *
   */
  template <class InputVector>
  void
  get_function_third_derivatives(
    const InputVector &                             fe_function,
    const ArrayView<const types::global_dof_index> &indices,
    std::vector<Tensor<3, spacedim, typename InputVector::value_type>>
      &third_derivatives) const;

  /**
   * 这个函数与上述get_function_third_derivatives()函数中的第一个函数的关系，与带有类似参数的get_function_values()函数与get_function_values()函数中的第一个函数的关系相同。更多信息见那里。
   * @dealiiRequiresUpdateFlags{update_3rd_derivatives}
   *
   */
  template <class InputVector>
  void
  get_function_third_derivatives(
    const InputVector &                             fe_function,
    const ArrayView<const types::global_dof_index> &indices,
    ArrayView<
      std::vector<Tensor<3, spacedim, typename InputVector::value_type>>>
               third_derivatives,
    const bool quadrature_points_fastest = false) const;
  //@}

  /// @name Cell degrees of freedom
  //@{

  /**
   * 返回一个对象，它可以被认为是一个包含从零（包括）到`dofs_per_cell`（不包括）所有索引的数组。这允许人们使用基于范围的
   * "for "循环来编写代码，如以下类型。
   * @code
   * FEValues<dim>      fe_values (...);
   * FullMatrix<double> cell_matrix (...);
   *
   * for (auto &cell : dof_handler.active_cell_iterators())
   *   {
   *     cell_matrix = 0;
   *     fe_values.reinit(cell);
   *     for (const auto q : fe_values.quadrature_point_indices())
   *       for (const auto i : fe_values.dof_indices())
   *         for (const auto j : fe_values.dof_indices())
   *           cell_matrix(i,j) += ...; // Do something for DoF indices (i,j)
   *                                    // at quadrature point q
   *   }
   * @endcode
   * 这里，我们在所有单元上的所有自由度上循环，`i'和`j'代表所有单元自由度的有效指数，由传递给`fe_values'的有限元定义。
   *
   */
  std_cxx20::ranges::iota_view<unsigned int, unsigned int>
  dof_indices() const;

  /**
   * 返回一个对象，可以认为是一个数组，包含从 @p
   * start_dof_index
   * （包括）到`dofs_per_cell`（包括）的所有索引。
   * 这允许人们使用基于范围的 "for
   * "循环来编写以下类型的代码。
   * @code
   * FEValues<dim>      fe_values (...);
   * FullMatrix<double> cell_matrix (...);
   *
   * for (auto &cell : dof_handler.active_cell_iterators())
   *   {
   *     cell_matrix = 0;
   *     fe_values.reinit(cell);
   *     for (const auto q : fe_values.quadrature_point_indices())
   *       for (const auto i : fe_values.dof_indices())
   *         for (const auto j : fe_values.dof_indices_starting_at(i))
   *           cell_matrix(i,j) += ...; // Do something for DoF indices (i,j)
   *                                    // at quadrature point q
   *   }
   * @endcode
   * 这里，我们在所有单元上的所有局部自由度上循环，`i`取所有单元自由度的有效指数，由传递给`fe_values`的有限元定义，而`j`取`i`范围的指定子集，从`i`本身开始，到单元自由度的数量为止。通过这种方式，我们可以构建刚度矩阵贡献的上半部分和对角线（假设它是对称的，并且只需要计算它的一半），例如。
   * @note 如果 @p start_dof_index
   * 等于单元格中的DoF数量，则返回的索引范围为空。
   *
   */
  std_cxx20::ranges::iota_view<unsigned int, unsigned int>
  dof_indices_starting_at(const unsigned int start_dof_index) const;

  /**
   * 返回一个对象，它可以被认为是一个包含从0（包括）到
   * @p end_dof_index
   * （包括）所有索引的数组。这允许人们使用基于范围的
   * "for "循环来编写以下类型的代码。
   * @code
   * FEValues<dim>      fe_values (...);
   * FullMatrix<double> cell_matrix (...);
   *
   * for (auto &cell : dof_handler.active_cell_iterators())
   *   {
   *     cell_matrix = 0;
   *     fe_values.reinit(cell);
   *     for (const auto q : fe_values.quadrature_point_indices())
   *       for (const auto i : fe_values.dof_indices())
   *         for (const auto j : fe_values.dof_indices_ending_at(i))
   *           cell_matrix(i,j) += ...; // Do something for DoF indices (i,j)
   *                                    // at quadrature point q
   *   }
   * @endcode
   * 这里，我们在所有单元上的所有局部自由度上循环，`i'代表所有单元自由度的有效指数，由传递给`fe_values'的有限元定义，而`j'代表`i'范围的指定子集，从零开始，到`i'本身结束。这样，我们可以构建刚度矩阵贡献的下半部分和对角线（假设它是对称的，并且只需要计算它的一半），例如。
   * @note 如果 @p end_dof_index
   * 等于零，则返回的索引范围为空。
   *
   */
  std_cxx20::ranges::iota_view<unsigned int, unsigned int>
  dof_indices_ending_at(const unsigned int end_dof_index) const;

  //@}

  /// @name Geometry of the cell
  //@{

  /**
   * 返回一个对象，可以认为是一个包含从零到`n_quadrature_points`所有索引的数组。这允许使用基于范围的`for'循环来编写以下类型的代码。
   * @code
   * FEValues<dim> fe_values (...);
   *
   * for (auto &cell : dof_handler.active_cell_iterators())
   *   {
   *     fe_values.reinit(cell);
   *     for (const auto q_point : fe_values.quadrature_point_indices())
   *       ... do something at the quadrature point ...
   *   }
   * @endcode
   * 这里，我们正在循环所有单元格上的所有正交点，`q_point`采用所有有效的正交点索引，由传递给`fe_values`的正交规则定义。
   * @see  CPP11
   *
   */
  std_cxx20::ranges::iota_view<unsigned int, unsigned int>
  quadrature_point_indices() const;

  /**
   * <tt>q</tt>第1个正交点在实空间的位置。
   * @dealiiRequiresUpdateFlags{update_quadrature_points}
   *
   */
  const Point<spacedim> &
  quadrature_point(const unsigned int q) const;

  /**
   * 返回实空间中正交点向量的参考。
   * @dealiiRequiresUpdateFlags{update_quadrature_points}
   *
   */
  const std::vector<Point<spacedim>> &
  get_quadrature_points() const;

  /**
   * 映射的正交点权重。如果这个对象指的是体积评价（即派生类是FEValues类型），那么这就是雅可比行列式乘以<tt>i</tt>第1个单位正交点的权重。
   * 对于表面评估（即类FEFaceValues或FESubfaceValues），它是映射的表面元素乘以正交点的权重。
   * 你可以把这个函数返回的数量看作是我们在这里通过正交实现的积分中的体积或表面元素
   * $dx, ds$ 。      @dealiiRequiresUpdateFlags{update_JxW_values}
   *
   */
  double
  JxW(const unsigned int quadrature_point) const;

  /**
   * 返回一个对存放JxW()返回值的数组的引用。
   *
   */
  const std::vector<double> &
  get_JxW_values() const;

  /**
   * 返回指定正交点上的变换的雅各布系数，即
   * $J_{ij}=dx_i/d\hat x_j$   @dealiiRequiresUpdateFlags{update_jacobians}
   * 。
   *
   */
  const DerivativeForm<1, dim, spacedim> &
  jacobian(const unsigned int quadrature_point) const;

  /**
   * 返回一个对存放jacobian()返回值的数组的引用。
   * @dealiiRequiresUpdateFlags{update_jacobians}
   *
   */
  const std::vector<DerivativeForm<1, dim, spacedim>> &
  get_jacobians() const;

  /**
   * 返回在指定的正交点，即 $G_{ijk}=dJ_{jk}/d\hat x_i$
   * ，从单位到实数单元转换的二阶导数，即雅各布式的一阶导数。
   * @dealiiRequiresUpdateFlags{update_jacobian_grads}
   *
   */
  const DerivativeForm<2, dim, spacedim> &
  jacobian_grad(const unsigned int quadrature_point) const;

  /**
   * 返回对存放jacobian_grads()返回值的数组的引用。
   * @dealiiRequiresUpdateFlags{update_jacobian_grads}
   *
   */
  const std::vector<DerivativeForm<2, dim, spacedim>> &
  get_jacobian_grads() const;

  /**
   * 返回从单位单元到实数单元转换的二阶导数，即雅各布式的一阶导数，在指定的正交点，向前推到实数单元坐标，即
   * $G_{ijk}=dJ_{iJ}/d\hat x_K (J_{jJ})^{-1} (J_{kK})^{-1}$  。
   * @dealiiRequiresUpdateFlags{update_jacobian_pushed_forward_grads}
   *
   */
  const Tensor<3, spacedim> &
  jacobian_pushed_forward_grad(const unsigned int quadrature_point) const;

  /**
   * 返回持有jacobian_pushed_forward_grads()返回值的数组的引用。
   * @dealiiRequiresUpdateFlags{update_jacobian_pushed_forward_grads}
   *
   */
  const std::vector<Tensor<3, spacedim>> &
  get_jacobian_pushed_forward_grads() const;

  /**
   * 返回在指定的正交点上，即 $G_{ijkl}=\frac{d^2J_{ij}}{d\hat x_k
   * d\hat x_l}$
   * ，从单位到实数单元转换的第三导数，即雅各布的第二导数。
   * @dealiiRequiresUpdateFlags{update_jacobian_2nd_derivatives}
   *
   */
  const DerivativeForm<3, dim, spacedim> &
  jacobian_2nd_derivative(const unsigned int quadrature_point) const;

  /**
   * 返回对存放jacobian_2nd_derivatives()返回值的数组的引用。
   * @dealiiRequiresUpdateFlags{update_jacobian_2nd_derivatives}
   *
   */
  const std::vector<DerivativeForm<3, dim, spacedim>> &
  get_jacobian_2nd_derivatives() const;

  /**
   * 返回从单位到实心单元转换的第三导数，即雅各布的第二导数，在指定的正交点，向前推到实心单元坐标，即
   * $G_{ijkl}=\frac{d^2J_{iJ}}{d\hat x_K d\hat x_L} (J_{jJ})^{-1}
   * (J_{kK})^{-1}(J_{lL})^{-1}$  。
   * @dealiiRequiresUpdateFlags{update_jacobian_pushed_forward_2nd_derivatives}
   *
   */
  const Tensor<4, spacedim> &
  jacobian_pushed_forward_2nd_derivative(
    const unsigned int quadrature_point) const;

  /**
   * 返回持有jacobian_pushed_forward_2nd_derivatives()返回值的数组的引用。
   * @dealiiRequiresUpdateFlags{update_jacobian_pushed_forward_2nd_derivatives}
   *
   */
  const std::vector<Tensor<4, spacedim>> &
  get_jacobian_pushed_forward_2nd_derivatives() const;

  /**
   * 返回在指定的正交点，即 $G_{ijklm}=\frac{d^2J_{ij}}{d\hat x_k
   * d\hat x_l d\hat x_m}$
   * ，从单位到实数单元转换的第四导数，即雅各布式的第三导数
   * 。      @dealiiRequiresUpdateFlags{update_jacobian_3rd_derivatives}
   *
   */
  const DerivativeForm<4, dim, spacedim> &
  jacobian_3rd_derivative(const unsigned int quadrature_point) const;

  /**
   * 返回对存放jacobian_3rd_derivatives()返回值的数组的引用。
   * @dealiiRequiresUpdateFlags{update_jacobian_3rd_derivatives}
   *
   */
  const std::vector<DerivativeForm<4, dim, spacedim>> &
  get_jacobian_3rd_derivatives() const;

  /**
   * 返回从单位到实心单元转换的第四导数，即雅各布的第三导数，在指定的正交点，向前推到实心单元坐标，即
   * $G_{ijklm}=\frac{d^3J_{iJ}}{d\hat x_K d\hat x_L d\hat x_M} (J_{jJ})^{-1}
   * (J_{kK})^{-1} (J_{lL})^{-1} (J_{mM})^{-1}$  。
   * @dealiiRequiresUpdateFlags{update_jacobian_pushed_forward_3rd_derivatives}
   *
   */
  const Tensor<5, spacedim> &
  jacobian_pushed_forward_3rd_derivative(
    const unsigned int quadrature_point) const;

  /**
   * 返回持有jacobian_pushed_forward_3rd_derivatives()返回值的数组的引用。
   * @dealiiRequiresUpdateFlags{update_jacobian_pushed_forward_2nd_derivatives}
   *
   */
  const std::vector<Tensor<5, spacedim>> &
  get_jacobian_pushed_forward_3rd_derivatives() const;

  /**
   * 返回指定正交点上的变换的逆雅各布系数，即
   * $J_{ij}=d\hat x_i/dx_j$
   * @dealiiRequiresUpdateFlags{update_inverse_jacobians}  。
   *
   */
  const DerivativeForm<1, spacedim, dim> &
  inverse_jacobian(const unsigned int quadrature_point) const;

  /**
   * 返回一个对存放inverse_jacobian()返回值的数组的引用。
   * @dealiiRequiresUpdateFlags{update_inverse_jacobians}
   *
   */
  const std::vector<DerivativeForm<1, spacedim, dim>> &
  get_inverse_jacobians() const;

  /**
   * 返回一个正交点的法向量。如果你为一个面调用这个函数（即，当使用FEFaceValues或FESubfaceValues对象时），那么这个函数返回面的<tt>i</tt>第1个正交点处的单元格的外向法向量。
   * 相反，如果你为一个一维的单元格调用这个函数（即，当使用`FEValues<dim,spacedim>`对象时，`spacedim>dim`），那么这个函数返回单元格的法向量
   *
   * - 换句话说，是对嵌入三角形的流形的法向量的一个近似值。在这种情况下，流形的法线方向当然有两个，这个函数返回由顶点编号引起的 "向上 "方向。    矢量的长度被规范化为1。      @dealiiRequiresUpdateFlags{update_normal_vectors}
   *
   */
  const Tensor<1, spacedim> &
  normal_vector(const unsigned int i) const;

  /**
   * 返回此对象所代表的所有正交点的法向量。关于法向量所代表的内容，请参见normal_vector()函数。
   * @dealiiRequiresUpdateFlags{update_normal_vectors}
   *
   */
  const std::vector<Tensor<1, spacedim>> &
  get_normal_vectors() const;

  //@}

  /// @name Extractors Methods to extract individual components
  //@{

  /**
   * 创建当前FEValues对象的一个视图，该视图代表可能是矢量值的有限元的特定标量分量。视图的概念在命名空间FEValuesViews的文档中解释，特别是在 @ref
   * vector_valued 模块中。
   *
   */
  const FEValuesViews::Scalar<dim, spacedim> &
  operator[](const FEValuesExtractors::Scalar &scalar) const;

  /**
   * 创建一个当前FEValues对象的视图，该视图代表了一组 <code>dim</code>
   * 矢量值有限元的标量分量（即矢量）。视图的概念在命名空间FEValuesViews的文档中解释，特别是在
   * @ref vector_valued 模块中。
   *
   */
  const FEValuesViews::Vector<dim, spacedim> &
  operator[](const FEValuesExtractors::Vector &vector) const;

  /**
   * 创建当前FEValues对象的视图，该视图代表了一组 <code>(dim*dim
   * + dim)/2</code>
   * 矢量值有限元的标量分量（即对称的二阶张量）。视图的概念在命名空间FEValuesViews的文档中解释，特别是在
   * @ref vector_valued 模块中。
   *
   */
  const FEValuesViews::SymmetricTensor<2, dim, spacedim> &
  operator[](const FEValuesExtractors::SymmetricTensor<2> &tensor) const;


  /**
   * 创建一个当前FEValues对象的视图，该视图表示一组 <code>(dim*dim)</code>
   * 矢量值有限元的标量分量（即二阶张量）。视图的概念在命名空间FEValuesViews的文档中解释，特别是在
   * @ref vector_valued 模块中。
   *
   */
  const FEValuesViews::Tensor<2, dim, spacedim> &
  operator[](const FEValuesExtractors::Tensor<2> &tensor) const;

  //@}

  /// @name Access to the raw data
  //@{

  /**
   * 对所选映射对象的常量引用。
   *
   */
  const Mapping<dim, spacedim> &
  get_mapping() const;

  /**
   * 对所选有限元对象的常数参考。
   *
   */
  const FiniteElement<dim, spacedim> &
  get_fe() const;

  /**
   * 返回为该对象设置的更新标志。
   *
   */
  UpdateFlags
  get_update_flags() const;

  /**
   * 返回到当前单元的三角形迭代器。
   *
   */
  const typename Triangulation<dim, spacedim>::cell_iterator
  get_cell() const;

  /**
   * 返回当前单元格与前一个单元格的关系。如果结果是
   * <tt>CellSimilarity::translation</tt>.
   * ，这允许重新使用一些单元格数据（如具有常数系数的方程的局部矩阵）。
   *
   */
  CellSimilarity::Similarity
  get_cell_similarity() const;

  /**
   * 确定这个对象的内存消耗（以字节为单位）的估计值。
   *
   */
  std::size_t
  memory_consumption() const;
  //@}


  /**
   * 如果要求FEValuesBase返回一个字段的值，而这个字段不是这个FEValuesBase的UpdateFlags所要求的，就会抛出这个异常。
   * @ingroup Exceptions
   *
   */
  DeclException1(
    ExcAccessToUninitializedField,
    std::string,
    << "You are requesting information from an FEValues/FEFaceValues/FESubfaceValues "
    << "object for which this kind of information has not been computed. What "
    << "information these objects compute is determined by the update_* flags you "
    << "pass to the constructor. Here, the operation you are attempting requires "
    << "the <" << arg1
    << "> flag to be set, but it was apparently not specified "
    << "upon construction.");

  /**
   * FEValues
   * FiniteElement和cell->get_dof_handler().get_fe()之间不匹配。
   * @ingroup Exceptions
   *
   */
  DeclExceptionMsg(
    ExcFEDontMatch,
    "The FiniteElement you provided to FEValues and the FiniteElement that belongs "
    "to the DoFHandler that provided the cell iterator do not match.");
  /**
   * 一个给定的形状函数不是原始的，但它需要是。
   * @ingroup Exceptions
   *
   */
  DeclException1(ExcShapeFunctionNotPrimitive,
                 int,
                 << "The shape function with index " << arg1
                 << " is not primitive, i.e. it is vector-valued and "
                 << "has more than one non-zero vector component. This "
                 << "function cannot be called for these shape functions. "
                 << "Maybe you want to use the same function with the "
                 << "_component suffix?");

  /**
   * 给定的FiniteElement不是原始元素，见
   * FiniteElement::is_primitive(). 。
   * @ingroup Exceptions
   *
   */
  DeclExceptionMsg(
    ExcFENotPrimitive,
    "The given FiniteElement is not a primitive element but the requested operation "
    "only works for those. See FiniteElement::is_primitive() for more information.");

protected:
  /**
   * FEValues类的对象需要存储一个指向当前单元的迭代器，以便能够在get_function_values()和各种函数中提取该单元上的自由度值。另一方面，这个类也应该适用于不同的迭代器，只要它们有相同的接口来提取自由度值（即，例如，它们需要有一个
   * @p get_interpolated_dof_values函数）。
   * 这就需要一个共同的迭代器类的基类，并使我们在这里需要的函数
   * @p virtual.
   * 另一方面，这是我们在库中唯一需要的地方，引入一个迭代器的基类并使一个函数虚化会惩罚
   * <em>  迭代器的所有  </em>
   * 用户，这些函数基本上是作为非常快速的存取函数。所以我们不想这样做。相反，我们在这里做的是让我们需要的函数变成虚拟的，只用于
   * <em>  这个类  </em>
   * 。这个想法是这样的：有一个共同的基类，它声明了一些纯虚拟的函数，对于每个可能的迭代器类型，我们有一个派生类，它将迭代器存储到单元格并实现这些函数。由于迭代器类具有相同的接口，我们可以使派生类成为模板，对迭代器类型进行模板化。
   * 这样一来，虚函数的使用就只限于这个类，其他迭代器的使用者就不用承担负面影响了。
   * @note  这个类是<a
   * href="https://www.artima.com/cppsource/type_erasure.html">type
   * erasure</a>设计模式的一个例子。
   *
   */
  class CellIteratorBase;

  /**
   * 向前声明源自CellIteratorBase的类。它们的定义和实现在.cc文件中给出。
   *
   */
  template <typename CI>
  class CellIterator;
  class TriaCellIterator;

  /**
   * 存储上次调用reinit()函数时选择的单元格。
   * 这对<tt>get_function_*</tt>函数以及提取器类中的同名函数是必要的。
   *
   */
  std::unique_ptr<const CellIteratorBase> present_cell;

  /**
   * 一个信号连接，我们用来确保每当三角结构因细化而发生变化时，我们会得到通知。我们需要知道这一点，因为它使所有的单元格迭代器失效，作为其中的一部分，我们在随后调用reinit()时保留了'present_cell'迭代器，以便计算单元格相似度。
   *
   */
  boost::signals2::connection tria_listener_refinement;

  /**
   * 一个信号连接，我们用它来确保每当三角结构因网格转换而发生变化时，我们都能得到通知。我们需要知道这一点，因为它使所有的单元格迭代器失效，作为其中的一部分，我们在后续调用reinit()时保留了'present_cell'迭代器，以便计算单元格的相似度。
   *
   */
  boost::signals2::connection tria_listener_mesh_transform;

  /**
   * 一个与三角结构相连的函数，以便在三角结构发生变化，迭代器随之失效时，将存储的'present_cell'迭代器重置为一个无效的迭代器。
   *
   */
  void
  invalidate_present_cell();

  /**
   * 这个函数被派生类中的各种 reinit()
   * 函数调用。给定参数所指示的单元格，测试我们是否必须丢弃之前存储的present_cell参数，因为这需要我们比较不同三角形的单元格。在检查这一切的时候，还要确保我们有tria_listener连接到我们将在调用此函数后立即设置present_cell的三角结构。
   *
   */
  void
  maybe_invalidate_previous_present_cell(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell);

  /**
   * 指向与该FEValues对象相关的映射对象的指针。
   *
   */
  const SmartPointer<const Mapping<dim, spacedim>, FEValuesBase<dim, spacedim>>
    mapping;

  /**
   * 一个指向映射内部数据对象的指针，从 Mapping::get_data(),
   * Mapping::get_face_data(), 或 Mapping::get_subface_data(). 获得。
   *
   */
  std::unique_ptr<typename Mapping<dim, spacedim>::InternalDataBase>
    mapping_data;

  /**
   * 一个对象， Mapping::fill_fe_values()
   * 和类似的函数将其输出放入其中。
   *
   */
  dealii::internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
    mapping_output;


  /**
   * 一个指向与此FEValues对象相关的有限元对象的指针。
   *
   */
  const SmartPointer<const FiniteElement<dim, spacedim>,
                     FEValuesBase<dim, spacedim>>
    fe;

  /**
   * 指向有限元内部数据对象的指针，从
   * FiniteElement::get_data(),  Mapping::get_face_data(), 或
   * FiniteElement::get_subface_data(). 获得。
   *
   */
  std::unique_ptr<typename FiniteElement<dim, spacedim>::InternalDataBase>
    fe_data;

  /**
   * 一个对象， FiniteElement::fill_fe_values()
   * 和类似的函数将其输出放入其中。
   *
   */
  dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim,
                                                                     spacedim>
    finite_element_output;


  /**
   * 交给FEValues构造器的原始更新标志。
   *
   */
  UpdateFlags update_flags;

  /**
   * 初始化一些更新标志。从派生类的 @p initialize
   * 函数中调用，这些函数又从它们的构造函数中调用。
   * 基本上，这个函数使用已经存储的有限元和映射对象找出需要设置的标志来计算用户想要的一切，正如通过作为参数传递的标志所表达的那样。
   *
   */
  UpdateFlags
  compute_update_flags(const UpdateFlags update_flags) const;

  /**
   * 一个枚举变量，可以存储当前单元与之前访问的单元的不同状态。如果需要，可以在这里检查额外的状态，并在重新启动时使用其中一个方法。
   *
   */
  CellSimilarity::Similarity cell_similarity;

  /**
   * 一个检查新单元是否与之前使用的单元相似的函数。然后，大量的数据可以被重新使用，例如实空间中的基函数导数，shape_grad。
   *
   */
  void
  check_cell_similarity(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell);

private:
  /**
   * 对所有可能的FEValuesViews对象进行缓存。
   *
   */
  dealii::internal::FEValuesViews::Cache<dim, spacedim> fe_values_views_cache;

  // Make the view classes friends of this class, since they access internal
  // data.
  template <int, int>
  friend class FEValuesViews::Scalar;
  template <int, int>
  friend class FEValuesViews::Vector;
  template <int, int, int>
  friend class FEValuesViews::SymmetricTensor;
  template <int, int, int>
  friend class FEValuesViews::Tensor;
};



/**
 * 在单元格的正交点上评估的有限元。
 * 这个函数实现了FEValuesBase的初始化程序，如果需要以单元格的正交点为单位的值。更多的文件请看这个类。
 *
 *
 * @ingroup feaccess
 *
 *
 */
template <int dim, int spacedim = dim>
class FEValues : public FEValuesBase<dim, spacedim>
{
public:
  /**
   * 我们对其进行积分的对象的尺寸。对于本类，这等于
   * <code>dim</code>  。
   *
   */
  static const unsigned int integral_dimension = dim;

  /**
   * 构造函数。从映射和有限元对象中获取单元独立数据，匹配正交规则和更新标志。
   *
   */
  FEValues(const Mapping<dim, spacedim> &      mapping,
           const FiniteElement<dim, spacedim> &fe,
           const Quadrature<dim> &             quadrature,
           const UpdateFlags                   update_flags);

  /**
   * 像上面的函数一样，但取一个正交规则的集合。
   * @note
   * 与FEFaceValues相反，我们要求集合中正交规则的数量为1。
   *
   */
  FEValues(const Mapping<dim, spacedim> &      mapping,
           const FiniteElement<dim, spacedim> &fe,
           const hp::QCollection<dim> &        quadrature,
           const UpdateFlags                   update_flags);

  /**
   * 构造函数。这个构造函数除了使对象隐含地使用 $Q_1$
   * 映射（即MappingQGeneric(1)类型的对象）外，与其他构造函数是等价的。
   *
   */
  FEValues(const FiniteElement<dim, spacedim> &fe,
           const Quadrature<dim> &             quadrature,
           const UpdateFlags                   update_flags);

  /**
   * 像上面的函数一样，但取一个正交规则的集合。
   * @note
   * 与FEFaceValues相反，我们要求集合中的正交规则的数量为1。
   *
   */
  FEValues(const FiniteElement<dim, spacedim> &fe,
           const hp::QCollection<dim> &        quadrature,
           const UpdateFlags                   update_flags);

  /**
   * 重新初始化类型为 "iterator into a DoFHandler object
   * "的给定单元的梯度、雅各比行列式等，以及与此对象相关的有限元。假设给定单元所使用的有限元也是这个FEValues对象所使用的有限元。
   *
   */
  template <bool level_dof_access>
  void
  reinit(
    const TriaIterator<DoFCellAccessor<dim, spacedim, level_dof_access>> &cell);

  /**
   * 重新初始化梯度、雅可比行列式等，用于给定类型为
   * "进入三角测量对象的迭代器
   * "的单元和给定的有限元。由于进入三角剖分的迭代器只传递单元的几何信息，而不传递可能与此单元相关的自由度的信息，所以如果需要自由度的信息，你将无法调用这一类的一些函数。这些函数首先是<tt>get_function_value/gradients/hessians/laplacians/third_derivatives</tt>函数。如果你想调用这些函数，你必须调用
   * @p
   * reinit变体，将迭代器带入DoFHandler或其他DoF处理程序类型对象。
   *
   */
  void
  reinit(const typename Triangulation<dim, spacedim>::cell_iterator &cell);

  /**
   * 返回对该对象所存储的正交公式副本的引用。
   *
   */
  const Quadrature<dim> &
  get_quadrature() const;

  /**
   * 确定此对象的内存消耗（以字节为单位）的估计值。
   *
   */
  std::size_t
  memory_consumption() const;

  /**
   * 返回对这个非常对象的引用。
   * 虽然看起来不是很有用，但是这个函数的存在是为了给
   * hp::FEValues
   * 类提供能力，在这种情况下，它提供了当前单元的FEValues对象（记住，对于hp-finite元素，实际使用的FE对象可能会在不同的单元之间变化，所以我们也需要不同的FEValues对象用于不同的单元；一旦你重新初始化了
   * hp::FEValues
   * ]对象，它就会为该单元格的FE检索FEValues对象，并通过与此相同的函数返回；因此，这里的这个函数只提供相同的接口，以便人们可以对FEValues和
   * hp::FEValues). 进行模板化。
   *
   */
  const FEValues<dim, spacedim> &
  get_present_fe_values() const;

private:
  /**
   * 在这里存储一份正交公式的副本。
   *
   */
  const Quadrature<dim> quadrature;

  /**
   * 做两个构造函数的共同工作。
   *
   */
  void
  initialize(const UpdateFlags update_flags);

  /**
   * reinit()函数只做需要了解迭代器类型的那部分工作。在设置完present_cell()后，它们会传递给这个函数，这个函数做真正的工作，而且与单元格迭代器的实际类型无关。
   *
   */
  void
  do_reinit();
};


/**
 * 将FEValuesBase的接口扩展到只有在评估单元格表面的东西时才有意义的值。所有在单元格内部的数据也都可以在这里得到。
 * 见FEValuesBase
 *
 *
 * @ingroup feaccess
 *
 *
 */
template <int dim, int spacedim = dim>
class FEFaceValuesBase : public FEValuesBase<dim, spacedim>
{
public:
  /**
   * 我们对其进行积分的对象的尺寸。对于本类，这等于
   * <code>dim-1</code>  。
   *
   */
  static const unsigned int integral_dimension = dim - 1;

  /**
   * 构造函数。调用基类的构造函数，用正确的尺寸设置本类的数组。
   * 实际上，填充这些数组是派生类的构造函数的职责。
   * @p n_faces_or_subfaces
   * 是这个对象所要存储的面或子面的数量。实际数量取决于派生类，对于FEFaceValues来说，它是<tt>2*dim</tt>，而对于FESubfaceValues类来说，它是<tt>2*dim*(1<<(dim-1))</tt>，即面的数量乘以每个面的子面的数量。
   *
   */
  FEFaceValuesBase(const unsigned int                  dofs_per_cell,
                   const UpdateFlags                   update_flags,
                   const Mapping<dim, spacedim> &      mapping,
                   const FiniteElement<dim, spacedim> &fe,
                   const Quadrature<dim - 1> &         quadrature);

  /**
   * 像上面的函数一样，但采取的是正交规则的集合。这允许给每个面分配一个不同的正交规则。在集合只包含一个面的正交规则的情况下，这个正交规则将用于所有面。
   *
   */
  FEFaceValuesBase(const unsigned int                  dofs_per_cell,
                   const UpdateFlags                   update_flags,
                   const Mapping<dim, spacedim> &      mapping,
                   const FiniteElement<dim, spacedim> &fe,
                   const hp::QCollection<dim - 1> &    quadrature);

  /**
   * 单元在<tt>i</tt>第1个正交点的变换的边界形式。 见
   * @ref GlossBoundaryForm  。
   * @dealiiRequiresUpdateFlags{update_boundary_forms}
   *
   */
  const Tensor<1, spacedim> &
  boundary_form(const unsigned int i) const;

  /**
   * 返回曲面映射的外向法向量乘以雅各布系数的列表。
   * @dealiiRequiresUpdateFlags{update_boundary_forms}
   *
   */
  const std::vector<Tensor<1, spacedim>> &
  get_boundary_forms() const;

  /**
   * 返回上次调用reinit()函数时选择的面的索引。
   *
   */
  unsigned int
  get_face_index() const;

  /**
   * 返回该对象所存储的正交公式副本的引用。
   *
   */
  const Quadrature<dim - 1> &
  get_quadrature() const;

  /**
   * 确定此对象的内存消耗（以字节为单位）的估计值。
   *
   */
  std::size_t
  memory_consumption() const;

protected:
  /**
   * 上次调用reinit()函数时选择的面的编号。
   *
   */
  unsigned int present_face_no;

  /**
   * 最后一次调用reinit()函数时选择的面的索引。
   *
   */
  unsigned int present_face_index;

  /**
   * 在这里存储一份正交公式的副本。
   *
   */
  const hp::QCollection<dim - 1> quadrature;
};



/**
 * 在一个面上的正交点中评估的有限元。
 * 这个类将FEFaceValuesBase的功能添加到FEValues中；更多的文档请看那里。
 * 由于有限元函数及其导数在单元格边界可能是不连续的，所以这个函数对一个网格面没有限制。但是，这些值从相邻的任何一个单元接近面的时候都有限制。
 *
 *
 * @ingroup feaccess
 *
 *
 */
template <int dim, int spacedim = dim>
class FEFaceValues : public FEFaceValuesBase<dim, spacedim>
{
public:
  /**
   * 该对象所处的尺寸。
   *
   */

  static const unsigned int dimension = dim;

  static const unsigned int space_dimension = spacedim;

  /**
   * 我们对其进行积分的对象的维度。对于本类，这等于
   * <code>dim-1</code>  。
   *
   */
  static const unsigned int integral_dimension = dim - 1;

  /**
   * 构造函数。从映射和有限元对象中获取单元独立数据，匹配正交规则和更新标志。
   *
   */
  FEFaceValues(const Mapping<dim, spacedim> &      mapping,
               const FiniteElement<dim, spacedim> &fe,
               const Quadrature<dim - 1> &         quadrature,
               const UpdateFlags                   update_flags);

  /**
   * 像上面的函数一样，但取一个正交规则的集合。这允许给每个面分配一个不同的正交规则。在集合只包含一个面的正交规则的情况下，这个正交规则将用于所有面。
   *
   */
  FEFaceValues(const Mapping<dim, spacedim> &      mapping,
               const FiniteElement<dim, spacedim> &fe,
               const hp::QCollection<dim - 1> &    quadrature,
               const UpdateFlags                   update_flags);

  /**
   * 构造函数。这个构造函数等同于其他的构造函数，只是它使对象隐含地使用
   * $Q_1$ 映射（即MappingQGeneric(1)类型的对象）。
   *
   */
  FEFaceValues(const FiniteElement<dim, spacedim> &fe,
               const Quadrature<dim - 1> &         quadrature,
               const UpdateFlags                   update_flags);

  /**
   * 像上面的函数一样，但采取的是正交规则的集合。这允许给每个面分配一个不同的正交规则。在集合只包含一个面的正交规则的情况下，这个正交规则将用于所有面。
   *
   */
  FEFaceValues(const FiniteElement<dim, spacedim> &fe,
               const hp::QCollection<dim - 1> &    quadrature,
               const UpdateFlags                   update_flags);

  /**
   * 重新初始化 @p face_no 的 @p cell
   * 面和给定的有限元的梯度、雅克比行列式等。
   *
   */
  template <bool level_dof_access>
  void
  reinit(
    const TriaIterator<DoFCellAccessor<dim, spacedim, level_dof_access>> &cell,
    const unsigned int face_no);

  /**
   * 重新初始化面 @p face 和单元 @p cell.
   * 的梯度、雅可比决定数等。
   * @note   @p face  必须是 @p cell's 面的迭代器之一。
   *
   */
  template <bool level_dof_access>
  void
  reinit(
    const TriaIterator<DoFCellAccessor<dim, spacedim, level_dof_access>> &cell,
    const typename Triangulation<dim, spacedim>::face_iterator &          face);

  /**
   * 重新初始化 "进入三角形对象的迭代器
   * "类型的给定单元上的给定面的梯度、雅可比行列式等，以及给定的有限元。由于进入三角剖分的迭代器只传递单元的几何信息，而不传递可能与此单元相关的自由度信息，如果需要自由度信息，你将无法调用该类的一些函数。这些函数首先是<tt>get_function_value/gradients/hessians/third_derivatives</tt>函数。如果你想调用这些函数，你必须调用
   * @p
   * reinit变体，将迭代器带入DoFHandler或其他DoF处理程序类型对象。
   *
   */
  void
  reinit(const typename Triangulation<dim, spacedim>::cell_iterator &cell,
         const unsigned int                                          face_no);

  /* 重新初始化 "进入三角形对象的迭代器 "类型的给定单元上的给定面的梯度、雅各比行列式等，以及给定的有限元。由于进入三角剖分的迭代器只传递单元的几何信息，而不传递可能与此单元相关的自由度的信息，所以如果需要自由度的信息，你将无法调用这一类的一些函数。这些函数首先是<tt>get_function_value/gradients/hessians/third_derivatives</tt>函数。如果你想调用这些函数，你必须调用 @p  reinit变体，将迭代器带入DoFHandler或其他DoF处理程序类型对象。   
*  @note   @p face 必须是 @p cell's 面的迭代器之一。 
*
*/
  void
  reinit(const typename Triangulation<dim, spacedim>::cell_iterator &cell,
         const typename Triangulation<dim, spacedim>::face_iterator &face);

  /**
   * 返回对这个非常对象的引用。
   * 虽然看起来不是很有用，但是这个函数的存在是为了给
   * hp::FEValues
   * 类提供能力，在这种情况下，它提供了当前单元的FEValues对象（请记住，对于HP-finite元素，实际使用的FE对象可能会在不同的单元之间变化，所以我们也需要不同的FEValues对象用于不同的单元；一旦你重新初始化了
   * hp::FEValues
   * ]对象，它就会为该单元格上的FE检索FEValues对象，并通过与此相同的函数返回；因此，这里的这个函数只提供相同的接口，以便人们可以对FEValues和
   * hp::FEValues). 进行模板化。
   *
   */
  const FEFaceValues<dim, spacedim> &
  get_present_fe_values() const;

private:
  /**
   * 做两个构造函数的共同工作。
   *
   */
  void
  initialize(const UpdateFlags update_flags);

  /**
   * reinit()函数只做需要了解迭代器类型的那部分工作。在设置完present_cell()之后，它们会传递给这个函数，这个函数会做真正的工作，而且与单元格迭代器的实际类型无关。
   *
   */
  void
  do_reinit(const unsigned int face_no);
};


/**
 * 在一个面的正交点上评估的有限元。
 * 这个类将FEFaceValuesBase的功能添加到FEValues中；更多的文档请看那里。
 * 该类用于位于细化边上的面。在这种情况下，邻近的单元被细化。为了能够计算内部和外部函数值之间的差异，邻近单元的细化必须在这个单元上模拟。这可以通过应用模拟细化的正交规则来实现。由此产生的数据字段被分割开来，以反映邻居的细化结构：一个子面的编号对应于邻居面的孩子的编号。
 *
 *
 * @ingroup feaccess
 *
 *
 */
template <int dim, int spacedim = dim>
class FESubfaceValues : public FEFaceValuesBase<dim, spacedim>
{
public:
  /**
   * 该对象所处的尺寸。
   *
   */
  static const unsigned int dimension = dim;

  /**
   * 该对象所处空间的尺寸。
   *
   */
  static const unsigned int space_dimension = spacedim;

  /**
   * 我们所整合的对象的维度。对于本类，这等于
   * <code>dim-1</code>  。
   *
   */
  static const unsigned int integral_dimension = dim - 1;

  /**
   * 构造函数。从映射和有限元对象中获取独立单元数据，匹配正交规则和更新标志。
   *
   */
  FESubfaceValues(const Mapping<dim, spacedim> &      mapping,
                  const FiniteElement<dim, spacedim> &fe,
                  const Quadrature<dim - 1> &         face_quadrature,
                  const UpdateFlags                   update_flags);

  /**
   * 像上面的函数一样，但取一个正交规则的集合。
   * @note
   * 与FEFaceValues相反，我们要求集合中正交规则的数量为1。
   *
   */
  FESubfaceValues(const Mapping<dim, spacedim> &      mapping,
                  const FiniteElement<dim, spacedim> &fe,
                  const hp::QCollection<dim - 1> &    face_quadrature,
                  const UpdateFlags                   update_flags);

  /**
   * 构造函数。这个构造函数除了使对象隐含地使用 $Q_1$
   * 映射（即MappingQGeneric(1)类型的对象）外，与其他构造函数是等价的。
   *
   */
  FESubfaceValues(const FiniteElement<dim, spacedim> &fe,
                  const Quadrature<dim - 1> &         face_quadrature,
                  const UpdateFlags                   update_flags);

  /**
   * 像上面的函数一样，但取一个正交规则的集合。
   * @note
   * 与FEFaceValues相反，我们要求集合中的正交规则的数量为1。
   *
   */
  FESubfaceValues(const FiniteElement<dim, spacedim> &fe,
                  const hp::QCollection<dim - 1> &    face_quadrature,
                  const UpdateFlags                   update_flags);

  /**
   * 重新初始化类型为 "iterator into a DoFHandler object
   * "的给定单元的梯度、雅各比行列式等，以及与该对象相关的有限元。假设给定单元所使用的有限元也是这个FESubfaceValues对象所使用的有限元。
   *
   */
  template <bool level_dof_access>
  void
  reinit(
    const TriaIterator<DoFCellAccessor<dim, spacedim, level_dof_access>> &cell,
    const unsigned int face_no,
    const unsigned int subface_no);

  /**
   * 替代的重新初始化函数，作为参数，接受面和子面的迭代器，而不是它们的数字。
   *
   */
  template <bool level_dof_access>
  void
  reinit(
    const TriaIterator<DoFCellAccessor<dim, spacedim, level_dof_access>> &cell,
    const typename Triangulation<dim, spacedim>::face_iterator &          face,
    const typename Triangulation<dim, spacedim>::face_iterator &subface);

  /**
   * 重新初始化给定子面的梯度、雅可比行列式等，这些梯度和雅可比行列式是在给定的
   * "进入三角形对象的迭代器
   * "和给定的有限元的单元上进行的。由于进入三角剖分的迭代器只传递单元的几何信息，而不传递可能与此单元相关的自由度的信息，如果需要自由度的信息，你将无法调用这一类的一些函数。这些函数首先是<tt>get_function_value/gradients/hessians/third_derivatives</tt>函数。如果你想调用这些函数，你必须调用
   * @p
   * reinit变体，将迭代器带入DoFHandler或其他DoF处理程序类型对象。
   *
   */
  void
  reinit(const typename Triangulation<dim, spacedim>::cell_iterator &cell,
         const unsigned int                                          face_no,
         const unsigned int subface_no);

  /**
   * 重新初始化给定子面的梯度、雅可比行列式等，这些梯度、雅可比行列式是在给定的
   * "进入三角形对象的迭代器
   * "类型的单元上，以及给定的有限元上。由于进入三角剖分的迭代器只传递单元的几何信息，而不传递可能与此单元相关的自由度的信息，如果需要自由度的信息，你将无法调用这一类的一些函数。这些函数首先是<tt>get_function_value/gradients/hessians/third_derivatives</tt>函数。如果你想调用这些函数，你必须调用
   * @p
   * reinit变体，将迭代器带入DoFHandler或其他DoF处理程序类型对象。
   * 这和前面的函数做的事情一样，但是把迭代器而不是数字作为参数。
   * @note   @p face  和  @p subface  必须对应于  @p cell.
   * 的一个面（以及该面的一个子面）。
   *
   */
  void
  reinit(const typename Triangulation<dim, spacedim>::cell_iterator &cell,
         const typename Triangulation<dim, spacedim>::face_iterator &face,
         const typename Triangulation<dim, spacedim>::face_iterator &subface);

  /**
   * 返回一个对这个对象的引用。
   * 虽然看起来不是很有用，但是这个函数的存在是为了给
   * hp::FEValues
   * 类提供能力，在这种情况下，它提供了当前单元的FEValues对象（记住，对于hp-finite元素，实际使用的FE对象可能在不同的单元中发生变化，所以我们也需要不同的FEValues对象用于不同的单元；一旦你重新初始化了
   * hp::FEValues
   * ]对象，它就会为该单元格上的FE检索FEValues对象，并通过与此相同的函数返回；因此，这里的这个函数只提供相同的接口，以便人们可以对FEValues和
   * hp::FEValues). 进行模板化。
   *
   */
  const FESubfaceValues<dim, spacedim> &
  get_present_fe_values() const;

  /**
   * @todo  Document this
   * @ingroup Exceptions
   *
   */
  DeclException0(ExcReinitCalledWithBoundaryFace);

  /**
   * @todo  记录这一点
   * @ingroup Exceptions
   *
   */
  DeclException0(ExcFaceHasNoSubfaces);

private:
  /**
   * 做两个构造函数的共同工作。
   *
   */
  void
  initialize(const UpdateFlags update_flags);

  /**
   * reinit()函数只做需要了解迭代器类型的那部分工作。在设置完present_cell()之后，它们会传递给这个函数，这个函数会做真正的工作，而且与单元格迭代器的实际类型无关。
   *
   */
  void
  do_reinit(const unsigned int face_no, const unsigned int subface_no);
};


#ifndef DOXYGEN


 /*------------------------ Inline functions: namespace FEValuesViews --------*/ 

namespace FEValuesViews
{
  template <int dim, int spacedim>
  inline typename Scalar<dim, spacedim>::value_type
  Scalar<dim, spacedim>::value(const unsigned int shape_function,
                               const unsigned int q_point) const
  {
    AssertIndexRange(shape_function, fe_values->fe->n_dofs_per_cell());
    Assert(
      fe_values->update_flags & update_values,
      ((typename FEValuesBase<dim, spacedim>::ExcAccessToUninitializedField(
        "update_values"))));

    // an adaptation of the FEValuesBase::shape_value_component function
    // except that here we know the component as fixed and we have
    // pre-computed and cached a bunch of information. See the comments there.
    if (shape_function_data[shape_function].is_nonzero_shape_function_component)
      return fe_values->finite_element_output.shape_values(
        shape_function_data[shape_function].row_index, q_point);
    else
      return 0;
  }



  template <int dim, int spacedim>
  inline typename Scalar<dim, spacedim>::gradient_type
  Scalar<dim, spacedim>::gradient(const unsigned int shape_function,
                                  const unsigned int q_point) const
  {
    AssertIndexRange(shape_function, fe_values->fe->n_dofs_per_cell());
    Assert(fe_values->update_flags & update_gradients,
           (typename FEValuesBase<dim, spacedim>::ExcAccessToUninitializedField(
             "update_gradients")));

    // an adaptation of the FEValuesBase::shape_grad_component
    // function except that here we know the component as fixed and we have
    // pre-computed and cached a bunch of information. See the comments there.
    if (shape_function_data[shape_function].is_nonzero_shape_function_component)
      return fe_values->finite_element_output
        .shape_gradients[shape_function_data[shape_function].row_index]
                        [q_point];
    else
      return gradient_type();
  }



  template <int dim, int spacedim>
  inline typename Scalar<dim, spacedim>::hessian_type
  Scalar<dim, spacedim>::hessian(const unsigned int shape_function,
                                 const unsigned int q_point) const
  {
    AssertIndexRange(shape_function, fe_values->fe->n_dofs_per_cell());
    Assert(fe_values->update_flags & update_hessians,
           (typename FEValuesBase<dim, spacedim>::ExcAccessToUninitializedField(
             "update_hessians")));

    // an adaptation of the FEValuesBase::shape_hessian_component
    // function except that here we know the component as fixed and we have
    // pre-computed and cached a bunch of information. See the comments there.
    if (shape_function_data[shape_function].is_nonzero_shape_function_component)
      return fe_values->finite_element_output
        .shape_hessians[shape_function_data[shape_function].row_index][q_point];
    else
      return hessian_type();
  }



  template <int dim, int spacedim>
  inline typename Scalar<dim, spacedim>::third_derivative_type
  Scalar<dim, spacedim>::third_derivative(const unsigned int shape_function,
                                          const unsigned int q_point) const
  {
    AssertIndexRange(shape_function, fe_values->fe->n_dofs_per_cell());
    Assert(fe_values->update_flags & update_3rd_derivatives,
           (typename FEValuesBase<dim, spacedim>::ExcAccessToUninitializedField(
             "update_3rd_derivatives")));

    // an adaptation of the FEValuesBase::shape_3rdderivative_component
    // function except that here we know the component as fixed and we have
    // pre-computed and cached a bunch of information. See the comments there.
    if (shape_function_data[shape_function].is_nonzero_shape_function_component)
      return fe_values->finite_element_output
        .shape_3rd_derivatives[shape_function_data[shape_function].row_index]
                              [q_point];
    else
      return third_derivative_type();
  }



  template <int dim, int spacedim>
  inline typename Vector<dim, spacedim>::value_type
  Vector<dim, spacedim>::value(const unsigned int shape_function,
                               const unsigned int q_point) const
  {
    AssertIndexRange(shape_function, fe_values->fe->n_dofs_per_cell());
    Assert(fe_values->update_flags & update_values,
           (typename FEValuesBase<dim, spacedim>::ExcAccessToUninitializedField(
             "update_values")));

    // same as for the scalar case except that we have one more index
    const int snc =
      shape_function_data[shape_function].single_nonzero_component;
    if (snc == -2)
      return value_type();
    else if (snc != -1)
      {
        value_type return_value;
        return_value[shape_function_data[shape_function]
                       .single_nonzero_component_index] =
          fe_values->finite_element_output.shape_values(snc, q_point);
        return return_value;
      }
    else
      {
        value_type return_value;
        for (unsigned int d = 0; d < dim; ++d)
          if (shape_function_data[shape_function]
                .is_nonzero_shape_function_component[d])
            return_value[d] = fe_values->finite_element_output.shape_values(
              shape_function_data[shape_function].row_index[d], q_point);

        return return_value;
      }
  }



  template <int dim, int spacedim>
  inline typename Vector<dim, spacedim>::gradient_type
  Vector<dim, spacedim>::gradient(const unsigned int shape_function,
                                  const unsigned int q_point) const
  {
    AssertIndexRange(shape_function, fe_values->fe->n_dofs_per_cell());
    Assert(fe_values->update_flags & update_gradients,
           (typename FEValuesBase<dim, spacedim>::ExcAccessToUninitializedField(
             "update_gradients")));

    // same as for the scalar case except that we have one more index
    const int snc =
      shape_function_data[shape_function].single_nonzero_component;
    if (snc == -2)
      return gradient_type();
    else if (snc != -1)
      {
        gradient_type return_value;
        return_value[shape_function_data[shape_function]
                       .single_nonzero_component_index] =
          fe_values->finite_element_output.shape_gradients[snc][q_point];
        return return_value;
      }
    else
      {
        gradient_type return_value;
        for (unsigned int d = 0; d < dim; ++d)
          if (shape_function_data[shape_function]
                .is_nonzero_shape_function_component[d])
            return_value[d] =
              fe_values->finite_element_output.shape_gradients
                [shape_function_data[shape_function].row_index[d]][q_point];

        return return_value;
      }
  }



  template <int dim, int spacedim>
  inline typename Vector<dim, spacedim>::divergence_type
  Vector<dim, spacedim>::divergence(const unsigned int shape_function,
                                    const unsigned int q_point) const
  {
    // this function works like in the case above
    AssertIndexRange(shape_function, fe_values->fe->n_dofs_per_cell());
    Assert(fe_values->update_flags & update_gradients,
           (typename FEValuesBase<dim, spacedim>::ExcAccessToUninitializedField(
             "update_gradients")));

    // same as for the scalar case except that we have one more index
    const int snc =
      shape_function_data[shape_function].single_nonzero_component;
    if (snc == -2)
      return divergence_type();
    else if (snc != -1)
      return fe_values->finite_element_output
        .shape_gradients[snc][q_point][shape_function_data[shape_function]
                                         .single_nonzero_component_index];
    else
      {
        divergence_type return_value = 0;
        for (unsigned int d = 0; d < dim; ++d)
          if (shape_function_data[shape_function]
                .is_nonzero_shape_function_component[d])
            return_value +=
              fe_values->finite_element_output.shape_gradients
                [shape_function_data[shape_function].row_index[d]][q_point][d];

        return return_value;
      }
  }



  template <int dim, int spacedim>
  inline typename Vector<dim, spacedim>::curl_type
  Vector<dim, spacedim>::curl(const unsigned int shape_function,
                              const unsigned int q_point) const
  {
    // this function works like in the case above

    AssertIndexRange(shape_function, fe_values->fe->n_dofs_per_cell());
    Assert(fe_values->update_flags & update_gradients,
           (typename FEValuesBase<dim, spacedim>::ExcAccessToUninitializedField(
             "update_gradients")));
    // same as for the scalar case except that we have one more index
    const int snc =
      shape_function_data[shape_function].single_nonzero_component;

    if (snc == -2)
      return curl_type();

    else
      switch (dim)
        {
          case 1:
            {
              Assert(false,
                     ExcMessage(
                       "Computing the curl in 1d is not a useful operation"));
              return curl_type();
            }

          case 2:
            {
              if (snc != -1)
                {
                  curl_type return_value;

                  // the single nonzero component can only be zero or one in 2d
                  if (shape_function_data[shape_function]
                        .single_nonzero_component_index == 0)
                    return_value[0] =
                      -1.0 * fe_values->finite_element_output
                               .shape_gradients[snc][q_point][1];
                  else
                    return_value[0] = fe_values->finite_element_output
                                        .shape_gradients[snc][q_point][0];

                  return return_value;
                }

              else
                {
                  curl_type return_value;

                  return_value[0] = 0.0;

                  if (shape_function_data[shape_function]
                        .is_nonzero_shape_function_component[0])
                    return_value[0] -=
                      fe_values->finite_element_output
                        .shape_gradients[shape_function_data[shape_function]
                                           .row_index[0]][q_point][1];

                  if (shape_function_data[shape_function]
                        .is_nonzero_shape_function_component[1])
                    return_value[0] +=
                      fe_values->finite_element_output
                        .shape_gradients[shape_function_data[shape_function]
                                           .row_index[1]][q_point][0];

                  return return_value;
                }
            }

          case 3:
            {
              if (snc != -1)
                {
                  curl_type return_value;

                  switch (shape_function_data[shape_function]
                            .single_nonzero_component_index)
                    {
                      case 0:
                        {
                          return_value[0] = 0;
                          return_value[1] = fe_values->finite_element_output
                                              .shape_gradients[snc][q_point][2];
                          return_value[2] =
                            -1.0 * fe_values->finite_element_output
                                     .shape_gradients[snc][q_point][1];
                          return return_value;
                        }

                      case 1:
                        {
                          return_value[0] =
                            -1.0 * fe_values->finite_element_output
                                     .shape_gradients[snc][q_point][2];
                          return_value[1] = 0;
                          return_value[2] = fe_values->finite_element_output
                                              .shape_gradients[snc][q_point][0];
                          return return_value;
                        }

                      default:
                        {
                          return_value[0] = fe_values->finite_element_output
                                              .shape_gradients[snc][q_point][1];
                          return_value[1] =
                            -1.0 * fe_values->finite_element_output
                                     .shape_gradients[snc][q_point][0];
                          return_value[2] = 0;
                          return return_value;
                        }
                    }
                }

              else
                {
                  curl_type return_value;

                  for (unsigned int i = 0; i < dim; ++i)
                    return_value[i] = 0.0;

                  if (shape_function_data[shape_function]
                        .is_nonzero_shape_function_component[0])
                    {
                      return_value[1] +=
                        fe_values->finite_element_output
                          .shape_gradients[shape_function_data[shape_function]
                                             .row_index[0]][q_point][2];
                      return_value[2] -=
                        fe_values->finite_element_output
                          .shape_gradients[shape_function_data[shape_function]
                                             .row_index[0]][q_point][1];
                    }

                  if (shape_function_data[shape_function]
                        .is_nonzero_shape_function_component[1])
                    {
                      return_value[0] -=
                        fe_values->finite_element_output
                          .shape_gradients[shape_function_data[shape_function]
                                             .row_index[1]][q_point][2];
                      return_value[2] +=
                        fe_values->finite_element_output
                          .shape_gradients[shape_function_data[shape_function]
                                             .row_index[1]][q_point][0];
                    }

                  if (shape_function_data[shape_function]
                        .is_nonzero_shape_function_component[2])
                    {
                      return_value[0] +=
                        fe_values->finite_element_output
                          .shape_gradients[shape_function_data[shape_function]
                                             .row_index[2]][q_point][1];
                      return_value[1] -=
                        fe_values->finite_element_output
                          .shape_gradients[shape_function_data[shape_function]
                                             .row_index[2]][q_point][0];
                    }

                  return return_value;
                }
            }
        }
    // should not end up here
    Assert(false, ExcInternalError());
    return curl_type();
  }



  template <int dim, int spacedim>
  inline typename Vector<dim, spacedim>::hessian_type
  Vector<dim, spacedim>::hessian(const unsigned int shape_function,
                                 const unsigned int q_point) const
  {
    // this function works like in the case above
    AssertIndexRange(shape_function, fe_values->fe->n_dofs_per_cell());
    Assert(fe_values->update_flags & update_hessians,
           (typename FEValuesBase<dim, spacedim>::ExcAccessToUninitializedField(
             "update_hessians")));

    // same as for the scalar case except that we have one more index
    const int snc =
      shape_function_data[shape_function].single_nonzero_component;
    if (snc == -2)
      return hessian_type();
    else if (snc != -1)
      {
        hessian_type return_value;
        return_value[shape_function_data[shape_function]
                       .single_nonzero_component_index] =
          fe_values->finite_element_output.shape_hessians[snc][q_point];
        return return_value;
      }
    else
      {
        hessian_type return_value;
        for (unsigned int d = 0; d < dim; ++d)
          if (shape_function_data[shape_function]
                .is_nonzero_shape_function_component[d])
            return_value[d] =
              fe_values->finite_element_output.shape_hessians
                [shape_function_data[shape_function].row_index[d]][q_point];

        return return_value;
      }
  }



  template <int dim, int spacedim>
  inline typename Vector<dim, spacedim>::third_derivative_type
  Vector<dim, spacedim>::third_derivative(const unsigned int shape_function,
                                          const unsigned int q_point) const
  {
    // this function works like in the case above
    AssertIndexRange(shape_function, fe_values->fe->n_dofs_per_cell());
    Assert(fe_values->update_flags & update_3rd_derivatives,
           (typename FEValuesBase<dim, spacedim>::ExcAccessToUninitializedField(
             "update_3rd_derivatives")));

    // same as for the scalar case except that we have one more index
    const int snc =
      shape_function_data[shape_function].single_nonzero_component;
    if (snc == -2)
      return third_derivative_type();
    else if (snc != -1)
      {
        third_derivative_type return_value;
        return_value[shape_function_data[shape_function]
                       .single_nonzero_component_index] =
          fe_values->finite_element_output.shape_3rd_derivatives[snc][q_point];
        return return_value;
      }
    else
      {
        third_derivative_type return_value;
        for (unsigned int d = 0; d < dim; ++d)
          if (shape_function_data[shape_function]
                .is_nonzero_shape_function_component[d])
            return_value[d] =
              fe_values->finite_element_output.shape_3rd_derivatives
                [shape_function_data[shape_function].row_index[d]][q_point];

        return return_value;
      }
  }



  namespace internal
  {
    /**
     * 返回一个张量的对称版本，该张量的第n行等于第二个参数，其他所有的行都等于0。
     *
     */
    inline dealii::SymmetricTensor<2, 1>
    symmetrize_single_row(const unsigned int n, const dealii::Tensor<1, 1> &t)
    {
      AssertIndexRange(n, 1);
      (void)n;

      return {{t[0]}};
    }



    inline dealii::SymmetricTensor<2, 2>
    symmetrize_single_row(const unsigned int n, const dealii::Tensor<1, 2> &t)
    {
      switch (n)
        {
          case 0:
            {
              return {{t[0], 0, t[1] / 2}};
            }
          case 1:
            {
              return {{0, t[1], t[0] / 2}};
            }
          default:
            {
              AssertIndexRange(n, 2);
              return {};
            }
        }
    }



    inline dealii::SymmetricTensor<2, 3>
    symmetrize_single_row(const unsigned int n, const dealii::Tensor<1, 3> &t)
    {
      switch (n)
        {
          case 0:
            {
              return {{t[0], 0, 0, t[1] / 2, t[2] / 2, 0}};
            }
          case 1:
            {
              return {{0, t[1], 0, t[0] / 2, 0, t[2] / 2}};
            }
          case 2:
            {
              return {{0, 0, t[2], 0, t[0] / 2, t[1] / 2}};
            }
          default:
            {
              AssertIndexRange(n, 3);
              return {};
            }
        }
    }
  } // namespace internal



  template <int dim, int spacedim>
  inline typename Vector<dim, spacedim>::symmetric_gradient_type
  Vector<dim, spacedim>::symmetric_gradient(const unsigned int shape_function,
                                            const unsigned int q_point) const
  {
    AssertIndexRange(shape_function, fe_values->fe->n_dofs_per_cell());
    Assert(fe_values->update_flags & update_gradients,
           (typename FEValuesBase<dim, spacedim>::ExcAccessToUninitializedField(
             "update_gradients")));

    // same as for the scalar case except that we have one more index
    const int snc =
      shape_function_data[shape_function].single_nonzero_component;
    if (snc == -2)
      return symmetric_gradient_type();
    else if (snc != -1)
      return internal::symmetrize_single_row(
        shape_function_data[shape_function].single_nonzero_component_index,
        fe_values->finite_element_output.shape_gradients[snc][q_point]);
    else
      {
        gradient_type return_value;
        for (unsigned int d = 0; d < dim; ++d)
          if (shape_function_data[shape_function]
                .is_nonzero_shape_function_component[d])
            return_value[d] =
              fe_values->finite_element_output.shape_gradients
                [shape_function_data[shape_function].row_index[d]][q_point];

        return symmetrize(return_value);
      }
  }



  template <int dim, int spacedim>
  inline typename SymmetricTensor<2, dim, spacedim>::value_type
  SymmetricTensor<2, dim, spacedim>::value(const unsigned int shape_function,
                                           const unsigned int q_point) const
  {
    AssertIndexRange(shape_function, fe_values->fe->n_dofs_per_cell());
    Assert(fe_values->update_flags & update_values,
           (typename FEValuesBase<dim, spacedim>::ExcAccessToUninitializedField(
             "update_values")));

    // similar to the vector case where we have more then one index and we need
    // to convert between unrolled and component indexing for tensors
    const int snc =
      shape_function_data[shape_function].single_nonzero_component;

    if (snc == -2)
      {
        // shape function is zero for the selected components
        return value_type();
      }
    else if (snc != -1)
      {
        value_type         return_value;
        const unsigned int comp =
          shape_function_data[shape_function].single_nonzero_component_index;
        return_value[value_type::unrolled_to_component_indices(comp)] =
          fe_values->finite_element_output.shape_values(snc, q_point);
        return return_value;
      }
    else
      {
        value_type return_value;
        for (unsigned int d = 0; d < value_type::n_independent_components; ++d)
          if (shape_function_data[shape_function]
                .is_nonzero_shape_function_component[d])
            return_value[value_type::unrolled_to_component_indices(d)] =
              fe_values->finite_element_output.shape_values(
                shape_function_data[shape_function].row_index[d], q_point);
        return return_value;
      }
  }



  template <int dim, int spacedim>
  inline typename SymmetricTensor<2, dim, spacedim>::divergence_type
  SymmetricTensor<2, dim, spacedim>::divergence(
    const unsigned int shape_function,
    const unsigned int q_point) const
  {
    AssertIndexRange(shape_function, fe_values->fe->n_dofs_per_cell());
    Assert(fe_values->update_flags & update_gradients,
           (typename FEValuesBase<dim, spacedim>::ExcAccessToUninitializedField(
             "update_gradients")));

    const int snc =
      shape_function_data[shape_function].single_nonzero_component;

    if (snc == -2)
      {
        // shape function is zero for the selected components
        return divergence_type();
      }
    else if (snc != -1)
      {
        // we have a single non-zero component when the symmetric tensor is
        // represented in unrolled form. this implies we potentially have
        // two non-zero components when represented in component form!  we
        // will only have one non-zero entry if the non-zero component lies on
        // the diagonal of the tensor.
        //
        // the divergence of a second-order tensor is a first order tensor.
        //
        // assume the second-order tensor is A with components A_{ij}.  then
        // A_{ij} = A_{ji} and there is only one (if diagonal) or two non-zero
        // entries in the tensorial representation.  define the
        // divergence as:
        // b_i \dealcoloneq \dfrac{\partial phi_{ij}}{\partial x_j}.
        // (which is incidentally also
        // b_j \dealcoloneq \dfrac{\partial phi_{ij}}{\partial x_i}).
        // In both cases, a sum is implied.
        //
        // Now, we know the nonzero component in unrolled form: it is indicated
        // by 'snc'. we can figure out which tensor components belong to this:
        const unsigned int comp =
          shape_function_data[shape_function].single_nonzero_component_index;
        const unsigned int ii =
          value_type::unrolled_to_component_indices(comp)[0];
        const unsigned int jj =
          value_type::unrolled_to_component_indices(comp)[1];

        // given the form of the divergence above, if ii=jj there is only a
        // single nonzero component of the full tensor and the gradient
        // equals
        // b_ii \dealcoloneq \dfrac{\partial phi_{ii,ii}}{\partial x_ii}.
        // all other entries of 'b' are zero
        //
        // on the other hand, if ii!=jj, then there are two nonzero entries in
        // the full tensor and
        // b_ii \dealcoloneq \dfrac{\partial phi_{ii,jj}}{\partial x_ii}.
        // b_jj \dealcoloneq \dfrac{\partial phi_{ii,jj}}{\partial x_jj}.
        // again, all other entries of 'b' are zero
        const dealii::Tensor<1, spacedim> &phi_grad =
          fe_values->finite_element_output.shape_gradients[snc][q_point];

        divergence_type return_value;
        return_value[ii] = phi_grad[jj];

        if (ii != jj)
          return_value[jj] = phi_grad[ii];

        return return_value;
      }
    else
      {
        Assert(false, ExcNotImplemented());
        divergence_type return_value;
        return return_value;
      }
  }



  template <int dim, int spacedim>
  inline typename Tensor<2, dim, spacedim>::value_type
  Tensor<2, dim, spacedim>::value(const unsigned int shape_function,
                                  const unsigned int q_point) const
  {
    AssertIndexRange(shape_function, fe_values->fe->n_dofs_per_cell());
    Assert(fe_values->update_flags & update_values,
           (typename FEValuesBase<dim, spacedim>::ExcAccessToUninitializedField(
             "update_values")));

    // similar to the vector case where we have more then one index and we need
    // to convert between unrolled and component indexing for tensors
    const int snc =
      shape_function_data[shape_function].single_nonzero_component;

    if (snc == -2)
      {
        // shape function is zero for the selected components
        return value_type();
      }
    else if (snc != -1)
      {
        value_type         return_value;
        const unsigned int comp =
          shape_function_data[shape_function].single_nonzero_component_index;
        const TableIndices<2> indices =
          dealii::Tensor<2, spacedim>::unrolled_to_component_indices(comp);
        return_value[indices] =
          fe_values->finite_element_output.shape_values(snc, q_point);
        return return_value;
      }
    else
      {
        value_type return_value;
        for (unsigned int d = 0; d < dim * dim; ++d)
          if (shape_function_data[shape_function]
                .is_nonzero_shape_function_component[d])
            {
              const TableIndices<2> indices =
                dealii::Tensor<2, spacedim>::unrolled_to_component_indices(d);
              return_value[indices] =
                fe_values->finite_element_output.shape_values(
                  shape_function_data[shape_function].row_index[d], q_point);
            }
        return return_value;
      }
  }



  template <int dim, int spacedim>
  inline typename Tensor<2, dim, spacedim>::divergence_type
  Tensor<2, dim, spacedim>::divergence(const unsigned int shape_function,
                                       const unsigned int q_point) const
  {
    AssertIndexRange(shape_function, fe_values->fe->n_dofs_per_cell());
    Assert(fe_values->update_flags & update_gradients,
           (typename FEValuesBase<dim, spacedim>::ExcAccessToUninitializedField(
             "update_gradients")));

    const int snc =
      shape_function_data[shape_function].single_nonzero_component;

    if (snc == -2)
      {
        // shape function is zero for the selected components
        return divergence_type();
      }
    else if (snc != -1)
      {
        // we have a single non-zero component when the tensor is
        // represented in unrolled form.
        //
        // the divergence of a second-order tensor is a first order tensor.
        //
        // assume the second-order tensor is A with components A_{ij},
        // then divergence is d_i := \frac{\partial A_{ij}}{\partial x_j}
        //
        // Now, we know the nonzero component in unrolled form: it is indicated
        // by 'snc'. we can figure out which tensor components belong to this:
        const unsigned int comp =
          shape_function_data[shape_function].single_nonzero_component_index;
        const TableIndices<2> indices =
          dealii::Tensor<2, spacedim>::unrolled_to_component_indices(comp);
        const unsigned int ii = indices[0];
        const unsigned int jj = indices[1];

        const dealii::Tensor<1, spacedim> &phi_grad =
          fe_values->finite_element_output.shape_gradients[snc][q_point];

        divergence_type return_value;
        // note that we contract \nabla from the right
        return_value[ii] = phi_grad[jj];

        return return_value;
      }
    else
      {
        Assert(false, ExcNotImplemented());
        divergence_type return_value;
        return return_value;
      }
  }



  template <int dim, int spacedim>
  inline typename Tensor<2, dim, spacedim>::gradient_type
  Tensor<2, dim, spacedim>::gradient(const unsigned int shape_function,
                                     const unsigned int q_point) const
  {
    AssertIndexRange(shape_function, fe_values->fe->n_dofs_per_cell());
    Assert(fe_values->update_flags & update_gradients,
           (typename FEValuesBase<dim, spacedim>::ExcAccessToUninitializedField(
             "update_gradients")));

    const int snc =
      shape_function_data[shape_function].single_nonzero_component;

    if (snc == -2)
      {
        // shape function is zero for the selected components
        return gradient_type();
      }
    else if (snc != -1)
      {
        // we have a single non-zero component when the tensor is
        // represented in unrolled form.
        //
        // the gradient of a second-order tensor is a third order tensor.
        //
        // assume the second-order tensor is A with components A_{ij},
        // then gradient is B_{ijk} := \frac{\partial A_{ij}}{\partial x_k}
        //
        // Now, we know the nonzero component in unrolled form: it is indicated
        // by 'snc'. we can figure out which tensor components belong to this:
        const unsigned int comp =
          shape_function_data[shape_function].single_nonzero_component_index;
        const TableIndices<2> indices =
          dealii::Tensor<2, spacedim>::unrolled_to_component_indices(comp);
        const unsigned int ii = indices[0];
        const unsigned int jj = indices[1];

        const dealii::Tensor<1, spacedim> &phi_grad =
          fe_values->finite_element_output.shape_gradients[snc][q_point];

        gradient_type return_value;
        return_value[ii][jj] = phi_grad;

        return return_value;
      }
    else
      {
        Assert(false, ExcNotImplemented());
        gradient_type return_value;
        return return_value;
      }
  }

} // namespace FEValuesViews



 /*---------------------- Inline functions: FEValuesBase ---------------------*/ 



template <int dim, int spacedim>
inline const FEValuesViews::Scalar<dim, spacedim> &FEValuesBase<dim, spacedim>::
                                                   operator[](const FEValuesExtractors::Scalar &scalar) const
{
  AssertIndexRange(scalar.component, fe_values_views_cache.scalars.size());

  return fe_values_views_cache.scalars[scalar.component];
}



template <int dim, int spacedim>
inline const FEValuesViews::Vector<dim, spacedim> &FEValuesBase<dim, spacedim>::
                                                   operator[](const FEValuesExtractors::Vector &vector) const
{
  AssertIndexRange(vector.first_vector_component,
                   fe_values_views_cache.vectors.size());

  return fe_values_views_cache.vectors[vector.first_vector_component];
}



template <int dim, int spacedim>
inline const FEValuesViews::SymmetricTensor<2, dim, spacedim> &
  FEValuesBase<dim, spacedim>::
  operator[](const FEValuesExtractors::SymmetricTensor<2> &tensor) const
{
  Assert(
    tensor.first_tensor_component <
      fe_values_views_cache.symmetric_second_order_tensors.size(),
    ExcIndexRange(tensor.first_tensor_component,
                  0,
                  fe_values_views_cache.symmetric_second_order_tensors.size()));

  return fe_values_views_cache
    .symmetric_second_order_tensors[tensor.first_tensor_component];
}



template <int dim, int spacedim>
inline const FEValuesViews::Tensor<2, dim, spacedim> &
  FEValuesBase<dim, spacedim>::
  operator[](const FEValuesExtractors::Tensor<2> &tensor) const
{
  AssertIndexRange(tensor.first_tensor_component,
                   fe_values_views_cache.second_order_tensors.size());

  return fe_values_views_cache
    .second_order_tensors[tensor.first_tensor_component];
}



template <int dim, int spacedim>
inline const double &
FEValuesBase<dim, spacedim>::shape_value(const unsigned int i,
                                         const unsigned int j) const
{
  AssertIndexRange(i, fe->n_dofs_per_cell());
  Assert(this->update_flags & update_values,
         ExcAccessToUninitializedField("update_values"));
  Assert(fe->is_primitive(i), ExcShapeFunctionNotPrimitive(i));
  Assert(present_cell.get() != nullptr,
         ExcMessage("FEValues object is not reinit'ed to any cell"));
  // if the entire FE is primitive,
  // then we can take a short-cut:
  if (fe->is_primitive())
    return this->finite_element_output.shape_values(i, j);
  else
    {
      // otherwise, use the mapping
      // between shape function
      // numbers and rows. note that
      // by the assertions above, we
      // know that this particular
      // shape function is primitive,
      // so we can call
      // system_to_component_index
      const unsigned int row =
        this->finite_element_output
          .shape_function_to_row_table[i * fe->n_components() +
                                       fe->system_to_component_index(i).first];
      return this->finite_element_output.shape_values(row, j);
    }
}



template <int dim, int spacedim>
inline double
FEValuesBase<dim, spacedim>::shape_value_component(
  const unsigned int i,
  const unsigned int j,
  const unsigned int component) const
{
  AssertIndexRange(i, fe->n_dofs_per_cell());
  Assert(this->update_flags & update_values,
         ExcAccessToUninitializedField("update_values"));
  AssertIndexRange(component, fe->n_components());
  Assert(present_cell.get() != nullptr,
         ExcMessage("FEValues object is not reinit'ed to any cell"));

  // check whether the shape function
  // is non-zero at all within
  // this component:
  if (fe->get_nonzero_components(i)[component] == false)
    return 0;

  // look up the right row in the
  // table and take the data from
  // there
  const unsigned int row =
    this->finite_element_output
      .shape_function_to_row_table[i * fe->n_components() + component];
  return this->finite_element_output.shape_values(row, j);
}



template <int dim, int spacedim>
inline const Tensor<1, spacedim> &
FEValuesBase<dim, spacedim>::shape_grad(const unsigned int i,
                                        const unsigned int j) const
{
  AssertIndexRange(i, fe->n_dofs_per_cell());
  Assert(this->update_flags & update_gradients,
         ExcAccessToUninitializedField("update_gradients"));
  Assert(fe->is_primitive(i), ExcShapeFunctionNotPrimitive(i));
  Assert(present_cell.get() != nullptr,
         ExcMessage("FEValues object is not reinit'ed to any cell"));
  // if the entire FE is primitive,
  // then we can take a short-cut:
  if (fe->is_primitive())
    return this->finite_element_output.shape_gradients[i][j];
  else
    {
      // otherwise, use the mapping
      // between shape function
      // numbers and rows. note that
      // by the assertions above, we
      // know that this particular
      // shape function is primitive,
      // so we can call
      // system_to_component_index
      const unsigned int row =
        this->finite_element_output
          .shape_function_to_row_table[i * fe->n_components() +
                                       fe->system_to_component_index(i).first];
      return this->finite_element_output.shape_gradients[row][j];
    }
}



template <int dim, int spacedim>
inline Tensor<1, spacedim>
FEValuesBase<dim, spacedim>::shape_grad_component(
  const unsigned int i,
  const unsigned int j,
  const unsigned int component) const
{
  AssertIndexRange(i, fe->n_dofs_per_cell());
  Assert(this->update_flags & update_gradients,
         ExcAccessToUninitializedField("update_gradients"));
  AssertIndexRange(component, fe->n_components());
  Assert(present_cell.get() != nullptr,
         ExcMessage("FEValues object is not reinit'ed to any cell"));
  // check whether the shape function
  // is non-zero at all within
  // this component:
  if (fe->get_nonzero_components(i)[component] == false)
    return Tensor<1, spacedim>();

  // look up the right row in the
  // table and take the data from
  // there
  const unsigned int row =
    this->finite_element_output
      .shape_function_to_row_table[i * fe->n_components() + component];
  return this->finite_element_output.shape_gradients[row][j];
}



template <int dim, int spacedim>
inline const Tensor<2, spacedim> &
FEValuesBase<dim, spacedim>::shape_hessian(const unsigned int i,
                                           const unsigned int j) const
{
  AssertIndexRange(i, fe->n_dofs_per_cell());
  Assert(this->update_flags & update_hessians,
         ExcAccessToUninitializedField("update_hessians"));
  Assert(fe->is_primitive(i), ExcShapeFunctionNotPrimitive(i));
  Assert(present_cell.get() != nullptr,
         ExcMessage("FEValues object is not reinit'ed to any cell"));
  // if the entire FE is primitive,
  // then we can take a short-cut:
  if (fe->is_primitive())
    return this->finite_element_output.shape_hessians[i][j];
  else
    {
      // otherwise, use the mapping
      // between shape function
      // numbers and rows. note that
      // by the assertions above, we
      // know that this particular
      // shape function is primitive,
      // so we can call
      // system_to_component_index
      const unsigned int row =
        this->finite_element_output
          .shape_function_to_row_table[i * fe->n_components() +
                                       fe->system_to_component_index(i).first];
      return this->finite_element_output.shape_hessians[row][j];
    }
}



template <int dim, int spacedim>
inline Tensor<2, spacedim>
FEValuesBase<dim, spacedim>::shape_hessian_component(
  const unsigned int i,
  const unsigned int j,
  const unsigned int component) const
{
  AssertIndexRange(i, fe->n_dofs_per_cell());
  Assert(this->update_flags & update_hessians,
         ExcAccessToUninitializedField("update_hessians"));
  AssertIndexRange(component, fe->n_components());
  Assert(present_cell.get() != nullptr,
         ExcMessage("FEValues object is not reinit'ed to any cell"));
  // check whether the shape function
  // is non-zero at all within
  // this component:
  if (fe->get_nonzero_components(i)[component] == false)
    return Tensor<2, spacedim>();

  // look up the right row in the
  // table and take the data from
  // there
  const unsigned int row =
    this->finite_element_output
      .shape_function_to_row_table[i * fe->n_components() + component];
  return this->finite_element_output.shape_hessians[row][j];
}



template <int dim, int spacedim>
inline const Tensor<3, spacedim> &
FEValuesBase<dim, spacedim>::shape_3rd_derivative(const unsigned int i,
                                                  const unsigned int j) const
{
  AssertIndexRange(i, fe->n_dofs_per_cell());
  Assert(this->update_flags & update_3rd_derivatives,
         ExcAccessToUninitializedField("update_3rd_derivatives"));
  Assert(fe->is_primitive(i), ExcShapeFunctionNotPrimitive(i));
  Assert(present_cell.get() != nullptr,
         ExcMessage("FEValues object is not reinit'ed to any cell"));
  // if the entire FE is primitive,
  // then we can take a short-cut:
  if (fe->is_primitive())
    return this->finite_element_output.shape_3rd_derivatives[i][j];
  else
    {
      // otherwise, use the mapping
      // between shape function
      // numbers and rows. note that
      // by the assertions above, we
      // know that this particular
      // shape function is primitive,
      // so we can call
      // system_to_component_index
      const unsigned int row =
        this->finite_element_output
          .shape_function_to_row_table[i * fe->n_components() +
                                       fe->system_to_component_index(i).first];
      return this->finite_element_output.shape_3rd_derivatives[row][j];
    }
}



template <int dim, int spacedim>
inline Tensor<3, spacedim>
FEValuesBase<dim, spacedim>::shape_3rd_derivative_component(
  const unsigned int i,
  const unsigned int j,
  const unsigned int component) const
{
  AssertIndexRange(i, fe->n_dofs_per_cell());
  Assert(this->update_flags & update_3rd_derivatives,
         ExcAccessToUninitializedField("update_3rd_derivatives"));
  AssertIndexRange(component, fe->n_components());
  Assert(present_cell.get() != nullptr,
         ExcMessage("FEValues object is not reinit'ed to any cell"));
  // check whether the shape function
  // is non-zero at all within
  // this component:
  if (fe->get_nonzero_components(i)[component] == false)
    return Tensor<3, spacedim>();

  // look up the right row in the
  // table and take the data from
  // there
  const unsigned int row =
    this->finite_element_output
      .shape_function_to_row_table[i * fe->n_components() + component];
  return this->finite_element_output.shape_3rd_derivatives[row][j];
}



template <int dim, int spacedim>
inline const FiniteElement<dim, spacedim> &
FEValuesBase<dim, spacedim>::get_fe() const
{
  return *fe;
}



template <int dim, int spacedim>
inline const Mapping<dim, spacedim> &
FEValuesBase<dim, spacedim>::get_mapping() const
{
  return *mapping;
}



template <int dim, int spacedim>
inline UpdateFlags
FEValuesBase<dim, spacedim>::get_update_flags() const
{
  return this->update_flags;
}



template <int dim, int spacedim>
inline const std::vector<Point<spacedim>> &
FEValuesBase<dim, spacedim>::get_quadrature_points() const
{
  Assert(this->update_flags & update_quadrature_points,
         ExcAccessToUninitializedField("update_quadrature_points"));
  Assert(present_cell.get() != nullptr,
         ExcMessage("FEValues object is not reinit'ed to any cell"));
  return this->mapping_output.quadrature_points;
}



template <int dim, int spacedim>
inline const std::vector<double> &
FEValuesBase<dim, spacedim>::get_JxW_values() const
{
  Assert(this->update_flags & update_JxW_values,
         ExcAccessToUninitializedField("update_JxW_values"));
  Assert(present_cell.get() != nullptr,
         ExcMessage("FEValues object is not reinit'ed to any cell"));
  return this->mapping_output.JxW_values;
}



template <int dim, int spacedim>
inline const std::vector<DerivativeForm<1, dim, spacedim>> &
FEValuesBase<dim, spacedim>::get_jacobians() const
{
  Assert(this->update_flags & update_jacobians,
         ExcAccessToUninitializedField("update_jacobians"));
  Assert(present_cell.get() != nullptr,
         ExcMessage("FEValues object is not reinit'ed to any cell"));
  return this->mapping_output.jacobians;
}



template <int dim, int spacedim>
inline const std::vector<DerivativeForm<2, dim, spacedim>> &
FEValuesBase<dim, spacedim>::get_jacobian_grads() const
{
  Assert(this->update_flags & update_jacobian_grads,
         ExcAccessToUninitializedField("update_jacobians_grads"));
  Assert(present_cell.get() != nullptr,
         ExcMessage("FEValues object is not reinit'ed to any cell"));
  return this->mapping_output.jacobian_grads;
}



template <int dim, int spacedim>
inline const Tensor<3, spacedim> &
FEValuesBase<dim, spacedim>::jacobian_pushed_forward_grad(
  const unsigned int i) const
{
  Assert(this->update_flags & update_jacobian_pushed_forward_grads,
         ExcAccessToUninitializedField("update_jacobian_pushed_forward_grads"));
  Assert(present_cell.get() != nullptr,
         ExcMessage("FEValues object is not reinit'ed to any cell"));
  return this->mapping_output.jacobian_pushed_forward_grads[i];
}



template <int dim, int spacedim>
inline const std::vector<Tensor<3, spacedim>> &
FEValuesBase<dim, spacedim>::get_jacobian_pushed_forward_grads() const
{
  Assert(this->update_flags & update_jacobian_pushed_forward_grads,
         ExcAccessToUninitializedField("update_jacobian_pushed_forward_grads"));
  Assert(present_cell.get() != nullptr,
         ExcMessage("FEValues object is not reinit'ed to any cell"));
  return this->mapping_output.jacobian_pushed_forward_grads;
}



template <int dim, int spacedim>
inline const DerivativeForm<3, dim, spacedim> &
FEValuesBase<dim, spacedim>::jacobian_2nd_derivative(const unsigned int i) const
{
  Assert(this->update_flags & update_jacobian_2nd_derivatives,
         ExcAccessToUninitializedField("update_jacobian_2nd_derivatives"));
  Assert(present_cell.get() != nullptr,
         ExcMessage("FEValues object is not reinit'ed to any cell"));
  return this->mapping_output.jacobian_2nd_derivatives[i];
}



template <int dim, int spacedim>
inline const std::vector<DerivativeForm<3, dim, spacedim>> &
FEValuesBase<dim, spacedim>::get_jacobian_2nd_derivatives() const
{
  Assert(this->update_flags & update_jacobian_2nd_derivatives,
         ExcAccessToUninitializedField("update_jacobian_2nd_derivatives"));
  Assert(present_cell.get() != nullptr,
         ExcMessage("FEValues object is not reinit'ed to any cell"));
  return this->mapping_output.jacobian_2nd_derivatives;
}



template <int dim, int spacedim>
inline const Tensor<4, spacedim> &
FEValuesBase<dim, spacedim>::jacobian_pushed_forward_2nd_derivative(
  const unsigned int i) const
{
  Assert(this->update_flags & update_jacobian_pushed_forward_2nd_derivatives,
         ExcAccessToUninitializedField(
           "update_jacobian_pushed_forward_2nd_derivatives"));
  Assert(present_cell.get() != nullptr,
         ExcMessage("FEValues object is not reinit'ed to any cell"));
  return this->mapping_output.jacobian_pushed_forward_2nd_derivatives[i];
}



template <int dim, int spacedim>
inline const std::vector<Tensor<4, spacedim>> &
FEValuesBase<dim, spacedim>::get_jacobian_pushed_forward_2nd_derivatives() const
{
  Assert(this->update_flags & update_jacobian_pushed_forward_2nd_derivatives,
         ExcAccessToUninitializedField(
           "update_jacobian_pushed_forward_2nd_derivatives"));
  Assert(present_cell.get() != nullptr,
         ExcMessage("FEValues object is not reinit'ed to any cell"));
  return this->mapping_output.jacobian_pushed_forward_2nd_derivatives;
}



template <int dim, int spacedim>
inline const DerivativeForm<4, dim, spacedim> &
FEValuesBase<dim, spacedim>::jacobian_3rd_derivative(const unsigned int i) const
{
  Assert(this->update_flags & update_jacobian_3rd_derivatives,
         ExcAccessToUninitializedField("update_jacobian_3rd_derivatives"));
  Assert(present_cell.get() != nullptr,
         ExcMessage("FEValues object is not reinit'ed to any cell"));
  return this->mapping_output.jacobian_3rd_derivatives[i];
}



template <int dim, int spacedim>
inline const std::vector<DerivativeForm<4, dim, spacedim>> &
FEValuesBase<dim, spacedim>::get_jacobian_3rd_derivatives() const
{
  Assert(this->update_flags & update_jacobian_3rd_derivatives,
         ExcAccessToUninitializedField("update_jacobian_3rd_derivatives"));
  Assert(present_cell.get() != nullptr,
         ExcMessage("FEValues object is not reinit'ed to any cell"));
  return this->mapping_output.jacobian_3rd_derivatives;
}



template <int dim, int spacedim>
inline const Tensor<5, spacedim> &
FEValuesBase<dim, spacedim>::jacobian_pushed_forward_3rd_derivative(
  const unsigned int i) const
{
  Assert(this->update_flags & update_jacobian_pushed_forward_3rd_derivatives,
         ExcAccessToUninitializedField(
           "update_jacobian_pushed_forward_3rd_derivatives"));
  Assert(present_cell.get() != nullptr,
         ExcMessage("FEValues object is not reinit'ed to any cell"));
  return this->mapping_output.jacobian_pushed_forward_3rd_derivatives[i];
}



template <int dim, int spacedim>
inline const std::vector<Tensor<5, spacedim>> &
FEValuesBase<dim, spacedim>::get_jacobian_pushed_forward_3rd_derivatives() const
{
  Assert(this->update_flags & update_jacobian_pushed_forward_3rd_derivatives,
         ExcAccessToUninitializedField(
           "update_jacobian_pushed_forward_3rd_derivatives"));
  Assert(present_cell.get() != nullptr,
         ExcMessage("FEValues object is not reinit'ed to any cell"));
  return this->mapping_output.jacobian_pushed_forward_3rd_derivatives;
}



template <int dim, int spacedim>
inline const std::vector<DerivativeForm<1, spacedim, dim>> &
FEValuesBase<dim, spacedim>::get_inverse_jacobians() const
{
  Assert(this->update_flags & update_inverse_jacobians,
         ExcAccessToUninitializedField("update_inverse_jacobians"));
  Assert(present_cell.get() != nullptr,
         ExcMessage("FEValues object is not reinit'ed to any cell"));
  return this->mapping_output.inverse_jacobians;
}



template <int dim, int spacedim>
inline std_cxx20::ranges::iota_view<unsigned int, unsigned int>
FEValuesBase<dim, spacedim>::dof_indices() const
{
  return {0U, dofs_per_cell};
}



template <int dim, int spacedim>
inline std_cxx20::ranges::iota_view<unsigned int, unsigned int>
FEValuesBase<dim, spacedim>::dof_indices_starting_at(
  const unsigned int start_dof_index) const
{
  Assert(start_dof_index <= dofs_per_cell,
         ExcIndexRange(start_dof_index, 0, dofs_per_cell + 1));
  return {start_dof_index, dofs_per_cell};
}



template <int dim, int spacedim>
inline std_cxx20::ranges::iota_view<unsigned int, unsigned int>
FEValuesBase<dim, spacedim>::dof_indices_ending_at(
  const unsigned int end_dof_index) const
{
  Assert(end_dof_index < dofs_per_cell,
         ExcIndexRange(end_dof_index, 0, dofs_per_cell));
  return {0U, end_dof_index + 1};
}



template <int dim, int spacedim>
inline std_cxx20::ranges::iota_view<unsigned int, unsigned int>
FEValuesBase<dim, spacedim>::quadrature_point_indices() const
{
  return {0U, n_quadrature_points};
}



template <int dim, int spacedim>
inline const Point<spacedim> &
FEValuesBase<dim, spacedim>::quadrature_point(const unsigned int i) const
{
  Assert(this->update_flags & update_quadrature_points,
         ExcAccessToUninitializedField("update_quadrature_points"));
  AssertIndexRange(i, this->mapping_output.quadrature_points.size());
  Assert(present_cell.get() != nullptr,
         ExcMessage("FEValues object is not reinit'ed to any cell"));

  return this->mapping_output.quadrature_points[i];
}



template <int dim, int spacedim>
inline double
FEValuesBase<dim, spacedim>::JxW(const unsigned int i) const
{
  Assert(this->update_flags & update_JxW_values,
         ExcAccessToUninitializedField("update_JxW_values"));
  AssertIndexRange(i, this->mapping_output.JxW_values.size());
  Assert(present_cell.get() != nullptr,
         ExcMessage("FEValues object is not reinit'ed to any cell"));

  return this->mapping_output.JxW_values[i];
}



template <int dim, int spacedim>
inline const DerivativeForm<1, dim, spacedim> &
FEValuesBase<dim, spacedim>::jacobian(const unsigned int i) const
{
  Assert(this->update_flags & update_jacobians,
         ExcAccessToUninitializedField("update_jacobians"));
  AssertIndexRange(i, this->mapping_output.jacobians.size());
  Assert(present_cell.get() != nullptr,
         ExcMessage("FEValues object is not reinit'ed to any cell"));

  return this->mapping_output.jacobians[i];
}



template <int dim, int spacedim>
inline const DerivativeForm<2, dim, spacedim> &
FEValuesBase<dim, spacedim>::jacobian_grad(const unsigned int i) const
{
  Assert(this->update_flags & update_jacobian_grads,
         ExcAccessToUninitializedField("update_jacobians_grads"));
  AssertIndexRange(i, this->mapping_output.jacobian_grads.size());
  Assert(present_cell.get() != nullptr,
         ExcMessage("FEValues object is not reinit'ed to any cell"));

  return this->mapping_output.jacobian_grads[i];
}



template <int dim, int spacedim>
inline const DerivativeForm<1, spacedim, dim> &
FEValuesBase<dim, spacedim>::inverse_jacobian(const unsigned int i) const
{
  Assert(this->update_flags & update_inverse_jacobians,
         ExcAccessToUninitializedField("update_inverse_jacobians"));
  AssertIndexRange(i, this->mapping_output.inverse_jacobians.size());
  Assert(present_cell.get() != nullptr,
         ExcMessage("FEValues object is not reinit'ed to any cell"));

  return this->mapping_output.inverse_jacobians[i];
}



template <int dim, int spacedim>
inline const Tensor<1, spacedim> &
FEValuesBase<dim, spacedim>::normal_vector(const unsigned int i) const
{
  Assert(this->update_flags & update_normal_vectors,
         (typename FEValuesBase<dim, spacedim>::ExcAccessToUninitializedField(
           "update_normal_vectors")));
  AssertIndexRange(i, this->mapping_output.normal_vectors.size());
  Assert(present_cell.get() != nullptr,
         ExcMessage("FEValues object is not reinit'ed to any cell"));

  return this->mapping_output.normal_vectors[i];
}



 /*--------------------- Inline functions: FEValues --------------------------*/ 


template <int dim, int spacedim>
inline const Quadrature<dim> &
FEValues<dim, spacedim>::get_quadrature() const
{
  return quadrature;
}



template <int dim, int spacedim>
inline const FEValues<dim, spacedim> &
FEValues<dim, spacedim>::get_present_fe_values() const
{
  return *this;
}


 /*---------------------- Inline functions: FEFaceValuesBase -----------------*/ 


template <int dim, int spacedim>
inline unsigned int
FEFaceValuesBase<dim, spacedim>::get_face_index() const
{
  return present_face_index;
}


 /*----------------------- Inline functions: FE*FaceValues -------------------*/ 

template <int dim, int spacedim>
inline const Quadrature<dim - 1> &
FEFaceValuesBase<dim, spacedim>::get_quadrature() const
{
  return quadrature[quadrature.size() == 1 ? 0 : present_face_no];
}



template <int dim, int spacedim>
inline const FEFaceValues<dim, spacedim> &
FEFaceValues<dim, spacedim>::get_present_fe_values() const
{
  return *this;
}



template <int dim, int spacedim>
inline const FESubfaceValues<dim, spacedim> &
FESubfaceValues<dim, spacedim>::get_present_fe_values() const
{
  return *this;
}



template <int dim, int spacedim>
inline const Tensor<1, spacedim> &
FEFaceValuesBase<dim, spacedim>::boundary_form(const unsigned int i) const
{
  AssertIndexRange(i, this->mapping_output.boundary_forms.size());
  Assert(this->update_flags & update_boundary_forms,
         (typename FEValuesBase<dim, spacedim>::ExcAccessToUninitializedField(
           "update_boundary_forms")));

  return this->mapping_output.boundary_forms[i];
}

#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif


