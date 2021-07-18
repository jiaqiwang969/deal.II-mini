//include/deal.II-translator/numerics/data_postprocessor_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2007 - 2020 by the deal.II authors
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

#ifndef dealii_data_postprocessor_h
#define dealii_data_postprocessor_h



#include <deal.II/base/config.h>

#include <deal.II/base/point.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/base/tensor.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_update_flags.h>

#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_component_interpretation.h>

#include <boost/any.hpp>

#include <string>
#include <vector>

DEAL_II_NAMESPACE_OPEN


/**
 * 一个用于数据结构的命名空间，这些数据结构将从DataOut传递给DataPostprocessor的成员函数。
 *
 *
 */
namespace DataPostprocessorInputs
{
  /**
   * 一个包含Scalar和Vector类共同元素的基类，这些元素被作为参数传递给
   * DataPostprocessor::evaluate_scalar_field() 和
   * DataPostprocessor::evaluate_vector_field().
   * 这个共同的基类提供了对正在评估解决方案的点的访问，以及其他一些领域，如下所述。
   * <h4>Normal vector access</h4>
   * 如果合适的话，也就是说，如果当前正在处理的对象是一个单元格的面，并且DataPostprocessor对象是从DataOutFaces或类似的类中调用的，那么当前对象也会在这些评价点上存储到生成输出的几何体的法向量。
   * 另一方面，如果解决方案是在一个单元格上被评估的，那么
   * @p normal_vectors 成员变量就不包含任何有用的东西。
   * <h4>Cell access</h4>
   * DataPostprocessor通常由DataOut或DataOutFaces等类调用，这些类以单元格为基础评估解决方案字段。因此，从DataPostprocessor（或DataPostprocessorScalar、DataPostprocessorVector或DataPostprocessorTensor）派生的类有时需要使用当前正在调查的单元。因此，DataOut和类似的类通过这个命名空间的类（特别是当前的基类）将它们当前正在处理的单元传递给DataPostprocessor。
   * 然而，情况并不那么简单。这是因为当前类（以及从它派生的类）只知道输出所处的空间维度。但这可能来自于很多方面。例如，如果我们是在三维空间，这可能是因为我们是在DoFHandler<3>或DoFHandler<2,3>上工作（也就是说，要么是三维网格，要么是嵌入三维空间的二维表面的二维网格）。另一种情况是像DataOutRotation或DataOutStack这样的类，那么
   * @p spacedim
   * 等于3可能意味着我们实际上是在做一个DoFHandler<2>。
   * 换句话说，仅仅因为我们知道当前类的 @p spacedim
   * 模板参数的值，并不意味着当前正在处理的单元格迭代器的数据类型是明显的。
   * 尽管如此，为了使单元格迭代器可以被访问，这个类使用一个
   * boost::any
   * 类型的对象来存储单元格迭代器。你可以把它看成是一个可以指向任何东西的无效指针。
   * 因此要使用被使用的东西需要用户知道被指向的东西的数据类型。
   * 为了使其发挥作用，DataOut和相关类在当前类型的对象中存储了一个单元格的表示。为了把它拿回来，你将使用get_cell()函数，该函数要求你说，作为一个模板参数，当前正在处理的单元格的尺寸。这是你在应用程序中通常拥有的知识：例如，如果你的应用程序在
   * @p dim  空间维度中运行，并且你目前正在使用 DataOut
   * 类，那么被处理的单元格具有数据类型
   * <code>DataOut<dim>::cell_iterator</code>  。
   * 因此，在后处理程序中，你可以调用<code>inputs.get_cell
   * @<dim@>
   * </code>。然而，由于技术原因，C++通常会要求你把它写成
   * <code>inputs.template get_cell@<dim@> </code>
   * ，因为我们在这里调用的成员函数要求我们明确提供模板参数。
   * 让我们考虑一个完整的后处理器的例子，从粘度 $\eta$
   * 和流体速度梯度 $\nabla u$ 计算应力 $\|\sigma\| = \|\eta \nabla
   * u\|$
   * 的流体规范，假设粘度是取决于单元的材料id的东西。这可以通过我们从DataPostprocessorScalar派生出来的一个类来完成，在这个类中，我们重载了
   * DataPostprocessor::evaluate_vector_field()
   * 函数，该函数接收速度的值和梯度（还有其他求解变量，如压力，但我们暂时忽略这些）。然后我们可以使用这样的代码。
   * @code
   * template <int dim>
   * class ComputeStress : public DataPostprocessorScalar<dim>
   * {
   *   public:
   *     ... // overload other necessary member variables
   *     virtual
   *     void
   *     evaluate_vector_field
   *     (const DataPostprocessorInputs::Vector<dim> &input_data,
   *      std::vector<Vector<double> > &computed_quantities) const override
   *     {
   *       const typename DoFHandler<dim>::cell_iterator current_cell =
   *         input_data.template get_cell<dim>();
   *       const viscosity = look_up_viscosity (current_cell->material_id());
   *
   *       for (unsigned int q=0; q<input_data.solution_gradients.size(); ++q)
   *         computed_quantities[q][0] =
   *           (viscosity input_data.solution_gradients[q]).norm();
   *     }
   * };
   * @endcode
   *
   *
   */
  template <int spacedim>
  struct CommonInputs
  {
    /**
     * 在我们生成图形输出的点上，对单元格面的法向数组进行评估。这个数组只被DataOutFaces类使用，而其他所有可以使用DataPostprocessor框架的类都留空。在DataOutFaces的情况下，这个数组包含了从单元内部看到的面的外向法向量。
     * 这个数组只有在用户派生的类重载了
     * DataPostprocessor::get_needed_update_flags(),
     * 并且函数返回（可能还有其他标志）
     * UpdateFlags::update_normal_vectors.
     * 时才会被填充。另外，从DataPostprocessorScalar、DataPostprocessorVector或DataPostprocessorTensor派生的类可以将这个标志传递给这三个类的构造器。
     *
     */
    std::vector<Tensor<1, spacedim>> normals;

    /**
     * 一个坐标数组，对应于我们在一个单元上生成图形输出的位置。
     * 这个数组只有在用户派生的类重载了
     * DataPostprocessor::get_needed_update_flags(),
     * 并且函数返回（可能还有其他标志）
     * UpdateFlags::update_quadrature_points.
     * 时才会被填充。另外，从DataPostprocessorScalar、DataPostprocessorVector或DataPostprocessorTensor派生的类可以将这个标志传递给这三个类的构造器。
     *
     */
    std::vector<Point<spacedim>> evaluation_points;

    /**
     * 设置当前用于评估DataPostprocessor对象被调用的数据的单元。
     * 这个函数通常不从用户空间调用，而是由DataOut和类似的类在创建对象时调用，然后传递给DataPostprocessor。
     *
     */
    template <int dim>
    void
    set_cell(const typename DoFHandler<dim, spacedim>::cell_iterator &cell);

    /**
     * 设置当前用于评估DataPostprocessor对象被调用的数据的单元。
     * 这个函数通常不从用户空间调用，而是由DataOut和类似的类在创建对象时调用，然后传递给DataPostprocessor。
     * @deprecated  使用带有dim模板参数的等效函数来代替。
     *
     */
    template <typename DoFHandlerType>
    DEAL_II_DEPRECATED void
    set_cell(const typename DoFHandlerType::cell_iterator &cell);

    /**
     * 查询我们当前产生图形输出的单元格。
     * 关于如何使用这个函数的例子，请参见当前类的文档。
     *
     */
    template <int dim>
    typename DoFHandler<dim, spacedim>::cell_iterator
    get_cell() const;

    /**
     * 查询我们当前产生图形输出的单元格。
     * 关于如何使用这个函数的例子，请参见当前类的文档。
     * @deprecated 使用带有dim模板参数的等效函数来代替。
     *
     */
    template <typename DoFHandlerType>
    DEAL_II_DEPRECATED typename DoFHandlerType::cell_iterator
    get_cell() const;

  private:
    /**
     * set_cell()存储单元格的地方。由于单元格迭代器的实际数据类型可以是许多不同的东西，接口在这里使用
     * boost::any
     * 。这使得set_cell()中的赋值很简单，但需要知道get_cell()中存储对象的数据类型。
     *
     */
    boost::any cell;
  };

  /**
   * 一个用于向 DataPostprocessor::evaluate_scalar_field().
   * 传递信息的结构，它包含了标量解变量在单元格或面的评估点的值和（如果要求）导数。然而，如果标量解是复值的，则不使用该类，因为在这种情况下，实部和虚部是分开处理的）。
   *
   * - 导致数据后处理程序的矢量值输入，然后被传递到 DataPostprocessor::evaluate_vector_field() 。    通过CommonInputs基类中的字段，该类也可以访问评估点的位置、法向量（如果合适），以及当前正在评估的单元数据（同样如果合适）。
   *
   */
  template <int spacedim>
  struct Scalar : public CommonInputs<spacedim>
  {
    /**
     * 每个评估点的（标量）解决方案的数值数组，用于从一个单元格、面或其他对象创建图形输出。
     *
     */
    std::vector<double> solution_values;

    /**
     * 在每个评估点的（标量）解决方案的梯度数组，用于从一个单元格、面或其他对象创建图形输出。
     * 这个数组只有在用户派生的类重载了
     * DataPostprocessor::get_needed_update_flags(),
     * 并且函数返回（可能还有其他标志）
     * UpdateFlags::update_gradients.
     * 时才会被填充。另外，从DataPostprocessorScalar、DataPostprocessorVector或DataPostprocessorTensor派生的类可以将这个标志传递给这三个类的构造器。
     *
     */
    std::vector<Tensor<1, spacedim>> solution_gradients;

    /**
     * 在每个评估点的（标量）解决方案的二阶导数数组，用于从一个单元格、面或其他对象创建图形输出。
     * 这个数组只有在用户派生的类重载了
     * DataPostprocessor::get_needed_update_flags(),
     * 并且函数返回（可能还有其他标志）
     * UpdateFlags::update_hessians.
     * 时才会被填充。另外，从DataPostprocessorScalar、DataPostprocessorVector或DataPostprocessorTensor派生的类可以将这个标志传递给这三个类的构造器。
     *
     */
    std::vector<Tensor<2, spacedim>> solution_hessians;
  };



  /**
   * 一个用于向 DataPostprocessor::evaluate_vector_field().
   * 传递信息的结构，它包含向量值解变量在单元格或面的评估点的值和（如果要求）导数。
   * 如果求解向量是复数值的，也可以使用这个类（在这种情况下，它是标量值还是向量值并不重要），因为在这种情况下，DataOut和相关类将求解向量的实部和虚部分开。在实践中，这意味着如果一个解向量有
   * $N$ 个向量分量（即，有 $N$
   * 个函数构成你所处理的PDE的解； $N$
   * 不是解向量的大小），那么如果解是实值的，下面的`solution_values`变量将是一个数组，有多少个条目就有多少个单元格上的评估点，每个条目是一个长度为
   * $N$ 的向量，代表在一个点上评估的 $N$
   * 解函数。另一方面，如果解是复值的（即传递给
   * DataOut::build_patches()
   * 的向量有复值的条目），那么这个类的`solution_values`成员变量将有每个评价点的
   * $2N$ 条目。这些条目中的第一个 $N$
   * 代表解决方案的实部，第二个 $N$
   * 条目对应于在评估点评估的解决方案的虚部。同样的布局被用于
   * "solution_gradients "和 "solution_hessians
   * "字段。首先是实部的梯度/黑森数，然后是虚部的所有梯度/黑森数。在DataPostprocessor类本身的文档中有更多关于这个主题的信息。
   * step-58
   * 提供了一个如何在复值情况下使用这个类的例子。
   * 通过CommonInputs基类中的字段，该类也可以访问评估点的位置、法向量（如果合适），以及当前正在评估的单元数据（同样如果合适）。
   *
   */
  template <int spacedim>
  struct Vector : public CommonInputs<spacedim>
  {
    /**
     * 用于从一个单元格、面或其他对象创建图形输出的每个评价点的向量值解决方案的数组。
     * 外向量在评估点上运行，而内向量在将被生成输出的有限元场的分量上运行。
     *
     */
    std::vector<dealii::Vector<double>> solution_values;

    /**
     * 在每个评估点的矢量值解决方案的梯度阵列，用于从一个单元格、面或其他对象创建图形输出。
     * 外向量在评估点上运行，而内向量在将被生成输出的有限元场的分量上运行。
     * 这个数组只有在用户派生的类重载了
     * DataPostprocessor::get_needed_update_flags(),
     * 并且函数返回（可能还有其他标志）
     * UpdateFlags::update_gradients.
     * 时才会被填充。另外，从DataPostprocessorScalar、DataPostprocessorVector或DataPostprocessorTensor派生的类可以将这个标志传递给这三个类的构造器。
     *
     */
    std::vector<std::vector<Tensor<1, spacedim>>> solution_gradients;

    /**
     * 矢量值解决方案在每个评估点的二阶导数数组，用于从一个单元格、面或其他对象创建图形输出。
     * 外向量在评估点上运行，而内向量在将被生成输出的有限元场的分量上运行。
     * 这个数组只有在用户派生的类重载了
     * DataPostprocessor::get_needed_update_flags(),
     * 并且函数返回（可能还有其他标志）
     * UpdateFlags::update_hessians.
     * 时才会被填充。另外，从DataPostprocessorScalar、DataPostprocessorVector或DataPostprocessorTensor派生的类可以将这个标志传递给这三个类的构造器。
     *
     */
    std::vector<std::vector<Tensor<2, spacedim>>> solution_hessians;
  };

} // namespace DataPostprocessorInputs


/**
 * 该类提供了一个接口，用于计算来自解决方案的派生量，然后可以使用诸如DataOut类的设施，以图形格式输出，用于可视化。
 * 对于FE解决方案的（图形）输出，人们经常希望包括派生量，这些派生量是由解决方案的值和可能的解决方案的第一和第二导数计算出来的。例如，在超声速流动计算中，根据速度和密度计算马赫数，或者如
 * step-29 和 step-58
 * 中所示，计算复值解的幅度（实际上是幅度的平方）。其他的用途在
 * step-32  和  step-33
 * 中显示。这个类提供了执行这种后处理的接口。给出解的值和导数，在我们想要生成输出的那些点上，这个类的函数可以被重载以计算新的数量。
 * 可以给 DataOut::add_data_vector()
 * 函数提供一个数据向量和一个从当前类派生出来的对象（对DataOutRotation和DataOutFaces也是如此）。这将导致
 * DataOut::build_patches()
 * 计算派生量，而不是使用数据向量提供的数据（通常是解向量）。注意，DataPostprocessor对象（即实际上是你的派生类的对象）必须活到DataOut对象被销毁为止，因为后者保持着一个指向前者的指针，如果指向的对象被销毁，而后者仍有一个指向它的指针，它就会抱怨。如果数据后处理器和DataOut对象都是一个函数的局部变量（例如，在
 * step-29
 * 中就是如此），那么你可以通过在DataOut变量之前声明数据后处理器变量来避免这个错误，因为对象的销毁顺序是相反的。
 * 为了不进行不必要的计算，DataPostprocessor必须提供计算派生量所需的输入数据的信息，即它是否需要所提供数据的值、一阶导数和/或二阶导数。与DataOutFaces对象结合使用的DataPostprocessor对象也可以要求提供每个点的法向量。需要哪些数据的信息必须通过虚拟函数get_need_update_flags()返回的UpdateFlags来提供。你有责任在计算派生量时只使用那些被更新的值。DataOut对象将在调用evaluate_scalar_field()或evaluate_vector_field()时提供对请求数据的引用（DataOut决定调用这两个函数中的哪一个，取决于所使用的有限元是只有一个，还是有多个矢量成分；注意，这只由所使用的有限元的成分数量决定，而不是由当前派生类计算的数据是标量还是矢量值）。
 * 此外，派生类必须实现get_names()函数，后者函数返回的输出变量的数量必须与前者返回的向量的大小相匹配。此外，这个数字还必须与计算量的数量相匹配，当然了。
 *
 *  <h3>Use in simpler cases</h3>
 * 从当前的类派生出来，可以实现非常普遍的后处理程序。例如，在
 * step-32
 * 程序中，我们实现了一个后处理器，它接收一个由速度、压力和温度（dim+2分量）组成的解决方案，并计算出各种输出量，其中一些是矢量值，一些是标量。另一方面，在
 * step-29
 * 中，我们实现了一个后处理器，只计算由双分量有限元给出的复数的大小。为此必须实现四个虚拟函数（evaluate_scalar_field()或evaluate_vector_field()、get_names()、get_update_flags()和get_data_component_interpretation()）似乎很愚蠢。
 * 为此，有三个类DataPostprocessorScalar、DataPostprocessorVector和DataPostprocessorTensor，如果输出量是一个单一的标量、一个单一的向量（这里用来指正好有
 * @p dim 个分量）或一个单一的张量（这里用来指正好有
 * <code>dim*dim</code>
 * 个分量），就可以使用它们。当使用这些类时，人们只需要编写一个构造函数，将输出变量的名称和更新标志传递给基类的构造函数，并重载实际计算结果的函数。
 * DataPostprocessorVector和DataPostprocessorTensor类的文档也包含了大量关于如何使用它们的例子。
 * step-29
 * 教程程序包含了一个使用DataPostprocessorScalar类的例子。
 *
 *  <h3>Complex-valued solutions</h3>
 * 有一些PDEs的解是复值的。例如， step-58 和 step-62
 * 解决的问题，其每一点的解都是由一个 `std::complex<double>`
 * 变量代表的复数组成。(  step-29
 * 也解决了这样的问题，但在那里我们选择用两个实值场来表示解。)
 * 在这种情况下，交给 DataOut::build_patches() 的向量是
 * `Vector<std::complex<double>>`,
 * 类型或基本等同于此的东西。这方面的问题，在DataOut本身的文档中也讨论过，最广泛使用的可视化文件格式（特别是VTK和VTU格式）实际上不能表示复数量。在这些数据文件中唯一可以存储的是实值量。
 * 因此，DataOut被迫将事物拆成实数和虚数部分，并将两者作为单独的量输出。这是由DataOut直接写入文件的数据的情况，但它也是首先通过DataPostprocessor对象（或其派生类的对象）的情况。这些对象所看到的都是实值的集合，即使底层的解决方案向量是复值的。
 * 所有这些都有两个含义。
 *
 *
 *
 * - 如果一个求解向量是复值的，那么这就导致在每个评估点至少有两个输入分量。因此， DataPostprocessor::evaluate_scalar_field() 函数永远不会被调用，即使底层有限元只有一个解分量。相反，DataOut将一直*调用 DataPostprocessor::evaluate_vector_field(). 。
 *
 *
 *
 * - 派生类中的 DataPostprocessor::evaluate_vector_field() 的实现必须理解解值是如何安排在他们作为输入收到的 DataPostprocessorInputs::Vector 对象中。  这里的规则是。如果有限元有 $N$ 矢量成分（包括 $N=1$ 的情况，即标量元素），那么复值解向量的输入将有 $2N$ 成分。这些分量首先包含所有求解分量的实部的值（或梯度，或Hessians），然后是所有求解分量的虚部的值（或梯度，或Hessians）。
 * step-58
 * 提供了一个例子，说明这个类（或者说，派生的DataPostprocessorScalar类）是如何在复值情况下使用。
 *
 *
 * @ingroup output
 *
 *
 */
template <int dim>
class DataPostprocessor : public Subscriptor
{
public:
  /**
   * 销毁器。这个函数实际上不做任何事情，但被标记为虚拟，以确保数据后处理器可以通过指向基类的指针被销毁。
   *
   */
  virtual ~DataPostprocessor() override = default;

  /**
   * 这是个主函数，实际执行后处理。第二个参数是对后处理数据的引用，它已经有了正确的大小，必须由这个函数来填充。
   * 该函数通过第一个参数获取所有评估点的解的值、梯度和高阶导数，以及其他数据，如单元格。并非这个参数的所有成员向量都会被填充数据
   *
   * - 事实上，只有当相应的标志被get_need_update_flags()函数（在用户的派生类中实现）的覆盖版本返回时，导数和其他数量才会包含有效的数据。  否则，这些向量将处于一个未指定的状态。    当被DataOut或类似类转换为图形数据的有限元字段代表标量数据时，即如果使用的有限元只有一个实值矢量分量，则调用此函数。
   *
   */
  virtual void
  evaluate_scalar_field(const DataPostprocessorInputs::Scalar<dim> &input_data,
                        std::vector<Vector<double>> &computed_quantities) const;

  /**
   * 与evaluate_scalar_field()函数相同，但当原始数据矢量代表矢量数据，即使用中的有限元有多个矢量分量时，该函数被调用。如果有限元是标量的，但是解向量是复值的，也会调用这个函数。如果要可视化的解向量是复值的（无论标量与否），那么输入数据首先包含解向量在每个评估点的所有实部，然后是所有虚部。
   *
   */
  virtual void
  evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                        std::vector<Vector<double>> &computed_quantities) const;

  /**
   * 返回描述计算量名称的字符串向量。
   *
   */
  virtual std::vector<std::string>
  get_names() const = 0;

  /**
   * 这个函数返回关于如何解释由一个以上数据集组成的输出文件的各个组成部分的信息。
   * 例如，如果一个人有一个2d的斯托克斯方程的有限元，代表分量（u,v,p），我们希望表明前两个，u和v，代表一个逻辑矢量，这样以后当我们生成图形输出时，我们可以把它们交给一个可视化程序，该程序将自动知道把它们作为一个矢量场来渲染，而不是作为两个独立的标量场。
   * 这个函数的默认实现返回一个值的向量
   * DataComponentInterpretation::component_is_scalar,
   * ，表示所有的输出组件都是独立的标量场。然而，如果一个派生类产生的数据表示向量，它可以返回一个包含数值的向量
   * DataComponentInterpretation::component_is_part_of_vector.
   * 在上面的例子中，对于(u,v,p)，人们会返回一个包含组件
   * (DataComponentInterpretation::component_is_part_of_vector,
   * DataComponentInterpretation::component_is_part_of_vector,
   * DataComponentInterpretation::component_is_scalar) 的向量。
   *
   */
  virtual std::vector<DataComponentInterpretation::DataComponentInterpretation>
  get_data_component_interpretation() const;

  /**
   * 返回，必须提供哪些数据来计算派生量。  这必须是 @p
   * update_values,  @p update_gradients,  @p update_hessians 和 @p
   * update_quadrature_points. 的组合。注意，标志 @p
   * update_quadrature_points 更新
   * DataPostprocessorInputs::CommonInputs::evaluation_points.
   * 如果DataPostprocessor要与DataOutFaces结合使用，你也可以通过
   * @p update_normal_vectors 标志要求更新法线。
   * 这些标志的描述可以在 dealii::UpdateFlags. 中找到。
   *
   */
  virtual UpdateFlags
  get_needed_update_flags() const = 0;
};



/**
 * 这个类为DataPostprocessor类所提供的功能提供了一个更简单的接口，如果人们想从传递给DataOut类的有限元场中只计算一个标量。对于这种特殊情况，
 * DataPostprocessor::get_data_component_interpretation()
 * 的返回值是很清楚的，我们将get_names()和get_need_update_flags()返回的值传递给构造函数，这样派生类就不必手工实现这些函数了。
 * 派生类所要做的就是实现一个构造函数，并重载
 * DataPostprocessor::evaluate_scalar_field() 或
 * DataPostprocessor::evaluate_vector_field()
 * ，这在DataPostprocessor类的文档中有所讨论。
 * 如何使用这个类的例子可以在 step-29
 * 中找到，在这种情况下，我们对计算一个复值溶液的大小（标量）感兴趣。在
 * step-29 中，解向量由独立的实部和虚部组成，而 step-58
 * 将解向量计算为具有复数项的向量，DataPostprocessorScalar类在那里被用来计算解的幅值和相位，那里的方式不同。
 * 关于如何使用密切相关的DataPostprocessorVector类的例子，可以在该类的文档中找到。对于DataPostprocessorTensor类也是如此。
 *
 *
 * @ingroup output
 *
 *
 */
template <int dim>
class DataPostprocessorScalar : public DataPostprocessor<dim>
{
public:
  /**
   * 构造函数。获取由当前类派生的单一标量变量的名称，以及计算该数量所需的更新标志。
   * @param  name
   * 这个类计算的标量变量的名称应该在图形输出文件中提供。
   * @param  update_flags 这必须是 @p update_values,   @p
   * update_gradients,   @p update_hessians  和  @p update_quadrature_points.
   * 的组合。注意，标志 @p update_quadrature_points 更新
   * DataPostprocessorInputs::CommonInputs::evaluation_points.
   * 如果DataPostprocessor要与DataOutFaces结合使用，你也可以通过
   * @p update_normal_vectors 标志要求更新法线。
   * 这些标志的描述可以在 dealii::UpdateFlags. 中找到。
   *
   */
  DataPostprocessorScalar(const std::string &name,
                          const UpdateFlags  update_flags);

  /**
   * 返回描述计算量名称的字符串向量。鉴于这个类的目的，这是一个单项的向量，等于给构造函数的名称。
   *
   */
  virtual std::vector<std::string>
  get_names() const override;

  /**
   * 这个函数返回关于如何解释由一个以上的数据集组成的输出文件的各个组成部分的信息。由于当前的类是为了用于单个标量结果变量，所以返回值显然是
   * DataComponentInterpretation::component_is_scalar. 。
   *
   */
  virtual std::vector<DataComponentInterpretation::DataComponentInterpretation>
  get_data_component_interpretation() const override;

  /**
   * 返回，必须提供哪些数据来计算派生量。
   * 这里返回的标志是传递给这个类的构造函数的标志。
   *
   */
  virtual UpdateFlags
  get_needed_update_flags() const override;

private:
  /**
   * 复制给这个类的构造函数的两个参数。
   *
   */
  const std::string name;
  const UpdateFlags update_flags;
};



/**
 * 该类为DataPostprocessor类所提供的功能提供了一个更简单的接口，以防人们只想从传递给DataOut类的有限元场中计算一个单一的矢量（定义为正好有
 * @p dim 分量）。对于这种特殊情况，
 * DataPostprocessor::get_data_component_interpretation()
 * 的返回值是很清楚的，我们将get_names()和get_need_update_flags()返回的值传递给构造函数，这样派生类就不必手工实现这些函数。
 * 派生类所要做的就是实现一个构造函数，并重载
 * DataPostprocessor::evaluate_scalar_field() 或
 * DataPostprocessor::evaluate_vector_field()
 * ，这在DataPostprocessor类的文档中有所讨论。
 * 与之密切相关的类DataPostprocessorScalar的使用方法的例子可以在
 * step-29
 * 中找到。关于如何使用DataPostprocessorTensor类的例子可以在该类的文档中找到。
 *
 *  <h3> An example </h3>
 * 人们想用后处理器做的一个常见的例子是，不仅要把解的值可视化，还要把梯度可视化。事实上，这正是
 * step-19
 * 所需要的，因此它几乎逐字逐句地使用下面的代码。为了简单起见，让我们假设你只有一个标量的解。事实上，因为它是现成的，让我们简单地使用
 * step-6
 * 求解器来产生这样一个标量解。梯度是一个矢量（正好有
 * @p dim
 * 分量），所以当前的类适合通过后处理来产生梯度。然后，下面的代码片段实现了可视化梯度所需的一切。
 *
 * @code
 * template <int dim>
 * class GradientPostprocessor : public DataPostprocessorVector<dim>
 * {
 * public:
 * GradientPostprocessor ()
 *   :
 *   // call the constructor of the base class. call the variable to
 *   // be output "grad_u" and make sure that DataOut provides us
 *   // with the gradients:
 *   DataPostprocessorVector<dim> ("grad_u",
 *                                 update_gradients)
 * {}
 *
 * virtual
 * void
 * evaluate_scalar_field
 * (const DataPostprocessorInputs::Scalar<dim> &input_data,
 *  std::vector<Vector<double> > &computed_quantities) const override
 * {
 *   // ensure that there really are as many output slots
 *   // as there are points at which DataOut provides the
 *   // gradients:
 *   AssertDimension (input_data.solution_gradients.size(),
 *                    computed_quantities.size());
 *
 *   // then loop over all of these inputs:
 *   for (unsigned int p=0; p<input_data.solution_gradients.size(); ++p)
 *     {
 *       // ensure that each output slot has exactly 'dim'
 *       // components (as should be expected, given that we
 *       // want to create vector-valued outputs), and copy the
 *       // gradients of the solution at the evaluation points
 *       // into the output slots:
 *       AssertDimension (computed_quantities[p].size(), dim);
 *       for (unsigned int d=0; d<dim; ++d)
 *         computed_quantities[p][d]
 *           = input_data.solution_gradients[p][d];
 *     }
 * }
 * };
 * @endcode
 * 唯一需要的是在该示例程序的 @p Step6 类的 @p run()
 * 函数中，为 DataOut::add_vector()
 * 的调用添加另一个输出。然后，相应的代码片断会是这样的（在这里我们也使用VTU作为文件格式来输出数据）。
 *
 * @code
 * GradientPostprocessor<dim> gradient_postprocessor;
 *
 * DataOut<dim> data_out;
 * data_out.attach_dof_handler (dof_handler);
 * data_out.add_data_vector (solution, "solution");
 * data_out.add_data_vector (solution, gradient_postprocessor);
 * data_out.build_patches ();
 *
 * std::ofstream output ("solution.vtu");
 * data_out.write_vtu (output);
 * @endcode
 *
 * 这将导致以下解决方案和梯度的输出（你可能想与 step-6
 * 的结果部分所示的解决方案进行比较；为了简单起见，目前的数据是在一个更粗的网格上产生的）。
*  @image html data_postprocessor_vector_0.png   @image html data_postprocessor_vector_1.png 。
 * 在第二张图片中，背景颜色对应于梯度矢量的大小，矢量字形对应于梯度本身。一开始，看到从每个顶点出发的多个向量，向不同的方向发展，可能会让人感到惊讶。但这是因为解决方案只是连续的：一般来说，梯度在边缘上是不连续的，所以从每个顶点出发的多个向量只是代表解决方案在每个相邻单元的不同梯度。
 * 上面的输出
 *
 * - 即解决方案的梯度 $\nabla u$ 。
 *
 * 如果把 step-6
 * 解释为解决一个稳态传热问题，则对应于温度梯度。它在域的中心部分非常小，因为在
 * step-6 中，我们正在解决一个方程，其系数 $a(\mathbf x)$
 * 在中心部分很大，在外部很小。这可以被认为是一种导热性好的材料，因此温度梯度很小。另一方面，"热通量
 * "与数量 $a(\mathbf x) \nabla u(\mathbf x)$
 * 相对应。对于该方程的解，通量在界面上应该是连续的。这很容易通过对后处理程序的以下修改得到验证。
 *
 * @code
 * template <int dim>
 * class HeatFluxPostprocessor : public DataPostprocessorVector<dim>
 * {
 * public:
 * HeatFluxPostprocessor ()
 *   :
 *   // like above, but now also make sure that DataOut provides
 *   // us with coordinates of the evaluation points:
 *   DataPostprocessorVector<dim> ("heat_flux",
 *                                 update_gradients | update_quadrature_points)
 * {}
 *
 * virtual
 * void
 * evaluate_scalar_field
 * (const DataPostprocessorInputs::Scalar<dim> &input_data,
 *  std::vector<Vector<double> > &computed_quantities) const override
 * {
 *   AssertDimension (input_data.solution_gradients.size(),
 *                    computed_quantities.size());
 *
 *   for (unsigned int p=0; p<input_data.solution_gradients.size(); ++p)
 *     {
 *       AssertDimension (computed_quantities[p].size(), dim);
 *       for (unsigned int d=0; d<dim; ++d)
 *         // like above, but also multiply the gradients with
 *         // the coefficient evaluated at the current point:
 *         computed_quantities[p][d]
 *           = coefficient (input_data.evaluation_points[p])
 *             input_data.solution_gradients[p][d];
 *     }
 * }
 * };
 * @endcode
 * 通过这个后处理器，我们可以得到以下的热通量图。
*  @image html data_postprocessor_vector_2.png
 * 如背景颜色所示，梯度乘以系数现在是一个连续函数。在界面周围有一些（大的）矢量，在那里系数跳跃（在圆盘中心到周边的一半距离处），似乎指向错误的方向；这是一个伪命题，因为在这些点上，解决方案有一个不连续的梯度，而且当前网格上的数值解决方案不能充分解决这个界面。然而，这对目前的讨论并不重要。
 *
 *  <h3> Extension to the gradients of vector-valued problems </h3>
 * 上面的例子使用了一个标量解和它的梯度作为例子。另一方面，人们可能想对一个矢量值位移场的梯度做类似的事情（比如位移场的应变或应力，像那些在
 * step-8 ,  step-17 ,  step-18 , 或 step-44
 * 中计算的）。在这种情况下，解决方案已经是矢量值的，应力是一个（对称的）张量。
 * deal.II目前不支持输出张量值的数量，但它们当然可以作为张量的标量值成分的集合来输出。这可以通过使用DataPostprocessorTensor类来实现。该类的文档包含一个例子。
 *
 *
 *
 * @ingroup output
 *
 */
template <int dim>
class DataPostprocessorVector : public DataPostprocessor<dim>
{
public:
  /**
   * 构造函数。获取由当前类派生的单一向量变量的名称，以及计算该数量所需的更新标志。
   * @param  name
   * 这个类计算的向量变量的名称应该在图形输出文件中提供。
   * @param  update_flags 这必须是 @p update_values,   @p
   * update_gradients,   @p update_hessians  和  @p update_quadrature_points.
   * 的组合。注意，标志 @p update_quadrature_points 更新
   * DataPostprocessorInputs::CommonInputs::evaluation_points.
   * 如果DataPostprocessor要与DataOutFaces结合使用，你也可以通过
   * @p update_normal_vectors 标志要求更新法向。
   * 这些标志的描述可以在 dealii::UpdateFlags. 中找到。
   *
   */
  DataPostprocessorVector(const std::string &name,
                          const UpdateFlags  update_flags);

  /**
   * 返回描述计算量名称的字符串向量。考虑到这个类的目的，这是一个有dim项的向量，都等于给构造函数的名称。
   *
   */
  virtual std::vector<std::string>
  get_names() const override;

  /**
   * 这个函数返回关于如何解释由一个以上的数据集组成的输出文件的各个组成部分的信息。由于当前的类是为了用于单个矢量结果变量，所以返回值显然是
   * DataComponentInterpretation::component_is_part 重复dim次。
   *
   */
  virtual std::vector<DataComponentInterpretation::DataComponentInterpretation>
  get_data_component_interpretation() const override;

  /**
   * 返回必须提供哪些数据来计算派生量。
   * 这里返回的标志是传递给这个类的构造函数的标志。
   *
   */
  virtual UpdateFlags
  get_needed_update_flags() const override;

private:
  /**
   * 复制给这个类的构造函数的两个参数。
   *
   */
  const std::string name;
  const UpdateFlags update_flags;
};



/**
 * 这个类为DataPostprocessor类所提供的功能提供了一个更简单的接口，如果人们想从传递给DataOut类的有限元场中只计算一个张量（定义为正好有
 * <code>dim*dim</code> 个分量）。
 * 对于这种情况，我们希望将所有这些分量作为张量值的一部分输出。不幸的是，以图形文件格式写入DataOut数据的各种后端（见DataOutBase命名空间，了解可以写入哪些格式）在当前不支持张量数据。事实上，提供语义信息的DataComponentInterpretation命名空间也不支持如何解释图形数据的单个组件。尽管如此，像DataPostprocessorScalar和DataPostprocessorVector一样，这个类有助于设置DataPostprocessor基类所要求的get_names()和get_need_update_flags()函数应该返回什么，因此当前类根据当前类的构造函数从进一步的派生类中收到的信息实现这些。
 * （为了将这个标量字段的集合可视化，然后将其解释为张量，我们必须（i）使用一个能够可视化张量的可视化程序，并且（ii）教它如何将标量字段重新组合成张量。在VisIt的情况下
 *
 * - 见https://wci.llnl.gov/simulation/computer-codes/visit/
 *
 * - 这是通过创建一个新的 "表达式 "来完成的：实质上，我们创建了一个变量，比如 "grad_u"，它是张量值，其值由表达式<code>{{grad_u_xx,grad_u_xy},{grad_u_yx, grad_u_yy}}</code>给出，其中引用的变量是标量场的名称，这里是由以下例子产生。然后VisIt能够将这个 "新 "变量可视化为一个张量）。)
 * 所有派生类要做的就是实现一个构造函数，并重载
 * DataPostprocessor::evaluate_scalar_field() 或
 * DataPostprocessor::evaluate_vector_field()
 * ，在DataPostprocessor类的文档中讨论过。
 * 与之密切相关的类DataPostprocessorScalar的使用方法的例子可以在
 * step-29
 * 中找到。关于如何使用DataPostprocessorVector类的例子可以在该类的文档中找到。
 *
 *  <h3> An example </h3>
 * 人们想用后处理器做的一个常见的例子是，不仅要将解的值可视化，还要将梯度可视化。这个类是为了张量值的输出，所以我们将从一个矢量值的解决方案开始：
 * step-8  的位移场。梯度是一个等级2的张量（正好有
 * <code>dim*dim</code>
 * 的分量），所以当前的类适合通过后处理产生梯度。然后，下面的代码片段实现了可视化梯度所需的一切。
 *
 * @code
 * template <int dim>
 * class GradientPostprocessor : public DataPostprocessorTensor<dim>
 * {
 * public:
 *   GradientPostprocessor ()
 *     :
 *     DataPostprocessorTensor<dim> ("grad_u",
 *                                   update_gradients)
 *   {}
 *
 *   virtual
 *   void
 *   evaluate_vector_field
 *   (const DataPostprocessorInputs::Vector<dim> &input_data,
 *    std::vector<Vector<double> > &computed_quantities) const override
 *   {
 *     // ensure that there really are as many output slots
 *     // as there are points at which DataOut provides the
 *     // gradients:
 *     AssertDimension (input_data.solution_gradients.size(),
 *                      computed_quantities.size());
 *
 *     for (unsigned int p=0; p<input_data.solution_gradients.size(); ++p)
 *       {
 *         // ensure that each output slot has exactly 'dim*dim'
 *         // components (as should be expected, given that we
 *         // want to create tensor-valued outputs), and copy the
 *         // gradients of the solution at the evaluation points
 *         // into the output slots:
 *         AssertDimension (computed_quantities[p].size(),
 *                          (Tensor<2,dim>::n_independent_components));
 *         for (unsigned int d=0; d<dim; ++d)
 *           for (unsigned int e=0; e<dim; ++e)
 *             computed_quantities[p][Tensor<2,dim>::component_to_unrolled_index(TableIndices<2>(d,e))]
 *               = input_data.solution_gradients[p][d][e];
 *       }
 *   }
 * };
 * @endcode
 * 这段代码中唯一棘手的部分是如何将应变张量的
 * <code>dim*dim</code> 元素排序到计算输出量的一个向量中去。
 *
 * - 换句话说，如何将张量的元素<i>unroll</i>放入矢量中。这由 Tensor::component_to_unrolled_index() 函数提供便利，该函数接收一对指定张量特定元素的索引，并返回一个向量索引，然后在上面的代码中用于填充 @p computed_quantities 数组。
 * 最后一件事是在该示例程序的 @p Step8 类的 @p output_results()
 * 函数中为 DataOut::add_vector()
 * 的调用添加另一个输出。然后，相应的代码片断会是这样的。
 *
 * @code
 *   GradientPostprocessor<dim> grad_u;
 *
 *   DataOut<dim> data_out;
 *   data_out.attach_dof_handler (dof_handler);
 *
 *   std::vector<DataComponentInterpretation::DataComponentInterpretation>
 *   data_component_interpretation
 *   (dim, DataComponentInterpretation::component_is_part_of_vector);
 *   data_out.add_data_vector (solution,
 *                             std::vector<std::string>(dim,"displacement"),
 *                             DataOut<dim>::type_dof_data,
 *                             data_component_interpretation);
 *   data_out.add_data_vector (solution, grad_u);
 *   data_out.build_patches ();
 *   data_out.write_vtk (output);
 * @endcode
 *
 * 这将导致以下位移场（即解决方案）和梯度的输出（你可能想与
 * step-8
 * 的结果部分显示的解决方案进行比较；为简单起见，当前数据是在均匀网格上生成的）。
*  @image html data_postprocessor_tensor_0.png   @image html data_postprocessor_tensor_1.png 。
 * 这些图片显示了平均每十个网格点上代表梯度张量的椭圆。你可能想通过阅读VisIt可视化程序的文档（见https://wci.llnl.gov/simulation/computer-codes/visit/）来了解张量究竟是如何被可视化的。
 * 在弹性中，人们往往对位移的梯度不感兴趣，而是对
 * "应变 "感兴趣，即梯度的对称版本  $\varepsilon=\frac 12
 * (\nabla u + \nabla u^T)$
 * 。这很容易通过以下的小修改来实现。
 *
 * @code
 * template <int dim>
 * class StrainPostprocessor : public DataPostprocessorTensor<dim>
 * {
 * public:
 *   StrainPostprocessor ()
 *     :
 *     DataPostprocessorTensor<dim> ("strain",
 *                                   update_gradients)
 *   {}
 *
 *   virtual
 *   void
 *   evaluate_vector_field
 *   (const DataPostprocessorInputs::Vector<dim> &input_data,
 *    std::vector<Vector<double> > &computed_quantities) const override
 *   {
 *     AssertDimension (input_data.solution_gradients.size(),
 *                      computed_quantities.size());
 *
 *     for (unsigned int p=0; p<input_data.solution_gradients.size(); ++p)
 *       {
 *         AssertDimension (computed_quantities[p].size(),
 *                          (Tensor<2,dim>::n_independent_components));
 *         for (unsigned int d=0; d<dim; ++d)
 *           for (unsigned int e=0; e<dim; ++e)
 *             computed_quantities[p][Tensor<2,dim>::component_to_unrolled_index(TableIndices<2>(d,e))]
 *               = (input_data.solution_gradients[p][d][e]
 *                  +
 *                  input_data.solution_gradients[p][e][d]) / 2;
 *       }
 *   }
 * };
 * @endcode
 *
 * 在 step-8 中使用这个类导致了以下的可视化。
*  @image html data_postprocessor_tensor_2.png
 * 鉴于输出应变很容易，写一个后处理程序来计算解场中的<i>stress</i>也不会很复杂，因为通过与应变-应力张量或在简单情况下与Lam&eacute;常数相乘，应力很容易从应变中计算出来。
 *
 *
 * @ingroup output
 *
 *
 */
template <int dim>
class DataPostprocessorTensor : public DataPostprocessor<dim>
{
public:
  /**
   * 构造函数。获取由当前类派生的单一矢量变量的名称，以及计算该数量所需的更新标志。
   * @param  name
   * 这个类计算的向量变量的名称应该在图形输出文件中提供。
   * @param  update_flags 这必须是 @p update_values,   @p
   * update_gradients,   @p update_hessians  和  @p update_quadrature_points.
   * 的组合。注意，标志 @p update_quadrature_points 更新
   * DataPostprocessorInputs::CommonInputs::evaluation_points.
   * 如果DataPostprocessor要与DataOutFaces结合使用，你也可以通过
   * @p update_normal_vectors 标志要求更新法向。
   * 这些标志的描述可以在 dealii::UpdateFlags. 中找到。
   *
   */
  DataPostprocessorTensor(const std::string &name,
                          const UpdateFlags  update_flags);

  /**
   * 返回描述计算量名称的字符串向量。考虑到这个类的目的，这是一个有dim项的向量，都等于给构造函数的名称。
   *
   */
  virtual std::vector<std::string>
  get_names() const override;

  /**
   * 这个函数返回关于如何解释由一个以上的数据集组成的输出文件的各个组成部分的信息。由于当前的类是为了用于单个矢量结果变量，所以返回值显然是
   * DataComponentInterpretation::component_is_part 重复dim次。
   *
   */
  virtual std::vector<DataComponentInterpretation::DataComponentInterpretation>
  get_data_component_interpretation() const override;

  /**
   * 返回必须提供哪些数据来计算派生数量。
   * 这里返回的标志是传递给这个类的构造函数的标志。
   *
   */
  virtual UpdateFlags
  get_needed_update_flags() const override;

private:
  /**
   * 复制给这个类的构造函数的两个参数。
   *
   */
  const std::string name;
  const UpdateFlags update_flags;
};



#ifndef DOXYGEN
// -------------------- template functions ----------------------

namespace DataPostprocessorInputs
{
  template <int spacedim>
  template <typename DoFHandlerType>
  void
  CommonInputs<spacedim>::set_cell(
    const typename DoFHandlerType::cell_iterator &new_cell)
  {
    return set_cell<DoFHandlerType::dimension>(new_cell);
  }



  template <int spacedim>
  template <int dim>
  void
  CommonInputs<spacedim>::set_cell(
    const typename DoFHandler<dim, spacedim>::cell_iterator &new_cell)
  {
    // see if we had previously already stored a cell that has the same
    // data type; if so, reuse the memory location and avoid calling 'new'
    // inside boost::any
    if (typename DoFHandler<dim, spacedim>::cell_iterator *storage_location =
          boost::any_cast<typename DoFHandler<dim, spacedim>::cell_iterator>(
            &cell))
      *storage_location = new_cell;
    else
      // if we had nothing stored before, or if we had stored a different
      // data type, just let boost::any replace things
      cell = new_cell;
  }



  template <int spacedim>
  template <typename DoFHandlerType>
  typename DoFHandlerType::cell_iterator
  CommonInputs<spacedim>::get_cell() const
  {
    return get_cell<DoFHandlerType::dimension>();
  }



  template <int spacedim>
  template <int dim>
  typename DoFHandler<dim, spacedim>::cell_iterator
  CommonInputs<spacedim>::get_cell() const
  {
    Assert(cell.empty() == false,
           ExcMessage(
             "You are trying to access the cell associated with a "
             "DataPostprocessorInputs::Scalar object for which no cell has "
             "been set."));
    Assert((boost::any_cast<typename DoFHandler<dim, spacedim>::cell_iterator>(
              &cell) != nullptr),
           ExcMessage(
             "You are trying to access the cell associated with a "
             "DataPostprocessorInputs::Scalar with a DoFHandler type that is "
             "different from the type with which it has been set. For example, "
             "if the cell for which output is currently being generated "
             "belongs to a DoFHandler<2, 3> object, then you can only call the "
             "current function with a template argument equal to "
             "DoFHandler<2, 3>, but not with any other class type or dimension "
             "template argument."));

    return boost::any_cast<typename DoFHandler<dim, spacedim>::cell_iterator>(
      cell);
  }
} // namespace DataPostprocessorInputs

#endif

DEAL_II_NAMESPACE_CLOSE

#endif


