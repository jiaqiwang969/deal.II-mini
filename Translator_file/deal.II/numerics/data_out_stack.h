//include/deal.II-translator/numerics/data_out_stack_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2020 by the deal.II authors
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

#ifndef dealii_data_out_stack_h
#define dealii_data_out_stack_h


#include <deal.II/base/config.h>

#include <deal.II/base/data_out_base.h>
#include <deal.II/base/smartpointer.h>

#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out_dof_data.h>

#include <string>
#include <vector>

DEAL_II_NAMESPACE_OPEN

// Forward declaration
#ifndef DOXYGEN
template <int dim, int spacedim>
class DoFHandler;
#endif

/**
 * @deprecated  使用DataOutStack<dim, spacedim>代替。
 *
 *
 */
template <int dim, int spacedim = dim, typename DoFHandlerType = void>
class DataOutStack;

#ifndef DOXYGEN
// prevent doxygen from complaining about potential recursive class relations
template <int dim, int spacedim, typename DoFHandlerType>
class DataOutStack : public DataOutStack<dim, spacedim, void>
{
public:
  DEAL_II_DEPRECATED
  DataOutStack()
    : DataOutStack<dim, spacedim, void>()
  {}
};
#endif // DOXYGEN

/**
 * 这个类用于将几个计算的输出堆叠到一个输出文件中，方法是将数据集在另一个与空间方向正交的坐标方向上堆叠。最常见的用途是将几个时间步骤的结果堆叠到一个时空输出文件中，或者例如将几个参数值的参数相关方程的解的结果连接到一起。该接口主要是以DataOut类为模型，更多的文件请看那里。
 * 我们将为一个时间相关的问题解释这个概念，但是可以用任何参数来代替时间。在我们的例子中，一个方程的解被计算到每个离散的时间水平。然后，这将被添加到本类的一个对象中，在所有的时间层被添加后，一个空间-时间图将以基类支持的任何输出格式被写入。在输出时，每个时间层上的（空间）解被扩展到时间方向上，写两次，一次是时间层本身，一次是等于时间层减去给定时间步长的时间。这两个副本被连接起来，形成一个时空板块，在时间上有恒定的值。
 * 由于时间上的片状常数输出，一般来说，写出的解在离散的时间水平上是不连续的，但在大多数情况下，输出仍然是足够的。未来可能会增加更复杂的时间内插值。
 *
 *  <h3>Example of Use</h3>
 * 下面的小例子将说明该类的不同使用步骤。假设使用的有限元由两个部分组成，
 * @p u 和 @p v, ，解向量被命名为 @p solution
 * ，并且计算了一个向量 @p error
 * ，其中包含每个空间单元的误差指标。
 * 请注意，与DataOut类不同的是，在第一次使用之前，有必要首先声明数据向量和组件的名称。这是因为在所有的时间层次上，都应该有相同的数据来产生合理的时间空间输出。生成的输出在每个空间和时间方向上都有两个分项，这适用于空间的二次有限元，例如。
 *
 *
 * @code
 * DataOutStack<dim> data_out_stack;
 *
 *                                // first declare the vectors
 *                                // to be used later
 * std::vector<std::string> solution_names;
 * solution_names.emplace_back ("u");
 * solution_names.emplace_back ("v");
 * data_out_stack.declare_data_vector (solution_names,
 *                                     DataOutStack<dim>::dof_vector);
 * data_out_stack.declare_data_vector ("error",
 *                                     DataOutStack<dim>::cell_vector);
 *
 *                                // now do computations
 * for (double parameter=0; ...)
 *   {
 *     DoFHandler<dim,spacedim> dof_handler;
 *     ...                        // compute something
 *
 *                                // now for output
 *     data_out_stack.new_parameter_value (parameter,
 *                                         delta_parameter);
 *     data_out_stack.attach_dof_handler (dof_handler);
 *     data_out_stack.add_data_vector (solution, solution_names);
 *     data_out_stack.add_data_vector (error, "error");
 *     data_out_stack.build_patches (2);
 *     data_out_stack.finish_parameter_value ();
 *   };
 * @endcode
 *
 *
 * @ingroup output
 *
 */
template <int dim, int spacedim>
class DataOutStack<dim, spacedim, void>
  : public DataOutInterface<dim + 1, spacedim + 1>
{
  static_assert(dim == spacedim, "Not implemented for dim != spacedim.");

public:
  /**
   * 补丁的尺寸参数。
   *
   */
  static constexpr int patch_dim      = dim + 1;
  static constexpr int patch_spacedim = spacedim + 1;

  /**
   * 声明该类中使用的两种类型的向量的数据类型。
   *
   */
  enum VectorType
  {
    /**
     * 数据描述每个单元的一个值。
     *
     */
    cell_vector,
    /**
     * 数据为每个DoF描述一个值。
     *
     */
    dof_vector
  };

  /**
   * 解构器。只声明使其成为 @p virtual. 。
   *
   */
  virtual ~DataOutStack() override = default;

  /**
   * 为一个特定的参数值开始下一组数据。参数 @p
   * parameter_step 表示间隔（向后方向，从 @p parameter_value)
   * 算起），输出将以参数方向扩展，即与空间方向正交。
   *
   */
  void
  new_parameter_value(const double parameter_value,
                      const double parameter_step);

  /**
   * 附加网格的DoF处理程序，以及与之前 @p new_parameter_value.
   * 设置的参数相关的数据
   * 这必须在为当前参数值添加数据向量之前发生。
   *
   */
  void
  attach_dof_handler(const DoFHandler<dim, spacedim> &dof_handler);

  /**
   * 申报一个数据向量。 @p vector_type
   * 参数决定了数据向量是否将被视为DoF或单元数据。
   * 如果DoFHandler目前使用的有限元（之前附加到这个对象）只有一个分量，因此只需要给出一个名称，则可以调用这个版本。
   *
   */
  void
  declare_data_vector(const std::string &name, const VectorType vector_type);

  /**
   * 声明一个数据矢量。 @p vector_type
   * 参数决定了该数据向量将被视为DoF或单元格数据。
   * 如果DoFHandler目前使用的有限元（之前附加到这个对象）有一个以上的分量，因此需要给出一个以上的名称，则必须调用这个版本。然而，如果有限元只有一个分量，你也可以用只包含一个元素的
   * <tt>std::vector@<std::string@></tt> 调用这个函数。
   *
   */
  void
  declare_data_vector(const std::vector<std::string> &name,
                      const VectorType                vector_type);


  /**
   * 为目前设定的参数值添加一个数据向量。
   * 如果DoFHandler目前使用的有限元（之前附加到这个对象）只有一个分量，因此只需要给出一个名称，则可以调用这个版本。
   * 如果 @p vec 是一个有多个分量的向量，该函数将通过在
   * @p name
   * 中添加下划线和每个分量的编号来为所有分量生成不同的名称。
   * 数据向量在第一次实际使用之前必须使用 @p
   * declare_data_vector函数注册。
   * 请注意，这个向量的副本会一直保存到下一次调用 @p
   * finish_parameter_value
   * 为止，所以如果你的内存不足，你可能想在所有涉及大矩阵的计算都已经完成之后再调用这个函数。
   *
   */
  template <typename number>
  void
  add_data_vector(const Vector<number> &vec, const std::string &name);

  /**
   * 为目前设定的参数值添加一个数据向量。
   * 如果DoFHandler目前使用的有限元（之前附加到这个对象上）有一个以上的分量，因此需要给出一个以上的名称，则必须调用这个版本。然而，如果有限元只有一个分量，你也可以用只包含一个元素的
   * <tt>std::vector@<std::string@></tt> 调用这个函数。
   * 在实际第一次使用之前，数据向量必须已经用 @p
   * declare_data_vector函数注册。
   * 请注意，这个向量的副本会一直保存到下一次调用 @p
   * finish_parameter_value
   * 为止，所以如果你的内存不足，你可能希望在所有涉及大矩阵的计算都已经完成之后再调用这个函数。
   *
   */
  template <typename number>
  void
  add_data_vector(const Vector<number> &          vec,
                  const std::vector<std::string> &names);

  /**
   * 这是这个类的中心函数，因为它建立了由基类的低级函数编写的补丁列表。从本质上讲，补丁是三角形和DoFHandler对象的每个单元上的数据的一些中间表示，然后可以用来以某种可视化程序可读的格式写入文件。
   * 你可以在这个类的一般文档中找到关于这个函数的使用概述。在这个类的基类DataOut_DoFData的文档中也提供了一个例子。
   * @param  n_subdivisions 参见 DataOut::build_patches()
   * 对这个参数的广泛描述。分数的数量在这个类所使用的类似时间的参数的方向上总是一个。
   *
   */
  void
  build_patches(const unsigned int n_subdivisions = 0);

  /**
   * 一旦 @p build_patches
   * 被调用，释放所有不再需要的数据，并且所有其他针对给定参数值的事务都完成。
   * 与 @p new_parameter_value. 相对应。
   *
   */
  void
  finish_parameter_value();

  /**
   * 确定这个对象的内存消耗（以字节为单位）的估计值。
   *
   */
  std::size_t
  memory_consumption() const;

  /**
   * 例外情况
   *
   */
  DeclException1(
    ExcVectorNotDeclared,
    std::string,
    << "The data vector for which the first component has the name " << arg1
    << " has not been added before.");
  /**
   * 异常情况
   *
   */
  DeclExceptionMsg(ExcDataNotCleared,
                   "You cannot start a new time/parameter step before calling "
                   "finish_parameter_value() on the previous step.");
  /**
   * 异常情况
   *
   */
  DeclExceptionMsg(
    ExcDataAlreadyAdded,
    "You cannot declare additional vectors after already calling "
    "build_patches(). All data vectors need to be declared "
    "before you call this function the first time.");
  /**
   * 异常情况
   *
   */
  DeclException1(ExcNameAlreadyUsed,
                 std::string,
                 << "You tried to declare a component of a data vector with "
                 << "the name <" << arg1
                 << ">, but that name is already used.");

private:
  /**
   * 目前的参数值。
   *
   */
  double parameter;

  /**
   * 目前的参数步骤，即接下来要写入的参数区间的长度。
   *
   */
  double parameter_step;

  /**
   * 对应于当前参数值的数据所使用的DoF处理程序。
   *
   */
  SmartPointer<const DoFHandler<dim, spacedim>,
               DataOutStack<dim, spacedim, void>>
    dof_handler;

  /**
   * 所有过去和现在的参数值数据集的补丁列表。
   *
   */
  std::vector<dealii::DataOutBase::Patch<patch_dim, patch_spacedim>> patches;

  /**
   * 保存当前参数值的数据向量（单元和dof数据）的结构。
   *
   */
  struct DataVector
  {
    /**
     * 数据向量。
     *
     */
    Vector<double> data;

    /**
     * 每个这样的数据集内的不同组件的名称。
     *
     */
    std::vector<std::string> names;

    /**
     * 确定这个对象的内存消耗（以字节为单位）的估计值。
     *
     */
    std::size_t
    memory_consumption() const;
  };

  /**
   * 列表中的DoF数据向量。
   *
   */
  std::vector<DataVector> dof_data;

  /**
   * 单元数据向量的列表。
   *
   */
  std::vector<DataVector> cell_data;

  /**
   * 这是一个函数，派生类通过这个函数将Patch结构形式的预处理数据（在基类DataOutBase中声明）传播给实际的输出函数。
   *
   */
  virtual const std::vector<dealii::DataOutBase::Patch<
    DataOutStack<dim, spacedim, void>::patch_dim,
    DataOutStack<dim, spacedim, void>::patch_spacedim>> &
  get_patches() const override;


  /**
   * 虚拟函数，基类的输出函数通过它获得数据集的名称。
   *
   */
  virtual std::vector<std::string>
  get_dataset_names() const override;
};


DEAL_II_NAMESPACE_CLOSE

#endif


