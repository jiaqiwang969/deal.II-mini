//include/deal.II-translator/numerics/data_out_dof_data_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2021 by the deal.II authors
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

#ifndef dealii_data_out_dof_data_h
#define dealii_data_out_dof_data_h



#include <deal.II/base/config.h>

#include <deal.II/base/data_out_base.h>
#include <deal.II/base/mg_level_object.h>
#include <deal.II/base/smartpointer.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/mapping.h>

#include <deal.II/grid/tria.h>

#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/fe_values.h>
#include <deal.II/hp/mapping_collection.h>
#include <deal.II/hp/q_collection.h>

#include <deal.II/numerics/data_component_interpretation.h>
#include <deal.II/numerics/data_postprocessor.h>

#include <memory>

DEAL_II_NAMESPACE_OPEN

namespace Exceptions
{
  /**
   * 异常的命名空间，这些异常在整个DataOut*类的集合中使用。
   * 类的集合中使用的异常命名空间。
   *
   */
  namespace DataOutImplementation
  {
    /**
     * 异常情况
     *
     */
    DeclException1(ExcInvalidNumberOfSubdivisions,
                   int,
                   << "The number of subdivisions per patch, " << arg1
                   << ", is not valid. It needs to be greater or equal to "
                      "one, or zero if you want it to be determined "
                      "automatically.");

    /**
     * 异常
     *
     */
    DeclExceptionMsg(ExcNoTriangulationSelected,
                     "For the operation you are attempting, you first need to "
                     "tell the DataOut or related object which DoFHandler or "
                     "triangulation you would like to work on.");

    /**
     * 异常情况
     *
     */
    DeclExceptionMsg(ExcNoDoFHandlerSelected,
                     "For the operation you are attempting, you first need to "
                     "tell the DataOut or related object which DoFHandler "
                     "you would like to work on.");

    /**
     * 异常情况
     *
     */
    DeclException3(ExcInvalidVectorSize,
                   int,
                   int,
                   int,
                   << "The vector has size " << arg1
                   << " but the DoFHandler object says that there are " << arg2
                   << " degrees of freedom and there are " << arg3
                   << " active cells. The size of your vector needs to be"
                   << " either equal to the number of degrees of freedom (when"
                   << " the data is of type type_dof_data), or equal to the"
                   << " number of active cells (when the data is of type "
                   << " type_cell_data).");
    /**
     * 异常情况
     *
     */
    DeclException2(
      ExcInvalidCharacter,
      std::string,
      size_t,
      << "Please use only the characters [a-zA-Z0-9_<>()] for" << std::endl
      << "description strings since some graphics formats will only accept these."
      << std::endl
      << "The string you gave was <" << arg1
      << ">, within which the invalid character is <" << arg1[arg2] << ">."
      << std::endl);
    /**
     * 异常情况
     *
     */
    DeclExceptionMsg(
      ExcOldDataStillPresent,
      "When attaching a triangulation or DoFHandler object, it is "
      "not allowed if old data vectors are still referenced. If "
      "you want to reuse an object of the current type, you first "
      "need to call the 'clear_data_vector()' function.");
    /**
     * 异常情况
     *
     */
    DeclException2(ExcInvalidNumberOfNames,
                   int,
                   int,
                   << "You have to give one name per component in your "
                   << "data vector. The number you gave was " << arg1
                   << ", but the number of components is " << arg2 << ".");
    /**
     * 异常情况
     *
     */
    DeclExceptionMsg(ExcIncompatibleDatasetNames,
                     "While merging sets of patches, the two sets to be merged "
                     "need to refer to data that agrees on the names of the "
                     "various variables represented. In other words, you "
                     "cannot merge sets of patches that originate from "
                     "entirely unrelated simulations.");
    /**
     * 异常情况
     *
     */
    DeclExceptionMsg(ExcIncompatiblePatchLists,
                     "While merging sets of patches, the two sets to be merged "
                     "need to refer to data that agrees on the number of "
                     "subdivisions and other properties. In other words, you "
                     "cannot merge sets of patches that originate from "
                     "entirely unrelated simulations.");

    DeclException2(ExcInvalidVectorDeclaration,
                   int,
                   std::string,
                   << "When declaring that a number of components in a data "
                   << "set to be output logically form a vector instead of "
                   << "simply a set of scalar fields, you need to specify "
                   << "this for all relevant components. Furthermore, "
                   << "vectors must always consist of exactly <dim> "
                   << "components. However, the vector component at "
                   << "position " << arg1 << " with name <" << arg2
                   << "> does not satisfy these conditions.");

    DeclException2(ExcInvalidTensorDeclaration,
                   int,
                   std::string,
                   << "When declaring that a number of components in a data "
                   << "set to be output logically form a tensor instead of "
                   << "simply a set of scalar fields, you need to specify "
                   << "this for all relevant components. Furthermore, "
                   << "tensors must always consist of exactly <dim*dim> "
                   << "components. However, the tensor component at "
                   << "position " << arg1 << " with name <" << arg2
                   << "> does not satisfy these conditions.");

  } // namespace DataOutImplementation
} // namespace Exceptions


namespace internal
{
  namespace DataOutImplementation
  {
    /**
     * DataEntry类抽象了用户可以附加到DataOut（和类似的）对象的向量的具体数据类型，并允许底层的DataOut函数查询解决方案向量的单个元素，而不必知道具体的向量类型。这避免了DataOut必须知道正在使用什么向量，但它的缺点是DataOut也不知道这些向量的底层标量类型。
     * 如果底层标量类型都代表实数（在数学意义上
     *
     * 即标量类型是 @p float,   @p double,
     * 等），那么这就不是一个问题。
     *
     * - DataOut只是接收单个向量组件的值作为 @p double 对象。另一方面，如果矢量类型使用 std::complex 标量类型，那么DataEntry为一个矢量条目返回 @p double 是不够的。
     *
     * - 我们需要为DataOut提供一种方法来查询实部和虚部，以便它们可以分别写入输出文件。        这个枚举允许DataOut告诉DataEntry函数它想查询矢量条目的哪个部分，也就是说，它想查询矢量条目的实部或虚部。
     *
     */
    enum class ComponentExtractor
    {
      real_part,
      imaginary_part
    };


    /**
     * 对于每个通过add_data_vector()函数添加的向量，我们需要跟踪它的一个指针，并在我们生成补丁时允许从中提取数据。不幸的是，我们需要对许多不同的向量类型做这件事。幸运的是，它们都有相同的接口。所以我们的方法是有一个基类，提供访问向量信息的函数，并有一个派生模板类，可以为每个向量类型实例化。
     * 由于向量都有相同的接口，这不是什么大问题，因为它们都可以使用相同的通用模板化代码。
     * @note  这个类是<a
     * href="https://www.artima.com/cppsource/type_erasure.html">type
     * erasure</a>设计模式的一个例子。
     *
     */
    template <int dim, int spacedim>
    class DataEntryBase
    {
    public:
      /**
       * 构造函数。给出一个矢量的各个组成部分的名称列表，以及它们作为标量或矢量数据的解释。这个构造函数假定不使用后处理程序。
       *
       */
      DataEntryBase(const DoFHandler<dim, spacedim> *dofs,
                    const std::vector<std::string> & names,
                    const std::vector<
                      DataComponentInterpretation::DataComponentInterpretation>
                      &data_component_interpretation);

      /**
       * 当要使用数据后处理程序时的构造函数。在这种情况下，名称和向量的声明将从后处理器中获得。
       *
       */
      DataEntryBase(const DoFHandler<dim, spacedim> *  dofs,
                    const DataPostprocessor<spacedim> *data_postprocessor);

      /**
       * 解构器变为虚拟。
       *
       */
      virtual ~DataEntryBase() = default;

      /**
       * 假设存储的向量是一个单元格向量，从其中提取给定的元素。
       *
       */
      virtual double
      get_cell_data_value(const unsigned int       cell_number,
                          const ComponentExtractor extract_component) const = 0;

      /**
       * 给定一个FEValuesBase对象，从我们实际存储的向量中提取当前单元格上的值。
       *
       */
      virtual void
      get_function_values(const FEValuesBase<dim, spacedim> &fe_patch_values,
                          const ComponentExtractor           extract_component,
                          std::vector<double> &patch_values) const = 0;

      /**
       * 给定一个FEValuesBase对象，从我们实际存储的向量中提取当前单元格上的值。这个函数的作用与上面的函数相同，但对矢量值的有限元而言。
       *
       */
      virtual void
      get_function_values(
        const FEValuesBase<dim, spacedim> &  fe_patch_values,
        const ComponentExtractor             extract_component,
        std::vector<dealii::Vector<double>> &patch_values_system) const = 0;

      /**
       * 给定一个FEValuesBase对象，从我们实际存储的向量中提取当前单元上的梯度。
       *
       */
      virtual void
      get_function_gradients(
        const FEValuesBase<dim, spacedim> &fe_patch_values,
        const ComponentExtractor           extract_component,
        std::vector<Tensor<1, spacedim>> & patch_gradients) const = 0;

      /**
       * 给定一个FEValuesBase对象，从我们实际存储的向量中提取当前单元格上的梯度。这个函数与上面的函数相同，但对矢量值的有限元而言。
       *
       */
      virtual void
      get_function_gradients(const FEValuesBase<dim, spacedim> &fe_patch_values,
                             const ComponentExtractor extract_component,
                             std::vector<std::vector<Tensor<1, spacedim>>>
                               &patch_gradients_system) const = 0;

      /**
       * 给定一个FEValuesBase对象，从我们实际存储的向量中提取当前单元上的二阶导数。
       *
       */
      virtual void
      get_function_hessians(
        const FEValuesBase<dim, spacedim> &fe_patch_values,
        const ComponentExtractor           extract_component,
        std::vector<Tensor<2, spacedim>> & patch_hessians) const = 0;

      /**
       * 给定一个FEValuesBase对象，从我们实际存储的向量中提取当前单元格上的二阶导数。这个函数的作用与上面的函数相同，但对矢量值的有限元而言。
       *
       */
      virtual void
      get_function_hessians(const FEValuesBase<dim, spacedim> &fe_patch_values,
                            const ComponentExtractor extract_component,
                            std::vector<std::vector<Tensor<2, spacedim>>>
                              &patch_hessians_system) const = 0;

      /**
       * 返回此对象（的派生类）所代表的数据是否代表复值（而不是实值）信息。
       *
       */
      virtual bool
      is_complex_valued() const = 0;

      /**
       * 清除对向量的所有引用。
       *
       */
      virtual void
      clear() = 0;

      /**
       * 确定这个对象的内存消耗（以字节为单位）的估计值。
       *
       */
      virtual std::size_t
      memory_consumption() const = 0;

      /**
       * 指向该向量所基于的DoFHandler对象的指针。
       *
       */
      SmartPointer<const DoFHandler<dim, spacedim>> dof_handler;

      /**
       * 这个数据向量的组成部分的名称。
       *
       */
      const std::vector<std::string> names;

      /**
       * 一个向量，对于当前数据集的n_output_variables变量中的每一个，表明它们是标量场、矢量场的一部分或任何其他支持的数据种类。
       *
       */
      const std::vector<
        DataComponentInterpretation::DataComponentInterpretation>
        data_component_interpretation;

      /**
       * 指向数据后处理对象的指针，该对象将被应用于该数据向量。
       *
       */
      SmartPointer<const dealii::DataPostprocessor<spacedim>> postprocessor;

      /**
       * 这个数据集提供的输出变量的数量（如果应用了DataPostprocessor，可以是向量值函数/数据向量中的组件数量，也可以是计算量的数量）。这个变量是通过<tt>names.size()</tt>确定的，因此等同于<tt>names.size()</tt>。
       *
       */
      unsigned int n_output_variables;
    };


    /**
     * 一个数据结构，在并行构建补丁时，在一个线程中保存所有需要的数据。这些数据结构是全局性的，而不是在每个单元上创建的，以避免线程中的内存分配。这是在WorkStream类的文档中讨论的AdditionalData类数据结构的基类。
     * <code>cell_to_patch_index_map</code>
     * 是一个数组，为索引<tt>[i][j]</tt>存储与索引为 @p j
     * 的单元格相关的补丁的编号，在级别 @p i.
     * 上，该信息是在生成补丁之前设置的，需要用来生成邻接信息。
     * 这个结构被几个DataOut*类使用，这些类从它那里派生出自己的ParallelData类，用于附加字段。
     *
     */
    template <int dim, int spacedim>
    struct ParallelDataBase
    {
      ParallelDataBase(
        const unsigned int               n_datasets,
        const unsigned int               n_subdivisions,
        const std::vector<unsigned int> &n_postprocessor_outputs,
        const Mapping<dim, spacedim> &   mapping,
        const std::vector<
          std::shared_ptr<dealii::hp::FECollection<dim, spacedim>>>
          &               finite_elements,
        const UpdateFlags update_flags,
        const bool        use_face_values);

      ParallelDataBase(
        const unsigned int               n_datasets,
        const unsigned int               n_subdivisions,
        const std::vector<unsigned int> &n_postprocessor_outputs,
        const dealii::hp::MappingCollection<dim, spacedim> &mapping,
        const std::vector<
          std::shared_ptr<dealii::hp::FECollection<dim, spacedim>>>
          &               finite_elements,
        const UpdateFlags update_flags,
        const bool        use_face_values);

      ParallelDataBase(const ParallelDataBase &data);

      void
      reinit_all_fe_values(
        std::vector<std::shared_ptr<DataEntryBase<dim, spacedim>>> &dof_data,
        const typename dealii::Triangulation<dim, spacedim>::cell_iterator
          &                cell,
        const unsigned int face = numbers::invalid_unsigned_int);

      const FEValuesBase<dim, spacedim> &
      get_present_fe_values(const unsigned int dataset) const;

      void
      resize_system_vectors(const unsigned int n_components);

      const unsigned int n_datasets;
      const unsigned int n_subdivisions;

      DataPostprocessorInputs::Scalar<spacedim>        patch_values_scalar;
      DataPostprocessorInputs::Vector<spacedim>        patch_values_system;
      std::vector<std::vector<dealii::Vector<double>>> postprocessed_values;

      const dealii::hp::MappingCollection<dim, spacedim> mapping_collection;
      const std::vector<
        std::shared_ptr<dealii::hp::FECollection<dim, spacedim>>>
                        finite_elements;
      const UpdateFlags update_flags;

      std::vector<std::shared_ptr<dealii::hp::FEValues<dim, spacedim>>>
        x_fe_values;
      std::vector<std::shared_ptr<dealii::hp::FEFaceValues<dim, spacedim>>>
        x_fe_face_values;
    };
  } // namespace DataOutImplementation
} // namespace internal


// TODO: Most of the documentation of DataOut_DoFData applies to DataOut.

/**
 * 这是一个抽象的类，它提供了从网格上的数据向量生成基类输出的补丁的功能。它允许将一个或多个指针附加到DoFHandler上，并附加节点和单元数据，表示网格上的功能，这些功能以后将以任何实现的数据格式写入。
 *
 *  <h3>User visible interface</h3>
 * 该类的用户可见界面允许用户以两种不同的方式指定数据。一种是让这个类知道一个DoFHandler对象，并添加数据向量，这些数据向量都对应于这个DoFHandler或网格单元，这些网格单元以后将以某种格式写入文件中。第二种方法是将一个DoFHandler对象与向量一起传递。这允许以一种整洁的方式设置来自不同DoFHandler的数据（当然，它们都需要基于相同的三角测量）。与其思考不同的函数，第一种的例子可能是最好的解释。
 *
 * @code
 * ...
 * ...   // compute solution, which contains nodal values
 * ...
 * ...   // compute error_estimator, which contains one value per cell
 *
 * std::vector<std::string> solution_names;
 * solution_names.emplace_back ("x-displacement");
 * solution_names.emplace_back ("y-displacement");
 *
 * DataOut<dim> data_out;
 * data_out.attach_dof_handler (dof_handler);
 * data_out.add_data_vector (solution, solution_names);
 * data_out.add_data_vector (error_estimator, "estimated_error");
 *
 * data_out.build_patches ();
 *
 * ofstream output_file ("output");
 * data_out.write_xxx (output_file);
 *
 * data_out.clear();
 * @endcode
 *
 * attach_dof_handler()告诉这个类，所有未来的操作都要通过DoFHandler对象和它所在的三角形来进行。然后我们添加解决方案向量和误差估计器；注意它们有不同的尺寸，因为解决方案是一个节点向量，这里由两个分量（"x-位移
 * "和
 * "y-位移"）组成，而误差估计器可能是一个容纳单元数据的向量。当附加一个数据向量时，你必须给向量的每个分量一个名字，这是通过一个<tt>vector<string>/tt>类型的对象作为第二个参数来完成的；如果向量中只有一个分量，例如我们像第二种情况那样添加单元数据，或者如果DoFHandler使用的有限元只有一个分量，那么你可以使用第二个add_data_vector()函数，它需要一个
 * @p string  ] 而不是<tt>vector<string></tt>。
 * add_data_vector()函数有额外的参数（有默认值），可以用来指定某些转换。特别是，它允许附加DataPostprocessor参数，以便在场将被评估的每一点上从数据向量中计算出派生信息，从而可以将其写入文件（例如，可以从密度和速度中计算出高超声速流的马赫数； step-29  ]也展示了一个例子）；另一个通过参数指定的默认值的信息是如何解释某些输出成分，即数据的每个成分在逻辑上是否是一个独立的标量场，或者其中一些成分在逻辑上是否形成一个矢量场（见 DataComponentInterpretation::DataComponentInterpretation 枚举，以及 @ref step_22  "  step-22  "
 * 教程程序）。
 * 由于内存消耗的原因，这个类不会通过add_data_vector()函数复制给它的向量。它只存储一个引用，所以你有责任确保数据向量存在足够长的时间。
 * 在添加完所有的数据向量后，你需要调用一个函数，该函数生成补丁（即一些中间数据表示），用于从存储的数据中输出。派生类将这个函数命名为build_patches()。最后，你以一种格式或其他方式将数据写入（）到一个文件中。
 * 在上面的例子中，使用了一个DataOut类型的对象，也就是一个派生类的对象。这是必要的，因为当前的类并没有提供实际生成补丁的方法，只是提供了存储和访问数据的帮助。任何真正的功能都在派生类中实现，如DataOut。
 * 请注意，这个类的基类，DataOutInterface提供了几个函数，以缓解运行时可确定的输出格式的编程（即你不需要通过调用上面例子中的
 * DataOutInterface::write_xxx
 * 来使用固定的格式，但你可以通过运行时参数来选择它，而不必自己写<tt>if
 * () ... else ...
 * </tt>条款），还有一些函数和类提供了通过为每个输出格式设置标志来控制输出外观的方法。
 *
 *  <h3>Information for derived classes</h3>
 * 这个类所缺乏的是一种从存储的数据和自由度信息中产生输出补丁的方法。因为这个任务通常是依赖于应用的，所以它被留给了派生类。例如，在许多应用中，可能希望将输出的深度限制在一定数量的细化水平上，并且仅以插值的方式将数据从较细的单元写入较粗的单元，以减少输出的数量。另外，在形成一个补丁时，可能希望在不同的单元上使用不同数量的细分，例如在不同的单元上完成对试验空间的不同多项式程度的划分。另外，输出不一定由每个单元的补丁组成，而可能是由面的补丁组成的，还有其他的东西。看看派生类在这方面有什么可能。
 * 由于这个原因，它留给派生类提供一个函数，名称通常是build_patches()或类似的，用来填充这个类的#patches数组。
 * 关于这个类的模板，它需要三个值：首先是三角形和DoF处理程序操作的空间维度，其次是补丁所代表的对象的维度。
 * 虽然在大多数情况下它们是相等的，但也有一些类不成立，例如，如果一个人输出利用原域的旋转对称性的计算结果（在这种情况下，输出的空间维度将比DoF处理程序的维度高一个，见DataOut_Rotation()类），或者我们可以设想，我们可以写一个类，只输出在域的切面上的解决方案，在这种情况下，输出的空间维度小于DoF处理程序的。最后一个模板参数表示嵌入补丁的空间维度；通常，这个维度与补丁本身的维度相同（这也是模板参数的默认值），但可能有一些情况并非如此。例如，在DataOut_Faces()类中，补丁是由三角形的面生成的。因此，补丁的维度比嵌入空间的维度少一个，在这种情况下，嵌入空间的维度等于三角形和DoF处理器的维度。然而，对于上述领域的切口，如果切口是直的，那么切口可以被嵌入到一个比三角形维度低一个维度的空间中，因此最后一个模板参数的值与第二个参数相同。
 *
 *
 * @ingroup output
 *
 *
 */
template <int dim,
          int patch_dim,
          int spacedim       = dim,
          int patch_spacedim = patch_dim>
class DataOut_DoFData : public DataOutInterface<patch_dim, patch_spacedim>
{
public:
  /**
   * 对所考虑的dof处理程序类的迭代器类型的类型定义。
   *
   */
  using cell_iterator = typename Triangulation<dim, spacedim>::cell_iterator;

public:
  /**
   * 描述给add_data_vector()的向量是什么的类型：在DoFHandler对象中每个自由度有一个条目的向量（如解决方案向量），或者在DoFHandler对象基础的三角形中每个单元有一个条目（如每个单元的误差数据）。值#type_automatic告诉add_data_vector()自己找出来（关于使用的方法，见add_data_vector()的文档）。
   *
   */
  enum DataVectorType
  {
    /**
     * 数据向量条目与自由度有关
     *
     */
    type_dof_data,

    /**
     * 数据向量条目为每个网格单元一个
     *
     */
    type_cell_data,

    /**
     * 自动找出
     *
     */
    type_automatic
  };

  /**
   * 构造函数
   *
   */
  DataOut_DoFData();

  /**
   * 解构器。
   *
   */
  virtual ~DataOut_DoFData() override;

  /**
   * 指定一个dof处理程序，用来提取几何数据以及节点和节点值之间的映射。如果所有添加的数据向量都补充了DoFHandler参数，那么这个调用就没有必要。
   * 这个调用是可选的：如果你用指定的DoFHandler对象添加数据向量，那么这就包含了生成输出所需的所有信息。
   *
   */
  void
  attach_dof_handler(const DoFHandler<dim, spacedim> &);

  /**
   * 指定一个用于提取几何数据和节点与节点值之间的映射的三角测量。
   * 这个调用是可选的：如果你用指定的DoFHandler对象添加数据向量，那么这就包含了生成输出所需的所有信息。
   * 当你只输出单元格向量而完全没有DoFHandler时，这个调用很有用，在这种情况下，它提供了几何图形。
   *
   */
  void
  attach_triangulation(const Triangulation<dim, spacedim> &);

  /**
   * 添加一个数据向量和它的名字。
   * 一个指向该向量的指针被存储，所以你必须确保该向量至少在你调用<tt>write_*</tt>函数时存在于该地址。
   * 假设该向量的分量与自由度处理程序中的自由度数量相同，在这种情况下，它被认为是一个存储节点数据的向量；或者大小可能是目前网格上活动单元的数量，在这种情况下，它被认为是一个单元数据向量。由于自由度和单元的数量通常不相等，函数可以自行决定给定哪种类型的向量。然而，在一些角落里，这种自动判断并不奏效。
   * 一个例子是如果你用片状常数元素计算，并且有一个标量解决方案，那么有多少个单元就有多少个自由度（尽管它们的编号可能不同）。
   * 另一种可能是，如果你有一个嵌入2D空间的1D网格，并且该网格由单元的封闭曲线组成；在这种情况下，有多少节点就有多少单元，当使用Q1元素时，你会有多少自由度就有多少单元。
   * 在这种情况下，你可以将函数的最后一个参数从默认值#type_automatic改为#type_dof_data或#type_cell_data，这取决于矢量代表什么。除了这种角落里的情况，你可以把参数留在默认值上，让函数决定向量本身的类型。
   * 如果它是一个持有DoF数据的向量，给出的名称应是底层有限元的每个分量。
   * 如果它是一个仅由一个子元素组成的有限元，那么下面还有一个函数，它接受一个单一的名字而不是一个名字的向量。
   * data_component_interpretation参数包含关于如何解释由多个数据集组成的输出文件的各个组成部分的信息。
   * 例如，如果一个人有一个2D的斯托克斯方程的有限元，代表组件（u,v,p），我们希望表明前两个，u和v，代表一个逻辑矢量，这样以后当我们生成图形输出时，我们可以把它们交给一个可视化程序，该程序将自动知道把它们作为一个矢量场来渲染，而不是作为两个独立的标量场。
   * 这个参数的默认值（即一个空的矢量）对应于一个值的矢量
   * DataComponentInterpretation::component_is_scalar,
   * ，表示所有的输出组件都是独立的标量场。然而，如果给定的数据向量代表逻辑向量，你可以传入一个包含数值
   * DataComponentInterpretation::component_is_part_of_vector. 的向量。
   * 在上面的例子中，人们会传入一个包含(u,v,p)的组件
   * (DataComponentInterpretation::component_is_part_of_vector,
   * DataComponentInterpretation::component_is_part_of_vector,
   * DataComponentInterpretation::component_is_scalar) 的向量。
   * 数据向量的名称只能包含字母、下划线和其他一些字符。请参考本类中声明的ExcInvalidCharacter异常，看看哪些字符是有效的，哪些是无效的。
   * @note
   * 矢量参数的实际类型可以是任何矢量类型，FEValues可以使用
   * FEValuesBase::get_function_values() 函数从单元格上提取数值。
   * @note
   * 当并行工作时，要写入的向量需要对本地拥有的单元上的所有自由度进行读取访问的幽灵化，详见
   * step-40 或 step-37
   * 教程程序，即可能需要调用data.update_ghost_values（）。
   *
   */
  template <class VectorType>
  void
  add_data_vector(
    const VectorType &              data,
    const std::vector<std::string> &names,
    const DataVectorType            type = type_automatic,
    const std::vector<DataComponentInterpretation::DataComponentInterpretation>
      &data_component_interpretation = std::vector<
        DataComponentInterpretation::DataComponentInterpretation>());

  /**
   * 这个函数是上面那个函数的缩写（关于各种参数的讨论见那里），旨在用于不是由子元素组成的有限元。在这种情况下，每个数据向量只需要给出一个名称，这就是这个函数的作用。它只是在将
   * @p name
   * 转换为字符串矢量后，将其参数转发给上面的另一个add_data_vector()函数。
   * 如果 @p data
   * 是一个有多个成分的向量，这个函数将通过在 @p name
   * 中附加下划线和每个成分的编号来为所有成分生成不同的名称。
   * 模板参数的实际类型可以是任何向量类型，FEValues可以使用
   * FEValuesBase::get_function_values() 函数在单元格中提取值。
   *
   */
  template <class VectorType>
  void
  add_data_vector(
    const VectorType &   data,
    const std::string &  name,
    const DataVectorType type = type_automatic,
    const std::vector<DataComponentInterpretation::DataComponentInterpretation>
      &data_component_interpretation = std::vector<
        DataComponentInterpretation::DataComponentInterpretation>());

  /**
   * 这个函数是上面那个函数的扩展（除了第一个参数，其他参数的讨论见那里），允许用自己的DoFHandler对象设置一个向量。这个DoFHandler需要与调用
   * @p  add_data_vector或 @p attach_dof_handler,
   * 分配的其他DoFHandler对象兼容，即所有的DoFHandler对象都需要基于相同的三角测量。这个函数允许你从描述不同解决方案组件的多个DoFHandler对象导出数据。在
   * step-61  中给出了一个使用此函数的例子。
   * 由于这个函数接受一个DoFHandler对象，因此自然地代表了dof数据，上面其他方法中出现的数据矢量类型参数就没有必要了。
   *
   */
  template <class VectorType>
  void
  add_data_vector(
    const DoFHandler<dim, spacedim> &dof_handler,
    const VectorType &               data,
    const std::vector<std::string> & names,
    const std::vector<DataComponentInterpretation::DataComponentInterpretation>
      &data_component_interpretation = std::vector<
        DataComponentInterpretation::DataComponentInterpretation>());


  /**
   * 这个函数是上述函数的缩写，只给出一个标量 @p
   * dof_handler 和一个数据名称。
   *
   */
  template <class VectorType>
  void
  add_data_vector(
    const DoFHandler<dim, spacedim> &dof_handler,
    const VectorType &               data,
    const std::string &              name,
    const std::vector<DataComponentInterpretation::DataComponentInterpretation>
      &data_component_interpretation = std::vector<
        DataComponentInterpretation::DataComponentInterpretation>());

  /**
   * 这个函数是上述函数的一个替代品，允许输出派生量而不是给定的数据。这种转换必须在一个派生自DataPostprocessor的类中完成。这个函数在
   * step-29  中使用。其他用途见  step-32  和  step-33  。
   * 这些派生量的名称由 @p
   * data_postprocessor参数提供。同样，其他add_data_vector()函数的data_component_interpretation参数也是由data_postprocessor参数提供的。由于只有类型为
   * @p type_dof_data
   * 的数据可以被转换，这个类型也是隐含地知道的，不需要给出。
   * @note
   * 矢量参数的实际类型可以是任何矢量类型，FEValues可以使用
   * FEValuesBase::get_function_values() 函数从单元格上提取数值。
   * @note
   * 数据后处理器对象（即实际上是你的派生类的对象）必须活到DataOut对象被销毁为止，因为后者保持着一个指向前者的指针，如果指向的对象被销毁，而后者仍有一个指向它的指针，就会抱怨。如果数据后处理器和DataOut对象都是一个函数的局部变量（例如，在
   * step-29
   * 中就是如此），那么你可以通过在DataOut变量之前声明数据后处理器变量来避免这个错误，因为对象的销毁顺序与声明顺序相反。
   *
   */
  template <class VectorType>
  void
  add_data_vector(const VectorType &                 data,
                  const DataPostprocessor<spacedim> &data_postprocessor);

  /**
   * 与上面的函数相同，但有一个DoFHandler对象，不需要与最初设置的DoFHandler重合。注意，后处理器只能从给定的DoFHandler和求解向量中读取数据，而不是其他求解向量或DoFHandler。
   *
   */
  template <class VectorType>
  void
  add_data_vector(const DoFHandler<dim, spacedim> &  dof_handler,
                  const VectorType &                 data,
                  const DataPostprocessor<spacedim> &data_postprocessor);

  /**
   * 添加一个多级数据向量。    这个函数将属于DoFHandler  @p
   * dof_handler的每一级的向量形式的向量值多级向量 @p data
   * 添加到图形输出中。这个函数通常与调用set_cell_selection()一起使用，该函数选择的是特定层次上的单元，而不是活动单元（默认）。
   * 矢量 @p data
   * 可以通过几种方式获得，例如在多网格循环期间或之后使用
   * Multigrid::solution 或 Multigrid::defect ，或者通过
   * MGTransferMatrixFree::interpolate_to_mg(). 对解进行插值处理 @p
   * names 和 @p data_component_interpretation
   * 与add_data_vector（）函数相同。
   *
   */
  template <class VectorType>
  void
  add_mg_data_vector(
    const DoFHandler<dim, spacedim> &dof_handler,
    const MGLevelObject<VectorType> &data,
    const std::vector<std::string> & names,
    const std::vector<DataComponentInterpretation::DataComponentInterpretation>
      &data_component_interpretation = std::vector<
        DataComponentInterpretation::DataComponentInterpretation>());

  /**
   * 上述函数的标量版本。
   *
   */
  template <class VectorType>
  void
  add_mg_data_vector(const DoFHandler<dim, spacedim> &dof_handler,
                     const MGLevelObject<VectorType> &data,
                     const std::string &              name);

  /**
   * 释放指向数据向量的指针。这允许输出一组新的向量，而无需再次提供DoF处理程序。因此，DataOut对象可以在代数背景下使用。注意，除了数据向量，已经计算的补丁也会被删除。
   *
   */
  void
  clear_data_vectors();

  /**
   * 释放指向所有输入数据元素的指针，即指向数据向量和DoF处理器对象的指针。当你调用了派生类的
   * @p build_patches
   * 函数时，这个函数可能很有用，因为此时补丁已经建立，不再需要输入数据，也不需要引用它。然后你就可以从主线程中分离出来输出补丁，而不需要再确保在输出线程结束之前，DoF处理程序对象和向量不能被删除。
   *
   */
  void
  clear_input_data_references();

  /**
   * 这个函数可以用来将使用作为参数的对象的 @p
   * build_patches
   * 函数创建的补丁合并到这个对象创建的补丁列表中。例如，如果有一个领域分解算法，其中每个块都由它自己的DoFHandler表示，但人们想同时输出所有块上的解决方案，这有时是很方便的。
   * 要做到这一点，给定的参数和这个对象需要有相同数量的输出向量，并且它们需要使用相同数量的每个补丁的细分。如果两个对象中的补丁在空间上有重叠，那么输出结果可能会看起来相当有趣。
   * 如果你在合并补丁后为这个对象调用build_patches()，之前的状态会被覆盖，而合并的补丁会丢失。
   * 第二个参数允许将第一个参数中传递的对象中的补丁的每个节点移动一定量。这对于生成一个区块集合的
   * "爆炸 "视图有时是有用的。
   * 如果这个对象或另一个对象还没有设置任何补丁，这个函数将会失败。
   *
   */
  template <int dim2, int spacedim2>
  void
  merge_patches(
    const DataOut_DoFData<dim2, patch_dim, spacedim2, patch_spacedim> &source,
    const Point<patch_spacedim> &shift = Point<patch_spacedim>());

  /**
   * @deprecated  用merge_patches()代替DoFHandlerType2模板。
   *
   */
  template <typename DoFHandlerType2>
  DEAL_II_DEPRECATED void
  merge_patches(const DataOut_DoFData<DoFHandlerType2::dimension,
                                      patch_dim,
                                      DoFHandlerType2::space_dimension,
                                      patch_spacedim> &source,
                const Point<patch_spacedim> &shift = Point<patch_spacedim>());

  /**
   * 释放指向数据向量和DoF处理程序的指针。你必须使用add_data_vector()函数重新设置所有的数据条目。Dof处理程序的指针也被清空，连同所有其他数据一起。
   * 实际上，这个函数把所有的东西都重设为一个处女状态。
   *
   */
  virtual void
  clear();

  /**
   * 确定这个对象的内存消耗（以字节为单位）的估计值。
   *
   */
  std::size_t
  memory_consumption() const;

protected:
  /**
   * 缩写有点冗长的Patch类的名字。
   *
   */
  using Patch = dealii::DataOutBase::Patch<patch_dim, patch_spacedim>;

  /**
   * 指向三角测量对象的指针。
   *
   */
  SmartPointer<const Triangulation<dim, spacedim>> triangulation;

  /**
   * 指向可选的处理程序对象的指针。
   *
   */
  SmartPointer<const DoFHandler<dim, spacedim>> dofs;

  /**
   * 带有每个自由度数值向量的数据元素的列表。
   *
   */
  std::vector<std::shared_ptr<
    internal::DataOutImplementation::DataEntryBase<dim, spacedim>>>
    dof_data;

  /**
   * 带有每个单元值向量的数据元素的列表。
   *
   */
  std::vector<std::shared_ptr<
    internal::DataOutImplementation::DataEntryBase<dim, spacedim>>>
    cell_data;

  /**
   * 这是一个补丁列表，每次调用build_patches()时都会创建一个补丁。这些补丁在基类的输出例程中使用。
   *
   */
  std::vector<Patch> patches;

  /**
   * %函数，基类的函数通过这个函数知道他们应该把哪些补丁写到文件中。
   *
   */
  virtual const std::vector<Patch> &
  get_patches() const override;

  /**
   * 虚拟函数，基类的输出函数通过它获得数据集的名称。
   *
   */
  virtual std::vector<std::string>
  get_dataset_names() const override;

  /**
   * 提取存储在dof_data对象中的有限元，包括一个FE_DGQ<dim>(0)的假对象，以防只使用三角法。
   *
   */
  std::vector<std::shared_ptr<dealii::hp::FECollection<dim, spacedim>>>
  get_fes() const;

  /**
   * 相关的 DataOutInterface::get_nonscalar_data_ranges()
   * 函数的重载。参见那里有更多的文档。
   *
   */
  virtual std::vector<
    std::tuple<unsigned int,
               unsigned int,
               std::string,
               DataComponentInterpretation::DataComponentInterpretation>>
  get_nonscalar_data_ranges() const override;

  // Make all template siblings friends. Needed for the merge_patches()
  // function.
  template <int, int, int, int>
  friend class DataOut_DoFData;

   /**
    */
  template <int, class>
  friend class MGDataOut;

private:
  /**
   * 由四个公共的add_data_vector方法调用的通用函数。
   *
   */
  template <class VectorType>
  void
  add_data_vector_internal(
    const DoFHandler<dim, spacedim> *dof_handler,
    const VectorType &               data,
    const std::vector<std::string> & names,
    const DataVectorType             type,
    const std::vector<DataComponentInterpretation::DataComponentInterpretation>
      &        data_component_interpretation,
    const bool deduce_output_names);
};



// -------------------- template and inline functions ------------------------
template <int dim, int patch_dim, int spacedim, int patch_spacedim>
template <typename VectorType>
void
DataOut_DoFData<dim, patch_dim, spacedim, patch_spacedim>::add_data_vector(
  const VectorType &   vec,
  const std::string &  name,
  const DataVectorType type,
  const std::vector<DataComponentInterpretation::DataComponentInterpretation>
    &data_component_interpretation)
{
  Assert(triangulation != nullptr,
         Exceptions::DataOutImplementation::ExcNoTriangulationSelected());
  std::vector<std::string> names(1, name);
  add_data_vector_internal(
    dofs, vec, names, type, data_component_interpretation, true);
}



template <int dim, int patch_dim, int spacedim, int patch_spacedim>
template <typename VectorType>
void
DataOut_DoFData<dim, patch_dim, spacedim, patch_spacedim>::add_data_vector(
  const VectorType &              vec,
  const std::vector<std::string> &names,
  const DataVectorType            type,
  const std::vector<DataComponentInterpretation::DataComponentInterpretation>
    &data_component_interpretation)
{
  Assert(triangulation != nullptr,
         Exceptions::DataOutImplementation::ExcNoTriangulationSelected());
  add_data_vector_internal(
    dofs, vec, names, type, data_component_interpretation, false);
}



template <int dim, int patch_dim, int spacedim, int patch_spacedim>
template <typename VectorType>
void
DataOut_DoFData<dim, patch_dim, spacedim, patch_spacedim>::add_data_vector(
  const DoFHandler<dim, spacedim> &dof_handler,
  const VectorType &               data,
  const std::string &              name,
  const std::vector<DataComponentInterpretation::DataComponentInterpretation>
    &data_component_interpretation)
{
  std::vector<std::string> names(1, name);
  add_data_vector_internal(&dof_handler,
                           data,
                           names,
                           type_dof_data,
                           data_component_interpretation,
                           true);
}



template <int dim, int patch_dim, int spacedim, int patch_spacedim>
template <typename VectorType>
void
DataOut_DoFData<dim, patch_dim, spacedim, patch_spacedim>::add_data_vector(
  const DoFHandler<dim, spacedim> &dof_handler,
  const VectorType &               data,
  const std::vector<std::string> & names,
  const std::vector<DataComponentInterpretation::DataComponentInterpretation>
    &data_component_interpretation)
{
  add_data_vector_internal(&dof_handler,
                           data,
                           names,
                           type_dof_data,
                           data_component_interpretation,
                           false);
}



template <int dim, int patch_dim, int spacedim, int patch_spacedim>
template <typename VectorType>
void
DataOut_DoFData<dim, patch_dim, spacedim, patch_spacedim>::add_data_vector(
  const VectorType &                 vec,
  const DataPostprocessor<spacedim> &data_postprocessor)
{
  Assert(dofs != nullptr,
         Exceptions::DataOutImplementation::ExcNoDoFHandlerSelected());
  add_data_vector(*dofs, vec, data_postprocessor);
}



template <int dim, int patch_dim, int spacedim, int patch_spacedim>
template <int dim2, int spacedim2>
void
DataOut_DoFData<dim, patch_dim, spacedim, patch_spacedim>::merge_patches(
  const DataOut_DoFData<dim2, patch_dim, spacedim2, patch_spacedim> &source,
  const Point<patch_spacedim> &                                      shift)
{
  const std::vector<Patch> &source_patches = source.get_patches();
  Assert((patches.size() != 0) && (source_patches.size() != 0),
         ExcMessage("When calling this function, both the current "
                    "object and the one being merged need to have a "
                    "nonzero number of patches associated with it. "
                    "Either you called this function on objects that "
                    "are empty, or you may have forgotten to call "
                    "the 'build_patches()' function."));
  // Check equality of component names
  Assert(get_dataset_names() == source.get_dataset_names(),
         Exceptions::DataOutImplementation::ExcIncompatibleDatasetNames());

  // Make sure patches are compatible. Ideally, we would check that all input
  // patches from both collections are all compatible, but we'll be content
  // with checking that just the first ones from both sources are.
  //
  // We check compatibility by testing that both sets of patches result
  // from the same number of subdivisions, and that they have the same
  // number of source vectors (they really should, since we already checked
  // that there are the same number of source components above, but you
  // never know). This implies that the data should have the same number of
  // columns. They should really have the same number of rows as well,
  // but depending on whether a patch has points included or not, the
  // number of rows may or may not include coordinates for the points,
  // and the comparison has to account for that because in each source
  // stream, the patches may include some that have points included.
  Assert(patches[0].n_subdivisions == source_patches[0].n_subdivisions,
         Exceptions::DataOutImplementation::ExcIncompatiblePatchLists());
  Assert(patches[0].data.n_cols() == source_patches[0].data.n_cols(),
         Exceptions::DataOutImplementation::ExcIncompatiblePatchLists());
  Assert((patches[0].data.n_rows() +
          (patches[0].points_are_available ? 0 : patch_spacedim)) ==
           (source_patches[0].data.n_rows() +
            (source_patches[0].points_are_available ? 0 : patch_spacedim)),
         Exceptions::DataOutImplementation::ExcIncompatiblePatchLists());

  // check equality of the vector data
  // specifications
  Assert(get_nonscalar_data_ranges().size() ==
           source.get_nonscalar_data_ranges().size(),
         ExcMessage("Both sources need to declare the same components "
                    "as vectors."));
  for (unsigned int i = 0; i < get_nonscalar_data_ranges().size(); ++i)
    {
      Assert(std::get<0>(get_nonscalar_data_ranges()[i]) ==
               std::get<0>(source.get_nonscalar_data_ranges()[i]),
             ExcMessage("Both sources need to declare the same components "
                        "as vectors."));
      Assert(std::get<1>(get_nonscalar_data_ranges()[i]) ==
               std::get<1>(source.get_nonscalar_data_ranges()[i]),
             ExcMessage("Both sources need to declare the same components "
                        "as vectors."));
      Assert(std::get<2>(get_nonscalar_data_ranges()[i]) ==
               std::get<2>(source.get_nonscalar_data_ranges()[i]),
             ExcMessage("Both sources need to declare the same components "
                        "as vectors."));
    }

  // merge patches. store old number
  // of elements, since we need to
  // adjust patch numbers, etc
  // afterwards
  const unsigned int old_n_patches = patches.size();
  patches.insert(patches.end(), source_patches.begin(), source_patches.end());

  // perform shift, if so desired
  if (shift != Point<patch_spacedim>())
    for (unsigned int i = old_n_patches; i < patches.size(); ++i)
      for (const unsigned int v : GeometryInfo<patch_dim>::vertex_indices())
        patches[i].vertices[v] += shift;

  // adjust patch numbers
  for (unsigned int i = old_n_patches; i < patches.size(); ++i)
    patches[i].patch_index += old_n_patches;

  // adjust patch neighbors
  for (unsigned int i = old_n_patches; i < patches.size(); ++i)
    for (const unsigned int n : GeometryInfo<patch_dim>::face_indices())
      if (patches[i].neighbors[n] != Patch::no_neighbor)
        patches[i].neighbors[n] += old_n_patches;
}



template <int dim, int patch_dim, int spacedim, int patch_spacedim>
template <typename DoFHandlerType2>
void
DataOut_DoFData<dim, patch_dim, spacedim, patch_spacedim>::merge_patches(
  const DataOut_DoFData<DoFHandlerType2::dimension,
                        patch_dim,
                        DoFHandlerType2::space_dimension,
                        patch_spacedim> &source,
  const Point<patch_spacedim> &          shift)
{
  this->merge_patches<DoFHandlerType2::dimension,
                      DoFHandlerType2::space_dimension>(source, shift);
}



namespace Legacy
{
  /**
   * @deprecated  使用没有DoFHandlerType模板的 dealii::DataOut_DoFData
   * 代替。
   *
   */
  template <typename DoFHandlerType,
            int patch_dim,
            int patch_space_dim = patch_dim>
  using DataOut_DoFData DEAL_II_DEPRECATED =
    dealii::DataOut_DoFData<DoFHandlerType::dimension,
                            patch_dim,
                            DoFHandlerType::space_dimension,
                            patch_space_dim>;
} // namespace Legacy


DEAL_II_NAMESPACE_CLOSE

#endif


