//include/deal.II-translator/numerics/point_value_history_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2020 by the deal.II authors
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


#ifndef dealii_point_value_history_h
#define dealii_point_value_history_h

#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/point.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/component_mask.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_tools.h>

#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_postprocessor.h>

#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

DEAL_II_NAMESPACE_OPEN

namespace internal
{
  namespace PointValueHistoryImplementation
  {
    /**
     * 一个存储所需数据的类，用于引用最接近一个请求点的支持点。
     *
     */
    template <int dim>
    class PointGeometryData
    {
    public:
      PointGeometryData(
        const Point<dim> &                          new_requested_location,
        const std::vector<Point<dim>> &             new_locations,
        const std::vector<types::global_dof_index> &new_sol_indices);
      Point<dim>                           requested_location;
      std::vector<Point<dim>>              support_point_locations;
      std::vector<types::global_dof_index> solution_indices;
    };
  } // namespace PointValueHistoryImplementation
} // namespace internal



/**
 * PointValueHistory解决了在网格上特定点绘制时间（或任何其他迭代过程）的解值图的开销。用户提前指定需要监测的解的点，并给每个需要记录的解向量一个记忆性的名字。然后，对于每一步，用户调用三个可用的
 * "评估域
 * "方法之一来存储每个时间步骤的数据，该类提取所要求的点的数据来存储。最后，一旦计算完成，用户可以要求生成输出文件；这些文件是Gnuplot格式，但基本上只是普通的文本，可以很容易地导入其他程序中，例如导入电子表格。
 * 用户可以存储与网格位置无关的额外变量，指定n_independent_variables。然后，该类期望在每个步骤中使用
 * @p push_back_independent.
 * 方法添加一个大小为n_independent_variables的 std::vector
 * ，这可用于例如记录外部输入，记录求解器性能数据，如求解步骤和求解器收敛前的步骤所需时间，保存计算的规范，或简单地保存时间、时间步骤数或非线性迭代数，以及从网格评估的数据。
 * 这三种 "评价场 "方法各自有不同的优点和缺点，使其适用于不同的情况。  <ol>   <li>  首先，不采取 @p  DataPostprocessor对象的 @p evaluate_field 版本选择离给定点最近的支持点（见 @ref GlossSupport  "术语表中的这个条目"）来提取数据。这使得每个时间步长需要运行的代码非常短，因为在网格上循环提取所需的dof_index可以在开始时只做一次。然而，这种方法不适合于不将dof分配到实际网格位置的有限元对象（即没有 @ref GlossSupport "支持点 "的FEs），或者使用自适应网格细化的情况。后者的原因是，离给定点最近的支持点的位置可能在网格细化时发生变化。如果三角剖分有任何变化，该类将抛出一个异常（尽管在网格细化时可以重新计算最近的支持点，但支持点的位置很可能会发生微小的变化，使得数据的解释变得困难，因此目前没有实现）。
 * <li>  其次， @p evaluate_field_at_requested_location 调用 @p
 * VectorTools::point_value
 * 来计算所要求的特定点的值。这个方法对任何被 @p
 * VectorTools::point_value.
 * 支持的FE都是有效的。具体来说，这个方法可以被使用自适应网格细化的代码所调用。
 * <li>  最后，该类提供了一个函数 @p evaluate_field
 * ，接收一个 @p
 * DataPostprocessor对象。该方法允许使用deal.II数据后处理程序来从解决方案中即时计算新的数量。这些值位于离请求点最近的正交点上。如果网格在两次调用之间被细化，这个点就会改变，所以在使用自适应细化的代码中使用这个方法时必须小心，但是由于输出是有意义的（在这个意义上，所选择的正交点被保证保持在同一附近，该类并不阻止在自适应代码中使用这个方法。如果网格发生了变化，该类在输出文件中提供警告。请注意，我们可以通过提供一个有更多点的正交公式来减少这个过程所带来的误差，但代价是要做更多的工作，因为此时最接近的正交点更接近于真正应该发生的评估点。(附带说明：为什么不立即在要求的点上进行评估？原因是这需要在每个单元上设置一个新的正交点对象，该对象只有一个与你真正想要的点的参考坐标相对应的点；然后用它初始化一个FEValues对象；然后在这个点评估解决方案；然后将结果交给DataPostprocessor对象。这一连串的事情是昂贵的
 *
 * - 这就是 VectorTools::point_value() 昂贵的原因。在我们想要评估解决方案的每个单元格上使用相同的正交公式，并且只需要初始化一次FEValue对象，这是一个更便宜的选择，当然，代价是只能得到一个近似的结果。)   </ol>
 * 该类根据所提供的助记符自动为存储的数据生成名称。方法
 * @p add_component_names 和 @p
 * add_independent_names允许用户在需要时提供名称列表来代替使用。
 * 以下是一个小的代码片段，显示了这个类的一个常见用法。
 *
 *
 * @code
 * #include <deal.II/numerics/point_value_history.h>
 * //....
 *
 * //... code to setup Triangulation, perform any refinement necessary
 * // and setup DoFHandler, sizing solution Vectors etc
 *
 * // just one independent value, which happens to be an input
 * unsigned int n_inputs = 1;
 *
 * // call the constructor
 * PointValueHistory<dim> node_monitor(dof_handler, n_inputs);
 *
 * // setup fields and points required
 * node_monitor.add_field_name("Solution");
 * std::vector <Point <dim> > point_vector(2);
 * point_vector[0] = Point <dim>(0, 0);
 * point_vector[1] = Point <dim>(0.25, 0);
 * node_monitor.add_points(point_vector); // multiple points at once
 * node_monitor.add_point(Point<dim>(1, 0.2)); // add a single point
 * node_monitor.close(); // close the class once the setup is complete
 * node_monitor.status(std::cout); // print out status to check if desired
 *
 * // ... more code ...
 *
 * // ... in an iterative loop ...
 * // double time, vector <double> with size 1 input_value,
 * // and Vector <double> solution calculated in the loop
 * node_monitor.start_new_dataset(time);
 * node_monitor.push_back_independent(input_value);
 * node_monitor.evaluate_field("Solution", solution);
 *
 * // ... end of iterative loop ...
 *
 * node_monitor.write_gnuplot("node"); // write out data files
 *
 * @endcode
 *
 *
 *
 */
template <int dim>
class PointValueHistory
{
public:
  /**
   * 提供一个不支持添加点或网格数据的精简的类实例。
   * 这可用于记录外部输入或记录求解器性能数据。
   *
   */
  PointValueHistory(const unsigned int n_independent_variables = 0);

  /**
   * 构造函数将该类链接到一个特定的 @p DoFHandler. ，该类从
   * @p DoFHandler
   * 中读取特定的数据并存储在内部以便快速访问（特别是请求点的最近邻居的dof指数），如果需要支持点的数据，该类对
   * @p DoFHandler 的改变是相当不容忍的。网格细化和 @p
   * DoFRenumbering 方法应该在 @p add_points
   * 方法被调用之前进行，而且自适应网格细化只被一些方法支持。
   * 用户可以通过使用n_independent_variables指定所需的数量来存储与网格位置无关的额外变量，并根据需要调用
   * @p push_back_independent  。
   * 这可用于记录外部输入或记录求解器性能数据。
   *
   */
  PointValueHistory(const DoFHandler<dim> &dof_handler,
                    const unsigned int     n_independent_variables = 0);

  /**
   * 复制构造函数。这个构造函数可以安全地与包含数据的
   * @p
   * PointValueHistory对象一起调用，但这可能是昂贵的，应予以避免。
   *
   */
  PointValueHistory(const PointValueHistory &point_value_history);

  /**
   * 赋值运算符。一旦类被关闭并添加了数据，这个赋值运算符就可以被安全地调用，但这主要是为了允许在类中声明的
   * @p PointValueHistory
   * 对象在类的后期被重新初始化而提供的。当对象包含数据时，使用赋值运算符可能会很昂贵。
   *
   */
  PointValueHistory &
  operator=(const PointValueHistory &point_value_history);

  /**
   * 解构器。
   *
   */
  ~PointValueHistory();

  /**
   * 在类中添加一个点。在网格中找到最接近该点的支持点（每个组件一个），并存储其详细信息，以便在调用
   * @p evaluate_field
   * 时使用。如果需要一个以上的点，宁可使用 @p add_points
   * 方法，因为这样可以将网格的迭代次数降到最低。
   *
   */
  void
  add_point(const Point<dim> &location);

  /**
   * 添加多个点到类中。在网格中找到最接近该点的支持点（每个组件一个），并储存其细节，以便在调用
   * @p evaluate_field
   * 时使用。如果需要一个以上的点，宁可调用这个方法，因为它比add_point方法更有效率，因为它使网格的迭代次数降到最低。这些点会按照它们在列表中出现的顺序被添加到内部数据库中，即使一个点被多次请求，所请求的点和所添加的点之间总是存在一对一的对应关系。
   *
   */
  void
  add_points(const std::vector<Point<dim>> &locations);



  /**
   * 将另一个助记符串（也就是 @p VectorType) 放入类中。
   * 这个方法为变量增加的存储空间等于component_mask中真值的数量。这也为已经在类中的点增加了额外的条目，所以
   * @p add_field_name 和 @p add_points 可以按任何顺序调用。
   *
   */
  void
  add_field_name(const std::string &  vector_name,
                 const ComponentMask &component_mask = ComponentMask());

  /**
   * 将另一个助记符串（也就是 @p VectorType) 放入类中。
   * 这个方法增加了n_components变量的存储空间。这也为已经在类中的点增加了额外的条目，因此
   * @p add_field_name和 @p add_points
   * 可以按任何顺序调用。这个方法生成一个 std::vector 0,
   * ..., n_components-1并调用前面的函数。
   *
   */
  void
  add_field_name(const std::string &vector_name,
                 const unsigned int n_components);

  /**
   * 为一个字段的每个组件提供可选的名称。如果提供的话，这些名字将被用来代替从字段名中产生的名字。
   *
   */
  void
  add_component_names(const std::string &             vector_name,
                      const std::vector<std::string> &component_names);

  /**
   * 为独立值提供可选的名称。如果提供的话，这些名字将代替
   * "Indep_... "使用。
   *
   */
  void
  add_independent_names(const std::vector<std::string> &independent_names);



  /**
   * 从提供的VectorType中提取存储点的值，并将它们添加到新的数据集中的vector_name中。在添加字段时提供的成分掩码被用来选择要提取的成分。如果使用
   * @p
   * DoFHandler，必须为每个数据集（时间步长、迭代等）的每个vector_name调用一个（且只有一个）evaluate_field方法，否则会发生
   * @p ExcDataLostSync 错误。
   *
   */
  template <class VectorType>
  void
  evaluate_field(const std::string &name, const VectorType &solution);


  /**
   * 使用提供的 @p DataPostprocessor
   * 对象计算数值，并将其添加到vector_name中的新数据集。添加字段时提供的component_mask用于选择要从
   * @p DataPostprocessor
   * 返回矢量中提取的成分。该方法需要处理一个字段名的向量，如果许多字段使用相同的
   * @p DataPostprocessor
   * 对象，该方法是首选，因为每个单元只定位一次。提供的正交对象用于矢量场的所有分量。虽然这个方法在网格发生变化时不会抛出异常。由于每次调用该函数时都会重新选择正交点，因此没有内部数据结构被废止）。
   * 然而，用户必须意识到，如果网格发生变化，所选择的点也会略有不同，这使得数据的解释更加困难。如果使用
   * @p DoFHandler
   * ，必须为每个数据集（时间步长、迭代等）的每个向量_名称调用一个（且只有一个）evaluate_field方法，否则会发生
   * @p ExcDataLostSync 错误。
   *
   */
  template <class VectorType>
  void
  evaluate_field(const std::vector<std::string> &names,
                 const VectorType &              solution,
                 const DataPostprocessor<dim> &  data_postprocessor,
                 const Quadrature<dim> &         quadrature);

  /**
   * 构建一个只包含vector_name的 std::vector  <std::string>
   * 并调用上述函数。如果多个字段使用同一个 @p
   * DataPostprocessor 对象，上述函数会更有效率。
   *
   */
  template <class VectorType>
  void
  evaluate_field(const std::string &           name,
                 const VectorType &            solution,
                 const DataPostprocessor<dim> &data_postprocessor,
                 const Quadrature<dim> &       quadrature);


  /**
   * 从提供的VectorType中提取实际要求的点的值，并将其添加到vector_name中的新数据集。不像其他的evaluate_field方法，这个方法不关心dof_handler是否被修改，因为它使用调用
   * @p VectorTools::point_value
   * 来提取数据。因此，如果只使用这个方法，该类与自适应细化完全兼容。添加字段时提供的
   * component_mask 被用来选择要提取的组件。如果使用 @p
   * DoFHandler，必须为每个数据集（时间步长、迭代等）的每个vector_name调用一个（并且只有一个）evaluate_field方法，否则会发生
   * @p ExcDataLostSync 错误。
   *
   */
  template <class VectorType>
  void
  evaluate_field_at_requested_location(const std::string &name,
                                       const VectorType & solution);


  /**
   * 将当前数据集的密钥添加到数据集中。虽然先调用这个方法是明智的，但是这个方法、
   * @p evaluate_field和 @p push_back_independent
   * 的顺序并不重要。然而重要的是，一个给定的数据集的所有数据被添加到每个数据集中，并且在开始一个新的数据集之前添加。这可以防止
   * @p ExcDataLostSync. 的发生。
   *
   */
  void
  start_new_dataset(const double key);

  /**
   * 如果已经设置了独立的值，这个方法会存储这些值。
   * 每个数据集只能调用一次，如果使用了独立值，必须对每个数据集都调用。如果不调用这个方法，就会出现
   * @p ExcDataLostSync 的异常。
   *
   */
  void
  push_back_independent(const std::vector<double> &independent_values);


  /**
   * 写出一系列名为base_name + "-00.gpl"、base_name + "-01.gpl
   * "等的.gpl文件。数据文件给出了关于支持点选择的位置和解释数据的信息。如果
   * @p n_indep ！=0，一个额外的文件 base_name + "_indep.gpl"
   * 包含关键和独立数据。文件名称与点被添加到类中的顺序一致。数据列的名称可以通过函数
   * @p  add_component_names和 @p add_independent_names.  提供。
   * 支持点的信息只有在dof_handler没有改变的情况下才有意义。
   * 因此，如果使用了自适应网格细化，则不应使用支持点数据。可选的参数postprocessor_locations用于在输出文件中添加postprocessor位置。如果需要，应该在dof_handler可用的情况下，通过调用get_postprocessor_locations获得数据。默认参数是一个空的字符串向量，并将抑制后处理程序位置的输出。
   *
   */
  void
  write_gnuplot(const std::string &            base_name,
                const std::vector<Point<dim>> &postprocessor_locations =
                  std::vector<Point<dim>>());


  /**
   * 返回一个 @p Vector
   * ，其中选择的点的索引用1标记。这个方法主要是用来测试和验证该类的工作是否正常。通过将这个向量传递给DataOut对象，用户可以验证
   * @p get_support_locations 返回的位置是否与 @p DataOut 从 @p
   * Vector
   * 返回的位置一致。下面的代码片断演示了如何做到这一点。
   * @code
   * // Make a DataOut object and attach the dof_handler
   * DataOut<dim> data_out;
   * data_out.attach_dof_handler(dof_handler);
   *
   * // Call the mark_locations method to get the vector with indices flagged
   * Vector<double> support_point_locations = node_monitor.mark_locations();
   *
   * // Add the vector to the data_out object and
   * // write out a file in the usual way
   * data_out.add_data_vector(support_point_locations, "Monitor_Locations");
   * data_out.build_patches(2);
   * std::ofstream output("locations.gpl");
   * data_out.write_gnuplot(output);
   * @endcode
   *
   *
   */
  Vector<double>
  mark_support_locations();

  /**
   * 存储由 @p add_point(s)方法选择的每个支持点的实际位置。
   * 这可以用来与要求的点进行比较，例如通过使用 @p
   * Point<dim>::distance
   * 函数。为了方便起见，位置被调整为正确的点数，方法是。
   *
   */
  void
  get_support_locations(std::vector<std::vector<Point<dim>>> &locations);

  /**
   * 存储由data_postprocessor使用的点的实际位置。
   * 这可以用来与要求的点进行比较，例如通过使用 @p
   * Point<dim>::distance
   * 函数。与support_locations不同，这些位置是在每次用后处理器调用evaluate_field方法时计算的。这个方法使用相同的算法，所以可以找到相同的点。为方便起见，位置被该方法调整为正确的点数。
   *
   */
  void
  get_postprocessor_locations(const Quadrature<dim> &  quadrature,
                              std::vector<Point<dim>> &locations);

  /**
   * 一旦数据集被添加到类中，要求添加额外的点会使数据解释不清楚。布尔值
   * @p closed
   * 定义了该类的一个状态，确保这种情况不会发生。额外的点或向量只能在类未关闭时被添加，在数据集被添加或写入文件之前，类必须被关闭。
   * @p PointValueHistory::get_support_locations 和 @p
   * PointValueHistory::status
   * 并不要求类被关闭。如果一个需要类打开或关闭的方法在错误的状态下被调用，就会抛出
   * @p ExcInvalidState 的异常。
   *
   */
  void
  close();


  /**
   * 删除此对象对上次创建类时使用的 @p DoFHandler 的锁。
   * 这个方法通常不需要被调用，但是如果 @p PointValueHistory
   * 类的寿命可能比它长的话，这个方法对于确保 @p
   * DoFHandler
   * 在超出范围之前被释放是很有用的。一旦这个方法被调用，大多数方法都会抛出一个
   * @p ExcInvalidState
   * 的异常，所以如果使用这个方法应该是最后一次调用这个类。
   *
   */
  void
  clear();

  /**
   * 打印关于该类的有用的调试信息，包括为每个点选择了哪些支持点以及存储的数据大小的细节。
   *
   */
  void
  status(std::ostream &out);


  /**
   * 检查内部数据大小以测试数据同步的损失。这通常用于
   * @p Assert 语句中的 @p ExcDataLostSync 例外。  如果 @p strict 是
   * @p false
   * ，如果所有的大小都在1以内（需要允许添加数据），这个方法返回
   * @p true ，与 @p strict = @p true 它们必须是完全相等。
   *
   */

  bool
  deep_check(const bool strict);

  /**
   * 异常情况
   *
   */
  DeclExceptionMsg(ExcNoIndependent,
                   "A call has been made to push_back_independent() when "
                   "no independent values were requested.");

  /**
   * 异常情况
   *
   */
  DeclExceptionMsg(
    ExcDataLostSync,
    "This error is thrown to indicate that the data sets appear to be out of "
    "sync. The class requires that the number of dataset keys is the same as "
    "the number of independent values sets and mesh linked value sets. The "
    "number of each of these is allowed to differ by one to allow new values "
    "to be added with out restricting the order the user choses to do so. "
    "Special cases of no FHandler and no independent values should not "
    "trigger this error.");


  /**
   * 异常情况
   *
   */
  DeclExceptionMsg(
    ExcDoFHandlerRequired,
    "A method which requires access to a @p DoFHandler to be meaningful has "
    "been called when have_dof_handler is false (most likely due to default "
    "constructor being called). Only independent variables may be logged with "
    "no DoFHandler.");

  /**
   * 异常情况
   *
   */
  DeclExceptionMsg(
    ExcDoFHandlerChanged,
    "The triangulation has been refined or coarsened in some way. This "
    "suggests that the internal DoF indices stored by the current "
    "object are no longer meaningful.");

private:
  /**
   * 存储键，在标线上的值。这通常是时间，但也可能是时间步长，迭代等。
   *
   */
  std::vector<double> dataset_key;

  /**
   * 不依赖于网格位置的值。
   *
   */
  std::vector<std::vector<double>> independent_values;

  /**
   * 保存一个向量，列出与独立值相关的组件名称。如果用户不提供名称，这将是一个空的向量。
   *
   */
  std::vector<std::string> indep_names;

  /**
   * 为每个助记符条目保存数据。 data_store: 助记符
   *
   * -> [point_0_components point_1_components ... point_n-1_components][key]
   * 这种格式有利于标量记忆法在矢量空间中的应用，因为标量记忆法在每个点上只有一个成分。矢量成分严格来说是FE.n_components()长。
   *
   */
  std::map<std::string, std::vector<std::vector<double>>> data_store;

  /**
   * 为每个助记符保存一个分量掩码。
   *
   */
  std::map<std::string, ComponentMask> component_mask;


  /**
   * 保存一个向量，列出与助记符相关的元件名称。如果用户不提供名称，这将是一个空的矢量。
   *
   */
  std::map<std::string, std::vector<std::string>> component_names_map;

  /**
   * 保存支持点的位置和其他网格信息。
   *
   */
  std::vector<internal::PointValueHistoryImplementation::PointGeometryData<dim>>
    point_geometry_data;


  /**
   * 用于强制执行某些方法的 @p closed 状态。
   *
   */
  bool closed;

  /**
   * 用于为某些方法强制执行 @p !cleared 状态。
   *
   */
  bool cleared;


  /**
   * 一个指向提供给构造函数的dof_handler的智能指针。这可以通过调用
   * @p clear(). 来释放。
   *
   */
  SmartPointer<const DoFHandler<dim>, PointValueHistory<dim>> dof_handler;


  /**
   * 变量，用于检查三角测量是否已经改变。如果改变了，某些数据就会过时（尤其是
   * PointGeometryData::solution_indices.  *
   */
  bool triangulation_changed;



  /**
   * 一个布尔值，记录该类是否用DoFHandler初始化。
   *
   */
  bool have_dof_handler;

  /**
   * 用来检测来自三角区的信号。
   *
   */
  boost::signals2::connection tria_listener;

  /**
   * 存储要求的独立变量的数量。
   *
   */
  unsigned int n_indep;


  /**
   * 一个函数，每当三角剖面被修改时，它将通过信号被触发。
   * 目前，它用于检查三角形是否发生了变化，使预先计算的值无效。
   *
   */
  void
  tria_change_listener();
};


DEAL_II_NAMESPACE_CLOSE
#endif  /* dealii_point_value_history_h */ 


