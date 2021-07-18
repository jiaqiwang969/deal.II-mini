//include/deal.II-translator/numerics/vector_tools_0.txt
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

#ifndef dealii_vector_tools_h
#define dealii_vector_tools_h


#include <deal.II/base/config.h>

#include <deal.II/numerics/vector_tools_boundary.h>
#include <deal.II/numerics/vector_tools_common.h>
#include <deal.II/numerics/vector_tools_constraints.h>
#include <deal.II/numerics/vector_tools_evaluate.h>
#include <deal.II/numerics/vector_tools_integrate_difference.h>
#include <deal.II/numerics/vector_tools_interpolate.h>
#include <deal.II/numerics/vector_tools_mean_value.h>
#include <deal.II/numerics/vector_tools_point_gradient.h>
#include <deal.II/numerics/vector_tools_point_value.h>
#include <deal.II/numerics/vector_tools_project.h>
#include <deal.II/numerics/vector_tools_rhs.h>


DEAL_II_NAMESPACE_OPEN

// TODO: Move documentation of functions to the functions!

/**
 * 提供一个命名空间，提供对向量的一些操作。其中包括标准向量的组合，有限元解和连续函数之差的积分，连续函数对有限元空间的插值和投影以及其他操作。
 *
 *
 * @note
 * 几乎所有的函数都有两个版本，一个需要明确的Mapping参数，一个不需要。第二个版本通常用一个隐含的
 * $Q_1$
 * 参数来调用第一个版本（即用一个MappingQGeneric(1)类型的参数）。如果你打算让你的代码使用一个不同于（双/三）线性的映射，那么你需要调用应该使用<b>with</b>映射参数的函数。
 *
 *  <h3>Description of operations</h3>
 * 这个方法集合提供了以下操作。  <ul>   <li>  插值：将向量中的每个自由度分配为作为参数给出的函数的值。这等于说，所产生的有限元函数（与输出向量同构）在所有试验函数的支持点都有精确的函数值。试算函数的支持点是其值等于1的点，例如，对于线性试算函数，支持点是一个元素的四个角。因此，这个函数依赖于一个假设，即使用有限元，其自由度是函数值（Lagrange元），而不是梯度、法向导数、二次导数等（Hermite元、Quintic Argyris元等）。
 * 要创建的向量的某些值被设置两次甚至更多，这似乎是不可避免的。原因是我们必须在所有单元格上循环，并获得位于其上的每个试用函数的函数值。这也适用于位于面和角上的函数，因此我们要访问不止一次。虽然设置向量中的值不是一个昂贵的操作，但考虑到必须调用一个虚拟函数，对给定函数的评估可能是昂贵的。
 * <li>
 * 投影：计算给定函数对有限元空间的<i>L</i><sup>2</sup>投影，即如果<i>f</i>
 * ]是要投影的函数，计算<i>V<sub>h</sub></i>中的<i>f<sub>h</sub></i>，使(<i>f<sub>h</sub></i>,
 * <i>v<sub>h</sub></i>)=(<i>f</i>,
 * <i>v<sub>h</sub></i>)对所有离散测试函数<i>v<sub>h</sub></i>。这是通过线性方程组<i>
 * M v = f</i>的解来实现的，其中<i>M</i>是质量矩阵 $m_{ij} =
 * \int_\Omega \phi_i(x) \phi_j(x) dx$ 和 $f_i = \int_\Omega f(x) \phi_i(x)
 * dx$  。那么解向量 $v$
 * 就是投影<i>f<sub>h</sub></i>的节点表示。project()函数在
 * step-21  和  step-23  教程中使用。
 * 为了得到正确的结果，可能需要正确处理边界条件。下面列出了一些可能需要这样做的情况。
 * 如果需要，可以通过<i>L</i><sup>2</sup>将给定函数的迹线投影到局限于域的边界的有限元空间来完成，然后利用这些信息，用
 * MatrixTools::apply_boundary_values()
 * 函数从整个域的质量矩阵中消除边界节点。函数的轨迹向边界的投影是通过
 * VectorTools::project_boundary_values()
 * （见下文）函数完成的，该函数是用边界函数地图
 * std::map<types::boundary_id,  const
 * Function<spacedim,number>*>调用的，其中所有边界指标从零到
 * numbers::internal_face_boundary_id-1   (numbers::internal_face_boundary_id
 * 用于其他目的，见三角类文档）指向要投影的函数。对边界的投影是使用给project()函数的边界上的第二个正交公式来进行的。第一个正交公式被用来计算右手边和质量矩阵的数字正交。
 * 通常不需要先对边界值进行投影，然后从全局方程组中消除它们。如果你想对投影函数的边界值实施特殊的限制，例如在时间相关的问题中，这可能是必要的：你可能想投影初始值，但需要与后来的边界值保持一致。由于后者是在每个时间步骤中投射到边界上的，所以我们有必要在将初始值投射到整个域之前，也投射初始值的边界值。
 * 很明显，这两种投影方案的结果是不同的。通常，当先投射到边界时，原始函数和投射到整个域上的差值的<i>L</i><sup>2</sup>-norm会更大（已观察到5个因子），而在边界上积分的误差的<i>L</i><sup>2</sup>-norm当然应该更少。如果不进行对边界的投影，相反的情况也应该成立。
 * 选择是否需要先向边界投影是通过传递给函数的<tt>project_to_boundary_first</tt>标志完成的。
 * 如果给出 @p false ，则忽略面的额外正交公式。
 * 你应该注意的是，如果不要求向边界投影，一个边界值为零的函数在投影后可能不会有边界值为零。对于这种特别重要的情况有一个标志，它告诉函数在各自的边界部分强制执行零边界值。由于强制零边界值也可以通过投影达到，但使用其他方法获得更经济，如果设置了
 * @p enforce_zero_boundary 标志， @p project_to_boundary_first
 * 标志将被忽略。
 * 线性系统的求解目前是用简单的CG方法完成的，没有预处理，也没有多重网格。这显然不是太有效，但在许多情况下是足够的，而且容易实现。这个细节在未来可能会改变。
 * <li>  创建右手边的向量。create_right_hand_side()
 * 函数计算向量  $f_i = \int_\Omega f(x) \phi_i(x) dx$
 * 。这与取右手边的 <tt>MatrixCreator::create_*</tt>
 * 函数所做的相同，但不需要组装矩阵。
 * <li>  创建点源的右手边向量。create_point_source_vector()
 * 函数计算向量  $F_i = \int_\Omega \delta(x-x_0) \phi_i(x) dx$  。
 * <li>
 * 创建边界右手向量。create_boundary_right_hand_side()函数计算向量
 * $f_i = \int_{\partial\Omega} g(x) \phi_i(x) dx$
 * 。这是在拉普拉斯方程或其他二阶算子中具有不均匀的诺伊曼边界值时，边界力的右手边贡献。这个函数还需要一个可选的参数，表示积分应扩展到边界的哪些部分。如果使用默认参数，它将应用于所有边界。
 * <li>  边界值的插值。 MatrixTools::apply_boundary_values()
 * 函数接收一个边界节点及其数值的列表。你可以通过使用interpolate_boundary_values()函数对边界函数进行内插来得到这样一个列表。要使用它，你必须指定一对边界指标（类型为
 * <tt>types::boundary_id</tt>;
 * 详见三角形类文档中的章节）和表示具有该边界指标的边界面上的节点的Dirichlet边界值的相应函数列表。
 * 通常，所有其他的边界条件，如不均匀的Neumann值或混合边界条件都在弱式计算中处理。因此，没有尝试将这些纳入矩阵和矢量的组装过程中。
 * 在这个函数中，边界值是插值的，也就是说，一个节点被赋予了边界函数的点值。在某些情况下，可能需要使用边界函数的L2投影或任何其他方法。为此，我们参考了下面的project_boundary_values()函数。
 * 你应该知道，边界函数可以在面的内部的节点上进行评估。然而，这些节点不需要在真正的边界上，而是在由单元格到实际单元格的映射所代表的边界的近似上。由于这种映射在大多数情况下不是面的精确映射，边界函数是在不在边界上的点上评估的，你应该确保返回的值在某种意义上是合理的。
 * 在1d中，情况有点不同，因为那里的面（即顶点）没有边界指标。我们假定，如果在边界函数列表中给出了边界指标0，那么左边的边界点将被插值，而右边的边界点则与地图中的边界指标1相关。然后，各自的边界函数在各自的边界点的地方被评估。
 * <li>
 * 边界值的投影。project_boundary_values()函数的作用类似于interpolate_boundary_values()函数，除了它不是通过内插得到边界节点的节点值，而是通过<i>L</i><sup>2</sup>函数的轨迹投影到边界上。
 * 投射发生在所有具有边界指标的边界部分，这些边界指标列在边界函数的地图
 * (std::map<types::boundary_id, const
 * Function<spacedim,number>*>中。这些边界部分可能是连续的，也可能不是。对于这些边界部分，使用
 * MatrixTools::create_boundary_mass_matrix()
 * 函数组装质量矩阵，以及适当的右手边。然后用一个简单的CG方法（没有预设条件）来解决所产生的方程组，在大多数情况下，这对目前的目的是足够的。
 * <li>  计算错误。函数 integrate_difference()
 * 可以计算给定的（连续）参考函数和不同规范的有限元解之间的误差。积分使用给定的正交公式进行，并假定给定的有限元对象等于用于计算解决方案的对象。
 * 结果被存储在一个向量中（命名为 @p difference),
 * ，其中每个条目等于一个单元上的差值的给定规范。当用
 * @p begin_active
 * 开始并用<tt>++</tt>操作符推广时，条目的顺序与 @p
 * cell_iterator 相同。
 * 这个数据，每个活动单元一个数字，可以通过
 * DataOut::add_data_vector
 * 函数直接传递给DataOut类来产生图形输出。最后，可以使用
 * VectorTools::compute_global_error(). 函数将
 * VectorTools::integrate_difference()
 * 中每个单元的输出插值到有限元场的结点上。
 * 目前，有可能从每个单元的差值中计算出以下数值。  @p
 * mean,   @p L1_norm,   @p L2_norm,   @p  Linfty_norm,  @p H1_seminorm  和
 * @p H1_norm,  见  VectorTools::NormType.
 * 对于平均差值，计算参考函数减去数值解，而不是反过来。
 * 某一单元格上的差值的无穷大准则返回正交公式参数给出的正交点上的差值的最大绝对值。这在某些情况下并不是一个很好的近似值，因为例如高斯正交公式并不评估单元格的端点或角点的差。在这种情况下，你可能想选择一个有更多正交点的正交公式，或者选择一个有另一种正交点分布的正交公式。你还应该考虑到有限元在某些点上的超融合特性：例如在一维中，标准的有限元方法是一种拼合方法，应该在节点点上返回精确值。因此，梯形规则应该总是返回一个消失的L-无穷大误差。相反，在二维中，最大的L-无穷大误差应该位于顶点或单元的中心，这使得使用辛普森正交法则是合理的。另一方面，在高斯积分点可能会出现超融合。这些例子不是作为经验法则，而是想说明，使用错误的正交公式可能会显示出明显错误的结果，应该注意选择正确的公式。
 * <i>H</i><sup>1</sup>半正则是<i>L</i><sup>2</sup>差值梯度的准则。全<i>H</i><sup>1</sup>规范的平方是seminorm的平方与<i>L</i><sup>2</sup>规范的平方之和。
 * 为了得到全局<i>L<sup>1</sup></i>的误差，你必须对 @p
 * difference, 中的条目进行求和，例如使用 Vector::l1_norm()
 * 函数。
 * 对于全局<i>L</i><sup>2</sup>差值，你必须将各条目的平方相加，然后取和的根，例如使用
 * Vector::l2_norm().
 * 这两个操作代表向量的<i>l</i><sub>1</sub>和<i>l</i><sub>2</sub>规范，但你不需要取每个条目的绝对值，因为单元格规范已经是正的了。
 * 要得到全局平均差，只需像上面那样将各元素相加。要得到
 * $L_\infty$ 准则，取向量元素的最大值，例如，使用
 * Vector::linfty_norm() 函数。
 * 对于全局<i>H</i><sup>1</sup>规范和半规范，与<i>L</i><sup>2</sup>规范的规则相同：计算细胞误差向量的<i>l</i><sub>2</sub>规范。
 * 注意，在一维的情况下，如果你要求一个需要计算梯度的规范，那么所提供的函数会自动沿曲线投影，并且只计算梯度的切向部分的差值，因为无论如何都没有关于梯度的法向部分的信息。  </ul>
 * 所有函数都使用上次自由度分布在三角形上时给DoFHandler对象的有限元。另外，如果需要访问描述边界确切形式的对象，则访问存储在三角剖分对象中的指针。
 *
 *
 * @note
 * 这个模板的实例化提供给一些矢量类型，特别是<code>Vector&lt;float&gt;,
 * Vector&lt;double&gt;, BlockVector&lt;float&gt;,
 * BlockVector&lt;double&gt;</code>；其他的可以在应用代码中生成（见手册中
 * @ref Instantiations 部分）。
 *
 *
 * @ingroup numerics
 *
 *
 */
namespace VectorTools
{}

DEAL_II_NAMESPACE_CLOSE

#endif


