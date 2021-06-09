//include/deal.II-translator/dofs/dof_handler_0.txt
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

#ifndef dealii_dof_handler_h
#define dealii_dof_handler_h



#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/function.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/iterator_range.h>
#include <deal.II/base/smartpointer.h>

#include <deal.II/distributed/tria_base.h>

#include <deal.II/dofs/block_info.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_faces.h>
#include <deal.II/dofs/dof_iterator_selector.h>
#include <deal.II/dofs/dof_levels.h>
#include <deal.II/dofs/number_cache.h>

#include <deal.II/hp/fe_collection.h>

#include <boost/serialization/split_member.hpp>

#include <map>
#include <memory>
#include <set>
#include <vector>

DEAL_II_NAMESPACE_OPEN

// Forward declarations
#ifndef DOXYGEN
template <int dim, int spacedim>
class FiniteElement;
template <int dim, int spacedim>
class Triangulation;

namespace internal
{
  namespace DoFHandlerImplementation
  {
    struct Implementation;

    namespace Policy
    {
      template <int dim, int spacedim>
      class PolicyBase;
      struct Implementation;
    } // namespace Policy
  }   // namespace DoFHandlerImplementation

  namespace DoFAccessorImplementation
  {
    struct Implementation;
  }

  namespace DoFCellAccessorImplementation
  {
    struct Implementation;
  }

  namespace hp
  {
    namespace DoFHandlerImplementation
    {
      struct Implementation;
    }
  } // namespace hp
} // namespace internal

namespace parallel
{
  namespace distributed
  {
    template <int dim, int spacedim, typename VectorType>
    class CellDataTransfer;
  }
} // namespace parallel
#endif

/**
 * 给定一个三角形和一个有限元的描述，这个类列举了三角形的所有顶点、边、面和单元的自由度。因此，它也为离散空间 $V_h$ 提供了一个<i>basis</i>，其元素是由FiniteElement对象在每个单元上定义的有限元函数。这个类满足了 @ref ConceptMeshType 的 "MeshType概念 "
 * 要求。 它首次在 step-2 教程程序中使用。
 * 对于每个顶点、直线、四边形等，这个类存储了一个生活在这个对象上的自由度指数的列表。这些指数指的是无约束的自由度，也就是说，有约束的自由度的编号与无约束的自由度的编号是一样的，只是后来被淘汰了。
 * 这导致全局向量和矩阵中的指数也是指所有的自由度，需要进行某种浓缩来限制方程组只包括无约束的自由度。指数的实际存储布局在
 * dealii::internal::DoFHandlerImplementation::DoFLevel 类文档中描述。
 * 该类提供了遍历所有单元的迭代器，与Triangulation类的方式基本相同。使用begin()和end()函数（和同伴，如begin_active()），可以获得遍历单元的迭代器，并查询自由度结构以及三角形的数据。这些迭代器建立在Triangulation类的迭代器之上，但与纯粹的三角形迭代器相比，提供了额外的自由度功能信息。自由度迭代器由<tt>++</tt>和<tt>--</tt>操作符呈现的顺序与遍历三角形的相应迭代器相同，该DoFHandler是在三角形上构建的。
 * 如果要解决曲面上的问题，必须使用<tt>spacedim</tt>参数。如果没有指定，这个参数的默认值是<tt>=dim</tt>，意味着我们要在一个维度等于它所嵌入的空间维度的域中解决问题。
 *
 *  <h3>Distribution of indices for degrees of freedom</h3>
 * 自由度（`dofs'）通过函数distribut_dofs()在给定的三角形上分布。它被传递给一个有限元对象，描述有多少自由度位于顶点、线等处。它逐个遍历三角形，如果尚未编号，则对该单元的自由度进行编号。对于非多重网格算法，只考虑活动单元。活动单元被定义为那些没有子单元的单元，也就是说，它们是最精细的单元。
 * 由于三角形的遍历是从最粗的活动层的单元开始，然后到更细的层次，因此最大的单元以及它们的边界线和顶点的系数最低，更细的单元的系数更高。
 * 这种编号方式意味着所得矩阵的带宽非常大，因此对于某些求解算法来说是非常不理想的。由于这个原因，DoFRenumbering类提供了几种算法来重新安排dof的编号。关于已实现的算法的讨论见这里。
 *
 *  <h3>Interaction with distributed meshes</h3>
 * 在构造时，这个类需要一个对三角形对象的引用。在大多数情况下，这将是对Triangulation类型对象的引用，即代表完全驻留在单个处理器上的三角形的类。然而，它也可以是 parallel::distributed::Triangulation
 * 类型（例如，见 step-32 、 step-40 ，特别是 @ref distributed
 * 模块），在这种情况下，DoFHandler对象将只管理本地拥有的和幽灵单元的自由度。这个过程对使用者来说是完全透明的。
 *
 *  <h3>User defined renumbering schemes</h3>
 * DoFRenumbering类提供了许多重新编号的方案，如Cuthill-McKee方案。基本上，该函数设置了一个数组，在这个数组中，我们为每个自由度存储了重新编号后这个自由度应该有的新索引。使用这个数组，本类的renumber_dofs()函数被调用，它实际上执行了从旧的DoF指数到数组中给出的指数的改变。然而，在某些情况下，用户可能想计算他们自己的重新编号顺序；在这种情况下，可以分配一个数组，每个自由度有一个元素，并在其中填入各自自由度应被分配的编号。例如，这个数字可以通过对下风方向的自由度的支持点进行排序来获得。
 * 然后用数组调用
 * <tt>renumber_dofs(vector<types::global_dof_index>)</tt>
 * 函数，将旧的自由度指数转换为新的自由度指数。
 *
 *  <h3>Serializing (loading or storing) DoFHandler objects</h3>
 * 像deal.II中的许多其他类一样，DoFHandler类可以使用BOOST的序列化设施将其内容流向一个存档。这样存储的数据以后可以再次从存档中检索，以恢复这个对象的内容。这个工具经常被用来将程序的状态保存到磁盘上，以便以后可能复活，通常是在长期运行的计算的检查点/重启策略的背景下，或者在不是很可靠的计算机上（例如，在非常大的集群上，个别节点偶尔会出现故障，然后导致整个MPI作业的失败）。
 * DoFHandler类这样做的模式与Triangulation类相似（见该类的一般文档中的章节）。特别是，load()函数并不完全恢复与之前使用save()函数存储的相同状态。相反，该函数假设你将数据加载到一个DoFHandler对象中，该对象已经与一个三角测量相关联，其内容与保存数据时使用的三角测量一致。同样，load()函数假定当前对象已经与一个有限元对象相关联，该对象与数据保存时的对象相匹配；后者可以通过在从序列化存档中重新加载数据之前使用同种有限元调用
 * DoFHandler::distribute_dofs() 实现。
 *
 *  <h3>hp-adaptive finite element methods</h3>
 * 该类不允许在所有单元上只使用一个特定的有限元，而是允许在每个单元上列举不同有限元的自由度。为此，我们给每个单元格分配一个
 * <code>active_fe_index</code>
 * ，表示在有限元集合中的哪个元素（由 hp::FECollection)
 * 类型的对象代表）是住在这个单元格上的。然后，该类列举了三角形的每个单元上与这些有限元相关的自由度，如果可能的话，还可以识别单元界面上的自由度，如果它们匹配的话。如果相邻的单元沿共同界面的自由度不立即匹配（例如，如果你有
 * $Q_2$ 和 $Q_3$
 * 元素在一个共同的面相遇），那么就需要计算约束，以确保网格上产生的有限元空间保持一致。
 * 使用这种类型的对象的整个过程在  step-27  中解释。这个类实现的许多算法在 @ref hp_paper  "hp-paper "
 * 中描述。
 *
 *  <h3>Active FE indices and their behavior under mesh refinement</h3>
 * 使用这个类的典型工作流程是创建一个网格，给每个活动单元分配一个活动FE索引，调用
 * DoFHandler::distribute_dofs(),
 * ，然后在这个有限元空间上组装一个线性系统并解决问题。
 * 活跃FE指数将在网格适应过程中自动从旧网格转移到新网格。未来的FE指数是用来确定网格适应后的活动FE指数，并用于为新网格准备旧网格上的数据。如果没有指定未来FE指数，则以有限元为准。
 * 特别是在适应过程中，以下规则适用。
 *
 *
 *
 * - 网格细化时，子单元会继承父单元的未来FE指数。
 *
 *
 *
 * - 当粗化单元时，父单元（现在处于活动状态）将被分配一个未来的FE指数，该指数由其（不再活动的）子单元决定，遵循FiniteElementDomination逻辑。在以前分配给前子女的元素集合中，我们选择一个由所有子女支配的元素作为父单元。如果没有找到，我们就在整个集合中挑选一个被所有前子女支配的最主要的元素。关于这个主题的进一步信息，请参见 hp::FECollection::find_dominated_fe_extended() 。
 * 在 hp::Refinement
 * 命名空间中，有自动hp-适应的策略，它将根据标准设置未来的FE指数。
 *
 *  <h3>Active FE indices and parallel meshes</h3> 当该类与
 * parallel::shared::Triangulation  或  parallel::distributed::Triangulation,
 * 一起使用时，你只能通过
 * <code>cell-@>set_active_fe_index(...)</code>
 * 这样的调用，在本地拥有的单元上设置活动FE指数。另一方面，不允许在幽灵或人造单元上设置活动FE指数。
 * 然而，幽灵单元确实获得了什么元素在其上处于活动状态的信息：每当你调用
 * DoFHandler::distribute_dofs(),
 * 时，所有参与并行网格的处理器都会以这样的方式交换信息，幽灵单元上的活动FE指数等于在拥有该特定幽灵单元的处理器上设置的活动FE指数。因此，人们可以在幽灵单元上<i>query</i>
 * @p active_fe_index ，只是不能用手去设置它。
 * 在人工单元上，没有关于那里使用的 @p active_fe_index 的信息。这是因为我们甚至不知道这些细胞是否存在，即使存在，目前的处理器也不知道关于它们的任何具体信息。更多信息见 @ref GlossArtificialCell "人工细胞的词汇表条目"
 * 。 在细化和粗化过程中，关于每个单元的 @p active_fe_index
 * 的信息将被自动转移。
 * 然而，在hp模式下使用带有DoFHandler的
 * parallel::distributed::Triangulation
 * ，在序列化过程中需要额外注意，因为没有关于活动FE指数的信息会被自动传输。这必须使用prepare_for_serialization_of_active_fe_indices()和deserialize_active_fe_indices()函数手动完成。前者必须在调用
 * parallel::distributed::Triangulation::save() 之前调用，后者需要在
 * parallel::distributed::Triangulation::load(). 之后运行。
 * 如果进一步的数据将通过 parallel::distributed::CellDataTransfer,
 * parallel::distributed::SolutionTransfer, 或 Particles::ParticleHandler
 * 类附加到三角形上，所有相应的准备和反序列化函数调用需要以相同顺序发生。更多信息请参考
 * parallel::distributed::SolutionTransfer 的文档。
 *
 *
 * @ingroup dofs
 *
 */
template <int dim, int spacedim = dim>
class DoFHandler : public Subscriptor
{
  using ActiveSelector =
    dealii::internal::DoFHandlerImplementation::Iterators<dim, spacedim, false>;
  using LevelSelector =
    dealii::internal::DoFHandlerImplementation::Iterators<dim, spacedim, true>;

public:
  /**
   * 一个别名，用于识别DoFHandler对象中的单元格迭代器。  迭代器的概念在 @ref Iterators "迭代器文档模块 "
   * 中有详细的讨论。    目前的别名在本质上与相应的
   * Triangulation::cell_accessor
   * 别名一样工作。然而，除了已经通过CellAccessor类提供的成员函数外，它还提供了DoFCellAccessor的成员函数。
   * @ingroup Iterators
   *
   */
  using cell_accessor = typename ActiveSelector::CellAccessor;

  /**
   * 一个别名，用于识别指向面的迭代器。  迭代器的概念在 @ref Iterators "迭代器文档模块 "
   * 中有详细的讨论。    当前的别名在本质上与相应的
   * Triangulation::face_accessor
   * 别名一样工作。然而，除了已经通过TriaAccessor类提供的成员函数外，它还提供了DoFAccessor的成员函数。
   * @ingroup Iterators
   *
   */
  using face_accessor = typename ActiveSelector::FaceAccessor;

  /**
   * 一个别名，定义了一个网格的（一维）线条的迭代器。在一维网格中，这些线是网格的单元，而在二维网格中，这些线是单元的面。
   * @ingroup Iterators
   *
   */
  using line_iterator = typename ActiveSelector::line_iterator;

  /**
   * 一个别名，允许在<i>active</i>线上迭代，即没有子节点的线的子集。在一维网格中，这些是网格的单元，而在二维网格中，线是单元的面。
   * 在二维或三维网格中，没有子节点的线（即活动线）是至少一个活动单元的一部分。每条这样的线还可能是与活动单元相邻的更粗的单元的线的子线。(这个较粗的邻居也是活动的)。
   * @ingroup Iterators
   *
   */
  using active_line_iterator = typename ActiveSelector::active_line_iterator;

  /**
   * 一个别名，定义了一个网格的（二维）四边形的迭代器。在二维网格中，这些是网格的单元，而在三维网格中，四边形是单元的面。
   * @ingroup Iterators
   *
   */
  using quad_iterator = typename ActiveSelector::quad_iterator;

  /**
   * 一个允许在<i>active</i>四边形上迭代的别名，即没有子节点的四边形子集。在二维网格中，这些是网格的单元，而在三维网格中，四边形是单元的面。
   * 在三维网格中，没有孩子的四边形（即活动四边形）是至少一个活动单元的面。每个这样的四边形还可能是与活动单元相邻的更粗的单元的四边形面的子。这个较粗的邻居也将是活动的）。
   * @ingroup Iterators
   *
   */
  using active_quad_iterator = typename ActiveSelector::active_quad_iterator;

  /**
   * 一个别名，定义了一个网格的（三维）六边形的迭代器。这个迭代器只有在三维网格中才有意义，在三维网格中六边形是网格的单元。
   * @ingroup Iterators
   *
   */
  using hex_iterator = typename ActiveSelector::hex_iterator;

  /**
   * 一个别名，允许遍历网格的<i>active</i>个六边形。
   * 这个迭代器只有在三维网格中才有意义，在三维网格中六边形是网格的单元。因此，在这些三维网格中，这个迭代器等同于
   * @p active_cell_iterator  别名。
   * @ingroup Iterators
   *
   */
  using active_hex_iterator = typename ActiveSelector::active_hex_iterator;

  /**
   * @ingroup Iterators
   *
   */
  using active_cell_iterator = typename ActiveSelector::active_cell_iterator;

  /**
   * @ingroup Iterators
   *
   */
  using cell_iterator = typename ActiveSelector::cell_iterator;

  /**
   * 一个别名，用于识别指向面的迭代器。  迭代器的概念在 @ref Iterators "迭代器文档模块 "
   * 中有详细的讨论。
   * 虽然这个别名的实际数据类型隐藏在几层（不幸的是必须的）间接因素后面，但它本质上是TriaIterator<DoFAccessor>。TriaIterator类的工作方式就像一个指向对象的指针，当你解除引用时，会产生一个DoFAccessor类型的对象。DoFAccessor又是一个可以用来查询面的DoF索引的类，但它也是从TriaAccessor派生出来的，因此可以用来查询几何属性，如面的顶点、面积等。
   * @ingroup Iterators
   *
   */
  using face_iterator = typename ActiveSelector::face_iterator;

  /**
   * 一个别名，用于识别指向活动面的迭代器，即指向没有子节点的面。活动面必须是至少一个活动单元的面。
   * 除了 "活动 "的限定，这个别名与 @p face_iterator
   * 的别名相同。特别是，取消引用这两个别名都会产生相同的对象。
   * @ingroup Iterators
   *
   */
  using active_face_iterator = typename ActiveSelector::active_face_iterator;

  using level_cell_accessor = typename LevelSelector::CellAccessor;
  using level_face_accessor = typename LevelSelector::FaceAccessor;

  using level_cell_iterator = typename LevelSelector::cell_iterator;
  using level_face_iterator = typename LevelSelector::face_iterator;


  /**
   * 让维度在函数模板中可用。
   *
   */
  static const unsigned int dimension = dim;

  /**
   * 使空间维度在函数模板中可用。
   *
   */
  static const unsigned int space_dimension = spacedim;

  /**
   * 在给定单元上使用的有限元的默认索引。
   *
   */
  static const unsigned int default_fe_index = 0;

  /**
   * 在给定单元上使用的有限元的无效索引。
   *
   */
  static const unsigned int invalid_fe_index = numbers::invalid_unsigned_int;

  /**
   * 我们存储活动FE索引的类型。
   *
   */
  using active_fe_index_type = unsigned short int;

  /**
   * 我们在CRS数据结构中存储偏移量的类型。
   *
   */
  using offset_type = unsigned int;

  /**
   * 无效的活动FE索引，它将被用作默认值，以确定未来的FE索引是否已经被设置。
   *
   */
  static const active_fe_index_type invalid_active_fe_index =
    static_cast<active_fe_index_type>(-1);

  /**
   * 标准构造函数，不初始化任何数据。用这个构造函数构造一个对象后，使用
   * reinit() 得到一个有效的 DoFHandler。
   *
   */
  DoFHandler();

  /**
   * 构造函数。以 @p tria 作为工作中的三角关系。
   *
   */
  explicit DoFHandler(const Triangulation<dim, spacedim> &tria);

  /**
   * 复制构造函数。DoFHandler对象很大，很昂贵。
   * 它们不应该被复制，特别是不应该被意外地复制，而应该被有意地构造。因此，这个构造函数被明确地从这个类的接口中移除。
   *
   */
  DoFHandler(const DoFHandler &) = delete;

  /**
   * 解构器。
   *
   */
  virtual ~DoFHandler() override;

  /**
   * 复制操作符。DoFHandler对象是大的和昂贵的。
   * 它们不应该被复制，特别是不应该被意外地复制，而应该被有意地构建。因此，这个操作符被明确地从这个类的接口中移除。
   *
   */
  DoFHandler &
  operator=(const DoFHandler &) = delete;

  /**
   * 给DoFHandler分配一个三角形和一个有限元素，并计算网格上的自由度分布。
   * @deprecated  使用 reinit() 和 distribute_dofs() 来代替。
   *
   */
  DEAL_II_DEPRECATED
  void
  initialize(const Triangulation<dim, spacedim> &tria,
             const FiniteElement<dim, spacedim> &fe);

  /**
   * 与上述相同，但需要一个 hp::FECollection 对象。
   * @deprecated  使用 reinit() 和 distribute_dofs() 来代替。
   *
   */
  DEAL_II_DEPRECATED
  void
  initialize(const Triangulation<dim, spacedim> &   tria,
             const hp::FECollection<dim, spacedim> &fe);

  /**
   * 给这个对象分配一个FiniteElement  @p fe  。
   * @note
   * 这个函数对作为参数的有限元进行复制，并将其作为成员变量存储。因此，可以编写如下代码
   * @code
   * dof_handler.set_fe(FE_Q<dim>(2));
   * @endcode
   * 然后你可以通过调用 DoFHandler::get_fe().
   * 来访问有限元。然而，通常更方便的做法是将命名的有限元对象作为成员变量保存在你的主类中，当你需要访问有限元的属性时直接引用它（如
   * FiniteElementData::dofs_per_cell).
   * 这就是所有教程程序的做法。      @warning
   * 这个函数只设置一个FiniteElement。自由度要么还没有被分配，要么是使用先前设置的元素进行分配。在这两种情况下，访问自由度将导致无效的结果。为了恢复一致性，请调用distribution_dofs()。
   * @deprecated  用distribution_dofs()代替。
   *
   */
  DEAL_II_DEPRECATED
  void
  set_fe(const FiniteElement<dim, spacedim> &fe);

  /**
   * 与上述相同，但需要一个 hp::FECollection 对象。
   * @deprecated  使用 distribute_dofs() 来代替。
   *
   */
  DEAL_II_DEPRECATED
  void
  set_fe(const hp::FECollection<dim, spacedim> &fe);

  /**
   * 穿过三角形，将所有活动单元的活动FE指数设置为 @p
   * active_fe_indices. 中给出的值。
   *
   */
  void
  set_active_fe_indices(const std::vector<unsigned int> &active_fe_indices);

  /**
   * 通过三角剖分，将所有活动单元的活动FE指数存储到向量
   * @p active_fe_indices. 中，如有必要，该向量将被调整大小。
   *
   */
  void
  get_active_fe_indices(std::vector<unsigned int> &active_fe_indices) const;

  /**
   * 给DoFHandler分配一个三角图。
   * 移除与之前的Triangulation对象的所有关联，并与新对象建立连接。所有关于以前的自由度的信息将被删除。激活hp-mode。
   *
   */
  void
  reinit(const Triangulation<dim, spacedim> &tria);

  /**
   * 通过三角剖分并 "分配
   * "给定有限元所需的自由度。"分布
   * "自由度包括分配内存来存储所有自由度的实体（如顶点、边、面等）的索引，然后列举所有自由度。换句话说，虽然网格和有限元对象本身只是定义了一个有限元空间
   * $V_h$
   * ，但分配自由度的过程确保了这个空间有一个基础，并且这个基础的形状函数是以可索引、可预测的方式列举的。
   * 网格上自由度的确切顺序，即有限元空间的基函数被列举的顺序，是deal.II作为一个实现细节处理的东西。一般来说，自由度的列举顺序与我们遍历单元的顺序相同，但你不应该依赖任何特定的编号。相反，如果你想要一个特定的顺序，可以使用命名空间DoFRenumbering中的函数。
   * 这个函数在 step-2 教程程序的介绍中首次讨论。
   * @note
   * 该函数对作为参数给定的有限元进行复制，并将其存储为成员变量，与上述函数set_fe()类似。
   *
   */
  void
  distribute_dofs(const FiniteElement<dim, spacedim> &fe);

  /**
   * 同上，但取一个 hp::FECollection 对象。
   *
   */
  void
  distribute_dofs(const hp::FECollection<dim, spacedim> &fe);

  /**
   * 为几何多网格在每个层面上分配层面自由度。在调用此函数之前，需要使用distribut_dofs()来分配活动的自由度。
   *
   */
  void
  distribute_mg_dofs();

  /**
   * 返回这个DoFHandler是否有hp能力。
   *
   */
  bool
  has_hp_capabilities() const;

  /**
   * 这个函数返回这个DoFHandler是否有分布在每个多网格层次上的DoF，或者换句话说，如果distribution_mg_dofs()已经被调用。
   *
   */
  bool
  has_level_dofs() const;

  /**
   * 这个函数返回这个DoFHandler是否有活动的DoF。这相当于询问(i)distribution_dofs()是否被调用，以及(ii)被分配自由度的有限元是否真的有自由度（例如FE_Nothing就不是这种情况）。
   * 如果这个对象是基于 parallel::distributed::Triangulation,
   * ，那么如果平行DoFHandler对象的<i>any</i>分区有任何自由度，当前函数返回true。换句话说，即使Triangulation在当前MPI进程上不拥有任何活动单元，但至少有一个进程拥有单元，而且至少这一个进程有任何自由度与之相关，该函数也会返回true。
   *
   */
  bool
  has_active_dofs() const;

  /**
   * 在用FESystem元素distribution_dofs()之后，全局和水平向量的块结构被存储在BlockInfo对象中，可以用block_info()访问。这个函数在同一对象的每个单元上初始化本地块结构。
   *
   */
  void
  initialize_local_block_info();

  /**
   * 清除此对象的所有数据。
   *
   */
  void
  clear();

  /**
   * 根据每个自由度的新DoF指数列表对自由度重新编号。
   * 这个函数是由DoFRenumbering函数中的函数在计算自由度指数的新排序后调用的。然而，它当然也可以从用户代码中调用。
   * @arg  new_number
   * 这个数组的大小必须等于当前处理器所拥有的自由度数量，也就是说，大小必须等于n_locally_owned_dofs()的返回值。如果只有一个处理器参与存储当前的网格，那么就等于总的自由度数，即n_dofs()的结果。这个数组的内容是由local_owned_dofs()返回的IndexSet中所列的每个自由度的新全局索引。在顺序网格的情况下，这意味着这个数组是当前网格上每个自由度的新索引列表。如果我们有一个
   * parallel::shared::Triangulation 或 parallel::distributed::Triangulation
   * 的底层DoFHandler对象，该数组是所有本地拥有的自由度的新索引列表，列举的顺序与当前本地拥有的DoF相同。换句话说，假设自由度
   * <code>i</code> 是当前本地拥有的，那么
   * <code>new_numbers[locally_owned_dofs().index_within_set(i)]</code>
   * 返回新的全局自由度索引 <code>i</code>
   * 。由于在顺序的情况下，local_owned_dofs()的IndexSet是完整的，在只有一个处理器参与网格的情况下，后一种对数组内容的约定可以简化为前一种。
   * @note  虽然从上面可以看出，知道并行计算中本地拥有的自由度的<i>number</i>是重新编号下的不变量可能会令人惊讶，即使与这些本地拥有的自由度相关的<i>indices</i>不是。从根本上说，这个不变性的存在是因为<i>decision</i>一个自由度是否为本地所有，与该自由度的（新旧）索引无关。事实上，如果自由度在本地拥有的单元上，而不是在相邻单元具有较低 @ref GlossSubdomainId "子域id "
   * 的单元之间的界面上，那么自由度就是本地拥有的。由于这两个条件都与与DoF相关的索引无关，一个本地拥有的自由度在重新编号后也将是本地拥有的。
   * 另一方面，诸如本地拥有的DoF的索引集是否形成一个连续的范围（即Local_owned_dofs()是否返回
   * IndexSet::is_contiguous() 返回 @p true)
   * 的IndexSet对象）等属性当然会受到这里进行的确切重新编号的影响。例如，虽然在
   * distribute_dofs() 中对 DoF
   * 指数进行的初始编号产生了一个连续的编号，但由
   * DoFRenumbering::component_wise()
   * 执行的重新编号通常不会产生连续的本地所有 DoF
   * 指数。
   *
   */
  void
  renumber_dofs(const std::vector<types::global_dof_index> &new_numbers);

  /**
   * 与上述功能相同，但对多网格层次结构中的单层自由度进行重新编号。
   *
   */
  void
  renumber_dofs(const unsigned int                          level,
                const std::vector<types::global_dof_index> &new_numbers);

  /**
   * 返回给定三角结构中一个自由度与给定有限元可能耦合的最大自由度数。
   * 这是系统矩阵中每行的最大条目数；因此这一信息可以在构建SparsityPattern对象时使用。
   * 返回的数字并不是真正的最大数字，而是基于有限元和在一个顶点相遇的最大单元数的估计。这个数字对受限矩阵也是有效的。
   * 耦合数的确定可以通过简单的画图来完成。在这个函数的实现中可以找到一个例子。
   * @note
   * 这个函数最常被用来确定稀疏模式的最大行长度。不幸的是，虽然这个函数返回的估计值在1d和2d中相当准确，但在3d中往往明显过高，导致SparsityPattern类在某些情况下分配太多的内存。除非有人能够改进目前的三维函数，否则对于这些情况，我们没有什么办法。解决这个问题的典型方法是使用一个中间压缩的稀疏模式，只在需要时分配内存。请参考
   * step-2 和 step-11
   * 的例子程序，了解如何做到这一点。这个问题在  @ref
   * Sparsity  模块的文档中也有讨论。
   *
   */
  unsigned int
  max_couplings_between_dofs() const;

  /**
   * 返回位于边界上的自由度的数量，边界上的另一个自由度可以与之耦合。
   * 这个数字与max_couplings_between_dofs()相同，少一个维度。
   * @note
   * 关于这个函数的性能，与max_couplings_per_dofs()相同。可以考虑用一个动态稀疏模式类来代替（见
   * @ref Sparsity  ）。
   *
   */
  unsigned int
  max_couplings_between_boundary_dofs() const;

   /*--------------------------------------*/ 

  /**
   * @name  细胞迭代器函数
   *
   */

  /*  
     * @{    
*
*/

  /**
   * 迭代到第一层使用的单元格  @p level.  。
   *
   */
  cell_iterator
  begin(const unsigned int level = 0) const;

  /**
   * 到第一层活动单元的迭代器  @p level.
   * 如果给定的层不包含任何活动单元（即该层的所有单元都被进一步细化），则该函数返回
   * <code>end_active(level)</code> ，因此，此类的循环是
   * @code
   * for (cell=dof_handler.begin_active(level);
   *      cell!=dof_handler.end_active(level);
   *      ++cell)
   *   {
   *     ...
   *   }
   * @endcode
   * 迭代次数为零，如果这一层没有活动单元，可能会出现这种情况。
   *
   */
  active_cell_iterator
  begin_active(const unsigned int level = 0) const;

  /**
   * 过去结束的迭代器；这个迭代器用于比较过去结束或开始前状态的迭代器。
   *
   */
  cell_iterator
  end() const;

  /**
   * 返回一个迭代器，它是第一个不在给定级别上的迭代器。如果
   * @p level 是最后一个层次，那么这将返回<tt>end()</tt>。
   *
   */
  cell_iterator
  end(const unsigned int level) const;

  /**
   * 返回一个活动的迭代器，它是不在给定级别上的第一个活动迭代器。如果
   * @p level 是最后一层，那么这个返回<tt>end()</tt>。
   *
   */
  active_cell_iterator
  end_active(const unsigned int level) const;

  /**
   * 迭代器到第 @p level.
   * 层的第一个使用的单元格，这将返回一个level_cell_iterator，当dof_indices()被调用时，会返回level
   * dofs。
   *
   */
  level_cell_iterator
  begin_mg(const unsigned int level = 0) const;

  /**
   * 遍历最后一个单元的迭代器  @p level.
   * 当调用dof_indices()时，这将返回一个level_cell_iterator，返回level
   * dofs。
   *
   */
  level_cell_iterator
  end_mg(const unsigned int level) const;

  /**
   * 过去的迭代器；这个迭代器用于比较具有过去或开始前状态的迭代器。
   *
   */
  level_cell_iterator
  end_mg() const;

  /**
   * @name  返回迭代器范围的单元格迭代器函数
   *
   */

  /**
   * 返回一个迭代器范围，该范围包含了组成这个DoFHandler的所有单元格（活动或不活动）。这样的范围对于初始化基于范围的for循环是很有用的，正如C++11所支持的那样。
   * 参见active_cell_iterators()文档中的例子。      @return
   * 半开范围  <code>[this->begin(), this->end())</code>
   * @ingroup CPP11
   *
   */
  IteratorRange<cell_iterator>
  cell_iterators() const;

  /**
   * #Range-based_for_loop">here</a>）。一个例子是，如果没有基于范围的for循环，人们往往会写出如下的代码。
   * @code
   * DoFHandler<dim> dof_handler;
   * ...
   * typename DoFHandler<dim>::active_cell_iterator
   *   cell = dof_handler.begin_active(),
   *   endc = dof_handler.end();
   * for (; cell!=endc; ++cell)
   *   {
   *     fe_values.reinit (cell);
   *     ...do the local integration on 'cell'...;
   *   }
   * @endcode
   * 使用C++11的基于范围的for循环，现在这完全等同于以下的代码。
   * @code
   * DoFHandler<dim> dof_handler;
   * ...
   * for (const auto &cell : dof_handler.active_cell_iterators())
   *   {
   *     fe_values.reinit (cell);
   *     ...do the local integration on 'cell'...;
   *   }
   * @endcode
   * @return  半开范围<code>[this->begin_active(), this->end())</code>。
   * @ingroup CPP11
   *
   */
  IteratorRange<active_cell_iterator>
  active_cell_iterators() const;

  /**
   * 返回一个迭代器范围，该范围包含了组成这个DoFHandler的所有单元格（无论是否激活）的水平单元格的形式。这样的范围对于初始化基于范围的for循环是很有用的，正如C++11所支持的那样。
   * 参见active_cell_iterators()文档中的例子。      @return
   * 半开范围<code>[this->begin_mg(), this->end_mg())</code>。
   * @ingroup CPP11
   *
   */
  IteratorRange<level_cell_iterator>
  mg_cell_iterators() const;

  /**
   * 返回一个迭代器范围，该范围包含所有构成该DoFHandler的给定级别的单元格（无论是否激活）。这样的范围对于初始化基于范围的for循环是很有用的，正如C++11所支持的那样。
   * 参见active_cell_iterators()文档中的例子。      @param[in]  level
   * 这个三角结构的细化层次中的一个给定的层次。
   * @return  半开放范围<code>[this->begin(level),
   * this->end(level))</code>  @pre  level必须小于this->n_levels()。
   * @ingroup CPP11
   *
   */
  IteratorRange<cell_iterator>
  cell_iterators_on_level(const unsigned int level) const;

  /**
   * 返回一个迭代器范围，该范围包含所有构成此DoFHandler的给定级别的活动单元。这样的范围对于初始化基于范围的for循环是很有用的，正如C++11所支持的那样。
   * 参见active_cell_iterators()文档中的例子。      @param[in]  level
   * 这个三角结构的细化层次中的一个给定的层次。
   * @return  半开放范围<code>[this->begin_active(level),
   * this->end(level))</code>  @pre  level必须小于this->n_levels()。
   * @ingroup CPP11
   *
   */
  IteratorRange<active_cell_iterator>
  active_cell_iterators_on_level(const unsigned int level) const;

  /**
   * 返回一个迭代器范围，该范围包含所有构成该DoFHandler的给定级别的单元格（无论是否激活）的level-cell形式。这样的范围对于初始化基于范围的for循环是很有用的，正如C++11所支持的那样。
   * 参见active_cell_iterators()文档中的例子。      @param[in]  level
   * 在这个三角结构的细化层次中的一个给定级别。
   * @return  半开范围<code>[this->begin_mg(level),
   * this->end_mg(level))</code>  @pre  level必须小于this->n_levels()。
   * @ingroup CPP11
   *
   */
  IteratorRange<level_cell_iterator>
  mg_cell_iterators_on_level(const unsigned int level) const;

  /*  @} .   
*
*/


   /*---------------------------------------*/ 


  /**
   * 返回全局自由度的数量。如果当前对象自己处理所有的自由度（即使你可能打算并行解决你的线性系统，如 step-17
   * 或 step-18
   * ），那么这个数字等于本地拥有的自由度数，因为这个对象不知道你想用它做什么，并认为它拥有它知道的每一个自由度。
   * 另一方面，如果这个对象在一个
   * parallel::distributed::Triangulation
   * 对象上操作，那么这个函数返回全局自由度的数量，在所有处理器上累积。
   * 在这两种情况下，返回的数字都包括那些被悬挂节点约束的自由度，见
   * @ref constraints  。
   * 从数学上讲，这个函数返回的数字等于有限元空间的尺寸（不考虑约束），对应于(i)它所定义的网格，和(ii)当前对象使用的有限元。当然，它也等于跨越这个空间的形状函数的数量。
   *
   */
  types::global_dof_index
  n_dofs() const;

  /**
   * 某一层次上的多层次自由度的（全局）数量。
   * 如果没有给这个层次分配自由度，则返回
   * numbers::invalid_dof_index.  否则返回这个层次的自由度数。
   *
   */
  types::global_dof_index
  n_dofs(const unsigned int level) const;

  /**
   * 返回位于边界上的本地拥有的自由度的数量。
   *
   */
  types::global_dof_index
  n_boundary_dofs() const;

  /**
   * 返回位于边界上有边界指标的部分的本地拥有的自由度数，这些边界指标列在给定的集合中。
   * 使用 @p map 而不是 @p set 的原因与
   * DoFTools::make_boundary_sparsity_pattern()
   * 的变体文件中描述的相同，该变体需要一个地图。
   * 然而，这个函数还有一个重载，需要一个 @p set
   * 的参数（见下文）。
   *
   */
  template <typename number>
  types::global_dof_index
  n_boundary_dofs(
    const std::map<types::boundary_id, const Function<spacedim, number> *>
      &boundary_ids) const;

  /**
   * 返回位于边界部分的自由度的数量，这些部分的边界指标列在给定的集合中。的
   *
   */
  types::global_dof_index
  n_boundary_dofs(const std::set<types::boundary_id> &boundary_ids) const;

  /**
   * 访问一个告知自由度处理程序块结构的对象。    如果在distribution_dofs()中使用了FESystem，自由度自然会被分成几个 @ref GlossBlock "块"
   * 。
   * 对于每个基本元素，出现的块数与它的倍数一样多。
   * 在distribut_dofs()的最后，每个块中的自由度被计算出来，并存储在BlockInfo对象中，可以在这里访问。如果你之前调用了distribution_mg_dofs()，那么在多网格层次结构的每一层都会做同样的事情。此外，每个单元上的块结构可以通过调用initialize_local_block_info()在此对象中生成。
   *
   */
  const BlockInfo &
  block_info() const;

  /**
   * 返回属于这个过程的自由度的数量。
   * 如果这是一个连续的DoFHandler，那么这个结果等于n_dofs()产生的结果。这里，"顺序
   * "意味着要么整个程序不使用MPI，要么使用MPI但只使用一个MPI进程，要么有多个MPI进程，但这个DoFHandler建立的Triangulation只在一个MPI进程上工作）。
   * 另一方面，如果我们在一个 parallel::distributed::Triangulation
   * 或 parallel::shared::Triangulation,
   * 上操作，那么它只包括当前处理器拥有的自由度。请注意，在这种情况下，这并不包括所有分布在当前处理器的网格图像上的自由度：特别是，这个处理器所拥有的单元和其他处理器所拥有的单元之间的界面上的一些自由度可能是他们的，而幽灵单元上的自由度也不一定包括在内。
   *
   */
  types::global_dof_index
  n_locally_owned_dofs() const;

  /**
   * 返回一个IndexSet，描述本地拥有的自由度集合，作为0...n_dofs()的一个子集。这个集合的元素数等于n_locally_owned_dofs()。
   *
   */
  const IndexSet &
  locally_owned_dofs() const;

  /**
   * 返回一个IndexSet，描述用于给定多网格级别的本地拥有的DoF集合，作为0...n_dofs(level)的一个子集。
   *
   */
  const IndexSet &
  locally_owned_mg_dofs(const unsigned int level) const;

  /**
   * 返回一个对该对象所使用的indexth有限元对象的常数引用。
   *
   */
  const FiniteElement<dim, spacedim> &
  get_fe(const unsigned int index = 0) const;

  /**
   * 返回对该对象所使用的有限元对象集合的常数引用。
   *
   */
  const hp::FECollection<dim, spacedim> &
  get_fe_collection() const;

  /**
   * 返回对该对象所依据的三角形的常数引用。
   *
   */
  const Triangulation<dim, spacedim> &
  get_triangulation() const;

  /**
   * 返回底层三角法所使用的MPI通信器。
   *
   */
  MPI_Comm
  get_communicator() const;

  /**
   * 每当考虑用 parallel::distributed::Triangulation
   * 作为底层三角结构进行序列化时，我们也需要考虑在所有活动单元上存储活动的FE指数。
   * 这个函数注册了这些指数，当
   * parallel::distributed::Triangulation::save()
   * 函数在底层三角结构上被调用时，这些指数将被存储。
   * @note  目前只对 parallel::distributed::Triangulation.
   * 类型的三角形实现，如果注册了不同的类型，将触发一个断言。
   * @see   parallel::distributed::SolutionTransfer
   * 的文档有关于序列化的进一步信息。
   *
   */
  void
  prepare_for_serialization_of_active_fe_indices();

  /**
   * 每当考虑用 parallel::distributed::Triangulation
   * 作为底层三角的序列化时，我们也需要考虑在所有活动单元上存储活动FE指数。
   * 这个函数反序列化并将之前存储的活动FE指数分配到所有活动单元上。
   * @note  目前只针对类型的三角形实现
   * parallel::distributed::Triangulation.
   * 如果注册了不同的类型，将触发一个断言。      @see
   * parallel::distributed::SolutionTransfer
   * 的文档有关于序列化的进一步信息。
   *
   */
  void
  deserialize_active_fe_indices();

  /**
   * 确定此对象的内存消耗（以字节为单位）的估计值。
   * 这个函数是虚拟的，因为dof
   * handler对象可能通过指向这个基类的指针来访问，尽管实际对象可能是一个派生类。
   *
   */
  virtual std::size_t
  memory_consumption() const;

  /**
   * 使用[BOOST序列化库](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html)将此对象的数据写入一个流中，以便进行序列化。
   *
   */
  template <class Archive>
  void
  save(Archive &ar, const unsigned int version) const;

  /**
   * 为了序列化的目的，使用[BOOST序列化库](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html)从流中读取此对象的数据。
   *
   */
  template <class Archive>
  void
  load(Archive &ar, const unsigned int version);

#ifdef DOXYGEN
  /**
   * 使用[BOOST序列化库](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html)从流中写入和读取此对象的数据，以达到序列化的目的。
   *
   */
  template <class Archive>
  void
  serialize(Archive &archive, const unsigned int version);
#else
  // This macro defines the serialize() method that is compatible with
  // the templated save() and load() method that have been implemented.
  BOOST_SERIALIZATION_SPLIT_MEMBER()
#endif

  /**
   * 异常情况
   *
   */
  DeclException0(ExcNoFESelected);
  /**
   * 异常情况
   * @ingroup Exceptions
   *
   */
  DeclException0(ExcInvalidBoundaryIndicator);
  /**
   * 异常情况
   * @ingroup Exceptions
   *
   */
  DeclException1(ExcInvalidLevel,
                 int,
                 << "The given level " << arg1
                 << " is not in the valid range!");
  /**
   * 异常情况
   * @ingroup Exceptions
   *
   */
  DeclException1(ExcNewNumbersNotConsecutive,
                 types::global_dof_index,
                 << "The given list of new dof indices is not consecutive: "
                 << "the index " << arg1 << " does not exist.");
  /**
   * 异常情况
   *
   */
  DeclException2(ExcInvalidFEIndex,
                 int,
                 int,
                 << "The mesh contains a cell with an active FE index of "
                 << arg1 << ", but the finite element collection only has "
                 << arg2 << " elements");

  /**
   * 当DoFHandler没有hp-capabilities时，某种功能没有意义时使用的异常。
   *
   */
  DeclExceptionMsg(ExcOnlyAvailableWithHP,
                   "The current function doesn't make sense when used with a "
                   "DoFHandler without hp-capabilities.");

  /**
   * 当DoFHandler有hp-capabilities时，某项功能没有实现时使用的异常。
   *
   */
  DeclExceptionMsg(ExcNotImplementedWithHP,
                   "The current function has not yet been implemented for a "
                   "DoFHandler with hp-capabilities.");

private:
  /**
   * 一个数据结构，用于存储与特定顶点相关的DoF指数。与单元不同，顶点生活在多网格层次结构的多个层面上；因此，我们需要为每个顶点存储它所处的每个层面的DoF指数。这个类就是这样做的。
   *
   */
  class MGVertexDoFs
  {
  public:
    /**
     * 构造函数。
     *
     */
    MGVertexDoFs();

    /**
     * 调用一个函数来分配必要的内存量，以存储该顶点在给定（包括）层次范围内的DoF的索引。
     *
     */
    void
    init(const unsigned int coarsest_level,
         const unsigned int finest_level,
         const unsigned int dofs_per_vertex);

    /**
     * 返回该结构存储数据的最粗的层次。
     *
     */
    unsigned int
    get_coarsest_level() const;

    /**
     * 返回该结构存储数据的最细层次。
     *
     */
    unsigned int
    get_finest_level() const;

    /**
     * 返回当前顶点存储的给定级别的 <code>dof_number</code>
     * 个自由度的索引。
     *
     */
    types::global_dof_index
    get_index(const unsigned int level,
              const unsigned int dof_number,
              const unsigned int dofs_per_vertex) const;

    /**
     * 为当前顶点存储的给定级别的 <code>dof_number</code>
     * 个自由度的索引设置为 <code>index</code>  。
     *
     */
    void
    set_index(const unsigned int            level,
              const unsigned int            dof_number,
              const unsigned int            dofs_per_vertex,
              const types::global_dof_index index);

  private:
    /**
     * 此对象存储自由度指数的最粗层次。
     *
     */
    unsigned int coarsest_level;

    /**
     * 该对象存储DoF指数的最细级别。
     *
     */
    unsigned int finest_level;

    /**
     * 一个指向数组的指针，在这个数组中，我们存储这个顶点存在的各个层次上的DoFs的索引。
     * 属于 @p level 的DoF的起始偏移由 <code>n_dofs_per_vertex()
     * (level-coarsest_level)</code> 给出。 因此， @p
     * n_dofs_per_vertex()
     * 必须作为一个参数传递给设置或读取索引的函数。
     *
     */
    std::unique_ptr<types::global_dof_index[]> indices;
  };

  /**
   * 每当底层三角结构通过h/p-refinement/coarsening和序列化发生变化时，单元格的活动FE索引就需要被转移。这个结构存储了该过程中所需要的所有临时信息。
   *
   */
  struct ActiveFEIndexTransfer
  {
    /**
     * 用于临时存储迭代器和未来持续存在的单元格的活动FE索引的容器。
     *
     */
    std::map<const cell_iterator, const unsigned int> persisting_cells_fe_index;

    /**
     * 容器用于临时存储将被精炼的单元格的迭代器和未来的活动FE索引。
     *
     */
    std::map<const cell_iterator, const unsigned int> refined_cells_fe_index;

    /**
     * 容器，用于临时存储粗化后将保留的父单元格的迭代器和未来有效的FE索引。
     *
     */
    std::map<const cell_iterator, const unsigned int> coarsened_cells_fe_index;

    /**
     * 容器，用于临时存储每个本地拥有的单元格的活动FE索引，以便在
     * parallel::distributed::Triangulation 对象之间转移。
     *
     */
    std::vector<unsigned int> active_fe_indices;

    /**
     * 帮助对象，在细化/粗化和序列化过程中转移
     * parallel::distributed::Triangulation
     * 对象上的所有活动FE指数。
     *
     */
    std::unique_ptr<
      parallel::distributed::
        CellDataTransfer<dim, spacedim, std::vector<unsigned int>>>
      cell_data_transfer;
  };

  /**
   * 一个包含块结构信息的对象。
   *
   */
  BlockInfo block_info_object;

  /**
   * 表示当前DoFHandler是否有hp-能力的布尔值。
   *
   */
  bool hp_capability_enabled;

  /**
   * 要工作的三角结构的地址。
   *
   */
  SmartPointer<const Triangulation<dim, spacedim>, DoFHandler<dim, spacedim>>
    tria;

  /**
   * 存储一个 hp::FECollection
   * 对象。如果在这个对象的初始化过程中只使用了一个FiniteElement，它就包含了（一个）FiniteElement。
   *
   */
  hp::FECollection<dim, spacedim> fe_collection;

  /**
   * 一个描述自由度应如何分配和重新编号的对象。
   *
   */
  std::unique_ptr<dealii::internal::DoFHandlerImplementation::Policy::
                    PolicyBase<dim, spacedim>>
    policy;

  /**
   * 一个包含各种数字的结构，这些数字表征了这个对象所工作的自由度。
   * 对于这个结构的大多数成员，在这个类中有一个访问函数来返回其值。
   *
   */
  dealii::internal::DoFHandlerImplementation::NumberCache number_cache;

  /**
   * 像number_cache一样的数据结构，但为每个多网格层次。
   *
   */
  std::vector<dealii::internal::DoFHandlerImplementation::NumberCache>
    mg_number_cache;

  /**
   * 缓存的所有活动单元的自由度指数。通过cell_dof_cache_ptr（CRS方案）来识别单元在向量中的适当位置。
   *
   */
  mutable std::vector<std::vector<types::global_dof_index>>
    cell_dof_cache_indices;

  /**
   * 指向cell_dof_cache_indices中活动单元的第一个缓存自由度的指针（通过级别和级别索引识别）。
   *
   */
  mutable std::vector<std::vector<offset_type>> cell_dof_cache_ptr;

  /**
   * 所有相关活动有限元的每个d+1几何对象（3D：顶点、线、四边形、六边形）的自由度指数。通过object_dof_ptr（CRS方案）识别适当的位置。
   *
   */
  mutable std::vector<std::array<std::vector<types::global_dof_index>, dim + 1>>
    object_dof_indices;

  /**
   * 指向所有相关活动有限元的几何对象的第一个缓存自由度的指针。
   * @note 在正常模式下，可以直接访问这个数据结构。
   * 在hp模式下，需要通过hp_object_fe_indices/hp_object_fe_ptr进行转接。
   *
   */
  mutable std::vector<std::array<std::vector<offset_type>, dim + 1>>
    object_dof_ptr;

  /**
   * 每个几何对象的有效FE指数。通过hp_object_fe_ptr(CRS方案)来识别一个单元在向量中的适当位置。
   *
   */
  mutable std::array<std::vector<active_fe_index_type>, dim + 1>
    hp_object_fe_indices;

  /**
   * 指向一个几何对象的第一个FE索引的指针。
   *
   */
  mutable std::array<std::vector<offset_type>, dim + 1> hp_object_fe_ptr;

  /**
   * 一个活动单元的活动FE索引（由级别和级别索引识别）。
   * 这个向量只在hp模式下使用。
   *
   */
  mutable std::vector<std::vector<active_fe_index_type>>
    hp_cell_active_fe_indices;

  /**
   * 一个活动单元的未来FE指数（由级别和级别指数识别）。
   * 这个向量只在hp模式下使用。
   *
   */
  mutable std::vector<std::vector<active_fe_index_type>>
    hp_cell_future_fe_indices;

  /**
   * 一个数组，用于存储位于顶点的水平自由度的索引。
   *
   */
  std::vector<MGVertexDoFs> mg_vertex_dofs;

  /**
   * 用于存储不同多网格层次的自由度编号的空间。
   *
   */
  std::vector<
    std::unique_ptr<dealii::internal::DoFHandlerImplementation::DoFLevel<dim>>>
    mg_levels;

  /**
   * 用于存储多网格背景下面的自由度数的空间。
   *
   */
  std::unique_ptr<dealii::internal::DoFHandlerImplementation::DoFFaces<dim>>
    mg_faces;

  /**
   * 我们把我们的数据结构嵌入到一个指针中，以控制所有与传输有关的数据只在实际传输过程中存在。
   *
   */
  std::unique_ptr<ActiveFEIndexTransfer> active_fe_index_transfer;

  /**
   * 一个连接的列表，这个对象用它连接到三角区，以获得三角区变化时的信息。
   *
   */
  std::vector<boost::signals2::connection> tria_listeners;

  /**
   * 这个对象连接到三角测量的连接的列表。当数据由于细化或重新分区而需要传输时，它们会被特别触发。只有在hp模式下才有效。
   *
   */
  std::vector<boost::signals2::connection> tria_listeners_for_transfer;

  /**
   * 释放所有用于非多重网格数据结构的内存。
   *
   */
  void
  clear_space();

  /**
   * 释放所有用于多栅格数据结构的内存。
   *
   */
  void
  clear_mg_space();

  /**
   * 返回指定对象的dof索引。
   *
   */
  template <int structdim>
  types::global_dof_index
  get_dof_index(const unsigned int obj_level,
                const unsigned int obj_index,
                const unsigned int fe_index,
                const unsigned int local_index) const;

  /**
   * 返回指定对象的dof索引。
   *
   */
  template <int structdim>
  void
  set_dof_index(const unsigned int            obj_level,
                const unsigned int            obj_index,
                const unsigned int            fe_index,
                const unsigned int            local_index,
                const types::global_dof_index global_index) const;

  /**
   * 设置DoFHandler策略。
   *
   */
  void
  setup_policy();

  /**
   * 设置与底层三角结构的信号的连接。
   *
   */
  void
  connect_to_triangulation_signals();

  /**
   * 为活动和未来的fe_indices创建默认表。
   * 活动指数用一个零指标初始化，这意味着fe[0]将被默认使用。未来的指数被初始化为一个无效的指标，这意味着默认情况下不安排任何p-适应。
   * 这个方法在构造时和底层三角结构被创建时被调用。这确保每个单元都有一个有效的活动和未来的fe_index。
   *
   */
  void
  create_active_fe_table();

  /**
   * 更新活动和未来fe_indices的表格。
   * 每当底层三角结构发生变化时（无论是通过适应还是反序列化），活动和未来的FE指数表将被调整为三角结构的当前结构。活动指数和未来指数的缺失值将被初始化为其默认值（参见create_active_fe_table()）。
   * 这个方法在细化后和反序列化后被调用。这可以确保每个单元格都有一个有效的活动和未来的fe_index。
   *
   */
  void
  update_active_fe_table();

  /**
   * 一个函数，它将在相关的三角函数或
   * parallel::shared::Triangulation
   * 被修改之前通过三角函数信号被触发。
   * 在细化发生之前，存储所有将被细化或粗化的单元的活动FE指数的函数，以便在细化之后可以再次设置这些指数。
   *
   */
  void
  pre_transfer_action();

  /**
   * 就在相关的三角化或 parallel::shared::Triangulation
   * 被修改后，将通过三角化信号触发的函数。
   * 恢复所有被细化或粗化的单元的活动FE指数的函数。
   *
   */
  void
  post_transfer_action();

  /**
   * 在相关的 parallel::distributed::Triangulation
   * 被修改之前，将通过三角测量信号触发的函数。
   * 将所有活动的FE指数存储在本地拥有的单元上的函数，以分配给所有参与的处理器。
   *
   */
  void
  pre_distributed_transfer_action();

  /**
   * 在相关的 parallel::distributed::Triangulation
   * 被修改后，就会通过三角信号触发一个函数。
   * 恢复本地拥有的单元上所有活动的FE指数的函数，这些单元已经被通信了。
   *
   */
  void
  post_distributed_transfer_action();


  // Make accessor objects friends.
  template <int, int, int, bool>
  friend class dealii::DoFAccessor;
  template <int, int, bool>
  friend class dealii::DoFCellAccessor;
  friend struct dealii::internal::DoFAccessorImplementation::Implementation;
  friend struct dealii::internal::DoFCellAccessorImplementation::Implementation;

  // Likewise for DoFLevel objects since they need to access the vertex dofs
  // in the functions that set and retrieve vertex dof indices.
  friend struct dealii::internal::DoFHandlerImplementation::Implementation;
  friend struct dealii::internal::hp::DoFHandlerImplementation::Implementation;
  friend struct dealii::internal::DoFHandlerImplementation::Policy::
    Implementation;

  // explicitly check for sensible template arguments, but not on windows
  // because MSVC creates bogus warnings during normal compilation
#ifndef DEAL_II_MSVC
  static_assert(dim <= spacedim,
                "The dimension <dim> of a DoFHandler must be less than or "
                "equal to the space dimension <spacedim> in which it lives.");
#endif
};

namespace internal
{
  namespace hp
  {
    namespace DoFHandlerImplementation
    {
      /**
       * 给定一个hp模式的DoFHandler对象，确保用户为本地拥有的单元设置的未来FE指数也被传达给所有其他相关单元。
       * 对于 parallel::shared::Triangulation
       * 对象，该信息同时分布在幽灵和人工单元上。
       * 在使用 parallel::distributed::Triangulation
       * 的情况下，指数只传达给幽灵单元。
       *
       */
      template <int dim, int spacedim>
      void
      communicate_future_fe_indices(DoFHandler<dim, spacedim> &dof_handler);

      /**
       * 返回来自整个 hp::FECollection
       * 的有限元的索引，该有限元被分配给 @p parent.
       * 的子单元的未来有限元所支配。
       * 我们在该单元的子单元上的未来有限元中找到相应的有限元。如果没有一个符合条件，我们将搜索范围扩大到整个
       * hp::FECollection,
       * ，这是描述最小的有限元空间的元素，包括分配给子单元的所有未来有限元。如果该函数根本无法找到有限元，将触发一个断言。
       * 通过这种方式，我们在hp-context中h-coarsening的情况下确定父单元的有限元。
       * @note
       * 这个函数只能在直接父单元上调用，即子单元都是活动的非活动单元。
       * @note  在 parallel::shared::Triangulation
       * 对象中，兄弟姐妹单元可以是幽灵单元，请确保未来的FE指数已经正确地与communication_future_fe_indices()沟通。否则，在不同的处理器上，结果可能不同。没有检查未来FE指数的一致性。
       *
       */
      template <int dim, int spacedim = dim>
      unsigned int
      dominated_future_fe_on_children(
        const typename DoFHandler<dim, spacedim>::cell_iterator &parent);

      /**
       * 异常情况
       *
       */
      DeclExceptionMsg(
        ExcNoDominatedFiniteElementOnChildren,
        "No FiniteElement has been found in your FECollection that is "
        "dominated by all children of a cell you are trying to coarsen!");
    } // namespace DoFHandlerImplementation
  }   // namespace hp
} // namespace internal

#ifndef DOXYGEN

 /* ----------------------- Inline functions ----------------------------------
 */ 


template <int dim, int spacedim>
inline bool
DoFHandler<dim, spacedim>::has_hp_capabilities() const
{
  return hp_capability_enabled;
}



template <int dim, int spacedim>
inline bool
DoFHandler<dim, spacedim>::has_level_dofs() const
{
  return mg_number_cache.size() > 0;
}



template <int dim, int spacedim>
inline bool
DoFHandler<dim, spacedim>::has_active_dofs() const
{
  return number_cache.n_global_dofs > 0;
}



template <int dim, int spacedim>
inline types::global_dof_index
DoFHandler<dim, spacedim>::n_dofs() const
{
  return number_cache.n_global_dofs;
}



template <int dim, int spacedim>
inline types::global_dof_index
DoFHandler<dim, spacedim>::n_dofs(const unsigned int level) const
{
  Assert(has_level_dofs(),
         ExcMessage(
           "n_dofs(level) can only be called after distribute_mg_dofs()"));
  Assert(level < mg_number_cache.size(), ExcInvalidLevel(level));
  return mg_number_cache[level].n_global_dofs;
}



template <int dim, int spacedim>
types::global_dof_index
DoFHandler<dim, spacedim>::n_locally_owned_dofs() const
{
  return number_cache.n_locally_owned_dofs;
}



template <int dim, int spacedim>
const IndexSet &
DoFHandler<dim, spacedim>::locally_owned_dofs() const
{
  return number_cache.locally_owned_dofs;
}



template <int dim, int spacedim>
const IndexSet &
DoFHandler<dim, spacedim>::locally_owned_mg_dofs(const unsigned int level) const
{
  Assert(level < this->get_triangulation().n_global_levels(),
         ExcMessage("The given level index exceeds the number of levels "
                    "present in the triangulation"));
  Assert(
    mg_number_cache.size() == this->get_triangulation().n_global_levels(),
    ExcMessage(
      "The level dofs are not set up properly! Did you call distribute_mg_dofs()?"));
  return mg_number_cache[level].locally_owned_dofs;
}



template <int dim, int spacedim>
inline const FiniteElement<dim, spacedim> &
DoFHandler<dim, spacedim>::get_fe(const unsigned int number) const
{
  Assert(fe_collection.size() > 0,
         ExcMessage("No finite element collection is associated with "
                    "this DoFHandler"));
  return fe_collection[number];
}



template <int dim, int spacedim>
inline const hp::FECollection<dim, spacedim> &
DoFHandler<dim, spacedim>::get_fe_collection() const
{
  return fe_collection;
}



template <int dim, int spacedim>
inline const Triangulation<dim, spacedim> &
DoFHandler<dim, spacedim>::get_triangulation() const
{
  Assert(tria != nullptr,
         ExcMessage("This DoFHandler object has not been associated "
                    "with a triangulation."));
  return *tria;
}



template <int dim, int spacedim>
inline MPI_Comm
DoFHandler<dim, spacedim>::get_communicator() const
{
  Assert(tria != nullptr,
         ExcMessage("This DoFHandler object has not been associated "
                    "with a triangulation."));
  return tria->get_communicator();
}



template <int dim, int spacedim>
inline const BlockInfo &
DoFHandler<dim, spacedim>::block_info() const
{
  Assert(this->hp_capability_enabled == false, ExcNotImplementedWithHP());

  return block_info_object;
}



template <int dim, int spacedim>
template <typename number>
types::global_dof_index
DoFHandler<dim, spacedim>::n_boundary_dofs(
  const std::map<types::boundary_id, const Function<spacedim, number> *>
    &boundary_ids) const
{
  Assert(!(dim == 2 && spacedim == 3) || this->hp_capability_enabled == false,
         ExcNotImplementedWithHP());

  // extract the set of boundary ids and forget about the function object
  // pointers
  std::set<types::boundary_id> boundary_ids_only;
  for (typename std::map<types::boundary_id,
                         const Function<spacedim, number> *>::const_iterator p =
         boundary_ids.begin();
       p != boundary_ids.end();
       ++p)
    boundary_ids_only.insert(p->first);

  // then just hand everything over to the other function that does the work
  return n_boundary_dofs(boundary_ids_only);
}



namespace internal
{
  /**
   * 返回一个字符串，代表给定参数的动态类型。
   * 这与typeid(...).name()的作用基本相同，但事实证明这在Intel
   * 13+上是坏的。    定义在dof_handler.cc中。
   *
   */
  template <int dim, int spacedim>
  std::string
  policy_to_string(const dealii::internal::DoFHandlerImplementation::Policy::
                     PolicyBase<dim, spacedim> &policy);
} // namespace internal



template <int dim, int spacedim>
template <class Archive>
void
DoFHandler<dim, spacedim>::save(Archive &ar, const unsigned int) const
{
  if (this->hp_capability_enabled)
    {
      ar & this->object_dof_indices;
      ar & this->object_dof_ptr;

      ar & this->cell_dof_cache_indices;
      ar & this->cell_dof_cache_ptr;

      ar & this->hp_cell_active_fe_indices;
      ar & this->hp_cell_future_fe_indices;

      ar &hp_object_fe_ptr;
      ar &hp_object_fe_indices;

      ar &number_cache;

      ar &mg_number_cache;

      // write out the number of triangulation cells and later check during
      // loading that this number is indeed correct; same with something that
      // identifies the policy
      const unsigned int n_cells = this->tria->n_cells();
      std::string        policy_name =
        dealii::internal::policy_to_string(*this->policy);

      ar &n_cells &policy_name;
    }
  else
    {
      ar & this->block_info_object;
      ar &number_cache;

      ar & this->object_dof_indices;
      ar & this->object_dof_ptr;

      ar & this->cell_dof_cache_indices;
      ar & this->cell_dof_cache_ptr;

      // write out the number of triangulation cells and later check during
      // loading that this number is indeed correct; same with something that
      // identifies the FE and the policy
      unsigned int n_cells     = this->tria->n_cells();
      std::string  fe_name     = this->get_fe(0).get_name();
      std::string  policy_name = internal::policy_to_string(*this->policy);

      ar &n_cells &fe_name &policy_name;
    }
}



template <int dim, int spacedim>
template <class Archive>
void
DoFHandler<dim, spacedim>::load(Archive &ar, const unsigned int)
{
  if (this->hp_capability_enabled)
    {
      ar & this->object_dof_indices;
      ar & this->object_dof_ptr;

      ar & this->cell_dof_cache_indices;
      ar & this->cell_dof_cache_ptr;

      ar & this->hp_cell_active_fe_indices;
      ar & this->hp_cell_future_fe_indices;

      ar &hp_object_fe_ptr;
      ar &hp_object_fe_indices;

      ar &number_cache;

      ar &mg_number_cache;

      // these are the checks that correspond to the last block in the save()
      // function
      unsigned int n_cells;
      std::string  policy_name;

      ar &n_cells &policy_name;

      AssertThrow(
        n_cells == this->tria->n_cells(),
        ExcMessage(
          "The object being loaded into does not match the triangulation "
          "that has been stored previously."));
      AssertThrow(
        policy_name == dealii::internal::policy_to_string(*this->policy),
        ExcMessage("The policy currently associated with this DoFHandler (" +
                   dealii::internal::policy_to_string(*this->policy) +
                   ") does not match the one that was associated with the "
                   "DoFHandler previously stored (" +
                   policy_name + ")."));
    }
  else
    {
      ar & this->block_info_object;
      ar &number_cache;

      object_dof_indices.clear();

      object_dof_ptr.clear();

      ar & this->object_dof_indices;
      ar & this->object_dof_ptr;

      ar & this->cell_dof_cache_indices;
      ar & this->cell_dof_cache_ptr;

      // these are the checks that correspond to the last block in the save()
      // function
      unsigned int n_cells;
      std::string  fe_name;
      std::string  policy_name;

      ar &n_cells &fe_name &policy_name;

      AssertThrow(
        n_cells == this->tria->n_cells(),
        ExcMessage(
          "The object being loaded into does not match the triangulation "
          "that has been stored previously."));
      AssertThrow(
        fe_name == this->get_fe(0).get_name(),
        ExcMessage(
          "The finite element associated with this DoFHandler does not match "
          "the one that was associated with the DoFHandler previously stored."));
      AssertThrow(policy_name == internal::policy_to_string(*this->policy),
                  ExcMessage(
                    "The policy currently associated with this DoFHandler (" +
                    internal::policy_to_string(*this->policy) +
                    ") does not match the one that was associated with the "
                    "DoFHandler previously stored (" +
                    policy_name + ")."));
    }
}



template <int dim, int spacedim>
inline types::global_dof_index
DoFHandler<dim, spacedim>::MGVertexDoFs::get_index(
  const unsigned int level,
  const unsigned int dof_number,
  const unsigned int dofs_per_vertex) const
{
  Assert((level >= coarsest_level) && (level <= finest_level),
         ExcInvalidLevel(level));
  return indices[dofs_per_vertex * (level - coarsest_level) + dof_number];
}



template <int dim, int spacedim>
inline void
DoFHandler<dim, spacedim>::MGVertexDoFs::set_index(
  const unsigned int            level,
  const unsigned int            dof_number,
  const unsigned int            dofs_per_vertex,
  const types::global_dof_index index)
{
  Assert((level >= coarsest_level) && (level <= finest_level),
         ExcInvalidLevel(level));
  indices[dofs_per_vertex * (level - coarsest_level) + dof_number] = index;
}



extern template class DoFHandler<1, 1>;
extern template class DoFHandler<1, 2>;
extern template class DoFHandler<1, 3>;
extern template class DoFHandler<2, 2>;
extern template class DoFHandler<2, 3>;
extern template class DoFHandler<3, 3>;


#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif


