//include/deal.II-translator/numerics/solution_transfer_0.txt
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

#ifndef dealii_solution_transfer_h
#  define dealii_solution_transfer_h


 /*----------------------------   solutiontransfer.h     ----------------------*/ 


#  include <deal.II/base/config.h>

#  include <deal.II/base/exceptions.h>
#  include <deal.II/base/smartpointer.h>

#  include <deal.II/dofs/dof_handler.h>

#  include <deal.II/lac/vector.h>

#  include <vector>

DEAL_II_NAMESPACE_OPEN

/**
 * 该类实现了离散FE函数（例如：一个求解向量）从一个网格到另一个网格的转移，该网格是通过一个细化和/或粗化步骤从第一个网格得到的。在插值过程中，向量被重新初始化为新的大小，并填充了插值的数值。这个类在
 * step-15 、 step-26 、 step-31 和 step-33
 * 的教程程序中使用。这个类的一个版本可以在平行三角形上工作，可作为
 * parallel::distributed::SolutionTransfer. 。 <h3>Usage</h3>
 * 该类以两种不同的方式实现算法。  <ul>   <li>  如果网格只被细化（即没有单元格被粗化），那么使用 @p SolutionTransfer 如下。
 *
 * @code
 * SolutionTransfer<dim, Vector<double> > soltrans(*dof_handler);
 *
 * // flag some cells for refinement, e.g.
 * GridRefinement::refine_and_coarsen_fixed_fraction(*tria,
 *                                                 error_indicators,
 *                                                 0.3,
 *                                                 0);
 * // prepare the triangulation for refinement,
 * tria->prepare_coarsening_and_refinement();
 *
 * // tell the SolutionTransfer object that we intend to do pure refinement,
 * soltrans.prepare_for_pure_refinement();
 *
 * // actually execute the refinement,
 * tria->execute_coarsening_and_refinement();
 *
 * // and redistribute dofs.
 * dof_handler->distribute_dofs (fe);
 * @endcode
 *
 * 然后继续做
 *
 * @code
 * // take a copy of the solution vector
 * Vector<double> solution_old(solution);
 *
 * // resize solution vector to the correct size, as the @p refine_interpolate
 * // function requires the vectors to be of right sizes
 * solution.reinit(dof_handler->n_dofs());
 *
 * // and finally interpolate
 * soltrans.refine_interpolate(solution_old, solution);
 * @endcode
 *
 * 尽管允许多次调用 @p refine_interpolate
 * 函数，例如用于插值几个解向量，但存在以下同时插值几个函数的可能性。
 *
 * @code
 * std::vector<Vector<double> > solutions_old(n_vectors, Vector<double> (n));
 * ...
 * std::vector<Vector<double> > solutions(n_vectors, Vector<double> (n));
 * soltrans.refine_interpolate(solutions_old, solutions);
 * @endcode
 * 这在几个教程程序中使用，例如  step-31  。
 * <li>  如果网格中的单元格将被粗化，那么使用 @p
 * SolutionTransfer，如下。
 *
 * @code
 * SolutionTransfer<dim, Vector<double> > soltrans(*dof_handler);
 *
 * // flag some cells for refinement and coarsening, e.g.
 * GridRefinement::refine_and_coarsen_fixed_fraction(*tria,
 *                                                 error_indicators,
 *                                                 0.3,
 *                                                 0.05);
 *
 * // prepare the triangulation,
 * tria->prepare_coarsening_and_refinement();
 *
 * // prepare the SolutionTransfer object for coarsening and refinement and give
 * // the solution vector that we intend to interpolate later,
 * soltrans.prepare_for_coarsening_and_refinement(solution);
 *
 * // actually execute the refinement,
 * tria->execute_coarsening_and_refinement ();
 *
 * // redistribute dofs,
 * dof_handler->distribute_dofs (fe);
 *
 * // and interpolate the solution
 * Vector<double> interpolate_solution(dof_handler->n_dofs());
 * soltrans.interpolate(solution, interpolated_solution);
 * @endcode
 *
 * 如果网格被划分到几个MPI进程中，那么需要注意的是，旧的解决方案必须被复制到一个也能提供对本地相关DoF值（插值过程所需的这些值）的访问。
 *
 * @code
 * // Create initial indexsets pertaining to the grid before refinement
 * IndexSet locally_owned_dofs, locally_relevant_dofs;
 * locally_owned_dofs = dof_handler.locally_owned_dofs();
 * DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);
 *
 * // The solution vector only knows about locally owned DoFs
 * TrilinosWrappers::MPI::Vector solution;
 * solution.reinit(locally_owned_dofs,
 *               mpi_communicator);
 * ...
 * // Transfer solution to vector that provides access to locally relevant DoFs
 * TrilinosWrappers::MPI::Vector old_solution;
 * old_solution.reinit(locally_owned_dofs,
 *                   locally_relevant_dofs,
 *                   mpi_communicator);
 * old_solution = solution;
 * ...
 * // Refine grid
 * // Recreate locally_owned_dofs and locally_relevant_dofs index sets
 * ...
 * solution.reinit(locally_owned_dofs, mpi_communicator);
 * soltrans.refine_interpolate(old_solution, solution);
 * @endcode
 *
 * 对函数<code>interpolate (const VectorType &in, VectorType &out)</code>的多次调用是不允许的。通过使用<tt>void interpolate (const vector<VectorType> &all_in, vector<VectorType> &all_out) const</tt>，并使用各自的 @p  prepare_for_coarsening_and_refinement函数，在实际细化和粗化三角测量之前，将几个矢量作为输入，可以一步完成几个函数的插值（见这里）。  </ul>
 * 对于删除 @p SolutionTransfer
 * 中所有存储的数据并重新初始化，使用<tt>clear()</tt>函数。
 * 模板参数 @p VectorType 表示你要传输的数据容器的类型。
 *
 *  <h3>Interpolating in the presence of hanging nodes and boundary
 * values</h3>
 * 插值到新的网格是一个局部操作，也就是说，它只插值到新的网格。如果新的网格有悬空节点，你将得到一个不满足悬空节点约束的解。边界值也是如此：内插后的解将只是旧解在边界上的内插，而这可能满足也可能不满足新引入的边界节点的边界值。
 * 因此，你可能不得不在插值后应用悬挂节点或边界值约束。
 * step-15 和 step-26 有处理这个问题的例子。
 *
 *  <h3>Implementation</h3>
 * <ul>   <li>  只用细化的解转移。假设我们在当前（原始）网格上得到了一个解向量。这个向量的每个条目都属于离散化的一个DoF。如果我们现在细化网格，那么对 DoFHandler::distribute_dofs() 的调用将至少改变一些DoF指数。因此，我们需要在细化前存储所有活动单元的DoF指数。每个活动单元的指针被用来指向该单元的这些DoF指数的向量。这是由prepare_for_pure_refinement()完成的。
 * 在函数<tt>refine_interpolate(in,out)</tt>中，在每个设置了指针的单元上（即在原始网格中处于活动状态的单元），我们现在可以通过使用存储的DoF指数访问该单元上的解向量
 * @p in 的局部值。这些局部值被内插并设置到向量 @p out
 * 中，该向量是在细化网格上内插的离散函数 @p 的末端。
 * <tt>refine_interpolate(in,out)</tt>函数可以为原始网格上任意多的离散函数（解向量）多次调用。
 * <li>  带有粗化和细化的解转移。在调用
 * Triangulation::prepare_coarsening_and_refinement
 * 后，一个（父）单元的所有子单元的粗化标志或者没有被设置。在粗化
 * (Triangulation::execute_coarsening_and_refinement)
 * 时，不再需要的单元将从三角图中删除。
 * 对于从（将要粗化的）子单元到其父单元的插值，需要子单元。因此，在这些子单元被粗化（和删除！）之前，需要进行这种内插和存储我们想要内插的每个离散函数的内插值的工作。同样，每个相关单元的指针被设置为指向这些值（见下文）。此外，不被粗化的单元的DoF指数需要根据纯细化的解决方案转移来存储（参见这里）。所有这些都由<tt>prepare_for_coarsening_and_refinement(all_in)</tt>执行，其中<tt>vector<VectorType>
 * all_in</tt>包括所有要插值到新网格的离散函数。
 * 由于我们需要两种不同的指针（<tt>vector<unsigned
 * int></tt>用于Dof指数，<tt>vector<VectorType></tt>用于内插的DoF值），我们使用包括这两个指针的
 * @p Pointerstruct ，每个单元的指针指向这些 @p
 * 指针结构。在每个单元中，每次只使用两个不同的指针中的一个，因此我们可以在一个时候将<tt>void指针</tt>作为<tt>vector<unsigned
 * int></tt>使用，在另一个时候作为<tt>vector<VectorType></tt>使用，但是在两者之间使用这个
 * @p
 * Pointerstruct使得这些指针的使用更加安全，并提供更好的可能性来扩展它们的使用。
 * 在<tt>interpolate(all_in, all_out)</tt>中，细化的单元是根据纯细化时的解决方案转移来处理的。此外，在每个被粗化的单元上（因此以前是父单元）， @p all_out 中离散函数的值被设置为存储的局部插值，由于 @p Pointerstruct 中的'vector<VectorType>'指针被该单元的指针所指向，所以可以访问这些值。很明显，<tt>interpolate(all_in, all_out)</tt>只能用之前作为<tt>prepare_for_coarsening_and_refinement(all_in)</tt>函数参数的<tt>向量<VectorType> all_in</tt>来调用。因此，<tt>interpolate(all_in, all_out)</tt>可以（与<tt>refine_interpolate(in, out)</tt>相反）只被调用一次。  </ul>
 *
 *  <h3>Interaction with hanging nodes</h3>
 * 这个类尽力在新网格上表示旧网格上存在的有限元函数，但这可能导致新网格上的函数在悬挂节点处不再符合要求的情况。为此，考虑一个两次细化的网格的情况，开始时只有一个方形单元（即我们现在有16个单元）。再考虑到我们将其中的4个单元粗化到第一个细化水平。在这种情况下，如果我们使用
 * $Q_1$ 元素，我们最终得到的网格将如下所示。
*  @image html hanging_nodes.png ""
 * 从旧网格插值到新网格的过程将意味着有限元函数的值在所有保持原样的单元（即精细单元）上不会改变，但在右上方的粗单元上，顶点上的四个值是通过从其以前的子单元向下插值而得到。
 * 如果原始函数不是线性的，这意味着被标记的悬空节点将保留它们的旧值，一般来说，这不会导致沿相应边缘的连续函数。换句话说，在
 * SolutionTransfer::interpolate()
 * 之后得到的解向量不满足悬空节点约束：它对应于点式插值，但不对应于插值<i>onto
 * the new finite element space that contains constraints from hanging
 * nodes</i>。
 * 这是否是一个你需要担心的问题，取决于你的应用。当然，这种情况很容易纠正，在转移后将
 * AffineConstraints::distribute()
 * 应用于你的解决方案向量，使用在新的DoFHandler对象上计算的约束对象（如果你有悬挂节点，你可能需要创建这个对象）。例如，这也是在
 * step-15 中所做的。
 *
 *
 * @note
 * 这种情况只有在你做粗加工时才会发生。如果所有单元保持原样或被细化，那么
 * SolutionTransfer::interpolate()
 * 就会计算出一个新的节点值矢量，但所代表的函数当然是完全相同的，因为旧的有限元空间是新的有限元空间的一个子空间。因此，如果旧函数是符合要求的（即满足悬挂节点约束），那么新函数也是符合要求的，没有必要调用
 * AffineConstraints::distribute().  。
 *
 *  <h3>Implementation in the context of hp-finite elements</h3>
 * 在具有hp-capabilities的DoFHandler的情况下，没有任何东西定义哪些属于与DoFHandler相关的
 * hp::FECollection
 * 的有限元，应该被考虑在不活动的单元上（即，有子单元）。这是因为自由度只分配给活动单元，事实上，不允许使用
 * DoFAccessor::set_active_fe_index().
 * 在非活动单元上设置活动FE索引。
 * 因此，如果，例如，几个单元被粗化掉，应该发生什么是不完全自然的。然后这个类实现了以下算法。
 *
 *
 *
 * 如果一个单元被细化，那么在细化之前，解向量的值将从活动有限元的空间内插到未来有限元的空间内。然后，这些值被分配到细化后的子单元的有限元空间。如果旧单元使用Q2空间，而子单元使用Q1空间，这可能会丢失信息；如果母单元使用Q1空间，而子单元是Q2空间，信息可能会被延长。
 *
 *
 *
 * - 如果单元格要被粗化，那么子单元格的值将使用子单元格未来有限元空间中最大的空间内插到母单元格，这将被识别为遵循FiniteElementDomination逻辑的最小主导元素（参考 hp::FECollection::find_dominated_fe_extended() 了解更多信息）。例如，如果一个单元的子单元使用Q1、Q2和Q3空间，那么来自子单元的值将被内插到母单元的Q3空间中。在细化之后，母单元上的这个Q3函数会被内插到用户为这个单元选择的空间中（在这个例子中，如果用户在细化之后和调用 DoFHandler::distribute_dofs()). 之前为不同的空间设置了活动FE索引，那么这个空间可能与Q3不同。
 *
 *
 * @note
 * 在HP精简的背景下，如果单元被粗化或者某些单元上的多项式度数被降低，那么旧的有限元空间就不是新空间的子空间，你可能会遇到与上面讨论的悬挂节点相同的情况。你可能要考虑对通过转移解决方案得到的向量调用
 * AffineConstraints::distribute() 。
 *
 *
 * @ingroup numerics
 *
 *
 */
template <int dim, typename VectorType = Vector<double>, int spacedim = dim>
class SolutionTransfer
{
public:
  /**
   * 构造函数，将当前的DoFHandler作为参数。
   *
   */
  SolutionTransfer(const DoFHandler<dim, spacedim> &dof);

  /**
   * 解除函数
   *
   */
  ~SolutionTransfer();

  /**
   * 在调用构造函数后直接将该类重设为它的状态。
   *
   */
  void
  clear();

  /**
   * 为 @p SolutionTransfer
   * 的纯细化做准备。它存储了每个单元的dof指数。调用此函数后，只允许调用
   * @p  refine_interpolate函数。
   *
   */
  void
  prepare_for_pure_refinement();

  /**
   * 为 @p SolutionTransfer
   * 的粗化和细化做准备。它存储了每个单元格的道夫指数，并将
   * @p all_in
   * 中的向量的道夫值存储在每个将被粗化的单元格中。  @p
   * all_in
   * 包括所有将被内插到新（细化和/或粗化）网格的向量。
   *
   */
  void
  prepare_for_coarsening_and_refinement(const std::vector<VectorType> &all_in);

  /**
   * 与前面的函数相同，但只对一个离散函数进行内插。
   *
   */
  void
  prepare_for_coarsening_and_refinement(const VectorType &in);

  /**
   * 这个函数将细化前网格上的离散函数 @p in,
   * 插值到细化后的网格上的函数 @p out
   * ，后者是一个矢量。它假定向量具有正确的尺寸（即<tt>in.size()==n_dofs_old</tt>,
   * <tt>out.size()==n_dofs_refined</tt>）。    只有当 @p
   * prepare_for_pure_refinement
   * 被调用并且细化被执行之前，才允许调用此函数。允许多次调用此函数。例如，用于插值几个函数。
   *
   */
  void
  refine_interpolate(const VectorType &in, VectorType &out) const;

  /**
   * 这个函数将存储在 @p
   * all_in中的离散函数插值到细化和/或粗化的网格上。它假定
   * @p all_in 中的向量表示与 @p all_in
   * 中的向量相同，作为<tt>prepare_for_refinement_and_coarsening(all_in)/tt>的参数。然而，内部没有办法验证这一点，所以这里要小心。
   * 只有在第一次  Triangulation::prepare_coarsening_and_refinement,
   * 第二次  @p   SolutionTransfer::prepare_for_coarsening_and_refinement,
   * 然后第三次  Triangulation::execute_coarsening_and_refinement
   * 被调用之前，才允许调用此函数。
   * 不允许多次调用该函数。几个函数的插值可以在一个步骤中进行。
   * 输出向量的数量被假定为与输入向量的数量相同。另外，假设输出向量的大小是正确的（
   * @p n_dofs_refined).  否则将抛出一个断言。
   *
   */
  void
  interpolate(const std::vector<VectorType> &all_in,
              std::vector<VectorType> &      all_out) const;

  /**
   * 和前面的函数一样。它只对一个函数进行插值。它假定向量具有正确的大小（即<tt>in.size()==n_dofs_old</tt>,
   * <tt>out.size()==n_dofs_refined</tt>）。
   * 不允许多次调用此函数。通过使用<tt>interpolate (all_in,
   * all_out)</tt>，可以在一个步骤中执行多个函数的插值。
   *
   */
  void
  interpolate(const VectorType &in, VectorType &out) const;

  /**
   * 确定这个对象的内存消耗（以字节为单位）的估计值。
   *
   */
  std::size_t
  memory_consumption() const;

  /**
   * 异常情况
   *
   */
  DeclExceptionMsg(ExcNotPrepared,
                   "You are attempting an operation for which this object is "
                   "not prepared. This may be because you either did not call "
                   "one of the prepare_*() functions at all, or because you "
                   "called the wrong one for the operation you are currently "
                   "attempting.");

  /**
   * 异常情况
   *
   */
  DeclExceptionMsg(
    ExcAlreadyPrepForRef,
    "You are attempting to call one of the prepare_*() functions "
    "of this object to prepare it for an operation for which it "
    "is already prepared. Specifically, the object was "
    "previously prepared for pure refinement.");

  /**
   * 异常情况
   *
   */
  DeclExceptionMsg(
    ExcAlreadyPrepForCoarseAndRef,
    "You are attempting to call one of the prepare_*() functions "
    "of this object to prepare it for an operation for which it "
    "is already prepared. Specifically, the object was "
    "previously prepared for both coarsening and refinement.");

private:
  /**
   * 指向自由度处理程序的指针，以便与之合作。
   *
   */
  SmartPointer<const DoFHandler<dim, spacedim>,
               SolutionTransfer<dim, VectorType, spacedim>>
    dof_handler;

  /**
   * 存储细化和/或粗化前的自由度数量。
   *
   */
  types::global_dof_index n_dofs_old;

  /**
   * 对 @p PreparationState 的声明，表示 @p SolutionTransfer:
   * 的三种可能状态：准备进行 "纯细化"，准备进行
   * "粗化和细化 "或未准备。
   *
   */
  enum PreparationState
  {
    /**
     * 溶液转移器还没有准备好。
     *
     */
    none,
    /**
     * 该SolutionTransfer是为纯粹的细化准备的。
     *
     */
    pure_refinement,
    /**
     * SolutionTransfer是为粗化和细化而准备的。
     *
     */
    coarsening_and_refinement
  };

  /**
   * 各个变量的定义。
   *
   */
  PreparationState prepared_for;


  /**
   * 用于 @p prepare_for_refining （当然也用于 @p
   * repare_for_refining_and_coarsening），并存储所有将被细化的单元格的dof指数。
   *
   */
  std::vector<std::vector<types::global_dof_index>> indices_on_cell;

  /**
   * 所有的单元格数据（dof指数和dof值）都应该可以从每个单元格中访问。由于每个单元只有一个
   * @p user_pointer,
   * ，所以需要将数据的多个指针打包到一个结构中。请注意，在我们的案例中，每个单元都需要<tt>向量<unsigned
   * int>索引</tt>（如果该单元将被细化）或<tt>向量<double>
   * dof_values</tt>（如果该单元的孩子将被删除），因此一个
   * @p
   * user_pointer应该足够了，但为了允许一些错误检查和保护用户不犯用户错误，
   * @p user_pointer  将被这个结构 "乘以"。
   *
   */
  struct Pointerstruct
  {
    Pointerstruct()
      : indices_ptr(nullptr)
      , dof_values_ptr(nullptr)
      , active_fe_index(0)
    {}
    Pointerstruct(std::vector<types::global_dof_index> *indices_ptr_in,
                  const unsigned int                    active_fe_index_in = 0)
      : indices_ptr(indices_ptr_in)
      , dof_values_ptr(nullptr)
      , active_fe_index(active_fe_index_in)
    {}
    Pointerstruct(
      std::vector<Vector<typename VectorType::value_type>> *dof_values_ptr_in,
      const unsigned int active_fe_index_in = 0)
      : indices_ptr(nullptr)
      , dof_values_ptr(dof_values_ptr_in)
      , active_fe_index(active_fe_index_in)
    {}
    std::size_t
    memory_consumption() const;

    std::vector<types::global_dof_index> *                indices_ptr;
    std::vector<Vector<typename VectorType::value_type>> *dof_values_ptr;
    unsigned int                                          active_fe_index;
  };

  /**
   * 从单元格的级别和索引到 @p Pointerstructs
   * 的地图映射（参见那里）。这个映射使得在这个对象中保留所有传输解决方案所需的信息成为可能，而不是为此目的使用Triangulation的用户指针。
   *
   */
  std::map<std::pair<unsigned int, unsigned int>, Pointerstruct> cell_map;

  /**
   * 用于 @p prepare_for_refining_and_coarsening
   * 所有将被粗化的单元的内插dof值将被存储在这个向量中。
   *
   */
  std::vector<std::vector<Vector<typename VectorType::value_type>>>
    dof_values_on_cell;
};

namespace Legacy
{
  /**
   * @deprecated  使用没有DoFHandlerType模板的 dealii::SolutionTransfer
   * 代替。
   *
   */
  template <int dim,
            typename VectorType     = Vector<double>,
            typename DoFHandlerType = DoFHandler<dim>>
  using SolutionTransfer DEAL_II_DEPRECATED =
    dealii::SolutionTransfer<dim, VectorType, DoFHandlerType::space_dimension>;
} // namespace Legacy


DEAL_II_NAMESPACE_CLOSE

#endif // dealii_solutiontransfer_h
 /*---------------------------- solutiontransfer.h ---------------------------*/ 


