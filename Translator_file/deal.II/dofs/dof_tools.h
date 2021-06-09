//include/deal.II-translator/dofs/dof_tools_0.txt
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

#ifndef dealii_dof_tools_h
#define dealii_dof_tools_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/point.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/component_mask.h>

#include <deal.II/hp/dof_handler.h>

#include <deal.II/lac/affine_constraints.h>

#include <map>
#include <ostream>
#include <set>
#include <vector>


DEAL_II_NAMESPACE_OPEN

// Forward declarations
#ifndef DOXYGEN
class BlockMask;
template <int dim, typename RangeNumberType>
class Function;
template <int dim, int spacedim>
class FiniteElement;
namespace hp
{
  template <int dim, int spacedim>
  class MappingCollection;
  template <int dim, int spacedim>
  class FECollection;
} // namespace hp
template <class MeshType>
class InterGridMap;
template <int dim, int spacedim>
class Mapping;
class SparsityPattern;
template <int dim, class T>
class Table;
template <typename Number>
class Vector;

namespace GridTools
{
  template <typename CellIterator>
  struct PeriodicFacePair;
}

namespace DoFTools
{
  namespace internal
  {
    /* make_flux_sparsity_pattern的face_has_flux_coupling参数的默认值。在此定义（而不是使用参数列表中的默认lambda），以避免gcc中的一个bug，即同一个lambda被多次定义。   
*
*/
    template <int dim, int spacedim>
    inline bool
    always_couple_on_faces(
      const typename DoFHandler<dim, spacedim>::active_cell_iterator &,
      const unsigned int)
    {
      return true;
    }
  } // namespace internal
} // namespace DoFTools

#endif

/**
 * 这是一个对自由度进行操作的函数集合，并对自由度的数量进行操作。成员函数的文档将提供更多的信息，但是对于存在于多个版本的函数，在这个全局文档中有一些章节说明了一些共同点。
 * <h3>Setting up sparsity patterns</h3>
 * 当组装系统矩阵时，条目通常是 $a_{ij} = a(\phi_i, \phi_j)$
 * 的形式，其中 $a$
 * 是一个双线性函数，通常是一个积分。因此，当使用稀疏矩阵时，我们只需要为那些
 * $a_{ij}$ 保留空间，它们是非零的，这等于说基函数 $\phi_i$
 * 和 $\phi_j$
 * 的支持有一个非空交点。由于基函数的支持只约束在它们所在的单元或与之相邻的单元上，为了确定稀疏模式，只需在所有单元上循环，并将每个单元上的所有基函数与该单元上的所有其他基函数连接起来。
 * 可能有一些有限元，其单元上的所有基函数并不相互连接，但由于作者不知道有这种情况发生，所以没有使用这种情况。
 *
 *  <h3>DoF numberings on boundaries</h3>
 * 当把函数的轨迹投射到边界或部分边界时，我们需要建立只作用于位于边界上的自由度的矩阵和向量，而不是作用于所有自由度。我们可以通过简单地建立所有内部自由度的条目为零的矩阵来做到这一点，但是这样的矩阵总是有很大的等级缺陷，而且在工作中不是很实用。
 * 在这种情况下，需要的是对边界自由度进行编号，也就是说，我们应该列举所有位于边界的自由度，并排除所有其他（内部）自由度。map_dof_to_boundary_indices()函数正是这样做的：它提供一个矢量，其条目数与整个域上的自由度一样多，每个条目是边界编号中的数字，如果自由度不在边界上，则为
 * numbers::invalid_dof_index 。
 * 有了这个向量，对于任何给定的自由度，可以在那些位于边界上的自由度中得到一个唯一的数字；或者，如果你的自由度在域的内部，结果就是
 * numbers::invalid_dof_index.
 * 我们需要这个映射，例如，在边界上建立质量矩阵（关于这一点，见make_boundary_sparsity_pattern()函数，下面相应部分，以及MatrixCreator命名空间文档）。
 * 实际上，有两个map_dof_to_boundary_indices()函数，一个产生所有边界自由度的编号，另一个只产生边界部分的编号，即那些边界指标被列在给函数的指标集中的部分。后一种情况是需要的，例如，我们只想投射边界的Dirichlet部分的边界值。然后，你给函数一个边界指标的列表，指的是要进行投影的Dirichlet部分。你想投射的边界部分不需要是连续的；但是，不保证每个边界部分的指数是连续的，也就是说，不同部分的自由度指数可能是混合的。
 * 边界上的自由度但不在指定的边界部分之一的自由度被赋予指数
 * numbers::invalid_dof_index,
 * ，就像它们在内部一样。如果没有给出边界指标，或者如果一个单元格的任何面都没有包含在给定列表中的边界指标，那么新指数的向量仅由
 * numbers::invalid_dof_index. 组成。
 * （作为一个附带说明，对于角落的情况。边界上的自由度是什么，这个问题不是那么容易。
 * 它实际上应该是一个自由度，其各自的基础函数在边界上有非零值。至少对于拉格朗日元素来说，这个定义等于说形状函数的离点，或者deal.II所说的支持点，即函数承担其名义值的点（对于拉格朗日元素来说，这是函数值为1的点），是位于边界上。我们并不直接检查这一点，这个标准是通过有限元类给出的信息来定义的：有限元类定义了每个顶点、每条线等的基函数数量，基函数是根据这些信息来编号的；一个基函数被认为是在一个单元的面上（如果该单元在边界上，也就在边界上），根据它属于一个顶点、线等，但不属于单元的内部。有限元使用相同的单元编号，因此我们可以说，如果一个自由度被编号为线上的一个道夫，我们就认为它位于线上。偏离点究竟在哪里，是有限元的秘密（好吧，你可以问，但我们在这里不做），在此不作讨论。
 *
 *  <h3>Setting up sparsity patterns for boundary matrices</h3>
 * 在某些情况下，人们只想处理位于边界上的DoF。例如，一个应用是，如果不是插值非同质的边界值，而是想投影它们。为此，我们需要两样东西：一种识别位于边界上（部分）的节点的方法，以及一种只用边界上的自由度建立矩阵的方法（即更小的矩阵，在其中我们甚至不建立大的零块，因为大多数自由度在域的边界上没有支持）。这些任务中的第一部分由map_dof_to_boundary_indices()函数完成（如上所述）。
 * 第二部分要求我们首先为边界节点之间的耦合建立一个稀疏模式，然后实际建立这个矩阵的组成部分。虽然实际计算这些小边界矩阵的条目在MatrixCreator命名空间中讨论，但创建稀疏模式是由create_boundary_sparsity_pattern()函数完成的。对于它的工作，它需要有一个所有这些自由度的编号，这些自由度在我们感兴趣的边界的那些部分。你可以从map_dof_to_boundary_indices()函数中得到这个数字。然后它建立对应于
 * $\int_\Gamma \varphi_{b2d(i)} \varphi_{b2d(j)} dx$
 * 这样的积分的稀疏模式，其中 $i$ 和 $j$ 是矩阵的索引，
 * $b2d(i)$ 是位于边界上的自由度的全局DoF编号（即 $b2d$
 * 是map_dof_to_boundary_indices() 函数返回的映射的逆值）。
 *
 *
 *
 * @ingroup dofs
 *
 *
 */
namespace DoFTools
{
  /**
   * 某些<tt>make_*_pattern</tt>函数在表格中使用的标志，用于描述解的两个分量是否在对应于单元项或面项的双线性形式中耦合。一个使用这些标志的例子在
   * step-46  的介绍中显示。
   * 在下面对各个元素的描述中，请记住这些标志是作为大小为
   * FiniteElement::n_components 乘以 FiniteElement::n_components
   * 的表格的元素使用的，其中每个元素表示两个组件是否耦合。
   *
   */
  enum Coupling
  {
    /**
     * 两个组件不耦合。
     *
     */
    none,
    /**
     * 两个组件是耦合的。
     *
     */
    always,
    /**
     * 只有当两个组件的形状函数在一个给定的面上都不为零时，它们才会耦合。这个标志只在计算单元格面上的积分时使用，例如，在
     * DoFTools::make_flux_sparsity_pattern(). 中 使用 Coupling::always
     * 在一般情况下，梯度等发生在面上的积分。
     *
     */
    nonzero
  };

  /**
   * @name  DoF耦合  @{
   *
   */

  /**
   * 将一个耦合表从用户友好的按组件组织映射到按块组织。
   * 返回的向量将被初始化为该函数中的正确长度。
   *
   */
  template <int dim, int spacedim>
  void
  convert_couplings_to_blocks(const DoFHandler<dim, spacedim> &dof_handler,
                              const Table<2, Coupling> &table_by_component,
                              std::vector<Table<2, Coupling>> &tables_by_block);

  /**
   * 给定一个有限元和它的向量分量如何相互耦合的表，计算并返回一个描述各个形状函数如何相互耦合的表。
   *
   */
  template <int dim, int spacedim>
  Table<2, Coupling>
  dof_couplings_from_component_couplings(
    const FiniteElement<dim, spacedim> &fe,
    const Table<2, Coupling> &          component_couplings);

  /**
   * 与上述有限元集合的函数相同，返回一个表格的集合。
   * 该函数目前对 DoFTools::Couplings::nonzero 的处理与
   * DoFTools::Couplings::always 相同。
   *
   */
  template <int dim, int spacedim>
  std::vector<Table<2, Coupling>>
  dof_couplings_from_component_couplings(
    const hp::FECollection<dim, spacedim> &fe,
    const Table<2, Coupling> &             component_couplings);
  /**
   * @}
   *
   */

  /**
   * @name  稀疏模式生成  
     * @{ 
   *
   */

  /**
   * 计算建立在给定 @p dof_handler
   * 上的矩阵的哪些条目可能是非零的，并创建一个代表这些非零位置的稀疏模式对象。
   * 这个函数通过<i>simulating</i>计算全局系统矩阵中非零项的可能位置，在实际组装矩阵的过程中，人们会将这些条目写入全局系统矩阵。为此，该函数假设每个有限元基函数只有在其自由度与该单元的内部、面、边或顶点相关时，才是该单元的非零值。
   * 因此，从两个具有（全局）指数 $i$ 和 $j$ 的基函数
   * $\varphi_i$ 和 $\varphi_j$ 计算出来的矩阵条目 $A_{ij}$
   * （例如，使用双线性形式 $A_{ij}=a(\varphi_i,\varphi_j)$
   * ）只有在这些形状函数对应于至少一个共同单元上定义的自由度时才可能是非零。因此，这个函数只是在所有单元中循环，找出所有自由度的全局指数，并假定所有与这些指数相联系的矩阵条目将导致一个非零矩阵条目。然后，这些将被添加到稀疏模式中。
   * 由于这个生成稀疏性模式的过程没有考虑到以后要解决的方程，所以产生的稀疏性模式是对称的。
   * 这种算法对每个单元上的形状函数不加区分，也就是说，它只是将一个单元上的所有自由度与一个单元上的所有其他自由度进行耦合。这通常是一种情况，而且总是一种安全的假设。然而，如果你对运算符的结构有所了解，知道它不会将某些形状函数与某些测试函数耦合，那么你可以通过调用下面描述的当前函数的变体来获得更稀疏的稀疏模式，该变体允许指定哪些向量分量与其他向量分量耦合。
   * 上面描述的方法基于这样的假设：自由度之间的耦合只发生在形状函数至少在一个单元上重叠的情况下。这是最常见的涉及保形元素的有限元公式的情况。然而，对于诸如非连续Galerkin有限元方法这样的公式，双线性形式包含了单元之间的界面条款，这些条款将生活在一个单元上的形状函数与生活在相邻单元上的形状函数相耦合。当前函数不会看到这些耦合，因此不会在稀疏模式中分配条目。然后，你会在矩阵组装过程中遇到麻烦，因为你试图向矩阵条目中写入疏散模式中没有分配到的空间。
   * 这可以通过调用 DoFTools::make_flux_sparsity_pattern()
   * 函数来避免，该函数考虑了相邻单元上自由度之间的耦合。
   * 在其他情况下，双线性形式包含非局部项，例如在处理积分方程时。这些情况需要不同的方法来建立稀疏模式，这取决于问题的确切表述。那么，你必须自己做这件事。
   * @param[in]  dof_handler
   * 描述哪些自由度存在于哪些单元的DoFHandler对象。
   * @param[out]  sparsity_pattern 要填入条目的稀疏性模式。
   * @param[in]  约束
   * 上述生成条目的过程完全是每个单元的局部。因此，稀疏性模式没有规定只有在消除悬空节点或其他约束时才会被写入的矩阵条目。它们必须通过后续调用
   * AffineConstraints::condense().
   * 来处理。另外，在创建稀疏模式时，自由度的约束已经被考虑在内。为此，将AffineConstraints对象作为第三个参数传递给当前函数。这样就不需要调用
   * AffineConstraints::condense() 了。这个过程在 step-6 、 step-27
   * 和其他教程程序中都有解释。      @param[in]
   * keep_constrained_dofs
   * 如果约束条件已经通过传递AffineConstraints对象在该函数中得到了处理，那么如果这些条目在实际装配该疏散模式的矩阵时也不会被写入，那么就可以放弃疏散模式中的一些非对角线条目。具体来说，当使用
   * AffineConstraints::distribute_local_to_global(),
   * 的装配方法时，没有条目会被写入那些对应于受限自由度的矩阵行或列中。在这种情况下，你可以将参数
   * @p  keep_constrained_dofs设置为 @p false
   * 来避免在稀疏模式中分配这些条目。      @param[in]
   * subdomain_id
   * 如果指定的话，疏散模式只建立在子域_id等于给定参数的单元上。这在矩阵和稀疏模式（例如
   * TrilinosWrappers::SparsityPattern)
   * ）可能是分布式的，并且不是每个MPI进程都需要构建整个稀疏模式的并行环境中很有用；在这种情况下，如果每个进程只构建与它负责的子域_id相对应的那部分稀疏模式就足够了。这个特征在
   * step-32 中使用。（对于 parallel::distributed::Triangulation
   * 类型的对象通常不需要这个参数，因为当前函数无论如何只在本地拥有的单元上循环；因此，这个参数通常只在你想把子域_id用于指示哪个处理器拥有一个单元以外的事情时才有意义，例如一个单元属于域的哪个几何成分）。
   * @note  稀疏模式的实际类型可以是SparsityPattern、DynamicSparsityPattern、BlockSparsityPattern、BlockDynamicSparsityPattern，或者其他任何满足类似要求的类。假设疏散模式的大小与自由度的数量相匹配，并且如果疏散模式是 "静态 "
   * 的，则有足够的未使用的非零条目来填充疏散模式（关于这意味着什么的更多信息，请参见
   * @ref Sparsity
   * ）。这个函数生成的非零条目被添加到对象以前可能的内容中，也就是说，以前添加的条目不会被删除。
   * @note
   * 如果稀疏模式由SparsityPattern类型的对象表示（而不是例如DynamicSparsityPattern），你需要记住在生成模式后使用
   * SparsityPattern::compress() 。
   * @ingroup constraints
   *
   */
  template <int dim,
            int spacedim,
            typename SparsityPatternType,
            typename number = double>
  void
  make_sparsity_pattern(
    const DoFHandler<dim, spacedim> &dof_handler,
    SparsityPatternType &            sparsity_pattern,
    const AffineConstraints<number> &constraints = AffineConstraints<number>(),
    const bool                       keep_constrained_dofs = true,
    const types::subdomain_id subdomain_id = numbers::invalid_subdomain_id);

  /**
   * 计算建立在给定的 @p dof_handler
   * 上的矩阵的哪些条目可能是非零的，并创建一个代表这些非零位置的稀疏模式对象。
   * 这个函数是以前的make_sparsity_pattern()函数的一个简单的变化（关于所有常用参数的描述见那里），但它为矢量有限元提供了功能，允许更具体地确定哪些变量在哪个方程中耦合。
   * 例如，如果你想解决斯托克斯方程。
   * @f{align*}{
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   * -\Delta \mathbf u + \nabla p &= 0,\\ \text{div}\ u &= 0
   * @f}
   * *在两个空间维度上，使用稳定的Q2/Q1混合元素（使用FESystem类），那么你不希望所有的自由度在每个方程中耦合。更具体地说，在第一个方程中，只有 $u_x$ 和 $p$ 出现；在第二个方程中，只有 $u_y$ 和 $p$ 出现；而在第三个方程中，只有 $u_x$ 和 $u_y$ 出现。(注意，这个讨论只谈及解变量的矢量分量和不同的方程，而与自由度无关，事实上也与任何一种离散化无关)。我们可以用以下的 "耦合 "模式来描述这一点。    @f[
   * \left[
   * \begin{array}{ccc}
   * 1 & 0 & 1 \\
   * 0 & 1 & 1 \\
   * 1 & 1 & 0
   * \end{array}
   * \right]
   * @f] 其中 "1 "表示两个变量（即FES系统的矢量分量）在各自的方程中耦合，而 "0 "表示没有耦合。  这些零意味着在通过标准的有限元公式进行离散化时，我们将不会向矩阵中写入条目，例如，将压力测试函数与压力形状函数耦合（与上述其他零类似）。那么为矩阵中的这些条目和稀疏模式分配内存就是一种浪费，你可以通过创建一个像上面那样的掩码来避免这种情况，该掩码向计算稀疏模式的（当前）函数描述这一点。如上所述，上面显示的掩码是指组成FES系统的组件，而不是自由度或形状函数。    这个函数被设计成通过 @p couplings 参数接受耦合模式，如上图所示，该参数包含#Coupling类型的值。它就像前面的函数一样建立矩阵结构，但如果不是由耦合模式指定的，就不创建矩阵元素。如果耦合是对称的，那么产生的稀疏模式也将是对称的。    如果使用中的有限元的一些或全部形状函数在一个以上的分量中是非零的（用交易二的话说：它们是 @ref GlossPrimitive "非原始有限元"），就会有一个复杂的情况。  在这种情况下，采取对应于第一个非零分量的耦合元素，并忽略该分量的其他耦合元素。
   * @ingroup constraints
   *
   */
  template <int dim,
            int spacedim,
            typename SparsityPatternType,
            typename number = double>
  void
  make_sparsity_pattern(
    const DoFHandler<dim, spacedim> &dof_handler,
    const Table<2, Coupling> &       coupling,
    SparsityPatternType &            sparsity_pattern,
    const AffineConstraints<number> &constraints = AffineConstraints<number>(),
    const bool                       keep_constrained_dofs = true,
    const types::subdomain_id subdomain_id = numbers::invalid_subdomain_id);

  /**
   * 构建一个稀疏模式，允许在两个不同但相关的网格上耦合自由度。
   * 这个想法是，如果两个给定的DoFHandler对象对应于两个不同的网格（并且可能对应于在这些单元上使用的不同的有限元），但是如果它们所基于的两个三角形是通过分层细化从同一个粗网格导出的，那么人们可能会设置一个问题，想用一个网格的形状函数来测试另一个网格的形状函数。特别是，这意味着来自第一个网格上的单元的形状函数要与第二个网格上位于相应单元的形状函数进行测试；这种对应关系是IntergridMap类可以确定的。
   * 这个函数然后构建一个稀疏模式，其中代表行的自由度来自第一个给定的DoFHandler，而对应列的自由度则来自第二个DoFHandler。
   *
   */
  template <int dim, int spacedim, typename SparsityPatternType>
  void
  make_sparsity_pattern(const DoFHandler<dim, spacedim> &dof_row,
                        const DoFHandler<dim, spacedim> &dof_col,
                        SparsityPatternType &            sparsity);

  /**
   * 计算建立在给定 @p dof_handler
   * 上的矩阵的哪些条目可能是非零的，并创建一个代表这些非零位置的稀疏模式对象。这个函数是上述make_sparsity_pattern()函数的一个变体，它假定你想用来生成矩阵的双线性形式也包含了单元格之间<i>faces</i>的积分项（即，它包含了单元格之间的
   * "通量"，解释了这个函数的名称）。
   * 这个函数对非连续加尔金方法很有用，标准的make_sparsity_pattern()函数只会为一个单元上的所有自由度与同一单元上的所有其他自由度耦合创建非零条目；然而，在DG方法中，每个单元上的所有或部分自由度也与通过共同面连接到当前单元的其他单元的自由度耦合。当前函数还创建了由这些额外耦合产生的矩阵中的非零条目。换句话说，与make_sparsity_pattern()所做的工作相比，这个函数计算了一个严格的非零项的超级集合。
   * @param[in]  dof_handler
   * 描述哪些自由度存在于哪些单元上的DoFHandler对象。
   * @param[out]  sparsity_pattern 要填入的稀疏度模式。
   * @note  稀疏度模式的实际类型可以是SparsityPattern、DynamicSparsityPattern、BlockSparsityPattern、BlockDynamicSparsityPattern或任何其他满足类似要求的类。假设疏散模式的大小与自由度的数量相匹配，并且如果疏散模式是 "静态 "
   * 的，则有足够的未使用的非零条目来填充疏散模式（关于这意味着什么的更多信息，请参见
   * @ref Sparsity
   * ）。这个函数生成的非零条目被添加到对象以前可能的内容中，也就是说，以前添加的条目不会被删除。
   * @note
   * 如果稀疏模式由SparsityPattern类型的对象表示（而不是例如DynamicSparsityPattern），你需要记住在生成模式后使用
   * SparsityPattern::compress() 。
   * @ingroup constraints
   *
   */
  template <int dim, int spacedim, typename SparsityPatternType>
  void
  make_flux_sparsity_pattern(const DoFHandler<dim, spacedim> &dof_handler,
                             SparsityPatternType &            sparsity_pattern);

  /**
   * 这个函数的作用与其他make_flux_sparsity_pattern()函数基本相同，但允许指定一些额外的参数。这些参数的含义与上述第一个make_sparsity_pattern()函数中讨论的相同。
   * @ingroup constraints
   *
   */
  template <int dim,
            int spacedim,
            typename SparsityPatternType,
            typename number>
  void
  make_flux_sparsity_pattern(
    const DoFHandler<dim, spacedim> &dof_handler,
    SparsityPatternType &            sparsity_pattern,
    const AffineConstraints<number> &constraints,
    const bool                       keep_constrained_dofs = true,
    const types::subdomain_id subdomain_id = numbers::invalid_subdomain_id);


  /**
   * 这个函数的作用与另一个make_flux_sparsity_pattern()函数基本相同，但允许指定耦合矩阵，说明在你离散化的每个方程中，解变量的哪些成分是耦合的。这与上面第二个make_sparsity_pattern()函数中讨论的完全类似。
   * 事实上，这个函数需要两个这样的掩码，一个是描述哪些变量在构成双线性方程的单元积分中相互耦合，哪些变量在面积分中相互耦合。如果你把只由1组成的掩码传递给这两个掩码，那么你将得到与你调用上述make_sparsity_pattern()函数中的第一个相同的稀疏性模式。通过将这些掩码中的一些条目设置为0，你可以得到一个更稀疏的稀疏模式。
   * @ingroup constraints
   *
   */
  template <int dim, int spacedim, typename SparsityPatternType>
  void
  make_flux_sparsity_pattern(
    const DoFHandler<dim, spacedim> &dof,
    SparsityPatternType &            sparsity,
    const Table<2, Coupling> &       cell_integrals_mask,
    const Table<2, Coupling> &       face_integrals_mask,
    const types::subdomain_id subdomain_id = numbers::invalid_subdomain_id);


  /**
   * 这个函数与之前的make_flux_sparsity_pattern()函数的功能基本相同，但允许应用AffineConstraints对象。这对于有限元的某些部分是连续的，某些部分是不连续的情况很有用，允许对连续部分施加约束，同时也建立不连续部分所需的通量项。
   * 可选的 @p face_has_flux_coupling
   * 可以用来指定在哪些面上发生通量耦合。这允许在使用双线性形式时创建一个更稀疏的模式，即通量项只出现在三角结构中的一个子集的面上。默认情况下，通量耦合被添加到所有内部面。
   * @p face_has_flux_coupling
   * 应该是一个函数，它接收一个active_cell_iterator和一个面的索引，如果该面有一个通量耦合，应该返回true。当使用
   * ::dealii::DoFHandler 时，我们可以，比如说，使用
   * @code
   * auto face_has_flux_coupling =
   *  [](const typename DoFHandler<dim>::active_cell_iterator &cell,
   *     const unsigned int                                    face_index) {
   *    const Point<dim> &face_center = cell->face(face_index)->center();
   *    return 0 < face_center[0];
   *  };
   * @endcode
   *
   *
   */
  template <int dim,
            int spacedim,
            typename SparsityPatternType,
            typename number>
  void
  make_flux_sparsity_pattern(
    const DoFHandler<dim, spacedim> &dof,
    SparsityPatternType &            sparsity,
    const AffineConstraints<number> &constraints,
    const bool                       keep_constrained_dofs,
    const Table<2, Coupling> &       couplings,
    const Table<2, Coupling> &       face_couplings,
    const types::subdomain_id        subdomain_id,
    const std::function<
      bool(const typename DoFHandler<dim, spacedim>::active_cell_iterator &,
           const unsigned int)> &face_has_flux_coupling =
      &internal::always_couple_on_faces<dim, spacedim>);

  /**
   * 创建边界矩阵的稀疏模式。更多信息请参见该类的一般文档。
   * 该函数基本上做了其他make_sparsity_pattern()函数所做的事情，但假定用于建立矩阵的双线性形式不包括域的积分，而只包括域的边界上的积分。
   *
   */
  template <int dim, int spacedim, typename SparsityPatternType>
  void
  make_boundary_sparsity_pattern(
    const DoFHandler<dim, spacedim> &           dof,
    const std::vector<types::global_dof_index> &dof_to_boundary_mapping,
    SparsityPatternType &                       sparsity_pattern);

  /**
   * 这个函数是之前make_boundary_sparsity_pattern()函数的一个变体，我们假设将产生矩阵的边界积分只延伸到边界的那些部分，这些部分的边界指标在这个函数的
   * @p boundary_ids 参数中列出。
   * 这个函数本来可以通过传递一个 @p set
   * 的边界_id数字来写。然而，整个deal.II中处理边界指标的大多数函数都采取边界指标和相应的边界函数的映射，即一个
   * std::map<types::boundary_id,  const
   * Function<spacedim,number>*>参数。相应地，这个函数也是这样做的，尽管实际的边界函数在这里被忽略了。
   * 因此，如果你没有任何这样的边界函数，只要用你想要的边界指标创建一个地图，并将函数指针设置为空指针）。
   *
   */
  template <int dim,
            int spacedim,
            typename SparsityPatternType,
            typename number>
  void
  make_boundary_sparsity_pattern(
    const DoFHandler<dim, spacedim> &dof,
    const std::map<types::boundary_id, const Function<spacedim, number> *>
      &                                         boundary_ids,
    const std::vector<types::global_dof_index> &dof_to_boundary_mapping,
    SparsityPatternType &                       sparsity);

  /**
   * @}
   *
   */

  /**
   * @name  悬空节点和其他约束  
     * @{ 
   *
   */

  /**
   * @ingroup constraints
   *
   */
  template <int dim, int spacedim, typename number>
  void
  make_hanging_node_constraints(const DoFHandler<dim, spacedim> &dof_handler,
                                AffineConstraints<number> &      constraints);

  /**
   * 当问题中的不同变量被离散在不同的网格上时，这个函数被使用，其中一个网格严格地比另一个更粗。一个例子是优化问题，控制变量通常在比状态变量更粗的网格上离散。
   * 该函数的结果可以在数学上表述如下。让 ${\cal T}_0$ 和
   * ${\cal T}_1$ 是两个网格，其中 ${\cal T}_1$ 是由 ${\cal T}_0$
   * 严格通过细化或不考虑 ${\cal T}_0$
   * 的单元而得到的。在两者上使用相同的有限元，有与这些网格相关的函数空间
   * ${\cal V}_0$ 和 ${\cal V}_1$ 。那么每个函数 $v_0 \in {\cal V}_0$
   * 当然也可以在 ${\cal V}_1$ 中精确表示，因为通过构造
   * ${\cal V}_0 \subset {\cal V}_1$  。  然而，并不是 ${\cal V}_1$
   * 中的每个函数都能被表达为 ${\cal V}_0$
   * 的形状函数的线性组合。可以表示的函数位于 ${\cal V}_1$
   * 的一个同质子空间（即 ${\cal V}_0$
   * ，当然），这个子空间可以用 $CV=0$
   * 形式的线性约束来表示，其中 $V$ 是函数 $v\in {\cal V}_1$
   * 的结点值的向量。换句话说，每个同时满足  $v_h\in {\cal
   * V}_0$  的函数  $v_h=\sum_j V_j \varphi_j^{(1)} \in {\cal V}_1$
   * 都自动满足  $CV=0$
   * 。这个函数以AffineConstraints对象的形式计算矩阵 $C$ 。
   * 这些约束的构造如下：对于粗网格上的每个自由度（即形状函数），我们计算它在细网格上的表现，即细网格上的形状函数的线性组合看起来如何与粗网格上的形状函数相像。从这些信息中，我们可以计算出如果一个线性方程在细网格上的解在粗网格上可以表示出来，那么这些约束就必须成立。如何计算这些约束条件的确切算法相当复杂，最好是通过阅读源代码来理解，其中包含许多注释。
   * 这个函数的用法如下：它接受两个DoF处理程序作为参数，第一个是指粗网格，第二个是指细网格。在这两者上，一个有限元由DoF处理程序对象表示，通常会有几个向量分量，可能属于不同的基元。因此，这个函数的第二和第四个参数说明粗网格上的哪个矢量分量应被用来限制细网格上的所述分量。两个网格上的各个分量所使用的有限元必须是相同的。一个例子可以说明这一点：考虑一个优化问题，控制
   * $q$ 在粗网格上离散，状态变量 $u$
   * （和相应的拉格朗日乘子 $\lambda$
   * ）在细网格上离散。它们分别使用片状常数不连续、连续线性和连续线性元素进行离散。在粗网格上只有参数
   * $q$
   * 被表示，因此粗网格上的DoFHandler对象只表示一个变量，用片状常数不连续元素离散。那么，粗网格上表示矢量分量的参数将是零（唯一可能的选择，因为粗网格上的变量是标量）。如果细网格FES系统中变量的排序是
   * $u, q, \lambda$
   * ，那么对应于矢量分量的函数的第四个参数将是1（对应于变量
   * $q$ ；0将是 $u$ ，2将是 $\lambda$ ）。
   * 该函数还需要一个IntergridMap类型的对象，代表如何从粗网格单元到细网格上的相应单元。原则上，这个对象可以由函数本身从两个DoFHandler对象中生成，但由于它在使用不同网格的程序中可能是可用的，所以函数只是把它作为一个参数。
   * 计算出的约束被输入一个AffineConstraints类型的变量中；之前的内容不会被删除。
   *
   */
  template <int dim, int spacedim>
  void
  compute_intergrid_constraints(
    const DoFHandler<dim, spacedim> &              coarse_grid,
    const unsigned int                             coarse_component,
    const DoFHandler<dim, spacedim> &              fine_grid,
    const unsigned int                             fine_component,
    const InterGridMap<DoFHandler<dim, spacedim>> &coarse_to_fine_grid_map,
    AffineConstraints<double> &                    constraints);


  /**
   * 这个函数生成一个矩阵，当一个元素数与粗网格上该分量的自由度相同的数据向量乘以这个矩阵时，我们得到一个元素数与细网格上全局自由度相同的向量。细网格上有限元场的其他矢量分量的所有元素都不被触及。
   * 细网格的三角化可以是分布式的。当并行调用时，每个进程必须有一份粗网格的副本。在这种情况下，函数返回本地拥有的一组单元的转移表示。
   * 这个函数的输出是一种压缩格式，可以用来构造相应的稀疏转移矩阵。
   *
   */
  template <int dim, int spacedim>
  void
  compute_intergrid_transfer_representation(
    const DoFHandler<dim, spacedim> &              coarse_grid,
    const unsigned int                             coarse_component,
    const DoFHandler<dim, spacedim> &              fine_grid,
    const unsigned int                             fine_component,
    const InterGridMap<DoFHandler<dim, spacedim>> &coarse_to_fine_grid_map,
    std::vector<std::map<types::global_dof_index, float>>
      &transfer_representation);

  /**
   * @}
   *
   */


  /**
   * @name  周期性的边界条件  
     * @{ 
   *
   */

  /**
   * 将周期性边界条件引起的（代数）约束插入AffineConstraints对象 @p
   * constraints.  给定一对不一定活动的边界面 @p face_1 和 @p
   * face_2，这个函数将与 @p face_1
   * 描述的边界相关的所有DoF约束到 @p face_2.
   * 描述的边界的各自DoF 更确切地说。    如果 @p face_1 和 @p
   * face_2 都是活动面，它将 @p 面_1的DoFs添加到 @p constraints
   * 的约束DoFs列表中，并添加条目将它们约束到 @p
   * 面_2的相应DoFs值。这发生在一个纯粹的代数层面上，意味着
   * @p face_1
   * 上具有（局部面）索引<tt>i</tt>的全局DoF被约束到 @p
   * face_2
   * 上具有（局部面）索引<tt>i</tt>的DoF（可能被纠正方向，见下文）。
   * 否则，如果 @p face_1 和 @p face_2
   * 不是活动面，这个函数就会递归地在 @p face_1 和 @p face_2.
   * 的子代上循环，如果两个面中只有一个是活动的，那么我们就会递归地迭代非活动面的子代，并确保精炼面的解函数与非精炼面的解函数相等，这与我们在不同的精炼单元聚集的地方执行悬挂节点约束一样的方式。然而，与悬挂节点不同的是，我们并不强制要求域的两边只能有一个细化级别的差异，你希望是周期性的）。
   * 这个程序只约束那些还没有被约束的DoF。如果这个例程遇到一个已经被约束的DoF（例如被Dirichlet边界条件约束），约束的旧设置（条目被约束的DoF，不均匀性）被保留，什么也不会发生。
   * @p component_mask 中的标志（见 @ref GlossComponentMask
   * ）表示有限元空间的哪些部分应受到周期性边界条件的约束。如果它与默认值相同，则所有分量都受到约束。如果它与默认值不同，则假定条目数等于有限元的分量数。这可以用来强制执行方程组中只有一个变量的周期性。
   * @p face_orientation,   @p face_flip 和 @p face_rotation
   * 描述了在匹配和约束DoF之前应该应用到 @p face_1
   * 的方向。这与给定面在各自单元中的实际方向无关（对于边界面总是默认的），而是你想看到周期性被强制执行的方式。例如，通过使用这些标志，你可以执行
   * $u(0,y)=u(1,1-y)$
   * 那种条件（即莫比乌斯带），或者在三维中执行扭曲的环。更确切地说，这些标志是以如下方式匹配局部面的DoF指数。
   * 在2d中。<tt>face_orientation</tt>必须总是<tt>true</tt>，<tt>face_rotation</tt>总是<tt>false</tt>，face_flip具有<tt>line_flip</tt>的意义；这意味着，例如对于<tt>Q1</tt>。
   * @code
   *
   * face_orientation = true, face_flip = false, face_rotation = false:
   *
   *   face1:           face2:
   *
   *   1                1
   *   |        <-->    |
   *   0                0
   *
   *   Resulting constraints: 0 <-> 0, 1 <-> 1
   *
   *   (Numbers denote local face DoF indices.)
   *
   *
   * face_orientation = true, face_flip = true, face_rotation = false:
   *
   *   face1:           face2:
   *
   *   0                1
   *   |        <-->    |
   *   1                0
   *
   *   Resulting constraints: 1 <-> 0, 0 <-> 1
   * @endcode
   * 同样地，对于Q1在3D中的情况也是如此。
   * @code
   *
   * face_orientation = true, face_flip = false, face_rotation = false:
   *
   *   face1:           face2:
   *
   *   2
   *
   * - 3            2
   *
   * - 3
   *   |   |    <-->    |   |
   *   0
   *
   * - 1            0
   *
   * - 1
   *
   *   Resulting constraints: 0 <-> 0, 1 <-> 1, 2 <-> 2, 3 <-> 3
   *
   *   (Numbers denote local face DoF indices.)
   *
   *
   * face_orientation = false, face_flip = false, face_rotation = false:
   *
   *   face1:           face2:
   *
   *   1
   *
   * - 3            2
   *
   * - 3
   *   |   |    <-->    |   |
   *   0
   *
   * - 2            0
   *
   * - 1
   *
   *   Resulting constraints: 0 <-> 0, 2 <-> 1, 1 <-> 2, 3 <-> 3
   *
   *
   * face_orientation = true, face_flip = true, face_rotation = false:
   *
   *   face1:           face2:
   *
   *   1
   *
   * - 0            2
   *
   * - 3
   *   |   |    <-->    |   |
   *   3
   *
   * - 2            0
   *
   * - 1
   *
   *   Resulting constraints: 3 <-> 0, 2 <-> 1, 1 <-> 2, 0 <-> 3
   *
   *
   * face_orientation = true, face_flip = false, face_rotation = true
   *
   *   face1:           face2:
   *
   *   0
   *
   * - 2            2
   *
   * - 3
   *   |   |    <-->    |   |
   *   1
   *
   * - 3            0
   *
   * - 1
   *
   *   Resulting constraints: 1 <-> 0, 3 <-> 1, 0 <-> 2, 2 <-> 3
   *
   * and any combination of that...
   * @endcode
   * 可以指定一个矩阵 @p matrix 和一个 std::vector  @p first_vector_components，描述在约束到 @p face_2. 的DoF之前， @p face_1的DoF应该如何被修改，这里有两种声明。如果 std::vector   @p  first_vector_components是非空的，该矩阵被解释为一个 @p  dim  $\times$   @p dim  旋转矩阵，它被应用于FESystem的 @p first_vector_components 中列出的所有矢量值块。如果 @p first_vector_components为空，该矩阵将被解释为大小为no_face_dofs  $\times$  no_face_dofs的插值矩阵。    这个函数确保身份约束不会在 @p constraints.  @p periodicity_factor 中产生循环，可以用来实现 $\psi(\mathbf{r})=e^{-i\mathbf{k}\cdot\mathbf{r}}u(\mathbf{r})$ 形式的布洛赫周期条件（又称相移周期条件），其中 $u$ 是与晶格相同的周期性， $\mathbf{k}$ 是波矢，见[https://en.wikipedia.org/wiki/Bloch_wave]（https://en.wikipedia.org/wiki/Bloch_wave）。  在 @p face_2 处的解等于 @p face_1 乘以 @p periodicity_factor. 处的解。例如，如果 @p face_1 处的解是 $\psi(0)$ ， $\mathbf{d}$ 是 @p face_2, 上的相应点，那么 @p face_2 处的解应该是 $\psi(d) = \psi(0)e^{-i \mathbf{k}\cdot \mathbf{d}}$  。这个条件可以用  $\mathrm{periodicity\_factor}=e^{-i \mathbf{k}\cdot \mathbf{d}}$  来实现。    详细信息可参见  @ref GlossPeriodicConstraints  "关于周期性边界条件的词汇条目"
   * 。
   *
   */
  template <typename FaceIterator, typename number>
  void
  make_periodicity_constraints(
    const FaceIterator &                         face_1,
    const typename identity<FaceIterator>::type &face_2,
    AffineConstraints<number> &                  constraints,
    const ComponentMask &            component_mask   = ComponentMask(),
    const bool                       face_orientation = true,
    const bool                       face_flip        = false,
    const bool                       face_rotation    = false,
    const FullMatrix<double> &       matrix           = FullMatrix<double>(),
    const std::vector<unsigned int> &first_vector_components =
      std::vector<unsigned int>(),
    const number periodicity_factor = 1.);



  /**
   * 将周期性边界条件引起的（代数）约束插入AffineConstraints对象
   * @p constraints.
   * 这是上述make_periodicity_constraints()的低级变量的主要高级接口。它接受一个
   * std::vector   @p periodic_faces
   * 作为参数，并在每个条目上应用上面的make_periodicity_constraints()。
   * @p periodic_faces  可以通过 GridTools::collect_periodic_faces.
   * 创建。
   * @note  对于建立在 parallel::distributed::Triangulation 对象上的DoFHandler对象， parallel::distributed::Triangulation::add_periodicity 在调用此函数之前必须被调用。      @see   @ref GlossPeriodicConstraints  "关于周期性边界条件的词汇条目 "
   * 和 step-45 的进一步信息。
   *
   */
  template <int dim, int spacedim, typename number>
  void
  make_periodicity_constraints(
    const std::vector<GridTools::PeriodicFacePair<
      typename DoFHandler<dim, spacedim>::cell_iterator>> &periodic_faces,
    AffineConstraints<number> &                            constraints,
    const ComponentMask &            component_mask = ComponentMask(),
    const std::vector<unsigned int> &first_vector_components =
      std::vector<unsigned int>(),
    const number periodicity_factor = 1.);

  /**
   * 与上述相同。      @deprecated
   * 使用以dim和spacedim为模板参数的函数。
   *
   */
  template <typename DoFHandlerType, typename number>
  DEAL_II_DEPRECATED void
  make_periodicity_constraints(
    const std::vector<
      GridTools::PeriodicFacePair<typename DoFHandlerType::cell_iterator>>
      &                              periodic_faces,
    AffineConstraints<number> &      constraints,
    const ComponentMask &            component_mask = ComponentMask(),
    const std::vector<unsigned int> &first_vector_components =
      std::vector<unsigned int>(),
    const number periodicity_factor = 1.);



  /**
   * 将周期性边界条件引起的（代数）约束插入AffineConstraints
   * @p constraints.
   * 这个函数作为make_periodicity_constraints()函数的一个高级接口。
   * 定义一个 "第一 "边界为所有边界面，其边界ID为 @p
   * b_id1，一个 "第二 "边界由所有属于 @p  b_id2的面组成。
   * 这个函数试图在orthogonal_equality()的帮助下将所有属于第一条边界的面与属于第二条边界的面相匹配。更确切地说，坐标只在
   * @p direction 部分不同的面被识别。
   * 如果这个匹配是成功的，它将所有与 "第一
   * "边界相关的DoFs约束到 "第二
   * "边界的相应DoFs，并尊重两个面的相对方向。
   * @note  这个函数是一个方便的封装器。它在内部调用 GridTools::collect_periodic_faces() 的参数，并将输出结果反馈给上述make_periodicity_constraints()变量。如果你需要更多的功能，直接使用 GridTools::collect_periodic_faces() 。      @see   @ref GlossPeriodicConstraints  "关于周期性边界条件的词汇条目"
   * ，以获得更多信息。
   *
   */
  template <int dim, int spacedim, typename number>
  void
  make_periodicity_constraints(
    const DoFHandler<dim, spacedim> &dof_handler,
    const types::boundary_id         b_id1,
    const types::boundary_id         b_id2,
    const unsigned int               direction,
    AffineConstraints<number> &      constraints,
    const ComponentMask &            component_mask     = ComponentMask(),
    const number                     periodicity_factor = 1.);



  /**
   * @note  这个版本的make_periodicity_constraints将不会在单元格不在 @ref GlossFaceOrientation "标准方向 "
   * 的网格上工作。
   * @note  这个函数是一个方便的封装器。它在内部调用 GridTools::collect_periodic_faces() 中提供的参数，并将输出结果反馈给上述make_periodicity_constraints()变量。如果你需要更多的功能，请直接使用 GridTools::collect_periodic_faces() 。      @see   @ref GlossPeriodicConstraints  "关于周期性边界条件的词汇条目"
   * ，以获得更多信息。
   *
   */
  template <int dim, int spacedim, typename number>
  void
  make_periodicity_constraints(
    const DoFHandler<dim, spacedim> &dof_handler,
    const types::boundary_id         b_id,
    const unsigned int               direction,
    AffineConstraints<number> &      constraints,
    const ComponentMask &            component_mask     = ComponentMask(),
    const number                     periodicity_factor = 1.);

  /**
   * @}
   *
   */

  /**
   * @name  识别具有特殊性质的自由度子集  
     * @{ 
   *
   */

  /**
   * 返回一个IndexSet，描述所有将被接口约束的自由度，即所有悬挂节点。
   * 在 parallel::shared::Triangulation 或
   * parallel::distributed::Triangulation
   * 的情况下，只考虑本地相关的道夫。
   *
   */
  template <int dim, int spacedim>
  IndexSet
  extract_hanging_node_dofs(const DoFHandler<dim, spacedim> &dof_handler);

  /**
   * 提取属于矢量值有限元的某些矢量分量的自由度的（本地拥有的）指数。 @p component_mask 定义了要从DoFHandler @p dof. 中提取FES系统或矢量值元素的哪些组件或块，然后输出对象中的条目对应于属于这些组件的自由度。    如果所考虑的有限元不是原始的，即它的一些或全部形状函数在一个以上的矢量分量中是非零的（例如，对于FE_Nedelec或FE_RaviartThomas元素来说，这一点是成立的），那么形状函数不能与单个矢量分量相关联。  在这种情况下，如果 <em> 这个元素的一个 </em> 形状向量分量在 @p component_mask 中被标记（见 @ref GlossComponentMask ），那么这相当于选择 <em> 与这个非原始基础元素对应的所有 </em> 向量分量。      @param[in]  dof_handler 其列举的自由度将被此函数过滤的DoFHandler。    @param[in]  component_mask 一个说明你要选择哪些组件的掩码。该掩码的大小必须与 @p dof_handler. 所使用的FiniteElement中的构件数量相匹配。更多信息请参见 @ref GlossComponentMask "构件掩码的词汇表条目"
   * 。    @return
   * 一个IndexSet对象，它将确切地包含那些(i)对应于上述掩码所选择的自由度的条目，以及(ii)是本地拥有的条目。索引集的大小等于全局自由度的数量。请注意，产生的对象总是
   * DoFHandler::locally_owned_dofs() 返回的一个子集。
   *
   */
  template <int dim, int spacedim>
  IndexSet
  extract_dofs(const DoFHandler<dim, spacedim> &dof_handler,
               const ComponentMask &            component_mask);

  /**
   * 这个函数等同于上面的 DoFTools::extract_dofs() 函数，除了不是根据组件（见 @ref GlossComponent ）而是根据它们是否是特定块的一部分（见 @ref GlossBlock ）来选择提取哪些自由度。  因此，第二个参数不是一个ComponentMask，而是一个BlockMask对象。      @param[in]  dof_handler 其列举的自由度将被这个函数过滤的DoFHandler。    @param[in]  block_mask 一个说明你想选择哪些块的掩码。该掩码的大小必须与 @p dof_handler. 所使用的FiniteElement中的块数相匹配。更多信息请参见 @ref GlossBlockMask "关于块掩码的词汇表条目"
   * 。    @return
   * 一个IndexSet对象，它将确切地包含那些(i)对应于由上述掩码选择的自由度的条目，以及(ii)本地拥有的条目。索引集的大小等于全局自由度的数量。请注意，产生的对象总是
   * DoFHandler::locally_owned_dofs() 返回的一个子集。
   *
   */
  template <int dim, int spacedim>
  IndexSet
  extract_dofs(const DoFHandler<dim, spacedim> &dof_handler,
               const BlockMask &                block_mask);

  /**
   * 对于多网格自由度编号的一个层次，做与相应的extract_dofs()函数相同的事情。
   *
   */
  template <int dim, int spacedim>
  void
  extract_level_dofs(const unsigned int               level,
                     const DoFHandler<dim, spacedim> &dof,
                     const ComponentMask &            component_mask,
                     std::vector<bool> &              selected_dofs);

  /**
   * 对一个多网格DoF编号的一个层次，做与相应的extract_dofs()函数相同的事情。
   *
   */
  template <int dim, int spacedim>
  void
  extract_level_dofs(const unsigned int               level,
                     const DoFHandler<dim, spacedim> &dof,
                     const BlockMask &                component_mask,
                     std::vector<bool> &              selected_dofs);

  /**
   * 提取所有在边界的自由度，并属于解决方案的指定组件。该函数在最后一个非缺省值参数中返回其结果，如果一个自由度在边界并属于所选分量之一，该参数包含 @p true ，否则包含 @p false 。    通过指定 @p boundary_ids 变量，你可以选择自由度所在的面必须有哪些边界指标才能被提取出来。如果它是一个空列表，那么所有的边界指标都被接受。     @p component_mask 的大小（见 @ref GlossComponentMask ）应等于 @p 自由度所使用的有限元中的组件数。 @p selected_dofs 的大小应等于<tt>dof_handler.n_dofs()</tt>。这个数组以前的内容会被覆盖掉。    使用通常的惯例，如果一个形状函数在一个以上的分量中是非零的（即它是非正则的），那么就使用分量掩码中对应于第一个非零分量的元素。  掩码中对应于后面分量的元素被忽略。      @deprecated  这个函数对建立在 parallel::distributed::Triangulation 对象上的DoFHandler对象不起作用。原因是输出参数 @p selected_dofs 的长度必须等于<i>all</i>全局自由度。因此，这不能扩展到非常大的问题，这也是该函数被废弃的原因。如果你需要这个函数的功能来进行平行三角计算，那么你需要使用另一个 DoFTools::extract_boundary_dofs() 函数，它通过IndexSet对象返回信息。      @param[in]  dof_handler 描述哪个自由度在哪个单元上的对象。    @param[in]  component_mask 表示应考虑的有限元的矢量分量的掩码（也见 @ref GlossComponentMask ）。    @param[out]  selected_dofs 一个被返回的布尔运算的向量，对于这个向量，如果相应的索引是位于边界上的自由度（并且对应于被选择的向量分量和边界指标，取决于 @p component_mask 和 @p boundary_ids参数的值），则该元素将是 @p true 。    @param[in]  boundary_ids 如果为空，该函数提取边界所有部分的自由度指数。如果它是一个非空的列表，那么这个函数只考虑具有这个参数中所列边界指标的边界面。      @see   @ref GlossBoundaryIndicator  "关于边界指标的词汇条目"
   *
   */
  template <int dim, int spacedim>
  DEAL_II_DEPRECATED void
  extract_boundary_dofs(const DoFHandler<dim, spacedim> &   dof_handler,
                        const ComponentMask &               component_mask,
                        std::vector<bool> &                 selected_dofs,
                        const std::set<types::boundary_id> &boundary_ids = {});

  /**
   * 提取在边界上的所有自由度，并属于解决方案的指定组件。该函数以IndexSet的形式返回其结果，IndexSet包含与这些选定的自由度相对应的条目，也就是说，这些自由度位于边界并属于选定的成分之一。
   * 通过指定 @p boundary_ids
   * 变量，你可以选择要提取的自由度所在的面必须有哪些边界指标。如果它是一个空列表（默认），那么所有的边界指标都被接受。
   * 这个功能在  step-11  和  step-15  中使用，例如。
   *
   */
  template <int dim, int spacedim>
  IndexSet
  extract_boundary_dofs(const DoFHandler<dim, spacedim> &dof_handler,
                        const ComponentMask &component_mask = ComponentMask(),
                        const std::set<types::boundary_id> &boundary_ids = {});

  /**
   * 与前一个函数相同，只是它通过第三个参数返回其信息。
   * @deprecated  用前面的函数代替。
   *
   */
  template <int dim, int spacedim>
  DEAL_II_DEPRECATED void
  extract_boundary_dofs(const DoFHandler<dim, spacedim> &   dof_handler,
                        const ComponentMask &               component_mask,
                        IndexSet &                          selected_dofs,
                        const std::set<types::boundary_id> &boundary_ids = {});

  /**
   * 这个函数与extract_boundary_dofs()函数类似，但它提取那些形状函数在所选边界的至少一部分为非零的自由度。对于连续元素，这正是自由度定义在边界面上的形状函数的集合。另一方面，如果使用的有限元是不连续元，所有的自由度都定义在单元内部，因此没有一个是边界自由度。  然而，其中有几个自由度的形状函数在边界上是不为零的。因此，这个函数提取了所有那些 FiniteElement::has_support_on_face 函数说它在所选边界部分的任何面上都是非零的。      @see   @ref GlossBoundaryIndicator  "关于边界指标的词汇条目"
   *
   */
  template <int dim, int spacedim>
  void
  extract_dofs_with_support_on_boundary(
    const DoFHandler<dim, spacedim> &   dof_handler,
    const ComponentMask &               component_mask,
    std::vector<bool> &                 selected_dofs,
    const std::set<types::boundary_id> &boundary_ids =
      std::set<types::boundary_id>());

  /**
   提取形状函数的所有索引，使其支持完全包含在 @p predicate 为 <code>true</code> 的单元格中。  结果以IndexSet的形式返回。    考虑以下的FE空间，其中谓词对域的左半部的所有单元格返回 <code>true</code> 。      @image html extract_dofs_with_support_contained_within.png  这个函数将返回这些单元格上的所有DoF指数的联合，减去DoF 11, 13, 2和0；结果将是<code>[9,10], 12, [14,38]</code>。在上图中，返回的DoFs被红线隔开，从本质上讲，这个函数回答的问题如下。  给定一个带有相关DoFs的子域，这些DoFs中允许非零的最大子集是什么，以便在调用 AffineConstraints::distribute() 后，得到的解向量将只在给定的域内有支持。这里， @p constraints 是包含悬挂节点约束的AffineConstraints容器。    在 parallel::distributed::Triangulation 的情况下， @p predicate 将只为本地拥有的和幽灵单元调用。产生的索引集可能包含与本地拥有的或幽灵单元相关的DoF，但不为当前MPI核所拥有。
   *
   */
  template <int dim, int spacedim, typename number = double>
  IndexSet
  extract_dofs_with_support_contained_within(
    const DoFHandler<dim, spacedim> &dof_handler,
    const std::function<
      bool(const typename DoFHandler<dim, spacedim>::active_cell_iterator &)>
      &                              predicate,
    const AffineConstraints<number> &constraints = AffineConstraints<number>());

  /**
   * 提取一个向量，代表DoFHandler对<tt>component_mask</tt>（见 @ref
   * GlossComponentMask ）选择的组件的恒定模式。
   * 离散化的恒定模式是所选分量上的拉普拉斯算子的无效空间，并应用诺伊曼边界条件。在使用类
   * TrilinosWrappers::PreconditionAMG.
   * 时，无效空间是获得良好的AMG预处理的必要成分，因为ML
   * AMG包只对各自矩阵的代数属性工作，它没有机会检测矩阵是来自标量还是矢量值问题。
   * 然而，一个近乎无效的空间正好提供了所需的关于矩阵中矢量分量位置的信息。空空间（或者说，常数模式）是由给定的DoFHandler底层的有限元提供的，对于大多数元素，空空间将由与<tt>component_mask</tt>（见
   * @ref GlossComponentMask
   * ）中的真实参数一样多的向量组成，每个向量在一个向量分量中为1，在所有其他分量中为0。
   * 然而，例如FE_DGP的常数函数的表示是不同的（每个元素上的第一个分量为一，其他所有分量为零），有些标量元素甚至可能有两个常数模式（FE_Q_DG0）。因此，我们将这个对象存储在一个向量中，其中外向量包含DoFHandler上的实际恒定模式的集合。每个内向量有多少个分量，就有多少个所选分量中的（本地拥有的）自由度。请注意，任何与这个无效空间相关的矩阵都必须使用相同的<tt>component_mask</tt>参数来构建，因为自由度的编号是相对于所选的道夫而言的，而不是相对于所有道夫而言的。
   * 这个程序的主要原因是使用AMG预处理程序的空空间。
   *
   */
  template <int dim, int spacedim>
  void
  extract_constant_modes(const DoFHandler<dim, spacedim> &dof_handler,
                         const ComponentMask &            component_mask,
                         std::vector<std::vector<bool>> & constant_modes);
  //@}

  /**
   * @name  并行化和域分解 
     * @{ 
   *
   */
  /**
   * 标记所有在给定子域id的单元上的自由度。请注意，面的自由度可以属于不同子域id的单元，所以对于不同的子域id来说，被标记的自由度集并不相互排斥。
   * 如果你想得到自由度与子域的唯一关联，请使用 @p
   * get_subdomain_association 函数。
   *
   */
  template <int dim, int spacedim>
  void
  extract_subdomain_dofs(const DoFHandler<dim, spacedim> &dof_handler,
                         const types::subdomain_id        subdomain_id,
                         std::vector<bool> &              selected_dofs);

  /**
   * 提取在当前DoFHandler上有效的全局DoF指数集合。对于普通的DoFHandler来说，这些都是DoF指数，但是对于建立在 parallel::distributed::Triangulation 上的DoFHandler对象来说，这个集合是 DoFHandler::locally_owned_dofs() 的超集，包含了所有住在本地拥有的单元上的DoF指数（包括在与幽灵单元的接口上）。然而，它不包含专门定义在幽灵或人工单元上的自由度指数（见 @ref GlossArtificialCell "词汇表"
   * ）。
   * 这个函数识别的自由度等于从Dof_indices_with_subdomain_association()函数中获得的自由度，当调用本地拥有的子域ID时。
   *
   */
  template <int dim, int spacedim>
  void
  extract_locally_active_dofs(const DoFHandler<dim, spacedim> &dof_handler,
                              IndexSet &                       dof_set);

  /**
   * 与上述函数相同，但适用于某个（多网格）层次。
   * 这个函数返回在给定层次上所有本地拥有的单元（包括与幽灵单元的接口）上的所有DoF指数。
   *
   */
  template <int dim, int spacedim>
  void
  extract_locally_active_level_dofs(
    const DoFHandler<dim, spacedim> &dof_handler,
    IndexSet &                       dof_set,
    const unsigned int               level);

  /**
   * 提取在当前DoFHandler上活动的全局DoF指数集合。对于普通的DoFHandler，这些都是DoF指数，但是对于建立在 parallel::distributed::Triangulation 上的DoFHandler对象，这个集合是 DoFHandler::locally_owned_dofs() 和所有鬼魂单元上的DoF指数的联合。实质上，它是所有非人造单元上的DoF指数（见 @ref GlossArtificialCell "术语表"
   * ）。
   *
   */
  template <int dim, int spacedim>
  void
  extract_locally_relevant_dofs(const DoFHandler<dim, spacedim> &dof_handler,
                                IndexSet &                       dof_set);


  /**
   * 为掩码内的每个组件提取本地拥有的DoF指数集，这些指数为当前处理器所拥有。对于被掩码禁用的组件，会返回一个空的IndexSet。对于建立在顺序三角形上的标量DoFHandler，返回的向量包含一个包含所有DoF指数的完整IndexSet。如果掩码包含所有组件（这也对应于默认值），那么返回的索引集的联合就相当于
   * DoFHandler::locally_owned_dofs() 的返回值。
   *
   */
  template <int dim, int spacedim>
  std::vector<IndexSet>
  locally_owned_dofs_per_component(
    const DoFHandler<dim, spacedim> &dof_handler,
    const ComponentMask &            components = ComponentMask());

  /**
   * 对于每个处理器，确定本地拥有的自由度集合为一个IndexSet。然后这个函数返回一个索引集的向量，其中向量的大小等于参与自由度处理对象的MPI进程的数量。
   * 该函数可用于 dealii::Triangulation 或
   * parallel::shared::Triangulation. 类型的对象，但对
   * parallel::distributed::Triangulation
   * 类型的对象不起作用，因为对于这样的三角形，我们没有关于三角形的所有单元的本地可用信息，因此不能对其他处理器本地拥有的单元上的自由度有任何明确的说法。
   *
   */
  template <int dim, int spacedim>
  std::vector<IndexSet>
  locally_owned_dofs_per_subdomain(
    const DoFHandler<dim, spacedim> &dof_handler);

  /**
   * 对于每个处理器，确定本地相关自由度的集合为IndexSet。然后这个函数返回一个索引集的向量，其中向量的大小等于参与自由度处理对象的MPI进程的数量。
   * 该函数可用于 dealii::Triangulation 或
   * parallel::shared::Triangulation. 类型的对象，但对
   * parallel::distributed::Triangulation
   * 类型的对象不起作用，因为对于这样的三角形，我们没有关于三角形的所有单元的本地可用信息，因此不能对其他处理器本地拥有的单元上的自由度有任何明确的说法。
   *
   */
  template <int dim, int spacedim>
  std::vector<IndexSet>
  locally_relevant_dofs_per_subdomain(
    const DoFHandler<dim, spacedim> &dof_handler);


  /**
   * 与extract_locally_relevant_dofs()相同，但对于给定的 @p level.
   * 的多网格DoFs。
   *
   */
  template <int dim, int spacedim>
  void
  extract_locally_relevant_level_dofs(
    const DoFHandler<dim, spacedim> &dof_handler,
    const unsigned int               level,
    IndexSet &                       dof_set);


  /**
   * 对于每个自由度，在输出数组中返回它属于哪个子域（由<tt>cell->subdomain_id()</tt>函数给出）。在调用这个函数时，输出数组应该已经有了合适的大小。
   * 请注意，与面、边和顶点相关的自由度如果位于分区的边界上，可能与多个子域相关。在这种情况下，我们将它们分配给具有较小子域ID的进程。这可能会导致分区中自由度的数量不同，即使单元格的数量是完全等分的。虽然这是令人遗憾的，但在实践中这并不是一个问题，因为只要分区的数量保持不变，当我们细化网格时，分区边界上的自由度数量是渐进式消失的。
   * 这个函数返回每个DoF与一个子域的关联。如果你正在寻找每个
   * @em
   * 单元与一个子域的关联，可以查询<tt>cell->subdomain_id()</tt>函数，或者使用
   * <tt>GridTools::get_subdomain_association</tt> 函数。
   * 请注意，这个函数对于建立在
   * parallel::distributed::Triangulation
   * 上的DoFHandler对象的用途值得怀疑，因为在这种情况下，MPI进程对各个自由度的所有权是由DoF
   * handler对象控制的，而不是基于某种几何算法与子域id相结合。特别是，这个命名空间中的函数所识别的与子域相关的自由度与DoFHandler类所识别的自由度不一样。
   *
   */
  template <int dim, int spacedim>
  void
  get_subdomain_association(const DoFHandler<dim, spacedim> & dof_handler,
                            std::vector<types::subdomain_id> &subdomain);

  /**
   * 计算有多少个自由度与给定的 @p subdomain
   * 索引唯一相关。
   * 请注意，可能有一些罕见的情况，即具有给定 @p subdomain
   * 索引的单元格存在，但它的自由度实际上都没有与之相关。在这种情况下，返回值将为零。
   * 如果没有具有给定 @p subdomain
   * 索引的单元格，该函数将产生一个异常。
   * 该函数返回与一个子域相关的DoFs数量。
   * 如果你要寻找与这个子域相关的 @em 单元，请使用
   * <tt>GridTools::count_cells_with_subdomain_association</tt> 函数。
   * 注意这个函数对于建立在 parallel::distributed::Triangulation
   * 上的DoFHandler对象的用途值得怀疑，因为在这种情况下，MPI进程对单个自由度的所有权是由DoF
   * handler对象控制的，而不是基于一些与子域id相关的几何算法。特别是，这个命名空间中的函数所识别的与子域相关的自由度与DoFHandler类所识别的自由度不一样。
   *
   */
  template <int dim, int spacedim>
  unsigned int
  count_dofs_with_subdomain_association(
    const DoFHandler<dim, spacedim> &dof_handler,
    const types::subdomain_id        subdomain);

  /**
   * 计算有多少个自由度与给定的 @p subdomain
   * 索引唯一相关。
   * 这个函数的作用与前一个函数相同，只是它将结果在DoFHandler对象所使用的有限元的向量分量之间进行分割。因此，最后一个参数（其长度必须等于矢量分量的数量）将存储每个矢量分量的多少个自由度与给定的子域相关。
   * 注意这个函数对于建立在 parallel::distributed::Triangulation
   * 上的DoFHandler对象的用途是值得怀疑的，因为在这种情况下，MPI进程对各个自由度的所有权是由DoF
   * handler对象控制的，而不是基于一些与子域id相关的几何算法。特别是，这个命名空间中的函数所识别的与子域相关的自由度与DoFHandler类所识别的自由度不一样。
   *
   */
  template <int dim, int spacedim>
  void
  count_dofs_with_subdomain_association(
    const DoFHandler<dim, spacedim> &dof_handler,
    const types::subdomain_id        subdomain,
    std::vector<unsigned int> &      n_dofs_on_subdomain);

  /**
   * 返回一组索引，表示生活在给定子域上的自由度，即当前处理器拥有的单元上的自由度。请注意，这包括这个子域
   * "拥有
   * "的自由度（即get_subdomain_association()返回的值等于这里给出的子域，并且被
   * DoFHandler::locally_owned_dofs()
   * 函数选中的自由度），也包括所有位于给定子域和其他子域之间边界上的自由度。从本质上讲，位于子域之间边界的自由度将出现在这个函数返回的多个子域的索引集中。
   * 注意这个函数对于建立在 parallel::distributed::Triangulation
   * 上的DoFHandler对象的用途是值得怀疑的，因为在这种情况下，MPI进程对各个自由度的所有权是由DoF
   * handler对象控制的，而不是基于一些与子域id相关的几何算法。特别是，这个命名空间中的函数所识别的与子域相关的自由度与DoFHandler类所识别的自由度不一样。
   *
   */
  template <int dim, int spacedim>
  IndexSet
  dof_indices_with_subdomain_association(
    const DoFHandler<dim, spacedim> &dof_handler,
    const types::subdomain_id        subdomain);
  // @}
  /**
   * @name  细胞斑块上的DoF指数 为小块的细胞斑块创建包含大量自由度的结构。由此产生的对象可用于RelaxationBlockSOR和相关类，以实现Schwarz预处理和平滑器，其中子域仅由少量单元组成。
   *
   */
  //@{

  /**
   * 返回由参数描述的一组单元（即补丁）上的自由度集合。
   * 补丁通常用于定义误差估计器，这些估计器需要解决网格中每个单元周围补丁上的局部问题。你可以使用
   * GridTools::get_patch_around_cell().
   * 得到一个构成给定单元周围补丁的单元列表。虽然
   * DoFTools::count_dofs_on_patch()
   * 可以用来确定这些局部问题的大小，这样就可以组装局部系统，然后进行求解，但仍然有必要提供一个住在补丁上的自由度的全局索引和局部枚举之间的映射。这个函数通过返回住在补丁上的自由度的集合来提供这样一个局部列举。
   * 由于这个集合是以 std::vector,
   * 的形式返回的，我们也可以把它看成是一个映射
   * @code
   * i
   *
   * -> global_dof_index
   * @endcode
   * 其中 <code>i</code>
   * 是返回向量的索引（即补丁上一个自由度的<i>local</i>索引），
   * <code>global_dof_index</code>
   * 是位于补丁上的自由度的全局索引。返回的数组大小等于
   * DoFTools::count_dofs_on_patch(). 。
   * @note
   * 返回的数组是按全局自由度索引排序的。因此，如果我们认为这个数组的索引是本地DoF索引，那么产生的本地系统就保留了全局系统的块状结构。
   * @param  补丁 一个DoFHandler<dim,  spacedim>::active_cell_iterator
   * @return
   * 位于补丁上的那些全局自由度的列表，如上定义。
   * @note
   * 在并行分布式计算的背景下，只有在本地拥有的单元周围的补丁上调用这个函数才有意义。这是因为本地拥有的单元的邻居要么是本地拥有的单元，要么是幽灵单元。对于这两种情况，我们知道这些单元实际上是完整的、平行的三角形的真实单元。我们还可以查询这些单元的自由度。换句话说，这个函数只有在补丁中的所有单元都是本地拥有的或者是幽灵单元的情况下才能工作。
   *
   */
  template <int dim, int spacedim>
  std::vector<types::global_dof_index>
  get_dofs_on_patch(
    const std::vector<typename DoFHandler<dim, spacedim>::active_cell_iterator>
      &patch);

  /**
   * 与上述相同。      @deprecated
   * 使用以dim和spacedim为模板参数的函数。
   *
   */
  template <typename DoFHandlerType>
  DEAL_II_DEPRECATED std::vector<types::global_dof_index>
                     get_dofs_on_patch(
                       const std::vector<typename DoFHandlerType::active_cell_iterator> &patch);

  /**
   * 创建一个稀疏模式，它列出了与给定层次上每个单元相关的自由度。这种模式可以在RelaxationBlock类中作为加法和乘法施瓦茨方法的块列表。
   * 该模式中的行指数是通过三角法的一个层次进行标准迭代而得到的单元格指数。对于一个
   * parallel::distributed::Triangulation,
   * 来说，只有本地拥有的单元被输入。
   * 稀疏模式在这个函数中被调整为包含与给定级别上本地拥有的单元格一样多的行，与该级别上的自由度一样多的列。
   * <tt>selected_dofs</tt>是一个由单元上的局部自由度索引的向量。如果它被使用，只有这些自由度被输入到块列表中被选择。例如，这允许排除组件或边界上的道夫。
   *
   */
  template <int dim, int spacedim>
  void
  make_cell_patches(SparsityPattern &                block_list,
                    const DoFHandler<dim, spacedim> &dof_handler,
                    const unsigned int               level,
                    const std::vector<bool> &        selected_dofs = {},
                    const types::global_dof_index    offset        = 0);

  /**
   * 创建一个入射矩阵，对于多级DoFHandler的某一层的每一个顶点，标志着哪些自由度与相邻单元相关。这个数据结构是一个矩阵，有多少行就有多少顶点，有多少列就有多少自由度，条目是真还是假。这个数据结构由一个SparsityPattern对象方便地表示。
   * 在进入这个函数时，稀疏性模式可能是空的，将被重新初始化为正确的大小。
   * 该函数有一些布尔参数（列在下面）控制生成补丁的细节。默认设置是Arnold-Falk-Winther类型的平滑器，用于具有基本边界条件的发散和曲率符合的有限元。其他应用也是可能的，特别是改变<tt>boundary_patches</tt>用于非基本边界条件。
   * 这个函数返回<tt>vertex_mapping</tt>，它包含从顶点索引到<tt>block_list</tt>块索引的映射。对于没有导致顶点补丁的顶点，<tt>vertex_mapping</tt>中的条目包含值<tt>invalid_unsigned_int</tt>。如果<tt>invert_vertex_mapping</tt>被设置为<tt>true</tt>，那么<tt>vertex_mapping</tt>将被倒置，这样它就包含了从块索引到相应顶点索引的映射。
   * @arg  <tt>block_list</tt>：将存储补丁的SparsityPattern。
   * @arg
   * <tt>dof_handler</tt>：提供拓扑结构操作的多级dof处理程序。
   * @arg
   * <tt>interior_dofs_only</tt>：对于一个顶点周围的每个单元补丁，只收集该补丁的内部自由度，而不考虑该补丁边界上的自由度。例如，这是Arnold-Falk-Winther类型的平滑器的设置。
   * @arg
   * <tt>boundary_patches</tt>：包括域的边界顶点周围的补丁。如果不包括，将只生成内部顶点周围的补丁。
   * @arg
   * <tt>level_boundary_patches</tt>：对朝向更粗的单元的细化边也是如此。
   * @arg
   * <tt>single_cell_patches</tt>：如果不为真，包含单个单元的补丁会被消除。
   * @arg
   * <tt>invert_vertex_mapping</tt>：如果为真，那么返回值包含每个块的一个顶点索引；如果为假，那么返回值包含每个顶点的一个块索引或<tt>invalid_unsigned_int</tt>。
   *
   */
  template <int dim, int spacedim>
  std::vector<unsigned int>
  make_vertex_patches(SparsityPattern &                block_list,
                      const DoFHandler<dim, spacedim> &dof_handler,
                      const unsigned int               level,
                      const bool                       interior_dofs_only,
                      const bool                       boundary_patches = false,
                      const bool level_boundary_patches                 = false,
                      const bool single_cell_patches                    = false,
                      const bool invert_vertex_mapping = false);

  /**
   * 与上述相同，但允许单独排除块上的边界道夫。
   * 如果你想使用，例如，Taylor
   * Hood元素，这很有帮助，因为它允许你不包括速度块在补丁上的边界DoFs，同时也允许你包括压力块的边界DoFs。
   * 对于顶点周围的每个单元补丁，如果 @p exclude_boundary_dofs
   * 的BlockMask中对应块的布尔值为false，则收集该补丁的所有内部自由度并忽略该补丁边界上的自由度。
   *
   */
  template <int dim, int spacedim>
  std::vector<unsigned int>
  make_vertex_patches(SparsityPattern &                block_list,
                      const DoFHandler<dim, spacedim> &dof_handler,
                      const unsigned int               level,
                      const BlockMask &exclude_boundary_dofs  = BlockMask(),
                      const bool       boundary_patches       = false,
                      const bool       level_boundary_patches = false,
                      const bool       single_cell_patches    = false,
                      const bool       invert_vertex_mapping  = false);

  /**
   * 创建一个入射矩阵，对于多级DoFHandler的某一层的每一个单元，都标志着哪些自由度与这个单元的子代相关。这个数据结构可以方便地用SparsityPattern对象表示。
   * 因此，该函数创建了一个稀疏模式，在每一行（行对应于该层的单元）列出与该单元的子单元相关的自由度。这里使用的自由度指数是多级层次结构中的一级自由度指数，也就是说，它们可能与本身并不活跃的子单元有关。进入这个函数时，稀疏模式可能是空的，将被重新初始化为正确的大小。
   * 该函数有一些布尔参数（列在下面）控制生成补丁的细节。默认设置是Arnold-Falk-Winther类型的平滑器，用于具有基本边界条件的发散和曲率符合的有限元。其他应用也是可能的，特别是改变<tt>boundary_dofs</tt>用于非基本边界条件。
   * @arg  <tt>block_list</tt>：将存储补丁的SparsityPattern。
   * @arg
   * <tt>dof_handler</tt>：提供所操作的拓扑结构的多级dof处理器。
   * @arg
   * <tt>interior_dofs_only</tt>：对于顶点周围的每个单元补丁，只收集该补丁的内部自由度，而忽略该补丁边界上的自由度。例如，这就是Arnold-Falk-Winther类型的平滑器的设置。
   * @arg  <tt>boundary_dofs</tt>:
   * 包括自由度，这些自由度将被<tt>interior_dofs_only</tt>排除，但位于域的边界上，因此需要平滑。如果<tt>interior_dofs_only</tt>是假的，这个参数就没有影响。
   *
   */
  template <int dim, int spacedim>
  void
  make_child_patches(SparsityPattern &                block_list,
                     const DoFHandler<dim, spacedim> &dof_handler,
                     const unsigned int               level,
                     const bool                       interior_dofs_only,
                     const bool                       boundary_dofs = false);

  /**
   * 创建一个只有一个补丁的块列表，它又包含了给定层次上的所有自由度。
   * 这个函数主要是对make_child_patches()和make_vertex_patches()等函数在第0层的一个闭合，这些函数可能会产生一个空的补丁列表。
   * @arg  <tt>block_list</tt>: 补丁将被存储到的SparsityPattern。
   * @arg
   * <tt>dof_handler</tt>：提供拓扑结构操作的多级dof处理程序。
   * @arg  <tt>level</tt> 用于建立列表的网格级别。      @arg
   * <tt>interior_dofs_only</tt>:
   * 如果为真，排除域的边界上的自由度。
   *
   */
  template <int dim, int spacedim>
  void
  make_single_patch(SparsityPattern &                block_list,
                    const DoFHandler<dim, spacedim> &dof_handler,
                    const unsigned int               level,
                    const bool interior_dofs_only = false);

  /**
   * @}
   *
   */
  /**
   * @name  计算自由度和相关函数  
     * @{ 
   *
   */

  /**
   * 计算在总数中，有多少自由度属于每个组件。如果有限元的构件数是一个（即你只有一个标量变量），那么这个构件中的数字显然等于自由度总数。
   * 否则，所有组件中的自由度之和需要等于总数量。
   * 然而，如果有限元不是原始的，即它的一些或全部形状函数在一个以上的矢量分量中是非零的，那么最后一句话就不成立了。例如，这适用于Nedelec或Raviart-Thomas元素。在这种情况下，一个自由度在每个分量中都被计算为非零，因此上述的总和大于自由度的总数。
   * 这种行为可以通过可选的参数<tt>vector_valued_once</tt>来关闭。如果这是<tt>true</tt>，非原始向量值元素的成分数只收集在第一个成分中。所有其他分量的计数将为零。
   * 额外的可选参数 @p target_component
   * 允许对组件进行重新排序和分组。为此，它包含了每个组件的组件编号，它应该被计算为。如果多次输入相同的号码，就会把几个组件归为同一个。这个参数的应用之一是当你想形成块状矩阵和向量，但又想把几个分量打包到同一个块中时（例如，当你有
   * @p dim
   * 速度和一个压力时，要把所有速度放到一个块中，而把压力放到另一个块中）。
   * 结果在 @p dofs_per_component. 中返回。注意， @p
   * dofs_per_component的大小需要足以容纳 @p target_component.
   * 中指定的所有索引。如果不是这样，会抛出一个断言。
   * 没有被target_components锁定的索引将不被触及。
   *
   */
  template <int dim, int spacedim>
  std::vector<types::global_dof_index>
  count_dofs_per_fe_component(
    const DoFHandler<dim, spacedim> &dof_handler,
    const bool                       vector_valued_once = false,
    const std::vector<unsigned int> &target_component   = {});

  /**
   * 计算每个块中的自由度。这个函数类似于count_dofs_per_component()，不同的是，计数是按块进行的。详见术语表中的 @ref GlossBlock  "块"
   * 。在调用这个函数之前，再次假设向量具有正确的大小。如果不是这样，就会抛出一个断言。
   * 这个函数在  step-22  ,  step-31  , 和  step-32
   * 教程中使用，还有其他一些程序。      @pre
   * dofs_per_block变量具有与dof_handler参数所使用的有限元的块数相同的组件，或者与target_blocks参数中列举的块数相同（如果给出的话）。
   *
   */
  template <int dim, int spacedim>
  std::vector<types::global_dof_index>
  count_dofs_per_fe_block(const DoFHandler<dim, spacedim> &dof,
                          const std::vector<unsigned int> &target_block =
                            std::vector<unsigned int>());

  /**
   * 对于DoFHandler的每个活动单元，提取活动的有限元索引并填充作为第二个参数的矢量。这个向量被认为具有与活动单元相同数量的条目。
   * 对于没有hp-capabilities作为第一个参数的DoFHandler对象，返回的向量将只由0组成，表示所有单元使用相同的有限元。在hp模式下，这些值可能是不同的，但是。
   *
   */
  template <int dim, int spacedim>
  void
  get_active_fe_indices(const DoFHandler<dim, spacedim> &dof_handler,
                        std::vector<unsigned int> &      active_fe_indices);

  /**
   * 计算参数所描述的一组单元（即一个补丁）上有多少自由度。
   * 补丁通常用于定义误差估计器，这些估计器需要解决网格中每个单元周围补丁上的局部问题。你可以使用
   * GridTools::get_patch_around_cell().
   * 得到一个围绕给定单元的补丁的单元列表，这个函数在设置用于解决单元周围补丁的局部问题的线性系统的大小时非常有用。然后，函数
   * DoFTools::get_dofs_on_patch()
   * 将有助于建立全局自由度和局部自由度之间的联系。
   * @param  patch 一个DoFHandler<dim,
   * spacedim>类型的对象内的单元的集合  @return
   * 与这个补丁的单元相关的自由度数。
   * @note
   * 在并行分布式计算的背景下，只有在本地拥有的单元格周围的补丁上调用这个函数才有意义。这是因为本地拥有的单元的邻居要么是本地拥有的单元，要么是幽灵单元。对于这两种情况，我们知道这些单元实际上是完整的、平行的三角形的真实单元。我们还可以查询这些单元的自由度。换句话说，这个函数只有在补丁中的所有单元都是本地拥有的或者是幽灵单元的情况下才能工作。
   *
   */
  template <int dim, int spacedim>
  unsigned int
  count_dofs_on_patch(
    const std::vector<typename DoFHandler<dim, spacedim>::active_cell_iterator>
      &patch);

  /**
   * 与上述相同。      @deprecated
   * 使用以dim和spacedim为模板参数的函数。
   *
   */
  template <typename DoFHandlerType>
  DEAL_II_DEPRECATED unsigned int
  count_dofs_on_patch(
    const std::vector<typename DoFHandlerType::active_cell_iterator> &patch);

  /**
   * @}
   *
   */

  /**
   * @name  返回不同DoF映射的函数  
     * @{ 
   *
   */

  /**
   * 创建一个从自由度指数到该自由度在边界上的指数的映射。在此操作之后，<tt>mapping[dof]</tt>给出边界上自由度列表中全局编号为
   * @p dof 的自由度的索引。
   * 如果要求的自由度不在边界上，则<tt>mapping[dof]</tt>的值为
   * numbers::invalid_dof_index.
   * 该函数主要用于从试验函数中设置边界上的矩阵和向量，而矩阵和向量使用边界本地的试验函数的编号。
   * @p mapping 的先前内容被删除。
   *
   */
  template <int dim, int spacedim>
  void
  map_dof_to_boundary_indices(const DoFHandler<dim, spacedim> &     dof_handler,
                              std::vector<types::global_dof_index> &mapping);

  /**
   * 与之前的函数相同，只是只考虑边界的那些部分，对于这些部分的边界指标列在第二个参数中。    更多信息请参见本类的一般文档。      @see   @ref GlossBoundaryIndicator  "关于边界指标的词汇条目"
   *
   */
  template <int dim, int spacedim>
  void
  map_dof_to_boundary_indices(const DoFHandler<dim, spacedim> &   dof_handler,
                              const std::set<types::boundary_id> &boundary_ids,
                              std::vector<types::global_dof_index> &mapping);

  /**
   * 返回该DoF处理对象处理的所有自由度的支持点（见此 @ref GlossSupport "术语条目"
   * ）的列表。当然，这个函数只有在DoF处理对象使用的有限元对象实际提供支持点时才起作用，即没有边缘元素或类似的东西。否则，就会抛出一个异常。
   * @pre  给定的数组的长度必须与自由度的元素数量相同。
   * @note
   * 这个函数的前提条件是输出参数的大小必须等于自由度的总数，这使得这个函数不适合给定的DoFHandler对象派生自
   * parallel::TriangulationBase 对象的情况（或任何派生自
   * parallel::TriangulationBase).
   * 的类，因此，如果用这样的DoFHandler调用，这个函数将产生一个错误。
   * @param[in]  映射 从参考单元到定义DoF的实际单元的映射。
   * @param[in]  dof_handler
   * 描述哪个DoF指数在三角结构的哪个单元上的对象。
   * @param[in,out]  support_points
   * 存储实空间坐标中斗室的相应位置的向量。这个对象以前的内容在这个函数中被删除。
   * @param[in]  mask
   * 一个可选的分量掩码，用于限制从中提取支持点的分量。
   *
   */
  template <int dim, int spacedim>
  void
  map_dofs_to_support_points(const Mapping<dim, spacedim> &   mapping,
                             const DoFHandler<dim, spacedim> &dof_handler,
                             std::vector<Point<spacedim>> &   support_points,
                             const ComponentMask &mask = ComponentMask());

  /**
   * 与前面的函数相同，但用于hp-case。
   *
   */
  template <int dim, int spacedim>
  void
  map_dofs_to_support_points(
    const dealii::hp::MappingCollection<dim, spacedim> &mapping,
    const DoFHandler<dim, spacedim> &                   dof_handler,
    std::vector<Point<spacedim>> &                      support_points,
    const ComponentMask &                               mask = ComponentMask());

  /**
   *
   */
  template <int dim, int spacedim>
  void
  map_dofs_to_support_points(
    const Mapping<dim, spacedim> &                      mapping,
    const DoFHandler<dim, spacedim> &                   dof_handler,
    std::map<types::global_dof_index, Point<spacedim>> &support_points,
    const ComponentMask &                               mask = ComponentMask());

  /**
   * 与前面的函数相同，但用于hp-case。
   *
   */
  template <int dim, int spacedim>
  void
  map_dofs_to_support_points(
    const dealii::hp::MappingCollection<dim, spacedim> &mapping,
    const DoFHandler<dim, spacedim> &                   dof_handler,
    std::map<types::global_dof_index, Point<spacedim>> &support_points,
    const ComponentMask &                               mask = ComponentMask());


  /**
   * 这是一个与上面那个相反的函数。它生成一个地图，其中键是自由度的支持点，而值是自由度指数。关于支持点的定义，请看这个 @ref GlossSupport "词汇表条目"
   * 。
   * 由于在点的空间中没有自然的顺序（除了1d的情况），你必须提供一个带有明确指定的比较器对象的地图。因此，这个函数在比较器对象上被模板化。
   * 在这个函数中，地图对象的先前内容被删除。
   * 就像上面的函数一样，假定这里使用的有限元实际上支持其所有组件的支持点的概念。
   * @todo
   * 这个函数应该生成一个多图，而不仅仅是一个地图，因为几个道夫可能位于同一个支持点。目前，只有map_dofs_to_support_points()为每个点返回的地图中的最后一个值将被返回。
   *
   */
  template <int dim, int spacedim, class Comp>
  void
  map_support_points_to_dofs(
    const Mapping<dim, spacedim> &   mapping,
    const DoFHandler<dim, spacedim> &dof_handler,
    std::map<Point<spacedim>, types::global_dof_index, Comp>
      &point_to_index_map);
  /**
   * @}
   *
   */

  /**
   * @name  杂项  
     * @{ 
   *
   */

  /**
   * 取一个住在单元上的值的向量（例如，每个单元的误差），并以这样的方式将其分配到道夫上，从而产生一个有限元场，然后可以进一步处理，例如，用于输出。你应该注意到，所产生的场在悬挂的节点上将不是连续的。
   * 然而，这可以通过在矢量完全组装后调用为该DoFHandler对象创建的AffineConstraints对象的适当
   * @p 分布函数来轻松安排。    假设 @p cell_data
   * 中的元素数等于活动单元的数量， @p dof_data
   * 中的元素数等于<tt>dof_handler.n_dofs()</tt>。
   * 注意，输入向量可以是任何数据类型的向量，只要它可以转换为
   * @p double.
   * 输出向量，作为DoF处理程序上的数据向量，总是由 @p
   * double. 类型的元素组成。
   * 如果这个DoFHandler使用的有限元由一个以上的分量组成，你需要指定输出向量中的哪个分量应该用来存储有限元场；默认是0（如果有限元只由一个分量组成则不允许有其他值）。矢量的所有其他分量保持不动，即它们的内容不被改变。
   * 如果所使用的有限元的形状函数在一个以上的向量分量中是非零的（用deal.II的话说：它们是非正则的），则不能使用这个函数。
   *
   */
  template <int dim, int spacedim, typename Number>
  void
  distribute_cell_to_dof_vector(const DoFHandler<dim, spacedim> &dof_handler,
                                const Vector<Number> &           cell_data,
                                Vector<double> &                 dof_data,
                                const unsigned int               component = 0);


  /**
   * 根据给定的地图，用点数据生成gnuplot可读的文本输出
   * @p support_points.
   * 对于每个支持点的位置，生成一个包含地图上所有DoF列表的字符串标签。
   * 该地图可以通过调用map_dofs_to_support_points()来生成，对于可视化未知数的位置和全局编号非常有用。
   * 输出中每一行的格式的例子是。
   * @code
   * x [y] [z] "dof1, dof2"
   * @endcode
   * 其中x、y和z（只存在于相应的维度）是支持点的坐标，后面是一串DoF编号。
   * 带标签的点可以在gnuplot中作如下图示。
   * @code
   * plot "./points.gpl" using 1:2:3 with labels point offset 1,1
   * @endcode
   * 例子（这也包括用GridOut单独编写的网格）。    <p
   * ALIGN="center">
   @image html support_point_dofs1.png
   @image html support_point_dofs2.png
   * </p>
   * 要在单个gnuplot文件中生成网格和支撑点信息，请使用类似的代码
   * @code
   * std::ofstream out("gnuplot.gpl");
   * out << "plot '-' using 1:2 with lines, "
   *   << "'-' with labels point pt 2 offset 1,1"
   *   << std::endl;
   * GridOut().write_gnuplot (triangulation, out);
   * out << "e" << std::endl;
   *
   * std::map<types::global_dof_index, Point<dim> > support_points;
   * DoFTools::map_dofs_to_support_points (MappingQ1<dim>(),
   *                                     dof_handler,
   *                                     support_points);
   * DoFTools::write_gnuplot_dof_support_point_info(out,
   *                                              support_points);
   * out << "e" << std::endl;
   * @endcode
   * 并在gnuplot中执行以下命令。
   * @code
   * load "gnuplot.gpl"
   * @endcode
   * 或者，以下gnuplot脚本在命令行中以<tt>gnuplot
   * gnuplot.gpl</tt>的形式执行时将生成一个png文件。
   * @code
   * std::ofstream out("gnuplot.gpl");
   *
   * out << "set terminal png size 400,410 enhanced font \"Helvetica,8\"\n"
   *   << "set output \"output.png\"\n"
   *   << "set size square\n"
   *   << "set view equal xy\n"
   *   << "unset xtics\n"
   *   << "unset ytics\n"
   *   << "unset grid\n"
   *   << "unset border\n"
   *   << "plot '-' using 1:2 with lines notitle, "
   *   << "'-' with labels point pt 2 offset 1,1 notitle"
   *   << std::endl;
   * GridOut().write_gnuplot (triangulation, out);
   * out << "e" << std::endl;
   *
   * std::map<types::global_dof_index, Point<dim> > support_points;
   * DoFTools::map_dofs_to_support_points (MappingQ1<dim>(),
   *                                     dof_handler,
   *                                     support_points);
   * DoFTools::write_gnuplot_dof_support_point_info(out,
   *                                              support_points);
   * out << "e" << std::endl;
   * @endcode
   *
   *
   */
  template <int spacedim>
  void
  write_gnuplot_dof_support_point_info(
    std::ostream &                                            out,
    const std::map<types::global_dof_index, Point<spacedim>> &support_points);


  /**
   * 为 @p zero_boundary_constraints 添加约束，对应于在给定的边界指标上强制执行零边界条件。    这个函数约束了边界给定部分的所有自由度。    在 step-36 中使用了这个函数的一个变体，参数不同。      @param  dof 要工作的DoFHandler。    @param  boundary_id 应该被计算约束的那部分边界的指标。如果这个数字等于 numbers::invalid_boundary_id ，那么该域的所有边界都将被处理。    @param  zero_boundary_constraints 约束对象，约束条件将被写入其中。由于零边界值而产生的新约束将被简单地添加，保留之前存在的任何其他约束。然而，这只有在该对象以前的内容由不在这里处理的边界上的自由度的约束组成时才有效。如果以前有位于边界上的自由度的约束，那么这将构成冲突。参见 @ref constraints 模块，以处理个别自由度上存在冲突约束的情况。    @param  component_mask 一个可选的组件掩码，将这个函数的功能限制在一个FES系统的子集上。对于非 @ref GlossPrimitive "原始 "
   * 形状函数，任何属于形状函数的自由度都会受到影响，其中至少有一个非零分量受到分量屏蔽的影响（见
   * @ref GlossComponentMask  ）。
   * 如果省略这个参数，有限元中所有在边界上有自由度的分量将被考虑。
   * @ingroup constraints   @see   @ref GlossBoundaryIndicator  "关于边界指标的词汇条目"
   *
   */
  template <int dim, int spacedim, typename number>
  void
  make_zero_boundary_constraints(
    const DoFHandler<dim, spacedim> &dof,
    const types::boundary_id         boundary_id,
    AffineConstraints<number> &      zero_boundary_constraints,
    const ComponentMask &            component_mask = ComponentMask());

  /**
   * 与前一个函数相同，只是对边界的所有部分进行处理，而不仅仅是那些有特定边界指标的部分。那么这个函数就相当于以
   * numbers::invalid_boundary_id 为第二个参数调用前一个函数。
   * 这个函数在  step-36  中使用，例如。
   * @ingroup constraints
   *
   */
  template <int dim, int spacedim, typename number>
  void
  make_zero_boundary_constraints(
    const DoFHandler<dim, spacedim> &dof,
    AffineConstraints<number> &      zero_boundary_constraints,
    const ComponentMask &            component_mask = ComponentMask());

  /**
   * @}
   *
   */

  /**
   * @name  异常情况  
     * @{ 
   *
   */

  /**
   * @todo  编写说明
   *
   */
  DeclException0(ExcFiniteElementsDontMatch);
  /**
   * @todo  撰写描述
   * @ingroup Exceptions
   *
   */
  DeclException0(ExcGridNotCoarser);
  /**
   * @todo  编写描述 异常情况
   * @ingroup Exceptions
   *
   */
  DeclException0(ExcGridsDontMatch);
  /**
   * DoFHandler没有用有限元进行初始化。请先调用
   * DoFHandler::distribute_dofs() 。
   * @ingroup Exceptions
   *
   */
  DeclException0(ExcNoFESelected);
  /**
   * @todo  编写说明
   * @ingroup Exceptions
   *
   */
  DeclException0(ExcInvalidBoundaryIndicator);
  /**
   *
   */
} // namespace DoFTools



 /* ------------------------- inline functions -------------- */ 

#ifndef DOXYGEN

namespace DoFTools
{
  /**
   * 操作员计算出两个中的最大耦合度。      @relatesalso
   * DoFTools
   *
   */
  inline Coupling
  operator|=(Coupling &c1, const Coupling c2)
  {
    if (c2 == always)
      c1 = always;
    else if (c1 != always && c2 == nonzero)
      return c1 = nonzero;
    return c1;
  }


  /**
   * 计算两个人中的最大耦合度的操作者。      @relatesalso
   * DoFTools
   *
   */
  inline Coupling
  operator|(const Coupling c1, const Coupling c2)
  {
    if (c1 == always || c2 == always)
      return always;
    if (c1 == nonzero || c2 == nonzero)
      return nonzero;
    return none;
  }


  // ---------------------- inline and template functions --------------------

  template <int dim, int spacedim, class Comp>
  void
  map_support_points_to_dofs(
    const Mapping<dim, spacedim> &   mapping,
    const DoFHandler<dim, spacedim> &dof_handler,
    std::map<Point<spacedim>, types::global_dof_index, Comp>
      &point_to_index_map)
  {
    // let the checking of arguments be
    // done by the function first
    // called
    std::vector<Point<spacedim>> support_points(dof_handler.n_dofs());
    map_dofs_to_support_points(mapping, dof_handler, support_points);
    // now copy over the results of the
    // previous function into the
    // output arg
    point_to_index_map.clear();
    for (types::global_dof_index i = 0; i < dof_handler.n_dofs(); ++i)
      point_to_index_map[support_points[i]] = i;
  }



  template <typename DoFHandlerType, typename number>
  inline void
  make_periodicity_constraints(
    const std::vector<
      GridTools::PeriodicFacePair<typename DoFHandlerType::cell_iterator>>
      &                              periodic_faces,
    AffineConstraints<number> &      constraints,
    const ComponentMask &            component_mask,
    const std::vector<unsigned int> &first_vector_components,
    const number                     periodicity_factor)
  {
    make_periodicity_constraints<DoFHandlerType::dimension,
                                 DoFHandlerType::space_dimension>(
      periodic_faces,
      constraints,
      component_mask,
      first_vector_components,
      periodicity_factor);
  }



  template <typename DoFHandlerType>
  inline std::vector<types::global_dof_index>
  get_dofs_on_patch(
    const std::vector<typename DoFHandlerType::active_cell_iterator> &patch)
  {
    return get_dofs_on_patch<DoFHandlerType::dimension,
                             DoFHandlerType::space_dimension>(patch);
  }



  template <typename DoFHandlerType>
  inline unsigned int
  count_dofs_on_patch(
    const std::vector<typename DoFHandlerType::active_cell_iterator> &patch)
  {
    return count_dofs_on_patch<DoFHandlerType::dimension,
                               DoFHandlerType::space_dimension>(patch);
  }
} // namespace DoFTools

#endif

DEAL_II_NAMESPACE_CLOSE

#endif


