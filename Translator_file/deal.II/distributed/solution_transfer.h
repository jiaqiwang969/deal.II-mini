//include/deal.II-translator/distributed/solution_transfer_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2021 by the deal.II authors
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

#ifndef dealii_distributed_solution_transfer_h
#define dealii_distributed_solution_transfer_h

#include <deal.II/base/config.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <vector>


DEAL_II_NAMESPACE_OPEN

namespace parallel
{
  namespace distributed
  {
    /**
     * 在细化和/或粗化分布式网格的同时，通过内插法转移离散的FE函数（如解决方案矢量），并处理必要的通信。
     * @note
     * 需要注意的是，如果你同时使用一个以上的SolutionTransfer对象，那么对prepare_*()和interpolate()/deserialize()的调用需要以相同的顺序进行。
     * <h3>Note on ghost elements</h3>
     * 在并行计算中，PETSc或Trilinos矢量可能包含幽灵元素，也可能不包含。对于用prepare_for_coarsening_and_refinement()或prepare_for_serialization()读入信息，你需要提供有重影元素的向量，这样就可以读到所有本地活动元素。另一方面，有重影的向量通常是不可写的，所以对于调用interpolate()或deserialize()，你需要提供没有重影元素的分布式向量。更确切地说，在插值过程中，当前的算法会写进所有局部活动的自由度。
     * <h3>Transferring a solution</h3>
     * 这里VectorType是你喜欢的向量类型，例如
     * PETScWrappers::MPI::Vector,   TrilinosWrappers::MPI::Vector,
     * 或相应的块向量。
     * @code
     * SolutionTransfer<dim, VectorType> soltrans(dof_handler);
     * // flag some cells for refinement and coarsening, e.g.
     * GridRefinement::refine_and_coarsen_fixed_fraction(tria,
     *                                                 error_indicators,
     *                                                 0.3,
     *                                                 0.05);
     *
     * // prepare the triangulation,
     * tria.prepare_coarsening_and_refinement();
     *
     * // prepare the SolutionTransfer object for coarsening and refinement
     * // and give the solution vector that we intend to interpolate later,
     * soltrans.prepare_for_coarsening_and_refinement(solution);
     *
     * // actually execute the refinement,
     * tria.execute_coarsening_and_refinement ();
     *
     * // redistribute dofs,
     * dof_handler.distribute_dofs (fe);
     *
     * // and interpolate the solution
     * VectorType interpolated_solution;
     *
     * //create VectorType in the right size here
     * soltrans.interpolate(interpolated_solution);
     * @endcode
     * 由于网格是分布式的，需要注意的是，旧的解决方案（）必须被复制到一个也能提供对本地相关DoF值的访问（这些值在插值过程中需要）。
     * @code
     * // Create initial indexsets pertaining to the grid before refinement
     * IndexSet locally_owned_dofs, locally_relevant_dofs;
     * locally_owned_dofs = dof_handler.locally_owned_dofs();
     * DoFTools::extract_locally_relevant_dofs(dof_handler,
     * locally_relevant_dofs);
     *
     * // The solution vector only knows about locally owned DoFs
     * TrilinosWrappers::MPI::Vector solution;
     * solution.reinit(locally_owned_dofs,
     *               mpi_communicator);
     * ...
     * // Transfer solution to vector that provides access to
     * // locally relevant DoFs
     * TrilinosWrappers::MPI::Vector old_solution;
     * old_solution.reinit(locally_owned_dofs,
     *                   locally_relevant_dofs,
     *                   mpi_communicator);
     * old_solution = solution;
     *
     * // Initialize SolutionTransfer object
     * SolutionTransfer<dim, VectorType> soltrans(dof_handler);
     * soltrans.prepare_for_coarsening_and_refinement(old_solution);
     * ...
     * // Refine grid
     * // Recreate locally_owned_dofs and locally_relevant_dofs index sets
     * ...
     * solution.reinit(locally_owned_dofs, mpi_communicator);
     * soltrans.interpolate(solution);
     * @endcode
     * 与PETSc和Trilinos矢量不同， LinearAlgebra::distributed::Vector
     * 允许写进重影元素。
     * 对于一个重影向量，插值步骤可以通过以下方式完成
     * @code
     * interpolated_solution.zero_out_ghost_values();
     * soltrans.interpolate(interpolated_solution);
     * interpolated_solution.update_ghost_values();
     * @endcode
     * <h3>Use for Serialization</h3>
     * 这个类可以用来序列化和随后反序列化一个带有求解向量的分布式网格到一个文件。如果你使用一个以上的DoFHandler，从而使用一个以上的SolutionTransfer对象，它们需要按照相同的顺序进行序列化和反序列化。
     * 如果向量有本地相关的DoF，序列化工作如下。
     * @code
     * parallel::distributed::SolutionTransfer<dim,VectorType>
     * sol_trans(dof_handler);
     * sol_trans.prepare_for_serialization (vector);
     *
     * triangulation.save(filename);
     * @endcode
     * 对于反序列化，向量需要是一个分布式向量（没有鬼魂元素）。
     * @code
     * //[create coarse mesh...]
     * triangulation.load(filename);
     *
     * parallel::distributed::SolutionTransfer<dim,VectorType>
     * sol_trans(dof_handler);
     * sol_trans.deserialize (distributed_vector);
     * @endcode
     * <h3>Note on usage with DoFHandler with hp-capabilities</h3>
     * 由于具有hp-capabilities的DoFHandler对象上的数据与许多不同的FiniteElement对象相关联，每个单元的数据都必须用其相应的`future_fe_index`来处理。此外，如果涉及到细化，数据将以其`future_fe_index`打包在父单元上，随后以相同的索引在子单元上解包。
     * 如果单元被粗化为一个，数据将被打包在子单元上，以其共同子空间中最不占优势的有限元，并在父单元上以这个特定的有限元进行解包（更多信息请参考
     * hp::FECollection::find_dominated_fe_extended()  ）。
     * 在细化过程中转移解决方案的工作方式与非hp情况下完全相同。然而，当考虑序列化时，我们还必须在一个额外的步骤中存储活动FE指数。下面提供了一个代码片段，演示了使用
     * parallel::distributed::SolutionTransfer
     * 类与具有hp-capabilities的DoFHandler对象进行序列化。这里VectorType是你喜欢的向量类型，例如
     * PETScWrappers::MPI::Vector,   TrilinosWrappers::MPI::Vector,
     * 或相应的块向量。
     * 如果向量有本地相关的DoF，序列化工作如下。
     * @code
     * parallel::distributed::
     * SolutionTransfer<dim, VectorType, DoFHandler<dim,spacedim>>
     *   sol_trans(hp_dof_handler);
     *
     * hp_dof_handler.prepare_for_serialization_of_active_fe_indices();
     * sol_trans.prepare_for_serialization(vector);
     *
     * triangulation.save(filename);
     * @endcode
     * 对于反序列化，向量需要是一个分布式向量（没有鬼魂元素）。
     * @code
     * //[create coarse mesh...]
     * triangulation.load(filename);
     *
     * hp::FECollection<dim,spacedim> fe_collection;
     * //[prepare identical fe_collection...]
     *
     * DoFHandler<dim,spacedim> hp_dof_handler(triangulation);
     * // We need to introduce our dof_handler to the fe_collection
     * // before setting all active FE indices.
     * hp_dof_handler.deserialize_active_fe_indices();
     * hp_dof_handler.distribute_dofs(fe_collection);
     *
     * parallel::distributed::
     * SolutionTransfer<dim,VectorType,DoFHandler<dim,spacedim>>
     *   sol_trans(hp_dof_handler);
     * sol_trans.deserialize(distributed_vector);
     * @endcode
     * <h3>Interaction with hanging nodes</h3>
     * 从本质上讲，这个类实现了与 dealii::SolutionTransfer
     * 相同的步骤（尽管实现完全是独立的）。因此，这个类会发生与
     * dealii::SolutionTransfer.
     * 相同的挂起节点和粗化的问题，见那里的扩展讨论。
     * @ingroup distributed
     *
     */
    template <int dim, typename VectorType, int spacedim = dim>
    class SolutionTransfer
    {
    public:
      /**
       * 构造函数。              @param[in]  dof
       * 所有操作都将发生在它上面的DoFHandler。
       * 在这个构造函数被调用的时候，DoFHandler仍然指向有关细化发生之前的三角结构。
       *
       */
      SolutionTransfer(const DoFHandler<dim, spacedim> &dof);

      /**
       * 解构器。
       *
       */
      ~SolutionTransfer() = default;

      /**
       * 为粗化和细化的当前对象做准备。它存储每个单元的dof指数，并在每个将被粗化的单元中存储
       * @p all_in 中向量的dof值。  @p all_in
       * 包括所有将被内插到新的（细化和/或粗化的）网格上的向量。
       *
       */
      void
      prepare_for_coarsening_and_refinement(
        const std::vector<const VectorType *> &all_in);

      /**
       * 与前面的函数相同，但只对一个离散函数进行内插。
       *
       */
      void
      prepare_for_coarsening_and_refinement(const VectorType &in);

      /**
       * 在细化或粗化网格之前，将先前存储在此对象中的数据插值到当前的单元格组中。对提供给prepare_for_coarsening_and_refinement()的每个向量进行插值，并将结果写入给定的向量集中。
       *
       */
      void
      interpolate(std::vector<VectorType *> &all_out);

      /**
       * 与前一个函数相同。它只对一个函数进行插值。它假定向量有正确的大小（即<tt>in.size()==n_dofs_old</tt>,
       * <tt>out.size()==n_dofs_refined</tt>）。
       * 不允许多次调用此函数。通过使用<tt>interpolate (all_in,
       * all_out)</tt>，可以在一个步骤中执行多个函数的插值。
       *
       */
      void
      interpolate(VectorType &out);

      /**
       * 准备对给定的矢量进行序列化。序列化由
       * Triangulation::save().
       * 完成，给定的向量需要本地活动DoF的所有信息（必须是重影的）。更多信息请参见该类的文档。
       *
       */
      void
      prepare_for_serialization(const VectorType &in);

      /**
       * 和上面的函数一样，只是针对一个向量的列表。
       *
       */
      void
      prepare_for_serialization(const std::vector<const VectorType *> &all_in);

      /**
       * 执行给定向量的反序列化。这需要在调用
       * Triangulation::load().
       * 后进行，给定的向量必须是一个完全分布的向量，没有重影元素。更多信息请参见该类的文档。
       *
       */
      void
      deserialize(VectorType &in);


      /**
       * 和上面的函数一样，只是针对一个向量列表。
       *
       */
      void
      deserialize(std::vector<VectorType *> &all_in);

    private:
      /**
       * 指向自由度处理程序的指针，用于工作。
       *
       */
      SmartPointer<const DoFHandler<dim, spacedim>,
                   SolutionTransfer<dim, VectorType, spacedim>>
        dof_handler;

      /**
       * 一个存储指向所有我们应该从旧网格复制到新网格的向量的指针。
       *
       */
      std::vector<const VectorType *> input_vectors;

      /**
       * Triangulation分配给这个对象的句柄，我们可以用它来访问我们的内存偏移和我们的打包函数。
       *
       */
      unsigned int handle;

      /**
       * 一个回调函数，用于将当前网格上的数据打包成对象，以后可以在细化、粗化和重新划分后进行检索。
       *
       */
      std::vector<char>
      pack_callback(
        const typename Triangulation<dim, spacedim>::cell_iterator &cell,
        const typename Triangulation<dim, spacedim>::CellStatus     status);

      /**
       * 一个回调函数，用来解开当前网格上的数据，这些数据在细化、粗化和重新划分之前，已经被打包在网格上了。
       *
       */
      void
      unpack_callback(
        const typename Triangulation<dim, spacedim>::cell_iterator &cell,
        const typename Triangulation<dim, spacedim>::CellStatus     status,
        const boost::iterator_range<std::vector<char>::const_iterator>
          &                        data_range,
        std::vector<VectorType *> &all_out);


      /**
       * 将pack_callback()函数注册到已经分配给DoFHandler类成员的
       * parallel::distributed::Triangulation ，并存储返回的句柄。
       *
       */
      void
      register_data_attach();
    };
  } // namespace distributed
} // namespace parallel

namespace Legacy
{
  namespace parallel
  {
    namespace distributed
    {
      /**
       * @deprecated  使用没有DoFHandlerType模板的
       * dealii::parallel::distributed::SolutionTransfer 代替。
       *
       */
      template <int dim,
                typename VectorType,
                typename DoFHandlerType = DoFHandler<dim>>
      using SolutionTransfer DEAL_II_DEPRECATED =
        dealii::parallel::distributed::
          SolutionTransfer<dim, VectorType, DoFHandlerType::space_dimension>;
    } // namespace distributed
  }   // namespace parallel
} // namespace Legacy


DEAL_II_NAMESPACE_CLOSE

#endif


