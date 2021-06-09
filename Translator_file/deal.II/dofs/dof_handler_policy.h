//include/deal.II-translator/dofs/dof_handler_policy_0.txt
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

#ifndef dealii_dof_handler_policy_h
#  define dealii_dof_handler_policy_h



#  include <deal.II/base/config.h>

#  include <deal.II/base/exceptions.h>
#  include <deal.II/base/template_constraints.h>

#  include <deal.II/dofs/dof_renumbering.h>
#  include <deal.II/dofs/dof_tools.h>

#  include <map>
#  include <set>
#  include <vector>

DEAL_II_NAMESPACE_OPEN

// Forward declaration
#  ifndef DOXYGEN
template <int, int>
class DoFHandler;
#  endif

namespace internal
{
  namespace DoFHandlerImplementation
  {
    struct NumberCache;

    /**
     * 一个命名空间，我们在其中定义描述如何分配和重新编号自由度的类。
     *
     */
    namespace Policy
    {
      struct Implementation;

      /**
       * 一个实现 DoFHandler::distribute_dofs 和
       * DoFHandler::renumber_dofs 函数应如何工作的策略的类。
       *
       */
      template <int dim, int spacedim>
      class PolicyBase
      {
      public:
        /**
         * 解构器。
         *
         */
        virtual ~PolicyBase() = default;

        /**
         * 在与该策略对象相关的DoFHandler对象上分配自由度。参数是对DoFHandler对象的NumberCache的引用。该函数可以修改它以使DoFHandler相关的函数在政策类中调用时正常工作。更新后的NumberCache被写入该参数中。
         *
         */
        virtual NumberCache
        distribute_dofs() const = 0;

        /**
         * 在与此策略对象相关的DoFHandler的每一层上分布多网格Dofs。返回所有层次的数字缓存的向量。
         *
         */
        virtual std::vector<NumberCache>
        distribute_mg_dofs() const = 0;

        /**
         * 按照第一个参数指定的自由度重新编号。
         * 在重新编号后，为DoFHandler返回一个更新的NumberCache。
         *
         */
        virtual NumberCache
        renumber_dofs(
          const std::vector<types::global_dof_index> &new_numbers) const = 0;

        /**
         * 对多网格层次结构中的一个层次的多级自由度进行重新编号。第二个参数指定新的DoF索引集。
         * 在重新编号后，为DoFHandler的指定层次返回一个更新的NumberCache。
         *
         */
        virtual NumberCache
        renumber_mg_dofs(
          const unsigned int                          level,
          const std::vector<types::global_dof_index> &new_numbers) const = 0;
      };


      /**
       * 该类实现了顺序操作的默认策略，即针对所有单元都得到自由度的情况。
       *
       */
      template <int dim, int spacedim>
      class Sequential : public PolicyBase<dim, spacedim>
      {
      public:
        /**
         * 构造函数。          @param  dof_handler
         * 这个策略类应该在其上工作的DoFHandler对象。
         *
         */
        Sequential(DoFHandler<dim, spacedim> &dof_handler);

        // documentation is inherited
        virtual NumberCache
        distribute_dofs() const override;

        // documentation is inherited
        virtual std::vector<NumberCache>
        distribute_mg_dofs() const override;

        // documentation is inherited
        virtual NumberCache
        renumber_dofs(const std::vector<types::global_dof_index> &new_numbers)
          const override;

        // documentation is inherited
        virtual NumberCache
        renumber_mg_dofs(const unsigned int level,
                         const std::vector<types::global_dof_index>
                           &new_numbers) const override;

      protected:
        /**
         * 这个策略对象赖以工作的DoFHandler对象。
         *
         */
        SmartPointer<DoFHandler<dim, spacedim>> dof_handler;
      };



      /**
       * 当我们使用 parallel::shared::Triangulation
       * 对象时，这个类实现了操作的策略。
       *
       */
      template <int dim, int spacedim>
      class ParallelShared : public PolicyBase<dim, spacedim>
      {
      public:
        /**
         * 构造函数。          @param  dof_handler
         * 这个策略类应该在其上工作的DoFHandler对象。
         *
         */
        ParallelShared(DoFHandler<dim, spacedim> &dof_handler);

        /**
         * 在作为第一个参数的对象上分配自由度。
         * 在分配时，自由度按子域重新编号，number_cache.n_locally_owned_dofs_per_processor[i]和number_cache.locally_owned_dofs被一致更新。
         *
         */
        virtual NumberCache
        distribute_dofs() const override;

        /**
         * 这个函数还没有实现。
         *
         */
        virtual std::vector<NumberCache>
        distribute_mg_dofs() const override;

        /**
         * 按照第一个参数指定的自由度重新编号。
         * 输入参数 @p new_numbers
         * 可以有和全局自由度一样多的条目（即dof_handler.n_dofs()）或者dof_handler.local_owned_dofs().n_elements()。因此，它可以利用为
         * parallel::distributed 情况下实施的重新编号函数。
         *
         */
        virtual NumberCache
        renumber_dofs(const std::vector<types::global_dof_index> &new_numbers)
          const override;

        // documentation is inherited
        virtual NumberCache
        renumber_mg_dofs(const unsigned int level,
                         const std::vector<types::global_dof_index>
                           &new_numbers) const override;

      private:
        /**
         * 这个策略对象所工作的DoFHandler对象。
         *
         */
        SmartPointer<DoFHandler<dim, spacedim>> dof_handler;
      };


      /**
       * 这个类实现了我们使用 parallel::DistributedTriangulationBase
       * 对象时的操作策略。
       *
       */
      template <int dim, int spacedim>
      class ParallelDistributed : public PolicyBase<dim, spacedim>
      {
      public:
        /**
         * 构造函数。          @param  dof_handler
         * 这个策略类应该在其上工作的DoFHandler对象。
         *
         */
        ParallelDistributed(DoFHandler<dim, spacedim> &dof_handler);

        // documentation is inherited
        virtual NumberCache
        distribute_dofs() const override;

        // documentation is inherited
        virtual std::vector<NumberCache>
        distribute_mg_dofs() const override;

        // documentation is inherited
        virtual NumberCache
        renumber_dofs(const std::vector<types::global_dof_index> &new_numbers)
          const override;

        // documentation is inherited
        virtual NumberCache
        renumber_mg_dofs(const unsigned int level,
                         const std::vector<types::global_dof_index>
                           &new_numbers) const override;

      private:
        /**
         * 这个策略对象赖以工作的DoFHandler对象。
         *
         */
        SmartPointer<DoFHandler<dim, spacedim>> dof_handler;
      };
    } // namespace Policy
  }   // namespace DoFHandlerImplementation
} // namespace internal



DEAL_II_NAMESPACE_CLOSE

#endif
 /*--------------------------   dof_handler_policy.h -------------------------*/ 


