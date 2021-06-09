//include/deal.II-translator/meshworker/functional_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2020 by the deal.II authors
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


#ifndef dealii_mesh_worker_functional_h
#define dealii_mesh_worker_functional_h

#include <deal.II/base/config.h>

#include <deal.II/algorithms/any_data.h>

#include <deal.II/base/mg_level_object.h>
#include <deal.II/base/smartpointer.h>

#include <deal.II/lac/block_vector.h>

#include <deal.II/meshworker/dof_info.h>

#include <deal.II/multigrid/mg_constrained_dofs.h>


DEAL_II_NAMESPACE_OPEN

namespace MeshWorker
{
  namespace Assembler
  {
    /**
     * 将对一个函数的局部贡献组装成全局函数的类。
     * @ingroup MeshWorker
     *
     */
    template <typename number = double>
    class Functional
    {
    public:
      /**
       * 初始化本地数据以存储函数。数字<tt>n</tt>是要计算的函数的数量。
       *
       */
      void
      initialize(const unsigned int n);
      /**
       * 初始化以后用于装配的DoFInfo对象中的本地数据。
       * 如果 <code>!face</code> ， @p info
       * 对象指的是一个单元，否则指的是一个内部或边界面。
       *
       */
      template <class DOFINFO>
      void
      initialize_info(DOFINFO &info, bool face);

      /**
       * 将局部值组装成全局向量。
       *
       */
      template <class DOFINFO>
      void
      assemble(const DOFINFO &info);

      /**
       * 将两个局部值集合到全局向量中。
       *
       */
      template <class DOFINFO>
      void
      assemble(const DOFINFO &info1, const DOFINFO &info2);

      /**
       * 在#results中的第i个条目的值。
       *
       */
      number
      operator()(const unsigned int i) const;

    private:
      /**
       * 将结果添加到其中的值。
       *
       */
      std::vector<double> results;
    };

    /**
     * 计算一个或几个函数的单元和面的贡献，通常用于误差估计。对于一个给定的单元或面，结果被储存在哪个组件中的信息是由其user_index变量传递的。因此，在使用这个类之前，你需要确保适当地设置这些变量。
     * @ingroup MeshWorker
     *
     */
    template <typename number = double>
    class CellsAndFaces
    {
    public:
      /**
       * 构造函数。初始化成员变量。
       *
       */
      CellsAndFaces();

      /**
       * 初始化函数，指定 @p results
       * 向量以及是否应单独收集脸部数据。              @p
       * results 应该包含两个名为 "单元格 "和 "面孔
       * "的块向量（后者只有在 @p separate_faces
       * 为真时才有）。在这两个中的每一个中，每个块都应该有相同的大小，并且大到足以容纳它所使用的循环所覆盖的单元格和面孔中设置的所有用户指数。通常，对于估计器来说，这分别是
       * Triangulation::n_active_cells() 和 Triangulation::n_faces(), 。
       * 使用BlockVector似乎很麻烦，但它允许我们同时组装几个函数，每个区块都有一个。误差估计的典型情况是在每个向量中只是有一个块。
       *
       */
      void
      initialize(AnyData &results, bool separate_faces = true);

      /**
       * 初始化以后用于组装的DoFInfo对象中的本地数据。
       * 如果 <code>!face</code> ， @p info
       * 对象指的是一个单元，否则指的是一个内部或边界面。
       *
       */
      template <class DOFINFO>
      void
      initialize_info(DOFINFO &info, bool face) const;

      /**
       * 将局部值组装成全局向量。
       *
       */
      template <class DOFINFO>
      void
      assemble(const DOFINFO &info);

      /**
       * 将两个局部值集合到全局向量中。
       *
       */
      template <class DOFINFO>
      void
      assemble(const DOFINFO &info1, const DOFINFO &info2);

      /**
       * @p results. 中第i项的值
       *
       */
      number
      operator()(const unsigned int i) const;

    private:
      AnyData results;
      bool    separate_faces;
    };
    //----------------------------------------------------------------------//

    template <typename number>
    inline void
    Functional<number>::initialize(const unsigned int n)
    {
      results.resize(n);
      std::fill(results.begin(), results.end(), 0.);
    }


    template <typename number>
    template <class DOFINFO>
    inline void
    Functional<number>::initialize_info(DOFINFO &info, bool)
    {
      info.initialize_numbers(results.size());
    }


    template <typename number>
    template <class DOFINFO>
    inline void
    Functional<number>::assemble(const DOFINFO &info)
    {
      for (unsigned int i = 0; i < results.size(); ++i)
        results[i] += info.value(i);
    }


    template <typename number>
    template <class DOFINFO>
    inline void
    Functional<number>::assemble(const DOFINFO &info1, const DOFINFO &info2)
    {
      for (unsigned int i = 0; i < results.size(); ++i)
        {
          results[i] += info1.value(i);
          results[i] += info2.value(i);
        }
    }


    template <typename number>
    inline number
    Functional<number>::operator()(const unsigned int i) const
    {
      AssertIndexRange(i, results.size());
      return results[i];
    }

    //----------------------------------------------------------------------//

    template <typename number>
    inline CellsAndFaces<number>::CellsAndFaces()
      : separate_faces(true)
    {}



    template <typename number>
    inline void
    CellsAndFaces<number>::initialize(AnyData &r, bool sep)
    {
      Assert(r.name(0) == "cells", AnyData::ExcNameMismatch(0, "cells"));
      if (sep)
        {
          Assert(r.name(1) == "faces", AnyData::ExcNameMismatch(1, "faces"));
          AssertDimension(r.entry<BlockVector<double> *>(0)->n_blocks(),
                          r.entry<BlockVector<double> *>(1)->n_blocks());
        }

      results        = r;
      separate_faces = sep;
    }

    template <typename number>
    template <class DOFINFO>
    inline void
    CellsAndFaces<number>::initialize_info(DOFINFO &info, bool) const
    {
      info.initialize_numbers(
        results.entry<BlockVector<double> *>(0)->n_blocks());
    }


    template <typename number>
    template <class DOFINFO>
    inline void
    CellsAndFaces<number>::assemble(const DOFINFO &info)
    {
      BlockVector<double> *v;
      if (separate_faces && info.face_number != numbers::invalid_unsigned_int)
        v = results.entry<BlockVector<double> *>(1);
      else
        v = results.entry<BlockVector<double> *>(0);

      for (unsigned int i = 0; i < info.n_values(); ++i)
        v->block(i)(info.cell->user_index()) += info.value(i);
    }


    template <typename number>
    template <class DOFINFO>
    inline void
    CellsAndFaces<number>::assemble(const DOFINFO &info1, const DOFINFO &info2)
    {
      for (unsigned int i = 0; i < info1.n_values(); ++i)
        {
          if (separate_faces)
            {
              BlockVector<double> *v1 = results.entry<BlockVector<double> *>(1);
              const double         J  = info1.value(i) + info2.value(i);
              v1->block(i)(info1.face->user_index()) += J;
              if (info2.face != info1.face)
                v1->block(i)(info2.face->user_index()) += J;
            }
          else
            {
              BlockVector<double> *v0 = results.entry<BlockVector<double> *>(0);
              v0->block(i)(info1.cell->user_index()) += .5 * info1.value(i);
              v0->block(i)(info2.cell->user_index()) += .5 * info2.value(i);
            }
        }
    }
  } // namespace Assembler
} // namespace MeshWorker

DEAL_II_NAMESPACE_CLOSE

#endif


