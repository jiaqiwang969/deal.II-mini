//include/deal.II-translator/meshworker/simple_0.txt
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


#ifndef dealii_mesh_worker_simple_h
#define dealii_mesh_worker_simple_h

#include <deal.II/base/config.h>

#include <deal.II/algorithms/any_data.h>

#include <deal.II/base/mg_level_object.h>
#include <deal.II/base/smartpointer.h>

#include <deal.II/lac/block_vector.h>

#include <deal.II/meshworker/dof_info.h>
#include <deal.II/meshworker/functional.h>

#include <deal.II/multigrid/mg_constrained_dofs.h>

/* 包含 MeshWorker::Assembler::MatrixSimple,  MeshWorker::Assembler::MGMatrixSimple,  MeshWorker::Assembler::ResidualSimple, 和 MeshWorker::Assembler::SystemSimple. 类的页眉。

 
*
*/

DEAL_II_NAMESPACE_OPEN

namespace MeshWorker
{
  namespace Assembler
  {
    /**
     * 组装残差，没有块状结构。
     * 这个汇编器类的数据结构是每个单元上的一个简单向量，条目从零到
     * FiniteElementData::dofs_per_cell
     * ，以及一个简单的全局向量，条目编号从零到
     * DoFHandler::n_dofs().
     * ，不需要BlockInfo，全局向量可以是任何类型的向量，通过<tt>operator()
     * (unsigned int)</tt>有元素访问。
     * @ingroup MeshWorker
     *
     */
    template <typename VectorType>
    class ResidualSimple
    {
    public:
      /**
       * 用一个AnyData对象初始化，该对象保存着组装的结果。
       * 组装目前写入<tt>results</tt>的第一个向量中。
       *
       */
      void
      initialize(AnyData &results);

      /**
       * 初始化约束。
       *
       */
      void
      initialize(
        const AffineConstraints<typename VectorType::value_type> &constraints);

      /**
       * 初始化以后用于装配的DoFInfo对象中的本地数据。
       * 如果 @p info 对象指的是一个单元格，如果
       * <code>!face</code> ，则指的是一个内部或边界面。
       *
       */
      template <class DOFINFO>
      void
      initialize_info(DOFINFO &info, bool face) const;

      /**
       * 将局部残差组合成全局残差。
       * 数值被添加到之前的内容中。如果约束被激活，则使用
       * AffineConstraints::distribute_local_to_global() 。
       *
       */
      template <class DOFINFO>
      void
      assemble(const DOFINFO &info);

      /**
       * 将两个局部残差组合成全局残差。
       *
       */
      template <class DOFINFO>
      void
      assemble(const DOFINFO &info1, const DOFINFO &info2);

    protected:
      /**
       * 由assemble()填充的全局残差向量。
       *
       */
      AnyData residuals;

      /**
       * 一个指向包含约束的对象的指针。
       *
       */
      SmartPointer<const AffineConstraints<typename VectorType::value_type>,
                   ResidualSimple<VectorType>>
        constraints;
    };


    /**
     * 将本地矩阵组装成一个单一的全局矩阵或与同一DoFHandler相关的几个全局矩阵。如果这些全局矩阵有块状结构，则不使用该结构，而是使用全局自由度的编号。
     * 在用一个SparseMatrix对象（或另一个提供相同功能的矩阵
     * SparseMatrix::add())
     * 或一个这样的向量）初始化后，这个类可以在
     * MeshWorker::loop()
     * 中使用，将单元和面的矩阵组合成全局矩阵。
     * 如果在初始化过程中提供了一个AffineConstraints，这个矩阵将被用来
     * (AffineConstraints::distribute_local_to_global(),
     * 准确地说）将局部矩阵输入到全局稀疏矩阵。
     * 汇编器可以处理两种不同类型的局部数据。首先，默认情况下，明显的选择是取一个单一的局部矩阵，其尺寸等于单元的自由度数。
     * 或者，可以在DoFInfo中初始化一个局部块结构。
     * 在这之后，本地数据将被排列成一个由n乘n的FullMatrix块组成的数组（n是DoFInfo中DoFHandler使用的FES系统中的块数），这些块在DoFInfo中以列索引最快的方式排序。如果矩阵被初始化为几个矩阵的向量，并且使用了本地块结构，那么LocalResults中的前n<sup>2</sup>个矩阵将被用于该向量中的第一个矩阵，第二组n<sup>2</sup>个矩阵，以此类推。
     * @ingroup MeshWorker
     *
     */
    template <typename MatrixType>
    class MatrixSimple
    {
    public:
      /**
       * 构造函数，初始化#阈值，它限制了可以进入矩阵的小数字。
       *
       */
      MatrixSimple(double threshold = 1.e-12);

      /**
       * 存储结果矩阵，以便以后进行组装。
       *
       */
      void
      initialize(MatrixType &m);

      /**
       * 存储几个结果矩阵供以后组装。
       *
       */
      void
      initialize(std::vector<MatrixType> &m);

      /**
       * 初始化约束。在用一个有效的AffineConstraints对象调用此函数后，函数
       * AffineConstraints::distribute_local_to_global()
       * 将被assemble()用来将单元和面矩阵分配到一个全局稀疏矩阵中。
       *
       */
      void
      initialize(
        const AffineConstraints<typename MatrixType::value_type> &constraints);

      /**
       * 初始化以后用于装配的DoFInfo对象中的局部数据。
       * 如果 @p info 对象指的是一个单元格，如果
       * <code>!face</code> ，则指的是一个内部或边界面。
       *
       */
      template <class DOFINFO>
      void
      initialize_info(DOFINFO &info, bool face) const;

      /**
       * 将与单个单元相关的局部矩阵组装成全局矩阵。
       *
       */
      template <class DOFINFO>
      void
      assemble(const DOFINFO &info);

      /**
       * 将 @p info1 和 @p info2
       * 对象中与一个内部面相关的所有局部矩阵组装到全局矩阵中。
       *
       */
      template <class DOFINFO>
      void
      assemble(const DOFINFO &info1, const DOFINFO &info2);

    protected:
      /**
       * 正在组装的全局矩阵的向量。
       *
       */
      std::vector<SmartPointer<MatrixType, MatrixSimple<MatrixType>>> matrix;

      /**
       * 将被输入全局矩阵的最小正数。所有更小的绝对值将被视为零，不会被装配。
       *
       */
      const double threshold;

    private:
      /**
       * 将单个矩阵 <code>M</code> 装配到向量#matrix中
       * <code>index</code> 处的元素。
       *
       */
      void
      assemble(const FullMatrix<double> &                  M,
               const unsigned int                          index,
               const std::vector<types::global_dof_index> &i1,
               const std::vector<types::global_dof_index> &i2);

      /**
       * 一个指向包含约束的对象的指针。
       *
       */
      SmartPointer<const AffineConstraints<typename MatrixType::value_type>,
                   MatrixSimple<MatrixType>>
        constraints;
    };


    /**
     * 将局部矩阵组装成水平矩阵，而不使用块结构。
     * @todo
     * 用局部细化和连续元素组装水平矩阵所需的矩阵结构缺失。
     * @ingroup MeshWorker
     *
     */
    template <typename MatrixType>
    class MGMatrixSimple
    {
    public:
      /**
       * 构造函数，初始化#threshold，它限制了可以输入矩阵的小数字。
       *
       */
      MGMatrixSimple(double threshold = 1.e-12);

      /**
       * 存储结果矩阵，以便以后进行组装。
       *
       */
      void
      initialize(MGLevelObject<MatrixType> &m);

      /**
       * 初始化多级约束。
       *
       */
      void
      initialize(const MGConstrainedDoFs &mg_constrained_dofs);

      /**
       * 初始化矩阵#flux_up和#flux_down，用于不连续Galerkin方法的局部细化。
       *
       */
      void
      initialize_fluxes(MGLevelObject<MatrixType> &flux_up,
                        MGLevelObject<MatrixType> &flux_down);

      /**
       * 初始化矩阵#interface_in和#interface_out，用于连续Galerkin方法的局部细化。
       *
       */
      void
      initialize_interfaces(MGLevelObject<MatrixType> &interface_in,
                            MGLevelObject<MatrixType> &interface_out);
      /**
       * 初始化DoFInfo对象中的局部数据，以后用于装配。
       * 如果 @p info
       * 对象指的是一个单元，否则指的是一个内部或边界面。
       *
       */
      template <class DOFINFO>
      void
      initialize_info(DOFINFO &info, bool face) const;

      /**
       * 将矩阵 DoFInfo::M1[0] 组装成全局矩阵。
       *
       */
      template <class DOFINFO>
      void
      assemble(const DOFINFO &info);

      /**
       * 将 @p info1 和 @p info2
       * 对象中的两个局部矩阵组装成全局矩阵。
       *
       */
      template <class DOFINFO>
      void
      assemble(const DOFINFO &info1, const DOFINFO &info2);

    private:
      /**
       * 将单个矩阵组装成全局矩阵。
       *
       */
      void
      assemble(MatrixType &                                G,
               const FullMatrix<double> &                  M,
               const std::vector<types::global_dof_index> &i1,
               const std::vector<types::global_dof_index> &i2);

      /**
       * 将一个单一的矩阵组装成一个全局矩阵。
       *
       */
      void
      assemble(MatrixType &                                G,
               const FullMatrix<double> &                  M,
               const std::vector<types::global_dof_index> &i1,
               const std::vector<types::global_dof_index> &i2,
               const unsigned int                          level);

      /**
       * 将一个单一的矩阵组合成一个全局矩阵。
       *
       */

      void
      assemble_up(MatrixType &                                G,
                  const FullMatrix<double> &                  M,
                  const std::vector<types::global_dof_index> &i1,
                  const std::vector<types::global_dof_index> &i2,
                  const unsigned int level = numbers::invalid_unsigned_int);
      /**
       * 将一个单一的矩阵组合成一个全局矩阵。
       *
       */

      void
      assemble_down(MatrixType &                                G,
                    const FullMatrix<double> &                  M,
                    const std::vector<types::global_dof_index> &i1,
                    const std::vector<types::global_dof_index> &i2,
                    const unsigned int level = numbers::invalid_unsigned_int);

      /**
       * 将一个单一的矩阵组合成一个全局矩阵。
       *
       */

      void
      assemble_in(MatrixType &                                G,
                  const FullMatrix<double> &                  M,
                  const std::vector<types::global_dof_index> &i1,
                  const std::vector<types::global_dof_index> &i2,
                  const unsigned int level = numbers::invalid_unsigned_int);

      /**
       * 将一个单一的矩阵组合成一个全局矩阵。
       *
       */

      void
      assemble_out(MatrixType &                                G,
                   const FullMatrix<double> &                  M,
                   const std::vector<types::global_dof_index> &i1,
                   const std::vector<types::global_dof_index> &i2,
                   const unsigned int level = numbers::invalid_unsigned_int);

      /**
       * 正在装配的全局矩阵。
       *
       */
      SmartPointer<MGLevelObject<MatrixType>, MGMatrixSimple<MatrixType>>
        matrix;

      /**
       * 用于细化边缘的面通量项的矩阵，将粗的耦合到细的。
       *
       */
      SmartPointer<MGLevelObject<MatrixType>, MGMatrixSimple<MatrixType>>
        flux_up;

      /**
       * 用于细化边缘的面通量项的矩阵，将细化与粗化相耦合。
       *
       */
      SmartPointer<MGLevelObject<MatrixType>, MGMatrixSimple<MatrixType>>
        flux_down;

      /**
       * 用于整个细化边缘的连续元素的面贡献的矩阵，耦合粗到细。
       *
       */
      SmartPointer<MGLevelObject<MatrixType>, MGMatrixSimple<MatrixType>>
        interface_in;

      /**
       * 用于整个细化边缘的连续元素的面贡献的矩阵，将细化与粗化相耦合。
       *
       */
      SmartPointer<MGLevelObject<MatrixType>, MGMatrixSimple<MatrixType>>
        interface_out;
      /**
       * 一个指向包含约束的对象的指针。
       *
       */
      SmartPointer<const MGConstrainedDoFs, MGMatrixSimple<MatrixType>>
        mg_constrained_dofs;

      /**
       * 将被输入全局矩阵的最小正数。所有更小的绝对值将被视为零，不会被集合。
       *
       */
      const double threshold;
    };


    /**
     * 一次性组装一个简单的矩阵和一个简单的右手边。我们使用MatrixSimple和ResidualSimple的组合来实现这一点。单元和面操作者应该在LocalResults中填充矩阵和向量对象，这个类将把它们组装成矩阵和向量对象。
     * @ingroup MeshWorker
     *
     */
    template <typename MatrixType, typename VectorType>
    class SystemSimple : private MatrixSimple<MatrixType>,
                         private ResidualSimple<VectorType>
    {
    public:
      /**
       * 构造函数在MatrixSimple中设置阈值。
       *
       */
      SystemSimple(double threshold = 1.e-12);

      /**
       * 存储两个对象的数据被组装成。
       *
       */
      void
      initialize(MatrixType &m, VectorType &rhs);

      /**
       * 初始化约束。在用一个有效的AffineConstraints对象调用此函数后，函数
       * AffineConstraints::distribute_local_to_global()
       * 将被assemble()用来将单元和面矩阵分配到一个全局稀疏矩阵中。
       *
       */
      void
      initialize(
        const AffineConstraints<typename VectorType::value_type> &constraints);

      /**
       * 初始化以后用于装配的DoFInfo对象中的局部数据。
       * 如果 @p info
       * 对象指的是一个单元，否则指的是一个内部或边界面。
       *
       */
      template <class DOFINFO>
      void
      initialize_info(DOFINFO &info, bool face) const;

      /**
       * 将矩阵 DoFInfo::M1[0] 组装成全局矩阵。
       *
       */
      template <class DOFINFO>
      void
      assemble(const DOFINFO &info);

      /**
       * 将 @p info1 和 @p info2
       * 对象中的两个局部矩阵组装成全局矩阵。
       *
       */
      template <class DOFINFO>
      void
      assemble(const DOFINFO &info1, const DOFINFO &info2);

    private:
      /**
       * 将单个矩阵 <code>M</code> 装配到向量#matrix中
       * <code>index</code> 处的元素中。
       *
       */
      void
      assemble(const FullMatrix<double> &                  M,
               const Vector<double> &                      vector,
               const unsigned int                          index,
               const std::vector<types::global_dof_index> &indices);

      void
      assemble(const FullMatrix<double> &                  M,
               const Vector<double> &                      vector,
               const unsigned int                          index,
               const std::vector<types::global_dof_index> &i1,
               const std::vector<types::global_dof_index> &i2);
    };


    //----------------------------------------------------------------------//

    template <typename VectorType>
    inline void
    ResidualSimple<VectorType>::initialize(AnyData &results)
    {
      residuals = results;
    }



    template <typename VectorType>
    inline void
    ResidualSimple<VectorType>::initialize(
      const AffineConstraints<typename VectorType::value_type> &c)
    {
      constraints = &c;
    }



    template <typename VectorType>
    template <class DOFINFO>
    inline void
    ResidualSimple<VectorType>::initialize_info(DOFINFO &info, bool) const
    {
      info.initialize_vectors(residuals.size());
    }



    template <typename VectorType>
    template <class DOFINFO>
    inline void
    ResidualSimple<VectorType>::assemble(const DOFINFO &info)
    {
      for (unsigned int k = 0; k < residuals.size(); ++k)
        {
          VectorType *v = residuals.entry<VectorType *>(k);
          for (unsigned int i = 0; i != info.vector(k).n_blocks(); ++i)
            {
              const std::vector<types::global_dof_index> &ldi =
                info.vector(k).n_blocks() == 1 ? info.indices :
                                                 info.indices_by_block[i];

              if (constraints != nullptr)
                constraints->distribute_local_to_global(info.vector(k).block(i),
                                                        ldi,
                                                        *v);
              else
                v->add(ldi, info.vector(k).block(i));
            }
        }
    }

    template <typename VectorType>
    template <class DOFINFO>
    inline void
    ResidualSimple<VectorType>::assemble(const DOFINFO &info1,
                                         const DOFINFO &info2)
    {
      assemble(info1);
      assemble(info2);
    }


    //----------------------------------------------------------------------//

    template <typename MatrixType>
    inline MatrixSimple<MatrixType>::MatrixSimple(double threshold)
      : threshold(threshold)
    {}


    template <typename MatrixType>
    inline void
    MatrixSimple<MatrixType>::initialize(MatrixType &m)
    {
      matrix.resize(1);
      matrix[0] = &m;
    }


    template <typename MatrixType>
    inline void
    MatrixSimple<MatrixType>::initialize(std::vector<MatrixType> &m)
    {
      matrix.resize(m.size());
      for (unsigned int i = 0; i < m.size(); ++i)
        matrix[i] = &m[i];
    }


    template <typename MatrixType>
    inline void
    MatrixSimple<MatrixType>::initialize(
      const AffineConstraints<typename MatrixType::value_type> &c)
    {
      constraints = &c;
    }


    template <typename MatrixType>
    template <class DOFINFO>
    inline void
    MatrixSimple<MatrixType>::initialize_info(DOFINFO &info, bool face) const
    {
      Assert(matrix.size() != 0, ExcNotInitialized());

      const unsigned int n = info.indices_by_block.size();

      if (n == 0)
        info.initialize_matrices(matrix.size(), face);
      else
        {
          info.initialize_matrices(matrix.size() * n * n, face);
          unsigned int k = 0;
          for (unsigned int m = 0; m < matrix.size(); ++m)
            for (unsigned int i = 0; i < n; ++i)
              for (unsigned int j = 0; j < n; ++j, ++k)
                {
                  info.matrix(k, false).row    = i;
                  info.matrix(k, false).column = j;
                  if (face)
                    {
                      info.matrix(k, true).row    = i;
                      info.matrix(k, true).column = j;
                    }
                }
        }
    }



    template <typename MatrixType>
    inline void
    MatrixSimple<MatrixType>::assemble(
      const FullMatrix<double> &                  M,
      const unsigned int                          index,
      const std::vector<types::global_dof_index> &i1,
      const std::vector<types::global_dof_index> &i2)
    {
      AssertDimension(M.m(), i1.size());
      AssertDimension(M.n(), i2.size());

      if (constraints == nullptr)
        {
          for (unsigned int j = 0; j < i1.size(); ++j)
            for (unsigned int k = 0; k < i2.size(); ++k)
              if (std::fabs(M(j, k)) >= threshold)
                matrix[index]->add(i1[j], i2[k], M(j, k));
        }
      else
        constraints->distribute_local_to_global(M, i1, i2, *matrix[index]);
    }


    template <typename MatrixType>
    template <class DOFINFO>
    inline void
    MatrixSimple<MatrixType>::assemble(const DOFINFO &info)
    {
      Assert(!info.level_cell, ExcMessage("Cell may not access level dofs"));
      const unsigned int n = info.indices_by_block.size();

      if (n == 0)
        for (unsigned int m = 0; m < matrix.size(); ++m)
          assemble(info.matrix(m, false).matrix, m, info.indices, info.indices);
      else
        {
          for (unsigned int m = 0; m < matrix.size(); ++m)
            for (unsigned int k = 0; k < n * n; ++k)
              {
                assemble(
                  info.matrix(k + m * n * n, false).matrix,
                  m,
                  info.indices_by_block[info.matrix(k + m * n * n, false).row],
                  info.indices_by_block[info.matrix(k + m * n * n, false)
                                          .column]);
              }
        }
    }


    template <typename MatrixType>
    template <class DOFINFO>
    inline void
    MatrixSimple<MatrixType>::assemble(const DOFINFO &info1,
                                       const DOFINFO &info2)
    {
      Assert(!info1.level_cell, ExcMessage("Cell may not access level dofs"));
      Assert(!info2.level_cell, ExcMessage("Cell may not access level dofs"));
      AssertDimension(info1.indices_by_block.size(),
                      info2.indices_by_block.size());

      const unsigned int n = info1.indices_by_block.size();

      if (n == 0)
        {
          for (unsigned int m = 0; m < matrix.size(); ++m)
            {
              assemble(info1.matrix(m, false).matrix,
                       m,
                       info1.indices,
                       info1.indices);
              assemble(info1.matrix(m, true).matrix,
                       m,
                       info1.indices,
                       info2.indices);
              assemble(info2.matrix(m, false).matrix,
                       m,
                       info2.indices,
                       info2.indices);
              assemble(info2.matrix(m, true).matrix,
                       m,
                       info2.indices,
                       info1.indices);
            }
        }
      else
        {
          for (unsigned int m = 0; m < matrix.size(); ++m)
            for (unsigned int k = 0; k < n * n; ++k)
              {
                const unsigned int row = info1.matrix(k + m * n * n, false).row;
                const unsigned int column =
                  info1.matrix(k + m * n * n, false).column;

                assemble(info1.matrix(k + m * n * n, false).matrix,
                         m,
                         info1.indices_by_block[row],
                         info1.indices_by_block[column]);
                assemble(info1.matrix(k + m * n * n, true).matrix,
                         m,
                         info1.indices_by_block[row],
                         info2.indices_by_block[column]);
                assemble(info2.matrix(k + m * n * n, false).matrix,
                         m,
                         info2.indices_by_block[row],
                         info2.indices_by_block[column]);
                assemble(info2.matrix(k + m * n * n, true).matrix,
                         m,
                         info2.indices_by_block[row],
                         info1.indices_by_block[column]);
              }
        }
    }


    //----------------------------------------------------------------------//

    template <typename MatrixType>
    inline MGMatrixSimple<MatrixType>::MGMatrixSimple(double threshold)
      : threshold(threshold)
    {}


    template <typename MatrixType>
    inline void
    MGMatrixSimple<MatrixType>::initialize(MGLevelObject<MatrixType> &m)
    {
      matrix = &m;
    }

    template <typename MatrixType>
    inline void
    MGMatrixSimple<MatrixType>::initialize(const MGConstrainedDoFs &c)
    {
      mg_constrained_dofs = &c;
    }


    template <typename MatrixType>
    inline void
    MGMatrixSimple<MatrixType>::initialize_fluxes(
      MGLevelObject<MatrixType> &up,
      MGLevelObject<MatrixType> &down)
    {
      flux_up   = &up;
      flux_down = &down;
    }


    template <typename MatrixType>
    inline void
    MGMatrixSimple<MatrixType>::initialize_interfaces(
      MGLevelObject<MatrixType> &in,
      MGLevelObject<MatrixType> &out)
    {
      interface_in  = &in;
      interface_out = &out;
    }


    template <typename MatrixType>
    template <class DOFINFO>
    inline void
    MGMatrixSimple<MatrixType>::initialize_info(DOFINFO &info, bool face) const
    {
      const unsigned int n = info.indices_by_block.size();

      if (n == 0)
        info.initialize_matrices(1, face);
      else
        {
          info.initialize_matrices(n * n, face);
          unsigned int k = 0;
          for (unsigned int i = 0; i < n; ++i)
            for (unsigned int j = 0; j < n; ++j, ++k)
              {
                info.matrix(k, false).row    = i;
                info.matrix(k, false).column = j;
                if (face)
                  {
                    info.matrix(k, true).row    = i;
                    info.matrix(k, true).column = j;
                  }
              }
        }
    }


    template <typename MatrixType>
    inline void
    MGMatrixSimple<MatrixType>::assemble(
      MatrixType &                                G,
      const FullMatrix<double> &                  M,
      const std::vector<types::global_dof_index> &i1,
      const std::vector<types::global_dof_index> &i2)
    {
      AssertDimension(M.m(), i1.size());
      AssertDimension(M.n(), i2.size());
      Assert(mg_constrained_dofs == 0, ExcInternalError());
      // TODO: Possibly remove this function all together

      for (unsigned int j = 0; j < i1.size(); ++j)
        for (unsigned int k = 0; k < i2.size(); ++k)
          if (std::fabs(M(j, k)) >= threshold)
            G.add(i1[j], i2[k], M(j, k));
    }


    template <typename MatrixType>
    inline void
    MGMatrixSimple<MatrixType>::assemble(
      MatrixType &                                G,
      const FullMatrix<double> &                  M,
      const std::vector<types::global_dof_index> &i1,
      const std::vector<types::global_dof_index> &i2,
      const unsigned int                          level)
    {
      AssertDimension(M.m(), i1.size());
      AssertDimension(M.n(), i2.size());

      if (mg_constrained_dofs == nullptr)
        {
          for (unsigned int j = 0; j < i1.size(); ++j)
            for (unsigned int k = 0; k < i2.size(); ++k)
              if (std::fabs(M(j, k)) >= threshold)
                G.add(i1[j], i2[k], M(j, k));
        }
      else
        {
          for (unsigned int j = 0; j < i1.size(); ++j)
            for (unsigned int k = 0; k < i2.size(); ++k)
              {
                // Only enter the local values into the global matrix,
                //  if the value is larger than the threshold
                if (std::fabs(M(j, k)) < threshold)
                  continue;

                // Do not enter, if either the row or the column
                // corresponds to an index on the refinement edge. The
                // level problems are solved with homogeneous
                // Dirichlet boundary conditions, therefore we
                // eliminate these rows and columns. The corresponding
                // matrix entries are entered by assemble_in() and
                // assemble_out().
                if (mg_constrained_dofs->at_refinement_edge(level, i1[j]) ||
                    mg_constrained_dofs->at_refinement_edge(level, i2[k]))
                  continue;

                // At the boundary, only enter the term on the
                // diagonal, but not the coupling terms
                if ((mg_constrained_dofs->is_boundary_index(level, i1[j]) ||
                     mg_constrained_dofs->is_boundary_index(level, i2[k])) &&
                    (i1[j] != i2[k]))
                  continue;

                G.add(i1[j], i2[k], M(j, k));
              }
        }
    }


    template <typename MatrixType>
    inline void
    MGMatrixSimple<MatrixType>::assemble_up(
      MatrixType &                                G,
      const FullMatrix<double> &                  M,
      const std::vector<types::global_dof_index> &i1,
      const std::vector<types::global_dof_index> &i2,
      const unsigned int                          level)
    {
      AssertDimension(M.n(), i1.size());
      AssertDimension(M.m(), i2.size());

      if (mg_constrained_dofs == nullptr)
        {
          for (unsigned int j = 0; j < i1.size(); ++j)
            for (unsigned int k = 0; k < i2.size(); ++k)
              if (std::fabs(M(k, j)) >= threshold)
                G.add(i1[j], i2[k], M(k, j));
        }
      else
        {
          for (unsigned int j = 0; j < i1.size(); ++j)
            for (unsigned int k = 0; k < i2.size(); ++k)
              if (std::fabs(M(k, j)) >= threshold)
                if (!mg_constrained_dofs->at_refinement_edge(level, i2[k]))
                  G.add(i1[j], i2[k], M(k, j));
        }
    }

    template <typename MatrixType>
    inline void
    MGMatrixSimple<MatrixType>::assemble_down(
      MatrixType &                                G,
      const FullMatrix<double> &                  M,
      const std::vector<types::global_dof_index> &i1,
      const std::vector<types::global_dof_index> &i2,
      const unsigned int                          level)
    {
      AssertDimension(M.m(), i1.size());
      AssertDimension(M.n(), i2.size());

      if (mg_constrained_dofs == nullptr)
        {
          for (unsigned int j = 0; j < i1.size(); ++j)
            for (unsigned int k = 0; k < i2.size(); ++k)
              if (std::fabs(M(j, k)) >= threshold)
                G.add(i1[j], i2[k], M(j, k));
        }
      else
        {
          for (unsigned int j = 0; j < i1.size(); ++j)
            for (unsigned int k = 0; k < i2.size(); ++k)
              if (std::fabs(M(j, k)) >= threshold)
                if (!mg_constrained_dofs->at_refinement_edge(level, i2[k]))
                  G.add(i1[j], i2[k], M(j, k));
        }
    }

    template <typename MatrixType>
    inline void
    MGMatrixSimple<MatrixType>::assemble_in(
      MatrixType &                                G,
      const FullMatrix<double> &                  M,
      const std::vector<types::global_dof_index> &i1,
      const std::vector<types::global_dof_index> &i2,
      const unsigned int                          level)
    {
      AssertDimension(M.m(), i1.size());
      AssertDimension(M.n(), i2.size());
      Assert(mg_constrained_dofs != nullptr, ExcInternalError());

      for (unsigned int j = 0; j < i1.size(); ++j)
        for (unsigned int k = 0; k < i2.size(); ++k)
          if (std::fabs(M(j, k)) >= threshold)
            // Enter values into matrix only if j corresponds to a
            // degree of freedom on the refinement edge, k does
            // not, and both are not on the boundary. This is part
            // the difference between the complete matrix with no
            // boundary condition at the refinement edge and
            // the matrix assembled above by assemble().

            // Thus the logic is: enter the row if it is
            // constrained by hanging node constraints (actually,
            // the whole refinement edge), but not if it is
            // constrained by a boundary constraint.
            if (mg_constrained_dofs->at_refinement_edge(level, i1[j]) &&
                !mg_constrained_dofs->at_refinement_edge(level, i2[k]))
              {
                if ((!mg_constrained_dofs->is_boundary_index(level, i1[j]) &&
                     !mg_constrained_dofs->is_boundary_index(level, i2[k])) ||
                    (mg_constrained_dofs->is_boundary_index(level, i1[j]) &&
                     mg_constrained_dofs->is_boundary_index(level, i2[k]) &&
                     i1[j] == i2[k]))
                  G.add(i1[j], i2[k], M(j, k));
              }
    }


    template <typename MatrixType>
    inline void
    MGMatrixSimple<MatrixType>::assemble_out(
      MatrixType &                                G,
      const FullMatrix<double> &                  M,
      const std::vector<types::global_dof_index> &i1,
      const std::vector<types::global_dof_index> &i2,
      const unsigned int                          level)
    {
      AssertDimension(M.n(), i1.size());
      AssertDimension(M.m(), i2.size());
      Assert(mg_constrained_dofs != nullptr, ExcInternalError());

      for (unsigned int j = 0; j < i1.size(); ++j)
        for (unsigned int k = 0; k < i2.size(); ++k)
          if (std::fabs(M(k, j)) >= threshold)
            if (mg_constrained_dofs->at_refinement_edge(level, i1[j]) &&
                !mg_constrained_dofs->at_refinement_edge(level, i2[k]))
              {
                if ((!mg_constrained_dofs->is_boundary_index(level, i1[j]) &&
                     !mg_constrained_dofs->is_boundary_index(level, i2[k])) ||
                    (mg_constrained_dofs->is_boundary_index(level, i1[j]) &&
                     mg_constrained_dofs->is_boundary_index(level, i2[k]) &&
                     i1[j] == i2[k]))
                  G.add(i1[j], i2[k], M(k, j));
              }
    }


    template <typename MatrixType>
    template <class DOFINFO>
    inline void
    MGMatrixSimple<MatrixType>::assemble(const DOFINFO &info)
    {
      Assert(info.level_cell, ExcMessage("Cell must access level dofs"));
      const unsigned int level = info.cell->level();

      if (info.indices_by_block.size() == 0)
        {
          assemble((*matrix)[level],
                   info.matrix(0, false).matrix,
                   info.indices,
                   info.indices,
                   level);
          if (mg_constrained_dofs != nullptr)
            {
              assemble_in((*interface_in)[level],
                          info.matrix(0, false).matrix,
                          info.indices,
                          info.indices,
                          level);
              assemble_out((*interface_out)[level],
                           info.matrix(0, false).matrix,
                           info.indices,
                           info.indices,
                           level);
            }
        }
      else
        for (unsigned int k = 0; k < info.n_matrices(); ++k)
          {
            const unsigned int row    = info.matrix(k, false).row;
            const unsigned int column = info.matrix(k, false).column;

            assemble((*matrix)[level],
                     info.matrix(k, false).matrix,
                     info.indices_by_block[row],
                     info.indices_by_block[column],
                     level);

            if (mg_constrained_dofs != nullptr)
              {
                assemble_in((*interface_in)[level],
                            info.matrix(k, false).matrix,
                            info.indices_by_block[row],
                            info.indices_by_block[column],
                            level);
                assemble_out((*interface_out)[level],
                             info.matrix(k, false).matrix,
                             info.indices_by_block[column],
                             info.indices_by_block[row],
                             level);
              }
          }
    }


    template <typename MatrixType>
    template <class DOFINFO>
    inline void
    MGMatrixSimple<MatrixType>::assemble(const DOFINFO &info1,
                                         const DOFINFO &info2)
    {
      Assert(info1.level_cell, ExcMessage("Cell must access level dofs"));
      Assert(info2.level_cell, ExcMessage("Cell must access level dofs"));
      const unsigned int level1 = info1.cell->level();
      const unsigned int level2 = info2.cell->level();

      if (info1.indices_by_block.size() == 0)
        {
          if (level1 == level2)
            {
              assemble((*matrix)[level1],
                       info1.matrix(0, false).matrix,
                       info1.indices,
                       info1.indices,
                       level1);
              assemble((*matrix)[level1],
                       info1.matrix(0, true).matrix,
                       info1.indices,
                       info2.indices,
                       level1);
              assemble((*matrix)[level1],
                       info2.matrix(0, false).matrix,
                       info2.indices,
                       info2.indices,
                       level1);
              assemble((*matrix)[level1],
                       info2.matrix(0, true).matrix,
                       info2.indices,
                       info1.indices,
                       level1);
            }
          else
            {
              Assert(level1 > level2, ExcInternalError());
              // Do not add info2.M1,
              // which is done by
              // the coarser cell
              assemble((*matrix)[level1],
                       info1.matrix(0, false).matrix,
                       info1.indices,
                       info1.indices,
                       level1);
              if (level1 > 0)
                {
                  assemble_up((*flux_up)[level1],
                              info1.matrix(0, true).matrix,
                              info2.indices,
                              info1.indices,
                              level1);
                  assemble_down((*flux_down)[level1],
                                info2.matrix(0, true).matrix,
                                info2.indices,
                                info1.indices,
                                level1);
                }
            }
        }
      else
        for (unsigned int k = 0; k < info1.n_matrices(); ++k)
          {
            const unsigned int row    = info1.matrix(k, false).row;
            const unsigned int column = info1.matrix(k, false).column;

            if (level1 == level2)
              {
                assemble((*matrix)[level1],
                         info1.matrix(k, false).matrix,
                         info1.indices_by_block[row],
                         info1.indices_by_block[column],
                         level1);
                assemble((*matrix)[level1],
                         info1.matrix(k, true).matrix,
                         info1.indices_by_block[row],
                         info2.indices_by_block[column],
                         level1);
                assemble((*matrix)[level1],
                         info2.matrix(k, false).matrix,
                         info2.indices_by_block[row],
                         info2.indices_by_block[column],
                         level1);
                assemble((*matrix)[level1],
                         info2.matrix(k, true).matrix,
                         info2.indices_by_block[row],
                         info1.indices_by_block[column],
                         level1);
              }
            else
              {
                Assert(level1 > level2, ExcInternalError());
                // Do not add info2.M1,
                // which is done by
                // the coarser cell
                assemble((*matrix)[level1],
                         info1.matrix(k, false).matrix,
                         info1.indices_by_block[row],
                         info1.indices_by_block[column],
                         level1);
                if (level1 > 0)
                  {
                    assemble_up((*flux_up)[level1],
                                info1.matrix(k, true).matrix,
                                info2.indices_by_block[column],
                                info1.indices_by_block[row],
                                level1);
                    assemble_down((*flux_down)[level1],
                                  info2.matrix(k, true).matrix,
                                  info2.indices_by_block[row],
                                  info1.indices_by_block[column],
                                  level1);
                  }
              }
          }
    }

    //----------------------------------------------------------------------//

    template <typename MatrixType, typename VectorType>
    SystemSimple<MatrixType, VectorType>::SystemSimple(double t)
      : MatrixSimple<MatrixType>(t)
    {}


    template <typename MatrixType, typename VectorType>
    inline void
    SystemSimple<MatrixType, VectorType>::initialize(MatrixType &m,
                                                     VectorType &rhs)
    {
      AnyData     data;
      VectorType *p = &rhs;
      data.add(p, "right hand side");

      MatrixSimple<MatrixType>::initialize(m);
      ResidualSimple<VectorType>::initialize(data);
    }

    template <typename MatrixType, typename VectorType>
    inline void
    SystemSimple<MatrixType, VectorType>::initialize(
      const AffineConstraints<typename VectorType::value_type> &c)
    {
      ResidualSimple<VectorType>::initialize(c);
    }


    template <typename MatrixType, typename VectorType>
    template <class DOFINFO>
    inline void
    SystemSimple<MatrixType, VectorType>::initialize_info(DOFINFO &info,
                                                          bool     face) const
    {
      MatrixSimple<MatrixType>::initialize_info(info, face);
      ResidualSimple<VectorType>::initialize_info(info, face);
    }

    template <typename MatrixType, typename VectorType>
    inline void
    SystemSimple<MatrixType, VectorType>::assemble(
      const FullMatrix<double> &                  M,
      const Vector<double> &                      vector,
      const unsigned int                          index,
      const std::vector<types::global_dof_index> &indices)
    {
      AssertDimension(M.m(), indices.size());
      AssertDimension(M.n(), indices.size());

      AnyData     residuals = ResidualSimple<VectorType>::residuals;
      VectorType *v         = residuals.entry<VectorType *>(index);

      if (ResidualSimple<VectorType>::constraints == nullptr)
        {
          for (unsigned int i = 0; i < indices.size(); ++i)
            (*v)(indices[i]) += vector(i);

          for (unsigned int j = 0; j < indices.size(); ++j)
            for (unsigned int k = 0; k < indices.size(); ++k)
              if (std::fabs(M(j, k)) >= MatrixSimple<MatrixType>::threshold)
                MatrixSimple<MatrixType>::matrix[index]->add(indices[j],
                                                             indices[k],
                                                             M(j, k));
        }
      else
        {
          ResidualSimple<VectorType>::constraints->distribute_local_to_global(
            M,
            vector,
            indices,
            *MatrixSimple<MatrixType>::matrix[index],
            *v,
            true);
        }
    }

    template <typename MatrixType, typename VectorType>
    inline void
    SystemSimple<MatrixType, VectorType>::assemble(
      const FullMatrix<double> &                  M,
      const Vector<double> &                      vector,
      const unsigned int                          index,
      const std::vector<types::global_dof_index> &i1,
      const std::vector<types::global_dof_index> &i2)
    {
      AssertDimension(M.m(), i1.size());
      AssertDimension(M.n(), i2.size());

      AnyData     residuals = ResidualSimple<VectorType>::residuals;
      VectorType *v         = residuals.entry<VectorType *>(index);

      if (ResidualSimple<VectorType>::constraints == nullptr)
        {
          for (unsigned int j = 0; j < i1.size(); ++j)
            for (unsigned int k = 0; k < i2.size(); ++k)
              if (std::fabs(M(j, k)) >= MatrixSimple<MatrixType>::threshold)
                MatrixSimple<MatrixType>::matrix[index]->add(i1[j],
                                                             i2[k],
                                                             M(j, k));
        }
      else
        {
          ResidualSimple<VectorType>::constraints->distribute_local_to_global(
            vector, i1, i2, *v, M, false);
          ResidualSimple<VectorType>::constraints->distribute_local_to_global(
            M, i1, i2, *MatrixSimple<MatrixType>::matrix[index]);
        }
    }


    template <typename MatrixType, typename VectorType>
    template <class DOFINFO>
    inline void
    SystemSimple<MatrixType, VectorType>::assemble(const DOFINFO &info)
    {
      AssertDimension(MatrixSimple<MatrixType>::matrix.size(),
                      ResidualSimple<VectorType>::residuals.size());
      Assert(!info.level_cell, ExcMessage("Cell may not access level dofs"));
      const unsigned int n = info.indices_by_block.size();

      if (n == 0)
        {
          for (unsigned int m = 0; m < MatrixSimple<MatrixType>::matrix.size();
               ++m)
            assemble(info.matrix(m, false).matrix,
                     info.vector(m).block(0),
                     m,
                     info.indices);
        }
      else
        {
          for (unsigned int m = 0; m < MatrixSimple<MatrixType>::matrix.size();
               ++m)
            for (unsigned int k = 0; k < n * n; ++k)
              {
                const unsigned int row = info.matrix(k + m * n * n, false).row;
                const unsigned int column =
                  info.matrix(k + m * n * n, false).column;

                if (row == column)
                  assemble(info.matrix(k + m * n * n, false).matrix,
                           info.vector(m).block(row),
                           m,
                           info.indices_by_block[row]);
                else
                  assemble(info.matrix(k + m * n * n, false).matrix,
                           info.vector(m).block(row),
                           m,
                           info.indices_by_block[row],
                           info.indices_by_block[column]);
              }
        }
    }


    template <typename MatrixType, typename VectorType>
    template <class DOFINFO>
    inline void
    SystemSimple<MatrixType, VectorType>::assemble(const DOFINFO &info1,
                                                   const DOFINFO &info2)
    {
      Assert(!info1.level_cell, ExcMessage("Cell may not access level dofs"));
      Assert(!info2.level_cell, ExcMessage("Cell may not access level dofs"));
      AssertDimension(info1.indices_by_block.size(),
                      info2.indices_by_block.size());

      const unsigned int n = info1.indices_by_block.size();

      if (n == 0)
        {
          for (unsigned int m = 0; m < MatrixSimple<MatrixType>::matrix.size();
               ++m)
            {
              assemble(info1.matrix(m, false).matrix,
                       info1.vector(m).block(0),
                       m,
                       info1.indices);
              assemble(info1.matrix(m, true).matrix,
                       info1.vector(m).block(0),
                       m,
                       info1.indices,
                       info2.indices);
              assemble(info2.matrix(m, false).matrix,
                       info2.vector(m).block(0),
                       m,
                       info2.indices);
              assemble(info2.matrix(m, true).matrix,
                       info2.vector(m).block(0),
                       m,
                       info2.indices,
                       info1.indices);
            }
        }
      else
        {
          for (unsigned int m = 0; m < MatrixSimple<MatrixType>::matrix.size();
               ++m)
            for (unsigned int k = 0; k < n * n; ++k)
              {
                const unsigned int row = info1.matrix(k + m * n * n, false).row;
                const unsigned int column =
                  info1.matrix(k + m * n * n, false).column;

                if (row == column)
                  {
                    assemble(info1.matrix(k + m * n * n, false).matrix,
                             info1.vector(m).block(row),
                             m,
                             info1.indices_by_block[row]);
                    assemble(info2.matrix(k + m * n * n, false).matrix,
                             info2.vector(m).block(row),
                             m,
                             info2.indices_by_block[row]);
                  }
                else
                  {
                    assemble(info1.matrix(k + m * n * n, false).matrix,
                             info1.vector(m).block(row),
                             m,
                             info1.indices_by_block[row],
                             info1.indices_by_block[column]);
                    assemble(info2.matrix(k + m * n * n, false).matrix,
                             info2.vector(m).block(row),
                             m,
                             info2.indices_by_block[row],
                             info2.indices_by_block[column]);
                  }
                assemble(info1.matrix(k + m * n * n, true).matrix,
                         info1.vector(m).block(row),
                         m,
                         info1.indices_by_block[row],
                         info2.indices_by_block[column]);
                assemble(info2.matrix(k + m * n * n, true).matrix,
                         info2.vector(m).block(row),
                         m,
                         info2.indices_by_block[row],
                         info1.indices_by_block[column]);
              }
        }
    }
  } // namespace Assembler
} // namespace MeshWorker

DEAL_II_NAMESPACE_CLOSE

#endif


