//include/deal.II-translator/meshworker/assembler_0.txt
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


#ifndef dealii_mesh_worker_assembler_h
#define dealii_mesh_worker_assembler_h

#include <deal.II/base/config.h>

#include <deal.II/algorithms/any_data.h>

#include <deal.II/base/mg_level_object.h>
#include <deal.II/base/smartpointer.h>

#include <deal.II/lac/block_vector.h>

#include <deal.II/meshworker/dof_info.h>
#include <deal.II/meshworker/functional.h>
#include <deal.II/meshworker/simple.h>

#include <deal.II/multigrid/mg_constrained_dofs.h>


DEAL_II_NAMESPACE_OPEN

namespace MeshWorker
{
  /**
   * 含有对象的命名空间，可用于将单元格和面的计算数据组装成全局对象。这可以达到从收集单元和面的贡献的总误差估计到组装矩阵和多级矩阵的目的。
   * <h3>Data models</h3>
   * 从这个命名空间中选择的类决定了哪种数据模型被使用。
   * 对于局部和全局对象，我们可以选择两种模型。
   * <h4>The comprehensive data model</h4>
   * 这是由FESystem类建立的结构。在全球范围内，这意味着，数据被集合到一个残差向量和一个矩阵中。这些对象可能是块状向量和块状矩阵，但组装的过程忽略了这个事实。
   * 同样地，只有一个单元向量和单元矩阵，分别由FES系统的所有自由度来索引。
   * 在建立单元矩阵时，有必要区分系统的不同组成部分，并为每对单元选择正确的运算符。
   * <h4>The blocked data
   * model</h4>在这里，所有的块都是单独处理的（尽管在其他地方使用FESystem是为了其方便）。例如，没有组装区块矩阵，而是一个区块的列表，以后可以通过BlockMatrixArray来组合。在本地，这意味着，一个系统的每个矩阵块都是单独生成的，并被组装成相应的全局块。
   * 如果全局系统中每个块的矩阵数量不同，这种方法是有利的。例如，Oseen问题的块状预处理需要3个压力矩阵，但只有一个发散和一个速度的平流-扩散算子。
   * 此外，这种方法能够从每个方程和耦合算子的构建块中构建一个方程系统。
   * 然而，由于必须为每个基本元素创建一个单独的FEValues对象，所以先验地不太清楚哪种数据模型更有效。
   * @ingroup MeshWorker
   *
   */
  namespace Assembler
  {
    /**
     * 将局部残差组装成全局残差。
     * 全局残差被期望为一个FEVectors对象。本地残差是块向量。
     * 根据BlockInfo对象是否用 BlockInfo::initialize_local(),
     * 初始化，在本地使用全面或块数据模型。
     * 在块模型中，本地向量的每个块都对应于系统的单个块对这个单元的限制（见
     * @ref GlossBlock ）。
     * 因此，这个局部块的大小是FES系统中相应基元的自由度数。
     * @todo  综合模型目前没有实现。
     * @ingroup MeshWorker
     *
     */
    template <typename VectorType>
    class ResidualLocalBlocksToGlobalBlocks
    {
    public:
      /**
       * 将BlockInfo和矩阵指针复制到本地变量中。
       *
       */
      void
      initialize(const BlockInfo *block_info, AnyData &residuals);

      /**
       * 初始化约束。
       *
       */
      void
      initialize(
        const AffineConstraints<typename VectorType::value_type> &constraints);

      /**
       * 初始化以后用于装配的DoFInfo对象中的局部数据。
       * 如果 <code>!face</code> ， @p info
       * 对象指的是一个单元，否则指的是一个内部或边界面。
       *
       */
      template <class DOFINFO>
      void
      initialize_info(DOFINFO &info, bool face) const;


      /**
       * 将局部残差组装成全局残差。
       *
       */
      template <class DOFINFO>
      void
      assemble(const DOFINFO &info);

      /**
       * 将两个局部残差集合到全局残差中。
       *
       */
      template <class DOFINFO>
      void
      assemble(const DOFINFO &info1, const DOFINFO &info2);

    private:
      /**
       * 将一个局部残差装配到全局中。
       *
       */
      void
      assemble(VectorType &                                global,
               const BlockVector<double> &                 local,
               const std::vector<types::global_dof_index> &dof);

      /**
       * 全局向量，以AnyData容器的指针形式存储。
       *
       */
      AnyData residuals;

      /**
       * 一个指向包含块结构的对象的指针。
       *
       */
      SmartPointer<const BlockInfo,
                   ResidualLocalBlocksToGlobalBlocks<VectorType>>
        block_info;

      /**
       * 一个指向包含约束的对象的指针。
       *
       */
      SmartPointer<const AffineConstraints<typename VectorType::value_type>,
                   ResidualLocalBlocksToGlobalBlocks<VectorType>>
        constraints;
    };


    /**
     * 一个帮助类，将本地矩阵组装成全局矩阵。
     * 全局矩阵应该是MatrixBlock对象的一个向量，每个对象都包含一个矩阵对象，其函数对应于
     * SparseMatrix::add()
     * ，以及该矩阵在一个块系统中代表的块行和块列的信息。
     * 本地矩阵被期望为一个类似的MatrixBlock对象的向量，但包含一个FullMatrix。
     * 与ResidualLocalBlocksToGlobalBlocks一样，BlockInfo对象的初始化决定了是使用综合数据模型还是块模型。
     * 在综合模型中，每个LocalMatrixBlocks的坐标（0,0）和尺寸等于FES系统的自由度数。
     * 在综合模型中，每个块都有自己的块坐标，其大小取决于相关的
     * FESystem::base_element().
     * 这些块可以单独生成，并将由这个对象组装成正确的矩阵块。
     * @ingroup MeshWorker
     *
     */
    template <typename MatrixType, typename number = double>
    class MatrixLocalBlocksToGlobalBlocks
    {
    public:
      /**
       * 构造函数，初始化#阈值，它限制了可以输入矩阵的小数字。
       *
       */
      MatrixLocalBlocksToGlobalBlocks(double threshold = 1.e-12);

      /**
       * 将BlockInfo和矩阵指针复制到局部变量中，并初始化单元矩阵向量。
       *
       */
      void
      initialize(const BlockInfo *              block_info,
                 MatrixBlockVector<MatrixType> &matrices);

      /**
       * 初始化约束。
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
       * 将局部矩阵组装成全局矩阵。
       *
       */
      template <class DOFINFO>
      void
      assemble(const DOFINFO &info);

      /**
       * 将所有局部矩阵集合到全局矩阵中。
       *
       */
      template <class DOFINFO>
      void
      assemble(const DOFINFO &info1, const DOFINFO &info2);

    private:
      /**
       * 将单个局部矩阵组合成全局矩阵。
       *
       */
      void
      assemble(MatrixBlock<MatrixType> &                   global,
               const FullMatrix<number> &                  local,
               const unsigned int                          block_row,
               const unsigned int                          block_col,
               const std::vector<types::global_dof_index> &dof1,
               const std::vector<types::global_dof_index> &dof2);

      /**
       * 全局矩阵，以指针向量的形式存储。
       *
       */
      SmartPointer<MatrixBlockVector<MatrixType>,
                   MatrixLocalBlocksToGlobalBlocks<MatrixType, number>>
        matrices;

      /**
       * 一个指向包含块结构的对象的指针。
       *
       */
      SmartPointer<const BlockInfo,
                   MatrixLocalBlocksToGlobalBlocks<MatrixType, number>>
        block_info;
      /**
       * 一个指向包含约束的对象的指针。
       *
       */

      SmartPointer<const AffineConstraints<typename MatrixType::value_type>,
                   MatrixLocalBlocksToGlobalBlocks<MatrixType, number>>
        constraints;

      /**
       * 将被输入全局矩阵的最小的正数。所有更小的绝对值将被视为零，不会被集合。
       *
       */
      const double threshold;
    };

    /**
     * 一个帮助类，将局部矩阵组装成全局多级矩阵。这个类是MatrixLocalBlocksToGlobalBlocks的多级等价物，该类的文档在很大程度上适用于此。
     * 全局矩阵被期望为指向MatrixBlock对象的一个向量，每个对象都包含一个MGLevelObject，其中的矩阵具有对应于
     * SparseMatrix::add()
     * 的函数，以及该矩阵在块系统中代表的块行和块列的信息。
     * 本地矩阵是一个类似MatrixBlock对象的向量，但包含一个FullMatrix。
     * 如果发生局部细化，多重网格方法需要更多的矩阵，两个用于连续元素，另外两个如果在界面上计算数值通量。第二组矩阵可以通过initialize_edge_flux()添加。一旦加入，所有参与矩阵的贡献将从单元和面的矩阵中自动集合起来。
     * @ingroup MeshWorker
     *
     */
    template <typename MatrixType, typename number = double>
    class MGMatrixLocalBlocksToGlobalBlocks
    {
    public:
      using MatrixPtrVector = MGMatrixBlockVector<MatrixType>;
      using MatrixPtrVectorPtr =
        SmartPointer<MatrixPtrVector,
                     MGMatrixLocalBlocksToGlobalBlocks<MatrixType, number>>;

      /**
       * 构造函数，初始化#threshold，它限制了可以输入矩阵的小数字。
       *
       */
      MGMatrixLocalBlocksToGlobalBlocks(double threshold = 1.e-12);

      /**
       * 将BlockInfo和矩阵指针复制到本地变量，并初始化单元格矩阵向量。
       *
       */
      void
      initialize(const BlockInfo *block_info, MatrixPtrVector &matrices);

      /**
       * 初始化多级约束。
       *
       */
      void
      initialize(const MGConstrainedDoFs &mg_constrained_dofs);

      /**
       * 局部细化网格上的多网格方法需要额外的矩阵。
       * 对于不连续的Galerkin方法，这些是跨越细化边缘的两个通量矩阵，由该方法设置。
       *
       */
      void
      initialize_edge_flux(MatrixPtrVector &up, MatrixPtrVector &down);

      /**
       * 局部细化网格上的多栅方法需要额外的矩阵。
       * 对于不连续的Galerkin方法，这些是跨越细化边缘的两个通量矩阵，由该方法设置。
       *
       */
      void
      initialize_interfaces(MatrixPtrVector &interface_in,
                            MatrixPtrVector &interface_out);
      /**
       * 初始化以后用于装配的DoFInfo对象中的局部数据。
       * 如果 <code>!face</code> ， @p info
       * 对象指的是一个单元，否则指的是一个内部或边界面。
       *
       */
      template <class DOFINFO>
      void
      initialize_info(DOFINFO &info, bool face) const;


      /**
       * 将局部矩阵组装成全局矩阵。
       *
       */
      template <class DOFINFO>
      void
      assemble(const DOFINFO &info);

      /**
       * 将所有局部矩阵集合到全局矩阵中。
       *
       */
      template <class DOFINFO>
      void
      assemble(const DOFINFO &info1, const DOFINFO &info2);

    private:
      /**
       * 将单个局部矩阵组合成全局矩阵。
       *
       */
      void
      assemble(MatrixType &                                global,
               const FullMatrix<number> &                  local,
               const unsigned int                          block_row,
               const unsigned int                          block_col,
               const std::vector<types::global_dof_index> &dof1,
               const std::vector<types::global_dof_index> &dof2,
               const unsigned int                          level1,
               const unsigned int                          level2,
               bool                                        transpose = false);

      /**
       * 将一个单一的局部矩阵组合成一个全局矩阵。
       *
       */
      void
      assemble_fluxes(MatrixType &                                global,
                      const FullMatrix<number> &                  local,
                      const unsigned int                          block_row,
                      const unsigned int                          block_col,
                      const std::vector<types::global_dof_index> &dof1,
                      const std::vector<types::global_dof_index> &dof2,
                      const unsigned int                          level1,
                      const unsigned int                          level2);

      /**
       * 将一个单一的本地矩阵组合成一个全局矩阵。
       *
       */
      void
      assemble_up(MatrixType &                                global,
                  const FullMatrix<number> &                  local,
                  const unsigned int                          block_row,
                  const unsigned int                          block_col,
                  const std::vector<types::global_dof_index> &dof1,
                  const std::vector<types::global_dof_index> &dof2,
                  const unsigned int                          level1,
                  const unsigned int                          level2);

      /**
       * 将一个单一的本地矩阵组合成一个全局矩阵。
       *
       */
      void
      assemble_down(MatrixType &                                global,
                    const FullMatrix<number> &                  local,
                    const unsigned int                          block_row,
                    const unsigned int                          block_col,
                    const std::vector<types::global_dof_index> &dof1,
                    const std::vector<types::global_dof_index> &dof2,
                    const unsigned int                          level1,
                    const unsigned int                          level2);

      /**
       * 将一个单一的本地矩阵组合成一个全局矩阵。
       *
       */
      void
      assemble_in(MatrixType &                                global,
                  const FullMatrix<number> &                  local,
                  const unsigned int                          block_row,
                  const unsigned int                          block_col,
                  const std::vector<types::global_dof_index> &dof1,
                  const std::vector<types::global_dof_index> &dof2,
                  const unsigned int                          level1,
                  const unsigned int                          level2);

      /**
       * 将一个单一的本地矩阵组合成一个全局矩阵。
       *
       */
      void
      assemble_out(MatrixType &                                global,
                   const FullMatrix<number> &                  local,
                   const unsigned int                          block_row,
                   const unsigned int                          block_col,
                   const std::vector<types::global_dof_index> &dof1,
                   const std::vector<types::global_dof_index> &dof2,
                   const unsigned int                          level1,
                   const unsigned int                          level2);

      /**
       * 水平矩阵，存储为一个指向性的向量。
       *
       */
      MatrixPtrVectorPtr matrices;

      /**
       * 细化边缘的细化水平和粗化水平之间的通量矩阵。
       *
       */
      MatrixPtrVectorPtr flux_down;

      /**
       * 细化边缘的粗级和细级之间的通量矩阵。
       *
       */
      MatrixPtrVectorPtr flux_up;

      /**
       * 细化边缘处精细级和粗化级之间的界面矩阵。
       *
       */
      MatrixPtrVectorPtr interface_out;

      /**
       * 细化边缘处的粗级和细级之间的界面矩阵。
       *
       */
      MatrixPtrVectorPtr interface_in;

      /**
       * 一个指向包含块结构的对象的指针。
       *
       */
      SmartPointer<const BlockInfo,
                   MGMatrixLocalBlocksToGlobalBlocks<MatrixType, number>>
        block_info;

      /**
       * 一个指向包含约束条件的对象的指针。
       *
       */
      SmartPointer<const MGConstrainedDoFs,
                   MGMatrixLocalBlocksToGlobalBlocks<MatrixType, number>>
        mg_constrained_dofs;


      /**
       * 将被输入全局矩阵的最小的正数。所有更小的绝对值将被视为零，不会被集合。
       *
       */
      const double threshold;
    };

    //----------------------------------------------------------------------//

    template <typename VectorType>
    inline void
    ResidualLocalBlocksToGlobalBlocks<VectorType>::initialize(
      const BlockInfo *b,
      AnyData &        m)
    {
      block_info = b;
      residuals  = m;
    }

    template <typename VectorType>
    inline void
    ResidualLocalBlocksToGlobalBlocks<VectorType>::initialize(
      const AffineConstraints<typename VectorType::value_type> &c)
    {
      constraints = &c;
    }


    template <typename VectorType>
    template <class DOFINFO>
    inline void
    ResidualLocalBlocksToGlobalBlocks<VectorType>::initialize_info(
      DOFINFO &info,
      bool) const
    {
      info.initialize_vectors(residuals.size());
    }

    template <typename VectorType>
    inline void
    ResidualLocalBlocksToGlobalBlocks<VectorType>::assemble(
      VectorType &                                global,
      const BlockVector<double> &                 local,
      const std::vector<types::global_dof_index> &dof)
    {
      if (constraints == 0)
        {
          for (unsigned int b = 0; b < local.n_blocks(); ++b)
            for (unsigned int j = 0; j < local.block(b).size(); ++j)
              {
                // The coordinates of
                // the current entry in
                // DoFHandler
                // numbering, which
                // differs from the
                // block-wise local
                // numbering we use in
                // our local vectors
                const unsigned int jcell =
                  this->block_info->local().local_to_global(b, j);
                global(dof[jcell]) += local.block(b)(j);
              }
        }
      else
        constraints->distribute_local_to_global(local, dof, global);
    }


    template <typename VectorType>
    template <class DOFINFO>
    inline void
    ResidualLocalBlocksToGlobalBlocks<VectorType>::assemble(const DOFINFO &info)
    {
      for (unsigned int i = 0; i < residuals.size(); ++i)
        assemble(*(residuals.entry<VectorType>(i)),
                 info.vector(i),
                 info.indices);
    }


    template <typename VectorType>
    template <class DOFINFO>
    inline void
    ResidualLocalBlocksToGlobalBlocks<VectorType>::assemble(
      const DOFINFO &info1,
      const DOFINFO &info2)
    {
      for (unsigned int i = 0; i < residuals.size(); ++i)
        {
          assemble(*(residuals.entry<VectorType>(i)),
                   info1.vector(i),
                   info1.indices);
          assemble(*(residuals.entry<VectorType>(i)),
                   info2.vector(i),
                   info2.indices);
        }
    }


    //----------------------------------------------------------------------//

    template <typename MatrixType, typename number>
    inline MatrixLocalBlocksToGlobalBlocks<MatrixType, number>::
      MatrixLocalBlocksToGlobalBlocks(double threshold)
      : threshold(threshold)
    {}


    template <typename MatrixType, typename number>
    inline void
    MatrixLocalBlocksToGlobalBlocks<MatrixType, number>::initialize(
      const BlockInfo *              b,
      MatrixBlockVector<MatrixType> &m)
    {
      block_info = b;
      matrices   = &m;
    }



    template <typename MatrixType, typename number>
    inline void
    MatrixLocalBlocksToGlobalBlocks<MatrixType, number>::initialize(
      const AffineConstraints<typename MatrixType::value_type> &c)
    {
      constraints = &c;
    }



    template <typename MatrixType, typename number>
    template <class DOFINFO>
    inline void
    MatrixLocalBlocksToGlobalBlocks<MatrixType, number>::initialize_info(
      DOFINFO &info,
      bool     face) const
    {
      info.initialize_matrices(*matrices, face);
    }



    template <typename MatrixType, typename number>
    inline void
    MatrixLocalBlocksToGlobalBlocks<MatrixType, number>::assemble(
      MatrixBlock<MatrixType> &                   global,
      const FullMatrix<number> &                  local,
      const unsigned int                          block_row,
      const unsigned int                          block_col,
      const std::vector<types::global_dof_index> &dof1,
      const std::vector<types::global_dof_index> &dof2)
    {
      if (constraints == nullptr)
        {
          for (unsigned int j = 0; j < local.n_rows(); ++j)
            for (unsigned int k = 0; k < local.n_cols(); ++k)
              if (std::fabs(local(j, k)) >= threshold)
                {
                  // The coordinates of
                  // the current entry in
                  // DoFHandler
                  // numbering, which
                  // differs from the
                  // block-wise local
                  // numbering we use in
                  // our local matrices
                  const unsigned int jcell =
                    this->block_info->local().local_to_global(block_row, j);
                  const unsigned int kcell =
                    this->block_info->local().local_to_global(block_col, k);

                  global.add(dof1[jcell], dof2[kcell], local(j, k));
                }
        }
      else
        {
          const BlockIndices &                 bi = this->block_info->local();
          std::vector<types::global_dof_index> sliced_row_indices(
            bi.block_size(block_row));
          for (unsigned int i = 0; i < sliced_row_indices.size(); ++i)
            sliced_row_indices[i] = dof1[bi.block_start(block_row) + i];

          std::vector<types::global_dof_index> sliced_col_indices(
            bi.block_size(block_col));
          for (unsigned int i = 0; i < sliced_col_indices.size(); ++i)
            sliced_col_indices[i] = dof2[bi.block_start(block_col) + i];

          constraints->distribute_local_to_global(local,
                                                  sliced_row_indices,
                                                  sliced_col_indices,
                                                  global);
        }
    }


    template <typename MatrixType, typename number>
    template <class DOFINFO>
    inline void
    MatrixLocalBlocksToGlobalBlocks<MatrixType, number>::assemble(
      const DOFINFO &info)
    {
      for (unsigned int i = 0; i < matrices->size(); ++i)
        {
          // Row and column index of
          // the block we are dealing with
          const types::global_dof_index row = matrices->block(i).row;
          const types::global_dof_index col = matrices->block(i).column;

          assemble(matrices->block(i),
                   info.matrix(i, false).matrix,
                   row,
                   col,
                   info.indices,
                   info.indices);
        }
    }


    template <typename MatrixType, typename number>
    template <class DOFINFO>
    inline void
    MatrixLocalBlocksToGlobalBlocks<MatrixType, number>::assemble(
      const DOFINFO &info1,
      const DOFINFO &info2)
    {
      for (unsigned int i = 0; i < matrices->size(); ++i)
        {
          // Row and column index of
          // the block we are dealing with
          const types::global_dof_index row = matrices->block(i).row;
          const types::global_dof_index col = matrices->block(i).column;

          assemble(matrices->block(i),
                   info1.matrix(i, false).matrix,
                   row,
                   col,
                   info1.indices,
                   info1.indices);
          assemble(matrices->block(i),
                   info1.matrix(i, true).matrix,
                   row,
                   col,
                   info1.indices,
                   info2.indices);
          assemble(matrices->block(i),
                   info2.matrix(i, false).matrix,
                   row,
                   col,
                   info2.indices,
                   info2.indices);
          assemble(matrices->block(i),
                   info2.matrix(i, true).matrix,
                   row,
                   col,
                   info2.indices,
                   info1.indices);
        }
    }


    // ----------------------------------------------------------------------//

    template <typename MatrixType, typename number>
    inline MGMatrixLocalBlocksToGlobalBlocks<MatrixType, number>::
      MGMatrixLocalBlocksToGlobalBlocks(double threshold)
      : threshold(threshold)
    {}


    template <typename MatrixType, typename number>
    inline void
    MGMatrixLocalBlocksToGlobalBlocks<MatrixType, number>::initialize(
      const BlockInfo *b,
      MatrixPtrVector &m)
    {
      block_info = b;
      AssertDimension(block_info->local().size(), block_info->global().size());
      matrices = &m;
    }


    template <typename MatrixType, typename number>
    inline void
    MGMatrixLocalBlocksToGlobalBlocks<MatrixType, number>::initialize(
      const MGConstrainedDoFs &mg_c)
    {
      mg_constrained_dofs = &mg_c;
    }


    template <typename MatrixType, typename number>
    template <class DOFINFO>
    inline void
    MGMatrixLocalBlocksToGlobalBlocks<MatrixType, number>::initialize_info(
      DOFINFO &info,
      bool     face) const
    {
      info.initialize_matrices(*matrices, face);
    }



    template <typename MatrixType, typename number>
    inline void
    MGMatrixLocalBlocksToGlobalBlocks<MatrixType, number>::initialize_edge_flux(
      MatrixPtrVector &up,
      MatrixPtrVector &down)
    {
      flux_up   = up;
      flux_down = down;
    }


    template <typename MatrixType, typename number>
    inline void
    MGMatrixLocalBlocksToGlobalBlocks<MatrixType, number>::
      initialize_interfaces(MatrixPtrVector &in, MatrixPtrVector &out)
    {
      interface_in  = in;
      interface_out = out;
    }


    template <typename MatrixType, typename number>
    inline void
    MGMatrixLocalBlocksToGlobalBlocks<MatrixType, number>::assemble(
      MatrixType &                                global,
      const FullMatrix<number> &                  local,
      const unsigned int                          block_row,
      const unsigned int                          block_col,
      const std::vector<types::global_dof_index> &dof1,
      const std::vector<types::global_dof_index> &dof2,
      const unsigned int                          level1,
      const unsigned int                          level2,
      bool                                        transpose)
    {
      for (unsigned int j = 0; j < local.n_rows(); ++j)
        for (unsigned int k = 0; k < local.n_cols(); ++k)
          if (std::fabs(local(j, k)) >= threshold)
            {
              // The coordinates of
              // the current entry in
              // DoFHandler
              // numbering, which
              // differs from the
              // block-wise local
              // numbering we use in
              // our local matrices
              const unsigned int jcell =
                this->block_info->local().local_to_global(block_row, j);
              const unsigned int kcell =
                this->block_info->local().local_to_global(block_col, k);

              // The global dof
              // indices to assemble
              // in. Since we may
              // have face matrices
              // coupling two
              // different cells, we
              // provide two sets of
              // dof indices.
              const unsigned int jglobal = this->block_info->level(level1)
                                             .global_to_local(dof1[jcell])
                                             .second;
              const unsigned int kglobal = this->block_info->level(level2)
                                             .global_to_local(dof2[kcell])
                                             .second;

              if (mg_constrained_dofs == 0)
                {
                  if (transpose)
                    global.add(kglobal, jglobal, local(j, k));
                  else
                    global.add(jglobal, kglobal, local(j, k));
                }
              else
                {
                  if (!mg_constrained_dofs->at_refinement_edge(level1,
                                                               jglobal) &&
                      !mg_constrained_dofs->at_refinement_edge(level2, kglobal))
                    {
                      if (mg_constrained_dofs->set_boundary_values())
                        {
                          if ((!mg_constrained_dofs->is_boundary_index(
                                 level1, jglobal) &&
                               !mg_constrained_dofs->is_boundary_index(
                                 level2, kglobal)) ||
                              (mg_constrained_dofs->is_boundary_index(
                                 level1, jglobal) &&
                               mg_constrained_dofs->is_boundary_index(
                                 level2, kglobal) &&
                               jglobal == kglobal))
                            {
                              if (transpose)
                                global.add(kglobal, jglobal, local(j, k));
                              else
                                global.add(jglobal, kglobal, local(j, k));
                            }
                        }
                      else
                        {
                          if (transpose)
                            global.add(kglobal, jglobal, local(j, k));
                          else
                            global.add(jglobal, kglobal, local(j, k));
                        }
                    }
                }
            }
    }


    template <typename MatrixType, typename number>
    inline void
    MGMatrixLocalBlocksToGlobalBlocks<MatrixType, number>::assemble_fluxes(
      MatrixType &                                global,
      const FullMatrix<number> &                  local,
      const unsigned int                          block_row,
      const unsigned int                          block_col,
      const std::vector<types::global_dof_index> &dof1,
      const std::vector<types::global_dof_index> &dof2,
      const unsigned int                          level1,
      const unsigned int                          level2)
    {
      for (unsigned int j = 0; j < local.n_rows(); ++j)
        for (unsigned int k = 0; k < local.n_cols(); ++k)
          if (std::fabs(local(j, k)) >= threshold)
            {
              // The coordinates of
              // the current entry in
              // DoFHandler
              // numbering, which
              // differs from the
              // block-wise local
              // numbering we use in
              // our local matrices
              const unsigned int jcell =
                this->block_info->local().local_to_global(block_row, j);
              const unsigned int kcell =
                this->block_info->local().local_to_global(block_col, k);

              // The global dof
              // indices to assemble
              // in. Since we may
              // have face matrices
              // coupling two
              // different cells, we
              // provide two sets of
              // dof indices.
              const unsigned int jglobal = this->block_info->level(level1)
                                             .global_to_local(dof1[jcell])
                                             .second;
              const unsigned int kglobal = this->block_info->level(level2)
                                             .global_to_local(dof2[kcell])
                                             .second;

              if (mg_constrained_dofs == 0)
                global.add(jglobal, kglobal, local(j, k));
              else
                {
                  if (!mg_constrained_dofs->non_refinement_edge_index(
                        level1, jglobal) &&
                      !mg_constrained_dofs->non_refinement_edge_index(level2,
                                                                      kglobal))
                    {
                      if (!mg_constrained_dofs->at_refinement_edge(level1,
                                                                   jglobal) &&
                          !mg_constrained_dofs->at_refinement_edge(level2,
                                                                   kglobal))
                        global.add(jglobal, kglobal, local(j, k));
                    }
                }
            }
    }

    template <typename MatrixType, typename number>
    inline void
    MGMatrixLocalBlocksToGlobalBlocks<MatrixType, number>::assemble_up(
      MatrixType &                                global,
      const FullMatrix<number> &                  local,
      const unsigned int                          block_row,
      const unsigned int                          block_col,
      const std::vector<types::global_dof_index> &dof1,
      const std::vector<types::global_dof_index> &dof2,
      const unsigned int                          level1,
      const unsigned int                          level2)
    {
      for (unsigned int j = 0; j < local.n_rows(); ++j)
        for (unsigned int k = 0; k < local.n_cols(); ++k)
          if (std::fabs(local(j, k)) >= threshold)
            {
              // The coordinates of
              // the current entry in
              // DoFHandler
              // numbering, which
              // differs from the
              // block-wise local
              // numbering we use in
              // our local matrices
              const unsigned int jcell =
                this->block_info->local().local_to_global(block_row, j);
              const unsigned int kcell =
                this->block_info->local().local_to_global(block_col, k);

              // The global dof
              // indices to assemble
              // in. Since we may
              // have face matrices
              // coupling two
              // different cells, we
              // provide two sets of
              // dof indices.
              const unsigned int jglobal = this->block_info->level(level1)
                                             .global_to_local(dof1[jcell])
                                             .second;
              const unsigned int kglobal = this->block_info->level(level2)
                                             .global_to_local(dof2[kcell])
                                             .second;

              if (mg_constrained_dofs == 0)
                global.add(jglobal, kglobal, local(j, k));
              else
                {
                  if (!mg_constrained_dofs->non_refinement_edge_index(
                        level1, jglobal) &&
                      !mg_constrained_dofs->non_refinement_edge_index(level2,
                                                                      kglobal))
                    {
                      if (!mg_constrained_dofs->at_refinement_edge(level1,
                                                                   jglobal) &&
                          !mg_constrained_dofs->at_refinement_edge(level2,
                                                                   kglobal))
                        global.add(jglobal, kglobal, local(j, k));
                    }
                }
            }
    }

    template <typename MatrixType, typename number>
    inline void
    MGMatrixLocalBlocksToGlobalBlocks<MatrixType, number>::assemble_down(
      MatrixType &                                global,
      const FullMatrix<number> &                  local,
      const unsigned int                          block_row,
      const unsigned int                          block_col,
      const std::vector<types::global_dof_index> &dof1,
      const std::vector<types::global_dof_index> &dof2,
      const unsigned int                          level1,
      const unsigned int                          level2)
    {
      for (unsigned int j = 0; j < local.n_rows(); ++j)
        for (unsigned int k = 0; k < local.n_cols(); ++k)
          if (std::fabs(local(k, j)) >= threshold)
            {
              // The coordinates of
              // the current entry in
              // DoFHandler
              // numbering, which
              // differs from the
              // block-wise local
              // numbering we use in
              // our local matrices
              const unsigned int jcell =
                this->block_info->local().local_to_global(block_row, j);
              const unsigned int kcell =
                this->block_info->local().local_to_global(block_col, k);

              // The global dof
              // indices to assemble
              // in. Since we may
              // have face matrices
              // coupling two
              // different cells, we
              // provide two sets of
              // dof indices.
              const unsigned int jglobal = this->block_info->level(level1)
                                             .global_to_local(dof1[jcell])
                                             .second;
              const unsigned int kglobal = this->block_info->level(level2)
                                             .global_to_local(dof2[kcell])
                                             .second;

              if (mg_constrained_dofs == 0)
                global.add(jglobal, kglobal, local(k, j));
              else
                {
                  if (!mg_constrained_dofs->non_refinement_edge_index(
                        level1, jglobal) &&
                      !mg_constrained_dofs->non_refinement_edge_index(level2,
                                                                      kglobal))
                    {
                      if (!mg_constrained_dofs->at_refinement_edge(level1,
                                                                   jglobal) &&
                          !mg_constrained_dofs->at_refinement_edge(level2,
                                                                   kglobal))
                        global.add(jglobal, kglobal, local(k, j));
                    }
                }
            }
    }

    template <typename MatrixType, typename number>
    inline void
    MGMatrixLocalBlocksToGlobalBlocks<MatrixType, number>::assemble_in(
      MatrixType &                                global,
      const FullMatrix<number> &                  local,
      const unsigned int                          block_row,
      const unsigned int                          block_col,
      const std::vector<types::global_dof_index> &dof1,
      const std::vector<types::global_dof_index> &dof2,
      const unsigned int                          level1,
      const unsigned int                          level2)
    {
      //      AssertDimension(local.n(), dof1.size());
      //      AssertDimension(local.m(), dof2.size());

      for (unsigned int j = 0; j < local.n_rows(); ++j)
        for (unsigned int k = 0; k < local.n_cols(); ++k)
          if (std::fabs(local(j, k)) >= threshold)
            {
              // The coordinates of
              // the current entry in
              // DoFHandler
              // numbering, which
              // differs from the
              // block-wise local
              // numbering we use in
              // our local matrices
              const unsigned int jcell =
                this->block_info->local().local_to_global(block_row, j);
              const unsigned int kcell =
                this->block_info->local().local_to_global(block_col, k);

              // The global dof
              // indices to assemble
              // in. Since we may
              // have face matrices
              // coupling two
              // different cells, we
              // provide two sets of
              // dof indices.
              const unsigned int jglobal = this->block_info->level(level1)
                                             .global_to_local(dof1[jcell])
                                             .second;
              const unsigned int kglobal = this->block_info->level(level2)
                                             .global_to_local(dof2[kcell])
                                             .second;

              if (mg_constrained_dofs == 0)
                global.add(jglobal, kglobal, local(j, k));
              else
                {
                  if (mg_constrained_dofs->at_refinement_edge(level1,
                                                              jglobal) &&
                      !mg_constrained_dofs->at_refinement_edge(level2, kglobal))
                    {
                      if (mg_constrained_dofs->set_boundary_values())
                        {
                          if ((!mg_constrained_dofs->is_boundary_index(
                                 level1, jglobal) &&
                               !mg_constrained_dofs->is_boundary_index(
                                 level2, kglobal)) ||
                              (mg_constrained_dofs->is_boundary_index(
                                 level1, jglobal) &&
                               mg_constrained_dofs->is_boundary_index(
                                 level2, kglobal) &&
                               jglobal == kglobal))
                            global.add(jglobal, kglobal, local(j, k));
                        }
                      else
                        global.add(jglobal, kglobal, local(j, k));
                    }
                }
            }
    }

    template <typename MatrixType, typename number>
    inline void
    MGMatrixLocalBlocksToGlobalBlocks<MatrixType, number>::assemble_out(
      MatrixType &                                global,
      const FullMatrix<number> &                  local,
      const unsigned int                          block_row,
      const unsigned int                          block_col,
      const std::vector<types::global_dof_index> &dof1,
      const std::vector<types::global_dof_index> &dof2,
      const unsigned int                          level1,
      const unsigned int                          level2)
    {
      //      AssertDimension(local.n(), dof1.size());
      //      AssertDimension(local.m(), dof2.size());

      for (unsigned int j = 0; j < local.n_rows(); ++j)
        for (unsigned int k = 0; k < local.n_cols(); ++k)
          if (std::fabs(local(k, j)) >= threshold)
            {
              // The coordinates of
              // the current entry in
              // DoFHandler
              // numbering, which
              // differs from the
              // block-wise local
              // numbering we use in
              // our local matrices
              const unsigned int jcell =
                this->block_info->local().local_to_global(block_row, j);
              const unsigned int kcell =
                this->block_info->local().local_to_global(block_col, k);

              // The global dof
              // indices to assemble
              // in. Since we may
              // have face matrices
              // coupling two
              // different cells, we
              // provide two sets of
              // dof indices.
              const unsigned int jglobal = this->block_info->level(level1)
                                             .global_to_local(dof1[jcell])
                                             .second;
              const unsigned int kglobal = this->block_info->level(level2)
                                             .global_to_local(dof2[kcell])
                                             .second;

              if (mg_constrained_dofs == 0)
                global.add(jglobal, kglobal, local(k, j));
              else
                {
                  if (mg_constrained_dofs->at_refinement_edge(level1,
                                                              jglobal) &&
                      !mg_constrained_dofs->at_refinement_edge(level2, kglobal))
                    {
                      if (mg_constrained_dofs->set_boundary_values())
                        {
                          if ((!mg_constrained_dofs->is_boundary_index(
                                 level1, jglobal) &&
                               !mg_constrained_dofs->is_boundary_index(
                                 level2, kglobal)) ||
                              (mg_constrained_dofs->is_boundary_index(
                                 level1, jglobal) &&
                               mg_constrained_dofs->is_boundary_index(
                                 level2, kglobal) &&
                               jglobal == kglobal))
                            global.add(jglobal, kglobal, local(k, j));
                        }
                      else
                        global.add(jglobal, kglobal, local(k, j));
                    }
                }
            }
    }


    template <typename MatrixType, typename number>
    template <class DOFINFO>
    inline void
    MGMatrixLocalBlocksToGlobalBlocks<MatrixType, number>::assemble(
      const DOFINFO &info)
    {
      const unsigned int level = info.cell->level();

      for (unsigned int i = 0; i < matrices->size(); ++i)
        {
          // Row and column index of
          // the block we are dealing with
          const unsigned int row = matrices->block(i)[level].row;
          const unsigned int col = matrices->block(i)[level].column;

          assemble(matrices->block(i)[level].matrix,
                   info.matrix(i, false).matrix,
                   row,
                   col,
                   info.indices,
                   info.indices,
                   level,
                   level);
          if (mg_constrained_dofs != 0)
            {
              if (interface_in != 0)
                assemble_in(interface_in->block(i)[level],
                            info.matrix(i, false).matrix,
                            row,
                            col,
                            info.indices,
                            info.indices,
                            level,
                            level);
              if (interface_out != 0)
                assemble_in(interface_out->block(i)[level],
                            info.matrix(i, false).matrix,
                            row,
                            col,
                            info.indices,
                            info.indices,
                            level,
                            level);

              assemble_in(matrices->block_in(i)[level],
                          info.matrix(i, false).matrix,
                          row,
                          col,
                          info.indices,
                          info.indices,
                          level,
                          level);
              assemble_out(matrices->block_out(i)[level],
                           info.matrix(i, false).matrix,
                           row,
                           col,
                           info.indices,
                           info.indices,
                           level,
                           level);
            }
        }
    }


    template <typename MatrixType, typename number>
    template <class DOFINFO>
    inline void
    MGMatrixLocalBlocksToGlobalBlocks<MatrixType, number>::assemble(
      const DOFINFO &info1,
      const DOFINFO &info2)
    {
      const unsigned int level1 = info1.cell->level();
      const unsigned int level2 = info2.cell->level();

      for (unsigned int i = 0; i < matrices->size(); ++i)
        {
          MGLevelObject<MatrixBlock<MatrixType>> &o = matrices->block(i);

          // Row and column index of
          // the block we are dealing with
          const unsigned int row = o[level1].row;
          const unsigned int col = o[level1].column;

          if (level1 == level2)
            {
              if (mg_constrained_dofs == 0)
                {
                  assemble(o[level1].matrix,
                           info1.matrix(i, false).matrix,
                           row,
                           col,
                           info1.indices,
                           info1.indices,
                           level1,
                           level1);
                  assemble(o[level1].matrix,
                           info1.matrix(i, true).matrix,
                           row,
                           col,
                           info1.indices,
                           info2.indices,
                           level1,
                           level2);
                  assemble(o[level1].matrix,
                           info2.matrix(i, false).matrix,
                           row,
                           col,
                           info2.indices,
                           info2.indices,
                           level2,
                           level2);
                  assemble(o[level1].matrix,
                           info2.matrix(i, true).matrix,
                           row,
                           col,
                           info2.indices,
                           info1.indices,
                           level2,
                           level1);
                }
              else
                {
                  assemble_fluxes(o[level1].matrix,
                                  info1.matrix(i, false).matrix,
                                  row,
                                  col,
                                  info1.indices,
                                  info1.indices,
                                  level1,
                                  level1);
                  assemble_fluxes(o[level1].matrix,
                                  info1.matrix(i, true).matrix,
                                  row,
                                  col,
                                  info1.indices,
                                  info2.indices,
                                  level1,
                                  level2);
                  assemble_fluxes(o[level1].matrix,
                                  info2.matrix(i, false).matrix,
                                  row,
                                  col,
                                  info2.indices,
                                  info2.indices,
                                  level2,
                                  level2);
                  assemble_fluxes(o[level1].matrix,
                                  info2.matrix(i, true).matrix,
                                  row,
                                  col,
                                  info2.indices,
                                  info1.indices,
                                  level2,
                                  level1);
                }
            }
          else
            {
              Assert(level1 > level2, ExcNotImplemented());
              if (flux_up->size() != 0)
                {
                  // Do not add M22,
                  // which is done by
                  // the coarser cell
                  assemble_fluxes(o[level1].matrix,
                                  info1.matrix(i, false).matrix,
                                  row,
                                  col,
                                  info1.indices,
                                  info1.indices,
                                  level1,
                                  level1);
                  assemble_up(flux_up->block(i)[level1].matrix,
                              info1.matrix(i, true).matrix,
                              row,
                              col,
                              info1.indices,
                              info2.indices,
                              level1,
                              level2);
                  assemble_down(flux_down->block(i)[level1].matrix,
                                info2.matrix(i, true).matrix,
                                row,
                                col,
                                info2.indices,
                                info1.indices,
                                level2,
                                level1);
                }
            }
        }
    }
  } // namespace Assembler
} // namespace MeshWorker

DEAL_II_NAMESPACE_CLOSE

#endif


