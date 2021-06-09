//include/deal.II-translator/A-headers/global_dof_index_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2013 - 2018 by the deal.II authors
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


/**
 * @page GlobalDoFIndex When to use types::global_dof_index instead of unsigned int
 * deal.II可以被配置为使用64位的自由度指数，而不是通常的无符号整数，在目前大多数系统上默认为32位。这是必要的，因为我们希望能够解决超过40亿个未知数的问题（32位无符号整数所能表示的极限）。同时，我们不想不分青红皂白地将deal.II中的所有整数替换成64位版本，因为这将增加许多地方的内存使用，我们在这些地方表示的数量肯定不会超过40亿。
 * 我们为这些指数定义的数据类型是 <code>\#ifdef</code>
 * ，以保持代码库的大部分不受 types::global_dof_index.
 * 的影响。如果deal.II被正常配置，这种类型是 <code>unsigned
 * int</code> ，但如果提供正确的标志，可以切换到
 * <code>unsigned long long int</code>
 * （见ReadMe文件）。本页旨在澄清何时必须使用
 * types::global_dof_index ，何时可以使用普通无符号整数。 <dl>
 * <dt class="glossary">  @anchor  GlobalDoFIndexBlockIndices
 * <b>BlockIndices</b></dt>  <dd>
 * 块的数量是一个无符号的int，因为这个数字预计会很低，即小于四亿。然而，块的大小是一个
 * types::global_dof_index ，因为每个块可以是任意的大。   </dd>
 * <dt class=" glossary">  @anchor  GlobalDoFIndexCell <b>Cell</b></dt>  <dd>
 * 单元的ID是不唯一的。不同细化程度的单元和/或不同处理器上的单元可以有相同的ID。因此，所有与单元相关的数据都可以是无符号的int，因为在一个处理器上，一个网格级别，肯定不会有超过40亿的单元。
 * </dd> <dt class=" glossary">  @anchor  GlobalDoFIndexDoFHandler
 * <b>DoFHandler</b></dt>  <dd>
 * 每个自由度的ID在并行计算中是唯一的。因此，自由度是
 * types::global_dof_index.   </dd>  。 <dt class=" glossary">  @anchor
 * GlobalDoFIndexFullMatrix <b>FullMatrix</b></dt>  <dd>
 * 行和列的数量是 types::global_dof_index
 * ，即使不期望有人会创建一个有这么多条目的FullMatrix。然而，AffineConstraints类的一些方法是以矩阵类型为模板的，因此，FullMatrix的大小必须与SparseMatrix的大小为同一类型。
 * </dd> <dt class=" glossary">  @anchor  GlobalDoFIndexSparseMatrix
 * <b>SparseMatrix</b></dt>  <dd>
 * SparseMatrix的大小可以是任意大的，可以想象，在单个节点上有足够的内存，可以生成一个超过40亿行或列的矩阵。因此，采用了
 * types::global_dof_index
 * 。然而，即使对于我们现在可以解决的大型复杂问题，期望稀疏矩阵中的非零条目数超过40亿是不合理的。因此，我们仍然使用无符号int，例如，
 * SparsityPattern::row_lengths 和类似的函数。   </dd> </dl>
 *
 */


