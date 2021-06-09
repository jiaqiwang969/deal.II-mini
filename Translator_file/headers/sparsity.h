//include/deal.II-translator/A-headers/sparsity_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2015 by the deal.II authors
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
 *    @defgroup Sparsity Sparsity patterns
 * 几乎所有的有限元公式都会导致矩阵的
 * "稀疏"，即每行的非零元素数量(i)与矩阵的整体大小相比相对较小，并且(ii)被一个固定的数字所限制，如果网格被细化，这个数字不会增长。在这种情况下，不存储矩阵的<i>all</i>元素，而只存储那些实际（或可能）为非零的元素会更有效。这需要为每一行存储非零项的列索引（我们称之为
 * "稀疏模式"），以及这些非零项的实际值。在实践中，有时会出现一些非零值实际上为零的情况。稀疏模式和稀疏矩阵只打算为<i>may</i>非零的条目提供空间，而且是在我们还不知道这些条目最终会有什么值的时候这样做；如果一个系数或单元格恰好有特定的值，它们可能会有一个零值）。)
 * 在deal.II中，稀疏模式通常与实际的稀疏矩阵分开（除了SparseMatrixEZ类和一些来自外部库接口的类，如PETSc）。原因是人们经常有几个共享相同稀疏模式的矩阵；例子包括时间步进方案所需的刚度和质量矩阵，或者广义特征值问题的左手和右手矩阵。因此，如果它们中的每一个都必须单独存储它们的稀疏性模式，那将是一种浪费。
 * 因此，deal.II具有矩阵类所建立的稀疏模式类。有两组主要的稀疏模式类，如下所述。
 *
 *  <h4>"Static" sparsity patterns</h4>
 * deal.II中的主要稀疏矩阵类，SparseMatrix，只为每个矩阵条目存储一个值，但不存储这些条目的位置。为此，它依赖于从与该矩阵相关的稀疏模式对象中获得的信息。这个稀疏性模式对象必须是SparsityPattern类型的。
 * 因为矩阵是大的对象，而且改变它们的成本相对较高，所以SparsityPattern对象的构建分为两个阶段：首先，在一个
 * "动态
 * "阶段，分配期望在其上构建的矩阵有非零条目的位置；在第二个
 * "静态 "阶段，这些非零位置的表示被 "压缩
 * "成通常的压缩稀疏行（CSR）格式。在这之后，不能再添加新的非零位置。只有在压缩之后才能将稀疏模式与矩阵联系起来，因为后者需要前者的高效压缩数据格式。在动态阶段建立一个稀疏模式经常发生在
 * DoFTools::make_sparsity_pattern()
 * 函数中。虽然这看起来是一个限制，但首先建立一个稀疏模式，然后只在先前分配的位置写入矩阵，这通常不是一个重要的问题，因为在有限元代码中，通常很清楚矩阵中哪些元素可能是非零的，哪些肯定是零的。
 * 这种两阶段生成稀疏模式的优点是，当它实际用于矩阵时，有一个非常有效的格式。特别是，条目的位置被存储在一个线性数组中，允许快速访问，对具有深层次缓存的现代CPU类型很友好。因此，静态SparsityPattern类是deal.II的主SparseMatrix类可以工作的唯一对象。
 * 静态稀疏模式的主要缺点是，它们的有效构造需要合理地猜测每一行最大可能有多少条目。在实际构建过程中，例如在
 * DoFTools::make_sparsity_pattern()
 * 函数中，最多只能分配到之前所说的那么多条目。这是一个问题，因为通常很难估计每行的最大条目数。因此，一个常见的策略是首先建立和中间的稀疏模式，在构建稀疏模式的过程中使用效率较低的存储方案，然后直接复制到静态的、压缩的形式。大多数教程程序都是这样做的，从
 * step-2 开始（也可参见，例如 step-11 ， step-18 ，和 step-27
 * 教程程序）。
 *
 *  <h4>"Dynamic" or "compressed" sparsity patterns</h4>
 * 如上所述，要获得一个稀疏模式的每一行的最大条目数的良好估计往往很复杂。因此，任何试图用不好的估计来分配一个普通的SparsityPattern的做法都需要大量的内存，几乎所有的内存都不会被使用，在压缩时被取消分配。
 * 为了避免这种情况，deal.II包含一个名为DynamicSparsityPattern的
 * "动态 "或 "压缩
 * "稀疏性模式，它只分配必要的内存来容纳当前添加的条目。虽然这比上面提到的最坏情况下的行为节省了很多内存，但它需要使用效率较低的存储方案来插入元素，而且频繁分配内存往往也需要大量计算时间。然而，避免过多的内存分配的权衡是无法避免的。
 * 该类通常以如下方式使用
 * @verbatim
 * DynamicSparsityPattern dsp (dof_handler.n_dofs());
 * DoFTools::make_sparsity_pattern (dof_handler,
 *                                dsp);
 * constraints.condense (dsp);
 *
 * SparsityPattern final_sparsity_pattern;
 * final_sparsity_pattern.copy_from (dsp);
 * @endverbatim
 *
 * 中间的、压缩的稀疏模式被直接复制到最终静态模式的
 * "压缩 "形式中。 <h4>Dynamic block sparsity patterns</h4>
 * BlockDynamicSparsityPattern类实现了一个用于构造块矩阵的动态稀疏模式数组。更多信息见文档和
 * step-22 。
 *
 * @ingroup LAC
 *
 */


