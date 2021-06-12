//include/deal.II-translator/A-headers/mg_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2020 by the deal.II authors
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
 * @defgroup mg Multilevel support 
 * 与多栅格算法有关的类。
 * 实现多网格方案的主要类是Multigrid，其功能是
 * Multigrid::cycle().  它使用以下抽象类来执行多网格循环。
 * <ol>   <li>  MGMatrixBase包含水平矩阵，在 mg::Matrix   <li>  MGCoarseGridBase是最粗层次上的求解器。   <li>  MGSmootherBase在每个层次上进行平滑处理。   <li>  MGTransferBase组织层次间的转移。   </ol>
 * 此外，还有一个PreconditionMG类，它是Multigrid的一个封装器，具有deal.II
 * @ref
 * Preconditioners的标准接口。PreconditionMG也使用继承自MGTransferBase的类，例如MGTransferPrebuilt，它使用
 * MGTransferPrebuilt::copy_to_mg() 和
 * MGTransferPrebuilt::copy_from_mg_add(),
 * ，在全局向量和水平向量之间进行转移。
 * 最后，我们有几个辅助类，即MGLevelObject，它在每个层次上存储一个对象*。
 * 关于如何使用这一功能，请参见 step-16  、 step-16  b和
 * step-39  示例程序。 <h3>Multigrid and hanging nodes</h3>
 * 在自适应细化网格上使用多网格方法比使用常规细化方法涉及更多的基础设施。首先，为了保持最佳的复杂度，我们需要决定如何在每个层次上进行平滑处理。为此，我们必须在多级分解的意义上定义什么是级。
 * 首先，我们定义多网格意义上的一个层次是由网格层次结构中某一层次的所有单元构成的。因此，某一级的平滑只限于由这一级或更细的单元组成的子域。这通常被称为局部平滑。这种定义的优点是，多网格方案的层次矩阵可以通过遍历某一层次的所有单元而轻松组装起来，而且这些层次矩阵不包含悬空节点。
 * 这种分解的缺点是，我们需要额外的矩阵来处理细化边上出现的问题。此外，根据方法是连续的（因此在细化边缘有自由度）还是不连续的（在细化边缘采用通量矩阵），处理方法是不同的。虽然这些矩阵很小，但我们必须把它们集合起来，并通知多棱镜方法。
 *
 */

/**
 * 这个命名空间包含了在我们知道在局部细化和块系统的背景下需要什么之后，对多层次支持的重新实现。
 *
 * @ingroup mg
 *
 */
namespace mg
{
}


