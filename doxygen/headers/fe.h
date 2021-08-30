//include/deal.II-translator/A-headers/fe_0.txt
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
 *    @defgroup feall Finite elements
 * 所有与形状函数和对形状函数的访问有关的类。
 * 这涉及到有限元的实际值。关于自由度的编号请参考  @ref
 * dofs  的模块。
 * 本模块的类和函数分为几个子组，在上面列出的各自子模块中讨论。此外，FETools类提供了提供有限元素、元素间变换等信息的函数。
 * 从总体上看，这个模块的各个部分与库中的各种其他部分相互作用。
 * @dot digraph G { graph[rankdir="TB",bgcolor="transparent"];
 *
 * node [fontname="FreeSans",fontsize=15, shape=box,height=0.2,width=0.4,
 * color="black", fillcolor="white", style="filled"]; edge [color="black",
 * weight=10];
 *
 * tria       [label="Triangulation",    URL="\ref grid"]; fe
 * [label="Finite elements",    URL="\ref feall", fillcolor="deepskyblue"];
 * mapping    [label="Mapping",          URL="\ref mapping"]; quadrature
 * [label="Quadrature",       URL="\ref Quadrature"]; dh
 * [label="DoFHandler",       URL="\ref dofs"]; fevalues   [label="FEValues",
 * URL="\ref feaccess"]; systems    [label="Linear systems",   URL="\ref
 * LAC"]; solvers    [label="Linear solvers",   URL="\ref Solvers"]; output
 * [label="Graphical output", URL="\ref output"]; manifold
 * [label="Manifold",         URL="\ref manifold"];
 *
 * tria
 *
 * -> dh              [color="black",style="solid"]; fe
 *
 * -> dh                [color="black",style="solid"]; fe
 *
 * -> fevalues          [color="black",style="solid"]; mapping
 *
 * -> fevalues     [color="black",style="solid"]; quadrature
 *
 * -> fevalues  [color="black",style="solid"]; dh
 *
 * -> systems           [color="black",style="solid"]; fevalues
 *
 * -> systems     [color="black",style="solid"]; systems
 *
 * -> solvers      [color="black",style="solid"]; solvers
 *
 * -> output       [color="black",style="solid"]; manifold
 *
 * -> tria        [color="black",style="solid"]; manifold
 *
 * -> mapping     [color="black",style="solid"];
 *
 * { rank=same mapping
 *
 * -> quadrature [dir="none", color="transparent"]; quadrature
 *
 * -> fe      [dir="none", color="transparent"]; fe
 *
 * -> tria            [dir="none", color="transparent"]; }
 *
 * } @enddot
 *
 */


/**
 *    @defgroup febase Base classes
 * 这个子模块的成员描述了有限元类的实现机制，而没有实际实现一个具体的元素。例如，FiniteElement基类声明了派生类如果要描述一个有限元空间所必须实现的虚拟函数。同样地，FiniteElementData持有描述有限元特征的某些数值的变量，例如每个顶点、线或面的自由度数。
 * 另一方面，像FE_Poly和FE_PolyTensor这样的类是高级抽象。它们描述了建立在单元格上的形状函数的多项式描述之上的有限元。从它们派生出来的类只需要提供一个特定的多项式的描述，而有限元就是从这个描述中建立的。例如，实现通常的拉格朗日元素的FE_Q类使用FE_Poly基类来生成有限元，为它提供一组拉格朗日插值多项式，对应于插值点的等距细分。
 * 最后，FESystem类用于处理矢量值问题。在这里，人们可能想把一些标量（或者也是矢量值）基元耦合在一起，形成矢量值算子的联合有限元。例如，对于三维Navier-Stokes流动，人们可能希望用三个Q1元素来表示速度的三个分量，用一个片状常数Q0元素来表示压力。FESystem类可以用来将这四个基本元素组合成一个具有4个矢量分量的单一矢量值元素。
 * step-8 、 step-17 和 step-18
 * 教程程序介绍了该类在矢量值弹性（Lam&eacute;）方程的使用。
 * step-20
 * 讨论了一个混合拉普拉斯离散化，也使用了矢量值元素。
 *
 * @ingroup feall
 *
 */


/**
 *    @defgroup feaccess Finite element access/FEValues classes
 * 本模块中的类在人们想要组装矩阵或矢量时使用。它们将有限元、正交对象和映射联系起来：有限元类描述单元格上的有限元空间（即单位线段、正方形或立方体<tt>[0,1]^d</tt>），正交类描述正交点的位置和它们的权重，映射类描述如何将一个点从单元格映射到实数单元并返回。由于积分发生在实单元上的正交点，需要知道它们的位置以及这些点的有限元形状函数的值和梯度。FEValues类可以协调获得这些信息。对于面的积分（例如边界上的积分，或者单元间的界面），FEFaceValues类提供了与FEValues类对单元类似的功能。最后，FESubfaceValues类提供了在面的一部分进行积分的可能性，如果相邻的单元被细化，并且目前的单元只与相邻的单元共享其面的一部分。如果使用矢量值元素，FEValues和相关的类允许访问所有的矢量组件；如果想挑选单个组件，有一些提取器类可以使这个任务更简单，如 @ref
 * vector_valued 模块中所述。
 * 这一组的最后一个成员，UpdateFlags枚举，是作为一种优化使用的：与其让FEValues类计算与一个单元上给定的有限元有关的每一个可能的数据，你必须预先指定你真正感兴趣的信息。UpdateFlags枚举用于提供符号名称，表示您希望FEValues类计算的内容。
 * 从 step-3 开始，所有这些类都在 @ref Tutorial 的 "教程程序 "
 * 中使用，并且在那里有详细的描述。
 * FEValues类和朋友们的实际工作很复杂，因为它必须是通用的，但又是有效的。 @ref
 * UpdateFlags 的页面试图对其工作原理做一个概述。
 * 从总体上看，这个模块的各个部分与库的其他部分相互作用。
 * @dot digraph G { graph[rankdir="TB",bgcolor="transparent"];
 *
 * node [fontname="FreeSans",fontsize=15, shape=box,height=0.2,width=0.4,
 * color="black", fillcolor="white", style="filled"]; edge [color="black",
 * weight=10];
 *
 * tria       [label="Triangulation",    URL="\ref grid"]; fe
 * [label="Finite elements",    URL="\ref feall"]; mapping
 * [label="Mapping",          URL="\ref mapping"]; quadrature
 * [label="Quadrature",       URL="\ref Quadrature"]; dh
 * [label="DoFHandler",       URL="\ref dofs"]; fevalues   [label="FEValues",
 * URL="\ref feaccess", fillcolor="deepskyblue"]; systems    [label="Linear
 * systems",   URL="\ref LAC"]; solvers    [label="Linear solvers",
 * URL="\ref Solvers"]; output     [label="Graphical output", URL="\ref
 * output"]; manifold   [label="Manifold",         URL="\ref manifold"];
 *
 * tria
 *
 * -> dh              [color="black",style="solid"]; fe
 *
 * -> dh                [color="black",style="solid"]; fe
 *
 * -> fevalues          [color="black",style="solid"]; mapping
 *
 * -> fevalues     [color="black",style="solid"]; quadrature
 *
 * -> fevalues  [color="black",style="solid"]; dh
 *
 * -> systems           [color="black",style="solid"]; fevalues
 *
 * -> systems     [color="black",style="solid"]; systems
 *
 * -> solvers      [color="black",style="solid"]; solvers
 *
 * -> output       [color="black",style="solid"]; manifold
 *
 * -> tria        [color="black",style="solid"]; manifold
 *
 * -> mapping     [color="black",style="solid"];
 *
 * { rank=same mapping
 *
 * -> quadrature [dir="none", color="transparent"]; quadrature
 *
 * -> fe      [dir="none", color="transparent"]; fe
 *
 * -> tria            [dir="none", color="transparent"]; }
 *
 * } @enddot
 * @ingroup feall
 *
 */


/**
 *    @defgroup fe Finite element space descriptions
 * 这里的类描述了有限元空间，如最简单的Q1（双/三线）空间，以及高阶拉格朗日空间Qp，但也有更专业的空间，如Nedelec或Raviart-Thomas空间。具体的实现是由抽象的FiniteElement基类派生的。
 * 本质上，这些类必须实现的函数提供了查询单元格上某一点的形状函数的值或导数的能力。为了在整合矩阵和右手条目中发挥作用，我们必须有能力将这些形状函数和梯度映射到实际单元中。这是由Mapping基类（见 @ref
 * mapping ）和FEValues类（见 @ref feaccess
 * ）共同派生的类来完成的。 <h3>Vector-valued finite elements</h3>
 * deal.II提供了两种不同类型的向量值元素。首先，有一组真正的矢量元素，通常通过以下事实来区分，即每个矢量分量由一组不同的各向异性多项式组成。这些元素通常与微分形式相关。目前，它们是
 * <ul>   <li>  FE_ABF  <li>  FE_BDM, FE_DGBDM  <li>  FE_Nedelec, FE_DGNedelec  <li>  FE_RaviartThomas, FE_DGRaviartThomas  </ul>  。
 * 另外，deal.II提供了一种机制，可以从现有的标量或矢量元素中创建一个矢量元素。FESystem类负责这个工作：它本身并不描述形状函数，而是从其他有限元对象中组装一个矢量值的有限元。这个功能在
 * step-8 ,  step-17 和之后的其他教程程序中有所描述。
 *
 *
 * @note
 * FE_PolyTensor类提供了对矢量值元素的实现支持。通常情况下，一个新的向量元素应该派生自这个类。
 * <h3>Discontinuous Galerkin</h3>
 * 对于每个符合任何弱可微函数空间的有限元，如<i>H<sup>1</sup></i>或<i>H<sup>curl</sup></i>，我们可以通过简单地将顶点、边或面的所有自由度分配到单元的内部来定义一个类似的DG空间。这要从拓扑学的角度来理解。这种自由度的插值算子仍然会在边界上。  虽然没有这样做，但我们提供了很多这样的元素，加上那些没有符合要求的对应元素，比如FE_DGP。以下是当前DG元素的列表。   <ul>   <li>  标量。FE_DGP, FE_DGQ  <li>  标量，不同的形状函数。FE_DGPMonomial, FE_DGPNonparametric, FE_DGQArbitraryNodes  <li>  矢量值的。  FE_DGBDM, FE_DGNedelec, FE_DGRaviartThomas  </ul>
 *
 *
 * @note
 * FE_DGVector类支持向量值DG元素的实现，其方式是只需要提供向量多项式空间。由此派生的实际类只需要实现一个构造函数和
 * FiniteElement::get_name().  。
 *
 * @ingroup feall
 *
 */


/**
 *    @defgroup mapping Mappings between reference and real cell
 * 本模块中的类用于从单位坐标映射到真实单元格的坐标。最常见的是，人们使用MappingQ1类，它提供了一个Q1（双线/三线）映射（即对通常的Q1元素来说是一个等价的映射）。然而，还有其他一些实现高阶映射的类，以提供曲线元素。这些在
 * step-11 和 step-12 的教程程序中讨论。
 * MappingQ1Eulerian类是对MappingQ1类的扩展，它接受一个描述域的每个位置的位移场的向量。这在欧拉计算中使用，不需要在每个时间步长后实际移动顶点。
 * 此外，MappingC1类提供了一个计算域的边界，这个边界不仅是弯曲的，而且在边界上两个单元之间的界面上有一个连续的导数。
 * 最后，MappingCartesian类是对砖形且边缘与坐标轴平行的元素的优化。
 * 从总体上看，这个模块的各个部分与库的其他各种部分相互作用。
 * @dot digraph G { graph[rankdir="TB",bgcolor="transparent"];
 *
 * node [fontname="FreeSans",fontsize=15, shape=box,height=0.2,width=0.4,
 * color="black", fillcolor="white", style="filled"]; edge [color="black",
 * weight=10];
 *
 * tria       [label="Triangulation",    URL="\ref grid"]; fe
 * [label="Finite elements",    URL="\ref feall"]; mapping
 * [label="Mapping",          URL="\ref mapping", fillcolor="deepskyblue"];
 * quadrature [label="Quadrature",       URL="\ref Quadrature"]; dh
 * [label="DoFHandler",       URL="\ref dofs"]; fevalues   [label="FEValues",
 * URL="\ref feaccess"]; systems    [label="Linear systems",   URL="\ref
 * LAC"]; solvers    [label="Linear solvers",   URL="\ref Solvers"]; output
 * [label="Graphical output", URL="\ref output"]; manifold
 * [label="Manifold",         URL="\ref manifold"];
 *
 * tria
 *
 * -> dh              [color="black",style="solid"]; fe
 *
 * -> dh                [color="black",style="solid"]; fe
 *
 * -> fevalues          [color="black",style="solid"]; mapping
 *
 * -> fevalues     [color="black",style="solid"]; quadrature
 *
 * -> fevalues  [color="black",style="solid"]; dh
 *
 * -> systems           [color="black",style="solid"]; fevalues
 *
 * -> systems     [color="black",style="solid"]; systems
 *
 * -> solvers      [color="black",style="solid"]; solvers
 *
 * -> output       [color="black",style="solid"]; manifold
 *
 * -> tria        [color="black",style="solid"]; manifold
 *
 * -> mapping     [color="black",style="solid"];
 *
 * { rank=same mapping
 *
 * -> quadrature [dir="none", color="transparent"]; quadrature
 *
 * -> fe      [dir="none", color="transparent"]; fe
 *
 * -> tria            [dir="none", color="transparent"]; }
 *
 * } @enddot
 * @ingroup feall
 *
 */


