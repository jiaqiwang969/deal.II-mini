//include/deal.II-translator/A-headers/main_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2020 by the deal.II authors
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
 * @mainpage
 * 这是deal.II类和函数文档的主要起始页。关于其他方面的文档，如构建系统，可以在其他地方找到。此外，还有<a
 * href="Tutorial.html">Tutorial programs on the use of the library</a>。
 * deal.II库中的许多类可以被分组为模块（见<a
 * href="modules.html">Modules
 * page</a>或本页面顶部菜单中的相应条目）。这些模块围绕着任何有限元程序的构建块而形成。下面的点击图给出了deal.II中主要类组的交互方式的概要，下面有更详细的描述（灰色方框表示可选的外部库的子集，灰色椭圆表示可选的外部应用程序的子集，deal.II可以与之交互）。
 *
* @dot
 digraph G
{
  graph[rankdir="TB",bgcolor="transparent"];
  node [fontname="FreeSans",fontsize=15,
        shape=record,height=0.2,width=0.4,
        color="black", fillcolor="white", style="filled"];
  edge [color="black", weight=10];
  tria       [label="Triangulation",    URL="\ref grid"];
  fe         [label="Finite elements",    URL="\ref feall"];
  mapping    [label="Mapping",          URL="\ref mapping"];
  quadrature [label="Quadrature",       URL="\ref Quadrature"];
  dh         [label="DoFHandler",       URL="\ref dofs"];
  fevalues   [label="FEValues",         URL="\ref feaccess"];
  systems    [label="Linear systems",   URL="\ref LAC"];
  solvers    [label="Linear solvers",   URL="\ref Solvers"];
  output     [label="Graphical output", URL="\ref output"];
  manifold   [label="Manifold",         URL="\ref manifold"];
  tria -> dh              [color="black",style="solid"];
  fe -> dh                [color="black",style="solid"];
  fe -> fevalues          [color="black",style="solid"];
  mapping -> fevalues     [color="black",style="solid"];
  quadrature -> fevalues  [color="black",style="solid"];
  dh -> systems           [color="black",style="solid"];
  fevalues -> systems     [color="black",style="solid"];
  systems -> solvers      [color="black",style="solid"];
  solvers -> output       [color="black",style="solid"];
  manifold -> tria        [color="black",style="solid"];
  manifold -> mapping     [color="black",style="solid"];
  node [fontname="FreeSans",fontsize=12,
        shape=record,height=0.2,width=0.4,
        color="gray55", fontcolor="gray55", fillcolor="white", style="filled"];
  edge [color="gray55", weight=1];
  opencascade [label="OpenCASCADE"];
  subgraph linalglibs {
    rank="same";
    petsc       [label="PETSc",    URL="\ref PETScWrappers"];
    trilinos    [label="Trilinos", URL="\ref TrilinosWrappers"];
    cuda        [label="CUDA",     URL="\ref CUDAWrappers"];
  }
  umfpack     [label="UMFPACK"];
  petsc -> systems        [dir="none"];
  petsc -> solvers        [dir="none"];
  trilinos -> systems     [dir="none"];
  trilinos -> solvers     [dir="none"];
  cuda -> systems         [dir="none"];
  cuda -> solvers         [dir="none"];
  umfpack -> solvers      [dir="none"];
  opencascade -> manifold [dir="none"];
  node [fontname="FreeSans",fontsize=12,
        shape=ellipse,height=0.2,width=0.4,
        color="gray55", fontcolor="gray55", fillcolor="white", style="filled"];
  edge [color="gray55", weight=1];
  gmsh        [label="gmsh", URL="\ref Gmsh"];
  visit       [label="VisIt"]
  paraview    [label="ParaView"]
  gmsh -> tria       [dir="none"];
  output -> visit    [dir="none"];
  output -> paraview [dir="none"];
}
 * @enddot

 * 这些组在教程程序中都有涉及，在 step-3
 * 中首先概述了它们的组合方式。下面是这个组的分类指南，以及与每个组有关的文档链接。
 * <ol>   <li>  <b>%Triangulation</b>。三角形是单元及其低维边界对象的集合。单元是参考超立方体[0,1]<sup>dim</sup>在 @ref mapping 模块中的适当映射下的图像。
 * 三角化存储了网格的几何和拓扑属性：单元如何连接以及它们的顶点在哪里。三角化不知道任何关于你可能想在这个网格上使用的有限元的信息，三角化甚至不知道任何关于其单元形状的信息：在2D中它只知道一个单元有4个面（线）和4个顶点（在3D中它有6个面（四边形）、12条线和8个顶点），但其他的一切都由映射类来定义。
 * 三角形的属性和数据几乎都是通过所有单元的循环来查询的，可能还会查询每个单元的所有面。因此，关于网格的大部分知识都隐藏在
 * @em
 * 迭代器的后面，即类似指针的结构，可以从一个单元迭代到下一个单元，并且可以询问它目前指向的单元的信息。
 * 描述三角形和单元格的类位于 @ref
 * grid 模块中，并有相关文档。迭代器在 @ref Iterators
 * 模块中描述。
 * <li>
 * <b>%Manifold</b>:矩阵描述了单元格的形状，更广泛地说，描述了要解决方程的领域的几何形状。它们使用微分几何的语言。更多信息可以在
 * @ref manifold  中找到。
 * <li>  <b>Finite Element</b>
 * 。有限元类描述了定义在单元格上的有限元空间的属性。这包括，例如，有多少自由度位于顶点，在线上，或在单元的内部。除此之外，有限元类当然还必须提供单元格上各点的单个形状函数的值和梯度。
 * 有限元类在 @ref
 * feall 模块中描述。
 * <li>
 * <b>%Quadrature</b>:与有限元一样，正交对象是在单元格上定义的。它们只描述单元格上正交点的位置，以及其上正交点的权重。
 * 描述特定正交公式的类的文档可以在 @ref
 * Quadrature 模块中找到。
 * <li>  <b>%DoFHandler</b>:
 * %DoFHandler对象是三角形和有限元的汇合点：有限元类描述了每个顶点、线条或单元需要多少自由度，DoFHandler类分配了这个空间，使三角形的每个顶点、线条或单元有正确的数量。它也给它们一个全局的编号。
 * 一个不同的观点是这样的。网格和有限元描述了我们寻求离散解的有限维空间
 * $V_h$
 * 的抽象属性，而%DoFHandler类列举了这个空间的具体基础，因此我们可以通过一个有序的系数集
 * $U_j$ 来表示离散解。
 * 就像三角化对象一样，对DoFHandlers的大多数操作都是通过在所有单元上循环，并对每个单元或其中的一个子集进行操作。因此，这两个类的接口相当相似：它们允许获得第一个和最后一个单元（或面，或线等）的迭代器，并通过这些迭代器提供信息。从这些迭代器中可以得到的信息是已经可以从三角形迭代器中得到的几何和拓扑信息（它们实际上是派生类），以及诸如当前单元上的自由度的全局数字。我们也可以要求一个迭代器从一个数据向量中提取与当前单元格上的自由度相对应的值，该向量存储了与三角形相关的所有自由度的值。
 * 值得注意的是，就像三角剖分一样，DoFHandler类不知道任何关于从单元格到其单个单元格的映射。它也不知道对应于它所管理的自由度的形状函数：它所知道的是，例如，每个顶点有2个自由度，每个单元格内部有4个。除了它们存在的事实之外，它们的具体细节与DoFHandler类无关。
 * DoFHandler类及其关联物在 @ref
 * dofs模块中描述。此外，还有一些专门的版本可以处理多级和hp-discretizations。这些都在
 * @ref mg 和 @ref hp
 * 模块中描述。有限元方法经常意味着对自由度的约束，例如对悬挂节点或适用边界条件的节点的约束；处理这种约束在
 * @ref constraints 模块中描述。
 * <li>
 * <b>%Mapping</b>。有限元程序的下一步是，人们希望利用有限元的形状函数和正交规则定义的正交点，计算三角形的每个单元上的矩阵和右手条目或其他数量。为此，有必要将形状函数、正交点和正交权重从单元格映射到三角形的每个单元。这不是由映射和派生类直接完成的，而是由映射和派生类来促进的：它们描述了如何将点从单位空间映射到实空间并返回，以及提供这个导数的梯度和雅各布决定因素。
 * 这些类都在 @ref
 * mapping 模块中描述。
 * <li>
 * <b>%FEValues</b>。下一步是实际取一个有限元，在映射到实数单元时，在正交公式定义的点上评估其形状函数及其梯度。这就是FEValues类和兄弟姐妹的领域：从某种意义上说，它们提供了一个有限元函数空间的点式视图。
 * 这似乎有局限性：在数学分析中，我们总是以单元或单元面的积分来写公式，涉及到有限元形状函数。因此，人们会认为有必要将有限元空间描述为连续空间。然而，在实践中，这是没有必要的：在实际计算中，所有的积分都被使用正交公式的近似值所取代，因此真正需要的是在域内有限数量的给定位置评估形状函数的能力。FEValues类正是提供了这种信息。给定有限元、正交和映射对象，它们计算连续函数空间（相对于离散的，而不是相对于不连续的）对离散的点的限制。
 * 有许多对象可以做到这一点。FEValues用于对单元格进行评估，FEFaceValues用于对单元格的面进行评估，FESubfaceValues用于对单元格的部分面进行评估。所有这些类都在 @ref
 * feaccess 模块中描述。
 * <li>  <b>Linear
 * Systems</b>。如果知道如何使用FEValues和朋友们评估单个单元上的形状函数的值和梯度，并且知道如何使用DoFHandler迭代器获得单元上自由度的全局数，那么下一步就是使用问题的双线性形式来组合线性系统的系统矩阵（和右手边）。然后我们将从这个线性系统中确定我们问题的解决方案。
 * 要做到这一点，我们需要有存储和管理矩阵和向量条目的类。deal.II为此提供了一整套的类，以及与其他提供类似功能的软件包的接口。这方面的文档可以在 @ref
 * LAC 模块中找到。
 * <li>
 * <b>Linear
 * Solvers</b>:为了确定一个有限维度的线性方程组的解，人们需要线性求解器。在有限元应用中，它们经常是迭代式的，但有时也可能要使用直接或稀疏的直接求解器。deal.II有相当多的此类求解器。它们被记录在
 * @ref Solvers 模块中。
 * <li>
 * <b>Output</b>。最后，一旦在给定的三角形上获得了有限元问题的解，人们往往希望用可视化程序对其进行后处理。这个库本身并不做可视化处理，而是生成各种图形格式的输出文件，这些格式可以被广泛使用的可视化工具所理解。
 * 在 @ref 输出模块中给出了关于这样做的类的描述。   </ol>
 * 此外，deal.II还有一些超越这里所列举的类组。它们涉及到上面介绍的层次结构中更细化的概念，或者涉及到诸如处理输入和输出这样的切身问题，这些问题不一定专门针对有限元程序，但也出现在那里。这些类都列在类和命名空间的视图中，可以从本页顶部的菜单栏中找到，并且也被分组为自己的模块（见本页顶部的<a
 * href="modules.html">Modules link</a>）。
 * 我们提供了Doxygen标签文件，供那些想把应用程序的文档直接链接到deal.II在线文档的用户使用。该标签文件在<a
 * href="../deal.tag"><code>deal.tag</code></a>。对于deal.II的每个版本，它都驻留在Doxygen参考文档的正上方的目录中。为了使用这个标签文件，你必须把它下载到一个Doxygen可以找到的地方。之后，在你的Doxygen选项文件中找到
 * <code>TAGFILES</code> 这个键，然后写一些类似于<pre> TAGFILES =
 * deal.tag=http://www.dealii.org/X.Y.Z/doxygen/deal.II </pre> 其中
 * <code>X.Y.Z</code>
 * 指的是你要链接的版本。请确保你使用匹配的标签文件。理论上，你也可以针对deal.II的发展中的修订版进行链接，但你必须担心，如果deal.II的结构发生变化，你的链接可能会变得无效。
 *
 */


