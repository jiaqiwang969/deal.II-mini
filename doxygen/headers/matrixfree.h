//include/deal.II-translator/A-headers/matrixfree_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2020 by the deal.II authors
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
 *    @defgroup matrixfree Matrix-free infrastructure
 * 本模块描述了deal.II中的无矩阵基础设施。以下是deal.II中主要的类群与无矩阵基础设施互动的概要，可点击图表，下面有更详细的描述。
 *
 *   @dot digraph G { graph[rankdir="TB",bgcolor="transparent"]; node
 * [shape=box,fontname="FreeSans",fontsize=15, height=0.2,width=0.4,
 * color="black", fillcolor="white", style="filled"]; edge [color="black",
 * weight=10]; subgraph base { rank="same"; tria       [label="Triangulation",
 * URL="\ref grid"]; fe         [label="FiniteElement",    URL="\ref feall"];
 * mapping    [label="Mapping",          URL="\ref mapping"]; quadrature
 * [label="Quadrature",       URL="\ref Quadrature"]; } dh
 * [label="DoFHandler",       URL="\ref dofs"]; simd     [label="SIMD",
 * fontname="FreeSans",fontsize=12, height=0.2,width=0.4, color="gray",
 * fontcolor="gray", fillcolor="white", style="filled"]; fevalues
 * [label="FEEvaluation", fillcolor="deepskyblue"]; mf
 * [label="MatrixFree loops", fillcolor="deepskyblue"]; cuda
 * [label="CUDA",     URL="\ref CUDAWrappers",
 * fontname="FreeSans",fontsize=12, height=0.2,width=0.4, color="gray",
 * fontcolor="gray", fillcolor="white", style="filled"]; tbb
 * [label="TBB", fontname="FreeSans",fontsize=12, height=0.2,width=0.4,
 * color="gray", fontcolor="gray", fillcolor="white", style="filled"];
 * {rank=same simd
 *
 * -> fevalues        [dir="none", color="transparent"]; fevalues
 *
 * -> mf          [dir="none", color="transparent"]; mf
 *
 * -> cuda              [dir="none", color="transparent"]; cuda
 *
 * -> tbb             [dir="none", color="transparent"]; } subgraph sol {
 * rank="same"; solvers [label="Solvers",   URL="\ref Solvers",
 * fillcolor="deepskyblue"]; gmg     [label="Geometric Multigrid",
 * fontname="FreeSans",fontsize=12, height=0.2,width=0.4, color="black",
 * fontcolor="black", fillcolor="white", style="dashed"]; } output
 * [label="Graphical output", URL="\ref output"]; manifold
 * [label="Manifold",         URL="\ref manifold"]; tria
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
 * -> mf                [color="black",style="solid"]; mf
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
 * -> mapping     [color="black",style="solid"]; node
 * [fontname="FreeSans",fontsize=12, shape=record,height=0.2,width=0.4,
 * color="gray", fontcolor="gray", fillcolor="white", style="filled"]; edge
 * [color="gray", weight=1]; opencascade [label="OpenCASCADE"]; subgraph
 * misclibs { systems    [label="Operators", fillcolor="deepskyblue"]; }
 * opencascade
 *
 * -> manifold [dir="none"];
 *
 * node [fontname="FreeSans",fontsize=12, shape=ellipse,height=0.2,width=0.4,
 * color="gray", fontcolor="gray", fillcolor="white", style="filled"]; edge
 * [color="gray", weight=1]; gmsh        [label="gmsh", URL="\ref Gmsh"];
 * visit       [label="VisIt"] paraview    [label="ParaView"] gmsh
 *
 * -> tria       [dir="none"]; output
 *
 * -> visit    [dir="none"]; output
 *
 * -> paraview [dir="none"]; } @enddot
 * 从本质上讲，FEEvaluation类在MatrixFree的数据存储之上提供的框架是一个专门的运算符评估框架。它目前只与库中提供的具有特殊结构的元素子集兼容，即那些基础可以被描述为一维多项式的张量乘积的元素。这为向量条目和正交点的值或梯度之间的有效转换提供了机会，这种技术被称为和因子化。这种技术起源于谱元界，由Orszag在1980年的工作开始。虽然这种技术最初不过是一种装配向量（或矩阵）的特殊技术，比通用车辆FEValues更快，但它的效率使得在迭代求解器中使用这些积分设施直接评估矩阵-向量产品成为可能，而不是先装配一个矩阵，然后用这个矩阵做矩阵-向量产品。这一步最初是非直观的，与许多人在数学和计算机科学教育中所接受的教育相违背，包括大多数deal.II的开发者，因为一次又一次地重新计算积分，而不是使用预先计算的数据，似乎是一种浪费。然而，正如教程程序
 * step-37 、 step-48 、 step-59 、 step-64 和 step-67
 * 所示，这些概念在现代计算机架构上通常比传统算法更出色。
 * 有利于无矩阵计算的两个主要原因如下。   <ol>   <li>  无矩阵方法跳过了大的全局稀疏矩阵的存储，并在飞行中计算基础的弱形式。由于内存传输，即从RAM内存中读取数据的速度，是基于矩阵的计算的瓶颈，而不是使用这些数据所做的实际算术，因此，一个读取数据较少的无矩阵评估，即使做了较多的计算，也会有优势。这个概念是建立在计算机架构的一个趋势之上的，这个趋势最好的描述就是<i>memory wall</i>，说的是计算性能比内存性能增长得更快。因此，一定程度的算术运算基本上是免费的，而且这个份额在过去20年里变得更大。除了在显式时间积分中的经典使用外，它还使这种激进的算法转换成为迭代求解器的矩阵-向量积的无矩阵实现。当然，实现必须是高效的，而且不能有多余的计算量，以使其在总体上获胜。deal.II库使用SIMD矢量化和基于多项式程度模板的高度优化内核来实现这一目标。为了提供一个视角，二次元FE_Q的稀疏矩阵-向量乘积曾经与2005-2007年左右设计的处理器（如奔腾4或AMD Opteron Barcelona，每个芯片有2-4个内核）上的无矩阵实现速度相当。到2018年，无矩阵评估的速度大约是8倍（在英特尔Skylake服务器上测量，14个核心）。   <li>  无矩阵方法随着度数的增加，每个自由度的复杂度更好，这是由于和因子化的原因。对于无矩阵方案来说，每个自由度的工作随着 $\mathcal O(k)$ 度的增加而增加，而对于基于矩阵的方法来说，它随着 $\mathcal O(k^d)$ 度的增加而增加。这使高阶方案具有优势。在无矩阵评估中一个特别好的特点是 $\mathcal O(1)$ 项经常占主导地位，所以看起来高阶方法在评估时间上和低阶方法一样快，当他们有相同的自由度数量时。对于deal.II中的实现，最佳吞吐量通常是在3到6的多项式度数之间实现的。   </ol>
 * 总而言之，无矩阵计算是高阶元素（高阶意味着除了线性形状函数之外的一切）的方法，并用于显式时间步进（
 * step-48
 * ）或迭代求解器，其中也可以用无矩阵的方式进行预处理，正如
 * step-37 和 step-59 教程程序中所展示的。 <h3>The matrix-free
 * evaluation infrastructure</h3>
 * 顶层接口由FEEvaluation类提供，其中还包含了对不同用例的广泛描述。
 * <h4>The FEEvaluation class hierarchy</h4>
 * FEEvaluation类派生于FEEvaluationAccess类，后者又继承于FEEvaluationBase。FEEvaluation类本身不仅对维度、分量的数量和数字类型（如双数或浮点数）进行模板化，而且对多项式程度和每个空间方向上的正交点数量进行模板化。这些信息用于将和分解中的循环长度传递给各自的内核（见`tensor_product_kernels.h`和`evaluation_kernels.h`）并确保最佳效率。所有访问向量或提供访问单个正交点的数据字段的方法都继承自FEEvaluationAccess。
 * FEEvaluationAccess类的动机是允许根据分量的多少，对内插解字段的值和梯度访问进行专业化。而基类FEEvaluationBase将梯度作为一个
 * "张量<1,n_components,张量<1,dim,矢量数组<Number>>"返回，外张量经过组件，内张量经过梯度的`dim'组件。对于一个标量场，即`n_components=1'，我们可以跳过外张量，简单地使用`Tensor<1,dim,VectorizedArray<Number>'作为梯度类型。同样地，对于一个`n_components=dim'的系统，梯度的适当格式是`Tensor<2,dim,VectorizedArray<Number>'。
 * <h4>The FEFaceEvaluation class</h4>
 * 面积分，如连续有限元中的不均匀诺伊曼条件或一大类非连续Galerkin方案，除了单元积分外，还需要评估面的正交点上的量。面评价的设施大多与FEEvaluation共享，即FEFaceEvaluation也继承自FEEvaluationAccess。所有关于自由度和形状函数的数据字段都可以重复使用，后者是因为所有信息都由一维形状数据组成。然而，关于映射数据，由于数据是
 * "structdim=dim-1
 * "的，所以使用了一个特殊化。因此，FEEvaluationAccess和FEEvaluationBase被赋予一个模板参数`is_face`，以分别持有指向单元和面的映射信息的指针。除了用
 * FEEvaluationAccess::get_value() 访问函数值或用
 * FEEvaluationAccess::get_gradient(),
 * 访问梯度外，面评估器还可以用
 * FEEvaluationAccess::get_normal_vector()
 * 访问法向量和一个专门的字段
 * FEEvaluationAccess::get_normal_derivative(),
 * ，返回解场对面的法向导数。这个量被计算为梯度（在实空间）乘以法向量。梯度和法向量的组合是许多（简单的）二阶椭圆方程的典型特征，例如用内部惩罚法对拉普拉斯进行离散化。如果不需要单独的梯度，联合操作大大减少了数据访问，因为每个正交点只需要`dim`数据条目`normal
 * Jacobian`，而单独访问时需要`dim^2`字段的Jacobian和`dim`字段的normal。
 * 计算面积分的一个重要优化是考虑必须访问的矢量数据量，以评估面的积分。例如，想想FE_DGQ的情况，即拉格朗日多项式，其部分节点在元素边界上。对于函数值的评估，只有
 * $(k+1)^{d-1}$ 自由度通过非零基函数做出贡献，而其余的
 * $(k+1)^d$
 * 基函数在该边界上评估为零。由于矢量访问是无矩阵计算的瓶颈之一，对矢量的访问应该限制在有趣的条目上。为了实现这种设置，方法
 * FEFaceEvaluation::gather_evaluate() （和
 * FEFaceEvaluation::integrate_scatter()
 * 的积分等价物）将矢量访问与插值到正交点相结合。存在两种特殊情况，包括前面提到的
 * "非零 "值的情况，它被存储为字段
 * internal::MatrixFreeFunctions::ShapeInfo::nodal_at_cell_boundaries.
 * 。对于在一个面上只有选定数量的基函数的值和第一导数评估为非零的情况，也可以有类似的属性。相关的元素类型是FE_DGQHermite，决定存储在属性
 * internal::MatrixFreeFunctions::tensor_symmetric_hermite.
 * 中，是否可以使用这样的优化内核的决定是在
 * FEFaceEvaluation::gather_evaluate() 和
 * FEFaceEvaluation::integrate_scatter().
 * 中自动做出的。为每个积分任务做这个决定似乎效率不高，但最终这是一个单一的`if`语句（条件跳跃），对于现代CPU来说很容易预测，因为决定在一个积分循环中总是相同。(人们只需支付一定程度上增加的编译时间，因为编译器需要为所有路径生成代码，虽然)。
 * <h3>The data storage through the MatrixFree class</h3>
 * 由FEEvaluation和FEFaceEvaluation执行的任务可以分成三类。<i>index
 * access into vectors</i>，<i>evaluation and integration on the unit
 * cell</i>，和<i>operation on quadrature points including the geometry
 * evaluation</i>。这种分割反映在MatrixFree所包含的主要数据字段上，分别用
 * internal::MatrixFreeFunctions::DoFInfo,
 * internal::MatrixFreeFunctions::ShapeInfo, 和
 * internal::MatrixFreeFunctions::MappingInfo
 * 表示这三个类别。它们的设计原则和内部布局在以下几个小节中描述。
 * 所有这些数据结构坚持的主要界面是，集成任务被分解成一个单元或面的范围，人们可以通过一个整数索引来索引。单元积分、面内积分和边界积分的整数范围的信息是由类
 * internal::MatrixFreeFunctions::TaskInfo,
 * 使用数据字段`cell_partition_data`、`face_partition_data`和`boundary_partition_data`提供。这个类还包含了用于使用线程并行调度任务的索引子范围的信息，以及在`{cell,face,boundary}_partition_data`内对索引范围的分组，用于交错单元和面的积分，这样对单元和面的积分的向量项的访问就会重新使用已经在缓存中的数据。
 * <h4>Index storage: the internal::MatrixFreeFunctions::DoFInfo struct</h4>
 * DoFInfo类的主要目的是提供矢量访问函数 FEEvaluationBase::read_dof_values() 和 FEEvaluationBase::distribute_local_to_global(). 所消耗的索引，这些索引布置如下。   <ol>   <li>  指数存储在MPI本地索引空间中，以实现直接的数组访问，而不是将全局索引转换为本地索引。后者绝对不利于性能。 </li>   <li>  指数被存储在一个叫做 internal::MatrixFreeFunctions::DoFInfo::dof_indices, 的字段中，这是一个长索引数组。以<i>cell index</i>为单位的访问粒度由辅助字段 internal::MatrixFreeFunctions::DoFInfo::row_starts 控制，它类似于压缩矩阵存储中的行开始索引。该方案支持可变长度，因为我们支持hp-adaptivity和由于主索引阵列中包含的约束而产生的索引间接性。由于单元格上的矢量化，访问粒度最初会以<i>cell batches</i>为单位。然而，我们必须能够同时访问单个单元，例如，对于面的积分，面的批次一般与单元的批次不同，因此访问不是线性的。此外，如果我们为每个单独的组件提供一个<i>start index</i>，那么对多组件系统的支持就变得透明了。因此，`row_starts`字段的长度为 `n_cell_batches()*VectorizedArray<Number>::%size()*n_components`.   </li>   <li>  在一个多基元的系统中，组件之间的转换由四个变量控制  <ol>   <li>`std::vector<unsigned  int> n_components`（每个基元的组件），  </li>   <li>`std::vector<unsigned  int> start_components`（从基元到唯一组件编号的转换），  </li>   <li>`std::vector<unsigned  int> component_to_base_index`（从唯一元件编号到基数索引的翻译），以及  </li>   <li>`std::vector<std::vector<unsigned  int>> component_dof_indices_offset`（特定元件的自由度范围在一个单元上的全部自由度列表中的偏移）。 </li>   </ol>   </li> 。
 * <li>  在hp-adaptive计算中提取FE指数的信息。 </li>   <li>  关于 "第一次访问 "特定向量条目的信息，该信息用于在第一次访问目标向量之前不久将其清零的 MatrixFree::loop 。这被用来避免向整个向量写零，从而破坏了数据位置性。 </li>   </ol> .
 * internal::MatrixFreeFunctions::DoFInfo 中数据结构的设置是在
 * internal::MatrixFreeFunctions::DoFInfo::read_dof_indices,
 * 中完成的，在这里我们首先假设了一个非常一般的有限元布局，无论是连续的还是不连续的元素，并且我们解决了由于悬挂节点而产生的约束。这个初始步骤是在单元的原始排序中完成的。在后面的阶段，这些单元一般会被重新排列，以反映我们在最终循环中通过单元的顺序，我们也会在DoF指数中寻找可以利用的模式，如单元内连续的指数范围。这种重新排序是为了实现与MPI的通信和计算的重叠（如果启用的话），并在单元格上形成更好的具有矢量化的批次组。指数的数据存储在这个最终的顺序中是线性的，安排在
 * internal::MatrixFreeFunctions::DoFInfo::reorder_cells. 。
 * 因为存储索引的数据量是不可忽略的，所以对于携带更多结构的特殊配置，减少数据量是值得的。一个例子是FE_DGQ的情况，每个单元的一个索引就足以描述其所有的自由度，其他的则是连续的顺序。类
 * internal::MatrixFreeFunctions::DoFInfo 包含一个特殊的向量数组
 * internal::MatrixFreeFunctions::DoFInfo::dof_indices_contiguous
 * ，每个单元包含一个数字。由于单元和面的积分使用不同的访问模式，而且这种特殊情况下的数据很小，我们最好存储3个这样的向量，一个用于装饰为`内部'的面（索引0），一个用于装饰为`外部'的面（索引1），一个用于单元（索引2），而不是通过
 * internal::MatrixFreeFunctions::FaceInfo. 使用定向。
 * DoFInfo中有一系列额外的特殊存储格式可用。关于在deal.II中实现的选项及其动机，我们参考该结构的文档
 * internal::MatrixFreeFunctions::DoFInfo::IndexStorageVariants 。
 * 最后，DoFInfo类还持有一个共享指针，描述向量的并行分区。由于
 * Utilities::MPI::Partitioner, 的限制，传递给 MatrixFree::reinit()
 * 函数的单个DoFHandler对象内的索引必须在每个MPI进程中是连续的，也就是说，本地范围必须最多包括一个块。除了基本的分区器，该类还提供了一组更严格的索引集，只涉及所有鬼魂索引的一个子集，被添加到向量的鬼魂范围。这些交换模式被设计为与通过
 * internal::MatrixFreeFunctions::ShapeInfo::nodal_at_cell_boundaries
 * 的减少索引访问相结合，例如。
 * MatrixFree类支持多个DoFHandler对象，以传递给
 * MatrixFree::reinit()
 * 函数。对于这些DoFHandler对象中的每一个，都会创建一个单独的
 * internal::MatrixFreeFunctions::DoFInfo
 * 对象。在MatrixFree中，我们存储了一个 `std::vector` 的
 * internal::MatrixFreeFunctions::DoFInfo 对象来说明这一事实。
 * <h4>The internal::MatrixFreeFunctions::ShapeInfo structure</h4>
 * 一维形状函数在一维正交点上的评估被存储在
 * internal::MatrixFreeFunctions::ShapeInfo. 类中
 * 更确切地说，我们保存了所有的函数值、梯度和豫备值。此外，面的形状函数的值和导数，即单位区间的0和1点，也被存储。对于悬空节点上的面积分，相邻两个单元中较粗的单元必须对数值进行插值，而不是插值到完整的正交点，而只是插值到一个子面（评估点要么按比例调整为[0，1/2]，要么为[1/2，1]）。这种情况由数据字段`values_within_subface`、`gradients_within_subface`和`hessians_within_subface`处理。这个数据结构也会检查形状函数相对于参考单元中心的对称性（在这种情况下，会应用所谓的偶数变换，进一步减少计算量）。
 * <h4>The internal::MatrixFreeFunctions::MappingInfo structure</h4>
 * 评估的几何信息存储在 internal::MatrixFreeFunctions::MappingInfo.
 * 类中 与 internal::MatrixFreeFunctions::DoFInfo
 * 类类似，在一个MatrixFree实例中可以有多种变体，在这种情况下是基于多个正交公式。此外，单元格和面的单独数据被存储。由于涉及到更多的逻辑，而且字段之间有协同作用，字段的
 * `std::vector` 被保存在 internal::MatrixFreeFunctions::MappingInfo.
 * 中。单个字段是 internal::MatrixFreeFunctions::MappingInfoStorage
 * 类型的，保存有反雅各布系数、JxW值、法向量、法向量乘反雅各布系数的数组（对于
 * FEEvaluationAccess::get_normal_derivative()),
 * 实空间的正交点，以及参考元素上的正交点。我们使用一个辅助索引数组，指向每个单元的数据起点，即雅各布、JxW值和法向量的`data_index_offsets`字段，以及正交点的`quadrature_point_offsets`。这种偏移使HP-adaptivity的字段长度可变，类似于DoFInfo的做法，但它也使我们称之为<i>geometry
 * compression</i>的东西得以实现。为了减少数据访问，我们检测单元格的简单几何形状，其中雅各布系数在一个单元格内是恒定的，或者也是跨单元格的，使用
 * internal::MatrixFreeFunctions::MappingInfo::cell_type:  。
 * <ol>   <li>  笛卡尔单元是指雅各布系数为对角线且在单元的每个正交点上相同的单元。每个单元只需要存储一个字段。由于单元内的相似性，我们还检查了当前处理器上所有单元的雅各布系数相同的其他单元批。这可以进一步减少内存访问。由于一般情况下的JxW值存储的是雅各布系数乘以正交权重，但我们只想为笛卡尔单元保留一个字段，所以我们在笛卡尔情况下误用了<i>JxW</i>这个名字，只存储了雅各布系数的行列式，没有正交权重。因此，在 FEEvaluationBase::submit_value() 和类似的情况下，我们需要注意，因为我们仍然必须乘以权重。   <li>  平行单元在整个单元内有恒定的雅各布系数，所以每个单元只需要存储一个字段。由于单元内的相似性，我们还检查了当前处理器上所有单元的雅各布系数相同的其他单元批。由于一般情况下的JxW值存储的是雅各布系数乘以正交权重，但我们只想为一个仿生单元保留一个字段，所以我们在仿生情况下滥用了<i>JxW</i>这个名字，就像在笛卡尔情况下一样，只存储雅各布系数的行列式，而不存储正交权重。因此，在 FEEvaluationBase::submit_value() 和类似的情况下，我们需要注意，因为我们仍然必须乘以权重。   <li> 在面孔上，我们可以有这样的特殊情况：当JxW值不同时，法向量在所有正交点都是相同的。这种情况适用于平坦的面孔。为了减少数据访问，我们在 internal::MatrixFreeFunctions::GeometryType. 中保留了这个作为压缩索引的第三个选项。与笛卡尔和仿射的情况相反，在数组中只保留一个字段，扁平面为所有正交点保留一个单独的条目（保留一个索引字段`data_index_offsets`），但只访问第一个字。   <li>  一般类型的索引是指没有发现压缩的单元或面。在这种情况下，我们也不寻找在一个以上的单元上找到相同图案的机会，尽管这种情况可能存在，比如对于挤压的网格。这种搜索操作是基于使用自定义浮点比较器`FPArrayComparator`将数据插入到 `std::map` 中，当每个单元使用单一数据字段时，效率足够高。然而，如果对所有单元的所有正交点（有许多不同的情况）进行，它将是相当昂贵的。   </ol>
 * internal::MatrixFreeFunctions::MappingInfo
 * 的实现被分成了单元格和面的部分，所以这两个部分可以很容易地被分开持有。让代码读起来有点尴尬的是，我们需要从FEValues对象中完成的原始标量评估中把几个对象批在一起，我们需要识别重复的数据字段，而且我们需要用
 * `std::map` 对笛卡尔和仿射的情况在几个单元中定义压缩。
 * internal::MatrixFreeFunctions::MappingInfo
 * 的数据计算部分除了明显的MPI并行化外，还通过任务来并行化。每个处理器在一个子范围内计算信息，然后数据最终被复制到一个单一的组合数据域中。
 * <h3>Identification and parallelization of face integrals</h3>
 * 目前MatrixFree中的面积分方案为所有的面建立了一个独立的任务列表，而不是明确地通过一个单元的`2*dim`面。这样做的好处是一个面的所有信息只需要处理一次。典型的DG方法计算的数字通量是保守的，也就是说，从面的两边看都是一样的，无论什么信息离开一个单元都必须正好再次进入邻居的单元。有了这个方案，它们必须只被计算一次。同时，这也确保了几何信息也必须只被加载一次。一个可能的缺点是，基于面的独立编号方法使得基于线程的并行变得比基于单元的方法复杂得多，因为在这种方法中，只有当前单元的信息被写入，而邻居的信息只被读取）。
 * 由于面是独立于单元的，它们得到了自己的矢量布局。从一批面孔来看，不管是什么连续的一批细胞都会被交织在一起，这是面孔的本质（在这里我们只把细胞内有相同面孔索引的面孔放在一起，以此类推）。脸部循环的设置，在文件`脸部_设置_内部.h`中完成，试图提供至少部分类似于单元格补丁的脸部批次，以增加数据的定位性。沿着这些思路，在典型的
 * MatrixFree::loop
 * 情况下，面的工作也与单元的工作交错进行，也就是说，函数调用中返回的`cell_range'和`face_range'参数通常很短。
 * 因为所有来自两边的积分都是一次性执行的，所以出现了一个问题，即在子域边界的两个处理器中哪一个被分配到一个面。本模块的作者进行了大量的实验，发现应用于自由度存储的方案，即把所有可能重叠的项目分配给一个处理器，是相当不平衡的，面的数量最多有20%的差异。为了提高性能，在`face_setup_internal.h`中实现了一个平衡的方案，将每对处理器之间的所有接口分成两块，一块由一个处理器完成，另一块由另一个处理器完成。尽管这增加了通过MPI发送的消息的数量，但这是值得的，因为负载变得更加平衡。另外，当本地问题的规模为100,000
 * DoFs（三维）时，消息相当大，约为5-50kB。在这种消息大小下，延迟通常小于吞吐量。
 * 脸部数据默认不被初始化，但必须由
 * MatrixFree::AdditionalData,
 * 中的脸部更新标志触发，即`mapping_update_flags_inner_faces`或`mapping_update_flags_boundary_faces`设置为与`update_default`不同的值。
 * <h3>Invoking MatrixFree::loop</h3>
 * MatrixFree类支持两种类型的实体的循环。第一种，从2012年开始在deal.II主分支上可用，是只执行单元积分，使用三个`cell_loop`函数中的一个，该函数需要一个指向单元操作的函数指针。第二种设置，在2018年引入，是一个循环，也可以执行面和/或边界积分，简单地称为`loop`。这需要三个函数指针，分别解决单元工作、内面工作和边界面工作。
 * 除了以适当的方式安排工作外，该循环还执行两个任务。   <ol>   <li>  对`src`和`dst`向量进行数据交换，分别调用`update_ghost_values()`和 `compress(VectorOperation::add)`, 。如果各自的标志 MatrixFree::AdditionalData::overlap_communication_computation 被设置为 "true"（默认），那么交换可以以异步的方式进行，与不需要远程处理器数据的单元的工作重叠。   <li>  使用相应的标志将`dst`向量归零。在循环内这样做的好处是，循环知道向量中的哪些条目（首先）被单元格和面的循环中的一些子ranges触及。因此，它可以将向量逐个归零，以确保我们不需要两次访问向量条目（一次归零，一次添加贡献）。这似乎是一个微小的优化，但事实上运算符的评估可以非常快，简单地将一个向量归零就可以花费运算符评估时间的20%左右，所以这确实是值得努力的由于这个参数有一些实验性，DoFInfo类保留了一个静态变量 internal::MatrixFreeFunctions::DoFInfo::chunk_size_zero_vector ，在这里可以调整（如果有人认为其他东西会更好，例如因为未来的计算机看起来和2018年推出这个参数时不同）。   </ol>
 * 最后， MatrixFree::loop
 * 函数还需要一个参数来传递面积分的数据访问类型，由结构
 * MatrixFree::DataAccessOnFaces,
 * 描述，以减少处理器之间需要交换的数据量。不幸的是，目前还没有办法将这些信息传达给
 * MatrixFree::loop
 * ，这些信息在FEFaceEvaluation内部通过评估类型（值和/或梯度）和底层形状函数的组合获得，以避免在第二个地方手动设置这类信息。
 *
 */


