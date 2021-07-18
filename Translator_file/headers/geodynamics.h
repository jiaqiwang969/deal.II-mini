//include/deal.II-translator/A-headers/geodynamics_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2020 by the deal.II authors
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
 * @defgroup geodynamics The geodynamics demonstration suite
 * deal.II的 @ref Tutorial "教程 "
 * 包含一组程序，它们共同构成了地球动力学示范套件。这些程序的想法是利用地球动力学的应用来演示高级有限元软件的技术，也就是对固体地球过程的研究。通过这样做，这些程序应该为解决实际地球动力学问题的更专业的专用程序提供一个基础，例如作为研究生或博士后工作的一部分。下面将对这些计划的动机进行更深入的讨论。
 * 目前，地球动力学测试套件包含以下程序。
 * -  step-8  : 弹性
 *
 * -  step-18  : 一个%平行弹性求解器
 *
 * -  step-20  : 多孔介质流
 *
 * -  step-21  : 通过多孔介质的多相流动
 *
 * -  step-22  :斯托克斯流
 *
 * -  step-31  : 热对流（Boussinesq流动
 * -  step-32  ：用于地幔对流的%平行布西尼克解算器
 * 其中一些程序是根据加州理工学院的合同开发的，得到了美国国家科学基金会EAR-0426271号奖的支持，这是资助<a
 * target="_top" href="http://www.geodynamics.org">Computational
 * Infrastructure in
 * Geodynamics</a>计划的第一笔拨款。接受者Wolfgang
 * Bangerth对这一支持来源表示衷心感谢。
 *
 *  <h3>Rationale</h3>
 * 自适应网格细化（AMR）长期以来一直被认为是一项关键技术，它可以帮助精确和有效地解决一些地球动力学应用的数值问题。它在地球动力学界已经讨论了好几年，并且自CIG成立以来一直是其任务清单上的一个议题。然而，到目前为止，在这个方向上发生的事情还比较少。直到最近，才有了在地球动力学中使用AMR的尝试。CIG于2007年10月在Boulder举办了一次关于AMR技术的研讨会；George
 * Biros、Omar Ghattas、Mike
 * Gurnis和ShijieZhong小组之间的合作目前正在开发一个%并行的自适应地幔对流处理器；deal.II的一些主要开发者最终开发了用于模拟地幔对流的<a
 * href="https://aspect.geodynamics.org">ASPECT
 * code</a>，现在已经是相当成熟和广泛使用的代码。
 * AMR技术在地球动力学中应用缓慢的原因之一是最初的障碍比较大：代码必须提供数据结构和算法来处理自适应网格，有限元必须能够处理悬挂节点，等等。要做到这一点，在足够的通用性下，有限元程序要增加几万行代码，对于普通学生来说，在学位论文的时间范围内是无法做到的。另一方面，有一些库提供了基础代码，支持AMR的应用可以在此基础上迅速建立。当然，deal.II正是提供了这种基础。
 * 地球动力学测试套件的目标是为与地球动力学相关的各种主题编写程序。继续保持现有教程程序的风格
 *
 * - 一个广泛的介绍，解释一个应用的背景和形式，以及在其解决方案中使用的数值方案的概念；整个代码中的详细评论，解释实施细节；以及一个显示数值结果的部分
 *
 * - 我们打算将所产生的程序作为解决模型问题的有据可查的应用程序来提供。特别是，它们旨在实现以下目标。 <ul>   <li>  <i>Starting points:</i> 现有的deal.II教程已被证明是一个很好的起点，可供研究生和研究人员快速开发自己的应用程序。通过提供已经接近目标应用的程序，通常可以很快获得第一个结果，既保持了开发过程中最初的热情，也允许将研究时间用于实现特定的应用行为，而不是将几个月的工作用于支持AMR的基本基础代码。
 * 支持这一观点的事实是，尽管有<a
 * href="https://www.dealii.org/publications.html">more than 1,000
 * publications</a>介绍了用deal.II获得的结果，但我们知道只有相对较少的应用是用deal.II从头开始建立的；所有其他的应用都是从某个教程程序的修改开始的。
 * <li>  <i>Training:</i>
 * 我们建议编写的教程程序将为学生和研究人员提供当前数值技术的参考实现，如AMR、高阶元素、复杂的线性和非线性求解器、稳定技术等。提供这些作为其他人进一步开发的起点，也将有助于实现在现代数值算法方面培训新一代地球动力学家的目标。
 * <li>  <i>Extending equations and formulations:</i>
 * 在deal.II中，用另一个方程扩展一组方程是相当简单的，例如，一个额外的平流量作为右手边或在一个系数中进入现有方程。由于应用通常使用封锁的矩阵，而不是用一个大矩阵代替所有东西的方法，所以为增强方程找到合适的线性求解器也不复杂。因此，deal.II是一个很好的工具，可以尝试更复杂的问题公式，或更完整的模型及其对解的准确性的影响。
 * <li>  <i>Rapid prototyping and benchmarking:</i>
 * deal.II提供了许多可互换的组件，允许快速建立有限元种类和顺序、稳定技术或线性求解器的原型。例如，通常只需要改变几行代码就可以用高阶元素取代低阶元素。通过这种方式，尝试高阶元素、不同的块消除求解器或不同的稳定技术变得相对简单。反过来，这可能有助于在计算求解时间和数值解的准确性方面对应用进行基准测试。
 * 本模块中的应用将已经过正确性的基准测试。现有的教程程序通常采用更简单而不是更复杂的求解器方案进行阐述，但经常建议采用更复杂的方案，包括在附录中提示如何实现这些方案。
 * <li>  <i>Try algorithms:</i> deal.II的快速原型能力也可能有助于在deal.II适用的程序规模上确定最佳算法，然后在可以在更大规模的机器上运行的专用程序中实现这种特定的算法（没有能力轻易改变它）。例如，一个建立在deal.II上的小型地幔对流代码可以用来确定二阶元素是否对这个目的有用（例如，见 step-31 中所示的结果）。如果是这样，那么人们就可以在更大的代码中使用这种知识，比如上面提到的ASPECT代码。 </ul>
 *
 *
 */


