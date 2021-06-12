// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2021 by the deal.II authors
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
 * @page Tutorial Tutorial programs
 *
 * New to deal.II? You might want to start with tutorial Step-1 and work
 * your way up to Step-5. At that point you can explore what features you
 * are interested in and look at the large collection of programs listed
 * below.
 *
 * The deal.II tutorial contains a collection of programs, each more or
 * less built atop of previous ones, which demonstrate various aspects of
 * the library. Each such example has the following structure:
 * <ol>
 *   <li> <b>Introduction:</b> What the program does, including
 *        the mathematical model, and
 *        what programming techniques are new.
 *   <li> <b>The commented program:</b> An extensively documented listing of the
 *        source code.
 *   <li> <b>Results:</b> The output of the program, with comments and
 *        interpretation.
 *   <li> <b>The plain program:</b> The source code stripped of
 *        all comments.
 * </ol>
 * You can browse the available tutorial programs
 * <ol>
 *   <li> as <b><a href="#graph">a graph</a></b> that shows how the major
 *      concepts of each tutorial programs builds on previous ones (though each
 *      program may also use minor pieces from other programs not specifically
*       connected in the graph).
 *   <li> as <b><a href="#list">a list</a></b> that provides a short
 *     synopsis of each program.
 *   <li> or <b><a href="#topic">grouped by topic</a></b>.
 * </ol>
 *
 * The programs are in the <code>examples/</code> directory of your local
 * deal.II installation. After compiling the library itself, if you go into
 * one of the tutorial directories, you can configure the program by typing
 * <code>cmake .</code>, build it via <code>make</code> and run it using
 * <code>make run</code>. The latter command also compiles the program if
 * that has not already been done. The CMakeLists.txt files in the
 * different directories are based on the
 * <a href="../../users/cmake_user.html#cmakeauto" target="_top">autopilot
 * style CMakeLists.txt example</a>.
 *
 * @note Some of the tutorial programs also jointly form
 *   the <a href="../../doxygen/deal.II/group__geodynamics.html">geodynamics
 *   demonstration suite</a>. More, often more complex but less well documented,
 *   deal.II-based programs than the ones that form the tutorial can also be
 *   found in the @ref CodeGallery .
 *
 *
 * <a name="graph"></a>
 * @anchor TutorialConnectionGraph
 * <h3>Connections between tutorial programs</h3>
 *
 * The following graph shows the connections between tutorial programs and
 * how their major components build on each other.
 * Click on any of the boxes to go to one of the programs. If you hover
 * your mouse pointer over a box, a brief description of the program
 * should appear.
 * @dot
digraph StepsMap
{
  overlap=false;
  edge [fontname="FreeSans",
        fontsize="10",
        labelfontname="FreeSans",
        labelfontsize="10",
        color="black",
        style="solid"];
  node [fontname="FreeSans",
        fontsize="10",
        shape="rectangle",
        height=0.2,
        width=0.4,
        color="black",
        fillcolor="white",
        style="filled"];
  Step1 [label="1", URL="\ref step_1", tooltip="Creating a mesh, refining it, writing it to a file.",height=.8,width=.8,shape="octagon",fillcolor="green"];
  Step10 [label="10", URL="\ref step_10", tooltip="Higher order mappings.",height=.35,width=.35,fillcolor="orange"];
  Step11 [label="11", URL="\ref step_11", tooltip="Higher order mappings. Dealing with constraints.",height=.35,width=.35,fillcolor="orange"];
  Step12 [label="12", URL="\ref step_12", tooltip="The Discontinuous Galerkin method for a linear advection problem.",height=.35,width=.35,fillcolor="orange"];
  Step12b [label="12b", URL="\ref step_12b", tooltip="Discontinuous Galerkin for linear advection, MeshWorker version.",height=.35,width=.35,fillcolor="orange"];
  Step13 [label="13", URL="\ref step_13", tooltip="Modularity. Software design.",height=.35,width=.35,fillcolor="orange"];
  Step14 [label="14", URL="\ref step_14", tooltip="Duality based error estimates. Adaptivity.",height=.35,width=.35,fillcolor="orange"];
  Step15 [label="15", URL="\ref step_15", tooltip="A nonlinear elliptic problem. Newton's method.",height=.35,width=.35,fillcolor="orange"];
  Step16 [label="16", URL="\ref step_16", tooltip="Multigrid on adaptive meshes.",height=.35,width=.35,fillcolor="orange"];
  Step16b [label="16b", URL="\ref step_16b", tooltip="MeshWorker for multigrid on adaptive meshes.",height=.35,width=.35,fillcolor="orange"];
  Step17 [label="17", URL="\ref step_17", tooltip="Parallel computing using MPI, using PETSc.",height=.25,width=.25,fillcolor="lightblue"];
  Step18 [label="18", URL="\ref step_18", tooltip="Quasistatic elasticity. More parallel computing.",height=.25,width=.25,fillcolor="lightblue"];
  Step19 [label="19", URL="\ref step_19", tooltip="Coupling particles to the solution of partial differential equations.",height=.35,width=.35,fillcolor="orange"];
  Step2 [label="2", URL="\ref step_2", tooltip="Assigning degrees of freedom to a grid.",height=.8,width=.8,shape="octagon",fillcolor="green"];
  Step20 [label="20", URL="\ref step_20", tooltip="Mixed finite elements for the mixed Laplacian. Block solvers.",height=.25,width=.25,fillcolor="yellow2"];
  Step21 [label="21", URL="\ref step_21", tooltip="Two-phase flow in porous media.",height=.25,width=.25,fillcolor="yellow2"];
  Step22 [label="22", URL="\ref step_22", tooltip="The Stokes equation on adaptive meshes.",height=.25,width=.25,fillcolor="yellow2"];
  Step23 [label="23", URL="\ref step_23", tooltip="Time dependent problems. The wave equation.",height=.25,width=.25,fillcolor="dodgerblue1"];
  Step24 [label="24", URL="\ref step_24", tooltip="The wave equation with absorbing boundary conditions. Extracting point values.",height=.25,width=.25,fillcolor="dodgerblue1"];
  Step25 [label="25", URL="\ref step_25", tooltip="The nonlinear sine-Gordon soliton equation.",height=.25,width=.25,fillcolor="dodgerblue1"];
  Step26 [label="26", URL="\ref step_26", tooltip="The heat equation. Time dependent meshes.",height=.25,width=.25,fillcolor="dodgerblue1"];
  Step27 [label="27", URL="\ref step_27", tooltip="Using the hp-finite element method for an elliptic problem.",height=.35,width=.35,fillcolor="orange"];
  Step28 [label="28", URL="\ref step_28", tooltip="Handling multiple meshes at the same time. Neutron transport.",height=.35,width=.35,fillcolor="orange"];
  Step29 [label="29", URL="\ref step_29", tooltip="A complex-valued Helmholtz equation. Sparse direct solvers.",height=.35,width=.35,fillcolor="orange"];
  Step3 [label="3", URL="\ref step_3", tooltip="Solving Poisson's equation.",height=.8,width=.8,shape="octagon",fillcolor="green"];
  Step30 [label="30", URL="\ref step_30", tooltip="Anisotropic refinement for DG methods.",height=.35,width=.35,fillcolor="orange"];
  Step31 [label="31", URL="\ref step_31", tooltip="Boussinesq flow for thermal convection.",height=.25,width=.25,fillcolor="yellow2"];
  Step32 [label="32", URL="\ref step_32", tooltip="A parallel Boussinesq flow solver for thermal convection in the earth mantle.",height=.25,width=.25,fillcolor="yellow2"];
  Step33 [label="33", URL="\ref step_33", tooltip="Hyperbolic conservation laws: the Euler equations of gas dynamics.",height=.25,width=.25,fillcolor="yellow2"];
  Step34 [label="34", URL="\ref step_34", tooltip="Boundary element methods for potential flow.",height=.25,width=.25,fillcolor="yellow2"];
  Step35 [label="35", URL="\ref step_35", tooltip="A projection solver for the Navier-Stokes equations.",height=.25,width=.25,fillcolor="yellow2"];
  Step36 [label="36", URL="\ref step_36", tooltip="Finding eigenvalues of the Schrödinger equation.",height=.35,width=.35,fillcolor="orange"];
  Step37 [label="37", URL="\ref step_37", tooltip="Matrix-free methods. Multigrid. Fast assembly techniques.",height=.35,width=.35,fillcolor="orange"];
  Step38 [label="38", URL="\ref step_38", tooltip="Solving the Laplace-Beltrami operator on a surface.",height=.35,width=.35,fillcolor="orange"];
  Step39 [label="39", URL="\ref step_39", tooltip="Interior Penalty for the Laplace equation. Adaptive refinement. Multigrid.",height=.35,width=.35,fillcolor="orange"];
  Step4 [label="4", URL="\ref step_4", tooltip="Dimension independent programming. Boundary conditions. Functions.",height=.8,width=.8,shape="octagon",fillcolor="green"];
  Step40 [label="40", URL="\ref step_40", tooltip="Solving the Laplace equation on adaptive meshes on thousands of processors.",height=.8,width=.8,shape="octagon",fillcolor="green"];
  Step41 [label="41", URL="\ref step_41", tooltip="Solving the obstacle problem: a variational inequality.",height=.25,width=.25,fillcolor="lightblue"];
  Step42 [label="42", URL="\ref step_42", tooltip="An adaptive, 3d solver for an elasto-plastic contact problem.",height=.25,width=.25,fillcolor="lightblue"];
  Step43 [label="43", URL="\ref step_43", tooltip="Efficient ways to solve two-phase flow problems on adaptive meshes in 2d and 3d.",height=.25,width=.25,fillcolor="yellow2"];
  Step44 [label="44", URL="\ref step_44", tooltip="Quasi-static finite-strain elasticity.",height=.25,width=.25,fillcolor="lightblue"];
  Step45 [label="45", URL="\ref step_45", tooltip="Periodic boundary conditions.",height=.35,width=.35,fillcolor="orange"];
  Step46 [label="46", URL="\ref step_46", tooltip="Coupling different physical models (flow, elasticity) in different parts of the domain.",height=.35,width=.35,fillcolor="orange"];
  Step47 [label="47", URL="\ref step_47", tooltip="Solving the fourth-order biharmonic equation.",height=.35,width=.35,fillcolor="orange"];
  Step48 [label="48", URL="\ref step_48", tooltip="Parallelization via MPI. The wave equation, in linear and nonlinear variants. Mass lumping. Fast assembly techniques.",height=.25,width=.25,fillcolor="dodgerblue1"];
  Step49 [label="49", URL="\ref step_49", tooltip="How to create and modify meshes.",height=.35,width=.35,fillcolor="orange"];
  Step5 [label="5", URL="\ref step_5", tooltip="Reading a grid from disk. Computations on successively refined grids. Variable coefficients. Assertions.",height=.8,width=.8,shape="octagon",fillcolor="green"];
  Step50 [label="50", URL="\ref step_50", tooltip="Multigrid on adaptive meshes distributed in parallel.",height=.35,width=.35,fillcolor="orange"];
  Step51 [label="51", URL="\ref step_51", tooltip="The convection-diffusion equation. Hybridizable discontinuous Galerkin methods. Face elements.",height=.35,width=.35,fillcolor="orange"];
  Step52 [label="52", URL="\ref step_52", tooltip="Time-dependent diffusion equation. Time stepping methods.",height=.25,width=.25,fillcolor="dodgerblue1"];
  Step53 [label="53", URL="\ref step_53", tooltip="Geometry: Dealing with deformed domains.",height=.35,width=.35,fillcolor="orange"];
  Step54 [label="54", URL="\ref step_54", tooltip="Geometry: Using industry standard IGES files as boundary descriptors.",height=.35,width=.35,fillcolor="orange"];
  Step55 [label="55", URL="\ref step_55", tooltip="Solving the Stokes problem in parallel.",height=.25,width=.25,fillcolor="yellow2"];
  Step56 [label="56", URL="\ref step_56", tooltip="Geometric multigrid preconditioners for the Stokes problem.",height=.35,width=.35,fillcolor="orange"];
  Step57 [label="57", URL="\ref step_57", tooltip="The Navier Stokes equations via Newton's iteration.",height=.25,width=.25,fillcolor="yellow2"];
  Step58 [label="58", URL="\ref step_58", tooltip="The complex-valued Nonlinear Schrödinger Equation.",height=.25,width=.25,fillcolor="dodgerblue1"];
  Step59 [label="59", URL="\ref step_59", tooltip="Matrix-free methods. Multigrid. Fast assembly techniques.",height=.35,width=.35,fillcolor="orange"];
  Step6 [label="6", URL="\ref step_6", tooltip="Adaptive local refinement. Higher order elements",height=.8,width=.8,shape="octagon",fillcolor="green"];
  Step60 [label="60", URL="\ref step_60", tooltip="The fictitious domain method using distributed Lagrange multipliers.",height=.35,width=.35,fillcolor="orange"];
  Step61 [label="61", URL="\ref step_61", tooltip="The Weak Galerkin method applied to the Poisson equation.",height=.35,width=.35,fillcolor="orange"];
  Step62 [label="62", URL="\ref step_62", tooltip="Elastic equation in the frequency domain. Calculating the transmission and resonance frequency of a phononic structure.",height=.35,width=.35,fillcolor="orange"];
  Step63 [label="63", URL="\ref step_63", tooltip="Block smoothers for Geometric Multigrid.",height=.35,width=.35,fillcolor="orange"];
  Step64 [label="64", URL="\ref step_64", tooltip="Matrix-free methods using CUDA and MPI.",height=.35,width=.35,fillcolor="orange"];
  Step65 [label="65", URL="\ref step_65", tooltip="Geometry: Working efficiently with expensive manifolds.",height=.35,width=.35,fillcolor="orange"];
  Step66 [label="66", URL="\ref step_66", tooltip="Parallel matrix-free geometric multigrid for a nonlinear problem.",height=.35,width=.35,fillcolor="orange"];
  Step67 [label="67", URL="\ref step_67", tooltip="An explicit time integrator for the Euler equations with matrix-free implementation.",height=.25,width=.25,fillcolor="yellow2"];
  Step68 [label="68", URL="\ref step_68", tooltip="A particle tracking problem using an analytically defined velocity field.",height=.35,width=.35,fillcolor="orange"];
  Step69 [label="69", URL="\ref step_69", tooltip="Hyperbolic conservation laws: a first-order guaranteed maximum wavespeed method for the compressible Euler equations.",height=.25,width=.25,fillcolor="yellow2"];
  Step7 [label="7", URL="\ref step_7", tooltip="Helmholtz equation. Computing errors. Boundary integrals.",height=.35,width=.35,fillcolor="orange"];
  Step70 [label="70", URL="\ref step_70", tooltip="A fluid structure interaction problem, using a penalty term.",height=.25,width=.25,fillcolor="yellow2"];
  Step71 [label="71", URL="\ref step_71", tooltip="Coupled constitutive modeling using automatic and symbolic differentiation",height=.35,width=.35,fillcolor="orange"];
  Step72 [label="72", URL="\ref step_72", tooltip="Newton's method for a nonlinear elliptic problem. Automatic differentiation",height=.35,width=.35,fillcolor="orange"];
  Step74 [label="74", URL="\ref step_74", tooltip="Symmetric interior penalty Galerkin for Poisson's equation. ",height=.35,width=.35,fillcolor="orange"];
  Step75 [label="75", URL="\ref step_75", tooltip="Solving the Laplace equation on hp-adaptive meshes with a MatrixFree hybrid",height=.35,width=.35,fillcolor="orange"];
  Step76 [label="76", URL="\ref step_76", tooltip="An alternative implementation for an explicit time integrator for the Euler equations with matrix-free implementation.",height=.25,width=.25,fillcolor="yellow2"];
  Step77 [label="77", URL="\ref step_77", tooltip="Newton's method for a nonlinear elliptic problem. Uses %SUNDIALS'",height=.35,width=.35,fillcolor="orange"];
  Step78 [label="78", URL="\ref step_78", tooltip="Black-Scholes equation for stock options.",height=.35,width=.35,fillcolor="orange"];
  Step79 [label="79", URL="\ref step_79", tooltip="Topology optimization of elastic media.",height=.25,width=.25,fillcolor="lightblue"];
  Step8 [label="8", URL="\ref step_8", tooltip="Systems of PDE. Elasticity. Tensors.",height=.35,width=.35,fillcolor="orange"];
  Step9 [label="9", URL="\ref step_9", tooltip="Advection equation. Multithreading. Refinement criteria.",height=.35,width=.35,fillcolor="orange"];
  Step7 -> Step10 [color="orange",weight=5,];
  Step10 -> Step11 [color="orange",weight=5,];
  Step7 -> Step12 [color="orange",weight=5,];
  Step16 -> Step12b [color="orange",weight=5,];
  Step7 -> Step12b [color="orange",weight=5,];
  Step39 -> Step12b [color="orange",weight=5,];
  Step6 -> Step13 [];
  Step13 -> Step14 [color="orange",weight=5,];
  Step6 -> Step15 [];
  Step6 -> Step16 [];
  Step16 -> Step16b [color="orange",weight=5,];
  Step8 -> Step17 [];
  Step17 -> Step18 [color="lightblue",weight=5,];
  Step6 -> Step19 [];
  Step1 -> Step2 [color="green",weight=100,];
  Step4 -> Step20 [];
  Step20 -> Step21 [color="yellow2",weight=5,];
  Step6 -> Step22 [];
  Step21 -> Step22 [color="yellow2",weight=5,];
  Step4 -> Step23 [];
  Step23 -> Step24 [color="dodgerblue1",weight=5,];
  Step24 -> Step25 [color="dodgerblue1",weight=5,];
  Step6 -> Step26 [];
  Step6 -> Step27 [];
  Step6 -> Step28 [];
  Step4 -> Step29 [];
  Step2 -> Step3 [color="green",weight=100,];
  Step12 -> Step30 [color="orange",weight=5,];
  Step22 -> Step31 [color="yellow2",weight=5,];
  Step31 -> Step32 [color="yellow2",weight=5,];
  Step55 -> Step32 [color="yellow2",weight=5,];
  Step12 -> Step33 [];
  Step71 -> Step33 [];
  Step4 -> Step34 [];
  Step22 -> Step35 [color="yellow2",weight=5,];
  Step4 -> Step36 [];
  Step16 -> Step37 [color="orange",weight=5,];
  Step40 -> Step37 [];
  Step34 -> Step38 [];
  Step12b -> Step39 [color="orange",weight=5,];
  Step3 -> Step4 [color="green",weight=100,];
  Step40a [style="invis"];  Step40b [style="invis"];  Step40c [style="invis"];  Step6 -> Step40a [style="invis"];  Step40a -> Step40b [style="invis"];  Step40b -> Step40c [style="invis"];  Step40c -> Step40 [style="invis"];  Step6 -> Step40 [weight=100,color="green"];  Step15 -> Step41 [];
  Step41 -> Step42 [color="lightblue",weight=5,];
  Step40 -> Step42 [];
  Step31 -> Step43 [color="yellow2",weight=5,];
  Step18 -> Step44 [color="lightblue",weight=5,];
  Step6 -> Step45 [];
  Step8 -> Step46 [color="orange",weight=5,];
  Step22 -> Step46 [];
  Step27 -> Step46 [color="orange",weight=5,];
  Step12 -> Step47 [color="orange",weight=5,];
  Step25 -> Step48 [color="dodgerblue1",weight=5,];
  Step37 -> Step48 [];
  Step1 -> Step49 [];
  Step4 -> Step5 [color="green",weight=100,];
  Step16 -> Step50 [color="orange",weight=5,];
  Step37 -> Step50 [color="orange",weight=5,];
  Step7 -> Step51 [color="orange",weight=5,];
  Step9 -> Step51 [color="orange",weight=5,];
  Step61 -> Step51 [color="orange",weight=5,];
  Step26 -> Step52 [color="dodgerblue1",weight=5,];
  Step49 -> Step53 [color="orange",weight=5,];
  Step53 -> Step54 [color="orange",weight=5,];
  Step40 -> Step55 [];
  Step22 -> Step55 [color="yellow2",weight=5,];
  Step16 -> Step56 [color="orange",weight=5,];
  Step22 -> Step56 [];
  Step15 -> Step57 [];
  Step22 -> Step57 [color="yellow2",weight=5,];
  Step26 -> Step58 [color="dodgerblue1",weight=5,];
  Step29 -> Step58 [];
  Step37 -> Step59 [color="orange",weight=5,];
  Step5 -> Step6 [color="green",weight=100,];
  Step6 -> Step60 [];
  Step51 -> Step61 [color="orange",weight=5,];
  Step8 -> Step62 [color="orange",weight=5,];
  Step40 -> Step62 [];
  Step16 -> Step63 [color="orange",weight=5,];
  Step7 -> Step64 [color="orange",weight=5,];
  Step37 -> Step64 [color="orange",weight=5,];
  Step49 -> Step65 [color="orange",weight=5,];
  Step15 -> Step66 [color="orange",weight=5,];
  Step37 -> Step66 [color="orange",weight=5,];
  Step33 -> Step67 [color="yellow2",weight=5,];
  Step48 -> Step67 [];
  Step59 -> Step67 [];
  Step19 -> Step68 [color="orange",weight=5,];
  Step33 -> Step69 [color="yellow2",weight=5,];
  Step40 -> Step69 [];
  Step6 -> Step7 [];
  Step19 -> Step70 [];
  Step32 -> Step70 [color="yellow2",weight=5,];
  Step60 -> Step70 [];
  Step71 -> Step72 [color="orange",weight=5,];
  Step15 -> Step72 [color="orange",weight=5,];
  Step12 -> Step74 [color="orange",weight=5,];
  Step27 -> Step75 [color="orange",weight=5,];
  Step37 -> Step75 [color="orange",weight=5,];
  Step40 -> Step75 [];
  Step67 -> Step76 [color="yellow2",weight=5,];
  Step15 -> Step77 [color="orange",weight=5,];
  Step26 -> Step78 [];
  Step8 -> Step79 [];
  Step15 -> Step79 [];
  Step6 -> Step8 [];
  Step6 -> Step9 [];
}
 * @enddot
 *
 * <b>Legend:</b><br />
 * @dot
graph StepsDescription
{
  overlap=false;
  edge [fontname="FreeSans",
        fontsize="10",
        labelfontname="FreeSans",
        labelfontsize="10",
        color="black",
        style="solid"];
  node [fontname="FreeSans",
        fontsize="10",
        shape="rectangle",
        height=0.2,
        width=0.4,
        color="black",
        fillcolor="white",
        style="filled"];
  code_gallery [label="" ,height=.08,width=.125,shape="circle", fillcolor="white"];
  fake_code_gallery [label="Code gallery", shape=plaintext];
  code_gallery -- fake_code_gallery [style=dotted, arrowhead=odot, arrowsize=1];
  time_dependent [label="" ,height=.25,width=.25, fillcolor="dodgerblue1"];
  fake_time_dependent [label="Time dependent problems", shape=plaintext];
  time_dependent -- fake_time_dependent [style=dotted, arrowhead=odot, arrowsize=1];
  unfinished [label="" ,height=.25,width=.25,style="dashed", fillcolor="white"];
  fake_unfinished [label="Unfinished codes", shape=plaintext];
  unfinished -- fake_unfinished [style=dotted, arrowhead=odot, arrowsize=1];
  solids [label="" ,height=.25,width=.25, fillcolor="lightblue"];
  fake_solids [label="Solid mechanics", shape=plaintext];
  solids -- fake_solids [style=dotted, arrowhead=odot, arrowsize=1];
  fluids [label="" ,height=.25,width=.25, fillcolor="yellow2"];
  fake_fluids [label="Fluid dynamics", shape=plaintext];
  fluids -- fake_fluids [style=dotted, arrowhead=odot, arrowsize=1];
  basic [label="" ,height=.8,width=.8,shape="octagon", fillcolor="green"];
  fake_basic [label="Basic techniques", shape=plaintext];
  basic -- fake_basic [style=dotted, arrowhead=odot, arrowsize=1];
  techniques [label="" ,height=.35,width=.35, fillcolor="orange"];
  fake_techniques [label="Advanced techniques", shape=plaintext];
  techniques -- fake_techniques [style=dotted, arrowhead=odot, arrowsize=1];
  basic -- techniques [style=invis];
  techniques -- fluids [style=invis];
  fluids -- solids [style=invis];
  solids -- time_dependent [style=invis];
  time_dependent -- unfinished [style=invis];
  unfinished -- code_gallery [style=invis];
  {rank=same; basic, techniques, fluids, solids, time_dependent, unfinished, code_gallery}}
 * @enddot
 *
 * <a name="list"></a>
 * <h3>Tutorial programs listed by number</h3>
 *
 * <table align="center" width="90%">
 *   <tr valign="top">
 *       <td width="100px">step-1</td>
 *       <td> Creating a grid. A simple way to write it to a file.
 *       <br/> Keywords: Triangulation, GridGenerator::hyper_cube(),
 *       GridGenerator::hyper_shell(), GridOut,
 *       Triangulation::execute_coarsening_and_refinement()
 *       </td></tr>
 *
 *   <tr valign="top">
 *       <td>step-2</td>
 *       <td> Associate degrees of freedom to
 *       each vertex and compute the resulting sparsity pattern of
 *       matrices. Show that renumbering reduces the bandwidth of
 *       matrices significantly, i.e. clusters nonzero entries around the
 *       diagonal.
 *       <br/> Keywords: FE_Q, DynamicSparsityPattern,
 *       DoFTools::make_sparsity_pattern(), DoFHandler::distribute_dofs(),
 *       DoFRenumbering, SparsityPattern
 *       </td></tr>
 *
 *   <tr valign="top">
 *       <td>step-3</td>
 *       <td> Actually solve Laplace's
 *       problem. Object-orientation. Assembling matrices and
 *       vectors. Boundary values.
 *       <br/> Keywords: FEValues, VectorTools::interpolate_boundary_values(),
 *       MatrixTools::apply_boundary_values(), SolverCG, Vector,
 *       SparseMatrix, DataOut
 *       </td></tr>
 *
 *   <tr valign="top">
 *       <td>step-4</td>
 *       <td> This example is programmed in a
 *       way that it is independent of the dimension for which we want to
 *       solve Laplace's equation; we will solve the equation in 2D and
 *       3D, although the program is exactly the same. Non-constant right
 *       hand side function. Non-homogeneous boundary values.
 *       <br/> Keywords: VectorTools::point_value(),
 *       VectorTools::compute_mean_value()
 *       </td></tr>
 *
 *   <tr valign="top">
 *       <td>step-5</td>
 *       <td> Computations on successively
 *       refined grids. Reading a grid from disk. Some optimizations.
 *       Using assertions. Non-constant coefficient in
 *       the elliptic operator (yielding the extended Poisson
 *       equation). Preconditioning the CG solver for the
 *       linear system of equations.
 *       <br/> Keywords: PreconditionSSOR, GridIn, SphericalManifold
 *       </td></tr>
 *
 *   <tr valign="top">
 *       <td>step-6</td>
 *       <td> Adaptive local
 *       refinement. Handling of hanging nodes. Higher order elements.
 *       Catching exceptions in the <code>main</code> function.
 *       <br/> Keywords: DoFTools::make_hanging_node_constraints(),
 *       AffineConstraints::distribute_local_to_global(), KellyErrorEstimator,
 *       GridRefinement::refine_and_coarsen_fixed_number()
 *       </td></tr>
 *
 *   <tr valign="top">
 *       <td>step-7</td>
 *       <td> Helmholtz
 *       equation. Non-homogeneous Neumann boundary conditions and
 *       boundary integrals. Verification of correctness of computed
 *       solutions. Computing the error between exact and numerical
 *       solution and output of the data in tables. Using counted pointers.
 *       <br/> Keywords: FEFaceValues, VectorTools::integrate_difference(),
 *       VectorTools::compute_global_error(), TableHandler
 *       </td></tr>
 *
 *   <tr valign="top">
 *       <td>step-8</td>
 *       <td> The elasticity equations will be
 *       solved instead of Laplace's equation. The solution is
 *       vector-valued and the equations form a system with as many
 *       equations as the dimension of the space in which it is posed.
 *       <br/> Keywords: FESystem
 *       </td></tr>
 *
 *   <tr valign="top">
 *       <td>step-9</td>
 *       <td> Linear advection equation, assembling
 *       the system of equations in parallel using multi-threading,
 *       implementing a refinement criterion based on a finite difference
 *       approximation of the gradient.
 *       <br/> Keywords: TensorFunction, WorkStream::run(), SolverGMRES
 *       </td></tr>
 *
 *   <tr valign="top">
 *       <td>step-10</td>
 *       <td> Higher order mappings. Do not
 *       solve equations, but rather compute the value of pi to high
 *       accuracy.
 *       <br/> Keywords: MappingQ, FE_Nothing, ConvergenceTable, GridOut, FEFaceValues
 *       </td></tr>
 *
 *   <tr valign="top">
 *       <td>step-11</td>
 *       <td> Solving a Laplace problem with
 *       higher order mappings. Using mean value constraints and
 *       intermediate representations of sparsity patterns.
 *       <br/> Keywords: AffineConstraints, DoFTools::extract_boundary_dofs(),
 *       TableHandler
 *       </td></tr>
 *
 *   <tr valign="top">
 *       <td>step-12</td>
 *       <td> Discontinuous Galerkin methods for linear advection problems.
 *       <br/> Keywords: FEInterfaceValues, MeshWorker::mesh_loop(),
 *       DoFTools::make_flux_sparsity_pattern()
 *       </td></tr>
 *
 *   <tr valign="top">
 *       <td>step-13</td>
 *       <td> Software design questions and
 *       how to write a modular, extensible finite element program.
 *       </td></tr>
 *
 *   <tr valign="top">
 *       <td>step-14</td>
 *       <td> Duality based error estimators,
 *       more strategies to write a modular, extensible finite element
 *       program.
 *       <br/> Keywords: KellyErrorEstimator
 *       </td></tr>
 *
 *   <tr valign="top">
 *       <td>step-15</td>
 *       <td> A nonlinear elliptic problem: The minimal surface equation.
 *       Newton's method. Transferring a solution across mesh refinement.
 *       <br/> Keywords: SolutionTransfer
 *       </td></tr>
 *
 *   <tr valign="top">
 *       <td>step-16</td>
 *       <td> Multigrid preconditioning of the Laplace equation on adaptive
 *       meshes.
 *       <br/> Keywords: Multigrid, PreconditionMG, mg::Matrix,
 *       MGTransferPrebuilt, MeshWorker::mesh_loop(), MGLevelObject,
 *       MGConstrainedDoFs
 *       </td></tr>
 *
 *   <tr valign="top">
 *       <td>step-16b</td>
 *       <td> A variant of step-16 but with MeshWorker for assembly: Multigrid
 *       preconditioning of the Laplace equation on adaptive meshes.
 *       </td></tr>
 *
 *   <tr valign="top">
 *       <td>step-17</td>
 *       <td> Using PETSc for linear algebra; running
 *       in parallel on clusters of computers linked together by MPI.
 *       <br/> Keywords: PETScWrappers::MPI::SparseMatrix, ConditionalOStream,
 *       PETScWrappers::PreconditionBlockJacobi
 *       </td></tr>
 *
 *   <tr valign="top">
 *       <td>step-18</td>
 *       <td> A time dependent problem; using a much
 *       simplified version of implementing elasticity; moving meshes; handling
 *       large scale output of parallel programs. Simple implicit (backward
 *       Euler) time stepping.
 *       <br/> Keywords: parallel::shared::Triangulation,
 *       DataOutInterface::write_vtu_with_pvtu_record()
 *       </td></tr>
 *
 *   <tr valign="top">
 *       <td>step-19</td>
 *       <td> Coupling particles to the solution of partial differential equations.
 *       <br/> Keywords: Particles
 *       </td></tr>
 *
 *   <tr valign="top">
 *       <td>step-20</td>
 *       <td> Mixed finite elements. Using block
 *       matrices and block vectors to define more complicated solvers and
 *       preconditioners working on the Schur complement.
 *       <br/> Keywords: FEValuesExtractors, LinearOperator, TensorFunction,
 *       FE_RaviartThomas
 *       </td></tr>
 *
 *   <tr valign="top">
 *       <td>step-21</td>
 *       <td> The time dependent two-phase flow in
 *       porous media. Extensions of mixed Laplace discretizations. More
 *       complicated block solvers. Simple explicit (forward Euler) time
 *       stepping.
 *       <br/> Keywords: TensorFunction, FE_RaviartThomas,
 *       VectorTools::project(), DiscreteTime
 *       </td></tr>
 *
 *   <tr valign="top">
 *       <td>step-22</td>
 *       <td> Solving the Stokes equations of slow fluid flow on adaptive
 *       meshes. More on Schur complement solvers. Advanced use of the
 *       AffineConstraints class.
 *       <br/> Keywords: AffineConstraints,
 *       VectorTools::compute_no_normal_flux_constraints(), SparseILU,
 *       SparseDirectUMFPACK, BlockDynamicSparsityPattern
 *       </td></tr>
 *
 *   <tr valign="top">
 *       <td>step-23</td>
 *       <td> Finally a "real" time dependent problem, the wave equation.
 *       Fractional time stepping (explicit, fully implicit and Crank-Nicholson
 *       method).
 *       <br/> Keywords: MatrixCreator, VectorTools::project()
 *       </td></tr>
 *
 *   <tr valign="top">
 *       <td>step-24</td>
 *       <td> A variant of step-23 with absorbing
 *       boundary conditions, and extracting practically useful data.
 *       Implicit time stepping.
 *       <br/> Keywords: VectorTools::point_value()
 *       </td></tr>
 *
 *   <tr valign="top">
 *       <td>step-25</td>
 *       <td> The sine-Gordon
 *       soliton equation, which is a nonlinear variant of the time
 *       dependent wave equation covered in step-23 and step-24.
 *       Fractional time stepping.
 *       <br/> Keywords: FunctionTime, VectorTools::integrate_difference()
 *       </td></tr>
 *
 *   <tr valign="top">
 *       <td>step-26</td>
 *       <td> The heat equation, solved on a mesh that is adapted
 *       every few time steps. Fractional time stepping.
 *       <br/> Keywords: KellyErrorEstimator, SolutionTransfer,
 *       VectorTools::interpolate(), VectorTools::create_right_hand_side()
 *       </td></tr>
 *
 *   <tr valign="top">
 *       <td>step-27</td>
 *       <td> The hp-finite element method.
 *       <br/> Keywords: hp::FECollection, hp::QCollection, hp::Refinement,
 *       FESeries::Fourier, Triangulation::create_triangulation()
 *       </td></tr>
 *
 *   <tr valign="top">
 *       <td>step-28</td>
 *       <td> Multiple grids for solving a multigroup diffusion equation
 *       in nuclear physics simulating a nuclear reactor core.
 *       </td></tr>
 *
 *   <tr valign="top">
 *       <td>step-29</td>
 *       <td> Solving a complex-valued Helmholtz equation. Sparse direct
 *       solvers. Dealing with parameter files.  </td></tr>
 *
 *   <tr valign="top">
 *       <td>step-30</td>
 *       <td> Anisotropic refinement for DG finite element methods.
 *       </td></tr>
 *
 *   <tr valign="top">
 *       <td>step-31</td>
 *       <td> Time-dependent Stokes flow driven by temperature
 *       differences in a fluid. Adaptive meshes that change between time
 *       steps. Implicit/explicit time stepping.
 *       </td></tr>
 *
 *   <tr valign="top">
 *       <td>step-32</td>
 *       <td> A massively parallel solver for time-dependent Stokes flow driven
 *       by temperature differences in a fluid. Adapting methods for real-world
 *       equations. Implicit/explicit time stepping.
 *       </td></tr>
 *
 *   <tr valign="top">
 *       <td>step-33</td>
 *       <td> A nonlinear hyperbolic conservation law: The Euler equations of
 *       compressible gas dynamics. Fractional time stepping.
 *       </td></tr>
 *
 *   <tr valign="top">
 *       <td>step-34</td>
 *       <td> Boundary element methods (BEM) of low order: Exterior irrotational
 *       flow. The ParsedFunction class.
 *       </td></tr>
 *
 *   <tr valign="top">
 *       <td>step-35</td>
 *       <td> A projection solver for the Navier&ndash;Stokes equations.
 *       </td></tr>
 *
 *   <tr valign="top">
 *       <td>step-36</td>
 *       <td> Using SLEPc for linear algebra; solving an eigenspectrum
 *       problem. The Schrödinger wave equation.
 *       </td></tr>
 *
 *   <tr valign="top">
 *       <td>step-37</td>
 *       <td> Solving a Poisson problem with a multilevel preconditioner without
 *       explicitly storing the matrix (a matrix-free method) in a massively
 *       parallel context.
 *       </td></tr>
 *
 *   <tr valign="top">
 *       <td>step-38</td>
 *       <td>Solving the Laplace-Beltrami equation on curved manifolds embedded
 *       in higher dimensional spaces.
 *       </td></tr>
 *
 *   <tr valign="top">
 *       <td>step-39</td>
 *       <td> Solving Poisson's equation once more, this time with the
 *       interior penalty method, one of the discontinuous Galerkin
 *       methods developed for this problem. Error estimator, adaptive
 *       meshes, and multigrid preconditioner, all using the MeshWorker
 *       framework.
 *       </td></tr>
 *
 *   <tr valign="top">
 *       <td>step-40</td>
 *       <td> Techniques for the massively parallel solution of the Laplace
 *       equation (up to 10,000s of processors).
 *       </td></tr>
 *
 *   <tr valign="top">
 *       <td>step-41</td>
 *       <td> Solving the obstacle problem, a variational inequality.
 *       </td></tr>
 *
 *   <tr valign="top">
 *       <td>step-42</td>
 *       <td> A solver for an elasto-plastic contact problem, running on
 *       parallel machines.
 *       </td></tr>
 *
 *   <tr valign="top">
 *       <td>step-43</td>
 *       <td> Advanced techniques for the simulation of porous media flow.
 *       Explicit time stepping.
 *       </td></tr>
 *
 *   <tr valign="top">
 *       <td>step-44</td>
 *       <td> Finite strain hyperelasticity based on a three-field formulation.
 *       Implicit time stepping.
 *       <br/> Keywords: CellDataStorage, FEValuesExtractors, WorkStream::run,
 *       BlockSparseMatrix, BlockVector, ComponentSelectFunction,
 *       Physics::Elasticity, FullMatrix::extract_submatrix_from(),
 *       FullMatrix::scatter_matrix_to(), LinearOperator, SolverSelector,
 *       PreconditionSelector, ReductionControl, MappingQEulerian
 *       </td></tr>
 *
 *   <tr valign="top">
 *       <td>step-45</td>
 *       <td> Periodic boundary conditions.
 *       <br/> Keywords: GridTools::collect_periodic_faces(),
 *       GridTools::PeriodicFacePair, Triangulation::add_periodicity()
 *       </td></tr>
 *
 *   <tr valign="top">
 *       <td>step-46</td>
 *       <td> Coupling different kinds of equations in different parts of the domain.
 *       </td></tr>
 *
 *   <tr valign="top">
 *       <td>step-47</td>
 *        <td> Solving the fourth-order biharmonic equation using the $C^0$
 *        Interior Penalty (C0IP) method.
 *        <br/> Keywords: FEInterfaceValues
 *       </td></tr>
 *
 *   <tr valign="top">
 *       <td>step-48</td>
 *       <td> Explicit time stepping for the Sine&ndash;Gordon equation based on
 *       a diagonal mass matrix. Efficient implementation of (nonlinear) finite
 *       element operators.
 *       </td></tr>
 *
 *   <tr valign="top">
 *       <td>step-49</td>
 *       <td> Advanced mesh creation and manipulation techniques.
 *       </td></tr>
 *
 *   <tr valign="top">
 *       <td>step-50</td>
 *       <td> Geometric multigrid on adaptive meshes distributed in parallel.
 *       <br/> Keywords: Multigrid, MGLevelObject, MGConstrainedDoFs, IndexSet, MGTools, PreconditionMG, MatrixFree, FEInterfaceValues, MeshWorker::mesh_loop()
 *       </td></tr>
 *
 *   <tr valign="top">
 *       <td>step-51</td>
 *       <td> Solving the convection-diffusion equation with a hybridizable
 *       discontinuous Galerkin method using face elements.
 *       </td></tr>
 *
 *   <tr valign="top">
 *       <td>step-52</td>
 *       <td> Solving the time dependent neutron diffusion equation using
 *       Runge-Kutta methods. Explicit and implicit time stepping.
 *       </td></tr>
 *
 *   <tr valign="top">
 *       <td>step-53</td>
 *       <td> Describing the geometry of complex domains and curved boundaries.
 *       </td></tr>
 *
 *   <tr valign="top">
 *       <td>step-54</td>
 *       <td> Using CAD files to describe the boundary of your domain.
 *       <br/> Keywords: Manifold, OpenCASCADE::read_IGES(),
 *       OpenCASCADE::NormalProjectionBoundary
 *       </td></tr>
 *
 *   <tr valign="top">
 *       <td>step-55</td>
 *       <td> Solving the Stokes problem in parallel.
 *       </td></tr>
 *
 *   <tr valign="top">
 *       <td>step-56</td>
 *       <td> Geometric Multigrid for Stokes.
 *       </td></tr>
 *
 *   <tr valign="top">
 *       <td>step-57</td>
 *       <td> Incompressible, stationary Navier Stokes equations.
 *       </td></tr>
 *
 *   <tr valign="top">
 *       <td>step-58</td>
 *       <td> The nonlinear Schr&ouml;dinger equation.
 *       </td></tr>
 *
 *   <tr valign="top">
 *       <td>step-59</td>
 *       <td> Solving a Poisson problem discretized with an interior penalty DG
 *       method and a multilevel preconditioner in a matrix-free fashion using
 *       a massively parallel implementation.
 *       </td></tr>
 *
 *   <tr valign="top">
 *       <td>step-60</td>
 *       <td> Distributed Lagrange multipliers for the solution of
 *       Poisson problems in complex domains with constraints defined
 *       on non-matching grids.
 *       </td></tr>
 *
 *   <tr valign="top">
 *       <td>step-61</td>
 *       <td> Solving the Poisson problem with the "weak Galerkin" finite element
 *       method.
 *       </td></tr>
 *
 *   <tr valign="top">
 *       <td>step-62</td>
 *       <td> Resonance frequency and bandgap of a phononic crystal. Elastic
 *       wave equation in the frequency domain with Perfectly Matched Layer
 *       boundary conditions. Parallelization via MUMPS and MPI.
 *       </td></tr>
 *
 *   <tr valign="top">
 *       <td>step-63</td>
 *       <td>Block smoothers for geometric multigrid. A scalar convection
 *       diffusion equation is solved with different additive or
 *       multiplicative multigrid smoothers.
 *       <br/> Keywords: Multigrid, MeshWorker::mesh_loop(),
 *       MGSmootherPrecondition, RelaxationBlock, DoFRenumbering::downstream()
 *       </td></tr>
 *
 *   <tr valign="top">
 *       <td>step-64</td>
 *       <td> Solving a Helmholtz problem using matrix-free methods on the GPU
 *       with MPI parallelization.
 *       </td></tr>
 *
 *   <tr valign="top">
 *       <td>step-65</td>
 *        <td> The TransfiniteInterpolationManifold and MappingQCache classes for
 *        advanced manifold operations.
 *       </td></tr>
 *
 *   <tr valign="top">
 *       <td>step-66</td>
 *       <td> A matrix-free geometric multigrid solver for a nonlinear problem.
 *       <br/> Keywords: MatrixFree, Multigrid,
 *       MGTransferMatrixFree::interpolate_to_mg(), MatrixFreeTools::compute_diagonal(),
 *       TransfiniteInterpolationManifold, MappingQGeneric
 *       </td></tr>
 *
 *   <tr valign="top">
 *       <td>step-67</td>
 *        <td> Solving the Euler equations of compressible gas dynamics with an
 *        explicit time integrator and high-order discontinuous Galerkin
 *        methods based on matrix-free implementations.
 *       </td></tr>
 *
 *   <tr valign="top">
 *       <td>step-68</td>
 *        <td> Simulation of the motion of massless tracer particles in a vortical flow.
 *             Parallel simulation of the advection of particles with load balancing.
 *        <br/> Keywords: Particles
 *       </td></tr>
 *
 *   <tr valign="top">
 *       <td>step-69</td>
 *        <td> Hyperbolic conservation laws: a first-order guaranteed maximum
 *        wavespeed method for the compressible Euler equations. Explicit time
 *        stepping.
 *       </td></tr>
 *
 *   <tr valign="top">
 *       <td>step-70</td>
 *       <td> A fluid structure interaction problem on fully distributed
 *       non-matching grids, using penalty methods, and a coupling constructed
 *       through a ParticleHandler object.
 *       <br/> Keywords: Particles
 *       </td></tr>
 *
 *   <tr valign="top">
 *       <td>step-71</td>
 *       <td> Constitutive modelling: a demonstration of how automatic and symbolic 
 *       differentiation can be used to rapidly implement a complex coupled constitutive
 *       law.
 *       <br/> Keywords: Automatic differentiation, Symbolic differentiation, 
 *       Constitutive modelling
 *       </td></tr>
 *
 *   <tr valign="top">
 *       <td>step-72</td>
 *       <td> A nonlinear elliptic problem: The minimal surface equation.
 *       Newton's method. Using automatic differentiation to linearize the 
 *       residual, or to compute the full linear system from an energy
 *       functional.
 *       <br/> Keywords: Automatic differentiation, MeshWorker::mesh_loop()
 *       </td></tr>
 *
 *   <tr valign="top">
 *       <td>step-74</td>
 *       <td> The Symmetric interior penalty Galerkin (SIPG) method for Poisson's equation.
 *       <br/> Keywords: MeshWorker::mesh_loop(), FEInterfaceValues, ConvergenceTable
 *       </td></tr>
 *
 *   <tr valign="top">
 *       <td>step-75</td>
 *       <td> Solving the Laplace equation on hp-adaptive meshes with MatrixFree
 *       methods on thousands of processors.
 *       <br/> Keywords: parallel::distributed::Triangulation, hp::Refinement, MatrixFree
 *       </td></tr>
 *
 *   <tr valign="top">
 *       <td>step-76</td>
 *       <td> Like step-67, but for demonstrating MPI-3.0 shared-memory features.
 *       <br/> Keywords: MatrixFree, MatrixFree::cell_loop()
 *       </td></tr>
 *
 *   <tr valign="top">
 *       <td>step-77</td>
 *       <td> The nonlinear minimal surface equation revisited.
 *       Interfacing with SUNDIALS' KINSOL nonlinear solver.
 *       </td></tr>
 *
 *   <tr valign="top">
 *       <td>step-78</td>
 *       <td> Solves the Black-Scholes equation for options pricing in 1-D.
 *       </td></tr>
 *
 *   <tr valign="top">
 *       <td>step-79</td>
 *       <td> A topology optimization program for elastic media using the solid
 *       isotropic material penalization (SIMP) formulation.
 *       </td></tr>
 *
 * </table>
 *
 *
 * <a name="topic"></a>
 * <h3>Tutorial programs grouped by topics</h3>
 *
 * <h4><b>Basic techniques</b></h4>
 * <table align="center" width="90%">
 *
 *   <tr valign="top">
 *     <td width="400px"> Creating a grid. A simple way to write it to a file
 *     <td>step-1</td>
 *     </td>
 *   </tr>
 *
 *   <tr valign="top">
 *     <td> Degrees of freedom
 *     <td>step-2</td>
 *     </td>
 *   </tr>
 *
 *   <tr valign="top">
 *     <td> Solving the Laplace equation
 *     <td>step-3</td>
 *     </td>
 *   </tr>
 *
 *   <tr valign="top">
 *     <td> Dimension independent programming, non-zero data
 *     <td>step-4</td>
 *     </td>
 *   </tr>
 *
 *   <tr valign="top">
 *     <td> Computing on uniformly refined meshes
 *     <td>step-5</td>
 *     </td>
 *   </tr>
 *
 *   <tr valign="top">
 *     <td> Adaptivity
 *     <td>step-6, step-26</td>
 *     </td>
 *   </tr>
 *
 *   <tr valign="top">
 *     <td> Evaluating errors
 *     <td>step-7</td>
 *     </td>
 *
 *   <tr valign="top">
 *     <td> Nonlinear problems, Newton's method
 *     </td>
 *     <td>step-15</td>
 *   </tr>
 *
 * </table>
 * <h4><b>Advanced techniques</b></h4>
 * <table align="center" width="90%">
 *
 *   <tr valign="top">
 *     <td width="400px"> Multithreading
 *     </td>
 *     <td>
 *       step-9,
 *       step-28,
 *       step-32,
 *       step-44,
 *       step-48,
 *       step-51,
 *       step-69
 *     </td>
 *   </tr>
 *
 *   <tr valign="top">
 *     <td> Block solvers and preconditioners
 *     </td>
 *     <td>
 *       step-20,
 *       step-21,
 *       step-22,
 *       step-31,
 *       step-32,
 *       step-43,
 *       step-44,
 *       step-55,
 *       step-56,
 *       step-57,
 *       step-60,
 *       step-70
 *     </td>
 *   </tr>
 *
 *   <tr valign="top">
 *     <td> Using Trilinos
 *     </td>
 *     <td>
 *       step-31,
 *       step-32,
 *       step-33,
 *       step-40,
 *       step-41,
 *       step-42,
 *       step-43,
 *       step-50,
 *       step-55,
 *       step-71,
 *       step-72,
 *       step-75
 *     </td>
 *   </tr>
 *
 *   <tr valign="top">
 *     <td> Parallelization via PETSc and MPI
 *     </td>
 *     <td>
 *       step-17,
 *       step-18,
 *       step-40,
 *       step-50,
 *       step-55
 *     </td>
 *   </tr>
 *
 *   <tr valign="top">
 *     <td> Parallelization via Trilinos and MPI
 *     </td>
 *     <td>
 *       step-32,
 *       step-40,
 *       step-42,
 *       step-50,
 *       step-55,
 *       step-70
 *     </td>
 *   </tr>
 *
 *   <tr valign="top">
 *     <td> Parallelization via MUMPS and MPI
 *     </td>
 *     <td>
 *       step-62
 *     </td>
 *   </tr>
 *
 *   <tr valign="top">
 *     <td> Parallelization via CUDA and MPI
 *     </td>
 *     <td>
 *       step-64
 *     </td>
 *   </tr>
 *
 *   <tr valign="top">
 *     <td> Parallelization on very large numbers of processors
 *     </td>
 *     <td>
 *       step-32,
 *       step-37,
 *       step-40,
 *       step-42,
 *       step-50,
 *       step-55,
 *       step-59,
 *       step-66,
 *       step-67,
 *       step-69,
 *       step-70,
 *       step-75
 *     </td>
 *   </tr>
 *
 *   <tr valign="top">
 *     <td> Input parameter handling
 *     </td>
 *     <td>
 *       step-28,
 *       step-29,
 *       step-32,
 *       step-33,
 *       step-34,
 *       step-35,
 *       step-36,
 *       step-42,
 *       step-44,
 *       step-60,
 *       step-62,
 *       step-69,
 *       step-70,
 *       step-71,
 *       step-72
 *     </td>
 *   </tr>
 *
 *   <tr valign="top">
 *     <td> Higher order mappings
 *     </td>
 *     <td>
 *       step-10,
 *       step-11,
 *       step-32,
 *       step-60,
 *       step-65,
 *       step-66,
 *       step-67
 *     </td>
 *   </tr>
 *
 *   <tr valign="top">
 *     <td> Error indicators and estimators
 *     </td>
 *     <td>
 *       step-6,
 *       step-9,
 *       step-14,
 *       step-39,
 *       step-50,
 *       step-74
 *     </td>
 *   </tr>
 *
 *   <tr valign="top">
 *     <td> Transferring solutions across mesh refinement
 *     </td>
 *     <td>
 *       step-15,
 *       step-28,
 *       step-31,
 *       step-32,
 *       step-33,
 *       step-42,
 *       step-43,
 *       step-57,
 *       step-70
 *     </td>
 *   </tr>
 *
 *   <tr valign="top">
 *     <td> Discontinuous Galerkin methods
 *     </td>
 *     <td>
 *       step-12,
 *       step-21,
 *       step-39,
 *       step-46,
 *       step-47,
 *       step-51,
 *       step-59,
 *       step-61,
 *       step-67,
 *       step-74
 *     </td>
 *   </tr>
 *
 *   <tr valign="top">
 *     <td> hp-finite elements
 *     </td>
 *     <td>
 *       step-27,
 *       step-46,
 *       step-75
 *     </td>
 *   </tr>
 *
 *   <tr valign="top">
 *     <td> Anisotropic refinement for DG finite element methods
 *     </td>
 *     <td>step-30</td>
 *   </tr>
 *
 *   <tr valign="top">
 *     <td> Computing Jacobians from residuals, automatic and symbolic differentiation
 *     </td>
 *     <td>
 *       step-33,
 *       step-71,
 *       step-72
 *     </td>
 *   </tr>
 *
 *   <tr valign="top">
 *     <td> Operator splitting
 *     </td>
 *     <td>step-21, step-31, step-32, step-58</td>
 *   </tr>
 *
 *   <tr valign="top">
 *     <td> Boundary element methods, curved manifolds
 *     </td>
 *     <td>
 *       step-32,
 *       step-34,
 *       step-38,
 *       step-53,
 *       step-54,
 *       step-65
 *     </td>
 *   </tr>
 *
 *   <tr valign="top">
 *     <td> Periodic boundary conditions
 *     </td>
 *     <td>
 *       step-45,
 *       step-59
 *     </td>
 *   </tr>
 *
 *   <tr valign="top">
 *     <td> Matrix-free methods with sum factorization
 *     </td>
 *     <td>
 *       step-37,
 *       step-48,
 *       step-59,
 *       step-64,
 *       step-66,
 *       step-67,
 *       step-75,
 *       step-76
 *     </td>
 *   </tr>
 *
 *   <tr valign="top">
 *     <td> Advanced meshes and geometries
 *     </td>
 *     <td>
 *       step-49,
 *       step-53,
 *       step-54,
 *       step-65
 *     </td>
 *   </tr>
 *
 *   <tr valign="top">
 *     <td> Non matching algorithms
 *     </td>
 *     <td>
 *       step-60,
 *       step-70
 *     </td>
 *   </tr>
 *
 *   <tr valign="top">
 *     <td> HDF5 and Python
 *     </td>
 *     <td>
 *       step-62
 *     </td>
 *   </tr>
 *
 *   <tr valign="top">
 *     <td> Particles
 *     </td>
 *     <td>
 *     step-19,
 *     step-68,
 *     step-70
 *     </td>
 *   </tr>
 *
 * </table>
 * <h4><b>Linear solvers</b></h4>
 * <table align="center" width="90%">
 *
 *   <tr valign="top">
 *     <td width="400px"> Conjugate Gradient solver
 *     </td>
 *     <td>step-3</td>
 *   </tr>
 *
 *   <tr valign="top">
 *     <td> Preconditioned CG solver
 *     </td>
 *     <td>step-5</td>
 *   </tr>
 *
 *   <tr valign="top">
 *     <td> BiCGStab
 *     </td>
 *     <td>step-9</td>
 *   </tr>
 *
 *   <tr valign="top">
 *     <td> Direct solvers
 *     </td>
 *     <td>
 *       step-29,
 *       step-44,
 *       step-47,
 *       step-46,
 *       step-58,
 *       step-62
 *     </td>
 *   </tr>
 *
 *   <tr valign="top">
 *     <td> Multilevel preconditioners
 *     </td>
 *     <td>
 *       step-16,
 *       step-16b
 *       step-31,
 *       step-32,
 *       step-37,
 *       step-39,
 *       step-41,
 *       step-42,
 *       step-43,
 *       step-50,
 *       step-56,
 *       step-59,
 *       step-63,
 *       step-66,
 *       step-75
 *     </td>
 *   </tr>
 *
 *   <tr valign="top">
 *     <td> Parallel solvers
 *     </td>
 *     <td>
 *       step-17,
 *       step-18,
 *       step-32,
 *       step-37,
 *       step-40,
 *       step-42,
 *       step-50,
 *       step-55,
 *       step-59,
 *       step-66,
 *       step-75
 *     </td>
 *   </tr>
 *
 *   <tr valign="top">
 *     <td> Block and Schur complement solvers
 *     </td>
 *     <td>
 *       step-20,
 *       step-21,
 *       step-22,
 *       step-31,
 *       step-32,
 *       step-43,
 *       step-55,
 *       step-56,
 *       step-57,
 *       step-60,
 *       step-70
 *     </td>
 *   </tr>
 *
 *   <tr valign="top">
 *     <td> Decoupled projection solvers
 *     </td>
 *     <td>step-35</td>
 *   </tr>
 *
 *   <tr valign="top">
 *     <td> Linear Newton systems from nonlinear equations
 *     </td>
 *     <td>
 *       step-33,
 *       step-41,
 *       step-42,
 *       step-44,
 *       step-57
 *     </td>
 *   </tr>
 *
 *   <tr valign="top">
 *     <td> Eigenvalue solvers
 *     </td>
 *     <td>step-36</td>
 *   </tr>
 *
 *   <tr valign="top">
 *     <td> Linear operators
 *     </td>
 *     <td>
 *       step-44,
 *       step-60,
 *       step-70
 *     </td>
 *   </tr>
 *
 * </table>
 *
 *
 * <h4><b>Other equations</b></h4>
 * <table align="center" width="90%">
 *
 *   <tr valign="top">
 *     <td width="400px"> Helmholtz equation
 *     </td>
 *     <td>
 *       step-7,
 *       step-29,
 *       step-62,
 *       step-64
 *     </td>
 *   </tr>
 *
 *   <tr valign="top">
 *     <td> Elasticity and elasto-plasticity equations
 *     </td>
 *     <td>
 *       step-8,
 *       step-42,
 *       step-46,
 *       step-62
 *     </td>
 *   </tr>
 *
 *   <tr valign="top">
 *     <td> Heat equation
 *     </td>
 *     <td>
 *       step-26
 *     </td>
 *   </tr>
 *
 *   <tr valign="top">
 *     <td> Minimal surface equation
 *     </td>
 *     <td>
 *       step-15
 *     </td>
 *   </tr>
 *
 *   <tr valign="top">
 *     <td> Quasi-static elasticity equations
 *     </td>
 *     <td>
 *       step-18,
 *       step-44
 *     </td>
 *   </tr>
 *
 *   <tr valign="top">
 *     <td> Transport (advection) equations
 *     </td>
 *     <td>step-9,
 *         step-21,
 *         step-31,
 *         step-32,
 *         step-43,
 *         step-51
 *     </td>
 *   </tr>
 *
 *   <tr valign="top">
 *     <td> The nonlinear hyperbolic Euler system of compressible gas dynamics
 *     </td>
 *     <td>
 *       step-33,
 *       step-67,
 *       step-69
 *     </td>
 *   </tr>
 *
 *   <tr valign="top">
 *     <td> Mixed Laplace, Darcy, Porous media
 *     </td>
 *     <td>
 *       step-20,
 *       step-21,
 *       step-43
 *     </td>
 *   </tr>
 *
 *   <tr valign="top">
 *     <td> Stokes and incompressible Navier-Stokes flow
 *     </td>
 *     <td>
 *       step-22,
 *       step-31,
 *       step-32,
 *       step-35,
 *       step-46,
 *       step-55,
 *       step-56,
 *       step-57,
 *       step-70
 *     </td>
 *   </tr>
 *
 *   <tr valign="top">
 *     <td> The wave equation, in linear and nonlinear variants
 *     </td>
 *     <td>
 *       step-23,
 *       step-24,
 *       step-25,
 *       step-48,
 *       step-58
 *     </td>
 *   </tr>
 *
 *   <tr valign="top">
 *     <td> A multigroup diffusion problem in neutron transport
 *     </td>
 *     <td>step-28</td>
 *   </tr>
 *
 *   <tr valign="top">
 *     <td> Irrotational flow
 *     </td>
 *     <td>step-34</td>
 *   </tr>
 *
 *   <tr valign="top">
 *     <td> An eigenspectrum problem
 *     </td>
 *     <td>step-36</td>
 *   </tr>
 *
 *   <tr valign="top">
 *     <td> Fourth-order biharmonic equation
 *     </td>
 *     <td>step-47</td>
 *   </tr>
 *
 *   <tr valign="top">
 *     <td> The obstacle problem, a variational inequality
 *     </td>
 *     <td>
 *       step-41,
 *       step-42
 *     </td>
 *   </tr>
 *
 *   <tr valign="top">
 *     <td> The nonlinear Schr&ouml;dinger equation
 *     </td>
 *     <td>step-58</td>
 *   </tr>
 *
 *   <tr valign="top">
 *     <td> Coupling different equations in different parts of the domain
 *     </td>
 *     <td>step-46</td>
 *   </tr>
 *
 * </table>
 *
 *
 *
 * <h4><b>%Vector problems</b></h4>
 * <table align="center" width="90%">
 *
 *   <tr valign="top">
 *     <td width="400px"> Elasticity and elasto-plasticity equations
 *     </td>
 *     <td>
 *       step-8,
 *       step-42
 *     </td>
 *   </tr>
 *
 *   <tr valign="top">
 *     <td> Mixed Laplace
 *     </td>
 *     <td>step-20</td>
 *   </tr>
 *
 *   <tr valign="top">
 *     <td> Mixed Laplace plus an advection equation
 *     </td>
 *     <td>step-21,
 *         step-43
 *     </td>
 *   </tr>
 *
 *   <tr valign="top">
 *     <td> Incompressible Stokes and Navier-Stokes flow
 *     </td>
 *     <td>step-22,
 *         step-31,
 *         step-32,
 *         step-35,
 *         step-55,
 *         step-56,
 *         step-57,
 *         step-70
 *    </td>
 *   </tr>
 *
 *   <tr valign="top">
 *     <td> A complex-valued Helmholtz problem
 *     </td>
 *     <td>step-29</td>
 *   </tr>
 *
 *   <tr valign="top">
 *     <td> The Euler equations of compressible gas dynamics
 *     </td>
 *     <td>
 *       step-33,
 *       step-67,
 *       step-69
 *     </td>
 *   </tr>
 *
 *   <tr valign="top">
 *     <td> Coupling different equations in different parts of the domain
 *     <td>step-46</td>
 *   </tr>
 *
 * </table>
 *
 *
 *
 * <h4><b>Time dependent problems</b></h4>
 * <table align="center" width="90%">
 *
 *   <tr valign="top">
 *     <td> The heat equation
 *     </td>
 *     <td>step-26
 *     </td>
 *   </tr>
 *
 *   <tr valign="top">
 *     <td width="400px"> Quasi-static elasticity
 *     </td>
 *     <td>
 *      step-18,
 *      step-44
 *     </td>
 *   </tr>
 *
 *   <tr valign="top">
 *     <td> Porous media flow
 *     </td>
 *     <td>step-21,
 *         step-43
 *     </td>
 *   </tr>
 *
 *   <tr valign="top">
 *     <td> The wave equation, in linear and nonlinear variants
 *     </td>
 *     <td>step-23,
 *         step-24,
 *         step-25,
 *         step-48,
 *         step-58
 *     </td>
 *   </tr>
 *
 *   <tr valign="top">
 *     <td> Time dependent Stokes flow driven by buoyancy
 *     </td>
 *     <td>step-31,
 *         step-32
 *     </td>
 *   </tr>
 *
 *   <tr valign="top">
 *     <td> The Euler equations of compressible gas dynamics
 *     </td>
 *     <td>
 *       step-33,
 *       step-67,
 *       step-69
 *     </td>
 *   </tr>
 *
 *   <tr valign="top">
 *     <td> The nonlinear Schr&ouml;dinger equation
 *     </td>
 *     <td>step-58</td>
 *   </tr>
 *
 *   <tr valign="top">
 *     <td> Time dependent neutron diffusion equation
 *     </td>
 *     <td>step-52</td>
 *   </tr>
 *
 *   <tr valign="top">
 *     <td> Time dependent fluid structure interaction problems
 *     </td>
 *     <td>step-70</td>
 *   </tr>
 * </table>
 */
