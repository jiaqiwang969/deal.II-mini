// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2021 by the deal.II authors
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
@page recent_changes Changes since the last release

<p>
This is the list of changes made since the last release of deal.II.
All entries are signed with the names of the authors.
</p>
<!-- ----------- INCOMPATIBILITIES ----------------- -->

<a name="incompatible"></a>
<h3 style="color:red">Incompatibilities</h3>

<p style="color:red">
Following are a few modifications to the library that unfortunately
are incompatible with previous versions of the library, but which we
deem necessary for the future maintainability of the
library. Unfortunately, some of these changes will require
modifications to application programs. We apologize for the
inconvenience this causes.
</p>
<ol>

 <li>
  Removed: The deprecated specializations of VectorTools::integrate_difference()
  and VectorTools::project_boundary_values_curl_conforming() have been removed.
  <br>
  (Daniel Arndt, 2021/06/04)
 </li>

 <li>
  Removed: The deprecated overloads of 
  FETools::lexicographic_to_hierarchic_numbering() and
  FETools::hierarchic_to_lexicographic_numbering() have been removed.
  <br>
  (Daniel Arndt, 2021/06/03)
 </li>

 <li>
  Removed: The deprecated member variable VectorizedArray::n_array_elements has
  been removed.
  <br>
  (Daniel Arndt, 2021/06/03)
 </li>

 <li>
  Removed: The deprecated member function CellAccessor::active() has been removed.
  <br>
  (Daniel Arndt, 2021/06/03)
 </li>

 <li>
  Removed: The deprecated member functions DataOut::first_cell() and DataOut::next_cell()
  have been removed.
  <br>
  (Daniel Arndt, 2021/06/02)
 </li>

 <li>
  Removed: The deprecated build_matrices() member functions
  in various multigrid transfer classes have been removed.
  <br>
  (Daniel Arndt, 2021/06/01)
 </li>

 <li>
  Removed: The deprecated versions of DoFTools::count_dofs_per_component()
  and DoFTools::count_dofs_per_block() have been removed.
  <br>
  (Daniel Arndt, 2021/06/01)
 </li>

 <li>
  Removed: The deprecated overloads for DoFTools::extract_hanging_node_dofs(),
  and DoFTools::extract_dofs() have been removed.
  <br>
  (Daniel Arndt, 2021/05/28)
 </li>

 <li>
  Removed: Deprecated parallel::CellWeights member functions have been removed.
  <br>
  (Daniel Arndt, 2021/05/27)
 </li>

 <li>
  Removed: The overloads in CUDAWrappers::FEEvaluation that take a local dof index
  or a quadrature point as arguement have been removed.
  Use the ones that don't use these arguments in their interface instead.
  <br>
  (Daniel Arndt, 2021/05/26)
 </li>

 <li>
  Removed: The deprecated parallel::Triangulation class has been removed.
  Use parallel::TriangulationBase instead.
  <br>
  (Daniel Arndt, 2021/05/26)
 </li>

 <li>
  Removed: The deprecated member field 
  MatrixFree::AdditionalData::level_mg_handler.
  <br>
  (Peter Munch, 2021/05/25)
 </li>

 <li>
  Removed: The deprecated
  parallel::distributed::CellDataTransfer::CoarseningStrategies struct has been removed.
  Use AdaptationStrategies::Coarsening instead.
  <br>
  (Daniel Arndt, 2021/05/25)
 </li>

 <li>
  Removed: The deprecated member functions 
  DoFHandler::locally_owned_dofs_per_processor(), 
  DoFHandler::n_locally_owned_dofs_per_processor(), and 
  DoFHandler::n_locally_owned_mg_dofs_per_processor() have been removed.
  <br>
  (Daniel Arndt, 2021/05/24)
 </li>

 <li>
  New: The incompatibilities list has been rotated.
  <br>
  (Matthias Maier, 2021/05/22)
 </li>

 <li>
  Deprecated: The template arguments of the following classes have been
  changed to avoid the legacy argument `DoFHandlerType`:
  <ul>
    <li> `SolutionTransfer<dim, VectorType, DoFHandlerType> -> SolutionTransfer<dim, VectorType, spacedim>`
    <li> `parallel::distributed::SolutionTransfer<dim, VectorType, DoFHandlerType> -> parallel::distributed::SolutionTransfer<dim, VectorType, spacedim>`
    <li> `Functions::FEFieldFunction<dim, DoFHandlerType, VectorType> -> Functions::FEFieldFunction<dim, VectorType, spacedim>`
    <li> `DataOut<dim, DoFHandlerType> -> DataOut<dim, spacedim>`
    <li> `DataOut_DoFData<DoFHandlerType, patch_dim, patch_space_dim> -> DataOut_DoFData<dim, patch_dim, spacedim, patch_spacedim>`
    <li> `DataOutFaces<dim, DoFHandlerType> -> DataOutFaces<dim, spacedim>`
    <li> `DataOutRotation<dim, DoFHandlerType> -> DataOutRotation<dim, spacedim>`
  </ul>
  Please change your code accordingly.
  <br>
  If for some reason, you need a code that is compatible with deal.II
  9.3 and the subsequent release, a Legacy namespace has been introduced
  with aliases of these classes with their original interface. You can
  make the following substitutions to your code for each of the affected
  classes:
  <ul>
    <li>X &rarr; Legacy::X
  </ul>
  To perform this substitution automatically, you may use a *search and
  replace* script like the following made for the *bash* shell:
  @code{.sh}
  classes=(SolutionTransfer parallel::distributed::SolutionTransfer Functions::FEFieldFunction DataOut DataOut_DoFData DataOutFaces DataOutRotation)
  for c in \${classes[@]}; do
    find /path/to/your/code -type f -exec sed -i -E "/(\w\${c}|\${c}[^<]|Particles::\${c}|distributed::\${c}|^\s*(\/\/|\*))/! s/\${c}/Legacy::\${c}/g" {} \;
  done
  @endcode
  (Marc Fehling, 2020/11/21)
 </li>

</ol>


<!-- ----------- GENERAL IMPROVEMENTS ----------------- -->

<a name="general"></a>
<h3>General</h3>
<ol>

 <li>
  New: The major changes list has been rotated.
  <br>
  (Matthias Maier, 2021/05/22)
 </li>

</ol>


<!-- ----------- SPECIFIC IMPROVEMENTS ----------------- -->

<a name="specific"></a>
<h3>Specific improvements</h3>
<ol>

 <li>
  New: The new function Utilities::MPI::reduce() allows to reduce 
  arbitrary types with a user-specified binary operation.
  <br>
  (Peter Munch, 2021/05/27)
 </li>

 <li>
  New: Fixed a bug in clear user data for standard Triangulation.
  <br>
  (Nicola Giuliani, 2021/05/25)
 </li>

 <li>
  New: The minor changes list has been rotated.
  <br>
  (Matthias Maier, 2021/05/22)
 </li>

</ol>

*/
