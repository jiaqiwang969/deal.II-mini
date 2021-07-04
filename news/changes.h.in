// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2018 by the deal.II authors
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
  Deprecated: The QTrapez class, poorly named because the proper English
  term is "trapezoidal quadrature rule", has been renamed to QTrapezoid,
  and the class with the old name has been deprecated.
  <br>
  (Wolfgang Bangerth, 2020/09/28)
 </li>

 <li>
  Replaced: Python wrapper for 'merge_triangulations' with a more generic equivalent.
  <br>
  (Alexander Grayver, 2020/09/01)
 </li>

 <li>
  Removed: CUDA 9, 10.0, and 10.1 are not supported anymore.
  <br>
  (Bruno Turcksin, 2020/08/05)
 </li>

 <li>
  Changed: The template arguments of the classes DoFAccessor and DoFCellAccessor have changed.
  The template argument DoFHandlerType has been replaced by dimension and
  space dimension.  
  <br>
  (Peter Munch, 2020/06/24)
 </li>

 <li>
  Deprecated: The functions MatrixFree::reinit(), which take 
  a vector of hp::DoFHandlers, have been deprecated. Users are asked
  to provide vectors of DoFhandlers, which may contain hp::DoFHandlers. This is
  possible now since hp::DoFHandler is deriving from DoFHandler.  
  <br>
  (Peter Munch, 2020/06/03)
 </li>

 <li>
  Removed: The deprecated class MGCoarseGridLACIteration has been removed.
  <br>
  (Daniel Arndt, 2020/06/12)
 </li>

 <li>
  Removed: The header file `deal.II/grid/tria_object.h` has been
  removed. It was only used for internal purposes.
  <br>
  (Wolfgang Bangerth, 2020/06/05)
 </li>

 <li>
  Changed: The binary representation of the Triangulation and DoFHandler classes 
  created by the function `save()` with release 9.2 cannot be read anymore with 
  `load()` due to major internal changes
  of these classes. This change also affects, i.a., the functions
  `GridIn::read_vtu()` and `GridOut::write_vtu()`.
  <br>
  (Peter Munch, 2020/06/03)
 </li>

 <li>
  Removed: The deprecated bindings to the legacy NETCDF C++ library have been
  removed.
  <br>
  (David Wells, 2020/05/27)
 </li>

 <li>
  Removed: The deprecated bindings to nanoflann have been removed.
  <br>
  (David Wells, 2020/05/27)
 </li>

 <li>
  Changed: The ThreadLocalStorage class has been reimplemented with C++14 STL
  primitives and does not depend on the TBB library any more. With that the
  obscure ThreadLocalStorage::get_implementation() function that exposed the
  underlying TBB container has been removed.
  <br>
  (Matthias Maier, 2020/05/23)
 </li>

 <li>
  Removed: The Threads::Task class had an `operator==` that allowed
  comparing objects of this type for equality. This operator has been
  removed. If you want to store objects of this kind in a collection
  that requires this kind of operator (say, `std::set`), then you
  probably can't do so any more in a reasonable way. However, this is
  exactly what the Threads::TaskGroup class is there for.
  <br>
  (Wolfgang Bangerth, 2020/05/26)
 </li>

 <li>
  Deprecated: The class hp::DoFHandler has been deprecated, since the DoFHandler
  has been extended with its functionalities.
  <br>
  (Peter Munch, 2020/05/23)
 </li>

 <li>
  Removed: The following preprocessor definitions have been removed from
  config.h.in: DEAL_II_NOEXCEPT, DEAL_II_USE_MT_POSIX,
  DEAL_II_USE_MT_POSIX_NO_BARRIERS
  <br>
  (Matthias Maier, 2020/05/23)
 </li>

 <li>
  Removed: The deprecated classes Threads::Mutex::ScopedLock,
  Threads::ConditionVariable, and deprecated functions
  Threads::Mutex::acquire(), Threads::Mutex::release(),
  Threads::n_existing_threads(), Threads::this_thread_id() have been removed.
  <br>
  (Matthias Maier, 2020/05/22)
 </li>

 <li>
  Removed: The <code>DEAL_II_WITH_CXX14</code> and
  <code>DEAL_II_WITH_CXX17</code> configuration options have been removed.
  The library will now be compiled with the default C++ standard enabled by
  the compiler. This is (as of May 2020) C++14 for all compilers. If you want
  to override that behavior, please set the C++ standard directly for example
  by configuring with <code>-DDEAL_II_CXX_FLAGS="-std=c++17"</code>, or by
  setting the environement variable <code>CXXFLAGS="-std=c++17"</code>.
  <br>
  (Matthias Maier, 2020/05/21)
 </li>

 <li>
  Updated: deal.II now requires a compiler with enabled C++14 support.
  <br>
  (Matthias Maier, 2020/05/21)
 </li>

 <li>
  Changed: The polynomial space template argument has been removed from
  FE_Poly and FE_PolyTensor.
  <br>
  (Graham Harper, Daniel Arndt, 2020/05/21)
 </li>

 <li>
  Removed: The deprecated class Threads::PosixThreadBarrier has been
  removed.
  <br>
  (Wolfgang Bangerth, 2020/04/21)
 </li>

</ol>


<!-- ----------- GENERAL IMPROVEMENTS ----------------- -->

<a name="general"></a>
<h3>General</h3>

<ol>

 <li>
  New: Tutorial example (step-68) showcasing parallel simulation of the advection
  of particles including load balancing.
  <br>
  (Bruno Blais, Toni El Geitani Nehme, Rene Gassm
  ller, Peter Munch, 2020/05/23)
 </li>

 <li>
  New: The step-19 tutorial program shows how to use particle methods.
  <br>
  (Wolfgang Bangerth, Rene Gassmoeller, Peter Munch, 2020/09/15)
 </li>

 <li>
  New: The multithreading framework in the library has been completely revamped:
  Intel TBB has been replaced by taskflow, which is bundled within the
  library. Multithreading support is now unconditionally enabled and we removed
  the CMake option DEAL_II_WITH_THREADS. The tasking framework can be controlled
  using MultithreadInfo and the DEAL_II_NUM_THREADS environment variable.
  <br>
  (Wolfgang Bangerth, Timo Heister, Matthias Maier, 2020/05/25)
 </li>

 <li>
  Changed: The internal data structures of DoFHandler have been modified to use 
  compressed row storage, enabling it to also handle hp::DoFHandler functionalities.
  Currently, the user can choose between the normal mode and the hp mode during
  calling the constructur. Please note that the multigrid functionalities are only 
  available during normal mode.
  <br>
  (Peter Munch, 2020/05/23)
 </li>

 <li>
  List rotated: The list of major changes is now empty.
  <br>
  (Matthias Maier, 2020/05/12)
 </li>

</ol>


<!-- ----------- SPECIFIC IMPROVEMENTS ----------------- -->

<a name="specific"></a>
<h3>Specific improvements</h3>

<ol>

 <li>
  Improved: MatrixFree now also works for hp in MPI-parallelized
  programs.
  <br>
  (Marc Fehling, Katharina Kormann, Martin Kronbichler, Peter Munch, 2020/10/20)
 </li>

 <li>
  New: step-9 uses the "streamline-upwind Petrov-Galerkin" method, but
  does not make any attempt at explaining what this method is or why it
  might be named like this. This has been rectified: The introduction
  now has a long section that explains the origin of the method and its name.
  <br>
  (Wolfgang Bangerth, 2020/10/10)
 </li>

 <li>
  New: Mapping::transform_points_real_to_unit_cell() can compute the operation
  of Mapping::transform_real_to_unit_cell() on many points simultaneously, which
  can be much faster for MappingQGeneric and derived classes that involve
  expensive operations to compute the support points of the mapping.
  <br>
  (Martin Kronbichler, 2020/10/07)
 </li>

 <li>
  New: SSP_THIRD_ORDER is added to the namespace TimeStepping to 
  implement the explicit third order Strong Stability Preserving (SSP) Runge-Kutta method, 
  which is also called the third order Total Variation Diminishing (TVD) Runge-Kutta method, see @cite gottlieb2001strong.
  <br>
  (Jiaqi Zhang, 2020/10/05)
 </li>

 <li>
  New: Helper function DoFCellAccessor::dominated_future_fe_on_children()
  to clean up code for hp-coarsening.
  <br>
  (Marc Fehling, 2020/10/05)
 </li>

 <li>
  New: GridTools::affine_cell_approximation() returns a matrix <i>A</i> and
  offset vector <i>b</i> that describe a least-squares fit of an affine
  approximation to a set of vertices of a cell.
  <br>
  (Martin Kronbichler, 2020/10/04)
 </li>

 <li>
  New: Helper functions CellAccessor::child_iterators() and
  DoFCellAccessor::child_iterators() which return iterators to children of
  a cell via `cell->child_iterators()`.
  <br>
  (Marc Fehling, 2020/10/03)
 </li>

 <li>
  New: CellId has a new constructor to create it from a std::string.
  <br>
  (Timo Heister, 2020/10/05)
 </li>

 <li>
  Improved: MappingQGeneric::transform_real_to_unit_cell() has been made much
  faster by directly working with the tensor product form of the mapping shape
  functions and avoiding many unnecessary memory allocations. The main cost is
  now MappingQGeneric::compute_mapping_support_points(), which can be made fast
  with MappingQCache, for example.
  <br>
  (Martin Kronbichler, 2020/09/30)
 </li>

 <li>
  New: The function BlockSparsityPattern::print_svg() outputs a block
  sparsity pattern in SVG format.
  <br>
  (Wolfgang Bangerth, 2020/09/25)
 </li>

 <li>
  Changed: step-29 no longer uses the `deallog` variable to generate
  output, but instead directly writes to `std::cout`.
  <br>
  (Wolfgang Bangerth, 2020/09/23)
 </li>

 <li>
  New: The classes FEEvaluation and FEFaceEvaluation with template parameter -1
  for the polynomial degree is now based on pre-compiled templated code for
  polynomial degrees between 1 and 6. This allows for fast execution of
  matrix-free evaluation for run-time polynomial degrees. The generated
  instantiations are controlled by
  `include/deal.II/matrix_free/evaluation_template_factory.templates.h` and can
  be pre-compiled for additional degrees in user code.
  <br>
  (Martin Kronbichler, Peter Munch, 2020/09/21)
 </li>

 <li>
  Fixed: Our cmake scripts forgot to install some of the files that are
  part of the code gallery. This is now fixed.
  <br>
  (Wolfgang Bangerth, 2020/09/17)
 </li>

 <li>
  Fixed: DoFTools::extract_dofs() return an IndexSet as result used to have a
  quadratic complexity in the number of extracted indices. This is now fixed.
  <br>
  (Martin Kronbichler, 2020/08/11)
 </li>

 <li>
  Added: method for returning list of all triangulation cells.
  <br>
  (Alexander Grayver, 2020/09/03)
 </li>

 <li>
  Added: python wrapper for GridTools::replicate_triangulation,
  more general version of the GridTools::merge_triangulations is
  implemented
  <br>
  (Alexander Grayver, 2020/09/01)
 </li>

 <li>
  New: The methods FEEvaluation::gather_evaluate(),
  FEEFaceEvaluation::gather_evaluate(), FEEvaluation::integrate_scatter() and
  FEfaceEvaluation::integrate_scatter() can now also accept block vectors.
  <br>
  (Peter Munch, Magdalena Schreter, Martin Kronbichler, 2020/08/31)
 </li>

 <li>
  New: A particle collection can now be copied into a new ParticleHandler object using
  the new ParticleHandler::copy_from function.
  <br>
  (Rene Gassmoeller, 2020/08/28)
 </li>

 <li>
  Added: DataOut now supports HDF5 file format with simplex meshes.
  <br>
  (Pasquale Claudio Africa, 2020/08/27)
 </li>

 <li>
  Updated: GridTools::transform() now works with simplex meshes.
  <br>
  (Pasquale Claudio Africa, 2020/08/26)
 </li>

 <li>
  Improved: The macro DEAL_II_PICKUP_TESTS can now also be run with on optional
  input parameter that can be used to manipulate the folder name of tests during
  ctest.
  <br>
  (Peter Munch, 2020/07/23)
 </li>

 <li>
  New: The test suite can now also be run with .json files.
  <br>
  (Peter Munch, 2020/08/17)
 </li>

 <li>
  Improved: The definitions of the templated functions of the HDF5 interface are
  now in hdf5.h, therefore it is possible to use additional template combinations.
  The instantiations are no longer necessary, therefore they have been removed.
  <br>
  (Daniel Garcia-Sanchez, 2020/08/14)
 </li>

 <li>
  New: MGTools::make_sparsity_pattern() can now take an optional
  AffineConstraints argument to add the effect of, e.g., periodic boundary
  conditions.
  <br>
  (Martin Kronbichler, 2020/08/11)
 </li>

 <li>
  Fixed: The DataPostprocessorTensor class erroneously announced that
  the components of its output are a bunch of scalars when, of course,
  the whole point of the class was to output things as one tensor
  field. This is now fixed.
  <br>
  (Wolfgang Bangerth, 2020/08/06)
 </li>

 <li>
  New: GridIn::read_vtk() now supports triangular and tetrahedral meshes.
  <br>
  (Peter Munch, 2020/07/23)
 </li>

 <li>
  New: GridIn::read_msh() now supports triangular and tetrahedral meshes.
  <br>
  (Daniel Paukner, 2020/07/20)
 </li>

 <li>
  New: The method hp::FEFaceValues::reinit() can now also accept face iterators.
  <br>
  (Peter Munch, 2020/07/16)
 </li>

 <li>
  New: The class hp::MappingCollection has a new constructor. This constructor creates 
  a MappingCollection from one or more mapping objects passed to the constructor.
  <br>
  (Peter Munch, 2020/07/15)
 </li>

 <li>
  Fixed: MeshWorker::mesh_loop() did not work on 1d with refined grids. This is now fixed.
  <br>
  (Luca Heltai, 2020/07/08)
 </li>

 <li>
  Added the functions CUDAWrappers::MatrixFree::get_vector_partitioner() and
  CUDAWrappers::MatrixFree::get_dof_handler()
  <br>
  (Bruno Turcksin, 2020/07/06)
 </li>

 <li>
  New: The new class Simplex::ScalarPolynomial provides polynomials defined on 
  simplices.
  <br>
  (Peter Munch, 2020/07/02)
 </li>

 <li>
  Fixed: In parallel hp-adaptive applications,
  DoFHandler::distribute_dofs() no longer fails to enumerate degrees of
  freedom on ghost interfaces if continuous finite elements do not
  dominate each other.
  <br>
  (Marc Fehling, 2020/07/03)
 </li>

 <li>
  New: A new quadrature rule for simplex geometric entities has been added.
  <br>
  (Peter Munch, 2020/07/02)
 </li>

 <li>
  New: Geometric objects of a Triangulation are assigned a ReferenceCell::Type. The value
  can be queried via TriaAccessor::reference_cell_type().
  <br>
  (Peter Munch, 2020/06/30)
 </li>

 <li>
  New: The class ArrayView can now also be constructed from std::array.
  <br>
  (Peter Munch, 2020/06/29)
 </li>

 <li>
  New: BoundingBox::real_to_unit() and BoundingBox::unit_to_real() allow one to
  apply both the direct and the inverse transformation that are needed to map the
  unit bounding box to the current box, and viceversa.
  <br>
  (Luca Heltai, 2020/06/29)
 </li>

 <li>
  New: The member function DiscreteTime::set_next_step_size() is added.
  <br>
  (Reza Rastak, 2020/06/27)
 </li>

 <li>
  New: There is now a constructor for class Tensor that takes
  an initializer from an object of type ArrayView.
  <br>
  (Wolfgang Bangerth, 2020/06/27)
 </li>

 <li>
  Fixed: The class parallel::distributed::SolutionTransfer can now
  also handle FE_Nothing.
  <br>
  (Dominic Soldner, Peter Munch, 2020/06/24)
 </li>

 <li>
  Fixed: FEInterfaceValues now works also for codim one and two. Instantiated
  also DoFTools::make_flux_sparsity_pattern() for codim one and two.
  <br>
  (Luca Heltai, 2020/06/24)
 </li>

 <li>
  New: The class TriaAccessor provides now the capability to query
  the number of vertices, lines, and faces (with n_vertices(), 
  n_lines(), n_faces(), vertex_indices(), 
  line_indices(), face_indices()). The new methods can be used as an 
  alternative to the methods in GeometryInfo.
  <br>
  (Peter Munch, 2020/06/23)
 </li>

 <li>
  New: Added FEInterfaceValues to MeshWorker::ScratchData.
  <br>
  (Luca Heltai, 2020/06/23)
 </li>

 <li>
  Fixed: The ParticleHandler::insert_particles() function forgot to
  associate particles with the common property pool. Consequently,
  setting properties on particles added to a ParticleHandler this way
  led to an error.
  <br>
  (Andrew Davis, Wolfgang Bangerth, 2020/06/23)
 </li>

 <li>
  New: Particles::Particle and Particles::ParticleAccessor can now be used as
  indexable in boost::rtree objects. 
  <br> (Luca Heltai, 2020/06/15)
 </li>

 <li>
  New: The function Particles::ParticleHandler::add_global_particles() now takes 
  another optional argument, that allows one to set ids arbitrarily. Moreover,
  now the numbering of the ids is correct also if we call the method more than
  one time. Newly added particles, if ids are not specified, now correctly get
  the first available ids. 
  Added a new version of Particles::ParticleHandler::add_global_particles() that
  takes a vector of Particles::Particle objects instead of just their positions.
  This can be used in conjunction with the signal
  Particles::ParticleHandler::Signals::particle_lost() to reinsert
  Particles::Particle objects that went out of the locally owned and ghost cells.
  <br> (Luca Heltai, 2020/06/11)
 </li>

 <li>
  Fixed: The AffineConstraints class had a bug where, if deal.II was
  compiled without threading, the class used a global variable. If a
  user program used threads nonetheless, this global variable led to
  race conditions. This is now fixed.
  <br>
  (Wolfgang Bangerth, 2020/06/11)
 </li>

 <li>
  Fixed: Fix a bug where CUDAWrappers::MatrixFree::cell_loop() would set the
  destination vector to zero if the partitioner of the MatrixFree object was
  different from the partitioner of the source or destination vector.
  <br>
  (Bruno Turcksin, 2020/06/10)
 </li>

 <li>
  New: Applying user constraints 
  before prolongation in MGTransferPrebuilt.
  <br>
  (Julian Roth and Conrad Clevenger, 2020/06/05)
 </li>

 <li>
  New: FEEvaluation::evaluate(), FEEvaluation::integrate(),
  FEFaceEvaluation::evaluate() and FEFaceEvaluation::integrate() now take an
  argument of type EvaluationFlags::EvaluationFlags to determine which of the
  values, gradients or hessians should be evaluated to quadrature points or
  integrated, allowing much more expressive programs than the previous list of
  bools. The evaluation flags can be combined with `operator|`, similarly to
  UpdateFlags for FEValues.
  <br>
  (Timo Heister, 2020/06/05)
 </li>

 <li>
  Changed: The vertices in CellData are now stored in form of a std::vector
  instead of C-style array.
  <br>
  (Peter Munch, 2020/05/31)
 </li>

 <li>
  Improved: The efficiency of the assembly of step-62 has been improved and now
  it is 7 times faster.
  <br>
  (Daniel Garcia-Sanchez, 2020/05/31)
 </li>

 <li>
  Fixed: Fix a bug where only one CUDAWrappers::MatrixFree object was valid at a
  given time. There is now a variable CUDAWrappers::mf_n_concurrent_objects in
  base/cuda_size.h that controls the maximum number of concurrent objects. The
  default value is five.
  <br>
  (Bruno Turcksin, 2020/05/29)
 </li>

 <li>
  New: Add multigrid transfer operators for distributed polynomial and 
  global coarsening.
  <br>
  (Peter Munch, Laura Prieto Saavedra, 2020/05/29)
 </li>

 <li>
  Improved: step-28 now uses tasks instead of threads.
  <br>
  (David Wells, 2020/05/28)
 </li>

 <li>
  New: The class Particles::DataOut can now output particle properties as
  scalars, vectors, or tensors, depending on the arguments handed over to the
  Particles::DataOut::build_patches() function.
  <br>
  (Rene Gassmoeller, 2020/05/27)
 </li>

 <li>
  New: When executing a task on a separate thread, if that task ends
  with throwing an exception instead of returing a value, then this
  exception will be obtained when you wait for the task using
  Threads::Task::join() or Threads::Task::return_value().
  <br>
  (Wolfgang Bangerth, 2020/05/27)
 </li>

 <li>
  New: GridTools::Cache::get_locally_owned_cell_bounding_boxes_rtree() extracts
  a tree of bounding boxes covering the locally owned cells of a triangulation.
  This can be used in conjunction with GridTools::Cache::get_covering_rtree() to
  make approximate geometrical queries on who owns what spatial region.
  <br>
  (Luca Heltai, 2020/05/26)
 </li>

 <li>
  New: pack_rtree_of_indices() and IndexableGetterFromIndices allow to construct an
  RTree object that stores indices to existing containers of indexable objects.
  <br>
  (Luca Heltai, 2020/05/24)
 </li>

 <li>
  New: The class ParticleHandler now provides a signal 'signals.particle_lost'
  that is triggered whenever a particles can not be associated with a cell while
  calling its function sort_particles_into_subdomains_and_cells().
  <br>
  (Rene Gassmoeller, 2020/05/25)
 </li>

 <li>
  Bugfix: hp::Refinement::choose_p_over_h() now works in parallel.
  <br>
  (Marc Fehling, 2020/05/22)
 </li>

 <li>
  New: There is now a second overload for
  Particles::ParticleAccessor::set_properties() that takes an ArrayView
  as argument.
  <br>
  (Wolfgang Bangerth, 2020/05/22)
 </li>

 <li>
  Removed: All headers under <code>base/std_cxx11/</code> have been removed.
  <br>
  (David Wells, 2020/05/21)
 </li>

 <li>
  Changed: In many other classes that are templated on `dim` and
  `spacedim`, the second template argument `spacedim` had a default
  value equal to `dim`. Particles::DataOut did not, but now does.
  <br>
  (Wolfgang Bangerth, 2020/05/21)
 </li>

 <li>
  New: A new BoundingBoxDataOut class is available, to output a collection
  of objects that can be converted to BoundingBox objects in graphical format.
  <br>
  (Luca Heltai, 2020/05/19)
 </li>

 <li>
  List rotated: The list of minor changes is now empty.
  <br>
  (Matthias Maier, 2020/05/12)
 </li>

</ol>

*/
