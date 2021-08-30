

/* --------------------------------------------------------------------- 
 * 
 * Copyright (C) 1999 - 2021 by the deal.II authors 
 * 
 * This file is part of the deal.II library. 
 * 
 * The deal.II library is free software; you can use it, redistribute 
 * it, and/or modify it under the terms of the GNU Lesser General 
 * Public License as published by the Free Software Foundation; either 
 * version 2.1 of the License, or (at your option) any later version. 
 * The full text of the license can be found in the file LICENSE.md at 
 * the top level directory of deal.II. 
 * 
 * --------------------------------------------------------------------- 
 * 
 * Author: Wolfgang Bangerth, University of Heidelberg, 1999 
 */ 




#include <deal.II/grid/tria.h> 
#include <deal.II/grid/grid_generator.h> 


#include <deal.II/dofs/dof_handler.h> 



#include <deal.II/fe/fe_q.h> 


#include <deal.II/dofs/dof_tools.h> 


#include <deal.II/lac/sparse_matrix.h> 


#include <deal.II/lac/dynamic_sparsity_pattern.h> 


#include <deal.II/dofs/dof_renumbering.h> 


#include <fstream> 


using namespace dealii; 


void make_grid(Triangulation<2> &triangulation) 
{ 
  const Point<2> center(1, 0); 
  const double   inner_radius = 0.5, outer_radius = 1.0; 
  GridGenerator::hyper_shell( 
    triangulation, center, inner_radius, outer_radius, 5); 

  for (unsigned int step = 0; step < 3; ++step) 
    { 
      for (auto &cell : triangulation.active_cell_iterators()) 
        for (const auto v : cell->vertex_indices()) 
          { 
            const double distance_from_center = 
              center.distance(cell->vertex(v)); 

            if (std::fabs(distance_from_center - inner_radius) <= 
                1e-6 * inner_radius) 
              { 
                cell->set_refine_flag(); 
                break; 
              } 
          } 

      triangulation.execute_coarsening_and_refinement(); 
    } 
} 




void distribute_dofs(DoFHandler<2> &dof_handler) 
{ 
  const FE_Q<2> finite_element(1); 
  dof_handler.distribute_dofs(finite_element); 




  DynamicSparsityPattern dynamic_sparsity_pattern(dof_handler.n_dofs(), 
                                                  dof_handler.n_dofs()); 


  DoFTools::make_sparsity_pattern(dof_handler, dynamic_sparsity_pattern); 


  SparsityPattern sparsity_pattern; 
  sparsity_pattern.copy_from(dynamic_sparsity_pattern); 


  std::ofstream out("sparsity_pattern1.svg"); 
  sparsity_pattern.print_svg(out); 



} 





void renumber_dofs(DoFHandler<2> &dof_handler) 
{ 
  DoFRenumbering::Cuthill_McKee(dof_handler); 

  DynamicSparsityPattern dynamic_sparsity_pattern(dof_handler.n_dofs(), 
                                                  dof_handler.n_dofs()); 
  DoFTools::make_sparsity_pattern(dof_handler, dynamic_sparsity_pattern); 

  SparsityPattern sparsity_pattern; 
  sparsity_pattern.copy_from(dynamic_sparsity_pattern); 

  std::ofstream out("sparsity_pattern2.svg"); 
  sparsity_pattern.print_svg(out); 
} 





int main() 
{ 
  Triangulation<2> triangulation; 
  make_grid(triangulation); 

  DoFHandler<2> dof_handler(triangulation); 

  distribute_dofs(dof_handler); 
  renumber_dofs(dof_handler); 
} 


