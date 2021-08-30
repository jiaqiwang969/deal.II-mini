

/* --------------------------------------------------------------------- 
 * 
 * Copyright (C) 2001 - 2021 by the deal.II authors 
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
 * Authors: Wolfgang Bangerth, Ralf Hartmann, University of Heidelberg, 2001 
 */ 




#include <deal.II/base/quadrature_lib.h> 
#include <deal.II/base/convergence_table.h> 
#include <deal.II/grid/grid_generator.h> 
#include <deal.II/grid/manifold_lib.h> 
#include <deal.II/grid/tria.h> 
#include <deal.II/grid/grid_out.h> 
#include <deal.II/dofs/dof_handler.h> 
#include <deal.II/fe/fe_values.h> 


#include <deal.II/fe/fe_nothing.h> 


#include <deal.II/fe/mapping_q.h> 


#include <iostream> 
#include <fstream> 
#include <cmath> 


namespace Step10 
{ 
  using namespace dealii; 


  const long double pi = 3.141592653589793238462643L; 



  template <int dim> 
  void gnuplot_output() 
  { 
    std::cout << "Output of grids into gnuplot files:" << std::endl 
              << "===================================" << std::endl; 

    Triangulation<dim> triangulation; 
    GridGenerator::hyper_ball(triangulation); 


    for (unsigned int refinement = 0; refinement < 2; ++refinement) 
      { 
        std::cout << "Refinement level: " << refinement << std::endl; 

        std::string filename_base = "ball_" + std::to_string(refinement); 

        for (unsigned int degree = 1; degree < 4; ++degree) 
          { 
            std::cout << "Degree = " << degree << std::endl; 


            const MappingQ<dim> mapping(degree); 



            GridOut               grid_out; 
            GridOutFlags::Gnuplot gnuplot_flags(false, 60); 
            grid_out.set_flags(gnuplot_flags); 


            std::string filename = 
              filename_base + "_mapping_q_" + std::to_string(degree) + ".dat"; 
            std::ofstream gnuplot_file(filename); 


            grid_out.write_gnuplot(triangulation, gnuplot_file, &mapping); 
          } 
        std::cout << std::endl; 


        triangulation.refine_global(); 
      } 
  } 


  template <int dim> 
  void compute_pi_by_area() 
  { 
    std::cout << "Computation of Pi by the area:" << std::endl 
              << "==============================" << std::endl; 


    const QGauss<dim> quadrature(4); 


    for (unsigned int degree = 1; degree < 5; ++degree) 
      { 
        std::cout << "Degree = " << degree << std::endl; 


        Triangulation<dim> triangulation; 
        GridGenerator::hyper_ball(triangulation); 

        const MappingQ<dim> mapping(degree); 


        const FE_Nothing<dim> fe; 


        DoFHandler<dim> dof_handler(triangulation); 



        FEValues<dim> fe_values(mapping, fe, quadrature, update_JxW_values); 


        ConvergenceTable table; 


        for (unsigned int refinement = 0; refinement < 6; 
             ++refinement, triangulation.refine_global(1)) 
          { 


            table.add_value("cells", triangulation.n_active_cells()); 


            dof_handler.distribute_dofs(fe); 


            long double area = 0; 


            for (const auto &cell : dof_handler.active_cell_iterators()) 
              { 
                fe_values.reinit(cell); 
                for (unsigned int i = 0; i < fe_values.n_quadrature_points; ++i) 
                  area += static_cast<long double>(fe_values.JxW(i)); 
              } 


            table.add_value("eval.pi", static_cast<double>(area)); 
            table.add_value("error", static_cast<double>(std::fabs(area - pi))); 
          } 


        table.omit_column_from_convergence_rate_evaluation("cells"); 
        table.omit_column_from_convergence_rate_evaluation("eval.pi"); 
        table.evaluate_all_convergence_rates( 
                                    ConvergenceTable::reduction_rate_log2); 


        table.set_precision("eval.pi", 16); 
        table.set_scientific("error", true); 

        table.write_text(std::cout); 

        std::cout << std::endl; 
      } 
  } 


  template <int dim> 
  void compute_pi_by_perimeter() 
  { 
    std::cout << "Computation of Pi by the perimeter:" << std::endl 
              << "===================================" << std::endl; 


    const QGauss<dim - 1> quadrature(4); 


    for (unsigned int degree = 1; degree < 5; ++degree) 
      { 
        std::cout << "Degree = " << degree << std::endl; 
        Triangulation<dim> triangulation; 
        GridGenerator::hyper_ball(triangulation); 

        const MappingQ<dim>   mapping(degree); 
        const FE_Nothing<dim> fe; 

        DoFHandler<dim> dof_handler(triangulation); 


        FEFaceValues<dim> fe_face_values(mapping, 
                                         fe, 
                                         quadrature, 
                                         update_JxW_values); 
        ConvergenceTable  table; 

        for (unsigned int refinement = 0; refinement < 6; 
             ++refinement, triangulation.refine_global(1)) 
          { 
            table.add_value("cells", triangulation.n_active_cells()); 

            dof_handler.distribute_dofs(fe); 


            long double perimeter = 0; 
            for (const auto &cell : dof_handler.active_cell_iterators()) 
              for (const auto &face : cell->face_iterators()) 
                if (face->at_boundary()) 
                  { 


                    fe_face_values.reinit(cell, face); 
                    for (unsigned int i = 0; 
                         i < fe_face_values.n_quadrature_points; 
                         ++i) 
                      perimeter += 
                        static_cast<long double>(fe_face_values.JxW(i)); 
                  } 


            table.add_value("eval.pi", static_cast<double>(perimeter / 2.0L)); 
            table.add_value( 
              "error", static_cast<double>(std::fabs(perimeter / 2.0L - pi))); 
          } 


        table.omit_column_from_convergence_rate_evaluation("cells"); 
        table.omit_column_from_convergence_rate_evaluation("eval.pi"); 
        table.evaluate_all_convergence_rates( 
          ConvergenceTable::reduction_rate_log2); 

        table.set_precision("eval.pi", 16); 
        table.set_scientific("error", true); 

        table.write_text(std::cout); 

        std::cout << std::endl; 
      } 
  } 
} // namespace Step10 


int main() 
{ 
  try 
    { 
      std::cout.precision(16); 

      const unsigned int dim = 2; 

      Step10::gnuplot_output<dim>(); 

      Step10::compute_pi_by_area<dim>(); 
      Step10::compute_pi_by_perimeter<dim>(); 
    } 
  catch (std::exception &exc) 
    { 
      std::cerr << std::endl 
                << std::endl 
                << "----------------------------------------------------" 
                << std::endl; 
      std::cerr << "Exception on processing: " << std::endl 
                << exc.what() << std::endl 
                << "Aborting!" << std::endl 
                << "----------------------------------------------------" 
                << std::endl; 

      return 1; 
    } 
  catch (...) 
    { 
      std::cerr << std::endl 
                << std::endl 
                << "----------------------------------------------------" 
                << std::endl; 
      std::cerr << "Unknown exception!" << std::endl 
                << "Aborting!" << std::endl 
                << "----------------------------------------------------" 
                << std::endl; 
      return 1; 
    } 

  return 0; 
} 


