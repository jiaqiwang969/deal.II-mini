

/* --------------------------------------------------------------------- 
 * 
 * Copyright (C) 2009 - 2021 by the deal.II authors 
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
 *  Authors: Andrea Mola, Luca Heltai, 2014 
 */ 




#include <deal.II/grid/tria.h> 
#include <deal.II/grid/grid_generator.h> 
#include <deal.II/grid/grid_in.h> 
#include <deal.II/grid/grid_out.h> 
#include <deal.II/numerics/data_out.h> 
#include <deal.II/numerics/vector_tools.h> 


#include <deal.II/opencascade/manifold_lib.h> 
#include <deal.II/opencascade/utilities.h> 


#include <cmath> 
#include <iostream> 
#include <fstream> 
#include <string> 


namespace Step54 
{ 
  using namespace dealii; 




  class TriangulationOnCAD 
  { 
  public: 
    enum ProjectionType 
    { 
      NormalProjection       = 0, 
      DirectionalProjection  = 1, 
      NormalToMeshProjection = 2 
    }; 

    TriangulationOnCAD( 
      const std::string &  initial_mesh_filename, 
      const std::string &  cad_file_name, 
      const std::string &  output_filename, 
      const ProjectionType surface_projection_kind = NormalProjection); 

    void run(); 

  private: 
    void read_domain(); 

    void refine_mesh(); 

    void output_results(const unsigned int cycle); 

    Triangulation<2, 3> tria; 

    const std::string initial_mesh_filename; 
    const std::string cad_file_name; 
    const std::string output_filename; 

    const ProjectionType surface_projection_kind; 
  }; 


  TriangulationOnCAD::TriangulationOnCAD( 
    const std::string &  initial_mesh_filename, 
    const std::string &  cad_file_name, 
    const std::string &  output_filename, 
    const ProjectionType surface_projection_kind) 
    : initial_mesh_filename(initial_mesh_filename) 
    , cad_file_name(cad_file_name) 
    , output_filename(output_filename) 
    , surface_projection_kind(surface_projection_kind) 
  {} 





  void TriangulationOnCAD::read_domain() 
  { 
    TopoDS_Shape bow_surface = OpenCASCADE::read_IGES(cad_file_name, 1e-3); 



    const double tolerance = OpenCASCADE::get_shape_tolerance(bow_surface) * 5; 


    std::vector<TopoDS_Compound>  compounds; 
    std::vector<TopoDS_CompSolid> compsolids; 
    std::vector<TopoDS_Solid>     solids; 
    std::vector<TopoDS_Shell>     shells; 
    std::vector<TopoDS_Wire>      wires; 

    OpenCASCADE::extract_compound_shapes( 
      bow_surface, compounds, compsolids, solids, shells, wires); 


    std::ifstream in; 

    in.open(initial_mesh_filename); 

    GridIn<2, 3> gi; 
    gi.attach_triangulation(tria); 
    gi.read_vtk(in); 


    output_results(0); 


    Triangulation<2, 3>::active_cell_iterator cell = tria.begin_active(); 
    cell->set_manifold_id(1); 

    for (const auto &face : cell->face_iterators()) 
      face->set_manifold_id(2); 




    Assert( 
      wires.size() > 0, 
      ExcMessage( 
        "I could not find any wire in the CAD file you gave me. Bailing out.")); 

    OpenCASCADE::ArclengthProjectionLineManifold<2, 3> line_projector( 
      wires[0], tolerance); 

    tria.set_manifold(2, line_projector); 


    switch (surface_projection_kind) 
      { 
        case NormalProjection: 
          { 
            OpenCASCADE::NormalProjectionManifold<2, 3> normal_projector( 
              bow_surface, tolerance); 
            tria.set_manifold(1, normal_projector); 

            break; 
          } 

        case DirectionalProjection: 
          { 
            OpenCASCADE::DirectionalProjectionManifold<2, 3> 
              directional_projector(bow_surface, 
                                    Point<3>(0.0, 1.0, 0.0), 
                                    tolerance); 
            tria.set_manifold(1, directional_projector); 

            break; 
          } 


        case NormalToMeshProjection: 
          { 
            OpenCASCADE::NormalToMeshProjectionManifold<2, 3> 
              normal_to_mesh_projector(bow_surface, tolerance); 
            tria.set_manifold(1, normal_to_mesh_projector); 

            break; 
          } 


        default: 
          AssertThrow(false, ExcInternalError()); 
      } 
  } 



  void TriangulationOnCAD::refine_mesh() 
  { 
    tria.refine_global(1); 
  } 



  void TriangulationOnCAD::output_results(const unsigned int cycle) 
  { 
    const std::string filename = 
      (output_filename + "_" + Utilities::int_to_string(cycle) + ".vtk"); 
    std::ofstream logfile(filename); 
    GridOut       grid_out; 
    grid_out.write_vtk(tria, logfile); 
  } 


  void TriangulationOnCAD::run() 
  { 
    read_domain(); 

    const unsigned int n_cycles = 5; 
    for (unsigned int cycle = 0; cycle < n_cycles; ++cycle) 
      { 
        refine_mesh(); 
        output_results(cycle + 1); 
      } 
  } 
} // namespace Step54 


int main() 
{ 
  try 
    { 
      using namespace Step54; 

      const std::string in_mesh_filename = "input/initial_mesh_3d.vtk"; 
      const std::string cad_file_name    = "input/DTMB-5415_bulbous_bow.iges"; 

      std::cout << "----------------------------------------------------------" 
                << std::endl; 
      std::cout << "Testing projection in direction normal to CAD surface" 
                << std::endl; 
      std::cout << "----------------------------------------------------------" 
                << std::endl; 
      std::string        out_mesh_filename = ("3d_mesh_normal_projection"); 
      TriangulationOnCAD tria_on_cad_norm(in_mesh_filename, 
                                          cad_file_name, 
                                          out_mesh_filename, 
                                          TriangulationOnCAD::NormalProjection); 
      tria_on_cad_norm.run(); 
      std::cout << "----------------------------------------------------------" 
                << std::endl;
      std::cout << std::endl; 
      std::cout << std::endl; 

      std::cout << "----------------------------------------------------------" 
                << std::endl; 
      std::cout << "Testing projection in y-axis direction" << std::endl; 
      std::cout << "----------------------------------------------------------" 
                << std::endl; 
      out_mesh_filename = ("3d_mesh_directional_projection"); 
      TriangulationOnCAD tria_on_cad_dir( 
        in_mesh_filename, 
        cad_file_name, 
        out_mesh_filename, 
        TriangulationOnCAD::DirectionalProjection);  
      tria_on_cad_dir.run();  
      std::cout << "----------------------------------------------------------"  
                << std::endl; 
      std::cout << std::endl; 
      std::cout << std::endl; 

      std::cout << "----------------------------------------------------------" 
                << std::endl; 
      std::cout << "Testing projection in direction normal to mesh elements" 
                << std::endl; 
      std::cout << "----------------------------------------------------------" 
                << std::endl; 
      out_mesh_filename = ("3d_mesh_normal_to_mesh_projection"); 
      TriangulationOnCAD tria_on_cad_norm_to_mesh( 
        in_mesh_filename, 
        cad_file_name, 
        out_mesh_filename, 
        TriangulationOnCAD::NormalToMeshProjection); 
      tria_on_cad_norm_to_mesh.run(); 
      std::cout << "----------------------------------------------------------" 
                << std::endl; 
      std::cout << std::endl; 
      std::cout << std::endl; 
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

