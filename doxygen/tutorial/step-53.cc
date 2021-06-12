

/* --------------------------------------------------------------------- 
 * 
 * Copyright (C) 2014 - 2020 by the deal.II authors 
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
 * Authors: Wolfgang Bangerth, Texas A&M University, 2014 
 *          Luca Heltai, SISSA, 2014 
 *          D. Sarah Stamps, MIT, 2014 
 */ 




#include <deal.II/grid/tria.h> 
#include <deal.II/grid/grid_generator.h> 
#include <deal.II/grid/grid_out.h> 
#include <deal.II/grid/manifold.h> 
#include <deal.II/grid/grid_tools.h> 


#include <deal.II/base/function_lib.h> 

#include <boost/iostreams/filtering_stream.hpp> 
#include <boost/iostreams/filter/gzip.hpp> 
#include <boost/iostreams/device/file.hpp> 

#include <fstream> 
#include <iostream> 
#include <memory> 


namespace Step53 
{ 
  using namespace dealii; 



  class AfricaTopography 
  { 
  public: 
    AfricaTopography(); 

    double value(const double lon, const double lat) const; 

  private: 
    const Functions::InterpolatedUniformGridData<2> topography_data; 

    static std::vector<double> get_data(); 
  }; 


  AfricaTopography::AfricaTopography() 
    : topography_data({{std::make_pair(-6.983333, 11.966667), 
                        std::make_pair(25, 35.95)}}, 
                      {{379, 219}}, 
                      Table<2, double>(380, 220, get_data().begin())) 
  {} 

  double AfricaTopography::value(const double lon, const double lat) const 
  { 
    return topography_data.value( 
      Point<2>(-lat * 180 / numbers::PI, lon * 180 / numbers::PI)); 
  } 



  std::vector<double> AfricaTopography::get_data() 
  { 
    std::vector<double> data; 


    boost::iostreams::filtering_istream in; 
    in.push(boost::iostreams::basic_gzip_decompressor<>()); 
    in.push(boost::iostreams::file_source("topography.txt.gz")); 

    for (unsigned int line = 0; line < 83600; ++line) 
      { 
        try 
          { 
            double lat, lon, elevation; 
            in >> lat >> lon >> elevation; 

            data.push_back(elevation); 
          } 
        catch (...) 
          { 
            AssertThrow(false, 
                        ExcMessage("Could not read all 83,600 data points " 
                                   "from the file <topography.txt.gz>!")); 
          } 
      } 

    return data; 
  } 


  class AfricaGeometry : public ChartManifold<3, 3> 
  { 
  public: 
    virtual Point<3> pull_back(const Point<3> &space_point) const override; 

    virtual Point<3> push_forward(const Point<3> &chart_point) const override; 

    virtual std::unique_ptr<Manifold<3, 3>> clone() const override; 

  private: 
    static const double R; 
    static const double ellipticity; 

    const AfricaTopography topography; 

    Point<3> push_forward_wgs84(const Point<3> &phi_theta_d) const; 
    Point<3> pull_back_wgs84(const Point<3> &x) const; 

    Point<3> push_forward_topo(const Point<3> &phi_theta_d_hat) const; 
    Point<3> pull_back_topo(const Point<3> &phi_theta_d) const; 
  }; 

  const double AfricaGeometry::R           = 6378137; 
  const double AfricaGeometry::ellipticity = 8.1819190842622e-2; 


  Point<3> AfricaGeometry::pull_back(const Point<3> &space_point) const 
  { 
    return pull_back_topo(pull_back_wgs84(space_point)); 
  } 

  Point<3> AfricaGeometry::push_forward(const Point<3> &chart_point) const 
  { 
    return push_forward_wgs84(push_forward_topo(chart_point)); 
  } 


  std::unique_ptr<Manifold<3, 3>> AfricaGeometry::clone() const 
  { 
    return std::make_unique<AfricaGeometry>(); 
  } 


  Point<3> AfricaGeometry::push_forward_wgs84(const Point<3> &phi_theta_d) const 
  { 
    const double phi   = phi_theta_d[0]; 
    const double theta = phi_theta_d[1]; 
    const double d     = phi_theta_d[2]; 

    const double R_bar = R / std::sqrt(1 - (ellipticity * ellipticity * 
                                            std::sin(theta) * std::sin(theta))); 

    return {(R_bar + d) * std::cos(phi) * std::cos(theta), 
            (R_bar + d) * std::sin(phi) * std::cos(theta), 
            ((1 - ellipticity * ellipticity) * R_bar + d) * std::sin(theta)}; 
  } 

  Point<3> AfricaGeometry::pull_back_wgs84(const Point<3> &x) const 
  { 
    const double b   = std::sqrt(R * R * (1 - ellipticity * ellipticity)); 
    const double ep  = std::sqrt((R * R - b * b) / (b * b)); 
    const double p   = std::sqrt(x(0) * x(0) + x(1) * x(1)); 
    const double th  = std::atan2(R * x(2), b * p); 
    const double phi = std::atan2(x(1), x(0)); 
    const double theta = 
      std::atan2(x(2) + ep * ep * b * std::pow(std::sin(th), 3), 
                 (p - 
                  (ellipticity * ellipticity * R * std::pow(std::cos(th), 3)))); 
    const double R_bar = 
      R / (std::sqrt(1 - ellipticity * ellipticity * std::sin(theta) * 
                           std::sin(theta))); 
    const double R_plus_d = p / std::cos(theta); 

    Point<3> phi_theta_d; 
    if (phi < 0) 
      phi_theta_d[0] = phi + 2 * numbers::PI; 
    else if (phi > 2 * numbers::PI) 
      phi_theta_d[0] = phi - 2 * numbers::PI; 
    else 
      phi_theta_d[0] = phi; 
    phi_theta_d[1] = theta; 
    phi_theta_d[2] = R_plus_d - R_bar; 
    return phi_theta_d; 
  } 


  Point<3> 
  AfricaGeometry::push_forward_topo(const Point<3> &phi_theta_d_hat) const 
  { 
    const double d_hat = phi_theta_d_hat[2]; 
    const double h = topography.value(phi_theta_d_hat[0], phi_theta_d_hat[1]); 
    const double d = d_hat + (d_hat + 500000) / 500000 * h; 
    return {phi_theta_d_hat[0], phi_theta_d_hat[1], d}; 
  } 

  Point<3> AfricaGeometry::pull_back_topo(const Point<3> &phi_theta_d) const 
  { 
    const double d     = phi_theta_d[2]; 
    const double h     = topography.value(phi_theta_d[0], phi_theta_d[1]); 
    const double d_hat = 500000 * (d - h) / (500000 + h); 
    return {phi_theta_d[0], phi_theta_d[1], d_hat}; 
  } 



  void run() 
  { 
    AfricaGeometry   geometry; 
    Triangulation<3> triangulation; 

    { 
      const Point<3> corner_points[2] = { 
        Point<3>(26 * numbers::PI / 180, -10 * numbers::PI / 180, -500000), 
        Point<3>(35 * numbers::PI / 180, 5 * numbers::PI / 180, 0)}; 
      std::vector<unsigned int> subdivisions(3); 
      subdivisions[0] = 1; 
      subdivisions[1] = 2; 
      subdivisions[2] = 1; 
      GridGenerator::subdivided_hyper_rectangle( 
        triangulation, subdivisions, corner_points[0], corner_points[1], true); 

      GridTools::transform( 
        [&geometry](const Point<3> &chart_point) { 
          return geometry.push_forward(chart_point); 
        }, 
        triangulation); 
    } 


    triangulation.set_manifold(0, geometry); 
    for (const auto &cell : triangulation.active_cell_iterators()) 
      cell->set_all_manifold_ids(0); 


    for (unsigned int i = 0; i < 6; ++i) 
      { 
        for (const auto &cell : triangulation.active_cell_iterators()) 
          for (const auto &face : cell->face_iterators()) 
            if (face->boundary_id() == 5) 
              { 
                cell->set_refine_flag(); 
                break; 
              } 
        triangulation.execute_coarsening_and_refinement(); 

        std::cout << "Refinement step " << i + 1 << ": " 
                  << triangulation.n_active_cells() << " cells, " 
                  << GridTools::minimal_cell_diameter(triangulation) / 1000 
                  << "km minimal cell diameter" << std::endl; 
      } 


    const std::string filename = "mesh.vtu"; 
    std::ofstream     out(filename); 
    GridOut           grid_out; 
    grid_out.write_vtu(triangulation, out); 
  } 
} // namespace Step53 



int main() 
{ 
  try 
    { 
      Step53::run(); 
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
} 


