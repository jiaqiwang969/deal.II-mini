

/* --------------------------------------------------------------------- 
 * 
 * Copyright (C) 2021 by the deal.II authors 
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
 * Author: Jean-Paul Pelteret, 2021 
 */ 




#include <deal.II/algorithms/general_data_storage.h> 


#include <deal.II/base/discrete_time.h> 
#include <deal.II/base/numbers.h> 
#include <deal.II/base/parameter_acceptor.h> 
#include <deal.II/base/symmetric_tensor.h> 
#include <deal.II/base/tensor.h> 
#include <deal.II/base/timer.h> 
#include <deal.II/base/utilities.h> 


#include <deal.II/physics/transformations.h> 
#include <deal.II/physics/elasticity/kinematics.h> 
#include <deal.II/physics/elasticity/standard_tensors.h> 


#include <deal.II/differentiation/ad.h> 
#include <deal.II/differentiation/sd.h> 


#include <fstream> 


namespace Step71 
{ 
  using namespace dealii; 




  namespace SimpleExample 
  { 



    template <typename NumberType> 
    NumberType f(const NumberType &x, const NumberType &y) 
    { 
      return std::cos(y / x); 
    } 


    double df_dx(const double x, const double y); 

    double df_dy(const double x, const double y); 

    double d2f_dx_dx(const double x, const double y); 

    double d2f_dx_dy(const double x, const double y); 

    double d2f_dy_dx(const double x, const double y); 


    double d2f_dy_dy(const double x, const double y); 


    void 
    run_and_verify_ad(const double x, const double y, const double tol = 1e-12) 
    { 


      constexpr unsigned int                     dim = 1; 
      constexpr Differentiation::AD::NumberTypes ADTypeCode = 
        Differentiation::AD::NumberTypes::sacado_dfad_dfad; 
      using ADHelper = 
        Differentiation::AD::ScalarFunction<dim, ADTypeCode, double>; 


      constexpr unsigned int n_independent_variables = 2; 


      ADHelper ad_helper(n_independent_variables); 
      using ADNumberType = typename ADHelper::ad_type; 


      ad_helper.register_independent_variables({x, y}); 


      const std::vector<ADNumberType> independent_variables_ad = 
        ad_helper.get_sensitive_variables(); 
      const ADNumberType &x_ad = independent_variables_ad[0]; 
      const ADNumberType &y_ad = independent_variables_ad[1]; 


      const ADNumberType f_ad = f(x_ad, y_ad); 











      ad_helper.register_dependent_variable(f_ad); 


      Vector<double>     Df(ad_helper.n_dependent_variables()); 
      FullMatrix<double> D2f(ad_helper.n_dependent_variables(), 
                             ad_helper.n_independent_variables()); 


      const double computed_f = ad_helper.compute_value(); 
      ad_helper.compute_gradient(Df); 
      ad_helper.compute_hessian(D2f); 


      AssertThrow(std::abs(f(x, y) - computed_f) < tol, 
                  ExcMessage(std::string("Incorrect value computed for f. ") + 
                             std::string("Hand-calculated value: ") + 
                             Utilities::to_string(f(x, y)) + 
                             std::string(" ; ") + 
                             std::string("Value computed by AD: ") + 
                             Utilities::to_string(computed_f))); 


      const double computed_df_dx = Df[0]; 
      const double computed_df_dy = Df[1]; 

      AssertThrow(std::abs(df_dx(x, y) - computed_df_dx) < tol, 
                  ExcMessage( 
                    std::string("Incorrect value computed for df/dx. ") + 
                    std::string("Hand-calculated value: ") + 
                    Utilities::to_string(df_dx(x, y)) + std::string(" ; ") + 
                    std::string("Value computed by AD: ") + 
                    Utilities::to_string(computed_df_dx))); 
      AssertThrow(std::abs(df_dy(x, y) - computed_df_dy) < tol, 
                  ExcMessage( 
                    std::string("Incorrect value computed for df/dy. ") + 
                    std::string("Hand-calculated value: ") + 
                    Utilities::to_string(df_dy(x, y)) + std::string(" ; ") + 
                    std::string("Value computed by AD: ") + 
                    Utilities::to_string(computed_df_dy))); 


      const double computed_d2f_dx_dx = D2f[0][0]; 
      const double computed_d2f_dx_dy = D2f[0][1]; 
      const double computed_d2f_dy_dx = D2f[1][0]; 
      const double computed_d2f_dy_dy = D2f[1][1]; 

      AssertThrow(std::abs(d2f_dx_dx(x, y) - computed_d2f_dx_dx) < tol, 
                  ExcMessage( 
                    std::string("Incorrect value computed for d2f/dx_dx. ") + 
                    std::string("Hand-calculated value: ") + 
                    Utilities::to_string(d2f_dx_dx(x, y)) + std::string(" ; ") + 
                    std::string("Value computed by AD: ") + 
                    Utilities::to_string(computed_d2f_dx_dx))); 
      AssertThrow(std::abs(d2f_dx_dy(x, y) - computed_d2f_dx_dy) < tol, 
                  ExcMessage( 
                    std::string("Incorrect value computed for d2f/dx_dy. ") + 
                    std::string("Hand-calculated value: ") + 
                    Utilities::to_string(d2f_dx_dy(x, y)) + std::string(" ; ") + 
                    std::string("Value computed by AD: ") + 
                    Utilities::to_string(computed_d2f_dx_dy))); 
      AssertThrow(std::abs(d2f_dy_dx(x, y) - computed_d2f_dy_dx) < tol, 
                  ExcMessage( 
                    std::string("Incorrect value computed for d2f/dy_dx. ") + 
                    std::string("Hand-calculated value: ") + 
                    Utilities::to_string(d2f_dy_dx(x, y)) + std::string(" ; ") + 
                    std::string("Value computed by AD: ") + 
                    Utilities::to_string(computed_d2f_dy_dx))); 
      AssertThrow(std::abs(d2f_dy_dy(x, y) - computed_d2f_dy_dy) < tol, 
                  ExcMessage( 
                    std::string("Incorrect value computed for d2f/dy_dy. ") + 
                    std::string("Hand-calculated value: ") + 
                    Utilities::to_string(d2f_dy_dy(x, y)) + std::string(" ; ") + 
                    std::string("Value computed by AD: ") + 
                    Utilities::to_string(computed_d2f_dy_dy))); 
    } 





    double df_dx(const double x, const double y) 
    { 
      Assert(x != 0.0, ExcDivideByZero()); 
      return y * std::sin(y / x) / (x * x); 
    } 
    double df_dy(const double x, const double y) 
    { 
      return -std::sin(y / x) / x; 
    } 


    double d2f_dx_dx(const double x, const double y) 
    { 
      return -y * (2 * x * std::sin(y / x) + y * std::cos(y / x)) / 
             (x * x * x * x); 
    } 
    double d2f_dx_dy(const double x, const double y) 
    { 
      return (x * std::sin(y / x) + y * std::cos(y / x)) / (x * x * x); 
    } 

    double d2f_dy_dx(const double x, const double y) 
    { 
      return (x * std::sin(y / x) + y * std::cos(y / x)) / (x * x * x); 
    } 
    double d2f_dy_dy(const double x, const double y) 
    { 
      return -(std::cos(y / x)) / (x * x); 
    } 





    void 
    run_and_verify_sd(const double x, const double y, const double tol = 1e-12) 
    { 


      const Differentiation::SD::Expression x_sd("x"); 
      const Differentiation::SD::Expression y_sd("y"); 


      const Differentiation::SD::Expression f_sd = f(x_sd, y_sd); 






      const Differentiation::SD::Expression df_dx_sd = f_sd.differentiate(x_sd); 
      const Differentiation::SD::Expression df_dy_sd = f_sd.differentiate(y_sd); 


      const Differentiation::SD::Expression d2f_dx_dx_sd = 
        df_dx_sd.differentiate(x_sd); 
      const Differentiation::SD::Expression d2f_dx_dy_sd = 
        df_dx_sd.differentiate(y_sd); 
      const Differentiation::SD::Expression d2f_dy_dx_sd = 
        df_dy_sd.differentiate(x_sd); 
      const Differentiation::SD::Expression d2f_dy_dy_sd = 
        df_dy_sd.differentiate(y_sd); 



      const Differentiation::SD::types::substitution_map substitution_map = 
        Differentiation::SD::make_substitution_map( 
          std::pair<Differentiation::SD::Expression, double>{x_sd, x}, 
          std::pair<Differentiation::SD::Expression, double>{y_sd, y}); 




      const double computed_f = 
        f_sd.substitute_and_evaluate<double>(substitution_map); 

      AssertThrow(std::abs(f(x, y) - computed_f) < tol, 
                  ExcMessage(std::string("Incorrect value computed for f. ") + 
                             std::string("Hand-calculated value: ") + 
                             Utilities::to_string(f(x, y)) + 
                             std::string(" ; ") + 
                             std::string("Value computed by AD: ") + 
                             Utilities::to_string(computed_f))); 


      const double computed_df_dx = 
        df_dx_sd.substitute_and_evaluate<double>(substitution_map); 
      const double computed_df_dy = 
        df_dy_sd.substitute_and_evaluate<double>(substitution_map); 

      AssertThrow(std::abs(df_dx(x, y) - computed_df_dx) < tol, 
                  ExcMessage( 
                    std::string("Incorrect value computed for df/dx. ") + 
                    std::string("Hand-calculated value: ") + 
                    Utilities::to_string(df_dx(x, y)) + std::string(" ; ") + 
                    std::string("Value computed by AD: ") + 
                    Utilities::to_string(computed_df_dx))); 
      AssertThrow(std::abs(df_dy(x, y) - computed_df_dy) < tol, 
 
 
 
 
 
                    Utilities::to_string(computed_df_dy))); 


      const double computed_d2f_dx_dx = 
        d2f_dx_dx_sd.substitute_and_evaluate<double>(substitution_map); 
      const double computed_d2f_dx_dy = 
        d2f_dx_dy_sd.substitute_and_evaluate<double>(substitution_map); 
      const double computed_d2f_dy_dx = 
        d2f_dy_dx_sd.substitute_and_evaluate<double>(substitution_map); 
      const double computed_d2f_dy_dy = 
        d2f_dy_dy_sd.substitute_and_evaluate<double>(substitution_map); 

      AssertThrow(std::abs(d2f_dx_dx(x, y) - computed_d2f_dx_dx) < tol, 
                  ExcMessage( 
                    std::string("Incorrect value computed for d2f/dx_dx. ") + 
                    std::string("Hand-calculated value: ") + 
                    Utilities::to_string(d2f_dx_dx(x, y)) + std::string(" ; ") + 
                    std::string("Value computed by SD: ") + 
                    Utilities::to_string(computed_d2f_dx_dx))); 
      AssertThrow(std::abs(d2f_dx_dy(x, y) - computed_d2f_dx_dy) < tol, 
                  ExcMessage( 
                    std::string("Incorrect value computed for d2f/dx_dy. ") + 
                    std::string("Hand-calculated value: ") + 
                    Utilities::to_string(d2f_dx_dy(x, y)) + std::string(" ; ") + 
                    std::string("Value computed by SD: ") + 
                    Utilities::to_string(computed_d2f_dx_dy))); 
      AssertThrow(std::abs(d2f_dy_dx(x, y) - computed_d2f_dy_dx) < tol, 
                  ExcMessage( 
                    std::string("Incorrect value computed for d2f/dy_dx. ") + 
                    std::string("Hand-calculated value: ") + 
                    Utilities::to_string(d2f_dy_dx(x, y)) + std::string(" ; ") + 
                    std::string("Value computed by SD: ") + 
                    Utilities::to_string(computed_d2f_dy_dx))); 
      AssertThrow(std::abs(d2f_dy_dy(x, y) - computed_d2f_dy_dy) < tol, 
                  ExcMessage( 
                    std::string("Incorrect value computed for d2f/dy_dy. ") + 
                    std::string("Hand-calculated value: ") + 
                    Utilities::to_string(d2f_dy_dy(x, y)) + std::string(" ; ") + 
                    std::string("Value computed by SD: ") + 
                    Utilities::to_string(computed_d2f_dy_dy))); 
    } 


    void run() 
    { 
      const double x = 1.23; 
      const double y = 0.91; 

      std::cout << "Simple example using automatic differentiation..." 
                << std::endl; 
      run_and_verify_ad(x, y); 
      std::cout << "... all calculations are correct!" << std::endl; 

      std::cout << "Simple example using symbolic differentiation." 
                << std::endl; 
      run_and_verify_sd(x, y); 
      std::cout << "... all calculations are correct!" << std::endl; 
    } 

  } // namespace SimpleExample 



  namespace CoupledConstitutiveLaws 
  { 








    class ConstitutiveParameters : public ParameterAcceptor 
    { 
    public: 
      ConstitutiveParameters(); 

      double mu_e       = 30.0e3; 
      double mu_e_inf   = 250.0e3; 
      double mu_e_h_sat = 212.2e3; 
      double nu_e       = 0.49; 






      double mu_v       = 20.0e3; 
      double mu_v_inf   = 35.0e3; 
      double mu_v_h_sat = 92.84e3; 
      double tau_v      = 0.6; 


      double mu_r = 6.0; 

      bool initialized = false; 
    }; 


    ConstitutiveParameters::ConstitutiveParameters() 
      : ParameterAcceptor("/Coupled Constitutive Laws/Constitutive Parameters/") 
    { 
      add_parameter("Elastic shear modulus", mu_e); 
      add_parameter("Elastic shear modulus at magnetic saturation", mu_e_inf); 
      add_parameter( 
        "Saturation magnetic field strength for elastic shear modulus", 
        mu_e_h_sat); 
      add_parameter("Poisson ratio", nu_e); 

      add_parameter("Viscoelastic shear modulus", mu_v); 
      add_parameter("Viscoelastic shear modulus at magnetic saturation", 
                    mu_v_inf); 
      add_parameter( 
        "Saturation magnetic field strength for viscoelastic shear modulus", 
        mu_v_h_sat); 
      add_parameter("Characteristic relaxation time", tau_v); 

      add_parameter("Relative magnetic permeability", mu_r); 

      parse_parameters_call_back.connect([&]() { initialized = true; }); 
    } 



    template <int dim> 
    class Coupled_Magnetomechanical_Constitutive_Law_Base 
    { 
    public: 
      Coupled_Magnetomechanical_Constitutive_Law_Base( 
        const ConstitutiveParameters &constitutive_parameters); 


      virtual void update_internal_data(const SymmetricTensor<2, dim> &C, 
                                        const Tensor<1, dim> &         H, 
                                        const DiscreteTime &time) = 0; 



      virtual double get_psi() const = 0; 



      virtual Tensor<1, dim> get_B() const = 0; 

      virtual SymmetricTensor<2, dim> get_S() const = 0; 





      virtual SymmetricTensor<2, dim> get_DD() const = 0; 

      virtual Tensor<3, dim> get_PP() const = 0; 

      virtual SymmetricTensor<4, dim> get_HH() const = 0; 


      virtual void update_end_of_timestep() 
      {} 


与材料的弹性响应有关的参数依次是：//。







    protected: 
      const ConstitutiveParameters &constitutive_parameters; 

      double get_mu_e() const; 

      double get_mu_e_inf() const; 

      double get_mu_e_h_sat() const; 

      double get_nu_e() const; 

      double get_lambda_e() const; 

      double get_kappa_e() const; 





粘弹性剪切模量的饱和磁场强度，和//--特征松弛时间。

      double get_mu_v() const; 

      double get_mu_v_inf() const; 

      double get_mu_v_h_sat() const; 

      double get_tau_v() const; 





      double get_mu_r() const; 

      constexpr double get_mu_0() const; 
      double           get_delta_t(const DiscreteTime &time) const; 
    }; 


    template <int dim> 
    Coupled_Magnetomechanical_Constitutive_Law_Base<dim>:: 
      Coupled_Magnetomechanical_Constitutive_Law_Base( 
        const ConstitutiveParameters &constitutive_parameters) 
      : constitutive_parameters(constitutive_parameters) 
    { 
      Assert(get_kappa_e() > 0, ExcInternalError()); 
    } 

    template <int dim> 
    double 
    Coupled_Magnetomechanical_Constitutive_Law_Base<dim>::get_mu_e() const 
    { 
      return constitutive_parameters.mu_e; 
    } 

    template <int dim> 
    double 
    Coupled_Magnetomechanical_Constitutive_Law_Base<dim>::get_mu_e_inf() const 
    { 
      return constitutive_parameters.mu_e_inf; 
    } 

    template <int dim> 
    double 
    Coupled_Magnetomechanical_Constitutive_Law_Base<dim>::get_mu_e_h_sat() const 
    { 
      return constitutive_parameters.mu_e_h_sat; 
    } 

    template <int dim> 
    double 
    Coupled_Magnetomechanical_Constitutive_Law_Base<dim>::get_nu_e() const 
    { 
      return constitutive_parameters.nu_e; 
    } 

    template <int dim> 
    double 
    Coupled_Magnetomechanical_Constitutive_Law_Base<dim>::get_lambda_e() const 
    { 
      return 2.0 * get_mu_e() * get_nu_e() / (1.0 - 2.0 * get_nu_e()); 
    } 

    template <int dim> 
    double 
    Coupled_Magnetomechanical_Constitutive_Law_Base<dim>::get_kappa_e() const 
    { 
      return (2.0 * get_mu_e() * (1.0 + get_nu_e())) / 
             (3.0 * (1.0 - 2.0 * get_nu_e())); 
    } 

    template <int dim> 
    double 
    Coupled_Magnetomechanical_Constitutive_Law_Base<dim>::get_mu_v() const 
    { 
      return constitutive_parameters.mu_v; 
    } 

    template <int dim> 
    double 
    Coupled_Magnetomechanical_Constitutive_Law_Base<dim>::get_mu_v_inf() const 
    { 
      return constitutive_parameters.mu_v_inf; 
    } 

    template <int dim> 
    double 
    Coupled_Magnetomechanical_Constitutive_Law_Base<dim>::get_mu_v_h_sat() const 
    { 
      return constitutive_parameters.mu_v_h_sat; 
    } 

    template <int dim> 
    double 
    Coupled_Magnetomechanical_Constitutive_Law_Base<dim>::get_tau_v() const 
    { 
      return constitutive_parameters.tau_v; 
    } 

    template <int dim> 
    double 
    Coupled_Magnetomechanical_Constitutive_Law_Base<dim>::get_mu_r() const 
    { 
      return constitutive_parameters.mu_r; 
    } 

    template <int dim> 
    constexpr double 
    Coupled_Magnetomechanical_Constitutive_Law_Base<dim>::get_mu_0() const 
    { 
      return 4.0 * numbers::PI * 1e-7; 
    } 

    template <int dim> 
    double Coupled_Magnetomechanical_Constitutive_Law_Base<dim>::get_delta_t( 
      const DiscreteTime &time) const 
    { 
      return time.get_previous_step_size(); 
    } 



    template <int dim, Differentiation::AD::NumberTypes ADTypeCode> 
    class Magnetoelastic_Constitutive_Law_AD final 
      : public Coupled_Magnetomechanical_Constitutive_Law_Base<dim> 
    { 
      using ADHelper = 
        Differentiation::AD::ScalarFunction<dim, ADTypeCode, double>; 
      using ADNumberType = typename ADHelper::ad_type; 

    public: 
      Magnetoelastic_Constitutive_Law_AD( 
        const ConstitutiveParameters &constitutive_parameters); 


      virtual void update_internal_data(const SymmetricTensor<2, dim> &C, 
                                        const Tensor<1, dim> &         H, 
                                        const DiscreteTime &) override; 

      virtual double get_psi() const override; 

      virtual Tensor<1, dim> get_B() const override; 

      virtual SymmetricTensor<2, dim> get_S() const override; 

      virtual SymmetricTensor<2, dim> get_DD() const override; 

      virtual Tensor<3, dim> get_PP() const override; 

      virtual SymmetricTensor<4, dim> get_HH() const override; 


    private: 
      const FEValuesExtractors::Vector             H_components; 
      const FEValuesExtractors::SymmetricTensor<2> C_components; 


      ADHelper ad_helper; 


      double             psi; 
      Vector<double>     Dpsi; 
      FullMatrix<double> D2psi; 
    }; 


    template <int dim, Differentiation::AD::NumberTypes ADTypeCode> 
    Magnetoelastic_Constitutive_Law_AD<dim, ADTypeCode>:: 
      Magnetoelastic_Constitutive_Law_AD( 
        const ConstitutiveParameters &constitutive_parameters) 
      : Coupled_Magnetomechanical_Constitutive_Law_Base<dim>( 
          constitutive_parameters) 
      , H_components(0) 
      , C_components(Tensor<1, dim>::n_independent_components) 
      , ad_helper(Tensor<1, dim>::n_independent_components + 
                  SymmetricTensor<2, dim>::n_independent_components) 
      , psi(0.0) 
      , Dpsi(ad_helper.n_independent_variables()) 
      , D2psi(ad_helper.n_independent_variables(), 
              ad_helper.n_independent_variables()) 
    {} 


    template <int dim, Differentiation::AD::NumberTypes ADTypeCode> 
    void 
    Magnetoelastic_Constitutive_Law_AD<dim, ADTypeCode>::update_internal_data( 
      const SymmetricTensor<2, dim> &C, 
      const Tensor<1, dim> &         H, 
      const DiscreteTime &) 
    { 
      Assert(determinant(C) > 0, ExcInternalError()); 


      ad_helper.reset(); 


      ad_helper.register_independent_variable(H, H_components); 
      ad_helper.register_independent_variable(C, C_components); 


      const Tensor<1, dim, ADNumberType> H_ad = 
        ad_helper.get_sensitive_variables(H_components); 
      const SymmetricTensor<2, dim, ADNumberType> C_ad = 
        ad_helper.get_sensitive_variables(C_components); 


      const ADNumberType det_F_ad = std::sqrt(determinant(C_ad)); 
      const SymmetricTensor<2, dim, ADNumberType> C_inv_ad = invert(C_ad); 
      AssertThrow(det_F_ad > ADNumberType(0.0), 
                  ExcMessage("Volumetric Jacobian must be positive.")); 


      const ADNumberType f_mu_e_ad = 
        1.0 + (this->get_mu_e_inf() / this->get_mu_e() - 1.0) * 
                std::tanh((2.0 * H_ad * H_ad) / 
                          (this->get_mu_e_h_sat() * this->get_mu_e_h_sat())); 


      const ADNumberType psi_ad = 
        0.5 * this->get_mu_e() * f_mu_e_ad * 
          (trace(C_ad) - dim - 2.0 * std::log(det_F_ad))                 // 
        + this->get_lambda_e() * std::log(det_F_ad) * std::log(det_F_ad) // 
        - 0.5 * this->get_mu_0() * this->get_mu_r() * det_F_ad * 
            (H_ad * C_inv_ad * H_ad); // 

      ad_helper.register_dependent_variable(psi_ad); 


      psi = ad_helper.compute_value(); 
      ad_helper.compute_gradient(Dpsi); 
      ad_helper.compute_hessian(D2psi); 
    } 


    template <int dim, Differentiation::AD::NumberTypes ADTypeCode> 
    double Magnetoelastic_Constitutive_Law_AD<dim, ADTypeCode>::get_psi() const 
    { 
      return psi; 
    } 

    template <int dim, Differentiation::AD::NumberTypes ADTypeCode> 
    Tensor<1, dim> 
    Magnetoelastic_Constitutive_Law_AD<dim, ADTypeCode>::get_B() const 
    { 
      const Tensor<1, dim> dpsi_dH = 
        ad_helper.extract_gradient_component(Dpsi, H_components); 
      return -dpsi_dH; 
    } 

    template <int dim, Differentiation::AD::NumberTypes ADTypeCode> 
    SymmetricTensor<2, dim> 
    Magnetoelastic_Constitutive_Law_AD<dim, ADTypeCode>::get_S() const 
    { 
      const SymmetricTensor<2, dim> dpsi_dC = 
        ad_helper.extract_gradient_component(Dpsi, C_components); 
      return 2.0 * dpsi_dC; 
    } 

    template <int dim, Differentiation::AD::NumberTypes ADTypeCode> 
    SymmetricTensor<2, dim> 
    Magnetoelastic_Constitutive_Law_AD<dim, ADTypeCode>::get_DD() const 
    { 
      const Tensor<2, dim> dpsi_dH_dH = 
        ad_helper.extract_hessian_component(D2psi, H_components, H_components); 
      return -symmetrize(dpsi_dH_dH); 
    } 


    template <int dim, Differentiation::AD::NumberTypes ADTypeCode> 
    Tensor<3, dim> 
    Magnetoelastic_Constitutive_Law_AD<dim, ADTypeCode>::get_PP() const 
    { 
      const Tensor<3, dim> dpsi_dC_dH = 
        ad_helper.extract_hessian_component(D2psi, C_components, H_components); 
      return -2.0 * dpsi_dC_dH; 
    } 

    template <int dim, Differentiation::AD::NumberTypes ADTypeCode> 
    SymmetricTensor<4, dim> 
    Magnetoelastic_Constitutive_Law_AD<dim, ADTypeCode>::get_HH() const 
    { 
      const SymmetricTensor<4, dim> dpsi_dC_dC = 
        ad_helper.extract_hessian_component(D2psi, C_components, C_components); 
      return 4.0 * dpsi_dC_dC; 
    } 




    template <int dim> 
    class Magnetoviscoelastic_Constitutive_Law_SD final 
      : public Coupled_Magnetomechanical_Constitutive_Law_Base<dim> 
    { 
    public: 
      Magnetoviscoelastic_Constitutive_Law_SD( 
        const ConstitutiveParameters &               constitutive_parameters, 
        const Differentiation::SD::OptimizerType     optimizer_type, 
        const Differentiation::SD::OptimizationFlags optimization_flags); 


      virtual void update_internal_data(const SymmetricTensor<2, dim> &C, 
                                        const Tensor<1, dim> &         H, 
                                        const DiscreteTime &time) override; 

      virtual double get_psi() const override; 

      virtual Tensor<1, dim> get_B() const override; 

      virtual SymmetricTensor<2, dim> get_S() const override; 

      virtual SymmetricTensor<2, dim> get_DD() const override; 

      virtual Tensor<3, dim> get_PP() const override; 

      virtual SymmetricTensor<4, dim> get_HH() const override; 


      virtual void update_end_of_timestep() override; 





    private: 
      SymmetricTensor<2, dim> Q_t; 
      SymmetricTensor<2, dim> Q_t1; 



      const Differentiation::SD::Expression mu_e_sd; 
      const Differentiation::SD::Expression mu_e_inf_sd; 
      const Differentiation::SD::Expression mu_e_h_sat_sd; 
      const Differentiation::SD::Expression lambda_e_sd; 
      const Differentiation::SD::Expression mu_v_sd; 
      const Differentiation::SD::Expression mu_v_inf_sd; 
      const Differentiation::SD::Expression mu_v_h_sat_sd; 
      const Differentiation::SD::Expression tau_v_sd; 
      const Differentiation::SD::Expression delta_t_sd; 
      const Differentiation::SD::Expression mu_r_sd; 


      const Tensor<1, dim, Differentiation::SD::Expression>          H_sd; 
      const SymmetricTensor<2, dim, Differentiation::SD::Expression> C_sd; 


      const SymmetricTensor<2, dim, Differentiation::SD::Expression> Q_t_sd; 
      const SymmetricTensor<2, dim, Differentiation::SD::Expression> Q_t1_sd; 


      Differentiation::SD::Expression                          psi_sd; 
      Tensor<1, dim, Differentiation::SD::Expression>          B_sd; 
      SymmetricTensor<2, dim, Differentiation::SD::Expression> S_sd; 
      SymmetricTensor<2, dim, Differentiation::SD::Expression> BB_sd; 
      Tensor<3, dim, Differentiation::SD::Expression>          PP_sd; 
      SymmetricTensor<4, dim, Differentiation::SD::Expression> HH_sd; 




      Differentiation::SD::BatchOptimizer<double> optimizer; 


      Differentiation::SD::types::substitution_map 
      make_substitution_map(const SymmetricTensor<2, dim> &C, 
                            const Tensor<1, dim> &         H, 
                            const double                   delta_t) const; 

      void initialize_optimizer(); 
    }; 


    template <int dim> 
    Magnetoviscoelastic_Constitutive_Law_SD<dim>:: 
      Magnetoviscoelastic_Constitutive_Law_SD( 
        const ConstitutiveParameters &               constitutive_parameters, 
        const Differentiation::SD::OptimizerType     optimizer_type, 
        const Differentiation::SD::OptimizationFlags optimization_flags) 
      : Coupled_Magnetomechanical_Constitutive_Law_Base<dim>( 
          constitutive_parameters) 
      , Q_t(Physics::Elasticity::StandardTensors<dim>::I) 
      , Q_t1(Physics::Elasticity::StandardTensors<dim>::I) 
      , mu_e_sd("mu_e") 
      , mu_e_inf_sd("mu_e_inf") 
      , mu_e_h_sat_sd("mu_e_h_sat") 
      , lambda_e_sd("lambda_e") 
      , mu_v_sd("mu_v") 
      , mu_v_inf_sd("mu_v_inf") 
      , mu_v_h_sat_sd("mu_v_h_sat") 
      , tau_v_sd("tau_v") 
      , delta_t_sd("delta_t") 
      , mu_r_sd("mu_r") 
      , H_sd(Differentiation::SD::make_vector_of_symbols<dim>("H")) 
      , C_sd(Differentiation::SD::make_symmetric_tensor_of_symbols<2, dim>("C")) 
      , Q_t_sd( 
          Differentiation::SD::make_symmetric_tensor_of_symbols<2, dim>("Q_t")) 
      , Q_t1_sd( 
          Differentiation::SD::make_symmetric_tensor_of_symbols<2, dim>("Q_t1")) 
      , optimizer(optimizer_type, optimization_flags) 
    { 
      initialize_optimizer(); 
    } 






    template <int dim> 
    Differentiation::SD::types::substitution_map 
    Magnetoviscoelastic_Constitutive_Law_SD<dim>::make_substitution_map( 
      const SymmetricTensor<2, dim> &C, 
      const Tensor<1, dim> &         H, 
      const double                   delta_t) const 
    { 
      return Differentiation::SD::make_substitution_map( 
        std::make_pair(mu_e_sd, this->get_mu_e()), 
        std::make_pair(mu_e_inf_sd, this->get_mu_e_inf()), 
        std::make_pair(mu_e_h_sat_sd, this->get_mu_e_h_sat()), 
        std::make_pair(lambda_e_sd, this->get_lambda_e()), 
        std::make_pair(mu_v_sd, this->get_mu_v()), 
        std::make_pair(mu_v_inf_sd, this->get_mu_v_inf()), 
        std::make_pair(mu_v_h_sat_sd, this->get_mu_v_h_sat()), 
        std::make_pair(tau_v_sd, this->get_tau_v()), 
        std::make_pair(delta_t_sd, delta_t), 
        std::make_pair(mu_r_sd, this->get_mu_r()), 
        std::make_pair(H_sd, H), 
        std::make_pair(C_sd, C), 
        std::make_pair(Q_t_sd, Q_t), 
        std::make_pair(Q_t1_sd, Q_t1)); 
    } 



    template <int dim> 
    void Magnetoviscoelastic_Constitutive_Law_SD<dim>::initialize_optimizer() 
    { 
      const Differentiation::SD::Expression det_F_sd = 
        std::sqrt(determinant(C_sd)); 
      const SymmetricTensor<2, dim, Differentiation::SD::Expression> C_inv_sd = 
        invert(C_sd); 


      const Differentiation::SD::Expression f_mu_e_sd = 
        1.0 + 
        (mu_e_inf_sd / mu_e_sd - 1.0) * 
          std::tanh((2.0 * H_sd * H_sd) / (mu_e_h_sat_sd * mu_e_h_sat_sd)); 

      const Differentiation::SD::Expression psi_ME_sd = 
        0.5 * mu_e_sd * f_mu_e_sd * 
          (trace(C_sd) - dim - 2.0 * std::log(det_F_sd)) + 
        lambda_e_sd * std::log(det_F_sd) * std::log(det_F_sd) - 
        0.5 * this->get_mu_0() * mu_r_sd * det_F_sd * (H_sd * C_inv_sd * H_sd); 


      const Differentiation::SD::Expression f_mu_v_sd = 
        1.0 + 
        (mu_v_inf_sd / mu_v_sd - 1.0) * 
          std::tanh((2.0 * H_sd * H_sd) / (mu_v_h_sat_sd * mu_v_h_sat_sd)); 

      const Differentiation::SD::Expression psi_MVE_sd = 
        0.5 * mu_v_sd * f_mu_v_sd * 
        (Q_t_sd * (std::pow(det_F_sd, -2.0 / dim) * C_sd) - dim - 
         std::log(determinant(Q_t_sd))); 


      psi_sd = psi_ME_sd + psi_MVE_sd; 




      B_sd = -Differentiation::SD::differentiate(psi_sd, H_sd); 
      S_sd = 2.0 * Differentiation::SD::differentiate(psi_sd, C_sd); 



      const SymmetricTensor<2, dim, Differentiation::SD::Expression> 
        Q_t_sd_explicit = 
          (1.0 / (1.0 + delta_t_sd / tau_v_sd)) * 
          (Q_t1_sd + 
           (delta_t_sd / tau_v_sd * std::pow(det_F_sd, 2.0 / dim) * C_inv_sd)); 

      const Differentiation::SD::types::substitution_map 
        substitution_map_explicit = Differentiation::SD::make_substitution_map( 
          std::make_pair(Q_t_sd, Q_t_sd_explicit)); 


      BB_sd = symmetrize(Differentiation::SD::differentiate( 
        Differentiation::SD::substitute(B_sd, substitution_map_explicit), 
        H_sd)); 
      PP_sd = -Differentiation::SD::differentiate( 
        Differentiation::SD::substitute(S_sd, substitution_map_explicit), H_sd); 
      HH_sd = 
        2.0 * 
        Differentiation::SD::differentiate( 
          Differentiation::SD::substitute(S_sd, substitution_map_explicit), 
          C_sd); 



      optimizer.register_symbols( 
        Differentiation::SD::Utilities::extract_symbols( 
          make_substitution_map({}, {}, 0))); 


      optimizer.register_functions(psi_sd, B_sd, S_sd, BB_sd, PP_sd, HH_sd); 


      optimizer.optimize(); 
    } 


    template <int dim> 
    void Magnetoviscoelastic_Constitutive_Law_SD<dim>::update_internal_data( 
      const SymmetricTensor<2, dim> &C, 
      const Tensor<1, dim> &         H, 
      const DiscreteTime &           time) 
    { 


      const double delta_t = this->get_delta_t(time); 

      const double                  det_F = std::sqrt(determinant(C)); 
      const SymmetricTensor<2, dim> C_inv = invert(C); 
      AssertThrow(det_F > 0.0, 
                  ExcMessage("Volumetric Jacobian must be positive.")); 


      Q_t = (1.0 / (1.0 + delta_t / this->get_tau_v())) * 
            (Q_t1 + (delta_t / this->get_tau_v()) * std::pow(det_F, 2.0 / dim) * 
                      C_inv); 


      const auto substitution_map = make_substitution_map(C, H, delta_t); 


      optimizer.substitute(substitution_map); 
    } 


    template <int dim> 
    double Magnetoviscoelastic_Constitutive_Law_SD<dim>::get_psi() const 
    { 
      return optimizer.evaluate(psi_sd); 
    } 

    template <int dim> 
    Tensor<1, dim> Magnetoviscoelastic_Constitutive_Law_SD<dim>::get_B() const 
    { 
      return optimizer.evaluate(B_sd); 
    } 

    template <int dim> 
    SymmetricTensor<2, dim> 
    Magnetoviscoelastic_Constitutive_Law_SD<dim>::get_S() const 
    { 
      return optimizer.evaluate(S_sd); 
    } 

    template <int dim> 
    SymmetricTensor<2, dim> 
    Magnetoviscoelastic_Constitutive_Law_SD<dim>::get_DD() const 
    { 
      return optimizer.evaluate(BB_sd); 
    } 

    template <int dim> 
    Tensor<3, dim> Magnetoviscoelastic_Constitutive_Law_SD<dim>::get_PP() const 
    { 
      return optimizer.evaluate(PP_sd); 
    } 

    template <int dim> 
    SymmetricTensor<4, dim> 
    Magnetoviscoelastic_Constitutive_Law_SD<dim>::get_HH() const 
    { 
      return optimizer.evaluate(HH_sd); 
    } 


    template <int dim> 
    void Magnetoviscoelastic_Constitutive_Law_SD<dim>::update_end_of_timestep() 
    { 
      Q_t1 = Q_t; 
    } 



















    template <int dim> 
    class Magnetoelastic_Constitutive_Law final 
      : public Coupled_Magnetomechanical_Constitutive_Law_Base<dim> 
    { 
    public: 
      Magnetoelastic_Constitutive_Law( 
        const ConstitutiveParameters &constitutive_parameters); 

      virtual void update_internal_data(const SymmetricTensor<2, dim> &C, 
                                        const Tensor<1, dim> &         H, 
                                        const DiscreteTime &) override; 

      virtual double get_psi() const override; 

      virtual Tensor<1, dim> get_B() const override; 

      virtual SymmetricTensor<2, dim> get_S() const override; 

      virtual SymmetricTensor<2, dim> get_DD() const override; 

      virtual Tensor<3, dim> get_PP() const override; 

      virtual SymmetricTensor<4, dim> get_HH() const override; 

    private: 
      double                  psi; 
      Tensor<1, dim>          B; 
      SymmetricTensor<2, dim> S; 
      SymmetricTensor<2, dim> BB; 
      Tensor<3, dim>          PP; 
      SymmetricTensor<4, dim> HH; 
    }; 

    template <int dim> 
    Magnetoelastic_Constitutive_Law<dim>::Magnetoelastic_Constitutive_Law( 
      const ConstitutiveParameters &constitutive_parameters) 
      : Coupled_Magnetomechanical_Constitutive_Law_Base<dim>( 
          constitutive_parameters) 
      , psi(0.0) 
    {} 



    template <int dim> 
    void Magnetoelastic_Constitutive_Law<dim>::update_internal_data( 
      const SymmetricTensor<2, dim> &C, 
      const Tensor<1, dim> &         H, 
      const DiscreteTime &) 
    { 
      const double                  det_F = std::sqrt(determinant(C)); 
      const SymmetricTensor<2, dim> C_inv = invert(C); 
      AssertThrow(det_F > 0.0, 
                  ExcMessage("Volumetric Jacobian must be positive.")); 


      const double two_h_dot_h_div_h_sat_squ = 
        (2.0 * H * H) / (this->get_mu_e_h_sat() * this->get_mu_e_h_sat()); 
      const double tanh_two_h_dot_h_div_h_sat_squ = 
        std::tanh(two_h_dot_h_div_h_sat_squ); 

      const double f_mu_e = 
        1.0 + (this->get_mu_e_inf() / this->get_mu_e() - 1.0) * 
                tanh_two_h_dot_h_div_h_sat_squ; 


      const double dtanh_two_h_dot_h_div_h_sat_squ = 
        std::pow(1.0 / std::cosh(two_h_dot_h_div_h_sat_squ), 2.0); 
      const Tensor<1, dim> dtwo_h_dot_h_div_h_sat_squ_dH = 
        2.0 * 2.0 / (this->get_mu_e_h_sat() * this->get_mu_e_h_sat()) * H; 

      const Tensor<1, dim> df_mu_e_dH = 
        (this->get_mu_e_inf() / this->get_mu_e() - 1.0) * 
        (dtanh_two_h_dot_h_div_h_sat_squ * dtwo_h_dot_h_div_h_sat_squ_dH); 


      const double d2tanh_two_h_dot_h_div_h_sat_squ = 
        -2.0 * tanh_two_h_dot_h_div_h_sat_squ * dtanh_two_h_dot_h_div_h_sat_squ; 
      const SymmetricTensor<2, dim> d2two_h_dot_h_div_h_sat_squ_dH_dH = 
        2.0 * 2.0 / (this->get_mu_e_h_sat() * this->get_mu_e_h_sat()) * 
        Physics::Elasticity::StandardTensors<dim>::I; 

      const SymmetricTensor<2, dim> d2f_mu_e_dH_dH = 
        (this->get_mu_e_inf() / this->get_mu_e() - 1.0) * 
        (d2tanh_two_h_dot_h_div_h_sat_squ * 
           symmetrize(outer_product(dtwo_h_dot_h_div_h_sat_squ_dH, 
                                    dtwo_h_dot_h_div_h_sat_squ_dH)) + 
         dtanh_two_h_dot_h_div_h_sat_squ * d2two_h_dot_h_div_h_sat_squ_dH_dH); 


      const double         log_det_F         = std::log(det_F); 
      const double         tr_C              = trace(C); 
      const Tensor<1, dim> C_inv_dot_H       = C_inv * H; 
      const double         H_dot_C_inv_dot_H = H * C_inv_dot_H; 


      const SymmetricTensor<2, dim> d_tr_C_dC = 
        Physics::Elasticity::StandardTensors<dim>::I; 
      const SymmetricTensor<2, dim> ddet_F_dC     = 0.5 * det_F * C_inv; 
      const SymmetricTensor<2, dim> dlog_det_F_dC = 0.5 * C_inv; 

      const Tensor<1, dim> dH_dot_C_inv_dot_H_dH = 2.0 * C_inv_dot_H; 

      SymmetricTensor<4, dim> dC_inv_dC; 
      for (unsigned int A = 0; A < dim; ++A) 
        for (unsigned int B = A; B < dim; ++B) 
          for (unsigned int C = 0; C < dim; ++C) 
            for (unsigned int D = C; D < dim; ++D) 
              dC_inv_dC[A][B][C][D] -=               // 
                0.5 * (C_inv[A][C] * C_inv[B][D]     // 
                       + C_inv[A][D] * C_inv[B][C]); // 

      const SymmetricTensor<2, dim> dH_dot_C_inv_dot_H_dC = 
        -symmetrize(outer_product(C_inv_dot_H, C_inv_dot_H)); 


      const SymmetricTensor<4, dim> d2log_det_F_dC_dC = 0.5 * dC_inv_dC; 

      const SymmetricTensor<4, dim> d2det_F_dC_dC = 
        0.5 * (outer_product(C_inv, ddet_F_dC) + det_F * dC_inv_dC); 

      const SymmetricTensor<2, dim> d2H_dot_C_inv_dot_H_dH_dH = 2.0 * C_inv; 

      Tensor<3, dim> d2H_dot_C_inv_dot_H_dC_dH; 
      for (unsigned int A = 0; A < dim; ++A) 
        for (unsigned int B = 0; B < dim; ++B) 
          for (unsigned int C = 0; C < dim; ++C) 
            d2H_dot_C_inv_dot_H_dC_dH[A][B][C] -= 
              C_inv[A][C] * C_inv_dot_H[B] + // 
              C_inv_dot_H[A] * C_inv[B][C];  // 

      SymmetricTensor<4, dim> d2H_dot_C_inv_dot_H_dC_dC; 
      for (unsigned int A = 0; A < dim; ++A) 
        for (unsigned int B = A; B < dim; ++B) 
          for (unsigned int C = 0; C < dim; ++C) 
            for (unsigned int D = C; D < dim; ++D) 
              d2H_dot_C_inv_dot_H_dC_dC[A][B][C][D] += 
                0.5 * (C_inv_dot_H[A] * C_inv_dot_H[C] * C_inv[B][D] + 
                       C_inv_dot_H[A] * C_inv_dot_H[D] * C_inv[B][C] + 
                       C_inv_dot_H[B] * C_inv_dot_H[C] * C_inv[A][D] + 
                       C_inv_dot_H[B] * C_inv_dot_H[D] * C_inv[A][C]); 


      psi = 
        (0.5 * this->get_mu_e() * f_mu_e) * 
          (tr_C - dim - 2.0 * std::log(det_F)) + 
        this->get_lambda_e() * (std::log(det_F) * std::log(det_F)) - 
        (0.5 * this->get_mu_0() * this->get_mu_r()) * det_F * (H * C_inv * H); 


      B = -(0.5 * this->get_mu_e() * (tr_C - dim - 2.0 * log_det_F)) * 
            df_mu_e_dH // 
          + 0.5 * this->get_mu_0() * this->get_mu_r() * det_F * 
              dH_dot_C_inv_dot_H_dH; // 

      S = 2.0 * (0.5 * this->get_mu_e() * f_mu_e) *                        // 
            (d_tr_C_dC - 2.0 * dlog_det_F_dC)                              // 
          + 2.0 * this->get_lambda_e() * (2.0 * log_det_F * dlog_det_F_dC) // 
          - 2.0 * (0.5 * this->get_mu_0() * this->get_mu_r()) *            // 
              (H_dot_C_inv_dot_H * ddet_F_dC                               // 
               + det_F * dH_dot_C_inv_dot_H_dC);                           // 


      BB = -(0.5 * this->get_mu_e() * (tr_C - dim - 2.0 * log_det_F)) * // 
             d2f_mu_e_dH_dH                                             // 
           + 0.5 * this->get_mu_0() * this->get_mu_r() * det_F * 
               d2H_dot_C_inv_dot_H_dH_dH; // 

      PP = -2.0 * (0.5 * this->get_mu_e()) *                                  // 
             outer_product(Tensor<2, dim>(d_tr_C_dC - 2.0 * dlog_det_F_dC),   // 
                           df_mu_e_dH)                                        // 
           +                                                                  // 
           2.0 * (0.5 * this->get_mu_0() * this->get_mu_r()) *                // 
             (outer_product(Tensor<2, dim>(ddet_F_dC), dH_dot_C_inv_dot_H_dH) // 
              + det_F * d2H_dot_C_inv_dot_H_dC_dH);                           // 

      HH = 
        4.0 * (0.5 * this->get_mu_e() * f_mu_e) * (-2.0 * d2log_det_F_dC_dC) // 
        + 4.0 * this->get_lambda_e() *                                       // 
            (2.0 * outer_product(dlog_det_F_dC, dlog_det_F_dC)               // 
             + 2.0 * log_det_F * d2log_det_F_dC_dC)                          // 
        - 4.0 * (0.5 * this->get_mu_0() * this->get_mu_r()) *                // 
            (H_dot_C_inv_dot_H * d2det_F_dC_dC                               // 
             + outer_product(ddet_F_dC, dH_dot_C_inv_dot_H_dC)               // 
             + outer_product(dH_dot_C_inv_dot_H_dC, ddet_F_dC)               // 
             + det_F * d2H_dot_C_inv_dot_H_dC_dC);                           // 
    } 

    template <int dim> 
    double Magnetoelastic_Constitutive_Law<dim>::get_psi() const 
    { 
      return psi; 
    } 

    template <int dim> 
    Tensor<1, dim> Magnetoelastic_Constitutive_Law<dim>::get_B() const 
    { 
      return B; 
    } 

    template <int dim> 
    SymmetricTensor<2, dim> Magnetoelastic_Constitutive_Law<dim>::get_S() const 
    { 
      return S; 
    } 

    template <int dim> 
    SymmetricTensor<2, dim> Magnetoelastic_Constitutive_Law<dim>::get_DD() const 
    { 
      return BB; 
    } 

    template <int dim> 
    Tensor<3, dim> Magnetoelastic_Constitutive_Law<dim>::get_PP() const 
    { 
      return PP; 
    } 

    template <int dim> 
    SymmetricTensor<4, dim> Magnetoelastic_Constitutive_Law<dim>::get_HH() const 
    { 
      return HH; 
    } 









    template <int dim> 
    class Magnetoviscoelastic_Constitutive_Law final 
      : public Coupled_Magnetomechanical_Constitutive_Law_Base<dim> 
    { 
    public: 
      Magnetoviscoelastic_Constitutive_Law( 
        const ConstitutiveParameters &constitutive_parameters); 

      virtual void update_internal_data(const SymmetricTensor<2, dim> &C, 
                                        const Tensor<1, dim> &         H, 
                                        const DiscreteTime &time) override; 

      virtual double get_psi() const override; 

      virtual Tensor<1, dim> get_B() const override; 

      virtual SymmetricTensor<2, dim> get_S() const override; 

      virtual SymmetricTensor<2, dim> get_DD() const override; 

      virtual Tensor<3, dim> get_PP() const override; 

      virtual SymmetricTensor<4, dim> get_HH() const override; 

      virtual void update_end_of_timestep() override; 

    private: 
      SymmetricTensor<2, dim> Q_t; 
      SymmetricTensor<2, dim> Q_t1; 

      double                  psi; 
      Tensor<1, dim>          B; 
      SymmetricTensor<2, dim> S; 
      SymmetricTensor<2, dim> BB; 
      Tensor<3, dim>          PP; 
      SymmetricTensor<4, dim> HH; 


      mutable GeneralDataStorage cache; 


      void set_primary_variables(const SymmetricTensor<2, dim> &C, 
                                 const Tensor<1, dim> &         H) const; 

      void update_internal_variable(const DiscreteTime &time); 



      const Tensor<1, dim> &get_H() const; 

      const SymmetricTensor<2, dim> &get_C() const; 


      double get_two_h_dot_h_div_h_sat_squ(const double mu_h_sat) const; 

      double get_tanh_two_h_dot_h_div_h_sat_squ(const double mu_h_sat) const; 

      double get_f_mu(const double mu, 
                      const double mu_inf, 
                      const double mu_h_sat) const; 


      double get_dtanh_two_h_dot_h_div_h_sat_squ(const double mu_h_sat) const; 

      Tensor<1, dim> 
      get_dtwo_h_dot_h_div_h_sat_squ_dH(const double mu_h_sat) const; 

      Tensor<1, dim> get_df_mu_dH(const double mu, 
                                  const double mu_inf, 
                                  const double mu_h_sat) const; 


      double get_d2tanh_two_h_dot_h_div_h_sat_squ(const double mu_h_sat) const; 

      SymmetricTensor<2, dim> 
      get_d2two_h_dot_h_div_h_sat_squ_dH_dH(const double mu_h_sat) const; 

      SymmetricTensor<2, dim> get_d2f_mu_dH_dH(const double mu, 
                                               const double mu_inf, 
                                               const double mu_h_sat) const; 


      const double &get_det_F() const; 

      const SymmetricTensor<2, dim> &get_C_inv() const; 

      const double &get_log_det_F() const; 

 

      const Tensor<1, dim> &get_C_inv_dot_H() const; 

 


      const SymmetricTensor<4, dim> &get_dC_inv_dC() const; 

      const SymmetricTensor<2, dim> &get_d_tr_C_dC() const; 

      const SymmetricTensor<2, dim> &get_ddet_F_dC() const; 

      const SymmetricTensor<2, dim> &get_dlog_det_F_dC() const; 

 

 


      const SymmetricTensor<4, dim> & 
      get_dQ_t_dC(const DiscreteTime &time) const; 


      const SymmetricTensor<4, dim> &get_d2log_det_F_dC_dC() const; 

      const SymmetricTensor<4, dim> &get_d2det_F_dC_dC() const; 

 

      const Tensor<3, dim> &get_d2H_dot_C_inv_dot_H_dC_dH() const; 

      const SymmetricTensor<4, dim> &get_d2H_dot_C_inv_dot_H_dC_dC() const; 
    }; 

    template <int dim> 
    Magnetoviscoelastic_Constitutive_Law< 
      dim>::Magnetoviscoelastic_Constitutive_Law(const ConstitutiveParameters 
                                                   &constitutive_parameters) 
      : Coupled_Magnetomechanical_Constitutive_Law_Base<dim>( 
          constitutive_parameters) 
      , Q_t(Physics::Elasticity::StandardTensors<dim>::I) 
      , Q_t1(Physics::Elasticity::StandardTensors<dim>::I) 
      , psi(0.0) 
    {} 

    template <int dim> 
    void Magnetoviscoelastic_Constitutive_Law<dim>::update_internal_data( 
      const SymmetricTensor<2, dim> &C, 
      const Tensor<1, dim> &         H, 
      const DiscreteTime &           time) 
    { 


      set_primary_variables(C, H); 
      update_internal_variable(time); 


 
                                     this->get_mu_e_inf(), 
                                     this->get_mu_e_h_sat()); 

      const double f_mu_v = get_f_mu(this->get_mu_v(), 
                                     this->get_mu_v_inf(), 
                                     this->get_mu_v_h_sat()); 


      const Tensor<1, dim> df_mu_e_dH = get_df_mu_dH(this->get_mu_e(), 
                                                     this->get_mu_e_inf(), 
                                                     this->get_mu_e_h_sat()); 

      const Tensor<1, dim> df_mu_v_dH = get_df_mu_dH(this->get_mu_v(), 
                                                     this->get_mu_v_inf(), 
                                                     this->get_mu_v_h_sat()); 


      const SymmetricTensor<2, dim> d2f_mu_e_dH_dH = 
        get_d2f_mu_dH_dH(this->get_mu_e(), 
                         this->get_mu_e_inf(), 
                         this->get_mu_e_h_sat()); 

      const SymmetricTensor<2, dim> d2f_mu_v_dH_dH = 
        get_d2f_mu_dH_dH(this->get_mu_v(), 
                         this->get_mu_v_inf(), 
                         this->get_mu_v_h_sat()); 


      const double &                 det_F = get_det_F(); 
      const SymmetricTensor<2, dim> &C_inv = get_C_inv(); 

      const double &log_det_F         = get_log_det_F(); 
      const double &tr_C              = get_trace_C(); 
      const double &H_dot_C_inv_dot_H = get_H_dot_C_inv_dot_H(); 


      const SymmetricTensor<2, dim> &d_tr_C_dC     = get_d_tr_C_dC(); 
      const SymmetricTensor<2, dim> &ddet_F_dC     = get_ddet_F_dC(); 
      const SymmetricTensor<2, dim> &dlog_det_F_dC = get_dlog_det_F_dC(); 

      const SymmetricTensor<4, dim> &dQ_t_dC = get_dQ_t_dC(time); 

      const Tensor<1, dim> &dH_dot_C_inv_dot_H_dH = get_dH_dot_C_inv_dot_H_dH(); 

      const SymmetricTensor<2, dim> &dH_dot_C_inv_dot_H_dC = 
        get_dH_dot_C_inv_dot_H_dC(); 


      const SymmetricTensor<4, dim> &d2log_det_F_dC_dC = 
        get_d2log_det_F_dC_dC(); 

      const SymmetricTensor<4, dim> &d2det_F_dC_dC = get_d2det_F_dC_dC(); 

      const SymmetricTensor<2, dim> &d2H_dot_C_inv_dot_H_dH_dH = 
        get_d2H_dot_C_inv_dot_H_dH_dH(); 

      const Tensor<3, dim> &d2H_dot_C_inv_dot_H_dC_dH = 
        get_d2H_dot_C_inv_dot_H_dC_dH(); 

      const SymmetricTensor<4, dim> &d2H_dot_C_inv_dot_H_dC_dC = 
        get_d2H_dot_C_inv_dot_H_dC_dC(); 







      psi = (0.5 * this->get_mu_e() * f_mu_e) * 
              (tr_C - dim - 2.0 * std::log(det_F)) + 
            this->get_lambda_e() * (std::log(det_F) * std::log(det_F)); 
      psi += (0.5 * this->get_mu_v() * f_mu_v) * 
             (Q_t * (std::pow(det_F, -2.0 / dim) * C) - dim - 
              std::log(determinant(Q_t))); 
      psi -= 
        (0.5 * this->get_mu_0() * this->get_mu_r()) * det_F * (H * C_inv * H); 


      B = 
        -(0.5 * this->get_mu_e() * (tr_C - dim - 2.0 * log_det_F)) * df_mu_e_dH; 
      B -= (0.5 * this->get_mu_v()) * 
           (Q_t * (std::pow(det_F, -2.0 / dim) * C) - dim - 
            std::log(determinant(Q_t))) * 
           df_mu_v_dH; 
      B += 0.5 * this->get_mu_0() * this->get_mu_r() * det_F * 
           dH_dot_C_inv_dot_H_dH; 

      S = 2.0 * (0.5 * this->get_mu_e() * f_mu_e) *                         // 
            (d_tr_C_dC - 2.0 * dlog_det_F_dC)                               // 
          + 2.0 * this->get_lambda_e() * (2.0 * log_det_F * dlog_det_F_dC); // 
      S += 2.0 * (0.5 * this->get_mu_v() * f_mu_v) * 
           ((Q_t * C) * 
              ((-2.0 / dim) * std::pow(det_F, -2.0 / dim - 1.0) * ddet_F_dC) + 
            std::pow(det_F, -2.0 / dim) * Q_t);                // dC/dC = II 
      S -= 2.0 * (0.5 * this->get_mu_0() * this->get_mu_r()) * // 
           (H_dot_C_inv_dot_H * ddet_F_dC                      // 
            + det_F * dH_dot_C_inv_dot_H_dC);                  // 


      BB = -(0.5 * this->get_mu_e() * (tr_C - dim - 2.0 * log_det_F)) * 
           d2f_mu_e_dH_dH; 
      BB -= (0.5 * this->get_mu_v()) * 
            (Q_t * (std::pow(det_F, -2.0 / dim) * C) - dim - 
             std::log(determinant(Q_t))) * 
            d2f_mu_v_dH_dH; 
      BB += 0.5 * this->get_mu_0() * this->get_mu_r() * det_F * 
            d2H_dot_C_inv_dot_H_dH_dH; 

      PP = -2.0 * (0.5 * this->get_mu_e()) * 
           outer_product(Tensor<2, dim>(d_tr_C_dC - 2.0 * dlog_det_F_dC), 
                         df_mu_e_dH); 
      PP -= 2.0 * (0.5 * this->get_mu_v()) * 
            outer_product(Tensor<2, dim>((Q_t * C) * 
                                           ((-2.0 / dim) * 
                                            std::pow(det_F, -2.0 / dim - 1.0) * 
                                            ddet_F_dC) + 
                                         std::pow(det_F, -2.0 / dim) * Q_t), 
                          df_mu_v_dH); 
      PP += 2.0 * (0.5 * this->get_mu_0() * this->get_mu_r()) * 
            (outer_product(Tensor<2, dim>(ddet_F_dC), dH_dot_C_inv_dot_H_dH) + 
             det_F * d2H_dot_C_inv_dot_H_dC_dH); 

      HH = 
        4.0 * (0.5 * this->get_mu_e() * f_mu_e) * (-2.0 * d2log_det_F_dC_dC) // 
        + 4.0 * this->get_lambda_e() *                                       // 
            (2.0 * outer_product(dlog_det_F_dC, dlog_det_F_dC)               // 
             + 2.0 * log_det_F * d2log_det_F_dC_dC);                         // 
      HH += 4.0 * (0.5 * this->get_mu_v() * f_mu_v) * 
            (outer_product((-2.0 / dim) * std::pow(det_F, -2.0 / dim - 1.0) * 
                             ddet_F_dC, 
                           C * dQ_t_dC + Q_t) + 
             (Q_t * C) * 
               (outer_product(ddet_F_dC, 
                              (-2.0 / dim) * (-2.0 / dim - 1.0) * 
                                std::pow(det_F, -2.0 / dim - 2.0) * ddet_F_dC) + 
                ((-2.0 / dim) * std::pow(det_F, -2.0 / dim - 1.0) * 
                 d2det_F_dC_dC)) + 
             outer_product(Q_t, 
                           (-2.0 / dim) * std::pow(det_F, -2.0 / dim - 1.0) * 
                             ddet_F_dC) + 
             std::pow(det_F, -2.0 / dim) * dQ_t_dC); 
      HH -= 4.0 * (0.5 * this->get_mu_0() * this->get_mu_r()) * // 
            (H_dot_C_inv_dot_H * d2det_F_dC_dC                  // 
             + outer_product(ddet_F_dC, dH_dot_C_inv_dot_H_dC)  // 
             + outer_product(dH_dot_C_inv_dot_H_dC, ddet_F_dC)  // 
             + det_F * d2H_dot_C_inv_dot_H_dC_dC);              // 


      cache.reset(); 
    } 

    template <int dim> 
    double Magnetoviscoelastic_Constitutive_Law<dim>::get_psi() const 
    { 
      return psi; 
    } 

    template <int dim> 
    Tensor<1, dim> Magnetoviscoelastic_Constitutive_Law<dim>::get_B() const 
    { 
      return B; 
    } 

    template <int dim> 
    SymmetricTensor<2, dim> 
    Magnetoviscoelastic_Constitutive_Law<dim>::get_S() const 
    { 
      return S; 
    } 

    template <int dim> 
    SymmetricTensor<2, dim> 
    Magnetoviscoelastic_Constitutive_Law<dim>::get_DD() const 
    { 
      return BB; 
    } 

    template <int dim> 
    Tensor<3, dim> Magnetoviscoelastic_Constitutive_Law<dim>::get_PP() const 
    { 
      return PP; 
    } 

    template <int dim> 
    SymmetricTensor<4, dim> 
    Magnetoviscoelastic_Constitutive_Law<dim>::get_HH() const 
    { 
      return HH; 
    } 

    template <int dim> 
    void Magnetoviscoelastic_Constitutive_Law<dim>::update_end_of_timestep() 
    { 
      Q_t1 = Q_t; 
    } 

    template <int dim> 
    void Magnetoviscoelastic_Constitutive_Law<dim>::update_internal_variable( 
      const DiscreteTime &time) 
    { 
      const double delta_t = this->get_delta_t(time); 

      Q_t = (1.0 / (1.0 + delta_t / this->get_tau_v())) * 
            (Q_t1 + (delta_t / this->get_tau_v()) * 
                      std::pow(get_det_F(), 2.0 / dim) * get_C_inv()); 
    } 


    template <int dim> 
    double 
    Magnetoviscoelastic_Constitutive_Law<dim>::get_two_h_dot_h_div_h_sat_squ( 
      const double mu_h_sat) const 
    { 
      const Tensor<1, dim> &H = get_H(); 
      return (2.0 * H * H) / (mu_h_sat * mu_h_sat); 
    } 

    template <int dim> 
    double Magnetoviscoelastic_Constitutive_Law< 
      dim>::get_tanh_two_h_dot_h_div_h_sat_squ(const double mu_h_sat) const 
    { 
      return std::tanh(get_two_h_dot_h_div_h_sat_squ(mu_h_sat)); 
    } 


 
    double Magnetoviscoelastic_Constitutive_Law<dim>::get_f_mu( 
      const double mu, 
      const double mu_inf, 
      const double mu_h_sat) const 
    { 
      return 1.0 + 
             (mu_inf / mu - 1.0) * get_tanh_two_h_dot_h_div_h_sat_squ(mu_h_sat); 
    } 


    template <int dim> 
    double Magnetoviscoelastic_Constitutive_Law< 
      dim>::get_dtanh_two_h_dot_h_div_h_sat_squ(const double mu_h_sat) const 
    { 
      return std::pow(1.0 / std::cosh(get_two_h_dot_h_div_h_sat_squ(mu_h_sat)), 
                      2.0); 
    } 

    template <int dim> 
    Tensor<1, dim> Magnetoviscoelastic_Constitutive_Law< 
      dim>::get_dtwo_h_dot_h_div_h_sat_squ_dH(const double mu_h_sat) const 
    { 
      return 2.0 * 2.0 / (mu_h_sat * mu_h_sat) * get_H(); 
    } 

    template <int dim> 
    Tensor<1, dim> Magnetoviscoelastic_Constitutive_Law<dim>::get_df_mu_dH( 
      const double mu, 
      const double mu_inf, 
      const double mu_h_sat) const 
  template <int dim> 
      return (mu_inf / mu - 1.0) * 
             (get_dtanh_two_h_dot_h_div_h_sat_squ(mu_h_sat) * 
              get_dtwo_h_dot_h_div_h_sat_squ_dH(mu_h_sat)); 
    } 

    template <int dim> 
    double Magnetoviscoelastic_Constitutive_Law< 
      dim>::get_d2tanh_two_h_dot_h_div_h_sat_squ(const double mu_h_sat) const 
    { 
      return -2.0 * get_tanh_two_h_dot_h_div_h_sat_squ(mu_h_sat) * 
             get_dtanh_two_h_dot_h_div_h_sat_squ(mu_h_sat); 
    } 

    template <int dim> 
    SymmetricTensor<2, dim> Magnetoviscoelastic_Constitutive_Law< 
      dim>::get_d2two_h_dot_h_div_h_sat_squ_dH_dH(const double mu_h_sat) const 
    { 
      return 2.0 * 2.0 / (mu_h_sat * mu_h_sat) * 
             Physics::Elasticity::StandardTensors<dim>::I; 
    } 

    template <int dim> 
    SymmetricTensor<2, dim> 
    Magnetoviscoelastic_Constitutive_Law<dim>::get_d2f_mu_dH_dH( 
      const double mu, 
      const double mu_inf, 
      const double mu_h_sat) const 
    { 
      return (mu_inf / mu - 1.0) * 
             (get_d2tanh_two_h_dot_h_div_h_sat_squ(mu_h_sat) * 
                symmetrize( 
                  outer_product(get_dtwo_h_dot_h_div_h_sat_squ_dH(mu_h_sat), 
                                get_dtwo_h_dot_h_div_h_sat_squ_dH(mu_h_sat))) + 
              get_dtanh_two_h_dot_h_div_h_sat_squ(mu_h_sat) * 
                get_d2two_h_dot_h_div_h_sat_squ_dH_dH(mu_h_sat)); 
    } 


    template <int dim> 
    void Magnetoviscoelastic_Constitutive_Law<dim>::set_primary_variables( 
      const SymmetricTensor<2, dim> &C, 
      const Tensor<1, dim> &         H) const 
    { 


      const std::string name_H("H"); 
      Assert(!cache.stores_object_with_name(name_H), 
             ExcMessage( 
               "The primary variable has already been added to the cache.")); 
      cache.add_unique_copy(name_H, H); 


      const std::string name_C("C"); 
      Assert(!cache.stores_object_with_name(name_C), 
             ExcMessage( 
               "The primary variable has already been added to the cache.")); 
      cache.add_unique_copy(name_C, C); 
    } 


    template <int dim> 
    const Tensor<1, dim> & 
    Magnetoviscoelastic_Constitutive_Law<dim>::get_H() const 
    { 
      const std::string name("H"); 
      Assert(cache.stores_object_with_name(name), 
             ExcMessage("Primary variables must be added to the cache.")); 
      return cache.template get_object_with_name<Tensor<1, dim>>(name); 
    } 

    template <int dim> 
    const SymmetricTensor<2, dim> & 
    Magnetoviscoelastic_Constitutive_Law<dim>::get_C() const 
    { 
      const std::string name("C"); 
      Assert(cache.stores_object_with_name(name), 
             ExcMessage("Primary variables must be added to the cache.")); 
      return cache.template get_object_with_name<SymmetricTensor<2, dim>>(name); 
    } 



    template <int dim> 
    const double &Magnetoviscoelastic_Constitutive_Law<dim>::get_det_F() const 
    { 
      const std::string name("det_F"); 
      if (cache.stores_object_with_name(name) == false) 
        { 
          const double det_F = std::sqrt(determinant(get_C())); 
          AssertThrow(det_F > 0.0, 
                      ExcMessage("Volumetric Jacobian must be positive.")); 
          cache.add_unique_copy(name, det_F); 
        } 

      return cache.template get_object_with_name<double>(name); 
    } 

    template <int dim> 
    const SymmetricTensor<2, dim> & 
    Magnetoviscoelastic_Constitutive_Law<dim>::get_C_inv() const 
    { 
      const std::string name("C_inv"); 
      if (cache.stores_object_with_name(name) == false) 
        { 
          cache.add_unique_copy(name, invert(get_C())); 
        } 

      return cache.template get_object_with_name<SymmetricTensor<2, dim>>(name); 
    } 

    template <int dim> 
    const double & 
    Magnetoviscoelastic_Constitutive_Law<dim>::get_log_det_F() const 
    { 
      const std::string name("log(det_F)"); 
      if (cache.stores_object_with_name(name) == false) 
        cache.add_unique_copy(name, std::log(get_det_F())); 

      return cache.template get_object_with_name<double>(name); 
    } 

    template <int dim> 
    const double &Magnetoviscoelastic_Constitutive_Law<dim>::get_trace_C() const 
    { 
      const std::string name("trace(C)"); 
      if (cache.stores_object_with_name(name) == false) 
        cache.add_unique_copy(name, trace(get_C())); 

      return cache.template get_object_with_name<double>(name); 
    } 

    template <int dim> 
    const Tensor<1, dim> & 
    Magnetoviscoelastic_Constitutive_Law<dim>::get_C_inv_dot_H() const 
    { 
      const std::string name("C_inv_dot_H"); 
      if (cache.stores_object_with_name(name) == false) 
        cache.add_unique_copy(name, get_C_inv() * get_H()); 

      return cache.template get_object_with_name<Tensor<1, dim>>(name); 
    } 

    template <int dim> 
    const double & 
    Magnetoviscoelastic_Constitutive_Law<dim>::get_H_dot_C_inv_dot_H() const 
    { 
      const std::string name("H_dot_C_inv_dot_H"); 
      if (cache.stores_object_with_name(name) == false) 
        cache.add_unique_copy(name, get_H() * get_C_inv_dot_H()); 

      return cache.template get_object_with_name<double>(name); 
    } 

    template <int dim> 
    const SymmetricTensor<4, dim> & 
    Magnetoviscoelastic_Constitutive_Law<dim>::get_dQ_t_dC( 
      const DiscreteTime &time) const 
    { 
      const std::string name("dQ_t_dC"); 
      if (cache.stores_object_with_name(name) == false) 
        { 
          const double  delta_t = this->get_delta_t(time); 
          const double &det_F   = get_det_F(); 

          const SymmetricTensor<4, dim> dQ_t_dC = 
            (1.0 / (1.0 + delta_t / this->get_tau_v())) * 
            (delta_t / this->get_tau_v()) * 
            ((2.0 / dim) * std::pow(det_F, 2.0 / dim - 1.0) * 
               outer_product(get_C_inv(), get_ddet_F_dC()) + 
             std::pow(det_F, 2.0 / dim) * get_dC_inv_dC()); 

          cache.add_unique_copy(name, dQ_t_dC); 
        } 

      return cache.template get_object_with_name<SymmetricTensor<4, dim>>(name); 
    } 

    template <int dim> 
    const SymmetricTensor<4, dim> & 
    Magnetoviscoelastic_Constitutive_Law<dim>::get_dC_inv_dC() const 
    { 
      const std::string name("dC_inv_dC"); 
      if (cache.stores_object_with_name(name) == false) 
        { 
          const SymmetricTensor<2, dim> &C_inv = get_C_inv(); 
          SymmetricTensor<4, dim>        dC_inv_dC; 

          for (unsigned int A = 0; A < dim; ++A) 
            for (unsigned int B = A; B < dim; ++B) 
              for (unsigned int C = 0; C < dim; ++C) 
                for (unsigned int D = C; D < dim; ++D) 
                  dC_inv_dC[A][B][C][D] -=               // 
                    0.5 * (C_inv[A][C] * C_inv[B][D]     // 
                           + C_inv[A][D] * C_inv[B][C]); // 

          cache.add_unique_copy(name, dC_inv_dC); 
        } 

      return cache.template get_object_with_name<SymmetricTensor<4, dim>>(name); 
    } 

    template <int dim> 
    const SymmetricTensor<2, dim> & 
    Magnetoviscoelastic_Constitutive_Law<dim>::get_d_tr_C_dC() const 
    { 
      const std::string name("d_tr_C_dC"); 
      if (cache.stores_object_with_name(name) == false) 
        cache.add_unique_copy(name, 
                              Physics::Elasticity::StandardTensors<dim>::I); 

      return cache.template get_object_with_name<SymmetricTensor<2, dim>>(name); 
    } 

    template <int dim> 
    const SymmetricTensor<2, dim> & 
    Magnetoviscoelastic_Constitutive_Law<dim>::get_ddet_F_dC() const 
    { 
      const std::string name("ddet_F_dC"); 
      if (cache.stores_object_with_name(name) == false) 
        cache.add_unique_copy(name, 0.5 * get_det_F() * get_C_inv()); 

      return cache.template get_object_with_name<SymmetricTensor<2, dim>>(name); 
    } 

    template <int dim> 
    const SymmetricTensor<2, dim> & 
    Magnetoviscoelastic_Constitutive_Law<dim>::get_dlog_det_F_dC() const 
    { 
      const std::string name("dlog_det_F_dC"); 
      if (cache.stores_object_with_name(name) == false) 
        cache.add_unique_copy(name, 0.5 * get_C_inv()); 

      return cache.template get_object_with_name<SymmetricTensor<2, dim>>(name); 
    } 

    template <int dim> 
    const Tensor<1, dim> & 
    Magnetoviscoelastic_Constitutive_Law<dim>::get_dH_dot_C_inv_dot_H_dH() const 
    { 
      const std::string name("dH_dot_C_inv_dot_H_dH"); 
      if (cache.stores_object_with_name(name) == false) 
        cache.add_unique_copy(name, 2.0 * get_C_inv_dot_H()); 

      return cache.template get_object_with_name<Tensor<1, dim>>(name); 
    } 

    template <int dim> 
    const SymmetricTensor<2, dim> & 
    Magnetoviscoelastic_Constitutive_Law<dim>::get_dH_dot_C_inv_dot_H_dC() const 
    { 
      const std::string name("dH_dot_C_inv_dot_H_dC"); 
      if (cache.stores_object_with_name(name) == false) 
        { 
          const Tensor<1, dim> C_inv_dot_H = get_C_inv_dot_H(); 
          cache.add_unique_copy( 
            name, -symmetrize(outer_product(C_inv_dot_H, C_inv_dot_H))); 
        } 

      return cache.template get_object_with_name<SymmetricTensor<2, dim>>(name); 
    } 

    template <int dim> 
    const SymmetricTensor<4, dim> & 
    Magnetoviscoelastic_Constitutive_Law<dim>::get_d2log_det_F_dC_dC() const 
    { 
      const std::string name("d2log_det_F_dC_dC"); 
      if (cache.stores_object_with_name(name) == false) 
        cache.add_unique_copy(name, 0.5 * get_dC_inv_dC()); 

      return cache.template get_object_with_name<SymmetricTensor<4, dim>>(name); 
    } 

    template <int dim> 
    const SymmetricTensor<4, dim> & 
    Magnetoviscoelastic_Constitutive_Law<dim>::get_d2det_F_dC_dC() const 
    { 
      const std::string name("d2det_F_dC_dC"); 
      if (cache.stores_object_with_name(name) == false) 
        cache.add_unique_copy(name, 
                              0.5 * 
                                (outer_product(get_C_inv(), get_ddet_F_dC()) + 
                                 get_det_F() * get_dC_inv_dC())); 

      return cache.template get_object_with_name<SymmetricTensor<4, dim>>(name); 
    } 

    template <int dim> 
    const SymmetricTensor<2, dim> & 
    Magnetoviscoelastic_Constitutive_Law<dim>::get_d2H_dot_C_inv_dot_H_dH_dH() 
      const 
    { 
      const std::string name("d2H_dot_C_inv_dot_H_dH_dH"); 
      if (cache.stores_object_with_name(name) == false) 
        cache.add_unique_copy(name, 2.0 * get_C_inv()); 

      return cache.template get_object_with_name<SymmetricTensor<2, dim>>(name); 
    } 

    template <int dim> 
    const Tensor<3, dim> & 
    Magnetoviscoelastic_Constitutive_Law<dim>::get_d2H_dot_C_inv_dot_H_dC_dH() 
      const 
    { 
      const std::string name("d2H_dot_C_inv_dot_H_dC_dH"); 
      if (cache.stores_object_with_name(name) == false) 
        { 
          const Tensor<1, dim> &         C_inv_dot_H = get_C_inv_dot_H(); 
          const SymmetricTensor<2, dim> &C_inv       = get_C_inv(); 

          Tensor<3, dim> d2H_dot_C_inv_dot_H_dC_dH; 
          for (unsigned int A = 0; A < dim; ++A) 
            for (unsigned int B = 0; B < dim; ++B) 
              for (unsigned int C = 0; C < dim; ++C) 
                d2H_dot_C_inv_dot_H_dC_dH[A][B][C] -= 
                  C_inv[A][C] * C_inv_dot_H[B] + // 
                  C_inv_dot_H[A] * C_inv[B][C];  // 

          cache.add_unique_copy(name, d2H_dot_C_inv_dot_H_dC_dH); 
        } 

      return cache.template get_object_with_name<Tensor<3, dim>>(name); 
    } 

    template <int dim> 
    const SymmetricTensor<4, dim> & 
    Magnetoviscoelastic_Constitutive_Law<dim>::get_d2H_dot_C_inv_dot_H_dC_dC() 
      const 
    { 
      const std::string name("d2H_dot_C_inv_dot_H_dC_dC"); 
      if (cache.stores_object_with_name(name) == false) 
        { 
          const Tensor<1, dim> &         C_inv_dot_H = get_C_inv_dot_H(); 
          const SymmetricTensor<2, dim> &C_inv       = get_C_inv(); 

          SymmetricTensor<4, dim> d2H_dot_C_inv_dot_H_dC_dC; 
          for (unsigned int A = 0; A < dim; ++A) 
            for (unsigned int B = A; B < dim; ++B) 
              for (unsigned int C = 0; C < dim; ++C) 
                for (unsigned int D = C; D < dim; ++D) 
                  d2H_dot_C_inv_dot_H_dC_dC[A][B][C][D] += 
                    0.5 * (C_inv_dot_H[A] * C_inv_dot_H[C] * C_inv[B][D] + 
                           C_inv_dot_H[A] * C_inv_dot_H[D] * C_inv[B][C] + 
                           C_inv_dot_H[B] * C_inv_dot_H[C] * C_inv[A][D] + 
                           C_inv_dot_H[B] * C_inv_dot_H[D] * C_inv[A][C]); 

          cache.add_unique_copy(name, d2H_dot_C_inv_dot_H_dC_dC); 
        } 

      return cache.template get_object_with_name<SymmetricTensor<4, dim>>(name); 
    } 


    class RheologicalExperimentParameters : public ParameterAcceptor 
    { 
    public: 
      RheologicalExperimentParameters(); 


      double sample_radius = 0.01; 
      double sample_height = 0.001; 





      double lambda_2 = 0.95; 
      double gamma_12 = 0.05; 
      double H_2      = 60.0e3; 





      double       frequency         = 1.0 / (2.0 * numbers::PI); 
      unsigned int n_cycles          = 5; 
      unsigned int n_steps_per_cycle = 2500; 


      bool        output_data_to_file = true; 
      std::string output_filename_rd = 
        "experimental_results-rate_dependent.csv"; 
      std::string output_filename_ri = 
        "experimental_results-rate_independent.csv"; 


      double start_time() const; 

      double end_time() const; 

      double delta_t() const; 


      Tensor<1, 3> get_H(const double time) const; 

      Tensor<2, 3> get_F(const double time) const; 


      bool print_status(const int step_number) const; 

      bool initialized = false; 
    }; 

    RheologicalExperimentParameters::RheologicalExperimentParameters() 
      : ParameterAcceptor("/Coupled Constitutive Laws/Rheological Experiment/") 
    { 
      add_parameter("Experimental sample radius", sample_radius); 
      add_parameter("Experimental sample radius", sample_height); 

      add_parameter("Axial stretch", lambda_2); 
      add_parameter("Shear strain amplitude", gamma_12); 
      add_parameter("Axial magnetic field strength", H_2); 

      add_parameter("Frequency", frequency); 
      add_parameter("Number of loading cycles", n_cycles); 
      add_parameter("Discretisation for each cycle", n_steps_per_cycle); 

      add_parameter("Output experimental results to file", output_data_to_file); 
      add_parameter("Output file name (rate dependent constitutive law)", 
                    output_filename_rd); 
      add_parameter("Output file name (rate independent constitutive law)", 
                    output_filename_ri); 

      parse_parameters_call_back.connect([&]() -> void { initialized = true; }); 
    } 

    double RheologicalExperimentParameters::start_time() const 
    { 
      return 0.0; 
    } 

    double RheologicalExperimentParameters::end_time() const 
    { 
      return n_cycles / frequency; 
    } 

    double RheologicalExperimentParameters::delta_t() const 
    { 
      return (end_time() - start_time()) / (n_steps_per_cycle * n_cycles); 
    } 

    bool 
    RheologicalExperimentParameters::print_status(const int step_number) const 
    { 
      return (step_number % (n_cycles * n_steps_per_cycle / 100)) == 0; 
    } 


    Tensor<1, 3> RheologicalExperimentParameters::get_H(const double) const 
    { 
      return Tensor<1, 3>({0.0, 0.0, H_2}); 
    } 


    Tensor<2, 3> RheologicalExperimentParameters::get_F(const double time) const 
    { 
      AssertThrow((sample_radius > 0.0 && sample_height > 0.0), 
                  ExcMessage("Non-physical sample dimensions")); 
      AssertThrow(lambda_2 > 0.0, 
                  ExcMessage("Non-physical applied axial stretch")); 

      const double sqrt_lambda_2     = std::sqrt(lambda_2); 
      const double inv_sqrt_lambda_2 = 1.0 / sqrt_lambda_2; 

      const double alpha_max = 
        std::atan(std::tan(gamma_12) * sample_height / 
                  sample_radius); // Small strain approximation 
      const double A       = sample_radius * alpha_max; 
      const double w       = 2.0 * numbers::PI * frequency; // in rad /s 
      const double gamma_t = A * std::sin(w * time); 
      const double tau_t = 
        gamma_t / 
        (sample_radius * sample_height); // Torsion angle per unit length 
      const double alpha_t = tau_t * lambda_2 * sample_height; 

      Tensor<2, 3> F; 
      F[0][0] = inv_sqrt_lambda_2 * std::cos(alpha_t); 
      F[0][1] = -inv_sqrt_lambda_2 * std::sin(alpha_t); 
      F[0][2] = -tau_t * sample_radius * sqrt_lambda_2 * std::sin(alpha_t); 
      F[1][0] = inv_sqrt_lambda_2 * std::sin(alpha_t); 
      F[1][1] = inv_sqrt_lambda_2 * std::cos(alpha_t); 
      F[1][2] = tau_t * sample_radius * sqrt_lambda_2 * std::cos(alpha_t); 
      F[2][0] = 0.0; 
      F[2][1] = 0.0; 
      F[2][2] = lambda_2; 

      AssertThrow((F[0][0] > 0) && (F[1][1] > 0) && (F[2][2] > 0), 
                  ExcMessage("Non-physical deformation gradient component.")); 
      AssertThrow(std::abs(determinant(F) - 1.0) < 1e-6, 
                  ExcMessage("Volumetric Jacobian is not equal to unity.")); 

      return F; 
    } 


    template <int dim> 
    void run_rheological_experiment( 
      const RheologicalExperimentParameters &experimental_parameters, 
      Coupled_Magnetomechanical_Constitutive_Law_Base<dim> 
        &material_hand_calculated, 
      Coupled_Magnetomechanical_Constitutive_Law_Base<dim> 
        &               material_assisted_computation, 
      TimerOutput &     timer, 
      const std::string filename) 
    { 


      const auto check_material_class_results = 
        []( 
          const Coupled_Magnetomechanical_Constitutive_Law_Base<dim> &to_verify, 
          const Coupled_Magnetomechanical_Constitutive_Law_Base<dim> &blessed, 
          const double tol = 1e-6) { 
          (void)to_verify; 
          (void)blessed; 
          (void)tol; 

          Assert(std::abs(blessed.get_psi() - to_verify.get_psi()) < tol, 
                 ExcMessage("No match for psi. Error: " + 
                            Utilities::to_string(std::abs( 
                              blessed.get_psi() - to_verify.get_psi())))); 

          Assert((blessed.get_B() - to_verify.get_B()).norm() < tol, 
                 ExcMessage("No match for B. Error: " + 
                            Utilities::to_string( 
                              (blessed.get_B() - to_verify.get_B()).norm()))); 
          Assert((blessed.get_S() - to_verify.get_S()).norm() < tol, 
                 ExcMessage("No match for S. Error: " + 
                            Utilities::to_string( 
                              (blessed.get_S() - to_verify.get_S()).norm()))); 

          Assert((blessed.get_DD() - to_verify.get_DD()).norm() < tol, 
                 ExcMessage("No match for BB. Error: " + 
                            Utilities::to_string( 
                              (blessed.get_DD() - to_verify.get_DD()).norm()))); 
          Assert((blessed.get_PP() - to_verify.get_PP()).norm() < tol, 
                 ExcMessage("No match for PP. Error: " + 
                            Utilities::to_string( 
                              (blessed.get_PP() - to_verify.get_PP()).norm()))); 
          Assert((blessed.get_HH() - to_verify.get_HH()).norm() < tol, 
                 ExcMessage("No match for HH. Error: " + 
                            Utilities::to_string( 
                              (blessed.get_HH() - to_verify.get_HH()).norm()))); 
        }; 


      std::ostringstream stream; 
      stream 
        << "Time;Axial magnetic field strength [A/m];Axial magnetic induction [T];Shear strain [%];Shear stress [Pa]\n"; 


      for (DiscreteTime time(experimental_parameters.start_time(), 
                             experimental_parameters.end_time() + 
                               experimental_parameters.delta_t(), 
                             experimental_parameters.delta_t()); 
           time.is_at_end() == false; 
           time.advance_time()) 
        { 
          if (experimental_parameters.print_status(time.get_step_number())) 
            std::cout << "Timestep = " << time.get_step_number() 
                      << " @ time = " << time.get_current_time() << "s." 
                      << std::endl; 


          const Tensor<1, dim> H = 
            experimental_parameters.get_H(time.get_current_time()); 
          const Tensor<2, dim> F = 
            experimental_parameters.get_F(time.get_current_time()); 
          const SymmetricTensor<2, dim> C = 
            Physics::Elasticity::Kinematics::C(F); 


          { 
            TimerOutput::Scope timer_section(timer, "Hand calculated"); 
            material_hand_calculated.update_internal_data(C, H, time); 
            material_hand_calculated.update_end_of_timestep(); 
          } 

          { 
            TimerOutput::Scope timer_section(timer, "Assisted computation"); 
            material_assisted_computation.update_internal_data(C, H, time); 
            material_assisted_computation.update_end_of_timestep(); 
          } 


          check_material_class_results(material_hand_calculated, 
                                       material_assisted_computation); 

          if (experimental_parameters.output_data_to_file) 
            { 


              const Tensor<1, dim> h = 
                Physics::Transformations::Covariant::push_forward(H, F); 
              const Tensor<1, dim> b = 
                Physics::Transformations::Piola::push_forward( 
                  material_hand_calculated.get_B(), F); 
              const SymmetricTensor<2, dim> sigma = 
                Physics::Transformations::Piola::push_forward( 
                  material_hand_calculated.get_S(), F); 
              stream << time.get_current_time() << ";" << h[2] << ";" << b[2] 
                     << ";" << F[1][2] * 100.0 << ";" << sigma[1][2] << "\n"; 
            } 
        } 


      if (experimental_parameters.output_data_to_file) 
        { 
          std::ofstream output(filename); 
          output << stream.str(); 
        } 
    } 


    void run(int argc, char *argv[]) 
    { 
      using namespace dealii; 

      constexpr unsigned int dim = 3; 

      const ConstitutiveParameters          constitutive_parameters; 
      const RheologicalExperimentParameters experimental_parameters; 

      std::string parameter_file; 
      if (argc > 1) 
        parameter_file = argv[1]; 
      else 
        parameter_file = "parameters.prm"; 
      ParameterAcceptor::initialize(parameter_file, "used_parameters.prm"); 


      { 
        TimerOutput timer(std::cout, 
                          TimerOutput::summary, 
                          TimerOutput::wall_times); 
        std::cout 
          << "Coupled magnetoelastic constitutive law using automatic differentiation." 
          << std::endl; 

        constexpr Differentiation::AD::NumberTypes ADTypeCode = 
          Differentiation::AD::NumberTypes::sacado_dfad_dfad; 

        Magnetoelastic_Constitutive_Law<dim> material(constitutive_parameters); 
        Magnetoelastic_Constitutive_Law_AD<dim, ADTypeCode> material_ad( 
          constitutive_parameters); 

        run_rheological_experiment(experimental_parameters, 
                                   material, 
                                   material_ad, 
                                   timer, 
                                   experimental_parameters.output_filename_ri); 

        std::cout << "... all calculations are correct!" << std::endl; 
      } 


      { 
        TimerOutput timer(std::cout, 
                          TimerOutput::summary, 
                          TimerOutput::wall_times); 
        std::cout 
          << "Coupled magneto-viscoelastic constitutive law using symbolic differentiation." 
          << std::endl; 

#ifdef DEAL_II_SYMENGINE_WITH_LLVM 
        std::cout << "Using LLVM optimizer." << std::endl; 
        constexpr Differentiation::SD::OptimizerType optimizer_type = 
          Differentiation::SD::OptimizerType::llvm; 
        constexpr Differentiation::SD::OptimizationFlags optimization_flags = 
          Differentiation::SD::OptimizationFlags::optimize_all; 
#else 
        std::cout << "Using lambda optimizer." << std::endl; 
        constexpr Differentiation::SD::OptimizerType optimizer_type = 
          Differentiation::SD::OptimizerType::lambda; 
        constexpr Differentiation::SD::OptimizationFlags optimization_flags = 
          Differentiation::SD::OptimizationFlags::optimize_cse; 
#endif 

        Magnetoviscoelastic_Constitutive_Law<dim> material( 
          constitutive_parameters); 

        timer.enter_subsection("Initialize symbolic CL"); 
        Magnetoviscoelastic_Constitutive_Law_SD<dim> material_sd( 
          constitutive_parameters, optimizer_type, optimization_flags); 
        timer.leave_subsection(); 

        run_rheological_experiment(experimental_parameters, 
                                   material, 
                                   material_sd, 
                                   timer, 
                                   experimental_parameters.output_filename_rd); 

        std::cout << "... all calculations are correct!" << std::endl; 
      } 
    } 

  } // namespace CoupledConstitutiveLaws 

} // namespace Step71 


int main(int argc, char *argv[]) 
{ 
  Step71::SimpleExample::run(); 
  Step71::CoupledConstitutiveLaws::run(argc, argv); 

  return 0; 
} 

