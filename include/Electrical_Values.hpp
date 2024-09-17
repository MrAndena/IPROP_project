#ifndef ELECTRICAL_VALUES_HPP
#define ELECTRICAL_VALUES_HPP

#include <deal.II/base/function.h>   // We inherit from the abstract template class Function
#include <deal.II/base/point.h>      // We need the Point class
#include "data_struct.hpp"
#include <cmath>

// This header file contains the declaration of four template classes that characterize the electrical problem,
// in particular we have:
// DopingValues -- PotentialValues -- ElectronInitialValues -- IonInitialValues
// These are four template classes, all of them inherit from the template class Function<dim> from deal ii and override the same 
// method "value". These template classes are used in the methods of the template class Problem.


// brief description of Function<dim>
/*
  In the deal.II library, Function<dim> is a template class representing a mathematical function defined on a
  domain with a specific dimension dim.
  This class provides a framework for defining and evaluating functions that are commonly used in finite element methods.
  It allows users to define functions symbolically or procedurally and provides methods for evaluating these functions at
  specified points in space.
  The class provides a virtual method: virtual double value(const Point<dim> &p, const unsigned int component = 0) const 
  for evaluating the function at a given point p in space. This method is expected to be overridden to define the specific
  behavior of the function.
  The optional component parameter allows for computing a specific component of the function's value if it's a vector or
  tensor function.
  Thus the main idea behinf Function<dim> is to override the method value in order to generate a specific math function
  for our proposals.
*/
using namespace dealii; 

// -------------------------------------- DOPING VALUES ------------------------------------------------------------------------

template <int dim>
class DopingValues : public Function<dim> {

public:
     
  DopingValues(const data_struct & d) : Function<dim>(), L(d.geometrical_parameters.L),
                                           A(d.drift_diffusion.physical_parameters.A),
                                           D(d.drift_diffusion.physical_parameters.D){}

  virtual double value(const Point<dim> & p, const unsigned int component = 0) const override; 

private:

  double L;
  double A;
  double D;

};

// ------------------------------------- POTENTIAL VALUES ------------------------------------------------------------------------

template <int dim>
class PotentialValues : public Function<dim> {

public:

  PotentialValues(const data_struct & d) : Function<dim>(), L(d.geometrical_parameters.L),
                                           V_TH(d.drift_diffusion.physical_parameters.V_TH),
                                           A(d.drift_diffusion.physical_parameters.A),
                                           D(d.drift_diffusion.physical_parameters.D),
                                           ni(d.drift_diffusion.physical_parameters.ni) {}

  virtual double value(const Point<dim> & p, const unsigned int component = 0) const override;

private:

  double A;
  double D;
  double L;
  double V_TH;
  double ni;

};

// ------------------------------------- ELECTRON INITIAL VALUES ---------------------------------------------------------------
template <int dim>
class ElectronInitialValues : public Function<dim> {

public:

  ElectronInitialValues(const data_struct & d) : Function<dim>(), L(d.geometrical_parameters.L)  {

      double A = d.drift_diffusion.physical_parameters.A;
      double D = d.drift_diffusion.physical_parameters.D;
      double ni = d.drift_diffusion.physical_parameters.ni;

      const double P2 = A/2. + std::sqrt(A*A + 4.*ni*ni)/2.; // [m^-3] Holes density on boundary 2  

      N1 = D/2. + std::sqrt(D*D + 4.*ni*ni)/2.; // [m^-3] Electron density on boundary 1
      N2 = ni*ni/P2; // [m^-3] Electron density on boundary 2
  }

  virtual double value(const Point<dim> & p, const unsigned int component = 0) const override;

private:
  
  double N1;
  double N2;
  double L;

};

// ------------------------------------- ION INITIAL VALUES --------------------------------------------------------------------

template <int dim>
class HoleInitialValues : public Function<dim> {

public:
  // Constructor that accepts the constants struct
  HoleInitialValues(const data_struct & d) 
    : Function<dim>(), L(d.geometrical_parameters.L) {
      
      double A = d.drift_diffusion.physical_parameters.A;
      double D = d.drift_diffusion.physical_parameters.D;
      double ni = d.drift_diffusion.physical_parameters.ni;

      const double N1 = D/2. + std::sqrt(D*D + 4.*ni*ni)/2.; // [m^-3] Electron density on boundary 1

      P2 = A/2. + std::sqrt(A*A + 4.*ni*ni)/2.; // [m^-3] Holes density on boundary 2   
      P1 = ni*ni/N1; // [m^-3] Holes density on boundary 1

    }

  // Override the value function
  virtual double value(const Point<dim> & p, const unsigned int component = 0) const override;

private:
  // Member variables to hold P1,P2,L
  double P1;
  double P2;
  double L;
};

#include "Electrical_Values_impl.hpp" // we refer to this file for the definitions of the templates

#endif // ELECTRICAL_VALUES_HPP
