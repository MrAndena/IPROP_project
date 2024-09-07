#ifndef NACA_POINTS_GENERATOR_HPP
#define NACA_POINTS_GENERATOR_HPP

#include "GmshPoint.hpp"
#include "data_struct.hpp"
#include <vector>
#include <iostream>
#include <fstream> 
#include <string>
#include <cmath>

class NACA_points_generator{

   private:
     
     data_struct m_data;

   public:

     //CONSTRUCTOR
     NACA_points_generator(const data_struct& d): m_data(d) {};

     //METHODS
     std::vector<double> GCL_nodes()  const;              //this method compute the coordinates of GCL nodes for the airfoil profile
     std::vector<GmshPoint> compute_profile() const;       //this method compute the points that made the airfoil profile
     void write_profile() const;    //this method writes in the output file the points that compose the airfoil

};

#endif //NACA_POINTS_GENERATOR_HPP