#include "NACA_points_generator.hpp"

std::vector<double> NACA_points_generator::GCL_nodes() const {

    // we define as a lambda function the GCL function
    double a = 0.0;
    double b = this->m_data.geometrical_parameters.chord_length;
    int n = 119; //ritorna n+1 points


    auto F = [a,b,n](double k) -> double {

    return (a+b)/2 - (b-a)/2*std::cos(k*M_PI/n);

    };
    
    // create and reserve memory for the vector that stores the coordinates
    std::vector<double> x_coord;
    x_coord.reserve(n); //ritorna n+1 punti
    
    // fill the vector
    for(int k=0; k<=n; ++k){

      x_coord.push_back(F(k));

    }

    return x_coord;

}

//----------------------------------------------------------------------------------------------------------------------------------------------------------

std::vector<GmshPoint> NACA_points_generator::compute_profile() const {


     // you want to retrive the coordinates of the points that build the airfoil; the front point of the airfoil is located in (0,0).

    // create vectors of double to store the coordinates, one for the x coordinates, other two for the upper
    // and the lower profile of the NACA respectively and finally the vector that stores the weights of the mesh ref 
    std::vector<double> x_coord;
    std::vector<double> up_coord;
    std::vector<double> low_coord;
    
    int n = 120;

    // reserve memory for the vectors
    x_coord.reserve(n);
    up_coord.reserve(n);
    low_coord.reserve(n);


    // The interval is (first,second). first is zero and second is equal to the length of the chord
    //we compute the x_coordinates exploiting the GCL_nodes function
    x_coord = GCL_nodes();

    // we compute now the the coordinates of both the profiles
    double t = this->m_data.geometrical_parameters.NACA_digits; //in the struct we save the last two digit, but we need the ratio of the thickness
    t = t/100;


    //lambda function that specify the profile of the NACA airfoil
    // NB: in this equation x=1 ( trailing edge of the airfoil) the thiockenss in not quit zero. if a zero thickness is required
    // modify the last coefficint : from -0.1015 to -0.1036 will result in a small change in the overall shap
    // we implemented this new version for computatinal purposis.

    double c = this->m_data.geometrical_parameters.chord_length;

    auto F_y = [t, c](double x) -> double {
    double x_norm = x / c;  // Normalize x to the range [0, 1]
    return c * 5 * t * (0.2969 * std::sqrt(x_norm) 
                      - 0.1260 * x_norm 
                      - 0.3516 * std::pow(x_norm, 2) 
                      + 0.2843 * std::pow(x_norm, 3) 
                      - 0.1036 * std::pow(x_norm, 4));
    };


    
    //we now compute the y coordinates of the upper and lower profile
    for(size_t i = 0; i<n ; ++i){

        up_coord[i] = F_y( x_coord[i]);
        low_coord[i] = -1*up_coord[i];

    }


    //finally we create a vector of Points to return, firstly we add the points of the upper profile and then the points of the lower part
    //NB: in order to have less problem with the creation of the Splines, we insert the lower points starting from the last one in vector "low_coord"

    std::vector<GmshPoint> Points;
    Points.reserve(n -2);

    for(size_t i = 0; i<n; ++i){

        GmshPoint temp( x_coord[i] ,  up_coord[i] ,0.0, 1.0);    
        Points.push_back(temp);

    }

    for(size_t i = n -2; i>0; --i){

        GmshPoint temp( x_coord[i] , low_coord[i] ,0.0, 1.0);
        Points.push_back(temp);

    }

    return Points;


}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

void NACA_points_generator::write_profile() const {
    
    //generate the output file .txt
    std::ofstream file("NACA_points.txt");

    if (!file.is_open()) {
        std::cerr << "Failed to open the file for writing."  << std::endl;
        return;
    }


    //first we retrive the points that we need
    std::vector<GmshPoint> airfoil_points = this->compute_profile();    
    
    //we then write all the points with the following sintax (remeber that in the get_local_mesh you retrive only the weigths, you 
    //have to multiply by the mesh refinement that you want)
    for(size_t i = 0; i<airfoil_points.size(); ++i){
        
        file << "Point("<<airfoil_points[i].get_tag()<<") = {"<<airfoil_points[i].get_x()<<", "<<airfoil_points[i].get_y()<<", "<<airfoil_points[i].get_z()<<", "<<airfoil_points[i].get_local_mesh_ref()<<"};"<<std::endl;

    }
    
    file.close();

    std::cout << "All the data are written in NACA_points.txt" << std::endl;

    return;





}