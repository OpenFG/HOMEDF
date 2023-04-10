#ifndef DISTRIBUTION_FUNCTION_H
#define DISTRIBUTION_FUNCTION_H
/////////////////////////////////////////////////////////////////////////
#include<math.h>
#include<sstream>
#include<iostream>
#include<algorithm> 
#include"Eigen/Dense"
#include"Output_Messages.h"
#include"Integration_Helper.h"

namespace homedf {

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
template<typename Moment_type>
class Distribution_Function : public Moment_type{

 public:

  /////////////////////////////////////////////////////////////////////////
  /// \brief Using the function W that is in the Moment_type class.
  using Moment_type::W;
  
  /////////////////////////////////////////////////////////////////////////
  /// \brief Defines the solution vector type.
  using Vector_type = typename Moment_type::Vector_type;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Defines the matrix type.
  using Matrix_type = typename Moment_type::Matrix_type;
  
  /////////////////////////////////////////////////////////////////////////
  /// \brief Default constructor, use compiler-generated version.
  Distribution_Function() = default;
  
  /////////////////////////////////////////////////////////////////////////
  /// \brief Copy constructor, use compiler-generated version.
  Distribution_Function(const Distribution_Function&) = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Move constructor, use compiler-generated version.
  Distribution_Function(Distribution_Function&&) = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Copy assignment, use compiler-generated version.
  Distribution_Function& operator=(const Distribution_Function&) = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Move assignment, use compiler-generated version.
  Distribution_Function& operator=(Distribution_Function&&) = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Constructor that takes the vector of alphas.
  Distribution_Function(double m_in, Vector_type alphas_in) : m(m_in),
                                                              a(alphas_in){}

  /////////////////////////////////////////////////////////////////////////
  /// \brief Evalultates the Distribution function.
  auto operator()(double v) const {
    return F(v);
  } 

  /////////////////////////////////////////////////////////////////////////
  /// \brief Mass of the molecule
  double m;
  
  /////////////////////////////////////////////////////////////////////////
  /// \brief Vector of all the alphas.
  Vector_type a;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Domain to be integrated in vx.
  double v_min = -16.0;
  double v_max =  16.0;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Dimensionalize the alphas based on a dimensional solution vector.
  auto Dimensionalize_Alphas(const Vector_type& U){
    return Moment_type::Dimensionalize_Alphas(U, a);
  }

  /////////////////////////////////////////////////////////////////////////
  /// \brief Function that calculates the Eigenvalues of the Jacobian.
  auto Eigen_Solver(double m, int num_integration_points); 

  
  
 private:
  
  /////////////////////////////////////////////////////////////////////////
  /// \brief Returns the value of the distribution function at one point.
  double F (double v) const {
    return exp(a.dot(W(v)));
  }
  

};
 
}//namespace homedf

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
//               Inline member functions
//

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
template<typename Moment_type>
auto homedf::Distribution_Function<Moment_type>::Eigen_Solver(double m,
							      int num_integration_points){

  Output_Begin_Status("Starting Eigen Solver to solve for Eigenvalues");
  
  Matrix_type H    = Integrate_H(*this, num_integration_points);
  Matrix_type J    = Integrate_J(*this, num_integration_points);
  Matrix_type dFdU = J*H.inverse();
  
  Output_Detail_Status("Solving Eigenvalues");

  Eigen::EigenSolver<Eigen::MatrixXd> es(dFdU);
  auto eigenvalues = es.eigenvalues().eval();
  
  if(isnan(a.norm())){
    eigenvalues.fill(NAN);
    Output_End_Status("The Eigenvalues could not be found");
  } else {
    Output_End_Status("The Eigenvalues were found Successfully");
  }

  return eigenvalues;
}

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

#endif //#ifndef DISTRIBUTION_FUNCTION_H
