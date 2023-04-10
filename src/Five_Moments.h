#ifndef FIVE_MOMENTS_H
#define FIVE_MOMENTS_H
/////////////////////////////////////////////////////////////////////////
#include<math.h>
#include<sstream>
#include<iostream>
#include<algorithm> 
#include"Eigen/Dense"
#include"Output_Messages.h"

namespace homedf {

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
class Five_Moments{

 public:
  /////////////////////////////////////////////////////////////////////////
  /// \brief Defines the solution vector type.
  using Vector_type = Eigen::Matrix<double,5,1>;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Defines the matrix type.
  using Matrix_type = Eigen::Matrix<double,5,5>;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Default constructor, use compiler-generated version.
  Five_Moments() = default;
  
  /////////////////////////////////////////////////////////////////////////
  /// \brief Copy constructor, use compiler-generated version.
  Five_Moments(const Five_Moments&) = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Move constructor, use compiler-generated version.
  Five_Moments(Five_Moments&&) = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Copy assignment, use compiler-generated version.
  Five_Moments& operator=(const Five_Moments&) = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Move assignment, use compiler-generated version.
  Five_Moments& operator=(Five_Moments&&) = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Returns vector of the velocity weight.
  static auto W (double v){
    
    Vector_type w;
    w(0) = 1.0;
    w(1) = v;
    w(2) = v*v;
    w(3) = v*v*v;
    w(4) = v*v*v*v;

    return w;
  }

  /////////////////////////////////////////////////////////////////////////
  /// \brief Funtion that returns the solution vector from the primitive variable.
  static auto Get_Solution_Vector_from_Primitive(const Vector_type& V);

  /////////////////////////////////////////////////////////////////////////
  /// \brief Funtion that returns the primitive variable from the solution vector.
  static auto Get_Primitive_from_Solution_Vector(const Vector_type& U);

  /////////////////////////////////////////////////////////////////////////
  /// \brief Funtion that returns the non-dimensional solution vector.
  static auto Non_Dimensionalize_Solution_Vector(const Vector_type& U);

  /////////////////////////////////////////////////////////////////////////
  /// \brief Funtion that returns the non-dimensional primitive vector.
  static auto Non_Dimensionalize_Primitive_Vector(const Vector_type& V);

  /////////////////////////////////////////////////////////////////////////
  /// \brief Funtion that returns the dimensional free coefficients.
  auto Dimensionalize_Alphas(const Vector_type& U, Vector_type a);

  /////////////////////////////////////////////////////////////////////////
  /// \brief Funtion that returns r for a specfic value of sigma.
  static double r_from_sigma(const Vector_type& V, double sigma);

  /////////////////////////////////////////////////////////////////////////
  /// \brief Funtion that returns sigma for a specfic value of r.
  static double sigma_from_r(const Vector_type& V);


};

}//namespace homedf

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
//               Inline member functions
//

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
auto homedf::Five_Moments::Get_Solution_Vector_from_Primitive(const Vector_type& V){

  double rho = V(0);
  double u   = V(1);
  double p   = V(2);
  double q   = V(3);
  double r   = V(4);

  Vector_type U;
  U(0) = rho;
  U(1) = rho*u;
  U(2) = rho*u*u + p;
  U(3) = rho*u*u*u + 3.0*u*p + q;
  U(4) = rho*u*u*u*u + 6.0*u*u*p + 4.0*u*q + r;

  return U;
  
}

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
auto homedf::Five_Moments::Get_Primitive_from_Solution_Vector(const Vector_type& U){

  double rho = U(0);
  double u   = U(1)/U(0);
  double p   = U(2) - rho*u*u;
  double q   = U(3) - rho*u*u*u - 3.0*u*p;
  double r   = U(4) - rho*u*u*u*u - 6.0*u*u*p - 4.0*u*q;

  Vector_type V;
  V(0) = rho;
  V(1) = u;
  V(2) = p;
  V(3) = q;
  V(4) = r;

  return V;
  
}

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
auto homedf::Five_Moments::Non_Dimensionalize_Solution_Vector(const Vector_type& U){

  double ND_factor = sqrt(U(2)/U(0));

  Vector_type U_ND;
  U_ND(0) = U(0)/U(0);
  U_ND(1) = U(1)/U(0)/ND_factor;
  U_ND(2) = U(2)/U(0)/ND_factor/ND_factor;
  U_ND(3) = U(3)/U(0)/ND_factor/ND_factor/ND_factor;
  U_ND(4) = U(4)/U(0)/ND_factor/ND_factor/ND_factor/ND_factor;
  
  return U_ND;
  
}

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
auto homedf::Five_Moments::Non_Dimensionalize_Primitive_Vector(const Vector_type& V){

  double ND_factor = sqrt((V(0)*V(1)*V(1) + V(2))/V(0));

  Vector_type V_ND;
  V_ND(0) = V(0)/V(0);
  V_ND(1) = V(1)/ND_factor;
  V_ND(2) = V(2)/V(0)/ND_factor/ND_factor;
  V_ND(3) = V(3)/V(0)/ND_factor/ND_factor/ND_factor;
  V_ND(4) = V(4)/V(0)/ND_factor/ND_factor/ND_factor/ND_factor;
  
  return V_ND;
  
}

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
auto homedf::Five_Moments::Dimensionalize_Alphas(const Vector_type& U,
						 Vector_type a){

  double NDF = sqrt(U(2)/U(0));

  a(0) += log(U(0)/NDF);
  a(1) *= 1.0/NDF;
  a(2) *= 1.0/NDF/NDF;
  a(3) *= 1.0/NDF/NDF/NDF;
  a(4) *= 1.0/NDF/NDF/NDF/NDF;

  return a;
  
}

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
double homedf::Five_Moments::r_from_sigma(const Vector_type& V, double sigma){

  double min_sigma = 1e-5;

  if(sigma < min_sigma){
    sigma = min_sigma;
  }
  
  double rho = V(0);
  double p   = V(2);
  double q   = V(3);

  double r = 1.0/sigma*q*q/p + (3.0-2.0*sigma)*p*p/rho;

  if(!std::isfinite(r)){
    throw std::invalid_argument( "Cannot find r for this value of sigma" );
  }
  
  return r;
  
}

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
double homedf::Five_Moments::sigma_from_r(const Vector_type& V){

  double rho = V(0);
  double p   = V(2);
  double q   = V(3);
  double r   = V(4);

  double A = 3.0*p*p-rho*r;
  
  double sigma = (A + sqrt(A*A + 8.0*rho*p*q*q))/(4.0*p*p);

  return sigma;
}

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

  
#endif //#ifndef FIVE_MOMENTS_H
