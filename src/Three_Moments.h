#ifndef THREE_MOMENTS_H
#define THREE_MOMENTS_H
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
class Three_Moments{

 public:
  /////////////////////////////////////////////////////////////////////////
  /// \brief Defines the solution vector type.
  using Vector_type = Eigen::Matrix<double,3,1>;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Defines the matrix type.
  using Matrix_type = Eigen::Matrix<double,3,3>;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Default constructor, use compiler-generated version.
  Three_Moments() = default;
  
  /////////////////////////////////////////////////////////////////////////
  /// \brief Copy constructor, use compiler-generated version.
  Three_Moments(const Three_Moments&) = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Move constructor, use compiler-generated version.
  Three_Moments(Three_Moments&&) = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Copy assignment, use compiler-generated version.
  Three_Moments& operator=(const Three_Moments&) = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Move assignment, use compiler-generated version.
  Three_Moments& operator=(Three_Moments&&) = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Returns vector of the velocity weight.
  static auto W (double v){
    
    Vector_type w;
    w(0) = 1.0;
    w(1) = v;
    w(2) = v*v;
    
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


};

}//namespace homedf

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
//               Inline member functions
//

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
auto homedf::Three_Moments::Get_Solution_Vector_from_Primitive(const Vector_type& V){

  double rho = V(0);
  double u   = V(1);
  double p   = V(2);
  
  Vector_type U;
  U(0) = rho;
  U(1) = rho*u;
  U(2) = rho*u*u + p;
  
  return U;
  
}

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
auto homedf::Three_Moments::Get_Primitive_from_Solution_Vector(const Vector_type& U){

  double rho = U(0);
  double u   = U(1)/U(0);
  double p   = U(2) - rho*u*u;
  
  Vector_type V;
  V(0) = rho;
  V(1) = u;
  V(2) = p;
  
  return V;
  
}

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
auto homedf::Three_Moments::Non_Dimensionalize_Solution_Vector(const Vector_type& U){

  double ND_factor = sqrt(U(2)/U(0));

  Vector_type U_ND;
  U_ND(0) = U(0)/U(0);
  U_ND(1) = U(1)/U(0)/ND_factor;
  U_ND(2) = U(2)/U(0)/ND_factor/ND_factor;
  
  return U_ND;
  
}

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
auto homedf::Three_Moments::Non_Dimensionalize_Primitive_Vector(const Vector_type& V){

  double ND_factor = sqrt((V(0)*V(1)*V(1) + V(2))/V(0));

  Vector_type V_ND;
  V_ND(0) = V(0)/V(0);
  V_ND(1) = V(1)/ND_factor;
  V_ND(2) = V(2)/V(0)/ND_factor/ND_factor;
  
  return V_ND;
  
}

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
auto homedf::Three_Moments::Dimensionalize_Alphas(const Vector_type& U,
  						  Vector_type a){

  double NDF = sqrt(U(2)/U(0));
  
  a(0) += log(U(0)/NDF);
  a(1) *= 1.0/NDF;
  a(2) *= 1.0/NDF/NDF;

  return a;
  
}

  
#endif //#ifndef THREE_MOMENTS_H
