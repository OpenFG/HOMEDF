#ifndef INTEGRATION_HELPER_H
#define INTEGRATION_HELPER_H
/////////////////////////////////////////////////////////////////////////
#include<cmath>
#include<math.h>
#include<iostream> 
#include"Eigen/Dense"
#include"Output_Messages.h"

namespace homedf {

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
template<typename Func_type>
auto Integrate_Solution_Vector(Func_type& F,
			       int num_points_v) {

  using Vector_type = typename Func_type::Vector_type;
  
  Output_Begin_Status("Starting Integration process for the Solution Vector");
  
  auto Integrand = [&] (const double& v) -> Vector_type {
    return F.m*F.W(v)*F(v);
  };

  Output_Detail_Status("Integrating the Solution Vector");
  
  Vector_type U = gaussian_integrate<5>(Integrand, F.v_min, F.v_max, num_points_v);

  Output_End_Status("The Solution Vector was Integrated Successfully");
  
  return U;
}

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
template<typename Func_type>
auto Integrate_Solution_Vector_prime(Func_type& F,
				     int num_points_v) {

  using Matrix_type = typename Func_type::Matrix_type;
  
  auto Integrand = [&] (const double& v) -> Matrix_type {
    return F.m*F.W(v)*F.W(v).transpose()*F(v);
  };
    
  Output_Detail_Status("Integrating the Solution Vector prime");

  Matrix_type U = gaussian_integrate<5>(Integrand, F.v_min, F.v_max, num_points_v);

  return U;
}

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
template<typename Func_type>
auto Integrate_Flux_Vector(Func_type& F,
			   int num_points_v) {

  using Vector_type = typename Func_type::Vector_type;
  
  Output_Begin_Status("Starting Integration process for the Flux Vector");
  
  auto Integrand = [&] (const double& v) -> Vector_type {
    return F.m*v*F.W(v)*F(v);
  };

  Output_Detail_Status("Integrating the Flux Vector");
  
  Vector_type Flux = gaussian_integrate<5>(Integrand, F.v_min, F.v_max, num_points_v);

  Output_End_Status("The Flux Vector was Integrated Successfully");
  
  return Flux;
}

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
template<typename Func_type>
auto Integrate_H(Func_type& F,
		 int num_points_v) {

  using Matrix_type = typename Func_type::Matrix_type;
  
  auto Integrand = [&] (const double& v) -> Matrix_type {
    return F.m*F.W(v)*F.W(v).transpose()*F(v);
  };
    
  Output_Detail_Status("Integrating H");

  Matrix_type U = gaussian_integrate<5>(Integrand, F.v_min, F.v_max, num_points_v);

  return U;
}

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
template<typename Func_type>
auto Integrate_J(Func_type& F,
		 int num_points_v) {

  using Matrix_type = typename Func_type::Matrix_type;
  
  auto Integrand = [&] (const double& v) -> Matrix_type {
    return F.m*v*F.W(v)*F.W(v).transpose()*F(v);
  };
    
  Output_Detail_Status("Integrating J");

  Matrix_type U = gaussian_integrate<5>(Integrand, F.v_min, F.v_max, num_points_v);

  return U;
}
  
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
template<int n, typename Func_type>
auto Integrate_Velocity_Moment(Func_type& F,
			       int num_points_v) {

  Output_Begin_Status("Starting Integration process for the Velocity Moment");

  auto Integrand = [&] (const double& v){
    return F.m*pow(v,n)*F(v);
  };
  
  Output_Detail_Status("Integrating the Velocity Moment");

  double M = gaussian_integrate<5>(Integrand, F.v_min, F.v_max, num_points_v);
  
  Output_End_Status("The Velocity Moment was Integrated Successfully");

  return M;
}

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
template<int n, typename Func_type>
auto Integrate_Random_Moment(Func_type& F,
			     int num_points_v) {

  Output_Begin_Status("Starting Integration process for the Random Moment");

  using Vector_type = typename Func_type::Vector_type;
  
  auto Integrand_U = [&] (const double& v) -> Vector_type {
    return F.m*F.W(v)*F(v);
  };
  
  Vector_type U = gaussian_integrate<5>(Integrand_U, F.v_min, F.v_max, num_points_v);
  
  double u = U(1)/U(0);
  
  auto Integrand = [&] (const double& v){
    return F.m*pow((v-u),n)*F(v);
  };
  
  Output_Detail_Status("Integrating the Random Moment");

  double M = gaussian_integrate<5>(Integrand, F.v_min, F.v_max, num_points_v);
  
  Output_End_Status("The Random Moment was Integrated Successfully");

  return M;
}
 
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////


}//namespace homedf


#endif //#ifndef INTEGRATION_HELPER_H
