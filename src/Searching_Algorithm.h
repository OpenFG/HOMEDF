#ifndef SEARCHING_ALGORITHM_H
#define SEARCHING_ALGORITHM_H
/////////////////////////////////////////////////////////////////////////
#include<cmath>
#include<math.h>
#include<iostream> 
#include"Eigen/Dense"
#include"Distribution_Function.h"
#include"Output_Messages.h"
#include"Domain_Adjustment.h"
#include"Integration_Helper.h"

namespace homedf {

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////
/// \brief Function that finds the alphas corresponding to the gaussian distribution.
template<typename Vector_type>
auto Get_Initial_Guess_for_Alphas(double m,
				  Vector_type& V){

  double rho = V(0);
  double u   = V(1);
  double p   = V(2);

  Vector_type alpha;
  alpha.fill(0.0);
  alpha(0) = log(rho/m*sqrt(rho/2.0/M_PI/p)) - rho/2.0/p*u*u;
  alpha(1) = rho/p*u;
  alpha(2) = -rho/2.0/p;

  return alpha;
}
  
/////////////////////////////////////////////////////////////////////////
/// \brief Function that solves for the alphas.
template<typename Vector_type, typename Func_type>
void Find_Alphas(Vector_type& U_target,
		 Func_type& F,
		 int num_points){

  using Matrix_type = typename Func_type::Matrix_type;
  
  Matrix_type d2Jda2       = Integrate_Solution_Vector_prime(F, num_points);
  Vector_type dJda         = d2Jda2.col(0) - U_target;
  Matrix_type d2Jda2_new;
  Vector_type dJda_new;
  Vector_type delta_a;
  Output_Detail_Status("The Norm of the Difference between U and U_target is: "+std::to_string(dJda.norm())+"\n");

  
  int max_num_newton_steps        = 500;
  int num_newton_steps            = 0;
  int max_num_newton_steps_resize = 10;
  int num_newton_steps_resize     = 0;
  double tol = 1e-9;

  while ((dJda.array().abs() > tol).any()){

    Output_Detail_Status("Updating the Alphas");
    
    delta_a = d2Jda2.inverse()*dJda;
    F.a = F.a - delta_a;

    d2Jda2_new = Integrate_Solution_Vector_prime(F, num_points);
    dJda_new = d2Jda2_new.col(0) - U_target;
    
    //Backtracking 
    num_newton_steps_resize = 0;
    while(!(dJda_new.norm() < dJda.norm()) && num_newton_steps_resize < max_num_newton_steps_resize){
      Output_Detail_Status("The step was too large, Taking a smaller step");
      delta_a = 0.5*delta_a;
      F.a += delta_a;
      d2Jda2_new = Integrate_Solution_Vector_prime(F, num_points);
      dJda_new = d2Jda2_new.col(0) - U_target;
      ++num_newton_steps_resize;
    }

    d2Jda2 = d2Jda2_new;
    dJda = dJda_new;

    ++num_newton_steps;

    //Returns NaN if could not converge
    if(num_newton_steps >= max_num_newton_steps || num_newton_steps_resize >= max_num_newton_steps_resize){
      Output_Detail_Status("Could not converge, problem aborted");
      Output_End_Status("The Alphas were not found for U_target, could not converge");
      F.a.fill(NAN);
      return;
    }
    
    Output_Detail_Status("The Norm of the Difference between U and U_target is: "
			 +std::to_string(dJda.norm())+"\n"); 
    
  }

  Output_Detail_Status("Alphas were found");
  
}

/////////////////////////////////////////////////////////////////////////
/// \brief Function that will keep solving for alphas until the interation doman is large enough.
template<typename Vector_type, typename Func_type>
auto Advance_F(double m,
	       Vector_type& U_target,
	       int num_points,
	       Func_type& F){

  Output_Begin_Status("Starting the Search Algorithm for the Alphas");
  
  //Initialize Distribution function with initial guess
  Check_and_Adjust_Domain(F, num_points);

  //Store initial alphas
  Vector_type initial_a = F.a;
  bool Domain_Large_Enough;

  //Solves entropy maximization problem
  Find_Alphas(U_target, F, num_points);

  if(isnan(F.a.norm())){
    return;
  }

  Output_Detail_Status("Verifying that the domain was big enough");
  Check_and_Adjust_Domain(F, num_points, &Domain_Large_Enough, 1);
  
  //Re-solves entropy maximization problem if domain was not big enough
  int num_domain_enlargement = 0;
 
  while(Domain_Large_Enough == false){

    if(num_domain_enlargement > 100){
      Output_Detail_Status("The maximum domain enlargement was exceeded");
      Output_End_Status("The Alphas were not found for U_target, domain was too small");
      F.a.fill(NAN);
      return;
      //throw std::runtime_error( "Domain is too small even after many enlargement" );
    }
    

    Output_Detail_Status("The Domain was too small resolving alphas");
    F.a = initial_a;
    Find_Alphas(U_target, F, num_points);

    if(isnan(F.a.norm())){
      return;
    }
    
    Output_Detail_Status("Verifying that the domain was big enough");
    Check_and_Adjust_Domain(F, num_points, &Domain_Large_Enough, 1);
    ++num_domain_enlargement;

  }
  
  Output_End_Status("The Alphas were found for U_target");
  
}



/////////////////////////////////////////////////////////////////////////
/// \brief Function that searches for alphas from scratch.
template<typename Moment_type, typename Vector_type>
auto Search_Algorithm(double m,
		      Vector_type& V,
		      int num_integration_points,
		      int num_steps = 1){

  auto a_in = homedf::Get_Initial_Guess_for_Alphas(m, V);
  auto F = Distribution_Function<Moment_type>(m, a_in);
  Vector_type U = F.Get_Solution_Vector_from_Primitive(V);

  Advance_F(m, U, num_integration_points, F);
  
  
  return F;
}

/////////////////////////////////////////////////////////////////////////
/// \brief Function that searches for alphas from previously obtained F.
template<typename Moment_type, typename Vector_type, typename Func_type>
auto Search_Algorithm(double m,
		      Vector_type& V,
		      int num_integration_points,
		      Func_type F){
  
  Vector_type U = F.Get_Solution_Vector_from_Primitive(V);
  Advance_F(m, U, num_integration_points, F);
  
  return F;
}

/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

}//namespace homedf


#endif //#ifndef SEARCHING_ALGORITHM_H
