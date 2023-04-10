#ifndef DOMAIN_ADJUSTMENT_H
#define DOMAIN_ADJUSTMENT_H
/////////////////////////////////////////////////////////////////////////
#include<math.h>
#include<iostream> 
#include"Eigen/Dense"
#include"Output_Messages.h"

namespace homedf {

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
template<typename Func_type>
void Check_and_Adjust_Domain(Func_type& F,
			     int num_points,
			     bool *Status=0,
			     int max_num_enlargement=100) {

  Output_Detail_Status("Checking and Adjusting Domain");

  double v_domain;
  
  double tol = 1e-3;

  bool domain_is_good_now = false;
  bool enlarge_pos;
  bool enlarge_neg;
  
  int num_enlargement = 0;

  while (domain_is_good_now == false && num_enlargement < max_num_enlargement) {

    enlarge_pos = false;
    enlarge_neg = false;
    
    
    if (enlarge_neg == false && F.m*pow(F.v_min*F.v_min, 4.0)*F(F.v_min) > tol){
      //enlarge in x-
      enlarge_neg = true;
    }
    if (enlarge_pos == false && F.m*pow(F.v_max*F.v_max, 4.0)*F(F.v_max) > tol){
      //enlarge in x+
      enlarge_pos = true;
    }
	 

    v_domain = F.v_max - F.v_min;
    

    if (enlarge_neg == true){
      F.v_min -= 0.5*v_domain;
    }
    if (enlarge_pos == true){
      F.v_max += 0.5*v_domain;
    }
    

    ++num_enlargement;
    
    if (enlarge_pos == false && enlarge_neg == false){
      domain_is_good_now = true;
    }
  }
  
  
  if(num_enlargement > 1 || domain_is_good_now == false){
    std::ostringstream os;
    os<<"["<<F.v_min<<", "<<F.v_max<<"] ";
    Output_Detail_Status("Domain was enlarged to: " + os.str());
    if(Status){
      *Status = false;
    }
  } else {
    Output_Detail_Status("Domain was large enough");
    if(Status){
      *Status = true;
    }
  }
  
}


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////



}//namespace homedf


#endif //#ifndef DOMAIN_ADJUSTMENT_H
