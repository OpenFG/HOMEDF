#ifndef OUTPUT_MESSAGES_H
#define OUTPUT_MESSAGES_H
/////////////////////////////////////////////////////////////////////////
#include<iostream>
#include"Eigen/Dense"

namespace homedf {

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////
/// \brief Function that outputs the header.
void Output_Header(){

  std::cout<<"\033[1;32m---------------------------------------------------------------------\033[0m"<<std::endl;
  std::cout<<"\033[1;32m---------------------------------------------------------------------\033[0m"<<std::endl;
  std::cout<<"\033[1;32m|  __    __   ________   __    __   ________   ______     ________  |\033[0m"<<std::endl;
  std::cout<<"\033[1;32m| |  |  |  | |        | |  \\  /  | |        | |      \\   |        | |\033[0m"<<std::endl;
  std::cout<<"\033[1;32m| |  |  |  | |   __   | |   \\/   | |   _____| |   _   \\  |   _____| |\033[0m"<<std::endl;
  std::cout<<"\033[1;32m| |  |__|  | |  |  |  | |        | |  |___    |  | \\   | |  |___    |\033[0m"<<std::endl;
  std::cout<<"\033[1;32m| |   __   | |  |  |  | |  |\\/|  | |   ___|   |  |  |  | |   ___|   |\033[0m"<<std::endl;
  std::cout<<"\033[1;32m| |  |  |  | |  |__|  | |  |  |  | |  |_____  |  |_/   | |  |       |\033[0m"<<std::endl;
  std::cout<<"\033[1;32m| |  |  |  | |        | |  |  |  | |        | |       /  |  |       |\033[0m"<<std::endl;
  std::cout<<"\033[1;32m| |__|  |__| |________| |__|  |__| |________| |______/   |__|       |\033[0m"<<std::endl;
  std::cout<<"\033[1;32m---------------------------------------------------------------------\033[0m"<<std::endl;
  std::cout<<"\033[1;32m---------------------------------------------------------------------\033[0m"<<std::endl;

}
  
/////////////////////////////////////////////////////////////////////////
/// \brief Function that outputs status messages.
void Output_Begin_Status(std::string message){
  
  std::cout<<"\033[1;32m[HOMEDF] \033[1;36m  ---->"+message+"<---- \033[0m"<<std::endl;

}

/////////////////////////////////////////////////////////////////////////
/// \brief Function that outputs status messages.
void Output_End_Status(std::string message){
  
  std::cout<<"\033[1;32m[HOMEDF] \033[1;36m  ---->"+message+"<---- \033[0m\n\n"<<std::endl;

}

/////////////////////////////////////////////////////////////////////////
/// \brief Function that outputs detail status messages.
void Output_Detail_Status(std::string message){
  
  std::cout<<"\033[1;32m[HOMEDF] \033[0;36m       "+message+" \033[0m"<<std::endl;

}

/////////////////////////////////////////////////////////////////////////
/// \brief Function that outputs vectors cleanly.
template<typename Vector_type>
void Output_Vector(std::string name, const Vector_type& U){

  Eigen::IOFormat fmt(10, 0, ", ", "\n", std::string( 16, ' ' )+"[", "]");
  std::cout<<"\033[1;35m"+std::string( 16, ' ' )+name+" =\n"<<U.format(fmt)<<"\033[0m\n"<<std::endl;

}

/////////////////////////////////////////////////////////////////////////
/// \brief Function that outputs scalar cleanly.
void Output_Scalar(std::string name, const double x){

  std::cout<<"\033[1;35m"+std::string( 16, ' ' )+name+" = "<<x<<"\033[0m\n"<<std::endl;

}




}//namespace homedf


#endif //#ifndef OUTPUT_MESSAGES_H
