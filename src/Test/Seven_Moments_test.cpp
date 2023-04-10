#include<iostream>
#include<iomanip>
#include<string>
#include<sstream>
#include<fstream>
#include<cmath>
#include<math.h>
#include"Eigen/Dense"
#include"gtest/gtest.h"
#include"../Integration.h"
#include"../Integration_Helper.h"
#include"../Distribution_Function.h"
#include"../Searching_Algorithm.h"
#include"../Output_Messages.h"
#include"../Three_Moments.h"
#include"../Five_Moments.h"
#include"../Seven_Moments.h"


/////////////////////////////////////////////////////////////////////////
/// Test for the right weight vector
TEST(Seven_Moments, W){

  double v  = 6.3;
  double W0 = 1.0;
  double W1 = v;
  double W2 = v*v;
  double W3 = v*v*v;
  double W4 = v*v*v*v;
  double W5 = v*v*v*v*v;
  double W6 = v*v*v*v*v*v;
  
  auto W = homedf::Seven_Moments::W(v);

  EXPECT_NEAR(W0, W(0), 1.0e-10);
  EXPECT_NEAR(W1, W(1), 1.0e-10);
  EXPECT_NEAR(W2, W(2), 1.0e-10);
  EXPECT_NEAR(W3, W(3), 1.0e-10);
  EXPECT_NEAR(W4, W(4), 1.0e-10);
  EXPECT_NEAR(W5, W(5), 1.0e-10);
  EXPECT_NEAR(W6, W(6), 1.0e-10);

}

/////////////////////////////////////////////////////////////////////////
/// Test for V -> U conversion
TEST(Seven_Moments, V_to_U){

  homedf::Seven_Moments::Vector_type V;
  V(0) = 1.3;
  V(1) = 1.25;
  V(2) = 2.65;
  V(3) = 1.2;
  V(4) = 6.9;
  V(5) = 6.08313;
  V(6) = 20.9793;
  
  homedf::Seven_Moments::Vector_type U_EXACT;
  U_EXACT(0) = 1.3;
  U_EXACT(1) = 1.625;
  U_EXACT(2) = 4.68125;
  U_EXACT(3) = 13.6765624999;
  U_EXACT(4) = 40.9175781249;
  U_EXACT(5) = 123.683227656;
  U_EXACT(6) = 377.2015298828125;
  
  auto U = homedf::Seven_Moments::Get_Solution_Vector_from_Primitive(V);

  EXPECT_NEAR(U_EXACT(0), U(0), 1.0e-8);
  EXPECT_NEAR(U_EXACT(1), U(1), 1.0e-8);
  EXPECT_NEAR(U_EXACT(2), U(2), 1.0e-8);
  EXPECT_NEAR(U_EXACT(3), U(3), 1.0e-8);
  EXPECT_NEAR(U_EXACT(4), U(4), 1.0e-8);
  EXPECT_NEAR(U_EXACT(5), U(5), 1.0e-8);
  EXPECT_NEAR(U_EXACT(6), U(6), 1.0e-8);

}

/////////////////////////////////////////////////////////////////////////
/// Test for U -> V conversion
TEST(Seven_Moments, U_to_V){

  homedf::Seven_Moments::Vector_type U;
  U(0) = 1.3;
  U(1) = 1.625;
  U(2) = 4.68125;
  U(3) = 13.6765624999;
  U(4) = 40.9175781249;
  U(5) = 123.683227656;
  U(6) = 377.2015298828125;
  
  homedf::Seven_Moments::Vector_type V_EXACT;
  V_EXACT(0) = 1.3;
  V_EXACT(1) = 1.25;
  V_EXACT(2) = 2.65;
  V_EXACT(3) = 1.2;
  V_EXACT(4) = 6.9;
  V_EXACT(5) = 6.08313;
  V_EXACT(6) = 20.9793;

  auto V = homedf::Seven_Moments::Get_Primitive_from_Solution_Vector(U);

  EXPECT_NEAR(V_EXACT(0), V(0), 1.0e-8);
  EXPECT_NEAR(V_EXACT(1), V(1), 1.0e-8);
  EXPECT_NEAR(V_EXACT(2), V(2), 1.0e-8);
  EXPECT_NEAR(V_EXACT(3), V(3), 1.0e-8);
  EXPECT_NEAR(V_EXACT(4), V(4), 1.0e-8);
  EXPECT_NEAR(V_EXACT(5), V(5), 1.0e-8);
  EXPECT_NEAR(V_EXACT(6), V(6), 1.0e-8);

}

/////////////////////////////////////////////////////////////////////////
/// Test for non-dimensionalization of U
TEST(Seven_Moments, U_ND){

  homedf::Seven_Moments::Vector_type U;
  U(0) = 1.3;
  U(1) = 1.625;
  U(2) = 4.68125;
  U(3) = 13.6765624999;
  U(4) = 40.9175781249;
  U(5) = 123.683227656;
  U(6) = 377.2015298828125;
  
  homedf::Seven_Moments::Vector_type U_ND_EXACT;
  U_ND_EXACT(0) = 1.0;
  U_ND_EXACT(1) = 0.65871988167203011;
  U_ND_EXACT(2) = 1.0;
  U_ND_EXACT(3) = 1.5395928235602456;
  U_ND_EXACT(4) = 2.4273343541217773;
  U_ND_EXACT(5) = 3.8665288142510428;
  U_ND_EXACT(6) = 6.2140485921659128;

  auto U_ND = homedf::Seven_Moments::Non_Dimensionalize_Solution_Vector(U);

  EXPECT_NEAR(U_ND_EXACT(0), U_ND(0), 1.0e-8);
  EXPECT_NEAR(U_ND_EXACT(1), U_ND(1), 1.0e-8);
  EXPECT_NEAR(U_ND_EXACT(2), U_ND(2), 1.0e-8);
  EXPECT_NEAR(U_ND_EXACT(3), U_ND(3), 1.0e-8);
  EXPECT_NEAR(U_ND_EXACT(4), U_ND(4), 1.0e-8);
  EXPECT_NEAR(U_ND_EXACT(5), U_ND(5), 1.0e-8);
  EXPECT_NEAR(U_ND_EXACT(6), U_ND(6), 1.0e-8);

}

/////////////////////////////////////////////////////////////////////////
/// Test for non-dimensionalization of V
TEST(Seven_Moments, V_ND){

  homedf::Seven_Moments::Vector_type V;
  V(0) = 1.2;
  V(1) = 0.15;
  V(2) = 1.1;
  V(3) = 0.2;
  V(4) = 3.5;
  V(5) = 0.1;
  V(6) = 13.0;
  
  homedf::Seven_Moments::Vector_type V_ND_EXACT;
  V_ND_EXACT(0) = 1.0;
  V_ND_EXACT(1) = 0.1547818111;
  V_ND_EXACT(2) = 0.9760425909;
  V_ND_EXACT(3) = 0.1831195636;
  V_ND_EXACT(4) = 3.306750732;
  V_ND_EXACT(5) = 0.09749045088;
  V_ND_EXACT(6) = 13.07778208;

  auto V_ND = homedf::Seven_Moments::Non_Dimensionalize_Primitive_Vector(V);

  EXPECT_NEAR(V_ND_EXACT(0), V_ND(0), 1.0e-8);
  EXPECT_NEAR(V_ND_EXACT(1), V_ND(1), 1.0e-8);
  EXPECT_NEAR(V_ND_EXACT(2), V_ND(2), 1.0e-8);
  EXPECT_NEAR(V_ND_EXACT(3), V_ND(3), 1.0e-8);
  EXPECT_NEAR(V_ND_EXACT(4), V_ND(4), 1.0e-8);
  EXPECT_NEAR(V_ND_EXACT(5), V_ND(5), 1.0e-8);
  EXPECT_NEAR(V_ND_EXACT(6), V_ND(6), 1.0e-8);

}

/////////////////////////////////////////////////////////////////////////
/// Test for dimensionalization of alphas
TEST(Seven_Moments, alpha){

  homedf::Seven_Moments::Vector_type V;
  V(0) = 1.2;
  V(1) = 0.15;
  V(2) = 1.1;
  V(3) = 0.2;
  V(4) = 3.5;
  V(5) = 0.1;
  V(6) = 13.0;
  
  auto U = homedf::Seven_Moments::Get_Solution_Vector_from_Primitive(V);

  homedf::Seven_Moments::Vector_type a_ND;
  a_ND(0) = -0.109853769;
  a_ND(1) = 0.007400686626;
  a_ND(2) = -5.2167344;
  a_ND(3) = 0.00384952724;
  a_ND(4) = 2.317463195;
  a_ND(5) = 0.009991441719;
  a_ND(6) = -0.2749650342;

  homedf::Seven_Moments::Vector_type a_EXACT;
  a_EXACT(0) = 0.1038489487;
  a_EXACT(1) = 0.007636611196;
  a_EXACT(2) = -5.554641774;
  a_EXACT(3) = 0.004229542489;
  a_EXACT(4) = 2.627407925;
  a_EXACT(5) = 0.0116888419;
  a_EXACT(6) = -0.3319322582;

  homedf::Seven_Moments Seven_Moments_Object;
  auto a = Seven_Moments_Object.Dimensionalize_Alphas(U, a_ND);
  
  EXPECT_NEAR(a_EXACT(0), a(0), 1.0e-8);
  EXPECT_NEAR(a_EXACT(1), a(1), 1.0e-8);
  EXPECT_NEAR(a_EXACT(2), a(2), 1.0e-8);
  EXPECT_NEAR(a_EXACT(3), a(3), 1.0e-8);
  EXPECT_NEAR(a_EXACT(4), a(4), 1.0e-8);
  EXPECT_NEAR(a_EXACT(5), a(5), 1.0e-8);
  EXPECT_NEAR(a_EXACT(6), a(6), 1.0e-8);

}

/////////////////////////////////////////////////////////////////////////
/// Test for sigma -> r
TEST(Seven_Moments, sigma_to_r){

  homedf::Seven_Moments::Vector_type V;
  V(0) = 1.2;
  V(1) = 0.15;
  V(2) = 1.1;
  V(3) = 0.2;
  V(4) = 3.5;
  V(5) = 0.1;
  V(6) = 13.0;

  double sigma   = 0.0608399;
  double r_EXACT = 3.5000000894;
  double r       = homedf::Seven_Moments::r_from_sigma(V, sigma);
  
  EXPECT_NEAR(r_EXACT, r, 1.0e-8);

}

/////////////////////////////////////////////////////////////////////////
/// Test for r -> sigma
TEST(Seven_Moments, r_to_sigma){

  homedf::Seven_Moments::Vector_type V;
  V(0) = 1.2;
  V(1) = 0.15;
  V(2) = 1.1;
  V(3) = 0.2;
  V(4) = 3.5;
  V(5) = 0.1;
  V(6) = 13.0;

  double sigma_EXACT = 0.0608399;
  double sigma = homedf::Seven_Moments::sigma_from_r(V);
  
  EXPECT_NEAR(sigma_EXACT, sigma, 1.0e-8);

}
