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
TEST(Five_Moments, W){

  double v  = 6.3;
  double W0 = 1.0;
  double W1 = v;
  double W2 = v*v;
  double W3 = v*v*v;
  double W4 = v*v*v*v;
  
  auto W = homedf::Five_Moments::W(v);

  EXPECT_NEAR(W0, W(0), 1.0e-10);
  EXPECT_NEAR(W1, W(1), 1.0e-10);
  EXPECT_NEAR(W2, W(2), 1.0e-10);
  EXPECT_NEAR(W3, W(3), 1.0e-10);
  EXPECT_NEAR(W4, W(4), 1.0e-10);
  
}

/////////////////////////////////////////////////////////////////////////
/// Test for V -> U conversion
TEST(Five_Moments, V_to_U){

  homedf::Five_Moments::Vector_type V;
  V(0) = 1.3;
  V(1) = 12.5;
  V(2) = 2.65;
  V(3) = 1.2;
  V(4) = 6.9;

  homedf::Five_Moments::Vector_type U_EXACT;
  U_EXACT(0) = 1.3;
  U_EXACT(1) = 16.25;
  U_EXACT(2) = 205.775;
  U_EXACT(3) = 2639.63749999999;
  U_EXACT(4) = 34289.55625;
  
  auto U = homedf::Five_Moments::Get_Solution_Vector_from_Primitive(V);

  EXPECT_NEAR(U_EXACT(0), U(0), 1.0e-8);
  EXPECT_NEAR(U_EXACT(1), U(1), 1.0e-8);
  EXPECT_NEAR(U_EXACT(2), U(2), 1.0e-8);
  EXPECT_NEAR(U_EXACT(3), U(3), 1.0e-8);
  EXPECT_NEAR(U_EXACT(4), U(4), 1.0e-8);

}

/////////////////////////////////////////////////////////////////////////
/// Test for U -> V conversion
TEST(Five_Moments, U_to_V){

  homedf::Five_Moments::Vector_type U;
  U(0) = 1.3;
  U(1) = 16.25;
  U(2) = 205.775;
  U(3) = 2639.63749999999;
  U(4) = 34289.55625;
  
  homedf::Five_Moments::Vector_type V_EXACT;
  V_EXACT(0) = 1.3;
  V_EXACT(1) = 12.5;
  V_EXACT(2) = 2.65;
  V_EXACT(3) = 1.2;
  V_EXACT(4) = 6.9;

  auto V = homedf::Five_Moments::Get_Primitive_from_Solution_Vector(U);

  EXPECT_NEAR(V_EXACT(0), V(0), 1.0e-8);
  EXPECT_NEAR(V_EXACT(1), V(1), 1.0e-8);
  EXPECT_NEAR(V_EXACT(2), V(2), 1.0e-8);
  EXPECT_NEAR(V_EXACT(3), V(3), 1.0e-8);
  EXPECT_NEAR(V_EXACT(4), V(4), 1.0e-8);

}

/////////////////////////////////////////////////////////////////////////
/// Test for non-dimensionalization of U
TEST(Five_Moments, U_ND){

  homedf::Five_Moments::Vector_type U;
  U(0) = 1.3;
  U(1) = 16.25;
  U(2) = 205.775;
  U(3) = 2639.63749999999;
  U(4) = 34289.55625;
  
  homedf::Five_Moments::Vector_type U_ND_EXACT;
  U_ND_EXACT(0) = 1.0;
  U_ND_EXACT(1) = 0.9935400628039;
  U_ND_EXACT(2) = 1.0;
  U_ND_EXACT(3) = 1.0195934812410;
  U_ND_EXACT(4) = 1.0527372649315;

  auto U_ND = homedf::Five_Moments::Non_Dimensionalize_Solution_Vector(U);

  EXPECT_NEAR(U_ND_EXACT(0), U_ND(0), 1.0e-8);
  EXPECT_NEAR(U_ND_EXACT(1), U_ND(1), 1.0e-8);
  EXPECT_NEAR(U_ND_EXACT(2), U_ND(2), 1.0e-8);
  EXPECT_NEAR(U_ND_EXACT(3), U_ND(3), 1.0e-8);
  EXPECT_NEAR(U_ND_EXACT(4), U_ND(4), 1.0e-8);

}

/////////////////////////////////////////////////////////////////////////
/// Test for non-dimensionalization of V
TEST(Five_Moments, V_ND){

  homedf::Five_Moments::Vector_type V;
  V(0) = 1.2;
  V(1) = 0.15;
  V(2) = 1.1;
  V(3) = 0.2;
  V(4) = 3.5;
  
  homedf::Five_Moments::Vector_type V_ND_EXACT;
  V_ND_EXACT(0) = 1.0;
  V_ND_EXACT(1) = 0.1547818111;
  V_ND_EXACT(2) = 0.9760425909;
  V_ND_EXACT(3) = 0.1831195636;
  V_ND_EXACT(4) = 3.306750732;

  auto V_ND = homedf::Five_Moments::Non_Dimensionalize_Primitive_Vector(V);

  EXPECT_NEAR(V_ND_EXACT(0), V_ND(0), 1.0e-8);
  EXPECT_NEAR(V_ND_EXACT(1), V_ND(1), 1.0e-8);
  EXPECT_NEAR(V_ND_EXACT(2), V_ND(2), 1.0e-8);
  EXPECT_NEAR(V_ND_EXACT(3), V_ND(3), 1.0e-8);
  EXPECT_NEAR(V_ND_EXACT(4), V_ND(4), 1.0e-8);

}

/////////////////////////////////////////////////////////////////////////
/// Test for dimensionalization of alphas
TEST(Five_Moments, alpha){

  homedf::Five_Moments::Vector_type V;
  V(0) = 1.3;
  V(1) = 1.25;
  V(2) = 2.65;
  V(3) = 1.2;
  V(4) = 6.9;

  auto U = homedf::Five_Moments::Get_Solution_Vector_from_Primitive(V);

  homedf::Five_Moments::Vector_type a_ND;
  a_ND(0) = 0.4304344733;
  a_ND(1) = -0.6065065751;
  a_ND(2) = -23.93543421;
  a_ND(3) = 30.6342280;
  a_ND(4) = -9.704300356;

  homedf::Five_Moments::Vector_type a_EXACT;
  a_EXACT(0) = 0.05219828583;
  a_EXACT(1) = -0.3196143515;
  a_EXACT(2) = -6.646956363;
  a_EXACT(3) = 4.483108159;
  a_EXACT(4) = -0.748388774;

  homedf::Five_Moments Five_Moments_Object;
  auto a = Five_Moments_Object.Dimensionalize_Alphas(U, a_ND);
  
  EXPECT_NEAR(a_EXACT(0), a(0), 1.0e-8);
  EXPECT_NEAR(a_EXACT(1), a(1), 1.0e-8);
  EXPECT_NEAR(a_EXACT(2), a(2), 1.0e-8);
  EXPECT_NEAR(a_EXACT(3), a(3), 1.0e-8);
  EXPECT_NEAR(a_EXACT(4), a(4), 1.0e-8);

}

/////////////////////////////////////////////////////////////////////////
/// Test for sigma -> r
TEST(Five_Moments, sigma_to_r){

  homedf::Five_Moments::Vector_type V;
  V(0) = 1.3;
  V(1) = 1.25;
  V(2) = 2.65;
  V(3) = 1.2;
  V(4) = 6.9;

  double sigma   = 0.91623346459;
  double r_EXACT = 6.9;
  double r       = homedf::Five_Moments::r_from_sigma(V, sigma);
  
  EXPECT_NEAR(r_EXACT, r, 1.0e-8);

}

/////////////////////////////////////////////////////////////////////////
/// Test for r -> sigma
TEST(Five_Moments, r_to_sigma){

  homedf::Five_Moments::Vector_type V;
  V(0) = 1.3;
  V(1) = 1.25;
  V(2) = 2.65;
  V(3) = 1.2;
  V(4) = 6.9;

  double sigma_EXACT = 0.91623346459;
  double sigma = homedf::Five_Moments::sigma_from_r(V);
  
  EXPECT_NEAR(sigma_EXACT, sigma, 1.0e-8);

}
