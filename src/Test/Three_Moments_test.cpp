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
/// Test for the weight vector
TEST(Three_Moments, W){

  double v  = 6.3;
  double W0 = 1.0;
  double W1 = v;
  double W2 = v*v;
  
  auto W = homedf::Three_Moments::W(v);

  EXPECT_NEAR(W0, W(0), 1.0e-10);
  EXPECT_NEAR(W1, W(1), 1.0e-10);
  EXPECT_NEAR(W2, W(2), 1.0e-10);
  
}

/////////////////////////////////////////////////////////////////////////
/// Test for V -> U conversion
TEST(Three_Moments, V_to_U){

  homedf::Three_Moments::Vector_type V;
  V(0) = 1.3;
  V(1) = 12.5;
  V(2) = 2.65;

  homedf::Three_Moments::Vector_type U_EXACT;
  U_EXACT(0) = 1.3;
  U_EXACT(1) = 16.25;
  U_EXACT(2) = 205.775;
  
  auto U = homedf::Three_Moments::Get_Solution_Vector_from_Primitive(V);

  EXPECT_NEAR(U_EXACT(0), U(0), 1.0e-10);
  EXPECT_NEAR(U_EXACT(1), U(1), 1.0e-10);
  EXPECT_NEAR(U_EXACT(2), U(2), 1.0e-10);

}

/////////////////////////////////////////////////////////////////////////
/// Test for U -> V conversion
TEST(Three_Moments, U_to_V){

  homedf::Three_Moments::Vector_type U;
  U(0) = 1.3;
  U(1) = 16.25;
  U(2) = 205.775;

  homedf::Three_Moments::Vector_type V_EXACT;
  V_EXACT(0) = 1.3;
  V_EXACT(1) = 12.5;
  V_EXACT(2) = 2.65;

  auto V = homedf::Three_Moments::Get_Primitive_from_Solution_Vector(U);

  EXPECT_NEAR(V_EXACT(0), V(0), 1.0e-10);
  EXPECT_NEAR(V_EXACT(1), V(1), 1.0e-10);
  EXPECT_NEAR(V_EXACT(2), V(2), 1.0e-10);

}

/////////////////////////////////////////////////////////////////////////
/// Test for non-dimensionalization of U
TEST(Three_Moments, U_ND){

  homedf::Three_Moments::Vector_type U;
  U(0) = 1.3;
  U(1) = 16.25;
  U(2) = 205.775;

  homedf::Three_Moments::Vector_type U_ND_EXACT;
  U_ND_EXACT(0) = 1.0;
  U_ND_EXACT(1) = 0.9935400628039;
  U_ND_EXACT(2) = 1.0;

  auto U_ND = homedf::Three_Moments::Non_Dimensionalize_Solution_Vector(U);

  EXPECT_NEAR(U_ND_EXACT(0), U_ND(0), 1.0e-10);
  EXPECT_NEAR(U_ND_EXACT(1), U_ND(1), 1.0e-10);
  EXPECT_NEAR(U_ND_EXACT(2), U_ND(2), 1.0e-10);

}

/////////////////////////////////////////////////////////////////////////
/// Test for non-dimensionalization of V
TEST(Three_Moments, V_ND){

  homedf::Three_Moments::Vector_type V;
  V(0) = 1.2;
  V(1) = 0.15;
  V(2) = 1.1;

  homedf::Three_Moments::Vector_type V_ND_EXACT;
  V_ND_EXACT(0) = 1.0;
  V_ND_EXACT(1) = 0.1547818111;
  V_ND_EXACT(2) = 0.9760425909;

  auto V_ND = homedf::Three_Moments::Non_Dimensionalize_Primitive_Vector(V);

  EXPECT_NEAR(V_ND_EXACT(0), V_ND(0), 1.0e-10);
  EXPECT_NEAR(V_ND_EXACT(1), V_ND(1), 1.0e-10);
  EXPECT_NEAR(V_ND_EXACT(2), V_ND(2), 1.0e-10);

}

/////////////////////////////////////////////////////////////////////////
/// Test for dimensionalization of alphas
TEST(Three_Moments, alpha){

  homedf::Three_Moments::Vector_type U;
  U(0) = 1.3;
  U(1) = 16.25;
  U(2) = 205.775;

  homedf::Three_Moments::Vector_type a_ND;
  a_ND(0) = -37.06829838;
  a_ND(1) = 77.14932318;
  a_ND(2) = -38.8254717;

  homedf::Three_Moments::Vector_type a_EXACT;
  a_EXACT(0) = -39.33814365;
  a_EXACT(1) = 6.132075472;
  a_EXACT(2) = -0.2452830189;

  homedf::Three_Moments Three_Moments_Object;
  auto a = Three_Moments_Object.Dimensionalize_Alphas(U, a_ND);
  
  EXPECT_NEAR(a_EXACT(0), a(0), 1.0e-8);
  EXPECT_NEAR(a_EXACT(1), a(1), 1.0e-8);
  EXPECT_NEAR(a_EXACT(2), a(2), 1.0e-8);

}
