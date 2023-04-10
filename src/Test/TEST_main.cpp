#include<iostream>
#include"gtest/gtest.h"
#include"Three_Moments_test.cpp"
#include"Five_Moments_test.cpp"
#include"Seven_Moments_test.cpp"
#include"Integration_test.cpp"
#include"Searching_Algorithm_test.cpp"
#include"Distribution_Function_test.cpp"

int main(int argc, char* argv[]) {

  /////////////////////////////////////////////////////////////////////
  // Initialize google test
  testing::InitGoogleTest(&argc, argv);

  return RUN_ALL_TESTS();
}
