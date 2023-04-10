#ifndef INTEGRATION_H
#define INTEGRATION_H
/////////////////////////////////////////////////////////////////////////
#include<array>
#include<cmath>

namespace homedf {

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
template<int quad_points, typename Point_type>
class Gaussian_Quadrature_Rule;

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
template<typename Point_type>
class Gaussian_Quadrature_Rule<1, Point_type> {
public:
  Gaussian_Quadrature_Rule(const Point_type& dp) :
    weight(dp) {}

  template<typename Func_type>
  auto operator()(const Func_type& f, const Point_type& p) const -> decltype(f(p)) {
    return f(p)*weight;
  }

private:
  const Point_type weight;
};

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
template<typename Point_type>
class Gaussian_Quadrature_Rule<2, Point_type> {
public:
  Gaussian_Quadrature_Rule(const Point_type& dp) :
    offset(dp*0.5/sqrt(3.0)),
    weight(dp*0.5) {}

  template<typename Func_type>
  auto operator()(const Func_type& f, const Point_type& p) const -> decltype(f(p)) {
    return (f(p-offset)+f(p+offset))*weight;
  }

private:
  const Point_type offset;
  const Point_type weight;
};

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
template<typename Point_type>
class Gaussian_Quadrature_Rule<3, Point_type> {
public:
  Gaussian_Quadrature_Rule(const Point_type& dp) :
    offset1(dp*0.5*sqrt(3.0/5.0)),
    weight0(dp*4.0/9.0),
    weight1(dp*5.0/18.0) {}

  template<typename Func_type>
  auto operator()(const Func_type& f, const Point_type& p) const -> decltype(f(p)) {
    return f(p)*weight0 + (f(p-offset1)+f(p+offset1))*weight1;
  }

private:
  const Point_type offset1;
  const Point_type weight0;
  const Point_type weight1;
};

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
template<typename Point_type>
class Gaussian_Quadrature_Rule<4, Point_type> {
public:

  Gaussian_Quadrature_Rule(const Point_type& dp) :
    offset0(dp*0.5*sqrt((3.0/7.0)-(2.0/7.0)*sqrt(6.0/5.0))),
    offset1(dp*0.5*sqrt((3.0/7.0)+(2.0/7.0)*sqrt(6.0/5.0))),
    weight0(dp*(18.0+sqrt(30.0))/72.0),
    weight1(dp*(18.0-sqrt(30.0))/72.0) {}

  template<typename Func_type>
  auto operator()(const Func_type& f, const Point_type& p) const -> decltype(f(p)) {
    return (f(p-offset1)+f(p+offset1))*weight1 + (f(p-offset0)+f(p+offset0))*weight0;
  }

private:
  const Point_type offset0;
  const Point_type offset1;
  const Point_type weight0;
  const Point_type weight1;
};

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
template<typename Point_type>
class Gaussian_Quadrature_Rule<5, Point_type> {
public:

  Gaussian_Quadrature_Rule(const Point_type& dp) :
    offset1(dp*sqrt(5.0-2.0*sqrt(10.0/7.0))/6.0),
    offset2(dp*sqrt(5.0+2.0*sqrt(10.0/7.0))/6.0),
    weight0(dp*64.0/225.0),
    weight1(dp*(322.0+13.0*sqrt(70.0))/1800.0),
    weight2(dp*(322.0-13.0*sqrt(70.0))/1800.0) {}

  template<typename Func_type>
  auto operator()(const Func_type& f, const Point_type& p) const -> decltype(f(p)) {
    return (f(p-offset2)+f(p+offset2))*weight2
         + (f(p-offset1)+f(p+offset1))*weight1
         + f(p)*weight0;
  }

private:
  const Point_type offset1;
  const Point_type offset2;
  const Point_type weight0;
  const Point_type weight1;
  const Point_type weight2;
};
 
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
template<int quad_points, typename Func_type, typename Point_type>
auto gaussian_integrate(const Func_type& func,
			const Point_type& p0,
			const Point_type& p1,
			const int& num_points) {
  const auto delta_p = (p1-p0)/num_points;
  
  const auto quad_rule = Gaussian_Quadrature_Rule<quad_points, Point_type>(delta_p);

  decltype(func(p0)) sum = quad_rule(func,p0 + 0.5*delta_p)-quad_rule(func,p0 + 0.5*delta_p);
  
  int num_thread = omp_get_max_threads();
  
#pragma omp parallel
  {
    auto my_sum = sum;
    auto p = p0;
    
#pragma omp for
    for(int i = 0; i < num_points; ++i) {
      p   = (i+0.5)*delta_p + p0;
      my_sum += quad_rule(func,p);
    }

    for(int i = 0; i < num_thread; ++i){
      if(i == omp_get_thread_num()){
	sum += my_sum;
      }
#pragma omp barrier
    }    
  }
  
  return sum;
}
//////////////////////////////////////////////////////////////////////////////////////////
  
}// namespace homedf

#endif //#ifndef INTEGRATION_H
