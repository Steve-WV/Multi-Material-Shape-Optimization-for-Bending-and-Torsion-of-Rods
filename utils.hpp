#ifndef __UTILS_HPP__
#define __UTILS_HPP__

#include <assert.h>
#include <algorithm>
#include <math.h>

//const int NumIntPts = 1;
//const double GaussWeights[NumIntPts] = {1.0};
//const double GaussPoints[NumIntPts][4] = {{1.0/3.0, 1.0/3.0, 1.0/3.0}};

const double pi = 3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628;

const int NumIntPts = 3;
const double GaussWeights[NumIntPts] = { 2.0/6.0, 2.0/6.0, 2.0/6.0 };
const double GaussPoints[NumIntPts][3] = 
	{ {1.0/6.0, 1.0/6.0, 4.0/6.0}, {1.0/6.0, 4.0/6.0, 1.0/6.0}, {4.0/6.0, 1.0/6.0, 1.0/6.0}, };

inline void die() { assert(false); }

template <class type>
inline type
min( type a, type b, type c ) {
  const type m = std::min( a, b );
  return std::min( m, c );
}

template <class type>
inline type
max( type a, type b, type c ) {
  const type m = std::max( a, b );
  return std::max( m, c );
}

template <class type>
inline type
sqr( type a ) {
   return a*a;
}
template <class type>
inline type 
cube(type a) {
return a*a*a;
}


template <class T>
struct Delete : public std::unary_function<T*, void>{
  void operator()( T* x ) const { delete x; }
};



#endif
