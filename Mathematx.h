#include<cmath> 
#include<limits>

#ifndef MATHEMATX_H_
#define MATHEMATX_H_ 

using std::numeric_limits; 

const long double PI = acos(-1.0) ; 
const long double mach_eps = numeric_limits<double>::epsilon() ; 



bool diff(long double a, long double b, long double eps = mach_eps ) { 
	return (fabs(a-b) > eps) ; 
}

bool close(double a , double b, double eps = mach_eps)  {
	return !diff(a, b , eps) ; 
}

long double atan_proper(long double x, long double y) { 
	long double theta ; 
	if ( close(x, 0) ) {
		if (y > 0) theta =  PI/2.0 ; 
		else if (y < 0) theta =  3*PI/2.0 ; 

	} else if (close(y,0)) {
		if (x > 0 ) theta = 0.0 ; 
		else if ( x < 0) theta = PI ; 

	} else { // both x and y nonzero 
		theta = atan(y/x ) ; 
		if(x > 0 and y > 0 ) {
			theta = theta ; 
		} else if (x > 0 && y< 0) {
			theta = theta + 2*PI ; 
		} else {
			theta = theta + PI ; 
		}
	}
	return theta ; 
}

int sgn(long double val) { 
	return (val > 0 ) - (val < 0) ; 
}

long double roundf(long double val , long double precision = 15){

	return ceil(val*pow(10,precision))/pow(10,precision) ; 
}

long double  isqrt(long double  x) { return 1./sqrt(x) ;  }

#endif 
