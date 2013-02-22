#include <iostream>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/math/common_factor_rt.hpp>
#include <cmath>


 using namespace boost::numeric::ublas;


int main() {
	
	/*    We will first implement the power method for finding
	 *    the largest eigenvalue
	 * 
	 *    
	 * 
	 */
		vector<double> xnp(3);
		vector<double> xn(3);
      xn(0) = 1; xn(1) = 1; xn(2)=1;
			double mod=sqrt(xn(0)*xn(0)+xn(1)*xn(1)+xn(2)*xn(2));
			xn=xn/mod;
			double scaling=0;
			
      matrix<double> A(3,3);
      A(0,0) = 1; A(0,1) = 2; A(0,2) = 0;
      A(1,0) = -2; A(1,1) = 1; A(1,2) = 2;
      A(2,0) = 1; A(2,1) = 3; A(2,2) = 1;
      
			std::cout << xn << std::endl;
			
			for (int i=0;i<=7;i++) {
				xnp = prod(A, xn);
				scaling=xnp(0);
				
				for (int j=1;j<=2;j++) {
					if ( fabs(xnp(j))>fabs(scaling) ) scaling = xnp(j);
				}
				
				//mod=sqrt(xnp(0)*xnp(0)+xnp(1)*xnp(1)+xnp(2)*xnp(2));
				//xnp=xnp/mod;
				xnp=xnp/scaling;
				
				
				
				
				std::cout << " M " << xn << " = " << xnp <<  " approx lambda = " << scaling << std::endl;
      	
				
				
				xn=xnp;
			}
			
			/* Rayleigh Quotiant for finding eigenvalue of dominant
			 * eigenvector.
			 * 
			 * Since A x= lambda x   -->    A (x.x) = lambda x.x
			 * 
			 * Therefore A (x.x)/ (x.x) = lamdba
			 * 
			 * So lambda = (Ax).x / (x.x)
			 */
			 
			 double lambda = inner_prod(prod(A,xnp),xnp)/inner_prod(xnp,xnp);
			 
			
			std::cout << "lambda = " << lambda << std::endl;
			
			
			
			
 
 
      return 0;
	  
	
	
	return 0;
	
}
