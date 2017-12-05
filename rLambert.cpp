// rLambert.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "stdlib.h"
#include "math.h"
#include "cmath"

double LambertW(const double z) {
  int i; 
  const double eps=4.0e-16, em1=0.3678794411714423215955237701614608; 
  double p,e,t,w;
  if (z<-em1) { 
    fprintf(stderr,"LambertW: bad argument %g, exiting.\n",z); return 0; 
  }
  if (0.0==z) return 0.0;
  if (z<-em1+1e-4) { // series near -em1 in sqrt(q)
    double q=z+em1,r=sqrt(q),q2=q*q,q3=q2*q;
    return 
     -1.0
     +2.331643981597124203363536062168*r
     -1.812187885639363490240191647568*q
     +1.936631114492359755363277457668*r*q
     -2.353551201881614516821543561516*q2
     +3.066858901050631912893148922704*r*q2
     -4.175335600258177138854984177460*q3
     +5.858023729874774148815053846119*r*q3
     -8.401032217523977370984161688514*q3*q;  // error approx 1e-16
  }
  /* initial approx for iteration... */
  if (z<1.0) { /* series near 0 */
    p=sqrt(2.0*(2.7182818284590452353602874713526625*z+1.0));
    w=-1.0+p*(1.0+p*(-0.333333333333333333333+p*0.152777777777777777777777)); 
  } else 
    w=log(z)-log(log(z)); /* asymptotic */
  for (i=0; i<10; i++) { /* Halley iteration */
    e=exp(w); 
    t=w*e-z;
    p=w+1.0;
    t/=e*p-0.5*(p+1.0)*t/p; 
    w-=t;
    if (fabs(t)<eps*(1.0+fabs(w))) return w; /* rel-abs error */
  }
  /* should never get here */
  fprintf(stderr,"LambertW: No convergence at z=%g, exiting.\n",z); 
  return 0;
}
double LambertW1(const double z) {
  int i; 
  const double eps=4.0e-16, em1=0.3678794411714423215955237701614608; 
  double p=1.0,e,t,w,l1,l2;
  if (z<-em1 || z>=0.0) { 
    fprintf(stderr,"LambertW1: bad argument %g, exiting.\n",z); return 0; 
  }
  /* initial approx for iteration... */
  if (z<-1e-6) { /* series about -1/e */
    p=-sqrt(2.0*(2.7182818284590452353602874713526625*z+1.0));
    w=-1.0+p*(1.0+p*(-0.333333333333333333333+p*0.152777777777777777777777)); 
  } else { /* asymptotic near zero */
    l1=log(-z);
    l2=log(-l1);
    w=l1-l2+l2/l1;
  }
  if (fabs(p)<1e-4) return w;
  for (i=0; i<10; i++) { /* Halley iteration */
    e=exp(w); 
    t=w*e-z;
    p=w+1.0;
    t/=e*p-0.5*(p+1.0)*t/p; 
    w-=t;
    if (fabs(t)<eps*(1.0+fabs(w))) return w; /* rel-abs error */
  }
  /* should never get here */
  fprintf(stderr,"LambertW1: No convergence at z=%g, exiting.\n",z);
  return 0;
}
double F ( const double x, const double r )
{
	return ( x * exp( x ) + r * x );
}

double DF ( const double x, const double r )
{
	return ( (1. + x) * exp( x ) + r );
}

double DDF ( const double x )
{
	return ( (2. + x) * exp( x ) );
}


double rLambert( const double x, const double r, const char prec = 6 )
{
	if ( x == 0. ) return 0.; // W(x,r)=0 always
	if ( r == 0. )
	{
		if ( x < 0 ) { printf( "%s\n", "There is one additional branch:" );
		printf( "W-1(%f, 0) = %.*f\n", x, prec, LambertW1(x) ); }
		return LambertW( x );
	}

	const double e = exp( 1. );
	const double em2 = exp ( -2. );

	double w, wprev; // intermediate values in the Halley method
	const unsigned int maxiter = 20; // max number of iterations. This eliminates problems with convergence
	unsigned int iter; // iteration counter
	

	// If r >= exp(-2), there is just one branch of W(x,r):
	if ( r >= em2 )
	{
		if ( (r == em2) && (x == -4. * em2) ) return -2.; // At this x W(x,em2) is not differentiable, so
		// Halley's method does not work. But it can be calculated that W(x,em2) = -2.
		// If x is close but not equal to -4*em2, there is no problem.
		// For example, if x= -4. * exp( -2. ) +/- .00000000001,
		// the program still gives the correct result

		// begin the iteration up to prec precision
		// initial value for the Halley method
		if ( x > 1. )  w = log( x ) - log( log( x ) );
		if ( x < -1.) w = 1. / r * x;
		if ( abs(x) <= 1 ) w = 0.;
		wprev = w;

		iter = 0;
		do
		{
			wprev = w;
			w -= 2.*((F(w,r)-x) * DF(w,r)) /
				(2.*DF(w,r)*DF(w,r) - (F(w,r)-x)*DDF(w));
			iter++;
		} while ( (abs( w - wprev ) > pow( 10., -prec )) && iter < maxiter );
		return w;
	} // if ( r >= em2 )
	
	//here comes the calculation when W(x,r) has three branches (i.e. 0 < r < e^(-2))
	else if ( 0 < r && r < em2 )
	{
		double alpha  = LambertW1( -r * e ) - 1.; // left  branch point
		double beta   = LambertW ( -r * e ) - 1.; // right branch point
		double falpha = F( alpha, r );
		double fbeta  = F( beta,  r );

		if ( x < fbeta ) // leftmost branch
		{
			if ( x < -40.54374295204823 ) return x / r; // because x*e^x < 1E-16 as x < -40.5437...
			wprev = w = x / r;

			iter = 0;
			do
			{
				wprev = w;
				w -= 2.*((F(w,r)-x) * DF(w,r)) /
					(2.*DF(w,r)*DF(w,r) - (F(w,r)-x)*DDF(w));
				iter++;
			} while ( (abs( w - wprev ) > pow( 10., -prec )) && iter < maxiter );
			return w;
		} // if ( x < fbeta )

		if ( x >= fbeta && x <= falpha ) // leftmost and inner branches
		{
			double wm2 = x / r; // initial value for the leftmost branch W-2(x,r)
			double wm1 = -3.; // initial value for the inner branch W-1(x,r)
			double wm0 = -1.;  // initial value for the rightmost branch W0(x,r)
			
			// Halley iteration for the leftmost branch
			iter = 0;
			do
			{
				wprev = wm2;
				wm2 -= 2.*((F(wm2,r)-x) * DF(wm2,r)) /
					(2.*DF(wm2,r)*DF(wm2,r) - (F(wm2,r)-x)*DDF(wm2));
				iter++;
			} while ( (abs( wm2 - wprev ) > pow( 10., -prec )) && iter < maxiter );
			printf( "%s\n", "There are two additional branches:" );
			printf( "W-2(%f, %f) = %.*f\n", x, r, prec, wm2 );

			// Halley iteration for the inner branch
			iter = 0;
			do
			{
				wprev = wm1;
				wm1 -= 2.*((F(wm1,r)-x) * DF(wm1,r)) /
					(2.*DF(wm1,r)*DF(wm1,r) - (F(wm1,r)-x)*DDF(wm1));
				iter++;
			} while ( (abs( wm1 - wprev ) > pow( 10., -prec )) && iter < maxiter );
			printf( "W-1(%f, %f) = %.*f\n", x, r, prec, wm1 );

			// Halley iteration for the rightmost branch
			iter = 0;
			do
			{
				wprev = wm0;
				wm0 -= 2.*((F(wm0,r)-x) * DF(wm0,r)) /
					(2.*DF(wm0,r)*DF(wm0,r) - (F(wm0,r)-x)*DDF(wm0));
				iter++;
			} while ( (abs( wm0 - wprev ) > pow( 10., -prec )) && iter < maxiter );
			return wm0;
		} // if ( x >= fbeta && x <= falpha ) 

		if ( x > falpha ) // rightmost branch
		{
			wprev = w = ( x > 1. ) ? log ( x ) - log( log ( x ) ) : 0.; // initial value

			iter = 0;
			do
			{
				wprev = w;
				w -= 2.*((F(w,r)-x) * DF(w,r)) /
					(2.*DF(w,r)*DF(w,r) - (F(w,r)-x)*DDF(w));
				iter++;
			} while ( (abs( w - wprev ) > pow( 10., -prec )) && iter < maxiter );
			return w;
		} // if ( x > falpha )
			
	} // else if
	else if ( r < 0 ) // two branches separated by W(-r*e)-1
	{
		// minimum of F: W(-r*e)-1, zeros: 0 and log(-r)
		if ( x < F( LambertW( -r * e ) - 1., r ) )
		{
			printf( "%s\n", "First argument of LambertW(x,r) is out of domain" );
			return NULL;
		}
		if ( x == log( -r ) ) return 0;
		
		if ( x < 0 ) // Two initial values, one less than the minimum of f, the other is one greater
		{
			// left branch
			wprev = w = LambertW( -r * e ) - 2;
	
			iter = 0;
			do
			{
				wprev = w;
				w -= 2.*((F(w,r)-x) * DF(w,r)) /
					(2.*DF(w,r)*DF(w,r) - (F(w,r)-x)*DDF(w));
				iter++;
			} while ( (abs( w - wprev ) > pow( 10., -prec )) && iter < maxiter );
			printf( "%s\n", "There is one additional branch:" );
			printf( "W-1(%f, %f) = %.*f\n", x, r, prec, w );

			// right branch
			wprev = w = LambertW( -r * e );
			iter = 0;
			do
			{
				wprev = w;
				w -= 2.*((F(w,r)-x) * DF(w,r)) /
					(2.*DF(w,r)*DF(w,r) - (F(w,r)-x)*DDF(w));
				iter++;
			} while ( (abs( w - wprev ) > pow( 10., -prec )) && iter < maxiter );
			return w;
			}

		if ( x > 0 )
		{
			double lzero = r > -1 ? log(-r) : 0;
			double rzero = r > -1 ? 0 : log(-r);

			// left branch
			wprev = w = lzero - 1;
			iter = 0;
			do
			{
				wprev = w;
				w -= 2.*((F(w,r)-x) * DF(w,r)) /
					(2.*DF(w,r)*DF(w,r) - (F(w,r)-x)*DDF(w));
				iter++;
			} while ( (abs( w - wprev ) > pow( 10., -prec )) && iter < maxiter );
			printf( "%s\n", "There is one additional branch:" );
			printf( "W-1(%f, %f) = %.*f\n", x, r, prec, w );

			// right branch
			wprev = w = rzero + 1;
			iter = 0;
			do
			{
				wprev = w;
				w -= 2.*((F(w,r)-x) * DF(w,r)) /
					(2.*DF(w,r)*DF(w,r) - (F(w,r)-x)*DDF(w));
				iter++;
			} while ( (abs( w - wprev ) > pow( 10., -prec )) && iter < maxiter );
			return w;
		} // if ( x < 0 )
		} // else if ( r < 0 )
	/* should never get here */
	fprintf(stderr, "LambertW1: No convergence at x=%g, exiting.\n", x);
	return 0;
} // rLambertW

int main()
{
	char prec = 9;

	double r = 2.7;
	double x = 1.35;

	printf( "W(%6.4f, %6.4f)= %.*f\n", x, r, prec, rLambert( x, r, prec ) );

	return 0;
}

