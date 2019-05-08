/*
   Copy from matlab expint.m
   yang.yang@nanosic.com
*/
#include<float.h>
#include<math.h>
#define Inf DBL_MAX
double polyval(double x, double*coeff,int len);
double  eps(double val);
double expint(double x){
//function y = expint(x)
//%EXPINT Exponential integral function.
//%   Y = EXPINT(X) is the exponential integral function for each
//%   element of X.  The exponential integral is defined as:
//%
//%   EXPINT(x) = integral from x to Inf of (exp(-t)/t) dt, for x > 0.
//%   
//%   By analytic continuation, EXPINT is a scalar-valued function in
//%   the complex plane cut along the negative real axis.
//%
//%   Another common definition of the exponential integral function is
//%   the Cauchy principal value integral from -Inf to X of (exp(t)/t)
//%   dt, for positive X.  This is denoted as Ei(x). The relationships
//%   between EXPINT(x) and Ei(x) are as follows:
//%
//%       EXPINT(-x+i*0) = -Ei(x) - i*pi, for real x > 0
//%       Ei(x) = REAL(-EXPINT(-x)), for real x > 0
//%
//%   Class support for input X:
//%      float: double, single
//
//%   Copyright 1984-2011 The MathWorks, Inc. 
//
//%       For elements of X in [-38,2], EXPINT uses a series expansion
//%       representation (equation 5.1.11 from Abramowitz & Stegun).
//%       For all other elements of X, EXPINT uses a continued fraction
//%       representation (equation 5.1.22 in A&S).
//%       References:
//%         [1] M. Abramowitz and I.A. Stegun, "Handbook of Mathematical
//%         Functions", Dover Publications, 1965, Ch. 5.
//
//siz = size(x);
//y = zeros(numel(x),1,superiorfloat(x));
//
//% make input a column vector
//x = x(:);
//
//% figure out which algorithm to use by evaluating interpolating polynomial 
//% at real(z)
//
//p = [-3.602693626336023e-09 -4.819538452140960e-07 -2.569498322115933e-05 ...
//     -6.973790859534190e-04 -1.019573529845792e-02 -7.811863559248197e-02 ...
//     -3.012432892762715e-01 -7.773807325735529e-01  8.267661952366478e+00];

   static double p[9] =  {-3.602693626336023e-09, -4.819538452140960e-07, -2.569498322115933e-05, 
                       -6.973790859534190e-04, -1.019573529845792e-02, -7.811863559248197e-02,
                       -3.012432892762715e-01, -7.773807325735529e-01,  8.267661952366478e+00};
    int n; double am2,am1,bm2,bm1; double alpha,a ,al,b,beta;
     double polyv,y;
     double xk, yk, pterm, term,f,oldf ; 
     int j;
//polyv = polyval(p,real(x));

     polyv = polyval(x, p, 9);
//
//% series expansion
//
//k = find( abs(imag(x)) <= polyv );
   
//if ~isempty(k)
//
//   %initialization
//   egamma=0.57721566490153286061;
//   xk = x(k);
//   yk = -egamma - log(xk);
//   j = 1;
//   pterm = xk;
//   term = xk;
//
//   while any(abs(term) > (eps(yk)))
//      yk = yk + term;
//      j = j + 1;
//      pterm = -xk.*pterm/j;
//      term = pterm/j;
//   end % end of the while loop
// 
//   y(k) = yk;
//end
   if(0<=polyv){
       double egamma = 0.57721566490153286061;
       
       xk = x;
       yk = -egamma - log(xk);
	   j = 1;
       pterm = xk;
       term = xk;
 
     while(fabs(term)> eps(yk)){
         yk = yk+term;
         j=j+1;
         pterm = -xk*pterm/j;
         term = pterm/j;
     }

       y = yk;
   }
 
//
//% continued fraction
//
//k = find( abs(imag(x)) > polyv );
//if ~isempty(k)
//   %   note: am1, bm1 corresponds to A(j-1), B(j-1) of recursion formulae
//   %         am2, bm2 corresponds to A(j-2), B(j-2) of recursion formulae
//   %         a,b      corresponds to A(j), B(j) of recursion formulae
//
//   n = 1; % we're calculating E1(x)
//
//   % initialization
//   xk = x(k);
//   am2 = zeros(size(xk));
//   bm2 = ones(size(xk));
//   am1 = ones(size(xk));
//   bm1 = xk;
//   f = am1 ./ bm1;
//   oldf = Inf(size(xk));
//   j = 2;
//
//   while any(abs(f-oldf) > (100*eps(class(f)).*abs(f)))
//       % calculate the coefficients of the recursion formulas for j even
//       alpha = n-1+(j/2); % note: beta= 1
//   
//       %calculate A(j), B(j), and f(j)
//       a = am1 + alpha * am2;
//       b = bm1 + alpha * bm2;
//   
//       % save new normalized variables for next pass through the loop
//       %  note: normalization to avoid overflow or underflow
//       am2 = am1 ./ b;
//       bm2 = bm1 ./ b;
//       am1 = a ./ b;
//       bm1 = 1;
//   
//       f = am1;
//       j = j+1;
//   
//       % calculate the coefficients for j odd
//       alpha = (j-1)/2;
//       beta = xk;
//       a = beta .* am1 + alpha * am2;
//       b = beta .* bm1 + alpha * bm2;
//       am2 = am1 ./ b;
//       bm2 = bm1 ./ b;
//       am1 = a ./ b;
//       bm1 = 1;
//       oldf = f;
//       f = am1;
//       j = j+1;
//   
//   end  % end of the while loop
//    
//   y(k)= exp(-xk) .* f - 1i*pi*((real(xk)<0)&(imag(xk)==0)); 
//end
  
   if(polyv<0){
       n = 1; 
	   xk = x;
	   am2 = 0; bm2 = 0;  am1 = 1;  
       bm1 = xk;  f = am1 / bm1;
       oldf = Inf;
       j = 2;

       if(abs(f-oldf)> (100*2.2204e-16 *abs(f))){  //  eps(class(f)) = eps('double') 浮点相对精度

         
            //%calculate A
             alpha = n-1+(j/2); //% note: beta= 1
             //%calculate A(j), B(j), and f(j)
             a = am1 + alpha * am2;
             b = bm1 + alpha * bm2;
             //% save new normalized variables for next pass through the loop
             //%  note: normalization to avoid overflow or underflow
             am2 = am1 / b;
             bm2 = bm1 / b;
             am1 = a / b;
             bm1 = 1;
             f = am1;
             j = j+1;
                          // calculate 
            // calculate the coefficients for j odd
             alpha = (j-1)/2;
             beta = xk;
             a = beta * am1 + alpha * am2;
             b = beta * bm1 + alpha * bm2;
             am2 = am1 / b;
             bm2 = bm1 / b;
             am1 = a / b;
             bm1 = 1;
             oldf = f;
             f = am1;
             j = j+1;
       }

       y = exp(-xk) * f;// - 1i*pi*(((xk)<0)&(imag(xk)==0)); 


       
   }


//
//y = reshape(y,siz);
    return y;
}

/*
   polyval from matlab function 'polyval'.

   y = x^8 + x^7 + x^6+x^5+x^4+x^3+x^2+x^1 + 1;
*/

double polyval(double x, double*coeff,int len){
    int i,j;
    double temp1, temp2, rev;
	double x_2 =  x*x;
    temp1 = x;	
    temp2 = x*x;

	rev = coeff[8];

    rev += temp1 * coeff[7] + temp2 *coeff[6];  temp1 *=x_2; temp2 *=x_2;
    rev += temp1 * coeff[5] + temp2 *coeff[4];  temp1 *=x_2; temp2 *=x_2;
    rev += temp1 * coeff[3] + temp2 *coeff[2];  temp1 *=x_2; temp2 *=x_2;
    rev += temp1 * coeff[1] + temp2 *coeff[0];
    return rev;
}


#define eps1  2.2204e-16
double  eps(double val){
   
    double rev,temp;
    temp = 1.0;
    val = fabs(val);

    rev = eps1;

    if(val>1.0){         
         while(temp> val){
             temp*=2;
             rev*=2;
         }
    }
    else{
          while( temp>val){
             temp/=2;
             rev/=2;
         }        
    }

 
	return rev;

}