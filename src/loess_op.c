#include <R.h>
#include <Rmath.h>

double loess_op( int *nn, 
                 int *ddegree, 
                 int *sspan, 
                 int *m, 
                 int *nn_m, 
                 double *result
                 );
inline double my_abs(double x);

double loess_op( int *nn, 
                 int *ddegree, 
                 int *sspan, 
                 int *m, 
                 int *nn_m, 
                 double *result
) {
   int n, n_m, degree, span, span2, span3, idx, offset, ll, rr, ok;
   int i, j, k, ii, jj;
   double *x, *w, *xw, *x2w, *x3w, max_dist, r, tmp1, tmp2;
   int start, end, tmp;

   n = *(nn);
   n_m = *(nn_m);
   degree = *(ddegree);
   span =  *(sspan);

   x = (double *) R_alloc(span, sizeof(double));
   w = (double *) R_alloc(span, sizeof(double));
   xw = (double *) R_alloc(span, sizeof(double));
   x2w = (double *) R_alloc(span, sizeof(double));
   x3w = (double *) R_alloc(span, sizeof(double));
   
   // variables for storing determinant intermediate values
   double a, b, c, d, e, a1, b1, c1, a11, b11, c11, det;
   
   span3 = span;
   if(span > n) span = n;
   span2 = (span - 1)/2;

   // want to start storing results at index 0, corresponding to the lowest m
   offset = *(m);
      
   // loop through all values of m
   for(i = 0; i < n_m; i++)
   {      
      start = *(m + i) - span2;
      end = *(m + i) + span2;
      max_dist = (double)span2;
      if(span == n) {
         start = 1;
         end = n;
         if(*(m + i) - 1 > n - *(m + i)) {
            max_dist = (double)(*(m + i) - 1);
         } else {
            max_dist = (double)(n - *(m + i));
         }
         max_dist = max_dist + (double)(span3-n)/2.0;
         // Rprintf("%f\n", max_dist);
      } else if(start < 1) {
         end = end - start + 1;
         start = 1;
         max_dist = (double)(end - *(m + i));
   	} else if(end > n) {
         start = start - (end - n);
         end = n;
         max_dist = (double)(*(m + i) - start);
   	}
      idx = (*(m + i)) - offset;

      // Rprintf("%d, %d, %d, %d, %f, %d\n", n, span, start, end, max_dist, *(m + i));

      // here's how to speed up the middle part:
      // this only works if it has been previously calculated
      // looks like it really doesn't make that big of a difference
      // if(start == 2) {
      //    tmp = *(m + i);
      //    while(*(m + i) <= n - tmp) {
      //       for(j=0; j < span; j++) {
      //          *(result + idx*n + start - 1 + j) = *(result + (idx - 1)*n + start - 2 + j);    
      //       }
      //       start = start + 1;
      //       idx = idx + 1;
      //       i = i + 1;
      //       // Rprintf("%d, %d, %d\n", start, idx, i);
      //    }
      //    end = start + span;
      // }

      a = 0.0;

      // get weights, x, and a
      for(j = 0; j < span; j++) {
         *(w + j) = 0.0;
         *(x + j) = (double)(start - *(m + i) + j);

         r = my_abs(*(x + j));
         // tricube - hard coded instead of using pow()
         tmp1 = r/max_dist;
         tmp2 = 1.0 - tmp1*tmp1*tmp1;
         *(w + j) = tmp2 * tmp2 * tmp2;                  

         a = a + *(w + j);
      }
      if(degree==0) {
         //TODO: make sure a is not zero
         a1 = 1/a;
         for(j=0; j < span; j++) {
            *(result + idx*n + start - 1 + j) = *(w + j) * a1;    
         }
      } else {
      
         // get xw, x2w, b, c for degree 1 or 2
         b = 0.0;
         c = 0.0;
         for(j = 0; j < span; j++) {
            *(xw + j) = *(x + j) * *(w + j);
            *(x2w + j) = *(x + j) * *(xw + j);
            b = b + *(xw + j);
            c = c + *(x2w + j);
         }
         if(degree==1) {
            // TODO: make sure this is 0
            det = 1/(a * c - b * b);
            a1 = c * det;
            b1 = -b * det;
            for(j=0; j < span; j++) {
               *(result + idx*n + start - 1 + j) = *(w + j) * a1 + *(xw + j) * b1;
            }
         } else {
            // TODO: make sure degree>2 cannot be specified (and <0 for that matter)
            // get x3w, d, and e for degree 2
            d = 0.0;
            e = 0.0;
            for(j = 0; j < span; j++) {
               *(x3w + j) = *(x + j) * *(x2w + j);
               d = d + *(x3w + j);
               e = e + *(x3w + j) * *(x + j);
            }
            a1 = e * c - d * d;
            b1 = c * d - e * b;
            c1 = b * d - c * c;
            // TODO: make sure this is 0
            det = 1/(a*a1 + b*b1 + c*c1);
            a11 = a1*det;
            b11 = b1*det;
            c11 = c1*det;
            for(j=0; j < span; j++) {
               *(result + idx*n + start - 1 + j) = *(w + j) * a11 + *(xw + j) * b11 + *(x2w + j) * c11;
            }
         }
      }
   }
}      

inline double my_abs(double x)
{
return (x > 0) ? x : -x;
}

