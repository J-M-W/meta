#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double mueller(double* r, int deriv, int derivCoord)
{
    /* Evaluates the deriv-th derivative of the Mueller potential with respect to x (if derivCoord is 0) or y (if derivCoord is 1) */
	int i;
	double V=0, expArg=0;

	double A[] = {-200,-100,-170,15};
	double a[]  = {-1, -1, -6.5, 0.7};
	double b[]  = {0, 0, 11, 0.6};
	double c[] = {-10, -10, -6.5, 0.7};
	double x0[] = {1, 0, -0.5, -1};
	double y0[] = {0, 0.5, 1.5, 1};

    if(deriv==0)
    {
        for(i=0; i<4; i++)
        {
            //printf("V = %lf\n", V);
            expArg=0;
            expArg += a[i]*(r[0]-x0[i])*(r[0]-x0[i]);
            expArg += b[i]*(r[0]-x0[i])*(r[1]-y0[i]);
            expArg += c[i]*(r[1]-y0[i])*(r[1]-y0[i]);
            V += A[i]*exp(expArg);
        }
        //printf("V = %lf\n", V);
        return V;
    }
    else if(deriv==1)
    {
        /* dV/dx */
        if(derivCoord==0)
        {
            for(i=0; i<4; i++)
            {
                expArg=0;
                expArg += a[i]*(r[0]-x0[i])*(r[0]-x0[i]);
                expArg += b[i]*(r[0]-x0[i])*(r[1]-y0[i]);
                expArg += c[i]*(r[1]-y0[i])*(r[1]-y0[i]);
                V += A[i]*exp(expArg)*(2*a[i]*(r[0]-x0[i]) + b[i]*(r[1]-y0[i]));
                return V;
            }
        }
        /* dV/dy */
        else if(derivCoord==1)
        {
            for(i=0; i<4; i++)
            {
                expArg=0;
                expArg += a[i]*(r[0]-x0[i])*(r[0]-x0[i]);
                expArg += b[i]*(r[0]-x0[i])*(r[1]-y0[i]);
                expArg += c[i]*(r[1]-y0[i])*(r[1]-y0[i]);
                V += A[i]*exp(expArg)*(2*c[i]*(r[1]-y0[i]) + b[i]*(r[0]-x0[i]));
                return V;
            }
        }
        else
        {
            printf("derivCoord must be 0 or 1, not %d. Exiting.\n", derivCoord);
            exit(EXIT_FAILURE);
        }
    }
    else
    {
        printf("Analytic Mueller function derivative of order %d implemented. Exiting.\n", deriv);
        exit(EXIT_FAILURE);
    }

}
