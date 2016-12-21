/* ============================================================================
* Function for numerically solving Ordinary Differential Equations
*
* - From numerical recipies in C, see Press et al. 1992.
*
* NOTES:
* This is quite a tricky setup. You have to pass the function to the odeint
* function which then passes it down two more levels to the rkqs and rkck funcs.
* To pass any parameters these have to be wrapped up and passed, adjusting all
* the expected parameter lists in the *derivs and *rkqs funcs. Note the order
* is obviously very important
*
* NB.I've refactored this function so that all the allocation/free calls have
* been extracted to speed things up as it was eating 40-50% of the time
* spent. Instead we now pass nrutil which contains all the needed vars
*
*
* AUTHOR:
*   Martin De Kauwe
*
* DATE:
*   23.11.2016
*
* =========================================================================== */


#include <math.h>
#define NRANSI
#include "nrutil.h"
#define MAXSTP 10000
#define TINY 1.0e-30

#include "odeint.h"

/*
I'm declaring these within the function...so this means we can't interface
with the results globally as numerical recipes intended. That was probably
a horrible idea anyway. We should pass them via the odeint func if we want
them. Frankly I can't be arsed passing N more vars, so we're not going to do
that
*/

/*extern int kmax,kount;
extern double *xp,**yp,dxsav;*/



void odeint(double ystart[], int nvar, double x1, double x2, double eps,
            double h1, double hmin, int *nok, int *nbad,
            double aa, double bb, double cc, double dd, double ee, nrutil *nr,
            void (*derivs)(double, double [], double [], double, double,
                           double, double, double),
            void (*rkqs)(double [], double [], int, double *, double, double,
                         double [], double *, double *, double, double, double,
                         double, double, nrutil *nr,
                         void (*)(double, double [], double [], double, double,
                                  double, double, double))) {


	int      kmax, kount;
	int      nstp, i;
	double   dxsav, xsav, x, hnext, hdid, h;
	//double  *yscal=NULL, *y=NULL, *dydx=NULL, *xp=NULL;
	//double **yp=NULL;

	/* initialising this within the func, which means this isn't generic*/
	//kmax = 100;

	//xp = dvector(1, kmax);
	//yp = dmatrix(1,nvar,1,kmax);
	//yscal=dvector(1,nvar);
	//y=dvector(1,nvar);
	//dydx=dvector(1,nvar);


    dxsav = (x2 - x1) / 20.0;
    x=x1;
	h=SIGN(h1, x2-x1);
	*nok = (*nbad) = kount = 0;

    //printf("* %lf\n", ystart[1]);
	for (i=1;i<=nvar;i++) nr->y[i]=ystart[i];
	if (kmax > 0) xsav=x-dxsav*2.0;
	for (nstp=1;nstp<=MAXSTP;nstp++) {
		(*derivs)(x,nr->y,nr->dydx, aa, bb, cc, dd, ee);
		for (i=1;i<=nvar;i++)
			nr->yscal[i]=fabs(nr->y[i])+fabs(nr->dydx[i]*h)+TINY;
		if (kmax > 0 && kount < kmax-1 && fabs(x-xsav) > fabs(dxsav)) {
			nr->xp[++kount]=x;
			for (i=1;i<=nvar;i++) nr->yp[i][kount]=nr->y[i];
			xsav=x;
		}

		if ((x+h-x2)*(x+h-x1) > 0.0) h=x2-x;
		(*rkqs)(nr->y,nr->dydx,nvar,&x,h,eps,nr->yscal,&hdid,&hnext,
                aa, bb, cc, dd, ee, nr,
				derivs);
		if (hdid == h) ++(*nok); else ++(*nbad);
		if ((x-x2)*(x2-x1) >= 0.0) {

			for (i=1;i<=nvar;i++) ystart[i]=nr->y[i];

			if (kmax) {
				nr->xp[++kount]=x;
				for (i=1;i<=nvar;i++) nr->yp[i][kount]=nr->y[i];
			}

			//free_dvector(dydx,1,nvar);
			//free_dvector(y,1,nvar);
			//free_dvector(yscal,1,nvar);
            //free_dvector(xp,1,kmax);
            //free_dmatrix(yp,1,nvar,1,kmax);

            //printf("** %lf\n", ystart[1]);
			return;
		}
		if (fabs(hnext) <= hmin) nrerror("Step size too small in odeint");
		h=hnext;
	}
	nrerror("Too many steps in routine odeint");


}
#undef MAXSTP
#undef TINY
#undef NRANSI
