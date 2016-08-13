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
*
* AUTHOR:
*   Martin De Kauwe
*
* DATE:
*   11.08.2016
*
* =========================================================================== */


#include <math.h>
#include "gday.h"
#define NRANSI
#include "nrutil.h"
#define MAXSTP 10000
#define TINY 1.0e-30

#include "rkqs.h"

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
			double aa, double bb, double cc, double dd, double ee,
			void (*derivs)(double, double [], double [], double, double,
						   double, double, double),
            void (*rkqs)(double [], double [], int, double *, double, double,
				         double [], double *, double *, double, double, double,
						 double, double,
						 void (*)(double, double [], double [], double, double,
							      double, double, double))) {

	int     nstp, i;
	int     kmax, kount;
	double  xsav,x,hnext,hdid,h;
	double *yscal,*y,*dydx;
	double *xp, **yp, dxsav;

	/* initialising this within the func, which means this isn't generic */
    kmax = 100;

    xp = dvector(1, kmax);
    yp = dmatrix(1,nvar,1,kmax);
    dxsav = (x2 - x1) / 20.0;
	yscal=dvector(1,nvar);
	y=dvector(1,nvar);
	dydx=dvector(1,nvar);
	x=x1;
	h=SIGN(h1,x2-x1);
	*nok = (*nbad) = kount = 0;

	for (i=1;i<=nvar;i++) y[i]=ystart[i];
	if (kmax > 0) xsav=x-dxsav*2.0;
	for (nstp=1;nstp<=MAXSTP;nstp++) {
		(*derivs)(x,y,dydx, aa, bb, cc, dd, ee);
		for (i=1;i<=nvar;i++)
			yscal[i]=fabs(y[i])+fabs(dydx[i]*h)+TINY;
		if (kmax > 0 && kount < kmax-1 && fabs(x-xsav) > fabs(dxsav)) {
			xp[++kount]=x;
			for (i=1;i<=nvar;i++) yp[i][kount]=y[i];
			xsav=x;
		}

		if ((x+h-x2)*(x+h-x1) > 0.0) h=x2-x;
		(*rkqs)(y,dydx,nvar,&x,h,eps,yscal,&hdid,&hnext, aa, bb, cc, dd, ee,
				derivs);
		if (hdid == h) ++(*nok); else ++(*nbad);
		if ((x-x2)*(x2-x1) >= 0.0) {

			for (i=1;i<=nvar;i++) ystart[i]=y[i];

			if (kmax) {
				xp[++kount]=x;
				for (i=1;i<=nvar;i++) yp[i][kount]=y[i];
			}

			free_dvector(dydx,1,nvar);
			free_dvector(y,1,nvar);
			free_dvector(yscal,1,nvar);
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
