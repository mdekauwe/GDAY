/* ============================================================================
* The Cash-Karp Runge-Kutta Step
*
* - From numerical recipies in C, see Press et al. 1992.
*
* NOTES:
* I've refactored this function so that all the allocation/free calls have
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

#define NRANSI
#include "nrutil.h"
#include "rkck.h"


void rkck(double y[], double dydx[], int n, double x, double h, double yout[],
	      double yerr[], double aa, double bb, double cc, double dd, double ee,
		  nrutil *nr,
	      void (*derivs)(double, double [], double [], double, double, double,
		  				 double, double))
{
	int           i;
	static double a2=0.2,a3=0.3,a4=0.6,a5=1.0,a6=0.875,b21=0.2,
				  b31=3.0/40.0,b32=9.0/40.0,b41=0.3,b42 = -0.9,b43=1.2,
				  b51 = -11.0/54.0, b52=2.5,b53 = -70.0/27.0,b54=35.0/27.0,
				  b61=1631.0/55296.0,b62=175.0/512.0,b63=575.0/13824.0,
				  b64=44275.0/110592.0,b65=253.0/4096.0,c1=37.0/378.0,
				  c3=250.0/621.0,c4=125.0/594.0,c6=512.0/1771.0,
				  dc5 = -277.0/14336.0;
	double        dc1=c1-2825.0/27648.0,dc3=c3-18575.0/48384.0,
		          dc4=c4-13525.0/55296.0,dc6=c6-0.25;
	//double       *ak2, *ak3, *ak4, *ak5, *ak6, *ytemp;


	//ak2=dvector(1,n);
	//ak3=dvector(1,n);
	//ak4=dvector(1,n);
	//ak5=dvector(1,n);
	//ak6=dvector(1,n);
	//ytemp=dvector(1,n);

	for (i=1;i<=n;i++)
		nr->ytemp[i]=y[i]+b21*h*dydx[i];
	(*derivs)(x+a2*h,nr->ytemp,nr->ak2, aa, bb, cc, dd, ee);
	for (i=1;i<=n;i++)
		nr->ytemp[i]=y[i]+h*(b31*dydx[i]+b32*nr->ak2[i]);
	(*derivs)(x+a3*h,nr->ytemp,nr->ak3, aa, bb, cc, dd, ee);
	for (i=1;i<=n;i++)
		nr->ytemp[i]=y[i]+h*(b41*dydx[i]+b42*nr->ak2[i]+b43*nr->ak3[i]);
	(*derivs)(x+a4*h,nr->ytemp,nr->ak4, aa, bb, cc, dd, ee);
	for (i=1;i<=n;i++)
		nr->ytemp[i]=y[i]+h*(b51*dydx[i]+b52*nr->ak2[i]+b53*nr->ak3[i]+b54*nr->ak4[i]);
	(*derivs)(x+a5*h,nr->ytemp,nr->ak5, aa, bb, cc, dd, ee);
	for (i=1;i<=n;i++)
		nr->ytemp[i]=y[i]+h*(b61*dydx[i]+b62*nr->ak2[i]+b63*nr->ak3[i]+b64*nr->ak4[i]+b65*nr->ak5[i]);
	(*derivs)(x+a6*h,nr->ytemp,nr->ak6, aa, bb, cc, dd, ee);
	for (i=1;i<=n;i++)
		yout[i]=y[i]+h*(c1*dydx[i]+c3*nr->ak3[i]+c4*nr->ak4[i]+c6*nr->ak6[i]);
	for (i=1;i<=n;i++)
		yerr[i]=h*(dc1*dydx[i]+dc3*nr->ak3[i]+dc4*nr->ak4[i]+dc5*nr->ak5[i]+dc6*nr->ak6[i]);

	//free_dvector(ytemp, 1, n);
	//free_dvector(ak6, 1, n);
	//free_dvector(ak5, 1, n);
	//free_dvector(ak4, 1, n);
	//free_dvector(ak3, 1, n);
	//free_dvector(ak2, 1, n);
}
#undef NRANSI
