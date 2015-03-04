/* ============================================================================
* Ross's Optimal rooting depth model.
*
*   Optimisation hypothesis = maximise total N uptake (Utot) by determining
*   optimal distribution of root mass (R(z)) with soil depth and optimal
*   root depth (D, m)
*
*    References:
*    ----------
*    * McMurtire, R. et al (2012) Increased nitrogen-uptake efficiency at
*      elevated  CO2 explained by an hypothesis of optimal root function.
*      Ecology and Evolution, 2, 1235--1250
*
* AUTHOR:
*   Martin De Kauwe
*
* DATE:
*   17.02.2015
*
* =========================================================================== */
#include "mate.h"

#include "optimal_root_model.h"

void calc_opt_root_depth(double d0, double r0, double top_soil_depth,
                         double rtoti, double nsupply, double depth_guess,
                         double *root_depth, double *nuptake, double *rabove) {
    /*

        Parameters:
        -----------
        d0 : float
            Length scale for exponential decline of Umax(z)
        r0 : float
            root C at half-maximum N uptake (kg C/m3)
        top_soil_depth : float
            depth of soil assumed by G'DAY, note Ross comment about 30 cm (units=m)
            [email]
        rtoti : float
            Initial fine root C mass -> from G'DAY
        nsupply : float
            daily net N mineralisation in top soil layer from G'DAY
        depth_guess : float
            Initial guess at the rooting depth, used as the first point in the
            root depth optimisation scheme [m].

        Returns:
        --------
        root_depth : float
            rooting depth [m]
        nuptake : float
            N uptake from roots [gN m-2 yr-1]
        rabove : float

    */
    double depth;
    /* Determine maximum rooting depth for model for a value root C */
    depth = estimate_max_root_depth(rtoti, depth_guess, r0, d0);
    *root_depth = depth;
    /* Optimised plant N uptake */
    *nuptake = calc_plant_nuptake(depth, nsupply, d0, top_soil_depth);

    /* G'DAY requires root litter input to the top 30 cm of soil, so
       return the roots above this depth */
    printf("%f\n", depth);
    *rabove = calculate_root_mass_above_depth(rtoti, depth, r0, d0,
                                              top_soil_depth);

    return;

}


double estimate_max_root_depth(double rtoti, double depth_guess, double r0,
                               double d0) {
    /* Determing the maximum rooting depth through solving Eqn. B6. for
    rooting depth

    Parameters:
    -----------
    rtoti : float
        Initial fine root root C mass [from G'DAY]
    depth_guess : float
        initial starting guess at the root depth [m]

    Returns:
    --------
    rooting_depth : float
        optimised rooting depth [m]

    */
    double (*fPtr)(double, double, double, double) = &rtot_wrapper;
    double (*fprimePtr)(double, double, double, double) = &rtot_derivative;
    double root_depth;

    root_depth = newton(fPtr, fprimePtr, depth_guess, rtoti, r0, d0);

    return (root_depth);
}


double rtot_wrapper(double dmax, double rtoti, double r0, double d0) {
    /* Wrapper method that calls rtot. Need to subtract rtoti because we
    are optimising the depth (Dmax) but the rtot estimate has to match the
    rtot from GDAY, i.e. diff = 0.0

    Parameters:
    -----------
    rtoti : float
        Initial fine root root C mass [from G'DAY]
    rtot : function
        call rtot function to estimate a value of rtot given a rooting
        depth iteration

    Returns:
    --------
    val  : float
        A optimised rooting depth iteration
    */

    return (rtot(dmax, r0, d0) - rtoti);
}

double rtot(double dmax, double r0, double d0) {
    /* Estimate the total root biomass per unit ground area, i.e. the
    integral of root mass per unit soil volume, R(z), over depth z from the
    soil surface to the maximim rooting depth, dmax. (Eqn 8, in McM 2012)

    Parameters:
    -----------
    dmax : float
        Rooting depth [m]
    r0 : float
        Root C at half-max N uptake.
    d0 : float
        Length scale for exponential decline of Umax(z)

    Returns:
    --------
    rtot : float
        Total root C mass given a rooting depth

    */
    return (r0 * (2.0 * d0 * (exp(dmax / (2.0 * d0)) - 1.0) - dmax));

}

double rtot_derivative(double dmax, double rtoti, double r0, double d0) {
    /* Derivative of maximum root depth equation, rtot

    Parameters:
    -----------
    dmax : float
        Rooting depth [m]
    rtoti : float
        Initial fine root root C mass [from G'DAY]
    r0 : float
        Root C at half-max N uptake.
    d0 : float
        Length scale for exponential decline of Umax(z)

    Returns:
    --------
    val : float
        derivative of rtot

    */
    return (r0 * (exp(0.5 * dmax / d0) - 1.0) );

}

double calculate_root_mass_above_depth(double rtoti, double root_depth,
                                       double r0, double d0,
                                       double top_soil_depth) {
    /* Estimate cumulative root mass above depth, 30 cm for the G'DAY model

    Parameters
    ----------
    rtoti : float
        Initial fine root root C mass [from G'DAY]
    root_depth : float
        model rooting depth (m)

    Returns
    -------
    val : float
        cumulative root C mass above soil depth assumed by G'DAY model, 30cm
    */
    double arg1, arg2;

    arg1 = rtoti + 2.0 * r0 * d0 + root_depth * r0;
    arg2 = 1.0 - exp(-top_soil_depth / (2.0 * d0));
    return (arg1 * arg2 - r0 * top_soil_depth);
}

double calc_plant_nuptake(double root_depth, double nsupply, double d0,
                          double top_soil_depth) {
    /* Plant N uptake (Utot) as a func of maximum rooting depth

    This is the alternative eqn from McM word document

    Parameters
    ----------
    root_depth : float
        max rooting depth [m]
    z : float
        incremental depth provided by integration func
    nsupply : float
        soil N supply rate to plants per day [N/m2]
    top_soil_depth : float
        Depth of soil assumed by G'DAY model [m]

    Returns
    -------
    nuptake : float
        plant N uptake
    */
    double Umax, arg;

    Umax = calc_umax(nsupply, top_soil_depth, d0);
    arg = 1.0 - exp(-root_depth / (2.0 * d0));

    return (Umax * (arg * arg));
}


double calc_umax(double nsupply, double top_soil_depth, double d0) {
    /* Calculate potential N uptake integrated over all soil depths

    Parameters
    ----------
    nsupply : float
        N supply rate to a specified soil depth (probably 30 cm)

    Returns
    -------
    Umax : float
        potential N uptake integrated over all soil depths
    */
    return (nsupply / (1.0 - exp(-top_soil_depth / d0)));
}

double calc_net_n_uptake(double nuptake, double Nr, double rabove,
                         double root_lifespan, double rootn) {
    /* Calculate total potential N uptake

    Parameters
    ----------
    nuptake : float
        plant N uptake [g N m-2 day-]
    Nr : float
        nitrogen concentration of the fine roots (N:C ratio) [kg Dm-1]
    rabove : float
        root mass above depth (top soil in gday)
    root_lifespan : float
        root lifespan [year]

    Returns
    -------
    Umax : float
        potential N uptake integrated over all soil depths

    */
    return ( nuptake - (rootn * rabove / root_lifespan) );
}


double newton(double (*func)(double, double, double, double),
              double (*fprime)(double, double, double, double), double x0,
              double arg1, double arg2, double arg3) {
    /* Newton-Raphson: finds a zero of the func, given an inital guess

    Parameters
    ----------
    f : function
        The function whose zero is wanted.
    x0 : float
        An initial guess.
    fprime : function
        The derivative of the function
    args : tuple, optional
        Extra arguments to be used in the function call.
    tol : float, optional
        The allowable error of the zero value.
    maxiter : int, optional
        Maximum number of iterations.

    Returns
    -------
    val : float
        Estimated location where function is zero

    */
    int    iter, maxiter = 250;
    double tol = 1E-6;
    double x;

    for (iter = 0; iter < maxiter; iter++) {

        x = (x0 - (func(x0, arg1, arg2, arg3) / fprime(x0, arg1, arg2, arg3)));
        if (fabs(x - x0) < tol) {
            return x;
        }
        x0 = x;
    }
    fprintf(stderr, "Minimum not found!!\n");
    exit(EXIT_FAILURE);
}
