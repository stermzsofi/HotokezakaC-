#ifndef CONSTANTS
#define CONSTANTS

namespace constants
{
    const double H0 = 67.4;                                         //km/s/Mpc
    const double H0_correct_Unit = H0 * 1.0227121650537078e-06;     //Mpc/Myr/Mpc
    const double one_per_H0 = 1.0/H0_correct_Unit;
    const double Omega_m = 0.315;                                   //matter density parameter
    const double Omega_lambda = 0.6847;                             //dark energy density parameter
    const double Omega_r = 3e-4;                                    //radiation density parameter
    const double My_Err = 1e-7;                                     //Error limit to check the equality of two number
};
#endif