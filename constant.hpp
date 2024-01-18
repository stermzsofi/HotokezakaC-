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
    const double cm_in_kpc = 3.240779289469756e-22;                 //cm in kpc
    const double cm_m3_in_kpc_m3 = 2.937998946027291e+64;           //cm^-3 in kpc^-3
    const double Today_Pu_abundance = 6e-17 * cm_m3_in_kpc_m3;      //today pu abundance from measuremets in kpc^-3
};
#endif