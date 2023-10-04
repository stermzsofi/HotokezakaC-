#include "hotokezaka.hpp"

double Fx_for_timeintegral::calc_Fx(double x)
{
    double denominator = constants::Omega_m /x + constants::Omega_r/(x*x) + constants::Omega_lambda * x * x;
    double res = constants::one_per_H0 * 1.0/(std::sqrt(denominator));
    return res; 
}

create_file_for_interpolation::create_file_for_interpolation()
{
    integral = Quad_Trapezoidal(fx_time, 0,1,1e-3);
}

void create_file_for_interpolation::save_data_to_file()
{
    std::ofstream myout;
    myout.open("interpolate_table.dat",std::ios::out | std::ios::binary);
    double last_integral_value = 0;
    double a_under_integral_limit = 0.0;
    double a_upper_integral_limit = 0.0;
    std::cout << a_under_integral_limit << "\t" << last_integral_value << std::endl;
    while(std::abs(1.0-a_upper_integral_limit) > constants::My_Err)
    {
        a_upper_integral_limit += a_step;
        integral.set_xstart(a_under_integral_limit);
        if(a_upper_integral_limit > 1.0)
        {
            a_upper_integral_limit = 1.0;
            integral.set_xend(a_upper_integral_limit);
        }
        else
        {
            //a_upper_integral_limit = a_under_integral_limit + a_step;
            integral.set_xend(a_upper_integral_limit);
        }
        last_integral_value += integral.qtrap();
        std::cout << a_upper_integral_limit << "\t" << last_integral_value << std::endl;
        a_under_integral_limit += a_step;
    }
}