#include "hotokezaka.hpp"



/*double rate_to_rate_density(double R, double r_Sun, double rd, double zd)
{
    double rho_star, M_star, R_rate_density;
    rho_star = std::exp(-r_Sun/rd) / (2.0*zd); 
        //It also include a sigma_d,0 multiplier, but it will falls out
    M_star = 2.0 * M_PI * rd * rd;
        //It also include a sigma_d,0 multiplier, but it will falls out
    R_rate_density = rho_star * R / M_star;
    return R_rate_density;
}

double rate_density_to_rate(double R_density, double r_Sun, double rd, double zd)
{
    double rho_star, M_star, R_rate;
    rho_star = std::exp(-r_Sun/rd) / (2.0*zd); 
        //It also include a sigma_d,0 multiplier, but it will falls out
    M_star = 2.0 * M_PI * rd * rd;
        //It also include a sigma_d,0 multiplier, but it will falls out
    R_rate = R_density * M_star / rho_star;
}*/

Calculated_Numbers_Based_on_read_in_parameters::Calculated_Numbers_Based_on_read_in_parameters(Read_In_Parameters& params) : param(params)
{
    //These also include a sigma_d,0 multiplier, but it will falls out
    M_star = 2.0 * M_PI * param.rd * param.rd;
    rho_star = std::exp(-param.r_Sun/param.rd);
}
double Calculated_Numbers_Based_on_read_in_parameters::rate_density_to_rate(double R_density)
{
    double R_rate = R_density * M_star / rho_star;
    return R_rate;
}

double Calculated_Numbers_Based_on_read_in_parameters::rate_to_rate_density(double R)
{
    double R_rate_density = rho_star * R / M_star;
    return R_rate_density;
}

void Calculated_Numbers_Based_on_read_in_parameters::calculate_number_of_events()
{
    //integral class inicializálása
    //qtrap meghívása
    //eredmény egészre kerekítve beírása number_of_events változóba
}

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
    myout.write((char*)(&a_under_integral_limit),sizeof(a_under_integral_limit));
    myout.write((char*)(&last_integral_value),sizeof(last_integral_value));
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
        myout.write((char*)(&a_upper_integral_limit),sizeof(a_upper_integral_limit));
        myout.write((char*)(&last_integral_value),sizeof(last_integral_value));
        a_under_integral_limit += a_step;
    }
    myout.close();
}

void Rate_density::set_r0_rate_density_with_rate_density(double r0)
{
    r0_rate_density = r0;
    calc_normalize_factor();
}

WandermanRate::WandermanRate(double r0)
{
    r0_rate_density = r0;
    normalize_factor = r0_rate_density / std::exp(-0.9/0.39);
}

void WandermanRate::calc_normalize_factor()
{
    normalize_factor = r0_rate_density / std::exp(-0.9/0.39);
}

double WandermanRate::get_rate_density_at_z(double z)
{
    double rate_density;
    if(z <= 0.9)
    {
        rate_density = normalize_factor * std::exp((z-0.9)/0.39);
    }
    else
    {
        rate_density = normalize_factor * std::exp(-(z-0.9)/0.26);
    }
    return rate_density;
}

double WandermanRate::calc_Fx(double x)
{
   return get_rate_density_at_z(x);
}