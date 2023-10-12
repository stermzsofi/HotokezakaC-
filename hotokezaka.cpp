#include "hotokezaka.hpp"


Calculated_Numbers_Based_on_read_in_parameters::Calculated_Numbers_Based_on_read_in_parameters(Read_In_Parameters& params) : param(params)
{
    //These also include a sigma_d,0 multiplier, but it will falls out
    M_star = 2.0 * M_PI * param.rd * param.rd;
    rho_star = std::exp(-param.r_Sun/param.rd);
    init();
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


void Calculated_Numbers_Based_on_read_in_parameters::init()
{
    //először: létezik-e fájl idő->z interpolációhoz. ha nem, csináljunk egyet, aztán olvassuk be.
    if(!std::filesystem::exists("interpolate_table.dat"))
    {
        create_file_for_interpolation file;
        file = create_file_for_interpolation();
        file.save_data_to_file();
    }
    time_to_z_interpolator.init("interpolate_table.dat");
    //+1: store the time of the universe
    Time_of_the_Universe = time_to_z_interpolator.genericInterpolationX(1.0);
    //second: calculate the number of events
        //simulation_Time -> z_simulation
        double z_simulation = time_to_z(param.simulation_Time);
        //integral initialize (set rate, ranges)
        //integral = Quad_Trapezoidal(*param.rate_function, 0.0, z_simulation);
        integral.set_newfx(*param.rate_function);
        //run qtrap and calculate the number_of_events_d
        //calculate the number_of_events
    //third: calculate cumulative_distribution function
        //first: create two vectors
        //calculate the integral of rate/number_of_events_d and store on the two vectors
        //initialize the cumulative_distribution_interpolator
    
}

double Calculated_Numbers_Based_on_read_in_parameters::time_to_z(double t)
{
    //t: this is the lookback time
    //but I need: t_real. First: need to calculate t_real
    double t_real = Time_of_the_Universe - t;
    //this will give back the a value
    double a = time_to_z_interpolator.genericInterpolationY(t_real);
    //z can be calculate like this:
    double z = 1.0/a - 1.0;
    return z;
}

double Calculated_Numbers_Based_on_read_in_parameters::z_to_time(double z)
{
    //I will need 'a' value for the interpolation
    double a = 1.0/(z + 1.0);
    //The interpolation will give back the real time
    double t_real = time_to_z_interpolator.genericInterpolationX(a);
    //the lookback time based on it:
    double t = Time_of_the_Universe - t_real;
    return t;
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

/*double WandermanRate::calc_Fx(double x)
{
   return get_rate_density_at_z(x);
}*/