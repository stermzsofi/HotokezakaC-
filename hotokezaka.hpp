#ifndef HOTOKEZAKA
#define HOTOKEZAKA

#include <iostream>
#include <fstream>
#include <string>
#include <memory>
#include <filesystem>
#include <algorithm>
#include <random>
//#include <boost/math/distributions/laplace.hpp>
#include <boost/random/laplace_distribution.hpp>
//#include <boost/math>
#include "Trapezoidal_rule/trapezoidal.hpp"
#include "Linear_interpol/linear_interpol.hpp"
#include "constant.hpp"

class Time_redshift
{
    public:
        Time_redshift();
        double time_to_z(double t);
        double z_to_time(double z);
        double time_of_universe();
        //Interpolator time_to_z_interpolator;
    private:
        Interpolator time_to_z_interpolator;
        double Time_of_the_Universe;
};

/*************************************************************/
/********DIFFERENT RATE DENSITY FUNCTIONS ********************/
/*************************************************************/
class Rate_density : public Fx
{
    public:
        //virtual double get_rate_density_at_t(double t) = 0;
        virtual double get_rate_density_at_z(double z) = 0;
        virtual double get_rate_density_at_t(double t) = 0;
        void set_r0_rate_density_with_rate_density(double r0);
        virtual void calc_normalize_factor() = 0;
        double calc_Fx(double x){ return get_rate_density_at_t(x);}
    protected:
        double r0_rate_density;
        double normalize_factor;
};

class WandermanRate : public Rate_density
{
    public:
        WandermanRate(Time_redshift& t_z) : time_z(t_z){};
        //WandermanRate(double r0);
        double get_rate_density_at_z(double z);
        double get_rate_density_at_t(double t);
        //double calc_Fx(double x);
        void calc_normalize_factor();
        //void set_r0_rate_density_with_rate_density(double r0);
        //void set_r0_rate_density_with_rate(double r0);
    private:
        Time_redshift& time_z;
        //double r0_rate_density;
        //double normalize_factor;    //the rate function at z=0 must return r0 
        //valószínűleg őt is át lesz érdemes vinni Rate_density functionba + ott kell egy virtual 
        // calc_normalize_factor fv is
};

class HopkinsRate : public Rate_density
{
    public:
        HopkinsRate(Time_redshift& t_z);
        double get_rate_density_at_z(double z);
        double get_rate_density_at_t(double t);
        void calc_normalize_factor();
    private:
        Time_redshift& time_z;
        double a,b,c,d,h;
};

struct Read_In_Parameters
{
    //in inputEventGenerator.in at initial code
    double h_scale;
    double r_Sun;
    double width;
    double rd;
    double simulation_Time;
    double r0_rate;
    double sampleDt;
    //in inputHotokezakaGalaxy.in at initial code
    double alpha;
    double vt;
    //I think we also will need the scaleheight of the disk
    double zd;
    //in tauList.in at initial code
    double tau;     //later could be a vector for more than one tau values
    //in intputEventGeneratator.in at initial code
    int number_of_runs;
    std::string read_in_rate_function;
    std::unique_ptr<Rate_density> rate_function; //unique_ptr kell valószínűleg

    //for time redshift interpolations
    Time_redshift time_z;

    void read_parameter_file(); //Read the values based on parameter file and run init function
    void init();
        //- create rate_function based in read_in_rate_function string
        // - run the set_r0_rate_density_with_rate_density function -> 
        //       not here, in Calculated_Numbers_Based_on_read_in_parameters class
};



class Calculated_Numbers_Based_on_read_in_parameters
{
    public:
        Calculated_Numbers_Based_on_read_in_parameters(Read_In_Parameters& params/*, Time_redshift& t_z*/);
        double rate_to_rate_density(double R);
        double rate_density_to_rate(double R_density);
        void calculate_number_of_events();
        //double time_to_z(double t);
        //double z_to_time(double z);
        //void init_cumulative_distribution();
        double interpolate_random_number_to_time(double myrand);
        int get_number_of_events() {return number_of_events;}
        void init();
        Read_In_Parameters& param;
        double D;       //D=alpha*(vt/7)*(H/0.2)
        //in Myr
        double taumix;  //taumix = 300 (R0/10)^-2/5 * (alpha/0.1)^-3/5*(vt/7)^-3/5*(H/0.2)^-3/5
        double Ni;
    private:
        //Read_In_Parameters& param;
        Quad_Trapezoidal integral;
        //Quad_Trapezoidal integrator_for_cumulative_distribution;
        //Interpolator time_to_z_interpolator;
        //Time_redshift& time_z;
        Interpolator cumulative_distribution_interpolator;
        double M_star;
        double rho_star;
        double Time_of_the_Universe;
        double number_of_events_d;
        int number_of_events;
};

class Fx_for_timeintegral : public Fx
{
    public:
        double calc_Fx(double x);
};

class create_file_for_interpolation
{
    public:
        create_file_for_interpolation();
        void save_data_to_file();
    private:
        Fx_for_timeintegral fx_time;
        Quad_Trapezoidal integral;
        double a_step = 1e-4;   //ezt lehet érdemes lenne majd nagyobbra állítani memória miatt
                                //1e-3 also enough 
};

struct randomEvent
{
    double time;
    double x;
    double y;
    double z;
    double distance;
};

class Create_events_and_calc_number_density
{
    public:
        Create_events_and_calc_number_density(Calculated_Numbers_Based_on_read_in_parameters& calc_par);
        void allEvent_number_densities();
    private:
        //Read_In_Parameters& param;
        //std::random_device rd;
        //std::mt19937 mt;
        std::uniform_real_distribution<double> rand_number_0_1;
        boost::random::laplace_distribution<double> rand_laplace;
        Calculated_Numbers_Based_on_read_in_parameters& calculated_parameters;
        double delta_t_crit; //delta_t_crit = ((8*pi*H*D)^2)/((4*pi*D)^3)
                                //if delta_t >= delta_t_crit, than 8*pi*D*deltaj will be smaller
                                //other case: (4*pi*D*deltat)^3/2 will be smaller
        std::vector<double> sampling_time_points;
        std::vector<double> number_densites;
        std::vector<double> median_number_densities;
        randomEvent create_random_event();
        void calc_number_density_for_an_event();
        double calc_Kj(double delta_time);
        double const_for_Kj_1;
        double const_for_Kj_2;
        
};

#endif