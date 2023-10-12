#ifndef HOTOKEZAKA
#define HOTOKEZAKA

#include <iostream>
#include <fstream>
#include <string>
#include <memory>
#include "Trapezoidal_rule/trapezoidal.hpp"
#include "Linear_interpol/linear_interpol.hpp"
#include "constant.hpp"

//I need include it to some class later
//double rate_to_rate_density(double R, double r_Sun, double rd, double zd);
//double rate_density_to_rate(double R_density, double r_Sun, double rd, double zd);


class Rate_density
{
    public:
        //virtual double get_rate_density_at_t(double t) = 0;
        virtual double get_rate_density_at_z(double z) = 0;
        void set_r0_rate_density_with_rate_density(double r0);
        virtual void calc_normalize_factor() = 0;
    protected:
        double r0_rate_density;
        double normalize_factor;
};

class WandermanRate : public Rate_density, public Fx
{
    public:
        WandermanRate(double r0);
        //double get_rate_density_at_t(double t);
        double get_rate_density_at_z(double z);
        double calc_Fx(double x);
        void calc_normalize_factor();
        //void set_r0_rate_density_with_rate_density(double r0);
        //void set_r0_rate_density_with_rate(double r0);
    private:
        //double r0_rate_density;
        //double normalize_factor;    //the rate function at z=0 must return r0 
        //valószínűleg őt is át lesz érdemes vinni Rate_density functionba + ott kell egy virtual 
        // calc_normalize_factor fv is
};

struct Read_In_Parameters
{
    double h_scale;
    double r_Sun;
    double width;
    double rd;
    double simulation_Time;
    double r0_rate;
    double sampleDt;
    //I think we also will need the scaleheight of the disk
    double zd;
    int number_of_runs;
    std::string read_in_rate_function;
    std::unique_ptr<Rate_density> rate_function; //unique_ptr kell valószínűleg

    void read_parameter_file(); //Read the values based on parameter file and run init function
    void init();    //ez valószínű nem is fog kelleni, már read közben létre lehet hozni rate_function-t.
        //- create rate_function based in read_in_rate_function string
        // - run the set_r0_rate_density_with_rate_density function;

    //double rate_to_rate_density(double R);
    //double rate_density_to_rate(double R_density);
};

class Calculated_Numbers_Based_on_read_in_parameters
{
    public:
        Calculated_Numbers_Based_on_read_in_parameters(Read_In_Parameters& params);
        double rate_to_rate_density(double R);
        double rate_density_to_rate(double R_density);
        void calculate_number_of_events();
        double time_to_z(double t);
        double z_to_time(double z);
        void init_cumulative_distribution();
        double interpolate_random_number_to_time(double myrand);
        void init();
    private:
        Read_In_Parameters& param;
        Quad_Trapezoidal integral;
        Interpolator time_to_z_interpolator;
        Interpolator cumulative_distribution_interpolator;
        double M_star;
        double rho_star;
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
};


class Create_events_and_calc_number_density
{

};

#endif