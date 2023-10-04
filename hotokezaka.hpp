#ifndef HOTOKEZAKA
#define HOTOKEZAKA

#include <iostream>
#include <fstream>
#include "Trapezoidal_rule/trapezoidal.hpp"
#include "constant.hpp"

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
        double a_step = 1e-4;
};

class Rate
{
    public:
        virtual double get_rate(double t) = 0;
};

class WandermanRate : public Rate, public Fx
{
    public:
        double get_rate(double t);
        double calc_Fx(double x);
};

#endif