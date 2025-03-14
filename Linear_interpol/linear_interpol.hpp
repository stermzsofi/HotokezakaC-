#ifndef LINEAR_INTERPOL
#define LINEAR_INTERPOL
#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <tuple>
#include <string>
#include <algorithm>

/*A struct, contain x and the corresponding Fx. Will use it at interpolation as a vector*/
struct Table_x_Fx
{
    double x;
    double Fx;
    Table_x_Fx(double _x, double _Fx) : x(_x), Fx(_Fx) {}
};


class Interpolator
{
    private:
        std::vector<Table_x_Fx> Tablenew;   //The datas
        double interpolate_linear(double x, double x1, double x0, double y1, double y0);    //It made the interpolation
    public:
        void init(std::string InFileName);  //uploaded Tablenew data vector from datafile, set kszimax
        void init(std::vector<double> x, std::vector<double> fx);   //uploaded Tablenew data vector from vectors
        double genericInterpolationX(double x);  //general version of interpolation
        double genericInterpolationY(double Fx);
};

bool isEqual(double x, double y, double eps);
#endif