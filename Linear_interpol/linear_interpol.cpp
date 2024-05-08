#include "linear_interpol.hpp"

bool isEqual(double x, double y, double eps)
{
    return std::abs(x-y) < eps;
}

void Interpolator::init(std::string InFileName)
{
    std::ifstream Infile;
    Infile.open(InFileName,std::ios::in | std::ios::binary);
    if(!Infile.is_open()) throw std::runtime_error("ERROR! "+InFileName+" cannot be opened.");
    double x, Fx;
    Infile.seekg(0);
    Infile.read((char*)&x,sizeof(x));
    Infile.read((char*)&Fx,sizeof(Fx));
    Tablenew.push_back(Table_x_Fx{x,Fx});
    while(!Infile.eof())
    {
        Infile.read((char*)&x,sizeof(x));
        Infile.read((char*)&Fx,sizeof(Fx));
        Tablenew.push_back(Table_x_Fx{x,Fx});
    }
    Infile.close();
}

double Interpolator::interpolate_linear(double x, double x1, double x0, double y1, double y0)
{
    return y0 + (x - x0) * (y1 - y0 ) / (x1 - x0);
}


double Interpolator::genericInterpolationX(double x)
{
    auto x_2 = std::find_if(Tablenew.begin(), Tablenew.end(), [x](Table_x_Fx now){return x <= now.x;});
    if (x_2 == Tablenew.end()) return std::prev(x_2,1)->Fx;
    auto x_1 = std::prev(x_2,1);
    double F_x = interpolate_linear(x, x_2 ->x,x_1->x, x_2->Fx, x_1->Fx );
    return F_x;
}

double Interpolator::genericInterpolationY(double Fx)
{
    auto Fx_2 = std::find_if(Tablenew.begin(), Tablenew.end(), [Fx](Table_x_Fx now){return Fx <= now.Fx;});
    if (Fx_2 == Tablenew.end()) return std::prev(Fx_2, 1) -> Fx;
    auto Fx_1 = std::prev(Fx_2,1);
    double x = interpolate_linear(Fx, Fx_2->Fx, Fx_1-> Fx, Fx_2->x, Fx_1 -> x);
    return x;
}


void Interpolator::init(std::vector<double> x, std::vector<double> fx)
{
    for( size_t i =0; i< fx.size(); i++) //this is evil: no check for sizes
    {
        this->Tablenew.push_back(Table_x_Fx{x[i],fx[i]});
    }
}