#include "trapezoidal.hpp"

/*
    **********

        Test classes: for test the operation of Trapezoidal_rule

    **********
*/


/**Functions for Line class*/
    /*Constructor*/
    Line::Line(double slope, double yintercept)
    {
        a = slope;
        b = yintercept;
    }

    /*calc_Fx function*/
    double Line::calc_Fx(double x)
    {
        double res;
        res = a*x + b;
        return res;
    }

/**Functions for Squarerootmy class**/
    /*Constructor*/
    Squarerootmy::Squarerootmy(double eltol)
    {
        a = eltol;
    }

    /*calc_Fx function*/
    double Squarerootmy::calc_Fx(double x)
    {
        double res;
        res = sqrt(x)+a;
        return res;
    }



    /*Power_function class functions*/
    Power_function::Power_function(double a, double b, double c)
    {
        A = a;
        B = b;
        C = c;
    }

    double Power_function::calc_Fx(double x)
    {
        double res;
        res = A * pow(x,B) + C;
        return res;
    }

    Quad_Trapezoidal::Quad_Trapezoidal(Fx& _fv, double _x_start, double _x_end, double _eps) : fv(&_fv), x_start(_x_start), x_end(_x_end), eps(_eps)
    {
        n = 0;
    }

    double Quad_Trapezoidal::next()
    {
        double x, tnm, del;
        long double sum;
        int it, j;
        n++;
        if ( 1 == n )
        {
            s = 0.5 * (x_end - x_start) * (fv->calc_Fx(x_start) + fv->calc_Fx(x_end));
            return s;
        }
        else
        {
            for (it = 1, j = 1; j < n - 1; j++)
            {
                it <<=1;    //bit eltolása balra: ha jól gondolom, ezzel adja vissza, hogy hány részre van felolsztva a tartomány
            }
            tnm = it;
            del = (x_end - x_start)/tnm;
            x = x_start + 0.5 * del;
            for(sum = 0.0, j = 0; j < it; j++, x += del)
            {
                sum += fv->calc_Fx(x);
            }
            s = 0.5 * (s + (x_end - x_start) * sum/tnm);
            return s;
        }
    }

    double Quad_Trapezoidal::qtrap()
    {
        const int JMAX = 20;        //max number of steps: 2^20
        double olds = 0.0;
        for (int j = 0; j < JMAX; j++)
        {
            next();
            if (!std::isfinite(s))
            {
                n = 0;
                //std::cout << j << "\t" << olds << std::endl;
                throw("WARNING: S is nan!");
            }  
            if(j > 3)
            {
                if (std::abs(s - olds) < eps * std::abs(olds))
                {
                    n = 0;
                    //std::cout << "The needed number of j is: " << j << " corresponds to " << std::pow(2,j) << " step. The error is: " << std::abs(s)/std::abs(olds) << std::endl;
                    return s;
                }
                //std::cout << "The error is: " << std::abs(s)/std::abs(olds) << std::endl;
            }
            //if(j < JMAX - 1) 
            olds = s;
        }
        n = 0;
        throw std::runtime_error("ERROR! Too many steps in routine qtrap.");
        //std::cout << "Too many steps in routine qtrap. The error is: " << std::abs(s)/std::abs(olds) << std::endl;
        //return s;
    }

    /*
        **********

            Trapezoidal_rule class functions: 
            It will calculate the integral with trapezoidal rule

        **********
    */

    /*Constructors*/
    Trapezoidal_rule::Trapezoidal_rule(Fx& fx, double xstart, double xend, double d_x, int n)
    {
        fv = &fx;
        x_start = xstart;
        x_end = xend;
        dx = d_x;
        partions = n;
    }

    Trapezoidal_rule::Trapezoidal_rule(Fx& fx, double xstart, double xend, int n)
    {
        fv = &fx;
        x_start = xstart;
        x_end = xend;
        partions = n;
    }

    void Trapezoidal_rule::set_ranges(double xstart, double xend)
    {
        x_start = xstart;
        x_end = xend;
    }

    void Trapezoidal_rule::set_xstart(double xstart)
    {
        x_start = xstart;
    }

    void Trapezoidal_rule::set_xend(double xend)
    {
        x_end = xend;
    }

    void Trapezoidal_rule::set_partions(int n)
    {
        partions = n;
    }

    void Trapezoidal_rule::set_newfx(Fx& fx)
    {
        fv = &fx;
    }

    /**Integrate part*/
    
    //This function will divide the intervall "partions" number of equivalent parts and return with the integral value
    double Trapezoidal_rule::calc_integral()
    {
        double integral = 0;
        double h = (x_end - x_start) / partions;
        double x0, x1, int_akt;
        for(int i = 1; i<= partions; i++)
        {
            x0 = x_start + (i-1)*h;
            x1 = x_start + i*h;
            int_akt = h*(fv ->calc_Fx(x1)+fv ->calc_Fx(x0))/2.0;
            integral += int_akt;
        }
        return integral;
    }

    //This function will goes through the intervall with dx steps and return with the vaule of the integral
    double Trapezoidal_rule::calc_integral_dx()
    {
        double integral = 0;
        double x0, x1, int_akt,x,h;
        for(x = x_start; x <= x_end; x += dx)
        {
            x0 = x;
            x1 = x + dx;
            h=dx;
            if(x1 > x_end)
            {
                x1 = x_end;
                h=x_end-x0;
            }
            int_akt = h*(fv ->calc_Fx(x1)+fv ->calc_Fx(x0))/2.0;
            integral += int_akt;
        }
        /*if(x < x_end + dx)
        {
            x0 = x - dx;
            x1 = x_end;
            double h = x1-x0;
            int_akt = h*(fv ->calc_Fx(x1)+fv ->calc_Fx(x0))/2.0;
            integral += int_akt;
        }*/
        return integral;
    }

    Newthon_Raphson::Newthon_Raphson(Fx& _fx, dFx& _dfx)
    {
        fx = &_fx;
        dfx = &_dfx;
        x0 = 1.0;
        eps = 1e-5;
        nmax = 1000;
    }

    double Newthon_Raphson::Newthon_Raphson_method()
    {
        double xn = x0;
        double xnp1=xn;
        for(int n = 0; n < nmax; n++)
        {
            xn = xnp1;
            xnp1 = xn - fx->calc_Fx(xn)/dfx->calc_dFx(xn);
            if(std::abs(xn-xnp1) < eps) return xnp1;
        }
        throw std::runtime_error("ERROR! Couldn't calculate Rc, too many steps needed!");
    }

    double Newthon_Raphson::Newthon_Raphson_safe()
    {
        double xh, xl;  //for lower and higher boundary
        double fl = fx->calc_Fx(x1);    
        double fh = fx->calc_Fx(x2);
        if((fl > 0 && fh > 0) || (fl < 0 && fh < 0))    //if in this boundary fx couldn't reach 0
        {
            std::cout << "xl = " << x1  << "\tfl = " << fl << "\txh = " << x2 << "\tfh = " << fh << std::endl;
            throw std::runtime_error("ERROR! Safe rutine boundary not appropriate!");
        }
        if(fl == 0) return x1; 
        if(fh == 0) return x2;
        if(fl < 0.0)
        {
            xl = x1;
            xh = x2;
        }
        else
        {
            xl = x2;
            xh = x1;
        }
        double xn = 0.5*(x1 + x2);  //find the center of the boundary
        double dxold = std::abs(x2-x1);     //last before and last step
        double dx = dxold;
        double f = fx->calc_Fx(xn);
        double df = dfx->calc_dFx(xn);
        for(int n = 0; n < nmax; n++)
        {
            if( (((xn - xh)*df - f) * ((xn - xl) * df - f)) > 0 || (std::abs(2.0*f) > std::abs(dxold * df)) )
            {
                dxold = dx;
                dx = 0.5 * (xh - xl);
                xn = xl + dx;
                if(xl == xn) return xn;
            } 
            else
            {
                dxold = dx;
                dx = f/df;
                double xnmin1 = xn;
                xn -= dx;
                if(xnmin1 == xn) return xn;
            }
            if(std::abs(dx) < eps) return xn;
            f = fx->calc_Fx(xn);
            df = dfx->calc_dFx(xn);
            if(f < 0.0)
            {
                xl = xn;
            }
            else
            {
                xh = xn;
            }
        }
        throw std::runtime_error("ERROR! Too many iteration step need in safe Newton-Raphson method!");
    }

    //Bad functions
    /*double Fx_for_Rc::calc_Fx(double x)
    {
        double f_Rc;
        f_Rc = B * std::pow(x,-5./6.) * 0.5 - (B * 0.5 * std::pow(rw,-2)* std::pow(x,7./6.)) - 1.0;
        return f_Rc;
    }

    double dFx_for_Rc::calc_dFx(double x)
    {
        double df_Rc;
        df_Rc = B * 5./12. * std::pow(x,-11./6.) - B * std::pow(rw, -2) * 7./12. * std::pow(x,1./6.);
        return df_Rc;
    }*/
