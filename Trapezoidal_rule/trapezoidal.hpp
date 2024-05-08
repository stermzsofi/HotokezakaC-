#ifndef TRAPEZOIDAL
#define TRAPEZOIDAL

#include <iostream>
#include <cmath>
#include <memory>


/*Fx class: this is a mother class. Children classes must be a calc_Fx function, which calculate the value of the function at x*/
class Fx
{
    public:
        virtual double calc_Fx(double x) = 0;
};

class dFx
{
	public:
		virtual double calc_dFx(double x) = 0;
};

	/*Some classes to the test*/
	/*Line class. It derive from fx. It definie a line, f(x)=a*x+b*/
	class Line : public Fx
	{
		public:
			Line(double, double);
			double calc_Fx(double x);
		protected:
			double a;
			double b;
	};


	/*Squareroot class. It derive from fx. It definie sqrt(x)+a*/
	class Squarerootmy : public Fx
	{
		public:
			Squarerootmy(double);
			double calc_Fx(double x);
		private:
			double a;
	};

	/*Power_function class. f(x)=A*x^B+C*/

	class Power_function : public Fx
	{
		public:
			Power_function(double, double, double);
			double calc_Fx(double x);
		private:
			double A;
			double B;
			double C;
	};

class Quadrature
{
	public:
		int n;
		virtual double next() = 0;
};

class Quad_Trapezoidal : public Quadrature
{
	public:
		Quad_Trapezoidal(){};
		Quad_Trapezoidal(Fx& _fv, double _x_start, double _x_end, double _eps);
		void set_eps(double _eps) {eps = _eps;}
		void set_xstart(double _x_start){x_start = _x_start;}
		void set_xend(double _x_end){x_end = _x_end;}
		void set_n(double _n){n = _n;}
		void set_newfx(Fx& fx){fv = &fx;}
		double next();
		double qtrap();
	private:
		Fx* fv;
		double x_start;
		double x_end;
		double s;
		double eps;
};


/*Trapezoidal rule class*/
	/* It need an Fx derived class
	 * and 2 double (the start of the integrate and the end of it
	 * An int variable, which set how many part do you want to divide the integrate (partions)
	 * We can set the start and the end of the integrate, so we don't need many classes for different ranges
	 * calc_integral will return the integral for the current setups
	 * */
class Trapezoidal_rule
{
    public:
        Trapezoidal_rule(Fx& fx, double xstart, double xend, double d_x, int n);
		Trapezoidal_rule(Fx& fx, double xstart, double xend, int n);
		Trapezoidal_rule(){}
        void set_ranges(double sxtart, double xend);
        void set_xstart(double xstart);
        void set_xend(double xend);
        void set_partions(int n);
		void set_newfx(Fx& fx);
        double calc_integral();
		double calc_integral_dx();
    private:
        Fx* fv;
        double x_start;
        double x_end;
		double dx;
        int partions;
};

class Newthon_Raphson
{
	public:
		Newthon_Raphson(Fx& _fx, dFx& _dfx);
		void set_x0(double _x0) {x0 = _x0;}
		void set_x1(double _x1) {x1 = _x1;}
		void set_x2(double _x2) {x2 = _x2;}
		void set_eps(double _eps) {eps = _eps;}
		void set_nmax(int _nmax) {nmax = _nmax;}
		double Newthon_Raphson_method();
		double Newthon_Raphson_safe();
	private:
		Fx* fx;
		dFx* dfx;
		double x0;
		double eps;
		double x1, x2;
		int nmax;
};

#endif
