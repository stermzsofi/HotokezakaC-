#include "hotokezaka.hpp"

bool IsDouble(const std::string& s)
{
    std::regex reg("^[+-]?[[:digit:]]+\\.?[[:digit:]]*[eE][+-]?[[:digit:]]+$|^[+-]?[[:digit:]]+\\.?[[:digit:]]*$");
    return std::regex_match(s,reg);
}

//initial values for all parameters
Read_In_Parameters::Read_In_Parameters()
{
    
}


void Read_In_Parameters::read_parameter_file()
{
    /*Open parameters.txt file*/
    std::ifstream input;
    input.open("parameters.txt");
    /*lines: it is a string vector. It will contain the lines of parameters.txt file*/
    std::vector<std::string> lines;
    /*Put the parameters.txt file lines to lines, string type vector*/
    if(input.is_open())
    {
        /*line: actual line of the file*/
        std::string line = "START";
        while(!input.eof())
        {
            std::getline(input, line);
            if(input.eof()) break;
            lines.push_back(line);
        }
    }
    else
    {
        throw std::runtime_error("ERROR! Couldn't open parameters.txt file\n");
    }
    /*One line contain two things: first is the name and second is the value*/
    std::string name;
    std::string value;
    for(auto actual_line : lines)
    {
        std::istringstream act_line_stream(actual_line);
        std::getline(act_line_stream, name, '\t');
        std::getline(act_line_stream, value, '\t');
        if("h_scale" == name)
        {
            h_scale = stod(value);
        }
        else if("r_Sun" == name)
        {
            r_Sun = stod(value);
        }
        else if("width" == name)
        {
            width = stod(value);
        }
        else if("rd" == name)
        {
            rd = stod(value);
        }
        else if("simulation_Time" == name)
        {
            simulation_Time = stod(value);
        }
        else if("r0_rate" == name)
        {
            r0_rate = stod(value);
        }
        else if("sampleDt" == name)
        {
            sampleDt = stod(value);
        }
        else if("alpha" == name)
        {
            alpha = stod(value);
        }
        else if("vt" == name)
        {
            vt = stod(value);
        }
        else if("tau" == name)
        {
            tau = stod(value);
        }
        else if("number_of_runs" == name)
        {
            number_of_runs = stoi(value);
        }
        else if("rate_function" == name)
        {
            read_in_rate_function = value;
        }
        else if("Ni_calc_method" == name)
        {
            ni_calculation_method = value;
        }
        else if("element_initial_prod_ratio" == name)
        {
            element_initial_production_ratio = stod(value);
        }
        else if("output" == name)
        {
            out_file = value;
        }
        else
        {
            throw std::runtime_error("ERROR! Error with parameter file line: "+name +" " +value);
        }
    }
    init();
}

void Read_In_Parameters::init()
{
    if(read_in_rate_function == "Wanderman")
    {
        rate_function.reset(new WandermanRate(time_z));
    }
    else if(read_in_rate_function == "Hopkins")
    {
        rate_function.reset(new HopkinsRate(time_z));
    }
    else
    {
        throw std::runtime_error("Unknown rate function "+read_in_rate_function);
    }
}

/************************************************
 * **********************************************
           Time-redshift exchange functions
*************************************************
************************************************/

Time_redshift::Time_redshift()
{
    //if interpolation file not exist, create it
    if(!std::filesystem::exists("interpolate_table.dat"))
    {
        create_file_for_interpolation file;
        file = create_file_for_interpolation();
        file.save_data_to_file();
    }
    //open interpolation file
    time_to_z_interpolator.init("interpolate_table.dat");
    //The time of the universe is the time at a=1
    Time_of_the_Universe = time_to_z_interpolator.genericInterpolationX(1.0);
}

//The time of the universe is the time at a=1
double Time_redshift::time_of_universe()
{
    return Time_of_the_Universe;
}

double Time_redshift::time_to_z(double t)
{
    //t: this is the lookback time
    //but I need: t_real. First: need to calculate t_real
    double t_real = Time_of_the_Universe - t;
    //this will give back the 'a' value
    double a = time_to_z_interpolator.genericInterpolationY(t_real);
    //z can be calculate like this:
    double z = 1.0/a - 1.0;
    return z;
}

double Time_redshift::z_to_time(double z)
{
    //I will need 'a' value for the interpolation
    double a = 1.0/(z + 1.0);
    //The interpolation will give back the real time
    double t_real = time_to_z_interpolator.genericInterpolationX(a);
    //the lookback time based on it:
    double t = Time_of_the_Universe - t_real;
    return t;
}

/********************************************************************
 ********************************************************************
    Calculated numbers based in read in parameters
*********************************************************************
********************************************************************/

Calculated_Numbers_Based_on_read_in_parameters::Calculated_Numbers_Based_on_read_in_parameters(Read_In_Parameters& params/*, Time_redshift& t_z*/) : param(params)/*, time_z(t_z)*/
{
    //These also include a sigma_d,0 multiplier, but it will falls out
    M_star = 2.0 * M_PI * param.rd * param.rd;
    //rho_star: it was integrated from z=-inf to z=inf, it pronounces the 2*zd
    rho_star = std::exp(-param.r_Sun/param.rd);/* / (2.0 * param.zd);*/
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
    //initialize integral class
    integral.set_eps(1e-3);
    integral.set_newfx(*param.rate_function);
    integral.set_xstart(0.0);
    integral.set_xend(param.simulation_Time);
    integral.set_n(0);
    
    //run qtrap and calculate the number_of_events_d
    std::cout << "Calculating the number of events" << std::endl;
    number_of_events_d = integral.qtrap();
    //integral resulted number of events/kpc^3
    //now: need multiple it with the volume of the investigated part
    number_of_events_d *= 2.0 * M_PI * param.r_Sun * param.width;
    
    //calculate the number_of_events (int)
    number_of_events = int(number_of_events_d);
    std::cout << "The number of events: " << number_of_events_d << std::endl;
}


void Calculated_Numbers_Based_on_read_in_parameters::init()
{
    //set r0 rate density at rate function
    param.rate_function->set_r0_rate_density_with_rate_density(rate_to_rate_density(param.r0_rate));
    //set Ni calculation method based on the string in param
    if("Hotokezaka" == param.ni_calculation_method)
    {
        Ni_calc.reset(new Ni_Hotokezaka(*this));
    }
    else if("Measurement" == param.ni_calculation_method)
    {
        Ni_calc.reset(new Ni_Measurement(*this));
    }
    else if(IsDouble(param.ni_calculation_method))
    {
        Ni_calc.reset(new Ni_Read_In(*this));
    }
    else
    {
        throw std::runtime_error("Unknown Ni calculation method: "+ param.ni_calculation_method);
    }

    //calculate taumix for r0 and Ni
    taumix = 300.0 * std::pow(param.r0_rate/10.0, -0.4) * std::pow(param.alpha/0.1,-0.6) * 
             std::pow(param.vt/7.0, -0.6) * std::pow(param.h_scale/0.2,-0.6);
    Ni = Ni_calc->calculate_Ni();
    
    //1e-3: convert 1/Gyr to 1/Myr
    D = param.alpha * (param.vt / 7.0) * (param.h_scale / 0.2) * 1e-3;      //kpc^2/Myr
    Time_of_the_Universe = param.time_z.time_of_universe();
    //second: calculate the number of events
    calculate_number_of_events();
        
    //third: calculate cumulative_distribution function
        //first: create two vectors
        std::vector<double> cumulative_dist_values;
        std::vector<double> time_values_for_cumulative_dist;
        double sum, t_last;
        t_last = 0.0;
        sum = 0.0;
        std::cout << "Started to create cumulative distribution function" << std::endl;
        time_values_for_cumulative_dist.push_back(0.0);
        cumulative_dist_values.push_back(0.0);
        //calculate the integral of rate/number_of_events_d and store on the two vectors
        for(double t = 0.6; t <= param.simulation_Time; t += 0.6)
        {
            integral.set_xstart(t_last);
            integral.set_xend(t);
            sum += integral.qtrap()/number_of_events_d * 2.0 * M_PI * param.r_Sun * param.width;
            cumulative_dist_values.push_back(sum);
            time_values_for_cumulative_dist.push_back(t);
            t_last = t;
        }
        //initialize the cumulative_distribution_interpolator
        cumulative_distribution_interpolator.init(cumulative_dist_values, time_values_for_cumulative_dist);
        std::cout << "Cumulative distribution function done" << std::endl;
}

double Calculated_Numbers_Based_on_read_in_parameters::interpolate_random_number_to_time(double myrand)
{
    return cumulative_distribution_interpolator.genericInterpolationX(myrand);
}


/**************************************************************
 * ************************************************************
 *          Fx function to create interploate file of
 *          the exchange between time and redshift
 *          and the creation of this file
 * ************************************************************
**************************************************************/
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
            integral.set_xend(a_upper_integral_limit);
        }
        last_integral_value += integral.qtrap();
        myout.write((char*)(&a_upper_integral_limit),sizeof(a_upper_integral_limit));
        myout.write((char*)(&last_integral_value),sizeof(last_integral_value));
        a_under_integral_limit += a_step;
    }
    myout.close();
}

/*****************************************
 *****************************************
        Rate density functions
******************************************
*****************************************/

void Rate_density::set_r0_rate_density_with_rate_density(double r0)
{
    r0_rate_density = r0;
    calc_normalize_factor();
}

    /************Wanderman Rate************/

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

    double WandermanRate::get_rate_density_at_t(double t)
    {
        double z = time_z.time_to_z(t);
        return get_rate_density_at_z(z);
    }

    /********************Hopkins rate**************/

    HopkinsRate::HopkinsRate(Time_redshift& t_z) : time_z(t_z)
    {
        a = 0.0170;
        b = 0.13;
        c = 3.3;
        d = 5.3;
        h = constants::H0 / 100.0;
    }

    double HopkinsRate::get_rate_density_at_z(double z)
    {
        double rate_density;
        rate_density = (a + b*z) * h/(1.0+std::pow(z/c,d)) * normalize_factor;
        return rate_density;
    }

    double HopkinsRate::get_rate_density_at_t(double t)
    {
        double z = time_z.time_to_z(t);
        return get_rate_density_at_z(z);
    }

    void HopkinsRate::calc_normalize_factor()
    {
        normalize_factor = r0_rate_density/(a*h);
    }

/***Ni calculation functions*/

    double Ni_Read_In::calculate_Ni()
    {
        //return calc_param.Ni;
        return stod(calc_param.param.ni_calculation_method);
    }

    double Ni_Hotokezaka::calculate_Ni()
    {
        return 3.9540983606557375e+49;
    }

    void Ni_Measurement::calculate_N_Pu()
    {
        //std::cout << constants::Today_Pu_abundance << "\t" << calc_param.taumix << std::endl;
        N_Pu = constants::Today_Pu_abundance / (calc_param.rate_to_rate_density(calc_param.param.r0_rate) * constants::tau_Pu * std::exp(-calc_param.taumix/(2.0*constants::tau_Pu)));
    }

    double Ni_Measurement::calculate_Ni()
    {
        calculate_N_Pu();
        double res = calc_param.param.element_initial_production_ratio * N_Pu * 1.0/(calc_param.param.Pu_initial_production_ratio);
        //std::cout << N_Pu << "\t" << res << std::endl;
        return res;
    }

Create_events_and_calc_number_density::Create_events_and_calc_number_density(Calculated_Numbers_Based_on_read_in_parameters& calc_par) : calculated_parameters(calc_par)
{
    const_for_Kj_1 = 4.0 * M_PI * calculated_parameters.D;
    const_for_Kj_2 = 8.0 * M_PI * calculated_parameters.param.h_scale * calculated_parameters.D;
    rand_number_0_1 = std::uniform_real_distribution<double>(0.0,1.0);
    rand_laplace = boost::random::laplace_distribution<double>(0.0, calculated_parameters.param.h_scale);
    double time;
    for(time = 0.0; time <= calculated_parameters.param.simulation_Time; time += calculated_parameters.param.sampleDt)
    {
        sampling_time_points.push_back(time);
        number_densites.push_back(0.0);
        median_number_densities.push_back(0.0);
    }
    if(sampling_time_points[sampling_time_points.size()-1] != calculated_parameters.param.simulation_Time)
    {
        sampling_time_points.push_back(calculated_parameters.param.simulation_Time);
        number_densites.push_back(0.0);
        median_number_densities.push_back(0.0);
    }
}

randomEvent Create_events_and_calc_number_density::create_random_event()
{
    randomEvent event;
    std::random_device rd;
    std::mt19937 mt(rd());
    double circumf = calculated_parameters.param.r_Sun * 2.0 * M_PI;
    event.time = calculated_parameters.interpolate_random_number_to_time(rand_number_0_1(mt));
    event.x = circumf * rand_number_0_1(mt) - 0.5 * circumf;
    event.y = calculated_parameters.param.width * rand_number_0_1(mt) - 0.5 * calculated_parameters.param.width;
    event.z = rand_laplace(mt);
    event.distance = std::sqrt((event.x*event.x + event.y*event.y + event.z * event.z));
    return event;
}

double Create_events_and_calc_number_density::calc_Kj(double delta_time)
{
    double left, right, Kj;
    left = std::pow(const_for_Kj_1*delta_time, 1.5);
    right = const_for_Kj_2 * delta_time;
    Kj = std::min(left, right);
    return Kj;
}

void Create_events_and_calc_number_density::calc_number_density_for_an_event()
{
    //first: need an event
    randomEvent event = create_random_event();
    //some usable variable
    double const_at_exp = - event.distance*event.distance/(4.0 * calculated_parameters.D);
    double delta_tj, number_density;
    //second: from the first sampling_time_points see all, if it < event.time, need calculating, else break
    for(long unsigned int i = 0; i < sampling_time_points.size(); i++)
    {
        if(sampling_time_points[i] <= event.time)
        {
            delta_tj = event.time - sampling_time_points[i];
            number_density = calculated_parameters.Ni/calc_Kj(delta_tj) * std::exp(const_at_exp/delta_tj - delta_tj/calculated_parameters.param.tau);
            number_densites[i] += number_density;
        }
        else
        {
            break;
        }
    }
}

void Create_events_and_calc_number_density::allEvent_number_densities()
{
    std::ofstream res;
    res.open(calculated_parameters.param.out_file);
    for(long unsigned int i = 0; i < sampling_time_points.size(); i++)
    {
        res << sampling_time_points[i] << "\t";
    }
    res << std::endl;
    std::cout << "Started to calculate median density" << std::endl;
    for(long unsigned int i = 0; i < median_number_densities.size(); i++)
    {
        double rate_dens_i = calculated_parameters.param.rate_function->get_rate_density_at_t(sampling_time_points[i]);
        double rate_i = calculated_parameters.rate_density_to_rate(rate_dens_i);
        double taumix = 300.0 * std::pow(rate_i/10.0, -0.4) * std::pow(calculated_parameters.param.alpha/0.1,-0.6) * 
             std::pow(calculated_parameters.param.vt/7.0, -0.6) * std::pow(calculated_parameters.param.h_scale/0.2,-0.6);
        median_number_densities[i] = calculated_parameters.Ni * rate_dens_i * 
                calculated_parameters.param.tau * std::exp(-taumix/(2.0*calculated_parameters.param.tau));
    }
    for(long unsigned int i = 0; i < median_number_densities.size(); i++)
    {
        res << median_number_densities[i] << "\t";
    }
    res << std::endl;
    std::cout << "Started to create and calculate random events" << std::endl;
    for(int i = 0; i < calculated_parameters.param.number_of_runs; i++)
    {
        //set all number_densities to 0
        std::fill(number_densites.begin(), number_densites.end(), 0.0);
        for(int j = 0; j < calculated_parameters.get_number_of_events(); j++)
        {
            calc_number_density_for_an_event();
        }
        //write to output file
        for(long unsigned int j = 0; j < number_densites.size(); j++)
        {
            res << number_densites[j] << "\t";
        }
        res << std::endl;
    }
    res.close();
}