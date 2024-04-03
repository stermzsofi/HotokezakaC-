#include "hotokezaka.hpp"


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
        else if("zd" == name)
        {
            zd = stod(value);
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
    Time_of_the_Universe = time_to_z_interpolator.genericInterpolationX(1.0);
}

//The time of the universe is the time at a=1
double Time_redshift::time_of_universe()
{
    //return time_to_z_interpolator.genericInterpolationX(1.0);
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

Calculated_Numbers_Based_on_read_in_parameters::Calculated_Numbers_Based_on_read_in_parameters(Read_In_Parameters& params/*, Time_redshift& t_z*/) : param(params)/*, time_z(t_z)*/
{
    //rho_star: it was integrated from z=-inf to z=inf, it pronounces the 2*zd
    //These also include a sigma_d,0 multiplier, but it will falls out
    M_star = 2.0 * M_PI * param.rd * param.rd;
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
    //integral class inicializálása
    //qtrap meghívása
    //eredmény egészre kerekítve beírása number_of_events változóba
}


void Calculated_Numbers_Based_on_read_in_parameters::init()
{
    //set r0 rate density at rate function
    param.rate_function->set_r0_rate_density_with_rate_density(rate_to_rate_density(param.r0_rate));
    //param.rate_function->set_r0_rate_density_with_rate_density(param.r0_rate);
    Ni = 1.0;       //later need to replace some calculations
    //Ni = (Pu/U)0 * NU = 0.4 * (XU * Mej * Msun)/(238 * 1.66e-24 g) = 
    // = 0.4 * YU * Mej * 1.1978e57  = 0.4 * 7.3e-13 * 1.1978e57 * 0.1/R0
    //Ni = 0.4 * 7.3e-13 * 1.1978e57 * 0.1/param.r0_rate; 
    
    
    taumix = 300.0 * std::pow(param.r0_rate/10.0, -0.4) * std::pow(param.alpha/0.1,-0.6) * 
             std::pow(param.vt/7.0, -0.6) * std::pow(param.h_scale/0.2,-0.6);
    //Ni = <ni>m / (rate_density_0 * taui * exp(-taumix/(2*taui)))
    //<ni>m comes from measurements: today Pu value is 6e-17 cm^-3, but we need kpc^-3
    Ni = constants::Today_Pu_abundance / (rate_to_rate_density(param.r0_rate) * param.tau * std::exp(-taumix/(2.0*param.tau)));
    Ni = 3.9540983606557375e+49;
    //Ni = 3.540983606557377e+50;
    //1e-3: convert 1/Gyr to 1/Myr
    // dimension : kpc^2/Myr
    D = param.alpha * (param.vt / 7.0) * (param.h_scale / 0.2) * 1e-3;
    //+1: store the time of the universe
    //Time_of_the_Universe = time_to_z_interpolator.genericInterpolationX(1.0);
    Time_of_the_Universe = param.time_z.time_of_universe();
    //second: calculate the number of events
        //simulation_Time -> z_simulation
        //double z_simulation = param.time_z.time_to_z(param.simulation_Time);
        //integral initialize (set rate, ranges)
        //integral = Quad_Trapezoidal(*param.rate_function, 0.0, z_simulation);
        integral.set_eps(1e-3);
        integral.set_newfx(*param.rate_function);
        integral.set_xstart(0.0);
        //integral.set_xend(z_simulation);
        integral.set_xend(param.simulation_Time);
        integral.set_n(0);
        //run qtrap and calculate the number_of_events_d
        std::cout << "Calculating the number of events" << std::endl;
        std::cout << param.rate_function->calc_Fx(300.0) << std::endl;
        number_of_events_d = integral.qtrap();
        //integral resulted number of events/kpc^3
        //now: need multiple it with the volume of the investigated part
        number_of_events_d *= 2.0 * M_PI * param.r_Sun * param.width;
        //number_of_events_d = rate_density_to_rate(number_of_events_d);
        //calculate the number_of_events
        number_of_events = int(number_of_events_d);
        std::cout << "The number of events: " << number_of_events_d << std::endl;
    //third: calculate cumulative_distribution function
        //first: create two vectors
        std::vector<double> cumulative_dist_values;
        std::vector<double> time_values_for_cumulative_dist;
        double z, z_last, sum, t_last;
        //int i = 0;
        z_last = 0.0;
        t_last = 0.0;
        sum = 0.0;
        //ezt itt nagyon át kell gondolni: rate vagy rate_density, hogyan is lesz belőle
        //eloszlás fv
        std::cout << "Started to create cumulative distribution function" << std::endl;
        time_values_for_cumulative_dist.push_back(0.0);
        cumulative_dist_values.push_back(0.0);
        //calculate the integral of rate/number_of_events_d and store on the two vectors
        for(double t = 0.6; t <= param.simulation_Time; t += 0.6)
        {
            //z = param.time_z.time_to_z(t);
            //integral.set_xstart(z_last);
            integral.set_xstart(t_last);
            //integral.set_xend(z);
            integral.set_xend(t);
            sum += integral.qtrap()/number_of_events_d * 2.0 * M_PI * param.r_Sun * param.width;
            cumulative_dist_values.push_back(sum);
            time_values_for_cumulative_dist.push_back(t);
            //i++;
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

/*double Calculated_Numbers_Based_on_read_in_parameters::time_to_z(double t)
{
    //t: this is the lookback time
    //but I need: t_real. First: need to calculate t_real
    double t_real = Time_of_the_Universe - t;
    //this will give back the a value
    double a = time_to_z_interpolator.genericInterpolationY(t_real);
    //z can be calculate like this:
    double z = 1.0/a - 1.0;
    return z;
}*/

/*double Calculated_Numbers_Based_on_read_in_parameters::z_to_time(double z)
{
    //I will need 'a' value for the interpolation
    double a = 1.0/(z + 1.0);
    //The interpolation will give back the real time
    double t_real = time_to_z_interpolator.genericInterpolationX(a);
    //the lookback time based on it:
    double t = Time_of_the_Universe - t_real;
    return t;
}*/
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
    //std::cout << a_under_integral_limit << "\t" << last_integral_value << std::endl;
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
        //std::cout << a_upper_integral_limit << "\t" << last_integral_value << std::endl;
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
    std::cout << "r0_rate_density = " << r0_rate_density << " norm factor = " << normalize_factor << std::endl; 
}

/*WandermanRate::WandermanRate(double r0)
{
    r0_rate_density = r0;
    normalize_factor = r0_rate_density / std::exp(-0.9/0.39);
}*/

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

Create_events_and_calc_number_density::Create_events_and_calc_number_density(Calculated_Numbers_Based_on_read_in_parameters& calc_par) : calculated_parameters(calc_par)
{
    //std::random_device rd;
    //mt.seed(std::random_device());
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
    //number_densites.reserve(sampling_time_points.size());
    //std::fill(number_densites.begin(), number_densites.end(), 0.0);
}

randomEvent Create_events_and_calc_number_density::create_random_event()
{
    randomEvent event;
    std::random_device rd;
    std::mt19937 mt(rd());
    double circumf = calculated_parameters.param.r_Sun * 2.0 * M_PI;
    //std::uniform_real_distribution<double> rand_number_0_1(0.0,1.0);
    //boost::math::laplace_distribution<double> rand_laplace(0,10);
    event.time = calculated_parameters.interpolate_random_number_to_time(rand_number_0_1(mt));
    event.x = circumf * rand_number_0_1(mt) - 0.5 * circumf;
    event.y = calculated_parameters.param.width * rand_number_0_1(mt) - 0.5 * calculated_parameters.param.width;
    //double temp = rand_number_0_1(mt) * calculated_parameters.param.zd - 0.5 * calculated_parameters.param.zd; //majd ki kell cserélni laplace-ra
    event.z = rand_laplace(mt);
    event.distance = std::sqrt((event.x*event.x + event.y*event.y + event.z * event.z));
    /*for(double i = 0; i < 1; i+=0.2)
    {
        double temp = boost::math::cdf(rand_laplace,i);
    }*/
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
    //std::cout << event.distance << "\t" << event.time << std::endl;
    double const_at_exp = - event.distance*event.distance/(4.0 * calculated_parameters.D);
    double delta_tj, number_density;
    //second: from the first sampling_time_points see all, if it < event.time, need calculating, else break
        //calculate the limit, where need to calculate
        //int range_limit = na jó, majd később
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
    /*std::cout << "Started to create and calculate random events" << std::endl;
    //create number_of_events random event and calculate the number densities
    for(int i = 0; i < calculated_parameters.get_number_of_events(); i++)
    {
        calc_number_density_for_an_event();
    }
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
    //write out the results to a file
    std::ofstream res;
    res.open("results.dat");
    res << "#time (BP)\tnumber_density (1/kpc^3)\tmedian number density (1/kpc^3)" << std::endl;
    for(long unsigned int i = 0; i < sampling_time_points.size(); i++)
    {
        res << sampling_time_points[i] << "\t" << number_densites[i] << "\t" << median_number_densities[i] << std::endl;
    }
    res.close();*/
    //Create random events and calculate the number density number_of_events times
    //Output file: first line: times
    //Second line: median value
    //After: random examples

    std::ofstream res;
    res.open("results.dat");
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
        //std::cout << "Write" << std::endl;
        for(long unsigned int j = 0; j < number_densites.size(); j++)
        {
            res << number_densites[j] << "\t";
        }
        res << std::endl;
    }
    res.close();
}

/*double WandermanRate::calc_Fx(double x)
{
   return get_rate_density_at_z(x);
}*/