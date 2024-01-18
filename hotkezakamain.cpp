
#include "hotokezaka.hpp"


int main(int argc, char* argv[])
{
    /*if(!std::filesystem::exists("interpolate_table.dat"))
    {
        create_file_for_interpolation file;
        file = create_file_for_interpolation();
        file.save_data_to_file();
    }
    Interpolator interpolate_a_time;
    interpolate_a_time.init("interpolate_table.dat");*/
    try{
        Read_In_Parameters readInParam;
        readInParam.read_parameter_file();

        Calculated_Numbers_Based_on_read_in_parameters calcParams(readInParam);
        //calcParams.init();

        /*for(double i = 0.0; i < 6000.0; i+=100.0)
        {
            std::cout << i << "\t" << calcParams.time_to_z(i) << std::endl;
        }*/

        Create_events_and_calc_number_density calc(calcParams);
        calc.allEvent_number_densities();
    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << '\n';
    }
    return 0;
}