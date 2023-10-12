
#include "hotokezaka.hpp"


int main(int argc, char* argv[])
{
    if(!std::filesystem::exists("interpolate_table.dat"))
    {
        create_file_for_interpolation file;
        file = create_file_for_interpolation();
        file.save_data_to_file();
    }
    Interpolator interpolate_a_time;
    interpolate_a_time.init("interpolate_table.dat");
    
    return 0;
}