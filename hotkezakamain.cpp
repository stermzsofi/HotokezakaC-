#include "hotokezaka.hpp"

int main(int argc, char* argv[])
{
    create_file_for_interpolation file;
    file = create_file_for_interpolation();
    file.save_data_to_file();
    return 0;
}