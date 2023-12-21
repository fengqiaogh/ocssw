#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <netcdf>

using namespace std;
using namespace netCDF;
using namespace netCDF::exceptions;

bool is_NCDF(const string filepath) {
    try {
        NcFile nc(filepath, NcFile::read);
        nc.close();
        return true;
    } catch (NcException const & e) {
        return false;
    }
}

vector<string> ncfiles(const char* filepath) {
    /* Returns a string vector of NetCDF filenames, listed in filepath.
       If filepath is a NetCDF file, vector will contain only that file.
    */
    vector<string> filelist;

    if (is_NCDF(filepath)) {
        filelist.push_back(filepath);
    } else {
        string line;
        ifstream infile(filepath);
        if (infile) {
            while (getline(infile,line)) {
                if (is_NCDF(line))
                    filelist.push_back(line);
            }
            infile.close();
        }
    }

    return(filelist);
}
