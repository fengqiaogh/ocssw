/*******************************************************************************
 *
 * NAME: RsViirs.cpp
 *
 * DESCRIPTION: Running the program on a granule will modify the data
 * in-place ..
 *
 * Based on a technique described in:
 * Remote Sensing 2016, 8(1), 79
 * Article Titled "Improved VIIRS and MODIS SST Imagery"
 * By Irina Gladkova, Alexander Ignatov, Fazlul Shahriar, Yury Kihai,
 * Don Hillger and Boris Petrenko
 *
 *
 *  Created on: February 17, 2019
 *      Author: Sam Anderson
 *
 *******************************************************************************/


#include "resam_viirs/RsViirs.h"

#include <netcdf>
#include <DDOptions.h>

using namespace std;
using namespace netCDF;
using namespace netCDF::exceptions;

// column break points
short SBP[NBREAKS+1] = {
    0, 5, 87, 170, 358, 567, 720, 850, 997, 1120, 1275, 1600
};

// relative row that the pixel comes from
short SORT_FIRST[NDETECTORS][NBREAKS] = {
    {0,   0,  0,  0,  0,  0,  0,  0,  0,  0, 0},
    {0,   0,  0,  0,  0,  0,  0,  0,  0,  0, 0},
    {0,   0,  0,  0,  0,  0,  0,  0,  0,  0, 0},
    {0,   0,  0,  0,  0,  0,  0,  0,  0,  0, 0},
    {0,   0,  0,  0,  0,  0,  0,  0,  0,  0, 0},
    {0,   0,  0,  0,  0,  0,  0,  0,  0,  0, 0},
    {0,   0,  0,  0,  0,  0,  0,  0,  0,  0, 0},
    {0,   0,  0,  0,  0,  0,  0,  0,  0,  0, 0},
    {+8, +8, +8,  0,  0,  0,  0,  0,  0,  0, 0},
    {+8, -1, -1, +7,  0,  0,  0,  0,  0,  0, 0},
    {-2, +7, +7, -1, +6,  0,  0,  0,  0,  0, 0},
    {+7, -2, -2, +6, -1, +5,  0,  0,  0,  0, 0},
    {-3, +6, +6, -2, +5, -1, +4,  0,  0,  0, 0},
    {+6, -3, -3, +5, -2, +4, -1, +3,  0,  0, 0},
    {-4, +5, +5, -3, +4, -2, +3, -1, +2,  0, 0},
    {+5, -4, -4, +4, -3, +3, -2, +2, -1, +1, 0},
};

// relative row that the pixel comes from
short SORT_MID[NDETECTORS][NBREAKS] = {
    {-5,  +4, +4, -4, +3, -3, +2, -2, +1, -1, 0},
    {+4,  -5, -5, +3, -4, +2, -3, +1, -2,  0, 0},
    {-6,  +3, +3, -5, +2, -4, +1, -3,  0,  0, 0},
    {+3,  -6, -6, +2, -5, +1, -4,  0,  0,  0, 0},
    {-7,  +2, +2, -6, +1, -5,  0,  0,  0,  0, 0},
    {+11, -7, -7, +1, -6,  0,  0,  0,  0,  0, 0},
    {+1,  +1, +1, -7,  0,  0,  0,  0,  0,  0, 0},
    {-9,  +9, -8,  0,  0,  0,  0,  0,  0,  0, 0},
    {+9,  -9, +8,  0,  0,  0,  0,  0,  0,  0, 0},
    {-1,  -1, -1, +7,  0,  0,  0,  0,  0,  0, 0},
    {-11, +7, +7, -1, +6,  0,  0,  0,  0,  0, 0},
    {+7,  -2, -2, +6, -1, +5,  0,  0,  0,  0, 0},
    {-3,  +6, +6, -2, +5, -1, +4,  0,  0,  0, 0},
    {+6,  -3, -3, +5, -2, +4, -1, +3,  0,  0, 0},
    {-4,  +5, +5, -3, +4, -2, +3, -1, +2,  0, 0},
    {+5,  -4, -4, +4, -3, +3, -2, +2, -1, +1, 0},
};

// relative row that the pixel comes from
short SORT_LAST[NDETECTORS][NBREAKS] = {
    {-5, +4, +4, -4, +3, -3, +2, -2, +1, -1, 0},
    {+4, -5, -5, +3, -4, +2, -3, +1, -2,  0, 0},
    {-6, +3, +3, -5, +2, -4, +1, -3,  0,  0, 0},
    {+3, -6, -6, +2, -5, +1, -4,  0,  0,  0, 0},
    {-7, +2, +2, -6, +1, -5,  0,  0,  0,  0, 0},
    {+2, -7, -7, +1, -6,  0,  0,  0,  0,  0, 0},
    {-8, +1, +1, -7,  0,  0,  0,  0,  0,  0, 0},
    {-8, -8, -8,  0,  0,  0,  0,  0,  0,  0, 0},
    {0,   0,  0,  0,  0,  0,  0,  0,  0,  0, 0},
    {0,   0,  0,  0,  0,  0,  0,  0,  0,  0, 0},
    {0,   0,  0,  0,  0,  0,  0,  0,  0,  0, 0},
    {0,   0,  0,  0,  0,  0,  0,  0,  0,  0, 0},
    {0,   0,  0,  0,  0,  0,  0,  0,  0,  0, 0},
    {0,   0,  0,  0,  0,  0,  0,  0,  0,  0, 0},
    {0,   0,  0,  0,  0,  0,  0,  0,  0,  0, 0},
    {0,   0,  0,  0,  0,  0,  0,  0,  0,  0, 0},
};


/**************************************************************************
 * NAME: RsViirs()
 *
 * DESCRIPTION: Class Constructor
 *
 *************************************************************************/

RsViirs::RsViirs() {
}

RsViirs::RsViirs( int lines, int pixels ) {
    lines_ = lines;
    pixels_ = pixels;
}

/**************************************************************************
 * NAME: ~RsViirs()
 *
 * DESCRIPTION: Class Destructor
 *
 *************************************************************************/

RsViirs::~RsViirs() {
}

/**************************************************************************
 * NAME: process()
 *
 * DESCRIPTION: Processes data and produces output values for granule
 *
 *************************************************************************/

int RsViirs::process() {
    int status = DT_SUCCESS;

    string filepath_l1b = get_option(INPUT_L1B);
    string filepath_geo = get_option(INPUT_GEO);

    if (filepath_l1b.empty() || filepath_geo.empty()) {
        return DT_FAIL;
    }

    NcFile* nc_io;
    try {
        nc_io = new NcFile(filepath_l1b, NcFile::read);
    } catch (NcException& e) {
        e.what();
        cerr << "RsInput:: Failure opening VIIRS L1B input file: "
                        + filepath_l1b << endl;
        return DT_FAIL;
    }
    NcDim line_dim = nc_io->getDim("number_of_lines");
    NcDim pixel_dim = nc_io->getDim("number_of_pixels");

    lines_ = line_dim.getSize();
    pixels_ = pixel_dim.getSize();
    delete nc_io;

    status = generate_sort_index();

    sio* sa = new sio;
    memset( sa, 0, sizeof(sio));
    usio* usa = new usio;
    memset( usa, 0, sizeof(usio));
    fio* fa = new fio;
    memset( fa, 0, sizeof(fio));

    try {
        nc_io = new NcFile(filepath_geo, NcFile::write);
    } catch (NcException& e) {
        e.what();
        cerr << "RsViirs:: Failure opening VIIRS Geolocation input file: "
                        + filepath_geo << endl;
        return DT_FAIL;
    }

    NcGroup nc_group = nc_io->getGroup("geolocation_data");
    NcVar nc_var = nc_group.getVar("quality_flag");
    nc_var.getVar(&sa->in[0][0]);
    resort(sa);
    nc_var.putVar(&sa->out[0][0]);
    nc_var = nc_group.getVar("height");
    nc_var.getVar(&sa->in[0][0]);
    resort(sa);
    nc_var.putVar(&sa->out[0][0]);
    nc_var = nc_group.getVar("range");
    nc_var.getVar(&sa->in[0][0]);
    resort(sa);
    nc_var.putVar(&sa->out[0][0]);
    nc_var = nc_group.getVar("sensor_azimuth");
    nc_var.getVar(&sa->in[0][0]);
    resort(sa);
    nc_var.putVar(&sa->out[0][0]);
    nc_var = nc_group.getVar("sensor_zenith");
    nc_var.getVar(&sa->in[0][0]);
    resort(sa);
    nc_var.putVar(&sa->out[0][0]);
    nc_var = nc_group.getVar("solar_azimuth");
    nc_var.getVar(&sa->in[0][0]);
    resort(sa);
    nc_var.putVar(&sa->out[0][0]);
    nc_var = nc_group.getVar("solar_zenith");
    nc_var.getVar(&sa->in[0][0]);
    resort(sa);
    nc_var.putVar(&sa->out[0][0]);
    nc_var = nc_group.getVar("latitude");
    nc_var.getVar(&fa->in[0][0]);
    resort(fa);
    nc_var.putVar(&fa->out[0][0]);
    nc_var = nc_group.getVar("longitude");
    nc_var.getVar(fa);
    resort(fa);
    nc_var.putVar(&fa->out[0][0]);
    delete nc_io;

    try {
        nc_io = new NcFile(filepath_l1b, NcFile::write);
    } catch (NcException& e) {
        e.what();
        cerr << "RsViirs:: Failure opening VIIRS L1B input file: "
                        + filepath_l1b << endl;
        return DT_FAIL;
    }
    nc_group = nc_io->getGroup("observation_data");
    for (int ib = 0; ib < NMBANDS; ib++) {
        char str[20];
        sprintf(str, "M%02d", ib + 1);
        nc_var = nc_group.getVar(string(str));
        nc_var.getVar(&usa->in[0][0]);
        resort(usa);
        fill_the_fills(usa);
        nc_var.putVar(&usa->out[0][0]);
        sprintf(str, "M%02d_quality_flags", ib + 1);
        nc_var = nc_group.getVar(string(str));
        nc_var.getVar(&sa->in[0][0]);
        resort(sa);
        nc_var.putVar(&sa->out[0][0]);
    }
    delete nc_io;

    delete sa;
    delete usa;
    delete fa;

    return status;
}

/**************************************************************************
 * NAME: generate_sort_index()
 *
 * DESCRIPTION: Generate a image of latitude sorting indices.
 *
 *************************************************************************/

int RsViirs::generate_sort_index() {

    int ND = NDETECTORS;
    for (int i = 0; i < NBREAKS; i++) {
        for (int ix = SBP[i]; ix < SBP[i+1]; ix++) {
            for (int iy = 0; iy < ND; iy++) {
                si_[iy][ix] = iy + SORT_FIRST[iy][i];
                si_[iy][pixels_-ix-1] = iy + SORT_FIRST[iy][i];
            }
            for (int iy = ND; iy < lines_ - ND; iy++) {
                si_[iy][ix] = iy + SORT_MID[iy % ND][i];
                si_[iy][pixels_-ix-1] = iy + SORT_MID[iy % ND][i];
            }
            for (int iy = lines_ - ND; iy < lines_; iy++) {
                si_[iy][ix] = iy + SORT_LAST[iy % ND][i];
                si_[iy][pixels_-ix-1] = iy + SORT_LAST[iy % ND][i];
            }
        }
    }
    return DT_SUCCESS;
}

/**************************************************************************
 * NAME: resort()
 *
 * DESCRIPTION: Returns the sorted image.
 * Si is the image of sort indices.
 *
 *************************************************************************/

int RsViirs::resort(sio* io) {
    for (int iy = 0; iy < lines_; iy++) {
        for (int ix = 0; ix < pixels_; ix++) {
            io->out[iy][ix] = io->in[si_[iy][ix]][ix];
        }
    }
    return DT_SUCCESS;
}

int RsViirs::resort(usio* io) {
    for (int iy = 0; iy < lines_; iy++) {
        for (int ix = 0; ix < pixels_; ix++) {
            io->out[iy][ix] = io->in[si_[iy][ix]][ix];
        }
    }
    return DT_SUCCESS;
}

int RsViirs::resort(fio* io) {
    for (int iy = 0; iy < lines_; iy++) {
        for (int ix = 0; ix < pixels_; ix++) {
            io->out[iy][ix] = io->in[si_[iy][ix]][ix];
        }
    }
    return DT_SUCCESS;
}

/**************************************************************************
 * NAME: fill_the_fills()
 *
 * DESCRIPTION: Replace fill values with average of above and below.
 *
 *************************************************************************/

int RsViirs::fill_the_fills(usio* io) {
    const unsigned short filltest = SOUB_UINT16_FILL;

    for (int iy = 1; iy < lines_ - 1; iy++) {
        for (int ix = 0; ix < pixels_; ix++) {
            if (io->out[iy][ix] >= filltest) {
                if (io->out[iy - 1][ix] < filltest &&
                    io->out[iy + 1][ix] < filltest) {
                    io->out[iy][ix] = io->out[iy - 1][ix] / 2 +
                                      io->out[iy + 1][ix] / 2;
                } else if (io->out[iy - 1][ix] < filltest) {
                    io->out[iy][ix] = io->out[iy - 1][ix];
                }  else if (io->out[iy + 1][ix] < filltest) {
                    io->out[iy][ix] = io->out[iy + 1][ix];
                }
            }
        }
    }
    for (int ix = 0; ix < pixels_; ix++) {
        if (io->out[0][ix] >= filltest) {
            if (io->out[1][ix] < filltest) {
                io->out[0][ix] = io->out[1][ix];
            }
        }
    }
    for (int ix = 0; ix < pixels_; ix++) {
        if (io->out[lines_-1][ix] >= filltest) {
            if (io->out[lines_-2][ix] < filltest) {
                io->out[lines_-1][ix] = io->out[lines_-2][ix];
            }
        }
    }

    return DT_SUCCESS;
}

int RsViirs::fill_the_fills(fio* io) {
    const unsigned short filltest = SOUB_UINT16_FILL;

    for (int iy = 1; iy < lines_ - 1; iy++) {
        for (int ix = 0; ix < pixels_; ix++) {
            if (io->out[iy][ix] < filltest) {
                if (io->out[iy - 1][ix] > filltest &&
                    io->out[iy + 1][ix] > filltest) {
                    io->out[iy][ix] = io->out[iy - 1][ix] / 2 +
                                      io->out[iy + 1][ix] / 2;
                } else if (io->out[iy - 1][ix] > filltest) {
                    io->out[iy][ix] = io->out[iy - 1][ix];
                }  else if (io->out[iy + 1][ix] > filltest) {
                    io->out[iy][ix] = io->out[iy + 1][ix];
                }
            }
        }
    }
    for (int ix = 0; ix < pixels_; ix++) {
        if (io->out[0][ix] < filltest) {
            if (io->out[1][ix] > filltest) {
                io->out[0][ix] = io->out[1][ix];
            }
        }
    }
    for (int ix = 0; ix < pixels_; ix++) {
        if (io->out[lines_-1][ix] < filltest) {
            if (io->out[lines_-2][ix] > filltest) {
                io->out[lines_-1][ix] = io->out[lines_-2][ix];
            }
        }
    }

    return DT_SUCCESS;
}

