#ifndef EXPAND_3D_HPP
#define EXPAND_3D_HPP
#define BAD_MAX 999.9

#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <string>
#include <boost/algorithm/string.hpp>
#include <netcdf>
#include "../l2bin/l2bin_input.h"
#include "genutils.h"
#include "readL2scan.h"


void set_mapper(const std::string &input_file_name, std::string &products_requested, std::vector<std::string> &products_requested_separated);
bool found_wv(int);
void set_prodname_3d_to_l2bin(const instr &input, const std::vector<std::string> &prodparam, l2_prod &l2_str, std::vector<std::string> &l2_prodname, std::vector<std::string> &l3_prodname, std::vector<int32_t> &thirdDimId, std::vector<float> &min_value, std::vector<float> &max_value);
bool set_l2_flags_use(const std::string &flagsuse);
#endif