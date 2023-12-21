#include "expand3D.hpp"
#include "productInfo.h"
#include "sensorInfo.h"
#include <limits>

namespace
{
    std::unordered_map<std::string, std::string> l3d_to_l2_prod_map;
    std::unordered_map<std::string, std::pair<int, int>> l3d_indexes;
    std::vector<int> wavelengths_3d;
    std::vector<int> wavelengths_all;
    std::unordered_set<std::string> wv_3d_products;
    std::unordered_set<std::string> wv_all_no_3d_products;
    std::unordered_set<std::string> original_l2_products;
    std::vector<int> wv_requested_3d;
    std::vector<int> wv_requested_all_other_no_3d;
    std::unordered_map<std::string, std::pair<float, float>> min_max_values;
    std::string platform, instrument;
    int sensorID;
    std::unordered_set<int> requested_wv;
    bool ini_3d{false};
    bool use_l2_flags{false};
}
bool found_wv(int wv)
{
    return requested_wv.count(wv) > 0;
}
std::vector<std::string> create_min_max_values(const std::vector<std::string> &products_requested_separated)
{
    std::vector<std::string> output_expanded_products;
    productInfo_t *info = allocateProductInfo();

    for (const auto &prod_with_min_max : products_requested_separated)
    {
        std::vector<std::string> prod_data;
        float max_val_preset, min_val_preset;
        boost::algorithm::split(prod_data, prod_with_min_max, boost::is_any_of("="), boost::token_compress_on); // products_requested}
        {
            if(findProductInfo(prod_data.at(0).c_str(), sensorID, info))
            {
                min_val_preset = info->validMin;
                max_val_preset = info->validMax;
            }
            else
            {
                min_val_preset = 0.0;
                max_val_preset =  std::numeric_limits<float>::max() ;
            }
        }
        if (prod_data.size() > 1)
        {
            output_expanded_products.push_back(prod_data.at(0));
            std::vector<std::string> min_max_sep;
            boost::algorithm::split(min_max_sep, prod_data.at(1), boost::is_any_of(":"), boost::token_compress_on); // products_requested}
            if (min_max_sep.size() == 1)
            {
                min_max_values[prod_data.at(0)] = {std::stof(min_max_sep.at(0)), max_val_preset};
            }
            else if (min_max_sep.size() == 2)
            {
                float max_val = min_max_sep.at(1).length() > 0 ? std::stof(min_max_sep.at(1)) : max_val_preset;
                float min_val = min_max_sep.at(0).length() > 0 ? std::stof(min_max_sep.at(0)) : min_val_preset;
                min_max_values[prod_data.at(0)] = {min_val, max_val};
            }
            else
            {
                std::cerr << "Can't parse the product " << prod_with_min_max << std::endl;
                exit(EXIT_FAILURE);
            }
        }
        else
        {
            output_expanded_products.push_back(prod_with_min_max);
            min_max_values[prod_with_min_max] = {min_val_preset, max_val_preset};
        }
    }
    // for(const auto & var :min_max_values )
    // {
    //     std::cout << var.first << " " << var.second.first << " " << var.second.second << std::endl; 
    // }
    freeProductInfo(info);
    return output_expanded_products;
}

/**
 * @brief 
 * The user can supply a 3D product in the form "Lt" or "Rrs" and it has to be expanded
 * @param products_requested_separated - list of all products. 
 */
void expand_l2_requested(std::vector<std::string> &products_requested_separated)
{
    std::unordered_set<std::string> remove;
    std::vector<std::string> add_3d_expanded;
    for (auto &prod : products_requested_separated)
    {
        if (wv_3d_products.count(prod) > 0)
        {
            remove.insert(prod);
            for (const auto &wv : wavelengths_3d)
            {
                auto exp3d_name = prod + "_" + std::to_string(wv);
                add_3d_expanded.push_back(exp3d_name);
                min_max_values[exp3d_name] = min_max_values.at(prod);
            }
        }
        if (wv_all_no_3d_products.count(prod) > 0)
        {
            remove.insert(prod);
            for (const auto &wv : wavelengths_all)
            {
                auto exp_all_name = prod + "_" + std::to_string(wv);
                add_3d_expanded.push_back(exp_all_name);
                min_max_values[exp_all_name] = min_max_values.at(prod);
            }
        }
    }

    if (!remove.empty())
    {
        for (const auto &rem : remove)
        {
            // std::cout << "All wavelength will be binned and sent to output for " << rem << std::endl;
            auto it = std::find(products_requested_separated.begin(), products_requested_separated.end(), rem);
            products_requested_separated.erase(it);
        }
        for (const auto &prod_exp_3d : add_3d_expanded)
        {
            products_requested_separated.push_back(prod_exp_3d);
        }
    }
}

// currently  works only for NC L2
/**
 * @brief
 *  readL2, openL2 uses the maps from this module. If the module is not initialized, the the function read and process an unxpanded L2 products
 * @param input_file_name path to nc file
 * @param products_requested products requested by the user, with delimeters
 * @param products_requested_separated list of the separated products
 */
void set_mapper(const std::string &input_file_name, std::string &products_requested, std::vector<std::string> &products_requested_separated)
{

    ini_3d = true;
    netCDF::NcFile l2_file(input_file_name, netCDF::NcFile::read);
    l2_file.getAtt("instrument").getValues(instrument);
    l2_file.getAtt("platform").getValues(platform);
    sensorID = instrumentPlatform2SensorId(instrument.c_str(),
                                           platform.c_str());
    {
        const auto grp = l2_file.getGroup("geophysical_data");
        const auto vars = grp.getVars();
        for (const auto &var : vars)
        {
            // std::cout << "Alex : Var geophysical_data name debug print " << var.first << std::endl;
            const auto dims = var.second.getDims();
            std::unordered_set<std::string> dim_names;
            for (const auto &dim : dims)
            {
                // std::cout << "Dim name " << dim.getName() << " " << dim.getSize() << std::endl;
                dim_names.insert(dim.getName());
            }
            // std::cout << "\n";
            if (dims.size() == 3 && dim_names.count("wavelength_3d"))
            {
                wv_3d_products.insert(var.first);
            }
            else if (dims.size() == 3)
            {
                wv_all_no_3d_products.insert(var.first);
            }
            original_l2_products.insert(var.first);
        }
    }
    {
        const auto grp = l2_file.getGroup("sensor_band_parameters");
        const auto vars = grp.getVars();
        for (const auto &var : vars)
        {
            const auto dims = var.second.getDims();
            size_t size = 1;
            for (const auto &dim : dims)
            {
                size *= dim.getSize();
            }
            // std::cout << "\n";
            if (var.first == "wavelength_3d")
            {
                wavelengths_3d = std::vector<int>(size, 0);
                var.second.getVar(wavelengths_3d.data());
            }
            if (var.first == "wavelength")
            {
                wavelengths_all = std::vector<int>(size, 0);
                var.second.getVar(wavelengths_all.data());
            }
        }
    }
    {
        // boost::erase_all(products_requested, " ");
        boost::algorithm::split(products_requested_separated, products_requested, boost::is_any_of(", "), boost::token_compress_on); // products_requested
        products_requested_separated = create_min_max_values(products_requested_separated);
        expand_l2_requested(products_requested_separated);
        {

            for (const auto &prod : products_requested_separated)
            {
                {
                    // std::cout << "A requested product " << prod << std::endl;
                }
                if (original_l2_products.count(prod))
                {
                    l3d_to_l2_prod_map[prod] = prod;
                    l3d_indexes[prod] = {0, 1};
                    continue;
                }
                {
                    std::vector<std::string> parsed_product;
                    boost::algorithm::split(parsed_product, prod, boost::is_any_of("_"));
                    if (parsed_product.size() != 2)
                    {
                        std::cerr << "Wrong product label : " << prod << std::endl;
                        for (size_t i = 0; i < parsed_product.size(); i++)
                        {
                            std::cout << i << " " << parsed_product.at(i) << std::endl;
                        }
                        exit(EXIT_FAILURE);
                    }
                    if (original_l2_products.count(parsed_product.at(0)))
                    {
                        const auto wv_str = parsed_product.at(1);
                        try
                        {
                            const int wv_num = std::stoi(wv_str);
                            int index = -1;
                            if (wv_3d_products.count(parsed_product.at(0)))
                            {
                                auto it = std::find(wavelengths_3d.begin(), wavelengths_3d.end(), wv_num);
                                if (it == wavelengths_3d.end())
                                {
                                    std::cerr << "Error : WAVELENGTH NOT FOUND: product " << prod << ", wv (nm) =  " << wv_num << std::endl;
                                    exit(EXIT_FAILURE);
                                }
                                index = it - wavelengths_3d.begin();
                                wv_requested_3d.push_back(wv_num);
                            }
                            else if (wv_all_no_3d_products.count(parsed_product.at(0)))
                            {
                                auto it = std::find(wavelengths_all.begin(), wavelengths_all.end(), wv_num);
                                if (it == wavelengths_3d.end())
                                {
                                    std::cerr << "Error : WAVELENGTH NOT FOUND: product " << prod << ", wv (nm) =  " << wv_num << std::endl;
                                    exit(EXIT_FAILURE);
                                }
                                index = it - wavelengths_all.begin();
                                wv_requested_all_other_no_3d.push_back(wv_num);
                            }

                            l3d_indexes[prod] = {index, index + 1};
                            l3d_to_l2_prod_map[prod] = parsed_product.at(0);
                        }
                        catch (const std::exception &e)
                        {
                            std::cerr << e.what() << '\n';
                            exit(EXIT_FAILURE);
                        }
                    }
                    else
                    {
                        std::cerr << "Product not found : " << prod << std::endl;
                        exit(EXIT_FAILURE);
                    }
                }
            }
            {
                for (const auto &wv : wv_requested_all_other_no_3d)
                {
                    requested_wv.insert(wv);
                }
                for (const auto &wv : wv_requested_3d)
                {
                    requested_wv.insert(wv);
                }
            }
        }
    }
    l2_file.close();
}

int32_t get_l2prod_index(const l2_prod &l2, const char *prodname)
{
    int32_t index;
    for (index = 0; index < l2.nprod; index++)
        if (strcmp(prodname, l2.prodname[index]) == 0)
            break;
    if (index == l2.nprod)
        index = -1;
    return index;
}

void set_prodname_3d_to_l2bin(const instr &input, const std::vector<std::string> &prodparam, l2_prod &l2_str, std::vector<std::string> &l2_prodname, std::vector<std::string> &l3_prodname, std::vector<int32_t> &thirdDimId, std::vector<float> &min_value, std::vector<float> &max_value)
{

    for (size_t iprod = 0; iprod < prodparam.size(); iprod++)
    {

        const auto prod_name = prodparam[iprod];
        int32_t l2_iprod = get_l2prod_index(l2_str, prod_name.c_str());
        if (l2_iprod == -1)
        {
            printf("-E- Product: %s was not found in the L2 file: %s\n", prod_name.c_str(), l2_str.filename);
            exit(EXIT_FAILURE);
        }
        int32_t tmpThirdDim = l2_str.thirdDim[l2_iprod];
        if (tmpThirdDim == 1)
        {
            l2_prodname.push_back(prod_name);
            l3_prodname.push_back(prod_name);
            thirdDimId.push_back(0);
            min_value.push_back(min_max_values.at(prod_name).first);
            max_value.push_back(min_max_values.at(prod_name).second);
        }
        else
        {
            std::cerr << "All products must be 2D. Error for " << prod_name << std::endl;
            exit(EXIT_FAILURE);
        }
    } /* iprod loop */
}

bool set_l2_flags_use(const std::string &flagsuse)
{
    // std::cout << "Flags to be used = " << flagsuse << " " << flagsuse.length() << std::endl;
    if (flagsuse.length() > 0)
        use_l2_flags = true;
    return use_l2_flags;
}
extern "C"
{
    void l3_l2_conversion(char **inp3, char **out2)
    {
        if (ini_3d)
        {
            *out2 = &l3d_to_l2_prod_map.at(*inp3)[0];
        }
        else
        {
            *out2 = *inp3;
        }
    }

    void l3_l2_index(char **inp3, int *start, int *count)
    {
        if (ini_3d)
        {
            *start = l3d_indexes.at(*inp3).first;
            *count = l3d_indexes.at(*inp3).second - l3d_indexes.at(*inp3).first;
        }
    }
    int get_set_flag()
    {
        return ini_3d ? 1 : 0;
    }
    int get_l2_flag_use()
    {
        auto ans = use_l2_flags ? 1 : 0;
        if (!ini_3d)
            return true;
        return ans;
    }
}
