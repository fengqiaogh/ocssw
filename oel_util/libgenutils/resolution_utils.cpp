#include <genutils.h>

#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/case_conv.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string.hpp>

bool is_digits(const std::string &str)
{
    return std::all_of(str.begin(), str.end(), ::isdigit);
}

double string2resolution(std::string resolutionStr) {
    double resolution = BAD_FLT;
    boost::trim(resolutionStr);
    boost::to_lower(resolutionStr);
    resolutionStr.erase(remove_if(resolutionStr.begin(), resolutionStr.end(), ::isspace), resolutionStr.end());

    if (resolutionStr.compare("90km") == 0)
        resolution = EARTH_CIRCUMFERENCE / 432.0;
    else if (resolutionStr.compare("36km") == 0)
        resolution = EARTH_CIRCUMFERENCE / 1080.0;
    else if (resolutionStr.compare("18km") == 0)
        resolution = EARTH_CIRCUMFERENCE / 2160.0;
    else if (resolutionStr.compare("9km") == 0)
        resolution = EARTH_CIRCUMFERENCE / 4320.0;
    else if (resolutionStr.compare("4km") == 0)
        resolution = EARTH_CIRCUMFERENCE / 8640.0;
    else if (resolutionStr.compare("2km") == 0)
        resolution = EARTH_CIRCUMFERENCE / 17280.0;
    else if (resolutionStr.compare("1km") == 0)
        resolution = EARTH_CIRCUMFERENCE / 34560.0;
    else if (resolutionStr.compare("hkm") == 0)
        resolution = EARTH_CIRCUMFERENCE / 69120.0;
    else if (resolutionStr.compare("qkm") == 0)
        resolution = EARTH_CIRCUMFERENCE / 138240.0;
    else if (resolutionStr.compare("hqkm") == 0)
        resolution = EARTH_CIRCUMFERENCE / 276480.0;
    else if (resolutionStr.compare("hhkm") == 0)
        resolution = EARTH_CIRCUMFERENCE / 552960.0;
    else if (resolutionStr.compare("smi") == 0)
        resolution = EARTH_CIRCUMFERENCE / 4096.0;
    else if (resolutionStr.compare("smi4") == 0)
        resolution = EARTH_CIRCUMFERENCE / 8192.0;
    else if (resolutionStr.compare("land") == 0)
        resolution = EARTH_CIRCUMFERENCE / 8640.0;
    else if (resolutionStr.compare("thirddeg") == 0)
        resolution = EARTH_CIRCUMFERENCE / 1080.0;
    else if (boost::ends_with(resolutionStr, "km")) {
        std::string val = resolutionStr.substr(0, resolutionStr.length() - 2);
        resolution = atof(val.c_str()) * 1000.0;
    }  else if (boost::ends_with(resolutionStr, "m")) {
        std::string val = resolutionStr.substr(0, resolutionStr.length() - 1);
        resolution = atof(val.c_str());
    } else if (boost::ends_with(resolutionStr, "deg")) {
        std::string val = resolutionStr.substr(0, resolutionStr.length() - 3);
        resolution = atof(val.c_str()) / 360.0 * EARTH_CIRCUMFERENCE;
    } else if (is_digits(resolutionStr)) {
        resolution = atof(resolutionStr.c_str());
    } else {
        resolution = BAD_FLT;
    }
    return resolution;
}

std::string resolution2string(double resolution) {
    std::string resolutionStr = "Unknown";
    if (resolution > 0) {
        std::string geoUnits = " m";
        if (resolution > 1000.0) {
            resolution = resolution / 1000.0;
            geoUnits = " km";
        } 
        resolutionStr = std::to_string(resolution) + geoUnits;
    }
    return resolutionStr;
}

extern "C" double str2resolution(char const * resolutionStr) {
    std::string str(resolutionStr);
    return string2resolution(str);
}

extern "C" const char* resolution2str(double resolution) {
    static std::string stringRes = resolution2string(resolution);
    return stringRes.c_str();
}

extern "C" double resolution2degrees(double resolution) {
    return resolution / EARTH_CIRCUMFERENCE * 360.0;
}

extern "C" double degrees2resolution(double degrees) {
    return degrees / 360.0 * EARTH_CIRCUMFERENCE;
}

