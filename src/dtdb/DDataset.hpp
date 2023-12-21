/***************************************************************************
 * NAME: DDataset.hpp
 *
 * DESCRIPTION: Template header file defining data structure, and  basic 
 * read and write functions
 *
 * Created on: July, 2020
 *     Author: Sam Anderson
 
 ***************************************************************************/

#ifndef INCLUDE_DDATASET_H_
#define INCLUDE_DDATASET_H_

#include <string>
#include <vector>
#include <netcdf>
#include <boost/array.hpp>
#include <boost/multi_array.hpp>
#include <boost/lexical_cast.hpp>

using namespace std;
using namespace netCDF;

constexpr float DFILL_FLOAT = -32767.0;
constexpr float DFILL_TEST = -999;
constexpr short DFILL_SHORT = -32767;
constexpr unsigned char  DFILL_UBYTE = 255;
constexpr short MAX_SHORT = 32765;

enum class dtype
{
  BYTE     = NC_BYTE,    //!< signed 1 byte integer
  CHAR     = NC_CHAR,    //!< ISO/ASCII character
  SHORT    = NC_SHORT,   //!< signed 2 byte integer
  INT      = NC_INT,     //!< signed 4 byte integer
  FLOAT    = NC_FLOAT,   //!< single precision floating point number
  DOUBLE   = NC_DOUBLE,  //!< double precision floating point number
  UBYTE    = NC_UBYTE,   //!< unsigned 1 byte int
  USHORT   = NC_USHORT,  //!< unsigned 2-byte int
  UINT     = NC_UINT,    //!< unsigned 4-byte int
  INT64    = NC_INT64,   //!< signed 8-byte int
  UINT64   = NC_UINT64,  //!< unsigned 8-byte int
  STRING   = NC_STRING   //!< string
};

template<typename U, size_t rank>
boost::array<U, rank> stl2boost(const vector<U> &vec)
{
	assert(vec.size() == rank);
	boost::array<U, rank> result;
	std::copy(vec.begin(), vec.end(), result.begin());
	return result;
};

class ddata
{
public:
	size_t rank;
	dtype type;
	vector<size_t> start;
	vector<size_t> count;
	bool bshort;
	const void *ptr;

	ddata(){
		rank = 1;
		type = dtype::DOUBLE;
		ptr = nullptr;
		bshort = false;
	};
	
	ddata(dtype itype,
		   vector<size_t> istart,
	      vector<size_t> icount,
		   void *iptr)
	{
		rank = icount.size();
		type = itype;
		start = istart;
		count = icount;
		ptr = iptr;
		bshort = false;
	};
	
	virtual ~ddata(){
		start.clear();
		count.clear();
	};
	
	void setShortFormat(bool bshrt) {bshort = bshrt;}
	void setType(dtype val) {type = val;}
	
	void putAtt(NcFile* ncfile, string name){
        switch(type) {
        case dtype::FLOAT:
            ncfile->putAtt(name, ncFloat, count[0], static_cast<const float*> (ptr));
            break;
        case dtype::INT:
        	ncfile->putAtt(name, ncInt, count[0], static_cast<const int*> (ptr));
            break;
        case dtype::SHORT:
        	ncfile->putAtt(name, ncShort, count[0], static_cast<const short*> (ptr));
            break;
        case dtype::UBYTE:
        	ncfile->putAtt(name, ncUbyte, count[0], static_cast<const unsigned char*> (ptr));
            break;
        case dtype::BYTE:
        	ncfile->putAtt(name, ncByte, count[0], static_cast<const char*> (ptr));
            break;
        case dtype::CHAR:
        	ncfile->putAtt(name, ncChar, count[0], static_cast<const char*> (ptr));
            break;
        case dtype::DOUBLE:
            ncfile->putAtt(name, ncDouble, count[0], static_cast<const double*> (ptr));
            break;
        case dtype::USHORT:
        	ncfile->putAtt(name, ncUshort, count[0], static_cast<const unsigned short*> (ptr));
            break;
        case dtype::UINT:
        	ncfile->putAtt(name, ncUint, count[0], static_cast<const unsigned int*> (ptr));
            break;
        case dtype::INT64:
        	ncfile->putAtt(name, ncInt64, count[0], static_cast<const long long*> (ptr));
            break;
        case dtype::UINT64:
        	ncfile->putAtt(name, ncUint64, count[0], static_cast<const unsigned long long*> (ptr));
            break;
        case dtype::STRING:
            ncfile->putAtt(name, static_cast<const char*> (ptr));
            break;
        }
    }
	
	void putAtt(NcVar& var, string name, string sval){
        switch(type) {
        case dtype::FLOAT:
            var.putAtt(name, ncFloat, boost::lexical_cast<float> (sval));
            break;
        case dtype::INT:
        	var.putAtt(name, ncInt, boost::lexical_cast<int> (sval));
            break;
        case dtype::SHORT:
        	var.putAtt(name, ncShort, boost::lexical_cast<short> (sval));
            break;
        case dtype::UBYTE:
        	var.putAtt(name, ncUbyte, boost::lexical_cast<unsigned int> (sval));
            break;
        case dtype::BYTE:
        	var.putAtt(name, ncByte, boost::lexical_cast<char> (sval));
            break;
        case dtype::CHAR:
        	var.putAtt(name, ncChar, boost::lexical_cast<char> (sval));
            break;
        case dtype::DOUBLE:
            var.putAtt(name, ncDouble, boost::lexical_cast<double> (sval));
            break;
        case dtype::USHORT:
        	var.putAtt(name, ncUshort, boost::lexical_cast<unsigned short> (sval));
            break;
        case dtype::UINT:
        	var.putAtt(name, ncUint, boost::lexical_cast<unsigned int> (sval));
            break;
        case dtype::INT64:
        	var.putAtt(name, ncInt64, boost::lexical_cast<long long> (sval));
            break;
        case dtype::UINT64:
        	var.putAtt(name, ncUint64, boost::lexical_cast<unsigned long long> (sval));
            break;
        case dtype::STRING:
            var.putAtt(name, sval);
            break;
        }
	}
	
	void setFill(NcVar& var, string sval){
		const string& str = sval;
        switch(type) {
        case dtype::FLOAT:
        	if(bshort) {
             	var.setFill(true, DFILL_SHORT);       		
        	} else {
	            var.setFill(true, DFILL_FLOAT);
         }
            break;
        case dtype::INT:
        	var.setFill(true, boost::lexical_cast<int> (str));
            break;
        case dtype::SHORT:
        	var.setFill(true, boost::lexical_cast<short> (str));
            break;
        case dtype::UBYTE:
        	var.setFill(true, boost::lexical_cast<unsigned int> (str));
            break;
        case dtype::BYTE:
        	var.setFill(true, boost::lexical_cast<char> (str));
            break;
        case dtype::CHAR:
        	var.setFill(true, boost::lexical_cast<char> (str));
            break;
        case dtype::DOUBLE:
        	var.setFill(true, boost::lexical_cast<double> (str));
            break;
        case dtype::USHORT:
        	var.setFill(true, boost::lexical_cast<unsigned short> (str));
            break;
        case dtype::UINT:
        	var.setFill(true, boost::lexical_cast<unsigned int> (str));
            break;
        case dtype::INT64:
        	var.setFill(true, boost::lexical_cast<long long> (str));
            break;
        case dtype::UINT64:
        	var.setFill(true, boost::lexical_cast<unsigned long long> (str));
            break;
        case dtype::STRING:
        	var.setFill(true, str);
            break;
        }
	}
	
	void putVar(NcVar& var){
        switch(type) {
        case dtype::FLOAT:
       		if(bshort && count.size()==2) {
					float scale, offset;
					std::map<std::string,NcVarAtt> matts = var.getAtts();
					if (matts.find("scale_factor") != matts.end()) {
						matts["add_offset"].getValues(&offset);
						matts["scale_factor"].getValues(&scale);
					}
					boost::multi_array<short, 2> spts;
					boost::array<size_t,2> bdimv = stl2boost<size_t,2>(count);
					spts.resize(bdimv);
					for (size_t i=0; i<spts.num_elements(); i++) {
						if (((float*)ptr)[i] <= DFILL_TEST) {
							((short*)spts.origin())[i] = DFILL_SHORT;
						} else {
						   float ftemp = (((float*)ptr)[i] - offset)/scale;
							((short*)spts.origin())[i] = 
								(ftemp > (float) MAX_SHORT) ? MAX_SHORT : (short) ftemp;
						}
					}
					var.putVar(start, count, spts.origin());
        		} else {
	            var.putVar(start, count, static_cast<const float*> (ptr));
            }
            break;
        case dtype::INT:
        		var.putVar(start, count, static_cast<const int*> (ptr));
            break;
        case dtype::SHORT:
        	var.putVar(start, count, static_cast<const short*> (ptr));
            break;
        case dtype::UBYTE:
        	var.putVar(start, count, static_cast<const unsigned char*> (ptr));
            break;
        case dtype::BYTE:
        	var.putVar(start, count, static_cast<const char*> (ptr));
            break;
        case dtype::CHAR:
        	var.putVar(start, count, static_cast<const char*> (ptr));
            break;
        case dtype::DOUBLE:
            var.putVar(start, count, static_cast<const double*> (ptr));
            break;
        case dtype::USHORT:
        	var.putVar(start, count, static_cast<const unsigned short*> (ptr));
            break;
        case dtype::UINT:
        	var.putVar(start, count, static_cast<const unsigned int*> (ptr));
            break;
        case dtype::INT64:
        	var.putVar(start, count, static_cast<const long long*> (ptr));
            break;
        case dtype::UINT64:
        	var.putVar(start, count, static_cast<const unsigned long long*> (ptr));
            break;
        case dtype::STRING:
            var.putVar(start, count, static_cast<const char*> (ptr));
            break;
        }
	}		
};

class ddstr: public ddata
{
public:
	string str;

	ddstr(){};
	~ddstr(){};

	ddstr(const string istr) {
		str = istr;
		type = dtype::STRING;
		start.push_back(0);
		count.push_back(1);
		ptr = str.c_str();
	}
};

template<typename T>
class ddval: public ddata
{
public:
	T val;

	ddval(){};
	~ddval(){};

	ddval<T>(dtype itype, const T ival) {
		val = ival;
		type = itype;
		start.push_back(0);
		count.push_back(1);
		ptr = &val;
	}
};

template<typename T, size_t ndim>
class ddma: public ddata
{
public:
	boost::multi_array<T, ndim> pts;

	ddma(){};
	~ddma(){};

	ddma<T,ndim>(dtype itype,
				 vector<size_t> istart,
				 vector<size_t> icount) {
		rank = ndim;
		type = itype;
		start = istart;
		count = icount;
		boost::array<size_t, ndim> bdimv = stl2boost<size_t,ndim>(count);
		pts.resize(bdimv);
		ptr = pts.origin();
	}
	
	ddma<T,ndim>(NcGroup& ncgrp,
			string name,
			dtype itype,
			vector<size_t> istart,
			vector<size_t> icount) 
	{
		NcVar var = ncgrp.getVar(name);
		if (!var.isNull()) {
			type = itype;
			rank = var.getDimCount();
			vector<NcDim> dims = var.getDims();
			size_t beg = (rank == ndim+1) ? 1 : 0;
			for (size_t i=beg; i<rank; i++){
				start.push_back(istart[i]);
				count.push_back(icount[i]);
			}
			boost::array<size_t, ndim> bdimv = stl2boost<size_t,ndim>(count);
			pts.resize(bdimv);
			ptr = pts.origin();
			var.getVar(istart, icount, (T*)ptr);
			
			T factor = 1.0;
			T offset = 0.0;
			T fillvalue = 0;
			std::map<std::string,NcVarAtt> matts = var.getAtts();
			if (matts.find("scale_factor") != matts.end()) {
				matts["scale_factor"].getValues(&factor);
				matts["add_offset"].getValues(&offset);
			}
	    	matts["_FillValue"].getValues(&fillvalue);
			for (size_t i=0; i<pts.num_elements(); i++) {
				if (((T*)ptr)[i] == fillvalue || ((T*)ptr)[i] == fillvalue-2) {
					((T*)ptr)[i] = (T)DFILL_FLOAT;
				} else {
					((T*)ptr)[i] = ((T*)ptr)[i]*factor + offset;
				}
			}
		} else {
			ptr = nullptr;
		}
	}
};

#endif
