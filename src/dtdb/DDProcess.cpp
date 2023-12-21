/*******************************************************************************
 *
 * NAME: DDProcess.cpp
 *
 * DESCRIPTION: Process granule using selected algorithm.
 *
 *  Created on: Aug 25, 2020
 *      Author: Sam Anderson
 *
 *******************************************************************************/

#include <boost/lexical_cast.hpp>
#include <DDataset.hpp>
#include <DDProcess.h>
#include <DDSensor.h>
#include <DDOptions.h>
#include <DDAncillary.h>
#include "darktarget/DtAlgorithm.h"
#include "deepblue/DbAlgorithm.h"

using namespace std;
using namespace pugi;

const vector<string> DDProcess::out_groups = { "geolocation", "ancillary",
		"observations", "geophysical_data", "statistics" };

const map<string,string> DDProcess::odps2dtdb = {
		{"cloud_mask", "cloud_mask"},
		{"land_water", "land_water"},
		{"relaz", "relative_azimuth"},
		{"solz", "solar_zenith"},
		{"senz", "sensor_zenith"},
		{"scattang", "scattang"},
		{"elev", "elevation"},
		{"windspeed", "windspeed"},
		{"windangle", "windangle"},
		{"water_vapor", "water_vapor"},
		{"pressure", "pressure"},
};

const map<string,dtype> DDProcess::input_names = {
		{"cloud_mask", dtype::UBYTE},
		{"land_water", dtype::UBYTE},
		{"relaz", dtype::FLOAT},
		{"solz", dtype::FLOAT},
		{"senz", dtype::FLOAT},
		{"elev", dtype::FLOAT},
		{"windspeed", dtype::FLOAT},
		{"windangle", dtype::FLOAT},
		{"water_vapor", dtype::FLOAT},
		{"pressure", dtype::FLOAT},
};

const map<string,dtype> DDProcess::dtdb_names = {
		{"quality", dtype::SHORT},
		{"aerosol_type", dtype::SHORT},
		{"l2_flags", dtype::UINT},
		{"fmf_550", dtype::FLOAT},
		{"angstrom", dtype::FLOAT},
		{"scattang", dtype::FLOAT},
};

const map<string,rhot_band> DDProcess::rhot_band_names = {
	    {"rhot_410", rhot_band::W410},           // m01
	    {"rhot_445", rhot_band::W445},           // m02
	    {"rhot_490", rhot_band::W490},           // m03
	    {"rhot_550", rhot_band::W550},           // m04
	    {"rhot_670", rhot_band::W670},           // m05
	    {"rhot_865", rhot_band::W865},           // m07
	    {"rhot_1240", rhot_band::W1240},         // m08
	    {"rhot_1380", rhot_band::W1380},         // m09
	    {"rhot_1610", rhot_band::W1610},         // m10
	    {"rhot_2250", rhot_band::W2250},         // m11
};

const map<string,aot_band> DDProcess::aot_band_names = {
	    {"aot_410", aot_band::W410},           // m01
	    {"aot_490", aot_band::W490},           // m03
	    {"aot_550", aot_band::W550},           // m04
	    {"aot_670", aot_band::W670},           // m05
	    {"aot_865", aot_band::W865},           // m07
	    {"aot_1240", aot_band::W1240},         // m08
	    {"aot_1610", aot_band::W1610},         // m10
	    {"aot_2250", aot_band::W2250},         // m11
};

const map<string,srf_band> DDProcess::srf_band_names = {
	    {"surface_410", srf_band::W410},           // m01
	    {"surface_490", srf_band::W490},           // m03
	    {"surface_670", srf_band::W670},           // m05
	    {"surface_2250", srf_band::W2250},         // m11
};

/**************************************************************************
 * NAME: DDProcess()
 *
 * DESCRIPTION: Class Constructor
 *
 *************************************************************************/
DDProcess::DDProcess()
{
    lines_ = 0.0;
    pixels_ = 0.0;
    alg_ = ALGO::NOTHING;
    instrument_ = SENSOR::NOTHING;
    psensor_ = nullptr;
    pancillary_ = nullptr;
    pl_ = nullptr;
    po_ = nullptr;
    bday_ = true;
    lprw_ = 1;
    bmaskglint_=true;
    bmaskcloud_=true;
    bmasksolz_=true;
    bmasksenz_=true;
    bgascorrect_=true;
}

/**************************************************************************
 * NAME: ~DDProcess()
 *
 * DESCRIPTION: Class Destructor
 *
 *************************************************************************/
DDProcess::~DDProcess()
{
    if (psensor_ != 0) {
        delete psensor_;
        psensor_ = 0L;
    }
    if (pancillary_ != 0) {
        delete pancillary_;
        pancillary_ = 0L;
    }
}

/**************************************************************************
 * NAME: initialize()
 *
 * DESCRIPTION: Initializes data and object classes for granule
 *
 *************************************************************************/
int DDProcess::initialize()
{
	string salg = get_option(ALGORITHM);
	if (salg == "darktarget") {
		std::cerr << "DDProcess:: Algorithm Dark Target" << std::endl;
		alg_ = ALGO::DARKTARGET;
		algStr_ = "dt";
        pl_ = new DtAlgLand();
        po_ = new DtAlgOcean();
	} else if (salg == "deepblue") {
		std::cerr << "DDProcess:: Algorithm Deep Blue" << std::endl;
		alg_ = ALGO::DEEPBLUE;
		algStr_ = "db";
        pl_ = new DbAlgLand();
	    po_ = new DbAlgOcean();
	} else {
		std::cerr << "DDProcess:: Invalid algorithm, exiting ..." << std::endl;
		return DTDB_FAIL;
	}
    if (query_granule() != DTDB_SUCCESS) {
		cerr << "DDProcess:: Query granule failure" << endl;
    	return DTDB_FAIL;
    } else if (!bday_) {
		cerr << "DDProcess:: Night granule detected" << endl;
    	return DTDB_FAIL;
    }
    switch (instrument_) {
    case SENSOR::POCI:
        psensor_ = new POCI();
        break;
    case SENSOR::VIIRS:
        psensor_ = new VIIRS();
        break;
    case SENSOR::NOTHING:
    default:
		cerr << "DDProcess:: Invalid sensor" << endl;
        return DTDB_FAIL;
        break;
    }
    pancillary_ = new DDAncillary();

	lprw_ = (size_t) get_option_int(INPUT_LINES_PER_RW);
	lprw_ = (lprw_<3) ? 3 : lprw_;

	vector<string> products = pl_->get_products();
	out_products.insert(out_products.end(), products.begin(), products.end());
	products = po_->get_products();
	out_products.insert(out_products.end(), products.begin(), products.end());

	return DTDB_SUCCESS;
}

/**************************************************************************
 * NAME: query_granule()
 *
 * DESCRIPTION: Read granule data essential for initialization
 *
 *************************************************************************/
int DDProcess::query_granule()
{
	string filepath = get_option(INPUT_L1B);

    if (!filepath.empty()) {
        NcFile* nc_input;
		try {
			nc_input = new NcFile(filepath, NcFile::read );
		}
		catch( NcException& e) {
			e.what();
	        cout << "DDProcess:: Failure opening L1B file" << endl;
			return DTDB_FAIL;
		}
		string str;
		multimap <string, NcGroupAtt> attributes = nc_input->getAtts(NcGroup::Current);
		multimap <string, NcGroupAtt>::iterator it;
        if ((it = attributes.find("instrument")) != attributes.end()) {
            (it->second).getValues(str);
        }
        if (str.find("OCI") != string::npos) {
            instrument_ =  SENSOR::POCI;
        } else if (str.find("VIIRS") != string::npos) {
            instrument_ = SENSOR::VIIRS;
        } else {
	        cout << "DDProcess:: Instrument not supported" << endl;
            return DTDB_FAIL;
        }
	    switch (instrument_) {
	    case SENSOR::VIIRS:
			if ((it = attributes.find("DayNightFlag")) != attributes.end()) {
				(it->second).getValues(str);
			} else if((it = attributes.find("day_night_flag")) != attributes.end()) {
					(it->second).getValues(str);
			}
			if (str != "Night") {
				bday_ = true;
			} else {
		        bday_ = false;
			}
	        lines_ = nc_input->getDim("number_of_lines").getSize();
	        pixels_ = nc_input->getDim("number_of_pixels").getSize();
	        title_ = "VIIRS Level-2 Data";
	        break;
	    case SENSOR::POCI:
        	bday_ = true;
	    	lines_ = nc_input->getDim("number_of_scans").getSize();
	    	pixels_ = nc_input->getDim("ccd_pixels").getSize();
	        title_ = "OCI Level-2 Data";
	        break;
	    default:
	        cout << "DDProcess:: Failure reading dimensions" << endl;
	        return DTDB_FAIL;
	    }
		delete nc_input;
	}  else {
        cout << "Invalid L1B filename ..." << endl;
        return DTDB_FAIL;
    }
	return DTDB_SUCCESS;
}

/**************************************************************************
 * NAME: process()
 *
 * DESCRIPTION: Processes data and produces output values for granule
 *
 *************************************************************************/
int DDProcess::process()
{
	map<string,ddata*> imap = psensor_->create( {0, 0}, {lines_, pixels_} );
    if (static_cast<ddval<int>*>(imap["status"])->val != DTDB_SUCCESS) {
        cerr << "DDProcess:: Create sensor failure at line " << endl;
        return DTDB_FAIL;
    }

    ddstr* dstr = new ddstr(title_);
    imap.insert({"title", dstr});
	ddma<unsigned char,2>* plw = static_cast<ddma<unsigned char,2>*>(imap["land_water"]);
	string parfilepath = get_option(INPUT_PAR);
    dstr = new ddstr(parfilepath);
    imap.insert({"config_file", dstr});
    ddval<unsigned char>* bval = new ddval<unsigned char>(dtype::UBYTE, get_bool(BOOL_GAS_CORRECTION));
    imap.insert({"bgascorrect", bval});
    bgascorrect_ = bval->val;
    bval = new ddval<unsigned char>(dtype::UBYTE, get_bool(BOOL_MASK_GLINT));
    imap.insert({"bmaskglint", bval});
    bmaskglint_ = bval->val;
    bval = new ddval<unsigned char>(dtype::UBYTE, get_bool(BOOL_MASK_CLOUD));
    imap.insert({"bmaskcloud", bval});
    bmaskcloud_ = bval->val;
    bval = new ddval<unsigned char>(dtype::UBYTE, get_bool(BOOL_MASK_SOLZ));
    imap.insert({"bmasksolz", bval});
    bmasksolz_ = bval->val;
    bval = new ddval<unsigned char>(dtype::UBYTE, get_bool(BOOL_MASK_SENZ));
    imap.insert({"bmasksenz", bval});
    bmasksenz_ = bval->val;
    ddval<float>* fval = new ddval<float>(dtype::FLOAT, get_option_float(INPUT_THRESH_SOLZ));
    imap.insert({"threshsolz", fval});
    fval = new ddval<float>(dtype::FLOAT, get_option_float(INPUT_THRESH_SENZ));
    imap.insert({"threshsenz", fval});
    fval = new ddval<float>(dtype::FLOAT, get_option_float(INPUT_THRESH_GLINT));
    imap.insert({"threshglint", fval});

    int status = pancillary_->initialize( (static_cast<ddval<int>*>(imap["start_hour"]))->val,
    		                          (static_cast<ddval<int>*>(imap["start_minute"]))->val);
    if (status != DTDB_SUCCESS) {
        cerr << "DDProcess:: Ancillary initialization failure" << endl;
        return DTDB_FAIL;
    }
    cerr << "DDProcess:: Initializing Land Process" << endl;
    if (pl_->initialize(imap) != DTDB_SUCCESS) {
        cerr << "DDProcess:: Land initialization failure" << endl;
        return DTDB_FAIL;
    }
    cerr << "DDProcess:: Initializing Ocean Process" << endl;
    if (po_->initialize(imap) != DTDB_SUCCESS) {
        cerr << "DDProcess:: Ocean initialization failure" << endl;
        return DTDB_FAIL;
    }

    cerr << "DDProcess:: Processing Granule" << endl;
    if (create_nc4(imap) != DTDB_SUCCESS) {
        cerr << "DDProcess:: Output file creation failure" << endl;
        return DTDB_FAIL;
    }

	for (auto &it : imap) {
		string name = (string) it.first;
		if (name != "land_water") {
			delete imap[name];
		}
	}

 	map<string, ddma<float, 2>*> flmap;
 	map<string, ddma<unsigned int, 2>*> uimap;
 	map<string, ddma<short, 2>*> shmap;
 	map<string, ddma<unsigned char, 2>*> ubmap;
 	map<string, ddata*> smap;
 	map<string, ddata*> amap;
 	map<string, ddata*> pmap;

	for (auto &it : dtdb_names) {
 		string name = (string) it.first;
 		if ( it.second == dtype::FLOAT ) {
 			ddma<float, 2>* ddfl =
 				new ddma<float, 2>(dtype::FLOAT, {0,0}, {lprw_, pixels_});
 			flmap.insert({ name, ddfl });
 		} else if ( it.second == dtype::UINT) {
 			ddma<unsigned int, 2>* ddin =
 				new ddma<unsigned int, 2>(dtype::UINT, {0,0}, {lprw_, pixels_});
 			uimap.insert({ name, ddin });
 		} else if ( it.second == dtype::SHORT ) {
 			ddma<short, 2>* ddsh =
 				new ddma<short, 2>(dtype::SHORT, {0,0}, {lprw_, pixels_});
 			shmap.insert({ name, ddsh });
 		} else if ( it.second == dtype::UBYTE ) {
			ddma<unsigned char, 2>* ddub =
				new ddma<unsigned char, 2>(dtype::UBYTE, {0,0}, {lprw_, pixels_});
			ubmap.insert({ name, ddub });
		}
 	}
	ddma<unsigned char, 2>* ddcm =
		new ddma<unsigned char, 2>(dtype::UBYTE, {0,0}, {lprw_, pixels_});
	ubmap.insert({ "cloud_mask", ddcm });

	for (auto &it : rhot_band_names) {
 		string name = (string) it.first;
 		ddma<float, 2>* ddfl =
 				new ddma<float, 2>(dtype::FLOAT, {0,0}, {lprw_, pixels_});
		string oname = "rhot_" + name.substr(1);
 		flmap.insert({ oname, ddfl });
 	}
	for (auto &it : aot_band_names) {
 		string name = (string) it.first;
 		ddma<float, 2>* ddfl =
 				new ddma<float, 2>(dtype::FLOAT, {0,0}, {lprw_, pixels_});
 		flmap.insert({ name, ddfl });
 	}
	for (auto &it : srf_band_names) {
 		string name = (string) it.first;
 		ddma<float, 2>* ddfl =
 				new ddma<float, 2>(dtype::FLOAT, {0,0}, {lprw_, pixels_});
 		flmap.insert({ name, ddfl });
 	}
    ddma<unsigned char, 2>* ddlw =
    		new ddma<unsigned char, 2>(dtype::UBYTE, {0,0}, {lprw_, pixels_});

// for special aot_380
	ddma<float, 2>* ddfl =
			new ddma<float, 2>(dtype::FLOAT, {0,0}, {lprw_, pixels_});
	flmap.insert({ "aot_380", ddfl });

 	size_t cline = 0;
 	while ( cline < lines_ ) {
		imap.clear();
		size_t sy = cline;
		size_t cy = lprw_;
		size_t left = lines_-cline;
		if (left < lprw_) {
			cy = left;
		}
		for (auto &it : smap) {
			delete smap[(string) it.first];
		}
		smap.clear();
		smap = psensor_->read( {sy, 0}, {cy, pixels_} );
	    if (static_cast<ddval<int>*>(smap["status"])->val != DTDB_SUCCESS) {
	        cerr << "DDProcess:: Read sensor failure at line " << sy << endl;
	    }
	    delete smap["status"];
	    smap.erase("status");
		imap.insert(smap.begin(), smap.end());
		for (auto &it : amap) {
			delete amap[(string) it.first];
		}
		amap.clear();
		amap = pancillary_->read(imap);
	    if (static_cast<ddval<int>*>(amap["status"])->val != DTDB_SUCCESS) {
	        cerr << "DDProcess:: Read ancillary failure at line " << sy << endl;
	    }
	    delete amap["status"];
	    amap.erase("status");
		imap.insert(amap.begin(), amap.end());
	    for ( size_t iy = 0; iy < cy; iy++ ) {
			for ( size_t ix = 0; ix < pixels_; ix++ ) {
				if (ix==5 && iy+sy==1714) {
//					bool testmode = true;
				}
				if (plw->pts[sy+iy][ix] == 0) {
					pmap = po_->process({iy,ix}, {cy,pixels_}, imap);
					ddlw->pts[iy][ix] = 0;
				    if (static_cast<ddval<int>*>(pmap["status"])->val != DTDB_SUCCESS) {
				        cerr << "DDProcess:: Ocean process failure at line " << iy+sy << " pixel " << ix << endl;
				    }
				} else if (plw->pts[sy+iy][ix] == 1) {
					pmap = pl_->process({iy,ix}, {cy,pixels_}, imap);
					ddlw->pts[iy][ix] = 1;
				    if (static_cast<ddval<int>*>(pmap["status"])->val != DTDB_SUCCESS) {
				        cerr << "DDProcess:: Land process failure at line " << iy+sy << " pixel " << ix << endl;
				    }
				}
				for (auto &it : flmap) {
					string name = (string) it.first;
					if (pmap.find(name) != pmap.end()) {
						flmap[name]->pts[iy][ix] =
							(static_cast<ddval<float>*>(pmap[name]))->val;
					}
				}
				for (auto &it : uimap) {
					string name = (string) it.first;
					if (pmap.find(name) != pmap.end()) {
						uimap[name]->pts[iy][ix] =
							(static_cast<ddval<unsigned int>*>(pmap[name]))->val;
					}
				}
				for (auto &it : shmap) {
					string name = (string) it.first;
					if (pmap.find(name) != pmap.end()) {
						shmap[name]->pts[iy][ix] =
							(static_cast<ddval<short>*>(pmap[name]))->val;
					}
				}
				for (auto &it : ubmap) {
					string name = (string) it.first;
					if (pmap.find(name) != pmap.end()) {
						ubmap[name]->pts[iy][ix] =
							(static_cast<ddval<short>*>(pmap[name]))->val;
					}
				}
				if (pmap.find("cloud_mask") != pmap.end()) {
					ubmap["cloud_mask"]->pts[iy][ix] =
							(static_cast<ddval<unsigned char>*>(pmap["cloud_mask"]))->val;
				}
				for (auto &it : pmap) {
					delete pmap[(string) it.first];
				}
				pmap.clear();
			}
		}
	    delete amap["cloud_mask"];
	    amap.erase("cloud_mask");
	    imap.erase("cloud_mask");
	 	for (auto &it : flmap) {
	 		string name = (string) it.first;
	 		flmap[name]->start = {cline, 0};
	 		flmap[name]->count = {cy,pixels_};
	 		imap.insert({name, flmap[name]});
	 	}
	 	for (auto &it : uimap) {
	 		string name = (string) it.first;
	 		uimap[name]->start = {cline, 0};
	 		uimap[name]->count = {cy,pixels_};
	 		imap.insert({name, uimap[name]});
	 	}
	 	for (auto &it : shmap) {
	 		string name = (string) it.first;
	 		shmap[name]->start = {cline, 0};
	 		shmap[name]->count = {cy,pixels_};
	 		imap.insert({name, shmap[name]});
	 	}
	 	for (auto &it : ubmap) {
	 		string name = (string) it.first;
	 		ubmap[name]->start = {cline, 0};
	 		ubmap[name]->count = {cy,pixels_};
	 		imap.insert({name, ubmap[name]});
	 	}
	 	ddlw->start = {cline,0};
	 	ddlw->count = {cy,pixels_};
	 	imap.insert({"land_water", ddlw});
		cerr << ".";
        status = write_nc4( imap );
    	cline += lprw_;
	}
	for (auto &it : imap) {
		delete imap[(string) it.first];
	}
	imap.clear();

	delete plw;
    delete pl_;
    delete po_;

    cerr << " DONE! " << endl;

    if (get_bool(BOOL_STATISTICS)) {
    	vector<string> dnames = {"fmf_550","angstrom","aot_410","aot_490",\
    			"aot_550","aot_670","aot_865","aot_1240","aot_1610","aot_2250"};
    	size_t window = 8;
        cerr << "DDProcess:: Generating statistics on decimated samples" << endl;
    	if (write_decimated(dnames, window) != DTDB_SUCCESS) {
            cerr << "DDProcess:: Write decimated statistics failure" << endl;
    	}
        cerr << "STATS DONE! " << endl;
    }

	return DTDB_SUCCESS;
}

/**************************************************************************
 * NAME: create_nc4()
 *
 * DESCRIPTION: Create NetCDF4 product file
 *
 *************************************************************************/
int DDProcess::create_nc4( map<string, ddata*> imap )
{
    vector<size_t> sdims, cdims;
    sdims = {0,0};
    cdims = {lines_, pixels_};
    ddata* odat = new ddata(dtype::FLOAT, sdims, cdims, nullptr);
// Create output file
    string output_filepath = get_option(OUTPUT_NC4);
    if (output_filepath.empty()) {
        cerr << "\nDDProcess:: Failure locating output file path.\n" << endl;
        return DTDB_FAIL;
    }
    NcFile* nc_output;
    try {
        nc_output = new NcFile( output_filepath, NcFile::replace );
    } catch( NcException& e) {
        e.what();
        cerr << "\nDDProcess:: Failure creating product file: " + output_filepath + "\n" << endl;
        return DTDB_FAIL;
    }
    try {
        string xfilepath = get_option(DTDB_XML);
        xml_document doc;
        xml_parse_result result = doc.load_file(xfilepath.c_str());
        if ( result.status !=  status_ok ) {
            cerr << "DDProcess:: Failure opening dtdb product configuration XML file " +
                    xfilepath << endl;
            return DTDB_FAIL;
        }
        xml_node dtdbnode = doc.child("products");
        for (xml_node attnode: dtdbnode.children()) {
            string sname = string(attnode.name());
            if (sname == "group" || sname == "metadata") continue;
            nc_output->putAtt(sname, attnode.child_value());
        }

        for (xml_node meta: dtdbnode.children("metadata"))
        {
            string mname = meta.attribute("name").value();
            if (imap.find(mname) != imap.end()) {
                ddata* dsp = imap[mname];
                dsp->putAtt(nc_output, mname);
            }
        }
        NcDim lines_dim = nc_output->addDim( "number_of_lines", lines_ );
        NcDim pixels_dim = nc_output->addDim( "pixels_per_line", pixels_ );
        NcDim all_bands_dim = nc_output->addDim( "number_of_bands", NTWL );
        dmap_.insert({lines_dim.getSize(), lines_dim});
        dmap_.insert({pixels_dim.getSize(), pixels_dim});
        dmap_.insert({all_bands_dim.getSize(), all_bands_dim});
        vector<NcDim> pdims;
        xfilepath = get_option(PRODUCT_XML);
        xml_document docx;
        result = docx.load_file(xfilepath.c_str());
        if ( result.status != status_ok ) {
            cerr << "DDProcess:: Failure opening product XML file " +
                    xfilepath << endl;
            return DTDB_FAIL;
        }
        xml_node prodnode = docx.child("products");
        for (xml_node group: dtdbnode.children("group"))  {
            string gname = group.attribute("name").value();
            if (!get_bool(gname)) continue;
            bool boutput = false;
            for(string str : out_groups) {
                if (str == gname) boutput = true;
            }
            NcGroup grp = nc_output->addGroup(gname);
            for (xml_node product: group.children("product")) {
                string pname = product.attribute("name").value();
                ddata* dsp;
                odat->setType(dtype::FLOAT);
                odat->setShortFormat(false);
                if (imap.find(pname) != imap.end()) {
                    dsp = imap[pname];
                } else if (input_names.find(pname) != input_names.end()) {
                	odat->setType(input_names.at(pname));
                	dsp = odat;
                } else if (dtdb_names.find(pname) != dtdb_names.end()) {
                	odat->setType(dtdb_names.at(pname));
                	dsp = odat;
                } else dsp = odat;
        		string oname = pname;
        		if(gname=="geophysical_data") {
            		auto rslt = find(out_products.begin(), out_products.end(), pname);
            		if (pname != "l2_flags") {
            			if (!get_bool(pname) || (rslt == out_products.end())) continue;
            			oname +=  "_" + algStr_;
        			}
        		}
                pdims.clear();
                for (size_t i=0; i<dsp->count.size(); i++) {
                    NcDim tdim = dmap_[dsp->count[i]];
                    pdims.push_back(tdim);
                }
                dtype nctype = dsp->type;
                xml_node productNode = prodnode.find_child_by_attribute("product", "name", pname.c_str());
                if (productNode) {
					if ( !strcmp(productNode.child("type").child_value(), "short") &&
					     dsp->type == dtype::FLOAT &&
						 get_bool(BOOL_SHORT_FORMAT))
					{
						nctype = dtype::SHORT;
						dsp->setShortFormat(true);
					}
                }
                NcVar var = grp.addVar(oname, NcType((nc_type)nctype), pdims);
       			vector<size_t> chunksizes = {lprw_, pixels_};
                var.setCompression(bShuffleFilter, bDeflateFilter, deflateLevel);
                var.setChunking(var.nc_CHUNKED, chunksizes);
                write_nc4_attributes(pname, prodnode, var, dsp);
                if (!boutput) dsp->putVar(var);
            }
        }
    } catch( NcException& e) {
        e.what();
        cerr << "\nDDProcess:: Failure initializing product file: " +
        		output_filepath + "\n" << endl;
        return DTDB_FAIL;
    }
	dmap_.clear();
    delete nc_output;

    return DTDB_SUCCESS;
}

/**************************************************************************
 * NAME: write_nc4()
 *
 * DESCRIPTION: Write a single line to output file.
 *
 *************************************************************************/
int  DDProcess::write_nc4( map<string, ddata*> omap )
{
    string output_filepath = get_option(OUTPUT_NC4);
    if (output_filepath.empty()) {
        cerr << "\nDDProcess:: Failure locating output file path.\n" << endl;
        return DTDB_FAIL;
    }
    NcFile* nc_output;
    try {
        nc_output = new NcFile( output_filepath, NcFile::write );
    }
    catch( NcException& e) {
        e.what();
        cerr << "\nDDProcess:: Failure writing to product file: " +
                output_filepath + "\n" << endl;
        return DTDB_FAIL;
    }

    string xfilepath = get_option(DTDB_XML);
    xml_document doc;
    xml_parse_result result = doc.load_file(xfilepath.c_str());
    if ( result.status !=  status_ok ) {
        cerr << "DDProcess::write_nc4: Failure opening XML product configuration file " +
                xfilepath << endl;
        return DTDB_FAIL;
    }
    xml_node dtdbnode = doc.child("products");
    xfilepath = get_option(PRODUCT_XML);
    xml_document docx;
    result = docx.load_file(xfilepath.c_str());
    if ( result.status != status_ok ) {
        cerr << "DDProcess::write_nc4: Failure opening product XML file " +
                xfilepath << endl;
        return DTDB_FAIL;
    }
    xml_node prodnode = docx.child("products");

    for (string gname : out_groups) {
        if (!get_bool(gname)) continue;
        try {
            xml_node group =
                 dtdbnode.find_child_by_attribute("group", "name", gname.c_str());
            NcGroup grp = nc_output->getGroup(gname);
            for (xml_node product: group.children("product")) {
                string pname = product.attribute("name").value();
        		string dname = pname;
                if (odps2dtdb.find(pname) != odps2dtdb.end()) {
                	dname = odps2dtdb.at(pname);
                }
        		string oname = pname;
        		if(gname=="geophysical_data") {
            		auto rslt = find(out_products.begin(), out_products.end(), pname);
            		if (pname != "l2_flags") {
            			if (!get_bool(pname) || (rslt == out_products.end())) continue;
            			oname +=  "_" + algStr_;
        			}
        		}
                if (omap.find(dname) != omap.end()) {
				    ddata* dsp = omap[dname];
	                xml_node pNode =
	                	prodnode.find_child_by_attribute("product", "name", pname.c_str());
	                if (pNode) {
						if ( !strcmp(pNode.child("type").child_value(), "short") &&
						     dsp->type == dtype::FLOAT &&  get_bool(BOOL_SHORT_FORMAT)) {
							dsp->setShortFormat(true);
						}
	                }
					NcVar var = grp.getVar(oname);
					if (!var.isNull()) {
						dsp->putVar(var);
					} else {
			            cerr << "\nDDProcess::write_nc4: Warning! Failed getVar for " +
			            		oname + " \n" << endl;
					}
                }
            }
        }
        catch( NcException& e) {
            e.what();
            cerr << "\nDDProcess::write_nc4 Failure writing line " +
            		output_filepath + "\n" << endl;
            return DTDB_FAIL;
        }
    }
    delete nc_output;

    return DTDB_SUCCESS;
}


/**************************************************************************
 * NAME: write_nc4_attributes()
 *
 * DESCRIPTION: Set dataset attributes from product.xml file.
 *
 *************************************************************************/
int DDProcess::write_nc4_attributes(string name, xml_node products,
		NcVar var, ddata* dsp)
{
	// set some defaults
	dsp->setFill(var, "-32767");

	if (name == "l2_flags") {
		vector<int> flag_masks = {1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024,
				2048, 4096, 8192, 16384, 32768, 65536, 131072, 262144, 524288,
				1048576, 2097152, 4194304, 8388608, 16777216, 33554432, 67108864,
				134217728, 268435456, 536870912, 1073741824};
		var.putAtt("flag_masks",ncInt, size_t(31), &flag_masks[0]);
		var.putAtt("flag_meanings","ATMFAIL LAND PRODWARN HIGLINT HILT HISATZEN COASTZ SPARE STRAYLIGHT CLDICE COCCOLITH TURBIDW HISOLZEN SPARE LOWLW CHLFAIL NAVWARN ABSAER SPARE MAXAERITER MODGLINT CHLWARN ATMWARN SPARE SEAICE NAVFAIL FILTER SPARE5 BOWTIEDEL HIPOL PRODFAIL");
	}

	// read attributes from product.xml file
    xml_node node;
    xml_node productNode = products.find_child_by_attribute("product", "name", name.c_str());
    if (!productNode) {
//		cout << "DDProcess::write_nc4_attributes(): product " << name << " not defined in product.xml" << endl;
        return DTDB_FAIL;
    }
    if ((node = productNode.child("standardName"))) {
        var.putAtt("standardName", node.child_value());
    }
    if ((node = productNode.child("units"))) {
        var.putAtt("units", node.child_value());
    } else {
		cout << "DDProcess::write_nc4_attributes(): units not defined in product.xml" << endl;
        return DTDB_FAIL;
    }
    if ((node = productNode.child("category"))) {
        var.putAtt("category", node.child_value());
    } else {
		cout << "DDProcess::write_nc4_attributes(): category not defined in product.xml" << endl;
        return DTDB_FAIL;
    }
    if ((node = productNode.child("displayScale"))) {
        var.putAtt("displayScale", node.child_value());
    }
    if ((node = productNode.child("reference"))) {
        var.putAtt("reference", node.child_value());
    }
    if ((node = productNode.child("comment"))) {
        var.putAtt("comment", node.child_value());
    }
    xml_node rangeNode = productNode.child("range");
    if (rangeNode) {
		if ((node = rangeNode.child("validMin")))
			dsp->putAtt(var, "valid_min", node.child_value());
		if ((node = rangeNode.child("validMax")))
			dsp->putAtt(var, "valid_max", node.child_value());
		if ((node = rangeNode.child("scaleFactor")))
			var.putAtt("scale_factor", ncFloat,
					boost::lexical_cast<float> (node.child_value()));
		if ((node = rangeNode.child("addOffset")))
			var.putAtt("add_offset", ncFloat,
					boost::lexical_cast<float> (node.child_value()));
    }
    xml_node algorithmNode =
    		productNode.find_child_by_attribute("algorithm", "name", algStr_.c_str());
    if (algorithmNode) {
		if ((node = algorithmNode.child("units"))) {
			var.putAtt("units", node.child_value());
		}
		if ((node = algorithmNode.child("fillValue"))) {
			dsp->setFill(var, node.child_value());
		}
		if ((node = algorithmNode.child("description"))) {
			var.putAtt("description", node.child_value());
		}
		if ((node = algorithmNode.child("reference"))) {
			var.putAtt("reference", node.child_value());
		}
		if ((node = algorithmNode.child("comment"))) {
			var.putAtt("comment", node.child_value());
		}
		rangeNode = algorithmNode.child("range");
		if (rangeNode) {
			if ((node = rangeNode.child("validMin")))
				dsp->putAtt(var, "valid_min", node.child_value());
			if ((node = rangeNode.child("validMax")))
				dsp->putAtt(var, "valid_max", node.child_value());
			if ((node = rangeNode.child("scaleFactor")))
				var.putAtt("scale_factor", ncFloat,
						boost::lexical_cast<float> (node.child_value()));
			if ((node = rangeNode.child("addOffset")))
				var.putAtt("add_offset", ncFloat,
						boost::lexical_cast<float> (node.child_value()));
		}
    }
    return DTDB_SUCCESS;
}


/**************************************************************************
 * NAME: write_decimated()
 *
 * DESCRIPTION: Compute and write decimated data sets.
 *
 *************************************************************************/
int DDProcess::write_decimated(vector<string> names, size_t& nwin)
{
	string filepath = get_option(OUTPUT_NC4);

	if (filepath.empty()) {
		cout << "DDAlgorithm:: Missing L2 file name for decimation" << endl;
		return DTDB_FAIL;
	}
	NcFile* nc_output;
	try {
		nc_output = new NcFile(filepath, NcFile::write);
	}
	catch( NcException& e) {
		e.what();
		cout << "DDAlgorithm:: Failure opening L2 file for decimation" << endl;
		return DTDB_FAIL;
	}
	size_t lines = nc_output->getDim("number_of_lines").getSize();
	size_t pixels = nc_output->getDim("pixels_per_line").getSize();
	vector<size_t> startp = {0,0};
	vector<size_t> countp = {lines, pixels};
	size_t dlines = lines/nwin+1;
	size_t dpixels = pixels/nwin+1;
	vector<size_t> startd = {0,0,0};
	vector<size_t> countd = {3, dlines, dpixels};

    NcGroup dgrp = nc_output->addGroup("statistics");
    NcDim sets_dim = nc_output->addDim( "number_of_views", countd[0] );
    NcDim lines_dim = nc_output->addDim( "number_of_decimated_lines", countd[1] );
    NcDim samples_dim = nc_output->addDim( "decimated_samples_per_line", countd[2] );
    vector<NcDim> dimsd = {sets_dim, lines_dim, samples_dim};
	ddma<float,3>* ddec = new ddma<float,3>(dtype::FLOAT, startd, countd);

	for (auto &it : names) {
		string name = (string) it;
		name +=  "_" + algStr_;
		NcGroup pgrp = nc_output->getGroup("geophysical_data");
		ddma<float,2>* ddat =
				new ddma<float,2>(pgrp, name, dtype::FLOAT, startp, countp);
		if (ddat->ptr == nullptr) {
			delete ddat;
			continue;
		}
		fill(ddec->pts.origin(), ddec->pts.origin()+ddec->pts.num_elements(),0);
        string strdec = name +"_"+ to_string(nwin) +"x"+ to_string(nwin);
   		NcVar vard = dgrp.addVar(strdec, NcType(NC_FLOAT), dimsd);
        if (vard.isNull()) {
    		cout << "DDProcess::write_decimated() Failure to create netcdf vars" << endl;
    		return DTDB_FAIL;
        }
        string lname = strdec + " Mean, Standard deviation, and sample counts on grid decimated by " + to_string(nwin);
        vard.putAtt("long_name", lname);
        ddec->setFill(vard, to_string(DFILL_FLOAT));
		vector<size_t> chunksizes = {countd[0], countd[1], countd[2]};
        vard.setCompression(bShuffleFilter, bDeflateFilter, deflateLevel);
        vard.setChunking(vard.nc_CHUNKED, chunksizes);

        for (size_t j=0; j<lines; j++) {
			for (size_t i=0; i<pixels; i++) {
				if (ddat->pts[j][i] < DFILL_TEST) continue;
				size_t iidx = i/nwin;
				size_t jidx = j/nwin;
				ddec->pts[0][jidx][iidx] += ddat->pts[j][i];
				ddec->pts[2][jidx][iidx]++;
			}
		}
		for (size_t j=0; j<dlines; j++) {
			for (size_t i=0; i<dpixels; i++) {
				if (ddec->pts[2][j][i] >9)  {
					ddec->pts[0][j][i] /= ddec->pts[2][j][i];
				} else {
					ddec->pts[0][j][i] = DFILL_FLOAT;
				}
			}
		}
        for (size_t j=0; j<lines; j++) {
			for (size_t i=0; i<pixels; i++) {
				if (ddat->pts[j][i] < DFILL_TEST) continue;
				size_t iidx = i/nwin;
				size_t jidx = j/nwin;
				ddec->pts[1][jidx][iidx] +=
						pow(ddat->pts[j][i]-ddec->pts[0][jidx][iidx],2.0);
			}
		}
		for (size_t j=0; j<dlines; j++) {
			for (size_t i=0; i<dpixels; i++) {
				if (ddec->pts[2][j][i] >9)  {
					ddec->pts[1][j][i] /= (ddec->pts[2][j][i]-1);
					ddec->pts[1][j][i] = sqrt(ddec->pts[1][j][i]);
				} else {
					ddec->pts[1][j][i] = DFILL_FLOAT;
				}
			}
		}
		ddec->putVar(vard);
		delete ddat;
	}
	delete ddec;
	delete nc_output;

	return DTDB_SUCCESS;
};


