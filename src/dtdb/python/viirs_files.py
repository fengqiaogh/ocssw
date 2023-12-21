'''
* NAME: viirs_files.py
*
* DESCRIPTION: VIIRS product and LUT object classes and supporting data structures.
* This is not an executable but is required for cal_viirs.py and obc_viirs.py.

* Created on July 15, 2016

* Author: Samuel Anderson
'''

import os
import sys

l1b_type = {"OBC_S":"OBC",
             "L1B-D":"L1B-DNB",
             "L1B-I":"L1B-IMG",
             "L1B-M":"L1B-MOD",
             "GEO-D":"GEO-DNB",
             "GEO-I":"GEO-IMG",
             "GEO-M_":"GEO-MOD",
             "L1A_S":"L1A"
             }

vcst_type = {"SVIMG_":"SVIMG",
             "SVMOD_":"SVMOD",
             "SVDNB_":"SVDNB",
             "IVCDB_":"IVCDB",
             "IVOBC_":"IVOBC",
             "GEOGRC":"GEOGRC",
             "GEODEG":"GEODEG",
             "GEOTC_":"GEOTC",
             "V-RVIR":"V-RVIRS"
             }

idps_type = {"SVI01":"VIIRS-I1-SDR",
             "SVI02":"VIIRS-I2-SDR",
             "SVI03":"VIIRS-I3-SDR",
             "SVI04":"VIIRS-I4-SDR",
             "SVI05":"VIIRS-I5-SDR",
             "SVM01":"VIIRS-M1-SDR",
             "SVM02":"VIIRS-M2-SDR",
             "SVM03":"VIIRS-M3-SDR",
             "SVM04":"VIIRS-M4-SDR",
             "SVM05":"VIIRS-M5-SDR",
             "SVM06":"VIIRS-M6-SDR",
             "SVM07":"VIIRS-M7-SDR",
             "SVM08":"VIIRS-M8-SDR",
             "SVM09":"VIIRS-M9-SDR",
             "SVM10":"VIIRS-M10-SDR",
             "SVM11":"VIIRS-M11-SDR",
             "SVM12":"VIIRS-M12-SDR",
             "SVM13":"VIIRS-M13-SDR",
             "SVM14":"VIIRS-M14-SDR",
             "SVM15":"VIIRS-M15-SDR",
             "SVM16":"VIIRS-M16-SDR",
             "SVM16":"VIIRS-M16-SDR",
             "SVDNB":"VIIRS-DNB-SDR",
             "IVOBC":"VIIRS-OBC-IP",
             "IVCDB":"VIIRS-DualGain-Cal-IP",
             "RNSCA":"VIIRS-Sci-RDR",
             "GDNBO":"VIIRS-DNB-SDR",
             "GMODO":"VIIRS-MOD-SDR",
             "GIMGO":"VIIRS-IMG-SDR"
             }

lut_type = {"TERRAIN-PATH":(0,"TERRAIN_PATH_LUN"),
            "LANDWATER-PATH":(0,"LANDWATER_PATH_LUN"),
            "CMNGEO-LEAPSEC":(0,"LEAPSEC_FILE_LUN"),
            "CMNGEO-PLANET-EPHEMERIS-LUT":(0,"JPL_EPHEMERIS_LUN"),            
            "CMNGEO-SAA-LUT":(0,"SAA_COEFFS_LUN"),
            "CMNGEO-POLAR-WANDER-LUT":(0,"POLAR_WANDER_LUN"),
            "CMNGEO-PARAM-LUT":(120,"CMNGEO_PARAM_LUT_LUN"),
            "VIIRS-SDR-GEO-DNB-PARAM-LUT":(6664,"GEOLOCATION_DNB_PARAMS_LUN"),
            "VIIRS-SDR-GEO-IMG-PARAM-LUT":(6664,"GEOLOCATION_IMG_PARAMS_LUN"),
            "VIIRS-SDR-GEO-MOD-PARAM-LUT":(6664,"GEOLOCATION_MOD_PARAMS_LUN"),
            "VIIRS-SDR-DG-ANOMALY-DN-LIMITS-LUT":(768,"SDR_DG_ANOMALY_DN_LIMITS_LUT_LUN"),
            "VIIRS-SDR-DNB-DN0-LUT":(1560576,"DNB_DN0_LUT_LUN"),
            "VIIRS-SDR-DNB-RVF-LUT":(520192,"DNB_RVS_LUT_LUN"),
            "VIIRS-SDR-DNB-STRAY-LIGHT-LUT":(88,"SDR_DNB_STRAY_LIGHT_LUT_LUN"),
            "VIIRS-SDR-DNB-FRAME-TO-ZONE-LUT":(16256,"DNB_FRAME_TO_ZONE_LUN"),
            "VIIRS-SDR-DNB-C-COEFFS-LUT":(73728,"DNB_C_COEFF_LUT_LUN"),
            "VIIRS-SDR-DNB-F-PREDICTED-LUT":(24592,"SDR_DNB_F_PREDICTED_LUT_LUN"),
            "VIIRS-SDR-F-PREDICTED-LUT":(101408,"SDR_F_PREDICTED_LUT_LUN"),
            "VIIRS-SDR-GAIN-LUT":(176,"SDR_GAIN_LUT_LUN"),
            "VIIRS-SDR-HAM-ER-LUT":(120024,"SDR_HAM_ER_LUT_LUN"),
            "VIIRS-SDR-RTA-ER-LUT":(120024,"SDR_RTA_ER_LUT_LUN"),
            "VIIRS-SDR-OBC-ER-LUT":(120024,"SDR_OBC_ER_LUT_LUN"),
            "VIIRS-SDR-OBC-RR-LUT":(120024,"SDR_OBC_RR_LUT_LUN"),
            "VIIRS-SDR-EBBT-LUT":(16800056,"SDR_EB_BT_LUT_LUN"),
            "VIIRS-SDR-TELE-COEFFS-LUT":(928,"SDR_TELEM_COEFFS_LUN"),
            "VIIRS-SDR-SOLAR-IRAD-LUT":(799216,"SDR_SOLAR_IRAD_LUT_LUN"),
            "VIIRS-SDR-RSR-LUT":(224056,"SDR_RSR_LUT_LUN"),
            "VIIRS-SDR-OBS-TO-PIXELS-LUT":(25216,"SDR_OBS_TO_PIXELS_LUT_LUN"),
            "VIIRS-SDR-RADIOMETRIC-PARAM-LUT":(220,"SDR_RADIOMETRIC_PARAMS_LUN"),
            "VIIRS-SDR-RADIOMETRIC-PARAM-V2-LUT":(232,"SDR_RADIOMETRIC_PARAMS_LUN"),
            "VIIRS-SDR-RADIOMETRIC-PARAM-V3-LUT":(236,"SDR_RADIOMETRIC_PARAMS_LUN"),
            "VIIRS-SDR-QA-LUT":(3812,"SDR_QA_LUT_LUN"),
            "VIIRS-SDR-EMISSIVE-LUT":(48,"SDR_EMISSIVE_LUT_LUN"),
            "VIIRS-SDR-REFLECTIVE-LUT":(8,"SDR_REFLECTIVE_LUT_LUN"),
            "VIIRS-SDR-RVF-LUT":(17533440,"SDR_RVS_LUT_LUN"),
            "VIIRS-SDR-BB-TEMP-COEFFS-LUT":(240,"SDR_BB_TEMP_CONSTANTS_LUN"),
            "VIIRS-SDR-COEFF-A-LUT":(901120,"SDR_COEFF_A_LUT_LUN"),
            "VIIRS-SDR-COEFF-B-LUT":(450560,"SDR_COEFF_B_LUT_LUN"),
            "VIIRS-SDR-DELTA-C-LUT":(4505920,"SDR_DELTA_C_LUT_LUN"),
            "VIIRS-SOLAR-DIFF-VOLT-LUT":(1280,"SOLAR_DIFF_VOLT_LUT_LUN"),
            "VIIRS-SOLAR-DIFF-AGG-HISTORY-AUX-AC":(184,"SOLAR_DIFF_HIST_AGG_AUX_LUN"),
            "VIIRS-SOLAR-DIFF-AGG-LUT":(1208,"SOLAR_DIFF_AGG_LUT_LUN"),
            "VIIRS-SOLAR-DIFF-PROC-COEFFS-LUT":(176,"SOLAR_DIFF_PROC_COEFFS_LUN"),
            "VIIRS-SOLAR-DIFF-REFL-LUT":(2296848,"SOLAR_DIFF_REFLECTANCE_LUT_LUN"),
            "VIIRS-SOLAR-DIFF-ROT-MATRIX-LUT":(72,"SOLAR_DIFF_ROT_MAT_LUT_LUN"),
            "VIIRS-SOLAR-DIFF-RVF-LUT":(7168,"SOLAR_DIFF_RVS_LUT_LUN"),
            "VIIRS-SOLAR-DIFF-SDSM-BRDF-LUT":(19520,"SOLAR_DIFF_SDSM_BRDF_LUT_LUN"),
            "VIIRS-SOLAR-DIFF-SDSM-TIME-LUT":(20,"SOLAR_DIFF_SDSM_TIME_LUT_LUN"),
            "VIIRS-SOLAR-DIFF-SDSM-TRANS-SCREEN-LUT":(16608,"SOLAR_DIFF_SDSM_TRN_SCR_LUT_LUN"),
            "VIIRS-SOLAR-DIFF-TRANS-SCREEN-LUT":(19520,"SOLAR_DIFF_TRANS_SCR_LUT_LUN"),
            "VIIRS-SDR-CAL-AUTOMATE-LUT":(21408,"SDR_CAL_AUTOMATE_LUT_LUN"),
            "VIIRS-SDR-ONBOARD-DNB-OFFSETS-LUT":(520192,"SDR_ONBOARD_DNB_OFFSETS_LUT_LUN"),
            "VIIRS-SDR-DNB-GAIN-RATIOS-LUT":(27648,"SDR_DNB_GAIN_RATIOS_LUT_LUN"),
            "VIIRS-SDR-DNB-LGS-GAINS-LUT":(27664,"SDR_DNB_LGS_GAINS_LUT_LUN"),
            "VIIRS-SDR-DNB-STRAY-LIGHT-CORRECTION-LUT":(487943852,"SDR_DNB_STRAY_LIGHT_CORRECTION_LUT_LUN"),
            "VIIRS-SDR-RELATIVE-SPECTRAL-RESPONSE-LUT":(240064,"SDR_RELATIVE_SPECTRAL_RESPONSE_LUT_LUN"),
            "VIIRS-RSBAUTOCAL-BRDF-RTA-VIEW-LUT":(538832,"RSBAUTOCAL_BRDF_RTA_VIEW_LUT_LUN"),
            "VIIRS-RSBAUTOCAL-BRDF-SDSM-VIEW-LUT":(158080,"RSBAUTOCAL_BRDF_SDSM_VIEW_LUT_LUN"),
            "VIIRS-RSBAUTOCAL-DNB-DARK-SIGNAL-AUTOMATE-LUT":(147472,"RSBAUTOCAL_DNB_DARK_SIGNAL_AUTOMATE_LUT_LUN"),
            "VIIRS-RSBAUTOCAL-DNB-GAIN-RATIOS-AUTOMATE-LUT":(41480,"RSBAUTOCAL_DNB_GAIN_RATIOS_AUTOMATE_LUT_LUN"),
            "VIIRS-RSBAUTOCAL-DNB-LGS-GAIN-AUTOMATE-LUT":(32264,"RSBAUTOCAL_DNB_LGS_GAIN_AUTOMATE_LUT_LUN"),
            "VIIRS-RSBAUTOCAL-DNB-MOON-ILLUMINATION-LUT":(175200,"RSBAUTOCAL_DNB_MOON_ILLUMINATION_LUT_LUN"),
            "VIIRS-RSBAUTOCAL-H-AUTOMATE-LUT":(508,"RSBAUTOCAL_H_AUTOMATE_LUT_LUN"),
            "VIIRS-RSBAUTOCAL-H-LUT":(624,"RSBAUTOCAL_H_LUT_LUN"),
            "VIIRS-RSBAUTOCAL-ROT-MATRIX-LUT":(72,"RSBAUTOCAL_ROT_MATRIX_LUT_LUN"),
            "VIIRS-RSBAUTOCAL-RSB-F-AUTOMATE-LUT":(50180,"RSBAUTOCAL_RSB_F_AUTOMATE_LUT_LUN"),
            "VIIRS-RSBAUTOCAL-RVF-LUT":(7680,"RSBAUTOCAL_RVF_LUT_LUN"),
            "VIIRS-RSBAUTOCAL-SDSM-SOLAR-SCREEN-TRANS-LUT":(106240,"RSBAUTOCAL_SDSM_SOLAR_SCREEN_TRANS_LUT_LUN"),
            "VIIRS-RSBAUTOCAL-VOLT-LUT":(1280,"RSBAUTOCAL_VOLT_LUT_LUN"),
            "VIIRS-RSBAUTOCAL-SDSM-TIME-LUT":(20,"RSBAUTOCAL_SDSM_TIME_LUT_LUN")
            }


iopath_type = {"Sci-RDR-LIST":("#VIIRS_SCIENCE_RDR_LUN","Sci_RDR_List"),
               "Sci-RDR-IN-PATH":("RDR_INPUT_PATH_LUN","Sci_RDR"),
               "OBCIP-IN-PATH":("OBC_IP_INPUT_PATH_LUN","OBCIP_INPUT"), 
               "IP-STAGING-PATH":("IP_STAGING_PATH_LUN","STAGING"), 
               "VRDR-OUT-PATH":("VERIFIED_RDR_PATH_LUN","VRDR"), 
               "GEO-OUT-PATH":("GEO_OUTPUT_PATH_LUN","GEO"), 
               "DNB-OUT-PATH":("DAYNIGHT_SDR_OUTPUT_PATH_LUN","SDR_DNB"),
               "IMG-OUT-PATH":("IMAGERY_SDR_OUTPUT_PATH_LUN","SDR_IMG"),
               "MOD-OUT-PATH":("MODERATE_SDR_OUTPUT_PATH_LUN","SDR_MOD"),
               "OBCIP-OUT-PATH":("OBC_IP_OUTPUT_PATH_LUN","OBCIP"),
               "DGIP-OUT-PATH":("DUAL_GAIN_IP_OUTPUT_PATH_LUN","DGIP"),
               "SD-OUT-PATH":("SD_OUTPUT_PATH_LUN","SD_HISTORY"),
               "CAL-HISTORY-PATH":("RSBAUTOCAL_HISTORY_AUX_LUN","HISTORY"),
               "OBCIP-HISTORY-PATH":("RSBAUTOCAL_OBCIP_HISTORY_AUX_LUN","HISTORY")
               }


                   
class l1_file:
    def __init__(self, arg1):
        self.path, self.filename = os.path.split(arg1)
        self.filepath = arg1
        self.type = self.parseL1Type(self.filename)
        self.date = self.parseL1Date(self.filename)
        self.time = self.parseL1Time(self.filename)
        self.corename = self.getCorename(self.filename)
        
    def display(self):
        print(self.filepath + "/" + self.filename)
        print(self.vtype, self.itype)
        print(self.date, self.time)

    def getStartYMD(self):
        y = int(self.date[0:4])
        m = int(self.date[4:6])
        d = int(self.date[6:8])
        return y, m, d

    def getStartHMS(self):
        s = self.time % 100
        m = (self.time - s)/100 % 100
        h = (self.time - m*100 - s)/10000
        return h, m, s

    def getEndHMS(self):
        s = self.endTime % 100
        m = (self.endTime - s)/100 % 100
        h = (self.endTime - m*100 - s)/10000
        return h, m, s
    
    def getCorename(self, str):
        pos = str.find(".L1B-M")
        if (pos > 0):
            return str[0:pos]
        pos = str.find(".L1A_SNPP")
        if (pos > 0):
            return str[0:pos]
        pos = str.find(".L1A_JPSS1")
        if (pos > 0):
            return str[0:pos]
        pos = str.find(".L1A.nc")
        if (pos > 0):
            return str[0:pos]        
        pos = str.find(".nc")
        if (pos > 0):
            return str[0:pos]        
        return ""
    
    def setOutputPath(self, arg):
        self.outputPath = arg
        if not os.path.isdir(arg):
            return False
        else:
            return True
               
    def getImgFilename(self):
        self.ImgFilename = self.corename + ".L1B-I_SNPP.nc"
        return self.ImgFilename
        
    def getModFilename(self):
        self.ModFilename = self.corename + ".L1B.nc"
        return self.ModFilename
        
    def getDnbFilename(self):
        self.DnbFilename = self.corename + ".L1B-D_SNPP.nc"
        return self.DnbFilename
        
    def getObcFilename(self):
        self.ObcFilename = self.corename + ".OBC_SNPP.nc"
        return self.ObcFilename
        
    def getCdgFilename(self):
        self.CdgFilename = self.corename + ".CDG_SNPP.nc"
        return self.CdgFilename
        
    def getGeoIFilename(self):
        self.GeoIFilename = self.corename + ".GEO-I_SNPP.nc"
        return self.GeoIFilename

    def getGeoMFilename(self):
        self.GeoMFilename = self.corename + ".GEO.nc"
        return self.GeoMFilename

    def getGeoDFilename(self):
        self.GeoDFilename = self.corename + ".GEO-D_SNPP.nc"
        return self.GeoDFilename
          
    def getDtFilename(self):
        self.DtFilename = self.corename + ".L2.AER_DT.nc"
        return self.DtFilename
          
    def getDbFilename(self):
        self.DbFilename = self.corename + ".L2.AER_DB.nc"
        return self.DbFilename
                    
    def parseL1Type(self, str):
        str_type = str[15:20]
        if str.find(".nc") < 0:
            return 0
        elif str.find("L1A") < 0:
            return "L1A"
        elif str_type not in l1b_type:
            return 0
        else: 
            return l1b_type[str_type]
    
    def parseL1Date(self, str):
        if str.startswith("OA"):
            str_date = str[2:9]
        elif str.startswith("VL1B"):
            str_date = str[13:19]
        elif str.startswith("PACE_OCI_SIM"):
            str_date = str[13:21]
        elif str.startswith("SNPP_VIIRS"):
            str_date = str[11:19]
        elif str.startswith("JPSS1_VIIRS") or str.startswith("JPSS2_VIIRS"):
            str_date = str[12:20]
        else:
            str_date = str[1:8]
        return str_date
                            
    def parseL1Time(self, str):
        if str.startswith("OA"):
            str_time = str[9:15]
        elif str.startswith("VL1B"):
            str_time = str[22:28]
        elif str.startswith("PACE_OCI_SIM"):
            str_time = str[22:28]
        elif str.startswith("SNPP_VIIRS"):
            str_time = str[20:26]
        elif str.startswith("JPSS1_VIIRS") or str.startswith("JPSS2_VIIRS"):
            str_time = str[21:27]
        else:
            str_time = str[8:14]
        return int(str_time)
 

class sdr_file:
    def __init__(self, arg1):
        self.path, self.filename = os.path.split(arg1)
        self.filepath = arg1
        self.vtype = self.parseVcstType(self.filename)
        self.itype = self.parseIdpsType(self.filename)
        self.date = self.parseSdrDate(self.filename)
        self.startTime = self.parseSdrStartTime(self.filename)
        self.endTime = self.parseSdrEndTime(self.filename)
    
    def display(self):
        print(self.filepath + "/" + self.filename)
        print(self.vtype, self.itype)
        print(self.date, self.startTime, self.endTime)

    def getStartHMS(self):
        s = self.startTime % 1000
        m = (self.startTime - s)/1000 % 100
        h = (self.startTime - m*1000 - s)/100000
        return h, m, s

    def getEndHMS(self):
        s = self.endTime % 1000
        m = (self.endTime - s)/1000 % 100
        h = (self.endTime - m*1000 - s)/100000
        return h, m, s

    def parseVcstType(self, str):
        if str.find("unscaled") > 0:
            return 0
        elif str.find(".h5") < 0:
            return 0
        elif str[:6] not in vcst_type:
            return 0
        else: 
            return vcst_type[str[:6]]
       
    def parseIdpsType(self, str):
        if str[-3:] != ".h5":
            return 0
        elif str[:5] not in idps_type:
            return 0
        else: 
            return idps_type[str[:5]]
      
    def parseSdrDate(self, str):
        pos = str.find("_d")
        return int(str[pos+2:pos+10])
                            
    def parseSdrStartTime(self, str):
        pos = str.find("_t")
        return int(str[pos+2:pos+9])
                            
    def parseSdrEndTime(self, str):
        pos = str.find("_e")
        return int(str[pos+2:pos+9])
                   
                   
class lut_file:
    def __init__(self, arg1,arg2):
        if os.path.isdir(arg2):
            self.path = arg2
            self.filename = ""
        else: 
            self.path, self.filename = os.path.split(arg2)
        self.filepath = arg2
        self.type = arg1
        self.date = parseLutDate(self.filename)
        self.lptype = lut_type[self.type][1]
    
    def display(self):
        print(self.filepath)
        print(self.type)
        print(self.date)

    def parseLutType(self, str):
        path,filename = os.path.split(str)
        pos = str.find("PolarWander")
        if pos > 0:
            return "CMNGEO-POLAR-WANDER-LUT"
        pos = str.find("VIIRS-RSBAUTOCAL-BRDF-SCREEN-TRANSMISSION-PRODUCT-RTA-VIEW-LUT")
        if pos >= 0:
            return "VIIRS-RSBAUTOCAL-BRDF-RTA-VIEW-LUT"
        pos = str.find("VIIRS-RSBAUTOCAL-BRDF-SCREEN-TRANSMISSION-PRODUCT-SDSM-VIEW-LUT")
        if pos >= 0:
            return "VIIRS-RSBAUTOCAL-BRDF-SDSM-VIEW-LUT"
        pos = filename.find("_npp_")
        if pos > 0:
            s = filename[:pos]
            if s in lut_type:
                return s
        return ""
       
    def parseLutDate(self, str):
        pos = str.find("PolarWander")
        if pos > 0:
            return int(str[pos+35:pos+43])
        pos = str.find("_npp_")
        if pos > 0:
            return int(str[pos+21:pos+29])
        return 0
     
                       
    def parseVcstLutDate(self, str):
        pos = str.find("_")
        if pos > 0:
            return int(str[pos+1:pos+9])
        return 0
 
                   
class io_path:
    def __init__(self, arg1, arg2):
        self.type = arg1
        self.path = arg2
        self.lptype = iopath_type[arg1][0]
    
    def display(self):
        print(self.type)
        print(self.path)


        