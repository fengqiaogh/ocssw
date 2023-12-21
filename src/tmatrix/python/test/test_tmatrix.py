"""
Sam Anderson
October 2017
"""

import unittest
import numpy as np
from pytmatrix.tmatrix import tmatrix


#some allowance for rounding errors etc
epsilon = 1e-7


def printUsage():
    print ""
    print "Usage: read_obc_cal_data BAND START_DATE END_DATE OUTPUT_DRECTORY"
    print "Example: read_obc_cal_data M6 20131130 20131202 /home/output_dir"
    
def getDaysInMonth(year,month):
    DaysInMonth = [ 0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 ]
    
    if (year == 2012 or year == 2016) and month == 2:
        return 29
    else:
        return DaysInMonth[month]

def getDateString(year,month,day):
    if day < 10:
        str_day = "0" + str(day)
    else:
        str_day = str(day)

    if month < 10:
        str_month = "0" + str(month)
    else:
        str_month = str(month)
        
    str_date = str(year) + str_month + str_day

    return str_date

def getExtensibleAtom(dataset_dtype):
    if dataset_dtype.name == "int8":
        extensible_atom = tables.Int8Atom()
    if dataset_dtype.name == "uint8":
        extensible_atom = tables.UInt8Atom()
    elif dataset_dtype.name == "int16":
        extensible_atom = tables.Int16Atom()
    elif dataset_dtype.name == "uint16":
        extensible_atom = tables.UInt16Atom()
    elif dataset_dtype.name == "int32":
        extensible_atom = tables.Int32Atom()
    elif dataset_dtype.name == "uint32":
        extensible_atom = tables.UInt32Atom()
    elif dataset_dtype.name == "int64":
        extensible_atom = tables.Int64Atom()
    elif dataset_dtype.name == "uint64":
        extensible_atom = tables.UInt64Atom()
    elif dataset_dtype.name == "float32":
        extensible_atom = tables.Float32Atom()
    elif dataset_dtype.name == "float64":
        extensible_atom = tables.Float64Atom()
        
    return extensible_atom
def run_tests():
    """Tests for the T-matrix code.

       Runs several tests that test the code. All tests should return ok.
       If they don't, please contact the author.
    """
    suite = unittest.TestLoader().loadTestsFromTestCase(TMatrixTests)
    unittest.TextTestRunner(verbosity=2).run(suite)
    
    
class TMatrixTests(unittest.TestCase):

    def test_reference(self):
        """
        Test basic table generation fortran code
        """
        tmr = tmatrix()
        (reff,veff,cext,csca,albedo,asym,f) = tmr.generate_tables()

        print "\nEffective radius = " + str(reff)
        print "Effective volume = " + str(veff)
        print "Extinction Cross Section = " + str(cext)
        print "Scattering Cross Section = " + str(csca)
        print "Albedo = " + str(albedo)
        print "Asymmetry = " + str(asym)
        print "Scatter = " + str(f)
        
    def test_netcdf(self, output):
        """
        Test basic table generation fortran code
        """
        tmr = tmatrix()
        (reff,veff,cext,csca,albedo,asym,f) = tmr.generate_tables()
        
        str_name = getDateString(2018,02,07)
        str_out_dir = "/accounts/ssander2/test/t_matrix_random" 
        out_filename = "T-MATRIX_" + str_name + ".nc"  
        out_file_path = str_out_dir + "/" + out_filename
              
        try:
            ncfileout = tables.openFile(out_file_path, "w", "Scattering matrices")
        except Exception, e:
            print "-- ERROR!  Invalid output filename entered --"
            sys.exit()
        
        extensible_float32 = tables.Float32Atom()
        effrad_list = []    
        effvol_list = []    
        extx_list = []    
        sctrx_list = []    
        albedo_list = []    
        asymm_list = []    
        scatter_list = []    
        
        tgroup = ncfileout.createGroup("/", "test1", 'test')
        
        effrad_list.append(ncfileout.createEArray(sv_group, 'Effective_Radius', extensible_float32, (0), "Effective Radius"))
        effrad_list.append(ncfileout.createEArray(sv_group, 'Effective Volume', extensible_float32, (0), "Effective Volume"))
        effrad_list.append(ncfileout.createEArray(sv_group, 'Extinction Cross Section', extensible_float32, (0), "Extinction Cross Section"))
        effrad_list.append(ncfileout.createEArray(sv_group, 'Scattering Cross Section', extensible_float32, (0), "Scattering Cross Section"))
        effrad_list.append(ncfileout.createEArray(sv_group, 'Asymmetry', extensible_float32, (0), "Asymmetry"))
        effrad_list.append(ncfileout.createEArray(sv_group, 'Scatter', extensible_float32, (0,tmr.angles,4,4), "Scatter"))
        


if __name__ == '__main__':
#    unittest.main()
    run_tests()
