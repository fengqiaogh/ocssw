#Creator: Emerson Sirk (emerson.a.sirk@nasa.gov)

#Last Edited: 1/14/22

from netCDF4 import Dataset
from create_images import *
from MDN.parameters import get_args
from MDN.utils import get_sensor_bands
import os

def create_copy(in_file, out_file, prod_select):
# This program takes an input l2 file with geophysical data and creates an output L2 file
# with the same global attributes and groups as the input. The geophysical data is replaced
# by the product data generated by a MDN neural network.
#in_file - the directory path of the input netcdf L2 file
#out_file - the directory path for the output netcdf L2 file
#prod_select - needs to be a string with products separated by commas. Options are chl, tss, cdom, pc.

	source = Dataset(in_file)
	copy = Dataset(out_file, mode='w')

	# Create the dimensions of the file
	for name, dim in source.dimensions.items():
		copy.createDimension(name, len(dim) if not dim.isunlimited() else None)

	# Copy the global attributes
	copy.setncatts({a:source.getncattr(a) for a in source.ncattrs()})

	#copy over other groups from input file to output file
	for gname in source.groups:

		#check to make sure it is not referencing input geophysical varibales
		if (gname == 'geophysical_data') == False:
			copy.createGroup(gname) #create new group in output file with same name

			#copy all variables in the group
			for vname in source[gname].variables:
				copy.createVariable('/'+gname+'/'+vname, source[gname][vname].datatype, source[gname][vname].dimensions)
				copy[gname][vname][:] = source[gname][vname][:]

	#create output data group
	copy.createGroup("geophysical_data")

	# Create the variables in the output data

	sensor = source.instrument
	kwargs = {
		'sensor'        : sensor,
		'product'       : 'chl,tss,cdom,pc', #all available products
		'sat_bands'     : True,
		'use_ratio'     : True,
		'use_excl_Rrs'  : True,
		'model_loc'    : os.environ['OCDATAROOT'] + '/aquaverse'
	}

	# Load the bands required for the given sensor
	req_bands = get_sensor_bands(sensor, get_args(**kwargs, use_cmdline = False))

	#calculate data
	input_data = source['geophysical_data']
	bands = sorted([int(k.replace('Rrs_', '')) for k in input_data.variables.keys() if 'Rrs_' in k])
	Rrs   = extract_data(input_data, bands, req_bands)

	# Generate product estimates - 'slices' contains the index of each product within 'products'
	products, slices = image_estimates(Rrs, **kwargs, use_cmdline = False)


	product_list = prod_select.split(',') #only user selected products
	unit_list = {'chl':'unit', 'tss':'unit', 'cdom':'unit', 'pc':'unit'}

	#assign variables from calculated data products
	for product_name in product_list:
		key = product_name
		idx = slices[key]
		#type = copy.createVLType(np.float64, product_name)
		type = np.float64
		var = copy.createVariable('/geophysical_data/'+product_name, type,("number_of_lines","pixels_per_line"))
		var.units = unit_list[key]
		var_data = products[...,idx]
		s = var_data.shape
		var_data = var_data.reshape((s[0],s[1]))
		var[:] = var_data		

	# Save the file
	copy.close()