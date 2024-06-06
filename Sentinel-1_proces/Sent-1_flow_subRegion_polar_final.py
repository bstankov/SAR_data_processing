#PROCESSING FLOW FOR SENTINEL-1

#   IMPORT LIBRARIES    ******************************
baseSNAP = 'C:/Program Files/snap/bin/gpt.exe'
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import os,sys,gc, glob,matplotlib
import numpy as np

sys.path.append(r'C:\U......conda/envs/mysnapy/Lib/snappy')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import snappy
from osgeo import gdal
from snappy import Product
from snappy import ProductIO
from snappy import ProductUtils
from snappy import WKTReader
from snappy import HashMap
from snappy import GPF
from snappy import jpy, ProgressMonitor

import matplotlib.image as mpimg #for plotting
# from termcolor import colored
from zipfile import ZipFile
from os.path import join
from glob import iglob
import pandas as pd

import subprocess
#change module setting
pd.options.display.max_colwidth = 80

polar=['VH','VV']
outpath = r'.....\bobo\New folder\sent_out'
###*************** 2. User defined functions    ****************
#*******************************************************************
def output_view(product, band, min_value_VV, max_value_VV, min_value_VH, max_value_VH):
    ''' Create visualization of processed Sentinel-1 SAR data

        Keyword arguments:
        product      --snappy GPF product's band to be vusulized
        band,         --list --> product's band to be visualize
        min_value_VV,  --int --> min value for color strech in VV band
        max_value_VV,   --int --> max value for color strech in VV band
        min_value_VH,   --int --> min value for color strech in VH band
        max_value_VH    --int --> max value for color strech in VH band

    '''
    band_data_list = []
    for i in band:
        band = product.getBand(i)
        w = band.getRasterWidth()
        h = band.getRasterHeight()
        band_data = np.zeros(w * h, np.float32)
        band.readPixels(0, 0, w, h, band_data)
        band_data.shape = h, w
        band_data_list.append(band_data)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 16))
    ax1.imshow(band_data_list[0], cmap='gray', vmin=min_value_VV, vmax=max_value_VV)
    ax1.set_title(output_bands[0])
    ax2.imshow(band_data_list[1], cmap='gray', vmin=min_value_VH, vmax=max_value_VH)
    ax1.set_title(output_bands[1])
    for ax in fig.get_axes():
        ax.label_outer()
def get_snap_info(operator):
    """
    Returns information about SNAP gpt operators and their parameters
    """
    op_spi = GPF.getDefaultInstance().getOperatorSpiRegistry().getOperatorSpi(operator)
    print('Op name:', op_spi.getOperatorDescriptor().getName())
    print('Op alias:', op_spi.getOperatorDescriptor().getAlias())
    param_Desc = op_spi.getOperatorDescriptor().getParameterDescriptors()
    for param in param_Desc:
        print(param.getName(), "or", param.getAlias())

#*******************************************************************

path = '......2019/IR/radar/PR/PR_South/'
s1meta = "manifest.safe"
os.chdir(path)
#list to store all products
products = []
pols = ["VH","VV"]
#name, sensing_mode, product_type, polarization, height, width, band_names = ([] for i in range(7))
#dates = []
for folder in os.listdir(path):
    gc.enable()
    print("folder",folder)
    output = path + folder
    timestamp = folder.split("_")[4]
    date = timestamp[:8]
    dates.append(date)
    s1prd = "%s/%s/%s.SAFE/%s" % (path, folder,folder, s1meta)
    print ("output:", output)
    print("s1prd:", s1prd)
    reader_p = ProductIO.getProductReader("SENTINEL-1")
    product_p = reader_p.readProductNodes(s1prd, None)
    products.append(product_p)
    print("Spath:",product_p)

    #Set target foledr and extract metadata
    ##THIS IS CRITICAL -DONOT FORGET WHOLE PATH, ENTIRE DIRECTORY (PATH) PLUS S1 DIRECT
    #Extract information about the Sentinel-1 GRD products:

for product in products:

    width = product.getSceneRasterWidth()
    height = product.getSceneRasterHeight()
    name = product.getName()
    band_names = product.getBandNames()
    print("Product: %s, %d x %d pixels" % (name, width, height))
    print("Bands:   %s" % (list(band_names)))
    print("name:",name)

#   *****************   3.SUBSETTING ***************************

HashMap = jpy.get_type('java.util.HashMap')
GPF.getDefaultInstance().getOperatorSpiRegistry().loadOperatorSpis()
#******************************# for georegion **********
# WKTReader = snappy.jpy.get_type('com.vividsolutions.jts.io.WKTReader') #".com to .org" is changed from SNAP 7 to SNAP 8.0 version
WKTReader = jpy.get_type('org.locationtech.jts.io.WKTReader')  #Note ".org"
#Geometry = jpy.get_type('org.locationtech.jts.geom.Geometry')
wkt = "POLYGON ((-117.35 57.02 , -117.1 57.03, -117.05 56.97, -117.3 56.97, -117.35 57.02 ))"
geom = WKTReader().read(wkt)

sub_parameters = HashMap()
sub_parameters.put('outputImageScaleInDb', False)
sub_parameters.put('copyMetadata', True)
sub_parameters.put('geoRegion', geom)

subsets = []  #keeps the list of all subsets objects

for product in products:
    subset = GPF.createProduct('Subset', sub_parameters, product)
    subsets.append(subset)
    print("output_bands:",list(subset.getBandNames()))
    # write file in linear units
    splits = product.getName().split("_")
    outsub = splits[0] + '_' +splits[4][:7] +'_sub'
   # snappy.ProductIO.writeProduct(subset, os.path.join(outpath, outsub), 'GeoTIFF')  #iskljuci

#********** 4   Apply orbit file    *********************
#****************************************************

aof_parameters = snappy.HashMap()
aof_parameters.put('Apply-Orbit-File', True)
aofs=[]
for sub in subsets:
    apply_orbit = snappy.GPF.createProduct('Apply-Orbit-File',aof_parameters,sub)
    aofs.append(apply_orbit)
    print("output_bands:",list(apply_orbit.getBandNames()))

#******    5. Thermal Noise Removeal  ************************
#**************************************************************

tnr = []  #bucket to store objects of TNR
TNR_parameters = snappy.HashMap()
TNR_parameters.put('removeThermalNoise', True)
for aof in aofs:
    thermal_noise = snappy.GPF.createProduct('ThermalNoiseRemoval',TNR_parameters,aof)
    tnr.append(thermal_noise)
#
#   **********  Step 6: Pre-processing - Calibration    *************
#**********************************************************

cal_parameters = snappy.HashMap()
cal_parameters.put('outputSigmaBand', True)
#cal_parameters.put('sourceBands','Intensity_VH', 'Intensity_VV')
cal_parameters.put('selectedPolarizations','VV')
cal_parameters.put('outputImageScaleInDb', False)

calibrates = []

for atnr in tnr:

    calibrate = GPF.createProduct('Calibration',cal_parameters, atnr)
    calibrates.append(calibrate)
    print("calib_bands:", list(calibrate.getBandNames()))
   # snappy.ProductIO.writeProduct(subset, os.path.join(outpath, "calibr"), 'GeoTIFF')
#
# ***********   #Step 7: Pre-processing - Speckle filtering ****************
#***********************************************************************
speckles = []

spc_parameters = snappy.HashMap()
spc_parameters.put('filter', 'Lee')
spc_parameters.put('filterSizeX',5)
spc_parameters.put('filterSizeY',5)
spc_parameters.put('dampingFactor',2)
spc_parameters.put('edgeThreshold',5000.0)
spc_parameters.put('estimateENL',True)
spc_parameters.put('enl',1.0)
#filtered = snappy.GPF.createProduct('Speckle-Filter', spc_parameters, calibrated)

for calibrate in calibrates:

    speckle = GPF.createProduct('Speckle-Filter', spc_parameters, calibrate)
    speckles.append(speckle)

#*****  8   Terrani   correction    ***********************
#****************************************************

mySRTMdem = '....ace/recon/SRTM_30m_2015/mosaic/srtm_30m.tif'
tcs=[]
#export separate polarizaiton
for p in polar:
    polarization = p
    tc_parameters = snappy.HashMap()
    tc_parameters.put('demResamplingMethod','NEAREST_NEIGHBOUR')
    tc_parameters.put('demName','External DEM')
    tc_parameters.put('externalDEMFile',mySRTMdem)
    tc_parameters.put('imgResamplingMethod','NEAREST_NEIGHBOUR')
    tc_parameters.put('pixelSpacingInMeter',10.0)
    tc_parameters.put('sourceBands','Sigma0_'+polarization)

    proj = '''PROJCS["NAD83 / UTM zone 11N",
    GEOGCS["NAD83",
    DATUM["North_American_Datum_1983",
    SPHEROID["GRS 1980",6378137,298.257222101,AUTHORITY["EPSG","7019"]],AUTHORITY["EPSG","6269"]],
    PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],
    UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4269"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],
    PROJECTION["Transverse_Mercator"],
    PARAMETER["latitude_of_origin",0],
    PARAMETER["central_meridian",-117],
    PARAMETER["scale_factor",0.9996],
    PARAMETER["false_easting",500000],
    PARAMETER["false_northing",0],
    AUTHORITY["EPSG","26911"],
    AXIS["Easting",EAST],
    AXIS["Northing",NORTH]]');'''

    tc_parameters.put('mapProjection',proj)
    tc_parameters.put('nodataValueAtSea',False)
    tc_parameters.put('saveSelectedSourceBand',True)


    for spec in speckles:
        terrain_correction = snappy.GPF.createProduct('Terrain-Correction', tc_parameters, spec)
        tcs.append(terrain_correction)
        splits_tc = spec.getName().split("_")

        outSent_tc = splits_tc[0] + '_' +splits_tc[1] +'_tc_'+ polarization
        print( 'splits_tc ',splits_tc )
        #print("tc_bands:", list(terrain_correction.getBandNames()))
        snappy.ProductIO.writeProduct(terrain_correction, os.path.join(outpath,outSent_tc), 'GeoTIFF')

        #********   9. cONVERSION FROM linear to dB
        #*****************************************************************

        dB_parameters = snappy.HashMap()
        lineartodb = GPF.createProduct('linearToFromdB', dB_parameters, terrain_correction)
        # write file in linear units
        snappy.ProductIO.writeProduct(lineartodb , os.path.join(outpath, outSent_tc.replace('tc','dB')), 'GeoTIFF')

print (" your output files are in: ", outpath  )
