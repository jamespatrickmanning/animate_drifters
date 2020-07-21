# -*- coding: utf-8 -*-
"""
Created on Mon Jul 13 09:48:00 2020

@author: JiM
routine to extract DOPPIO forecast tracks given mth, day, year
"""
import netCDF4
import datetime
import numpy as np
from matplotlib import pyplot as plt

# example for July 10, 2020
#get the data 
url='http://tds.marine.rutgers.edu/thredds/dodsC/floats/doppio_flt_20200710.nc?ocean_time[0:1:312],lon[0:1:312][0:1:31],lat[0:1:312][0:1:31],depth[0:1:312][0:1:31],temp[0:1:312][0:1:31],salt[0:1:312][0:1:31]'
nc=netCDF4.Dataset(url)
dlons=nc.variables['lon'][:].filled(np.nan)
dlats=nc.variables['lat'][:].filled(np.nan)
#dtime=nc.variables['ocean_time'] #seconds since Jan 1st 2016
#drho=nc.variables['salt']
#dtemp=nc.variables['temp']

# make ddatetime
#ddatetime=datetime.datetime(2016,1,1,0,0,0)+datetime.timedelta(seconds=dtime)
sh=np.shape(dlons)
for k in range(sh[1]): #  loop through all the drifters
    plt.plot(np.ma.getdata(dlons[:,k])[:],np.ma.getdata(dlats[:,k])[:])
