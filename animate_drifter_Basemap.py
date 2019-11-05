"""
Purpose: To create drifter animation from .dat file containing drifter lat, lot data
Original Author: Tanya Stoyanova (Cape Cod Community College Computer Science student)
Updated: 12 Dec 2018 by JiM where he:
   - stored in https://github.com/jamespatrickmanning/animate_drifters
   - added option for NCEP wind assuming user has downloaded the file needed
   - added option to color tracks according to drifter 'type'
   - added option to plot model vector fields using "uvmodel_function.uvmodel_plot" call assuming ctrl_uvmodel.csv is setup
   - considering redoing how the data is read in using pandas
Updated: 2 July 2019 by JiM to:
   - add "getgbox" function defining typical geographic areas

Note: Uses Linux "convert" so the final step might have to be done on a machine with that installed like NOVA. There are ways to do this in Python but that is documented in another file.
Note: Uses Basemap from mpl_toolkits
"""
 
from datetime import datetime, timedelta
import matplotlib
matplotlib.use('Agg')
from mpl_toolkits.basemap import Basemap
import csv
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import netCDF4 # for wind access
from math import sqrt
import uvmodel_function # needed to call "uvmodel_plot()"
#sys.path.append("mygit/modules")
#import basemap as bm # this is only needed when you might need more detail coastlines
from operator import itemgetter

#### HARDCODES ########
color_mode='type' # where type refers to drifter type (surf,drogue,etc) otherwise 'id'
area='BoF' # region to feed to "getgbox" function defining lat/lon box options include SNE, NorthShore, CCBAY, BoF
include_model_vectors='no' # runs the uvmodel routine inside loop
include_wind='no'
include_moorings='yes'
lat_w,lon_w=41.2,-71.6 # base of wind vector
ZOOM_STEPS = 5   #Number of frames until it reaches zoomed in state
dSIZE = 1.0  #marker size
dCoords = .1   #degrees to increase visible map area / margin around the data 
INPUT_FILENAME = '/net/pubweb_html/drifter/drift_ecohab_2019_3.dat' # you may need to download this to your machine (see nefsc.noaa.gov/drifter/drift_ecohab_2019_3.dat)
PATH_ANIM = ''#'c:\\Users\\Tanya\\eMOLT\\eMOLT\\src\\animations\\'
PATH_IMG = 'animations/frame'#'c:\\Users\\Tanya\\eMOLT\\eMOLT\\src\\images\\%03d.png'
PATH_FFMPEG = ''#c:\\ffmpeg\\'

def getgbox(area):
  # gets geographic box based on area
  if area=='SNE':
    gbox=[-70.,-65.,38.,42.] # for SNE
  elif area=='NorthShore':
    gbox=[-71.,-70.,42.,43.] # for north shore
  elif area=='CCBAY':
    gbox=[-70.75,-69.8,41.5,42.23] # CCBAY
  elif area=='BoF':
    gbox=[-68.5,-65.5,44.,45.25] # Bay of Fundy
  else:
    gbox=[]
  return gbox

def get_wind_ncep(starttime,endtime,lat,lon):
        #function get a time series of u & v wind m/s from pre-downloaded ncep file
        year=starttime.year
        url_input=""
        url_uwind=url_input+'uwnd.sig995.'+str(year)+'.nc'## where I downloaded these from ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis/surface/
        url_vwind=url_input+'vwnd.sig995.'+str(year)+'.nc'
        print url_vwind
        ncv=netCDF4.Dataset(url_vwind)
        ncu=netCDF4.Dataset(url_uwind)
        
        t=ncv.variables['time'][:]
        u_wind=ncu.variables['uwnd']
        v_wind=ncv.variables['vwnd']
        
        LAT=ncu.variables['lat'][:]
        LON=ncu.variables['lon'][:]
        
        for i in range(len(LON)):#transfer lon from (0,360) to (-180.180)
            if(LON[i]>180):
                LON[i]=-360+LON[i]
        ###########find index of nearest point###########
        index=[]
        d=[]
        for a in np.arange(len(LAT)):
            d1=[]
            for b in np.arange(len(LON)):
                d2=sqrt((LAT[a]-lat)*(LAT[a]-lat)+(LON[b]-lon)*(LON[b]-lon))
                d1.append(d2)
            d.append(d1) 
        #print np.argmin(d)/len(LON), np.argmin(d)%len(LON),np.hstack(d)[np.argmin(d)],d[np.argmin(d)/len(LON)][np.argmin(d)%len(LON)]
        index.append(np.argmin(d)/len(LON))#index of LAT
        index.append(np.argmin(d)%len(LON))#index of LON
        #print index
        cptime="%i,01,01,00,00"  %year
        cptimes=datetime.strptime(cptime, '%Y,%m,%d,%H,%M')
        moddate=[]
        for k in range(len(t)):
            moddate.append(cptimes+timedelta(hours=t[k]-t[0]))
        time_s=(starttime-cptimes).total_seconds()
        timeindex_s=int(time_s/60/60/6)# index of startime in model, NCEP data's time interval is 6 hours
        time_e=(endtime-cptimes).total_seconds()
        timeindex_e=int(time_e/60/60/6)# index of endtime in model

        u=[]
        v=[]
        u=u_wind[timeindex_s:timeindex_e,index[0],index[1]]#uwd[time,lat,lon]
        v=v_wind[timeindex_s:timeindex_e,index[0],index[1]]
        dtime=t[timeindex_s:timeindex_e]
        dtimes=[]
        for k in range(len(dtime)):
            times=starttime+timedelta(hours=(dtime[k]-dtime[0]))
            dtimes.append(times)   
        return dtimes,u,v


def selectUnique(myList):
    # order preserving
    noDupes = []
    [noDupes.append(i) for i in myList if not noDupes.count(i)]
    return noDupes
 
 
def yearday2date(my_day_number, theYear=2018):
   
    #find the ordinal number of the first day of the current year
    now = datetime.today()
    first_of_year = now.replace(year=theYear, month=1, day=1)
    first_ordinal = first_of_year.toordinal()
    #add the day of year that we are interested in and subtract 1
    my_day_ordinal = first_ordinal - 1 + my_day_number
    #convert the ordinal to date
    my_date = datetime.fromordinal(int(my_day_ordinal)) + timedelta(days=my_day_number % 1)
    #return in YYYY-MM-DD
    return my_date #my_date.isoformat() # or use strftime() to print it in different formats
 
def getLatLon(drift_id):
    lat = []
    lon = []
    if drift_id == 0:                                    
        for vert in verts:
            lon.append(float(vert[7].strip()))
            lat.append(float(vert[8].strip()))
            
        return lat, lon                                
    else:              
        for vert in verts:
            for v in vert:
                v.strip()
            if vert[0] == drift_id and (vert[12] < zEnd) and (vert[12] > zStart):
                lon.append(float(vert[7].strip()))
                lat.append(float(vert[8].strip()))
                                                        
        return lat, lon  
 
def selectMapCoords(lon, lat):

    llLon = round((min(lon) - dCoords), 3)
    print str(llLon)
    urLon = round((max(lon) + dCoords), 3)
    print str(urLon)
    llLat = round((min(lat) - dCoords), 3)
    print str(llLat)
    urLat = round((max(lat) + dCoords), 3)
    print str(urLat)
    return llLon, urLon, llLat, urLat  

def makeMap(llLon, urLon, llLat, urLat):
    #set up the map in a Equidistant Cylindrical projection
    m = Basemap(projection='cyl', llcrnrlat=llLat, urcrnrlat=urLat, \
                    llcrnrlon=llLon, urcrnrlon=urLon, resolution='h')
    return m

######## MAIN PROGRAM ######################### 

gbox=getgbox(area)# gets geographic box based on hardcoded "area='BoF', for example
verts, vert, rows = [], [], []
#read the data file
csv.register_dialect('spacedelimitedfixedwidth', delimiter=' ', skipinitialspace=True, quoting=csv.QUOTE_NONE)    
#with open('C:\\Users\\Tanya\\Desktop\\drift_2012\\' + INPUT_FILENAME, 'rb') as f:
with open(INPUT_FILENAME, 'rb') as f:
    dataReader = csv.reader(f, 'spacedelimitedfixedwidth')
    dataReader.next()
    try:
        for row in dataReader:
            if row not in verts:
                verts.append(row)
    except csv.Error, e:
        sys.exit()


#Delete duplicate data rows
verts = selectUnique(verts)   

lat, lon = [], []   
sub_verts, sub_ids, cur_ids, myIds, new_verts = [], [], [], [], []


colors = [
        "#FF0000", "#FFCC00", "#0000FF", "#00FF00", "#FF00FF", "#00FFFF", "#000000",
        "#800000", "#008000", "#000080", "#808000", "#800080", "#008080", "#808080",
        "#C00000", "#00C000", "#0000C0", "#C0C000", "#C000C0", "#00C0C0", "#C0C0C0",
        "#400000", "#004000", "#000040", "#404000", "#400040", "#004040", "#404040",
        "#200000", "#002000", "#000020", "#202000", "#200020", "#002020", "#202020",
        "#600000", "#006000", "#000060", "#606000", "#600060", "#006060", "#606060",
        "#A00000", "#00A000", "#0000A0", "#A0A000", "#A000A0", "#00A0A0", "#A0A0A0",
        "#E00000", "#00E000", "#0000E0", "#E0E000", "#E000E0", "#00E0E0", "#E0E0E0"]
    
allIds = [vert[0] for vert in verts]
uniqueIds = selectUnique(allIds)

#Figure out the year and create date 
years = []
prevYearday = 0

increasedYear = 0
year = int('20' + uniqueIds[0][0:2])
prevId = uniqueIds[0]

for vert in verts:
    if vert[0] != prevId and increasedYear == 1:
        increasedYear = 0
        #year -= 1 #JiM commented out
    if (float(vert[6]) - prevYearday) < 0 and vert[0] == prevId and increasedYear == 0:
        #year += 1 # JiM commented this out in Dec 2018
        increasedYear = 1
    vert.append(str(year))
    years.append(year)
    #print vert
    prevYearday = float(vert[6])
    prevId = vert[0]
       
d = 0
dates = []
for vert in verts:
    vert.append(yearday2date(float(vert[6]), years[d]))
    dates.append(yearday2date(float(vert[6]), years[d]))
    d = d + 1


lat, lon = getLatLon(0)
#sort data by time
verts.sort(key=lambda x: x[12]) 

#ask user to select parameters 
strIds = raw_input("If you want to see only specific drifters, enter their ids separated by space:")
myIds = strIds.split()
zID = raw_input("If you want to ZOOM specific drifter, enter id:")
zID.strip()
if zID:
    zStart = raw_input("Zoom start time: ")
    zStart = datetime.strptime(zStart, "%Y-%m-%d")
    zEnd = raw_input("Zoom end time: ")
    zEnd = datetime.strptime(zEnd, "%Y-%m-%d")
tailLength = timedelta(hours=int(raw_input("Please enter length of tail (in hours): ")))
tailLength = tailLength or datetime.timedelta(hours=6)
hrsBtwFrames = timedelta(hours=int(raw_input("Please enter hours between frames: ")))
#default start and end time
dft_start_time = min(dates)
dft_end_time = max(dates)
 
start_date = raw_input('Enter start time: [%s]' % dft_start_time)
if start_date:
    start_date = datetime.strptime(start_date, "%Y-%m-%d")
start_time = start_date or dft_start_time

end_date = raw_input('Enter end time: [%s]' % dft_end_time)
if end_date:
    end_date = datetime.strptime(end_date, "%Y-%m-%d")
end_time = end_date or dft_end_time

 
#If the user wants to see specific ids only, remove the unwanted data
if myIds:
    for myId in myIds:
        for vert in verts:
            if (vert[0].strip() == myId):
                new_verts.append(vert)
    verts = new_verts;


if zID:
    zLat, zLon = getLatLon(zID)
    xmin, xmax, ymin, ymax = selectMapCoords(zLon, zLat)

if include_wind=='yes':   
  [dtimes,u,v]=get_wind_ncep(start_time,end_time,lat_w,lon_w)

myColors = dict([(i, colors[uniqueIds.index(i)]) for i in uniqueIds])
#draw drifter track
i = 1
plotMade = False
st_time = start_time
mapChanged = 0
zoom = 0
zoom_verts = []

fig = plt.figure(1)
ax = fig.add_subplot(111)
if len(gbox)==0:
  llLon, urLon, llLat, urLat = selectMapCoords(lon, lat)
else:
  llLon, urLon, llLat, urLat=gbox
m = makeMap(llLon, urLon, llLat, urLat)
m.drawcoastlines()
m.fillcontinents(color='gray')
#m.drawrivers()
m.drawparallels(np.arange(int(llLat), int(urLat)+1,.1), labels=[1, 0, 0, 0])
m.drawmeridians(np.arange(int(llLon), int(urLon)+.1), labels=[0, 0, 0, 1])
#m.drawparallels(np.arange(int(llLat), int(urLat)+1,round((urLat-llLat)/3,1)), labels=[1, 0, 0, 0])
#m.drawmeridians(np.arange(int(llLon), int(urLon)+1,round((urLon-llLon)/3,2)),  labels=[0, 0, 0, 2])
count = ZOOM_STEPS
if include_wind=='yes':
  ax.quiver(lon_w,lat_w-0.3,10.0,0.0,color='red',scale=50.,zorder=2)
  ax.text(lon_w-.1,lat_w-0.35,'10 m/s (~20 knots)',color='red',zorder=2)
while st_time < end_time - tailLength:
    if include_wind=='yes':
      min_list= [abs(x-st_time) for x in dtimes] # find the nearest wind time to this time
      idex=min(enumerate(min_list), key=itemgetter(1))[0] # get the index of this case
      wind_arrow=ax.quiver(lon_w,lat_w,u[idex],v[idex],scale=50.0,color='red',zorder=2) # plot an arrow
    if include_moorings=='yes':
      plt.plot(-66.4856,44.7033,'*',markersize=20,color='k')
      plt.plot(-66.9450,44.3820,'*',markersize=20,color='k')
      plt.plot(-67.2950,44.5680,'*',markersize=20,color='k')
      plt.plot(-68.1448,44.1190,'*',markersize=20,color='k')
    print str(st_time)
    del sub_verts[:]
    del sub_ids[:]
    del cur_ids[:]
    for vert in verts:  
        if st_time < vert[12] < st_time + tailLength:
            sub_verts.append(vert)
            cur_ids.append(vert[0])
    
    #find unique ids in this time period        
    sub_ids = selectUnique(cur_ids) 
    #check if the drifter id that should be zoomed is among the sub_ids and within the specified time period
    for k in sub_ids:
        if zID == k.strip():
            if (st_time < zStart < (st_time + tailLength)) or (st_time > zStart and zEnd > st_time + tailLength):
                zoom = 1
      
  
    if zoom == 1:
        del sub_ids[:]        
        sub_ids.append(zID)
        del zoom_verts[:]
        zoom_verts = [vert for vert in sub_verts
                      if vert[0] == zID]
                        
        sub_verts[:] = zoom_verts
        
    #print "After checking for ZOOM sub_ids are:"
    theDate = st_time.strftime("%Y-%m-%d")
    res = []
    plt.title(theDate)
    #plt.xlabel('Longitude W', labelpad=20)
    #plt.ylabel('Latitude N', labelpad=35)
    plt.gca().set_autoscale_on(False)
    #plt.tight_layout()
    if sub_verts:          
        sub_verts.sort(key=lambda x: x[0])
        j = 0
        #print sub_ids
        for sub_id in sub_ids:
            #bm.basemap_region('sne')
            sub_lon = [float(sub_vert[7]) for sub_vert in sub_verts if sub_vert[0] == sub_ids[j]]
            sub_lat = [float(sub_vert[8]) for sub_vert in sub_verts if sub_vert[0] == sub_ids[j]]
            dep     = [float(sub_vert[9]) for sub_vert in sub_verts if sub_vert[0] == sub_ids[j]]
            sizes = [((k + 1) * dSIZE) for k in range (len(sub_lon))]
            lwid=[((k + 1) * dSIZE) for k in range (len(sub_lon))]
            if color_mode=='type': # color by drifter type
              if dep[0]==-1.0: # standard surface drifter
                  col='r'
              elif dep[0]==-0.3:# wooden turtle
                  col='b'
              elif dep[0]==-0.05:# lost drogue
                  col='y'
              else:
                  col='c' # drogue
            else:
              col=myColors[sub_id] # color by drifter id
            if sub_lat and sub_lon:
                ax.plot(sub_lon[-1],sub_lat[-1],'.', color=col,markersize=dSIZE*len(sub_lon)*2,markeredgecolor=col)
                for jj in range(len(sub_lat)-1):
                   ax.plot([sub_lon[jj],sub_lon[jj+1]],[sub_lat[jj],sub_lat[jj+1]],'-', color=col,linewidth=lwid[jj])
            del sub_lon[:] 
            del sub_lat[:]
            plotMade = True
            j = j + 1
  
    st_time = st_time + hrsBtwFrames 

    #change axis (zoom) 
    if zoom == 1 and mapChanged == 0:
        if count > 0:
            a = plt.gca()
            print "Before change: "
            print count
            print ax.get_xlim()
            print ax.get_ylim()
            x1, x2 = ax.get_xlim()
            y1, y2 = ax.get_ylim()
            ax.set_xlim([(x1 - float((llLon - xmin) / ZOOM_STEPS)), (x2 + float((urLon - xmax) / ZOOM_STEPS))])
            ax.set_ylim([(y1 - float((llLat - ymin) / ZOOM_STEPS)), (y2 - float((urLat - ymax) / ZOOM_STEPS))])
            print "After change: "
            print ax.get_xlim()
            print ax.get_ylim()
            count -= 1
            
        else:
            mapChanged = 1
    if zoom == 0 and mapChanged == 1:
        if count < ZOOM_STEPS:
            x1, x2 = ax.get_xlim()
            y1, y2 = ax.get_ylim()
            ax.set_xlim([(x1 + float((llLon - xmin) / 20)), (x2 - float((urLon - xmax) / 20))])
            ax.set_ylim([(y1 + float((llLat - ymin) / 20)), (y2 + float((urLat - ymax) / 20))])
            count += 1
    if include_model_vectors=='yes':
        Q=uvmodel_function.uvmodel_plot(ax,st_time)
    #save plot
    if plotMade == True:
        filename = PATH_IMG +str('%05d' % i) + '.png'
        i = i + 1
        plt.show()
        plt.savefig(filename, dpi=100)
        #print theDate
        plotMade = False
        #kk=get(ax.quiver)
        #del kk # deletes quiver
        if include_wind=='yes':
          wind_arrow.set_visible(False) # removes wind arrow
        if include_model_vectors=='yes':
          Q.set_visible(False) # removes model vectors
        del ax.lines[:]
        for r in res:
            r.remove()
        del res[:]
                
    zoom = 0  #clear zoom
    
     
#make animation 
try:
    # Make .flv
    #anim_name = INPUT_FILENAME[0:-4] + '.flv' #Change extension to create other type of files; Make sure that ffmpeg supports it
    anim_name = INPUT_FILENAME[0:-4] + '.gif'
    #cmd = PATH_FFMPEG + 'ffmpeg.exe -i ' + PATH_IMG + '  -r 15 -b 614000 ' + PATH_ANIM + anim_name
    cmd = 'convert -delay 10 -loop 0 '+ PATH_IMG + '*.png ' + anim_name
    os.system(cmd)
    print "Animation was created successfully."
 
except:
    print  "Could not create animation."  


print "You should probably delete the png files."
plt.close() 
 
