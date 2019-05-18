#!/usr/bin/python

# basic usr params
nodes=['eida-routing'] # list of nodes to search/download info. Can either a FDSN node
                       # or one of 'eida-routing' and 'iris-federator'
period='2018-03-01','2050-01-01' # selection period
p1=(21.6255,37.9624) # SE point of selection rectangle
p2=(22.7428,38.5516) # NW point of selection rectangle
level='response' # information layer to download (network,station,channel,response)
outName='wgoc_stations_eida.txt' # full name for destination file of results 
outType='STATIONXML'    # output type of file (CSS, KML, SACPZ, SHAPEFILE, STATIONTXT, STATIONXML)
                        # use 'CLASSTXT' for a custom text format
#############################################
#
#
#                  ~ getStations ~
#
# Fetch station metadata from FDSN in required
# format for selected area and time period from 
# selected nodes.
#
#
# Spingos I. (2019)
#############################################
# defs
VERSION = '1.0.0'
VERDATE = '15/04/2019'
## imports ##
import logging
from obspy import UTCDateTime
if "iris-federator" in nodes or "eida-routing" in nodes:
    from obspy.clients.fdsn import RoutingClient as Client
else:
    from obspy.clients.fdsn import Client
# Logging params
logging.basicConfig(
                    level=logging.DEBUG,
                    format='%(asctime)s %(levelname)s %(message)s'
                    )

# setup client #
for node in nodes:
    logging.info("Trying to get information from node '%s'" % node)
    clt=Client(node)
    inv=clt.get_stations(
                          starttime=UTCDateTime(period[0]),endtime=UTCDateTime(period[1]),
                          minlongitude=p1[0],minlatitude=p1[1],
                          maxlongitude=p2[0],maxlatitude=p2[1],
                          level=level
                          )
    try:
        sumInv += inv
    except NameError:
        sumInv = inv

# how many stations did we find?
# correct end times while iterating to avoid bugs with matching metadata
now=UTCDateTime()
Nnet=len(sumInv)
Nsta=0; Nchan=0
logging.info("Correcting end-times...")
for i,net in enumerate(sumInv):
    Nsta+=len(net)
    netEnd=net.end_date
    if not netEnd is None:
      if netEnd > now:
        sumInv[i].end_date=None
    for j,sta in enumerate(net):
        Nchan+=len(sta)
        staEnd=sta.end_date
        if not staEnd is None:
          if staEnd > now:
            sumInv[i][j].end_date=None
        for k,chan in enumerate(sta):
          chanEnd=chan.end_date
          if not chanEnd is None:
            if chanEnd > now:
              sumInv[i][j][k].end_date=None          

logging.info("Found %i networks, %i stations and %i channels" % (Nnet,Nsta,Nchan))

# write output
if outType is 'CLASSTXT':
    txtLines = []
    for net in sumInv:
        for sta in net:
            line = sta.code.ljust(5) + net.code.rjust(4) + \
                   '{:8.4f}'.format(sta.latitude) + \
                   '{:10.4f}'.format(sta.longitude) + \
                   '{:6.1f}'.format(sta.elevation) + '\n'
            txtLines.append(line)
    with open(outName,'w') as fid: 
        fid.writelines(txtLines)
else:
    sumInv.write(outName,outType)
logging.info("Saved results to '%s'" % outName)
try: # Py2/Py3 comp
    raw_input("Press ENTER to exit...")
except:
    input("Press ENTER to exit...")