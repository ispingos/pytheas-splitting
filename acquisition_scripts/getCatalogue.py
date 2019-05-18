#!/usr/bin/python

# basic usr params
nodes=['ISC']   # list of nodes to search/download info. Can either a FDSN node
			          # or one of 'eida-routing' and 'iris-federator'
period='2019-03-01','2019-03-31' # selection period
p1=(21.6255,37.9624) # SE point of selection rectangle
p2=(22.7428,38.5516) # NW point of selection rectangle
arrivals=False # select to download responses or not
outName='wgoc_events_example_ISC.xml' # full name for destination file of results 
outType='QUAKEML' # output type of file ( CMTSOLUTION, CNV, JSON, KML, NLLOC_OBS, 
                  #                       NORDIC, QUAKEML, SC3ML, SCARDEC, SHAPEFILE, ZMAP)
#############################################
#
#
#                  ~ getCatalogue ~
#
# Fetch event metadata from FDSN in required
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
    inv=clt.get_events(
                          starttime=UTCDateTime(period[0]),endtime=UTCDateTime(period[1]),
                          minlongitude=p1[0],minlatitude=p1[1],
                          maxlongitude=p2[0],maxlatitude=p2[1],
                          includearrivals=arrivals
                          )
    try:
        sumInv+=inv
    except NameError:
        sumInv=inv
# how many stations found?
logging.info("Found %i events" % len(sumInv))
# write output
sumInv.write(outName,outType)
logging.info("Saved results to '%s'" % outName)
try: # Py2/Py3 comp
    raw_input("Press ENTER to exit...")
except:
    input("Press ENTER to exit...")