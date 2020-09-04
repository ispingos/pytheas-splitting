#!/usr/bin/python

## input parameters ##
# catalogue of events. either a QuakeML or a simple text file with the following structure:
#   YYYY mm dd HH MM SS.ffff  LAT_in_DD LON_in_DD MAG DEPTH_in_km
catalogue = '../example/wgoc_events_example_NKUA.xml'
# which networks to download. can use '*' to get all network-station combinations.
networks = ['*']
# stations can be one of the following: 
#   - a list containing the names of the stations to download
#   - a StationXML file, from which all available stations will be acquired
#   - a radius (in km), given as simple float. All available stations in this 
#     radius around the epicenter will be downloaded.
#   - '*', a string containing only an asterisk. In this case ALL of the stations
#     in the relevant service will be acquired, so use with caution!

# stations = ['MALA','TRIZ','PSAR'] # json list of stations e.g. ['ST1','ST2']. 
                                    # If * is in the list, all available stations
                                    # will be downloaded.
# stations = "../example/wgoc_stations_eida.xml" # station xml file
stations = 20. # radius in km

# comma separated STRING of channels. * for all. can use widlcards such as 'HH?,HN?'
channels = "*" 
# comma separated STRING of SEED location codes. * for all
locations = "*"

## output parameters ##
# directory for saving data. the underlying structure will be based on the origin time
#   and will be of the format: destination/YYYY/mm/dd/YYYY-mm-dd-HH-MM-SS/
destination = '../example/data/'
# either 'MSEED' or 'SAC'
output = 'MSEED'

## other parameters ##
start = -60 # seconds relative to the pick time for start time
end = 60    # seconds relative to the pick time for end time
provider = 'NOA' # either 'eida-routing' or 'iris-federator' or a fdsn node
#############################################
#
#
#                  ~ getEvents ~
#
# Fetch event waveform data from FDSN in required
# format for the selected stations.
#
#
# Spingos I. (2019)
#############################################
# import modules
print("..:Loading modules...")
import os, logging, sys
from obspy import UTCDateTime, read_inventory, read_events
from obspy.geodetics.base import kilometer2degrees

SCRIPT_START=UTCDateTime()
# initial params
VERDATE="04/09/2020"
VERSION="3.2"
# Logging params
logging.basicConfig(
                    level=logging.DEBUG,
                    format='%(asctime)s %(levelname)s %(message)s'
                    )
if not destination.endswith("/"):
    destination+="/"
## functions ##
def getInvStations(inv):
    '''
    get station/networks from the inventory
    '''
    codeList=[]
    for net in inv:
        for sta in net:
            codeList.append("%s.%s" % (net.code,sta.code))
    codeList.sort()
    return codeList

def getServiceStations(client,event,network="*",station="*",channel="*",radius=None,lat=None,lon=None):
    '''
    get all the available stations for the event
    '''
    # if all networks or stations required, get a list
    logging.info("Searching for all networks/stations in network "+network+" for "+str(event))
    if lat is None or lon is None or radius is None:
        inv=client.get_stations(
                     network=network,station=station,channel=channel,
                     location='*',level='station',starttime=event-1,
                     endtime=event
                               )
    else:
        inv=client.get_stations(latitude=lat,longitude=lon,maxradius=kilometer2degrees(radius),level='station',
                                starttime=event-1,endtime=event)
    codeList=[]
    for net in inv:
        netcode=net.code
        for sta in net:
            stacode=sta.code
            code='.'.join((netcode,stacode))
            if code not in codeList:
                codeList.append(code)
    codeList.sort()
    logging.info("A total of "+str(len(codeList))+" potential stations for "+str(event)+" found.")
    return codeList

def makeEventFolders(events,destination):
    '''
    make the usual directory tree (year->month->day->events)
    '''
    pathDict={}
    for event in events:
        logging.info("Creating directory for "+str(event))
        year=str(event.year).zfill(4)
        month=str(event.month).zfill(2)
        day=str(event.day).zfill(2)
        # make successive directories
        fullPath=destination+year+"/"
        try:
            os.makedirs(fullPath)
        except FileExistsError:
            pass
        fullPath+=month+"/"
        try:
            os.makedirs(fullPath)
        except FileExistsError:
            pass
        fullPath+=day+"/"
        try:
            os.makedirs(fullPath)
        except FileExistsError:
            pass
        frmt="%Y-%m-%d-%H-%M-%S"
        fullPath+=event.strftime(frmt)+"/"
        try:
            os.makedirs(fullPath)
        except FileExistsError:
            logging.warning("Event folder already exists!")
        pathDict.update({str(event):fullPath})
    return pathDict


## prep ##
rad=False
# read inventory?
if not any((isinstance(stations,list),isinstance(stations,float),isinstance(stations,int))):
    inv=read_inventory(stations)
    codeList=getInvStations(inv)
# read the catalogue file
try:
    #-- test for quakeml
    with open(catalogue, 'r') as fid:
        test = fid.read()
    if not 'quakeml' in test:  # alternative should be a plain catalogue file, without any such strings
        raise IOError('File %s is not of QuakeML format' % catalogue)
    #
    evs=read_events(catalogue)
    # convert to text catalogue format
    temp = []
    for e in evs:
        line = e.preferred_origin().time.strftime("%Y %m %d %H %M %S.%f")[:-3]
        line += '{:8.4f}'.format(e.preferred_origin().latitude)
        line += '{:10.4f}'.format(e.preferred_origin().longitude)
        try:
            line += '{:6.1f}'.format(e.preferred_origin().depth/1000.) # QuakeML depth should
                                                                       # be in m
        except TypeError: # if depth is None, we'll get an 'unsupported operand' exception
            line += '{:6.1f}'.format(0)
        #line += '{:5.1f}'.format(e.preferred_magnitude().mag)
        line += '  1.0' # this the magnitude, which is not really used and since
                        # magnitude selection is not always easy automatically
                        # we use a dummy value
        temp.append(line)
except:  # if the file is not of QuakeML type, a simple text file will be read
    with open(catalogue,'r') as fid:
        temp=fid.readlines()
# read catalogue lines into a py container and make the directories
catList=[ UTCDateTime(x.strip()[:23]) for x in temp]
#catList.sort() # sort the origin times
latList=[]; lonList=[]
try:
    latList=[float(x[23:].split()[0]) for x in temp]
    lonList=[float(x[23:].split()[1]) for x in temp]
except:
    stations="*"
pathDict=makeEventFolders(catList,destination)
# get module and connect to client
if "iris-federator" == provider or "eida-routing" == provider:
    from obspy.clients.fdsn import RoutingClient as Client
    client=Client(provider)
else:
    from obspy.clients.fdsn import Client
    client=Client(provider)
# start iterations
for event in catList:
    if isinstance(stations,list):
        if '*' in stations and '*' in networks:
            codeList=getServiceStations(client,event)
        else:
            codeList=[]
            for net in networks:
                if '*' in stations:
                  try:
                    codeList+=getServiceStations(client,event,network=net)
                  except:
                    logging.exception("Could not access data for %s"%net)
                else:
                  for sta in stations:
                    code='.'.join((net,sta))
                    codeList.append(code)
            codeList.sort()
    elif isinstance(stations,float) or isinstance(stations,int):
        lat=latList[catList.index(event)]
        lon=lonList[catList.index(event)]
        rad=stations
        try:
            codeList=getServiceStations(client,event,lat=lat,lon=lon,radius=rad)
        except:
            logging.error("No available data")
            continue
    for code in codeList:
        net=code.split(".")[0]
        sta=code.split(".")[1]
        try:
            if not rad:
                logging.info("Attempting to grab "+code+" for "+str(event)+"...")
                wf=client.get_waveforms(
                                       network=net,station=sta,
                                       channel=channels,location=locations,
                                       starttime=event+start,
                                       endtime=event+end
                                       )
            else:
                logging.info("Attempting to grab data for "+str(event)+"...")
                wf=client.get_waveforms(
                                       network=",".join([x.split(".")[0] for x in codeList]),
                                       station=",".join([x.split(".")[1] for x in codeList]),
                                       channel=channels,location=locations,
                                       starttime=event+start,
                                       endtime=event+end
                                       )    
            # make sure they start/end at the same time
            if output == 'MSEED':
                # get unique station codes in mseed
                seedSta = list(set([tr.stats.station for tr in wf]))
                for stid in sorted(seedSta):
                    twf = wf.select(station=stid)
                    netid = twf[0].stats.network

                    name = '.'.join((
                                str(event.year),str(event.julday).zfill(3),
                                str(event.hour).zfill(2),str(event.minute).zfill(2),
                                str(event.second).zfill(2),
                                str(event.microsecond)[:4],
                                netid,stid,'mseed'
                            ))
                    path=pathDict[str(event)]+name
                    twf.write(path,'MSEED')
            elif output == 'SAC':
                wf.detrend('linear')
                wf.detrend('demean')
                wf.merge()
                for trace in wf:
                    loc=trace.stats.location
                    comp=trace.stats.channel
                    netid=trace.stats.network
                    stitd=trace.stats.station
                    name='.'.join((
                             str(event.year),str(event.julday).zfill(3),
                             str(event.hour).zfill(2),str(event.minute).zfill(2),
                             str(event.second).zfill(2),
                             str(event.microsecond)[:4],
                             netid,stitd,loc,comp,'SAC'
                            ))
                    path=pathDict[str(event)]+name                  
                    trace.write(path,'SAC')
            else:
                logging.error("Unknown format for output specified!")
                raise ValueError("Unknown output format")
        except KeyboardInterrupt:
            logging.error("User interrupt!")
            sys.exit()
        except:
            logging.exception("Couldn't download "+code+" for "+str(event))
        if rad:
            break

SCRIPT_END=UTCDateTime()
SCRIPT_DUR=str(SCRIPT_END-SCRIPT_START)
logging.info("events found: "+str(len(catList)))
logging.info("stations found: "+str(len(codeList)))
logging.info("time to download: "+SCRIPT_DUR)
try: # Py2/Py3 comp
    raw_input("Press ENTER to exit...")
except:
    input("Press ENTER to exit...")