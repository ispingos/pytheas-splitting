#!/usr/bin/python

## input parameters ##

#---- PATHS ----#
# catalogue of events. either a QuakeML or a simple text file with the following structure:
# YYYY mm dd HH MM SS.ffff  LAT_in_DD LON_in_DD MAG DEPTH_in_km 
catalogue = r"../example/wgoc_events_example_NKUA.xml"

# directory for saving data. the underlying structure will be based on the origin time
#   and will be of the format: destination/YYYY/mm/dd/YYYY-mm-dd-HH-MM-SS/
destination = r'../example/data'

#---- SELECTION ----#
# which networks to download. can use '*' to get all network-station combinations.
networks = ['*']

# stations can be one of the following: 
#   - a list containing the names of the stations to download
#   - a StationXML file, from which all available stations will be acquired
#   - a radius (in km), given as simple float. All available stations in this 
#     radius around the epicenter will be downloaded.
#   - '*', a string containing only an asterisk. In this case ALL of the stations
#     in the relevant service will be acquired, so use with caution!

# stations = ['MALA', 'TRIZ', 'PSAR'] # json list of stations e.g. ['ST1','ST2']. 
                                      # If * is in the list, all available stations
                                      # will be downloaded.

# stations = "../example/wgoc_stations_eida.xml" # station xml file

stations = 20. # radius in km

# comma separated STRING of channels. * for all. can use widlcards such as 'HH?,HN?'
channels = "*" 

# comma separated STRING of SEED location codes. * for all
locations = "*"

# Select the output type of the data. use 'disp', 'vel' or 'acc'. For no 
# instrument correction, use 'none'. If 'auto', declared accelerograph data 
# will be converted to velocity, while seismograph data won't be deconvoluted.
remove_response = 'none'

#---- OUTPUT ----#
# format of written waveform, either 'MSEED' or 'SAC'
output = 'MSEED'

#---- OTHER ----#
start = -10 # seconds relative to the pick time for start time
end = 30  # seconds relative to the pick time for end time
providers = ['RESIF', 'NOA']  # either 'eida-routing' or 'iris-federator' or a list of fdsn nodes

#############################################
#
#
#        ~ fdsn_get_waveform_data ~
#
# Fetch event waveform data from FDSN in required
# format for the selected stations.
#
#
# Spingos I. (2019-2021)
#############################################
# import modules
print("..:Loading modules...")
import os, logging, sys
from obspy import UTCDateTime, read_inventory, read_events
from obspy.geodetics.base import kilometer2degrees

SCRIPT_START = UTCDateTime()
# initial params
VERDATE="03/04/2021"
VERSION="3.4.1"
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
    code_list=[]
    for net in inv:
        for sta in net:
            code_list.append("%s.%s" % (net.code,sta.code))
    code_list.sort()
    return code_list

def getServiceStations(
    client, event, 
    network="*", station="*", location='*', channel="*",
    radius=None, lat=None, lon=None
    ):
    '''
    get all the available stations for the event
    '''
    # if all networks or stations required, get a list
    logging.info("%s: Searching %s.*.%s for %s " % (provider, network, channel, str(event)))
    if lat is None or lon is None or radius is None:
        inv = client.get_stations(
                     network=network,
                     station=station,
                     channel=channel,
                     location=location,
                     level='station',
                     starttime=event - 1,
                     endtime=event + 1
                               )
    else:
        inv = client.get_stations(
            network=network,
            location=location,
            channel=channel,
            latitude=lat, 
            longitude=lon, 
            maxradius=kilometer2degrees(radius),
            level='station',
            starttime=event - 1,
            endtime=event + 1
            )
    code_list = []
    for net in inv:
        netcode=net.code
        for sta in net:
            stacode = sta.code
            code = '.'.join((netcode,stacode))
            if code not in code_list:
                code_list.append(code)
    code_list.sort()
    logging.info("%s a total of %i potential stations for %s found" % (provider, len(code_list), str(event)))
    return code_list

## prep ##
rad=False
# read inventory?
if not any((isinstance(stations, list),isinstance(stations,float),isinstance(stations,int))):
    inv = read_inventory(stations)
    code_list = getInvStations(inv)
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
active_code_list = []
try:
    latList=[float(x[23:].split()[0]) for x in temp]
    lonList=[float(x[23:].split()[1]) for x in temp]
except:
    stations="*"

#-- attach response?
if output == 'none':
    attach_response_flag = False
else:
    attach_response_flag =True
for provider in providers:
    logging.info('-- connecting to data provider %s' % provider)
    # get module and connect to client
    if "iris-federator" == provider or "eida-routing" == provider:
        from obspy.clients.fdsn import RoutingClient as Client
        client=Client(provider)
    else:
        from obspy.clients.fdsn import Client
        client=Client(provider)
    # start iterations
    for i_ev, event in enumerate(catList):
        logging.info('%s: iterating event %s/%s' % (
            provider,
            '{:,}'.format(i_ev + 1),
            '{:,}'.format(len(catList))
            ))
        if isinstance(stations,list):
            if '*' in stations and '*' in networks:
                code_list = getServiceStations(
                    client,
                    event
                    )
            else:
                code_list = []
                for net in networks:
                    if '*' in stations:
                      try:
                        code_list+=getServiceStations(
                            client,
                            event,
                            network=net,
                            channel=channels
                            )
                      except:
                        logging.exception("%s: Could not access data for %s" % (provider, net))
                    else:
                      for sta in stations:
                        code='.'.join((net,sta))
                        code_list.append(code)
                code_list.sort()
        elif isinstance(stations,float) or isinstance(stations,int):
            lat = latList[catList.index(event)]
            lon = lonList[catList.index(event)]
            rad = stations
            try:
                code_list = getServiceStations(
                    client,
                    event,
                    lat=lat,
                    lon=lon,
                    radius=rad,
                    channel=channels
                    )
            except:
                logging.error("%s: No available data" % provider)
                continue
        for code in code_list:
            net = code.split(".")[0]
            sta = code.split(".")[1]
            try:
                if not rad:
                    logging.info("%s: attempting to fetch %s for %s" % (provider, code, str(event)))
                    wf = client.get_waveforms(
                                           network=net, station=sta,
                                           channel=channels, location=locations,
                                           starttime=event + start,
                                           endtime=event + end,
                                           attach_response=attach_response_flag
                                           )
                else:
                    str_codes = ','.join(code_list)
                    logging.info("%s: attempting to fetch '%s' for %s" % (provider, str_codes, str(event)))
                    wf = client.get_waveforms(
                                           network=",".join([x.split(".")[0] for x in code_list]),
                                           station=",".join([x.split(".")[1] for x in code_list]),
                                           channel=channels,
                                           location=locations,
                                           starttime=event + start,
                                           endtime=event + end,
                                           attach_response=attach_response_flag
                                           ) 
                # deconvolve?
                if remove_response.upper() in ['DISP', 'VEL', 'ACC']:
                    wf.remove_response(output=remove_response.upper())
                if remove_response.upper() == 'AUTO':
                    #-- acc convert to velocity
                    sel_chn = '?N?'
                    wf.select(channel=sel_chn).remove_response(output='VEL')
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
                        year = str(event.year).zfill(4)
                        month = str(event.month).zfill(2)
                        day = str(event.day).zfill(2)
                        frmt = "%Y-%m-%d-%H-%M-%S"
                        path_wf_dir = os.path.join(
                            destination,
                            year,
                            month, 
                            day, 
                            event.strftime(frmt)
                                                    )
                        try:
                            os.makedirs(path_wf_dir)
                        except FileExistsError:
                            pass
                        path = os.path.join(path_wf_dir, name)                  
                        twf.write(path,'MSEED')
                        if not code in active_code_list:
                            active_code_list.append(code)
                elif output == 'SAC':
                    try:
                        wf.detrend('linear')
                        wf.detrend('demean')
                    except:
                        logging.exception("%s: couldn't download %s for %s" % (provider, code, str(event)))
                    wf.merge(method=1, fill_value=0, interpolation_samples=0)
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
                        year = str(event.year).zfill(4)
                        month = str(event.month).zfill(2)
                        day = str(event.day).zfill(2)
                        frmt = "%Y-%m-%d-%H-%M-%S"
                        path_wf_dir = os.path.join(
                            destination, 
                            year,
                            month, 
                            day, 
                            event.strftime(frmt)
                                                    )
                        try:
                            os.makedirs(path_wf_dir)
                        except FileExistsError:
                            pass
                        path = os.path.join(path_wf_dir, name)                  
                        trace.write(path,'SAC')
                        if not code in active_code_list:
                            active_code_list.append(code)
                else:
                    logging.error("%s: unknown format for output specified!" % (provider))
                    raise ValueError("Unknown output format")
            except KeyboardInterrupt:
                logging.error("user interrupt!")
                sys.exit()
            except:
                logging.exception("%s: couldn't download %s for %s" % (provider, code, str(event)))
            if rad:
                break

SCRIPT_END = UTCDateTime()
SCRIPT_DUR = SCRIPT_END - SCRIPT_START
logging.info("events found: %s" % '{:,}'.format(len(catList)))
logging.info("stations found: %s" % '{:,}'.format(len(active_code_list)))
logging.info("time to download: %s" % '{:,.2f}'.format(SCRIPT_DUR))
logging.info('catalog: %s' % catalogue)
logging.info('stored data to %s' % destination)
try: # Py2/Py3 comp
    raw_input("Press ENTER to exit...")
except:
    input("Press ENTER to exit...")