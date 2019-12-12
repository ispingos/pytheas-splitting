# -*- coding: utf-8 -*-

"""
================================================================================
Pytheas: an open-source software solution for local shear-wave splitting studies
================================================================================

Pytheas is a tool that aims to introduce a new mentality in shear-wave splitting analysis 
from local recordings, incorporating manual, semi- and fully- automatic methods under one 
Graphical User Interface. Integrating databases and offering compatibility with popular data 
and metadata file formats, Pytheas is designed with the simplification of analysis in mind, 
while, at the same time, enhanching the effectiveness of processing and quality control of results.

Pytheas is released under the GNU GPLv3 license.

Authors: Spingos I. & Kaviris G. (c) 2019
Special thanks to Millas C. for testing the software and providing valuable feedback from the 
very early stages of this endeavor!

For any issues, comments or suggestions please contact us at ispingos@geol.uoa.gr or through GitHub 
at https://www.github.com/ispingos/pytheas-splitting

"""
####################################################################
#                                                                  #
# This module includes Database related functions.                 #
#                                                                  #
####################################################################

## imports
import os, logging
import sqlite3, io
from obspy import UTCDateTime
import numpy as np

## global variables ##
# default keys for each layer
defKeys={
            "event":("origin","latitude","longitude","depth","magnitude"),
            "station":("station","network","epicentral","azimuth","incidence","s_p",
                       "s_obs","s_theo","s_auto","orientation"),
            "method":("station","network","method","phi","td","pol","CC_NE","CC_FS","var_T","err_phi","err_td",
                      "N_contours","grade_score","grade","filter_min","filter_max",
                      "window_min","window_max","SNR","s_freq","comment",
                      "C_array","C_max","phi_test","td_test", # EV/RC related
                      "initial_clusters","initial_clusters_errors","calinski_score",
                      "clusters1","clusters2" # cluster analysis related
                       )
        }
# assign the SQL datatype related to each key
typeKeys={"origin":"DATETIME","latitude":"REAL","longitude":"REAL","depth":"REAL","magnitude":"REAL",
          "station":"TEXT","network":"TEXT","epicentral":"REAL","azimuth":"REAL","incidence":"REAL","SNR":"REAL",
          "s_obs":"DATETIME","s_theo":"DATETIME","s_auto":"DATETIME",
          "s_freq":"REAL","s_p":"REAL","orientation":"REAL",
          "method":"TEXT","phi":"REAL","td":"REAL","pol":"REAL","CC_NE":"REAL","CC_FS":"REAL",
          "var_T":"REAL","err_phi":"REAL",
          "err_td":"REAL","N_contours":"INTEGER","grade":"TEXT","grade_score":"REAL",
          "filter_min":"REAL","filter_max":"REAL",
          "window_min":"DATETIME","window_max":"DATETIME","comment":"TEXT",
          "C_array":"ARRAY","C_max":"REAL","phi_test":"ARRAY","td_test":"ARRAY", # EV/RC related
          "initial_clusters":"ARRAY","initial_clusters_errors":"ARRAY","calinski_score":"ARRAY",
          "clusters1":"ARRAY","clusters2":"ARRAY" # cluster analysis related
          }

## definitions ##
def adapt_array(arr):
    """
    Adapter function modified from 
    https://www.pythonforthelab.com/blog/storing-data-with-sqlite/#storing-numpy-arrays-into-databases
    
    :type arr: array-like
    :param arr: a numpy-array
    :returns: converted array to bytecode

    """
    out=io.BytesIO()
    np.save(out, arr)
    out.seek(0)
    return out.read()

def convert_array(text):
    """
    Converter function modified from
    https://www.pythonforthelab.com/blog/storing-data-with-sqlite/#storing-numpy-arrays-into-databases

    :type text: str
    :param text: bytecode string to convert to array
    :returns: the array obtained from the initial bytecode

    """
    out=io.BytesIO(text)
    out.seek(0)
    return np.load(out, allow_pickle=True)

def init(dbname):
    """
    initializes the database based on the
    provided catalogue and db output directory

    :type dbname: str
    :param dbname: path to the database file.
    :returns: (database connection, database cursor)

    """
    # check if db exists
    if os.path.exists(dbname):
        logging.info("%s already exists" % dbname)
        raise IOError("%s already exists" % dbname)
    # start it up
    logging.info("initializing db at: %s" % dbname)
    # arrays<->sqlite3
    sqlite3.register_adapter(np.ndarray,adapt_array)
    sqlite3.register_converter("ARRAY",convert_array)
    conn=sqlite3.connect(dbname,detect_types=sqlite3.PARSE_DECLTYPES)
    cur=conn.cursor()
    # make event table
    sql_com="CREATE TABLE event (evid TEXT, %s)" % ", ".join([x+" "+typeKeys[x] for x in defKeys['event']])
    cur.execute(sql_com)
    # make station table
    sql_com="CREATE TABLE station (evid TEXT, %s)" % ", ".join([x+" "+typeKeys[x] for x in defKeys['station']])
    cur.execute(sql_com)
    # make method table
    sql_com="CREATE TABLE method (evid TEXT, %s)" % ", ".join([x+" "+typeKeys[x] for x in defKeys['method']])
    cur.execute(sql_com)
    # finalize
    conn.commit()
    return conn,cur

def open(dbname):
    """
    connects to given database and creates the relevant cursor
    
    :type dbname: str
    :param dbname: path to the database file.
    :returns: (database connection, database cursor)

    """
    # arrays<->sqlite3
    sqlite3.register_adapter(type(np.zeros((1,1))), adapt_array)
    sqlite3.register_converter("ARRAY", convert_array)
    #
    logging.info("connecting to db at: %s" % dbname)
    conn=sqlite3.connect(dbname,detect_types=sqlite3.PARSE_DECLTYPES)
    cur=conn.cursor()  
    return conn,cur

def close(conn):
    """Close the connection to the database."""
    conn.close()

def addValues(conn,cur,splittingDict):
    """
    Add new values to the database.

    :type conn: :class: `~sqlite3.Connection`
    :param conn: the connection object, linked to the
        database.
    :type cur: :class: `~sqlite3.Cursor`
    :param cur: the cursor instance, linked to the
        database
    :type splittingDict: dict
    :param splittingDict: a dict containing all the splitting
        information to be added to the database.

    """
    ## get existing values
    # event ids
    selcom='SELECT evid FROM event'
    cur.execute(selcom)
    evids=cur.fetchall()    
    # station ids
    selcom='SELECT evid,station FROM station'
    cur.execute(selcom)
    stids=cur.fetchall()
    # method ids
    selcom='SELECT evid,station,method FROM method'
    cur.execute(selcom)
    mids=cur.fetchall()
    # add values for all three layers
    for evid in splittingDict:
        # update the event table
        eveDict=splittingDict[evid]
        if (evid,) not in evids:
            ecom='INSERT INTO event (%s) values (?' % ', '.join(["evid"]+[x for x in sorted(defKeys['event'])])
            ecom+=',?'*(len(defKeys['event']))+')'
            cur.execute(ecom,[evid]+[eveDict[x] for x in sorted(defKeys['event'])])
        else:
            ecom='UPDATE event SET %s WHERE evid=?' % \
                (
                  ', '.join(['%s = ?' % x for x in sorted(defKeys['event'])])
                )
            cur.execute(ecom,[eveDict[x] for x in sorted(defKeys['event'])]+[evid])
        # update the station table
        for skey in sorted(splittingDict[evid]):
            if skey in defKeys['event']: continue
            staDict=splittingDict[evid][skey]
            if (evid,skey) not in stids:
                scom='INSERT INTO station (%s) values (?' % ', '.join(["evid"]+[x for x in sorted(defKeys['station'])])
                scom+=',?'*(len(defKeys['station']))+')'
                cur.execute(scom,[evid]+[staDict[x] for x in sorted(defKeys['station'])])
            else:
                scom='UPDATE station SET %s WHERE evid=? AND station=?' % \
                    (
                      ', '.join(['%s = ?' % x for x in sorted(defKeys['station'])])
                    )
                cur.execute(scom,[staDict[x] for x in sorted(defKeys['station'])]+[evid,skey])
            # update the method table
            for mkey in sorted(staDict):
                if mkey in defKeys['station']: continue
                metDict=staDict[mkey]
                if metDict['grade'] in ['X']:
                    continue # skip X grades
                if (evid,skey,mkey) not in mids:
                    mcom='INSERT INTO method (%s) values (?' % ', '.join(["evid"]+[x for x in sorted(defKeys['method'])])
                    mcom+=',?'*(len(defKeys['method']))+')'
                    cur.execute(mcom,[evid]+[metDict[x] for x in sorted(defKeys['method'])])
                else:
                    mcom='UPDATE method SET %s WHERE evid=? AND station=? AND method=?' % \
                        (
                          ', '.join(['%s = ?' % x for x in sorted(defKeys['method'])])
                        )
                    cur.execute(mcom,[metDict[x] for x in sorted(defKeys['method'])]+[evid,skey,mkey])
    conn.commit() 
    logging.info("Successfully updated database.")

def find(cur,table,column,value):
    """
    Find event-station-method ids for the specified value in
    the specified column.

    :type cur: :class: `~sqlite3.Cursor`
    :param cur: the cursor instance, linked to the
        database.
    :type table: str
    :param table: the table name ('event','station' or 'method').
    :type column: str
    :param column: the column name.
    :type value: float, int or str
    :param value: the value to look for.
    :returns: (headers, values) from the retrieved data.

    """
    if table == "event":
        scom='SELECT evid FROM %s WHERE %s == ?' % (table,column)
    elif table == "station":
        scom='SELECT evid,station FROM %s WHERE %s == ?' % (table,column)
    elif table == "method":
        scom="SELECT evid,station,method FROM %s WHERE %s == ?" % (table,column)
    args=(value,)
    cur.execute(scom,args)
    hdr=[description[0] for description in cur.description]
    val=cur.fetchall()
    return hdr,val

def retrieve(cur,table,column,uid):
    """
    retrieve a specific value from a specific table. the uid is a string in the form:
    event_id/station/method

    :type cur: :class: `~sqlite3.Cursor`
    :param cur: the cursor instance, linked to the
        database.
    :type table: str
    :param table: the table name ('event','station' or 'method').
    :type column: str
    :param column: the column name.
    :type uid: str
    :param uid: a combination of ids in the form of 'event/station/method'.
    :returns: (headers, values) from the retrieved data.

    """
    uid=uid.replace("*",'%')
    uidList=uid.split("/")
    if len(uidList) == 1:
        evid=uidList[0]
        scom='SELECT %s FROM %s WHERE evid LIKE ?' % (column,table)
        args=(evid,)
    elif len(uidList) == 2:
        evid,station=uidList
        scom='SELECT %s FROM %s WHERE evid LIKE ? AND station LIKE ?' % (column,table)
        args=(evid,station)
    elif len(uidList) == 3:        
        evid,station,method=uidList
        scom='SELECT %s FROM %s WHERE evid LIKE ? AND station LIKE ? AND method LIKE ?' % (column,table)
        args=(evid,station,method)
    cur.execute(scom,args)
    hdr=[description[0] for description in cur.description]
    val=cur.fetchall()
    return hdr,val

def load(cur,evids="*",stids="*",meids="*"):
    """
    Load database into a splitting dictionary.

    :type cur: :class: `~sqlite3.Cursor`
    :param cur: the cursor instance, linked to the
        database.
    :type evids: list or str, optional
    :param eivds: list of event ids to load from the database, If '*'
        is specified, all events will be loaded. Defaults to '*'.
    :type stids: list or str, optional
    :param stids: list of station ids to load from the database, If '*'
        is specified, all stations will be loaded. Defaults to '*'.
    :type meids: list or str, optional
    :param meids: list of method ids to load from the database, If '*'
        is specified, all methods will be loaded. Defaults to '*'.        
    :returns: a dict containing all splitting information for the query.
    
    """
    # make check list
    chkList=(evids=="*",stids=="*",meids=="*")
    # set meids to all available
    if meids == "*":
        meids=["MAN","MEV","MME","MRC","AEV","AME","ARC"]
    # get event data
    if any(chkList):
        cur.execute('SELECT * FROM event')
    else:
        cur.execute(
                    'SELECT * FROM event WHERE evid IN (%s) AND evid IN \
                     (SELECT evid FROM method)' \
                    % (
                        ','.join(['"%s"' % x for x in evids])
                      )
                   )
    event_hdr=[description[0] for description in cur.description]
    event=cur.fetchall()
    if evids == "*" and not all((chkList[1],chkList[2])):
        evids=[x[event_hdr.index("evid")] for x in event]
    # get station data
    if any(chkList):
        cur.execute('SELECT * FROM station')
    elif evids != "*":
        cur.execute(
                    'SELECT * FROM station WHERE evid IN (%s) AND evid IN \
                     (SELECT evid FROM method)' \
                    % (
                        ','.join(['"%s"' % x for x in evids])
                      )
                   )    
    else:
        cur.execute('SELECT * FROM station WHERE evid IN (%s) AND station IN (%s) AND evid IN \
                     (SELECT evid FROM method)' \
                     % (
                        ','.join(['"%s"' % x for x in evids]),
                        ','.join(['"%s"' % x for x in stids])
                       )
                   )
    station_hdr=[description[0] for description in cur.description]
    station=cur.fetchall()
    if stids == "*" and not all((chkList[0],chkList[2])):
        stids=[x[station_hdr.index("station")] for x in station]
    # get method data
    if all(chkList):
        cur.execute('SELECT * FROM method')
    else:
        cur.execute('SELECT * FROM method WHERE evid IN (%s) AND station IN (%s) AND method IN (%s)' \
            % (
               ','.join(['"%s"' % x for x in evids]),
               ','.join(['"%s"' % x for x in stids]),
               ','.join(['"%s"' % x for x in meids])
               )
                    )
    method_hdr=[description[0] for description in cur.description]
    method=cur.fetchall()
    ##
    splittingDict={}
    for event_dat in event:
        evid=event_dat[event_hdr.index("evid")]
        splittingDict[evid]={x:y for x,y in zip(event_hdr,event_dat) if x != "evid"}
        splittingDict[evid]['origin']=UTCDateTime(splittingDict[evid]['origin']).datetime
    for station_dat in station:
        evid=station_dat[station_hdr.index("evid")]
        if (evids != "*") and (evid not in evids):
            continue
        stid=station_dat[station_hdr.index("station")]
        if (stids != "*") and (stid not in stids):
            continue
        splittingDict[evid][stid]={x:y for x,y in zip(station_hdr,station_dat) if x != "evid"}
        # correct np nans
        for x,y in splittingDict[evid][stid].items():
            if y == None:
                splittingDict[evid][stid][x]=np.nan
    for method_dat in method:
        evid=method_dat[method_hdr.index("evid")]
        stid=method_dat[method_hdr.index("station")]
        meid=method_dat[method_hdr.index("method")]
        try:
            splittingDict[evid][stid][meid]={x:y for x,y in zip(method_hdr,method_dat) if x != "evid"}
        except KeyError:
            logging.exception("Key missing")
            continue
        # correct np nans
        for x,y in splittingDict[evid][stid][meid].items():
            if type(y) is np.ndarray: 
                continue
            if y == None:
                splittingDict[evid][stid][meid][x]=np.nan
    ###
    logging.info("Found a total of %i events in database" % len(splittingDict))
    return splittingDict