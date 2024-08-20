#!/usr/bin/python
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

Authors: Spingos I. & Kaviris G. (c) 2019-2024
Special thanks to Millas C. for testing the software and providing valuable feedback from the 
very early stages of this endeavor!

For any issues, comments or suggestions please contact us at ispingos@geol.uoa.gr or through GitHub 
at https://www.github.com/ispingos/pytheas-splitting

"""
####################################################################
#                                                                  #
# This module contains the Clustering Thread                       #
#                                                                  #
####################################################################

## imports
import sys, logging, time
import pickle as pk
import numpy as np
from datetime import datetime
import time
import multiprocessing


#-- definitions
def clustering_station_wrapper(
        event, station, 
        dict_events, dict_stations, dict_picks_obs, 
        inventory, xml_inventory, vmod,
        cca_settings,
        tp_cnf, ca_cnf, general_cnf, grade_cnf,
        path_db,
        failed_pairs, path_full_events
                    ):
    """
    """
    #-- connect to db
    db_conn, db_cur = DB.open(path_db)
    #--
    uid = "%s/%s/A%s" % (event, station, cca_settings.sws_method)
    logging.info('================================================')
    logging.info('%s: STARTING' % uid)
    try:
        tic_pair = datetime.now()
        #- setup pair
        s_arr_taup = np.nan
        s_arr_auto = np.nan
        try:
            s_arr_obs = dict_picks_obs[station]
        except KeyError:
            s_arr_obs = np.nan
            if cca_settings.cca_pick_flag == 0:
                logging.warning('%s: SKIP_NO_SPICK_f0' % uid)
                return
        #-- start tests for possible skip
        #- check for already analyzed
        if cca_settings.skip_pairs:
            _, vals = DB.retrieve(db_cur, "method", "method", uid)
            if vals:
                logging.warning('%s: SKIP_PREV_PROC' % uid)
                return
        #- check for previously failed attempt?
        if cca_settings.skip_failed:
            if uid in failed_pairs:
                logging.warning('%s: SKIP_PREV_FAILED' % uid)
                return
        #- check if station in selected
        if station not in cca_settings.selected_stations:
            logging.warning('%s: SKIP_NOT_SELECT_STATION' % uid)
            return                                    
        #- check for incidence angle
        if tp_cnf.ainFlag:
            try:
                taupdat = tools.getTheorArrivals( # TODO: display secpndary phases too
                                        (
                                            dict_events["ORIGIN"],
                                            dict_events["LAT"],
                                            dict_events["LON"],
                                            dict_events["DEPTH"]
                                            ),
                                        (
                                            inventory[station]["latitude"],
                                            inventory[station]["longitude"]
                                            ),
                                        vmod
                                            )
                # grab measurements for the s arrival (not S)
                for reg in taupdat.values():
                    if reg["phase"] == "s": 
                        ain = reg["ain"]
                        dict_stations[station]['AN'] = ain
                        s_arr_taup = reg["time"]
                        break # can there be more than 1 "s" arrivals??
            except:
                logging.exception("Could not calculate TauP ain for %s/%s" % (event, station)) 
                s_arr_taup = np.nan
        else:
            ain = dict_stations[station]['AN']
        
        if dict_stations[station]['AN'] > cca_settings.max_ain:
            logging.warning('%s: SKIP_AIN' % uid)
            return
        # check for taup pick eligibility
        if cca_settings.cca_pick_flag == 3:
            if tools.isnone(s_arr_taup) and tools.isnone(s_arr_obs):
                logging.warning('%s: SKIP_NO_SPICK_f3' % uid)
                return
        #-- Read the waveforms for the selected event/station pair
        try:
            path_full_pair = os.path.join(
                path_full_events[event],
                '*.%s.*' % station
                    )
            stream, or_correction = tools.get_wave_record(
                                                path_full_pair,
                                                inventory=xml_inventory,
                                                pref_channels=general_cnf.chanPref,
                                                orient_corr=general_cnf.orientFlag
                                                )
            if not stream:
                logging.warning('%s: SKIP_NO_STREAM' % uid)
                return
            #-- trim?
            if general_cnf.trimFlag and s_arr_obs not in [stream[0].stats.starttime, 0, np.nan, None]:
                stream.trim(
                            starttime=s_arr_obs + general_cnf.trimStart,
                            endtime=s_arr_obs + general_cnf.trimEnd,
                            pad=True, nearest_sample=False, fill_value=0
                            )
            #-- sync waveforms in any case (trimmed or not)
            stream = tools.sync_waveforms(stream.copy())
            if abs(stream[1].stats.starttime - stream[2].stats.starttime) > stream[1].stats.delta:
                logging.warning('%s: SKIP_NO_SYNC' % uid)
                return
        except:
            logging.exception('%s: error in reading the Stream from file %s' % (uid, path_full_pair))
            return
        #-- get araic pick
        if cca_settings.cca_pick_flag in [1, 2]:
            try:
                _, s_pick_auto = tools.arPicker(stream, pytheas.pkCNF)
                s_arr_auto = stream[0].stats.starttime + s_pick_auto
            except:
                logging.exception('Could not estimate AR-AIC arrival')
                s_pick_auto = np.nan
                s_arr_auto = np.nan
        #-- setup other picks too
        try:
            s_pick_obs = s_arr_obs - stream[0].stats.starttime
        except:
            logging.exception('Could not estimate observed pick')
            s_pick_obs = np.nan
        try:
            s_pick_taup = s_arr_taup - stream[0].stats.starttime
        except:
            logging.exception('Could not estimate TauP pick')
            s_pick_taup = np.nan
        #-- assign picks and finish up
        s_pick_type = None
        s_pick = np.nan
        if cca_settings.cca_pick_flag == 0:  # obs-only
            s_pick = s_pick_obs
            s_pick_type = 'obs'
        elif cca_settings.cca_pick_flag == 1:  # obs->auto->taup
            for p, t in zip([s_pick_obs, s_pick_auto, s_pick_taup], ['obs', 'auto', 'taup']):
                if not tools.isnone(p):
                    s_pick = p
                    s_pick_type = t
                    break
        elif cca_settings.cca_pick_flag == 2:  # obs->auto
            for p, t in zip([s_pick_obs, s_pick_auto], ['obs', 'auto']):
                if not tools.isnone(p):
                    s_pick = p 
                    s_pick_type = t
                    break
        elif cca_settings.cca_pick_flag == 3:  # obs->taup
            for p, t in zip([s_pick_obs, s_pick_taup], ['obs', 'taup']):
                if not tools.isnone(p):
                    s_pick = p 
                    s_pick_type = t
                    break
        #-- finish up with the pick spaghetti
        try:
            s_arr = stream[0].stats.starttime + s_pick
        except ValueError:  # i.e. s_pick is np.nan
            s_arr = np.nan
            logging.warning('%s: SKIP_NO_SPICK' % uid)
            return
        logging.info('Will use arrival/pick/type: %s / %.3f / %s' % (s_arr, s_pick, s_pick_type))
        #-- filter the waveform (if req)
        if cca_settings.filter_auto_flag:
            b_idx, b_snr = tools.filtering_savage_et_al(
                            stream,
                            s_pick,
                            filter_cnf.filter_ranges,
                            general_cnf.snrStart, 
                            general_cnf.snrEnd
                                                        )
            filter_bounds = filter_cnf.filter_ranges[b_idx]
            # snr = b_snr  # redundant
        else:
            filter_bounds = cca_settings.filter_band
        if not np.isnan(filter_bounds[0]) and not np.isnan(filter_bounds[1]):
            logging.debug('Filtering in the %.2f / %.2f band' % (min(filter_bounds), max(filter_bounds)))
            stream.filter(
                type='bandpass',
                freqmin=min(filter_bounds),
                freqmax=max(filter_bounds),
                zerophase=True,
                corners=4
                        )
        else:
            logging.debug('No filter for this pair')
        #-- find the SNR and check
        snr = tools.calculate_snr(
                            stream, 
                            s_pick, 
                            general_cnf.snrStart, 
                            general_cnf.snrEnd
                            )
        logging.info('SNR is %.3f' % snr)
        if snr < cca_settings.min_snr:
            logging.warning('%s: SKIP_SNR' % uid)
            return
        #-- rotate to LQT or ZRT?
        baz = dict_stations[station]['BAZ']
        if cca_settings.rotate_to_lqt:
            ain = dict_stations[station]['AN']
            logging.info('Will rotate to LQT (ain: %.2f, baz: %.2f)' % (ain, baz))
        else:
            ain = 0.0
            logging.info('Will rotate to ZRT (baz: %.2f)' % baz)
        #-- add final method tag
        sws_method_tag = 'A%s' % cca_settings.sws_method
        ##############################################################
        #---- apply CA
        logging.info('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
        logging.info('Performing CA on %s ...' % uid)
        #-- calculate S-P
        try:
            tpts = s_arr - (dict_events['ORIGIN'] + dict_stations[station]['TOBSP'])
        except:
            logging.exception('Could not estimate S-P')
            tpts = np.inf
        if tpts > ca_cnf.tptsMax:
            tpts = ca_cnf.tptsMax
        logging.debug('Final S-P is %.4f s' % tpts)
        #-- get dominant period for shear-waves
        freq_dom_window = (
                int(np.floor(s_pick * stream[0].stats.sampling_rate)),
                int(np.ceil((s_pick + ca_cnf.specWindow) * stream[0].stats.sampling_rate))
                            )
        period_dom = 1. / tools.get_dom_freq(stream, freq_dom_window)
        # check whether period is within accepted bounds
        if period_dom < ca_cnf.minPeriod:
            period_dom = ca_cnf.minPeriod
        elif period_dom > ca_cnf.maxPeriod:
            period_dom = ca_cnf.maxPeriod
        logging.debug('Dominant S period is %.4f' % period_dom)
        #-- finalize ca settings
        ca_cnf.Tbeg0 = -tpts / 2.
        #ca_cnf.Tend0 = dPer # + ca_cnf.Tend0  # maybe consider this for the future?
        ca_cnf.Tend1 = (ca_cnf.multPeriod * period_dom)  # + ca_cnf.Tbeg0
        ca_cnf.Nbeg = int(np.ceil((ca_cnf.Tbeg1 - ca_cnf.Tbeg0) / ca_cnf.DTbeg))
        ca_cnf.Nend = int(np.ceil((ca_cnf.Tend1 - ca_cnf.Tend0) / ca_cnf.DTend))
        logging.debug("CA beg windows: [%.4f , %.4f]" % (ca_cnf.Tbeg0,ca_cnf.Tbeg1))
        logging.debug("CA end windows: [%.4f , %.4f]" % (ca_cnf.Tend0,ca_cnf.Tend1))
        logging.debug("N beg/end windows: %i x %i" % (ca_cnf.Nbeg,ca_cnf.Nend))
        #-- start the thread! our job here is done
        ca_thread = CA.clustering(
                            stream,
                            cca_settings.sws_method,
                            s_pick,
                            baz,
                            ain,
                            snr,
                            ca_cnf
                                        )
        ca_thread.run()
        # ca_thread.start()
        # while ca_thread.isRunning():
        #     time.sleep(0.1) # stability hack
        #     return
        toc_pair = datetime.now()
        elapsed_time_pair_dt = toc_pair - tic_pair
        elapsed_time_pair_s = elapsed_time_pair_dt.seconds + elapsed_time_pair_dt.microseconds / 10 ** 6
        logging.info('%s: COMPLETED' % uid)

        #-- make new splitting dictionary
        dict_splitting = {}
        dict_splitting[event] = tools.initSplittingDict("event")
        dict_splitting[event][station] = tools.initSplittingDict("station")
        dict_splitting[event][station][sws_method_tag] = tools.initSplittingDict("method")
        # update the event layer
        dict_splitting[event]["origin"] = dict_events["ORIGIN"].datetime
        dict_splitting[event]["latitude"] = float(dict_events["LAT"])
        dict_splitting[event]["longitude"] = float(dict_events["LON"])
        dict_splitting[event]["depth"] = float(dict_events["DEPTH"])
        dict_splitting[event]["magnitude"] = float(dict_events["MAG"])
        # update the station layer
        dict_splitting[event][station]["station"] = station
        dict_splitting[event][station]["network"] = dict_stations[station]["NET"]
        dict_splitting[event][station]["epicentral"] = float(dict_stations[station]["DIST"])
        dict_splitting[event][station]["azimuth"] = float(dict_stations[station]["BAZ"])
        dict_splitting[event][station]["incidence"] = float(dict_stations[station]["AN"])
        dict_splitting[event][station]["s_p"] = tpts
        dict_splitting[event][station]["orientation"] = or_correction
        # update arrivals and picks
        if not s_arr in [None, 0, np.nan, stream[0].stats.starttime]:
            dict_splitting[event][station]["s_obs"] = s_arr.datetime
        if not s_arr_taup in [None, 0, np.nan, stream[0].stats.starttime]:
            dict_splitting[event][station]["s_theo"] = s_arr_taup.datetime
        if not s_arr_auto in [None, 0, np.nan, stream[0].stats.starttime]:
            dict_splitting[event][station]["s_auto"] = s_arr_auto.datetime
        # update the method layer
        dict_splitting[event][station][sws_method_tag]["station"] = station
        dict_splitting[event][station][sws_method_tag]["method"] = sws_method_tag
        dict_splitting[event][station][sws_method_tag]["network"] = dict_stations[station]["NET"]
        dict_splitting[event][station][sws_method_tag]["SNR"] = snr
        dict_splitting[event][station][sws_method_tag]["s_freq"] = period_dom
        dict_splitting[event][station][sws_method_tag]["comment"] = 'elapsed time (s): %.3f' % elapsed_time_pair_s
        if ca_thread.excFlag:
            logging.error('%s: FAIL_CA_THREAD_ERROR' % uid)
            result_grade = 'F'
            dict_splitting[event][station][sws_method_tag]["grade"] = result_grade
        else:
            #-- store results
            result_phi = ca_thread.phi % 180
            result_td = ca_thread.dt
            result_errors = (ca_thread.sphi, ca_thread.sdt)
            result_snr = ca_thread.snr
            result_window_min, result_window_max = ca_thread.optWindow
            _, _, result_c_array, result_phi_test, result_td_test,\
            result_pol,_ ,_ ,result_c_e_max, \
            result_cc_fs, result_cc_ne, result_t_var, \
            result_n_contours = ca_thread.tempRes  
            result_grade_score, result_grade = tools.autoGrading(
                        (result_phi, result_pol, grade_cnf.polOff),
                        (result_errors[0], result_errors[1], result_snr, result_cc_fs, result_cc_ne),  # results
                        (grade_cnf.error_bounds[0],
                            grade_cnf.error_bounds[1],
                            grade_cnf.snr_bound,
                            grade_cnf.CC_FS_bound, grade_cnf.CC_NE_bound),  # bounds
                            grade_cnf.gradeDict
                                            )
            result_initial_clusters = ca_thread.initial
            result_initial_clusters_errors = ca_thread.initial_err
            result_calinski = ca_thread.calinski
            result_clusters1 = ca_thread.clusters1
            result_clusters2 = ca_thread.clusters2                                              
            #-- write results to dictionary
            dict_splitting[event][station][sws_method_tag]["phi"] = result_phi
            dict_splitting[event][station][sws_method_tag]["td"] = abs(result_td) * (10 ** 3)
            dict_splitting[event][station][sws_method_tag]["pol"] = result_pol % 180
            dict_splitting[event][station][sws_method_tag]["CC_FS"] = result_cc_fs
            dict_splitting[event][station][sws_method_tag]["CC_NE"] = result_cc_ne
            dict_splitting[event][station][sws_method_tag]["var_T"] = result_t_var
            dict_splitting[event][station][sws_method_tag]["err_phi"] = result_errors[0]
            dict_splitting[event][station][sws_method_tag]["err_td"] = result_errors[1] * (10 ** 3)
            dict_splitting[event][station][sws_method_tag]["N_contours"] = result_n_contours
            dict_splitting[event][station][sws_method_tag]["grade_score"] = result_grade_score
            dict_splitting[event][station][sws_method_tag]["grade"] = result_grade
            dict_splitting[event][station][sws_method_tag]["filter_min"] = min(filter_bounds)
            dict_splitting[event][station][sws_method_tag]["filter_max"] = max(filter_bounds)
            try:
                dict_splitting[event][station][sws_method_tag]["window_min"] = (stream[0].stats.starttime + result_window_min).datetime
                dict_splitting[event][station][sws_method_tag]["window_max"] = (stream[0].stats.starttime + result_window_max).datetime
            except:
                pass
            dict_splitting[event][station][sws_method_tag]["phi_test"] = np.asarray(result_phi_test)
            dict_splitting[event][station][sws_method_tag]["td_test"] = np.asarray(result_td_test)
            dict_splitting[event][station][sws_method_tag]["C_max"] = result_c_e_max
            dict_splitting[event][station][sws_method_tag]["C_array"] = np.asarray(result_c_array)
            dict_splitting[event][station][sws_method_tag]["initial_clusters"] = np.asarray(result_initial_clusters)
            dict_splitting[event][station][sws_method_tag]["initial_clusters_errors"] = np.asarray(result_initial_clusters_errors)
            dict_splitting[event][station][sws_method_tag]["calinski_score"] = np.asarray(result_calinski)
            dict_splitting[event][station][sws_method_tag]["clusters1"] = np.asarray(result_clusters1)
            dict_splitting[event][station][sws_method_tag]["clusters2"] = np.asarray(result_clusters2)    
        #-- store results to DB
        DB.addValues(db_conn, db_cur, dict_splitting)
        logging.info('%s added to databse' % uid)
        logging.info('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')                
    except:
        logging.exception('%s: FAIL_GENERAL_ERROR' % uid)
        toc_pair = datetime.now()
        elapsed_time_pair_dt = toc_pair - tic_pair
        elapsed_time_pair_s = elapsed_time_pair_dt.seconds + elapsed_time_pair_dt.microseconds / 10 ** 6
        logging.info('%s: COMPLETED' % uid)
        #-- make new splitting dictionary
        dict_splitting = {}
        dict_splitting[event] = tools.initSplittingDict("event")
        dict_splitting[event][station] = tools.initSplittingDict("station")
        dict_splitting[event][station][sws_method_tag] = tools.initSplittingDict("method")
        # update the event layer
        dict_splitting[event]["origin"] = dict_events["ORIGIN"].datetime
        dict_splitting[event]["latitude"] = float(dict_events["LAT"])
        dict_splitting[event]["longitude"] = float(dict_events["LON"])
        dict_splitting[event]["depth"] = float(dict_events["DEPTH"])
        dict_splitting[event]["magnitude"] = float(dict_events["MAG"])
        # update the station layer
        dict_splitting[event][station]["station"] = station
        dict_splitting[event][station]["network"] = dict_stations[station]["NET"]
        dict_splitting[event][station]["epicentral"] = float(dict_stations[station]["DIST"])
        dict_splitting[event][station]["azimuth"] = float(dict_stations[station]["BAZ"])
        dict_splitting[event][station]["incidence"] = float(dict_stations[station]["AN"])
        # update the method layer
        dict_splitting[event][station][sws_method_tag]["station"] = station
        dict_splitting[event][station][sws_method_tag]["method"] = sws_method_tag
        dict_splitting[event][station][sws_method_tag]["network"] = dict_stations[station]["NET"]
        dict_splitting[event][station][sws_method_tag]["comment"] = 'elapsed time (s): %.3f' % elapsed_time_pair_s
        dict_splitting[event][station][sws_method_tag]["grade"] = 'F'                    
        DB.addValues(db_conn, db_cur, dict_splitting)
        logging.info('%s added to databse' % uid)
        logging.info('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')