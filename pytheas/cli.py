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

Authors: Spingos I. & Kaviris G. (c) 2019-2021
Special thanks to Millas C. for testing the software and providing valuable feedback from the 
very early stages of this endeavor!

For any issues, comments or suggestions please contact us at ispingos@geol.uoa.gr or through GitHub 
at https://www.github.com/ispingos/pytheas-splitting

"""
#######################################################################
#                                                                     #
# This module includes functions/classes related to CLI functionality #
#                                                                     #
#######################################################################

import os
global _MODE_FILE
_MODE_FILE = __file__
# _MODE_FILE = 'lib/cli.py'
global _DIR_CWD
_DIR_CWD = os.path.dirname(os.path.realpath(_MODE_FILE))

## imports
import sys, logging, time
from logging.handlers import RotatingFileHandler
import pickle as pk
import numpy as np
from datetime import datetime
import time
import multiprocessing

#-- setup CCA logger
logger_cca_file = os.path.join(
    _DIR_CWD, 
    'logs', 
    'CCA_c%i_%s.log' % (int(sys.argv[1]), datetime.now().strftime('%Y%m%d_%H%M%S'))
    )
# logger_cca = ROOT_LOGGER.getLogger('logger_cca')
# logger_cca_fh = ROOT_LOGGER.FileHandler(logger_cca_file)
# logger_cca_fh.setLevel(ROOT_LOGGER.DEBUG)
# logger_cca_fh.setFormatter(ROOT_LOGGER.Formatter('%(asctime)s %(levelname)s %(message)s'))
# logger_cca.addHandler(logger_cca_fh)
# logger_cca_sh = ROOT_LOGGER.StreamHandler()
# logger_cca_sh.setLevel(ROOT_LOGGER.DEBUG)
# logger_cca_sh.setFormatter(ROOT_LOGGER.Formatter('%(asctime)s %(levelname)s %(message)s'))
# logger_cca.addHandler(logger_cca_sh)

# setup rotating files
core_str = '[p{:s}] '.format(sys.argv[1])
ROOT_LOGGER = logging.getLogger('my_root_logger')
ROOT_LOGGER.setLevel(logging.DEBUG)
root_handler = RotatingFileHandler(logger_cca_file, maxBytes=50 * 1000000,  # 50 MB
                                   backupCount=100)
root_handler.setFormatter(logging.Formatter(
     core_str + '%(asctime)s %(levelname)s %(message)s'
    ))
ROOT_LOGGER.addHandler(root_handler)
# setup logging to terminal
terminal = logging.StreamHandler()
terminal.setLevel(logging.DEBUG)
terminal.setFormatter(logging.Formatter(
    core_str +  '%(asctime)s %(levelname)s %(message)s'
))
ROOT_LOGGER.addHandler(terminal)

ROOT_LOGGER.info("Will save CCA related information to %s" % logger_cca_file)


# pytheas
from lib import eigenvalue as SC
from lib import clustering as CA
from lib import rotationcorrelation as RC
from lib import db_handler as DB
from lib import tools
from lib import parsers

#-- definitions
def clustering_station_wrapper(
        event, station, 
        dict_events, dict_stations, dict_picks_obs, 
        inventory, xml_inventory, vmod,
        cca_settings,
        tp_cnf, ca_cnf, general_cnf, grade_cnf, filter_cnf,
        path_db,
        failed_pairs, path_full_events
                    ):
    """
    """
    #-- connect to db
    db_conn, db_cur = DB.open(path_db)
    #--
    uid = "%s/%s/A%s" % (event, station, cca_settings.sws_method)
    ROOT_LOGGER.info('================================================')
    ROOT_LOGGER.info('%s: STARTING' % uid)
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
                ROOT_LOGGER.warning('%s: SKIP_NO_SPICK_f0' % uid)
                return
        #-- start tests for possible skip
        #- check for already analyzed
        if cca_settings.skip_pairs:
            _, vals = DB.retrieve(db_cur, "method", "method", uid)
            if vals:
                ROOT_LOGGER.warning('%s: SKIP_PREV_PROC' % uid)
                return
        #- check for previously failed attempt?
        if cca_settings.skip_failed:
            if uid in failed_pairs:
                ROOT_LOGGER.warning('%s: SKIP_PREV_FAILED' % uid)
                return
        #- check if station in selected
        if station not in cca_settings.selected_stations:
            ROOT_LOGGER.warning('%s: SKIP_NOT_SELECT_STATION' % uid)
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
                ROOT_LOGGER.exception("Could not calculate TauP ain for %s/%s" % (event, station)) 
                s_arr_taup = np.nan
        else:
            ain = dict_stations[station]['AN']
        
        if dict_stations[station]['AN'] > cca_settings.sws_window:
            ROOT_LOGGER.warning('%s: SKIP_AIN' % uid)
            return
        # check for taup pick eligibility
        if cca_settings.cca_pick_flag == 3:
            if tools.isnone(s_arr_taup) and tools.isnone(s_arr_obs):
                ROOT_LOGGER.warning('%s: SKIP_NO_SPICK_f3' % uid)
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
                ROOT_LOGGER.warning('%s: SKIP_NO_STREAM' % uid)
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
                ROOT_LOGGER.warning('%s: SKIP_NO_SYNC' % uid)
                return
        except:
            ROOT_LOGGER.exception('%s: error in reading the Stream from file %s' % (uid, path_full_pair))
            return
        #-- get araic pick
        if cca_settings.cca_pick_flag in [1, 2]:
            try:
                _, s_pick_auto = tools.arPicker(stream, pytheas.pkCNF)
                s_arr_auto = stream[0].stats.starttime + s_pick_auto
            except:
                ROOT_LOGGER.exception('Could not estimate AR-AIC arrival')
                s_pick_auto = np.nan
                s_arr_auto = np.nan
        #-- setup other picks too
        try:
            s_pick_obs = s_arr_obs - stream[0].stats.starttime
        except:
            ROOT_LOGGER.exception('Could not estimate observed pick')
            s_pick_obs = np.nan
        try:
            s_pick_taup = s_arr_taup - stream[0].stats.starttime
        except:
            ROOT_LOGGER.exception('Could not estimate TauP pick')
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
            ROOT_LOGGER.warning('%s: SKIP_NO_SPICK' % uid)
            return
        ROOT_LOGGER.info('Will use arrival/pick/type: %s / %.3f / %s' % (s_arr, s_pick, s_pick_type))
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
            ROOT_LOGGER.debug('Filtering in the %.2f / %.2f band' % (min(filter_bounds), max(filter_bounds)))
            stream.filter(
                type='bandpass',
                freqmin=min(filter_bounds),
                freqmax=max(filter_bounds),
                zerophase=True,
                corners=4
                        )
        else:
            ROOT_LOGGER.debug('No filter for this pair')
        #-- find the SNR and check
        snr = tools.calculate_snr(
                            stream, 
                            s_pick, 
                            general_cnf.snrStart, 
                            general_cnf.snrEnd
                            )
        ROOT_LOGGER.info('SNR is %.3f' % snr)
        if snr < cca_settings.snr_threshold:
            ROOT_LOGGER.warning('%s: SKIP_SNR' % uid)
            return
        #-- rotate to LQT or ZRT?
        baz = dict_stations[station]['BAZ']
        if cca_settings.rotate_to_lqt:
            ain = dict_stations[station]['AN']
            ROOT_LOGGER.info('Will rotate to LQT (ain: %.2f, baz: %.2f)' % (ain, baz))
        else:
            ain = 0.0
            ROOT_LOGGER.info('Will rotate to ZRT (baz: %.2f)' % baz)
        #-- add final method tag
        sws_method_tag = 'A%s' % cca_settings.sws_method
        ##############################################################
        #---- apply CA
        ROOT_LOGGER.info('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
        ROOT_LOGGER.info('Performing CA on %s ...' % uid)
        #-- calculate S-P
        try:
            tpts = s_arr - (dict_events['ORIGIN'] + dict_stations[station]['TOBSP'])
        except:
            ROOT_LOGGER.exception('Could not estimate S-P')
            tpts = np.inf
        if tpts > ca_cnf.tptsMax:
            tpts = ca_cnf.tptsMax
        ROOT_LOGGER.debug('Final S-P is %.4f s' % tpts)
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
        ROOT_LOGGER.debug('Dominant S period is %.4f' % period_dom)
        #-- finalize ca settings
        ca_cnf.Tbeg0 = -tpts / 2.
        #ca_cnf.Tend0 = dPer # + ca_cnf.Tend0  # maybe consider this for the future?
        ca_cnf.Tend1 = (ca_cnf.multPeriod * period_dom)  # + ca_cnf.Tbeg0
        ca_cnf.Nbeg = int(np.ceil((ca_cnf.Tbeg1 - ca_cnf.Tbeg0) / ca_cnf.DTbeg))
        ca_cnf.Nend = int(np.ceil((ca_cnf.Tend1 - ca_cnf.Tend0) / ca_cnf.DTend))
        ROOT_LOGGER.debug("CA beg windows: [%.4f , %.4f]" % (ca_cnf.Tbeg0,ca_cnf.Tbeg1))
        ROOT_LOGGER.debug("CA end windows: [%.4f , %.4f]" % (ca_cnf.Tend0,ca_cnf.Tend1))
        ROOT_LOGGER.debug("N beg/end windows: %i x %i" % (ca_cnf.Nbeg,ca_cnf.Nend))
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
        ROOT_LOGGER.info('%s: COMPLETED' % uid)

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
            ROOT_LOGGER.error('%s: FAIL_CA_THREAD_ERROR' % uid)
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
        ROOT_LOGGER.info('%s added to databse' % uid)
        ROOT_LOGGER.info('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')                
    except:
        ROOT_LOGGER.exception('%s: FAIL_GENERAL_ERROR' % uid)
        toc_pair = datetime.now()
        elapsed_time_pair_dt = toc_pair - tic_pair
        elapsed_time_pair_s = elapsed_time_pair_dt.seconds + elapsed_time_pair_dt.microseconds / 10 ** 6
        ROOT_LOGGER.info('%s: COMPLETED' % uid)
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
        ROOT_LOGGER.info('%s added to databse' % uid)
        ROOT_LOGGER.info('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')


#-- classes
class CLI():
    """The main CLI class"""

    def __init__(self, core_i):
        """initialization phase"""
        self.core_i = core_i
        ROOT_LOGGER.info('Running process #%i' % self.core_i)
        #-- load Pytheas configuration
        self.general_cnf = parsers.parseGeneralCnf(os.path.join(_DIR_CWD, 'etc', 'options', 'general.cnf'))
        self.pk_cnf = parsers.parsePickerCnf(os.path.join(_DIR_CWD, 'etc', 'options', 'picker.cnf'))
        self.ca_cnf = parsers.parseClusteringCnf(os.path.join(_DIR_CWD, 'etc', 'options', 'clustering.cnf'))
        self.tp_cnf = parsers.parseTaupCnf(os.path.join(_DIR_CWD, 'etc', 'options', 'taup.cnf'))
        self.grade_cnf = parsers.parseGradeCnf(os.path.join(_DIR_CWD, 'etc', 'options', 'grading.cnf'))
        self.filter_cnf = parsers.ParseFilterCnf(os.path.join(_DIR_CWD, 'etc', 'options', 'filters.cnf'))
        self.ca_cnf.maxTd = self.general_cnf.maxTd

        #-- load CCA settings
        self.cli_settings_dir = os.path.join(_DIR_CWD, 'etc', 'data_cli')
        self.cli_settings_path = os.path.join(self.cli_settings_dir, 'cca_settings_c%i.bin' % self.core_i)        

        try:
            self.cca_settings = pk.load(open(self.cli_settings_path, 'rb'))
            ROOT_LOGGER.info('Successfully read CCA parameters from file %s' % self.cli_settings_path)
        except FileNotFoundError:
            # self.cca_settings = tools.Dummy()  # maybe have defaults params in the future?
            ROOT_LOGGER.error('Could not read CCA parameters! Please make sure you made the required configurations through the GUI first.')
            return

        #-- load paths to data/catalog/db
        self.paths_settings_file = os.path.join(_DIR_CWD, 'etc', 'options', 'paths.exp')
        self.path_data, self.path_catalog, self.path_db = parsers.path_parser(self.paths_settings_file, mode='read')

        #-- ask user to run
        if self.cca_settings.cca_pick_flag == 0:
            self.ca_order_of_picks = ['Observed']
        elif self.cca_settings.cca_pick_flag == 1:
            self.ca_order_of_picks = ['Observed', 'AR-AIC', 'TauP']
        elif self.cca_settings.cca_pick_flag == 2:
            self.ca_order_of_picks = ['Observed', 'AR-AIC']
        elif self.cca_settings.cca_pick_flag == 3:
            self.ca_order_of_picks = ['Observed', 'TauP']
        user_msg = (
                '\n==============================================\n' + \
                'The following settings will be used in CCA:\n' + \
                'Data path: %s\n' % self.path_data + \
                'Catalog path: %s\n' % self.path_catalog + \
                'Database path: %s\n' % self.path_db + \
                'Method: %s\n' % self.cca_settings.sws_method + \
                '# of Stations: %s\n' % '{:,}'.format(len(self.cca_settings.selected_stations)) + \
                'Shear-wave Window: %.3f\n' % self.cca_settings.sws_window + \
                'Minimum SNR: %.3f\n' % self.cca_settings.snr_threshold + \
                'Filter: %s\n' % ('automatic\n' if self.cca_settings.filter_auto_flag else str(self.cca_settings.filter_band)) + \
                'Skip pairs already in DB: %s\n' % self.cca_settings.skip_pairs + \
                "Skip previously failed pairs: % s\n" % self.cca_settings.skip_failed + \
                'Order of picks: %s\n' % '->'.join(self.ca_order_of_picks)
                )
        # print(user_msg)
        # usr_ans = input('Do you want to continue (use s to show selected stations)? (y/n/s)\n')
        usr_ans = 'y'
        # usr_ans = 'y'
        if usr_ans.upper() in ['Y', 'YES', 'YEAH', 'YE']:
            pass
        elif usr_ans.upper() in ['S']:
            print('Selected stations are:\n%s\n' % ', '.join(self.cca_settings.selected_stations))
            usr_ans = input('Do you want to continue?? (y/n)\n')
            if usr_ans.upper() in ['Y', 'YES', 'YEAH', 'YE']:
                pass
            else:
                ROOT_LOGGER.warning('Aborting CCA processing per user request...')
                sys.exit()
        else:
            ROOT_LOGGER.warning('Aborting CCA processing per user request...')
            sys.exit()
        print('==============================================\n')

    def run(self):
        """"""
        #-- load station metadata
        try:
            self.inventory, self.xml_inventory = parsers.read_stations(self.tp_cnf.stations)
        except (OSError, FileNotFoundError):
            self.inventory, self.xml_inventory = parsers.read_stations(
                    os.path.join(os.path.split(_DIR_CWD)[0], self.tp_cnf.stations)
                    )
        #-- get event paths
        self.path_full_events, self.waveforms_num = tools.get_tree(self.path_data)
        self.events_num = len(self.path_full_events)
        #-- load event metadata from pytheas
        self.snapshot_pytheas = pk.load(
            open(
                os.path.join(self.cli_settings_dir, ('snapshot_pytheas_c%i.bin' % self.core_i)), 'rb')
            )
        
        #-- start
        tic = datetime.now()
        ROOT_LOGGER.info('Initiating CCA in CLI mode with the following settings')
        ROOT_LOGGER.info("Shear-wave Window: %.1f deg" % self.cca_settings.sws_window) 
        ROOT_LOGGER.info("Filter Band (Hz): %s" % ('automatic' if self.cca_settings.filter_auto_flag else str(self.cca_settings.filter_band)))
        ROOT_LOGGER.info('Order of picks is: %s' % '->'.join(self.ca_order_of_picks))
        ROOT_LOGGER.info("Process started @ %s" % tic)

        #-- prep
        # connect to the database
        self.db_conn, self.db_cur = DB.open(self.path_db)
        # find all failed attempts
        if self.cca_settings.skip_failed:
            _, uids = DB.find(self.db_cur, "method", "grade", "F")
            failed_pairs = ["%s/%s/%s" % x for x in uids]
            ROOT_LOGGER.debug('Found %s previously failed pairs to skip' % '{:,}'.format(len(failed_pairs)))
        else:
            failed_pairs = []
        # start iterating over events
        ROOT_LOGGER.debug('Starting iterating over events...')
        for i, event in enumerate(sorted(self.path_full_events)):
            ROOT_LOGGER.debug('Event code: %s (%s/%s)' % (
                                                        event, 
                                                        '{:,}'.format(i + 1),
                                                        '{:,}'.format(self.events_num)
                                                        ))
            #-- retrieve event/station/pick information
            try:
                dict_events = self.snapshot_pytheas['events'][event]
                dict_stations = self.snapshot_pytheas['stations'][event]
                dict_picks_obs = self.snapshot_pytheas['picks'][event]
                vmod = self.snapshot_pytheas['vmod']
            except KeyError:
                ROOT_LOGGER.exception('Could not find key of event %s in dicts' % event)
                continue
            ROOT_LOGGER.debug('Found %s stations & %s picks in metadata' % (
                        '{:,}'.format(len(dict_stations)),
                        '{:,}'.format(len(dict_picks_obs))
                                                                            ))
            #-- start iterating over stations
            ROOT_LOGGER.debug('Starting iterating over stations...')
            
            #-- setup multi
            # first, arguments list
            # m_args_list = ((
            #                 event, station, 
            #                 dict_events, dict_stations, dict_picks_obs, 
            #                 self.inventory, self.xml_inventory, vmod,
            #                 self.cca_settings,
            #                 self.tp_cnf, self.ca_cnf, self.general_cnf, self.grade_cnf,
            #                 self.path_db,
            #                 failed_pairs, self.path_full_events 
            #              ) for station in sorted(dict_stations) )
            # m_pool = multiprocessing.Pool(multiprocessing.cpu_count())
            # # m_pool = multiprocessing.pool.ThreadPool(multiprocessing.cpu_count())
            # m_pool.starmap(clustering_station_wrapper, m_args_list)
            # m_pool.close()
            # m_pool.join()
            
            for j, station in enumerate(sorted(dict_stations)):
                clustering_station_wrapper(
                            event, station, 
                            dict_events, dict_stations, dict_picks_obs, 
                            self.inventory, self.xml_inventory, vmod,
                            self.cca_settings,
                            self.tp_cnf, self.ca_cnf, self.general_cnf, self.grade_cnf,
                            self.filter_cnf,
                            self.path_db,
                            failed_pairs, self.path_full_events 
                                         )


        ROOT_LOGGER.info('Closing database...')
        DB.close(self.db_conn)   

        #-- show final stats
        toc = datetime.now()
        elapsed_time_process_dt = toc - tic
        elapsed_time_process_s = elapsed_time_process_dt.seconds + elapsed_time_process_dt.microseconds / 10 ** 6
        ROOT_LOGGER.info('Total elapsed time: %s' % time.strftime('%H:%M:%S', time.gmtime(elapsed_time_process_s)))
        ROOT_LOGGER.info('Finished processing catalogue file: %s' % self.path_catalog)

if __name__ == '__main__':
    cli = CLI(int(sys.argv[1]))
    cli.run()