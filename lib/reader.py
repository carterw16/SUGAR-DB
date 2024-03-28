"""Reads gridlabd or OpenDSS files.

Author(s): Naeem Turner-Bandele
Created Date: 06-30-2020
Updated Date: 10-18-2020
Email: nturnerb@cmu.edu
Status: Development

Reads GridLab-D or OpenDSS files using DiTTo's reader feature and then passes the output dictionaries to parser.py.

"""

import logging
import os
import time
import warnings

from scipy.sparse.linalg import MatrixRankWarning
from termcolor import colored

# Import the gridlabd and opendss readers
import ditto.readers.gridlabd.read as gld
import ditto.readers.opendss.read as dss
# Import dittos store module
from ditto.store import Store
from lib.get_root_directory import get_root_directory


def reader(case):
    """Reads the power systems case file.

        Reads the provided case path, searches for the power systems case file, and then uses DiTTo's reader to scan
        casefiles for the power systems elements, standardize the elements into a uniform format, and then outputs
        power grid elements as a dictionary.

        Args:
          case (string):
            A power systems test case.

        Returns:
          A dict with the power grid elements.

        """
    # initialize variables
    file_type = None
    case_file = None
    r = None

    # # File Types # #
    # 0 : gridlabd 1: opendss 2: CYME
    if 'gridlabd' in case:
        file_type = 0
    elif 'opendss' in case:
        file_type = 1
    elif 'cyme' in case:
        file_type = 2

    # Determine case path
    root_dir = get_root_directory()

    path_to_case = os.path.join(
        root_dir + os.path.sep + 'testcases' + os.path.sep,
        os.path.normpath(case))
    path_to_case = path_to_case.replace(
        '{seperator}validate'.format(seperator=os.path.sep), '')
    #################################
    #  STEP 1: CREATE AN EMPTY MODEL  #
    #################################
    logging.info(
        colored(
            '======================= Creating empty DiTTo model =======================',
            'white'))
    model = Store()
    start_reader_time = time.time()
    # Read file using gridlabd reader
    if file_type == 0:
        #################################
        #  STEP 2: READ FROM GRIDLAB-D  #
        #################################

        # Search for gld file and ignore recorders, schedules, and other irrelevant files
        with os.scandir(path_to_case) as it:
            for entry in it:
                entry_name = entry.name.lower()
                ignored = [
                    ".tmy", ".tmy2", "recorder", "schedule", "load", "loads"
                , "module", "player", ".csv",
                           ".json", ".txt", ".omd", "pv",  "equip", "vvo"]

                if not any(substr in entry_name for substr in ignored):
                    case_file = entry.name

        # Instantiate a Reader object
        r = gld.Reader(input_file=os.path.join(path_to_case, case_file))

        # Parse (can take some time for large systems...)
        logging.info(
            colored(
                '======================= Reading from Gridlab-D =======================',
                'white'))
    elif file_type == 1:
        ###############################
        #  STEP 2: READ FROM OPENDSS  #
        ###############################

        # Search for buscoords are in directory
        with os.scandir(path_to_case) as it:
            for entry in it:
                if 'master' in entry.name.lower() and entry.is_file():
                    case_file = entry.name

        # Instantiate a Reader object
        r = dss.Reader(master_file=os.path.join(path_to_case, case_file))

        # Parse (can take some time for large systems...)
        warnings.filterwarnings(
            "ignore",
            category=UserWarning)  # ignore the pandas not installed warning
        warnings.filterwarnings("always", category=MatrixRankWarning
                               )  # ignore the pandas not installed warning
        logging.info(
            colored(
                '======================= Reading from OpenDSS =======================',
                'white'))
    elif file_type == 2:
        ###############################
        #  STEP 2: READ FROM CYME  #
        ###############################
        pass
    else:
        logging.info("Invalid file type!")
        return
    end_reader_time = time.time()
    reader_time = end_reader_time - start_reader_time
    #################################
    #  STEP 3: PARSE MODEL  #
    #################################
    start_model_time = time.time()  # basic timer
    r.parse(model)
    # assign objects to casedata Dict based on key (k) and value(v)
    casedata = r.all_api_objects
    casedata["type"] = file_type
    end_model_time = time.time()
    modeler_time = end_model_time - start_model_time
    logging.info(
        colored(
            '================= DiTTo Model Parsed in {} seconds =================='
            .format(modeler_time), 'white'))
    simulation_stats = [reader_time, modeler_time]

    return casedata, simulation_stats
