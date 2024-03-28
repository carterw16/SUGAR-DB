from __future__ import absolute_import, division, print_function

import argparse
import logging
import logging.config as lcf
import os
import subprocess
from os import path
from pathlib import Path

import opendssdirect as dss

log_file_path = path.join(path.dirname(path.abspath(__file__)), 'logging.conf')
lcf.fileConfig(log_file_path)
# create logger
logger = logging.getLogger('runPowerflows')


def get_root_directory():
    if os.path.basename((os.path.abspath(os.getcwd()))) != 'SUGAR3':
        root_dir = os.path.abspath('.').replace(
            '{sep}{dirname}'.format(dirname=os.path.basename(
                os.path.abspath(os.getcwd())),
                                    sep=os.path.sep), '')
    else:
        root_dir = os.path.abspath('.')

    return root_dir


def run_opendss_powerflow(path_to_master, path_to_output):
    """
	Run opendss power flow on a given test case using the opendss python API.
	:param path_to_master:
	:type path_to_master:
	:param path_to_output:
	:type path_to_output:
	:return:
	:rtype:

	"""
    # run openDSS powerflow on case file
    result = dss.utils.run_command(
        'Redirect {}master.dss'.format(path_to_master))
    if result != '':
        logger.error('Unable to run openDSS powerflow on {}master.dss'.format(
            path_to_master))
        return
    else:
        logger.debug('opendss powerflow on {}master.dss successful'.format(
            path_to_master))
    # export openDSS voltage_profile
    dss.utils.run_command(
        'Export Voltages {}voltage_profile.csv'.format(path_to_output))


def run_gridlabd_powerflow(path_to_case, path_to_output):
    """
	Run gridlabd power flow on a given test case using gridlabd.
	:param path_to_case:
	:type path_to_case:
	:param path_to_output:
	:type path_to_output:
	:return:
	:rtype:
	"""

    # Search for gld file and ignore recorders, schedules, and other irrelevant files
    with os.scandir(path_to_case) as it:
        for entry in it:
            entry_name = entry.name.lower()
            if not entry_name.startswith(
                    'recorder') and not entry_name.startswith(
                        'schedule') and entry.is_file():
                case_file = entry.name

    # update case path with .glm file
    path_to_case = os.path.join(path_to_case, case_file)

    # run Gridlab-D powerflow on case file
    result = subprocess.run([
        'gridlabd', '-D',
        'outfile={}voltage_profile.csv'.format(path_to_output), '-q',
        path_to_case
    ],
                            capture_output=True,
                            encoding='utf-8',
                            universal_newlines=True)

    if result.returncode != 0:
        logger.error(
            'Unable to run gridlabd powerflow on {}'.format(path_to_case))
        return
    else:
        logger.debug(
            'Gridlab-D powerflow on {} successful'.format(path_to_case))


def run_cyme_powerflow():
    pass


def run_sugar_powerflow(path_to_case, path_to_output):
    """
	:param path_to_case:
	:type path_to_case:
	:param path_to_output:
	:type path_to_output:
	:return:
	:rtype:
	"""
    if path.basename((path.abspath(os.getcwd()))) != 'SUGAR3':
        root_dir = os.path.abspath('.').replace(
            '{sep}{dirname}'.format(dirname=path.basename(
                path.abspath(os.getcwd())),
                                    sep=os.path.sep), '')
    else:
        root_dir = os.path.abspath('.')

    path_to_SUGAR3 = os.path.abspath('{root}{sep}SUGAR3.py'.format(
        root=root_dir, sep=os.path.sep))
    path_to_settings = os.path.abspath(
        '{root}{sep}validate{sep}settings.json'.format(root=root_dir,
                                                       sep=os.path.sep))
    path_to_features = os.path.abspath(
        '{root}{sep}validate{sep}features.json'.format(root=root_dir,
                                                       sep=os.path.sep))

    result = subprocess.run([
        'python', path_to_SUGAR3, '--case', path_to_case, '--settings',
        path_to_settings, '--features', path_to_features, '--out',
        path_to_output
    ],
                            capture_output=True,
                            encoding='utf-8',
                            universal_newlines=True)
    if result.returncode != 0:
        logger.error(
            'Unable to run SUGAR3 powerflow on {}'.format(path_to_case))
    else:
        logger.debug('SUGAR3 powerflow on {} successful'.format(path_to_case))

    return result.returncode


def run_powerflows(case_names=None, powerflow_formats=None):
    """Runs opendss, gridlabd, cyme (not yet implemented), and SUGAR power flow analyses.
	This function is based on NREL DiTTo run_power_flow.py
	The input and output format is:
		- ./testcases/{format}/case_name
		- ./outputs/validation/case_name
	{format} refers to an opendss, gridlabd or cyme case
	**usage**: python run_powerflows -f ieee_4_node_stepdown_bal_Y-Y -f ieee_13node
	:return:
	:rtype:
	"""

    # If no cases are provided then run a set of standard cases
    if case_names is None:
        case_list = ['ieee_4node_stepdown_bal_Y-Y', 'ieee_13node']
    else:
        case_list = case_names

    # Possible formats for powerflow
    if powerflow_formats is None:
        formats_for_PF = ['opendss', 'gridlabd']
    else:
        formats_for_PF = powerflow_formats

    # Map power flow function to format name
    power_flow_functions = {
        'opendss': run_opendss_powerflow,
        'gridlabd': run_gridlabd_powerflow,
        'sugar': run_sugar_powerflow
    }

    root_dir = get_root_directory()

    returncode = 0
    for _format in formats_for_PF:
        for case in case_list:
            input_path = os.path.join(root_dir +
                                      '/testcases/{format_type}/{case}/'.format(
                                          format_type=_format, case=case))
            output_path = os.path.join(
                root_dir + '/output/validation/{format_type}/{case}/'.format(
                    format_type=_format, case=case))
            sugar_input_path = '{format_type}/{case}'.format(
                format_type=_format, case=case)
            sugar_output_path = os.path.join(
                root_dir + '/output/validation/sugar/{format_type}'.format(
                    format_type=_format))

            # if output path does not exist, create the path
            Path(output_path).mkdir(parents=True, exist_ok=True)
            Path(sugar_output_path).mkdir(parents=True, exist_ok=True)

            power_flow_functions[_format](input_path, output_path)
            returncode = power_flow_functions['sugar'](sugar_input_path,
                                                       sugar_output_path)

    return returncode


if __name__ == '__main__':
    # Parse the argumments
    cli_parser = argparse.ArgumentParser()

    # Case list
    cli_parser.add_argument('-c',
                            action='append',
                            dest='case_list',
                            default=None)
    # Format list
    cli_parser.add_argument(
        '-f',
        action='append',
        dest='format_list',
        default=None,
        help=
        'Provide a powerflow format (OpenDSS, Gridlabd-D) to compare to SUGAR.')

    parsed_args = cli_parser.parse_args()
    run_powerflows(parsed_args.case_list, parsed_args.format_list)
