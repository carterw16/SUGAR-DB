import argparse
import csv
import os

import numpy as np
import pandas as pd
from natsort import natsorted


def get_root_directory():
    if os.path.basename((os.path.abspath(os.getcwd()))) != 'SUGAR3':
        root_dir = os.path.abspath('.').replace(
            '{sep}{dirname}'.format(dirname=os.path.basename(
                os.path.abspath(os.getcwd())),
                                    sep=os.path.sep), '')
    else:
        root_dir = os.path.abspath('.')

    return root_dir


def compare_results(case_names=None,
                    powerflow_formats=None,
                    voltmag_tol=0.01,
                    voltang_tol=0.05,
                    append_file=False):

    # If no cases are provided, then compare the standard set of cases
    if case_names is None:
        case_list = ['ieee_4node_stepdown_bal_Y-Y', 'ieee_13node']
    else:
        case_list = case_names

    # If no specific powerflow comparison, then look at both of the defaults
    if powerflow_formats is None:
        formats_for_PF = ['opendss', 'gridlabd']
    else:
        formats_for_PF = powerflow_formats

    # get the root directory path
    root_dir = get_root_directory()

    # Open and if necessary create a compare_results.csv to analyze case results
    file_param = 'a' if append_file else 'w'
    with open('voltage_profile_comparison.csv', file_param,
              newline='\n') as csvfile:
        result_writer = csv.writer(csvfile, delimiter=',')
        if not append_file:
            result_writer.writerow(
                ['VOLTAGE COMPARISON BETWEEN SUGAR AND A GROUND TRUTH'])
            result_writer.writerow([
                'Case', 'Ground Truth', 'Max Voltage Magnitude Difference (%)',
                'Max Voltage Angle Difference (%)',
                'Node Magnitude Maximum Difference',
                'Node Angle Maximum Difference'
            ])
        for _format in formats_for_PF:
            for case_name in case_list:
                path_to_comparison = root_dir + '/output/validation/{format_type}/{case}'.format(
                    format_type=_format, case=case_name)
                with os.scandir(path_to_comparison) as it:
                    for entry in it:
                        if 'voltage' in entry.name.lower() and entry.is_file():
                            comparison_file = entry.name

                path_to_sugar = root_dir + '/output/validation/sugar/{format_type}/{case}'.format(
                    format_type=_format, case=case_name)

                with os.scandir(path_to_sugar) as it:
                    for entry in it:
                        if 'voltage' in entry.name.lower() and entry.is_file():
                            sugar_file = entry.name

                # Load the sugar dataframe and drop unnecessary columns
                sugar_df = pd.read_csv(os.path.join(path_to_sugar, sugar_file),
                                       skiprows=1)
                sugar_df.drop(columns=[
                    ' BasekV', ' phase1', ' phase2', ' phase3', ' pu1', ' pu2',
                    ' pu3', ' Unbalance (%)'
                ],
                              inplace=True)

                # Load the comparison dataframe, perform filtering, shifting if necessary, drop unnecessary columns,
                # and merge the comparison dataframe with the sugar results dataframe.
                if _format == 'opendss':
                    comparison_df = pd.read_csv(
                        os.path.join(path_to_comparison, comparison_file))
                    comparison_df.sort_values(by=['Bus'], inplace=True)

                    # create filters that will be used to shift columns
                    node1_filter2 = (comparison_df.iloc[:, 2] == 2
                                    )  # find phase 2 in Node1
                    node1_filter3 = (comparison_df.iloc[:, 2] == 3
                                    )  # find phase 3 in Node1

                    # find the rows that correspond to the columns which need shifting
                    misplaced_node1_phase3_idx = node1_filter3[
                        node1_filter3].index
                    misplaced_node1_phase2_idx = node1_filter2[
                        node1_filter2].index

                    # shift the rows
                    comparison_df.loc[misplaced_node1_phase2_idx,
                                      2:] = comparison_df.loc[
                                          misplaced_node1_phase2_idx,
                                          ' Node1':].shift(3,
                                                           axis=1,
                                                           fill_value=0)
                    comparison_df.loc[misplaced_node1_phase3_idx,
                                      2:] = comparison_df.loc[
                                          misplaced_node1_phase3_idx,
                                          ' Node1':].shift(6,
                                                           axis=1,
                                                           fill_value=0)

                    # after initial shifting, check if Node2 needs shifting
                    node2_filter3 = (comparison_df.iloc[:, 6] == 3
                                    )  # find phase 3 in Node2
                    misplaced_node2_phase3_idx = node2_filter3[
                        node2_filter3].index
                    comparison_df.loc[misplaced_node2_phase3_idx,
                                      6:] = comparison_df.loc[
                                          misplaced_node2_phase3_idx,
                                          ' Node2':].shift(3,
                                                           axis=1,
                                                           fill_value=0)
                    comparison_df.drop(columns=[
                        ' BasekV', ' Node1', ' Node2', ' Node3', ' pu1', ' pu2',
                        ' pu3'
                    ],
                                       inplace=True)
                    comparison_df['Bus'] = comparison_df['Bus'].str.lower()

                    # merge dataframes
                    Result_Merge = pd.merge(sugar_df,
                                            comparison_df,
                                            how='inner',
                                            on=['Bus'])

                else:
                    comparison_df = pd.read_csv(os.path.join(
                        path_to_comparison, comparison_file),
                                                skiprows=1)
                    comparison_df.node_name = comparison_df.node_name.astype(
                        'category')
                    comparison_df.node_name.cat.reorder_categories(natsorted(
                        comparison_df.node_name),
                                                                   inplace=True,
                                                                   ordered=True)
                    comparison_df.sort_values('node_name', inplace=True)
                    volt_angles = ['voltA_angle', 'voltB_angle', 'voltC_angle']
                    comparison_df[volt_angles] = comparison_df[
                        volt_angles].apply(lambda row: np.rad2deg(row))
                    sugar_df.rename(columns={"Bus": "node_name"}, inplace=True)

                    # Merge the two data frames
                    Result_Merge = pd.merge(sugar_df,
                                            comparison_df,
                                            how='inner',
                                            on=['node_name'])

                # Change negative angles to positive
                Result_Merge = Result_Merge.apply(
                    lambda x: x % 360 if x.name in list(Result_Merge.columns[
                        [2, 4, 6, 8, 10, 12]]) else x,
                    result_type='broadcast')

                # convert to matrix
                Result = Result_Merge.to_numpy()

                # Calculate differences
                sugar_results = Result[:, 1:7]
                comparison_results = Result[:, 7:]
                Result_Diff = np.abs(sugar_results - comparison_results)

                # Find the maximum error
                max_error = Result_Diff.max(axis=0)
                max_index = Result_Diff.argmax(axis=0)

                # Calculate

                # Find the relative error
                try:
                    max_VoltA_pu = max_error[0] / Result[max_index[0], 7]
                except ZeroDivisionError:
                    max_VoltA_pu = np.nan
                try:
                    max_VoltB_pu = max_error[2] / Result[max_index[2], 9]
                except ZeroDivisionError:
                    max_VoltB_pu = np.nan
                try:
                    max_VoltC_pu = max_error[4] / Result[max_index[4], 11]
                except ZeroDivisionError:
                    max_VoltC_pu = np.nan
                try:
                    max_angA_pu = max_error[1] / Result[max_index[1], 8]
                except ZeroDivisionError:
                    max_angA_pu = 0
                try:
                    max_angB_pu = max_error[3] / Result[max_index[3], 10]
                except ZeroDivisionError:
                    max_angB_pu = 0
                try:
                    max_angC_pu = max_error[5] / Result[max_index[5], 12]
                except:
                    max_angC_pu = 0

                max_voltmag_bus = np.argmax(
                    [max_error[0], max_error[2], max_error[4]])
                max_voltang_bus = np.argmax(
                    [max_error[1], max_error[3], max_error[5]])
                if _format == 'opendss':
                    max_err_V_node = Result_Merge.Bus[
                        max_index[max_voltmag_bus]]
                    max_err_ang_node = Result_Merge.Bus[
                        max_index[max_voltang_bus]]
                else:
                    max_err_V_node = Result_Merge.node_name[
                        max_index[max_voltmag_bus]]
                    max_err_ang_node = Result_Merge.node_name[
                        max_index[max_voltang_bus]]
                max_voltmag_pu = max(max_VoltA_pu, max_VoltB_pu, max_VoltC_pu)
                max_voltang_pu = max(max_angA_pu, max_angB_pu, max_angC_pu)
                max_voltmag_percent = np.around(max_voltmag_pu * 100, 5)
                max_voltang_percent = np.around(max_voltang_pu * 100, 5)
                # Result Writer
                result_writer.writerow([
                    case_name, _format, max_voltmag_percent,
                    max_voltang_percent, max_err_V_node, max_err_ang_node
                ])

    # save return code which is used only for the unit tests in validate.py
    if _format == 'opendss':
        return 0 if max_voltmag_pu <= voltmag_tol else 1
    else:
        return 0 if max_voltmag_pu <= voltmag_tol and max_voltang_pu <= voltang_tol else 1


if __name__ == '__main__':
    cli_parser = argparse.ArgumentParser(
        description='Process compare_results inputs.')
    cli_parser.add_argument('-c',
                            action='append',
                            dest='case_list',
                            default=None,
                            help='Provide a test case to compare to SUGAR.')
    cli_parser.add_argument(
        '-f',
        action='append',
        dest='format_list',
        default=None,
        help=
        'Provide a powerflow format (OpenDSS, Gridlabd-D) to compare to SUGAR.')
    cli_parser.add_argument('-mag_tol',
                            dest='voltmag_tol',
                            type=float,
                            default=0.01,
                            help='The acceptable voltage magnitude difference.')
    cli_parser.add_argument('-ang_tol',
                            dest='voltang_tol',
                            type=float,
                            default=0.05,
                            help='The acceptable voltage angle difference.')
    cli_parser.add_argument(
        '-a',
        dest='append_file',
        type=bool,
        default=False,
        help=
        'Decide if you would like the results written or appended with True or False.'
    )
    args = cli_parser.parse_args()
    compare_results(args.case_list, args.format_list, args.voltmag_tol,
                    args.voltang_tol, args.append_file)
