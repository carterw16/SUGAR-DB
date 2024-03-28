"""Unit testing for the SUGAR3 code.

Author(s): Naeem Turner-Bandele
Created Date: 07-18-2020
Updated Date: 07-23-2020
Email: nturnerb@cmu.edu
Status: Validated

Runs a series of unit tests using OpenDSS and GridLab-D as the ground truth for comparison.

  Typical usage example:
  python validate.py
  make test
"""

import logging
import unittest

from compare_results import compare_results
from run_powerflows import run_powerflows


class TestPowerFlows(unittest.TestCase):

    def test_sugar_gridlabd(self):
        """Runs gridlabd unit test.

		 Executes powerflow on a specified test case in gridlabd and in SUGAR. After, compares the results of the
		 files in order to certify that the outputs match.

		 Args:
			  None

		 Returns:
			  Certifies if the unit test was passed by asserting true or false.

		 Raises:
			  Will raise an error for any issues that occur during SUGAR3 powerflow.

		  """
        case_name = ['ieee_13node']
        powerflow_format = ['gridlabd']
        powerflow_flag = run_powerflows(case_name, powerflow_format)
        if powerflow_flag == 0:
            validate_flag = compare_results(case_name,
                                            powerflow_format,
                                            append_file=True)
            self.assertTrue(validate_flag == 0)
            self.assertFalse(validate_flag != 0)
        elif powerflow_flag == 1:
            logging.info('SUGAR3 powerflow on case %s failed to converge.' %
                         case_name)
            self.assertFalse(True)
        else:
            logging.info(
                'SUGAR3 powerflow on case %s failed due to unknown error. Use runSUGAR3.py to debug.'
                % case_name)
            self.assertFalse(True)

    # def test_sugar_opendss(self):
    #     """Runs openDSS unit test.
    #
    #      Executes powerflow on a specified test case in OpenDSS and in SUGAR. After, compares the results of the
    #      files in order to certify that the outputs match.
    #
    #      Args:
    #           None
    #
    #      Returns:
    #           Certifies if the unit test was passed by asserting true or false.
    #
    #      Raises:
    #           Will raise an error for any issues that occur during SUGAR3 powerflow.
    #
    #       """
    #     case_name = ['ckt7']
    #     powerflow_format = ['opendss']
    #     powerflow_flag = run_powerflows(case_name, powerflow_format)
    #     if powerflow_flag == 0:
    #         validate_flag = compare_results(case_name,
    #                                         powerflow_format,
    #                                         append_file=True)
    #         self.assertTrue(validate_flag == 0)
    #         self.assertFalse(validate_flag != 0)
    #     elif powerflow_flag == 1:
    #         logging.info('SUGAR3 powerflow on case %s failed to converge.' %
    #                      case_name)
    #         self.assertFalse(True)
    #     else:
    #         logging.info(
    #             'SUGAR3 powerflow on case %s failed due to unknown error. Use runSUGAR3.py to debug.'
    #             % case_name)
    #         self.assertFalse(True)


if __name__ == '__main__':
    unittest.main()
