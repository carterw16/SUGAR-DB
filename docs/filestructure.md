# File Structure
* [classes](../classes) contains the steady-state models and any controls
* [ditto](../ditto) is an external library used to read distribution test system case files.
    It is a lightly modified version of the NREL Distribution Transformation Tool.
* [docs](../docs) is where all guidelines and style conventions are located
* [.gitlab](../.gitlab) is for gitlab templates and other configuration files.
* [input](../input) is for files that are used as external inputs into the SUGAR3 solver
* [lib](../lib) contains all functions used to support the SUGAR3 solver including the parser
* [output](../output) is where all results go after a SUGAR3 simulation completes. Results
    for testcases get their own folder.
* [testcases](../testcases) contains all Gridlab-D, OpenDSS, and CYME cases available for simulation.
* [validate](../validate) has the SUGAR3 validation scripts and unit tests.