# SUGAR3
Three-Phase Power Flow + Optimization Solver. Designed for distribution systems but capable of solving transmission system cases.

## Steps to Run SUGAR3:
### 1. Install Requirements
    make init
### 2. Compile Files
    make compile
  You only need to run this once, unless changes are made to a Cython (*.pyx) file then you should rerun it.

### 3. Run SUGAR3
There are two ways to execute SUGAR3.

#### Via [runSUGAR3.py](runSUGAR3.py)
 * Open [runSUGAR3.py](runSUGAR3.py):
 * Select desired test case from testcases and format per instructions in [Test Cases](#test-cases)
 * Turn settings on/off
 * Turn features on/off
 * Add generator setpoints to "gen_setpoints.csv" if IBDGs enabled
 * Execute`python runSUGAR3.py`

### Via Command Line:
If you want to use the default settings run:
 
    python SUGAR3.py --case format/caseName

If settings, features, and the output location are important then create a [settings.json](/input/settings.json) and
 [features.json](/input/features.json
) file and execute:
```
python SUGAR3.py --case format/caseName --settings input/settings.json --features input/features.json --out
 [insert path_to_output here]
```

### 4. Check Outputs
If outputs are saved, they are in the [output](output) folder by default.

## Requirements
See [requirements.txt](requirements.txt) for a list of necessary packages.

## Developer Guide:
If you are contributing to this repository, the links below contain vital information.

* [File Structure](docs/filestructure.md)
* [Style Guide](docs/styleguide.md)
* [Contributing Guide](docs/contributing.md)


## Tests
Use the following command to run tests before committing any new code. Always run tests.
```
bash test_infeas_all_cases
```
This bash script will test your code against Gridlab-D results for six specific test cases. The percent error is below 0.2% between the present code base and the Gridlab-D solutions. If you commit any new code, please ensure that it does not cause additional error.

## Feature Status
  A list of features and their status is [here](docs/feature_status.md). Within issues, there is a list of feature requests that have not yet been implemented or are in progress.
    
## Test Cases:
  Test cases are located in [testcases](testcases). All cases are sorted by simulator type (e.g., OpenDSS, GridLab-D or
   CYME). To run a case, simply specify the format (e.g. gridlabd or opendss) and the name of the case. For example
   , to run a GridLab-D case:
   ```
       python SUGAR3.py --case gridlabd/ieee_8500node
   ``` 
   To run an OpenDSS case:
   ```
      python SUGAR3.py --case opendss/ieee_34node
   ```

   A detailed accounting of the status of all testcases is [here](testcases/testcasestatus.xlsx).
 
   #### Transmission System Cases 
   * Transmission test cases marked with "_T"
   * Transmission cases need their load factor (LF) reduced
   * A load factor (LF) between 0.6-0.1 usually works but varies by case

  ## Code Citation
  A. Pandey, N. Turner-Bandele, and E. Foster, _SUGAR-D: A Distribution Systems Analysis and Optimization Tool v2.02_, CMU ECE Pileggi Lab Group, Pittsburgh, Pennsylvania, 2021. 

