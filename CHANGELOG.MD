# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2.0.3] - 2021-06-14
### Added
- Incorporated faster revised homotopy ([Homotopy.py](classes/Homotopy.py)).
- Modified [setup.py](setup.py) to be more efficient and not consistently update.

## [2.0.2] - 2020-10-18
### Added
- Removed print statements and replaced them with logging. All logs will go to a "log/" directory.
- Cythonized fuse, switch, transformer, and line stamps ([FuseStamp.pyx](classes/FuseStamp.pyx), [SwitchesStamp.pyx
](classes/SwitchesStamp.pyx), [TransformerStamp.pyx](classes/TransformerStamp.pyx), 
[Linestamp.pyx](classes/LinesStamp.pyx))
- J1 EPRI OpenDSS test case

### Changed
- Fixed line, transformer, switch, and fuse current measurements. Previous current measurements impacted the final
 voltages and were implemented incorrectly.
See [CurrentMeasurements.py](classes/CurrentMeasurements.py), [Transformers.py](classes/Transformers.py), 
[Switches.py](classes/Switches.py), [Fuses.py](classes/Fuses.py) for more.
- Modified [Regulators.py](classes/Regulators.py) slightly in order to ensure only fixed taps are enabled for now
. Also, modified default regulator type depending on if the test case is Gridlab-D or OpenDSS.
- Cleaned up outputs. Fuse, Switch, and Transformer currents are now saved as outputs 
[save_output.py](lib/save_output.py).
- Parser modified to include new components (IBDGs) and case objects are now saved to a dictionary. 
See [parser.py](lib/parser.py).
- New object properties added to [retrieve_attributes.py](lib/retrieve_attributes.py).
- Removed unnecessary code and some global variables from [SUGAR3.py](SUGAR3.py).
- G and B homotopy values are now input parameters to provide to a settings.json or settings dictionary. See
 [runSUGAR3.py](runSUGAR3.py) for an example.
- DiTTo Gridlab-D reader modified to parse slack or swing node voltages.

## [2.0.1] - 2020-07-27
### Added
- Neutral phase stamp for reactors.
- Exponential Load model to [Loads.py](classes/Loads.py)
- Voltage Setpoints for OpenDSS test circuits.
- K1 EPRI OpenDSS test circuit.

### Changed
- If test case is OpenDSS, the voltage setpoint is now taken into consideration. Changes made to 
[Nodes.py](classes/Nodes.py) and [Generators.py](classes/Generators.py) illustrate this change.
- Changed OpenDSS unit test to ckt5.
- Updated [test case status](testcases/testcasestatus.xlsx).

## [2.0.0] - 2020-07-23
Officially on version 2.0.0 of the SUGAR three-phase solver. 2.0.0 chosen for clarity and lack of semantic
 versioning in the past.
### Added
- New parser created with use of DiTTo library. See [parser.py](lib/parser.py).
- Reactors are now an element available for use. See [reactors.py](classes/Reactors.py)
- Unit testing and the ability to validate results by comparing against ground truths. See [validate](validate).
- [Makefile](Makefile) to standardize common tasks.
- [Style Guide](docs/styleguide.md) to maintain order.
- [Contributing Guide](docs/contributing.md) to help lost souls figure out where to start.
- [Feature Status](docs/feature_status.md) which list things yet to be implemented.
- [File Structure](docs/filestructure.md) which highlights how the branch is organized.
- Auto formatting libraries to keep code consistent.
- [OpenDSS test cases](testcases/opendss).
- Test case status tracker. See [testcasestatus](testcases/testcasestatus.xlsx).

### Changed
- Capacitors, CurrentMeasurements, Fuses, Generators, Lines, Loads, Nodes, Reactors, Regulators, Switches
, Transformer, and TriplexLoads.
- [README](README.MD) updated.
- How outputs are saved and what information is collected. See [save_output.py](lib/save_output.py).

### Removed
- GridParser is dead. Long live GridParser! :)
- Miscellaneous test cases that held no value.
- Island checking. Will need to be fixed in a future update.
- Transformer Configuration and Line Configurations removed because they are no longer necessary.
