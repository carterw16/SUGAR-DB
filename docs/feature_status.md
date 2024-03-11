### Development (Work-In-Progress)
   * IBDG model and IBDG Control (Q-V, P-V, etc.)
   * Converting necessary files to Cython (ongoing)
    
### Not Yet Implemented
  * Variable Limiting
  * Optimization
  * Models:
    * PV Generator
    * Phase Shifting Transformers
    * Induction Generator
    * Storage
    * Houses
    * Water Heater
    * Regulators:
        * Open Delta
        * Closed Delta
    * Loads:
        * Exponential: voltage dependency of P and Q is defined by exponential
          parameters.
        * Constant P and Quadratic Q: the reactive power, Q, varies quadratically with
          the voltage (as a constant reactance) while the active power, P, is independent from the
          voltage, somewhat like a motor.
        * Constant P and Fixed Q: Q is a fixed value independent of time and voltage.
        * Constant P and Quadratic Q: Q varies with square of the voltage.
  * Parsing:
    * GridLab-D houses, inverters, diesel generators, and water heaters.
    * OpenDSS photovoltaics and storage.
  * Stamping only phases that are present instead of all phases
 
#### Needs Validating
 * Gmin Stepping

#### To Fix
   * Min and max bus voltage reports in the output summary file and during homotopy iterations
   * Current Measurements:
     * Loads
     * Delta-Connected components
