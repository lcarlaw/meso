# meso
This repo will create GR-readable placefiles of various parameters important to severe weather forecasting using data from the High Resolution Rapid Refresh (HRRR) or Rapid Refresh (RAP) models. The goal is to provide a very-near-realtime dataset that mimics the [SPC Mesoanalysis](https://www.spc.noaa.gov/exper/mesoanalysis/new/viewsector.php?sector=20#), but at a higher spatial resolution for use in GRAnalyst. Slight differences may be noted due to SPC's use of an objective analysis pass that incorporates recent surface observations, a step that is not performed with this code at this time.

This repo also has an archived mode, allowing the creation post-event reanalyses for use in local case studies.

## Code Execution Time Improvements Using Numba
The main overhaul was to accelerate the CPU-intensive thermodynamic calculations performed by SHARPpy using the [Python Numba](http://numba.pydata.org/) module. This is a non-trivial task, as several standard Python modules and code are not supported by Numba (use of `**kwargs` in function calls, masked numpy arrays, among other things). In addition, the nature of the "Just-in-time" compilation requires explicit `type` declarations within Python `Classes`. As a result, the original SHARPpy code was parsed out line-by-line to allow it to work with Numba and the `@njit` decorator, and some flexibility has certainly been lost here. In this current iteration of code, it's assumed that the input meteorological arrays are full without any missing/masked data.

The main overhead is due to the nature of "just-in-time" compilation whereby each of the "jitted" Python functions are translated and optimized to machine code. This is well worth it in this case due to the computationally-expensive lifting routines in SHARPpy, resulting in execution improvements which are orders of magnitude better than pure, interpreted, Python.

Here is how a few benchmarks compare run on a 2019 Macbook Pro with a 2.3 GHz Intel Core i9-9880H CPU (8 cores, 16 threads). The domain is 240 x 210 with a 13 km grid-spacing, discounting initial JIT times:

| Test Description      | Jitted? | SHARPpy Execution Time | Percent Change      |
| --------------------- | ------- | ---------------------- | ------------------- |
| Parallel (16 threads) | Yes     | 24.13s                 | N/A                 |
| Serial (1 thread)     | Yes     | 124.35s                | 415.34%             |
| Serial (1 thread)     | No      | 1959.84s               | 8022.01%            |

### To do:
- [ ] Build simple Cressman or barnes surface analysis?
- [ ] Create gridded "on-the-fly" storm motion to refine ESRH, deviant tornado, etc. calculations.
- [ ] Figure out GR's polygon fill rules: stripes on contour-filled plots?
- [X] Centralized hosting of `windicons.png`. Changes to `barbconfigs` in config file
- [ ] Add capability to output images (.tif, high-res .png) of mesoanalysis parameters & upper-air variables
- [X] Allow bundling of several placefiles together (i.e. MUCAPE and EBWD)
- [X] Improve download execution for archived THREDDS requests
- [X] Build in automated checks for hung processes in the `run.py` driver
- [X] Add better error logging to the download step

## Basic Setup Notes
The setup proceeds using Anaconda, as well as assuming a completely vanilla Python3 install.

### Creating the base environment
Run the following to create an Anaconda environment called `meso` with the required libraries:

```
conda env create -f environment.yml
```

### Dependencies and config file alterations
You will need working `wget` and `wgrib2` binaries on your filesystem. Add these to the `WGRIB2` and `WGET` variables in the `config.py` file. Specify a location for logfile output via the `LOG_DIR` variable, the default directory to download model data to in `MODEL_DIR`, and where to send the final placefiles in `OUTPUT_DIR`.

If you want to run this in realtime mode (see section below), specify the `PYTHON` variable to point to the particular anaconda `meso` environment (or whatever you named it) on your filesystem. If desired, locations to `WIND_ICONS` and `SHEAR1_ICONS` can also be specified in the config file.

#### Installing the latest WGRIB2 binary
Wgrib2 version 3.0.2 or higher is necessary for time interpolations and for decoding older versions of the RAP. The latest wgrib2 binary has an added flag called `new_grid_order` which is needed if you want to use this repository to read older RUC data stored on the NCEI THREDDS servers. Some of the older RAP/RUC grib files store the UGRD and VGRD entries in separate "blocks", and wgrib2 needs these to be paired together, one VGRD after a UGRD entry. The steps to install (at least on my 2019 Macbook Pro running 10.15.3 Catalina) were straightforward, although I needed a separate `gcc` install than the pre-packaged XCode version on my machine which was installed via [`homebrew`](https://brew.sh/). This may be different on your machine. If not using `brew`, the usual cautions of installing external libraries on your local machine apply.

```
brew install gcc@9
cd etc
tar -xvzf wgrib2.tgz.v3.0.2
cd grib2
export CC=/usr/local/bin/gcc-9
export CXX=c++
export FC=gfortran-9
make -j4
```

If `make` is successful, you should have an executable `wgrib2` file in the `etc/grib2/wgrib2/` directory.

If desired, you can then softlink this into the standard location on most file systems with a `sudo ln -s etc/grib2/wgrib2/wgrib2 /usr/local/bin`. Either way, update the `WGRIB2` variable in the `configs.py` file to point to the `WGRIB2` executable's location.

#### Set the perl scripts to executable
Before being able to cron the driver script, you may also have to type: `chmod +x IO/*.pl` within the parent meso directory to make the various perl scripts executable.

## Running in real time
We can avoid the use of the system's cron-scheduler and issues with system `$PATH` variables by using the `schedule` module in Python.

The `run.py` script is setup to activate the `get_data.py` module at :54 after each hour. In realtime mode, this downloads real time 13-km RAP data for all runs except at 00 and 12z when HRRR data are used (the 00 and 12z RAP cycles are delayed by ~30 minutes to assimilate some RAOBs). Both are on the native hybrid-sigma coordinate system to maximize the number of model levels near the surface to improve SHARPpy's lifting calculations. The HRRR runs are upscaled to the 13-km RAP horizontal grid-spacing. 1 and 2-hour forecasts are linearly interpolated in time to produce a 1.5 hour forecast to update meteorological parameters three times per hour.

Assuming all variable `PATHs` have been set correctly in the `config.py` file, to automate, type:

```
conda activate meso
python run.py
```

Important log files will be located in the `LOG_DIR` directory. These can all be monitored in-line with `tail -f *.log`. `process.py` will run after model data has been successfully downloaded, and will output placefiles in the `output` directory. These files will automatically time match in GR, with an update occurring at :15 and :45.

##  Adding parameters
Parameters to be output as placefiles are defined in the config file in the `SCALAR_PARAMS` and `VECTOR_PARAMS` dictionaries. The key-value pairs will be used for the placefile name and verbose info string, respectively. Base plot style specifications in the `contourconfigs` and `barbconfigs` dictionaries are overridden by individual entries in the `PLOTCONFIGS` dictionary with each key needing to match a key in either `SCALAR_PARAMS` or `VECTOR_PARAMS`.

Parameter calculations are performed in `calc.compute`. Numba does not allow easy passing of a dictionary registry containing pointers to parameter calculation functions, so each parameter must be specified as an if-block within function `worker`. Specific calculations must be added as a function in `sharptab.derived` and called in the associated if-block within the `worker` function. For information on how to add SHARPpy calculations, see the [SHARPpy Scripting documentation](https://sharppy.github.io/SHARPpy/scripting.html).

### Parameter example
Suppose we want to add 0-500 m Storm-Relative Helicity as a placefile. To do this:
1. Add an entry to the `SCALAR_PARAMS` dictionary in the `configs.py` file:
```
'srh500': '0-500 m Storm-Relative Helicity'
```
2. Add a function to the `calc.derived.py` file. For example:
```
@njit
def srh500(prof):
    RM5 = rm5(prof)
    srh = winds.helicity(prof, 0, 500, stu=RM5[0], stv=RM5[1])[0]
    return srh
```
3. Add an if-block to the `calc.compute.py` file:
```
if 'srh500' in SCALARS: d['srh500'][j,i] = derived.srh500(prof, eff_inflow)
```
4. Specify contouring configurations in the `PLOTCONFIGS` dictionary within the `configs.py` file. The dictionary key must match the key we used in step 1. If no entry is provided, plotting parameters will be specified via the `barbconfigs` or `contourconfigs` dictionaries.

## Specifying multi-parameter placefile bundles  
It is possible to load multiple placefiles with one entry by concatenating placefiles together. To turn this feature on, add entries to the `BUNDLES` dictionary in the `configs.py` file. The individual parameter names must match dictionary keys in the `SCALAR_PARAMS` or `VECTOR_PARAMS` dictionaries.

## Creating an archived case
### Download the model data
The `get_data.py` script can download archived 1-hour forecasts either from the NCEI THREDDS or Google Cloud servers. `get_data.py` will accept a single time or a time range.

For this example, we'll download HRRR data during the August 10th, 2020 Midwest Derecho:

```
python get_data.py -s 2020-08-10/17 -e 2020-08-10/23 -m HRRR
```

Archived native hybrid-sigma coordinate HRRR data will be downloaded into the `MODEL_DIR` directory specified in the config file (or to a specific directory with the -p flag) and upscaled to 13 km grid spacing (same as the RAP). The original data files (>500 MB) will be deleted from the file system.

The HRRR archive on the Google Cloud appears to go back to the 18z run on 2014-07-30, while RAP data is available back to the 00z run on 2021-02-22. The RAP/RUC archive on the THREDDS server goes back further, but you may notice more errors when downloading due to data response latencies during the web retrieval steps, as well as missing model runs.

#### Special Notes
If you request a considerable amount of data (numerous runs, HRRR data, etc.), you may need to increase the `TIMEOUT` variable in the `configs.py` file.

### Creating placefiles
```
python process.py -s 2020-08-10/17 -e 2020-08-10/23 -meso
```

You can view logs with `tail -f ./logs/*.log`. This will take a few minutes. When the scripts finish, text placefiles should be available in the `output` directory. These will be named with a trailing `YYYYmmddHH-YYYYmmddHH` corresponding to the valid times of the data within the placefiles. Data will automatically time-match in GR to the closest hour.

### Hodographs
![](https://github.com/lcarlaw/meso/blob/master/hodograph_example.png)

Usage to create hodographs is as follows:

```
python process.py [ -s START_TIME ] [ -e END_TIME ] [ -t TIME ] [ -hodo LAT_LON ]
                  [ -m STORM_MOTION ] [ -sw SFC_WIND ] [ -sr ]
```

* `START_TIME`: Initial model cycle time (this will be a 1-hr forecast) in the form `YYYY-mm-dd/HH`.
* `END_TIME`: Last model cycle time (this will be a 1-hr forecast) in the form `YYYY-mm-dd/HH`. Required if `START_TIME` has been declared.
* `TIME`: For a single model cycle time in the form `YYYY-mm-dd/HH`. Don't use with `-s` and `-e` specified.
* `LAT_LON`: Latitude/longitude pair for hodograph creation. Form is `LAT/LON`. Currently only accepts one point at a time.
* `STORM_MOTION`: The storm motion vector. It can take one of two forms. The first is either `BRM` for the Bunkers right-mover vector or `BLM` for the Bunkers left-mover vector. The second form is `DDD/SS`, where `DDD` is the direction the storm is coming from in degrees, and `SS` is the storm speed in knots. An example might be 240/35 (from the WSW at 35 kts). If the argument is not specified, the default is to use the Bunkers right-mover vector.
* `SFC_WIND`: The surface wind vector. Its form is the same as the `DDD/SS` form of the storm motion vector. If none is specified, the lowest-model level wind will be used.
* `-sr`: Storm-relative flag. If set, hodographs are altered to plot in a storm-relative sense, similar to [Cameron Nixon's work here](https://cameronnixonphotography.wordpress.com/research/the-storm-relative-hodograph/).

This plotting is done using modified code from Tim Supinie's [vad-plotter repo](https://github.com/tsupinie/vad-plotter).
