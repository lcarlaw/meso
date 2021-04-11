# meso
This repo will create GR-readable placefiles of various parameters important to severe weather forecasting using data from the High Resolution Rapid Refresh (HRRR) or Rapid Refresh (RAP) models. The goal is to provide a very-near-realtime dataset that mimics the [SPC Mesoanalysis](https://www.spc.noaa.gov/exper/mesoanalysis/new/viewsector.php?sector=20#), but at a higher spatial resolution for use in GRAnalyst, as well as the ability to easily produce post-event reanalyses for use in local case studies.

## Code Execution Time Improvements Using Numba
The main overhaul here was to massively accelerate the CPU-intensive thermodynamic calculations performed by SHARPpy (mainly from parcel lifting calculations) using the [Python Numba](http://numba.pydata.org/) module. This is a non-trivial task, as several standard Python modules and code are not supported by Numba. In addition, the nature of the "Just-in-time" compilation requires explicit `type` declarations within Python `Classes`. As a result, the original SHARPpy code had to be parsed out line-by-line to allow it to work with Numba and the `@njit` decorator, and some flexibility has certainly been lost here. The biggest issues were the lack of `**kwarg` support and numpy masked arrays. In this current iteration of code, it's assumed that the input meteorological arrays are full without any missing/masked data.

The main overhead--roughly 40-45 seconds for each run--is due to the nature of "just-in-time" compilation whereby each of the "jitted" Python functions are translated and optimized to machine code. This is well worth it, however, for the computation time improvements which are orders of magnitude better than pure, interpreted, Python.

Here is how a few benchmarks compare. The domain is 240 x 210 with a 13 km grid-spacing:

| Test Description      | Jitted? | Execution Time | Percent Improvement |
| --------------------- | ------- | -------------- | ------------------- |
| Parallel (16 threads) | Yes     | 91.21s         | N/A                 |
| Serial (1 thread)     | Yes     | 124.35s        | 36.34%              |
| Serial (1 thread)     | No      | 1959.84s       | 2048.71%            |

### To do:
- Build in automated checks for hung processes in the `run.py` driver
- Investigate Python Ray and cluster-computing
- Global Interpreter Lock (GIL) bottlenecks in numba?
- Support for CUDA/GPU-based calculations for further speed increases?  

## Basic Setup Notes
The setup here proceeds using Anaconda, as well as assuming a completely vanilla Python3 install. I've also edited my `~/.condarc` file to add conda-forge to the default channels.

### Exporting the base environment

```
conda env export --no-builds > environment.yml
```

### Creating the base environment
Run the following to create an Anaconda environment called `meso` with the required libraries:

```
conda env create -f environment.yml
```

### Dependencies and config files
You will need working `wget` and `wgrib2` binaries on your filesystem. Add these to the `WGRIB2` and `WGET` variables in the `config.py` file.

#### Installing the latest WGRIB2 binary
The latest version of wgrib2 has an added flag called `new_grid_order` which is necessary if you want to use this repository to read older RUC data stored on the NCEI THREDDS servers. The basic information here is that some of the older RAP/RUC grib files store the UGRD and VGRD entries separately, and wgrib2 needs these to be paired together, one VGRD after a UGRD entry. The steps to install (at least on my 2019 Macbook Pro running 10.15.3 Catalina) were straightforward, although I needed a separate `gcc` install than the pre-packaged XCode version on my machine which was installed via [`homebrew`](https://brew.sh/). This may be different on your machine.

```
brew install gcc@9
cd libs
tar -xvzf wgrib2.tgz.v3.0.2
cd grib2
export CC=/usr/local/bin/gcc-9
export CXX=/usr/local/bin/gcc++-9
export FC=gfortran
make
```

If `make` is successful, you should have an executable `wgrib2` binary in the `libs/grib2/wgrib2/` directory.

If desired, you can then softlink this into the standard location on most file systems with a `sudo ln -s libs/grib2/wgrib2/wgrib2 /usr/local/bin`. Either way, update the `WGRIB2` variable in the `configs.py` file to point to the `WGRIB2` binary location.

#### Set the perl scripts to executable
Before being able to cron the driver script, you may also have to type: `chmod +x IO/*.pl` within the parent meso directory.

## Running in real time
We can avoid the use of the system's cron-scheduler and issues with system `$PATH` variables by using the `schedule` module in Python.

The `run.py` script is setup to activate the `get_data.py` module at :56 after each hour to download real time HRRR data on the native hybrid-sigma coordinate system and then upscale to the 13-km RAP horizontal grid-spacing. The native coordinate system affords more accurate parcel calculations in SHARPpy over the smaller isobaric files. The NOMADS server is queried first, followed by the backup FTPPRD servers, and finally the GOOGLE cloud, which uploads HRRR data in near-real time. Once data is available on the local filesystem, `process.py` is called to produce GR-readable placefiles. To run, simply type:

```
conda activate meso
python run.py
```

The log file can be found in `logs/downloads.log`. You can leave this process running indefinitely.

## Creating an archived case
To-do.
