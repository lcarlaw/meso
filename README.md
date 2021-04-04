# meso
This repo will create GR-readable placefiles of various parameters important to severe weather forecasting using data from the High Resolution Rapid Refresh (HRRR) model. The goal is to provide a very-near-realtime dataset that mimics the [SPC Mesoanalysis](https://www.spc.noaa.gov/exper/mesoanalysis/new/viewsector.php?sector=20#), but at a higher spatial resolution for use in GRAnalyst. 

## Code Execution Time Improvements Using Numba
The main overhaul here was to massively accelerate the CPU-intensive thermodynamic calculations performed by SHARPpy (mainly from parcel lifting calculations) using the [Python Numba](http://numba.pydata.org/) module. This is a non-trivial task, as several standard Python modules and code are not supported by Numba. In addition, the nature of the "Just-in-time" compilation requires explicit `type` declarations within Python `Classes`. As a result, the original SHARPpy code had to be parsed out line-by-line to allow it to work with Numba and the njit decorator, and some flexibility has certainly been lost here. The biggest issues were the lack of `**kwarg` support and numpy masked arrays. In this current iteration of code, it's assumed that the input meteorological arrays are full without any missing/masked data.

The main overhead--roughly 30 seconds for each run--is due to the nature of "just-in-time" compilation whereby each of the "jitted" Python functions are translated and optimized to machine code. This is well worth it, however, for the compuation time improvements which are orders of magnitude better than pure python. 

### To do:
- Investigate Python Ray and cluster-computing
- Global Interpreter Lock (GIL) bottlenecks in numba?  

## Basic Setup Notes
The setup here proceeds using Anaconda, as well as assuming a completely vanilla Python3 install. I've also edited my `~/.condarc` file to add conda-forge to the default channels. 

### Creating the base environment
Run the following to create an Anaconda environment with the required libraries:

```
conda env create -f environment.yml
pip install geojsoncontour
```

### Dependencies and config files
You will need working `wget` and `wgrib2` binaries on your filesystem. Add these to the `WGRIB2` and `WGET` variables in the `config.py` file.

Before being able to cron the driver script, you may also have to type: `chmod +x IO/*.pl` within the parent meso directory. 
