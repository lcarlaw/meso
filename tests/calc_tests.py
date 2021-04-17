import sys
from pathlib import Path
file = Path(__file__).resolve()
parent, root = file.parent, file.parents[1]
sys.path.append(str(root))
import IO.read as read
import numpy as np
import pandas as pd
import sharppy
import sharppy.sharptab.profile as profile
import sharppy.sharptab.interp as interp
import sharppy.sharptab.winds as winds
import sharppy.sharptab.utils as utils
import sharppy.sharptab.params as params
import sharppy.sharptab.thermo as thermo

import sharptab.interp as interp

pressure = read.read_data('./rap.t14z.awp130pgrbf01.grib2.reduced')
native = read.read_data('./rap.t14z.awp130bgrbf01.grib2.reduced')

point = [[-97.97, 27.5]]
idx = interp.nearest_idx(point, pressure['lons'], pressure['lats'])[0]


j, i = idx
n_levs = np.shape(pressure['pres'][:,j,i])[0]

profile_pressure = profile.create_profile(pres=pd.unique(pressure['pres'][:,j,i]),
                                         hght=pd.unique(pressure['hght'][:,j,i]),
                                         tmpc=pd.unique(pressure['tmpc'][:,j,i]),
                                         dwpc=pd.unique(pressure['dwpc'][:,j,i]),
                                         wdir=pd.unique(pressure['wdir'][:,j,i]),
                                         wspd=pd.unique(pressure['wspd'][:,j,i]))
mupcl = params.parcelx(profile_pressure, flag=3)
print(mupcl.bplus)


profile_native = profile.create_profile(pres=native['pres'][:,j,i], hght=native['hght'][:,j,i],
                                      tmpc=native['tmpc'][:,j,i], dwpc=native['dwpc'][:,j,i],
                                      wdir=native['wdir'][:,j,i], wspd=native['wspd'][:,j,i])
mupcl = params.parcelx(profile_native, flag=3)
print(mupcl.bplus)
