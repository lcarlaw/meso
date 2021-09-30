import numpy as np
from scipy.ndimage import gaussian_filter
from configs import SCALAR_PARAMS, VECTOR_PARAMS, SIGMA

def filter(data):
    """
    Perform smoothing and data filtering

    Parameters:
    -----------
    data: list of dictionaries
        A list of dictionaries, with each entry corresponding to a model time step or
        forecast hour, containing the data to be filtered and smoothed

    Returns:
    --------
    data: list of dictionaries
        Filtered and smoothed data. Same form as input.
    """

    for t in range(len(data)):
        vars = list(data[t].keys())

        # Smooth the scalar parameter fields. SIGMA specified in the configs.py file.
        for v in vars:
            if v in SCALAR_PARAMS.keys():
                data[t][v] = gaussian_filter(data[t][v], sigma=SIGMA)

        # Filter the effective bulk wind difference (< 20 kts)
        idx = [i for i, s in enumerate(vars) if 'ebwd' in s]
        if len(idx) > 0:
            tmp = np.sqrt(data[t][vars[idx[0]]]**2 + data[t][vars[idx[0]]]**2)
            data[t][vars[idx[0]]] = np.where(tmp < 20, 0, data[t][vars[idx[0]]])
            data[t][vars[idx[1]]] = np.where(tmp < 20, 0, data[t][vars[idx[1]]])

        # Filter the 0-1 km shear (< 10 kts)
        idx = [i for i, s in enumerate(vars) if 'shr1' in s]
        if len(idx) > 0:
            tmp = np.sqrt(data[t][vars[idx[0]]]**2 + data[t][vars[idx[0]]]**2)
            data[t][vars[idx[0]]] = np.where(tmp < 10, 0, data[t][vars[idx[0]]])
            data[t][vars[idx[1]]] = np.where(tmp < 10, 0, data[t][vars[idx[1]]])

        # Filter the 0-3 km shear (< 20 kts)
        idx = [i for i, s in enumerate(vars) if 'shr3' in s]
        if len(idx) > 0:
            tmp = np.sqrt(data[t][vars[idx[0]]]**2 + data[t][vars[idx[0]]]**2)
            data[t][vars[idx[0]]] = np.where(tmp < 20, 0, data[t][vars[idx[0]]])
            data[t][vars[idx[1]]] = np.where(tmp < 20, 0, data[t][vars[idx[1]]])

    return data
