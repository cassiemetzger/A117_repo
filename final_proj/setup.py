import os 
print("Installing dependencies...")
os.system("pip install -r requirements.txt")
print("Done!")

import numpy as np
import matplotlib.pyplot as plt
import astropy 
from astropy.io import fits, ascii
from astropy.table import Table
import astroquery 
import glob
from astroquery.mast import Observations, Mast, MastMissions, MastMissionsClass, Catalogs
import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.xmatch import XMatch
from astropy.timeseries import BoxLeastSquares
from astroquery.sdss import SDSS
import fitsio
import pandas as pd

__all__ = [
    "np", 
    "plt", 
    "fits", 
    "ascii", 
    "Table"
    "Observations"
    "u", 
    "SkyCoord", 
    "XMatch", 
    "BoxLeastSquares", 
    "SDSS", 
    "pd"
]
# configure matplotlib settings
params = {'axes.labelsize': 12, 'xtick.labelsize': 12, 'ytick.labelsize': 12, 
          'text.usetex': False, 'lines.linewidth': 1,
          'axes.titlesize': 18, 'font.family': 'serif', 'font.size': 12}
plt.rcParams.update(params)

current_files = glob.glob("*")
if("allStar-dr17-synspec_rev1.fits" not in current_files): 
    print("Downloading APOGEE allStar...")
    Observations.download_file('mast:SDSS/apogee/allStar-dr17-synspec_rev1.fits')
else: 
    print("APOGEE allStar found!")
# function for downloading tess data 
def tess_data_download(num_obj = None, order_random = False, search_radius = 0.1): 
    transiting_exo = ascii.read('transiting_exoplanets.csv')
    os.system("mkdir tess_data")
    if(num_obj is not None): 
        obj_limit = num_obj
    else: 
        obj_limit = len(transiting_exo)
    if(order_random): 
        for i in range(obj_limit): 
            idx = np.random.randint(0, len(transiting_exo))
            row = transiting_exo[idx]
            ra = row['ra']
            dec = row['dec']
            observations = Observations.query_region(f"{ra} {dec}", radius=f"{search_radius} deg")
            obs_wanted = ((observations['dataproduct_type'] == 'timeseries') &
            (observations['obs_collection'] == 'TESS'))
            data_products = Observations.get_product_list(observations[obs_wanted])
            products_wanted = Observations.filter_products(data_products, 
                                                productSubGroupDescription=["LC"])
            obj_name = row['pl_name']
            obj_name=obj_name.replace(" ", "_")
            os.system(f"mkdir tess_data/{obj_name}")
            download = Observations.download_products(products_wanted)
            for file in download['Local Path']: 
                os.system(f"cp {file} ./tess_data/{obj_name}")
        os.system(f"rm -r mastDownload")
    else: 
        for i, row in enumerate(transiting_exo): 
            if(i < obj_limit): 
                ra = row['ra']
                dec = row['dec']
                observations = Observations.query_region(f"{ra} {dec}", radius=f"{search_radius} deg")
                obs_wanted = ((observations['dataproduct_type'] == 'timeseries') &
                (observations['obs_collection'] == 'TESS'))
                data_products = Observations.get_product_list(observations[obs_wanted])
                products_wanted = Observations.filter_products(data_products, 
                                                    productSubGroupDescription=["LC"])
                obj_name = row['pl_name']
                obj_name=obj_name.replace(" ", "_")
                os.system(f"mkdir tess_data/{obj_name}")
                download = Observations.download_products(products_wanted)
                for file in download['Local Path']: 
                    os.system(f"cp {file} ./tess_data/{obj_name}")
            else: 
                break
            os.system(f"rm -r mastDownload")
