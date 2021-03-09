from firefly import tess, firefly
from subprocess import run
from pandas import read_csv
import numpy as np
import signal
import ray
ray.init(num_cpus = 2)


class timeout:
    def __init__(self, seconds=1, error_message='Timeout'):
        self.seconds = seconds
        self.error_message = error_message
    def handle_timeout(self, signum, frame):
        raise TimeoutError(self.error_message)
    def __enter__(self):
        signal.signal(signal.SIGALRM, self.handle_timeout)
        signal.alarm(self.seconds)
    def __exit__(self, type, value, traceback):
        signal.alarm(0)
        
@ray.remote
def main(exoplanet):
    # Kill process in 5 days
    # max_runtime = 60 * 60 * 24 * 5
    # with timeout(seconds=max_runtime):
    print(f'Executing firefly on {exoplanet}')
    firefly(
        # Firefly Interface
        exoplanet,
        archive='nasa',
        curve_sample=1,
        clean=True,
        cache=False,
        auto=True,
        # MAST Search
        hlsp=['SPOC'],
        cadence=[120],
        bitmask='default',
        # TransitFit Variables
        walks=25,
        slices=5,
        cutoff=0.25,
        window=8,
        nlive=625,
        fit_ttv=False,
        detrending_list=[['nth order', 2]],
        dynesty_sample='rslice',
        fitting_mode='folded',
        limb_darkening_model='quadratic',
        ld_fit_method='coupled',
        max_batch_parameters=25,
        batch_overlap=2,
        dlogz=0.01,
        maxiter=None,
        maxcall=None,
        dynesty_bounding='multi',
        normalise=True,
        detrend=True,
        detrending_limits=[[-10,10]],
        # Plotting
        #plot=True,
        #marker_color='dimgray',
        #line_color='black',
        #bin_data=True,
        #binned_color='red',
        )
    print(f'Finished {exoplanet}')

# Set how many targets to run in parallel
# Define various lists to pass
# targets, all_targets, ttv_targets = tess(archive='nasa', survey=None)

# all_targets = [
#     'wasp100b', 'wasp126b', 'lhs1815b', 'kepler42c', 'wasp119b',
#     'wasp18b', 'hip65ab', 'l9859b', 'gj1252b', 'wasp62b',
#     'toi157b',
# ]
nasa = read_csv('nasa_tess_viable.csv')['Exoplanet'].values
spear = read_csv('firefly/data/spear.csv')['pl_name'].values
diff = np.setdiff1d(nasa, spear).tolist()
# Redefine list to use here
exoplanets = diff[::]
para = [main.remote(exoplanet) for i, exoplanet in enumerate(exoplanets)]
results = [ray.get(main.remote(exoplanet)) for i, exoplanet in enumerate(exoplanets)]
