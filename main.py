from firefly import tess, firefly
from multiprocessing import Pool
import sys
import os


class suppress_print():
    def __enter__(self):
        self.original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self.original_stdout
        

def main(exoplanet):
    firefly(
    # Firefly Interface
    exoplanet,
    archive='eu',
    curve_sample=1,
    email=False,
    to=['transitfit.server@gmail.com'],
    clean=True,
    cache=False,
    auto=True,
    # TransitFit Variables
    cutoff=0.25,
    window=2.5,
    nlive=1000,
    fit_ttv=False,
    detrending_list=[['nth order', 2]],
    dynesty_sample='rslice',
    fitting_mode='folded',
    limb_darkening_model='quadratic',
    ld_fit_method='coupled',
    max_batch_parameters=25,
    batch_overlap=2,
    dlogz=None,
    maxiter=None,
    maxcall=None,
    dynesty_bounding='multi',
    normalise=True,
    detrend=True
    )


processes = 1
targets, all_targets, ttv_targets = tess(archive='eu', survey='WASP')
all_targets = [
    'wasp100b', 'wasp126b', 'lhs1815b', 'kepler42c', 'wasp119b',
    'wasp18b', 'hip65ab', 'l9859b', 'gj1252b', 'wasp62b',
    'toi157b',
]

if __name__ == '__main__':
    with Pool(processes=processes) as pool:
        if processes==1:
            pool.map(main, all_targets)
        else:
            with suppress_print():
                pool.map(main, all_targets)
