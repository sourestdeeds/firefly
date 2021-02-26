from firefly import firefly
from sys import argv
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
    with suppress_print():
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
        detrending_limits=[[-10,10]]
        )


main(argv[1:])
