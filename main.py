from firefly import firefly
from sys import argv

def main(exoplanet):
    firefly(
        # Firefly Interface
        exoplanet,
        sigma=3,
        curve_sample=1, 
        email=False,
        to=['transitfit.server@gmail.com'], 
        clean=False,
        cache=False,
        auto=True,
        # TransitFit Variables
        cutoff=0.25,
        window=2,
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

main(argv[1:])
