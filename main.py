from firefly import firefly, tess_viable


def main():
    targets, all_targets = tess_viable(k=5)
    firefly(
        # Firefly Interface
        targets, 
        curve_sample=0.01, 
        email=True,
        to=['transitfit.server@gmail.com'], 
        clean=False,
        cache=True,
        # TransitFit Variables
        cutoff=0.25,
        window=2,
        nlive=500, 
        fit_ttv=False,
        detrending_list=[['nth order', 2]],
        dynesty_sample='rwalk', 
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

main()
