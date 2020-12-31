from firefly import firefly, tess_viable

def main():
    targets = tess_viable(k=10)
    firefly(
        # Firefly Interface
        targets, 
        archive='nasa', 
        curve_sample=0.01, 
        email=True,
        to=['transitfit.server@gmail.com'], 
        clean=False,
        # TransitFit Variables
        cutoff=0.25,
        window=2.5,
        nlive=300, 
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

main()
