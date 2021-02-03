from firefly import tess
from subprocess import run

def targets():
    targets, all_targets, ttv_targets = tess(k=25, survey='WASP')
    all_targets = ['lhs1815b', 'wasp119b', 'toi157b', 'kepler42c', 'hip65Ab',
                   'l9859b', 'gj1252b', 'wasp62b', 'hd213885b', 'qatar1b']
    for i, exoplanet in enumerate(all_targets):
        run(['python', 'main.py', exoplanet])

targets()

