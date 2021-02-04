from firefly import tess
from multiprocessing import Pool
from subprocess import run

def main(exoplanet):
    run(['python', 'main.py', f'{exoplanet}']
    )

# Set how many targets to run in parallel
processes = 2
# Define various lists to pass
targets, all_targets, ttv_targets = tess(archive='eu', survey='WASP')
all_targets = [
    'wasp100b', 'wasp126b', 'lhs1815b', 'kepler42c', 'wasp119b',
    'wasp18b', 'hip65ab', 'l9859b', 'gj1252b', 'wasp62b',
    'toi157b',
]
# Redefine list to use here
exoplanet = all_targets
if __name__ == '__main__':
    with Pool(processes=processes) as pool:
        pool.map(main, exoplanet)
