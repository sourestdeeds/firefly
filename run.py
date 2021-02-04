from firefly import tess
from multiprocessing import Pool
from subprocess import run
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
    run(['python', 'main.py', f'{exoplanet}']
    )

# Set how many targets to run in parallel
processes = 1
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
        if processes==1:
            pool.map(main, exoplanet)
        else:
            with suppress_print():
                pool.map(main, exoplanet)
