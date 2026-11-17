from os import listdir
from pathlib import Path
from pyshortcuts.utils import read_textfile
from feffdb import FeffDatabase
import csv
from feffdb.utils import parse_feffdat

FIRST_SHELL_TOL = 0.05

ATOMIC_NUMBER = {
    "C": 6, "N": 7, "O": 8, "F": 9, "S": 16, "Cl": 17,
    "Co": 27, "Ni": 28, "Cu": 29, "Se": 34, "Pt": 78,
}

ANION_FLOOR = {"O": 1.5, "S": 1.8, "N": 1.5, "Se": 2.0, "Cl": 2.0,
               "C": 1.5, "F": 1.5}


def _first_shell_cut(row):
    """Longest crystallographic first-shell contact for this structure, + tol."""
    raw = (row.get("shell") or "").strip()
    vals = []
    for part in raw.split(";"):
        part = part.strip()
        if part:
            try:
                vals.append(float(part))
            except ValueError:
                pass
    if not vals:
        try:
            vals = [float(row["mean_R"])]
        except (KeyError, TypeError, ValueError):
            return None
    return max(vals) + FIRST_SHELL_TOL


def firstshell_paths(run, scatterer_z=None, cut=None, floor=None):
    """
    find first shell paths for a Feff run. This may find multiple paths:
      -  single-scattering only
      -  in the first 5 paths listed
      -  up to R = 1.25*shortest path length
      -  total coordination number up to 16
      -  no paths close than 0.03 Ang.
    """
    filesdat = Path(run, 'files.dat')
    if not filesdat.exists():
        return [Path(run, 'feff0001.dat')]

    paths = []
    n_total = 0
    reff0, reff_last = 0., 0.
    if filesdat.exists():
        fdlines = read_textfile(filesdat).split('\n')
        header = True
        for i, line in enumerate(fdlines):
            if line.strip().startswith('----'):
                break
        ihead = i + 2

        for line in fdlines[ihead:]:
            try:
                fname, sig2, amp, deg, nlegs, reff = line.split()
            except:
                break
            deg = float(deg)
            nlegs = int(nlegs)
            reff = float(reff)
            if reff0 < 0.1:  # first path
                reff0, reff_last = reff, reff
                n_total = deg
                paths = [Path(run, fname)]
            else:
                if (nlegs > 2 or n_total > 16 or reff > 1.25*reff0
                    or len(paths) > 4):
                    continue
                if reff > reff_last + 0.03:
                    paths.append(Path(run, fname))
                    reff_last = reff
                    n_total + n_total + deg
    if scatterer_z is None:
        return paths
    kept = []
    for p in paths:
        try:
            dat = parse_feffdat(p.as_posix())
        except Exception:
            continue
        if int(dat['nleg']) != 2 or int(dat['scatterer']) != scatterer_z:
            continue
        if floor is not None and float(dat['reff']) < floor:
            continue
        if cut is not None and float(dat['reff']) > cut:
            continue
        kept.append(p)
    return kept or paths[:1]

def Add_Liao_Calcs():
    top_folder = 'Liao_Calcs'
    cif_folder = 'cif_files'
    feffdat_stub = 'FEFF_paths(output)_'

    for pair in listdir(top_folder):
        fpath = Path(top_folder, pair, f'{feffdat_stub}{pair}')
        cpath = Path(top_folder, pair, cif_folder)
        if not fpath.exists():
            continue
        scat = pair.split('-')[1]
        scatterer_z = ATOMIC_NUMBER.get(scat)
        floor = ANION_FLOOR.get(scat)
        cuts = {}
        csv_path = Path(top_folder, pair, 'feff_bins_summary_AMCSD.csv')
        if csv_path.exists():
            with open(csv_path, newline='') as fh:
                for row in csv.DictReader(fh):
                    cuts[row['id']] = _first_shell_cut(row)
        for drange in listdir(fpath):
            fpathd = Path(fpath, drange)
            if fpathd.is_dir():
                for run in listdir(fpathd):
                    cif_file = Path(cpath, drange, f'{run}.cif')
                    if cif_file.exists():
                        frun = Path(fpathd, run)
                        feff_inp = frun / 'feff.inp'
                        if feff_inp.exists():
                            print(feff_inp)
                            for feff_dat in firstshell_paths(frun, scatterer_z, cuts.get(run), floor):
                                feffdb.add_feffdat(feff_dat, feff_inp, cif_file=cif_file)
                                print(f'added {feff_dat}, {pair=} {drange=} {run=}')



def Add_LarixCalcs():
    top_folder = 'Larix'
    for calc_dir in listdir(top_folder):
        cpath = Path(top_folder, calc_dir)
        cif_file = None
        for fname in listdir(cpath):
            if fname.endswith('.cif'):
                cif_file = cpath / fname
        if cif_file is None:
            continue
        feff_inp = cpath / 'feff.inp'
        for feff_dat in firstshell_paths(cpath):
            feffdb.add_feffdat(feff_dat, feff_inp, cif_file=cif_file)
            print(f'added {feff_dat}, {cif_file} ')
##
if __name__ == '__main__':
    feffdb = FeffDatabase('./feffdat.db', create_ok=True)
    Add_LarixCalcs()
    Add_Liao_Calcs()
    feffdb.close()
