import numpy as np
from ase import Atoms
from ase.io import read, write
from ase.geometry import get_distances

L   = 48.0                       
tol = 2.0                        
rng = np.random.default_rng(42)  

molecules = [('PbI2.pdb', 200), ('MAI.pdb',  200), ('gbl.pdb',  112)]

box = Atoms(cell=[L, L, L], pbc=True)

for filename, n in molecules:
    template = read(filename)
    placed = 0
    while placed < n:
        mol = template.copy()

        # 1. random orientation: spin by a random angle about a random axis
        mol.rotate(rng.uniform(0, 360), rng.normal(size=3), center='COM')

        # 2. random position: move the centre of mass to a random point
        point = rng.uniform(0, L, 3)
        mol.translate(point - mol.get_center_of_mass())

        # 3. overlap check against everything already in the box
        if len(box) > 0:
            _, d = get_distances(mol.positions, box.positions, cell=box.cell, pbc=True)

            if d.min() < tol:
                continue            

        box += mol
        placed += 1

    print(f'{filename}: {placed} molecules placed')

box.wrap()
write('data.lmp', box, format='lammps-data', masses=True, specorder=['C', 'H', 'I', 'N', 'O', 'Pb'])
print(len(box), 'atoms written to data.lmp')