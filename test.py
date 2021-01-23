import ase
from ase.build import fcc100
from pdb import set_trace



def _Modify_poscar(poscar, n):
    with open(poscar) as file:
        fp = file.readlines()
    set_trace()
    fp.insert(n, fp[0])
    with open(poscar, 'w') as file:
        for line in fp:
            file.write(line)


if __name__ == '__main__':
    poscar = 'poscar3.vasp'
    n = 5
    _Modify_poscar(poscar, n)
