from pdb import set_trace
import os
from ase.build import cut
from ase.io import read, write
from ase.visualize import view
import numpy as np




def write_poscar(Atoms, path):
    with open(path, 'w') as file:
        n1 = np.size(np.where(Atoms.numbers == 79))
        n2 = np.size(np.where(Atoms.numbers == 29))
        file.write('Au%sCu%s\n' % (n1, n2))
        file.write('1.00000\n')
        cell = Atoms.cell
        if len(cell) == 1:
            cell = np.array([[cell[0], 0, 0], [0, cell[0], 0], [0, 0, cell[0]]])
        for line in cell:
            str_array=" ".join(map(str, line))
            file.write(str_array + '\n')
        file.write('Au Cu\n')
        file.write('%s %s\n' % (n1, n2))
        file.write('Direct\n')
        position = Atoms.positions @ np.linalg.inv(cell)
        pos_num = Atoms.numbers
        set_trace()
        pos_sym = []
        for i in pos_num:
            if i == 79:
                pos_sym.append('Au')
            else:
                pos_sym.append('Cu')
        for i, line in enumerate(position):
            str_array = " ".join(map(str, line))
            file.write(str_array + ' ' + pos_sym[i] + '\n')

        # file.write(cell)


def _Modify_poscar(poscar, n):
    with open(poscar) as file:
        fp = file.readlines()
    set_trace()
    fp.insert(n, fp[0])
    with open(poscar, 'w') as file:
        for line in fp:
            file.write(line)


""" 读取结构扩胞 """
pri_struct = read('POSCAR1.vasp')   # 读取原胞POSCAR结构
write(filename='poscar3.vasp', images=pri_struct, format='vasp')   # 默认format为文件名后缀
# # view(pri_struct)
# ase.io.vasp.write_vasp("POSCAR4x4x4.vasp", pri_struct*(4,4,4), label='444supercell',direct=True,sort=True)   # 扩胞操作
a1 = read('POSCAR4x4x4.vasp')


""" 切表面 """
# Cutsurface = cut(a1, c=(0, 0, 1), origo=(0, 0, 0), nlayers=4)    # 切（0,0,1）面,保留4层
# Cutsurface.center(vacuum=10, axis=2)   # 加真空层,axis=0是x方向,1是y方向,z是c方向
# view(Cutsurface)   # alt+x/y/z  切换视图的方向
# print(Cutsurface.cell)
# print('--------------')
# print(Cutsurface.positions)
# print('--------------')
# print(Cutsurface.numbers)
# write_poscar(Cutsurface, '001.vasp')


Cutsurface = cut(a1, c=(1, 1, 1), origo=(0, 0, 0), nlayers=4)       # 切（1,1,1）面,保留4层
Cutsurface.center(vacuum=10, axis=2)       # 加真空层
# view(Cutsurface)
# print(Cutsurface.cell)
# write_poscar(Cutsurface, '111.vasp')
# write(filename='poscar3.vasp', images=Cutsurface)


# Cutsurface = cut(a1, c=(0, 1, 1), origo=(0, 0, 0), nlayers=4)    # 切（1,1,1）面,保留4层
# Cutsurface.center(vacuum=10, axis=2)   # 加真空层
# # view(Cutsurface)
# write_poscar(Cutsurface, '011.vasp')
