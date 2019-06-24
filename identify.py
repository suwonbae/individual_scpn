# -*- coding: utf-8 -*-
"""
Created on Wed Mar  9 01:15:49 2016
This python script activates atoms so that they can be cross-linked.
Atom of type 2 is active atom.
Atom of type 3 is cross-linked atom.
@author: suwonbae
"""

import subprocess
import random
import sys

variables=[]
for arg in sys.argv[1:]:
    variables.append(arg)

monomer=int(variables[0])
target=int(variables[1])

file=open('data.cg'+str(monomer)+'scpn.txt','r')

percentnum=int(monomer*target*0.01)

flag=0
while flag==0:
    line=file.readline()
    if 'atoms' in line:
        atomnum=int(line.split()[0])
    if 'bonds' in line:
        bondnum=int(line.split()[0])
    if 'angles' in line:
        anglenum=int(line.split()[0])
    if 'dihedrals' in line:
        dihedralnum=int(line.split()[0])
    if 'atom types' in line:
        atomtypes=int(line.split()[0])
    if 'bond types' in line:
        bondtypes=int(line.split()[0])
    if 'angle types' in line:
        angletypes=int(line.split()[0])
    if 'dihedral types' in line:
        dihedraltypes=int(line.split()[0])
    if 'xlo xhi' in line:
        xlo, xhi=float(line.split()[0]), float(line.split()[1])
    if 'ylo yhi' in line:
        ylo, yhi=float(line.split()[0]), float(line.split()[1])     
    if 'zlo zhi' in line:
        zlo, zhi=float(line.split()[0]), float(line.split()[1])
        flag=1

masses=[]
atoms=[]
vels=[]
bonds=[]
angles=[]
dihedrals=[]
while flag==1:
    line=file.readline()
    if 'Masses' in line:
        line=file.readline()
        for i in range(atomtypes):
            line=file.readline()
            masses.append([int(line.split()[0]),float(line.split()[1])])
    if 'Atoms' in line:
        line=file.readline()
        for i in range(atomnum):
            line=file.readline()
            atoms.append([int(line.split()[0]),int(line.split()[1]),int(line.split()[2]),float(line.split()[3]),float(line.split()[4]),float(line.split()[5])])
    if 'Velocities' in line:
        line=file.readline()
        for i in range(atomnum):
            line=file.readline()
            vels.append([int(line.split()[0]),float(line.split()[1]),float(line.split()[2]),float(line.split()[3])])
    if 'Bonds' in line:
        line=file.readline()
        for i in range(bondnum):
            line=file.readline()
            bonds.append([int(line.split()[0]),int(line.split()[1]),int(line.split()[2]),int(line.split()[3])])
    if 'Angles' in line:
        line=file.readline()
        for i in range(anglenum):
            line=file.readline()
            angles.append([int(line.split()[0]),int(line.split()[1]),int(line.split()[2]),int(line.split()[3]),int(line.split()[4])])
    if 'Dihedrals' in line:
        line=file.readline()
        for i in range(dihedralnum):
            line=file.readline()
            dihedrals.append([int(line.split()[0]),int(line.split()[1]),int(line.split()[2]),int(line.split()[3]),int(line.split()[4]),int(line.split()[5])])
            flag=2

def getKey(item):
    return item[0]

atoms=sorted(atoms, key=getKey)

atomtypes+=2
masses.append([2, 14.02])
masses.append([3, 14.02])

pool=range(len(atoms))

randlist=[]
for i in range(percentnum):
    temp=random.sample(pool,1)[0]
    randlist.append(temp)
    if temp==0:
        pool.remove(temp)
        for j in range(4):
            if (temp+(j+1)) in pool:
                pool.remove(temp+(j+1))
    elif temp==pool[len(pool)-1]:
        pool.remove(temp)
        for j in range(4):
            if (temp-(j+1)) in pool:
                pool.remove(temp-(j+1))
    else:
        pool.remove(temp)
        for j in range(4):
            if (temp-(j+1)) in pool:
                pool.remove(temp-(j+1))
            if (temp+(j+1)) in pool:
                pool.remove(temp+(j+1))

for i in range(len(randlist)):
    atoms[randlist[i]][2]=2

f=open('data.prepared.txt','w')
f.write('LAMMPS data file\n\n')
f.write('%s atoms\n' % (atomnum))
f.write('%s atom types\n' % (atomtypes))
f.write('%s bonds\n' % (bondnum))
f.write('%s bond types\n' % (bondtypes))
f.write('%s angles\n' % (anglenum))
f.write('%s angle types\n' % (angletypes))
f.write('%s dihedrals\n' % (dihedralnum))
f.write('%s dihedral types\n\n' % (dihedraltypes))
f.write('10 extra bond per atom\n')
f.write('10 extra angle per atom\n')
f.write('10 extra dihedral per atom\n')
f.write('10 extra special per atom\n\n')
f.write('%s %s xlo xhi\n' % (xlo, xhi))
f.write('%s %s ylo yhi\n' % (ylo, yhi))
f.write('%s %s zlo zhi\n\n' % (zlo, zhi))
f.write('Masses\n\n')
for line in masses:
    f.write('%s %s\n' % (line[0], line[1]))
f.write('\n')
f.write('Atoms\n\n')
for line in atoms:
    f.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (line[0], line[1], line[2], line[3], line[4], line[5]))
f.write('\n')
f.write('Velocities\n\n')
for line in vels:
    f.write('%s\t%s\t%s\t%s\n' % (line[0], line[1], line[2], line[3]))
f.write('\n')
f.write('Bonds\n\n')
for line in bonds:
    f.write('%s\t%s\t%s\t%s\n' % (line[0], line[1], line[2], line[3]))
f.write('\n')
f.write('Angles\n\n')
for line in angles:
    f.write('%s\t%s\t%s\t%s\t%s\n' % (line[0], line[1], line[2], line[3], line[4]))
f.write('\n')
f.write('Dihedrals\n\n')
for line in dihedrals:
    f.write('%s\t%s\t%s\t%s\t%s\t%s' % (line[0], line[1], line[2], line[3], line[4], line[5]))
    if line[0]!=dihedralnum:
        f.write('\n')
f.close()

subprocess.Popen("cp data.prepared.txt data.input.txt", shell=True).wait()
subprocess.Popen("rm data.linked.txt", shell=True).wait()