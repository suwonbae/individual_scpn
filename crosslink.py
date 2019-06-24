# -*- coding: utf-8 -*-
"""
Created on Fri Jan  1 02:57:05 2016
This python script performs an algorithm which creates cross-links in a Coarse-Grained PE SCPN.
It has three parts.

1. The corresponding (to the desired cross-link percentage) number of beads become
active and the NVT dynamics simulation is performed for 5 ps.
 - A fully equilibrated linear CG PE SCPN model is processed by the python script 'identify.py'

2. Pairs of two active beads, both of which fall within the cature radius, are identified and
cross-links are created between the two active beads. After the cross-links are created,
the functionality of cross-linked active atoms is turned off.

3. As soon as cross-links are created, the spring constant and the equilibrium bond length are not
equal to the values at the equilibrium state. This part performs the multistep relaxation.
Through 5 steps, the parameters finally get to the real values.
Each step has a 50 ps NVT dynamics simulation.

###############################################################################
Parameters

ini_rc : Initial capture radius
k0 : Force constant of cross-link
r0 : Equilibrium bond length

target : Target percentage of cross-links as integer (ex. 5 % -> 5)

lammps : The full address of lammps executable

@author: suwonbae
"""

import numpy as np
import random as random
import math as math
import subprocess
import sys

#initial capture radius = ini_rc and cutoff radius for non-bonded potential
ini_rc=10.5
cutoff=10.5

#force constant and equilibrium bond length for new crosslinks
k0=350.0
r0=1.53

#lammps executable
lammps='ibrun lmp_knl -pk intel 0 -in'

#input for number of monomers and target percentage
variables=[]
for arg in sys.argv[1:]:
    variables.append(arg)

monomer=int(variables[0])
target=int(variables[1])

subprocess.Popen("python identify.py "+str(monomer)+" "+str(target), shell=True).wait()

file=open('in.crosslink.txt','w')
file.write('# Initialization\n')
file.write('units           real\n')
file.write('boundary        f f f\n')
file.write('atom_style      molecular\n')
file.write('read_data       data.input.txt extra/bond/types 2 extra/angle/types 1 extra/dihedral/types 1\n\n')
file.write('# Dreiding potential information\n')
file.write('neighbor        0.4 bin\n')
file.write('neigh_modify    delay 0 one 4000\n')
file.write('atom_modify     sort 0 0\n')
file.write('bond_style      harmonic\n')
file.write('bond_coeff      1 350.0000 1.53\n')
file.write('bond_coeff      2 %s %s\n' % (k0, r0))
file.write('bond_coeff      3 %s %s\n' % (k0, r0))
file.write('angle_style     harmonic\n')
file.write('angle_coeff     1 60.0000 109.50\n')
file.write('angle_coeff     2 60.0000 109.50\n')
file.write('dihedral_style  multi/harmonic\n')
file.write('dihedral_coeff  1 1.73 -4.49 0.776 6.99 0.0\n')
file.write('dihedral_coeff  2 1.73 -4.49 0.776 6.99 0.0\n')
file.write('pair_style      hybrid/overlay lj/cut 10.5\n')
file.write('pair_coeff      * * none\n')
file.write('pair_coeff      1 1 lj/cut 0.112 4.01\n')
file.write('pair_coeff      1 2 lj/cut 0.112 4.01\n')
file.write('pair_coeff      1 3 lj/cut 0.112 4.01\n')
file.write('pair_coeff      2 2 lj/cut 0.112 4.01\n')
file.write('pair_coeff      2 3 lj/cut 0.112 4.01\n')
file.write('pair_coeff      3 3 lj/cut 0.112 4.01\n\n')
file.write('comm_style      tiled\n\n')
file.write('group           chain type 1 2 3\n')
file.write('group           candi type 2\n')
file.write('compute         rg chain gyration\n')
file.write('compute         center chain com\n\n')
file.write('run             0\n')
file.write('variable        xlo_temp equal "bound(chain,xmin)-30"\n')
file.write('variable        xhi_temp equal "bound(chain,xmax)+30"\n')
file.write('variable        ylo_temp equal "bound(chain,ymin)-30"\n')
file.write('variable        yhi_temp equal "bound(chain,ymax)+30"\n')
file.write('variable        zlo_temp equal "bound(chain,zmin)-30"\n')
file.write('variable        zhi_temp equal "bound(chain,zmax)+30"\n')
file.write('variable        xlo_init equal ${xlo_temp}\n')
file.write('variable        xhi_init equal ${xhi_temp}\n')
file.write('variable        ylo_init equal ${ylo_temp}\n')
file.write('variable        yhi_init equal ${yhi_temp}\n')
file.write('variable        zlo_init equal ${zlo_temp}\n')
file.write('variable        zhi_init equal ${zhi_temp}\n')
file.write('change_box      all x final ${xlo_init} ${xhi_init} y final ${ylo_init} ${yhi_init} z final ${zlo_init} ${zhi_init} boundary f f f\n')
file.write('#####################################################\n')
file.write('# Equilibration (NVT at 300 K)\n')
file.write('reset_timestep  0\n')
file.write('fix             bal all balance 500 1.1 rcb out balancing\n')
file.write('fix             1 all nvt temp 300.0 300.0 100.0\n')
file.write('fix             2 all momentum 1 linear 1 1 1 angular\n')
file.write('thermo_style    custom step temp press etotal pe epair ebond eangle edihed lx ly lz c_rg c_center[1] c_center[2] c_center[3]\n')
file.write('thermo          100\n')
file.write('run_style       verlet\n')
file.write('run             1000000\n')
file.write('unfix           1\n')
file.write('unfix           2\n\n')
file.write('#####################################################\n')
file.write('# Equilibration (NVT at 300 K -> 500 K)\n')
file.write('reset_timestep  0\n')
file.write('velocity        all zero linear\n')
file.write('fix             1 all nvt temp 300.0 500.0 100.0\n')
file.write('fix             2 all momentum 1 linear 1 1 1 angular\n')
file.write('thermo_style    custom step temp press etotal pe epair ebond eangle edihed lx ly lz c_rg c_center[1] c_center[2] c_center[3]\n')
file.write('thermo          100\n')
file.write('run_style       verlet\n')
file.write('run             800000\n')
file.write('unfix           1\n')
file.write('unfix           2\n\n')
file.write('#####################################################\n')
file.write('# Equilibration (NVT at 500 K)\n')
file.write('reset_timestep  0\n')
file.write('velocity        all zero linear\n')
file.write('fix             1 all nvt temp 500.0 500.0 100.0\n')
file.write('fix             2 all momentum 1 linear 1 1 1 angular\n')
file.write('thermo_style    custom step temp press etotal pe epair ebond eangle edihed lx ly lz c_rg c_center[1] c_center[2] c_center[3]\n')
file.write('thermo          100\n')
file.write('run_style       verlet\n')
file.write('run             1000000\n')
file.write('unfix           1\n')
file.write('unfix           2\n')
file.write('unfix           bal\n\n')
file.write('write_data      data.input.txt\n\n')
file.write('print "All done"')
file.close()
subprocess.Popen(lammps+' in.crosslink.txt', shell=True).wait()

file=open('data.input.txt','r')
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
file.close()
def getKey(item):
    return item[0]

atoms=sorted(atoms, key=getKey)

active=0
linked=0
for i in range(len(atoms)):
    if atoms[i][2]==2:
        active+=1
    if atoms[i][2]==3:
        linked+=1

percentage=float(linked)/monomer

f=open('data.input.txt','w')
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

#endflag is initialized to 0
endflag=0
while (endflag==0):
    file=open('data.input.txt','r')
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
    file.close()
    
    def getKey(item):
        return item[0]
    
    atoms=sorted(atoms, key=getKey)
    
    #'fn_atoms' refers to atoms with a functionality
    fn_atoms=[]
    fn_atoms_num=0
    for i in range(atomnum):
        if atoms[i][2]==2:
            fn_atoms.append(atoms[i])
            fn_atoms_num+=1
    fn_atoms_memory=fn_atoms
    
    #'crosslink' is the list of new bonds
    crosslink=[]
    crosslink_num=0
    while fn_atoms_num > 0:
        c1=random.randrange(fn_atoms_num)
        ID=[row[0] for row in fn_atoms]
        molID=[row[1] for row in fn_atoms]
        atype=[row[2] for row in fn_atoms]
        coord=[[row[3] for row in fn_atoms],[row[4] for row in fn_atoms],[row[5] for row in fn_atoms]]
        
        d=[]
        for i in range(fn_atoms_num):
            distance=math.sqrt(math.pow(coord[0][i]-coord[0][c1],2)+math.pow(coord[1][i]-coord[1][c1],2)+math.pow(coord[2][i]-coord[2][c1],2))
            d.append(distance)
            
        candi=[]
        candi_d=[]
        for i in range(fn_atoms_num):
            if d[i]!=0 and d[i]<cutoff:
                candi.append(i)
                candi_d.append(d[i])
        
        if len(candi)==0:
            break
        else:
            rn=random.randrange(len(candi))
            c2=candi[rn]
            crosslink_num+=1
            
            rc=candi_d[rn]
            
            if c1 > c2:
                temp=c1
                c1=c2
                c2=temp
            
            if ID[c1] > ID[c2]:
                crosslink.append([bondnum+crosslink_num,3,ID[c2],ID[c1]])
            else:
                crosslink.append([bondnum+crosslink_num,3,ID[c1],ID[c2]])
                
            fn_atoms=fn_atoms[0:c1]+fn_atoms[c1+1:c2]+fn_atoms[c2+1:atomnum]
            fn_atoms_num=0
            
    if crosslink_num==0:
        active=0
        linked=0
        for i in range(len(atoms)):
            if atoms[i][2]==2:
                active+=1
            if atoms[i][2]==3:
                linked+=1
        
        percentage=float(linked)/monomer
        
        f=open('data.input.txt', 'w')
        f.write('LAMMPS data file, percentage = %s\n\n' % percentage)
        f.write('%s atoms\n' % atomnum)
        f.write('%s atom types\n' % atomtypes)
        f.write('%s bonds\n' % (bondnum))
        f.write('%s bond types\n' % bondtypes)
        f.write('%s angles\n' % (anglenum))
        f.write('%s angle types\n' % angletypes)
        f.write('%s dihedrals\n' % (dihedralnum))
        f.write('%s dihedral types\n\n' % dihedraltypes)
        f.write('%s %s xlo xhi\n' % (xlo, xhi))
        f.write('%s %s ylo yhi\n' % (ylo, yhi))
        f.write('%s %s zlo zhi\n\n' % (zlo, zhi))
        f.write('Masses\n\n')
        for line in masses:
            f.write('%s %s\n' % (int(line[0]), line[1]))
        f.write('\n')
        f.write('Atoms\n\n')
        for line in range(atomnum):
            f.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (atoms[line][0], atoms[line][1], atoms[line][2], atoms[line][3], atoms[line][4], atoms[line][5]))
        f.write('\n')
        f.write('Velocities\n\n')
        for line in vels:
            f.write('%s\t%s\t%s\t%s\n' % (line[0], line[1], line[2], line[3]))
        f.write('\n')
        f.write('Bonds\n\n')
        for line in range(len(bonds)):
            f.write('%s\t%s\t%s\t%s\n' % (bonds[line][0], bonds[line][1], bonds[line][2], bonds[line][3]))
        f.write('\n')
        f.write('Angles\n\n')
        for line in range(len(angles)):
            f.write('%s\t%s\t%s\t%s\t%s\n' % (angles[line][0], angles[line][1], angles[line][2], angles[line][3], angles[line][4]))
        f.write('\n')
        f.write('Dihedrals\n\n')
        for line in range(len(dihedrals)):
            f.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (dihedrals[line][0], dihedrals[line][1], dihedrals[line][2], dihedrals[line][3], dihedrals[line][4], dihedrals[line][5]))
        f.close()
        
        file=open('in.equil.txt','w')
        file.write('# Initialization\n')
        file.write('units           real\n')
        file.write('boundary        f f f\n')
        file.write('atom_style      molecular\n')
        file.write('read_data       data.input.txt\n\n')
        file.write('# Dreiding potential information\n')
        file.write('neighbor        0.4 bin\n')
        file.write('neigh_modify    delay 0 one 4000\n')
        file.write('atom_modify     sort 0 0\n')
        file.write('bond_style      harmonic\n')
        file.write('bond_coeff      1 350.0000 1.53\n')
        file.write('bond_coeff      2 %s %s\n' % (k0, r0))
        file.write('bond_coeff      3 %s %s\n' % (k0, r0))
        file.write('angle_style     harmonic\n')
        file.write('angle_coeff     1 60.0000 109.50\n')
        file.write('angle_coeff     2 60.0000 109.50\n')
        file.write('dihedral_style  multi/harmonic\n')
        file.write('dihedral_coeff  1 1.73 -4.49 0.776 6.99 0.0\n')
        file.write('dihedral_coeff  2 1.73 -4.49 0.776 6.99 0.0\n')
        file.write('pair_style      hybrid/overlay lj/cut 10.5\n')
        file.write('pair_coeff      * * none\n')
        file.write('pair_coeff      1 1 lj/cut 0.112 4.01\n')
        file.write('pair_coeff      1 2 lj/cut 0.112 4.01\n')
        file.write('pair_coeff      1 3 lj/cut 0.112 4.01\n')
        file.write('pair_coeff      2 2 lj/cut 0.112 4.01\n')
        file.write('pair_coeff      2 3 lj/cut 0.112 4.01\n')
        file.write('pair_coeff      3 3 lj/cut 0.112 4.01\n\n')
        file.write('comm_style      tiled\n\n')
        file.write('group           chain type 1 2 3\n')
        file.write('group           candi type 2\n')
        file.write('compute         rg chain gyration\n')
        file.write('compute         center chain com\n\n')
        file.write('#####################################################\n')
        file.write('# Equilibration (NVT at 500 K)\n')
        file.write('reset_timestep  0\n')
        file.write('fix             bal all balance 500 1.1 rcb out balancing\n')
        file.write('fix             1 all nvt temp 500.0 500.0 100.0\n')
        file.write('fix             2 all momentum 1 linear 1 1 1 angular\n')
        file.write('thermo_style    custom step temp press etotal pe epair ebond eangle edihed lx ly lz c_rg c_center[1] c_center[2] c_center[3]\n')
        file.write('thermo          100\n')
        file.write('run_style       verlet\n')
        file.write('run             100000\n')
        file.write('unfix           1\n')
        file.write('unfix           2\n')
        file.write('unfix           bal\n\n')
        file.write('write_data      data.linked.txt\n\n')
        file.write('print "All done"')
        file.close()
        subprocess.Popen(lammps+' in.equil.txt', shell=True).wait()

        file=open('data.linked.txt','r')
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
        file.close()
        
        def getKey(item):
            return item[0]
        
        atoms=sorted(atoms, key=getKey)
        
        active=0
        linked=0
        for i in range(len(atoms)):
            if atoms[i][2]==2:
                active+=1
            if atoms[i][2]==3:
                linked+=1
        
        for i in range(len(bonds)):
            if bonds[i][1]==3:
                bonds[i][1]=2
        
        percentage=float(linked)/monomer  
        
        f=open('data.input.txt','w')
        f.write('LAMMPS data file percentage = %s\n\n' % percentage)
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
    else:
        for i in range(crosslink_num):
            atoms[crosslink[i][2]-1][2]=3
            atoms[crosslink[i][3]-1][2]=3
        
        bonds_memory=bonds
        crosslink=np.asarray(crosslink)
        bonds=np.r_[bonds,crosslink]
        
        #builds new angles
        newangle=[]
        newangle_num=0
        for i in range(crosslink_num):
            center=crosslink[i][2]
            branch=crosslink[i][3]
            
            ct=[]
            for e in range(bondnum):
                if bonds_memory[e][2]==center:
                    ct.append(e)
            if len(ct)!=0:
                for j in range(len(ct)):
                    temp=[branch, center, bonds_memory[ct[j]][3]]
                    if len(temp)==3:
                        newangle_num+=1
                        newangle.append([anglenum+newangle_num,2,branch,center,bonds_memory[ct[j]][3]])
        
            ct=[]
            for e in range(bondnum):
                if bonds_memory[e][3]==center:
                    ct.append(e)
            if len(ct)!=0:
                for j in range(len(ct)):
                    temp=[branch, center, bonds_memory[ct[j]][2]]
                    if len(temp)==3:
                        newangle_num+=1
                        newangle.append([anglenum+newangle_num,2,branch,center,bonds_memory[ct[j]][2]])
        
        for i in range(crosslink_num):
            center=crosslink[i][3]
            branch=crosslink[i][2]
            
            ct=[]
            for e in range(bondnum):
                if bonds_memory[e][3]==center:
                    ct.append(e)
            if len(ct)!=0:
                for j in range(len(ct)):
                    temp=[branch, center, bonds_memory[ct[j]][2]]
                    if len(temp)==3:
                        newangle_num+=1
                        newangle.append([anglenum+newangle_num,2,branch,center,bonds_memory[ct[j]][2]])
        
            ct=[]
            for e in range(bondnum):
                if bonds_memory[e][2]==center:
                    ct.append(e)
            if len(ct)!=0:
                for j in range(len(ct)):
                    temp=[branch, center, bonds_memory[ct[j]][3]]
                    if len(temp)==3:
                        newangle_num+=1
                        newangle.append([anglenum+newangle_num,2,branch,center,bonds_memory[ct[j]][3]])
        
        angles_memory=angles
        newangle=np.asarray(newangle)
        angles=np.r_[angles,newangle]
        
        #builds new dihedrals
        newdihedral=[]
        newdihedral_num=0
        temp=[]
        for i in range(crosslink_num):
            center1=crosslink[i][2]
            center2=crosslink[i][3]
            
            for e in range(bondnum):
                if bonds_memory[e][2]==center1:
                    branch_new1=bonds[e][3]
                if bonds_memory[e][3]==center2:
                    branch_new2=bonds[e][2]
                elif bonds_memory[e][2]==center2:
                    branch_new2=bonds[e][3]
            comb=[branch_new2, center2, center1, branch_new1]
            if len(comb)==4:
                comb.sort()
                if not comb in temp:
                    temp.append(comb)
                    newdihedral_num+=1
                    newdihedral.append([dihedralnum+newdihedral_num,2,branch_new2,center2,center1,branch_new1])
            
            for e in range(bondnum):
                if bonds_memory[e][3]==center1:
                    branch_new1=bonds[e][2]
                if bonds_memory[e][2]==center2:
                    branch_new2=bonds[e][3]
                elif bonds_memory[e][3]==center2:
                    branch_new2=bonds[e][2]
            comb=[branch_new2, center2, center1, branch_new1]
            if len(comb)==4:
                comb.sort()
                if not comb in temp:
                    temp.append(comb)
                    newdihedral_num+=1
                    newdihedral.append([dihedralnum+newdihedral_num,2,branch_new2,center2,center1,branch_new1])
            
        for i in range(newangle_num):
            end1=newangle[i][2]
            center=newangle[i][3]
            end2=newangle[i][4]
    
            ct=[]
            for e in range(bondnum):
                if bonds_memory[e][2]==end1:
                    ct.append(e)
            if len(ct)!=0:
                for j in range(len(ct)):
                    if bonds_memory[ct[j]][3]!=center:
                        comb=[bonds_memory[ct[j]][3],end1,center,end2]
                        if len(comb)==4:
                            comb.sort()
                            if not comb in temp:
                                temp.append(comb)
                                newdihedral_num+=1
                                newdihedral.append([dihedralnum+newdihedral_num,2,bonds_memory[ct[j]][3],end1,center,end2])
    
            ct=[]
            for e in range(bondnum):
                if bonds_memory[e][3]==end1:
                    ct.append(e)
            if len(ct)!=0:
                for j in range(len(ct)):
                    if bonds_memory[ct[j]][2]!=center:
                        comb=[bonds_memory[ct[j]][2],end1,center,end2]
                        if len(comb)==4:
                            comb.sort()
                            if not comb in temp:
                                temp.append(comb)
                                newdihedral_num+=1
                                newdihedral.append([dihedralnum+newdihedral_num,2,bonds_memory[ct[j]][2],end1,center,end2])
                        
        for i in range(newangle_num):
            end1=newangle[i][4]
            center=newangle[i][3]
            end2=newangle[i][2]
    
            ct=[]
            for e in range(bondnum):
                if bonds_memory[e][3]==end1:
                    ct.append(e)
            if len(ct)!=0:
                for j in range(len(ct)):
                    if bonds_memory[ct[j]][2]!=center:
                        comb=[bonds_memory[ct[j]][2],end1,center,end2]
                        if len(comb)==4:
                            comb.sort()
                            if not comb in temp:
                                temp.append(comb)
                                newdihedral_num+=1
                                newdihedral.append([dihedralnum+newdihedral_num,2,bonds_memory[ct[j]][2],end1,center,end2])
    
            ct=[]
            for e in range(bondnum):
                if bonds_memory[e][2]==end1:
                    ct.append(e)
            if len(ct)!=0:
                for j in range(len(ct)):
                    if bonds_memory[ct[j]][3]!=center:
                        comb=[bonds_memory[ct[j]][3],end1,center,end2]
                        if len(comb)==4:
                            comb.sort()
                            if not comb in temp:
                                temp.append(comb)
                                newdihedral_num+=1
                                newdihedral.append([dihedralnum+newdihedral_num,2,bonds_memory[ct[j]][3],end1,center,end2])
        
        dihedrals_memory=dihedrals
        newdihedral=np.array(newdihedral)
        dihedrals=np.r_[dihedrals,newdihedral]
        
        f=open('data.input.txt', 'w')
        f.write('LAMMPS data file\n\n')
        f.write('%s atoms\n' % atomnum)
        f.write('%s atom types\n' % atomtypes)
        f.write('%s bonds\n' % (bondnum+crosslink_num))
        f.write('%s bond types\n' % bondtypes)
        f.write('%s angles\n' % (anglenum+newangle_num))
        f.write('%s angle types\n' % angletypes)
        f.write('%s dihedrals\n' % (dihedralnum+newdihedral_num))
        f.write('%s dihedral types\n\n' % dihedraltypes)
        f.write('%s %s xlo xhi\n' % (xlo, xhi))
        f.write('%s %s ylo yhi\n' % (ylo, yhi))
        f.write('%s %s zlo zhi\n\n' % (zlo, zhi))
        f.write('Masses\n\n')
        for line in masses:
            f.write('%s %s\n' % (int(line[0]), line[1]))
        f.write('\n')
        f.write('Atoms\n\n')
        for line in range(atomnum):
            f.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (atoms[line][0], atoms[line][1], atoms[line][2], atoms[line][3], atoms[line][4], atoms[line][5]))
        f.write('\n')
        f.write('Velocities\n\n')
        for line in vels:
            f.write('%s\t%s\t%s\t%s\n' % (line[0], line[1], line[2], line[3]))
        f.write('\n')
        f.write('Bonds\n\n')
        for line in range(len(bonds)):
            f.write('%s\t%s\t%s\t%s\n' % (bonds[line][0], bonds[line][1], bonds[line][2], bonds[line][3]))
        f.write('\n')
        f.write('Angles\n\n')
        for line in range(len(angles)):
            f.write('%s\t%s\t%s\t%s\t%s\n' % (angles[line][0], angles[line][1], angles[line][2], angles[line][3], angles[line][4]))
        f.write('\n')
        f.write('Dihedrals\n\n')
        for line in range(len(dihedrals)):
            f.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (dihedrals[line][0], dihedrals[line][1], dihedrals[line][2], dihedrals[line][3], dihedrals[line][4], dihedrals[line][5]))
        f.close()
    
        file=open('in.crosslink_temp.txt','w')
        file.write('# Initialization\n')
        file.write('units           real\n')
        file.write('boundary        f f f\n')
        file.write('atom_style      molecular\n')
        file.write('read_data       data.input.txt\n\n')
        file.write('# Dreiding potential information\n')
        file.write('neighbor        0.4 bin\n')
        file.write('neigh_modify    delay 0 one 4000\n')
        file.write('atom_modify     sort 0 0\n')
        file.write('bond_style      harmonic\n')
        file.write('bond_coeff      1 350.0000 1.53\n')
        file.write('bond_coeff      2 %s %s\n' % (k0, r0))
        file.write('bond_coeff      3 %s %s\n' % (k0/5, rc))
        file.write('angle_style     harmonic\n')
        file.write('angle_coeff     1 60.0000 109.50\n')
        file.write('angle_coeff     2 60.0000 109.50\n')
        file.write('dihedral_style  multi/harmonic\n')
        file.write('dihedral_coeff  1 1.73 -4.49 0.776 6.99 0.0\n')
        file.write('dihedral_coeff  2 1.73 -4.49 0.776 6.99 0.0\n')
        file.write('pair_style      hybrid/overlay lj/cut 10.5\n')
        file.write('pair_coeff      * * none\n')
        file.write('pair_coeff      1 1 lj/cut 0.112 4.01\n')
        file.write('pair_coeff      1 2 lj/cut 0.112 4.01\n')
        file.write('pair_coeff      1 3 lj/cut 0.112 4.01\n')
        file.write('pair_coeff      2 2 lj/cut 0.112 4.01\n')
        file.write('pair_coeff      2 3 lj/cut 0.112 4.01\n')
        file.write('pair_coeff      3 3 lj/cut 0.112 4.01\n\n')
        file.write('comm_style      tiled\n\n')
        file.write('group           chain type 1 2 3\n')
        file.write('group           candi type 2\n')
        file.write('compute         rg chain gyration\n')
        file.write('compute         center chain com\n\n')
        file.write('#####################################################\n')
        file.write('# Equilibration (NVT at 500 K)\n')
        file.write('reset_timestep  0\n')
	file.write('dump            1 all xyz 100000000 look_%s.xyz\n' % percentage)
        file.write('fix             bal all balance 500 1.1 rcb out balancing\n\n')
        file.write('# Multistep relaxation\n')
        for i in range(5):
            file.write('# step %s\n' % (i+1))
            file.write('bond_coeff      3 %s %s\n' % (k0*(i+1)/5, r0+(rc-r0)*(4-i)/4))
            file.write('fix             1 all nvt temp 500.0 500.0 100.0\n')
            file.write('fix             2 all momentum 1 linear 1 1 1 angular\n')
            file.write('thermo_style    custom step temp press etotal pe epair ebond eangle edihed lx ly lz c_rg c_center[1] c_center[2] c_center[3]\n')
            file.write('thermo          100\n')
            file.write('timestep        0.5\n')
            file.write('run_style       verlet\n')
            file.write('run             200000\n')
            file.write('unfix           1\n')
            file.write('unfix           2\n\n')
        file.write('undump          1\n')
        file.write('write_data      data.linked.txt\n\n')
        file.write('print "All done"')
        file.close()
        subprocess.Popen(lammps+' in.crosslink_temp.txt', shell=True).wait()

        file=open('data.linked.txt','r')
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
        file.close()
        
        def getKey(item):
            return item[0]
        
        atoms=sorted(atoms, key=getKey)
        
        active=0
        linked=0
        for i in range(len(atoms)):
            if atoms[i][2]==2:
                active+=1
            if atoms[i][2]==3:
                linked+=1
        
        for i in range(len(bonds)):
            if bonds[i][1]==3:
                bonds[i][1]=2
        
        percentage=float(linked)/monomer
        
        if ((percentage*100)>=(target-0.5)):
            endflag=1        
        
        f=open('data.input.txt','w')
        f.write('LAMMPS data file percentage = %s\n\n' % percentage)
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

file=open('in.anneal.txt','w')
file.write('# Initialization\n')
file.write('units           real\n')
file.write('boundary        f f f\n')
file.write('atom_style      molecular\n')
file.write('read_data       data.input.txt\n\n')
file.write('# Dreiding potential information\n')
file.write('neighbor        0.4 bin\n')
file.write('neigh_modify    delay 0 one 4000\n')
file.write('atom_modify     sort 0 0\n')
file.write('bond_style      harmonic\n')
file.write('bond_coeff      1 350.0000 1.53\n')
file.write('bond_coeff      2 %s %s\n' % (k0, r0))
file.write('bond_coeff      3 %s %s\n' % (k0/5, rc))
file.write('angle_style     harmonic\n')
file.write('angle_coeff     1 60.0000 109.50\n')
file.write('angle_coeff     2 60.0000 109.50\n')
file.write('dihedral_style  multi/harmonic\n')
file.write('dihedral_coeff  1 1.73 -4.49 0.776 6.99 0.0\n')
file.write('dihedral_coeff  2 1.73 -4.49 0.776 6.99 0.0\n')
file.write('pair_style      hybrid/overlay lj/cut 10.5\n')
file.write('pair_coeff      * * none\n')
file.write('pair_coeff      1 1 lj/cut 0.112 4.01\n')
file.write('pair_coeff      1 2 lj/cut 0.112 4.01\n')
file.write('pair_coeff      1 3 lj/cut 0.112 4.01\n')
file.write('pair_coeff      2 2 lj/cut 0.112 4.01\n')
file.write('pair_coeff      2 3 lj/cut 0.112 4.01\n')
file.write('pair_coeff      3 3 lj/cut 0.112 4.01\n\n')
file.write('comm_style      tiled\n\n')
file.write('group           chain type 1 2 3\n')
file.write('group           candi type 2\n')
file.write('compute         rg chain gyration\n')
file.write('compute         center chain com\n\n')
file.write('#####################################################\n')
file.write('# Simulated Annealing Stage 1 (NVT at 500 K)\n')
file.write('reset_timestep  0\n')
file.write('dump            1 all xyz 100000000 look_%s.xyz\n' % percentage)
file.write('fix             bal all balance 500 1.1 rcb out balancing\n')
file.write('fix             1 all nvt temp 500.0 500.0 100.0\n')
file.write('fix             2 all momentum 1 linear 1 1 1 angular\n')
file.write('thermo_style    custom step temp press etotal pe epair ebond eangle edihed lx ly lz c_rg c_center[1] c_center[2] c_center[3]\n')
file.write('thermo          100\n')
file.write('timestep        1\n')
file.write('run_style       verlet\n')
file.write('run             1000000\n')
file.write('unfix           1\n')
file.write('unfix           2\n')
file.write('undump          1\n\n')
file.write('#####################################################\n')
file.write('# Simulated Annealing Stage 2 (NVT at 500 K -> 300 K)\n')
file.write('fix             1 all nvt temp 500.0 300.0 100.0\n')
file.write('fix             2 all momentum 1 linear 1 1 1 angular\n')
file.write('thermo_style    custom step temp press etotal pe epair ebond eangle edihed lx ly lz c_rg c_center[1] c_center[2] c_center[3]\n')
file.write('thermo          100\n')
file.write('timestep        1\n')
file.write('run_style       verlet\n')
file.write('run             800000\n')
file.write('unfix           1\n')
file.write('unfix           2\n\n')
file.write('#####################################################\n')
file.write('# Simulated Annealing Stage 3 (NVT at 300 K)\n')
file.write('fix             1 all nvt temp 300.0 300.0 100.0\n')
file.write('fix             2 all momentum 1 linear 1 1 1 angular\n')
file.write('thermo_style    custom step temp press etotal pe epair ebond eangle edihed lx ly lz c_rg c_center[1] c_center[2] c_center[3]\n')
file.write('thermo          100\n')
file.write('timestep        1\n')
file.write('run_style       verlet\n')
file.write('run             1000000\n')
file.write('unfix           1\n')
file.write('unfix           2\n')
file.write('unfix           bal\n\n')
file.write('write_data      data.linked.txt\n\n')
file.write('print "All done"')
file.close()
subprocess.Popen(lammps+' in.anneal.txt', shell=True).wait()

file=open('data.linked.txt','r')
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
file.close()

def getKey(item):
    return item[0]

atoms=sorted(atoms, key=getKey)

active=0
linked=0
for i in range(len(atoms)):
    if atoms[i][2]==2:
        active+=1
    if atoms[i][2]==3:
        linked+=1

for i in range(len(bonds)):
    if bonds[i][1]==3:
        bonds[i][1]=2

percentage=float(linked)/monomer  

f=open('data.input.txt','w')
f.write('LAMMPS data file percentage = %s\n\n' % percentage)
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

subprocess.Popen('mv data.input.txt data.'+str(target)+'crosslinked.txt', shell=True).wait()
subprocess.Popen("rm data.linked.txt", shell=True).wait()
