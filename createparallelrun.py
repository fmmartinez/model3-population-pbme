#!/usr/bin/python
import os
import shutil

nproc = 10

#basis functions used
g = 1  
b = 1  
d = 1  

#delta
delta = 1

#trajectories per processor
tpp = 5000

dirs = []

l = []
m = []

if (nproc > 9):
	for i in range(0,10):
		name = 'map-0' + str(i)
		dirs.append(name)
	for i in range(10,nproc):
		name = 'map-'  + str(i)
		dirs.append(name)
else:
	for i in range(0,nproc):
		name = 'map-'  + str(i)
		dirs.append(name)

#generate folders
for i in range(0,nproc):
	if not(os.path.exists(dirs[i])):
			os.mkdir(dirs[i])

#original structure of md.in file in map00:
l.append('Np\tDELTA\tNOSC\tOME_MAX\t\r\n')
l.append('1\t' + str(delta) + 'd0\t20\t50\r\n')
l.append('NMCS\tNMDS\tseed\tDT\tLAMDA_D\r\n')
l.append('25000\t2000\t5\t5d-5\t10d0\r\n')
l.append('Eg\tEb\tEd\tmu\tE0\tbeta\tvib_omega\r\n')
l.append('0\t240\t240\t1\t0.9\t0.24\t37.7d0\r\n')
l.append('TIME_J\tTAU_J\tOMEGA_J\r\n')
l.append('0.15d0\t0.045d0\t260d0\r\n')
l.append('BATH(0:B EQ 1:T EQ)\tINIT\r\n')
l.append('0\t3\r\n')
l.append('basispercenter\tmapinG\tmapinB\tmapinD\tcount\r\n')
l.append('3\t' + str(g) + '\t' + str(b)+ '\t' + str(d) + '\t0\r\n')

#creating md.in files
for i in range(0,nproc):
	mdfile = open('./' + dirs[i] + '/md.in','w')
	
	seed = i*9 + 5

	l[3] = str(tpp) + '\t4000\t' + str(seed) + '\t5d-5\t10d0\r\n'

	for j in range(0,12):
		mdfile.write(l[j])

	mdfile.close()

#original structure of pbs fie
m.append('#!/bin/bash -l\n')
m.append('#PBS -r n\n')
m.append('#PBS -N pbme\n')
m.append('#PBS -l walltime=10:00:00\n')
m.append('#PBS -l procs=1\n')
m.append('#PBS -l mem=256mb\n')
m.append('#PBS -m bea\n')
m.append('#PBS -M fmmartin@ualberta.ca\n')
m.append('cd $PBS_O_WORKDIR\n')
m.append('module load compiler/intel/12.1\n')
m.append('time ./a.out < md.in\n')

#creating pbs files
for i in range(0,nproc):
	pbsfile = open('./' + dirs[i] + '/submit_cluster.pbs','w')
	
	m[2] = '#PBS -N p' + str(g) + str(b) + str(d) + '-' + str(i) + '\n'
	
	if ((i == 0)|(i == nproc-1)):
		m[6] = '#PBS -m bea\n'
	else:
		m[6] = '#PBS -m n\n'

	for j in range(0,11):
		pbsfile.write(m[j])

	pbsfile.close()


#copy executables
for i in range(0,nproc):
	shutil.copy2('a.out',dirs[i])
#	shutil.copy2('m_vib.mod',dirs[i])
