#!/usr/bin/env python
#md_output_onebyall.py
#Made by Jungmin Kim, SDlab

#find list to python
import glob
import os

print('='*60)
print('Made by Jkim. Made 180703 Ver 180703')

filelist_rdf=[] #List Declare
filelist_xyz=[]
command_rdf='plotrdf -a -np '
command_xyz='mdshak -r '
dirname_rdf='rdf_results'
dirname_xyz='xyz_results'

################ Dircheck % Mkdir ###############

if not os.path.isdir(dirname_rdf):
	os.mkdir(dirname_rdf)
if not os.path.isdir(dirname_xyz):
	os.mkdir(dirname_xyz)

################ filelist reading ###############

for file in glob.glob("*.mds"):
	filelist_xyz.append(file)
for file in glob.glob("*s_*.out"): # Read filelist_rdf
	filelist_rdf.append(file)

#################### plotrdf ####################

for ii in filelist_rdf:
	(temp,trash)=os.path.splitext(ii) # Remove extension
	command="%s%s.out > %s/%s.rdf" % (command_rdf,temp,dirname_rdf,temp)
	os.system(command)
	#print(command)
	print_temp="Current Operation is 'plotrdf' of '%s'" % (temp)
	print(print_temp)

#################### mdshak ####################

for ii in filelist_xyz:
	(temp,trash)=os.path.splitext(ii) # Remove extension
	command="%s%s.mds -f xyz -o %s/%s.xyz" % (command_xyz,temp,dirname_xyz,temp)
	os.system(command)
	#print(command)
	print_temp="Current Operation is 'mdshak' of '%s'" % (temp)
	print(print_temp)

print('='*60)

################# RDF_average.py ################

os.system("cd %s/" % dirname_rdf)
os.system("rdf_average_md_v2.py")
os.system("cd ../")
