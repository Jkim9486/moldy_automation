#!/usr/bin/env python
#RDF_average_Jkim.py
#<RDF_average_Hki.cpp> Code made by Dr. Hosung Ki, SDlab
#RDF_average_Jkim.py made by Jungmin Kim, SDlab

import glob
import os
import sys
#import numpy as np

#Making RDF configuration manually

config=[]
filelist={}
time_tag_dict={}

for file in glob.glob("rdf_results/*s_*.rdf"):
	if -1 != file.find("as"):
		time_unit=1e-18
		posi_time=file.find("as")
	elif -1 != file.find("fs"):
		time_unit=1e-15
		posi_time=file.find("fs")
	elif -1 != file.find("ps"):
		time_unit=1e-12
		posi_time=file.find("ps")
	elif -1 != file.find("ns"):
		time_unit=1e-9
		posi_time=file.find("ns")
	elif -1 != file.find("us"):
		time_unit=1e-6
		posi_time=file.find("us")
	elif -1 != file.find("ms"):
		time_unit=1e-3
		posi_time=file.find("ms")
	elif -1 != file.find("s"):
		time_unit=1
		posi_time=file.find("s")
	else:
		continue
	time_tag=file[12:posi_time]
	time=float(time_tag)*time_unit
	time_tag_dict[time]=file[12:posi_time+2]
	filelist[time]=file

sfilekey=sorted(filelist.keys())

weight_list=[]
for ii in range(0,len(sfilekey)):
	if ii==0:
		weight_list.append(format(sfilekey[ii],'.4g'))
	else:
		weight_list.append(format(sfilekey[ii]-sfilekey[ii-1],'.4g'))

time_tag_list=[]
for ii in range(0,len(sfilekey)):
	time_tag_list.append(time_tag_dict[sfilekey[ii]])

for ii in range(0,len(time_tag_list)):
	config.append([time_tag_list[ii],float(weight_list[ii])])

###################################################################################
filelist=[]
for file in glob.glob("rdf_results/*s_*.rdf"):
	filelist.append(file)

##len check
numoffile=len(filelist)
if not numoffile==len(config):
	print('''Error!! filelist and config don't match up''')
	sys.quit()


##Filename match up!!
config_counter=range(0,numoffile) ## config
for temp in filelist:
	for ii in config_counter:
		if config[ii][0] in temp:
			config[ii].append(temp)
##File data read
for ii in config_counter:
	with open(config[ii][2]) as data:
    		lines = data.readlines() 
	temp=[]
	for line in lines:
    		temp.append(line.split())
	data=[[0 for kk in range(len(temp[1])-1)]for jj in range(len(temp)-1)]
	for jj in range(1,len(temp)):
		for kk in range(1,len(temp[1])):
			data[jj-1][kk-1]=float(temp[jj][kk])
	config[ii].append(data) # i - 3 for data
	#config[ii].append(np.array(data)) # i - 4 for numpy data

##File pairname and q
with open(config[ii][2]) as data:
    	lines = data.readlines()
temp=[]
for line in lines:
    	temp.append(line.split())
q=[]
pair=[]

for ii in range(2,len(temp[0])):
	pair.append(temp[0][ii])

for ii in range(1,len(temp)):
	q.append(float(temp[ii][0]))

#print(q) Debug
#print(pairname) Debug

## Weight counter
weight=0
for ii in config_counter:
	weight=weight+config[ii][1]
#print(weight) Debug
## Weight check

#if not weight==weight_original:
	#print('''Error!! weight do not match up !!! ''')
	#exit()

## Average
numofq=len(config[0][3])
numofpair=len(config[0][3][0])
data_aver=[[0 for kk in range(numofpair)]for jj in range(numofq)]
data_var=[[0 for kk in range(numofpair)]for jj in range(numofq)]
data_std=[[0 for kk in range(numofpair)]for jj in range(numofq)]

for ii in config_counter:
	for jj in range(0,numofq):
		for kk in range(0,numofpair):
			data_aver[jj][kk]=data_aver[jj][kk]+config[ii][1]*config[ii][3][jj][kk]/weight

for ii in config_counter:
	for jj in range(0,numofq):
		for kk in range(0,numofpair):
			data_var[jj][kk]=data_var[jj][kk]+config[ii][1]*(config[ii][3][jj][kk]-data_aver[jj][kk])**2.0*numoffile/(numoffile-1)/weight

#print(data_var)

for jj in range(0,numofq):
		for kk in range(0,numofpair):
			data_std[jj][kk]=data_var[jj][kk]**(1.0/2.0)

#print(data_std)

#print(config)

##Make filename
posi=config[0][2].find('s_')
results=config[0][2][posi+2:]
results_filename='rdf_results/Average_%s' %(results)
results_stdname='rdf_results/Average_std_%s' %(results)
##File write
f=open(results_filename,'w')
f.write('#%7s' %('r'))
for ii in range(0,numofpair):
	f.write('%9s' %(pair[ii]))
for ii in range(0,numofq):
	f.write('\n')
	f.write('%8.3f' %(q[ii]))
	for jj in range(0,numofpair):
		f.write('%9.3f' %(data_aver[ii][jj]))
f.close()
##File write std
f=open(results_stdname,'w')
f.write('#%7s' %('r'))
for ii in range(0,numofpair):
	f.write('%9s' %(pair[ii]))
for ii in range(0,numofq):
	f.write('\n')
	f.write('%8.3f' %(q[ii]))
	for jj in range(0,numofpair):
		f.write('%9.3g' %(data_std[ii][jj]))
f.close()
##Print
print('='*60)
print('Made by Jkim. ver 180703')
for ii in config_counter:
	print("Averged rdf files : '%s' with weight : %.4g" %(config[ii][2],config[ii][1]))
print("Done. Result file is '%s'"%(results_filename))
print("Done. Result file is '%s'"%(results_stdname))
print('='*60)
