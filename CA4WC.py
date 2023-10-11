#!/usr/bin/env python
import numpy as np
import sys
from math import sqrt


fileName = sys.argv[1]
cutoff = 4.147
delta = 5.5
maxClusterSize = 10
numberOfOutputs = 10
clusterMembersWithTime = []

distAvgWithTime   = np.empty(shape=[0,maxClusterSize])
distAvgWithTimeII = np.empty(shape=[0,maxClusterSize])


for i in range(2,maxClusterSize+1):
	clusterIMembersWithTime = []
	clusterMembersWithTime.append(clusterIMembersWithTime)


def distance(x1,y1,z1,x2,y2,z2):
	dist = sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z2-z1)*(z2-z1))
	return dist
			
with open(fileName, 'r') as eqFile :
	dump = eqFile.readlines()
	
tsFound = False #Interested in all the frames
numAtoms = 0
numFrames = 0
xmin = xmax = ymin = ymax = zmin = zmax = 0.0
for i in range(len(dump)):
	line = dump[i]
	if "ITEM: TIMESTEP" in line:
		ts = int(dump[i+1])
		i = i + 1
		if (1) : # ts >= ts0 and ts <= tsf): # Remove the first N-steps as the equilibrium stage
			tsFound = True
			print("Working on the frame ", ts)
		else:
			tsFound = False
			print("Skipping the frame ", ts)
	elif "ITEM: NUMBER OF ATOMS" in line:
		numAtoms = int(dump[i+1])
		i = i + 1
	elif "ITEM: BOX BOUNDS" in line:
		line = dump[i+1].split()
		i = i + 1
		xmin = float(line[0])
		xmax = float(line[1])
		line = dump[i+1].split()
		i = i + 1
		ymin = float(line[0])
		ymax = float(line[1])
		line = dump[i+1].split()
		i = i + 1
		zmin = float(line[0])
		zmax = float(line[1])
	elif "ITEM: ATOMS" in line:
		if tsFound == True :
			tsFound = False
			ids = []
			xs = []
			ys = []
			zs = []
			for j in range(numAtoms):
				i = i + 1
				line   = dump[i].split()
				aid    = int(line[1])-1 # The Lammps has 1-based indexing system
				pType  = int(line[4])
				x = float(line[5])
				y = float(line[6])
				z = float(line[7])
				x = xmin + x*(xmax-xmin)
				y = ymin + y*(ymax-ymin)
				z = zmin + z*(zmax-zmin)
				if (pType == 2 or pType == 3): # Cation or Anion
					ids.append(aid)
					xs.append(x)
					ys.append(y)
					zs.append(z)
					
			clusters = []
			for aid in ids:
				clusters.append(aid)
				
			for i in range(len(xs)):
				for j in range(i+1,len(xs)):
					dist = distance(xs[i],ys[i],zs[i],xs[j],ys[j],zs[j])
					if (dist < cutoff):
						clusters[i] = min([clusters[i],clusters[j]])
						#clusters[i] = (clusters[i]+clusters[j])/2
						clusters[j] = clusters[i]
			
			clusters = np.array(clusters)
			clusterIds, clusterSizes = np.unique(clusters,return_counts=True)
			clusterIds   = np.array(clusterIds)
			clusterSizes = np.array(clusterSizes)
			
			
			distAvgTimeI = []        # in plane distance
			distAvgTimeI.append(ts)
			distAvgTimeII = []       # out of plane distance
			distAvgTimeII.append(ts)
			for i in range(2,maxClusterSize+1):
				index = np.where(clusterSizes == i)[0]
				clusterIdforI = clusterIds[index]
				
				for j in clusterIdforI:
					count = 0
					clusters = np.array(clusters)
					indxs = np.where(clusters == j)[0]
					for k in indxs:
						if (zs[k] >= 46.05-delta and zs[k] <= 60.88+delta):
							count = -1
							break
						elif (zs[k] >= 139.05-delta and zs[k] <= 153.856+delta):
							count = -1
							break
					if (count == 0):
						clusterIdforI = np.setdiff1d(clusterIdforI,[j])
				
				distAvg = 0.0
				distNum = 0
				distAvgII = 0.0
				distNumII = 0
				for j in clusterIdforI:
					indxs =  np.where(clusters == j)[0]
					for k in indxs:
						for l in indxs:
							distnc = distance(xs[k],ys[k],zs[k],xs[l],ys[l],zs[l])
							if (distnc < cutoff):
								dist      = distance(xs[k],ys[k],0.0,xs[l],ys[l],0.0)
								distII    = distance(0.0,0.0,zs[k],0.0,0.0,zs[l]) 
								distAvg   = distAvg   + dist
								distAvgII = distAvgII + distII
								distNum   = distNum   + 1
								distNumII = distNumII + 1
				if (distNum == 0):
					distAvg = 0
				else:
					distAvg = distAvg/distNum
				if (distNumII == 0):
					distAvgII = 0
				else:
					distAvgII = distAvgII/distNumII
					
					
				distAvgTimeI.append(distAvg)
				distAvgTimeII.append(distAvgII)
					
				for clusterId in clusterIdforI:
					clusterMembersWithTime[i-2].append(clusterId)
			
			
			distAvgTimeI      = np.array(distAvgTimeI)
			distAvgTimeII     = np.array(distAvgTimeII) 
			distAvgWithTime   = np.append(distAvgWithTime  ,[distAvgTimeI] ,axis=0)
			distAvgWithTimeII = np.append(distAvgWithTimeII,[distAvgTimeII],axis=0)		
			



clusterAgeMax = np.zeros(shape=[numberOfOutputs, maxClusterSize])

			
for i in range(2,maxClusterSize+1):
	print("Averaging cluster {}".format(i))
	clusterMembersWithTimeI = np.array(clusterMembersWithTime[i-2])
	clusterIds, clusterCounts = np.unique(clusterMembersWithTimeI,return_counts=True)
	clusterCountsInds = clusterCounts.argsort()
	clusterCounts = clusterCounts[clusterCountsInds]
	clusterIds    = clusterIds[clusterCountsInds]
	for j in range(numberOfOutputs):
		if (len(clusterCounts) > j):
			clusterAgeMax[j,i-2] = clusterCounts[-(j+1)] - 1
	
	
with open('Cluster-Age-Summary-delta=5.5-v02.csv', 'w') as filew:
	filew.write("cluster-size")
	for j in range(numberOfOutputs):
		filew.write(",cluster-age-(frames)-Max{}".format(str(j)))
	filew.write("\n")
	for i in range(2,maxClusterSize):
		filew.write(str(i))
		for j in range(numberOfOutputs):
			filew.write(","+str(clusterAgeMax[j,i-2]))
		filew.write("\n")
		
with open('InterparticleDistance-Delta=5.5_v02-InPlanevsOutPlane.csv', 'w') as filew:
	filew.write("time")
	for j in range(numberOfOutputs-1):
		filew.write(",interparticledistance-inplane-clusterSize={}#".format(str(j+2)))
	for j in range(numberOfOutputs-1):
		filew.write(",interparticledistance-outplane-clusterSize={}#".format(str(j+2)))
	filew.write("\n")
	
	for i in range(np.shape(distAvgWithTime)[0]):
		filew.write(str(distAvgWithTime[i,0]))
		for j in range(numberOfOutputs-1):
			filew.write(","+str(distAvgWithTime[i,j+1]))
		for j in range(numberOfOutputs-1):
			filew.write(","+str(distAvgWithTimeII[i,j+1]))
		filew.write("\n")
