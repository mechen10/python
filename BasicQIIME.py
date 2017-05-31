#!/bin/python

import os
import copy

def loadOTUTable(OTUFP): # Loads and converts to relative abundance
	# Make OTU Table
	os.system('biom convert -i ' + OTUFP + ' --to-tsv --header-key taxonomy -o OTU_Table_text.txt')
	# Load OTU Table
	biomOpen = open("OTU_Table_text.txt", 'r') # This is in Unix format
	biomTemp = []
	for i in biomOpen:
		tempLine = i.strip()
		tempLineSplit = tempLine.split("\t")
		biomTemp.append(tempLineSplit)
	biomOpen.close()
	# Get rid of last column
	OTUTable = {} # Now make dictionary
	taxaIDs = {} # Make taxonomy reference 
	for lineN in range(len(biomTemp)):
		if lineN == 0: # This is first line; skip this
			pass
		elif lineN == 1: # This is the header line
			headers = biomTemp[lineN][1:len(biomTemp[lineN])-1]
		else:
			OTUTable[str(biomTemp[lineN][0])] = {}
			taxaIDs[str(biomTemp[lineN][0])] = biomTemp[lineN][len(biomTemp[lineN])-1]
			for abund in range(len(biomTemp[lineN][1:])-1):
				OTUTable[str(biomTemp[lineN][0])][headers[abund]] = biomTemp[lineN][1:][abund]
	# Get total counts for each site
	totalCounts = {}
	for h in range(len(headers)):
		totalCounts[headers[h]] = 0
	for OTU in OTUTable:
		for h in range(len(headers)):
			totalCounts[headers[h]] += float(OTUTable[OTU][str(headers[h])])
	# Convert to relative abundance
	for OTU in OTUTable.keys():
		tempOTUlist = OTUTable[OTU].copy()
		for sites in OTUTable[OTU].keys():
			tempOTUlist[sites] = float(OTUTable[OTU][sites])/float(totalCounts[sites])
		OTUTable[OTU] = tempOTUlist
	os.system("rm OTU_Table_text.txt")
	return OTUTable,taxaIDs # Output is a 2-layer dictionary; first is OTU IDs and second is samples. Also, one-layer dict with taxa IDs
	
def loadMetadata(metadataFP):
	metadataOpen = open(metadataFP, 'U') # U is for 'Universal read'-- automatically turns into Unix LF
	metadataTemp = []
	for i in metadataOpen:
		lineTemp = i.strip()
		lineTempSplit = lineTemp.split("\t")
		metadataTemp.append(lineTempSplit)
	metadataOpen.close()
	metadata = {}
	metadataSites = []
	for lineN in range(len(metadataTemp)):
		if lineN == 0:
			headerList = metadataTemp[lineN]
			for headerName in metadataTemp[lineN]:
				metadata[headerName] = {}
		else:
			for i in range(1,len(metadataTemp[lineN])):
				metadataSites.append(metadataTemp[lineN][0])
				sortHeader = headerList[i]
				metadata[sortHeader][metadataTemp[lineN][0]] = metadataTemp[lineN][i]
	return metadata # output is 2-layer dictionary: first is Metadata and second is samples
	
def printTableFromDictionary(dictionary, output):
	toPrint = ''
	first = True
	for row in dictionary:
		if first == True:
			for column in dictionary[row]:
				toPrint += "\t" + str(column)
				first = False
			toPrint += "\n"
		toPrint += str(row)
		for column in dictionary[row]:
			toPrint += "\t" + str(dictionary[row][column])
		toPrint += "\n"
	open(str(output)+".txt", 'w').write(toPrint)
	print "DONE"
	
	
def getOTUSubset(OTUTable, dictSamples, dictOTUs):
	OTUTableSubsets = {}
	for treatment in allTreatmentList:
		colnames = allTreatmentList[treatment]
		rownames = ALLCORES[treatment]
		newOTUTable = {}
		for OTU in OTUTable:
			if OTU in rownames:
				newOTUTable[OTU] = {}
				for sample in OTUTable[OTU]:
					if sample in colnames:
						newOTUTable[OTU][str(sample)] = OTUTable[OTU][str(sample)]
		OTUTableSubsets[treatment] = newOTUTable
	return(OTUTableSubsets)
	
def removeMinOTUs(originalOTUTable, minThreshold):
	OTUTableInter = copy.deepcopy(originalOTUTable)
	toDelete = []
	for OTU in OTUTableInter:
		allValues = []
		for sample in OTUTableInter[OTU]:
			if OTUTableInter[OTU][sample] < minThreshold:
				OTUTableInter[OTU][sample] = 0
			allValues.append(OTUTableInter[OTU][sample])
		if all(item == 0 for item in allValues):
			toDelete.append(OTU)
	for i in toDelete:
		del OTUTableInter[i]
	return OTUTableInter
	
