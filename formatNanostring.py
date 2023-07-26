#!/usr/bin/python

'''
The MIT License (MIT)

Copyright (c) 2016 Charles Lin

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
'''

#131126_formatNanostring
#131126
#Charles Lin


#Description:

#This is a generic python template that has functions from utils.py imported and can be used on CFCE1



#================================================================================
#=============================DEPENDENCIES=======================================
#================================================================================

import sys
import os
import argparse
from collections import defaultdict
import pandas as pd

#print "Using python version %s" % sys.version






#================================================================================
#============================GLOBAL PARAMETERS===================================
#================================================================================

#add locations of files and global parameters in this section





#================================================================================
#===================================CLASSES======================================
#================================================================================

#user defined classes here

#================================================================================
#=================================FUNCTIONS======================================
#================================================================================


def formatFolder(folderName, create=False):

    '''
    makes sure a folder exists and if not makes it
    returns a bool for folder
    '''
    
    if folderName[-1] != '/':
        folderName +='/'

    try: 
        os.listdir(folderName)
        return folderName
    except OSError:
        print('folder %s does not exist' % (folderName))
        if create:
            os.mkdir(folderName)
            return folderName
    return False


def unParseTable(table, output, sep):
    '''
    turns a python list table into a delimited text table
    '''

    fh_out = open(output, 'w')
    if len(sep) == 0:
        for i in table:
            fh_out.write(str(i) + '\n')
    else:
        for line in table:
            fh_out.write(sep.join([str(x) for x in line]) + '\n')

    fh_out.close()




def parseRCC(rccFile):

    '''
    returns a dict keyed by refseq
    '''

    rccList = []

    rcc = open(rccFile,'r')

    for line in rcc:
        # a probe line!
        if ['Endogenous','Positive','Negative','Housekeeping'].count(line.split(',')[0]):
            print(line)
            line = line.rstrip().split(',')
            rccList.append(float(line[3]))

    rcc.close()


    print(rccList)

    return rccList

def parseRCCPlex(rccFile): 
    columnNum = str(int(rccFile.split('/')[-1].split('_')[-1].split('.')[0]))
    sampleNumToLetter = {'1':'A', '2':'B','3':'C','4':'D','5':'E','6':'F','7':'G','8':'H'}
    rccDict =defaultdict(dict)
    refseqToName = {}
    rcc = open(rccFile, 'r')
    for line in rcc: 
        if 'Endogenous' in line.split(',')[0]: 
            line = line.rstrip().split(',')
            sampleRow = sampleNumToLetter[line[0][-2]]
            rccDict["{}_{}".format(columnNum, sampleRow)][line[2]] = float(line[3])
            refseqToName[line[2]] = line[1]
        elif 'Housekeeping'in line.split(',')[0]: 
            line = line.rstrip().split(',')
            sampleRow = sampleNumToLetter[line[0][-2]]
            rccDict["{}_{}".format(columnNum, sampleRow)][line[2]] = float(line[3])
            refseqToName[line[2]] = line[1]
        elif ['Positive','Negative'].count(line.split(',')[0]): 
            line = line.rstrip().split(',')
            sampleRow=sampleNumToLetter[line[1].split('_')[-1]]
            rccDict["{}_{}".format(columnNum, sampleRow)][line[2]] = float(line[3])
            refseqToName[line[2]] = line[1].split('_')[0]
            
    rcc.close()

    rccDF = pd.DataFrame.from_dict(rccDict)
    rccDF['NAME'] = rccDF.apply(lambda row: refseqToName[row.name] , axis=1)

    return rccDF



def labelTable(nanoTable,rccFile):

    '''
    returns a dict keyed by refseq
    '''

    rccDict = {}

    rcc = open(rccFile,'r')

    for line in rcc:
        # a probe line!
        if ['Endogenous','Positive','Negative','Housekeeping'].count(line.split(',')[0]):
            
            line = line.rstrip().split(',')
            refID = line[2].split('.')[0]
            name = line[1]
            nanoTable.append([line[0],name,refID])



    return nanoTable



    


def formatNanostring(folderList,sampleList=[],output =''):

    '''
    formats nanostring into a reasonable table
    '''

    nanoTable =[]
    header = ['CODE_CLASS','NAME','ACCESSION']

    #first fill out the header and create a fileList

    fileList = []
    namesList = []

    folderList = [formatFolder(folderName,False) for folderName in folderList]

    nanoDict = {}

    
    for folderName in folderList:

        folderFileList = os.listdir(folderName)
        folderFileList = [x for x in folderFileList if x[0] != '.' and x.count('RCC') == 1]
        
        fileList += [folderName + x for x in folderFileList]
        header += folderFileList
        namesList += folderFileList
        
    #print fileList
    #parseRCC(fileList[0])
    #print namesList

    nanoTable = [header]

    #fill out the probeID and name columns

    nanoTable = labelTable(nanoTable,fileList[0])

    #print nanoTable

    for i in range(len(fileList)):

        filePath = fileList[i]
        fileName = namesList[i]
        #sampleName = sampleList[i]
        rccList = parseRCC(filePath)
        for j in range(len(rccList)):
            nanoTable[(j+1)].append(rccList[j])


    unParseTable(nanoTable,output,'\t')

    #now run the R code

    
def formatNanostringPlex(folderList,sampleList=[],output =''):
    nanoTable =[]
    header = ['CODE_CLASS','NAME','ACCESSION']

    #first fill out the header and create a fileList

    fileList = []
    namesList = []

    folderList = [formatFolder(folderName,False) for folderName in folderList]

    nanoDict = {}

    
    for folderName in folderList:

        folderFileList = os.listdir(folderName)
        folderFileList = [x for x in folderFileList if x[0] != '.' and x.count('RCC') == 1]
        
        fileList += [folderName + x for x in folderFileList]
        header += folderFileList
        namesList += folderFileList
        
    rccDFList = [parseRCCPlex(file) for file in fileList]
    plexDataDF = pd.concat(rccDFList, axis=1)
    plexDataDF = plexDataDF.loc[:,~plexDataDF.columns.duplicated()]

    return plexDataDF


#================================================================================
#===============================MAIN RUN=========================================
#================================================================================


#write the actual script here


def main():

    '''
    this is the main run function for the script
    all of the work should occur here, but no functions should be defined here
    '''





    #usage = "usage: %prog [options] -i [INPUT_FOLDER] -o [OUTPUT_FOLDER]"
    parser = argparse.ArgumentParser(usage='%(prog)s  -i [INPUT_FOLDER] -o [OUTPUT_FOLDER]')

    # required flags
    parser.add_argument("-i", "--input", dest="input", type=str,
                      help="Enter a folder containing RCC files to be processed", required=True)
    parser.add_argument("-o", "--output", dest="output", type=str,
                      help="Enter the output folder.", required=True)
    parser.add_argument('--plex', dest='plex', default=False, action='store_true')


    # optional
    parser.add_argument("-n", "--name", dest="name", type=str,
                      help="optional name for output table", required=False)


    args = parser.parse_args()

    if args.input and args.output:

        inputFolder = formatFolder(args.input,False)

        #first try to get a naming convention to set up the output        
        outputFolder = formatFolder(args.output,True)

        if args.name:
            outputPath = outputFolder + args.name + '.txt'
        else:
            #get the name from the inputFolder
            name = inputFolder.split('/')[-2]
            outputPath = outputFolder + name + '.txt'

    folderList = [inputFolder]
    if args.plex: 
        print("Using Plex-Set Formatting Rules")
        plexDF = formatNanostringPlex(folderList,[],outputPath)
        plexDF.to_csv("{}".format(outputPath), sep='\t')
    else:    
        formatNanostring(folderList,[],outputPath)

if __name__ == "__main__":
    main()
