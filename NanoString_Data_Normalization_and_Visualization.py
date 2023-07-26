from numpy.core.fromnumeric import mean
from plotly.subplots import make_subplots
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go

from collections import defaultdict
#sets the project directory to instruct program where files are located 
projectDir='c:/Users/brackr2/Documents/Brown_Lab/Experiments/TNF-a_IFN-g_Synergy_Project/NANO-String/RB83'
outputDir='c:/Users/brackr2/Documents/Brown_Lab/Experiments/TNF-a_IFN-g_Synergy_Project/NANO-String/RB83'
#obtains the plexDataFile using the project directory 
plexDataFile = "{}/20220526_6067-JB-PS3869-cart3_RCC.txt".format(projectDir)
#reads in the plexdatafile as a tab seperated file and uses the first column as the index 
plexDataDF = pd.read_csv(plexDataFile, sep='\t', index_col=0)

#NEED TO INCLUDE WAY TO DROP PLEX SETS IF NOT INCLUDED

#sets plexDataDF to all the values besides those in row H
plexDataDF = plexDataDF.loc[:, [sample for sample in plexDataDF.columns if '_C' or'_F' not in sample]]
drop_lst=[]
for i in plexDataDF:
    if 'C' in i:    
        del plexDataDF[i]
    elif 'F' in i:
        del plexDataDF[i]
    #elif '12' in i:
        #del plexDataDF[i]
    #elif '11' in i:
        #del plexDataDF[i]
#del plexDataDF['8_F']
del plexDataDF['12_G']
del plexDataDF['12_H']      
print(plexDataDF)


#drops the title of the index and renames it NAME for plexDataDF
plexDataDF.set_index('NAME', drop=True, inplace=True)
#reads the following file using columns 0 -12 and not including any rows past 
plexPlateDF = pd.read_csv("{}/RB83_Plate_Layout.csv".format(projectDir), header=None,
                         usecols=list(range(0,12)), skiprows=lambda x: x>7)
#Sets the index of PlexPlateDF to A-H
plexPlateDF.index=['A', 'B', 'D','E','G','H']
#names the columns 1-12
plexPlateDF.columns = [column+1 for column in plexPlateDF.columns]

###########################################################################################
#####code to do column wise normalization via geo-mean on pos controls 

#Reflects the DF diagonally so that columns in the original DF become rows and vice virsa
plexDataDF_T = plexDataDF.T
print(plexDataDF_T)

#creates the column plex_group by going through each row and taking the last value in the index(name) and
#assiging that value to the column plex group
plexDataDF_T['plex_group'] = plexDataDF_T.apply(lambda row: row.name[:-2], axis=1)
plexDataDF_T.to_csv('{}/plexDataDF_T.csv'.format(outputDir))
#defines the function posNormFactor which returns the value to which all values in the same column will be multipled by 
#to normalize all values in that column
def posNormFactor(group): 
    #resets multipleFactor to one so that each column is reset
    multipleFactor = 1
    #goes through each value in the Pos colum 
    print(group.POS)
    for control in group.POS: 
        #sets multipleFactor equal to control * multipleFactor
        multipleFactor = control*multipleFactor
    #n is set to number of samples in the same column
    n = len(group.POS)
    #We use geometic mean since this number will be used to normalize all samples in a given column
    geoMean = (multipleFactor)**(1/n)
    #geoMean is returned so that it can be used for normilization downstream
    return geoMean
#creates the variable normalizedGroups which is an empty array for right now    
normalizedGroups = []
lst_geoMean=[]
#for the index and plex group in plexDataDF which is now orderd so that plex_groups are in together and in numerical order
for name,group in plexDataDF_T.groupby('plex_group'): 
    #geoMean is set to the value of the geoMean calculated by the posNormFactor function for each column
    geoMean = posNormFactor(group)
    #creates list with all the geoMeans
    lst_geoMean.append(geoMean)
    #Calculates mean and std for export later in a .txt file 
    mean_1st_norm=np.mean(lst_geoMean)
    stdDev_1st_norm=np.std(lst_geoMean)
    #drops the column plex_group and uses lamda functions to multiply each value in the column by the geoMean of 
    #pos values in a given column
    colNormGroup = group.drop(['plex_group'], axis=1).apply(lambda row: row*(geoMean/row.POS), axis=1)
    #sets the column plex_group to the grouped plex group
    colNormGroup['plex_group'] = group.plex_group
    #appends colNormGroup to the df normalized Groups 
    normalizedGroups.append(colNormGroup)

#adds the values unique to normalized Groups to normPlexDataDF_T which should just be normalized pos values. 
#pos values in the same plex group i.e. all in 1A-1? should be the same 
normPlexDataDF_T = pd.concat(normalizedGroups)
normPlexDataDF_T.to_csv('{}/PosNormPlexDataDF_T.csv'.format(outputDir))
print(normalizedGroups)
print(normPlexDataDF_T)



###################################################################################################
####rowwise normalization by calibration sample (column 1A)

#Column 1 is our calibration sample so rowNormGroup is set to the values from the normPlexDataDF_T in plex group 1
#i.e. all rows with 1 in the name 1A-?
rowNormGroup = normPlexDataDF_T[normPlexDataDF_T.plex_group=='1']
#plex_group, POS, and NEG are dropped from rowNormGroup so that only genes remain in rowNormGroup_genes
rowNormGroup_genes = rowNormGroup.drop(['plex_group', 'POS', 'NEG'],axis=1)

#defines the function geneNormFactor which will calculate the normlization factor to be used for each respective gene 
#accross all plex sets
def geneNormFactor(sampleGeneCounts):
        x=1
        #resets the multpleFactor to 1
        multipleFactor = 1
        #Goes through each calibration sample value for each gene
        for sampleCounts in sampleGeneCounts: 
            #multipleFactor = the value of mulitpleFactor times each respective value for a given gene in the each of
            #the calibration sample
            multipleFactor = multipleFactor*sampleCounts
        #n is set to the length of sampleGeneCounts
        n = len(sampleGeneCounts)
        #We use geometic mean since this number will be multipled by all values for that respective gene in all other 
        #conditions
        geoMean = (multipleFactor)**(1/n)
        #returns geoMean for use in downstream analysis
        return geoMean
#defines the function rowNormFactors which will calculate the normalization factor to be used across each respective plex
#set

def rowNormFactors(group): 
    #obtains geneGeoMeans by calling the geneNormFactor function using he column of the plex group
    geneGeoMeans = group.apply(lambda column: geneNormFactor(column))
    #Calculates the rowScaleFactors by dividing 1 by the value of the gene in the calibration sample by the geomean of 
    #that gene giving the scale factor for each gene in a given row
    rowScaleFactors = 1/group*geneGeoMeans
    #returns the rowScaleFactor so that it can be utilized in downstream analysis
    return rowScaleFactors
#rowScaleFactors is set to the rowScaleFactors returned for the calibtation sample gene values
rowScaleFactors = rowNormFactors(rowNormGroup_genes)
#Calculates the mean and std for the row Scale factors will be used later to create .txt file 
row_norm_mean= rowScaleFactors.mean()
row_norm_std_dev=rowScaleFactors.std()
#scale_means is created as a new column in the df rowScaleFactors and the mean of the column is set to the scale_means cols
rowScaleFactors['scale_means'] = rowScaleFactors.mean(axis=1)
#scale_variance is created as a new column in the df rowScaleFactors and the variance of the col is set to the 
#scale_Variance col
rowScaleFactors['scale_variance'] = rowScaleFactors.var(axis=1)
#zip creates a tuple with the index of the rowscale factor as the first value and the means of the the scale factor as the
#second value. dict is then used to create this into a dictonrary with the index as the key and the mean as the value
rowScaleDict = dict(zip([sample[-1] for sample in rowScaleFactors.index.values.tolist()],rowScaleFactors.scale_means.values.tolist()))
#creates a new column row by pulling the last value off the index which is the row (A-?)
normPlexDataDF_T['row'] = [sample[-1] for sample in normPlexDataDF_T.index.values.tolist()]
#creates the empty list rowNormGroups
rowNormGroups = []
#for the index and plex group in normPlexDataDF_T grouped by the row (A-?) the following will be done
for name,group in normPlexDataDF_T.groupby('row'): 
    #the columns plex_group and row are absent from normGroup and the rest of the values are scaled by the proper
    #rowScaleDict key corresponding to the proper row i.e. row A gets multipled by the mean of row A scale factors
    normGroup = group.drop(['plex_group', 'row'], axis=1).multiply(rowScaleDict[name])
    print(normGroup)
    #the column plex_group is added back to the normGroup df 
    normGroup['plex_group'] = group.plex_group
    print(normGroup)
    #the column row is added back to the normGroup df 
    normGroup['row'] = group.row
    print(normGroup)
    #rowNormGroups assumes the value of normGroup
    rowNormGroups.append(normGroup)
#These changes are updated in the normPlexDataDF_T df
normPlexDataDF_T = pd.concat(rowNormGroups)
print(normPlexDataDF_T)

#Obtain the cond given the plate location starting with row and then column using normPlexDataDF_T

#defines the function get_row_condition which uses the the plexPlateDF layout to match the condition with its subsequent 
#values

#Need to go through and get geomean of each position in the row and then use this to get arthematic mean
#Divide geomean by arthematic to get scaling factor. Apply this across all values in that postion 


def get_row_cond(row,plexPlateDF):
    #print(plexPlateDF.loc['A',int(11)])
    #location is set to the name of the row i.e. 1A
    location=row.name
    #plate row is set to the last character of location i.e A
    plate_row=location[-1]
    #plate_col is set to the characters from index 0 to -1 the end ''.join is used to combine the characters if more than
    #one is presnt
    plate_col=''.join(location[0:-2])
    #checks to see if 'NaN' are present for anyreason 
    if plexPlateDF.loc[plate_row,int(plate_col)]=='nan':
        #if they are it is passed
        pass
    else:
        #otherwise the condition is set to the value of the found in the plexPlateDF at that coordinate
        cond=plexPlateDF.loc[plate_row,int(plate_col)]
        #the condtion is returned so that it can be used in downstream analysis
        return(cond)

with open ('{}/Norm_Output.txt'.format(projectDir),'w') as output:
    output.write('Column Norm Mean \t Column Norm Std Dev \n')
    output.write('{}\t{}\n'.format(mean_1st_norm,stdDev_1st_norm))
    output.write('Row Scale Factor Mean \n')
    output.write('{}\n'.format(row_norm_mean))
    output.write('Row Scale Factor Std Dev \n')
    output.write('{}\n'.format(row_norm_std_dev))
    output.close()
  
#the column condition is created by calling the get_row_cond on every value in a given row
normPlexDataDF_T['condition']=normPlexDataDF_T.apply(lambda row:get_row_cond(row, plexPlateDF), axis=1)
#normPlexDataGB is equal to normPlexDataDF_T grouped by condition
print(normPlexDataDF_T)
normPlexDataDF_T.to_csv('{}/normPlexDataDF_T.csv'.format(outputDir))


#Need to normalize my data by house keeper counts. First need to grab all counts in a given row. 
#Find geomean of all house keepers at every position then find average of all geo means 
#This is scaling factor that should be applied to all counts in that row. house_normalize does so
def house_normalize(row):
    #WANT TO REPLACE WITH COMMAND LINE INSTRUCTIONS TO GRAB HOUSEKEEPERS
    house_lst=['NOL7', 'GUSB','HPRT1']
    #INDEX ORDER IS REALY WEIRD THIS NEEDS TO BE SET AS I CAN NOT GET IT TO GO 1-12
    index_order=['1','10','11','12','2','3','4','5','6','7','8','9']
    forbidden_index=['12_G','12_H']
    #sets up variables for future use
    geo_lst=[]
    geo_dict={}
    ticker=0
    counter=0
    temp_counter=0
    #for plex in normPlexDataDF_T['plex_group']:    
    for plex in row['row']:
        print(len(row['row']))
        if len(row['row'])<12:
            index_order=['1','10','11','2','3','4','5','6','7','8','9']
        #geomean is set to 1 for later use 
        geomean=1
        #this lets us check if this is the first instance or not 
        if '{}_{}'.format(index_order[ticker],plex) in forbidden_index:
            print('{}_{}'.format(index_order[ticker],plex)) 
            pass         
        else:
            if counter==0:
            #temp_plex is utilized to ensure that no rows are misstankely utilized for analysis
                temp_plex=plex
            if temp_plex==plex:
                #for each gene in house_list the values will be mulitplied together
                for gene in house_lst:                               
                    print("{} = {}".format(gene,str(row[gene].loc['{}_{}'.format(index_order[ticker],plex)])))
                    geomean=geomean*float(row[gene].loc['{}_{}'.format(index_order[ticker],plex)])   
                #geomean is generated by raising geomean to 1/ length of house gene              
                geomean=geomean**(1/(len(house_lst)))
                print(geomean)
                #the subsequent geomean is then appended into geolist
                geo_lst.append(geomean)
                #geo_dict is created to pair index positions to the geomean at that position
                geo_dict[(str(index_order[ticker])+'_'+plex)]=geomean
                #values are reset / incrememnted to allow for next calculations 
                geomean=1
                counter+=1
        ticker+=1
    #artimetic mean is generated 
    print(mean(geo_lst))
    arithmetic_mean=mean(geo_lst)
    #returns arithemetic_mean and geo_dict
    print(geo_lst)
    return (arithmetic_mean, geo_dict)
plexset=['A', 'B', 'D','E','G','H']
pos_track=12
score=1 
HouseRowNormGroups=[]
houseNormGroupDF=pd.DataFrame
index_order=['1','10','11','12','2','3','4','5','6','7','8','9']
forbidden_index=['12_G','12_H']
#groups data to be passed by row letter
for name,group in normPlexDataDF_T.groupby('row'):
    #calls house normalize function to be performed for the row. Output data is in the format 
    #of arthmatic mean and dict of index paired to geomean of housekeepers for that given position
    print(group)
    group.to_csv('{}/group.csv'.format(outputDir))
    HouseNorm = house_normalize(group)
    #iterates through using counter to the lenght of the dictonary produced by house_normalize
    #this allows for the tailoring of how many wells are actually present in the sample 
    for counter in index_order:
        #this gets the data at a particular index and drops non float values 
        #allowing for the mulitplication by the index specific geo mean divided by the arthmatic mean
        #to produce scaling factor this is applied across the counts for each gene
        print("Group = {}".format(group))
        print("Postion= {}".format(str(counter)+'_'+name))
        if '{}_{}'.format(str(counter),name) in forbidden_index:
            pass
        else:       
            scaleFactor=(HouseNorm[0]/HouseNorm[1][str(counter)+'_'+name])
            print("Scale Factor = {}".format(scaleFactor))
            houseNormGroup=group.loc[str(counter)+'_'+name].drop(['plex_group','row','condition']).multiply(scaleFactor)
            print(houseNormGroup)
            #the resultant product is a series and is transformed to a data frame by using .to_frame()
            houseNormGroup=houseNormGroup.to_frame()
            #the axis are the inverse of what we need so the df is transformed diagonlly using .T
            houseNormGroup=houseNormGroup.T
            #plex_group is then appended to this df 
            houseNormGroup['plex_group'] = group.plex_group
            #the column row is added back to the normGroup df 
            houseNormGroup['row'] = group.row
            #rowNormGroups assumes the value of normGroup
            houseNormGroup['condition']=group.condition
            #assigns value of condition back to df
            HouseRowNormGroups.append(houseNormGroup)
#the df with all adjust values is then used to replace the old values in the master df
print(HouseRowNormGroups)
normPlexDataDF_T.groupby('row')
print(normPlexDataDF_T)
normPlexDataDF_T = pd.concat(HouseRowNormGroups)
normPlexDataDF_T.to_csv('{}/HouseNorm.csv'.format(outputDir))
normPlexDataGB=normPlexDataDF_T.groupby('condition')
#creates an empty default dictonary which will be filled in later
condition_metrics=defaultdict(dict)
#goes through each condition in  the normPlexDataDF_T df 
print(normPlexDataDF_T['condition'])
normPlexDataDF_T['condition'].to_csv('{}/condition.csv'.format(outputDir))
for cond in normPlexDataDF_T['condition']:
    print(cond)
    #creates a key within the specific cond called mean and is set to the value of the mean of values of each gene in a
    #given conditon   
    condition_metrics[cond]['Mean']=normPlexDataGB.get_group(cond).mean()
    #creates a key within the specific cond called StdDev and is set to the std dev of the values of each gene in a 
    #given cond
    condition_metrics[cond]['StdDev']=normPlexDataGB.get_group(cond).std()

###REPLACE THIS WITH A LIST THAT GETS APPENEDED FOR EVERYTHING IN NAMES
gene_list=['POS','NEG','SELE','CX3CL1','CXCL11','SOCS3','CXCL10','CXCL9','IL32','GUSB','NOL7','HPRT1']
ticker=0
ordered_mean_dict=defaultdict(dict)
ordered_std_dev_dict=defaultdict(dict)
for cond in condition_metrics:
    ticker=0
    for gene in gene_list:
        ordered_mean_dict[gene][cond]=condition_metrics[cond]['Mean'][ticker]
        ordered_std_dev_dict[gene][cond]=condition_metrics[cond]['StdDev'][ticker]
        ticker+=1

#IFN= index 1 TNF = index 0

def order_by_IFN(input_dict):
    ticker=0
    ticker2=0
    cond_list=[]
    counter=0
    for gene in input_dict:
        if ticker==0:
            for condition in input_dict[gene]:
                if condition=='NT':
                    pass
                elif condition=='SYN':
                    pass
                else:
                    conc=condition.split(';')
                    counter=0
                    if len(cond_list)==0:
                        cond_list.append(condition)
                        conc1=cond_list[0].split(';')              
                    elif float(conc[1])<=float(conc1[1]):
                        cond_list.insert(0,condition) 
                        conc1=cond_list[0].split(';')                      
                    else:
                        while counter<=len(cond_list):
                            conc2=cond_list[counter].split(';')
                            if float(conc[1])<=float(conc2[1]):
                                cond_list.insert(counter,condition)  
                                break
                            counter+=1                   
                            if counter==len(cond_list):
                                cond_list.append(condition)
                                break
                            
            ticker+=1
        else:
            break
    return(cond_list)
def sort_by_TNF(input_list):
    master_list=[]
    temp_list=[]
    switch=0
    ticker=0
    counter2=0
    for condition in input_list:
        counter=0        
        conc=condition.split(';')
        if ticker==0:
            temp_conc_TNF=conc[0]
            temp_conc_IFN=conc[1]
            ticker+=1
        if conc[1]!=temp_conc_IFN:
            master_list.append(temp_list)
            temp_list=[]
            temp_conc_TNF=conc[0]
            temp_conc_IFN=conc[1]
        if len(temp_list)==0:
            temp_conc_TNF=conc[0]
            temp_conc_IFN=conc[1]
            temp_list.append(condition)
        elif conc[0]<=temp_conc_TNF:
            temp_list.insert(0,condition)
            temp_conc_TNF=conc[0]
            temp_conc_IFN=conc[1]                  
        else:
            while counter<=len(temp_list):                
                if counter==len(temp_list):
                    temp_list.append(condition)
                    counter+=1
                    break
                conc2=temp_list[counter].split(';')
                if float(conc[0])<=float(conc2[0]):
                    temp_list.insert(counter,condition) 
                    break                   
                counter+=1
        counter2+=1
        if counter2==len(input_list):
            master_list.append(temp_list)
    return(master_list)


def generate_labels(input_list):
    tnf_conc=[]
    ifn_conc=[]
    for group in input_list:
        for cond in group:
            temp_cond=cond.split(';')
            if temp_cond[1] in tnf_conc:
                pass
            else:
                tnf_conc.append(temp_cond[1])
            if temp_cond[0] in ifn_conc:
                pass 
            else:
                ifn_conc.append(temp_cond[0])
    return (ifn_conc, tnf_conc)


def generate_predicted_values(ordered_mean_dict):
    cond_lst=(sort_by_TNF(order_by_IFN(ordered_mean_dict)))
    predicted_dict=defaultdict(dict)
    IFN_Cond=[]
    TNF_Cond=[]
    for gene in ordered_mean_dict:
        count=0
        for cond1 in cond_lst:
            for condition in cond1:
                cond=condition.split(';')
                if cond[0]==cond[1]:
                    pass
                elif len(condition)==3 and count==0:
                    pass
                elif cond[0] =='0' and len(cond[0])==1:
                    IFN_Cond.append(condition)
                elif cond[1]=='0'and len(cond[1])==1:
                    TNF_Cond.append(condition)
                count+=1
        for cond2 in TNF_Cond:
            temp_TNF_cond=cond2.split(';')
            for cond3 in IFN_Cond:
                temp_IFN_cond=cond3.split(';')
                if temp_TNF_cond[0]!=0 and temp_IFN_cond[1]!=0:
                    predicted_dict[gene][str(temp_TNF_cond[0])+';'+str(temp_IFN_cond[1])]=(ordered_mean_dict[gene][cond2]+ordered_mean_dict[gene][cond3])
    return (predicted_dict)

def generate_labels_for_expected_plot():
    cond_lst=(sort_by_TNF(order_by_IFN(ordered_mean_dict)))
    count=0
    TNF_list=[]
    IFN_list=[]
    for gene in ordered_mean_dict:
        if count==0:
            for cond1 in cond_lst:
                for condition in cond1:
                    cond=condition.split(';')
                    if len(cond[0])>1 and len(cond[1])>1:
                        if cond[1] in TNF_list:
                            pass
                        else:
                            TNF_list.append(cond[1])
                        if cond[0] in IFN_list:
                            pass 
                        else:
                            IFN_list.append(cond[0])
        else:
            break
        count+=1
    return (IFN_list,TNF_list)

def generate_data(input_dict, input_list, expected_dict):
    temp_data_lst=[]
    temp_expected_data_lst=[]
    master_list=[]
    master_expected_list=[]
    temp_dict={}
    expected_temp_dict={}
    for gene in input_dict:        
        temp_dict={}
        master_list=[]
        if gene =='POS' or gene =='NEG':
            pass 
        else:
            for cond in input_list:
                for conc in cond:
                    temp_data_lst.append(input_dict[gene][conc])
                    cond_check=conc.split(';') 
                    if '0' in conc:
                        if len(cond_check[0])==1 or len(cond_check[1])==1:
                            temp_expected_data_lst.append(input_dict[gene][conc]) 
                        else:                  
                            temp_expected_data_lst.append(expected_dict[gene][conc])
                    else:
                        temp_expected_data_lst.append(expected_dict[gene][conc])
                master_list.append(temp_data_lst)
                master_expected_list.append(temp_expected_data_lst)
                temp_data_lst=[]     
                temp_expected_data_lst=[]

            expected_temp_dict[gene]=master_expected_list
            master_expected_list=[]
            temp_dict[gene]=master_list

        if len(temp_dict)!=0 and len(expected_temp_dict)!=0:
            pass
    generate_heat_map(temp_dict,expected_temp_dict)   
def generate_data_no_expected(input_dict,input_list):
    temp_data_lst=[]
    master_list=[]
    temp_dict={}
    expected_temp_dict={}
    for gene in input_dict:        
        temp_dict={}
        master_list=[]
        if gene =='POS' or gene =='NEG':
            pass 
        else:
            for cond in input_list:
                for conc in cond:
                    temp_data_lst.append(input_dict[gene][conc])
                master_list.append(temp_data_lst)
                temp_data_lst=[]             
            temp_dict[gene]=master_list
        if len(temp_dict)!=0:
            pass
            generate_heat_map_no_expect(temp_dict)

def generate_heat_map(input_dict,expected_dict):
    labels=generate_labels(sort_by_TNF(order_by_IFN(ordered_mean_dict)))
    expected_labels=generate_labels_for_expected_plot()
    print(expected_labels[0])
    expected_labels[0].insert(0,0)
    expected_labels[1].insert(0,0)
    for gene in expected_dict:
        #data=input_dict[gene]
        expected_data=expected_dict[gene]
        #fig = make_subplots(rows=1, cols=2)
        #fig.add_traces(px.imshow(data,
        #        labels=dict(x="IFNg", y="TNFa", color="Gene Expression"),
        #        x=labels[0],
        #        y=labels[1],
        #        color_continuous_scale='Reds',
        #                     ))
        fig=(px.imshow(expected_data,
                labels=dict(x="IFNg", y="TNFa", color="Gene Expression"),
                x=expected_labels[0],
                y=expected_labels[1],
                color_continuous_scale='RdBu_R',))
        fig.update_xaxes(side="top")
        fig.update_layout(
        title='{}'.format(gene),
        xaxis_nticks=36)
        fig.write_image("{}/HeatMap_{}_Expected.pdf".format(outputDir,gene))
        fig.show() 
def generate_heat_map_no_expect(input_dict):
    labels=generate_labels(sort_by_TNF(order_by_IFN(ordered_mean_dict)))
    print(input_dict)
    for gene in input_dict:
        data=input_dict[gene]
        data = np.array(data, dtype=float)
        fig=px.imshow(data,
                labels=dict(x="IFNg", y="TNFa", color="Gene Expression"),
                x=labels[0],
                y=labels[1],
                color_continuous_scale='RdBu_R',)
        fig.update_xaxes(side="top")
        fig.update_layout(
        title='{}'.format(gene),
        xaxis_nticks=36)
        fig.write_image("{}/HeatMap_{}.pdf".format(outputDir,gene))
        fig.show() 

def generate_bar_plot(ordered_mean_dict,ordered_std_dev_dict):
    gene_counts=[]
    gene_stddev=[]
    condition_list=[]
    gene_df=pd.DataFrame()
    labels=sort_by_TNF(order_by_IFN(ordered_mean_dict))
    for gene in ordered_mean_dict:
        gene_counts=[]
        gene_stddev=[]
        condition_list=[]
        gene_df=gene_df.iloc[0:0]
        print(gene_df)
        if gene== 'POS':
            pass 
        elif gene =='NEG':
            pass
        else:
            for cond in ordered_mean_dict[gene]:
                gene_counts.append(ordered_mean_dict[gene][cond])
                gene_stddev.append(ordered_std_dev_dict[gene][cond])
                condition_list.append(cond)
            print(gene_df)
            gene_df['Counts']=gene_counts
            gene_df['Values']=condition_list
            print(gene_df)
            fig = px.bar(gene_df, x='Values', y='Counts',error_y=gene_stddev,title='{}'.format(gene))
            fig.write_image("{}/{}.png".format(outputDir,gene))
            fig.show()

generate_data(ordered_mean_dict, sort_by_TNF(order_by_IFN(ordered_mean_dict)),generate_predicted_values(ordered_mean_dict))
generate_data_no_expected(ordered_mean_dict, sort_by_TNF(order_by_IFN(ordered_mean_dict)))
generate_bar_plot(ordered_mean_dict,ordered_std_dev_dict)
