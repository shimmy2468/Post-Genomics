"""
Original, vcf2csv.py, writeCSV function author: Bofei

2023 Modified version author: Amanda
    Includes Info column from VCF, splitInfo, CaseID, and ParseCSQ functions
"""

import os
import pandas as pd
import vcf


################# Split the info column for the CSQ information ##################

def splitINFO(file): # file would have to be a dataframe
    csq = []
    for row in range(file.shape[0]):
        info1 = str(file["info"][row])
        start = info1.find("CSQ")
        end = info1.find("SOMATIC")
        if end > start:
            info = info1[start+6:end-3]
            info = info.replace("'", "")
            info = info.replace("[", "")
            info = info.replace("]", "")
            info = info.replace(" ", "")
            info = info.split(",")
        else:
            info = info1[start+6:-2]
            info = info.replace("'", "")
            info = info.replace("[", "")
            info = info.replace("]", "")
            info = info.replace(" ", "")
            info = info.split(",")

        csq.append(info)
    file.insert(file.shape[1], "CSQ", csq)
    return file



################### separate CSQ list into individual elements for new columns ######################
    # AML headers
AMLheader = "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|ALLELE_NUM|DISTANCE|STRAND|FLAGS|VARIANT_CLASS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|RefSeq|GENE_PHENO|SIFT|PolyPhen|DOMAINS|HGVS_OFFSET|GMAF|AFR_MAF|AMR_MAF|EAS_MAF|EUR_MAF|SAS_MAF|AA_MAF|EA_MAF|ExAC_MAF|ExAC_Adj_MAF|ExAC_AFR_MAF|ExAC_AMR_MAF|ExAC_EAS_MAF|ExAC_FIN_MAF|ExAC_NFE_MAF|ExAC_OTH_MAF|ExAC_SAS_MAF|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|ENTREZ|EVIDENCE"

def parseCSQ (file, header): #file needs to be a dataframe, header needs to be based on the information within the vcf
    newColheader = header.split("|")

    # creates new column headers based on inputted list where the entries are empty list that can be appended to:
    for col in newColheader:
        file[col] = [list() for x in range(len(file.index))]


    for row in range(file.shape[0]):
        info1 = str(file["CSQ"][row])
        info1 = info1.split(",")

        for ele in info1:
            count = 0
            ele = ele.split("|")

            while count < len(ele):
                file[newColheader[count]][row].append(ele[count])
                count += 1
    return file

def CaseIDs(vcffile):
    temp = vcffile.split(".")
    vcffile = open('%s'%(vcffile), 'r')
    vcf = vcffile.readlines()
    vcffile.close()

    for l in range(0, len(vcf)): #looking line by line in the inputted VCF file
        if vcf[l].startswith('##INDIVIDUAL'):

            ind1 = vcf[l].find("ID=")
            ind2 = vcf[l].find("NAME=")

            line = vcf[l]
            if ind1 == -1:
                print("Problem with: " + temp[0])
                return
            else:
                #print(line[ind+3:-2])
                ele1 = temp[0] #Patient ID/VCF file name
                ele2 = line[ind1+3:-2] #Case ID that will be matched with cBio portal
                ele3 = line[ind2+5:ind1-1]
                #print(row)

    return ele1, ele2, ele3

def writecsv(inputfile, header):
    vcf_reader = vcf.Reader(open(inputfile,'r'))  # read vcf file
    #use nested list to store variant information, one for tumor, one for normal
    normal_var=[]
    tumor_var=[]
    var_score = [28]
    refct = []
    altct = []
    refct1 = []
    altct1 = []

    for record in vcf_reader:
        curr_var = []
        left = []
        right = []
        ref_seq = []
        alt_seq = []
        info = []
        # Because the AML and ALL have a reverse order for which comes first in the VCF, the sample[0/1] needed to be changed in order to make sure
            # the tumor and normal variants were saved to the correct csv. Can be compared to the original PrCa code
        tumorAD = record.samples[1]["AD"]
        normalAD = record.samples[0]["AD"]
        ########################################### tumor_var #############################################
        #######check for 2 ALT options:
        if len(tumorAD) == 3:
            sum = tumorAD[0] + tumorAD[1] + tumorAD[2]
            x = tumorAD[0]
            y = tumorAD[1]
            z = tumorAD[2]

            """
            INSERT CODE HERE
            if two out of the three are not greater than 0.05 throw away variant

            Hint: insert conditional statement that uses the defined variables x, y, z, and sum
            that passes over variant if the above is not true.
            
            """

        ########only 1 base for the ALT sequence
        else:
            sum = tumorAD[0] + tumorAD[1]
            refcount = tumorAD[0]
            altcount = tumorAD[1]

            # both the REF and ALT have to be greater
            if (refcount / sum) > 0.05 and (altcount / sum) > 0.05:
                curr_var.append(record.CHROM)
                left.append(record.POS)
                right.append(record.POS + len(record.ALT))
                ref_seq.append(record.REF)
                alt_seq.append(str((record.ALT[0])))
                refct.append(refcount)
                altct.append(altcount)
                info.append(str(record.INFO))
                lists = curr_var + left + right + ref_seq + alt_seq + var_score + info
                tumor_var.append(lists)


            # throw away variant if not
            else:
                pass

        ############################# normal_var ##############################3
        curr_var1 = []
        left1 = []
        right1 = []
        ref_seq1 = []
        alt_seq1 = []
        info = []
        if len(normalAD) == 3:
            sum = normalAD[0] + normalAD[1] + normalAD[2]
            x = normalAD[0]
            y = normalAD[1]
            z = normalAD[2]

            """
            INSERT CODE HERE
            if two out of the three are not greater than 0.05 throw away variant

            Hint: insert conditional statement that uses the defined variables x, y, z, and sum
            that passes over variant if the above is not true.

            """



        ########only 1 base for the ALT sequence
        else:
            sum = normalAD[0] + normalAD[1]
            refcount = normalAD[0]
            altcount = normalAD[1]

            # both the REF and ALT have to be greater
            if (refcount / sum) > 0.05 and (altcount / sum) > 0.05:
                curr_var1.append(record.CHROM)
                left1.append(record.POS)
                right1.append(record.POS + len(record.ALT))
                ref_seq1.append(record.REF)
                alt_seq1.append(str((record.ALT[0])))
                info.append(str(record.INFO))
                refct1.append(refcount)
                altct1.append(altcount)

                lists = curr_var1 + left1 + right1 + ref_seq1 + alt_seq1 + var_score + info
                normal_var.append(lists)


            # throw away variant if not
            else:
                pass


    #transfer list into dataframe, it would be easier for following manipulation
    normal=pd.DataFrame(normal_var,columns=['chrom','left','right','ref_seq','var_seq1','var_score', 'info'])
    tumor=pd.DataFrame(tumor_var,columns=['chrom','left','right','ref_seq','var_seq1','var_score', 'info'])

    normal_seq2 = []
    tumor_seq2 = []

    normal.insert(5, "count1", altct1)
    normal.insert(6, "count2", refct1)
    for i in range(normal.shape[0]):

        if normal['ref_seq'][i] > normal['var_seq1'][i]:
            normal_seq2.append(normal['ref_seq'][i])

        else:
            normal_seq2.append(normal['var_seq1'][i])
            normal['var_seq1'][i] = normal['ref_seq'][i]
            # switching count values because the ref order is now var_seq1
            temp = normal['count1'][i]
            normal['count1'][i] = normal['count2'][i]
            normal['count2'][i] = temp

    normal.insert(5, "var_seq2", normal_seq2)  # inserts column with information given

    tumor.insert(5, "count1", altct)
    tumor.insert(6, "count2", refct)
    for j in range(tumor.shape[0]):
        if tumor['ref_seq'][j] > tumor['var_seq1'][j]:
            tumor_seq2.append(tumor['ref_seq'][j])

        else:
            tumor_seq2.append(tumor['var_seq1'][j])
            tumor['var_seq1'][j] = tumor['ref_seq'][j]
            # switching count values because the ref order is now var_seq1
            temp = tumor['count1'][j]
            tumor['count1'][j] = tumor['count2'][j]
            tumor['count2'][j] = temp

    tumor.insert(5, "var_seq2", tumor_seq2)

    #inserts and sets the var_index using a list of number within the range of the rows
    tumor.insert(0, "var_index", list(range(tumor.shape[0])))
    normal.insert(0, "var_index", list(range(normal.shape[0])))
    #name = "Patient" + str(counts)
    fil = inputfile.split(".")
    id = fil[0]

    ## inserts VCF ID as a column
    tumor.insert(tumor.shape[1], "VCF_ID", id)
    normal.insert(normal.shape[1], "VCF_ID", id)



    ################################
    # Inserts three columns, vcf file name, patient ID, and the alternate ID
    ele1, ele2, ele3 = CaseIDs(inputfile)

    tumor.insert(tumor.shape[1], "Patient_ID", ele3)
    tumor.insert(tumor.shape[1], "Case_ID", ele2)

    normal.insert(normal.shape[1], "Patient_ID", ele3)
    normal.insert(normal.shape[1], "Case_ID", ele2)

    ################################



    ################################
    # for VCF with info columns, use the previously defined functions to split this column
    # if info column contains no information, comment that code out below
    tumor = splitINFO(tumor)
    tumor = parseCSQ(tumor, header)

    normal = splitINFO(normal)
    normal = parseCSQ(normal, header)

    # Remove the Info and CSQ Column
    tumor = tumor.drop(columns=["info", "CSQ"])
    normal = normal.drop(columns=["info","CSQ"])
    ##################################



    ############### export the CSV file ###############
    tumor.to_csv(id + "_tumor.csv")
    normal.to_csv(id + "_normal.csv")
    return



## Make sure you change the directory to where your VCF files are saved on your device
# os.chdir("/Users/abataycan/AML/AMLFiles")
# files = os.listdir()
#
# for file in files:
#     writecsv(file, AMLheader)






