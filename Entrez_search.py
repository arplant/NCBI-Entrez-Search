# -*- coding: utf-8 -*-
"""
arplant

Many thanks to leadbot and the following sources:
https://biopython.org/docs/dev/api/Bio.Entrez.html
http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec145
https://github.com/gumption/Using_Biopython_Entrez/blob/master/Biopython_Tutorial_and_Cookbook_Chapter_9.ipynb

"""

#Overview
#Read query terms from file, retrieve fasta, write the fasta to an output, then process it

#from Bio import SeqIO
#from io import BytesIO
#import os
from Bio import Entrez
from Bio import SeqIO
import time
Entrez.email = 'your_email@gmail.com'
searchName = "LisH_CTLH" # Use appropriate name; this search collects LisH, CTLH and WD40 proteins

#Read input terms from file "input.txt" into a list
inFile = open("input.txt","r")
inputList = inFile.read().splitlines() #Split at linebreaks, does not include newline char in item
print(inputList)
inFile.close()

#Create a set to store unique IDs
#This is used if different search terms are being used to retrieve the same type of protein
combinedIdList = set()


#Loop through the query terms, searching the database and returning IDs to the combined ID list
for query in inputList:
    queryterm = query + "[protein] NOT partial[Description]"
    handle = Entrez.esearch(db="protein", term=queryterm, retmax = 10, retmode="fasta")
    record = Entrez.read(handle)
    handle.close()
    for id in record["IdList"]:
        combinedIdList.add(id)
        #print("Added " + id)


#The combined ID list can be used to retrieve sequence data, etc.
print(combinedIdList)


for entry in combinedIdList:
    handle = Entrez.efetch(db="protein", id=entry, rettype="fasta", retmode="text")


#Write the fasta for each ID to a file
proFasta = open(searchName + ".fa", "a")
for entry in combinedIdList:
    handle = Entrez.efetch(db="protein", id=entry, rettype="fasta", retmode="text")
    proFasta.write(handle.read())
    handle.close()
    time.sleep(0.3) # Manually control query rate; NCBI does not permit high query rates
proFasta.close()

sequences = []
for seq_entry in SeqIO.parse(searchName + ".fa", "fasta"):
    sequences.append(str(seq_entry.seq))
    print(seq_entry.id)





"""
#Nested loop to find IDs based on search terms, then write the associated FASTA sequences to separate files
#For each input term, search the database and recover a list of IDs, then write matching data to file
for query in inputList:
    queryterm = query + "[protein] AND Bacteria[Orgn] NOT partial[Description]"
    handle = Entrez.esearch(db="protein", term=queryterm, retmax = 2, retmode="fasta")
    record = Entrez.read(handle)
    handle.close()
    print(record["IdList"])
    
    proFasta = open(query + ".fa", "a")
    #Working script to write fasta results to file iteratively
    for protId in record["IdList"]:
        handle = Entrez.efetch(db="protein", id=protId, rettype="fasta", retmode="text")
        proFasta.write(handle.read())
        handle.close()
        combinedIdList.add(protId)
        time.sleep(1) #Manually control query rate
    proFasta.close()


print(combinedIdList)

"""




























"""
#esearch for genomes within nuccore that contain gene of interest, remove shotgun
queryterm='BcsZ[gene] AND Bacteria[Orgn] NOT partial[Description]'
handle=Entrez.esearch(db="nuccore", term=queryterm, retmax = 10, retmode="xml")
record=Entrez.read(handle)
handle.close()
print(record["Count"])
print(record.keys())
"""




"""
#Working example to return genbank format data from nucleotide database using ID
handle = Entrez.efetch(db="nucleotide", id="2469782787", rettype="gb", retmode="text")
print(handle.read())
handle.close()
"""




"""
idVal = 2461233154;

handle = Entrez.efetch(db="nuccore", id=idVal, rettype="gb", retmode="text")
record = SeqIO.read(handle, "genbank")
handle.close()
print(record.id)
"""





"""
# Test loop to iterate through dictionary and write IDs to file
with open('bcsZ_nuc_ID.txt', 'w') as f:
    for i in record.get(u'IdList'):
        #print (i)
        f.write(i + '\n')
        #time.sleep(0.2)
f.close()


nuccore_id = "420775483"
idLink = Entrez.elink(dbfrom="nuccore", id=nuccore_id, linkname="nuccore")
linkRecord = Entrez.read(idLink)
print(linkRecord)
idLink.close()


#print(linkRecord[0]["LinkSetDb"][0]["LinkName"])
#linked = [link["Id"] for link in record[0]["LinkSetDb"][0]["Link"]]
#"420775483"in linked

#lr = linkRecord[0]
#print("Link record: " + lr)
"""








"""
#Return single record summary based on Id
testId = Entrez.esummary(db="nuccore", id="NC_000913.3", retmode="xml")
testRecord = Entrez.parse(testId)
for record in testRecord:
    print(record['Title'])
"""



"""
handle = Entrez.efetch(db = "nucleotide", id = "EU490707", rettype = "fasta")
record = SeqIO.read( handle, "fasta" )
record
"""
