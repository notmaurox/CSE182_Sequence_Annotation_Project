import sys
import urllib,urllib2
import os
import subprocess
import json
from Bio import SeqIO
from Bio.ExPASy import ScanProsite

#Example annotated data
#E. coli protein ANK03648.1 with the annotation Dihydropteroate synthase
# current_directory = os.getcwd()
# final_directory = os.path.join(current_directory, r'rawProtienAnnotationData')
# if not os.path.exists(final_directory):
#    os.makedirs(final_directory)
os.system("mkdir rawProtienAnnotationData")
os.system("mkdir BLAST")
os.system("mkdir BLAST_match")
os.system("mkdir Pfam")
os.system("mkdir Pfam_match")
os.system("mkdir Prosite")
os.system("mkdir Prosite_match")


indexFile = {}

fasta_sequences = SeqIO.parse(open(sys.argv[1]),'fasta')
startindex = 701
stopindex = 801
seqindex = startindex
annotationFile = open("protienAnnotations.txt","w")
chengzeAnnotations = open("chengzeAnnotations.tsv","w")
chengzeSpeciesFile = open("chengzeSpecies.tsv","w")
harambaesAnnotations = open("harambaesAnnotations.tsv","w")
irisAnnotations = open("irisAnnotations.tsv","w")
header = "Seq ID\tBLAST on UniProt\tPfam\tProsite\tKEGG IDs\tGO IDs\tComments\n"
headerPersonal = "Seq ID\tBLAST_Match\tBLAST_Function\tPfam_Match\tPfam_Function\tProsite_Match\tComments\n"
irisAnnotations.write(header)
annotationFile.write(headerPersonal)
#assuming first seq = 1
fastasRead = 0
while seqindex < stopindex:
    for fasta in fasta_sequences:
            if fastasRead == seqindex-1:
                seqid, seq = fasta.id, str(fasta.seq)
                fastasRead += 1
                break
            fastasRead += 1
    seqid = seqid.replace("embl-cds:","")
    seqData = {}
    indexFile[seqid] = seqData
    seqDataLoc = "rawProtienAnnotationData/" + seqid
    makeSeqDir = "mkdir " + seqDataLoc
    os.system(makeSeqDir)
    queryFileLoc = seqDataLoc + "/"+seqid+"query.seq"
    queryFile = open(queryFileLoc,"w")
    queryFile.write(seq)
    queryFile.close()
    localQueryFile = open("query.seq","w")
    localQueryFile.write(seq)
    seqData["query"] = seq
    localQueryFile.close()
    #Search against UniProt using BLAST
    db = "uniprot_sprot.fasta"
    blastOutLoc = "BLAST" + "/"+seqid+"_raw.txt"
    cmd = "blastp -query "+queryFileLoc+" -db " + db + " -out "+blastOutLoc+" -outfmt 7"
    os.system(cmd)
    #Get best hit ID
    blastOut = open(blastOutLoc,"r")
    linesRead = 0
    matchID = ""
    blastOutString = ""
    for line in blastOut:
        blastOutString = blastOutString + line
        if linesRead >= 5:
            if line[8:11]=="sp|" and line[17] == "|":
                matchID = line[11:17]
                break
        linesRead+=1
    seqData["BLAST"] = blastOutString
    #Search uniprot for best hit ID
    if matchID != "":
        link = "http://www.uniprot.org/uniprot/"+ matchID +".txt"
        matchPage = urllib.urlopen(link)
        matchInfo = matchPage.read()
        uniprotMatchLoc = "BLAST_match" + "/" + seqid + "_match.txt"
        uniprotMatch = open(uniprotMatchLoc,"w")
        uniprotMatch.write(matchInfo)
        seqData["BLASTmatch"] = matchInfo
        uniprotMatch.close()
        #Pull data,
        uniprotMatchRead = open(uniprotMatchLoc,"r")
        uniprotResult = ""
        uniprotKeywords = ""
        uniprotOS = ""
        DElinesRead = 0
        uniprotFirstKeyword = ""
        keggIDs = ""
        goIDs = ""
        functionStart = False
        functionEnd = False
        uniprotFunction = ""
        for line in uniprotMatchRead:
            #Descriptions
            if line[0:2] == "DE":
                ##Pull alt names(older classification) and ref names(latest),
                if DElinesRead == 0:
                    currLine  = line[5:len(line)].lstrip().replace('RecName: Full=','').replace('AltName: Full=','')
                    currLineSplit = currLine.split("{")
                    uniprotFirstKeyword = currLineSplit[0]
                uniprotKeywords = uniprotKeywords + line[5:len(line)].lstrip().replace('RecName: Full=','RecName: ').replace('AltName: Full=','AltName:')
                uniprotResult = uniprotResult + line[5:len(line)].lstrip().replace('RecName: Full=','RecName: ').replace('AltName: Full=','AltName:')
                DElinesRead += 1
            #OrganismClassification
            if line[0:2] == "OC":
                uniprotResult = uniprotResult + line[5:len(line)].lstrip()
                uniprotOS = uniprotOS + line[5:len(line)].lstrip()
            #KEGG ID's
            if line[0:11] == "DR   KEGG; ":
                keggSPLIT = line.split(";")
                keggIDs = keggIDs + keggSPLIT[1] + ";"
            if line[0:9] == "DR   GO; ":
                goSPLIT = line.split(";")
                goIDs = goIDs + goSPLIT[1] + ";"
            if "-!-" in line and functionStart == True:
                functionEnd = True
            if "-!- FUNCTION:" in line:
                functionStart = True
                uniprotFunction = line.replace("CC   -!- FUNCTION:","")
            if functionStart == True and functionEnd == False and "CC       " in line:
                uniprotFunction = uniprotFunction + line.replace("CC       "," ")
    #print("BLAST Fx")
    #print(uniprotFunction)
    keggIDs = keggIDs.replace(" ","")
    goIDs = goIDs.replace(" ","")
    if uniprotResult == "":
        uniprotResult = "No Hits"
    #Search against Pfam using HMMER
    cmd2 = "curl -L -H 'Expect:' -H 'Accept:text/xml' -F hmmdb=pfam -F seq='<query.seq' http://www.ebi.ac.uk/Tools/hmmer/search/hmmscan"
    p = subprocess.Popen(cmd2, stdout=subprocess.PIPE, shell=True)
    out, err = p.communicate()
    result = out.split('>')
    output = ""
    for lin in result:
        if not lin.startswith('#'):
            output = output + " " + lin
    outPfamLoc =  "Pfam" + "/" + seqid + "_raw.txt"
    outPfam = open(outPfamLoc,"w")
    outPfam.write(output)
    seqData["Pfam"] = output
    outPfam.close()
    pfamResult = ""
    bestPfamID = ""
    bestPfamDesc = ""
    decRead = 0
    idRead = 0
    for i in range(len(output)):
        if output[i:i+5] == " acc=":
            info = ""
            i2 = i+5+1
            while output[i2] != output[i+5]:
                if idRead == 0:
                    bestPfamID = bestPfamID + output[i2]
                i2 += 1
            idRead += 1
        if output[i:i+6] == " desc=":
            info = ""
            i2 = i+6+1
            while output[i2] != output[i+6]:
                pfamResult = pfamResult + output[i2]
                if decRead == 0:
                    bestPfamDesc = bestPfamDesc + output[i2]
                i2 += 1
            pfamResult = pfamResult + "; "
            decRead += 1
    if pfamResult == "":
        pfamResult = "No Hits"
        pfamMatchData = "No Hits"
    #Get Pfam webpage
    if pfamResult != "No Hits":
        pfamMatchPageURL = "http://pfam.xfam.org/family/" + bestPfamID
        pfamMatchPage = urllib.urlopen(pfamMatchPageURL)
        pfamMatchInfo = pfamMatchPage.read()
        pfamMatchLoc = "Pfam_match" + "/" + seqid + "_match.txt"
        pfamMatchFile = open(pfamMatchLoc,"w")
        pfamMatchFile.write(pfamMatchInfo)
        pfamMatchFile.close()
        pfamMatchFileRead = open(pfamMatchLoc,"r")
        pfamDataFound = False
        h1FlagFound = False
        pfamMatchData = ""
        for line in pfamMatchFileRead:
            if "pfamData" in line:
                pfamDataFound = True
            if "</h1>" in line and pfamDataFound == True:
                h1FlagFound = True
            if "<p>" in line and pfamDataFound == True and h1FlagFound == True:
                 pfamMatchData = line
                 break
        pfamMatchData = pfamMatchData.replace("          <p>","")
        pfamMatchData = pfamMatchData.replace("</p>","")
        #print(pfamMatchData)
        #print("NEW SEQ")
        if pfamMatchData == "":
            pfamMatchData = "" #Not sure how to fill this space, probably best to leave blank to make it easiest for other teams. 
        #pfamMatchInfo stores the html page for the prosite match with given ID
        #Not sure what to do with it yet
    #Search against prosite
    handle = ScanProsite.scan(seq=seq)
    result = ScanProsite.read(handle)
    type(result)
    #get 'signature_ac': u'XXXXXXX'
    # print("PROSITE DATA")
    # print( result.n_match )
    # print( result[0] )
    resultSig = ""
    prositeResult = "No Hits"
    #Make file containing prosite output

    #+++++++++++++++++++++++++++++++++++
    if result.n_match > 0:
        #Make file containing prosite output
        prositeOutLoc =  "Prosite" + "/" + seqid + "_raw.txt"
        prositeData = open(prositeOutLoc,"w")
        prositeDataString = ""
        for i in range(result.n_match):
            prositeData.write(str(result[i])+"\n")
            prositeDataString = prositeDataString + str(result[i])+"\n"
        resultID = str(result[0])
        seqData["Prosite"] = prositeDataString
        prositeData.close()
        for i in range(len(resultID)):
            if resultID[i:i+18] == "'signature_ac': u'":
                resultSig = resultID[i+18:i+25]
                #print(resultSig)
        prositeLink = "http://prosite.expasy.org/" + resultSig + ".txt"
        prositematchPage = urllib.urlopen(prositeLink)
        prositematchInfo = prositematchPage.read()
        #print prositematchInfo
        prositeLoc =  "Prosite_match" + "/" + seqid + "_match.txt"
        prositeMatch = open(prositeLoc,"w")
        prositeMatch.write(prositematchInfo)
        seqData["Prositematch"] = prositematchInfo
        prositeMatch.close()
        prositeMatchRead = open(prositeLoc,"r")
        prositeResult = ""
        for line in prositeMatchRead:
            if line[0:2] == "DE":
                prositeResult = prositeResult + line[5:len(line)]

    if bestPfamID == "":
        bestPfamID = "No ID"
    if bestPfamDesc == "":
        bestPfamDesc = "No Hits"
    if resultSig == "":
        resultSig = "No ID"
    if keggIDs == "":
        keggIDs = "No IDs"
    if goIDs == "":
        goIDs = "No IDs"
    pfamResult = pfamResult.replace("\n","")
    prositeResult = prositeResult.replace("\n","")
    # #MAKE LIST OF ALL AVALIABE VARIABLES UP TO THIS POINT
    # #ID best BLAST hit on UniProt(SwissProt):
    # # - matchID
    # print("BLAST hit ID")
    # print(matchID)
    # #Keywords of best BLAST hit on UniProt:
    # # - uniprotKeywords
    # print("BLAST hit keywords")
    # print(uniprotKeywords)
    # #First Keyword of best BLAST hit on Uniprot (first DE line)
    # # - uniprotFirstKeyword
    # print("BLAST hit first keyword")
    # print(uniprotFirstKeyword)
    # #Organism classification of best BLAST hit on UniProt:
    # # - uniprotOS
    # print("BLAST hit organism specs")
    # print(uniprotOS)
    # #Combination of uniprotKeywords and uniprotOS
    # # - uniprotResult
    # print("BLAST keywords +  organism specs")
    # print(uniprotResult)
    # #Keywords from all hits to Pfam
    # # - pfamResult
    # print("Pfam keywords: ALL")
    # print(pfamResult)
    # #ID from hit to Pfam with best score, follow http://pfam.xfam.org/family/...
    # # - bestPfamID
    # print("BEST Pfam ID")
    # print(bestPfamID)
    # #Keywords from hit to Pfam with best score
    # # - bestPfamDesc
    # print("Pfam keywords: BEST")
    # print(bestPfamDesc)
    # #ID best Prosite match
    # # - resultSig
    # print("Prosite match ID")
    # print(resultSig)
    # #Keywords of best Prosite hit
    # # - prositeResult
    # print("Prosite keywords")
    # print(prositeResult)
    # #KEGG IDs
    # # - keggIDs
    # print("KEGG IDs")
    # print(keggIDs)
    # #Go IDs
    # # - goIDs
    # print("GO IDs")
    # print(goIDs)
    comment = "COMMENTxPLACEHOLDER"
    finalOut = seqid+"\t"+uniprotResult.replace('\n', ' ')+"\t"+pfamResult.replace('\n', ' ')+"\t"+prositeResult.replace('\n', ' ')
    #annotationFile.write( finalOut + "\n" )
    #Write chengze's file
    uniprotKeywordsClean = uniprotKeywords.replace('\n', ' ').replace('RecName: Full=','').replace('AltName: Full=','').replace('RecName: ','').replace('AltName: ','').replace('AltName:','').replace('RecName:','').replace('Short=','')
    chengzeOut = seqid.replace("embl-cds:","")+"\t"+matchID+"\t"+uniprotKeywordsClean+"\t"+bestPfamID.replace('\n', ' ')+"\t"+bestPfamDesc.replace('\n', ' ')+"\t"+resultSig.replace('\n', ' ')+"\t"+prositeResult+"\t"+keggIDs+"\t"+goIDs+"\t"+comment+"\t"+"w/e"
    chengzeOut.replace('\n',"")
    chengzeOut = chengzeOut.replace(";\t","\t").replace("; \t","\t")
    chengzeAnnotations.write( chengzeOut )
    chengzeSpeciesFile.write(seqid.replace("embl-cds:","")+"\t"+uniprotOS[0:len(uniprotOS)-2].replace("\n",""))
    #Write Haram-Baes file
    harambaesOut = seqid.replace("embl-cds:","")+"\t"+uniprotKeywordsClean+"\t"+pfamResult+"\t"+prositeResult+"\t"+keggIDs+"\t"+goIDs+"\t"+comment
    harambaesOut.replace("\n","")
    harambaesAnnotations.write(harambaesOut)
    #Wrtie Iris file
    irisOut = seqid.replace("embl-cds:","")+"\t"+uniprotFunction.replace("\n"," ")+"\t"+pfamMatchData.replace("\n"," ")+"\t"+prositeResult+"\t"+keggIDs+"\t"+goIDs+"\t"+comment
    irisAnnotations.write(irisOut)
    #Write MY FILE SMH, TREAT YO SELF
    #"Seq ID\tBLAST_Match\tBLAST_Function\tPfam_Match\tPfam_Function\tProsite_Match\tComments"
    myOut = seqid.replace("embl-cds:","")+"\t"+uniprotKeywords.replace('\n',"")+"\t"+uniprotFunction.replace("\n"," ")+"\t"+bestPfamDesc.replace("\n","")+"\t"+pfamMatchData.replace("\n"," ")+"\t"+prositeResult.replace("\n","")+"\t"+comment
    myOut = myOut.replace("\n"," ")
    annotationFile.write( myOut )
    if seqindex != stopindex-1:
        chengzeAnnotations.write("\n")
        chengzeSpeciesFile.write("\n")
        harambaesAnnotations.write("\n")
        irisAnnotations.write("\n")
        annotationFile.write("\n")
    seqindex += 1
annotationFile.close()
chengzeAnnotations.close()
chengzeSpeciesFile.close()
harambaesAnnotations.close()
irisAnnotations.close()

os.system("rm query.seq")

#Create Json file.
j = json.dumps(indexFile, indent=4)
f = open('rawProtienAnnotationDataIndexd.json', 'w')
print >> f, j
f.close()
