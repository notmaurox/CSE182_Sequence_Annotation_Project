commentedFile = open("CSE182ProtienAnnotationSheetNEW.tsv","r")
chengzeAnnotations = open("chengzeAnnotations.tsv","r")
#chengzeAnnotationsFinal = open("chengzeAnnotationsFINAL.tsv","w")
#ripAnnotations = open("irisAnnotations.tsv","r")
harambaesAnnotations = open("harambaesAnnotations.tsv","r")
harambaesAnnotationsFinal = open("harambaesAnnotationsFinal.tsv","w")
ripAnnotationsFinal = open("ripAnnotations.tsv","w")
linesRead = 0
harambeHeader ="Sequence\tBLAST\tPfam\tProsite\tALL_OTHER_ANNOTATIONS\tComments\n"
ripHeader = "ID\tBlast\tPfam\tProsite\tKEGG\tGO\tGo-Slim\tComments\n"
harambaesAnnotationsFinal.write(harambeHeader)
ripAnnotationsFinal.write(ripHeader)
for line in commentedFile:
    if linesRead > 0:
        lineC = chengzeAnnotations.readline()
        lineCSplit = lineC.split("\t")
        #lineI = irisAnnotations.readline()
        lineH = harambaesAnnotations.readline()
        splitLine = line.split("\t")
        harambeLine = splitLine[0]+"\t"+splitLine[2]+"\t"+splitLine[4]+"\t"
        harambeLine = harambeLine+splitLine[5]+"\t"+lineCSplit[8]+"\t"+splitLine[6]
        harambeLine = harambeLine.replace("No Hits","null")
        ripLine = splitLine[0]+"\t"+splitLine[2]+"\t"+splitLine[4]+"\t"+splitLine[5]+"\t"+lineCSplit[7]+"\t"+lineCSplit[8]+"\t"+""+"\t"+splitLine[6]
        ripLine = ripLine.replace("No Hits","")
        ripLine = ripLine.replace("No IDs","")
        harambaesAnnotationsFinal.write(harambeLine)
        ripAnnotationsFinal.write(ripLine)
    linesRead+=1
