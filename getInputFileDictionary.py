###
# Creates a dictionary out of the relevant GENRAY or CQL3D namelist file
# the top level keys are the various sections of the namelist, and variables are then a level below
# Example: genrayDict['plasma']['temp_scale(1)']
###

import numpy as np
#this method works by looking for '&'s, which denote where the sections are
def createInputFileDictionary(path):
    inputFile = open(path,'r')
    
    inputFileDict = {}
    
    inputLines = inputFile.readlines()
    
    variableValue = ''
    variableName = ''

    sectionName = ''

    for i in range(len(inputLines)):
        line = inputLines[i].strip()
        if len(line) == 0:
            continue

        if line[0] == '&':# and line[1:]!='end':
            #if this is not the first section, add the last variable to the dictionary before moving onto this new section
            if sectionName != '':
                addVariable(inputFileDict, sectionName, variableName, variableValue)
                variableValue = ''
                variableName = ''

            sectionName = line[1:]
            inputFileDict[sectionName] = {}
            continue 


        #"""       
        splitLine = line.split('=')
        
        #if this line is a continuation of a previous variable (such as dentab, for example)
        if len(splitLine) == 1:
            splitLine = splitLine[0].strip().replace(' ',',')
            variableValue = variableValue + splitLine + ','
            if i == len(inputLines) - 1:
                addVariable(inputFileDict, sectionName, variableName, variableValue)
        #if this is a new variable
        else:
            if variableName != '':
                addVariable(inputFileDict, sectionName, variableName, variableValue)
                variableValue = ''
                variableName = ''
            
            variableName = splitLine[0].strip(); variableValue = splitLine[1].strip().replace(' ',',') + ','
        #"""
            
    return inputFileDict

#adds a variable to the dictionary under the relevant section 
def addVariable(dictionary, sectionName, variableName, variableValue):
    #get rid of extra comma at the end
    variableValue = variableValue[:-1]
    splitVariable = np.array(variableValue.split(','))
    
    if len(splitVariable) > 1:
        #try to store the variable as a number
        try:
            dictionary[sectionName][variableName]=splitVariable.astype(np.float64)
        #if that didn't work, just store it as a string
        except Exception as exc:
            dictionary[sectionName][variableName] = variableValue
    else:
        try: 
            dictionary[sectionName][variableName]=float(variableValue)
        except:
            dictionary[sectionName][variableName]=variableValue
   
# main function
def getInputFileDictionary(gen_or_cql,targetDir = None):
    if targetDir is None:
        import getTargetInfo
        targetDir = getTargetInfo.getTargetDir()
    print(f'targetDir: {targetDir}')
    if gen_or_cql == 'genray':
        return createInputFileDictionary(f'{targetDir}/genray.in')
    elif gen_or_cql == 'cql3d':
        return createInputFileDictionary(f'{targetDir}/cqlinput')
        
