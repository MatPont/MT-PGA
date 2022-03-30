#!/usr/bin/env pvpython

from paraview.simple import *
import sys
import os

pt = 0.
eps1 = 5.
eps2 = 95.
eps3 = 90.
coef = 0.0
barySizeLimit = 0.0
nodePerTask = 32
isPD = False
deterministic = 1
noThreads = 1

openCDB = True
noGeod = 2
noPS = 4

mode = 0 # 0 : parallel _ 1 : embed _ 2 : data red.

# Filled
fileName = None
saveCoef = False 
saveOutput = False
saveRecTrees = False
dontSave = False

def getParametrizedName(fileName, pt, eps1, eps2, eps3, coef=None, isPD=False, noGeod=None, barySizeLimit=None):
  import os

  outFileName = os.path.basename(fileName)[:-4]
  outFileName += "_PT_"+str(pt)
  if isPD:
    outFileName += "_PD"
  else:
    outFileName += "_E1_"+str(eps1)
    outFileName += "_E2_"+str(eps2)
    outFileName += "_E3_"+str(eps3)
  outFileName += "" if coef is None else "_C_" + str(coef)
  outFileName += "" if noGeod is None else "_NG_" + str(noGeod)
  outFileName += "" if barySizeLimit is None or barySizeLimit == 0.0 else "_BSL_" + str(barySizeLimit)
  return outFileName

def cdbRead(fileName):
  # create a new 'TTK CinemaReader'
  tTKCinemaReader1 = TTKCinemaReader(DatabasePath=fileName)

  # create a new 'TTK CinemaQuery'
  tTKCinemaQuery1 = TTKCinemaQuery(InputTable=tTKCinemaReader1)

  # create a new 'TTK CinemaProductReader'
  tTKCinemaProductReader1 = TTKCinemaProductReader(Input=tTKCinemaQuery1)
  tTKCinemaProductReader1.AddFieldDataRecursively = 1

  # create a new 'TTK FlattenMultiBlock'
  tTKFlattenMultiBlock1 = TTKFlattenMultiBlock(Input=tTKCinemaProductReader1)

  return tTKFlattenMultiBlock1

def getGpcaMTOutFileName(port = 0):
  return dirGPCAMT + getParametrizedName(fileName, pt, eps1, eps2, eps3, coef=coef, isPD=isPD, barySizeLimit=barySizeLimit) + ("_bary.vtm" if port == 0 else "_coef.vtt" if port == 1 else "_vec.vtt" if port == 2 else "_corr.vtt")
  
def isGpcaMTOutFilesExist():
  exist = True
  for port in range(4):
    exist &= os.path.isfile(getGpcaMTOutFileName(port))
  return exist

def getEmbedOutFileName():
  paramName = getParametrizedName(fileName, pt, eps1, eps2, eps3, coef=coef, isPD=isPD, barySizeLimit=barySizeLimit)
  if deterministic:
    return dirEmbed + paramName + "_coef.csv"
  else:
    i = 0
    while True:
      outFileName = dirEmbed + "nonDeterministic/" + paramName + "_coef_" + str(i) + ".csv"
      if not os.path.isfile(outFileName):
        break
      i += 1
    return outFileName
    
def getDataRedOutFileName(isSecondInput = False):
  return dirDataRed + getParametrizedName(fileName, pt, eps1, eps2, eps3, coef=coef, isPD=isPD, noGeod=noGeod, barySizeLimit=barySizeLimit) + "_trees" + ("" if not isSecondInput else "2") + ".vtm"

def decodeAndSave(tTKMergeTreePrincipalGeodesics, isSecondInput = False):
    # create a new 'TTK MergeTreePrincipalGeodesicsDecoding'
    tTKMergeTreePrincipalGeodesicsDecoding1 = TTKMergeTreePrincipalGeodesicsDecoding(Barycenter=tTKMergeTreePrincipalGeodesics,
        Coefficients=OutputPort(tTKMergeTreePrincipalGeodesics,1),
        GeodesicsVector=OutputPort(tTKMergeTreePrincipalGeodesics,2),
        CorrelationMatrix=OutputPort(tTKMergeTreePrincipalGeodesics,3),
        InputTreesoptionnal=None)
    tTKMergeTreePrincipalGeodesicsDecoding1.OutputBarycenter = 0
    if isSecondInput:
      tTKMergeTreePrincipalGeodesicsDecoding1.ProcessSecondInput = 1
    tTKMergeTreePrincipalGeodesicsDecoding1.UpdatePipeline()
       
    if not dontSave:
      dataRedOutFileName = getDataRedOutFileName(isSecondInput) 
      writer = XMLMultiBlockDataWriter(Input=tTKMergeTreePrincipalGeodesicsDecoding1, FileName=dataRedOutFileName)
      writer.UpdatePipeline()
      print(dataRedOutFileName)

def gpcaMT(group1, group2=None, tTKMergeTreePrincipalGeodesics1=None):
  if tTKMergeTreePrincipalGeodesics1 is None:
    # create a new 'TTK MergeTreePrincipalGeodesics'
    tTKMergeTreePrincipalGeodesics1 = TTKMergeTreePrincipalGeodesics(Input=group1, OptionalInput=group2)
  tTKMergeTreePrincipalGeodesics1.NumberOfGeodesics = noGeod
  tTKMergeTreePrincipalGeodesics1.NumberOfProjectionIntervals = 16
  tTKMergeTreePrincipalGeodesics1.NumberOfProjectionSteps = noPS
  tTKMergeTreePrincipalGeodesics1.Deterministic = deterministic
  tTKMergeTreePrincipalGeodesics1.Epsilon1 = eps1
  tTKMergeTreePrincipalGeodesics1.Epsilon2 = eps2
  tTKMergeTreePrincipalGeodesics1.Epsilon3 = eps3
  tTKMergeTreePrincipalGeodesics1.PersistenceThreshold = pt
  tTKMergeTreePrincipalGeodesics1.NodePerTask = nodePerTask
  tTKMergeTreePrincipalGeodesics1.ParallelMode = 0
  tTKMergeTreePrincipalGeodesics1.NormalizedWasserstein = not isPD
  tTKMergeTreePrincipalGeodesics1.UseAugmentedEnergy = 0
  tTKMergeTreePrincipalGeodesics1.JoinSplitMixtureCoefficient = coef
  tTKMergeTreePrincipalGeodesics1.BarycenterSizeLimitPercent = barySizeLimit
  
  tTKMergeTreePrincipalGeodesics1.UseAllCores = 0
  tTKMergeTreePrincipalGeodesics1.ThreadNumber = noThreads
  
  tTKMergeTreePrincipalGeodesics1.UpdatePipeline()
  
  if saveCoef:
    if not dontSave:
      outFileName = getEmbedOutFileName()
      print(outFileName)
      writer = CSVWriter(Input=OutputPort(tTKMergeTreePrincipalGeodesics1, 1), FileName=outFileName)
      writer.UpdatePipeline()
    
    if saveOutput and not dontSave:
      for port in range(4):
        outFileName = getGpcaMTOutFileName(port)
        print(outFileName)
        if port == 0:
          writer = XMLMultiBlockDataWriter(Input=(tTKMergeTreePrincipalGeodesics1, port), FileName=outFileName)
        else:
          writer = XMLTableWriter(Input=OutputPort(tTKMergeTreePrincipalGeodesics1, port), FileName=outFileName)
        writer.UpdatePipeline()
    
  if saveRecTrees:
    decodeAndSave(tTKMergeTreePrincipalGeodesics1)
    
if __name__ == "__main__":

  #####################
  # Load parameters
  #####################
  dirPath = sys.argv[1]
  pt = float(sys.argv[2])
  eps1 = float(sys.argv[3])
  eps2 = float(sys.argv[4])
  eps3 = float(sys.argv[5])
  coef = float(sys.argv[6]) # 0: ST ; 1: JT
  barySizeLimit = float(sys.argv[7])
  isPD = isPD if len(sys.argv) < 9 else int(sys.argv[8])
  noThreads = noThreads if len(sys.argv) < 10 else int(sys.argv[9])
    
  print("persistence_threshold =", pt)
  print("epsilon               =", eps1)
  print("epsilon2              =", eps2)
  print("epsilon3              =", eps3)
  print("coef                  =", coef)
  print("barySizeLimit         =", barySizeLimit)
  print("noThreads             =", noThreads)
  
  if saveCoef:
    print("Save coef.")
  if saveOutput:
    print("Save output")
  if saveRecTrees:
    print("Save rec. trees")
  
  #####################
  # Load data and execute
  #####################
  # Get file path
  files = os.listdir(dirPath)
  
  stFilePath = dirPath+[a for a in files if a[-12:] == "ST_light.cdb"][0]
  dataST = cdbRead(stFilePath)
  
  fileName = "".join(os.path.basename(stFilePath).split("_ST_light"))

  gpcaMT(dataST)
    
    
    
    
