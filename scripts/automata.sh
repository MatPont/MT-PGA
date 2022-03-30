#!/bin/bash

if [ $# -ne 1 ]; then
  echo "Usage "$0 "numberOfThreads"
  exit
fi

prepath="../data/"

paths=(startingVortexGoodEnsemble_TA/ isabella_velocity_goodEnsemble_TA/ seaSurfaceHeightGoodEnsemble_TA/ vortexStreetGoodEnsemble2_TA/ particularEnsemble_TA/ cloud5_TA/ astroTurbulence_TA/ impactEnsemble3CTev_TA/ darkSky100S_TA/ volcanic2_TA/ astro3DTurbulence_TA/ earthquake2_TA/)
len=`expr ${#paths[@]} - 1`
for i in `seq 0 $len`; do
    path=${paths[$i]}
    paths[$i]=$prepath$path
done

defaultPt=0.25
pts=($defaultPt 1 1 $defaultPt $defaultPt 10 $defaultPt $defaultPt 10 $defaultPt 1 $defaultPt)
coefs=(0.5 0.5 0.75 0.5 0 0.5 0 0 0 0 0 0)
defaultEps1=5
eps1s=($defaultEps1 $defaultEps1 $defaultEps1 $defaultEps1 $defaultEps1 10 $defaultEps1 $defaultEps1 $defaultEps1 $defaultEps1 $defaultEps1 $defaultEps1)

GREEN='\033[0;32m'
NC='\033[0m' # No Color

eps2=95
eps3=90

noRun=1
barySizeLimit=20

isPD=$1
noThread=$2

for j in `seq 1 $noRun`; do
  len=`expr ${#paths[@]} - 1`
  for i in `seq 0 $len`; do
    path=${paths[$i]}
    pt=${pts[$i]}
    coef=${coefs[$i]}
    eps=${eps1s[$i]}
    
    coef=0 # only compute on split trees
    
    echo -e "${GREEN}_____________ gpcaMT.py $path $pt $eps $eps2 $eps3 $coef $barySizeLimit $isPD $noThreads _____________${NC}"
    pvpython gpcaMT.py $path $pt $eps $eps2 $eps3 $coef $barySizeLimit $isPD $noThreads 
  done
done
