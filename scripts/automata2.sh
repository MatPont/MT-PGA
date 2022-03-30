if [ $# -ne 1 ]; then
  echo "Usage "$0 "numberOfThreads"
  exit
fi

for isPD in 0 1; do
  for noThread in $1 1; do
  
    if [ $isPD -eq 0 ]; then
      if [ $noThread -eq 1 ]; then
        outFile="nohupTGGPP_PD_Seq_BSL20.out"
      else
        outFile="nohupTGGPP_PD_Para_BSL20.out"
      fi
    else
      if [ $noThread -eq 1 ]; then
        outFile="nohupTGGPP_Seq_BSL20.out"
      else
        outFile="nohupTGGPP_Para_BSL20.out"
      fi
    fi
  
    ./automata.sh $isPD $noThread > $outFile
  done  
done

