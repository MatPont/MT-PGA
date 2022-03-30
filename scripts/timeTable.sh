paths=(sciVis2018/impactEnsemble3CTev_TA/ sciVis2017/cloud5_TA/ sciVis2016/particularEnsemble_TA/ sciVis2015/darkSky100S_TA/ sciVis2014/volcanic2_TA/ sciVis2008/astroTurbulence_TA/ sciVis2008/astro3DTurbulence_TA/ sciVis2006/earthquake2_TA/ dataJules/isabella_velocity_goodEnsemble_TA/ dataJules/startingVortexGoodEnsemble_TA/ dataJules/seaSurfaceHeightGoodEnsemble_TA/ dataJules/vortexStreetGoodEnsemble2_TA/)

names=("Asteroid Impact \cite{scivis2018} (3D)" "Cloud processes \cite{scivis2017} (2D)" "Viscous fingering \cite{scivis2016} (3D)" "Dark matter \cite{scivis2015} (3D)" "Volcanic eruptions \cite{scivis2014} (2D)" "Ionization front \cite{scivis2008} (2D)" "Ionization front \cite{scivis2008} (3D)" "Earthquake \cite{scivis2006} (3D)" "Isabel \cite{scivis2004} (3D)" "Starting Vortex \cite{favelier2018} (2D)" "Sea Surface Height \cite{favelier2018} (2D)" "Vortex Street \cite{favelier2018} (2D)")

names=("Asteroid Impact (3D)" "Cloud processes (2D)" "Viscous fingering (3D)" "Dark matter (3D)" "Volcanic eruptions (2D)" "Ionization front (2D)" "Ionization front (3D)" "Earthquake (3D)" "Isabel (3D)" "Starting Vortex (2D)" "Sea Surface Height (2D)" "Vortex Street (2D)")

defaultPt=0.25
pts=($defaultPt 1 1 $defaultPt $defaultPt 10 $defaultPt $defaultPt 10 $defaultPt 1 $defaultPt)
coefs=(0.5 0.5 0.75 0.5 0 0.5 0 0 0 0 0 0)
defaultEps1=5
eps1s=($defaultEps1 $defaultEps1 $defaultEps1 $defaultEps1 $defaultEps1 10 $defaultEps1 $defaultEps1 $defaultEps1 $defaultEps1 $defaultEps1 $defaultEps1)

mtParaFile="nohupTGGPP_Para_BSL20.out"
mtSeqFile="nohupTGGPP_Seq_BSL20.out"
pdParaFile="nohupTGGPP_PD_Para_BSL20.out"
pdSeqFile="nohupTGGPP_PD_Seq_BSL20.out"

###############################################################################
noPairs=()
noData=()

len=`expr ${#paths[@]} - 1`
for i in `seq 0 $len`; do
  d=${paths[$i]}
  #noPair=`cat $mtParaFile | grep "_______\|trees" | grep $d -A 1 | tail -1 | cut -d\/ -f2 | cut -d, -f1`
  noPair=`cat $mtParaFile | grep "_______\|trees" | grep $d -A 1 | tail -1 | cut -d: -f2 | cut -d\/ -f1`
  noPairs[$i]=$noPair
  noDataI=`cat $mtParaFile | grep "_______\|trees" | grep $d -A 1 | tail -1 | cut -d] -f2 | cut -dt -f1`
  noData[$i]=$noDataI
  
  noPairs[$i]=
  noData[$i]=
done

###############################################################################
get_time_res(){
  res=$(./nohupTimeGetter.sh $1 $2 $3 $4 $5)
  echo $res | cut -d{ -f 2 | cut -d} -f 1
}

###############################################################################
# Create table
###############################################################################

echo -e "\n\n\n"

echo """\begin{table} 
  \caption{Computation times for computing 2 geodesics for all ensembles? both
for MT-PGA and PD-PGA?}
  \centering
  \scalebox{0.63}{
  \hspace{-0.15cm}
  %\makebox[\linewidth]{
    \begin{tabular}{|l|r|r||r|r|r||r|r|r|}
      \hline
      \rule{0pt}{2.25ex} \textbf{Dataset} & \$N$ & $|\branchtree|$ & 
\multicolumn{3}{c||}{PD-PGA} & \multicolumn{3}{c|}{MT-PGA}\\\\
                       &      &                 & 1 c. & 20 c. & Speedup & 1 c. & 20 c. & Speedup\\\\\

      \hline"""

len=`expr ${#paths[@]} - 1`
for i in `seq 0 $len`; do
  d=${paths[$i]}
  d=`echo $d | rev | cut -d\/ -f 1-2 | rev`

  resMtPara=$(get_time_res $mtParaFile $d $2 $3 $4)
  resMtSeq=$(get_time_res $mtSeqFile $d $2 $3 $4)
  mtSpeedup=$(LC_NUMERIC="en_US.UTF-8" printf %.2f $(echo "$resMtSeq / $resMtPara" | bc -l) )

  resPdPara=$(get_time_res $pdParaFile $d $2 $3 $4)  
  resPdSeq=$(get_time_res $pdSeqFile $d $2 $3 $4)
  pdSpeedup=$(LC_NUMERIC="en_US.UTF-8" printf %.2f $(echo "$resPdSeq / $resPdPara" | bc -l) )
  
  #data=$(echo $d | rev | cut -d\/ -f 2 | rev)
  
  echo -n "      "
  if [ $i -eq 0 ]; then
    echo -n "\rule{0pt}{2.25ex} "
  fi
  echo "${names[$i]} & ${noData[$i]} & ${noPairs[$i]} & $resPdSeq & $resPdPara & $pdSpeedup & $resMtSeq & $resMtPara & $mtSpeedup \\\\"
done

echo """      \hline
    \end{tabular}
  }
  \label{tab_timings}
  \vspace{-0.5cm}
\end{table}"""
