#!/bin/bash

. ./utils.sh

if [ $# -lt 1 ]; then
    echo "Usage: fileName [e1 e2 e3]"
    exit
fi

echo "=============" $1 "============="
for d in ${paths[@]}; do
    d=`echo $d | rev | cut -d\/ -f 1-2 | rev`
    echo -e "\n=======" $d "======="
    ./nohupTimeGetter.sh $1 $d $2 $3 $4
done



###############################################################################
# Create table
###############################################################################

echo -e "\n\n\n"

echo """\begin{table} 
    \centering
    %\scalebox{0.65}{
    \makebox[\linewidth]{%
    %\begin{tabular}{|l||r|r||r||r|r|r|}
    \begin{tabular}{|l||r|}
        \hline
        %\textbf{Dataset} & \$N$ & $|\branchtree|$ & $\wasserstein{2}$ \julien{(1 c.)} & $\wassersteinTree$ (1 c.) & $\wassersteinTree$ (20 c.) & Speedup \\\\
        \textbf{Dataset} & MT-PGA (20 c.) \\\\
        \hline"""
        
for d in ${paths[@]}; do
    d=`echo $d | rev | cut -d\/ -f 1-2 | rev`
    res=$(./nohupTimeGetter.sh $1 $d $2 $3 $4)
    res=$(echo $res | cut -d{ -f 2 | cut -d} -f 1)
    data=$(echo $d | rev | cut -d\/ -f 2 | rev)
    echo "          $data & $res \\\\"
done

echo """        \hline
    \end{tabular}
    }
    }
    \label{tab_timeSeq}
\end{table}"""
