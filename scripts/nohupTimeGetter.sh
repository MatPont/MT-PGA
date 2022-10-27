get_min_max_avg(){
    awk 'BEGIN {max=-999999;min=999999} {max=($1>max)?$1:max;min=($1<min)?$1:min;total += $1; count++} END { count=(count==0)?1:count; printf "\\textbf{%.2f} [%.2f ; %.2f] (%d)\n", total/count,min,max,count }' $1
}

if [ $# -lt 2 ]; then
    echo "Usage file dataset [e1 e2 e3]"
    exit
fi

e1=5
e2=95
e3=90
if [ $# -gt 2 ]; then
    e1=$3
    e2=$4
    e3=$5
fi

data=$(echo "$2" | sed -r 's/\//\\\//g') # replace / by \/

if [[ $data == @(*cloud*) ]]; then
  e1=10
fi

allText=$(sed -E -e '/'$data'.*'$e1'.'$e2'.'$e3'.* _/,/Total time/!d' $1)

echo "Total time" $(echo "$allText" | grep "Total time" | cut -d[ -f 6 | cut -ds -f 1 | get_min_max_avg)

echo "Assignment" $(echo "$allText" | grep "Assignment" | cut -d[ -f 6 | cut -ds -f 1 | get_min_max_avg)

echo "Projection" $(echo "$allText" | grep "Projection" | cut -d[ -f 6 | cut -ds -f 1 | get_min_max_avg)

echo "$allText" | grep "Cumulative explained T-Variance" | tail -1


###############################################################################
# Energy
###############################################################################

exit

dir=nohupEnergy
if [ ! -d $dir ]; then
  mkdir $dir
fi

name=$(echo $data | rev | cut -d/ -f 2 | rev | cut -d\\ -f 1)
file=${dir}/${name}_energy.txt

res="Energy, Prop. cost, Ortho. cost, Map. cost, Total\n"

noIteration=`echo "$allText" | grep "Iteration\|Total time" | grep time --max-count 1 -B 99999 | head -n -1 | wc -l`

for i in `seq 1 $noIteration`; do  
  energy=`echo "$allText" | grep "Energy\|Total time" | grep time --max-count 1 -B 99999 | head -n -1 | head -$i | tail -1 | cut -d= -f2`
  res=${res}${energy},
  
  propCost=`echo "$allText" | grep "Prop. cost\|Total time" | grep time --max-count 1 -B 99999 | head -n -1 | head -$i | tail -1 | cut -d= -f2`
  res=${res}${propCost},
  
  orthoCost=`echo "$allText" | grep "Ortho. cost\|Total time" | grep time --max-count 1 -B 99999 | head -n -1 | head -$i | tail -1 | cut -d= -f2`
  res=${res}${orthoCost},
  
  mapCost=`echo "$allText" | grep "Map. cost\|Total time" | grep time --max-count 1 -B 99999 | head -n -1 | head -$i | tail -1 | cut -d= -f2`
  res=${res}${mapCost},
  
  total=`echo "$allText" | grep "Total\|Total time" | grep time --max-count 1 -B 99999 | head -n -1 | head -$i | tail -1 | cut -d= -f2`
  res=${res}${total}"\n"
done

echo -e $res > $file

###############################################################################
# Variance
###############################################################################

file=${dir}/${name}_variances.txt

variances=`echo "$allText" | grep "Explained T-Variance" | cut -d: -f2 | cut -d% -f1`
variances=`echo -e "$variances" | head -2` # TODO get the number of geodesics
echo -e "Variances" > $file
echo -e "$variances" >> $file

###############################################################################

#echo "$allText" | grep "Total" | grep time --max-count 1 -B 99999 | cut -d= -f2 | head -n -1 > $file

#python3 energyPlot.py $file
