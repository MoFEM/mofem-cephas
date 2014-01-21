#!/bin/sh

# Make files for ploting load-disp. paths

# ./make_plot.sh LOG_FILE_NAME OUTPUT_FILE

LOGFILE=$1
BASEFILENAME=$2

# Node under force
nODEFORCE=141
#NODEFORCE=321

echo "0 0" | tee $BASEFILENAME && awk "/load_path_disp ent $NODEFORCE dim 1/ {print \$13,10*\$15}" $LOGFILE | tee -a $BASEFILENAME


# Nodes on the right/left side for crack mouth opening
NODE1=279
NODE1=280
#NODE1=1394
#NODE2=1395

awk '/load_path_disp ent 1 dim 2/ {print 10*$15,$13}' $LOGFILE | tee tmp1
awk "/load_path_disp ent $NODE1 dim 2/ {print \$13}" $LOGFILE | tee tmp2

awk '/load_path_disp ent 2 dim 2/ {print $13}' $LOGFILE | tee tmp3
awk "/load_path_disp ent $NODE2 dim 2/ {print \$13}" $LOGFILE | tee tmp4

paste -d '\t' tmp1 tmp2 tmp3 tmp4 | tee tmp5

echo "0 0" | tee $BASEFILENAME'_gap' | awk '{print ($2-$3+$4-$5)/2,$1}' tmp5 | tee -a $BASEFILENAME'_gap'



