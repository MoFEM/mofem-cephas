#!/bin/sh

######
# To run: ./create_graph_three_point_bending_SC.sh >/dev/null -o 1,2,3,4 -r 1,2
#####
rm -r Plot_files_SC
rm *.eps
mkdir Plot_files_SC

while test $# -gt 0; do
        case "$1" in
                -h|--help)
                        echo "create graph"
                        echo " "
                        echo "options:"
                        echo "-h, --help"
			echo "-o, --orders=1,2,3,4,5 (max is 5)\""
			echo "-r, --ref=0,1,2,3,4 ... (memory and time is your limit)\""
                        exit 0
                        ;;
		-o)
                        shift
                        if test $# -gt 0; then
                                export ORDERS=$1
                        else
                                echo "no griffith energy specified"
                                exit 1
                        fi
                        shift
                        ;;
                --orders*)
                        export ORDERS=`echo $1 | sed -e 's/^[^=]*=//g'`
                        shift
                        ;;
		-r)
                        shift
                        if test $# -gt 0; then
                                export REF=$1
                        else
                                echo "no griffith energy specified"
                                exit 1
                        fi
                        shift
                        ;;
                --ref*)
                        export REF=`echo $1 | sed -e 's/^[^=]*=//g'`
                        shift
                        ;;

	        esac
done

if [ -z ${ORDERS+x} ]; then 
  echo "Approximation order is unset, create_graph_three_point_bending --help"
  exit 1
fi
if [ -z ${REF+x} ]; then 
  echo "Refinement levels is unset, create_graph_three_point_bending --help"
  exit 1
fi



for order in `echo $ORDERS | sed -e 's/,/ /g'`
do
	for ref_level in `echo $REF | sed -e 's/,/ /g'`
	do
	
	#copy log file from each arc_length analysis and rename with order and refinement in name
	#cp arc_length_out_spatial_"$order"_"$ref_level"_*/log log_"$order"_"$ref_level"

	echo "0 0" | tee Plot_files_SC/load_path_CrackArea_vs_Lambda_"$order"_"$ref_level".dat && awk '/load_path:/ {print $4,-$6}' log_files_SC/log_"$order"_"$ref_level" | tee -a Plot_files_SC/load_path_CrackArea_vs_Lambda_"$order"_"$ref_level".dat
	#ent 276 is immidiately to the left of the crack, ent 279 is immediately to the right of the crack
	echo "0 0" | tee Plot_files_SC/load_path_DispY_vs_Lambda_"$order"_"$ref_level".dat && awk '/load_path_disp ent 279 dim 1/ {print -$13,-$15}' log_files_SC/log_"$order"_"$ref_level" | tee -a Plot_files_SC/load_path_DispY_vs_Lambda_"$order"_"$ref_level".dat
	done

done

