#!/bin/sh

######
# To run: ./create_graph_horizontal_crack.sh >/dev/null -o 1,2,3,4,5
#####

rm *dat *eps
while test $# -gt 0; do
        case "$1" in
                -h|--help)
                        echo "create graph"
                        echo " "
                        echo "options:"
                        echo "-h, --help"
			echo "-o, --orders=1,2,3,4,5 (max is 5)\"" 
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
	        esac
done

if [ -z ${ORDERS+x} ]; then 
  echo "Approximation order is unset, create graph --help"
  exit 1
fi

sed -n '/average griffith force order 1/,/average griffith force order 2/p' log_griffith_forces | tee order_1.txt
sed -n '/average griffith force order 2/,/average griffith force order 3/p' log_griffith_forces | tee order_2.txt
sed -n '/average griffith force order 3/,/average griffith force order 4/p' log_griffith_forces | tee order_3.txt
sed -n '/average griffith force order 4/,/average griffith force order 5/p' log_griffith_forces | tee order_4.txt
sed -n '/average griffith force order 5/,/average griffith force order 6/p' log_griffith_forces | tee order_5.txt


for order in `echo $ORDERS | sed -e 's/,/ /g'`
do	
	egrep '(average griffith force order|Problem ELASTIC_MECHANICS Nb. rows|average griffith force)' order_"$order".txt | tee order_"$order"_2.txt

	awk '{ printf "%s ", $0 }' order_"$order"_2.txt | tee order_"$order"_3.txt

	sed -e 's/average griffith force order /\'$'\n/g' order_"$order"_3.txt | tee order_"$order"_4.txt

	sed '/^$/d' order_"$order"_4.txt | tee order_"$order"_5.txt

	awk '{print "average griffith force order " $0}' order_"$order"_5.txt | tee input_file_"$order".dat

	rm  order_"$order"*.txt
done

gnuplot plot_file_horizontal_crack.in

open output.eps
