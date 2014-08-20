#!/bin/sh

#### 
# ./plot_graphs_PMMA.sh -o 1,2,3 -r 1,2,3,4 -a 45/60/75
####

while test $# -gt 0; do
        case "$1" in
                -h|--help)
			echo "you must be using gnuplot 4.4 or later"
                        echo "create graph"
                        echo " "
                        echo "options:"
                        echo "-h, --help"
			echo "-o, --orders=1,2,3,4,5 (max is 5)\""
			echo "-r, --ref=0,1,2,3,4 ... (memory and time is your limit)\""
			echo "-a, --angle=45,60 or 75\""
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
		-a)
			shift
			if test $# -gt 0; then
				export ANGLE=$1
			else
				echo "no angle specified"
				exit 1
			fi
			shift
			;;
	        esac
done

if [ -z ${ORDERS+x} ]; then 
  echo "Approximation order is unset, plot_graph --help"
  exit 1
fi
if [ -z ${REF+x} ]; then 
  echo "Refinement levels is unset, plot_graph --help"
  exit 1
fi
if [ -z ${ANGLE+x} ]; then
  echo "Angle of inclination of initial crack is unset, plot_graph -- help"
  exit 1
fi

#Remove previous plot files directory and create new directory
rm -r plot_files
mkdir plot_files

#Determine maximum and minimum approximation orders to be plotted
MAX_O=$(echo "$ORDERS" | sed -e 's/,/ /g' | sed -e 's/^.* //')
MIN_O=$(echo "$ORDERS" | sed -e 's/,/ /g' |  awk -F" " '{print $1}')

#Determine maximum and minimum refinement levels to be plotted
MAX_R=$(echo "$REF" | sed -e 's/,/ /g' | sed -e 's/^.* //')
MIN_R=$(echo "$REF" | sed -e 's/,/ /g' |  awk -F" " '{print $1}')

echo "Please wait, generating graphs..."

sed -n '/-my_order 1/,/-my_order 2/p' log_files/log_"$ANGLE" | tee log_files/log_1 >/dev/null
sed -n '/-my_order 2/,/-my_order 3/p' log_files/log_"$ANGLE" | tee log_files/log_2 >/dev/null
sed -n '/-my_order 3/,/-my_order 4/p' log_files/log_"$ANGLE" | tee log_files/log_3 >/dev/null

sed -n '/-my_ref 1/,/-my_ref 2/p' log_files/log_1 | tee log_files/log_1_1 >/dev/null
sed -n '/-my_ref 2/,/-my_ref 3/p' log_files/log_1 | tee log_files/log_1_2 >/dev/null
sed -n '/-my_ref 3/,/-my_ref 4/p' log_files/log_1 | tee log_files/log_1_3 >/dev/null
sed -n '/-my_ref 4/,/-my_ref 5/p' log_files/log_1 | tee log_files/log_1_4 >/dev/null
sed -n '/-my_ref 1/,/-my_ref 2/p' log_files/log_2 | tee log_files/log_2_1 >/dev/null
sed -n '/-my_ref 2/,/-my_ref 3/p' log_files/log_2 | tee log_files/log_2_2 >/dev/null
sed -n '/-my_ref 3/,/-my_ref 4/p' log_files/log_2 | tee log_files/log_2_3 >/dev/null
sed -n '/-my_ref 4/,/-my_ref 5/p' log_files/log_2 | tee log_files/log_2_4 >/dev/null
sed -n '/-my_ref 1/,/-my_ref 2/p' log_files/log_3 | tee log_files/log_3_1 >/dev/null
sed -n '/-my_ref 2/,/-my_ref 3/p' log_files/log_3 | tee log_files/log_3_2 >/dev/null
sed -n '/-my_ref 3/,/-my_ref 4/p' log_files/log_3 | tee log_files/log_3_3 >/dev/null
sed -n '/-my_ref 4/,/-my_ref 5/p' log_files/log_3 | tee log_files/log_3_4 >/dev/null

for order in `echo $ORDERS | sed -e 's/,/ /g'`
do
	for ref_level in `echo $REF | sed -e 's/,/ /g'`
	do
	# $7,$8,$9 are xyz coordinates, $11 is g1, $18 is g3, $13 is j.

	#45 degrees
	awk '/griffith force at ent/ {print $7,$8,$9,$11,$18,$13}' log_files/log_"$order"_"$ref_level" | tee -a plot_files/all_coords_and_g_"$order"_"$ref_level".dat >/dev/null
	awk '{if ($1 >= -1e-10) print $1,$2,$3,$4,$5,$6}' plot_files/all_coords_and_g_"$order"_"$ref_level".dat | tee -a plot_files/coords_and_g_"$order"_"$ref_level".dat >/dev/null
	sort -gk 1 plot_files/coords_and_g_"$order"_"$ref_level".dat | tee plot_files/sorted_coords_and_g_"$order"_"$ref_level".dat >/dev/null
	done
done


rm *eps
if [ $ANGLE == 45 ]; then
gnuplot << EOF
		#Set up gnuplot defining axis labels, key, etc
		set terminal postscript eps enhanced "Helvetica" 25 size 10,7
		set output "K2_K1_vs_x3_d.eps"
		set xlabel "x3/d"
		set ylabel "KII/KI"
		#set yrange [0.0:]
		set xrange [:1.2]
		set key above
		set key spacing 2
		set grid back lc rgb "black" linewidth 1
		set grid mxtics mytics ls 101 lc rgb "black" linewidth 5
		set title "PMMA beam with $ANGLE{/Symbol \260} notch"
		
		#K2 = sqrt(fabs((\$6*\$6)-(\$4*\$4)-(\$5*\$5)))

		plot for [i=$MIN_O:$MAX_O] 'plot_files/sorted_coords_and_g_'.i.'_1.dat' u ((sqrt(\$1*\$1+\$3*\$3))/(5/sin($ANGLE*pi/180))):((sqrt(abs((\$6*\$6)-(\$4*\$4)-(\$5*\$5)))/\$4)) w linespoints lt 1 pt 1 ps 2 linewidth 6 lc i title "$ANGLE{/Symbol \260} order ".i." h1",\
		for [i=$MIN_O:$MAX_O] 'plot_files/sorted_coords_and_g_'.i.'_2.dat' u ((sqrt(\$1*\$1+\$3*\$3))/(5/sin($ANGLE*pi/180))):((sqrt(abs((\$6*\$6)-(\$4*\$4)-(\$5*\$5)))/\$4)) w linespoints lt 1 pt 2 ps 2 linewidth 6 lc i+3 title "$ANGLE{/Symbol \260} order ".i." h2",\
		for [i=$MIN_O:$MAX_O] 'plot_files/sorted_coords_and_g_'.i.'_3.dat' u ((sqrt(\$1*\$1+\$3*\$3))/(5/sin($ANGLE*pi/180))):((sqrt(abs((\$6*\$6)-(\$4*\$4)-(\$5*\$5)))/\$4)) w linespoints lt 1 pt 3 ps 2 linewidth 6 lc i+6 title "$ANGLE{/Symbol \260} order ".i." h3",\
		for [i=$MIN_O:$MAX_O] 'plot_files/sorted_coords_and_g_'.i.'_4.dat' u ((sqrt(\$1*\$1+\$3*\$3))/(5/sin($ANGLE*pi/180))):((sqrt(abs((\$6*\$6)-(\$4*\$4)-(\$5*\$5)))/\$4)) w linespoints lt 1 pt 4 ps 2 linewidth 6 lc i+9 title "$ANGLE{/Symbol \260} order ".i." h4",\
		'SIF_graphs_K2_K1.txt' u (\$1):(\$2) w line lt 3 linewidth 10 lc rgb "black" title "Lazarus et al."
EOF
fi


if [ $ANGLE == 60 ]; then
gnuplot << EOF
		#Set up gnuplot defining axis labels, key, etc
		set terminal postscript eps enhanced "Helvetica" 25 size 10,7
		set output "K2_K1_vs_x3_d.eps"
		set xlabel "x3/d"
		set ylabel "KII/KI"
		#set yrange [0.0:]
		set xrange [:1.2]
		set key above
		set key spacing 2
		set grid back lc rgb "black" linewidth 1
		set grid mxtics mytics ls 101 lc rgb "black" linewidth 5
		set title "PMMA beam with $ANGLE{/Symbol \260} notch"
		
		#K2 = sqrt(fabs((\$6*\$6)-(\$4*\$4)-(\$5*\$5)))

		plot for [i=$MIN_O:$MAX_O] 'plot_files/sorted_coords_and_g_'.i.'_1.dat' u ((sqrt(\$1*\$1+\$3*\$3))/(5/sin($ANGLE*pi/180))):((sqrt(abs((\$6*\$6)-(\$4*\$4)-(\$5*\$5)))/\$4)) w linespoints lt 1 pt 1 ps 2 linewidth 6 lc i title "$ANGLE{/Symbol \260} order ".i." h1",\
		for [i=$MIN_O:$MAX_O] 'plot_files/sorted_coords_and_g_'.i.'_2.dat' u ((sqrt(\$1*\$1+\$3*\$3))/(5/sin($ANGLE*pi/180))):((sqrt(abs((\$6*\$6)-(\$4*\$4)-(\$5*\$5)))/\$4)) w linespoints lt 1 pt 2 ps 2 linewidth 6 lc i+3 title "$ANGLE{/Symbol \260} order ".i." h2",\
		for [i=$MIN_O:$MAX_O] 'plot_files/sorted_coords_and_g_'.i.'_3.dat' u ((sqrt(\$1*\$1+\$3*\$3))/(5/sin($ANGLE*pi/180))):((sqrt(abs((\$6*\$6)-(\$4*\$4)-(\$5*\$5)))/\$4)) w linespoints lt 1 pt 3 ps 2 linewidth 6 lc i+6 title "$ANGLE{/Symbol \260} order ".i." h3",\
		for [i=$MIN_O:$MAX_O] 'plot_files/sorted_coords_and_g_'.i.'_4.dat' u ((sqrt(\$1*\$1+\$3*\$3))/(5/sin($ANGLE*pi/180))):((sqrt(abs((\$6*\$6)-(\$4*\$4)-(\$5*\$5)))/\$4)) w linespoints lt 1 pt 4 ps 2 linewidth 6 lc i+9 title "$ANGLE{/Symbol \260} order ".i." h4",\
		'SIF_graphs_K2_K1.txt' u (\$3):(\$4) w line lt 3 linewidth 10 lc rgb "black" title "Lazarus et al."
EOF
fi


if [ $ANGLE == 75 ]; then
gnuplot << EOF
		#Set up gnuplot defining axis labels, key, etc
		set terminal postscript eps enhanced "Helvetica" 25 size 10,7
		set output "K2_K1_vs_x3_d.eps"
		set xlabel "x3/d"
		set ylabel "KII/KI"
		#set yrange [0.0:]
		set xrange [:1.2]
		set key above
		set key spacing 2
		set grid back lc rgb "black" linewidth 1
		set grid mxtics mytics ls 101 lc rgb "black" linewidth 5
		set title "PMMA beam with $ANGLE{/Symbol \260} notch"
		
		#K2 = sqrt(fabs((\$6*\$6)-(\$4*\$4)-(\$5*\$5)))

		plot for [i=$MIN_O:$MAX_O] 'plot_files/sorted_coords_and_g_'.i.'_1.dat' u ((sqrt(\$1*\$1+\$3*\$3))/(5/sin($ANGLE*pi/180))):((sqrt(abs((\$6*\$6)-(\$4*\$4)-(\$5*\$5)))/\$4)) w linespoints lt 1 pt 1 ps 2 linewidth 6 lc i title "$ANGLE{/Symbol \260} order ".i." h1",\
		for [i=$MIN_O:$MAX_O] 'plot_files/sorted_coords_and_g_'.i.'_2.dat' u ((sqrt(\$1*\$1+\$3*\$3))/(5/sin($ANGLE*pi/180))):((sqrt(abs((\$6*\$6)-(\$4*\$4)-(\$5*\$5)))/\$4)) w linespoints lt 1 pt 2 ps 2 linewidth 6 lc i+3 title "$ANGLE{/Symbol \260} order ".i." h2",\
		for [i=$MIN_O:$MAX_O] 'plot_files/sorted_coords_and_g_'.i.'_3.dat' u ((sqrt(\$1*\$1+\$3*\$3))/(5/sin($ANGLE*pi/180))):((sqrt(abs((\$6*\$6)-(\$4*\$4)-(\$5*\$5)))/\$4)) w linespoints lt 1 pt 3 ps 2 linewidth 6 lc i+6 title "$ANGLE{/Symbol \260} order ".i." h3",\
		for [i=$MIN_O:$MAX_O] 'plot_files/sorted_coords_and_g_'.i.'_4.dat' u ((sqrt(\$1*\$1+\$3*\$3))/(5/sin($ANGLE*pi/180))):((sqrt(abs((\$6*\$6)-(\$4*\$4)-(\$5*\$5)))/\$4)) w linespoints lt 1 pt 4 ps 2 linewidth 6 lc i+9 title "$ANGLE{/Symbol \260} order ".i." h4",\
		'SIF_graphs_K2_K1.txt' u (\$5):(\$6) w line lt 3 linewidth 10 lc rgb "black" title "Lazarus et al."
EOF
fi

#gnuplot << EOF
#		#Set up gnuplot defining axis labels, key, etc
#		set terminal postscript eps enhanced "Helvetica" 25 size 10,7
#		set output "K3_K1_vs_x3_d.eps"
#		set xlabel "x3/d"
#		set ylabel "KIII/KI"
#		#set yrange [0.0:]
#		set xrange [:1.2]
#		set key above
#		set key spacing 3
#		set grid back lc rgb "black" linewidth 1
#		set grid mxtics mytics ls 101 lc rgb "black" linewidth 5
#		set title "PMMA beam with $ANGLE degree notch"
#
#		plot for [i=$MIN_O:$MAX_O] 'plot_files/sorted_coords_and_g_'.i.'_1.dat' u ((sqrt(\$1*\$1+\$3*\$3))/(5/sin($ANGLE*pi/180))):((\$5/\$4)) w linespoints lt 1 pt 1 ps 2 linewidth 6 lc i title "$ANGLE degrees order ".i." h1",\
#		for [i=$MIN_O:$MAX_O] 'plot_files/sorted_coords_and_g_'.i.'_2.dat' u ((sqrt(\$1*\$1+\$3*\$3))/(5/sin($ANGLE*pi/180))):((\$5/\$4)) w linespoints lt 1 pt 2 ps 2 linewidth 6 lc i+3 title "$ANGLE degrees order ".i." h2",\
#		for [i=$MIN_O:$MAX_O] 'plot_files/sorted_coords_and_g_'.i.'_3.dat' u ((sqrt(\$1*\$1+\$3*\$3))/(5/sin($ANGLE*pi/180))):((\$5/\$4)) w linespoints lt 1 pt 3 ps 2 linewidth 6 lc i+6 title "$ANGLE degrees order ".i." h3",\
#		for [i=$MIN_O:$MAX_O] 'plot_files/sorted_coords_and_g_'.i.'_4.dat' u ((sqrt(\$1*\$1+\$3*\$3))/(5/sin($ANGLE*pi/180))):((\$5/\$4)) w linespoints lt 1 pt 4 ps 2 linewidth 6 lc i+9 title "$ANGLE degrees order ".i." h4",\
#		'SIF_graphs_K3_K1.txt' u (\$1):(\$2) w line lt 3 linewidth 10 lc rgb "black" title "Lazarus"
#EOF




echo "Finished!"

open *eps


