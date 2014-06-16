#!/bin/sh

#### NOTE:
#	1) you must be using gnuplot 4.4 or later
#	2) log files must have the naming convention log_ORDER_REFINEMENT
#	3) log files must be stored in a sub-directory called 'log_files'
#	4) ensure you enter the correct order(s), refinement level(s), entities, dimensions and scaling factor
#	5) entities can be found in the log file, search for "load_path:" and examine the next few lines
#	6) two options are available:
#		a) to generate the load vs crack area and load vs displacement plots for a single node, only one entity and one dimension are required, eg:
#			./plot_graphs.sh -o 1,2,3,4 -r 1,2 -ent 279 -dim 1 -s 1e4
#		b) to generate the load vs crack area and load vs crack mouth opening plots, two entities (either side of the crack) and all three dimensions are required eg:
#			./plot_graphs.sh -o 1,2,3,4 -r 1,2 -ent 1,785 -dim 0,1,2 -s 1e4
#	7) gnuplot settings can be changed below, see line 210 and after
#	8) after running the script, raw data files are located in a new sub-directory called "plot_files"
####

while test $# -gt 0; do
        case "$1" in
                -h|--help)
			echo "you must be using gnuplot 4.4 or later"
                        echo "create graph"
                        echo " "
                        echo "options:"
                        echo "-h, --help"
			echo "-s, --scaling_factor=SCALING_FACTOR"
			echo "-o, --orders=1,2,3,4,5 (max is 5)\""
			echo "-r, --ref=0,1,2,3,4 ... (memory and time is your limit)\""
			echo "-ent, --entity=162 or --entity=2,502\""
			echo "-dim, --dimension=0 for x, --dimension=1 for y, --dimension=2 for z"
                        exit 0
                        ;;
		-s)
			shift
			if test $# -gt 0; then
				export SCALING_FACTOR=$1
			else
				echo "no scaling factor specified"
				exit 1
			fi
			shift
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
		-ent)
                        shift
                        if test $# -gt 0; then
                                export ENTITY=$1
                        else
                                echo "no entity specified"
                                exit 1
                        fi
                        shift
                        ;;
		-dim)
                        shift
                        if test $# -gt 0; then
                                export DIMENSION=$1
                        else
                                echo "no dimension specified"
                                exit 1
                        fi
                        shift
                        ;;
	        esac
done

if [ -z ${SCALING_FACTOR+x} ] ; then
  echo "Scaling factor is unset, plot_graph --help"
  exit 1
fi
if [ -z ${ORDERS+x} ]; then 
  echo "Approximation order is unset, plot_graph --help"
  exit 1
fi
if [ -z ${REF+x} ]; then 
  echo "Refinement levels is unset, plot_graph --help"
  exit 1
fi
if [ -z ${ENTITY+x} ]; then 
  echo "Entity is unset, plot_graph --help"
  exit 1
fi
if [ -z ${DIMENSION+x} ]; then 
  echo "Dimension is unset, plot_graph --help"
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

#Determine entities to be plotted
ENT_A=$(echo "$ENTITY" | sed -e 's/,/ /g' |  awk -F" " '{print $1}')
ENT_B=$(echo "$ENTITY" | sed -e 's/,/ /g' | sed -e 's/^.* //')

#Determine dimensions to be plotted
DIM_A=$(echo "$DIMENSION" | sed -e 's/,/ /g' |  awk -F" " '{print $1}')
DIM_B=$(echo "$DIMENSION" | sed -e 's/,/ /g' |  awk -F" " '{print $2}')
DIM_C=$(echo "$DIMENSION" | sed -e 's/,/ /g' |  awk -F" " '{print $3}')

echo "Please wait, generating graphs..."

###NOTE###
#If only one entity is declared by the user, the script will only generate the crack area vs load factor and displacement vs load factor graphs
#If two entities are declared by the user, the script will assume these two entities are either side of a crack opening and generate a crack mouth opening vs load factor graph as well as crack area vs load factor and displacement vs load factor graphs
###NOTE###


if [ $ENT_A != $ENT_B ]; then
	#Loop over orders and loop over refinements to extract data from log files
	for order in `echo $ORDERS | sed -e 's/,/ /g'`
	do 
		for ref_level in `echo $REF | sed -e 's/,/ /g'`
		do
	
		#Extract crack area and load factor from log file and write to new file, suppressing output to terminal screen
		awk '/load_path:/ {print $4,$6}' log_files/log_"$order"_"$ref_level" | tee -a plot_files/load_path_CrackArea_vs_Lambda_"$order"_"$ref_level".dat >/dev/null
	
		#Extract displacement and load factor for entity A and dimension from log file and write to new file, suppressing output to terminal screen
		echo "0 0" | tee plot_files/load_path_DispX"$ENT_A"_vs_Lambda_"$order"_"$ref_level".dat >/dev/null && \
			awk -v ENTITY=$ENT_A -v DIM=$DIM_A '/load_path_disp ent [0-9]* dim [0]*/ {if($3==ENTITY && $5==DIM) print $13,$15}' log_files/log_"$order"_"$ref_level" | \
			tee -a plot_files/load_path_DispX"$ENT_A"_vs_Lambda_"$order"_"$ref_level".dat >/dev/null
		echo "0 0" | tee plot_files/load_path_DispY"$ENT_A"_vs_Lambda_"$order"_"$ref_level".dat >/dev/null && \
			awk -v ENTITY=$ENT_A -v DIM=$DIM_B '/load_path_disp ent [0-9]* dim [0]*/ {if($3==ENTITY && $5==DIM) print $13,$15}' log_files/log_"$order"_"$ref_level" | \
			tee -a plot_files/load_path_DispY"$ENT_A"_vs_Lambda_"$order"_"$ref_level".dat >/dev/null
		echo "0 0" | tee plot_files/load_path_DispZ"$ENT_A"_vs_Lambda_"$order"_"$ref_level".dat >/dev/null && \
			awk -v ENTITY=$ENT_A -v DIM=$DIM_C '/load_path_disp ent [0-9]* dim [0]*/ {if($3==ENTITY && $5==DIM) print $13,$15}' log_files/log_"$order"_"$ref_level" | \
			tee -a plot_files/load_path_DispZ"$ENT_A"_vs_Lambda_"$order"_"$ref_level".dat >/dev/null


		#Extract displacement and load factor for entity B and dimension from log file and write to new file, suppressing output to terminal screen
		echo "0 0" | tee plot_files/load_path_DispX"$ENT_B"_vs_Lambda_"$order"_"$ref_level".dat >/dev/null && \
			awk -v ENTITY=$ENT_B -v DIM=$DIM_A '/load_path_disp ent [0-9]* dim [0]*/ {if($3==ENTITY && $5==DIM) print $13,$15}' log_files/log_"$order"_"$ref_level" | \
			tee -a plot_files/load_path_DispX"$ENT_B"_vs_Lambda_"$order"_"$ref_level".dat >/dev/null
		echo "0 0" | tee plot_files/load_path_DispY"$ENT_B"_vs_Lambda_"$order"_"$ref_level".dat >/dev/null && \
			awk -v ENTITY=$ENT_B -v DIM=$DIM_B '/load_path_disp ent [0-9]* dim [0]*/ {if($3==ENTITY && $5==DIM) print $13,$15}' log_files/log_"$order"_"$ref_level" | \
			tee -a plot_files/load_path_DispY"$ENT_B"_vs_Lambda_"$order"_"$ref_level".dat >/dev/null
		echo "0 0" | tee plot_files/load_path_DispZ"$ENT_B"_vs_Lambda_"$order"_"$ref_level".dat >/dev/null && \
			awk -v ENTITY=$ENT_B -v DIM=$DIM_C '/load_path_disp ent [0-9]* dim [0]*/ {if($3==ENTITY && $5==DIM) print $13,$15}' log_files/log_"$order"_"$ref_level" | \
			tee -a plot_files/load_path_DispZ"$ENT_B"_vs_Lambda_"$order"_"$ref_level".dat >/dev/null

		cd plot_files
		paste load_path_DispX"$ENT_A"_vs_Lambda_"$order"_"$ref_level".dat load_path_DispY"$ENT_A"_vs_Lambda_"$order"_"$ref_level".dat  load_path_DispZ"$ENT_A"_vs_Lambda_"$order"_"$ref_level".dat > CMO_"$ENT_A"_vs_Lambda_"$order"_"$ref_level".dat
		awk '{print $1" "$3" "$5" "$6}'  CMO_"$ENT_A"_vs_Lambda_"$order"_"$ref_level".dat | tee CMO_"$ENT_A"_"$order"_"$ref_level".dat >/dev/null
		paste load_path_DispX"$ENT_B"_vs_Lambda_"$order"_"$ref_level".dat load_path_DispY"$ENT_B"_vs_Lambda_"$order"_"$ref_level".dat  load_path_DispZ"$ENT_B"_vs_Lambda_"$order"_"$ref_level".dat > CMO_"$ENT_B"_vs_Lambda_"$order"_"$ref_level".dat
		awk '{print $1" "$3" "$5" "$6}'  CMO_"$ENT_B"_vs_Lambda_"$order"_"$ref_level".dat | tee CMO_"$ENT_B"_"$order"_"$ref_level".dat >/dev/null

		paste CMO_"$ENT_A"_"$order"_"$ref_level".dat CMO_"$ENT_B"_"$order"_"$ref_level".dat > Disp_vs_Lambda_"$order"_"$ref_level".dat
		rm CMO_*.dat
		rm load_path_Disp*.dat
		cd ..
		done
	done
fi
if [ $ENT_A == $ENT_B ]; then
	#Loop over orders and loop over refinements to extract data from log files
	for order in `echo $ORDERS | sed -e 's/,/ /g'`
	do 
		for ref_level in `echo $REF | sed -e 's/,/ /g'`
		do
	
		#Extract crack area and load factor from log file and write to new file, suppressing output to terminal screen
		awk '/load_path:/ {print $4,$6}' log_files/log_"$order"_"$ref_level" | tee -a plot_files/load_path_CrackArea_vs_Lambda_"$order"_"$ref_level".dat >/dev/null
	
		#Extract displacement and load factor for the current (only) entity and dimension from log file and write to new file, suppressing output to terminal screen
		echo "0 0" | tee plot_files/load_path_DispY_vs_Lambda_"$order"_"$ref_level".dat >/dev/null && \
			awk -v ENTITY=$ENT_A -v DIM=$DIM_A '/load_path_disp ent [0-9]* dim [0]*/ {if($3==ENTITY && $5==DIM) print $13,$15}' log_files/log_"$order"_"$ref_level" | \
			tee -a plot_files/load_path_DispY_vs_Lambda_"$order"_"$ref_level".dat >/dev/null
		done

	done
fi

#Remove previous plot outputs
rm output*.eps
if [ $ENT_A != $ENT_B ]; then
	#Open gnuplot *Load-Crack mouth opening Plot*
	gnuplot << EOF
		#Set up gnuplot defining axis labels, key, etc
		set terminal postscript eps enhanced "Helvetica" 50 size 20,14
		set output "output_cmo.eps"
		set xlabel "Crack mouth opening norm (mm)"
		set ylabel "Load (N)"
		#set yrange [0.0:6000]
		set xrange [0.0:1e-3]
		set key above
		set key spacing 3
		set grid back lc rgb "black" linewidth 1
		set grid mxtics mytics ls 101 lc rgb "black" linewidth 5

		#(sqrt((\$1-\$5)**2+(\$2-\$6)**2+(\$3-\$7)**2))
	
		#Depending on the maximum refinement level assessed, for each order assessed, plot the load-crack mouth opening graph
		if ($MAX_R == 1) \
			plot for [i=$MIN_O:$MAX_O] 'plot_files/Disp_vs_Lambda_'.i.'_1.dat' u (sqrt((\$1-\$5)**2+(\$2-\$6)**2+(\$3-\$7)**2)):(abs(\$8*$SCALING_FACTOR)) w linespoints lt 5 pt 1 ps 2 linewidth 6 lc i title "order ".i." h1"
		if ($MAX_R == 2) \
			plot for [i=$MIN_O:$MAX_O] 'plot_files/Disp_vs_Lambda_'.i.'_1.dat' u (sqrt((\$1-\$5)**2+(\$2-\$6)**2+(\$3-\$7)**2)):(abs(\$8*$SCALING_FACTOR)) w linespoints lt 5 pt 1 ps 2 linewidth 6 lc i title "order ".i." h1",\
			for [i=$MIN_O:$MAX_O] 'plot_files/Disp_vs_Lambda_'.i.'_2.dat' u (sqrt((\$1-\$5)**2+(\$2-\$6)**2+(\$3-\$7)**2)):(abs(\$8*$SCALING_FACTOR)) w linespoints lt 1 pt 2 ps 2 linewidth 6 lc i title "order ".i." h2"
		if ($MAX_R == 3) \
			plot for [i=$MIN_O:$MAX_O] 'plot_files/Disp_vs_Lambda_'.i.'_1.dat' u (sqrt((\$1-\$5)**2+(\$2-\$6)**2+(\$3-\$7)**2)):(abs(\$8*$SCALING_FACTOR)) w linespoints lt 5 pt 1 ps 2 linewidth 6 lc i title "order ".i." h1",\
			for [i=$MIN_O:$MAX_O] 'plot_files/Disp_vs_Lambda_'.i.'_2.dat' u (sqrt((\$1-\$5)**2+(\$2-\$6)**2+(\$3-\$7)**2)):(abs(\$8*$SCALING_FACTOR)) w linespoints lt 1 pt 2 ps 2 linewidth 6 lc i title "order ".i." h2",\
			for [i=$MIN_O:$MAX_O] 'plot_files/Disp_vs_Lambda_'.i.'_3.dat' u (sqrt((\$1-\$5)**2+(\$2-\$6)**2+(\$3-\$7)**2)):(abs(\$8*$SCALING_FACTOR)) w linespoints lt 2 pt 3 ps 2 linewidth 6 lc i title "order ".i." h3"
		if ($MAX_R == 4) \
			plot for [i=$MIN_O:$MAX_O] 'plot_files/Disp_vs_Lambda_'.i.'_1.dat' u (sqrt((\$1-\$5)**2+(\$2-\$6)**2+(\$3-\$7)**2)):(abs(\$8*$SCALING_FACTOR)) w linespoints lt 5 pt 1 ps 2 linewidth 6 lc i title "order ".i." h1",\
			for [i=$MIN_O:$MAX_O] 'plot_files/Disp_vs_Lambda_'.i.'_2.dat' u (sqrt((\$1-\$5)**2+(\$2-\$6)**2+(\$3-\$7)**2)):(abs(\$8*$SCALING_FACTOR)) w linespoints lt 1 pt 2 ps 2 linewidth 6 lc i title "order ".i." h2",\
			for [i=$MIN_O:$MAX_O] 'plot_files/Disp_vs_Lambda_'.i.'_3.dat' u (sqrt((\$1-\$5)**2+(\$2-\$6)**2+(\$3-\$7)**2)):(abs(\$8*$SCALING_FACTOR)) w linespoints lt 2 pt 3 ps 2 linewidth 6 lc i title "order ".i." h3",\
			for [i=$MIN_O:$MAX_O] 'plot_files/Disp_vs_Lambda_'.i.'_4.dat' u (sqrt((\$1-\$5)**2+(\$2-\$6)**2+(\$3-\$7)**2)):(abs(\$8*$SCALING_FACTOR)) w linespoints lt 3 pt 4 ps 2 linewidth 6 lc i title "order ".i." h4"	
#Exit gunplot
EOF
fi

if [ $ENT_A == $ENT_B ]; then
	#Open gnuplot *Load-Displacement Plot*
	gnuplot << EOF
		#Set up gnuplot defining axis labels, key, etc
		set terminal postscript eps enhanced "Helvetica" 50 size 20,14
		set output "output_disp.eps"
		set xlabel "Displacement (mm)"
		set ylabel "Load (N)"
		#set yrange [0.0:6000]
		set xrange [0.0:1e-3]
		set key above
		set key spacing 3
		set grid back lc rgb "black" linewidth 1
		set grid mxtics mytics ls 101 lc rgb "black" linewidth 5
	
		#Depending on the maximum refinement level assessed, for each order assessed, plot the load-displacement graph
		if ($MAX_R == 1) \
			plot for [i=$MIN_O:$MAX_O] 'plot_files/load_path_DispY_vs_Lambda_'.i.'_1.dat' u (abs(\$1)):(abs(\$2*$SCALING_FACTOR)) w linespoints lt 5 pt 1 ps 2 linewidth 6 lc i title "order ".i." h1"
		if ($MAX_R == 2) \
			plot for [i=$MIN_O:$MAX_O] 'plot_files/load_path_DispY_vs_Lambda_'.i.'_1.dat' u (abs(\$1)):(abs(\$2*$SCALING_FACTOR)) w linespoints lt 5 pt 1 ps 2 linewidth 6 lc i title "order ".i." h1",\
			for [i=$MIN_O:$MAX_O] 'plot_files/load_path_DispY_vs_Lambda_'.i.'_2.dat' u (abs(\$1)):(abs(\$2*$SCALING_FACTOR)) w linespoints lt 1 pt 2 ps 2 linewidth 6 lc i title "order ".i." h2"
		if ($MAX_R == 3) \
			plot for [i=$MIN_O:$MAX_O] 'plot_files/load_path_DispY_vs_Lambda_'.i.'_1.dat' u (abs(\$1)):(abs(\$2*$SCALING_FACTOR)) w linespoints lt 5 pt 1 ps 2 linewidth 6 lc i title "order ".i." h1",\
			for [i=$MIN_O:$MAX_O] 'plot_files/load_path_DispY_vs_Lambda_'.i.'_2.dat' u (abs(\$1)):(abs(\$2*$SCALING_FACTOR)) w linespoints lt 1 pt 2 ps 2 linewidth 6 lc i title "order ".i." h2",\
			for [i=$MIN_O:$MAX_O] 'plot_files/load_path_DispY_vs_Lambda_'.i.'_3.dat' u (abs(\$1)):(abs(\$2*$SCALING_FACTOR)) w linespoints lt 2 pt 3 ps 2 linewidth 6 lc i title "order ".i." h3"
		if ($MAX_R == 4) \
			plot for [i=$MIN_O:$MAX_O] 'plot_files/load_path_DispY_vs_Lambda_'.i.'_1.dat' u (abs(\$1)):(abs(\$2*$SCALING_FACTOR)) w linespoints lt 5 pt 1 ps 2 linewidth 6 lc i title "order ".i." h1",\
			for [i=$MIN_O:$MAX_O] 'plot_files/load_path_DispY_vs_Lambda_'.i.'_2.dat' u (abs(\$1)):(abs(\$2*$SCALING_FACTOR)) w linespoints lt 1 pt 2 ps 2 linewidth 6 lc i title "order ".i." h2",\
			for [i=$MIN_O:$MAX_O] 'plot_files/load_path_DispY_vs_Lambda_'.i.'_3.dat' u (abs(\$1)):(abs(\$2*$SCALING_FACTOR)) w linespoints lt 2 pt 3 ps 2 linewidth 6 lc i title "order ".i." h3",\
			for [i=$MIN_O:$MAX_O] 'plot_files/load_path_DispY_vs_Lambda_'.i.'_4.dat' u (abs(\$1)):(abs(\$2*$SCALING_FACTOR)) w linespoints lt 3 pt 4 ps 2 linewidth 6 lc i title "order ".i." h4"		
#Exit gunplot
EOF
fi


#Open gnuplot *Load-Crack Area Plot*
gnuplot << EOF	
	#Set up gnuplot defining axis labels, key, etc
	set terminal postscript eps enhanced "Helvetica" 50 size 20,14
	set output "output_area.eps"
	set xlabel "Crack Area (mm^2)" 
	set ylabel "Load (N)"
	#set yrange [0.0:6000]
	#set xrange [:13000]
	set key above
	set key spacing 3
	set grid back lc rgb "black" linewidth 1
	set grid mxtics mytics ls 101 lc rgb "black" linewidth 5

	#Depending on the maximum refinement level assessed, for each order assessed, plot the load-crack area graph
	if ($MAX_R == 1) \
		plot for [i=$MIN_O:$MAX_O] 'plot_files/load_path_CrackArea_vs_Lambda_'.i.'_1.dat' u (abs(\$1)):(abs(\$2*$SCALING_FACTOR)) w linespoints lt 5 pt 1 ps 2 linewidth 6 lc i title "order ".i." h1"
	if ($MAX_R == 2) \
		plot for [i=$MIN_O:$MAX_O] 'plot_files/load_path_CrackArea_vs_Lambda_'.i.'_1.dat' u (abs(\$1)):(abs(\$2*$SCALING_FACTOR)) w linespoints lt 5 pt 1 ps 2 linewidth 6 lc i title "order ".i." h1",\
		for [i=$MIN_O:$MAX_O] 'plot_files/load_path_CrackArea_vs_Lambda_'.i.'_2.dat' u (abs(\$1)):(abs(\$2*$SCALING_FACTOR)) w linespoints lt 1 pt 2 ps 2 linewidth 6 lc i title "order ".i." h2"
	if ($MAX_R == 3) \
		plot for [i=$MIN_O:$MAX_O] 'plot_files/load_path_CrackArea_vs_Lambda_'.i.'_1.dat' u (abs(\$1)):(abs(\$2*$SCALING_FACTOR)) w linespoints lt 5 pt 1 ps 2 linewidth 6 lc i title "order ".i." h1",\
		for [i=$MIN_O:$MAX_O] 'plot_files/load_path_CrackArea_vs_Lambda_'.i.'_2.dat' u (abs(\$1)):(abs(\$2*$SCALING_FACTOR)) w linespoints lt 1 pt 2 ps 2 linewidth 6 lc i title "order ".i." h2",\
		for [i=$MIN_O:$MAX_O] 'plot_files/load_path_CrackArea_vs_Lambda_'.i.'_3.dat' u (abs(\$1)):(abs(\$2*$SCALING_FACTOR)) w linespoints lt 2 pt 3 ps 2 linewidth 6 lc i title "order ".i." h3"
	if ($MAX_R == 4) \
		plot for [i=$MIN_O:$MAX_O] 'plot_files/load_path_CrackArea_vs_Lambda_'.i.'_1.dat' u (abs(\$1)):(abs(\$2*$SCALING_FACTOR)) w linespoints lt 5 pt 1 ps 2 linewidth 6 lc i title "order ".i." h1",\
		for [i=$MIN_O:$MAX_O] 'plot_files/load_path_CrackArea_vs_Lambda_'.i.'_2.dat' u (abs(\$1)):(abs(\$2*$SCALING_FACTOR)) w linespoints lt 1 pt 2 ps 2 linewidth 6 lc i title "order ".i." h2",\
		for [i=$MIN_O:$MAX_O] 'plot_files/load_path_CrackArea_vs_Lambda_'.i.'_3.dat' u (abs(\$1)):(abs(\$2*$SCALING_FACTOR)) w linespoints lt 2 pt 3 ps 2 linewidth 6 lc i title "order ".i." h3",\
		for [i=$MIN_O:$MAX_O] 'plot_files/load_path_CrackArea_vs_Lambda_'.i.'_4.dat' u (abs(\$1)):(abs(\$2*$SCALING_FACTOR)) w linespoints lt 3 pt 4 ps 2 linewidth 6 lc i title "order ".i." h4"
#Exit gnuplot
EOF

echo "Finished!"

#Open gnuplot output files
open output*.eps

#End
