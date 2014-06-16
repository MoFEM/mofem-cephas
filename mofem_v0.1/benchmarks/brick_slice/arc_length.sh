#!/bin/sh

ARGS=$@

while test $# -gt 0; do
        case "$1" in
                -h|--help)
                        echo "crack propagation"
                        echo " "
                        echo "arc_length.sh [arguments]"
                        echo " "
                        echo "options:"
                        echo "-h, --help"
                        echo "-f, --file_name=FILE"
			echo "-a, --da=CRACK_AREA_INCREMENT"
			echo "-n, --nb_steps=NUMBER_OF_STEPS"
			echo "-g, --gc=GRIFFITH_ENERGY"
			echo "-d, --dir=DIR_NAME"
			echo "-p, --proc=NB_OF_PROCESSORS"
      			echo "-t, --rol=TOLERANCE"
			echo "-i. --its_d=DESIRED_NUMBER_OF_ITS"
			echo "-s, --spit_iter=NB_SPLIT"
			echo "--pc_factor_mat_solver_package=[superlu_dist,mumps]"
                        exit 0
                        ;;
                -f)
                        shift
                        if test $# -gt 0; then
                                export FILE=$1
                        else
                                echo "no file name specified"
                                exit 1
                        fi
                        shift
                        ;;
                --file_name*)
                        export FILE=`echo $1 | sed -e 's/^[^=]*=//g'`
                        shift
                        ;;
                -d)
                        shift
                        if test $# -gt 0; then
                                export DIR_NAME=$1
                        else
                                echo "no file name specified"
                                exit 1
                        fi
                        shift
                        ;;
                --dir*)
                        export DIR_NAME=`echo $1 | sed -e 's/^[^=]*=//g'`
                        shift
                        ;;

		-a)
                        shift
                        if test $# -gt 0; then
                                export DA=$1
                        else
                                echo "no crack area increment energy specified"
                                exit 1
                        fi
                        shift
                        ;;
                --da*)
                        export DA=`echo $1 | sed -e 's/^[^=]*=//g'`
                        shift
                        ;;
		-n)
                        shift
                        if test $# -gt 0; then
                                export NBSTEPS=$1
                        else
                                echo "no number of steps specified"
                                exit 1
                        fi
                        shift
                        ;;
                --nb_steps*)
                        export NBSTEPS=`echo $1 | sed -e 's/^[^=]*=//g'`
                        shift
                        ;;

		-g)
                        shift
                        if test $# -gt 0; then
                                export GC=$1
                        else
                                echo "no griffith energy specified"
                                exit 1
                        fi
                        shift
                        ;;
                --gc*)
                        export GC=`echo $1 | sed -e 's/^[^=]*=//g'`
                        shift
                        ;;

		-p)
                        shift
                        if test $# -gt 0; then
                                export NB_PROC=$1
                        else
                                echo "no processors number specified"
                                exit 1
                        fi
                        shift
                        ;;
                --proc*)
                        export NB_PROC=`echo $1 | sed -e 's/^[^=]*=//g'`
                        shift
                        ;;
		-t)
                        shift
                        if test $# -gt 0; then
                                export TOL=$1
                        else
                                echo "no tolerance specified"
                                exit 1
                        fi
                        shift
                        ;;
                --tol*)
                        export TOL=`echo $1 | sed -e 's/^[^=]*=//g'`
                        shift
                        ;;
		-s)
                        shift
                        if test $# -gt 0; then
                                export NBSPLIT=$1
                        else
                                echo "number of splits"
                                exit 1
                        fi
                        shift
                        ;;
                --split_iter*)
                        export NBSPLIT=`echo $1 | sed -e 's/^[^=]*=//g'`
                        shift
                        ;;
		-i)
                        shift
                        if test $# -gt 0; then
                                export ITSD=$1
                        else
                                echo "no desired number of iterations specified"
                                exit 1
                        fi
                        shift
                        ;;
                --its_d*)
                        export ITSD=`echo $1 | sed -e 's/^[^=]*=//g'`
                        shift
                        ;;
		--snes_linesearch_type*) 
                        export LINESEARCHTYPE=`echo $1 | sed -e 's/^[^=]*=//g'`
                        shift
                        ;;
		--pc_factor_mat_solver_package*) 
                        export PCFACTOR=`echo $1 | sed -e 's/^[^=]*=//g'`
                        shift
                        ;;
		--my_nb_sub_steps*) 
                        export NBSUBSTEPS=`echo $1 | sed -e 's/^[^=]*=//g'`
                        shift
                        ;;
                *)
                        break
                        ;;
        esac
done

if [ -z \$\{FILE+x\} ]; then 
  echo "Mesh file is unset, arc_length --help"
fi

if [ -z ${DA+x} ]; then 
  echo "Crack area increment is unset, arc_length --help"
  exit 1
fi

if [ -z ${NBSTEPS+x} ]; then 
  echo "Number of steps is unset, arc_length --help"
  exit 1
fi

if [ -z ${GC+x} ]; then 
  echo "Griffith energy is unset, arc_length --help"
  exit 1
fi

if [ -z ${DIR_NAME+x} ]; then
  DIR_SURFIX=`date "+DATE-%Y-%m-%d_TIME_%H_%M_%S"`
  DIR_PREFIX=`basename $FILE .h5m`
  export DIR_NAME=arc_length_"$DIR_PREFIX"_"$DIR_SURFIX"
fi

if [ ! -d $DIR_NAME ]; then
  mkdir $DIR_NAME
fi

if [ -z ${NB_PROC+x} ]; then 
  NB_PROC=4
fi

if [ -z ${TOL+x} ]; then 
  TOL=1e-9
fi

if [ -z ${ITSD+x} ]; then 
  ITSD=6
fi

if [ -z ${LINESEARCHTYPE+x} ]; then 
  LINESEARCHTYPE="basic"
fi

if [ -z ${NBSUBSTEPS+x} ]; then 
  NBSUBSTEPS=20
fi

if [ -z ${PCFACTOR+x} ]; then 
  PCFACTOR=superlu_dist
fi

if [ -z ${NBSPLIT+x} ]; then
  NBSPLIT=0
fi

BIN_PATH="/Users/rossmackenzie/mofem/examples/material_forces"
MPIRUN="/opt/build_for_gcc-mp-4.4/local/bin/mpirun -np $NB_PROC"

cp $FILE $DIR_NAME/
cd $DIR_NAME

echo "Line command ..." | tee log
echo $0 $ARGS | tee -a log
echo | tee -a log
echo "Satrt time: "`date` | tee -a log

COMMAND="$MPIRUN $BIN_PATH/arc_length_material_coupled \
  -my_file `basename $FILE` \
  -ksp_type fgmres -pc_type lu -pc_factor_mat_solver_package $PCFACTOR \
  -snes_type newtonls -snes_linesearch_type $LINESEARCHTYPE -snes_linesearch_monitor \
  -ksp_atol 1e-14 -ksp_rtol 1e-12 -ksp_max_it 500\
  -snes_max_it 25 -snes_atol 1e-9 -snes_rtol 1e-9 -snes_stol 0 \
  -snes_monitor -snes_converged_reason -ksp_monitor \
  -my_tol $TOL -my_its_d $ITSD \
  -my_gc $GC -my_alpha2 1e-3 -my_gamma 0. \
  -my_da $DA -my_load_steps $NBSTEPS \
  -my_nb_sub_steps $NBSUBSTEPS \
  -my_odd_face_split $NBSPLIT
  -mat_ascii_output_large \
  -on_error_attach_debugger \
  -log_summary"

echo | tee -a log
echo "Start calculations ..." | tee -a log
echo $COMMAND | tee -a log
echo | tee -a log

exec $COMMAND 2>&1 | tee -a log

echo | tee -a log
echo "End time: "`date` | tee -a log


if [ $? -ne 0 ]; then
  exit $?
fi


cd ..
