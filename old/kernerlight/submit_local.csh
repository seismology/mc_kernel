#!/bin/csh -f

echo "";echo "::::::::::::::::::::::: k e r n e r ::::::::::::::::::::::::::::"; echo ""
set homedir = $PWD

if ( -d $1) then
    echo " Run or directory" $1 "exists......."
    ls $1
    exit
endif

set datapath = 'Data'
set infopath = 'Info'

# determine number of simulations
set nsim = `grep "number of simulations" input_kernel.dat|awk '{print $1}'`
echo "Number of simulations: $nsim"

set fwdpath = `head -n 1 input_sem.dat  | cut -c 2- |awk '{print $1}'  `
set bwdpath = `head -n 2 input_sem.dat  | tail -n 1 | cut -c 2- |awk '{print $1}'  `

if ( $nsim > 1 ) then 
  set fwdpath =  `echo $fwdpath"MZZ"`
endif

echo "Forward simulation path: $fwdpath "
echo "Backward simulation path: $bwdpath"

if ( ! -d $fwdpath) then
  echo "Forward path does not exist"
  exit
endif

if ( ! -d $bwdpath) then
  echo "Backward path does not exist"
  exit
endif

cp -f $fwdpath/*.bm .

# only copy these files if changed, to avoid unnecessary compilation time:
set files=(mesh_params.h mesh_params_kernel.h Code/background_models.f90) 
foreach file ($files) 
    if ( -f $file ) then
        # if $file exists, check for differences:
        set changed = `diff $fwdpath/$file $file | wc -l`
        if ( $changed > 0 ) then 
            cp -f $fwdpath/$file .
        endif
    else
        # if $file does not exist (clean checkout)
        cp -f $fwdpath/$file .
    endif
end


if ( { make all } == 0 ) then
  echo "Compilation failed, please check the errors."
  exit
endif

mkdir $1
cd $1

mkdir $datapath
echo "creating $datapath"

mkdir $infopath
echo "creating $infopath"

mkdir Code
cp -p $homedir/*.f90 Code
cp -p $homedir/Makefile Code
cp -p $homedir/makemake.pl Code
cp $homedir/xkerner .
cp $homedir/input_*dat .
cp $homedir/inparam_cs .
cp $homedir/*.bm .
cp $homedir/mesh_params.h .
cp $homedir/mesh_params_kernel.h .
cp $homedir/param_snaps .
cp $homedir/input_taup_tt.csh .
rm -f input_timewindows_taup.dat

/bin/cp -f $fwdpath/Info/backgroundmodel_kmscale.dat0000 $homedir/$1/

/bin/cp -f $fwdpath/sourceparams.dat sourceparams_fwd.dat
/bin/cp -f $bwdpath/sourceparams.dat sourceparams_bwd.dat
/bin/cp -f $fwdpath/inparam inparam_fwd
/bin/cp -f $bwdpath/inparam inparam_bwd
/bin/cp -f $fwdpath/simulation.info simulation.info_fwd
/bin/cp -f $bwdpath/simulation.info simulation.info_bwd

if ( ${#argv} > 1 ) then
  set nodnum = $2
else
  set nodnum = `grep nproc_mesh mesh_params.h |awk '{print $6}'`
endif

echo " nproc= "$nodnum > Data/simulation.info
echo "preparing job on $nodnum nodes..."

#echo "# Sample PBS for parallel jobs" >runkerner_np$nodnum
#echo "#PBS -l nodes=$nodnum,walltime=400:00:00" >> runkerner_np$nodnum
#echo "ulimit -s unlimited" >> runkerner_np$nodnum
##echo "export P4_GLOBMEMSIZE=167772160000" >> runkerner_np$nodnum
#echo "cd $PWD" >> runkerner_np$nodnum
#echo "/usr/local/bin/mpiexec kerner > OUTPUT_$1" >> runkerner_np$nodnum

# ./input_taup_tt.csh

set outname = `echo $1 | sed 's/\// /g'|awk '{print NR, $NF}'|awk '{print $2}'`
set outputname = `echo "OUTPUT_$outname"`

mpiexec -np $nodnum xkerner > $outputname &
echo "Running in directory $1 "
echo "To check on the run:"
echo "tail -f $1/$outputname"

#qsub runkerner_np$nodnum
#qstat -n
