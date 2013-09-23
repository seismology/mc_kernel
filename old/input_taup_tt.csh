#!/bin/csh -f

set depth = `grep "source depth" simulation.info_fwd |awk '{print $1}' |sed 's/\./ /'|awk '{print $1}'`
set dist = `grep "theta/colatitude" input_kernel.dat |awk '{print $1}' |sed 's/\./ /'|awk '{print $1}'`
set period = `grep "filtering period" input_kernel.dat |awk '{print $1}' |sed 's/\./ /'|awk '{print $1}'`
echo "Distance: $dist , Depth: $depth , Period: $period "

rm -f taup_tt.dat taup_phase.dat taup_tt_phase.dat period.dat
if ( $dist <100 ) then 
    set phases = ( "P" "PP" "PcP" "S" "SS" "ScS" )
    foreach i (${phases})
        echo "computing traveltime for" $i
	taup_time -mod prem -h $depth -ph $i, -deg $dist | grep $depth | awk '{print $4}' |head -n 1  >> taup_tt.dat
	taup_time -mod prem -h $depth -ph $i, -deg $dist |grep $depth | awk '{print $3}' |head -n 1 >> taup_phase.dat
	echo $period >> period.dat
    end
else 
    set phases = ( "Pdiff" "PKiKP" "PP" "SKS" "SKIKS" "Sdiff" )
    foreach i (${phases})
        echo "computing traveltime for" $i
        set length_tt = ` taup_time -mod prem -h $depth -ph $i, -deg $dist | grep $depth  | awk '{print $4}'|wc -l `
        if ( $length_tt > 0 ) then
          echo $period >> period.dat
	  taup_time -mod prem -h $depth -ph $i, -deg $dist | grep $depth  | awk '{print $4}' |head -n 1 >> taup_tt.dat
          taup_time -mod prem -h $depth -ph $i, -deg $dist |grep $depth  | awk '{print $3}' |head -n 1 >> taup_phase.dat
        endif
    end
endif

paste taup_tt.dat taup_phase.dat period.dat> taup_tt_phase.dat
wc -l taup_tt_phase.dat |awk '{print $1}' > input_timewindows_taup.dat
cat taup_tt_phase.dat |awk '{print $1-0.75*$3, $1+0.75*$3,"   ", $2 }' >> input_timewindows_taup.dat


