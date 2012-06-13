#!/bin/bash
#This script reads the results from a stacking and modelling run to generate
#a table for making statistics plots

#The list variable is the fluxbinned list that the stack was run on
list=$1

#These lines read the flux bin edges from the filename
Imin=`ls  $list | grep -o [0-9]*[.][0-9][0-9][0-9] | tail -n 1`
Imax=`ls  $list | grep -o [0-9]*[.][0-9][0-9][0-9] | head -n 1`

#This line reads in the true P0 put into the simulation
Pinmed=`cat $list | ~stil/perl/onecol 4 | ~stil/perl/stat | grep "! median" |  grep -o [0-9]*[.][0-9Ee+-]*`

#The next lines read in data from the stacked images
noise=`basename $list`
stackI=`~/Code/pystack-multi/readheader.py I/noise_$list".dat*_median.fits" | grep "Median Stack Intensity" | grep -o [0-9][.][0-9Ee+-]*`
stackP=`~/Code/pystack-multi/readheader.py PI/$list"*_median.fits" | grep "Median Stack Intensity" | grep -o [0-9][.][0-9Ee+-]*`
P25=`~/Code/pystack-multi/readheader.py PI/$list"*_median.fits" | grep "Lowest Quartile Intensity" | grep -o [0-9][.][0-9Ee+-]*`
P75=`~/Code/pystack-multi/readheader.py PI/$list"*_median.fits" | grep "Highest Quartile Intensity" | grep -o [0-9][.][0-9Ee+-]*`
Nstack=`~/Code/pystack-multi/readheader.py PI/$list"*_median.fits" | grep "Number " | grep -o [0-9][0-9]*`
meanPoff=`~/Code/pystack-multi/readheader.py PI/$list"*_median.fits" | grep "Edge Median Noise" | grep -o [0-9][.][0-9Ee+-]*`
sigPoff=`~/Code/pystack-multi/readheader.py PI/$list"*_median.fits" | grep "Edge RMS Noise" | grep -o [0-9][.][0-9Ee+-]*`

#The next lines read in the background levels (with no source) of the model
modelPoff=`cat PI/noise_$list".dat_varp0_0.000000_tst" | ~stil/perl/onecol 3 | ~stil/perl/stat | grep "! median" | grep -o [0-9][.][0-9Ee+-]*`
sigPmodel=`cat PI/noise_$list".dat_varp0_0.000000_tst" | ~stil/perl/onecol 3 | ~stil/perl/stat | grep "! rms" | grep -o [0-9][.][0-9Ee+-]*`

#These variables are used to store the most likely P0 and the 16.5th and 83.5th
#percentile determined by the modelling
deltaML=99999999
delta16=99999999
delta83=99999999

#This loop find the most likely P0
for test in `ls PI/noise_$list*tst`
	do
		tempML=`cat $test | ~stil/perl/onecol 3 | ~stil/perl/percentile 0.165 | grep "Median" | grep -o [0-9][.][0-9Ee+-]*`
		temp16=`cat $test | ~stil/perl/onecol 3 | ~stil/perl/percentile 0.165 | grep "Upper percentile" | grep -o [0-9][.][0-9Ee+-]*`
		temp83=`cat $test | ~stil/perl/onecol 3 | ~stil/perl/percentile 0.165 | grep "Lower percentile" | grep -o [0-9][.][0-9Ee+-]*`
		if [ "`perl -e "if (abs($tempML-$stackP)<$deltaML){printf 1;}"`" == "1" ]
			then
			P0ML=`echo $test | grep -o [0-9][.][0-9Ee+-]*_tst | sed -e "s/\_tst//g"`
			deltaML=`perl -e "printf abs($tempML-$stackP)"`
			fi
		if [ "`perl -e "if (abs($temp16-$stackP)<$delta16){printf 1;}"`" == "1" ]
			then
			P016=`echo $test | grep -o [0-9][.][0-9Ee+-]*_tst | sed -e "s/\_tst//g"`
			delta16=`perl -e "printf abs($temp16-$stackP)"`
			fi
		if [ "`perl -e "if (abs($temp83-$stackP)<$delta83){printf 1;}"`" == "1" ]
			then
			P083=`echo $test | grep -o [0-9][.][0-9Ee+-]*_tst | sed -e "s/\_tst//g"`
			delta83=`perl -e "printf abs($temp83-$stackP)"`
			fi
	done

#These lines format the various output variables for printing
stackI=`echo $stackI | sed -e "s/e/\*10\^/" | bc -l`
stackI=`echo $stackI*1000 | bc -l`
stackP=`echo $stackP | sed -e "s/e/\*10\^/" | bc -l`
stackP=`echo $stackP*1000 | bc -l`
P25=`echo $P25 | sed -e "s/e/\*10\^/" | bc -l`
P25=`echo $P25*1000 | bc -l`
P75=`echo $P75 | sed -e "s/e/\*10\^/" | bc -l`
P75=`echo $P75*1000 | bc -l`
meanPoff=`echo $meanPoff | sed -e "s/e/\*10\^/" | bc -l`
meanPoff=`echo $meanPoff*1000 | bc -l`
sigPoff=`echo $sigPoff | sed -e "s/e/\*10\^/" | bc -l`
sigPoff=`echo $sigPoff*1000 | bc -l`
modelPoff=`echo $modelPoff | sed -e "s/e/\*10\^/" | bc -l`
modelPoff=`echo $modelPoff*1000 | bc -l`
sigPmodel=`echo $sigPmodel | sed -e "s/e/\*10\^/" | bc -l`
sigPmodel=`echo $sigPmodel*1000 | bc -l`
P0ML=`echo $P0ML*1000 | bc -l`
P016=`echo $P016*1000 | bc -l`
P083=`echo $P083*1000 | bc -l`

#This line prints the output variables in a formatted line
printf "%8.3f %8.3f %8.5f %8.5f %8.5f %8.5f %6d %8.5f %8.5f %8.5f %8.5e %8.3f %8.3f %8.3f %8.3f\n" "$Imax" "$Imin" "$stackI" "$stackP" "$P25" "$P75" "$Nstack" "$meanPoff" "$sigPoff" "$modelPoff" "$sigPmodel" "$P0ML" "$P016" "$P083" "$Pinmed"
