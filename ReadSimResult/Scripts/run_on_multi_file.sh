#!/bin/bash

inputdir=/eos/home-m/mikumar/HGCAl_Validation/Extended2026D86
pydir=$PWD/ReadSimResult/SimTrackAna/python
for i in `seq 3 9`
do
  echo Processing loop $i with file step2_${i}.root
  if [ -f $PWD/step2.root ] ; then
      rm $PWD/step2.root
  fi
  ln -s $inputdir/step2_${i}.root $PWD/step2.root 
  #cmsRun $pydir/CellHitSum_cfg.py #-n 4
  cmsRun $pydir/digiHit_cfg.py #-n 4
  mv geantoutput.root geantoutput_${i}.root
done

ls $PWD/geantoutput_*.root > /tmp/mikumar/geanoutputfile.txt
#source ~/scripts/addhisto_file.sh /tmp/idas/fl.txt
hadd merged.root geantoutput_*.root
#rm geantoutput_*.root
mv merged.root geantoutput_merged.root
