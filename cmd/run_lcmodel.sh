#!/bin/bash

echo

curdir=$(pwd)

echo Current Directory is: 
echo $curdir

DIR="LCModel/Block"
if [ -d "$DIR" ]; then
  # Take action if $DIR exists. #
  echo "fMRS block analysis detected ..."
  blockList=($( ls $DIR))

  x=0
  for n in ${blockList[@]}; do
	echo Analysing Block $n
	cd $curdir/LCModel/Block/$n
	for x in $(ls *.CONTROL); do
		echo $x
		~/.lcmodel/bin/lcmodel < $x
		fname="${x%%.*}"
		ps2pdf $fname.PS $fname.pdf
	done
  done
else
	echo "No fMRS block analysis detected ..."		
fi

# Analyse whole exp #
x=0
echo Analysing Data
cd $curdir/LCModel
	for x in $(ls *.CONTROL); do
	echo $x
	~/.lcmodel/bin/lcmodel < $x
	fname="${x%%.*}"
	ps2pdf $fname.PS $fname.pdf
done
