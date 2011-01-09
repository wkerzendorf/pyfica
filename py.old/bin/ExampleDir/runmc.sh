#!/bin/sh
#

ulimit -s unlimited

#replace "example" w/ specific dir

# Location of atomic data etc.
export MCDATA=$HOME/modelling/code/InputData
export OUTDATA=/$HOME/modelling/code/ExampleDir/OutputData

EXE=$HOME/modelling/code/bin/fica.exe
CODETAR=$HOME/modelling/code/bin/MCcode.tar.gz

#make sure we are in right dir
cd $OUTDATA
cd ..

if [ ! -x $EXE ]
then
  echo Executable not found!
  exit 1
fi

if [ -f $CODETAR ]
then
  cp $CODETAR .
else
  echo No archive found in $CODETAR.
fi


# copy old model back.

cd $OUTDATA

case $1  in
  '-r')
    echo reverting to previous model.
    for i in 'stst.dat' 'diagn.dat' 'yhea.dat' 'lhea.dat' 'spcp.dat' 'sica.dat' 'eica.dat' 'spct.dat' 'atmd.oud' 'ptfn.dat' 'taul.dat' 'sbib.dat'
    do
      echo $i
      cp $i.old $i
    done
    cd ..
    echo dica.dat, comp.ind
    cp dica.dat.old dica.dat
    cp comp.ind.old comp.ind
    echo Calling Gnuplot.
    gnuplot plot.gplt
    echo Finished.
    exit
    ;;
    *)
    ;;
esac


echo Saving old output...
for i in 'stst.dat' 'diagn.dat' 'yhea.dat' 'lhea.dat' 'spcp.dat' 'sica.dat' 'eica.dat' 'spct.dat' 'atmd.oud' 'ptfn.dat' 'taul.dat' 'sbib.dat'
do
  mv $i $i.old
done
echo Running new model

cd ..

#call the code
 time $EXE |tee out.out


echo Moving data to output direcotry. 
for i in 'stst.dat' 'diagn.dat' 'yhea.dat' 'lhea.dat' 'spcp.dat' 'sica.dat' 'eica.dat' 'spct.dat' 'atmd.oud' 'ptfn.dat' 'taul.dat' 'sbib.dat'
do
  mv $i $OUTDATA
done


echo Calling Gnuplot
cd $OUTDATA
cd ..
gnuplot plot.gplt

echo Finished.
