#!/bin/bash

if [ $1 == "-r" ] 
then 
  cp -R --preserve $2/* .
else 
  mkdir ./$1   
  cp -R --preserve OutputData ./$1/OutputData
  cp -R --preserve *.* ./$1
fi
