#!/bin/sh
#$ -S /bin/sh
#$ -N "job name"
#$ -cwd -j y
#$ -V
#$ -M mynetID@duke.edu -m e
/usr/local/bin/matlab -nodesktop -nodisplay -nosplash -nojvm < myfile.m > myfile.log