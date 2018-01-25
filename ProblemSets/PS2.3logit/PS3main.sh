#!/bin/sh
#$ -S /bin/sh
#$ -N "Integrate"
#$ -cwd -j y
#$ -V
#$ -M tmr17@duke.edu -m e
/usr/local/bin/matlab -nodesktop -nodisplay -nosplash -nojvm < PS3main.m > PS3main.log