#!/bin/sh
#$ -S /bin/sh
#$ -N "LeBron FTW"
#$ -cwd -j y
#$ -V
#$ -M tmr17@duke.edu -m e
/usr/local/bin/matlab -nodesktop -nodisplay -nosplash -nojvm < PS2main.m > PS2main.log