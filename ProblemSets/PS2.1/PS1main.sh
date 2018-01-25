#!/bin/sh
#$ -S /bin/sh
#$ -N "matsub"
#$ -cwd -j y
#$ -V
# -M tmr17@duke.edu -m e
/usr/local/bin/matlab -nodesktop -nodisplay -nosplash -nojvm < PS1main.m > PS1main.log
TMP=`mktemp -t ml`
trap "rm $TMP* 2>/dev/null" 0
echo 'done' > ml
mutt -s 'happy news' -a PS1main.log tyleransom@gmail.com < ml
exit
