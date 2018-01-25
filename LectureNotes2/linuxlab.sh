cd 
pwd
cd Teaching
cd ..
cd 
ls
mkdir Matlab_Module
cd Matlab_Module
mkdir Linux_Lab
cd Linux_Lab
touch linuxlab.m
mv linuxlab.m linuxlab1.m
cp -f linuxlab1.m linuxlab2.m
cd ..
cp -r Linux_Lab Linux_Lab1
ls -l
cd Linux_Lab1
rm linuxlab2.m
ls -a
cd ..
rm -r Linux_Lab1
ls
w
uname -a
whoami
finger tmr17
cd ~/Teaching/PS1.3
matsub ps3main.m ps3main.log
qsub -q all.q ps3main.sh
qstat
qstat -u "*"
qstat -f
qstat -u "*" -f
