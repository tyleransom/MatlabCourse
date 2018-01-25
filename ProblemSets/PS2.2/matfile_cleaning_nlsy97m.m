ID = data(:,1);
male = data(:,2);
AFQT = data(:,3);
Mhgc = data(:,4);
hgc = data(:,5:9);
hgc(600,3)=18;
hgc(600,4)=18;
hgc(600,5)=18;
hgc(5928,3)=20;
hgc(5928,4)=20;
hgc(5928,5)=20;
hgc(2316,3)=20;
hgc(2316,4)=20;
hgc(2316,5)=20;
hgc(78,5)=20;
hgc(2527,5)=20;
hgc(1233,:)=9*ones(1,5);
hgc([595 2584 5136],3:5) = 15*ones(3);
exper = data(:,10:14);
Diploma = data(:,15:19);
AA = data(:,20:24);
BA = data(:,25:29);
log_wage = data(:,30:34);
activity = data(:,35:39);
clear data
clear colheaders
clear textdata
save nlsy97m
