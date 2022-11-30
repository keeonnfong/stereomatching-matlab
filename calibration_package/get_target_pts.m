clear all
close all

cam={'L','R'};
n=10;


for j=1:2
    
name=strcat('dewarp_',cam{j},'.Cpt2');
fich=fopen(name);



data=textscan(fich,'%s%f%s%f%s%f%s%f%s%f%s%f%s', 'HeaderLines', 10,'Delimiter', '"','TreatAsEmpty', '~');

targ_pts_all=[cell2mat(data(:,2)),cell2mat(data(:,4)),cell2mat(data(:,6)),cell2mat(data(:,8)),cell2mat(data(:,10))];

for i=1:n
    
    name=strcat('Targ_pts_',num2str(i,'%06.f'),'.',cam{j},'.mat');
    targ_pts=targ_pts_all(targ_pts_all(:,end)==i,:);
    
    save(name,'targ_pts')
    
end
    
end