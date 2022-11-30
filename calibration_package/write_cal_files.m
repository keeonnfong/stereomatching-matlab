clear all
% close all

n=11; %number of planes

cam={'L','R'};
tifdir = dir('*.TIF');
for i=1:n
    for j=1:2
        
        name=strcat('tempGrid',num2str(i,'%06.f'),'.',cam{j},'.CPT');
        tifname=strcat('trimmed_CalPTV-DAN-WT',num2str(i,'%06d'),...
            '.T000.D000.P000.H000.',cam{j},'A.TIF');
        calImg = imread(tifname);
%         name=strcat('Targ_pts_',num2str(i,'%06.f'),'.',cam{j},'.mat');
        name_cal=strcat('Cal_C',int2str(j),'_',int2str(i),'.mat');
        bla=dlmread(name,',',9,0);
%         bla=load(name);
%         bla=bla.targ_pts;
        pimg=bla(:,1:2);%pimg(:,2)=abs(1600-pimg(:,2));
        pos2D=bla(:,3:4);%pos2D=fliplr(pos2D);
%         Tinv = cp2tform(pimg,pos2D,'polynomial',3);
%         T = cp2tform(pos2D,pimg,'polynomial',3);
        Tinv = cp2tform(pimg,pos2D,'projective');
        T = cp2tform(pos2D,pimg,'projective');

        aRoi=[min(pimg(:,1)),min(pimg(:,2));...
            max(pimg(:,1)),min(pimg(:,2));...
            max(pimg(:,1)),max(pimg(:,2));...
            min(pimg(:,1)),max(pimg(:,2))];
%         
%         figure; plot(pos2D(:,1),pos2D(:,2),'.');
%         
        if mod(j,2)==1
            subplot(121); cla;
            imagesc(flipud(calImg)); hold on;
            plot(pimg(:,1),(pimg(:,2)),'b.');
            hold on; plot(pimg(3,1),(pimg(3,2)),'r.');
            hold on; plot(pimg(2,1),(pimg(2,2)),'g.');
        else
            subplot(122); cla;
            imagesc(flipud(calImg)); hold on;
            plot(pimg(:,1),(pimg(:,2)),'b.');
            hold on; plot(pimg(3,1),(pimg(3,2)),'r.');
            hold on; plot(pimg(2,1),(pimg(2,2)),'g.');
        end
        axis tight; daspect([1 1 1]);
       
        save(name_cal,'aRoi', 'pimg', 'pos2D', 'T', 'Tinv');
    end
    if i==1
        pause();
    else
        pause(1);
    end
end