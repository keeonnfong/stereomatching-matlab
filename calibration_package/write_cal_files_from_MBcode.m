clear all
% close all

n=11; %number of planes

cam={'1','2'};
tifdir = dir('*.TIF');
ct = 1;
for i=5:15
    for j=1:2
        
        name=strcat('calib2D_',num2str(i),'_cam',cam{j},'.mat');
        tifname=strcat('CalibrationPlan_',num2str(i,'%02d'),...
            '_cam',cam{j},'.tif');
        calImg = imread(tifname);
        name_cal=strcat('Cal_C',int2str(j),'_',int2str(ct),'.mat');
        load(name,'pimg','pos3D');
        pos2D=pos3D(:,1:2); clear pos3D;
        Tinv = cp2tform(pimg,pos2D,'projective');
        T = cp2tform(pos2D,pimg,'projective');

        aRoi=[min(pimg(:,1)),min(pimg(:,2));...
            max(pimg(:,1)),min(pimg(:,2));...
            max(pimg(:,1)),max(pimg(:,2));...
            min(pimg(:,1)),max(pimg(:,2))];
         
        if mod(j,2)==1
            subplot(121); cla;
            imagesc((calImg)); hold on;
            plot(pimg(:,1),(pimg(:,2)),'bo');
            hold on; plot(pimg(3,1),(pimg(3,2)),'ro');
            hold on; plot(pimg(2,1),(pimg(2,2)),'go');
        else
            subplot(122); cla;
            imagesc((calImg)); hold on;
            plot(pimg(:,1),(pimg(:,2)),'bo');
            hold on; plot(pimg(3,1),(pimg(3,2)),'ro');
            hold on; plot(pimg(2,1),(pimg(2,2)),'go');
        end
        axis tight; daspect([1 1 1]);
       
        save(name_cal,'aRoi', 'pimg', 'pos2D', 'T', 'Tinv');
    end
%     if i==1
        ct = ct+1; pause();
%     else
%         pause(1);
%     end
end