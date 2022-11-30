%% calibration files to be used  
% Calibration files for camera 1 (13 in total)
% Cal_C1_1.mat ... Cal_C1_13.mat
% Calibration files for camera 2 (13 in total)
% Cal_C2_1.mat ... Cal_C2_13.mat

% each calibration file contains 
% pimg : positions in pix of dots detected in the image of the calibration mask
% pos2D : actual real_world coordonates of dots
% T direct transform: pos2D -> pimg obtained using cp2tform
% Tinv inverse transform: pimg -> pos2D obtained using cp2tform
% run planebyplane_calib.m to learn more on planebyplane calibration
demo_planebyplane_calib
clear

%% functions to be used to get the interpolant for each camera
% interppixel2line.m: is the main function
% pixel2line.m: subroutine
% fit3Dline.m: subroutine
% outputs have been saved in interpolants_7_16_4.mat
% this file can be loaded for stereomatching in next cell

% interpol = interppixel2line(NL,NC,NiL,NiC,calib_prefix,NZZ,ZZ)
% interpol is the interpolant used for stereomatching
% NL = number of lines per image (px)
% NC = number of columns per image (px)
% NiL = number of interpolating points in lines
% NiC = number of interpolating points in columns
% calib_prefix - Name of the calibration file '.mat' 
% containing the calibration transform T for each plane 
% NZZ - Number of the plane  
% ZZ - Position of the plane in the real world (here in mm)

% pixel2line.m
% subfunction of interppixel2line. not to be called by the user outside,this mfile
% compute equation of the ray of light for an ensemble of pixels
% for each couple of pixel coordinate (i,j), uses plane by plane calibration 
% transforms to compute corresponding real-world coordinates in each 
% calibration plane

% fit3Dline.m
% subfunction pixel2line.m not to called by the user outside,this mfile
% computes the ray of light from the real-world coordinates corresponding 

% running script
clear
NZ = 1:10;% interpolation using 7 planes 
ZZ = -4.5:1:4.5;% 1mm spacing between planes. Plane number 7 is at ZZ=0; 
% interpolant F1 for camera 1, size of original image 1088x2048 pixels 
tic;
% runs for about 3 minutes with these parameters
F1 = interppixel2line(1088/4,2048/4,272/4,512/4,'Cal_C1_',NZ,ZZ);
F2 = interppixel2line(1088/4,2048/4,272/4,512/4,'Cal_C2_',NZ,ZZ);
% parameters used to recompute F1 F2 in interpolants2019.mat
% runs for about 1h
% F1 = interppixel2line(1088,2048,272,512,'Cal_C1_',NZ,ZZ);
% F2 = interppixel2line(1088,2048,272,512,'Cal_C2_',NZ,ZZ);
toc;

save interpolants2019_new.mat
clear

%% building a ray of light using the interpolant
% F1,F2 interpolants for camera 1 and camera 2 obtained with interppixel2line.m
% a particle located (ypx_cam1,xpx_cam1) in image from camera 1
% same particle located (ypx_cam2,xpx_cam2) in image from camera 2
% [o1,x1]=constrline(ypx_cam1,xpx_cam1,F1);
% outpout: o1, x1 are two points belonging to the line for camera 1
% [o2,x2]=constrline(ypx_cam2,xpx_cam2,F2);
% outpout: o2, x2 are two points belonging to the line for camera 2

%% stereo-matching a particle using two rays of light 
% [xm,d1,d2]=lines_stereomatching(o1,x1,o2,x2);
% xm stereomatched particle position
% d1 distance from xm to the line from camera 1
% d2 distance from xm to the line from camera 2

%% demo of stereomatching using the interpolant
% we stereomatch a particle situated at [0 0 0] in the real world
% it corresponds to images of the mask at location 7
% positions are in Cal_C1_7.mat for camera 1
% positions are in Cal_C2_7.mat for camera 2

clear
load Cal_C1_7.mat;pimg1=pimg;pos2D_1=pos2D;
% find index of point closest to (0,0,0) in world coordinates in image 1
d=pos2D_1(:,1).^2+pos2D_1(:,2).^2;
[~,i1]=min(d);
%disp([pos2D_1(i1,2) pos2D_1(i1,1) 0])% should display [0 0 0]
clear T Tinv NZ ZZ aRoi pimg pos2D*

load Cal_C2_7.mat;pimg2=pimg;pos2D_2=pos2D;
d=pos2D_2(:,1).^2+pos2D_2(:,2).^2;
[~,i2]=min(d);
%disp([pos2D_2(i2,2) pos2D_2(i2,1) 0])% should display [0 0 0]
clear T Tinv NZ ZZ aRoi pimg pos2D*

% loading interpolant
load interpolants2019_new.mat
% computing ray of light for each detected point in images 1 and 2
[o1,x1]=constrline(pimg1(i1,2),pimg1(i1,1),F1);
[o2,x2]=constrline(pimg2(i2,2),pimg2(i2,1),F2);
% stereomatching the 2 rays of light
[xm,d1,d2]=lines_stereomatching(o1,x1,o2,x2);
err=sqrt(xm(1)^2+xm(2)^2+xm(3)^2);

disp('expected position = [0 0 0] mm');
disp(['distance to lines d=' num2str(d1) ' mm']);
disp(['stereomatched position = [' num2str(xm(1)) ' ' num2str(xm(2)) ' ' num2str(xm(3))  '] mm']);
disp(['error ' num2str(err) ' mm']);


