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
% demo_planebyplane_calib, careful this script could be misleading, and it 
% is not consistent in the coordinate system with this script.
%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%% Calibration files saving order, pos2D
% The files should be saved in a order consistent with script coordinate
% system, where the x-coordinate represents the number of rows, and the
% y-the number of lines. This may pose a problem if the local coordinate
% system of the camera or the experiment does not match the local
% coordinate system of the experiment.

% functions to be used to get the interpolant for each camera
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
NZ = 1:17;% interpolation using 11 planes 
ZZ = -7:1:9;% 1mm spacing between planes. Plane number 7 is at ZZ=0; 
% The issue is that interppixel2line should be feed with lines and rows
% not the other way around.

% interpolant F1 for camera 1, size of original image 2560x1600 pixels 
% Lines(1600) comes first, and rows (2560) come second
% change resolution

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% calculate interpolants
tic;
% runs for about 3 minutes with these parameters
F1 = interppixel2line(800/4,800/4,800/16,800/16,'Cal_C1_',NZ,ZZ);
F2 = interppixel2line(800/4,800/4,800/16,800/16,'Cal_C2_',NZ,ZZ);

toc; save interpolants.mat
% clear

%% calculate interpolants (HQ)
% runs for about 20 mins
tic;
F1 = interppixel2line(800,800,800/4,800/4,'Cal_C1_',NZ,ZZ);
F2 = interppixel2line(800,800,800/4,800/4,'Cal_C2_',NZ,ZZ);
toc;

save interpolants_HQ.mat
% clear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%% demo of stereomatching using the interpolant (old)
% we stereomatch a particle situated at [0 0 0] in the real world
% it corresponds to images of the mask at location 7
% positions are in Cal_C1_7.mat for camera 1
% positions are in Cal_C2_7.mat for camera 2
% error corresponds to the plate displacement
% close all;
for planeno = 6
disp('###############################')
load(['Cal_C1_',num2str(planeno),'.mat']);pimg1=pimg;pos2D_1=pos2D;
% find index of point closest to (0,0,0) in world coordinates in image 1
d=pos2D_1(:,1).^2+pos2D_1(:,2).^2;
%[~,i1]=min(d);
i1=1;
disp('expected position camera 1');

disp([pos2D_1(i1,1) pos2D_1(i1,2) ZZ(planeno)])% should display [0 0 0]
clear T Tinv aRoi pimg pos2D*

load(['Cal_C2_',num2str(planeno),'.mat']);pimg2=pimg;pos2D_2=pos2D;
d=pos2D_2(:,1).^2+pos2D_2(:,2).^2;
%[~,i2]=min(d);
i2=1;
disp('expected position camera 2');

disp([pos2D_2(i2,1) pos2D_2(i2,2) ZZ(planeno)])% should display [0 0 0]
ref=[[pos2D_2(i2,1) pos2D_2(i2,2) ZZ(planeno)]];
clear T Tinv aRoi pimg pos2D*

% loading interpolant
load interpolants_HQ.mat
% computing ray of light for each detected point in images 1 and 2
%[o1,x1]=constrline(pimg1(i1,2),pimg1(i1,1),F1);
%[o2,x2]=constrline(pimg2(i2,2),pimg2(i2,1),F2);

% y,x pixels
% check that function

[o1,x1]=constrline(pimg1(i1,2),pimg1(i1,1),F1);
[o2,x2]=constrline(pimg2(i2,2),pimg2(i2,1),F2);

% stereomatching the 2 rays of light
[xm,d1,d2]=lines_stereomatching(o1,x1,o2,x2);

err=sqrt((xm(1)-ref(1))^2+(xm(2)-ref(2))^2+(xm(3)-ref(3))^2);
errorlist(planeno)=err;
%disp('expected position = [0 0 0] mm');
disp(['distance to lines d=' num2str(d1) ' mm']);
disp(['stereomatched position = [' num2str(xm(1)) ' ' num2str(xm(2)) ' ' num2str(xm(3))  '] mm']);
disp(['error ' num2str(err) ' mm']);
end
plot(errorlist); xlabel('plane #');ylabel('error');
%% use particle_stereomatching.m to reconstruct calibration points in 3D
clear success_rate matching_rate
mkdir stereomatch-map
mkdir stereomatch-3D
d = dir('Cal_C1_*.mat'); D = numel(d)
for planeno = 8%1:D
load(['Cal_C1_',num2str(planeno),'.mat']);posLA = pimg(1:end,:);
load(['Cal_C2_',num2str(planeno),'.mat']);posRA = pimg(1:end,:); clear pos2D*
[Coord_3D,matching_rate(planeno),varargout] = particle_stereomatching(posLA,posRA,'interpolants_HQ',sqrt(800^2+800^2),[-7 9]);
% [Coord_3D,matching_rate,varargout] = particle_stereomatching(posLA,posRA,'interpolants_HQ',[150 40],[-20 20],[-5 5]);
success_rate(planeno) = sum(round(Coord_3D(:,3))==ZZ(planeno));%/size(Coord_3D,1)*100;
zlim([-10 10]); view(-38,9); %pause()

pause(0.5);
% figure(1); print(['stereomatch-map/',sprintf('%02d',planeno),'.png'],'-dpng');
% figure(2); print(['stereomatch-3D/',sprintf('%02d',planeno),'.png'],'-dpng');
end
figure(3); plot(ZZ,success_rate,'.-'), xlim([-10 10]);
ylabel('# of particles on correct plane');
xlabel('x (mm)');
goodplot([9;7],16,10);
% print(['stereomatch-3D/success_rate.png'],'-dpng');
%%
d = dir(['*',sprintf('%06d',planeno),'*.TIF']);
LA = imread(d(1).name); RA = imread((d(2).name));
figure(2); imagesc(LA); caxis([0 2^12]); daspect([1 1 1]); colormap(gray);
hold on; scatter(posLA(:,1),size(LA,1)-posLA(:,2),'r');
figure(3); imagesc(RA); caxis([0 2^12]);daspect([1 1 1]); colormap(gray);
hold on; scatter(posRA(:,1),size(RA,1)-posRA(:,2),'b');
hold on; figure(2);scatter(posRA(:,1),size(RA,1)-posRA(:,2),'b','filled');
figure(1);