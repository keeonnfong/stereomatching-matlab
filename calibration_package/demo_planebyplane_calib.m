% loads calibration file for plane 'plane_num' and camera 'cam_num'
%
% contains 
% pimg : positions in pix of dots detected in the image of the mask
% pos2D : actual real_world coordonates of dots
% T transform pos2D -> pimg obtained using cp2tform
% Tinv inverse transform pimg -> pos2D obtained using cp2tform
% 
% left subplot
% plots detected dots pimg in the image of the mask  
% left subplot
% actual positions pos2D of dots, and real_world coordinates of pimg
% computed using a projective Tform

clear
cam_num=1;% number of camera (1 or 2)
plane_num=11;% number of plane (between 1 and 13)
zpos=-5:1:5;
z_plane=zpos(plane_num);
disp(['z position of calibration plane: z_plane = ' num2str(z_plane) ' mm']);
fname=['Cal_C' int2str(cam_num) '_' int2str(plane_num) '.mat'];
load(fname);

% plots detected dots in the image of the mask
figure;
subplot(1,2,1)
plot(pimg(:,2),pimg(:,1),'bo');grid
xlabel('xpix (pix)')
ylabel('ypix (pix)')
title('detected dots in calibration image')
legend('detected dots')

% running projective calibration 
% computes Tform assuming projective transform
Tform = cp2tform(pimg,pos2D,'projective');
% computes real-world coordinates using pimg and Tform
[pimg_tform] = tformfwd(Tform,pimg);

% plots actual positions 
%figure;
subplot(1,2,2)
plot(pos2D(:,2),pos2D(:,1),'ks');grid
% adds real_world coordinates pimg_tform
hold on;
plot(pimg_tform(:,2),pimg_tform(:,1),'r.');
hold off
xlabel('X (mm)')
ylabel('Y (mm)')
title('real world coordinates')
legend('actual positions','pimg tform')




