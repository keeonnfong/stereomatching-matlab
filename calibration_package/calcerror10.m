% For different number of planes
% Load the interpolants F1=camera1, F2=camera2;
clear all
% close all
%tic;
%for N = 0:5
N=0;

%eval (['load interpolant/interpolant' num2str(N)]);
load 'interpolants_HQ.mat'
%
rrr=-5:1:5;
zr=rrr;

for i=1%:11
    %     ind=find(zr(i)==rrr);
    ind=i;
    % Camera 1
    %clear all
    %i=1;
    fnc1 = ['Cal_C1_' num2str(ind) '.mat'];
    load (fnc1,'Tinv','pos2D','pimg');
    Tinv1 = Tinv; clear Tinv
    pos2D1 = pos2D; clear pos2D
    pimg1 = pimg; clear pimg
    
    % Camera 2
    fnc2 = ['Cal_C2_' num2str(ind) '.mat'];
    load (fnc2,'Tinv','pos2D','pimg');
    Tinv2 = Tinv; clear Tinv
    pos2D2 = pos2D; clear pos2D
    pimg_ref = pimg; clear pimg
    
    
    [pos2D,ind1,ind2]=intersect(pos2D1,pos2D2,'rows');
    pimg1= pimg1(ind1,:);
    pimg_ref= pimg_ref(ind2,:);
    pos2D1=pos2D1(ind1,:);
    pos2D2=pos2D2(ind2,:);
    % Transform original centroids of camera1 to original images of camera2
    %im2= imread(['Cal_14_03_' num2str(i) '_Camera2_1.bmp']);
%    im2=imread('CalPTV-DAN-WT000001.T000.D000.P000.H000.RA.TIF');
    [pimg2(:,1), pimg2(:,2)] = tforminv(Tinv2, pos2D1(:,1), pos2D1(:,2));
    
%     figure; plot(pimg_ref(:,1),pimg_ref(:,2),'b.');
%     hold on; plot(pimg2(:,1),pimg2(:,2),'r.');
    
%     figure
%     imagesc(im2);colormap(gray); hold on;
%     plot(pimg2(:,1),pimg2(:,2),'+r');
    
    Ypxc1 = pimg1(:,2);
    Xpxc1 = pimg1(:,1);
    %
    Ypxc2 = pimg2(:,2);
    Xpxc2 = pimg2(:,1);
    
    [ o1, x1 ] = constrline( Ypxc1', Xpxc1', F1 );
    [ o2, x2 ] = constrline( Ypxc2', Xpxc2', F2 );
    
    for k=1:length(o1)
        [xm(:,k),d1(k),d2(k)] = lines_stereomatching(o1(:,k),x1(:,k),o2(:,k),x2(:,k));
    end
    
%        figure; plot3(xm(1,:),xm(2,:),xm(3,:),'.b'); daspect([1 1 1]);
    
    nan = find(isnan(xm(3,:))>0);
    xm(:,nan)=[];
    d1(nan)=[];
    d2=[];
    pos2D1(nan,:)=[];
    
 
    
    plan(i).x = xm(1,:);
    %plan(i).x = xm(1,xm(1,:)<50 & xm(1,:)>-50 & xm(2,:)<30 & xm(2,:)>-30);
    plan(i).y = xm(2,:);
    %plan(i).y = xm(2,xm(1,:)<50 & xm(1,:)>-50 & xm(2,:)<30 & xm(2,:)>-30);
    plan(i).z = xm(3,:);
    %plan(i).z = xm(3,xm(1,:)<50 & xm(1,:)>-50 & xm(2,:)<30 & xm(2,:)>-30);
    
    plan(i).xr=pos2D1(:,1)';
    %plan(i).xr=pos2D1(xm(1,:)<50 & xm(1,:)>-50 & xm(2,:)<30 & xm(2,:)>-30,1)';
    plan(i).yr=pos2D1(:,2)';
    %plan(i).yr=pos2D1(xm(1,:)<50 & xm(1,:)>-50 & xm(2,:)<30 & xm(2,:)>-30,2)';
    plan(i).zr=zr(i)*ones(1,size(pos2D1,1));
    %plan(i).zr=zr(i)*ones(1,size(pos2D1(xm(1,:)<50 & xm(1,:)>-50 & xm(2,:)<30 & xm(2,:)>-30,2),1));
    
    Nmax=13;
    cc=bone(Nmax+3);
    
    %plot3(plan(i).x, plan(i).z,plan(i).y,'.','markersize',10)
    %plot3(plan(i).x, plan(i).z,plan(i).y,'.','markersize',10,'color',cc(i,:))
    %hold on
    
    %
    %medZ(i) = mean(xm(3,:));
    %stdZ(i) = std(xm(3,:));
    %
    clear fnc1 Tinv1 pos2D1 pimg1 fnc2 Tinv2 pimg2 Ypxc1 Xpxc1 Ypxc2 Xpxc2 o1 o2 x1 x2 xm d1 d2 i k nan
end
%
X=[plan.x];
Y=[plan.y];
Z=[plan.z];
Xr=[plan.xr];
Yr=[plan.yr];
Zr=[plan.zr];

% Small volume
%X=[plan(6:8).x];
%Y=[plan(6:8).y];
%Z=[plan(6:8).z];
%Xr=[plan(6:8).xr];
%Yr=[plan(6:8).yr];
%Zr=[plan(6:8).zr];

dx=X-Xr;
dy=Y-Yr;
dz=Z-Zr;
d=sqrt(dx.^2+dy.^2+dz.^2);

[PDFx xx]=hist(dx,640);
[PDFy yy]=hist(dy,640);
[PDFz zz]=hist(dz,640);
[PDFd dd]=hist(d,640);

PDFx=PDFx/sum(PDFx)/(xx(2)-xx(1));
PDFy=PDFy/sum(PDFy)/(yy(2)-yy(1));
PDFz=PDFz/sum(PDFz)/(zz(2)-zz(1));
PDFd=PDFd/sum(PDFd)/(dd(2)-dd(1));

DX (N+1) = mean(abs(dx));
DY (N+1) = mean(abs(dy));
DZ (N+1) = mean(abs(dz));
D (N+1) = mean(abs(d));

%end

% PDF2D
%[H XX YY]=hist2d([d ; Z],16,12,[0 0.2]);
%figure;pcolor(XX,YY,H)
%toc;
%%
% Nmax=12;
% cc=winter(Nmax);
% 
% NN = 3:2:13
% plot(NN,fliplr(DX),'+-','color',cc(7,:))
% hold on
% plot(NN,fliplr(DY),'+-','color',cc(11,:))
% plot(NN,fliplr(DZ),'+-','color',cc(4,:))
% plot(NN,fliplr(D),'+-','color',cc(1,:))
% 

%%
% Contourf of the error
[xq,yq] = meshgrid(-55:1:55, -47:1:42);
vqn=zeros(size(xq,1),size(xq,2));

for i=1:10
    X=[plan(i).x];
    Y=[plan(i).y];
    Z=[plan(i).z];
    Xr=[plan(i).xr];
    Yr=[plan(i).yr];
    Zr=[plan(i).zr];
    
    dx=X-Xr;
    dy=Y-Yr;
    dz=Z-Zr;
    d=sqrt((dx.^2+dy.^2+dz.^2)/3);
    
    vq = griddata( plan(i).xr , plan(i).yr , d ,xq,yq);
    
    vqn = vqn + vq;
end

figure
contourf(xq,yq,vqn/10);
axis equal tight
h=colorbar;
h.Label.String = '\langle [(dx^2 + dy^2 + dz^2)/3]^{1/2}\rangle(mm)';


for i=1:10
    X=[plan(i).x];
    Y=[plan(i).y];
    Z=[plan(i).z];
    Xr=[plan(i).xr];
    Yr=[plan(i).yr];
    Zr=[plan(i).zr];
    
    dx=X-Xr;
    dy=Y-Yr;
    dz=Z-Zr;
    d=sqrt(dz.^2);
%     d=dz;
    
    vq = griddata( plan(i).xr , plan(i).yr , d ,xq,yq);
    
    vqn = vqn + vq;
end

figure
contourf(xq,yq,vqn/10);
axis equal tight
h=colorbar;
h.Label.String = '\langle [dz]\rangle(mm)';

