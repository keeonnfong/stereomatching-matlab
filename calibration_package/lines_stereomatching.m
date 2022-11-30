function [xm,d1,d2]=lines_stereomatching(o1,x1,o2,x2)

% [xm]=lines_stereomatching(o1,x1,o2,x2);
%
% this program computes the perpendicular line to two non parallel, non
% intersecting lines
% xm coordinates of the mid point
% D1 o1+t1*u1
% D2 o2+t2*u2
%
% DP line joining h1, h2 crossing points between DP,D1 and DP,D2
% xm midle of [h1 h2]
% input : coordinate of o1, x1, o2, x2 in world coordinate Rw
% output : xm=(h1+h2)/2;
% simple case
%
% o1=[-1;0;1];
% x1=[0;0;1];
% o2=[0;1;-1];
% x2=[0;0;-1];
% xm_th=[0;0;0];
%
% [xm,d1,d2]=lines_stereomatching(o1,x1,o2,x2);
%
% R.V. 30/20/2010

if size(o1,1)==3&size(o1,2)==1

% unit vector of D1
u=x1-o1;
u=u/norm(u);
% unit vector of D2
v=x2-o2;
v=v/norm(v);
% unit vector of DP
w=cross(u,v);
w=w/norm(w);

% we have h1 \in D1 et h2 \in D2
%
% h1=o1+lambda(1)*u;
% h2=o2+lambda(2)*v;
% vec(h1h2)=h2-h1=lambda(3)*w;
%
% we compute lambda(1),lambda(2),lambda(3)
% h2-h1=o2-o1+lambda(2)*v-lambda(1)*u=lambda(3)*w;
%
% then
% xm=(h1+h2)/2;
% xm=(o1+o2)/2+lambda(1)*u+lambda(2)*v;

M=zeros(3,3);
M(:,1)=u;
M(:,2)=-v;
M(:,3)=w;
invM=M^-1;
lambda=invM*(o2-o1);

xm=o1+o2+lambda(1)*u+lambda(2)*v;
xm=xm/2;
d1=norm(cross(xm-o1,u));
d2=norm(cross(xm-o2,v));

else
    xm=[];
    d1=[];
    d2=[];
    disp('vectors must be of size size(o1)=[3,1]');
end