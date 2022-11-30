function [Coord_3D,matching_rate,varargout] = particle_stereomatching(LA,RA,interpstr,npix,zlims)
% this function stereomatches particles from camera A and B 
% interpstr is the interpolants mat-file
% npix is the search area around the particle 
% searchbox is the searchbox around particle from camA in [x y] coordinates
% offset is the distance to offset the serachbox from the particle in camA
% in [x y] coordinates

% activate if running as funtion
load(interpstr);

% activate if running as script
% LA = posLA;
% RA = posRA;
% load interpolants_HQ.mat
% npix=400; % is the default
% searchbox = [150 40]; 
% offset = [-20 20];
% zlims = [-5 5];

%Camera L
n1=length(LA(:,1));
LA=[LA,linspace(1,n1,n1)'];
n2=length(RA(:,1));
RA=[RA,linspace(1,n2,n2)'];

Connect = [];%zeros(max([n1,n2]),2);
Coord_3D = [];%zeros(min([n1,n2]),3);
Calib=zlims;
% computing ray of light for each detected point in images 1 and 2
[o1,x1]=constrline(LA(:,2)',LA(:,1)',F1);
o1=[o1;linspace(1,n1,n1)];
x1=[x1;linspace(1,n1,n1)];
[o2,x2]=constrline(RA(:,2)',RA(:,1)',F2);
o2=[o2;linspace(1,n2,n2)];
x2=[x2;linspace(1,n2,n2)];

% for plotting
Fig1 = figure(1); clf; 
plot(LA(:,1),LA(:,2),'.r'); hold on;
plot(RA(:,1),RA(:,2),'xb'); 
daspect([1 1 1]); axis([0 800 0 800]); xlabel('X (px)'); ylabel('Y (px)');
% rectangle('Position',[LA(1,1)-searchbox(1)/2+offset(1) LA(1,2)-searchbox(2)/2+offset(2) searchbox(1:2)]);

    
%% stereomatching LA particles to RA
Connect = [];%zeros(max([n1,n2]),2);
Coord_3D = [];%zeros(min([n1,n2]),3);
ct = 1;
for k=1:size(LA,1) %number of particles in LA
    % region of 'interest' - points in RA found near the same pixel
    % location of the particle in LA
    ROI=RA(RA(:,1)>LA(k,1)-npix & RA(:,1)<LA(k,1)+npix & RA(:,2)>LA(k,2)-npix & RA(:,2)<LA(k,2)+npix,:);
%     ROI=RA(RA(:,1)>LA(k,1)-npix & RA(:,1)<LA(k,1) & RA(:,2)>LA(k,2)-npix/2 & RA(:,2)<LA(k,2)+npix/2,:);
    
    if ~isempty(ROI) % if particles are found
        posp=zeros(1,6);
        for j=1:length(ROI(:,1)) % for all particles found in the region of interest          
            % stereomatching the 2 rays of light
            [xm,d1,d2]=lines_stereomatching(o1(1:3,k),x1(1:3,k),o2(1:3,ROI(j,end)),x2(1:3,ROI(j,end)));
            posp = [posp;[xm',d1,d2,ROI(j,end)]];
        end
        % get a bunch of prospective particle matches (throw out first row)
        posp = posp(2:end,:); 
        % first, throw out matches where the z-pos lies outside the calibration range+5%
        posp=posp(bitand(posp(:,3)>(Calib(1)*1.05),posp(:,3)<(Calib(2)*1.05)),:);
        % next, choose the match that results in the closest dist. between rays of light that also meets the max distance (0.01)
        posp=posp(posp(:,4)<0.1,:); % max dist
%         posp=posp(posp(:,4)>0.001,:) % min dist?
        posp=posp(posp(:,4)==min(posp(1:end,4)),:);
        % finally, choose the final particle (if exist) and link to the corresponding particle in camL
        if ~isempty(posp)
            for ii = 1:size(posp,1)
                idx=posp(ii,6);
                Connect(ct,:)=[k,idx];
                Coord_3D(ct,:) = posp(ii,1:3);
                ct = ct+1;
            end
        end
    end
end

%% Stereomatch RA particles to LA
% find particles in RA that matches to mutliple particles in LA
% choose the one that is closest and eliminate the others

asdf = [unique(Connect(:,2)) groupcounts(Connect(:,2))];
asdf = asdf((asdf(:,2)>1),:);

for k=1:size(asdf,1)
    idx = find(Connect(:,2)==asdf(k,1));
    posp=[];
    for j=1:size(idx,1) 
        % stereomatching the 2 rays of light
        idxL = Connect(idx(j),1);
        idxR = Connect(idx(j),2);
        [xm,d1,d2]=lines_stereomatching(o2(1:3,idxR),x2(1:3,idxR),o1(1:3,idxL),x1(1:3,idxL));
        posp = [posp;[xm',d1,d2,idxL]];
    end
    idxL = posp(find(min(posp(2:end,4))),6); % the one with the closest ray
    for j=idx' % nan the rest 
        if Connect(j)~=idxL
            Connect(j,:) = [nan nan];
            Coord_3D(j,:) = [nan nan nan];
        end
    end
end
Connect = Connect(~isnan(Connect(:,1)),:);
Coord_3D = Coord_3D(~isnan(Coord_3D(:,1)),:);

%% matching rate
matching_rate = size(Coord_3D,1)/min(n1,n2);
disp(['% particles matched = ',num2str(matching_rate*100)]);

%% for plotting
for k=1:1:size(Connect,1)
    line([LA(Connect(k,1),1) RA(Connect(k,2),1)],[LA(Connect(k,1),2) RA(Connect(k,2),2)]);
end
figure(2); scatter3(Coord_3D(:,1),Coord_3D(:,2),Coord_3D(:,3),Coord_3D(:,1)*0+20,Coord_3D(:,3)*0,'filled');
 xlabel('X (mm)'); ylabel('Y (mm)'); zlabel('Z(mm)');

try varargout{1} = Connect; end
