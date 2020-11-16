function [err] = gp_element_sim(E_y)
%Periodically poled LiNbO3
rng(1)
addpath('util')
ShowPlots=true;
dxf_out=false;
type = 'len';
%% All units are in um
% Constants
c       = (2.99792458e8)*(1e6);%in um/sec 
lambda  = 1.0642;

% Crystal geometry, in um. The real crystal axes are: [x,y,z] (sim)=>[y,z,x] (crystal)
L_x=12e3; % Re-calculated later, small difference
L_y=1e3;
L_z=0.5e3;
% Simulation grid
dy=0.2;
dy_thresh=50;%in microns
dy_thresh_pix=dy_thresh/dy;% in pixels
y=-L_y/2:dy:L_y/2-dy;% Tranversial coordiantes
ny=length(y);

% Applied field, birefringence and main axes
% E_y = V_app/L_y*(1e6); % in kV/m
E_DC = [0;E_y;0]; % input field in V/m, positive poling
% No Field applied
[ n ,~ ] = refraction_field( E_DC*0 , lambda );
n_o=n(2);
n_e=n(3);
beta = 2*pi*(n_e-n_o)/lambda; % Birefringence coefficient for positive domain
% Positive domain
[ n_eig_pos,theta_pos ] = refraction_field( E_DC , lambda );
n_o_pos=n_eig_pos(2);
n_e_pos=n_eig_pos(3);
% Negative domain
[ n_eig_neg,theta_neg ] = refraction_field( -E_DC , lambda );
n_o_neg=n_eig_neg(2);
n_e_neg=n_eig_neg(3);
% Note - beta * dz << 1
% Recalculate so that length of crystal is an integer amount of wavelengths
L_FW=abs(2*pi/beta); % Find the positive domain FW length
L_HW =L_FW/2; % Length of half waveplate
L_x=round(L_x/L_FW)*L_FW; % Corrected length of crystal to accomodate an integer amount of wavelengths
n_FW=round(L_x/L_FW); % Number of whole waveplates making up the length of the crystal
disp(['Maximal relative phase gained in propagation: ',num2str(abs(deg2rad(theta_pos-theta_neg)*2*n_FW))])
%% Create  mask
mask_size=[n_FW ny];
crystal_mask=~create_mask(mask_size,dy_thresh_pix,type);
% Insert gap by request of HCP
needed_gap = 20/2;
gap_dy = ceil(needed_gap/dy);
crystal_mask = int8(padarray(crystal_mask,[0,2*gap_dy],1));
gap_mask = circshift(crystal_mask, [0,gap_dy]) - circshift(crystal_mask, [0,-gap_dy]);
crystal_mask(gap_mask~=0) = -1;
% Fix wrong area in deflector
if isequal(type,'def')
    if dy == 0.5
        crystal_mask(:,1960:end) = 0; 
    elseif dy == 0.2
        crystal_mask(:,5051:end) = 0; 
    else
        crystal_mask(:,10100:end) = 0; 
    end
end
ny=size(crystal_mask,2);
L_y = ny*dy;
y=-L_y/2:dy:L_y/2-dy;% Tranversial coordiantes
x=L_FW*1e-3*(1:n_FW);
% gap_mask = padarray(crystal_mask,[2,0],1,'pre') - padarray(crystal_mask,[2,0],1,'post');
% crystal_mask = int8(padarray(crystal_mask,[2,0],1,'pre'));
% crystal_mask(gap_mask==1) = -1;

n_FW=size(crystal_mask,1);
L_x=n_FW*L_FW;
if ShowPlots
    fshow(crystal_mask,y,x);title('Crystal mask profile');xlabel('y [\mum]');ylabel('x [mm]')
end
%% Input Gaussian
Ey_0=exp(-(y/(0.25*L_y)).^2)/sqrt(2);
Ez_0=Ey_0*1i;
[Ey,Ez]=rot_2d(Ey_0,Ez_0,theta_pos);% Postive domain is default
Ey=repmat(Ey,[n_FW,1]);
Ez=repmat(Ez,[n_FW,1]);
%% Error in duty cycle
err_arr = randn(n_FW-1, 2) / 3;
err_arr(abs(err_arr)>1) = 0;
err_arr = 1+(0.01 * err_arr);
if ShowPlots
    figure;histogram(err_arr(:))
    title('Domain length multiplicative factor')
end
err_arr_jones = err_arr';
err_arr_jones = err_arr_jones(:);
err_arr_jones = [ 0; 0; err_arr_jones];
%% Propagation utilities
fy=-0.5/dy:1/L_y:0.5/dy-1/L_y;
% Shifted TF for efficiency
TF_pos_o=conj(ifftshift(exp(1i*2*pi*L_HW*sqrt((n_o_pos/lambda)^2-fy.^2))));
TF_pos_e=conj(ifftshift(exp(1i*2*pi*L_HW*sqrt((n_e_pos/lambda)^2-fy.^2))));
TF_neg_o=conj(ifftshift(exp(1i*2*pi*L_HW*sqrt((n_o_neg/lambda)^2-fy.^2))));
TF_neg_e=conj(ifftshift(exp(1i*2*pi*L_HW*sqrt((n_e_neg/lambda)^2-fy.^2))));

prop_pos_o=@(In, err) ifft((TF_pos_o.^err).*fft(In));
prop_pos_e=@(In, err) ifft((TF_pos_e.^err).*fft(In));
prop_neg_o=@(In, err) ifft((TF_neg_o.^err).*fft(In));
prop_neg_e=@(In, err) ifft((TF_neg_e.^err).*fft(In));
%% Jones matrix method
PP1=zeros([2,1].*size(crystal_mask));
PP2=PP1;
PP1(1:2:end)=(crystal_mask==0);
PP2(2:2:end)=(crystal_mask==1);
PP = double(or(PP1,PP2));
base_mat_pos = [exp(1i*2*pi*n_o_pos/lambda*L_HW),  0;
                0,                                 exp(1i*2*pi*n_e_pos/lambda*L_HW)];
base_mat_neg = [exp(1i*2*pi*n_o_neg/lambda*L_HW),  0;
                0,                                 exp(1i*2*pi*n_e_neg/lambda*L_HW)];
rot_mat_pos = [cosd(theta_pos), sind(theta_pos);
           -sind(theta_pos),cosd(theta_pos)];
rot_mat_neg = [cosd(theta_neg), sind(theta_neg);
           -sind(theta_neg),cosd(theta_neg)];       
E_curr = [Ey_0;Ez_0];
E_jon_tot = zeros([2 size(PP)]);
for cidx=1:size(PP,1)
    poling_mask=PP(cidx,:,1);%mask along x axis    
    j_mat_pos = rot_mat_pos'*(base_mat_pos.^(err_arr_jones(cidx)))*rot_mat_pos;
    j_mat_neg = rot_mat_neg'*(base_mat_neg.^(err_arr_jones(cidx)))*rot_mat_neg;
    E_next_pos = j_mat_pos*E_curr;
    E_next_neg = j_mat_neg*E_curr;
    E_next = E_next_neg;
    E_next(:,poling_mask == 1) = E_next_pos(:,poling_mask == 1);
    E_curr = E_next;
    E_jon_tot(:, cidx, :) = E_curr;
end
EL_j=(E_curr(1,:)-1i*E_curr(2,:))/sqrt(2);
EL_jon_tot = (E_jon_tot(1,:,:)-1i*E_jon_tot(2,:,:))/sqrt(2);
angle_EL_j=unwrap(angle(EL_j));
figure;plot(y,angle_EL_j)
xlabel('y [\mum]');
ylabel('Radians');
title('Jones matrix method, LCP phase at crystal output');
figure;plot(y,abs(EL_j).^2)
xlabel('y [\mum]');
ylabel('Intensity');
title('Jones matrix method, LCP intensity at crystal output');
figure;imagesc(1:1714,y,abs(squeeze(EL_jon_tot)).^2);
xlabel('y [\mum]');
ylabel('x [\mum]');
title('Jones matrix method, LCP intensity along crystal');
%% Change of medium 
Ey_j_o = medium_change(E_curr(1,:),n_o,fy,lambda,false);
Ez_j_o = medium_change(E_curr(2,:),n_e,fy,lambda,true);
%% Propagate to focal plane
Ex_j_o=padarray(Ey_j_o,[0,length(E_curr(1,:))]);
Ey_j_o=padarray(Ez_j_o,[0,length(E_curr(2,:))]);
dy_new=y(2)-y(1);
L_new=length(Ex_j_o)*dy_new;
y_new = -0.5*L_new:dy_new:0.5*L_new-dy_new;
fy_fs=-0.5/dy_new:1/L_new:0.5/dy_new-1/L_new;
TF_fs=ifftshift(exp(-1i*2*pi*real(sqrt(1/(lambda^2)-fy_fs.^2))));
prop_fs=@(In,z) ifft((TF_fs.^z).*fft(In));
%% Propagate Jones method
delta_z=300;% In microns
prop_size=600; % Cycles to run outside crystal
Ex_j_o_op0=Ex_j_o;
Ex_j_o_op=repmat(Ex_j_o_op0,[prop_size,1]);
Ey_j_o_op0=Ey_j_o;
Ey_j_o_op=repmat(Ey_j_o_op0,[prop_size,1]);
reflection_decay_mask = abs(linspace(-1,1,size(Ex_j_o_op,2)));
reflection_decay_mask = exp(-(reflection_decay_mask/0.96).^50);
for idx=1:prop_size-1
    Ex_j_o_op(idx+1,:)=reflection_decay_mask.*prop_fs(Ex_j_o_op(idx,:),delta_z);
    Ey_j_o_op(idx+1,:)=reflection_decay_mask.*prop_fs(Ey_j_o_op(idx,:),delta_z);
end
if ShowPlots
    ashow(abs((Ex_j_o_op-1i*Ey_j_o_op)/sqrt(2)).^2,y_new,delta_z*1e-3*(1:prop_size));
    title('Jones matrix method, LCP propagation outside crystal');xlabel('y [\mum]');ylabel('x [mm]');
    ashow(abs(Ex_j_o_op).^2+abs(Ey_j_o_op).^2,y_new,delta_z*1e-3*(1:prop_size));
    title('Jones matrix method, Intensity outside crystal');xlabel('y [\mum]');ylabel('x [mm]');
end
%% Whole domain propagation
% Define first domain as positive
% Define x axes as ordinary
err=zeros(size(crystal_mask,3),1);
for pidx=1:size(crystal_mask,3)
    for cidx=1:(n_FW-1)
        poling_mask=crystal_mask(cidx,:,pidx);%mask along x axis
        % Principle domain - first positive domain then negative
        % Positive domain propagation along principle axes
        Ey_pos_n=prop_pos_o(Ey(cidx,:),err_arr(cidx,1));
        Ez_pos_n=prop_pos_e(Ez(cidx,:),err_arr(cidx,1));
        % Negative domain
        [Ey_neg,Ez_neg]=rot_2d(Ey_pos_n,Ez_pos_n,theta_neg-theta_pos);% Angle rotates from positive to negative
        % Propagate through negative domain to positive domain
        Ey_neg_n=prop_neg_o(Ey_neg,err_arr(cidx,2));
        Ez_neg_n=prop_neg_e(Ez_neg,err_arr(cidx,2));    
        [Ey(cidx+1,:),Ez(cidx+1,:)]=rot_2d(Ey_neg_n,Ez_neg_n,theta_pos-theta_neg);
        % Conjugated domain - first negative then positive
        [Ey_neg,Ez_neg]=rot_2d(Ey(cidx,:),Ez(cidx,:),theta_neg-theta_pos);
        Ey_neg_n=prop_neg_o(Ey_neg,err_arr(cidx,1));
        Ez_neg_n=prop_neg_e(Ez_neg,err_arr(cidx,1));
        [Ey_pos_n,Ez_pos_n]=rot_2d(Ey_neg_n,Ez_neg_n,theta_pos-theta_neg);
        Ey_pos_nn=prop_pos_o(Ey_pos_n,err_arr(cidx,2));
        Ez_pos_nn=prop_pos_e(Ez_pos_n,err_arr(cidx,2));
        % Twice positive - gap domain between to gratings according to HCP
        Ey_t_pos_nn=prop_pos_o(prop_pos_o(Ey(cidx,:),err_arr(cidx,1)),err_arr(cidx,2));
        Ez_t_pos_nn=prop_pos_e(prop_pos_e(Ez(cidx,:),err_arr(cidx,1)),err_arr(cidx,2));
        % Seperate to different regions
        Ey(cidx+1,poling_mask==1)=Ey_pos_nn(poling_mask==1);
        Ez(cidx+1,poling_mask==1)=Ez_pos_nn(poling_mask==1);        
        Ey(cidx+1,poling_mask==-1)=Ey_t_pos_nn(poling_mask==-1);
        Ez(cidx+1,poling_mask==-1)=Ez_t_pos_nn(poling_mask==-1);
        %EL=(Ex(cidx+1,:)-1i*Ey(cidx+1,:))/sqrt(2);
    end
end
[Ey,Ez]=rot_2d(Ey,Ez,-theta_pos);
%% Post process - find phase profile
% TODO - check if another step is necessary to shift from crystal to air
EL=(Ey-1i*Ez)/sqrt(2);
if ShowPlots
    ashow(abs(EL).^2,y,x);title('LCP propagation inside crystal');xlabel('y [\mum]');ylabel('x [mm]')
    LCP_inten_along = sum(abs(EL).^2,2);
    figure;plot(x(1:size(LCP_inten_along,1)),LCP_inten_along);title('LCP intensity along crystal');xlabel('x [mm]');ylabel('Intensity')
    ylim([0,1.1*max(LCP_inten_along(:))]);
    figure;plot(y,abs(EL(end,:).^2));title('LCP intensity at crystal output');xlabel('y [\mum]');ylabel('Intensity')
    angle_EL=unwrap(angle(EL(end,:)));
    figure;plot(y,angle_EL);title('LCP angle at inside crystal, output');xlabel('y [\mum]');ylabel('Intensity')
end
%% Change of medium 
Ey_o = medium_change(Ey(end,:),n_o,fy,lambda,false);
Ez_o = medium_change(Ez(end,:),n_e,fy,lambda,true);
EL_o=(Ey_o-1i*Ez_o)/sqrt(2);
angle_EL_o=unwrap(angle(EL_o));
if ShowPlots
    figure;plot(y,abs(EL_o).^2);title('LCP intensity outside crystal');xlabel('y [\mum]');ylabel('Intensity')
    figure;plot(y,angle_EL_o);title('LCP angle outside crystal');xlabel('y [\mum]');ylabel('Radians')
end
%% Propagate to focal plane
EL_o=padarray(EL_o,[0,length(EL_o)]);
L_new=length(EL_o)*dy;
y_new = -0.5*L_new:dy_new:0.5*L_new-dy_new;
fy_fs=-0.5/dy:1/L_new:0.5/dy-1/L_new;
TF_fs=ifftshift(exp(-1i*2*pi*real(sqrt((1/lambda^2-fy_fs.^2)))));
% TF_fs=ifftshift(exp(1i*pi*lambda*fy_fs.^2));
prop_fs=@(In,z) ifft((TF_fs.^z).*fft(In));
sqr_fac=polyfit(y,angle_EL_o,2);f_dist=pi/(lambda*sqr_fac(1));
I_f=abs(prop_fs(EL_o,f_dist)).^2;
%% Propagate
delta_z=300;% In microns
prop_size=600; % Cycles to run outside crystal
EL_op0=EL_o;
EL_op=repmat(EL_op0,[prop_size,1]);
reflection_decay_mask = abs(linspace(-1,1,size(EL_op,2)));
reflection_decay_mask = exp(-(reflection_decay_mask/0.96).^50);
for idx=1:prop_size-1
    EL_op(idx+1,:)=reflection_decay_mask.*prop_fs(EL_op(idx,:),delta_z);
end
if ShowPlots
    ashow(abs(EL_op).^2,y_new,delta_z*1e-3*(1:prop_size));title('LCP propagation outside crystal');
    xlabel('y [\mum]');ylabel('x [mm]')
end
%% Output DXF file
if dxf_out
    file_name = ['D:\NON-backup\PPLN_mask_' type '_s'];
    dxf_mask_out( crystal_mask ,file_name,L_HW,L_y,L_x,dy);
end
end