lambda = 1.064;
N = 500;
E_max = 1e6;
E_vec = linspace(0,E_max,N);
n_vec = zeros(3,N);
for idx=1:N
    [ n_eig,theta ] = refraction_field( [0;0;E_vec(idx)] , lambda );
    n_vec(:, idx) = n_eig;
end
figure;
plot(E_vec,n_vec(2,:)-n_vec(2,1));
title('n2(E)-n2(0)');
xlabel('E');
figure;
plot(E_vec,n_vec(3,:)-n_vec(3,1));
title('n3(E)-n3(0)');
xlabel('E');
delta_phi_regular = 2*pi/lambda * (n_vec(2,end)-n_vec(2,1));
[ ~,theta ] = refraction_field( [0;E_vec(end);0] , lambda );
beta = 2*pi*(n_vec(3,1)-n_vec(2,1))/lambda; % Birefringence coefficient for positive domain
L_FW=abs(2*pi/beta); % Find the positive domain FW length
delta_phi = 4*deg2rad(theta) / L_FW;