% Implication of the intensity variation over the atomic cloud. Calculated
% taking into account the variation of Omega_0 and s over the atomic cloud

clear all
close all

%% Parameters for 88Sr

c = 3e8;
kB = 1.38e-23;
m = 87.9*1.66e-27;          %atomic mass
lambda = 461e-9;            %Blue line
k = 2*pi/lambda;            %Wavevector
G_Sr = 2*pi*30.5e6;         %Blue line
I_sat = 40.5;               %Saturation intensity for the blue line, in mW/cm^2

%% Experimental parameters

%Experimental setup
h = 0.0042;                 %distance between the atomic cloud and the virtual mirror
theta_0 = 4.3*pi/180;       %incidence angle of the incident beam
theta_f = pi/(k*h*theta_0); %fringe period
N_point_theta = 300;        %to calculate the intensity as a function of theta and estimate the contrast
theta = linspace(theta_0-3*theta_f,theta_0+3*theta_f,N_point_theta);

L = 0.6;                    % distance between atoms and real mirror
tau_c = 2*L/c;              % 2 x travel time between atoms and real mirror

%Atomic cloud
sigma_x = 0.5e-3;           %sizes of the atomic cloud from the red MOT
sigma_y = 0.5e-3;
sigma_z = 0.5e-3;
b0 = 0.6;                   %On-resonance OD

%Probe beam
w0 = 2.1e-3;               %waist of the probe beam
z0 = pi*w0^2/lambda;
d_l = 0;%1.08e-3/2            %distance between the center of the cloud and the center of the laser incoming beam

%N_points_Omega = 30;
%Omega_0 = linspace(0.3,1.7,N_points_Omega)/tau_c;     % Rabi frequency at the center of the beam, for just one beam
% Omega_0 = linspace(1,2,N_points_Omega)/tau_c;     % Rabi frequency at the center of the beam, for just one beam
%s_0 = 8*Omega_0.^2/G_Sr^2;                          % Saturation parameter at maximum of standing wave, created by both incoming and reflected beam
%Power = I_sat*s_0./4*pi*(w0*100)^2/2;               % Total power in the incident beam, in mW

s_in = [0.4, 0.8, 1.2, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0];  % saturation parameter of incoming beam at the center of the cloud
s_0 = 4*s_in;
Omega_0 = sqrt(s_in/2)*G_Sr;
Power = I_sat*s_in*pi*(w0*100)^2/2;  
N_points_Omega = length(s_in);

% %% Section 1 - Reduction of the contrast due to the amplitude variation of the probe beam (Gaussian beam)
% % Calculation with an incident laser beam centered on the atomic cloud and
% % as if the amplitude was not tilted
% % Calculations done using the equation of I as a function of s, Omega_0 and Omega_M
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% r_max = 3*min(w0, max([sigma_x,sigma_y]));  % Integration limit
% z_max = 3*min(w0, sigma_z);  % Integration limit
% N_points_r = 50;                                    % Number of points for r to perform the integral
% r = linspace(0,r_max,N_points_r);
% N_points_z = 100;                                   % Number of points for z to perform the integral
% z = linspace(-z_max-h,z_max-h,N_points_z);
% N_points_zsw = 20;                                   % Number of points for z to perform the integral just on a period of the standing wave
% zsw = linspace(0,lambda*(N_points_zsw-1)/N_points_zsw,N_points_zsw);
% [R,Z,Zsw] = meshgrid(r,z,zsw);
% 
% rho = exp(-(R.^2/2/sigma_x^2 + (Z+h).^2/2/sigma_z^2));      % Atomic density profile
% R_sum = sum(R,3);
% rho_sum = sum(rho,3);
% 
% Intensities = zeros (N_points_Omega,N_point_theta);
% 
% for j = 1:N_points_Omega    
%     Omega = 2*Omega_0(j)*cos(k*(Z+Zsw)*cos(theta_0)).*exp(-R.^2/w0^2);
%     Omega_M = sqrt(Omega.^2-G_Sr^2/16);
%     s = 2*(abs(Omega)).^2/G_Sr^2;
%     
%     ss = s./(4*(1+s)).*(2./(1+s) + exp(-G_Sr*tau_c/2) + (s-1)./(s+1).*cos(Omega_M*tau_c)*exp(-3*G_Sr*tau_c/4) + G_Sr./(4*Omega_M).*(5*s-1)./(s+1).*sin(Omega_M*tau_c)*exp(-3*G_Sr*tau_c/4));
%     for i = 1:N_point_theta
%         %I(i) = sum(sum(R.*rho.*(s./(1+s) + 2.*cos(2*k*Z*cos(theta(i))).*ss)));  % Intensity as a function of theta (fringes)
%         %I(i) = sum(sum(R.*rho.*(s./(1+s) + 2*real(exp(-sqrt(-1)*2*k*Z*cos(theta(i))).*ss))));  % Intensity as a function of theta (fringes);
%         ss_avg = sum(s./(1+s) + 2*real(exp(-sqrt(-1)*2*k*(Z+Zsw)*cos(theta(i))).*ss),3);
%         I(i) = sum(sum(R_sum.*rho_sum.*ss_avg));  % Intensity as a function of theta (fringes);
%         
%     end
% %     figure(20)
% %     hold on
% %     plot(theta-theta_0,I/mean(I))
%     Intensities(j,:) = I;
%     Contraste_GaussianBeam(j) = (max(I)-min(I))/(max(I)+min(I))*2;     % Calculation of the fringe contrast
% end
% 
% fid = fopen( 'Fringes_centered_beam.txt', 'w' );
% for n=1:N_points_Omega
%     for m=1:N_point_theta
%         fprintf( fid, '%.1f  ', Intensities(n,m));
%     end
%     fprintf(fid,'\r\n');
% end
% fclose(fid);
% 
% figure (3)
% hold on
% plot(s_0/4,Contraste_GaussianBeam,'b');
% xlabel('\fontsize{16} s')
% ylabel('\fontsize{16} C')

%% Section 2 - Reduction of the contrast due to the amplitude variation of the probe beam (Gaussian beam) + attenuation
% Calculation with an incident laser beam centered on the atomic cloud and
% as if the amplitude was not tilted
% Calculations done using the equation of I as a function of s, Omega_0 and Omega_M
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r_max = 3*min(w0, max([sigma_x,sigma_y]));  % Integration limit
z_max = 3*min(w0, sigma_z);  % Integration limit
N_points_r = 50;                                    % Number of points for r to perform the integral
r = linspace(0,r_max,N_points_r);
N_points_z = 100;                                   % Number of points for z to perform the integral
z = linspace(-z_max-h,z_max-h,N_points_z);
dz = z(2)-z(1);
N_points_zsw = 20;                                   % Number of points for z to perform the integral just on a period of the standing wave
zsw = linspace(0,lambda*(N_points_zsw-1)/N_points_zsw,N_points_zsw);
[R,Z,Zsw] = meshgrid(r,z,zsw);

rho = exp(-(R.^2/2/sigma_x^2 + (Z+h).^2/2/sigma_z^2));      % Atomic density profile
R_sum = sum(R,3);
rho_sum = sum(rho,3);

rho_sum_BL = sum(rho, 1);  % Sum along dimension 2 → (Nx, 1, Nz)

disp(size(rho))
disp(size(rho_sum_BL))
% rho_perm = permute(rho,[2 1 3]);
% rho_sum_BL_perm = permute(rho_sum_BL,[3 2 1]);
% disp(size(rho_perm))
% disp(size(rho_sum_BL_perm))

rho_center = zeros(size(rho));
disp(size(rho_center))
for ii = 1: N_points_z
    rho_center(ii,:,:) = rho(ii,:,:)./rho_sum_BL(1,:,:)/dz*b0;  %Normalized for Beer-Lambert lawZ
end


disp(size(rho_center))

Intensities = zeros (N_points_Omega,N_point_theta);

for j = 1:N_points_Omega    
    %Omega versus Z
    s_0_Z = zeros(size(Z));
    s_0_Z(1,:,:) = 2*(Omega_0(j)/G_Sr)^2.*exp(-2*R(1,:,:).^2/w0^2);

    %Saturation variation for the incoming beam
    for i=1:(length(z)-1)
        ds_0 = -rho_center(i,:,:)./(1+s_0_Z(i,:,:)).*s_0_Z(i,:,:)*dz;
        s_0_Z(i+1,:,:) = s_0_Z(i,:,:)+ds_0;
    end
    
    %Saturation variation for the reflected beam
    Z_reflected = flip(Z, 2);
    s_0_Z_reflected = zeros(size(Z));
    s_0_Z_reflected(1,:,:) = s_0_Z(end,:,:);

    rho_center_reflected = flip(rho_center,2);

    for i=1:(length(z)-1)
        ds_0_reflected = -rho_center_reflected(i,:,:)./(1+s_0_Z_reflected(i,:,:)).*s_0_Z_reflected(i,:,:)*dz;
        s_0_Z_reflected(i+1,:,:) = s_0_Z_reflected(i,:,:)+ds_0_reflected;
    end
    
%     figure(100)
%     plot(s_0_Z(:,1,1))
%     hold on
%     plot(s_0_Z(:,50,1))
    
    %Interference of incoming and reflected beam
    Omega_Z = sqrt(s_0_Z/2)*G_Sr;
    Omega_Z_reflected = flip(sqrt(s_0_Z_reflected/2)*G_Sr,2);

    Omega_total = abs(Omega_Z.*exp(sqrt(-1)*k*(Z+Zsw)*cos(theta_0)) + Omega_Z_reflected.*exp(-sqrt(-1)*k*(Z+Zsw)*cos(theta_0)));
    
    %Next
    Omega = Omega_total;    %2*Omega_0(j)*cos(k*(Z+Zsw)*cos(theta_0)).*exp(-R.^2/w0^2);
    Omega_M = sqrt(Omega.^2-G_Sr^2/16);
    s = 2*(abs(Omega)).^2/G_Sr^2;
    
    ss = s./(4*(1+s)).*(2./(1+s) + exp(-G_Sr*tau_c/2) + (s-1)./(s+1).*cos(Omega_M*tau_c)*exp(-3*G_Sr*tau_c/4) + G_Sr./(4*Omega_M).*(5*s-1)./(s+1).*sin(Omega_M*tau_c)*exp(-3*G_Sr*tau_c/4));
    for i = 1:N_point_theta
        %I(i) = sum(sum(R.*rho.*(s./(1+s) + 2.*cos(2*k*Z*cos(theta(i))).*ss)));  % Intensity as a function of theta (fringes)
        %I(i) = sum(sum(R.*rho.*(s./(1+s) + 2*real(exp(-sqrt(-1)*2*k*Z*cos(theta(i))).*ss))));  % Intensity as a function of theta (fringes);
        ss_avg = sum(s./(1+s) + 2*real(exp(-sqrt(-1)*2*k*(Z+Zsw)*cos(theta(i))).*ss),3);
        I(i) = sum(sum(R_sum.*rho_sum.*ss_avg));  % Intensity as a function of theta (fringes);
        
    end
%     figure(20)
%     hold on
%     plot(theta-theta_0,I/mean(I))
%     Intensities(j,:) = I;
    Contraste_GaussianBeam_Attenuation(j) = (max(I)-min(I))/(max(I)+min(I))*2;     % Calculation of the fringe contrast
end

fid = fopen( 'Fringes_centered_beam.txt', 'w' );
for n=1:N_points_Omega
    for m=1:N_point_theta
        fprintf( fid, '%.1f  ', Intensities(n,m));
    end
    fprintf(fid,'\r\n');
end
fclose(fid);

figure (3)
hold on
plot(s_0/4,Contraste_GaussianBeam_Attenuation,'g');
xlabel('\fontsize{16} s')
ylabel('\fontsize{16} C')

fid = fopen( 'Contrast_Mathilde.txt', 'w' );
for n=1:N_points_Omega
    fprintf(fid, '%.3f  ', s_0(n)/4);
    fprintf(fid, '%.3f  ', Contraste_GaussianBeam_Attenuation(n));
    fprintf(fid,'\r\n');
end
fclose(fid);

% %% Section 3 - Plane waves + attenuation
% % Calculations done using the equation of I as a function of s, Omega_0 and Omega_M
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% r_max = 3*min(w0, max([sigma_x,sigma_y]));  % Integration limit
% z_max = 3*min(w0, sigma_z);  % Integration limit
% N_points_r = 50;                                    % Number of points for r to perform the integral
% r = linspace(0,r_max,N_points_r);
% N_points_z = 100;                                   % Number of points for z to perform the integral
% z = linspace(-z_max-h,z_max-h,N_points_z);
% dz = z(2)-z(1);
% N_points_zsw = 20;                                   % Number of points for z to perform the integral just on a period of the standing wave
% zsw = linspace(0,lambda*(N_points_zsw-1)/N_points_zsw,N_points_zsw);
% [R,Z,Zsw] = meshgrid(r,z,zsw);
% 
% rho = exp(-(R.^2 + (Z+h).^2)/2/sigma_x^2);      % Atomic density profile
% R_sum = sum(R,3);
% rho_sum = sum(rho,3);
% 
% rho_center = exp(-((z+h).^2)/2/sigma_z^2); %for R = 0
% rho_center = rho_center/sum(rho_center)/dz*b0;  %Normalized for Beer-Lambert law
% 
% Intensities = zeros(N_points_Omega,N_point_theta);
% 
% for j = 1:N_points_Omega 
%     %Omega versus Z
%     s_0_z = zeros(1,length(z));
%     s_0_z(1) = 2*(Omega_0(j)/G_Sr)^2;
% 
%     %Saturation variation for the incoming beam
%     for i=1:(length(z)-1)
%         ds_0 = -rho_center(i)/(1+s_0_z(i))*s_0_z(i)*dz;
%         s_0_z(i+1) = s_0_z(i)+ds_0;
%     end
%     
%     %Saturation variation for the reflected beam
%     z_reflected = flip(z);
%     s_0_z_reflected = zeros(1,length(z));
%     s_0_z_reflected(1) = s_0_z(end);
%     rho_center_reflected = flip(rho_center);
% 
%     for i=1:(length(z)-1)
%         ds_0_reflected = -rho_center_reflected(i)/(1+s_0_z_reflected(i))*s_0_z_reflected(i)*dz;
%         s_0_z_reflected(i+1) = s_0_z_reflected(i)+ds_0_reflected;
%     end
%     
%     %Interference of incoming and reflected beam
%     Omega_z = sqrt(s_0_z/2)*G_Sr;
%     [~, Omega_Z, ~] = meshgrid(r, Omega_z, zsw);
%     Omega_z_reflected = flip(sqrt(s_0_z_reflected/2)*G_Sr);
%     [~, Omega_Z_reflected, ~] = meshgrid(r, Omega_z_reflected, zsw);
% 
%     Omega_total = abs(Omega_Z.*exp(sqrt(-1)*k*(Z+Zsw)*cos(theta_0)) + Omega_Z_reflected.*exp(-sqrt(-1)*k*(Z+Zsw)*cos(theta_0)));
% 
%     %Next
%     Omega = Omega_total;    %2*Omega_0(j)*cos(k*(Z+Zsw)*cos(theta_0));
%     Omega_M = sqrt(Omega.^2-G_Sr^2/16);
%     s = 2*(abs(Omega)).^2/G_Sr^2;
%     
%     ss = s./(4*(1+s)).*(2./(1+s) + exp(-G_Sr*tau_c/2) + (s-1)./(s+1).*cos(Omega_M*tau_c)*exp(-3*G_Sr*tau_c/4) + G_Sr./(4*Omega_M).*(5*s-1)./(s+1).*sin(Omega_M*tau_c)*exp(-3*G_Sr*tau_c/4));
%     for i = 1:N_point_theta
%         %I(i) = sum(sum(R.*rho.*(s./(1+s) + 2.*cos(2*k*Z*cos(theta(i))).*ss)));  % Intensity as a function of theta (fringes)
%         ss_avg = sum(s./(1+s) + 2*real(exp(-sqrt(-1)*2*k*(Z+Zsw)*cos(theta(i))).*ss),3);    % integral just on a period of the standing wave
%         I(i) = sum(sum(R_sum.*rho_sum.*ss_avg));  % Intensity as a function of theta (fringes);
%     end
% %     figure(10)
% %     hold on
% %     Intensities(j,:) = I;
% %     plot(theta-theta_0,I/mean(I))
%     Contraste_PW_Attenuation(j) = (max(I)-min(I))/(max(I)+min(I))*2;     % Calculation of the fringe contrast
% end
% 
% figure (3)
% hold on
% plot(s_0/4,Contraste_PW_Attenuation,'r');
% xlabel('\fontsize{16} s')
% ylabel('\fontsize{16} C')

%% Section 4 - Reduction of the contrast due to the amplitude variation of the probe beam (Gaussian beam) + attenuation + Dirac
% Calculation with an incident laser beam centered on the atomic cloud and
% as if the amplitude was not tilted
% Calculations done using the equation of I as a function of s, Omega_0 and Omega_M
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r_max = 3*min(w0, max([sigma_x,sigma_y]));  % Integration limit
z_max = 3*min(w0, sigma_z);  % Integration limit
N_points_r = 50;                                    % Number of points for r to perform the integral
r = linspace(0,r_max,N_points_r);
N_points_z = 200;                                   % Number of points for z to perform the integral
z = linspace(-z_max-h,z_max-h,N_points_z);
dz = z(2)-z(1);
N_points_zsw = 20;                                   % Number of points for z to perform the integral just on a period of the standing wave
zsw = linspace(0,lambda*(N_points_zsw-1)/N_points_zsw,N_points_zsw);
[R,Z,Zsw] = meshgrid(r,z,zsw);

rho = exp(-(R.^2/2/sigma_x^2 + (Z+h).^2/2/sigma_z^2));      % Atomic density profile
R_sum = sum(R,3);
rho_sum = sum(rho,3);

rho_sum_BL = sum(rho, 1);                % Sum along dimension 2 → (Nx, 1, Nz)

rho_center = zeros(size(rho));
disp(size(rho_center))
for ii = 1: N_points_z
    rho_center(ii,:,:) = rho(ii,:,:)./rho_sum_BL(1,:,:)/dz*b0;  %Normalized for Beer-Lambert lawZ
end

Intensities = zeros (N_points_Omega,N_point_theta);

% for j = 1:N_points_Omega 
%     %Omega versus Z
%     s_0_z = zeros(1,length(z));
%     s_0_z(1) = 2*(Omega_0(j)/G_Sr)^2;
% 
%     %Saturation variation for the incoming beam
%     for i=1:(length(z)-1)
%         ds_0 = -rho_center(i)/(1+s_0_z(i))*s_0_z(i)*dz;
%         s_0_z(i+1) = s_0_z(i)+ds_0;
%     end
%     
%     %Saturation variation for the reflected beam
%     z_reflected = flip(z);
%     s_0_z_reflected = zeros(1,length(z));
%     s_0_z_reflected(1) = s_0_z(end);
%     rho_center_reflected = flip(rho_center);
% 
%     for i=1:(length(z)-1)
%         ds_0_reflected = -rho_center_reflected(i)/(1+s_0_z_reflected(i))*s_0_z_reflected(i)*dz;
%         s_0_z_reflected(i+1) = s_0_z_reflected(i)+ds_0_reflected;
%     end
%     
%     %Interference of incoming and reflected beam
%     Omega_z = sqrt(s_0_z/2)*G_Sr;
%     [~, Omega_Z, ~] = meshgrid(r, Omega_z, zsw);
%     Omega_z_reflected = flip(sqrt(s_0_z_reflected/2)*G_Sr);
%     [~, Omega_Z_reflected, ~] = meshgrid(r, Omega_z_reflected, zsw);
% 
%     Omega_total = abs(Omega_Z.*exp(sqrt(-1)*k*(Z+Zsw)*cos(theta_0)) + Omega_Z_reflected.*exp(-sqrt(-1)*k*(Z+Zsw)*cos(theta_0)));
% 
%     %Next
%     Omega = Omega_total;    %2*Omega_0(j)*cos(k*(Z+Zsw)*cos(theta_0));
%     Omega_M = sqrt(Omega.^2-G_Sr^2/16);
%     s = 2*(abs(Omega)).^2/G_Sr^2;
%     
%     ss = s./(4*(1+s)).*(2./(1+s) + exp(-G_Sr*tau_c/2) + (s-1)./(s+1).*cos(Omega_M*tau_c)*exp(-3*G_Sr*tau_c/4) + G_Sr./(4*Omega_M).*(5*s-1)./(s+1).*sin(Omega_M*tau_c)*exp(-3*G_Sr*tau_c/4));
%     for i = 1:N_point_theta
%         %I(i) = sum(sum(R.*rho.*(s./(1+s) + 2.*cos(2*k*Z*cos(theta(i))).*ss)));  % Intensity as a function of theta (fringes)
%         ss_avg = sum(s./(1+s) + 2*real(exp(-sqrt(-1)*2*k*(Z+Zsw)*cos(theta(i))).*ss),3);    % integral just on a period of the standing wave
%         I(i) = sum(sum(R_sum.*rho_sum.*ss_avg));  % Intensity as a function of theta (fringes);
%     end
% %     figure(10)
% %     hold on
% %     Intensities(j,:) = I;
% %     plot(theta-theta_0,I/mean(I))
%     Contraste_PW_Attenuation(j) = (max(I)-min(I))/(max(I)+min(I))*2;     % Calculation of the fringe contrast
% end

b_0_final = zeros(N_points_Omega);
b_0_final_reflected = zeros(N_points_Omega);

for j = 1:N_points_Omega    
    %Omega versus Z
    s_0_Z = zeros(size(Z));
    s_0_Z(1,:,:) = 2*(Omega_0(j)/G_Sr)^2.*exp(-2*R(1,:,:).^2/w0^2);

    %Saturation variation for the incoming beam
    for i=1:(length(z)-1)
        ds_0 = -rho_center(i,:,:)./(1+s_0_Z(i,:,:)).*s_0_Z(i,:,:)*dz;
        s_0_Z(i+1,:,:) = s_0_Z(i,:,:)+ds_0;
    end
    
    b_0_final(j) = -log(s_0_Z(end,N_points_r/2,N_points_zsw/2)/s_0_Z(1,N_points_r/2,N_points_zsw/2));
    
    
    %Saturation variation for the reflected beam
    Z_reflected = flip(Z, 2);
    s_0_Z_reflected = zeros(size(Z));
    s_0_Z_reflected(1,:,:) = s_0_Z(end,:,:);

    rho_center_reflected = flip(rho_center,2);

    for i=1:(length(z)-1)
        ds_0_reflected = -rho_center_reflected(i,:,:)./(1+s_0_Z_reflected(i,:,:)).*s_0_Z_reflected(i,:,:)*dz;
        s_0_Z_reflected(i+1,:,:) = s_0_Z_reflected(i,:,:)+ds_0_reflected;
    end
    
    b_0_final_reflected(j) = -log(s_0_Z_reflected(end,N_points_r/2,N_points_zsw/2)/s_0_Z_reflected(1,N_points_r/2,N_points_zsw/2));
    
    figure(100)
    plot(s_0_Z(:,N_points_r/2,N_points_zsw/2))
    hold on;
    
    %Interference of incoming and reflected beam
    Omega_Z = sqrt(s_0_Z/2)*G_Sr;
    Omega_Z_reflected = flip(sqrt(s_0_Z_reflected/2)*G_Sr,2);

    Omega_total = abs(Omega_Z.*exp(sqrt(-1)*k*(Z+Zsw)*cos(theta_0)) + Omega_Z_reflected.*exp(-sqrt(-1)*k*(Z+Zsw)*cos(theta_0)));
    
    %Next
    Omega = Omega_total;    %2*Omega_0(j)*cos(k*(Z+Zsw)*cos(theta_0)).*exp(-R.^2/w0^2);
    Omega_M = sqrt(Omega.^2-G_Sr^2/16);
    s = 2*(abs(Omega)).^2/G_Sr^2;
%     
%     ss = s./(4*(1+s)).*(2./(1+s) + exp(-G_Sr*tau_c/2) + (s-1)./(s+1).*cos(Omega_M*tau_c)*exp(-3*G_Sr*tau_c/4) + G_Sr./(4*Omega_M).*(5*s-1)./(s+1).*sin(Omega_M*tau_c)*exp(-3*G_Sr*tau_c/4));
%     for i = 1:N_point_theta
%         %I(i) = sum(sum(R.*rho.*(s./(1+s) + 2.*cos(2*k*Z*cos(theta(i))).*ss)));  % Intensity as a function of theta (fringes)
%         %I(i) = sum(sum(R.*rho.*(s./(1+s) + 2*real(exp(-sqrt(-1)*2*k*Z*cos(theta(i))).*ss))));  % Intensity as a function of theta (fringes);
%         ss_avg = sum(s./(1+s) + 2*real(exp(-sqrt(-1)*2*k*(Z+Zsw)*cos(theta(i))).*ss),3);
%         I(i) = sum(sum(R_sum.*rho_sum.*ss_avg));  % Intensity as a function of theta (fringes);    
%     end
    
    P_el = 1./(1+s);
    P_in = s./(1+s);
    g1_el = 1;
    g1_in = 1;
    g1_total = P_el.*g1_el + P_in.*g1_in;

    for i = 1:N_point_theta
%         ss_avg = sum(s./(1+s).*(1 + real(exp(-sqrt(-1)*2*k*(Z+Zsw)*cos(theta(i))).*g1_total)),3);
        ss_avg = sum(s./(1+s).*(1 + g1_total.*cos(2*k*(Z+Zsw)*cos(theta(i)))),3);
        I(i) = sum(sum(R_sum.*rho_sum.*ss_avg));  % Intensity as a function of theta (fringes);   
    end
%     figure(20)
%     hold on
%     plot(theta-theta_0,I/mean(I))
%     Intensities(j,:) = I;
    Contraste_GaussianBeam_Attenuation_Dirac(j) = (max(I)-min(I))/(max(I)+min(I))*2;     % Calculation of the fringe contrast
end

figure(101)
plot(s_0/4, b_0_final); hold on;
plot(s_0/4, b_0_final_reflected);

fid = fopen( 'Attenuation_Mathilde.txt', 'w' );
for n=1:N_points_Omega
    fprintf(fid, '%.3f  ', s_0/4);
    fprintf(fid, '%.3f  ', b_0_final);
    fprintf(fid, '%.3f  ', b_0_final_reflected);
    fprintf(fid,'\r\n');
end
fclose(fid);

fid = fopen( 'Fringes_centered_beam.txt', 'w' );
for n=1:N_points_Omega
    for m=1:N_point_theta
        fprintf( fid, '%.1f  ', Intensities(n,m));
    end
    fprintf(fid,'\r\n');
end
fclose(fid);

figure (3)
hold on
plot(s_0/4,Contraste_GaussianBeam_Attenuation_Dirac,'m');
xlabel('\fontsize{16} s')
ylabel('\fontsize{16} C')

fid = fopen( 'Contrast_Mathilde_Dirac.txt', 'w' );
for n=1:N_points_Omega
    fprintf(fid, '%.3f  ', s_0(n)/4);
    fprintf(fid, '%.3f  ', Contraste_GaussianBeam_Attenuation_Dirac(n));
    fprintf(fid,'\r\n');
end
fclose(fid);

%%
legend('Gaussian beam','Gaussian beam + attenuation', 'Plane Weve + attenuation')