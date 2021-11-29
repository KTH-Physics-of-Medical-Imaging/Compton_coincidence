%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Example script of maximum likelihood estimation
%Version: 211122 Christel Sundberg 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('Maximum_likelihood_parameters.mat')
%%
rng('shuffle');
maximum_likelihood = [];
E_gamma = [];
kVp = 120; %kVp x-ray tube voltage
eActVectorkeV = 15:0.1:kVp; %keV
sigma_p = photo(141:end)./(photo(141:end) + compton(141:end)); %Relative cross section, photo
sigma_c = compton(141:end)./(photo(141:end) + compton(141:end)); %Relative cross section, Compton
mu_original_E = 1:eActVectorkeV(2)-eActVectorkeV(1):kVp; 
sigma_p_original = photo./(photo + compton);
sigma_c_original = compton./(photo + compton);

energy_res = 0.5; %keV Energy resolution represented by a Gaussian with zero mean and standard deviation energy_res
spatial_res = [0.001 0.05 0.05]; %cm Spatial resolution represented by three Gaussians (x,y,z) with zero mean and standard deviation spatial_res in (x,y,z), respectively. 
%%
interaction_chain = Chain_2Compton; %Manually input which interaction chain should be analyzed
combinations = perms([1 2]); %All possible interaction orders, e.g. for 2 Compton interactions: combinations = perms([1 2]);
C_2Compton_observed = cell(1,length(interaction_chain)); %Saves data with errors. 
for i = 1:length(interaction_chain) 
    C_temp = interaction_chain{i};
    error_x_first = zeros(size(C_temp(:,2:4))); %Holds errors in position
    error_E_first = zeros(size(C_temp(:,1))); %Holds errors in energy
    for b = 1:length(C_temp(:,2)) %Obtain observed data by adding errors to simulated data
        error_x_first(b,:) = ([spatial_res(1)*sqrt(-2*log(rand()))*cos(2*pi*rand()) spatial_res(2)*sqrt(-2*log(rand()))*cos(2*pi*rand()) spatial_res(3)*sqrt(-2*log(rand()))*cos(2*pi*rand())]);
        error_E_first(b) = round(energy_res*sqrt(-2*log(rand()))*cos(2*pi*rand()),1);
        C_temp(b,1) = C_temp(b,1)+ error_E_first(b); %Error in energy
        if C_temp(b,1) <= 0 %Ensure energy is not below 0 keV
            C_temp(b,1) = 0.01;
        end
    end 
    C_temp(:,2:4) = C_temp(:,2:4) + error_x_first; %Error in position
    C_2Compton_observed{i} = C_temp; %Save observed data
    for k = 1:size(combinations,1) %Evaluate the likelihood for each possible interaction order
        C = [C_temp(combinations(k,:),:)];
        tic
        [maximum_likelihood(i,:,k), E_gamma(i,k)] = run_likelihood_estimation(C,S_tot,mu_original_E,mu_original,sigma_p_original,sigma_c_original,mu,sigma_c,sigma_p,energy_res,spatial_res,eActVectorkeV); %Output: likelihood function as a function of incident photon energy E_gamma along with the estimated incident photon energy E_gamma, obtained as the energy coordinate of the maximum of the likelihood function.
        toc
    end
end

function [maximum_likelihood, E_gamma] = run_likelihood_estimation(C,S_tot,mu_original_E,mu_original,sigma_p_original,sigma_c_original,mu,sigma_c,sigma_p,energy_res,spatial_res,eActVectorkeV)

C_temp = C;
likelihood = 0;

error_x =  zeros([length(C(:,2)),3]);
error_E = zeros([1000,length(C(:,2))]);
error_tempE = zeros(size(C(:,2)));
for a = 1:1500 %Number of samples in Monte Carlo integration   
    for b = 1:length(C(:,2)) %Number of interactions
        error_x(b,:) = ([spatial_res(1)*sqrt(-2*log(rand()))*cos(2*pi*rand()) spatial_res(2)*sqrt(-2*log(rand()))*cos(2*pi*rand()) spatial_res(3)*sqrt(-2*log(rand()))*cos(2*pi*rand())]); %Sample error in position
        error_tempE(b) = round(energy_res*sqrt(-2*log(rand()))*cos(2*pi*rand()),1); %Sample error in energy
    end 
    error_E(a,:) = error_tempE;
    interaction_energy = C_temp(:,1)+error_E(a,:)'; %Deposited energy
    interaction_position = [C_temp(:,2), C_temp(:,3), C_temp(:,4)]+error_x; %Interaction positions

    E_1 = interaction_energy(1);
    x_1 = interaction_position(1,:);
    x_2 = x_1;
    L_d = 8; %cm Detector depth
    dirac_x = 1; %Dummy
    temp_E = eActVectorkeV;
    exit_length = [x_1(1) x_1(2) L_d]; %Coordinate of detector edge for incident photon

    eActVectorkeV_temp = 0:eActVectorkeV(2)-eActVectorkeV(1):eActVectorkeV(end);
    dirac_E_temp(1,:) = 1/(energy_res*sqrt(2*pi))*exp(-(eActVectorkeV_temp-round(E_1,1)).^2/(2*energy_res^2)); %Gaussian representing energy resolution
    dirac_E_temp(1,:) = dirac_E_temp(1,:)/sum(dirac_E_temp(1,:));
    dirac_E(1,:) = dirac_E_temp(1,end-length(eActVectorkeV)+1:end); %Non-ideal case including limited energy resolution
    P_interacting = [];

    beer(1,:) = mu.*exp(-mu*x_1(3))./(1-exp(-mu*L_d)); %Probability of interacting in position x_1(3)
    P_escape = [];

    P_interacting(1,:) = (1-exp(-mu*L_d)); %Probability of interacting at, Beer Lambert
    likelihood_temp = dirac_x.*beer(1,:).*P_interacting(1,:);
    vector1 = x_1 - [x_1(1),x_1(2),0];
    last = 0;

    L_exit = [];
    L_exit(1) = L_d;
    KC1 =[];
    theta = [];
    for i = 1:length(C(:,1))-1 
        x_1 = x_2;
        x_2 = interaction_position(i+1,:);
        L_exit(i+1) = Escape_length(x_1,x_2,L_d);
        temp_E = temp_E - E_1;
        temp_E(temp_E < 1 ) = [];
        dirac_E_temp(i+1,:) = 1/(energy_res*sqrt(2*pi))*exp(-(eActVectorkeV_temp-sum(interaction_energy(1:i))-round(interaction_energy(i+1),1)).^2/(2*energy_res^2));
        dirac_E_temp(i+1,:) = dirac_E_temp(i+1,:)/sum(dirac_E_temp(i+1,:));
        dirac_E(i+1,:) = dirac_E_temp(i+1,end-length(eActVectorkeV)+1:end); %Including non-ideal energy resolution
        mu_start = find(round(mu_original_E,1) == round(temp_E(1),1));
        mu_end = find(round(mu_original_E,1) == round(temp_E(end),1));
        sigma_c(i+1,:) = [zeros([1 size(eActVectorkeV,2)-length(temp_E)]), sigma_c_original(mu_start:mu_end)]; %Relative interaction cross section, Compton
        sigma_p(i+1,:) = [zeros([1 size(eActVectorkeV,2)-length(temp_E)]), sigma_p_original(mu_start:mu_end)]; %Relative interaction cross section, photo
        mu_scattered = [zeros([1 size(eActVectorkeV,2)-length(temp_E)]), mu_original(mu_start:mu_end)]; %cm^-1 Attenuation coefficient of scattered photon  

        P_interacting(i+1,:) = (1-exp(-mu_scattered.*L_exit(i+1))); %Probability of interacting at all, Beer Lambert

        beer(i+1,:) = mu_scattered.*exp(-mu_scattered*norm(x_2-x_1))./(1-exp(-mu_scattered*L_exit(i+1))); %Probability of interaction position
        beer(i+1,isnan(beer(i+1,:))) = 0; %mu_scattered can be 0, remove NaN entries
        vector2 = x_2-x_1; %Direction of photon
        theta(i) = acos(dot(vector1,vector2)/(norm(vector1)*norm(vector2))); %Scattering angle theta

        KC1(i,:) = ComptonScattering(eActVectorkeV,interaction_energy(i),theta(i),S_tot,last,sum(interaction_energy(1:i-1))*(i>1)); %Probability of depositing energy interaction_energy(i) and scattering at angle theta(i)
        vector1 = vector2;
        E_1 = interaction_energy(i+1);

        likelihood_temp = likelihood_temp.*beer(i+1,:).*P_interacting(i+1,:).*sigma_c(i,:).*KC1(i,:);
    end
    temp_E = temp_E - E_1;
    temp_E(temp_E < 1 ) = [];
    mu_start = find(round(mu_original_E,1) == round(temp_E(1),1));
    mu_end = find(round(mu_original_E,1) == round(temp_E(end),1));
    mu_scattered = [zeros([1 size(eActVectorkeV,2)-length(temp_E)]), mu_original(mu_start:mu_end)]; %cm^-1 Attenuation coefficient of scattered photon         
    
    [P_escape, temp_thetas2] = P_escape_function(x_1,x_2,eActVectorkeV,interaction_energy(i+1),mu_scattered,S_tot,L_d,sum(interaction_energy(1:i)));
    KC1(i+1,:) = ComptonScattering(eActVectorkeV,interaction_energy(i+1),theta(i),S_tot,1,sum(interaction_energy(1:i)));
    likelihood_temp = likelihood_temp.*((P_escape.*KC1(i+1,:).*sigma_c(i+1,:)+sigma_p(i+1,:).*dirac_E(i+1,:)));
    likelihood = likelihood + likelihood_temp;
end
likelihood = likelihood/a;
[~, I] = max(likelihood);
maximum_likelihood = likelihood;
E_gamma = eActVectorkeV(I);
end

function L_exit = Escape_length(x_1,x_2,L_d)
x_direction = x_2 - x_1;
if x_direction(3) > 0 
    z_max = L_d; 
elseif x_direction(3) < 0 
    z_max = 0;
elseif x_direction(3) == 0
    z_max = 0;
end
if x_direction(2) == 0 
    x_2(2) = x_2(2)+1e-4;
end
if x_direction(1) == 0 
    x_2(1) = x_2(1)+1e-4;
end
y_exit = (z_max - x_1(3))*x_direction(2)/x_direction(3) + x_1(2);
x_exit = (z_max - x_1(3))*x_direction(1)/x_direction(3) + x_1(1);
L_exit = norm([x_exit y_exit z_max]-x_1);
end
function KC = ComptonScattering(eActVectorkeV,E,theta_scatter,S_tot,last,corr_E)
eActVectorkeV2 = eActVectorkeV; %keV
KC = zeros([size(eActVectorkeV2)]);
E_f_vector = 0:0.1:160;
dpi = 0.02;
theta = 0:dpi:pi;
m_e = 9.109e-31; %kg
c = 3e8; %m/s 
E_gamma = eActVectorkeV2;
E_prime = E_gamma-E;
if last == 1
    for i = 1:length(eActVectorkeV2)
        if E_prime(i) > 0 
            I = round((E_prime(i))*10+1);
            KC(i) = (sum((S_tot(:,I,i))));
        else 
            KC(i) = 0;
        end
    end 
elseif last == 0  
    I1 = find(theta >= theta_scatter,1);
    for i = 1:length(eActVectorkeV2)
        if E_prime(i) > 0 
            I2 = round((E_prime(i))*10+1);
            try
            KC(i) = S_tot(I1,I2,i);
            catch
                KC(i) = S_tot(end,I2,i);
            end
        else 
            KC(i) = 0;
        end
    end
end
KC_temp = zeros(size(eActVectorkeV2));
KC_temp((1+round(corr_E*10)):end) = KC(1:(end-round(corr_E*10)));
KC = KC_temp;
end
function [P_escape, temp_thetas2] = P_escape_function(x_1,x_2,eActVectorkeV,E_1,mu_2,S_tot,L_d,corr_E)
eActVectorkeV2 = eActVectorkeV;
u = (x_2-x_1); %Incident direction of photon
u = u/norm(u);
z = [0 0 1];
alpha = acos(dot(z,u)/(norm(z)*norm(u))); %Angle between z and u 
r = cross(u,z); %Orthogonal to u,z, and e 
P_escape = zeros(size(eActVectorkeV2));
m_e = 9.109e-31;
c = 3e8;
dphi = 0.04;
dpi = 0.02;
theta = 0:dpi:pi;
E_f_vector = 0:0.1:160;
E_f_vector = round(E_f_vector,1);
phi = 0:dphi:2*pi; %Rotate e with angle phi around u 
L_exit = zeros([size(phi,2),length(theta)]);
R = zeros([3*length(phi) 3]);
for i = 1:length(phi) 
	R(i*3+1-3:i*3,:) = [cos(phi(i))+u(1)^2*(1-cos(phi(i))) u(1)*u(2)*(1-cos(phi(i)))-u(3)*sin(phi(i)) u(1)*u(3)*(1-cos(phi(i)))+u(2)*sin(phi(i)); u(2)*u(1)*(1-cos(phi(i)))+u(3)*sin(phi(i)) cos(phi(i))+u(2)^2*(1-cos(phi(i))) u(2)*u(3)*(1-cos(phi(i)))-u(1)*sin(phi(i)); u(3)*u(1)*(1-cos(phi(i)))-u(2)*sin(phi(i)) u(3)*u(2)*(1-cos(phi(i)))+u(1)*sin(phi(i)) cos(phi(i))+u(3)^2*(1-cos(phi(i)))];    
end
for n = 1:length(theta)
    if theta(n) > alpha
    	beta = theta(n)-alpha;
    else 
        beta = alpha-theta(n); %Angle between e and z
    end
    if theta(n) > pi/2
        norm_e = norm(u)/cos(pi-theta(n));
    elseif theta(n) < pi/2
        norm_e = norm(u)/cos(theta(n));
    end    
    e_3 = norm_e*cos(beta);
    e_2 = (norm(u)*norm_e*cos(theta(n))-norm_e*cos(beta)*u(3))/(u(2)+u(1)*(-r(2))/r(1));
    e_1 = -r(2)/r(1)*e_2;
    e = [e_1 e_2 e_3]';
    e_phi = zeros([size(phi,2),3]);
    e_phi = R*e;
    e_phi = reshape(e_phi,[3, size(phi,2)])';
    x_3 = e_phi + x_2;
    x_first = x_2;
    x_second = x_3;
    x_gamma = x_first;
    i = 1;
    x_direction = x_second-x_first;
    
    for k = 1:size(x_direction,1)
        if x_direction(k,3) > 0 
            z_max = L_d; 
        elseif x_direction(k,3) < 0 
            z_max = 0;
        elseif x_direction(k,3) == 0
            z_max = 0;
        end
        y_exit = (z_max - x_first(3))*x_direction(k,2)/x_direction(k,3) + x_first(2);
        x_exit = (z_max - x_first(3))*x_direction(k,1)/x_direction(k,3) + x_first(1);
        L_exit(k,n) = sqrt((x_exit(i)-x_first(1)).^2+(y_exit(i)-x_first(2)).^2+(z_max-x_first(3)).^2);
    end
end
dphi_temp = dphi/(phi(end)+dphi);
input_temp = exp(-L_exit(:)*mu_2);
input_temp = reshape(input_temp,[size(L_exit) 1051]);
input_temp2 = sum(input_temp);
temp_thetas2 = 0;

for i = 1:length(eActVectorkeV)
    E_gamma = eActVectorkeV2(i);
    E_prime(i) = E_gamma-E_1-corr_E; %Scattered photon energy
    if E_prime(i) > 0 && floor(i)-round(corr_E*10) > 0 
        I = round((E_prime(i))*10+1);
        temp_thetas = (S_tot(:,I,(floor(i)-round(corr_E*10))));
        temp_thetas = temp_thetas/sum(temp_thetas);
        [~,I2] = max(temp_thetas);
        temp_thetas2(i) = theta(I2);
        P_escape(i) = input_temp2(:,:,i)*dphi_temp*temp_thetas;  
    else 
        P_escape(i) = 0;
    end
end
end
