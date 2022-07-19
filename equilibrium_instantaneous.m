function [] = biomech_analysis_equilibrium_mohammad(~)
%%%%%%%%%%%%%%%%%%%%%%%
% Mohammadhossein Ebrahimi 14.2.2018
%
% This script is to analyze EQUILIBRIUM and instantaneous moduli based on s
% multiple tress-relaxation protocol
% 
%
% Script loads the 4 step biomechanical testing data, finds the equilibrium
% (full relaxation) points, fits a line into the data and calculates the
% Youngs modulus for the samples.
% it also calucaltes instantaneous moduli from each loading step, fits a line and
% extracts initial and strain-dependent instantaneous moduli.
%
% some variables like relaxation time, strain, strain rate, the unit of force and displacement, 
% number of relaxation points etc, must be modified according to the measurement protocol
%
% This code also includes Correction of modulus based on Hayes correction factor.
%
% Mohammadhossein Ebrahimi 14.2.2018
%
%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%   EQUILIBRIUM MODULUS   %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


disp('Pick up the biomechanical testing data folder including both static and dynamic loading.')
folder = uigetdir(pwd, 'Pick up the biomechanical testing data folder including both static and dynamic loading.');

folder_stress_strain = strcat(folder, '\01_Stress_Strain');
folder_dynamic = strcat(folder, '\02_Dynamic');

% add path where choosen file is
cd(folder_stress_strain)

% Fetch the filenames from the folder:
d = dir();
files = {d(1:end,1).name}';
files(ismember(files,{'.','..'})) = [];
filename = files(1);
Biomech_data = importdata(char(filename));

% Load data from the raw-file:
data = Biomech_data.data;
information = Biomech_data.textdata;
Biomech_data.name = filename;

% Rename data to separate matrixes
time = (data(:,1));
displacement = (-1*data(:,2));
load1 = (-1*data(:,3));
bs = nanmean(load1(3:max((find(time<5)))));
load = load1-bs;    % Normalize data to the first measurement value

% Filter the measured signal with median filter:
window = 10;                    % window of the filter
fil = medfilt1(load,window);    % filtered signal


% calculate sampling frequency:
format long
df = diff(time);
dff = df(find(df~=0));
fs = floor(nanmean( 1./dff ));

% Apply relaxation time (seconds), this can be modified according to the measurement protocol
relax = 900;    % 900 sec = 15 min
% Apply how many steps are made in the test:
steps = 4; %this can be modified according to the measurement protocol
step = nan(steps+1,1);

% Times of the steps:
for kl = 1:steps+1;
    step(kl) = [relax*(kl-1)] + 5;
end

% Find the accurate points of the relaxation phases:
for ii = 1:steps+1;
    tt = find(step(ii) - mean(diff(time)) <= time & time <= step(ii) + mean(diff(time)));
    
    % Check it the threshold limits are too tight:
    if isempty(tt) == 1     
        
        tt = find(step(ii) - mean(diff(time)*2) <= time & time <= step(ii) + mean(diff(time))*2);
    
    end

    if size(tt,1) >= 1
        ts(ii) = ceil(mean(tt));

    else 
        ts(ii) = tt;
    end

        
end
    
ts = ts(:);
    
minimum = zeros(size(step,2),2);
win = 50;
for ii = 1:length(step);
    in_loc(ii) = ts(ii);
    in_val_meas(ii) = load(ts(ii));     % relaxation points from measured signal
    in_val(ii) = fil(ts(ii));           % relaxation points from filrered signal
    in_time(ii) = time(ts(ii));
    
    ind_relax(:,ii) = [ts(ii)-win:ts(ii)];
    
    in_val_mean(ii) = mean(fil(ind_relax(:,ii)));
    in_val_meas_mean(ii) = mean(load(ind_relax(:,ii)));
    minimum(ii,2) = in_time(ii);
    minimum(ii,1) = in_val_meas_mean(ii);   % relaxation points from measured signal
    minimum2(ii,1) = in_val_mean(ii);       % relaxation points from filtered signal

end


% Plot figure with all information from the data:
h1 = figure('units','normalized','position',[0.005 0.04 0.492 0.88]);
subplot 211, plot(time,displacement), title('Indenter displacement'), xlabel('time'), ylabel('Disp. (µm)')
subplot 212, plot(time,load), title('load'), xlabel('time'), ylabel('load (g)')
hold on, subplot 212, plot(minimum(:,2),minimum(:,1),'ro')


% fit line with "least squares"-method to the equilibrium points
% in a loop iteration such that last 4, 3 and 2 points are analyzed


r = (0.73e-3)/2;   % radius of the indenter
A = pi*r.^2;    % indented area



stress = (minimum(:,1).*9.81*1e-3)./A;  % equilibrium stresses applied to the sample (Pa)
stress = stress(:);

% Strains (according to the steps adapted previously):
strain = nan(steps+1,1);
for kl = 1:steps+1;
    strain(kl) = [(1-0.95^(kl-1))];
end
strain = strain(:);

% CALCULATE INSTANT MODULAE:
% Find the load peaks from the filtered data:
[pks, locs] = findpeaks(fil,'MINPEAKDISTANCE',mean(diff(fs*step))*0.9, 'NPEAKS', steps);

% Find the load peaks from the raw data:
ww = 5;
for kk = 1:steps;
    [mm(kk,1) mm(kk,2)] = max(load(locs(kk)-ww:locs(kk)+ww));
end

pks2 = mm(:,1);
locs2 = mm(:,2) + locs(1:steps,:) - ww - 1;

maximum = [load(locs) time(locs)];     % create matrix with load peaks and time points they appear

% Create matrix with the instant loads (load peaks - (minus) previous relaxation load)
instant = maximum(:,1) - minimum(1:end-1,1);    % relaxation points from measured signal
instant2 = maximum(:,1) - minimum2(1:end-1,1);  % relaxation points from filtered signal

h1 = figure('units','normalized','position',[0.005 0.04 0.492 0.88]);

subplot(4,3,1:3), plot(time,displacement), title('Indenter displacement'), xlabel('time'), ylabel('Disp. (µm)'), axis tight
subplot(4,3,4:6), plot(time,load), title('load'), xlabel('time'), ylabel('load (g)')
hold on, subplot(4,3,4:6), plot(time(ind_relax), load(ind_relax),'ro')
subplot(4,3,4:6), plot(time(locs2),pks2,'go')
subplot(4,3,7:9), plot(time,fil), title('FILTERED load'), xlabel('time'), ylabel('load (g)')
hold on, subplot(4,3,7:9), plot(time(ind_relax), fil(ind_relax),'ro')
subplot(4,3,7:9), plot(time(locs),pks,'go')

legend('Load', 'Relaxation points' , 'Load peaks','Location', 'Best')

disp(['Sample: ', char(files(1))])

% select manually the points which wanted to be analysed:
po = [1 5; 2 5; 2 4];

for kk = 1:3

i1 = po(kk,1); i2 = po(kk,2);
    
H = [minimum(i1:i2,2).^1 minimum(i1:i2,2).^0];
theta(:,kk) = H\minimum(i1:i2,1);

xx = linspace(minimum(i1,2),minimum(i2,2),1000);
yy = theta(1,kk)*xx + theta(2,kk);



% LEAST SQUARES FIT TO STRESS/STRAIN DATA:
HH = [strain(i1:i2).^1 strain(i1:i2).^0];         % observation model
theta_hat_ss(:,kk) = HH\stress(i1:i2);                 % slope and coefficient of the least squares solution
stress_hat = HH*theta_hat_ss(:,kk);                       % estimated values to stress
sg_hat = (norm(stress(i1:i2) - stress_hat).^2)/(size(HH,1)-size(HH,2));    % variance^2
Cth = sg_hat*inv(HH'*HH);                           % Covariance of the THETA value
mth(:,kk) = sqrt(diag(Cth));                              % Error of theta values


subplot(4,3,kk+9), plot(strain, stress,'ro')
hold on, subplot(4,3,kk+9), plot(strain(i1:i2), strain(i1:i2)*theta_hat_ss(1,kk) + theta_hat_ss(2,kk),'b--')
xlabel('Strain (%)'), ylabel('Stress [Pa]')
title(['Youngs modulus: ', num2str(theta_hat_ss(1,kk)/1e6), ' \pm ', num2str(mth(1,kk)/1e6), ' [Pa]'])
axis 'tight'


end

Eq_mod = theta_hat_ss(1,:)*1e-6;
Eq_mod = Eq_mod(:);

% The same calculations with Hayes et al 1972 - correction factor
v = 0.3;        % Poisson's ratio
% r = 0.001/2;     % indenter radius
thickness = str2double(strtok(cell2mat(Biomech_data.textdata(6)),'zero strain, µm: ')); % thickness of the sample [µm]
h = thickness*10^-6;

% Calculate KAPPA-value
kappa_eq = kappa_value(r, h, v);

% Corrected Equilibrium modulae (Hayes et al 1972, check also doctoral
% thesises of Petro Julkunen and Janne Mäkelä)
Es = ((1-v.^2)*pi*r*Eq_mod) ./ (2*kappa_eq*h);

% calculate Hayes corrected instant modulae:
v = 0.5;   % For instant modulus analysis, Poisson's ratio is assumed to be v = 0.5

% Kappa values:
kappa1 = kappa_value(r, h*0.95 , v);
kappa2 = kappa_value(r, h*0.95^2 , v);
kappa3 = kappa_value(r, h*0.95^3 , v);
kappa4 = kappa_value(r, h*0.95^4 , v);


kappa_instant = [kappa1 kappa2 kappa3 kappa4]';
% Corrected Instant modulae (Hayes et al 1972)
Es_instant = ((1-v.^2)*pi*r*instant) ./ (2*kappa_instant*h);    % relaxation points from measured signal
Es_instant2 = ((1-v.^2)*pi*r*instant2) ./ (2*kappa_instant*h);  % relaxation points from filtered signal



% % print additional information:
disp(' ')
disp('Uncorrected equilibrium moduli: ')
disp(['Points 1-5: ', num2str(theta_hat_ss(1,1)/1e6), ' ± ', num2str(mth(1,1)/1e6), ' [MPa]'])
disp(['Points 2-5: ', num2str(theta_hat_ss(1,2)/1e6), ' ± ', num2str(mth(1,2)/1e6), ' [MPa]'])
disp(['Points 2-4: ', num2str(theta_hat_ss(1,3)/1e6), ' ± ', num2str(mth(1,3)/1e6), ' [MPa]'])
disp(' ')

disp(' ')
disp('Corrected equilibrium moduli, aka Youngs moduli: ')
disp(['Points 1-5: ', num2str(Es(1)), ' [MPa]'])
disp(['Points 2-5: ', num2str(Es(2)), ' [MPa]'])
disp(['Points 2-4: ', num2str(Es(3)), ' [MPa]'])
disp(' ')
disp(['Sample thickness: ', num2str(thickness), ' [µm]'])
disp(' ')

disp(' ')
disp('Uncorrected instant moduli: ')
disp(['1st step: ', num2str(instant(1)), ' [MPa]'])
disp(['2nd step: ', num2str(instant(2)), ' [MPa]'])
disp(['3rd step: ', num2str(instant(3)), ' [MPa]'])
disp(['4rd step: ', num2str(instant(4)), ' [MPa]'])
disp(' ')

disp(' ')
disp('Hayes corrected instant moduli: ')
disp(['1st step: ', num2str(Es_instant(1)), ' [MPa]'])
disp(['2nd step: ', num2str(Es_instant(2)), ' [MPa]'])
disp(['3rd step: ', num2str(Es_instant(3)), ' [MPa]'])
disp(['4rd step: ', num2str(Es_instant(4)), ' [MPa]'])

disp(' ')

% LEAST SQUARES FIT TO STRESS/STRAIN DATA (FOR INSTANTANEOUS MODULUS):
HH = [strain(2:end).^1 strain(2:end).^0];         % observation model
theta_hat_imod(:,1) = HH\Es_instant;                 % slope and coefficient of the least squares solution
imod_hat = HH*theta_hat_imod(:,1);                       % estimated values to stress
sg_hat_imod = (norm(Es_instant - imod_hat).^2)/(size(HH,1)-size(HH,2));    % variance^2
Cth_imod = sg_hat_imod*inv(HH'*HH);                           % Covariance of the THETA value
mth_imod(:,1) = sqrt(diag(Cth_imod));                              % Error of theta values


strainvals = linspace(0,strain(end));
figure(3),clf
hold on
plot(strain(2:end), Es_instant,'ro')
plot(strainvals, strainvals*theta_hat_imod(1,1) + theta_hat_imod(2,1),'b--')
xlabel('Strain (%)'), ylabel('Modulus [MPa]')
title({['Strain-dependent instantaneous modulus: ', num2str(theta_hat_imod(1,1)), ' \pm ', num2str(mth_imod(1,1)), ' [MPa]'];[ 'Initial instantaneous modulus' num2str(theta_hat_imod(2,1)), ' \pm ', num2str(mth_imod(2,1)), ' [MPa]']})
hold off


%%%%% EXPORT DATA INTO .xlsx-file
caption = {};
caption(1,1) = filename; 
caption(1,2) = {['thickness (µm): ']}; caption(1,3) = {[num2str(thickness)]};


caption(2,2) = {'1st step'}; caption(2,3) = {'2nd step'};
caption(2,4) = {'3rd step'}; caption(2,5) = {'4th step'};

caption(3,1) = {'Relaxation stress (MPa) (3 steps)'};
caption(4,1) = {'Strain (%) (3 steps)'};

for ll = 1:steps;
    caption(3,1+ll) = {stress(ll+1)/1e6};
    caption(4,1+ll) = {strain(ll+1)};
end

caption(5,2) = {'Points 1-4'}; caption(5,3) = {'Points 2-4'}; caption(5,4) = {'Points 2-3'};
caption(6,1) = {'Uncorrected equilibrium modulus (MPa)'};
caption(7,1) = {'Corrected equilibrium modulus (MPa)'};


for ll = 1:size(po,1);
    caption(6,1+ll) = {Eq_mod(ll)};
    caption(7,1+ll) = {Es(ll)};
end

caption(8,1) = {'Kappa-value'};
caption(8,2) = {kappa_eq};   

caption(11,1) = {'Uncorrected instant modulus (MPa)'};
caption(12,1) = {'Corrected instant modulus (MPa)'};
caption(13,1) = {'Kappa value for instant modulus (each step separately)'};


caption(10,2) = {'1st load peak'}; caption(10,3) = {'2nd load peak'}; 
caption(10,4) = {'3rd load peak'}; caption(10,5) = {'4th load peak'};

caption(15,1) = {'Corrected instant strain-dependent modulus (MPa)'}; caption(15,2)={theta_hat_imod(1,1)}; caption(15,3) = {'±'}; caption(15,4)={mth_imod(1,1)};
caption(16,1) = {'Corrected instant initial modulus (MPa)'};  caption(16,2)={theta_hat_imod(2,1)}; caption(16,3) = {'±'}; caption(16,4)={mth_imod(2,1)};

for ll = 1:steps;
    caption(11,1+ll) = {instant(ll)};
    caption(12,1+ll) = {Es_instant(ll)};
    caption(13,1+ll) = {kappa_instant(ll)};
end




cd(folder)


% ADD INFO SHEET TO THE EXCEL-FILE:
info = ({'This excel-file contains some of the biomechanical analysis of ..... ' ; 
        '. 4 step relaxation tests were conducted and analysed. ' ;
        'In this analysis file, ONLY relaxation with points 1-5, 2-5 and 2-4 were done. ' ;
        'Poissons ratios are assumed to be: v=0.1 in the stress-relaxation test and v=0.5 in the dynamic and instant test. '
        });

xlswrite(strcat(char(filename), '.xlsx'), info , 'info') 

xlswrite(strcat(char(filename), '.xlsx'), caption , 'Eq. mod.')    

% Save image:
saveas(h1 , strcat(char(filename),'_eq') , 'fig')  % as a Matlab-file
saveas(h1 , strcat(char(filename),'_eq') , 'png')  % as a regular image

% Save into struct-format (to open with MATLAB later on):
Biomech_data_equilibrium.stress = stress;
Biomech_data_equilibrium.strain = strain;
Biomech_data_equilibrium.Eq_mod = theta_hat_ss(1,:);
Biomech_data_equilibrium.instant = instant;
Biomech_data_equilibrium.Es_instant = Es_instant;
Biomech_data_equilibrium.strain = strain;
Biomech_data_equilibrium.Eq_mod_error = mth(1,:);
Biomech_data_equilibrium.Kappa_factor_equilibrium = kappa_eq;
Biomech_data_equilibrium.Kappa_factor_instant = kappa_instant;
Biomech_data_equilibrium.Young_corrected = Es;
Biomech_data_equilibrium.Thickness = h;
Biomech_data_equilibrium.data = data;
Biomech_data_equilibrium.textdata = information;
Biomech_data_equilibrium.name = filename;


% save data and information to .mat-file
save(strcat(char(filename),'_eq.mat'),'-struct','Biomech_data_equilibrium');


end


%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % % % % % % % % % SEPARATE SUBFUNCTIONS % % % % % % % % % % % % %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%  CALCULATE KAPPA VALUE %%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [ kappa ] = kappa_value( radius, thickness , v )
% This function calculates tke KAPPA-value to the calculation of effective
% Young's modulus with respect to the sample thickness [m].
%
%   Effective Young's modulus:
%   Es = (1-v^2)*pi*a* E / 2*KAPPA*thickness
%   where, E = measured Young's modulus, a = radius of the indenter,
%   v = Poisson's ratio
%   thickness = sample thickness [m]
% 
%   Poisson's ratio is determined to be v = 0.1 in equilibrium modulus and 
%   v = 0.5 in dynamic modulus 
%
%
%
a = (radius);                         % indenter radius [m]


if v == 0.3
    ah = [0 : 0.2: 2]';
    % Hayes et al, 1972
    poisson = [1000 1207 1472 1784 2124 2480 2845 3214 3586 3960 4336]'.*10^-3; 
   
elseif v == 0.5
    ah = [0 : 0.2: 2]';
    % Hayes et al, 1972
    poisson = [1000 1268 1645 2129 2704 3359 4085 4878 5737 6659 7644]'.*10^-3;
else
    disp('Poissons coefficient you gave, is not v = 0.1 or v = 0.5. Please, give one of those. ')
    
end

spline_ah = 0:0.001:2;
%interpolate cubic spline for kappa values presented in Hayes et al. 1976 (or so)
spline_poisson = spline(ah,poisson,spline_ah);
%if a/h value is not excatly in the cubic spline data points, use linear
%interpolation for the kappa value between cubic spline data points 
%(error here is negligible due to the dense sampling in spline)
interp_poisson = interp1(spline_ah, spline_poisson, a/thickness);

%cubic spline approximated kappa
kappa = interp_poisson;

end

