% Script written for use with paper
% Written by Sam Wagner Aug 3, 2021
% INPUTS: 
%        - parameters in Table 1                         (Defaults given)
%        - .s2p file with simulated antenna coupling     (Default provided)
%        - .s2p file with simulated antenna transmission (Default provided)
% 
% .\utils contains support functions for this script.
%% housekeeping 
addpath('.\utils\');
clear all; close all; clc;                      % reset script
set(0,'defaultfigurecolor',      [1 1 1]);      % white figure backgrounds 
set(0,'defaultAxesFontSize',     20);           % font size
set(0,'defaultlinelinewidth',    2);            % line width
set(0,'defaultLegendInterpreter','latex');      % use latex
set(0,'defaultTextInterpreter',  'latex');      % use latex

%% Table 1 parameters & Input .s2p files
% fundamental GPR system parameters
V_pp_t  = 1;        % peak to peak voltage (Volts)
f_low   = 0.5e9;    % minimum operating frequency (Hz)
f_high  = 3e9;      % maximum operating frequency (Hz)
o_j     = 20e-12;   % system jitter or equivalent (s)
h       = 40e-2;    % expected antenna operating height (m)
d       = 30e-2;    % expected target depth (m)
d_max   = 1;        % maximum target depth (m)
G       = 6;        % antenna gain at center frequency (dB)
G_amp   = 14.5;     % receiver amplifier gain (dB)
P_1dB   = 19;       % receiver amplifier output 1dB compression point (dBm)
Z0      = 50;       % characteristic impedance
o       = 0;        % expected target radar cross section (dBsm)
er      = 5;        % expected ground relative permittivity
alpha   = 0.01;     % expected ground EM attenuation constant (Np/m)
dx      = 0;        % cross-track board-to-board separation (cm)
p_max   = 5.5;      % (UAV) - maximum payload
m_sys   = 1;        % (UAV) - expected weight of external system
N_ants  = 8;        % (UAV) - number of antennas total
V_noise_rms = 1e-3; % RMS noise of receiver 

% threshold parameters 
k_coupling      = 1;        % coupling slope threshold 
coupling_Zscore = 1.96;     % coupling slope Z score (Z_1-alpha/2)
k_distortion    = 0.9;      % pulse integrity threshold
k_RCS           = 1;        % RCS threshold
m_RCS           = 1;        % RCS bounce number 
k_ringing       = -20;      % max ringing specification (in dB)

% coupling and transmission filenames (s2p)
perform_post_design_spec = true;                            % want to calculate?
coupling_s2p_filename    = ".\data\coupling_20cm.s2p";      % file source for coupling
transmit_s2p_filename    = ".\data\transmission_1m.s2p";    % file source for transmit

% physical parameters
c = 3e8; % speed of light

%% 
% Pulse I/O Simulation Parameters
Fs  = (2*f_high)*10;            % sampling rate - 10x oversampling
Ts  = 1/Fs;                     % sampling period
N   = 1024;                     % number of points in t, f vectors
t   = (0:(N-1)).*Ts;            % time vector
f   = linspace(-Fs/2,Fs/2,N);   % frequency vector

% create an example input pulse based on f_lo, f_hi, and a FIR filter
BW_filter   = fir1(50, [f_low f_high]./(Fs/2),'bandpass');                  % FIR filter with N=50 at f_low,f_high
extra_delay = 4e-9;                                                         % extra time to delay pulse, for plotting consistency
H_filter    = freqz(BW_filter,1,f,Fs).*exp(-1j.*2.*pi.*f.*extra_delay);     % impulse response (f-domain)
x_t         = real(ifft(ifftshift(H_filter)));                              % IFFT for simulated input pulse
x_t         = x_t.*(V_pp_t/(max(x_t)-min(x_t)));                            % normalize to V_pp_t
X_f         = fftshift(fft(x_t));

% derived parameters
f_c     = (f_low + f_high)/2;   % center frequency
l_c     = c/f_c;                % center wavelength
G       = 10^(G/10);            % convert to linear
G_amp   = 10^(G_amp/10);        % convert to linear
o       = 10^(o/10);            % convert to linear
P_1dB   = 10^((P_1dB-30)/10);   % convert to Watts

%% Section 3 - GPR Antenna Coupling Specifications
% in this section, simulate a coupling signal and pass it through
% specifications to see if it violates

% calculate Vpp specification
max_coupling_Vpp = 2*sqrt(2)*sqrt(P_1dB*Z0/G_amp);

% calculate slope specification
max_coupling_slope = (k_coupling*V_pp_t/coupling_Zscore/o_j) * ...
    (G*l_c*sqrt(er)*sqrt(o))/((sqrt(er)+1)^2*(4*pi)^(3/2)*(h+d)^2) * ...
    exp(-alpha*d)/sqrt(2) ...
    *1e3/1e9;

if(perform_post_design_spec)
    % get simulated coupling signal
    Coupling_H  = get_tf_from_s2p(coupling_s2p_filename,f); % transfer function
    C_f         = Coupling_H.*X_f;                          % frequency-domain coupling signal
    c_t         = real(ifft(ifftshift(C_f)));               % time-domain coupling signal
    d_c_t       = gradient(c_t,Ts) .*1e3/1e9;               % coupling slope (mV/ns)
    time_gate   = (t-extra_delay) > (2*h/c);

    % find peaks to find vpp
    [peaks, locs] = findpeaks(abs(c_t),'MinPeakProminence',max(c_t)/10);
    [coupling_Vpp,max_ind] = max(abs(diff(c_t(locs))));   % find the maximum peak-to-peak voltage difference

    % find max coupling slope
    coupling_slope = max(abs(d_c_t(time_gate)));
    
    % does the antenna pass specification?
    passed_coupling_vpp_spec   = double(coupling_Vpp   <= max_coupling_Vpp);
    passed_coupling_slope_spec = double(coupling_slope <= max_coupling_slope);

    % send everything to a plot function
    plot_fig_1(t, c_t, locs, max_ind, d_c_t, time_gate, max_coupling_slope);
else
    passed_coupling_vpp_spec   = 2;
    passed_coupling_slope_spec = 2; 
end

%% Section 4 - GPR Antenna Distortion Requirements
% in this section, simulate a transmitted signal and pass it through
% specifications to see if it violates

% pulse integrity specification
if(perform_post_design_spec)
    Transmit_H  = get_tf_from_s2p(transmit_s2p_filename,f); % get voltage T.F.
    Y_f         = Transmit_H.*X_f;                          % calculate output spectrum
    y_t         = real(ifft(ifftshift(Y_f)));               % calculate output signal

    % find the cross-correlation. abs() taken to ensure flipping is OK.
    chi = max(abs(xcorr(x_t,y_t,'coeff')));

    % did it pass the distortion spec?
    passed_distortion_spec = double(chi > k_distortion);
    
    % plot the figure (in external function to save space)
    plot_fig_2(t,x_t,y_t,chi);
else
    passed_distortion_spec = 2;
end

% residual ringing specification
if(perform_post_design_spec)

    % find the FWTM (full width at a tenth-max) of input
    fwtm        = calculate_fwtm(t,x_t);

    [max_mag,max_ind] = max(abs(y_t));          % maximum of abs(y(t)) and argmax ind
    t_0         = t(max_ind) + fwtm/2;          % t_0 - time at which a perfectly transmitted pulse would approx. end
    R           = 20.*log10(abs(y_t)./max_mag); % Residual
    R(t < t_0)  = 0;                            % step function
    
    % we want to ignore the first falling edge, as it is a false alarm.
    % flip the residual and findpeaks() to only search after first "peak"
    R_flipped       = -R;
    [~,R_locs]      = findpeaks(R_flipped,'minpeakprominence',3);
    ind_first_pk    = R_locs(1);

    %does R pass specification?
    max_ringing_after_first_pk = max(R(ind_first_pk:end));
    passed_ringing_spec = double(max_ringing_after_first_pk < k_ringing);
    
    plot_fig_3(t,t_0,R,y_t,k_ringing);

else
    passed_ringing_spec = 2;
end

%% Section 5: GPR Antenna Gain, Beamwidth Requirements
% this section will generate a plot with horizontal resolution vs.
% beamwidth/gain
Nm      = 1024;                 % number of different "beamwidths"
m       = linspace(0.5,5,Nm);   % m used as parameter to sinc(m*theta)
Nt      = 1024;                 % number of theta (for simulation)
theta   = linspace(0,pi/2,Nt);  % vector from 0 to pi/2
G_0     = 4;                    % G0 - sets the absolute value of gain curve, contains efficiency- check Fig.2.

HPBWs   = zeros(Nm,1);  % contains the HPBWs
F       = zeros(Nt,Nm); % contains Nm pattern functions
tol     = 0.03;         % tolerance to zeroing out all but main beam

is_lt_noise_rms = zeros(Nm,1);          % logical vector used in finding theta_m
theta_ms        = zeros(Nm,1);          % stores values of theta_m
xa              = linspace(0,2+(2*h),512); % "search" vector to find theta_m 

% iterate over m to define the different gain patterns
for ii = 1:Nm
    mi  = m(ii);        
    
    % generate the pattern p (normalized by power integral)
    p   = abs(sinc(mi.*theta))...
        ./(2*pi.*trapz(theta, abs(sinc(mi.*theta)).*sin(theta)));
    
    % zero out everything but main beam to avoid unrealistic "choppiness"
    first_min_ind = find(p<tol,1,'first');
    p(first_min_ind:end)=1e-10.*max(p);
    
    % save p to the pattern matrix
    F(:,ii) = p;

    % calculate the HPBW of p (external function)
    HPBWs(ii) = find_hpbw(10*log10(p),theta);
end

G = F.*G_0; % un-normalize F to contain gain patterns

% sweep through xa to find the last theta_m
for ii = 1:numel(xa)
    % find theta_a and theta_g (angle in air, angle in ground) - ray based
    [theta_a,theta_g,~] = find_GPR_transmission_angles(er,h,xa(ii),0, 0, 0, d);

    % find the index which corresponds to theta_a (angle in air)
    [~,theta_ind] = min(abs(theta-theta_a));
    
    % find ground power reflectivity and transmissivity at angle theta (TE
    % example)
    R = ((cos(theta_a) - sqrt(er)*sqrt(1-(1/sqrt(er)*sin(theta_a) )))...
        /(cos(theta_a) + sqrt(er)*sqrt(1-(1/sqrt(er)*sin(theta_a))^2 )))^2;
    T = 1-R;

    % iterate over the possible patterns, finding the received voltage
    % using (6)
    V_rx_rms = zeros(Nm,1);
    for kk = 1:Nm
        Gi = G(theta_ind,kk);

        V_rx_rms(kk) = V_pp_t/2/sqrt(2)*(Gi*T*l_c*sqrt(o))...
            /((4*pi)^(3/2)*(h*sec(theta_a)+d*sec(theta_g))^2)...
            *exp(-alpha*d*sec(theta_g));
    end
    
    is_gt_noise_rms = V_rx_rms > V_noise_rms;               % logical vector-is V_rx_rms greater than noise?
    is_lt_noise_rms = is_lt_noise_rms | [~is_gt_noise_rms]; % has there been a less-than-noise value previously?

    % if greater than noise and have never experienced a less-than-noise event, save theta_g as theta_m.       
    theta_ms(is_gt_noise_rms & [~is_lt_noise_rms]) = theta_g; 
end

% calculate horizontal resolution based on theta_m
HRs = l_c./2./sqrt(er)./sin(theta_ms).*1e2;

[min_hpbw, max_hpbw] = find_optimal_beamwidths(HPBWs, HRs);

% plot the figure (in external function to save space)
plot_fig_4(HPBWs,HRs,max(G,[],1));

%% Section 6: Antenna RCS Specification

R = 1/(8*pi*h)*abs((1-sqrt(er))/(1+sqrt(er)));
oant_max = -20*log10(R) ...
         + 20/m_RCS*log10(k_RCS*(l_c*sqrt(er*o)*exp(-alpha*d))/((sqrt(er)+1)^2*(h+d)^2)) ...
         - 20.94/m_RCS;
     
%% Section 7: Size, Weight

maximum_width = l_c/4 - dx;
maximum_weight = (p_max-m_sys)/N_ants;
     
%% Generate console output with specifications

spec_pass_msgs = {"No Pass","Pass","N/A"};
fprintf("\n_________________________ Post-Design Specifications _________________________\n");
fprintf("Antenna Coupling Vpp:\n");
fprintf("\t Maximum: %8.2f V      | Calculated: %8.2f V     | Passed? %s\n", max_coupling_Vpp, coupling_Vpp, spec_pass_msgs{passed_coupling_vpp_spec+1});
fprintf("Antenna Coupling Slope:\n");
fprintf("\t Maximum: %8.2f  mV/ns | Calculated: %8.2f mV/ns | Passed? %s\n", max_coupling_slope, coupling_slope, spec_pass_msgs{passed_coupling_slope_spec+1});
fprintf("Pulse Integrity/Distortion:\n");
fprintf("\t Minimum: %8.2f        | Calculated: %8.3f       | Passed? %s\n", k_distortion, chi, spec_pass_msgs{passed_distortion_spec+1});
fprintf("Residual Ringing:\n");
fprintf("\t Minimum: %8.2f dB     | Calculated: %8.3f dB    | Passed? %s\n", k_ringing, max_ringing_after_first_pk, spec_pass_msgs{passed_ringing_spec+1});

fprintf("\n__________________________ Pre-Design Specifications _________________________\n");
fprintf("Beamwidth: See Figure 2. Rec. Range: [%8.2f, %8.2f] deg.\n",min_hpbw*180/pi,max_hpbw*180/pi);
fprintf("RCS (check if problematic):       %8.2f dBsm\n", oant_max);
fprintf("Size:      %8.2f cm\n", maximum_width*1e2);
fprintf("Weight:    %8.2f g\n\n", maximum_weight.*1e3);


%% custom functions - plotting figures
function plot_fig_1(t, c_t, locs, max_ind, d_c_t, time_gate, max_coupling_slope)
    first_nongated_time_ind = find(time_gate,1,'first');
    
    figure(); set(gcf,'position',[453   287   718   682]);
    subplot(211);
    hold on; grid on;
    plot(t.*1e9,c_t.*1e3,'k');
    plot(t(locs([max_ind max_ind+1])).*1e9, c_t(locs([max_ind max_ind+1])).*1e3,'rx');
    xlabel("Time (ns)");
    ylabel("Voltage (mV)");
    title("Peak-to-Peak Specification");
    legend("Coupling signal","Max. Vpp peaks");
    aa=gca;
    aa.GridColor=[0.1 0.1 0.1];
    aa.GridAlpha=0.3;
    
    subplot(212);
    hold on; grid on;
    plot(t(time_gate).*1e9,d_c_t(time_gate),'k');
    plot([min(t) max(t)].*1e9,max_coupling_slope.*[1 1],'r-');
    plot([min(t) max(t)].*1e9,max_coupling_slope.*[-1 -1],'r-');
    plot(t.*1e9,d_c_t,'k--');
    plot(t(first_nongated_time_ind).*1e9.*[1 1], max_coupling_slope.*[-1 1],'r-');
    aa=gca;
    aa.GridColor=[0.1 0.1 0.1];
    aa.GridAlpha=0.3;
    ylim(aa.YLim.*1.5);
    xlabel("Time (ns)");
    ylabel("Slope (mV/ns)");
    title("Slope specification");
    legend("Coupling slope (gated)","Maximum Slope Threshold");
end
function plot_fig_2(t,x_t,y_t,chi)
    figure(); set(gcf,'position',[312   488   806   420]);
    hold on; grid on; 
    plot(t.*1e9,normalize(x_t),'k');
    plot(t.*1e9,normalize(y_t),'r');
    legend("Input pulse","Output pulse ($$\chi$$ = "+num2str(chi)+")");
    xlabel("Time, ns");
    ylabel("Mag., Norm.");
    title("Distortion Specification");
    aa=gca;
    aa.GridColor=[0.1 0.1 0.1];
    aa.GridAlpha=0.3;
end


function plot_fig_3(t,t_0,R,y_t,k)
    [~,max_ind]=max(abs(y_t));
    t_max = t(max_ind);
    fwtm_2 = t_0 - t_max;
    
    figure(); 
    hold on; grid on;
    plot((t-t_0).*1e9,R,'k');
    plot([-2 3],k.*[1 1],'k--')
    p1 = patch([-2 3 3 -2],[-20 -20 0 0], [1 0 0]);
    p2 = patch([-2 3 3 -2],[-20 -20 -60 -60], [0 1 0]);
    
    p1.FaceAlpha=0.3;
    p2.FaceAlpha=0.3;
    xlabel("$$t-t_{0}$$ (ns)");
    ylabel("R (dB)");
    ylim([-60 0])
    xlim([-2 3])
    aa=gca;
    aa.GridColor=[0.1 0.1 0.1];
    aa.GridAlpha=0.3;
    

    figure();
    hold on; grid on;
    plot((t-t_0).*1e9,normalize(y_t),'k');
    plot(-fwtm_2.*[1 1].*1e9,[-1 1],'k--');
    plot(0.*[1 1].*1e9,[-1 1],'k--');
    xlabel("$$t-t_{0}$$ (ns)");
    xlim([-2 2])
    aa=gca;
    aa.GridColor=[0.1 0.1 0.1];
    aa.GridAlpha=0.3;
    
end

function plot_fig_4(HPBWs,HRs,maxGains)
    min_HR  = min(HRs);

    gainfit = fit(HPBWs.*180./pi,maxGains','power1');
    
    figure(); set(gcf,'position',[100, 421, 1000, 455]);
    subplot(121);
    hold on; grid on;
    plot(flip(HPBWs).*180./pi, flip(HRs),'k');
    min_HR_area = area([10,120],(min_HR+0.1*min_HR).*[1 1],min_HR);
    min_HR_area.FaceColor=[0,0,1];
    min_HR_area.FaceAlpha=0.3;
    xlim([10 120]);
    ylim([0 max(HRs)*3/4])
    xlabel("Beamwidth ($$^{\circ}$$)");
    ylabel("Resolution (cm)");
    aa=gca;
    aa.GridColor=[0.1 0.1 0.1];
    aa.GridAlpha=0.3;
    title("Horizontal Resolutions");
    
    subplot(122);
    hold on; grid on;
    plot(HPBWs.*180./pi,10*log10(maxGains),'k');
    plot(HPBWs.*180./pi,10.*log10(feval(gainfit,HPBWs.*180./pi)),'b');
    xlabel("Beamwidth ($$^{\circ}$$)");
    ylabel("Boresight Gain (dB)");
    legend('Numerical','Fit');
    aa=gca;
    aa.GridColor=[0.1 0.1 0.1];
    aa.GridAlpha=0.3;
    title("Gain vs. Beamwidth (Adjust $$G_0$$)");
end
