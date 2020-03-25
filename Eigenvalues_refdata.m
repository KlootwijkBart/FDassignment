% Load reference data
load matlab.mat

c      = 2.0569;	      % mean aerodynamic cord [m]
V0     = 76.87;          % true airspeed in the stationary flight condition [m/sec]

% Flight Data
FD_t   = flightdata.time.data;                 % time [s]
FD_h   = flightdata.Dadc1_alt.data*0.3048;     % altitude [m]
FD_cas = flightdata.Dadc1_cas.data*0.51444444; % calibrated air speed [m/s]
FD_tas = flightdata.Dadc1_tas.data*0.51444444; % true air speed [m/s]  
FD_TAT = flightdata.Dadc1_tat.data+273.15;     % total air temperature [K]

% Angle of attack
FD_aoa = flightdata.vane_AOA.data;             % angle of attack [deg]

% Euler angles
FD_th  = flightdata.Ahrs1_Pitch.data;          % pitch angle [deg]
FD_phi = flightdata.Ahrs1_Roll.data;           % roll angle [deg]  
   
% Rotation rates
FD_q   = flightdata.Ahrs1_bPitchRate.data;     % pitch rate [deg/s]
FD_p   = flightdata.Ahrs1_bRollRate.data;      % roll rate [deg/s]
FD_r   = flightdata.Ahrs1_bYawRate.data;       % yaw rate [deg/s]

% Control surface deflection
FD_de  = flightdata.delta_e.data;              % elevator deflection [deg]
FD_da  = flightdata.delta_a.data;              % aileron deflection [deg]
FD_dr  = flightdata.delta_r.data;              % rudder deflection [deg]

% Stick input
FD_fe  = flightdata.column_fe.data;            % stick force [N]
% FD_Se  = flightdata.column_Se.data;            % stick deflection [deg]

% indices for start and end of different eigenmodes, in the order shown
% below
% 1: Phugoid, 2: Short period, 3: A-periodic roll, 4: Dutch roll, 5: Dutch
% roll damped, 6: Spiral
idxstart = [32430,36340,35507,37210,37670,39200]-89;
idxend   = [34560,36480,35621,37380,37711,39411]-89;

% ask user input 
eigenmode = input(['',...
                   '\n1: Phugoid',...
                   '\n2: Short period',...
                   '\n3: A-periodic roll',...
                   '\n4: Dutch roll',...
                   '\n5: Dutch roll damped',...
                   '\n6: Spiral',...
                   '\n',...
                   '\nWhich eigenmode is to be simulated? ']);
%% Phugoid
if eigenmode == 1

idxstart = idxstart(eigenmode);
idxend   = idxend(eigenmode);

t_array = FD_t(idxstart:idxend);
q_array = FD_th(idxstart:idxend);

clf()
plot(t_array, q_array)

% Find peaks
[q_peaks_max, t_peaks_max] = findpeaks(q_array);
[q_peaks_min, t_peaks_min] = findpeaks(-1*q_array);
q_peaks_min = -1*q_peaks_min;

% Find middle of oscillation q_zero
q_diff = q_peaks_max(1) - interp1(t_peaks_min, q_peaks_min, t_peaks_max(1), 'linear');
q_zero = q_peaks_max(1) - q_diff/2;

% Normalize peaks around q_zero
q_peaks_max = q_peaks_max - q_zero;
q_peaks_min = q_peaks_min - q_zero;
q_peaks_min = -1*q_peaks_min; % to get array of positive amplitudes

% Concatenate and sort peaks to make array of amplitudes
q_peaks = cat(1, q_peaks_max, q_peaks_min);
t_peaks = cat(1, t_peaks_max, t_peaks_min);

[t_peaks, sortIdx] = sort(t_peaks);
q_peaks = q_peaks(sortIdx);

t_peaks = t_peaks/10 + t_array(1);

T_half = interp1(q_peaks, t_peaks, 0.5*q_peaks(1), 'linear', 'extrap') - t_peaks(1);
lambda_abs = log(0.5)*c/(T_half*V0);

Period_array = zeros(length(t_peaks)-1);
for i = 1:length(t_peaks)-1
    Period_i = t_peaks(i+1) - t_peaks(i);
    Period_array(i) = Period_i;
end

delta_array = zeros(length(t_peaks)-1);
for i = 1:length(q_peaks_max)-1
   delta_i = q_peaks_max(i+1)/q_peaks_max(i);
   delta_array(i) = delta_i; 
end

% Calculation period & imaginary part eta
Period = mean(Period_array);
Period = Period(1)*2;
eta = 2*pi/(Period);

% Logarithmic decrement delta & real part xi
delta = mean(delta_array);
delta = log(delta(1));
xi = delta/(Period);

lambda = xi + eta*1i
end

%% Short period
if eigenmode == 2
    idxstart = idxstart(eigenmode);
    idxend   = idxend(eigenmode);

    t_array = FD_t(idxstart:idxend);
    th_array = FD_th(idxstart:idxend);
    de_array = FD_de(idxstart:idxend);

    clf()
    plot(t_array, [th_array de_array])
end

%% Aperiodic roll
if eigenmode == 3
    idxstart = idxstart(eigenmode);
    idxend   = idxend(eigenmode);
    
    t_array = FD_t(idxstart:idxend);
    p_array = FD_p(idxstart:idxend);
    
    clf()
    plot(t_array, p_array)
    
    p_peaks = findpeaks(p_array);
    p_0 = p_peaks(1);
    p_first = p_array(1);
    p_halfamp = (p_0-p_first)/2 + p_first;
    
    p = p_first;
    i = 1;
    while p < p_halfamp
        i = i+1;
        p = p_array(i);
    end
    
    T_half = t_array(i)-t_array(1);
    lambda = log(0.5)/(T_half)
end

%% Dutch roll
if eigenmode == 4
    idxstart = idxstart(eigenmode);
    idxend   = idxend(eigenmode);
    
    t_array = FD_t(idxstart:idxend);
    r_array = FD_r(idxstart:idxend);
    dr_array = FD_dr(idxstart:idxend);
    
    % Find peaks
    [r_peaks_max, t_peaks_max] = findpeaks(r_array);
    [r_peaks_min, t_peaks_min] = findpeaks(-1*r_array);
    r_peaks_min = -1*r_peaks_min;

    % Find middle of oscillation q_zero
    r_diff = r_peaks_max(1) - interp1(t_peaks_min, r_peaks_min, t_peaks_max(1), 'linear');
    r_zero = r_peaks_max(1) - r_diff/2;

    % Normalize peaks around q_zero
    r_peaks_max = r_peaks_max - r_zero;
    r_peaks_min = r_peaks_min - r_zero;
    r_peaks_min = -1*r_peaks_min; % to get array of positive amplitudes

    % Concatenate and sort peaks to make array of amplitudes
    r_peaks = cat(1, r_peaks_max, r_peaks_min);
    t_peaks = cat(1, t_peaks_max, t_peaks_min);

    [t_peaks, sortIdx] = sort(t_peaks);
    r_peaks = r_peaks(sortIdx);

    t_peaks = t_peaks/10 + t_array(1);

    T_half = interp1(r_peaks, t_peaks, 0.5*r_peaks(1), 'linear', 'extrap') - t_peaks(1);
    lambda_abs = log(0.5)*c/(T_half*V0);

    Period_array = zeros(length(t_peaks)-1);
    for i = 1:length(t_peaks)-1
        Period_i = t_peaks(i+1) - t_peaks(i);
        Period_array(i) = Period_i;
    end

    delta_array = zeros(length(t_peaks)-1);
    for i = 1:length(r_peaks_max)-1
       delta_i = r_peaks_max(i+1)/r_peaks_max(i);
       delta_array(i) = delta_i; 
    end

    % Calculation period & imaginary part eta
    Period = mean(Period_array);
    Period = Period(1)*2;
    eta = 2*pi/(Period);

    % Logarithmic decrement delta & real part xi
    delta = mean(delta_array);
    delta = log(delta(1));
    xi = delta/(Period);

    lambda = xi + eta*1i
end

%% Spiral
if eigenmode == 6
    idxstart = idxstart(eigenmode);
    idxend   = idxend(eigenmode);
    
    t_array = FD_t(idxstart:idxend);
    p_array = FD_p(idxstart:idxend);
    
    clf()
    plot(t_array, p_array)
    
    % Find peaks
    [p_peaks_max, t_peaks_max] = findpeaks(p_array);
    [p_peaks_min, t_peaks_min] = findpeaks(-1*p_array);
    p_peaks_min = -1*p_peaks_min;

    % Find middle of oscillation q_zero
    p_diff = p_peaks_max(1) - interp1(t_peaks_max(1:2), p_peaks_max(1:2), t_peaks_min(1), 'linear');
    p_zero = p_peaks_max(1) - p_diff/2;

    % Normalize peaks around q_zero
    p_peaks_max = p_peaks_max - p_zero;
    p_peaks_min = p_peaks_min - p_zero;
    p_peaks_min = -1*p_peaks_min; % to get array of positive amplitudes

    % Concatenate and sort peaks to make array of amplitudes
    p_peaks = cat(1, p_peaks_max, p_peaks_min);
    t_peaks = cat(1, t_peaks_max, t_peaks_min);

    [t_peaks, sortIdx] = sort(t_peaks);
    p_peaks = p_peaks(sortIdx);

    t_peaks = t_peaks/10 + t_array(1);

    % Calculation period & imaginary part eta
    Period = t_peaks_max(2)-t_peaks_max(1);
    eta = 2*pi/(Period);

    % Logarithmic decrement delta & real part xi
    delta = log(p_peaks_max(2)/p_peaks_max(1));
    xi = delta/(Period);

    lambda = xi + eta*1i
end