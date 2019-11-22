clear all; close all; clc;

addpath([pwd '\Functions\'])


%%
data_init = load('initial_data');
data_t = load('MDF_N120');


dk = 4*pi/(852*1e-9)/52;
k_arr_exp = 0:dk:(length(data_init.N120_after) - 1)*dk;
k_arr_exp = k_arr_exp - mean(k_arr_exp);

nk_init_exp = 120*data_init.N120_after/trapz(k_arr_exp, data_init.N120_after);

nk_exp = 120*data_t.MDF_until30./trapz(k_arr_exp, data_t.MDF_until30')';

t_arr_exp = (0:size(nk_exp,1)-1) * 1e-3;


%% Simulation parameters

N           = 2^7;
M           = 2^6;

Nsteps      = 120;
stepOrder   = 2;

xmax_si     = 70*1e-6;
kmax_si     = 3*1e7;
tmax_si     = 12*1e-3;


%% Time-of-flight parameters

xmax_tof    = 700*1e-6;
tmax_tof    = 30*1e-3;


%% Experimental units

m_si        = 87*1.6605402e-27;
hbar_si     = 1.054571726e-34;
kB_si       = 1.38065e-23;
as_si       = 5.2e-9;

Natoms      = 120;
T_si        = 100*1e-9; 
omega_si    = 2*pi*83.3; % [s^-1]
omegaT_si   = 2*pi*31*1e3;
k_Bragg_si  = 4*pi/(852*1e-9); % [m^-1]

g1D_si      = 2*hbar_si*omegaT_si*as_si;
c_si        = g1D_si*m_si/hbar_si^2;


%% create arrays
x_array     = linspace(-xmax_si,xmax_si,M);
k_array     = linspace(-kmax_si,kmax_si,N);
t_array     = linspace(0, tmax_si, Nsteps+1);

dx          = x_array(2) - x_array(1);
x_array_tof = -xmax_tof:dx:xmax_tof;

dt          = t_array(2) - t_array(1);

%% Initialize state 
Vx      = @(t,x) 0.5*m_si*omega_si^2*x.^2;
dmudx   = @(t,x) -m_si*omega_si^2*x;

couplings = { @(t,x) 0              , @(t,x) c_si;
              []                    , []            ;
              @(t,x) dmudx(t,x)     , []            };

LLS         = LiebLinigerSolver_SI_test(omega_si, x_array, k_array, couplings);
mu0_fit     = LLS.fitAtomnumber(T_si, Vx, Natoms, true);

theta_init  = LLS.calcThermalState(T_si);

delta_k     = (k_array(2)-k_array(1));
theta_bragg = applyBraggPulse_test(theta_init, k_Bragg_si, delta_k, 0.1, 2.5, LLS);
% theta_bragg = applyBraggPulse(theta_init, k_Bragg_si, delta_k, LLS);


%% Initialize solver for TOF
ramp_fast   = @(t) heaviside( 1.8 - 1e3*t - eps ).*(1 - 0.5e3*t);
ramp_slow   = @(t) heaviside( 4.5 - 1e3*t - eps ).*(1 - 0.2e3*t);

dr_fast     = @(t) 0.5e3*heaviside( 1.8 - 1e3*t - eps );
dr_slow     = @(t) 0.2e3*heaviside( 4.5 - 1e3*t - eps );

couplings_tof1 = {  @(t,x) 0              , @(t,x) 0      ;
                    []                    , []            ;
                    []                    , []            };
                
couplings_tof2 = {  @(t,x) 0              , @(t,x)  ramp_fast(t)*c_si;
                    []                    , @(t,x)  dr_fast(t)*c_si;
                    []                    , []            };
                
couplings_tof3 = {  @(t,x) 0              , @(t,x)  ramp_slow(t)*c_si;
                    []                    , @(t,x)  dr_slow(t)*c_si;
                    []                    , []            };


                
                
LLS_tof     = LiebLinigerSolver_SI_test(1e5*omega_si, x_array_tof, k_array, couplings);
c_tof       = { couplings_tof1, couplings_tof2, couplings_tof3  };


skip = 2.5;
t_array_tof{1} = 0:skip*dt:tmax_tof;
t_array_tof{2} = [ 0:dt:1.8e-3  (1.8e-3 + skip*dt):skip*dt:tmax_tof ];
t_array_tof{3} = [ 0:dt:4.5e-3  (4.5e-3 + skip*dt):skip*dt:tmax_tof ];


figure
hold on
box on
plot(t_array_tof{2}*1e3, ramp_fast(t_array_tof{2}))
plot(t_array_tof{3}*1e3, ramp_slow(t_array_tof{3}))


%% Solve dynamics in cradle
theta_t     = LLS.propagateTheta(theta_bragg, t_array);


%% Calculate TOF profile at every 1 ms
nk_rapid = zeros(13, N);
nx_tof = zeros(13, length(x_array_tof), length(c_tof));
theta_tof = cell(13, length(c_tof));

for i = 1:13
    index =  1 + (i-1)*10;
    theta = theta_t{ index };
    
    % Calculate momentum space distribution given by rapidities
    nk_rapid(i,:) = LLS.calcMomentumDistr(theta, t_array(index));
    
    % Interpolate theta to TOF-grid
    pad = (length(x_array_tof) - M)/2;
    theta_temp = permute(double(theta), [1 5 2 3 4]);
    theta_temp = [ zeros(N,pad) , theta_temp , zeros(N,pad) ];
    theta_temp = permute( theta_temp, [1 5 3 4 2] );
    
    theta = GHDtensor( theta_temp );
    
    for j = 1:length(c_tof)
        LLS_tof.setCouplings( c_tof{j} );
        theta_free  = LLS_tof.propagateTheta(theta, t_array_tof{j});
        
        theta_tof{i,j} = theta_free{ end };
        nx_tof(i,:,j) = LLS_tof.calcCharges(theta_free{ end }, 0, t_array_tof{j}(end));
    end

end


%% Plot results

k_array_tof = m_si*x_array_tof/tmax_tof/hbar_si;
nk_tof      = nx_tof / (m_si/hbar_si/tmax_tof);


for i = 1:13
    figure
    hold on
    box on
    
    for j = 1:length(c_tof)
        plot(k_array_tof, nk_tof(i,:,j))
    end
    
    plot(k_array, nk_rapid(i,:), 'k--')
    plot(k_arr_exp, nk_exp(i,:), 'k')
    
    xlabel('k')
    ylabel('n(k)')
    
    title([ 't = ' num2str(i-1) 'ms'])
end



%% ------------------------------------------------
function theta_bragg = applyBraggPulse(theta, k_Bragg, dk, LLS)
    rho         = LLS.transform2rho(theta);
    rho_shift   = circshift( double(rho)/2, [round(k_Bragg/dk), 0]);
    rho_bragg   = GHDtensor( rho_shift + flipud(rho_shift) );
    theta_bragg = LLS.transform2theta(rho_bragg);
end


function theta_bragg = applyBraggPulse_test(theta, k_Bragg, dk, subamp, subwidth, LLS)
    rho         = LLS.transform2rho(theta);
    
    % Create central blob
    rho_center  = squeeze( subamp*double(rho) );
    center_amp  = sum(sum(rho_center));
    
    x_grid      = 1:size(rho_center, 2);
    y_grid      = 1:size(rho_center, 1);
    
    x_grid      = x_grid - mean(x_grid);
    y_grid      = y_grid - mean(y_grid);
    
    [X,Y]       = meshgrid( x_grid, y_grid' );
    [XX,YY]     = meshgrid( x_grid, subwidth*y_grid' );
    
    rho_center  = interp2( XX,YY,rho_center,X,Y, 'spline');
    rho_center  = rho_center * center_amp/sum(sum(rho_center));
    rho_center  = permute(rho_center, [1 5 3 4 2]);

    % Create bragg peaks
    rho_peak    = circshift( (1-subamp)*double(rho)/2, [round(k_Bragg/dk), 0]);
     
    rho_bragg   = GHDtensor( rho_peak + flipud(rho_peak) + rho_center );
    theta_bragg = LLS.transform2theta(rho_bragg);
end

