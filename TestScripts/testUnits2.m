clear all; close all;

addpath(['..' filesep 'Functions' filesep])

% Reproduces the results found in the paper:
% "Hydrodynamics of the interacting Bose gas in the Quantum Newton Cradle
% setup" (http://arxiv.org/abs/1711.00873)
% NOTE: Newton craddle paper uses velocity as rapidity, whereas code uses
% momentum as rapidity! 

% NOTE: There is a factor 2 difference when transforming from k to v due to
% m = 1/2 in my units ...


%% Simulation parameters

N           = 2^7;
M           = 2^6;

Nsteps      = 200;
stepOrder   = 2;

xmax_si     = 70*1e-6;
vmax_si     = 35*1e-3;
tmax_si     = 13*1e-3;

Options.autoDerivCoup = false;

%% Experimental units

m_si        = 87*1.6605402e-27;
hbar_si     = 1.054571726e-34;
kB_si       = 1.38065e-23;
as_si       = 5.2e-9;

Natoms      = 350;
c_si        = 12.3*1e6; % [ m^-1]
g1D_si      = c_si*hbar_si^2/m_si;
T_si        = 57*1e-9; 
omega_si    = 2*pi/0.013; % [s^-1]
omegaT_si   = g1D_si/(2*hbar_si*as_si);
k_Bragg_si  = 2*12.5*1e6; % [m^-1] ADDED FACTOR 2


%% create arrays
x_array     = linspace(-xmax_si,xmax_si,M);
v_array     = linspace(-vmax_si,vmax_si,N);
k_array     = m_si/hbar_si * v_array;
t_array     = linspace(0, tmax_si, Nsteps);


%% Initialize state and solve dynamics for Harmonic Oscillator
Vx_HO      = @(t,x) 0.5*m_si*omega_si^2 *x.^2;
dmux_HO    = @(t,x) -m_si*omega_si^2 *x;

couplings  = {  @(t,x) 0            , @(t,x) c_si   ;    % couplings
                []                  , []            ;    % d/dt
                dmux_HO             , []            };   % d/dx


LLS_HO      = LiebLinigerSolver_SI(omegaT_si, x_array, k_array, couplings, Options);
mu0_fit     = LLS_HO.fitAtomnumber(T_si, Vx_HO, Natoms, true);


theta_init  = LLS_HO.calcThermalState(T_si);

delta_k     = (k_array(2)-k_array(1));
theta_bragg = applyBraggPulse(theta_init, k_Bragg_si, delta_k, LLS_HO);

theta_t_HO  = LLS_HO.propagateTheta(theta_bragg, t_array);


%% Initialize state and solve dynamics for Anharmonic Oscillator
l_si        = 100*1e-6;
Vx_AHO      = @(t,x) 1/pi^2 *omega_si^2 *l_si^2 *m_si* (1 - cos(pi*x/l_si));
dmux_AHO    = @(t,x) -l_si/pi*omega_si^2 *m_si* sin(pi*x/l_si);

couplings{3,1} = dmux_AHO;

LLS_AHO     = LiebLinigerSolver_SI(omegaT_si, x_array, k_array, couplings, Options);
mu0_fit     = LLS_AHO.fitAtomnumber(T_si, Vx_AHO, Natoms, true);

theta_init  = LLS_AHO.calcThermalState(T_si);

delta_k     = (k_array(2)-k_array(1));
theta_bragg = applyBraggPulse(theta_init, k_Bragg_si, delta_k, LLS_AHO);


theta_t_AHO = LLS_AHO.propagateTheta(theta_bragg, t_array);


%% Calculate quantities
n_t_HO      = LLS_HO.calcCharges(theta_t_HO, 0, t_array) * 1e-6; % [um^-1]
rho_t_HO    = LLS_HO.transform2rho( theta_t_HO, t_array); 

n_t_AHO     = LLS_AHO.calcCharges(theta_t_AHO, 0, t_array) * 1e-6; % [um^-1]
rho_t_AHO   = LLS_AHO.transform2rho( theta_t_AHO, t_array);


%% Plot results

figure
for i = 1:6
    index = ceil( (length(t_array)-1)/5*(i-1) + 1 );
    
    rho_i_HO = 2*(hbar_si/m_si)*rho_t_HO{index}* 1e9; % [ms/um^2]
    rho_i_AHO = 2*(hbar_si/m_si)*rho_t_AHO{index}* 1e9; % [ms/um^2]
    
    % plot Harmonic quasiparticle distribution
    subplot(3,6,i)
    imagesc(x_array*1e6, v_array*1e3,squeeze(double(rho_i_HO)))
    set(gca,'YDir','normal')
    colormap(hot)
    caxis([0 0.36])
    ylabel('v [um/ms]')
    set(gca,'xticklabel',[])
    title(['t = ' num2str(round(t_array(index)/(2*pi/omega_si) ,2 )) 'T'])
    
    % plot Anharmonic quasiparticle distribution
    subplot(3,6,i+6)
    imagesc(x_array*1e6, v_array*1e3,squeeze(double(rho_i_AHO)))
    set(gca,'YDir','normal')
    colormap(hot)
    caxis([0 0.36])
    ylabel('v [um/ms]')
    set(gca,'xticklabel',[])

    % plot density profiles
    subplot(3,6,i+12)
    hold on
    box on
    plot(x_array*1e6, n_t_HO(:,index))
    plot(x_array*1e6, n_t_AHO(:,index))
    xlim([ -xmax_si*1e6, xmax_si*1e6])
    ylim([0 12])
    xlabel('x [um]')
    ylabel('n [um^-1]')
    xticks([ -60, -30 , 0 30, 60])
end




%% ------------------------------------------------
function theta_bragg = applyBraggPulse(theta, k_Bragg, dk, LLS)
    rho         = LLS.transform2rho(theta);
    rho_shift   = circshift( double(rho)/2, [round(k_Bragg/dk), 0]);
    rho_bragg   = GHDtensor( rho_shift + flipud(rho_shift) );
    theta_bragg = LLS.transform2theta(rho_bragg);
end