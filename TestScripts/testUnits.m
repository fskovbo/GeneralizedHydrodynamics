clear all; close all;

addpath('..\Functions\')

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

Nsteps      = 2*400;
stepOrder   = 2;

xmax_si     = 70*1e-6;
vmax_si     = 35*1e-3;
tmax_si     = 2*13*1e-3;

%% Experimental units

m_si        = 87*1.6605402e-27;
hbar_si     = 1.054571726e-34;
kB_si       = 1.38065e-23;
as_si       = 5.2e-9;

Natoms      = 350;
omega_si    = 2*pi/0.013; % [s^-1]
c_si        = 12.3*1e6; % [ m^-1]
g1D_si      = c_si*hbar_si^2/m_si;
T_si        = 57*1e-9; 
k_Bragg_si  = 2*12.5*1e6; % [m^-1] ADDED FACTOR 2

%% TBA units
Eg_si       = 0.5*m_si*g1D_si^2 /hbar_si^2;
Lg_si       = hbar_si^2 /(g1D_si*m_si);
t_si        = hbar_si/Eg_si;


%% create arrays
x_array_si  = linspace(-xmax_si,xmax_si,M);
v_array_si  = linspace(-vmax_si,vmax_si,N);
k_array_si  = m_si/hbar_si * v_array_si;
t_array_si  = linspace(0, tmax_si, Nsteps);


%% Convert to TBA units
x_array     = x_array_si/Lg_si;
k_array     = k_array_si*Lg_si;
t_array     = t_array_si/t_si;

T_tba       = T_si*kB_si/Eg_si;
k_Bragg_tba = k_Bragg_si*Lg_si;


%% Initialize state and solve dynamics for Harmonic Oscillator
Vx_HO      = @(t,x) 0.5*m_si*omega_si^2*(x*Lg_si).^2/Eg_si;

couplings.mu    = @(t,x) 0;
couplings.dmudx = @(t,x) -m_si*omega_si^2*(x*Lg_si^2)/Eg_si;
couplings.c     = @(t,x) c_si*Lg_si;
couplings.dcdt  = [];
couplings.dcdx  = [];
couplings.dmudt = [];


LLS_HO      = LiebLinigerSolver(x_array, k_array, couplings);
mu0_fit     = LLS_HO.fitAtomnumber(T_tba, Vx_HO, Natoms, true);

theta_init  = LLS_HO.calcThermalState(T_tba);

delta_k     = (k_array(2)-k_array(1));
theta_bragg = applyBraggPulse(theta_init, k_Bragg_tba, delta_k, LLS_HO);

theta_t_HO  = LLS_HO.propagateTheta(theta_bragg, t_array);


%% Initialize state and solve dynamics for Anharmonic Oscillator
l_si        = 100*1e-6;
Vx_AHO      = @(t,x) 1/pi^2 *omega_si^2 *l_si^2 *m_si* (1 - cos(pi*(x*Lg_si)/l_si))/Eg_si;
dVx_AHO     = @(t,x) l_si*Lg_si/pi*omega_si^2 *m_si* sin(pi*(x*Lg_si)/l_si)/Eg_si;

couplings.mu    = @(t,x) 0;
couplings.dmudx = @(t,x) -dVx_AHO(t,x);
couplings.c     = @(t,x) c_si*Lg_si;
couplings.dcdt  = [];
couplings.dcdx  = [];
couplings.dmudt = [];


LLS_AHO     = LiebLinigerSolver(x_array, k_array, couplings);
mu0_fit     = LLS_AHO.fitAtomnumber(T_tba, Vx_AHO, Natoms, true);

theta_init  = LLS_AHO.calcThermalState(T_tba);

delta_k     = (k_array(2)-k_array(1));
theta_bragg = applyBraggPulse(theta_init, k_Bragg_tba, delta_k, LLS_AHO);

theta_t_AHO = LLS_AHO.propagateTheta(theta_bragg, t_array);


%% Calculate quantities
n_t_HO      = LLS_HO.calcCharges(theta_t_HO, 0, t_array);
rho_t_HO    = LLS_HO.transform2rho( theta_t_HO, t_array);

n_t_AHO     = LLS_AHO.calcCharges(theta_t_AHO, 0, t_array);
rho_t_AHO   = LLS_AHO.transform2rho( theta_t_AHO, t_array);


%% Plot results
rho_t_si_HO = 2 * rho_t_HO/Lg_si/(m_si/hbar_si)*Lg_si * 1e9; % [ms/um^2] HAS AN EXTRA FACTOR 2!!!
n_t_si_HO   = n_t_HO/Lg_si *1e-6; % [um^-1]

rho_t_si_AHO= 2 * rho_t_AHO/Lg_si/(m_si/hbar_si)*Lg_si * 1e9; % [ms/um^2] HAS AN EXTRA FACTOR 2!!!
n_t_si_AHO  = n_t_AHO/Lg_si *1e-6; % [um^-1]


figure
for i = 1:6
    index = ceil( (length(t_array)-1)/5*(i-1) + 1 );
    
    % plot Harmonic quasiparticle distribution
    subplot(3,6,i)
    imagesc(x_array_si*1e6, v_array_si*1e3,squeeze(rho_t_si_HO(:,:,:,index)))
    set(gca,'YDir','normal')
    colormap(hot)
    caxis([0 0.36])
    ylabel('v [um/ms]')
    set(gca,'xticklabel',[])
    title(['t = ' num2str(round(t_array_si(index)/(2*pi/omega_si) ,2 )) 'T'])
    
    % plot Anharmonic quasiparticle distribution
    subplot(3,6,i+6)
    imagesc(x_array_si*1e6, v_array_si*1e3,squeeze(rho_t_si_AHO(:,:,:,index)))
    set(gca,'YDir','normal')
    colormap(hot)
    caxis([0 0.36])
    ylabel('v [um/ms]')
    set(gca,'xticklabel',[])

    % plot density profiles
    subplot(3,6,i+12)
    hold on
    box on
    plot(x_array_si*1e6, n_t_si_HO(:,index))
    plot(x_array_si*1e6, n_t_si_AHO(:,index))
    xlim([ -xmax_si*1e6, xmax_si*1e6])
    ylim([0 12])
    xlabel('x [um]')
    ylabel('n [um^-1]')
    xticks([ -60, -30 , 0 30, 60])
end



%% Run movie
figure
for i = 1:length(t_array)
    subplot(2,1,1)
    imagesc(x_array_si*1e6, v_array_si*1e3,squeeze(rho_t_si_AHO(:,:,:,i)))
    colormap(hot)
    caxis([0 0.36])
    ylabel('v [um/ms]')
    set(gca,'YDir','normal') 
    
    subplot(2,1,2)
    plot(x_array_si*1e6, n_t_si_AHO(:,i) )
    xlim([ -xmax_si*1e6, xmax_si*1e6])
    ylim([0 12])
    xlabel('x [um]')
    ylabel('n [um^-1]')
    xticks([ -60, -30 , 0 30, 60])

    suptitle(['t = ' num2str(round(t_array_si(i)/(2*pi/omega_si) ,2 )) 'T'])
    pause(0.02)
end


%% ------------------------------------------------
function theta_bragg = applyBraggPulse(theta, k_Bragg, dk, LLS)
    rho         = LLS.transform2rho(theta);
    rho_shift   = circshift(permute(rho/2, [3 2 1]), [0, round(k_Bragg/dk)]);
    rho_bragg   = permute(rho_shift + fliplr(rho_shift) , [3 2 1]);
    theta_bragg = LLS.transform2theta(rho_bragg);
end