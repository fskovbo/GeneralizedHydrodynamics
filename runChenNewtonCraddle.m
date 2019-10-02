clear all; close all;

addpath([pwd '\Functions\'])

%% Simulation parameters

N           = 2^7;
M           = 2^6;

Nsteps      = 200;
stepOrder   = 2;

xmax_si     = 50*1e-6;
kmax_si     = 3*1e7;
tmax_si     = 12*1e-3;

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


%% Initialize state 
Vx      = @(t,x) 0.5*m_si*omega_si^2*x.^2;
dmudx   = @(t,x) -m_si*omega_si^2*x;


couplings.mu    = @(t,x) 0;
couplings.dmudx = dmudx;
couplings.c     = @(t,x) c_si;
couplings.dcdt  = [];
couplings.dcdx  = [];
couplings.dmudt = [];


LLS         = LiebLinigerSolver_SI(omegaT_si, x_array, k_array, couplings);
mu0_fit     = LLS.fitAtomnumber(T_si, Vx, Natoms, true);

theta_init  = LLS.calcThermalState(T_si);

delta_k     = (k_array(2)-k_array(1));
theta_bragg = applyBraggPulse(theta_init, k_Bragg_si, delta_k, LLS);


%% Solve dynamics and calculate quantities
theta_t     = LLS.propagateTheta(theta_bragg, t_array);
n_t         = LLS.calcCharges(theta_t, 0, t_array) * 1e-6; % [um^-1]
rho_t       = LLS.transform2rho( theta_t, t_array); % [um um^-1]
nk_t        = squeeze(trapz(x_array, permute(rho_t, [3 2 1 4]))) * 1e6; % [um]

%% Plot results

figure
for i = 1:6
    index = ceil( (length(t_array)-1)/5*(i-1) + 1 );
    
    % plot Harmonic quasiparticle distribution
    subplot(3,6,i)
    imagesc(x_array*1e6, k_array*1e-6,squeeze(rho_t(:,:,:,index)))
    set(gca,'YDir','normal')
    colormap(hot)
    caxis([0 0.3])
    ylabel('k [um^-1]')
    set(gca,'xticklabel',[])
    title(['t = ' num2str(round(t_array(index)/(2*pi/omega_si) ,2 )) 'T'])

    % plot density profiles in real space
    subplot(3,6,i+6)
    hold on
    box on
    plot(x_array*1e6, n_t(:,index))
    xlim([ -xmax_si*1e6, xmax_si*1e6])
    ylim([0 10])
    xlabel('x [um]')
    ylabel('n [um^-1]')
    xticks([ -60, -30 , 0 30, 60])
    
    
    % plot density momentum space
    subplot(3,6,i+12)
    hold on
    box on
    plot(k_array*1e-6, nk_t(:,index))
    xlim([ -kmax_si*1e-6, kmax_si*1e-6])
%     ylim([0 10])
    xlabel('k [um^-1]')
    ylabel('n [um]')
%     xticks([ -60, -30 , 0 30, 60])
end



%%
figure
for i = 1:24
    index = ceil( (length(t_array)-1)/23*(i-1) + 1 );
    
    subplot(4,6,i)
    plot(k_array*1e-6, nk_t(:,index))
    xlim([ -kmax_si*1e-6, kmax_si*1e-6])
    ylim([0 10])
    xlabel('k [um^-1]')
    ylabel('n [um]')
end


figure
for i = 1:24
    index = ceil( (length(t_array)-1)/23*(i-1) + 1 );
    
    subplot(4,6,i)
    plot(x_array*1e6, n_t(:,index))
    xlim([ -xmax_si*1e6, xmax_si*1e6])
    ylim([0 10])
    xlabel('x [um]')
    ylabel('n [um^-1]')    
end


figure
for i = 1:24
    index = ceil( (length(t_array)-1)/23*(i-1) + 1 );
    
    subplot(4,6,i)
    imagesc(x_array*1e6, k_array*1e-6,squeeze(rho_t(:,:,:,index)))
    set(gca,'YDir','normal')
    colormap(hot)
    caxis([0 0.3])
    ylabel('k [um^-1]')
    xlabel('x [um]')
  
end


%% ------------------------------------------------
function theta_bragg = applyBraggPulse(theta, k_Bragg, dk, LLS)
    rho         = LLS.transform2rho(theta);
    rho_shift   = circshift(permute(rho/2, [3 2 1]), [0, round(k_Bragg/dk)]);
    rho_bragg   = permute(rho_shift + fliplr(rho_shift) , [3 2 1]);
    theta_bragg = LLS.transform2theta(rho_bragg);
end