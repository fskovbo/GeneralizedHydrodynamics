clear all; close all;

addpath('..\Functions\')

%% Define simulation parameters

N           = 2^7;
M           = 2^6;
dt          = 0.025;

kmax        = 8;
xmax        = 7;
tmax        = 3;

k_array     = linspace(-kmax, kmax, N);
x_array     = linspace(-xmax, xmax, M);
t_array     = linspace(0, tmax, tmax/dt+1);

stepOrder   = 2;
extrapFlag  = false;


%% Define physical couplings and temperature
couplings   = { @(t,x) 4 - 2*x.^2 , @(t,x) 1 };

T           = 2;


%% Initialize solver and generate initial state

LLS         = LiebLinigerSolver(x_array, k_array, couplings);
theta_init  = LLS.calcThermalState(T);

k_bragg     = 3.5;
delta_k     = k_array(2)-k_array(1);
theta_init  = applyBraggPulse(theta_init, k_bragg, delta_k, LLS);


%% Solve dynamics and calculate density

theta_t     = LLS.propagateTheta(theta_init, t_array);
n_t         = LLS.calcCharges(theta_t, 0);


%% Plot density carpet
figure
imagesc(x_array, t_array, n_t')
set(gca,'YDir','normal') 
xlabel('x')
ylabel('t')
caxis([0 1.7])
colormap(jet)
colorbar
title('Atomic Density')


%% Plot atomnumber
figure
plot(t_array, sum(n_t,1) )
xlabel('t')
ylabel('Atomnumber')


%% Plot theta snapshots
figure
for i = 1:12
    theta = double(theta_t{i});
    
    subplot(3,4,i)
    index = floor( (length(t_array)-1)/11*(i-1) + 1 );
    imagesc(x_array,k_array,squeeze(theta))
    colormap(hot)
    caxis([ 0 0.5])
    set(gca,'YDir','normal') 
    title(['t = ' num2str(t_array(index))])
end


%% Plot density snapshots
figure
for i = 1:12
   subplot(3,4,i)
   index = floor( (length(t_array)-1)/11*(i-1) + 1 );
   plot(x_array, n_t(:,index) )
   xlim([ -xmax, xmax])
   ylim([0 max( n_t(:) )])
   title(['t = ' num2str(t_array(index))])
end


%% Run movie
figure
for i = 1:length(t_array)
    theta = double(theta_t{i});
    
    subplot(2,1,1)
    imagesc(x_array,k_array,squeeze(theta))
    colormap(hot)
    caxis([ 0 0.5])
    ylabel('rapidity')
    set(gca,'YDir','normal') 
    
    subplot(2,1,2)
    plot(x_array, n_t(:,i) )
    xlim([ -xmax, xmax])
    ylim([0 max( n_t(:) )])
    xlabel('x')
    ylabel('density')
    
    suptitle(['t = ' num2str(t_array(i))])
    pause(0.05)
end


%%
function theta_bragg = applyBraggPulse(theta, k_Bragg, dk, LLS)
    rho         = LLS.transform2rho(theta);
    rho_shift   = circshift( double(rho)/2, [round(k_Bragg/dk), 0]);
    rho_bragg   = GHDtensor( rho_shift + flipud(rho_shift) );
    theta_bragg = LLS.transform2theta(rho_bragg);
end

