clear all; close all;

addpath('..\Functions\')

%% define parameters
N           = 2^6;
M           = 2^7;
Ntypes      = 1;
dt          = 0.005;

kmax        = pi/2;
xmax        = 1;

k_array     = linspace(-kmax, kmax, N);
x_array     = linspace(-xmax, xmax, M);
tcorr_array = [0.1, 0.2];%0.1:0.1:0.3;

stepOrder   = 2;
extrapFlag  = false;


T           = 0.25;


%% Define physical couplings
couplings.B         = @(t,x) -1; 
couplings.dBdx      = @(t,x) 0;
couplings.dBdt      = [];
couplings.Delta     = @(t,x) 1.5;
couplings.dDeltadt  = [];
couplings.dDeltadx  = [];


coup_init           = couplings;
coup_init.B         = @(t,x) -1 - 8*x.^2; 


%% Run simulation
c_idx           = [2, 2];
areCurrents     = [0, 0];

XXZ             = XXZchainSolver(x_array, k_array, couplings, Ntypes, stepOrder, extrapFlag);
[corrmat, theta, C1P] = XXZ.calcCorrelationMatrix(T, coup_init, tcorr_array, dt, c_idx, areCurrents);


%% PLOT
for i = 1:length(tcorr_array)
 
    figure
    
    subplot(3,2,[1 3])
    imagesc(x_array,x_array, corrmat(:,:,1,i))
    ylabel('x')
    xlabel('y')
    c1 = colorbar;
    caxis([0 1.9])
    c1.Location = 'northoutside';
    set(gca,'YDir','normal')
    title('direct')
    
    subplot(3,2,[2 4])
    imagesc(x_array,x_array, corrmat(:,:,2,i))
    ylabel('x')
    xlabel('y')
    c2 = colorbar;
    c2.Location = 'northoutside';
    set(gca,'YDir','normal')
    title('indirect')
    
    subplot(3,2,5)
    plot(x_array, C1P(:,1))
    xlim([-xmax, xmax])
    ylabel(['O(y, 0)'])
    xlabel('y')
    
    subplot(3,2,6)
    plot(x_array, C1P(:,1+i))
    xlim([-xmax, xmax])
    ylabel(['O(x, t = ' num2str(tcorr_array(i)) ')'])
    xlabel('x')
    
    suptitle(['<n(x, t = ' num2str(tcorr_array(i)) ') n(y, 0)>'])
    
end

%%
for i = 1:length(tcorr_array)
 
    figure
    
    imagesc(x_array,x_array, sum(corrmat(:,:,:,i),3))
    ylabel('x')
    xlabel('y')
    colorbar
    caxis([-0.1 1])
    set(gca,'YDir','normal')

    title(['<n(x, t = ' num2str(tcorr_array(i)) ') n(y, 0)>'])
    
end
