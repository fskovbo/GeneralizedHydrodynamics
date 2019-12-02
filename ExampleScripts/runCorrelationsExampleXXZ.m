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
% tcorr_array = [0.1, 0.2];%0.1:0.1:0.3;
tcorr_array = dt:dt:1;


T           = 0.25;


%% Define physical couplings
couplings       = { @(t,x) -1 , @(t,x) acosh(1.5) };

coup_init       = couplings;
coup_init{1}    = @(t,x) -1 - 8*x.^2; 


%% Run simulation
c_idx           = [0, 0];
areCurrents     = [0, 0];

XXZ             = XXZchainSolver(x_array, k_array, couplings, Ntypes);
% [corrmat, theta, C1P] = XXZ.calcCorrelationMatrix(T, coup_init, tcorr_array, dt, c_idx, areCurrents);

corrfunc = XXZ.calcCorrelationFunction(T, coup_init, tcorr_array, dt, c_idx, areCurrents);


%% PLOT
% for i = 1:length(tcorr_array)
%  
%     figure
%     
%     subplot(3,2,[1 3])
%     imagesc(x_array,x_array, corrmat(:,:,1,i))
%     ylabel('x')
%     xlabel('y')
%     c1 = colorbar;
%     caxis([0 0.7])
%     c1.Location = 'northoutside';
%     set(gca,'YDir','normal')
%     title('direct')
%     
%     subplot(3,2,[2 4])
%     imagesc(x_array,x_array, corrmat(:,:,2,i))
%     ylabel('x')
%     xlabel('y')
%     c2 = colorbar;
%     caxis([-0.01 0.05])
%     c2.Location = 'northoutside';
%     set(gca,'YDir','normal')
%     title('indirect')
%     
%     subplot(3,2,5)
%     plot(x_array, C1P(:,1))
%     xlim([-xmax, xmax])
%     ylabel(['O(y, 0)'])
%     xlabel('y')
%     
%     subplot(3,2,6)
%     plot(x_array, C1P(:,1+i))
%     xlim([-xmax, xmax])
%     ylabel(['O(x, t = ' num2str(tcorr_array(i)) ')'])
%     xlabel('x')
%     
%     suptitle(['<n(x, t = ' num2str(tcorr_array(i)) ') n(y, 0)>'])
%     
% end
% 
% %%
% for i = 1:length(tcorr_array)
%  
%     figure
%     
%     imagesc(x_array,x_array, sum(corrmat(:,:,:,i),3))
%     ylabel('x')
%     xlabel('y')
%     colorbar
%     caxis([-0.1 1])
%     set(gca,'YDir','normal')
% 
%     title(['<n(x, t = ' num2str(tcorr_array(i)) ') n(y, 0)>'])
%     
% end


%% PLOT correlation function
figure

subplot(1,3,1)
imagesc(x_array,tcorr_array,corrfunc(:,:,1)')
set(gca,'YDir','normal')
xlabel('x')
ylabel('t')
caxis([0 0.2])
title('direct')

subplot(1,3,2)
imagesc(x_array,tcorr_array,corrfunc(:,:,2)')
set(gca,'YDir','normal')
xlabel('x')
ylabel('t')
caxis([-0.01 0.01])
title('indirect')

subplot(1,3,3)
imagesc(x_array,tcorr_array,corrfunc(:,:,1)' + corrfunc(:,:,2)')
set(gca,'YDir','normal')
xlabel('x')
ylabel('t')
caxis([-0.01 0.2])
title('full')

suptitle(['<O(x, t) O(0, 0)>'])