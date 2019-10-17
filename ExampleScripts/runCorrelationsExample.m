clear all; close all;

addpath('..\Functions\')

%% define parameters

N       = 2^7;
M       = 2^6;

kmax    = 7;
xmax    = 5;

T       = 5;

dt      = 0.01;


k_array     = linspace(-kmax,kmax,N);
x_array     = linspace(-xmax,xmax,M);
tcorr_array = 0.2:0.2:0.4;


%% Run simulations
couplings   = { @(t,x) 5 , @(t,x) 1 };

LLS             = LiebLinigerSolver(x_array, k_array, couplings);

c_idx           = [0, 0];
areCurrents     = [0, 0];


c = 4;
b = 5/c^4;
a = 2*b*c^2;

coup_init       = couplings;
% coup_init.mu    = @(t,x) a*x.^2 - b*x.^4;
coup_init{1}    = @(t,x) 5 - 5*x.^2;

[corrmat, theta, C1P] = LLS.calcCorrelationMatrix(T, coup_init, tcorr_array, dt, c_idx, areCurrents);


%% PLOT
for i = 1:length(tcorr_array)
 
    figure
    
    
    subplot(3,2,[1 3])
    imagesc(x_array,x_array, corrmat(:,:,1,i))
    ylabel('x')
    xlabel('y')
    c1 = colorbar;
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
% 
% for i = 1:length(tcorr_array)
% %     q_0 = LLS.calcCharges(GHDcorrC.theta(:,:,:,1), 0, 0);
% %     q_t = LLS.calcCharges(GHDcorrC.theta(:,:,:,i+1), 0, 1);
%     
%     figure
%     
%     
%     imagesc(x_array,x_array, sum(corrmat(:,:,:,i), 3))
%     ylabel('x')
%     xlabel('y')
%     colorbar
%     set(gca,'YDir','normal')
%     
%     suptitle(['<n(x, t = ' num2str(tcorr_array(i)) ') n(y, 0)>'])
%     
% end


% for i = 1:length(tcorr_array)
%     figure
%     hold on
%     imagesc(x_array,x_array, corrmat(:,:,2,i))
%     line([-xmax, xmax],[-xmax, xmax],'Color','k')
%     ylabel('x')
%     xlabel('y')
%     axis tight
%     
%     set(gca,'YDir','normal')
%      
% end