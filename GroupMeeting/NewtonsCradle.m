clear all; close all;

addpath(['..' filesep 'Functions' filesep])

%% Define simulation parameters

N           = 2^7;
M           = 2^7;
dt          = 0.025;

kmax        = 13;
xmax        = 9;
tmax        = 15;

k_array     = linspace(-kmax, kmax, N);
x_array     = linspace(-xmax, xmax, M);
t_array     = linspace(0, tmax, tmax/dt+1);


%% Define physical couplings and temperature
couplings   = { @(t,x) 4 - 2*x.^2 , @(t,x) 1 };
coup_offset = { @(t,x) 4 - 2*(x-6).^2 , @(t,x) 1 };


T           = 2;


%% Initialize solver and generate initial state

LLS         = LiebLinigerSolver(x_array, k_array, couplings);
theta_init1 = LLS.calcThermalState(T);
theta_init2 = LLS.calcThermalState(T, coup_offset);

theta_init  = theta_init1 + theta_init2;



%% Solve dynamics and calculate density

theta_t     = LLS.propagateTheta(theta_init, t_array);
n_t         = LLS.calcCharges(theta_t, 0);


%% Plot density carpet
figure
imagesc(x_array, t_array, n_t')
set(gca,'YDir','normal') 
xlabel('x')
ylabel('t')
caxis([0 max(n_t(:))])
colormap(jet)
colorbar
title('Atomic Density')

saveas(gcf,'DensityCarpet.png')

%% Plot atomnumber
figure
plot(t_array, sum(n_t,1) )
xlabel('t')
ylabel('Atomnumber')


%% Plot theta snapshots
figure
for i = 1:12
    index = floor( (length(t_array)-1)/11*(i-1) + 1 );
    theta = double(theta_t{index});
    
    subplot(3,4,i)
    imagesc(x_array,k_array,squeeze(theta))
    colormap(hot)
    caxis([ 0 1])
    set(gca,'YDir','normal') 
    title(['t = ' num2str(t_array(index))])
end


%% Plot density snapshots
figure
for i = 1:12
   subplot(3,4,i)
   index = floor( (length(t_array)-1)/11*(i-1) + 1 );
   area(x_array, n_t(:,index) )
   xlim([ -xmax, xmax])
   ylim([0 max( n_t(:) )])
   title(['t = ' num2str(t_array(index))])
end


%% Run movie
figure
for i = 1:length(t_array)
    theta = double(theta_t{i});
    
    ax1 = subplot(2,1,1);
    imagesc(x_array,k_array,squeeze(theta))
    colormap(hot)
    caxis([ 0 1])
    ylabel('rapidity')
    set(gca,'YDir','normal') 

    
    ax2 = subplot(2,1,2);
    area(x_array, n_t(:,i) )
    xlim([ -xmax, xmax])
    ylim([0 max( n_t(:) )])
    xlabel('x')
    ylabel('density')
    
    samexaxis('xmt','on','ytac','join','yld',1,'YTickAntiClash')
%     suptitle(['t = ' num2str(t_array(i))])


    pause(0.05)
end


%% Make gifs
fig1 = figure;
% axis tight manual % this ensures that getframe() returns a consistent size
filename = 'density_short_dist6.gif';
for i = 1:91
    % Draw plot 
    area(x_array, n_t(:,i) )    
    xlabel('x')
    ylabel('density') 
    xlim([ -xmax, xmax])
    ylim([0 max( n_t(:) )])
    
    hold on
    plot(x_array, 0.1*x_array.^2, 'r-','Linewidth',2)
    hold off
    
    drawnow 
  
    % Capture the plot as an image 
    frame = getframe(fig1); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
    % Write to the GIF File 
    if i == 1 
        imwrite(imind,cm,filename,'gif','DelayTime', 0.03, 'Loopcount',inf); 
    else 
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime', 0.03); 
    end 
end


fig2 = figure;
% axis tight manual % this ensures that getframe() returns a consistent size
filename = 'density_long_dist6.gif';
for i = 1:length(t_array)
    % Draw plot 
    theta = double(theta_t{i});
    
    ax1 = subplot(2,1,1);
    imagesc(x_array,k_array,squeeze(theta))
    colormap(hot)
    caxis([ 0 1])
    ylabel('rapidity')
    set(gca,'YDir','normal') 
    title(['t = ' num2str(round(t_array(i),1) )],'Fontsize',14)
    
    ax2 = subplot(2,1,2);
    area(x_array, n_t(:,i) )
    xlim([ -xmax, xmax])
    ylim([0 max( n_t(:) )])
    xlabel('x')
    ylabel('density')
    
    hold on
    plot(x_array, 0.1*x_array.^2, 'r-','Linewidth',2)
    hold off

    
    samexaxis('xmt','on','ytac','join','yld',1,'YTickAntiClash')
    drawnow 
  
    % Capture the plot as an image 
    frame = getframe(fig2); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
    % Write to the GIF File 
    if i == 1 
        imwrite(imind,cm,filename,'gif','DelayTime', 0.03, 'Loopcount',inf); 
    else 
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime', 0.03); 
    end 
end


%%
function theta_bragg = applyBraggPulse(theta, k_Bragg, dk, LLS)
    rho         = LLS.transform2rho(theta);
    rho_shift   = circshift( double(rho)/2, [round(k_Bragg/dk), 0]);
    rho_bragg   = GHDtensor( rho_shift + flipud(rho_shift) );
    theta_bragg = LLS.transform2theta(rho_bragg);
end

