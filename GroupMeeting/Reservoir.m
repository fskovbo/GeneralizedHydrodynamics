clear all; close all;

addpath(['..' filesep 'Functions' filesep])

%% Define simulation parameters

N           = 2^7;
M           = 2^8;
dt          = 0.025;

kmax        = 9;
xmax        = 20;
tmax        = 20;

k_array     = linspace(-kmax, kmax, N);
x_array     = linspace(-xmax, xmax, M);
t_array     = linspace(0, tmax, tmax/dt+1);

Options.simplifySteps = 0;


%% Define physical couplings and temperature
boxTrap     = @(t,x) 14*(-tanh(2.5*(x+18.5)) + tanh(2.5*(x-18.5)) + 2);
barrier     = @(t,x) 25*exp(-x.^2 / 0.5^2);
 
% couplings   = { @(t,x) 4 - 14*(-tanh(2.5*(x+18.5)) + tanh(2.5*(x-18.5)) + 2) - 0.124*tanh(2*x) - (1-tanh(3*t))*25*exp(-x.^2 / 0.5^2)    , @(t,x) 1 };


couplings   = { @(t,x) 4 - 14*(-tanh(2.5*(x+18.5)) + tanh(2.5*(x-18.5)) + 2) , @(t,x) 1};

coup_init   = couplings;
% coup_init{1}= @(t,x) 4 - boxTrap(t,x) - barrier(t,x) - 0.124*tanh(2*x);
coup_init{1}= @(t,x) 4 - boxTrap(t,x) - 0.124*tanh(100*x);

T           = @(x) 2*heaviside(x) + 4*heaviside(-x);


% coup_init{1}= @(t,x) 4 - boxTrap(t,x) - 0.065*tanh(100*x);
% 
% T           = @(x) 2*heaviside(x) + 8*heaviside(-x);



%% Initialize solver and generate initial state

LLS         = LiebLinigerSolver(x_array, k_array, couplings, Options);
theta_init  = LLS.calcThermalState(T, coup_init);

n_0         = LLS.calcCharges(theta_init, 0);
% figure
% subplot(2,1,1)
% imagesc( squeeze( double(theta_init)) )
% subplot(2,1,2)
% plot(n_0)
% return

%% Solve dynamics and calculate density

theta_t     = LLS.propagateTheta(theta_init, t_array);
n_t         = LLS.calcCharges(theta_t, 0);
p_t         = LLS.calcCharges(theta_t, 1);
e_t         = LLS.calcCharges(theta_t, 2);


%% Plot density carpet
figure
imagesc(x_array, t_array, n_t')
set(gca,'YDir','normal') 
xlabel('x')
ylabel('t')
caxis([0 max( n_t(:))])
colormap(jet)
colorbar
title('Atomic Density')

saveas(gcf,'DensityCarpet_reservoir.png')

%% Plot momentum carpet
figure
imagesc(x_array, t_array, p_t')
set(gca,'YDir','normal') 
xlabel('x')
ylabel('t')
caxis([min( p_t(:)) max( p_t(:))])
colormap(jet)
colorbar
title('Momentum Density')

saveas(gcf,'MomentumCarpet_reservoir.png')

%% Plot energy carpet
figure
imagesc(x_array, t_array, e_t')
set(gca,'YDir','normal') 
xlabel('x')
ylabel('t')
caxis([min( e_t(:)) max( e_t(:))])
colormap(jet)
colorbar
title('Energy Density')


saveas(gcf,'EnergyCarpet_reservoir.png')

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
    caxis([ 0 1 ])
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
% figure
% for i = 1:length(t_array)
%     theta = double(theta_t{i});
%     
%     ax1 = subplot(2,1,1);
%     imagesc(x_array,k_array,squeeze(theta))
%     colormap(hot)
%     caxis([ 0 1])
%     ylabel('rapidity')
%     set(gca,'YDir','normal') 
% 
%     
%     ax2 = subplot(2,1,2);
%     area(x_array, n_t(:,i) )
%     xlim([ -xmax, xmax])
%     ylim([0 max( n_t(:) )])
%     xlabel('x')
%     ylabel('density')
%     
%     samexaxis('xmt','on','ytac','join','yld',1,'YTickAntiClash')
% %     suptitle(['t = ' num2str(t_array(i))])
% 
% 
%     pause(0.05)
% end

%% make gifs
fig1 = figure;
% axis tight manual % this ensures that getframe() returns a consistent size
filename = 'density_reservoir.gif';
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

    
    samexaxis('xmt','on','ytac','join','yld',1,'YTickAntiClash')
    drawnow 
  
    % Capture the plot as an image 
    frame = getframe(fig1); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
    % Write to the GIF File 
    if i == 1 
        imwrite(imind,cm,filename,'gif','DelayTime', 0.02, 'Loopcount',inf); 
    else 
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime', 0.02); 
    end 
end



