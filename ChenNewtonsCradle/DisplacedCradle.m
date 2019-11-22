clear all; close all; clc;

addpath([pwd '\Functions\'])


%% Simulation parameters

N           = 2^7;
M           = 2^7;

Nsteps      = 10*2*120;
stepOrder   = 2;

xmax_si     = 60*1e-6;
kmax_si     = 4*1e7;
tmax_si     = 10*12*1e-3;

%% Experimental units

m_si        = 87*1.6605402e-27;
hbar_si     = 1.054571726e-34;
kB_si       = 1.38065e-23;
as_si       = 5.2e-9;

Natoms      = 120;
T_si        = 100*1e-9; 
omega_si    = 2*pi*83.3; % [s^-1]
omegaT_si   = 2*pi*31*1e3;

g1D_si      = 2*hbar_si*omegaT_si*as_si;
c_si        = g1D_si*m_si/hbar_si^2;


%% create arrays
x_array     = linspace(-xmax_si,xmax_si,M);
k_array     = linspace(-kmax_si,kmax_si,N);
t_array     = linspace(0, tmax_si, Nsteps+1);


%% Initialize state 
offset  = 35*1e-6;    

Vx      = @(t,x) 0.5*m_si*omega_si^2*x.^2;
Vx_off  = @(t,x) 0.5*m_si*omega_si^2*(x-offset).^2;

dmudx   = @(t,x) -m_si*omega_si^2*x;


couplings = { @(t,x) 0              , @(t,x) c_si;
              []                    , []            ;
              @(t,x) dmudx(t,x)     , []            };


LLS         = LiebLinigerSolver_SI_test(omega_si, x_array, k_array, couplings);

mu0_fit     = LLS.fitAtomnumber(T_si, Vx_off, Natoms/2, true);
theta_init1 = LLS.calcThermalState(T_si);


couplings{1,1} = @(t,x) mu0_fit - Vx(t,x);
LLS.setCouplings( couplings );
theta_init2  = LLS.calcThermalState(T_si);


theta_init  = theta_init1 + theta_init2;


figure
subplot(1,3,1)
imagesc(squeeze(double(theta_init1)))

subplot(1,3,2)
imagesc(squeeze(double(theta_init2)))

subplot(1,3,3)
imagesc(squeeze(double(theta_init)))



%% Solve dynamics and calculate quantities
theta_t     = LLS.propagateTheta(theta_init, t_array);
n_t         = LLS.calcCharges(theta_t, 0, t_array);

%% Plot theta
figure

for i = 1:120
    index = 1 + 10*(i-1);
    
    subaxis(10,12,i,'SpacingHoriz',0.001,'SpacingVert',0.001) 
    imagesc( squeeze(double(theta_t{index})) )
    caxis([0 1])
    colormap(hot)
    set(gca,'xticklabel',[])
    set(gca,'yticklabel',[])
    set(gca,'YDir','normal')
    
end

%% Plot density carpet
figure
imagesc(x_array*1e6, t_array*1e3, n_t' * 1e-6)
xlabel('x [um]')
ylabel('t [ms]')
set(gca,'YDir','normal')

colormap(jet)
cb = colorbar;
ylabel(cb, 'atomic density [um\^(-1)]')
caxis([0 max(n_t(:)*1e-6)])

saveas(gcf,'DisplacedCradleDensityCarpet.png')

%% Run movie 
figure
for i = 1:length(t_array)
    theta = double(theta_t{i});
    
    ax1 = subplot(2,1,1);
    imagesc(x_array*1e6,k_array*1e-6,squeeze(theta))
    colormap(hot)
    caxis([ 0 1])
    ylabel('rapidity')
    set(gca,'YDir','normal') 

    
    ax2 = subplot(2,1,2);
    area(x_array*1e6, n_t(:,i)*1e-6 )
    xlim([ -xmax_si*1e6, xmax_si*1e6])
    ylim([0 max( n_t(:)*1e-6 )])
    xlabel('x [um]')
    ylabel('atomic density [um\^(-1)]')
    
    samexaxis('xmt','on','ytac','join','yld',1,'YTickAntiClash')
%     suptitle(['t = ' num2str(t_array(i))])


    pause(0.01)
end


%% Create gif

giffig = figure;
filename = 'DisplayedCradle.gif';
for i = 1:2:length(t_array)
    % Draw plot 
    theta = double(theta_t{i});
    
    ax1 = subplot(2,1,1);
    imagesc(x_array*1e6,k_array*1e-6,squeeze(theta))
    colormap(hot)
    caxis([ 0 1])
    ylabel('rapidity')
    set(gca,'YDir','normal') 
    title(['t = ' num2str(round(t_array(i)*1e3,1)) 'ms'],'Fontsize',14)
    
    ax2 = subplot(2,1,2);
    area(x_array*1e6, n_t(:,i)*1e-6 )
    xlim([ -xmax_si*1e6, xmax_si*1e6])
    ylim([0 max( n_t(:)*1e-6 )])
    xlabel('x [um]')
    ylabel('atomic density [um\^(-1)]')

    
    samexaxis('xmt','on','ytac','join','yld',1,'YTickAntiClash')
    drawnow 
  
    % Capture the plot as an image 
    frame = getframe(giffig); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
    % Write to the GIF File 
    if i == 1 
        imwrite(imind,cm,filename,'gif','DelayTime', 0.02, 'Loopcount',inf); 
    else 
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime', 0.02); 
    end 
end