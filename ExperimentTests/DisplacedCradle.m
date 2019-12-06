clear all; close all; clc;

addpath([pwd '\Functions\'])


%% Simulation parameters

N           = 2^7;
M           = 2^7;

Nsteps      = 10*2*120;
stepOrder   = 2;

xmax_si     = 150*1e-6;
kmax_si     = 1*1e7;
tmax_si     = 10*12*1e-3;

%% Experimental units

m_si        = 87*1.6605402e-27;
hbar_si     = 1.054571726e-34;
kB_si       = 1.38065e-23;
as_si       = 5.2e-9;

Natoms      = 4000;
T_si        = 60*1e-9; 
omega_si    = 2*pi*14; % [s^-1]
omegaT_si   = 2*pi*1.6*1e3;

g1D_si      = 2*hbar_si*omegaT_si*as_si;
c_si        = g1D_si*m_si/hbar_si^2;


%% create arrays
x_array     = linspace(-xmax_si,xmax_si,M);
k_array_i   = linspace(-kmax_si,kmax_si,N);
% [k_array_i,kw]=legzo(N, -kmax_si, kmax_si);
t_array     = linspace(0, tmax_si, Nsteps+1);


%% Initialize state 
Vx_harm = @(t,x) 0.5*m_si*omega_si^2*x.^2;
a = 3e-21;
b = 9e-14;
Vx_dw   = @(t,x) -0.25* a *x.^2 + 0.5* b *x.^4 + 1/32 * a^2/b;

figure
hold on
plot(x_array, Vx_harm(0,x_array))
plot(x_array, Vx_dw(0,x_array))


dmudx   = @(t,x) -m_si*omega_si^2*x;


couplings = { @(t,x) 0              , @(t,x) c_si;
              []                    , []            ;
              @(t,x) dmudx(t,x)     , []            };


LLS         = LiebLinigerSolver_SI(omegaT_si, x_array, k_array_i, k_array_i(2)-k_array_i(1), couplings);

mu0_fit     = LLS.fitAtomnumber(T_si, Vx_dw, Natoms, 1e-31, true);
theta_init = LLS.calcThermalState(T_si);


couplings{1,1} = @(t,x) mu0_fit - Vx_harm(t,x);



%%

% [k_array,kw]=legzo(2*N, -4*kmax_si, 4*kmax_si);
% LLS         = LiebLinigerSolver_SI(omegaT_si, x_array, k_array, kw, couplings);
k_array   = linspace(-2*kmax_si,2*kmax_si,2*N);

theta_mat = squeeze(double(theta_init));
theta_mat2 = padarray(theta_mat, [N/2, 0]);


figure
subplot(1,2,1)
imagesc(theta_mat)

subplot(1,2,2)
imagesc(theta_mat2)


theta_init = GHDtensor( permute(theta_mat2 , [1 5 3 4 2]) );


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

% %% Run movie 
% figure
% for i = 1:length(t_array)
%     theta = double(theta_t{i});
%     
%     ax1 = subplot(2,1,1);
%     imagesc(x_array*1e6,k_array*1e-6,squeeze(theta))
%     colormap(hot)
%     caxis([ 0 1])
%     ylabel('rapidity')
%     set(gca,'YDir','normal') 
% 
%     
%     ax2 = subplot(2,1,2);
%     area(x_array*1e6, n_t(:,i)*1e-6 )
%     xlim([ -xmax_si*1e6, xmax_si*1e6])
%     ylim([0 max( n_t(:)*1e-6 )])
%     xlabel('x [um]')
%     ylabel('atomic density [um\^(-1)]')
%     
%     samexaxis('xmt','on','ytac','join','yld',1,'YTickAntiClash')
% %     suptitle(['t = ' num2str(t_array(i))])
% 
% 
%     pause(0.01)
% end


