clear all; close all; clc;

addpath([pwd '\Functions\'])


%%
data_init   = load('initial_data');

dk          = 4*pi/(852*1e-9)/52;
k_arr_exp   = 0:dk:(length(data_init.N120_after) - 1)*dk;
k_arr_exp   = k_arr_exp - mean(k_arr_exp);

nk_init_b   = 120*data_init.N120_before/trapz(k_arr_exp, data_init.N120_before);
nk_init_a   = 120*data_init.N120_after/trapz(k_arr_exp, data_init.N120_after);


%% Simulation parameters

N           = 2^8;
M           = 2^7;

xmax_si     = 20*1e-6;
kmax_si     = 1*1e7;

%% Experimental units
m_si        = 87*1.6605402e-27;
hbar_si     = 1.054571726e-34;
kB_si       = 1.38065e-23;
as_si       = 5.2e-9;

Natoms      = 120;
omega_si    = 2*pi*83.3; % [s^-1]
omegaT_si   = 2*pi*31*1e3;
k_Bragg_si  = 4*pi/(852*1e-9); % [m^-1]

g1D_si      = 2*hbar_si*omegaT_si*as_si;
c_si        = g1D_si*m_si/hbar_si^2;

%% create arrays
x_array     = linspace(-xmax_si,xmax_si,M);
k_array     = linspace(-kmax_si,kmax_si,N);


%% Generate initial state
T_si        = [1 10 50 100 200] * 1e-9; 
nk_init     = zeros(N, length(T_si)+1);

Vx          = @(t,x) 0.5*m_si*omega_si^2*x.^2;
couplings   = { @(t,x) 0              , @(t,x) c_si };
LLS         = LiebLinigerSolver_SI_test(omega_si, x_array, k_array, couplings);

for i = 1:length(T_si)
    mu0_fit     = LLS.fitAtomnumber(T_si(i), Vx, Natoms, true);
    theta_init  = LLS.calcThermalState(T_si(i));
    nk_init(:,i)= LLS.calcMomentumDistr(theta_init, 0);
end
    
%%
rho_init        = calcZeroTempState(x_array, k_array, c_si, Vx(0,x_array), Natoms);
nk_init(:,end)  = trapz(x_array, rho_init');
%% plot initial state
figure
hold on
box on

plot(k_array*1e-6 , nk_init*1e6)
plot(k_arr_exp*1e-6, nk_init_b*1e6, 'k')

ylim([0 40])
xlabel('k [um^-1]')
ylabel('n [um^-1]')
xlim([ -kmax_si*1e-6, kmax_si*1e-6])


% figure
% imagesc( k_array*1e-6, x_array*1e6, rho_init )



%% ------------------------------------------------
function theta_bragg = applyBraggPulse(theta, k_Bragg, dk, LLS)
    rho         = LLS.transform2rho(theta);
    rho_shift   = circshift( double(rho)/2, [round(k_Bragg/dk), 0]);
    rho_bragg   = GHDtensor( rho_shift + flipud(rho_shift) );
    theta_bragg = LLS.transform2theta(rho_bragg);
end


function theta_bragg = applyBraggPulse_test(theta, k_Bragg, dk, subamp, subwidth, LLS)
    rho         = LLS.transform2rho(theta);
    
    % Create central blob
    rho_center  = squeeze( subamp*double(rho) );
    center_amp  = sum(sum(rho_center));
    
    x_grid      = 1:size(rho_center, 2);
    y_grid      = 1:size(rho_center, 1);
    
    x_grid      = x_grid - mean(x_grid);
    y_grid      = y_grid - mean(y_grid);
    
    [X,Y]       = meshgrid( x_grid, y_grid' );
    [XX,YY]     = meshgrid( x_grid, subwidth*y_grid' );
    
    rho_center  = interp2( XX,YY,rho_center,X,Y, 'spline');
    rho_center  = rho_center * center_amp/sum(sum(rho_center));
    rho_center  = permute(rho_center, [1 5 3 4 2]);

    % Create bragg peaks
    rho_peak    = circshift( (1-subamp)*double(rho)/2, [round(k_Bragg/dk), 0]);
     
    rho_bragg   = GHDtensor( rho_peak + flipud(rho_peak) + rho_center );
    theta_bragg = LLS.transform2theta(rho_bragg);
end

