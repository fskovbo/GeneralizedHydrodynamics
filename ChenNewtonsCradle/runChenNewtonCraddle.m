clear all; %close all; clc;

addpath([pwd '\Functions\'])


%%
% data_init = load('initial_data');
% data_t = load('MDF_N030');
% 
% 
% dk = 4*pi/(852*1e-9)/18;
% k_arr_exp = 0:dk:(length(data_init.N030_after) - 1)*dk;
% k_arr_exp = k_arr_exp - mean(k_arr_exp);
% 
% nk_init_exp = 30*data_init.N030_after/trapz(k_arr_exp, data_init.N030_after);
% 
% nk_exp = 30*data_t.MDF_until15/trapz(k_arr_exp, data_init.N030_after);
% 
% t_arr_exp = (0:size(nk_exp,1)-1) * 1e-3;
% 
% figure
% plot(k_arr_exp, nk_init_exp)
% 
% trapz(k_arr_exp, nk_init_exp)


%%
data_init = load('initial_data');
data_t = load('MDF_N120');


dk = 4*pi/(852*1e-9)/52;
k_arr_exp = 0:dk:(length(data_init.N120_after) - 1)*dk;
k_arr_exp = k_arr_exp - mean(k_arr_exp);

nk_init_exp = 120*data_init.N120_after/trapz(k_arr_exp, data_init.N120_after);

nk_exp = 120*data_t.MDF_until30/trapz(k_arr_exp, data_init.N120_after);

t_arr_exp = (0:size(nk_exp,1)-1) * 1e-3;

figure
plot(k_arr_exp, nk_init_exp)

trapz(k_arr_exp, nk_init_exp)


%% Simulation parameters

N           = 2^7;
M           = 2^7;

Nsteps      = 100*120;
stepOrder   = 2;

xmax_si     = 40*1e-6;
kmax_si     = 3*1e7;
tmax_si     = 100*12*1e-3;

%% Experimental units

m_si        = 87*1.6605402e-27;
hbar_si     = 1.054571726e-34;
kB_si       = 1.38065e-23;
as_si       = 5.2e-9;

Natoms      = 120;
T_si        = 100*1e-9; 
omega_si    = 1.0165*2*pi*83.3; % [s^-1]
% omega_si    = 1.014*2*pi*83.3; % [s^-1]
% omega_si    = 2*pi*83.3; % [s^-1]
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

couplings = { @(t,x) 0              , @(t,x) c_si;
              []                    , []            ;
              @(t,x) dmudx(t,x)     , []            };

% LLS         = LiebLinigerSolver_SI(omegaT_si, x_array, k_array, couplings);
LLS         = LiebLinigerSolver_SI_test(omega_si, x_array, k_array, couplings);
mu0_fit     = LLS.fitAtomnumber(T_si, Vx, Natoms, true);

theta_init  = LLS.calcThermalState(T_si);



couplings{1,1} = @(t,x) mu0_fit - Vx(t,x);
couplings{1,2} = @(t,x) c_si;


delta_k     = (k_array(2)-k_array(1));
theta_bragg = applyBraggPulse_test(theta_init, k_Bragg_si, delta_k, 0.2, 2.5, LLS);
% theta_bragg = applyBraggPulse(theta_init, k_Bragg_si, delta_k, LLS);

%% plot initial state
nk_init  = LLS.calcMomentumDistr(theta_bragg, 0);

figure
hold on
box on

plot(k_array*1e-6 , nk_init * 1e6)
plot(k_arr_exp*1e-6, nk_init_exp*1e6)

ylim([0 20])
xlabel('k [um^-1]')
ylabel('n [um^-1]')
xlim([ -kmax_si*1e-6, kmax_si*1e-6])


figure
imagesc(squeeze(double(theta_bragg)))


trapz( k_array, nk_init )
trapz(k_arr_exp, nk_init_exp)


%% Solve dynamics and calculate quantities
theta_t     = LLS.propagateTheta(theta_bragg, t_array);
% n_t         = LLS.calcCharges(theta_t, 0, t_array) * 1e-6; % [um^-1]
% rho_t       = LLS.transform2rho( theta_t, t_array); % [um um^-1]\


%% Plot results

% nk_rapid  = LLS.calcMomentumDistr(theta_t(1 + 10*((1:360) - 1) ), t_array(1 + 10*((1:360) - 1)));
% nk_rapid  = LLS.calcMomentumDistr(theta_t(1 + 10*((1:180) - 1) ), t_array(1 + 10*((1:180) - 1)));

% figure
% for i = 1:180
%     subaxis(15,12,i,'SpacingHoriz',0.001,'SpacingVert',0.001)
%     hold on
%     box on
%     plot(k_array*1e-6 , nk_rapid(:,i)*1e6, 'LineWidth', 1.5)
%     plot(k_arr_exp*1e-6 , nk_exp(i,:)*1e6, 'k')
%     
%     
%     line([-k_Bragg_si*1e-6  , -k_Bragg_si*1e-6],[0 20 ],'Color','black','LineStyle','--')
%     line([k_Bragg_si*1e-6 , k_Bragg_si*1e-6],[0 20 ],'Color','black','LineStyle','--')
%     
%     xlim([-30 30])
%     ylim([0 15])
%     set(gca,'xticklabel',[])
%     set(gca,'yticklabel',[])
% end
% 
% figure
% for i = 181:360
%     subaxis(15,12,i-180,'SpacingHoriz',0.001,'SpacingVert',0.001)
%     hold on
%     box on
%     plot(k_array*1e-6 , nk_rapid(:,i)*1e6, 'LineWidth', 1.5)
%     plot(k_arr_exp*1e-6 , nk_exp(i,:)*1e6, 'k')
%     
%     
%     line([-k_Bragg_si*1e-6  , -k_Bragg_si*1e-6],[0 20 ],'Color','black','LineStyle','--')
%     line([k_Bragg_si*1e-6 , k_Bragg_si*1e-6],[0 20 ],'Color','black','LineStyle','--')
%     
%     xlim([-30 30])
%     ylim([0 15])
%     set(gca,'xticklabel',[])
%     set(gca,'yticklabel',[])
% end
% 
% 
% %%
% LS_dist = zeros(1, 180);
% for i = 1:180
%     nke = nk_exp(i,:)*1e6;
%     nks = nk_rapid(:,i)*1e6;
%     
%     % try interpolating data to simultation grid
%     nke = interp1(k_arr_exp*1e-6, nke, k_array*1e-6);
%     
%     % compute square distance    
%     LS_dist(i) = sum( (nks - nke').^2 );
%     
% end
% 
% figure
% plot(0:179 , LS_dist)


%%
step = 1;

nk_rapid  = LLS.calcMomentumDistr(theta_t(1:step:end), t_array(1:step:end));

save_data.x_array = x_array;
save_data.k_array = k_array;
save_data.t_array = t_array(1:step:end);

save_data.nk_rapid = nk_rapid;
save_data.x_theta_t = theta_t(1:step:end);

save('GHD_T100_c1mid7.mat', 'save_data')

%%
% figure(69)
% for i = 1:1000
%     clf
%     index = i + (349-1)*10;
%     nk_rapid  = LLS.calcMomentumDistr(theta_t{index}, 0);
%     
%     
%     hold on
%     box on
%     
%     plot(k_array*1e-6 , nk_rapid*1e6, 'LineWidth', 1.5)
%     plot(k_arr_exp*1e-6 , nk_exp(349,:)*1e6, 'k')
%     
%     line([-k_Bragg_si*1e-6  , -k_Bragg_si*1e-6],[0 20 ],'Color','black','LineStyle','--')
%     line([k_Bragg_si*1e-6 , k_Bragg_si*1e-6],[0 20 ],'Color','black','LineStyle','--')
%     
%     xlim([-30 30])
%     ylim([0 15])
%     
%     title(num2str(i))
% end


%%
% figure
% for i = 1:length(t_array)
%     index_exp = floor(i/20) + 1;
%     
%     theta = double(theta_t{i});
%     nk_rapid  = LLS.calcMomentumDistr(theta_t{i}, 0);
%     
%     subplot(2,1,1)
%     imagesc(x_array*1e6,k_array*1e-6,squeeze(theta))
%     colormap(hot)
%     caxis([ 0 0.5])
%     ylabel('rapidity')
%     set(gca,'YDir','normal') 
%     
%     subplot(2,1,2)
%     hold on
%     box on
%     plot(k_arr_exp*1e-6 , nk_exp(index_exp,:)* 1e6 ,'k-', 'LineWidth', 1.5)
%     plot(k_array*1e-6, nk_rapid*1e6, 'LineWidth', 1.5)
%     ylim([0 20])
%     xlim([-20 20])
%     xlabel('k')
%     ylabel('density')
%     
%     suptitle(['t = ' num2str(t_array(i)*1e3 , 2) 'ms'])
%     pause(0.05)
% end


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

