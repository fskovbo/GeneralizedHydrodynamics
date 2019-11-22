clear all; close all;

addpath([pwd '\Functions\'])


%% Simulation parameters

N           = 2^7;
M           = 2^8;

Nsteps      = 2500;
stepOrder   = 2;

xmax_si     = 110*1e-6;
kmax_si     = 2*1e7;
tmax_si     = 200*1e-3;

%% Experimental units

m_si        = 87*1.6605402e-27;
hbar_si     = 1.054571726e-34;
kB_si       = 1.38065e-23;
as_si       = 5.2e-9;

Natoms      = 1000; 
omegaT_si   = 2*pi*1600;
boxL_si     = 200*1e-6;
boxH_si     = 3e-30;
barrierW_si = 0.5*1e-6;
barrierH_si = 1e-29;
barrierT_si = 25*1e-3;
PSF_si      = 3 * 1e-6;

g1D_si      = 2*hbar_si*omegaT_si*as_si;
c_si        = g1D_si*m_si/hbar_si^2;


%% create arrays
x_array     = linspace(-xmax_si,xmax_si,M);
k_array     = linspace(-kmax_si,kmax_si,N)/2;
t_array     = linspace(0, tmax_si, Nsteps+1);

%% Define physical couplings and temperature
boxTrap     = @(t,x) boxH_si*(-tanh(1/PSF_si*(x+boxL_si/2)) + tanh(1/PSF_si*(x-boxL_si/2)) + 2);
barrier     = @(t,x) (1 - tanh(2/barrierT_si * t)) * barrierH_si*(tanh(2/PSF_si*(x+barrierW_si/2)) - tanh(2/PSF_si*(x-barrierW_si/2)) );

dbox        = @(t,x) boxH_si/PSF_si*(-sech(1/PSF_si*(x+boxL_si/2)).^2 + sech(1/PSF_si*(x-boxL_si/2)).^2);
dbarrier    = @(t,x) (1 - tanh(2/barrierT_si * t)) * 2*barrierH_si/PSF_si*(sech(2/PSF_si*(x+barrierW_si/2)).^2 - sech(2/PSF_si*(x-barrierW_si/2)).^2 );


couplings   = { @(t,x) (boxTrap(t,x) + barrier(t,x)) , @(t,x) c_si   ; 
                []                                   , []            ; 
                @(t,x) -(dbox(t,x) + dbarrier(t,x))  , []            };

figure
hold on
plot(x_array, boxTrap(0, x_array))
plot(x_array, barrier(0, x_array))

figure
hold on
plot(x_array, dbox(0, x_array))
plot(x_array, dbarrier(0, x_array))


T_si        = 50*1e-9;   


%% Initialize solver and generate initial state

LLS         = LiebLinigerSolver_SI(omegaT_si, x_array, k_array, couplings);
mu0_fit     = LLS.fitAtomnumber(T_si,  [], Natoms, true);
theta_init  = LLS.calcThermalState(T_si);

%%
n_0         = LLS.calcCharges(theta_init, 0 ,0);

figure
subplot(2,1,1)
imagesc( squeeze( double(theta_init)) )
subplot(2,1,2)
plot(n_0)


%% Solve dynamics and calculate quantities
theta_t     = LLS.propagateTheta(theta_init, t_array);
n_t         = LLS.calcCharges(theta_t, 0, t_array) * 1e-6; % [um^-1]
rho_t       = LLS.transform2rho( theta_t, t_array); % [um um^-1]

%% Plot results

figure
for i = 1:10
    index   = ceil( (length(t_array)-1)/9*(i-1) + 1 );
    
    theta   = theta_t{index};
    
    % plot Harmonic quasiparticle distribution
    subplot(2,10,i)
    imagesc(x_array*1e6, k_array*1e-6,squeeze(double(theta)) )
    set(gca,'YDir','normal')
    colormap(hot)
    caxis([0 1])
    ylabel('k [um^-1]')
    set(gca,'xticklabel',[])
    title(['t = ' num2str(t_array(index)*1e3) 'ms'])

    % plot density profiles in real space
    subplot(2,10,i+10)
    hold on
    box on
    plot(x_array*1e6, n_t(:,index))
    xlim([ -xmax_si*1e6, xmax_si*1e6])
    ylim([0 1.1*max(n_t(:))])
    xlabel('x [um]')
    ylabel('n [um^-1]')
    
end


%% -----------------------------------------------------------------------

function out = convPSF(x_array, func)
    PSF     = 2.1 * 1e-6;
    gauss   =  exp( -x_array.^2 / (2*PSF^2) );
    out     = conv(func,gauss,'same');
end


function ddx = getSpatialDeriv(x_array, func)
    x_fine      = linspace( x_array(1), x_array(end), 8*(length(x_array)-1)+1  );
    [~,idx]     = ismembertol( x_array, x_fine, 1e-8);
    
    func_fine   = interp1( x_array, func, x_fine, 'spline' );
    deriv       = gradient(x_fine, func_fine);
    figure
    plot(x_fine, func_fine)
    figure
    plot(x_fine, deriv)
    ddx         = deriv(idx);
end
