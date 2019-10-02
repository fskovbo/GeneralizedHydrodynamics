% Test whether the characteristic u-function is correctly calculated, i.e.
% that theta_t(x, lambda) = theta_0( u( x,t,lambda), lambda)
% NOTE: assumes homogeneous evolution


clear all; close all;

%% define parameters

N       = 2^7;
M       = 2^6;

kmax    = 5;
xmax    = 7;

T       = 5;

stepOrder   = 2;
dt      = 0.025;


k_array     = linspace(-kmax,kmax,N);
x_array     = linspace(-xmax,xmax,M);
t_array     = 0:dt:0.4;

%% Setup state

mu          = @(t,x) 4 - 4*heaviside(-t+eps)*x.^2; 
dmudx       = @(t,x) 0;
c           = @(t,x) 1;
dcdt        = [];
dcdx        = [];

couplings.mu    = mu;
couplings.dmudx = dmudx;
couplings.c     = c;
couplings.dcdt  = dcdt;
couplings.dcdx  = dcdx;


LLS         = LiebLinigerSolver(x_array, k_array, couplings, stepOrder);
theta_init  = LLS.calcThermalState(T);


[theta_t, u_t] = LLS.propagateTheta(theta_init, t_array);

%% Test u - should follow theta(x,t,lambda) = theta( u(x,t,lambda), 0, lambda)
Nsamples    = 4;
t_idx       = floor(linspace(1,length(t_array),Nsamples));

for i = 1:Nsamples    
    u_int = permute(u_t(:,:,:,t_idx(i)), [3 2 1]);
    r_int = repmat(k_array, M, 1);

    % interpolates theta(x) to u(t,x,lambda)
    theta_u = interp2( k_array, x_array, permute(theta_init, [3 2 1]), r_int(:), u_int(:), 'spline');
    theta_u = reshape(theta_u, M, N);
    theta_u = permute(theta_u, [3 2 1]);
    
    figure
    subplot(2,2,1)
    imagesc(x_array, k_array, squeeze( theta_t(:,:,:,t_idx(i)) ))
    set(gca,'YDir','normal') 
    title('theta( x, t)')
    
    subplot(2,2,2)
    imagesc(x_array, k_array, squeeze( theta_u ))
    set(gca,'YDir','normal') 
    title('theta( u, 0)')
    
    subplot(2,2,3)
    imagesc(x_array, k_array, squeeze(u_t(:,:,:,t_idx(i))))
    set(gca,'YDir','normal') 
    title('u( x, t)')
    
    subplot(2,2,4)
    imagesc(x_array, k_array, squeeze( theta_t(:,:,:,t_idx(i)) - theta_u  ))
    set(gca,'YDir','normal') 
    title('difference')
    colorbar
    
    suptitle(['t = ' num2str(t_array(t_idx(i)))])
end