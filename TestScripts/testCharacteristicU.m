% Test whether the characteristic u-function is correctly calculated, i.e.
% that theta_t(x, lambda) = theta_0( u( x,t,lambda), lambda)
% NOTE: assumes homogeneous evolution


clear all; close all;

addpath('..\Functions\')

%% define parameters

N       = 2^8;
M       = 2^7;
Ntypes  = 1;

kmax    = 5;
xmax    = 7;

T       = 5;

dt      = 0.01;


k_array     = linspace(-kmax,kmax,N);
x_array     = linspace(-xmax,xmax,M);
t_array     = 0:dt:0.4;

%% Setup state
couplings   = { @(t,x) 4, @(t,x) 1 };
couplings_i = couplings;
couplings_i{1} = @(t,x) 4 - 4*(x-2).^2;


LLS         = LiebLinigerSolver(x_array, k_array, couplings);
theta_init  = LLS.calcThermalState(T, couplings_i);


[theta_t, u_t] = LLS.propagateTheta(theta_init, t_array);

%% Test u - should follow theta(x,t,lambda) = theta( u(x,t,lambda), 0, lambda)
Nsamples    = 4;
t_idx       = floor(linspace(1,length(t_array),Nsamples));

for i = 1:Nsamples
    for j = 1:Ntypes
        the_t   = theta_t{t_idx(i)};
        the_t   = squeeze( the_t(:,:,j,:,:) );

        uj_t    = u_t{t_idx(i)};
        uj_t    = squeeze( uj_t(:,:,j,:,:) );
        
        u_int   = permute( uj_t, [2 1]); % spatial index first
        the_int = permute(theta_init(:,:,j,:,:), [5 1 3 2 4]); % spatial index first
        r_int   = repmat(k_array, M, 1);

        % interpolates theta(x) to u(t,x,lambda)
        theta_u = interp2( k_array, x_array, the_int, r_int(:), u_int(:), 'spline');
        theta_u = reshape(theta_u, M, N);
        theta_u = permute(theta_u, [2 1]);

        figure
        subplot(2,2,1)
        imagesc(x_array, k_array, the_t )
        set(gca,'YDir','normal') 
        title('theta( x, t)')

        subplot(2,2,2)
        imagesc(x_array, k_array,  theta_u )
        set(gca,'YDir','normal') 
        title('theta( u, 0)')

        subplot(2,2,3)
        imagesc(x_array, k_array, uj_t )
        set(gca,'YDir','normal') 
        title('u( x, t)')

        subplot(2,2,4)
        imagesc(x_array, k_array, the_t - theta_u )
        set(gca,'YDir','normal') 
        title('difference')
        colorbar

        suptitle(['t = ' num2str(t_array(t_idx(i))) ', j = ' num2str(j)])
    end
end