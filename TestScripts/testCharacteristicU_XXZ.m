% Test whether the characteristic u-function is correctly calculated, i.e.
% that theta_t(x, lambda) = theta_0( u( x,t,lambda), lambda)
% NOTE: assumes homogeneous evolution


clear all; close all;

addpath('..\Functions\')

%% Define simulation parameters

N           = 2^7;
M           = 2^8;
Ntypes      = 2;
dt          = 0.01;

kmax        = pi/2;
xmax        = 1;
tmax        = 1;

k_array     = linspace(-kmax, kmax, N);
x_array     = linspace(-xmax, xmax, M);
t_array     = linspace(0, tmax, tmax/dt+2);


T           = 0.25;

%% Define physical couplings
couplings       = { @(t,x) -1 , @(t,x) acosh(1.5) };

coup_init       = couplings;
coup_init{1}    = @(t,x) -1 - 8*x.^2; 

%% Run simulation
XXZ         = XXZchainSolver(x_array, k_array, couplings, Ntypes);
theta_init  = XXZ.calcThermalState(T, coup_init);


[theta_t, u_t] = XXZ.propagateTheta(theta_init, t_array);

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


%% play movie

figure
for i = 1:length(t_array)
    dist = theta_t{i};
    
    subplot(1,2,1)
    imagesc(x_array,k_array,squeeze(dist(:,:,1,:,:)))
    colormap(hot)
    caxis([ 0 1])
    ylabel('rapidity')
    set(gca,'YDir','normal') 
    
    subplot(1,2,2)
    imagesc(x_array,k_array,squeeze(dist(:,:,2,:,:)))
    colormap(hot)
    caxis([ 0 0.02])
    ylabel('rapidity')
    set(gca,'YDir','normal') 
    
    suptitle(['t = ' num2str(t_array(i))])
    pause(0.05)
end