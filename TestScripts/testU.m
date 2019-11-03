% Test whether the characteristic u-function is correctly calculated, i.e.
% that theta_t(x, lambda) = theta_0( u( x,t,lambda), lambda)
% NOTE: assumes homogeneous evolution


clear all; close all;

addpath(['..' filesep 'Functions' filesep])

%% define parameters

N       = 2^7;
M       = 2^7;
Ntypes  = 1;

kmax    = 10;
xmax    = 10;

T       = 4;

dt      = 0.01;


k_array     = linspace(-kmax,kmax,N);
x_array     = linspace(-xmax,xmax,M);
t_array     = 0:dt:1;

%% Setup state
c = 4;
b = 5/c^4;
a = 2*b*c^2;

couplings   = { @(t,x) 5 - 1*x.^2, @(t,x) 1 };
couplings_i = { @(t,x) a*x.^2 - b*x.^4 , @(t,x) 1 };


LLS         = LiebLinigerSolver(x_array, k_array, couplings);
theta_init  = LLS.calcThermalState(T, couplings_i);


[theta_t, u_t, w_t] = LLS.propagateTheta(theta_init, t_array);

%% Test u - should follow theta(x,t,lambda) = theta( u(x,t,lambda), 0, lambda)
Nsamples    = 4;
t_idx       = floor(linspace(1,length(t_array),Nsamples));


for i = 1:Nsamples
    for j = 1:Ntypes
        the_t   = theta_t{t_idx(i)};
        the_t   = squeeze( the_t(:,:,j,:,:) );

        uj_t    = u_t{t_idx(i)};
        uj_t    = squeeze( uj_t(:,:,j,:,:) );
        
        wj_t    = w_t{t_idx(i)};
        wj_t    = squeeze( wj_t(:,:,j,:,:) );
        
        u_int   = permute( uj_t, [2 1]); % spatial index first
        w_int   = permute( wj_t, [2 1]); % spatial index first
        the_int = permute(theta_init(:,:,j,:,:), [5 1 3 2 4]); % spatial index first

        % interpolates theta(x) to u(t,x,lambda)
        theta_u = interp2( k_array, x_array, the_int, w_int(:), u_int(:), 'spline');
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

%         suptitle(['t = ' num2str(t_array(t_idx(i))) ', j = ' num2str(j)])
    end
end


%% trajectory plot
Umat = zeros(N,1,Ntypes,1,M,length(t_array));
Wmat = zeros(N,1,Ntypes,1,M,length(t_array));

for i = 1:length(t_array)
    Umat(:,:,:,:,:,i) = u_t{i};
    Wmat(:,:,:,:,:,i) = w_t{i};
end




rapids = 10:10:60;

for j = 1:Ntypes
    figure
    hold on
    plot(t_array, squeeze(Umat(rapids,1,1,1,48,:)) ,'-')
    plot(t_array, x_array(48) - 2*t_array.*k_array(rapids)' ,'k--')
    
end


positions = 10:10:60;

for j = 1:Ntypes
    figure
    hold on
    plot(t_array, squeeze(Umat(40,1,1,1,positions,:)) ,'-')
    plot(t_array, x_array(positions)' - 2*t_array.*k_array(40)' ,'k--')
    
end



%%

endpoints = [ 40 64; 70, 70; 100, 100];

for j = 1:Ntypes
    
    the_t   = theta_t{end};
    the_t   = squeeze( the_t(:,:,j,:,:) );

    
    col     = t_array;
    
    figure

    ax1 = axes;
    imagesc(x_array, k_array, the_t )
    set(gca,'YDir','normal') 
    
    xlabel('Space')
    ylabel('Rapidity')

    xlim([-xmax , xmax])
    ylim([-kmax , kmax])
    
    
    ax2 = axes;
    hold on
    
    for i = 1:size(endpoints,1)
        x       = flip(squeeze(Umat(endpoints(i,1),1,1,1,endpoints(i,2),:)));
        lambda  = flip(squeeze(Wmat(endpoints(i,1),1,1,1,endpoints(i,2),:)));

        scatter(x(:),lambda(:),[],col,'filled')
    end
    
    
    %Link them together
    linkaxes([ax1,ax2])
    %Hide the top axes
    ax2.Visible = 'off';
    ax2.XTick = [];
    ax2.YTick = [];
    %Give each one its own colormap
    colormap(ax1,'hot')
    colormap(ax2,'cool')
    
    
    c = colorbar;
    c.Label.String = 'Time';
    c.Label.FontSize = 12;
    
    pos = ax2.Position;
    set([ax1,ax2],'Position',pos);
end

