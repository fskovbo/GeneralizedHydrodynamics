clear all; close all;


addpath(['..' filesep 'Functions' filesep])

%% Define simulation parameters

N           = 2^7;
M           = 2^7;
dt          = 0.025;

kmax        = 3;
xmax        = 10;
tmax        = 6;

% k_array     = linspace(-kmax, kmax, N);
% kw = k_array(2) - k_array(1);
[k_array,kw]=legzo(N, -kmax, kmax);
x_array     = linspace(-xmax, xmax, M);
t_array     = linspace(0, tmax, tmax/dt+1);

options.extrapFlag = true;

%% Define physical couplings and temperature
% couplings  = { @(t,x) 1 - 1*x.^2 , @(t,x) 0.4 + 1.5*tanh(3*t) };
% T           = 0.5;

couplings  = { @(t,x) 0 , @(t,x) 1 };
T       = @(x) 1.5*heaviside(x) + 2*heaviside(-x);

%% Initialize state and solve dynamics

% Harmonic potential
tic
shG        = sinhGordonSolver(x_array, k_array, kw, couplings, options);
theta_init = shG.calcThermalState(T);
theta_t    = shG.propagateTheta(theta_init, t_array);
n_t        = shG.calcCharges(theta_t, 0, t_array);
toc


%% ------------ Plot results -------------------

% Initial state plot
figure
hold on
subplot(2,1,1)
imagesc(x_array,k_array, squeeze(double(theta_init)) )
colormap(hot)
caxis([ 0 1])
set(gca,'YDir','normal') 


subplot(2,1,2)
plot(x_array, n_t(:,1))
xlabel('x')
ylabel('n')
ylim([0 1])
xlim([-xmax, xmax])


xlim([-xmax, xmax])


% Stability plot
figure
hold on
box on
plot(t_array, sum(n_t,1))
xlabel('t')
ylabel('N')
xlim([0 tmax])

% find indices of specific times
plottimes       = [0, 0.1, 2, 6];
[~, t_idx]      = ismember(plottimes, t_array);



%% Harmonic Pontential
figure
for i = 1:4
    subplot(2,4, i)
    index = t_idx(i);
    imagesc(x_array,k_array, squeeze(double(theta_t{index})) ) 
    colormap(hot)
    caxis([ 0 1])
    set(gca,'YDir','normal') 
    title(['t = ' num2str(t_array(index))])
end

subplot(2,4,[5 6])
hold on
box on
for i = 1:4
    plot(x_array, n_t(:,t_idx(i)))
end
xlabel('x')
ylabel('n')
ylim([0 1])
xlim([-xmax, xmax])

subplot(2,4,[7 8])
box on
plot(t_array, n_t(M/2, :))
xlabel('t')
ylabel('n')
ylim([0 1])
xlim([0 tmax])
