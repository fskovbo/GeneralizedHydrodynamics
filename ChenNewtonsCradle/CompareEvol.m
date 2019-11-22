clear all; close all; clc;


%% Load experimental data
exp_data_init = load('initial_data');
exp_data_t = load('MDF_N120');


dk = 4*pi/(852*1e-9)/52;
k_arr_exp = 0:dk:(length(exp_data_init.N120_after) - 1)*dk;
k_arr_exp = k_arr_exp - mean(k_arr_exp);

nk_init_exp = 120*exp_data_init.N120_after/trapz(k_arr_exp, exp_data_init.N120_after);

nk_exp = permute(exp_data_t.MDF_full, [2 1]);
nk_exp = 120*nk_exp./trapz(k_arr_exp, nk_exp);


t_arr_exp = exp_data_t.time;


%% Load GHD data
disp('Loading simulation data ...')

% GHD_data_t1 = load('GHD_T100_c1'); % (c = 1)
% GHD_data_t2 = load('GHD_T100_c01'); % (c = 0.1)
% GHD_data_t3 = load('GHD_T100_c1mid2'); % (c = 1, stuff in middle)
% GHD_data_t4 = load('GHD_T100_c1mid3'); % (c = 1, less stuff in middle)

GHD_data_t1 = load('GHD_T100_c1mid6'); % (c = 1, less stuff in middle)
GHD_data_t2 = load('GHD_T100_c1mid7'); % (c = 1, less stuff in middle)

nk_ghd(:,:,1) = GHD_data_t1.save_data.nk_rapid;
nk_ghd(:,:,2) = GHD_data_t2.save_data.nk_rapid;
% nk_ghd(:,:,3) = GHD_data_t3.save_data.nk_rapid;
% nk_ghd(:,:,4) = GHD_data_t4.save_data.nk_rapid;


k_array(:,:,1) = GHD_data_t1.save_data.k_array;
k_array(:,:,2) = GHD_data_t2.save_data.k_array;
% k_array(:,:,3) = GHD_data_t3.save_data.k_array;
% k_array(:,:,4) = GHD_data_t4.save_data.k_array;


t_array = GHD_data_t1.save_data.t_array;

Nsets = size(k_array,3);
Nsteps = length(t_array);
N = length(k_array);

disp('... finished loading!')

%% Calculate average of 120 shots ~ 1 period (GHD)
avg_ghd = zeros( N , Nsteps-119, Nsets);

avg_t_ghd = zeros(1,Nsteps-119);

for i = 1:Nsteps-119
    range = i:(i+119);
    
    avg_t_ghd(i) = 1e3*mean(t_array(range));
    
    for j = 1:Nsets
        avg_ghd(:,i,j) = mean( nk_ghd(:, range, j) , 2 );
    end    
end


%% Calculate average of 12 shots ~ 1 period (Exp)
avg_exp = zeros( size(nk_exp) );
avg_t_exp = zeros(1,size(nk_exp,2));

% average first period (first 60 shots)
avg_exp(:, 1) = mean( nk_exp(:, 1:60) ,2 );
avg_t_exp(1) = mean( t_arr_exp(1:60) );


% average next 44 periods
ii = 1;
for i = 1:(44*12)
    range = 60 + (i:(i+11));
    
    t_diff = t_arr_exp(range(end)) - t_arr_exp(range(1));
    if t_diff == 11
        ii = ii + 1;
        avg_exp(:,ii) = mean( nk_exp(:, range) , 2 );
        avg_t_exp(ii) = mean( t_arr_exp(range) );
    else
        disp(['Warning: Detected jump of dt=' num2str(t_diff) ' at i=' num2str(i)])
    end
end

avg_exp = avg_exp(:,1:ii);
avg_t_exp = avg_t_exp(:,1:ii);



%% Fit experimental data with Gaussian
residuals_exp = zeros(1, length(avg_t_exp));

for i = 1:length(avg_t_exp)
    gfit = fit(k_arr_exp'*1e-6, avg_exp(:, i) *1e6 ,'gauss1');
    fitcurve = feval(gfit, k_arr_exp*1e-6);
    
    residuals_exp(i) = sum( (avg_exp(:, i) *1e6 - fitcurve).^2 );
end


%% Fit simulation data with Gaussian
residuals_ghd = zeros(Nsets, length(avg_t_ghd));

options = fitoptions('gauss1');
options.StartPoint = [3.5 , 0 , 20];

for i = 1:length(avg_t_ghd)
    for j = 1:size(avg_ghd,3)
        kmax        = k_array(1,end,j);
        dk          = k_array(1,2,j)-k_array(1,1,j);
        k_arr_fit   = -2*kmax:dk:2*kmax;
        
        avg_fit     = padarray(avg_ghd(:,i,j), N/2 -1 );
        gfit        = fit(k_arr_fit'*1e-6, avg_fit*1e6, 'gauss1', options);
        fitcurve    = feval(gfit, k_arr_fit*1e-6);
    
        residuals_ghd(j,i) = sum( (avg_fit*1e6 - fitcurve).^2 );
    end
end

%% Plot residuals

figure
hold on
box on
plot(avg_t_exp,residuals_exp,'ko')
plot(avg_t_ghd,residuals_ghd,'-')
set(gca,'yscale','log')
xlabel('t')
ylabel('deviation from Gaussian')


%% Calculate difference between experiment and simulation averages
sqdist = zeros( Nsets, length(avg_t_exp) );
sqdist_err = zeros( Nsets, length(avg_t_exp) , 2);


for i = 1:length(avg_t_exp)
    % Find closest matching time
    [~,index] = min(abs( avg_t_exp(i) - avg_t_ghd)); 
    
    for j = 1:Nsets
        % Interpolate experimental data so difference can be taken
        nke = interp1(k_arr_exp, avg_exp(:,i), k_array(:,:,j));
        
        sqdist(j,i) = sum( (avg_ghd(:,index,j)*1e6 - nke' *1e6).^2 );
        
        % Caculate "errorbars" given by the biggest differences within 10
        % timesteps
        range = (index-5):(index+5);
        range( range < 1) = 1;
        range( range > length(avg_t_ghd) ) = length(avg_t_ghd);
        
        diffs = sum( (avg_ghd(:,range,j)*1e6 - nke' *1e6).^2, 1 );
        sqdist_err(j,i,1) = min(diffs);
        sqdist_err(j,i,2) = max(diffs);
    end
end


%% Plot differences 
figure
hold on
box on

for j = 1:Nsets
    errorbar(avg_t_exp, sqdist(j,:), sqdist(j,:)-sqdist_err(j,:,1), sqdist_err(j,:,2)-sqdist(j,:)  )
end
xlabel('t')
ylabel('difference between averages')


%% Plot stability of simulation
atomnumber = zeros(Nsets, Nsteps);

for i = 1:Nsteps
   for j = 1:Nsets
       atomnumber(j,i) = trapz(k_array(:,:,j), nk_ghd(:,i,j) );
   end
end

figure
hold on
box on
plot(t_array, atomnumber)
xlabel('t')
ylabel('atomnumber')


%% plot last 30 periods of mid
figure(68)
for i = 1:180
    index_ghd = length(t_array) + 20*(i - 360);
    
    subaxis(15,12,i,'SpacingHoriz',0.001,'SpacingVert',0.001) 
    imagesc( squeeze(double(GHD_data_t4.save_data.x_theta_t{index_ghd})) )
    caxis([0 1])
    colormap(hot)
    set(gca,'xticklabel',[])
    set(gca,'yticklabel',[])
    set(gca,'YDir','normal')
    
end


%% Plot period average of experiment vs simulation
k_Bragg_si  = 4*pi/(852*1e-9); % [m^-1]
figure
for i = 1:45        
    ii = 10*(i-1) + 1;
    [~,index] = min(abs( avg_t_exp(ii) - avg_t_ghd));
    
    subplot(5,9,i)
    hold on
    box on
    
    for j = 1:Nsets
        plot(k_array(:,:,j)*1e-6,avg_ghd(:,index,j)*1e6,'LineWidth',1.5)
    end
    
    plot(k_arr_exp*1e-6,avg_exp(:, ii)*1e6,'k','LineWidth',1.5)
    
    line([-k_Bragg_si*1e-6  , -k_Bragg_si*1e-6],[0 20 ],'Color','black','LineStyle','--')
    line([k_Bragg_si*1e-6 , k_Bragg_si*1e-6],[0 20 ],'Color','black','LineStyle','--')
    
    ylim([0 5])
    xlim([-30 30])
    set(gca,'xticklabel',[])
    set(gca,'yticklabel',[])
    
    title(['t=' num2str(avg_t_exp(ii))])
    
end


%% Test Gaussian fit of simulation at specific times
% hej = figure(69);
% 
% for i = (length(avg_t_ghd) - 20):length(avg_t_ghd)
%     for j = 1
%         clf(hej)
%         hold on
%         box on
%         
%         kmax        = k_array(1,end,j);
%         dk          = k_array(1,2,j)-k_array(1,1,j);
%         k_arr_fit   = -2*kmax:dk:2*kmax;
%         
%         avg_fit     = padarray(avg_ghd(:,i,j), N/2 -1 );
%         gfit        = fit(k_arr_fit'*1e-6, avg_fit*1e6, 'gauss1', options);
%         fitcurve    = feval(gfit, k_arr_fit*1e-6);
%     
%         plot( k_arr_fit*1e-6, avg_fit*1e6 )
%         plot( k_arr_fit*1e-6, fitcurve )
%         title(['res = ' num2str(sum( (avg_fit *1e6 - fitcurve).^2 ))])
%     end
% end
