clear all; close; all clc;

data_init = load('initial_data');
data_t = load('MDF_N120');


dk = 4*pi/(852*1e-9)/52;
k_arr_exp = 0:dk:(length(data_init.N120_after) - 1)*dk;
k_arr_exp = k_arr_exp - mean(k_arr_exp);

nk_init_exp = 120*data_init.N120_after/trapz(k_arr_exp, data_init.N120_after);

nk_exp = 120*data_t.MDF_full/trapz(k_arr_exp, data_init.N120_after);

t_arr_exp = data_t.time;


k_Bragg_si  = 4*pi/(852*1e-9); % [m^-1]
figure
for i = 1:180
    subaxis(15,12,i,'SpacingHoriz',0.001,'SpacingVert',0.001)
    hold on
    box on
    plot(k_arr_exp*1e-6 , nk_exp(i,:)*1e6, 'k')
    
    
    line([-k_Bragg_si*1e-6  , -k_Bragg_si*1e-6],[0 20 ],'Color','black','LineStyle','--')
    line([k_Bragg_si*1e-6 , k_Bragg_si*1e-6],[0 20 ],'Color','black','LineStyle','--')
    
    xlim([-30 30])
    ylim([0 15])
    set(gca,'xticklabel',[])
    set(gca,'yticklabel',[])
   
    
end