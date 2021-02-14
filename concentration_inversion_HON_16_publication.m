%particle diameter                                        a    [m]
%the first isonification frequency                        f1    [Hz]
%the second isonification frequency                       f2    [Hz]
%first frequency designation (high,middle,low->1,2,3)     fb1   [integer]
%second frequency designation (high,middle,low->1,2,3)    fb2   [integer]
%sound velocity in water                                  c    [m/s]
%the particle density                                     rho_solid  [kg/m^3]
%particle scattering constant (theoretical)               ks   [-]
%the transducer constant                                  kt   [-]
%Fluid Temperature                                        T    [Celcius]
%Attenuation due to the particle                                alpha_s
%Outputs water attenuation constant                             alpha_w   [-]
%first Particle concentration from standard inversion           m1     [g/l]
%second Particle concentration from standard inversion          m2     [g/l]
%particle concentration form dual frequency inversion           m_dual [g/l]
%transducer radius                                        at    [m]

clear all
close all

  FREQ_RANGE = 1:3; %%select frequency range (3rd frequency excluded as fit found to be poor)
     CHANNEL_RANGE_DEFINED= [1 2 3 4 5 6 8]; %%select channel range 
     sample_range_DEFINED = 2:9;  %choose sample value range
        

     
     
     
     
              
RMS_DATA_V= 'RMS_DATA_V_HON_16.mat' ;  %SELECT .mat DATA FILE TO LOAD FROM
load(RMS_DATA_V)
load 'atten_coeff_hon_16.mat';
A= RMS_DATA_V_HON_16;            %select RMS data file loaded from data file       
atten_coeff_std=std(atten_coeff,0,1)./mean(atten_coeff,1);
a= 0.00003929;
f= [2000000, 2250000, 2500000];
at =   0.00318; 
c= 1474;
rho_solid = 2469.6;
x= 2*3.14.*f.*a/c;
T= 17;
alpha_w= 0.05641.*((f./1000000).^2).*exp(-T/27);

%particle scattering coefficient estimation%
sigma= (1-0.5.*exp(-(((x-1.5)./0.5).^2))) .* (1+0.4.*exp(-(((x-1.5)./3).^2))) .* (1-0.5.*exp(-(((x-5.9)./0.7).^2)));

form_function = (sigma.*(x.^2))./(1.17+0.95.*(x.^2));

ks= form_function./(sqrt(a.*rho_solid));

%formatting sample data dimensions for streamlining%
sample_range = zeros(max(FREQ_RANGE),max(CHANNEL_RANGE_DEFINED),length(sample_range_DEFINED));
     for int1 = FREQ_RANGE
         CHANNEL_RANGE(int1,:) = CHANNEL_RANGE_DEFINED;
         for int2 = CHANNEL_RANGE_DEFINED
             
          sample_range(int1,int2,:) = sample_range_DEFINED;       
                 
             
         end
     end
    




%setting up marker string for plots%
markers = [ "+" "^" "<" "v" "o" "s" "d" "*" "x" "+"];
atten_coeff = mean(atten_coeff,1);
correction=1;
CHANNEL_STRING="";

%%%%%%%%kt calcualtion starts here%%%%%%%%

for int1 = 1:9;

	for int2 = 1:8;

		for int3 = 1:3;

alpha_s = atten_coeff(int3)*sample_values(int1);
kt(int1, int2, int3, :) = ( squeeze((RMS_DATA_V_HON_16(int1, 1, int2, int3, :))).*(squeeze(r(int1, 1, int2, int3, :))) ) ./ (  ( sample_values(int1)^(1/2) ) .* ( exp( (-2*squeeze(r(int1, 1, int2, int3, :))).*(alpha_w(int3) + alpha_s) ) * ks(int3)));

		end
	end
end


%%%%%%Averageing starts here%%%%%%
%%%% ks concentration average is ks averaged across all probes, select concentrations and over a set distance range. Then seperated into each frequency %%%%
kt_refined = kt (2:3, :, :, 100:200);

kt_distance_average = mean(kt_refined, 4);

kt_conc_average = squeeze(mean(kt_distance_average, 1));

kt_probe_average = mean(kt_conc_average, 1);

%%Single Frequency concentration inversion%%

for int2 = FREQ_RANGE;   %%select frequency range
for CHANNEL_NUMBER = CHANNEL_RANGE(FREQ_RANGE(1),:);




fig = figure;
legend_direct=num2str( (sample_values(sample_range(FREQ_RANGE(1), CHANNEL_RANGE(FREQ_RANGE(1),1), :)))' ); 

    for int1= sample_range(FREQ_RANGE(1), CHANNEL_RANGE(FREQ_RANGE(1),1), :);

        
m_nominal=sample_values;
alpha_s(int1, int2) = m_nominal(int1)*atten_coeff(int2);   

%storing ks profiles for further analysis%
ks_alternate_stored(int1, 1, CHANNEL_NUMBER, int2, 1:312) = ( A(int1, 1, CHANNEL_NUMBER, int2, :) .* r(int1, 1, 1, int2, :) .* ( m_nominal(int1).^(-0.5) ) .* exp(2.*r(int1, 1, 1, int2, :).* (alpha_s(int1, int2) + alpha_w(int2)) ) ) ./ kt_conc_average(CHANNEL_NUMBER, int2) ;

ks_alternate = mean(ks_alternate_stored(int1, 1, CHANNEL_NUMBER, int2, 100:200), 5);
ks_alternate_stored_mean(int1, 1, CHANNEL_NUMBER, int2) = ks_alternate;

m1(int1, 1, CHANNEL_NUMBER, int2, :)= ( ( squeeze(A(int1, 1, CHANNEL_NUMBER, int2, :))...
    .*squeeze(r(int1, 1, 1, int2, :))...
    .*squeeze(Phi_exp(int1, 1, CHANNEL_NUMBER, int2, :))... 
    ./(ks_alternate*kt_conc_average(CHANNEL_NUMBER, int2)) ).^2)...
    .*exp(4.*squeeze(r(int1, 1, 1, int2, :))...
    .*(alpha_w(int2) + alpha_s(int1, int2))) ;  

hold on
plot(squeeze(m1(int1, 1, CHANNEL_NUMBER, int2, :)),squeeze(r(int1, 1, 1, int2, :)),char(markers(int1)));
CHANNEL_STRING(int1-1) = strjoin( {char(num2str(sample_values(int1),'%.1f')) 'g.l^{-1}'});

    end
    
h_legend=legend(CHANNEL_STRING');
title(['Single frequency profile at ' num2str(f(int2)/1000000) 'MHz without correction factor (Honite 16, probe' num2str(CHANNEL_NUMBER) ')'])
axis([0 170 0 0.3])
xlabel('Concentration (g/l)')
ylabel('r (m)')
%  saveas(fig, ['Single frequency profile using ' num2str(f(int2)) 'Hz for Honite 16 without correction factor for Honite 16 (probe' num2str(CHANNEL_NUMBER) ')'], 'png')

end    

end

%%Dual frequency concentration inversion%%
for int1 = 2:9
   
    for int2 = CHANNEL_RANGE_DEFINED;
       for int3 = FREQ_RANGE;
phi(int1,1,int2,int3,:)= sqrt( (((ks_alternate_stored_mean(int1, 1, int2, int3))*kt_conc_average(int2,int3)./(r(1,1,1,1,:).*Phi_exp(int1,1,int2,int3,:))).^2) .* exp(-4.*(r(1,1,1,1,:).*alpha_w(int3))) );
       end
    end
end

J = (RMS_DATA_V_HON_16.^2)./(phi.^2);
for int2 = CHANNEL_RANGE_DEFINED;
    for int3 = FREQ_RANGE;
        
        for int1 = 2:9
m_dual_1_2(int1,1,int2,1,:)= (J(int1,1,int2,1,:).^((1-(atten_coeff(1)/atten_coeff(2)))^-1)) .* (J(int1,1,int2,2,:).^((1-(atten_coeff(2)/atten_coeff(1)))^-1));
m_dual_1_3(int1,1,int2,1,:)= (J(int1,1,int2,1,:).^((1-(atten_coeff(1)/atten_coeff(3)))^-1)) .* (J(int1,1,int2,3,:).^((1-(atten_coeff(3)/atten_coeff(1)))^-1));
m_dual_2_3(int1,1,int2,1,:)= (J(int1,1,int2,2,:).^((1-(atten_coeff(2)/atten_coeff(3)))^-1)) .* (J(int1,1,int2,3,:).^((1-(atten_coeff(3)/atten_coeff(2)))^-1));
        end
    end
end

    for int2 = CHANNEL_RANGE_DEFINED;
        
        figure
        hold on
        legend_direct=num2str( (sample_values(sample_range(FREQ_RANGE(1), CHANNEL_RANGE(FREQ_RANGE(1),1), :)))' ); 

        for int1 = 2:9;
          plot(squeeze(m_dual_1_2(int1, 1, int2, 1, :)),squeeze(r(int1, 1, 1, 1, :)), char(markers(int1)));
             CHANNEL_STRING(int1-1) = strjoin( {char(num2str(sample_values(int1),'%.1f')) 'g.l^{-1}'});
        end
         h_legend=legend(CHANNEL_STRING');
         title(['Dual frequency profile for ' num2str(f(2)/1000000) ' and ' num2str(f(3)/1000000) 'MHz (Honite 16, probe' num2str(int2) ')'])
         axis([0 140 0 0.3])
         xlabel('Concentration (g.l^{-1})')
         ylabel('Distance from Transducer (m)')
         set(gca,'XAxisLocation','top','YAxisLocation','left','ydir','reverse','FontSize', 12) 
         set(legend,'FontSize',10)
    end
    
    
    for int2 = CHANNEL_RANGE_DEFINED;
        
        figure
        hold on
        legend_direct=num2str( (sample_values(sample_range(FREQ_RANGE(1), CHANNEL_RANGE(FREQ_RANGE(1),1), :)))' ); 

        for int1 = 2:9;
          plot(squeeze(m_dual_1_3(int1, 1, int2, 1, :)),squeeze(r(int1, 1, 1, 1, :)), char(markers(int1)));
             CHANNEL_STRING(int1-1) = strjoin( {char(num2str(sample_values(int1),'%.1f')) 'g.l^{-1}'});
        end
      h_legend=legend(CHANNEL_STRING');
         title(['Dual frequency profile for ' num2str(f(2)/1000000) ' and ' num2str(f(3)/1000000) 'MHz (Honite 16, probe' num2str(int2) ')'])
         axis([0 140 0 0.3])
         xlabel('Concentration (g.l^{-1})')
         ylabel('Distance from Transducer (m)')
         set(gca,'XAxisLocation','top','YAxisLocation','left','ydir','reverse','FontSize', 12) 
         set(legend,'FontSize',10)
    end
    
    
    for int2 = CHANNEL_RANGE_DEFINED;
        
        figure
        hold on
        legend_direct=num2str( (sample_values(sample_range(FREQ_RANGE(1), CHANNEL_RANGE(FREQ_RANGE(1),1), :)))' ); 

        for int1 = 2:9;
          plot(squeeze(m_dual_2_3(int1, 1, int2, 1, :)),squeeze(r(int1, 1, 1, 1, :)), char(markers(int1)));
             CHANNEL_STRING(int1-1) = strjoin( {char(num2str(sample_values(int1),'%.1f')) 'g.l^{-1}'});
        end
      h_legend=legend(CHANNEL_STRING');
         title(['Dual frequency profile for ' num2str(f(2)/1000000) ' and ' num2str(f(3)/1000000) 'MHz (Honite 16, probe' num2str(int2) ')'])
         axis([0 140 0 0.3])
         xlabel('Concentration (g.l^{-1})')
         ylabel('Distance from Transducer (m)')
         set(gca,'XAxisLocation','top','YAxisLocation','left','ydir','reverse','FontSize', 12) 
         set(legend,'FontSize',10)
    end
    
    
    %%statistical analysis of experimentally determined values%%
    ks_std = std(ks_alternate_stored(:,:,:,:,108:250), 0, 5);
    ks_std_norm = ks_std./ks_alternate_stored_mean;
    
    standard_deviation_1_2 = std(m_dual_1_2(2:9,1,:,1,108:250),0,5,'omitnan')./mean(m_dual_1_2(2:9,1,:,1,108:250),5);
    standard_deviation_1_3 = std(m_dual_1_3(2:9,1,:,1,108:250),0,5,'omitnan')./mean(m_dual_1_3(2:9,1,:,1,108:250),5);
    standard_deviation_2_3 = std(m_dual_2_3(2:9,1,:,1,108:250),0,5,'omitnan')./mean(m_dual_2_3(2:9,1,:,1,108:250),5);
    std_m_single = std(m1(2:9,1,:,1,108:250),0,5,'omitnan')./mean(m1(2:9,1,:,1,108:250),5);
    std_mean_pairings = [ (mean(nonzeros((mean((standard_deviation_1_2), 1))))) (mean(nonzeros((mean((standard_deviation_1_3), 1))))) (mean(nonzeros((mean((standard_deviation_2_3), 1))))) ];
    G_fit_std = sqrt((( Gfit-squeeze(G) )./Gfit ).^2);
    G_fit_std = mean(G_fit_std(:,:,:,108:250),4);
    l=squeeze(standard_deviation_1_3)./sample_values(2:9)';
    atten_error_percent = 1-exp(-2*atten_coeff_std.*r(1,1,1,1,108:250).*sample_values(2:9)');
    
    ATTEN_COEFF_RATIO(1) = atten_coeff(1)/atten_coeff(2);
    ATTEN_COEFF_RATIO(2) = atten_coeff(1)/atten_coeff(3);
    ATTEN_COEFF_RATIO(3) = atten_coeff(2)/atten_coeff(3);

    %%plotting distance averaged particle scattering coefficient as a function of concentration
    %%normalised to maximum coefficient value for each probe. Intended to convey
    %%percentage variation in particle scattering coefficient%%
    
figure
hold on
for int1 = 1:8
    for int2 = 1:3
        ks_alternate_stored_mean_norm(2:9,1,int1,int2)= ks_alternate_stored_mean(2:9,1,int1,int2)./(max(ks_alternate_stored_mean(2:9,1,int1,int2)));
        plot(sample_values(2:9)',squeeze(ks_alternate_stored_mean_norm(2:9,1,int1,int2)));
    end
end


%%Plotting 

for int1 = 1:8
   for int2=1:3
   figure
   hold on 
        for int3 = 2:9
           
  
        plot(squeeze(r(int3,1,int1,int2,1:312)),squeeze(ks_alternate_stored(int3,1,int1,int2,1:312)),char(markers(int3)) )
      

        
        end
    CHANNEL_STRING(int3-1) = strjoin( {char(num2str(sample_values(int3),'%.1f')) 'g.l^{-1}'});
    xlabel('Distance from Transducer (m)')
    ylabel('{\itk_{s}} (m.kg^{-1/2})')
    h_legend = legend(CHANNEL_STRING')
    legend boxoff
    axis([0 0.3 0 1.5])
    end
   
end

%%plotting particle scattering coefficient as a function of concentration
%%for probe 6%%
ks_overall_mean = mean(squeeze(mean(ks_alternate_stored_mean(2:9,:,:,:),3)),1);   
markers = ["+" "*" "x"; "^" "<" "v"; "s" "d" "p"];
RBG= [  0.6902    0.0627    0.1569; 1.0000    0.6314    0.2588; 0.1725    0.5098    0.149];

    for int1 = 6
        figure
        hold on
        for int2=1:3
            plot(sample_values(2:9),ks_alternate_stored_mean(2:9,1,int1,int2),[char(markers(2,int2))], 'Color', RBG(int2,:), 'MarkerFaceColor', RBG(int2,:))
        end
        xlabel ('Concentration (g.l^{-1})')
        ylabel ('k_{s} (m.kg^{-1/2})')
        legend('2 MHz', '2.25 MHz', '2.5 MHz')
        legend boxoff
        set(gca,'FontSize', 12)
        set(legend,'FontSize',12)
        axis([0 140 0 0.7])
    end 
    
