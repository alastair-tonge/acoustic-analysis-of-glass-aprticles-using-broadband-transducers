%Accepts the particle diameter                            a    [m]
%the first isonification frequency                        f1    [Hz]
%the second isonification frequency                       f2    [Hz]
%first frequency designation (high,middle,low->1,2,3)     fb1   [integer]
%second frequency designation (high,middle,low->1,2,3)    fb2   [integer]
%sound velocity in water                                  c    [m/s]
%the particle density                                     rho_solid  [kg/m^3]
%the first particle scattering constant (theoretical)     ks1   [-]
%the second particle scattering constant (theoretical)    ks2   [-]
%the transducer constant                                  kt_conc_average   [-]
%Fluid Temperature                                        T    [Celcius]
% First Bacscatter Voltage profiles from 'Analysis_HONxx_DC.m'   A     [RMS_DATA_V.mat]
% Second Bacscatter Voltage profiles from 'Analysis_HONxx_DC.m'  B     [RMS_DATA_V.mat]
% Attenuation due to the particle                                alpha_s
%Outputs water attenuation constant                             alpha_w   [-]
%first Particle concentration from standard inversion           m1     [g/l]
%second Particle concentration from standard inversion          m2     [g/l]
%particle concentration form dual frequency inversion           m_dual [g/l]
clear all
close all
%%kt_conc_average calculation based on first 3 concs only%%



     
     
     
     
              
RMS_DATA_V= 'RMS_DATA_V_HON_16.mat' ;  %SELECT .mat DATA FILE TO LOAD FROM
load(RMS_DATA_V)
load 'atten_coeff_hon_16.mat';
A= RMS_DATA_V_HON_16;            %select RMS data file loaded from data file       
r = 0.3/312:0.3/312:0.3;
r=repmat(r,1,1,9,8,3);
r=permute(r, [3 1 4 5 2]);
a= 0.000078737/2

f= [2000000, 2250000, 2500000];
at =   0.00318 %radius of active face of the transducer in meters
c= 1474;
rho_solid = 2469.6;
x= 2*3.14.*f.*a/c;

sigma= (1-0.5.*exp(-(((x-1.5)./0.5).^2))) .* (1+0.4.*exp(-(((x-1.5)./3).^2))) .* (1-0.5.*exp(-(((x-5.9)./0.7).^2)));

form_function = (sigma.*(x.^2))./(1.17+0.95.*(x.^2));

ks= form_function./(sqrt(a.*rho_solid));

kt_fit_range = 100:200;
sample_range = zeros(1);
%%%%%%%%%%CREATE A LINE ON THE PHI MODEL TO EXCLUDE DATA, PLOT IT WHERE YOU
%%%%%%%%%%WANT DAT AUP TO AND EXCLUDE VALUES < OR > THE VALUE OF THE LINE
%%%%%%%%%%MUCH EASIER TO DESCRIBE IN METHODOLOGY
  FREQ_RANGE = 1:3; %%select frequency range (3rd frequency excluded as fit found to be poor)
     CHANNEL_RANGE_DEFINED= [1 2 3 4 5 6 7 8]; %%select channel range (probes 4 ad 7 seemed to function fine for this experiment)
     sample_range_DEFINED = 2:9;  %choose sample value range
      sample_range = zeros(length(FREQ_RANGE), length(CHANNEL_RANGE_DEFINED), length(sample_range_DEFINED));
      
     for int1 = FREQ_RANGE
         CHANNEL_RANGE(int1,:) = CHANNEL_RANGE_DEFINED;
         for int2 = CHANNEL_RANGE_DEFINED
             
          sample_range(int1,int2,:) = sample_range_DEFINED;       
                 
             
         end
     end
    

T= 17;

alpha_w= 0.05641.*((f./1000000).^2).*exp(-T/27);

% atten_coeff = mean(atten_coeff,1);
correction=1;

                                             %%%%temp in celcius!!%%%%


% sample_values = [0 2.5 5.1 12.8 25.7 50 78.6 100 133.7 ];
sample_values = [0 2 3.8 9.9 19.2 43 70.1 88.6 127.1 ];
legend_direct=num2str(sample_values');



 markers = [ "+" "^" "<" "v" "o" "s" "d" "*" "x" "+"];
CHANNEL_STRING =" ";



%%%%NFC ERROR ANALYSIS STARTS HERE%%%%

ln_Phifit = squeeze(Gfit) - squeeze(G) ;
Phifit = exp(ln_Phifit);


f2 = repmat(f, 9,1,8,312);
f = permute(f2, [1,3,2,4]);
rn = pi*(at^2)./(c./f);

r_range = 15:60; %%range over which to take fit for Phifitted calculation
z = squeeze(r)./(rn);
Phioriginal = ( 1 + 1.35.*z + ((2.5.*z).^3.2) ) ./ ( 1.35.*z + ((2.5.*z).^3.2) );
z_range = z(sample_range(FREQ_RANGE(1), CHANNEL_RANGE(FREQ_RANGE(1),1), :),CHANNEL_RANGE(FREQ_RANGE(1),:),FREQ_RANGE,r_range);  %%define range of z values used for fit
Phifit_range=(Phifit(sample_range(FREQ_RANGE(1), CHANNEL_RANGE(FREQ_RANGE(1),1), :),CHANNEL_RANGE(FREQ_RANGE(1),:),FREQ_RANGE,r_range)); %%define range of Phifit values used for fit (must be same dims as z_range)

%%%This loop excludes selected data from being fitted in case only certain
%%%combinations of conentrations/probes/frequencies want to be omitted from
%%%the fit data (the points will be plotted still, just excluded in the
%%%calcualtion of fitting aprameters)


%%comment this out if no particular combinations are to be avoided%%

     FREQ_RANGE_exclude = 1:3; %%select frequency range 
     CHANNEL_RANGE_exclude= 1:8; %%select channel range 
     sample_range_exclude = 2:3;  %choose sample value range



     FREQ_RANGE_exclude = 1:3; %%select frequency range 
     CHANNEL_RANGE_exclude= [7]; %%select channel range 
     sample_range_exclude = 1:9;  %choose sample value range


   for int1 = FREQ_RANGE_exclude
         
         for int2 = CHANNEL_RANGE_exclude
             for int3 = sample_range_exclude
                 
          A1 = find(sample_range(FREQ_RANGE(1), CHANNEL_RANGE(1),:)==int3,1);   
          A2 = find(CHANNEL_RANGE(FREQ_RANGE(1),:)==int2,1);  
          A3 = find(FREQ_RANGE==int1,1); 
          Phifit_range(A1, A2, A3, :) = 0;
          z_range(A1, A2, A3, :) = 0;
          
             end
          end
   end 
   

Phifit_range(find(Phifit_range<1)) = 0;
z_range(find(Phifit_range<1)) = 0;

%%calling fit function for NFCF and plotting results%%
[fitresult, gof, opts] = create_fit_NFCF(nonzeros(z_range), nonzeros(Phifit_range));
a1234 = coeffvalues(fitresult);

Phifitted = ( a1234(4) + a1234(1).*z + ((a1234(2).*z).^a1234(3)) ) ./ ( a1234(1).*z + ((a1234(2).*z).^a1234(3)) );
Phioriginal = ( 1 + 1.35.*z + ((2.5.*z).^3.2) ) ./ ( 1.35.*z + ((2.5.*z).^3.2) );
G_new = squeeze(G) + log(Phifitted);

Phifit_stats(sample_range(FREQ_RANGE(1), CHANNEL_RANGE(FREQ_RANGE(1),1), :)) = gof;
R = (gof.rsquare);
n = sample_values(sample_range(FREQ_RANGE(1), CHANNEL_RANGE(FREQ_RANGE(1),1), :));      %%%sample plot range for dB vs distance graphs
b = CHANNEL_RANGE(FREQ_RANGE(1),:);     %%%channel plot range -- change this and first subplot dimension to be the same to plot specific probes%%%
d = FREQ_RANGE;     %%%frequency plot range

b1 = min(b) / min(b);
b2 = max(b) - ( min(b) - 1);
d1 = min(d) / min(d);
d2 = max(d) - (min(d) - 1);


l = jet( length(CHANNEL_RANGE(FREQ_RANGE(1),:)) ); %%GETS RGB VALUES FOR NUMBER OF PROBES SELECTED FOR REFERENCE WHEN PLOTTING
           l(max((sample_range(FREQ_RANGE(1), CHANNEL_RANGE(1),:)))+1, :) = [1 1 1]; %%SETS LAST VALUE OF L TO BE WHITE FOR FACELESS MARKER PLOTTING  
  FREQ_STRING = ["2 MHz " "2.25 MHZ " "2.5 MHz "];
  
  CHANNEL_STRING = "";
  markers = ["+" "*" "x"; "^" "<" "v"; "s" "d" "p"];
  model_data_plot_range = 20:50;
  fig = figure 
  hold on
      r_range_plot=1:length(r_range);     
      
for int3 = FREQ_RANGE
%         CHANNEL_RANGE(FREQ_RANGE(1),:) = 1:7;
  for int2 = CHANNEL_RANGE(FREQ_RANGE(1),:)   
    
%        sample_range(FREQ_RANGE(1), CHANNEL_RANGE(FREQ_RANGE(1),1), :)=1:8;     
           
    for int1 = sample_range(FREQ_RANGE(1), CHANNEL_RANGE(FREQ_RANGE(1),1), :)
      

            
            plot(squeeze(z_range(int1-1,int2,int3,r_range_plot)), squeeze(Phifit_range(int1-1,int2,int3,r_range_plot)), [char(markers(3,int3))], 'MarkerEdgeColor',[l(int2,:)], 'MarkerFaceColor',[l(int1,:)],'MarkerSize',15);
            
%             axis([0 0.3 -15 -5]);
            xlabel ('z (-)');
            ylabel ('\Psi_G (-)');
             
    int4 =  ( (find(FREQ_RANGE==int3,1)-1)*( (length(sample_range(FREQ_RANGE(1), CHANNEL_RANGE(FREQ_RANGE(1),1), :))) * (length(CHANNEL_RANGE(FREQ_RANGE(1),:))) ) ) +  ( (find(CHANNEL_RANGE(FREQ_RANGE(1),:)==int2,1)-1) * (length(sample_range(FREQ_RANGE(1), CHANNEL_RANGE(FREQ_RANGE(1),1), :))) ) + find(sample_range(FREQ_RANGE(1), CHANNEL_RANGE(FREQ_RANGE(1),1), :)==int1,1);
 
       
    CHANNEL_STRING(int4) = strjoin( {char(FREQ_STRING(int3)),  'probe',  char(string(int2)), char(string(sample_values(int1))), 'g/l'} );
       
    
      
     end
      
    
  end
end

 CHANNEL_STRING(int4+1)= "Phi fit curve";
             CHANNEL_STRING(int4+2)= "Phi fit curve - RMSE";
             CHANNEL_STRING(int4+3)= "Phi fit curve + RMSE";
             CHANNEL_STRING(int4+4)= "Original NFCF model";

 title(['\Psi_G fit curve for HON 16']);
 
 %%trimming data outside of one RMSE for re-fitting%%
 logic_Phifit_range_above = Phifit_range < Phifitted(2:9,:,:,r_range) - 2*(gof(1).rmse)  ;
 logic_Phifit_range_below =  Phifit_range > Phifitted(2:9,:,:,r_range) + 2*(gof(1).rmse);

Phifit_range(logic_Phifit_range_above) = 0;
Phifit_range(logic_Phifit_range_below) = 0;
z_range(logic_Phifit_range_above) = 0;
z_range(logic_Phifit_range_below) = 0;
 [fitresult, gof, opts] = create_fit_NFCF(nonzeros(z_range), nonzeros(Phifit_range));
a1234 = coeffvalues(fitresult);

Phifitted = ( a1234(4) + a1234(1).*z + ((a1234(2).*z).^a1234(3)) ) ./ ( a1234(1).*z + ((a1234(2).*z).^a1234(3)) )
Phioriginal = ( 1 + 1.35.*z + ((2.5.*z).^3.2) ) ./ ( 1.35.*z + ((2.5.*z).^3.2) );
G_new = squeeze(G) + log(Phifitted);



figure           
plot(squeeze(z(1,1,1,20:50)), squeeze(Phioriginal(2,1,1,20:50)), 'k');
hold on

    plot(squeeze(z(1,1,1,20:50)), squeeze(Phifitted(2,1,1,20:50)));
    plot(squeeze(z(1,1,1,20:50)), squeeze(Phifitted(2,1,1,20:50) - (gof(1).rmse)))
    plot(squeeze(z(1,1,1,20:50)), squeeze(Phifitted(2,1,1,20:50) + (gof(1).rmse)))



