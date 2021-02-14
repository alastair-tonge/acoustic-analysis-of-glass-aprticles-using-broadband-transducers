close all
clear all
load RMS_DATA_V_HON_16.mat
r = 0.3/312:0.3/312:0.3;
r=repmat(r,1,1,9,8,3);
r=permute(r, [3 1 4 5 2]);

db_data = 20*log10(RMS_DATA_V_HON_16);
G=log(r.*RMS_DATA_V_HON_16);
n = 2:9;      %%%sample plot range for dB vs distance graphs
b = 1:8;     %%%channel plot range -- change this and first subplot dimension to be the same to plot specific probes%%%
c = 1:3;     %%%frequency plot range

% load phifit_combined_hon.mat
noise_floor = 85;
noise = db_data(1,:,:,:,:);
SNR = db_data - noise ; 
SOS = 1474;
f= [2000000, 2250000, 2500000];
f2 = repmat(f, 8,1,1,9,312);
f = permute(f2, [4,3,1,2,5]);
at =   0.00318; %radius of active face of the transducer in meters
rn = pi*(at.^2)./(SOS./f);
z = r./rn;
% Phi_OG = ( 1 + 1.35.*z + ((2.5.*z).^3.2) ) ./ ( 1.35.*z + ((2.5.*z).^3.2) );
% Phi_exp = ( a1234(4) + a1234(1).*z + ((a1234(2).*z).^a1234(3)) ) ./ ( a1234(1).*z + ((a1234(2).*z).^a1234(3)) );


b1 = min(b) / min(b);
b2 = max(b) - ( min(b) - 1);
c1 = min(c) / min(c);
c2 = max(c) - (min(c) - 1);
%%following code calcualtes 'subplot index value' values to adust placement of plots onto the subplot based on desired
%%channel and frequency plot range (all sample values are plotted, adust 'n' in order to choose sample value range)

sample_values = [0 2 3.8 9.9 19.2 43 70.1 88.6 127.1 ];



for int2 = b1:b2;  %%%channel/probe plot range -- change this and subplot dimensions to plot specific probes%%%
    figure
   for int1 = n;    %%% sample plot range
  
        for int3 = c1:c2;   %%%frequency plot range
%                 int4 = ((max(c))*int2 + int3 - (max(c)));
%             
%             subplot(b2,c2,int4);
            hold on
            plot(squeeze(r(1,1,1,1,:)), squeeze(db_data(int1,1,b(int2),c(int3),:)), 'displayname', num2str(sample_values(int1-1)));
            xlabel ('Distance (m)');
            ylabel ('Backscatter Strength (dB)');
            title(['CH', num2str(b(int2)), 'FB', num2str(c(int3))]);
            legend([num2str(sample_values(n)')]);
            xlim([0 0.015])

        end
    end
end
%    
% %%%check G vs r graphs for appropriate range to take dG/dR fit!%%%
% %%%%fits a stright line over them for channels 1&2 (fit not shown)%%%%
% 
%  
%%%% Gfit a stright line over them for ALL DATA (fit not shown)%%%%


sample_range = 2:9;
channel_range = 1:8;
frequency_range = 1:3;
fit_length = 50;
fit_step = 1;
overall_fit_range_olympus = 40:220; %overall range to fir G with r, the linear fit function will search within this range for the most negative gradient fit


R = zeros (max(sample_range), max(channel_range), max(frequency_range)); 
Gfit = zeros (max(sample_range), max(channel_range), max(frequency_range));
fit_distance_range_stored = zeros (max(sample_range), max(channel_range), max(frequency_range));
max_range = zeros (max(sample_range), max(channel_range), max(frequency_range));
fit_distance_length_stored = zeros (max(sample_range), max(channel_range), max(frequency_range));

noise_floor_olympus = -86;

for int1 = sample_range
    for int2 = channel_range
   
                overall_fit_range = overall_fit_range_olympus;
                noise_floor_temp = noise_floor_olympus;
        
    
        for int3 = frequency_range
            
            if isempty(find(db_data(int1,1,int2,int3,overall_fit_range)<noise_floor_temp,1)) == 1
                max_range(int1,int2,int3) = length(overall_fit_range);
            else
                max_range(int1,int2,int3) = find(db_data(int1,1,int2,int3,overall_fit_range)<noise_floor_temp,1);
            end
            
overall_fit_range_temp = min(overall_fit_range):(max_range(int1,int2,int3) + min(overall_fit_range));
            if length(overall_fit_range_temp ) < fit_length
             fit_length_temp = length(overall_fit_range_temp);
            else
                fit_length_temp = fit_length;
            end    
    
     [fit_distance_range, R(int1,int2,int3)] = linear_fit_finder_min_grad( squeeze(r(1,1,1,1,overall_fit_range_temp)),  squeeze(G(int1,1,int2,int3,overall_fit_range_temp)),  fit_length_temp, fit_step );
P = polyfit(squeeze(r(1,1,1,1,fit_distance_range+(overall_fit_range_temp(1)-1))),squeeze(G(int1,1,int2,int3,fit_distance_range+(overall_fit_range_temp(1)-1))),1);
 Gfit(int1, int2, int3, 1:312) = (P(1).*squeeze(r(1,1,1,1,1:312)))+P(2);  
fits(int1,int2,int3)=P(1);
fit_distance_range_stored(int1,int2,int3) = fit_distance_range(1) +(overall_fit_range_temp(1)-1);
fit_distance_length_stored(int1,int2,int3) =    fit_length_temp;
%     yfit = (P(1).*squeeze(r(1,1,1,1,:)))+P(2);   
%     plot(squeeze(r(1,1,1,1,:)),yfit,'r-.');
        end
    end
end
 

%%%% plots G functions vs distance %%%%
figure
b = 1:8;     %%%channel plot range -- change this and first subplot dimension to be the same to plot specific probes%%%
c = 1:3;     %%%frequency plot range

%%following code calcualtes 'subplot index value' values to adust placement of plots onto the subplot based on desired
%%channel and frequency plot range (all sample values are plotted, adust 'n' in order to choose sample value range)
b1 = min(b) / min(b);
b2 = max(b) - ( min(b) - 1);
c1 = min(c) / min(c);
c2 = max(c) - (min(c) - 1);


             %%%subplot index range for sample number (same as sample number range)
for int2 = b1 : b2;   %%%subplot index range for selected channel plot range
   for int3 = c1 : c2;
        figure
        hold on
       for int1 = n;       
%            int4 = ((max(c))*int2 + int3 - (max(c)));
%             
%                 
%             subplot(max(b),max(c),int4);
           
            plot(squeeze(r(1,1,1,1,1:312)), squeeze(G(int1,1,b(int2),c(int3),1:312)));
             axis([0 0.3 -15 -5]);
            xlabel ('Distance (m)');
            ylabel ('G Function (-)');
            
    
%     plot(squeeze(r(1,1,1,1,:)),squeeze(Gfit(int1,b(int2),c(int3),1:312)),'k--');
       end
        title(['CH', num2str(b(int2)), 'FB', num2str(c(int3))]);
            legend([num2str(sample_values(n)')]);
            
    end
end


%%Plots dG/dr/dM for selected probes, chanels and concentrations
figure
hold on

n = 2:9;
b = [1 2 3 4 5 6 7 8];
c = 1:3;
fit_dGdr=zeros(max(b),max(c),2);
       

RBG= [  0.6902    0.0627    0.1569; 1.0000    0.6314    0.2588; 0.1725    0.5098    0.149]

 markers = [ "^" "<" "v"];
    for int2 = b;              %selects probe/ channel numbers to plot
        for int3=c;
        dgdrdm_plot(int3)=plot(sample_values(n)', squeeze(fits(n,int2,int3)),[char(markers(int3))], 'Color', RBG(int3,:), 'MarkerFaceColor', RBG(int3,:));
        Q(int3,:) = polyfit((sample_values(n)'), fits(n, int2, int3),1);
         fit_dGdr(int2,int3,:) = Q(int3,:);
         dgdrdm_fit_plot(int3)= plot(sample_values(n)',polyval(Q(int3,:),sample_values(n)'), '--', 'Color', get(dgdrdm_plot(int3),'Color'))
    yy(int2,int3,:)= polyval(Q(int3,:),sample_values(:))';
    xx(int3,n)= sample_values(n)';
 xlabel ('Concentration (g.l^{-1})')
     ylabel ('{\itdG/dr} (Np.m^{-1})')
     
          
     R_1= corrcoef(sample_values(n),squeeze(fits(n,int2,int3)));
     R_2= corrcoef(squeeze(fits(n,int2,int3)),squeeze(yy(int2,int3,n)));
     R_test(int2,int3) = R_2(2);
     R_dgdrdm(int2,int3) = R_1(2).^2;
        end
    end
 
    
     
    
    xx= squeeze(fits(:,6,:));
    zz= squeeze(yy(6,:,:))';
 
% zz = squeeze(fits(:,6,:));
 legend('2 MHz', '2.25 MHz', '2.5 MHz', '2 MHz Fit', '2.25 MHz Fit', '2.5 MHz Fit');
    legend boxoff;
    set(gca, 'FontSize', 18);
    set(dgdrdm_fit_plot(1:3), 'LineWidth', 2)
    set(dgdrdm_plot(1:3), 'MarkerSize', 8)
%     axis([0 140 -30 5])
    xticks(0:20:140)
    yticks(-15:5:5)
    
 
 

atten_coeff = -0.5*fit_dGdr(:,:,1);


save 'atten_coeff_hon_16.mat';
