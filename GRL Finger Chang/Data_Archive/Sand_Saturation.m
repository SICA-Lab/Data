%% load source data
clear; clc; close all;
load('SRC.mat');
Grid_x = [0:1:30];
Grid_y = [0:1:13];
Grid_z = [0:1:30];
Gx = numel(Grid_x);
Gy = numel(Grid_y);
Gz = numel(Grid_z);
Rcv = [15,138,15];

freq = [0.1e6:0.01e6:0.8e6];
N_fr = numel(freq);
c = 1500;

pt_fd = zeros(5,N_fr);
for nn = 1:5
    for ff = 1:N_fr
        src_sum = 0;
        k = 2*pi*freq(ff)/c;
        for ii = 1:Gz
            for jj = 1:Gy
                if jj == 1
                    SRC1 = -SRC;
                else
                    SRC1 = SRC;
                end
                for kk = 1:Gx
                    src = SRC1(nn,Gx*Gy*(ii-1) + Gy*(jj-1) + kk);
                    r = norm(Rcv - [Grid_x(kk), Grid_y(jj), Grid_z(ii)],2)*1e-3;
                    src_sum = src_sum + 1/r*exp(-j*k*r)*src;
                end
            end
        end
        pt_fd(nn,ff) = src_sum;
    end
end

T = [0:0.1:120]'*1e-6;
pt_td = cell(5,1);
pt_pp = zeros(5,1);
for nn = 1:5
    pt_td{nn,1} = real(IFFT_TA(pt_fd(nn,:).',freq,T));
    pt_pp(nn,1) = max(pt_td{nn,1}) - min(pt_td{nn,1});
end

%% load source strength
load('Q3.mat');

%% bcg
dir = pwd;
cur_folder = strcat(dir,'\Sand Saturation\bcg-Direct')
cd(cur_folder);
load('Meas.mat'); bcg=Meas;


%% dry
cur_folder = strcat(dir,'\Sand Saturation\0-Direct')
cd(cur_folder);

load('Meas.mat'); Meas(:,2)=Meas(:,2)-bcg(:,2);Meas01=Meas; 

%% 20
cur_folder = strcat(dir,'\Sand Saturation\20-Direct')
cd(cur_folder);
load('Meas.mat'); Meas(:,2)=Meas(:,2)-bcg(:,2);Meas21=Meas;

%% 40
cur_folder = strcat(dir,'\Sand Saturation\40-Direct')
cd(cur_folder);
load('Meas.mat'); Meas(:,2)=Meas(:,2)-bcg(:,2);Meas41=Meas;
%% 60
cur_folder = strcat(dir,'\Sand Saturation\60-Direct')
cd(cur_folder);
load('Meas.mat'); Meas(:,2)=Meas(:,2)-bcg(:,2);Meas61=Meas;
%% 80
cur_folder = strcat(dir,'\Sand Saturation\80-Direct')
cd(cur_folder);
load('Meas.mat'); Meas(:,2)=Meas(:,2)-bcg(:,2);Meas81=Meas;
%% 100
cur_folder = strcat(dir,'\Sand Saturation\100-Direct')
cd(cur_folder);
load('Meas.mat'); Meas(:,2)=Meas(:,2)-bcg(:,2);Meas101=Meas;
%% plot
close all; peak=[];
figure;plot(Meas01(:,1)*1e6,Meas01(:,2),'linewidth',2); title('Dry sand');xlabel('time(\mus)');ylabel('Amp(mV)');axis([50 170 -2.5e-3 2.5e-3]);set(gca,'fontsize',15);
figure;plot(Meas01(:,1)*1e6,Meas21(:,2),'linewidth',2); title('20% Saturation');xlabel('time(\mus)');ylabel('Amp(mV)');axis([50 170 -2.5e-3 2.5e-3]);set(gca,'fontsize',15);
figure;plot(Meas01(:,1)*1e6,Meas41(:,2),'linewidth',2); title('40% Saturation');xlabel('time(\mus)');ylabel('Amp(mV)');axis([50 170 -2.5e-3 2.5e-3]);set(gca,'fontsize',15);
figure;plot(Meas01(:,1)*1e6,Meas61(:,2),'linewidth',2); title('60% Saturation');xlabel('time(\mus)');ylabel('Amp(mV)');axis([50 170 -2.5e-3 2.5e-3]);set(gca,'fontsize',15);
figure;plot(Meas01(:,1)*1e6,Meas81(:,2),'linewidth',2); title('80% Saturation');xlabel('time(\mus)');ylabel('Amp(mV)');axis([50 170 -2.5e-3 2.5e-3]);set(gca,'fontsize',15);
figure;plot(Meas01(:,1)*1e6,Meas101(:,2),'linewidth',2); title('100% Saturation');xlabel('time(\mus)');ylabel('Amp(mV)');axis([50 170 -2.5e-3 2.5e-3]);set(gca,'fontsize',15);

peak=[max(Meas01(:,2))-min(Meas01(:,2)),max(Meas21(:,2))-min(Meas21(:,2)),max(Meas41(:,2))-min(Meas41(:,2)),...
    max(Meas61(:,2))-min(Meas61(:,2)),max(Meas81(:,2))-min(Meas81(:,2)),max(Meas101(:,2))-min(Meas101(:,2))];

%% compare
close all;
cd(dir);
% sat = [0 6 12 18 24 30]; 
sat = [0 20 40 60 80 100]; 
[axis,line1,line2] = plotyy(sat,peak,sat,[0;Q3]./max(Q3));
axes(axis(2)); hold on;
plot(sat,[0;pt_pp]/max(pt_pp),':','linewidth',2,'color','r'); 
axes(axis(1));
set(line1,'LineWidth',2);
set(line2,'LineWidth',2,'LineStyle','--');
xlabel('Saturation Level (%)');
ylabel(axis(1),'Experiment Results (V)')
ylabel(axis(2),'Normalized Simulation Results (a.u.)')
yticks(axis(1),[1.4e-3 3.2e-3 4.8e-3]);
set(axis(1),'YLim',[1.4e-3 4.8e-3]);
legend('P-P Amp','Source Strength','P-P Amp')
 set(axis(2),'YLim',[0 1])
set(gca,'FontSize',12)
set(axis(2),'FontSize',12)