% this code works for loading finger imaging data
% use hann TD smooth
%% load data
clear;clc;close all;
%%
dir = pwd;
begin_x = 65-12;
begin_y = 38-12;
end_x = 145+12;
end_y = 118+12;
gap_x = 4;
gap_y = 4;
pt_x = abs(end_x - begin_x)/gap_x + 1;
pt_y = abs(end_y - begin_y)/gap_y + 1;

Rcv_x = [begin_x:gap_x:end_x];
Rcv_y = [begin_y:gap_y:end_y];
Rcv_z = [0];
Rx = numel(Rcv_x);
Ry = numel(Rcv_y);
Rz = numel(Rcv_z);

Grid_x = [65-12:2:145+12];
Grid_y = [38-12:2:118+12];
Grid_z = [100:2:140]; %define z slice here;

Gx = numel(Grid_x);
Gy = numel(Grid_y);
Gz = numel(Grid_z);

c=1500; fmin=0.1e6; fmax=0.8e6; N_fr=71; freq=linspace(fmin,fmax,N_fr);
%% Reconstruct the sensing matrix
tic;
BGF=[];
for ry=1:Ry
    for rx=1:Rx
        d=[];
        Rcv_pt=[Rcv_x(rx),Rcv_y(ry),Rcv_z];
        for gz=1:Gz
            for gy=1:Gy
                for gx=1:Gx
                    Grid_pt=[Grid_x(gx),Grid_y(gy),Grid_z(gz)];
                    r=norm([Grid_pt-Rcv_pt],2);
                    d=[d,r];
                end
            end
        end
        d=d/1000; % unit mm to m
        lambda=c./freq;
        k=2*pi./lambda;
        BGF1=[];
        for f=1:N_fr
            P=1./d.*exp(j*k(f).*d);
            BGF1=[BGF1;P];
        end
        BGF=[BGF;BGF1];
    end
end
toc;

%% load marker data
% load background data
cur_folder = strcat(dir,'\Finger Imaging\Bcg');marker =1;
cd(cur_folder);
meas_bcg = cell(pt_y,pt_x);
pow_bcg = zeros(pt_y,pt_x); 
for ii = 1:pt_x
    for jj = 1:pt_y
        s1 = num2str(begin_x+(ii-1)*gap_x);
        s2 = '-';
        s3 = num2str(begin_y+(jj-1)*gap_y);
        load(strcat(s1,s2,s3));
        meas_bcg{jj,ii} = Data;
        if marker == 1
            load(strcat(s1,s2,s3,'_P'));
            pow_bcg(jj,ii) = P*1e-3;
        else
            pow_bcg(jj,ii) = 2;
        end
    end
end
cd(dir);
meas_bcg = BcgProcess(meas_bcg);
cur_folder = strcat(dir,'\Finger Imaging\Marker'); % load marker data
cd(cur_folder);
meas_cell_1 = cell(pt_y,pt_x);
pow_meas = zeros(pt_y,pt_x);
for ii = 1:pt_x
    for jj = 1:pt_y
        s1 = num2str(begin_x+(ii-1)*gap_x);
        s2 = '-';
        s3 = num2str(begin_y+(jj-1)*gap_y);
        load(strcat(s1,s2,s3));
        meas_cell_1{jj,ii} = Data;
        load(strcat(s1,s2,s3,'_P'));
        pow_meas(jj,ii) = P*1e-3;
    end
end
for ii = 1:pt_x
    for jj = 1:pt_y
        if pow_bcg(jj,ii) == 0
            if ii == 1 && pow_bcg(jj,ii+1) ~= 0
                pow_bcg(jj,ii) = pow_bcg(jj,ii+1);
            else
                pow_bcg(jj,ii) = pow_bcg(jj,ii-1);
            end
        end
        
        if pow_meas(jj,ii) == 0
            if ii == 1 && pow_meas(jj,ii+1) ~= 0
                pow_meas(jj,ii) = pow_meas(jj,ii+1);
            else
                pow_meas(jj,ii) = pow_meas(jj,ii-1);
            end
        end
        
    end
end
if nnz(pow_meas) ~= pt_x*pt_y || nnz(pow_bcg)~= pt_x*pt_y
    error('Power = 0');
end
T_diff = zeros(pt_y,pt_x);
for ii = 1:pt_x
    for jj = 1:pt_y
        Data1 = meas_cell_1{jj,ii};
        Data2 = meas_bcg{jj,ii};
        d1 = Data1(400:550,2);
        d2 = Data2(400:550,2);
        ind1 = find(d1>max(d1)*0.8);
        ind2 = find(d2>max(d2)*0.8);
        t1 = mean(Data1(ind1,1));
        t2 = mean(Data2(ind2,1));
        T_diff (jj,ii) = t2 - t1;
    end
end
G_diff_1 = round(T_diff/0.02e-6);
seg0 = max(max(abs(G_diff_1)));
for ii= 1:pt_x
    for jj  = 1:pt_y
        data = meas_bcg{jj,ii}(1+seg0:end-seg0,:);
        meas_bcg{jj,ii} = data;
        data = meas_cell_1{jj,ii}(1+seg0-G_diff_1(jj,ii):end-seg0-G_diff_1(jj,ii),:);
        meas_cell_1{jj,ii} = data;
        meas_cell_1{jj,ii}(:,1) = meas_bcg{jj,ii}(:,1);
    end
end
Han = hann(1200 - 2*seg0);
cd(dir);
meas_subtr_1 = meas_cell_1;% background removal of meas1
meas_subtr_2 = meas_cell_1;
for ii = 1:pt_x
    for jj = 1:pt_y
        meas_subtr_1{jj,ii}(:,2) = meas_cell_1{jj,ii}(:,2)./pow_meas(jj,ii) - meas_bcg{jj,ii}(:,2)./pow_bcg(jj,ii);
        meas_subtr_2{jj,ii}(:,2) = meas_subtr_1{jj,ii}(:,2).*Han;
    end
end
meas_0 = meas_subtr_2;


%% load finger 1 data
% load background data
cur_folder = strcat(dir,'\Finger Imaging\Marker');marker =1;
cd(cur_folder);
meas_bcg = cell(pt_y,pt_x);
pow_bcg = zeros(pt_y,pt_x); 
for ii = 1:pt_x
    for jj = 1:pt_y
        s1 = num2str(begin_x+(ii-1)*gap_x);
        s2 = '-';
        s3 = num2str(begin_y+(jj-1)*gap_y);
        load(strcat(s1,s2,s3));
        meas_bcg{jj,ii} = Data;
        if marker == 1
            load(strcat(s1,s2,s3,'_P'));
            pow_bcg(jj,ii) = P*1e-3;
        else
            pow_bcg(jj,ii) = 2;
        end
    end
end
cd(dir);

% load meas data
cur_folder = strcat(dir,'\Finger Imaging\Finger1')
cd(cur_folder);
meas_cell_1 = cell(pt_y,pt_x);
pow_meas = zeros(pt_y,pt_x);
for ii = 1:pt_x
    for jj = 1:pt_y
        s1 = num2str(begin_x+(ii-1)*gap_x);
        s2 = '-';
        s3 = num2str(begin_y+(jj-1)*gap_y);
        load(strcat(s1,s2,s3));
        meas_cell_1{jj,ii} = Data;
        load(strcat(s1,s2,s3,'_P'));
        pow_meas(jj,ii) = P*1e-3;
    end
end

for ii = 1:pt_x
    for jj = 1:pt_y
        if pow_bcg(jj,ii) == 0
            if ii == 1 && pow_bcg(jj,ii+1) ~= 0
                pow_bcg(jj,ii) = pow_bcg(jj,ii+1);
            else
                pow_bcg(jj,ii) = pow_bcg(jj,ii-1);
            end
        end
        
        if pow_meas(jj,ii) == 0
            if ii == 1 && pow_meas(jj,ii+1) ~= 0
                pow_meas(jj,ii) = pow_meas(jj,ii+1);
            else
                pow_meas(jj,ii) = pow_meas(jj,ii-1);
            end
        end
        
    end
end

if nnz(pow_meas) ~= pt_x*pt_y || nnz(pow_bcg)~= pt_x*pt_y
    error('Power = 0');
end

T_diff = zeros(pt_y,pt_x);
for ii = 1:pt_x
    for jj = 1:pt_y
        Data1 = meas_cell_1{jj,ii};
        Data2 = meas_bcg{jj,ii};
        d1 = Data1(400:550,2);
        d2 = Data2(400:550,2);
        ind1 = find(d1>max(d1)*0.8);
        ind2 = find(d2>max(d2)*0.8);
        t1 = mean(Data1(ind1,1));
        t2 = mean(Data2(ind2,1));
        T_diff (jj,ii) = t2 - t1;
    end
end
G_diff_1 = round(T_diff/0.02e-6);
seg1 = max(max(abs(G_diff_1)));
for ii= 1:pt_x
    for jj  = 1:pt_y
        data = meas_bcg{jj,ii}(1+seg1:end-seg1,:);
        meas_bcg{jj,ii} = data;
        data = meas_cell_1{jj,ii}(1+seg1-G_diff_1(jj,ii):end-seg1-G_diff_1(jj,ii),:);
        meas_cell_1{jj,ii} = data;
        meas_cell_1{jj,ii}(:,1) = meas_bcg{jj,ii}(:,1);
    end
end
Han = hann(1200 - 2*seg1);
cd(dir);
% background removal of meas1
meas_subtr_1 = meas_cell_1;
meas_subtr_2 = meas_cell_1;
for ii = 1:pt_x
    for jj = 1:pt_y
        meas_subtr_1{jj,ii}(:,2) = meas_cell_1{jj,ii}(:,2)./pow_meas(jj,ii) - meas_bcg{jj,ii}(:,2)./pow_bcg(jj,ii);
        meas_subtr_2{jj,ii}(:,2) = meas_subtr_1{jj,ii}(:,2).*Han;
    end
end
meas_1 = meas_subtr_2;

%% load finger 2 data
% load background data
cur_folder = strcat(dir,'\Finger Imaging\Marker');marker =1;
cd(cur_folder);
meas_bcg = cell(pt_y,pt_x);
pow_bcg = zeros(pt_y,pt_x); 
for ii = 1:pt_x
    for jj = 1:pt_y
        s1 = num2str(begin_x+(ii-1)*gap_x);
        s2 = '-';
        s3 = num2str(begin_y+(jj-1)*gap_y);
        load(strcat(s1,s2,s3));
        meas_bcg{jj,ii} = Data;
        if marker == 1
            load(strcat(s1,s2,s3,'_P'));
            pow_bcg(jj,ii) = P*1e-3;
        else
            pow_bcg(jj,ii) = 2;
        end
    end
end
cd(dir);

% load meas data
cur_folder = strcat(dir,'\Finger Imaging\Finger2')
cd(cur_folder);
meas_cell_1 = cell(pt_y,pt_x);
pow_meas = zeros(pt_y,pt_x);
for ii = 1:pt_x
    for jj = 1:pt_y
        s1 = num2str(begin_x+(ii-1)*gap_x);
        s2 = '-';
        s3 = num2str(begin_y+(jj-1)*gap_y);
        load(strcat(s1,s2,s3));
        meas_cell_1{jj,ii} = Data;
        load(strcat(s1,s2,s3,'_P'));
        pow_meas(jj,ii) = P*1e-3;
    end
end

for ii = 1:pt_x
    for jj = 1:pt_y
        if pow_bcg(jj,ii) == 0
            if ii == 1 && pow_bcg(jj,ii+1) ~= 0
                pow_bcg(jj,ii) = pow_bcg(jj,ii+1);
            else
                pow_bcg(jj,ii) = pow_bcg(jj,ii-1);
            end
        end
        
        if pow_meas(jj,ii) == 0
            if ii == 1 && pow_meas(jj,ii+1) ~= 0
                pow_meas(jj,ii) = pow_meas(jj,ii+1);
            else
                pow_meas(jj,ii) = pow_meas(jj,ii-1);
            end
        end
        
    end
end

if nnz(pow_meas) ~= pt_x*pt_y || nnz(pow_bcg)~= pt_x*pt_y
    error('Power = 0');
end

T_diff = zeros(pt_y,pt_x);
for ii = 1:pt_x
    for jj = 1:pt_y
        Data1 = meas_cell_1{jj,ii};
        Data2 = meas_bcg{jj,ii};
        d1 = Data1(400:550,2);
        d2 = Data2(400:550,2);
        ind1 = find(d1>max(d1)*0.8);
        ind2 = find(d2>max(d2)*0.8);
        t1 = mean(Data1(ind1,1));
        t2 = mean(Data2(ind2,1));
        T_diff (jj,ii) = t2 - t1;
    end
end
G_diff_1 = round(T_diff/0.02e-6);
seg4 = max(max(abs(G_diff_1)));
for ii= 1:pt_x
    for jj  = 1:pt_y
        data = meas_bcg{jj,ii}(1+seg4:end-seg4,:);
        meas_bcg{jj,ii} = data;
        data = meas_cell_1{jj,ii}(1+seg4-G_diff_1(jj,ii):end-seg4-G_diff_1(jj,ii),:);
        meas_cell_1{jj,ii} = data;
        meas_cell_1{jj,ii}(:,1) = meas_bcg{jj,ii}(:,1);
    end
end
Han = hann(1200 - 2*seg4);
cd(dir);
% background removal of meas1
meas_subtr_1 = meas_cell_1;
meas_subtr_2 = meas_cell_1;
for ii = 1:pt_x
    for jj = 1:pt_y
        meas_subtr_1{jj,ii}(:,2) = meas_cell_1{jj,ii}(:,2)./pow_meas(jj,ii) - meas_bcg{jj,ii}(:,2)./pow_bcg(jj,ii);
        meas_subtr_2{jj,ii}(:,2) = meas_subtr_1{jj,ii}(:,2).*Han;
    end
end
meas_2 = meas_subtr_2;

%% Conduct FFT of measured data
% FFT of marker
cd(dir);
N_rcv = pt_x*pt_y;

% define excitation signal
Signal=[[0:0.1:120]',zeros(1201,1)];
Signal(1:10,2)=1;
Signal(:,1) = Signal(:,1)*1e-6;
Signal=FFT1(Signal,freq);

meas = [];
for jj = 1:pt_y
    for ii = 1:pt_x
        meas = [meas;[meas_0{jj,ii}(:,1),meas_0{jj,ii}(:,2)]];
    end
end
YY=[];
for ii=1:N_rcv
    Y=FFT1(meas((ii-1)*(1200-2*seg0)+1:ii*(1200-2*seg0),:),freq);
    Y=Y./Signal;
    YY=[YY,Y];
end
YY0=YY';

% FFT of finger 1
meas = [];
for jj = 1:pt_y
    for ii = 1:pt_x
        meas = [meas;[meas_1{jj,ii}(:,1),meas_1{jj,ii}(:,2)]];
    end
end
YY=[];
for ii=1:N_rcv
    Y=FFT1(meas((ii-1)*(1200-2*seg1)+1:ii*(1200-2*seg1),:),freq);
    Y=Y./Signal;
    YY=[YY,Y];
end
YY1=YY';

% FFT of finger 2
meas = [];
for jj = 1:pt_y
    for ii = 1:pt_x
        meas = [meas;[meas_2{jj,ii}(:,1),meas_2{jj,ii}(:,2)]];
    end
end
YY=[];
for ii=1:N_rcv
    Y=FFT1(meas((ii-1)*(1200-2*seg4)+1:ii*(1200-2*seg4),:),freq);
    Y=Y./Signal;
    YY=[YY,Y];
end
YY2=YY';


