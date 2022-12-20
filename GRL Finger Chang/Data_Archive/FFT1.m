function [data]=FFT1(BP,freq)
%derive FFT value of background pressure
F=[];
T=BP(:,1);
N=length(T); % # of data 
N_fr=length(freq);
for i=1:N_fr
    aux=BP(:,2).*exp(-j*2*pi*freq(i)*T);
%     for s=1:N
%         aux(s)=BP(s,2)*exp(-j*2*pi*f(i)*t(s));
%     end
    F(i)=trapz(T,aux);
end
data=F;
% figure;
% plot(f,abs(F))