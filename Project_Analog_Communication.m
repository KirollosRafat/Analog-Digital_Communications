##Analog_Communication and Frequency Division Multiplexing
clc
clear all
close all

fs = 100;#sampling frequency
df = 0.01;#sampling frequency step
fm = 1;#message m(t) frequency
T = 100;#simulation time
ts = 1/fs;#sampling time
N = ceil(T/ts);#Number of steps

t = -(N*ts/2):ts:((N-1)*ts/2);
x = zeros(size(t));

## Define the piece-wise function
x(t > -2 & t <= -1)= t(t > -2 & t <= -1) + 2;
x(t >= -1 & t <= 1)= 1;
x(t >= 1 & t < 2 )= 2 - t(t >= 1 & t < 2 );

##plotting x(t)
figure (1)
plot(t,x)
xlim([-2 2])
ylim([0 1])
grid on
xlabel("Time(sec)","fontsize",13);
ylabel("x(t)","fontsize",13);
title("Message Siganl in time domain");
box off

if(rem(N,2)==0)
  f = - (0.5*fs) : df : (0.5*fs-df) ;
else
  f = - (0.5*fs-0.5*df) : df : (0.5*fs-0.5*df) ;
end

##X(f) using analytical method
X_analytical = 3*(sinc(f).*sinc(3*f));

##X(f) using fft built-in function
X = fftshift(fft(x))*ts;# dt non perodic

##plotting the results from the analytical and built-in functions on the same graph
figure(2)
subplot(211)
plot(f,abs(X_analytical),'r','linewidth',1,f,abs(X),'b')
grid on
xlabel("Frequency(Hz)","fontsize",13)
ylabel("|X(f)|","fontsize",13)
legend("Analytical with absolute","Using fft","fontsize",15)
box off
subplot(212)
plot(f,X_analytical,'r','linewidth',1.2,f,abs(X),'b')
grid on
xlabel("Frequency(Hz)","fontsize",13)
ylabel("|X(f)|","fontsize",13)
legend("Analytical","Using fft","fontsize",15)
box off


##Calculate the BandWidth of message x
Energy_x = sum(abs(X).^2)*df;#sum(abs(x).^2)*ts
Index = find(f==0);
Energy_acc = 0;
for index_x = Index:length(f)
  Energy_acc = df*abs(X(index_x)).^2 + Energy_acc;
  if(Energy_acc>=0.95*0.5*Energy_x)
    BW1 = f(index_x);
    break
  end
end

##Applying LPF to x(t)
LPF1 = zeros(size(f));
LPF1(abs(f) < 1)=1;

LPF2 = zeros(size(f));
LPF2(abs(f) < 0.3)=1;

x_after_LPF1 = (ifft(ifftshift(X.*LPF1))/ts); ## 1/ts as it is non-perodic signal #
x_after_LPF2 = (ifft(ifftshift(X.*LPF2))/ts);

figure (3)
subplot(311)
plot(t,x,'r')
xlim([-2 2])
ylim([0 1])
grid on
xlabel("Time(sec)","fontsize",13)
ylabel("x(t)","fontsize",13)
title("Message siganl before LPF")
box off

subplot(312)
plot(t,x_after_LPF1,'b')
xlim([-2 2])
ylim([0 1.2])
grid on
xlabel("Time(sec)","fontsize",13)
ylabel("x(t)","fontsize",13)
title("Message siganl after LPF with bandwidth = 1 Hz")
box off

subplot(313)
plot(t,x_after_LPF2,'c')
xlim([-2 2])
ylim([0 1.2])
grid on
xlabel("Time(sec)","fontsize",13)
ylabel("x(t)","fontsize",13)
title("Message siganl after LPF with bandwidth = 0.3 Hz")
box off

##New message siganl m(t)
m = cos(2*pi*fm*t);


##plotting the message siganl in time domain
figure(4)
plot(t,m,'m','linewidth',1.2)
xlim([0 6])
ylim([-1 1])
grid on
title("Sinusoidal message signal in time domain")
xlabel("Time(sec)","fontsize",13)
ylabel("m(t)","fontsize",13)
box off

m(t > 6) = 0;
m(t < 0) = 0;

##M(f) using Analytical analysis
y1 = fftshift(fft(cos(2*pi*fm*t)))/N;
rec_pulse = zeros(size(t));
rec_pulse(t >= 0 & t <= 6)=1;
y2 =fftshift(fft(rec_pulse))*ts;
M_analytical = conv(y1,y2,'same');



##M(f) using fft built-in function
M = fftshift(fft(m))*ts;

##Both plottings on the same graph
figure (5)
plot(f,abs(M_analytical),'r','linewidth',1.2,f,abs(M),'b')
grid on
xlabel("Frequency(Hz)","fontsize",13)
ylabel("|M(f)|","fontsize",13)
legend("Analytical","fft","fontsize",15)
box off


##Calculate the BandWidth of m(t)
Energy_m = sum(abs(M).^2)*df;
Index_m = find(f==0);
Energy_Acc = 0;
for indexx = Index_m : length(f)
  Energy_Acc = df*abs(M(indexx)).^2 + Energy_Acc;
  if(Energy_Acc>=0.95*0.5*Energy_m)%% 0.5 because we are calculting accumulated power in the USB only while the Total power = Power of USB + Power of LSB
    BW2 = f(indexx);
    break
 end
end


#### FDM scheme #####
BW2 = floor(BW2);
fc1 = 20; #carrier frequency used in DSB modulation
Guard_BW = 2.5; #required Guard BandWidth

## Modulating the first carrier with x_after_LPF1
c1 = cos(2*pi*fc1*t);
s1 = x_after_LPF1 .*c1;#DSB modulation
S1= fftshift(fft(s1))*ts;

# We chose to modulate the second carrier with USB
#To find the frequency of the carrier used in SSB modulation:
fc2 = ceil(fc1 + BW1 + Guard_BW + BW2);

#Modulating the second carrier with SSB scheme
c2 = cos(2*pi*fc2*t);
s2 = m .* c2;
#Band pass filter to pick the USB only
BPF = zeros(size(f));
BPF(abs(f) > fc2 & abs(f) < (fc2+3*BW2) )=1 ;


S2_DSB = fftshift(fft(s2))*ts;
S2 = S2_DSB .* BPF;
s2 = ifft(ifftshift(S2))/ts;

##step 11
s = s1 + s2;
S=fftshift(fft(s))*ts;

##Plotting the sum of the two signals in time domain
figure(6)
plot(t,s,'g','linewidth',1.3)
grid on
xlim([-2 6.4])
ylim([-1.5 1.6])
title('The sum representation of the two message signals in time domain')
xlabel('Time(sec)','fontsize',13)
ylabel('s(t)','fontsize',13)
box off

##Plotting the sum of the two signals in frequency domain
figure(7)
plot(f,abs(S),'k','linewidth',0.9)
grid on
title('Frequency representation')
xlabel('Frequency(Hz)','fontsize',13)
ylabel('|s(f)|','fontsize',13)
box off


##Demodulation of s1(t)
g1 = s1 .* cos(2*pi*fc1*t);
G1 = fftshift(fft(g1))*ts;
filter_1 = zeros(size(f));
filter_1(abs(f) < 1)=1;
X_Demodulated = G1 .* filter_1;
x_demodulated = ifft(ifftshift(X_Demodulated))/ts;



##Demodulation of s2(t)
g2 = s2 .* cos(2*pi*fc2*t);
G2 = fftshift(fft(g2))*ts;
filter_2 = zeros(size(f));
filter_2(abs(f) < (3*BW2))=1;
M_demodulated = G2.*filter_2;
m_demodulated = ifft(ifftshift(M_demodulated))/ts;


figure(8)
plot(t,(x_demodulated)/max(x_demodulated),'r','linewidth',1.3,t,x_after_LPF1/max(x_after_LPF1),'b')
grid on
xlim([-2 2])
ylim([0 1])
title('Normalized plotting of x(t) in time domain')
xlabel('Time(sec)','fontsize',13)
ylabel('x(t)','fontsize',13)
legend("Message after demodulation","Original Message",'fontsize',14)
box off

figure(9)
plot(t,(m_demodulated)/max(m_demodulated),'r','linewidth',1.5,t,m/max(m),'b')
grid on
xlim([.5 5.5])
ylim([-1 1])
title('Normalized plotting of m(t) in time domain')
xlabel('Time(sec)','fontsize',13)
ylabel('m(t)','fontsize',13)
legend("Message after demodulation","Original Message",'fontsize',14)
box off







