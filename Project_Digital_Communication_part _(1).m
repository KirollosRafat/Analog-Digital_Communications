##Two Selections of Line Codes: Unipolar NRZ and Manchester
clc
clear all
close all

##Main Program.
num_bits = 64;#Number of bits required
bits_generated = generate_random_bits(num_bits);#generate random bits by the program code itself.
bit_rate = 100; #A bitrate of 10 bit/sec


## Applying Unipolar NRZ ##
[t1,waveform_1] = Unipolar_NRZ(bits_generated,bit_rate);
figure (1)
plot(t1,waveform_1,'r','linewidth',0.8)
grid on
title('Unipolar Non-return to Zero line code in time domain')
xlim([0 .64])
ylim([0 1])
xlabel('Time(sec)','fontsize',14)
ylabel('Amplitude','fontsize',14)
box off

##Applying Manchester line code
[t2,waveform_2] = Manchester(bits_generated,bit_rate);
figure(2)
plot(t2, waveform_2,'b','linewidth',0.8);
grid on
xlim([0 0.64])
ylim([-2 2])
xlabel('Time(sec)','fontsize',14)
ylabel('Amplitude','fontsize',14)
title('Manchester line code in time domain');
box off


##Frequency Domain for UniPolar NRZ
fft_waveform_1 = fft(waveform_1)/bit_rate;
mag_fft_1 = abs(fft_waveform_1(1:floor(length(fft_waveform_1)/2)+1));
f1 = linspace(0, 2*bit_rate, length(mag_fft_1));
figure(3)
plot(f1, mag_fft_1,'r');
grid on
xlim([0 2*bit_rate])
xlabel('Frequency(Hz)');
ylabel('Magnitude(dB)');
title('Spectrum of Unipoloar NRZ line code');
box off

##Frequency Domain for Manchester line code
fft_waveform_2 = fft(waveform_2)/bit_rate;
mag_fft_2 = abs(fft_waveform_2(1:floor(length(fft_waveform_2)/2)+1));
f2= linspace(0,2*bit_rate,length(mag_fft_2));
figure(4)
plot(f2,mag_fft_2,'b');
grid on
xlim([0 2*bit_rate])
xlabel('Frequency(Hz)');
ylabel('Magnitude(dB)');
title('Spectrum of Manchester line code');
box off





