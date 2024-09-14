##ASK Modulation and Demodulation Techniques
clc
clear all
close all

#Parameters
number_of_bits = 64;
bit_rate = 100; # bit rate (bits per second)
bit_duration = 1/bit_rate;
sampling_rate = 10*bit_rate;
sampling_time = 1/sampling_rate;
fc = 5*bit_rate; ##Carrier Frequency




##Use User Defined function to genrate stream of bits
bits = generate_random_bits(number_of_bits);

##initiate Time Vector.
t = 0:sampling_time:number_of_bits*bit_duration-sampling_time;

##initiate Frequncy Vector.
f = linspace(-sampling_rate/2, sampling_rate/2, length(t));


#Generate carrier signal at the transmitter
carrier = cos(2*pi*fc*t);

#Adjust bits length to match time vector
bits = repelem(bits, bit_duration * sampling_rate);

#ASK modulation
modulated_signal = carrier .* bits ;

Modulated_signal = fft(modulated_signal)*bit_duration;#fft for modulated signal

figure(1);
plot(t,modulated_signal,'k')
grid on
title('Modulated ASK Signal')
xlabel('Time (s)');
ylabel('Amplitude');
box off
figure(2)
plot(f,Modulated_signal,'b')
grid on
title('Frequency Representation');
xlabel('Frequency(Hz)');
ylabel('Amplitude');
box off



##Phase shifts at reciever
phases = [0, 30, 60, 90]; #in degrees #We added zero degree as reference for comparsion as at zero degree synchronization between the transmitter and osiclliator is done

#Demodulate signal (each iteration corresponds to one of phase shifts
for phase_index = 1:length(phases)

    phase_shift = phases(phase_index);
    carrier_at_reciever = cos(2*pi*fc*t + deg2rad(phase_shift));
    demodulated_signal = modulated_signal .* carrier_at_reciever;

    #Low pass filter to select the siganl centerd at zero
    LPF = zeros(size(f));
    LPF = abs(f < bit_rate);
    filtered_demodulated_signal = demodulated_signal .* LPF;
    Demodulated_signal = fftshift(fft(filtered_demodulated_signal))*bit_duration;

    # Plot temporal and spectrum of the signal at the reciever.
    figure;
    subplot(211)
    plot(t,demodulated_signal,'k');
    grid on
    title(['Demodulated Signal with Phase Shift = ', num2str(phase_shift), ' degrees']);
    xlabel('Time (s)');
    ylabel('Amplitude');
    box off
    subplot(212);
    plot(f,abs(Demodulated_signal),'b');
    grid on
    title(['Spectrum of Demodulated Signal with Phase Shift = ', num2str(phase_shift), ' degrees']);
    xlabel('Frequency (Hz)');
    ylabel('Amplitude');
    box off
end
## It is observed that as the phase angle changes from 0 to 90 degrees, Distorion happens to the reterived signal till it  reaches maximum distortion at 90 degree phase.



