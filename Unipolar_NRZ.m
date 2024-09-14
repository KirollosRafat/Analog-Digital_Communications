#Function for generating Unipolar non-return to zero line coding
function [t,x] = Unipolar_NRZ(bits, bitrate)
T = length(bits)/bitrate;
n = 10;
N = n*length(bits);
dt = T/N;
t = 0:dt:T;
x = zeros(length(t));

for i = 0:length(bits)-1
  if bits(i+1) == 1
    x(i*n+1:(i+1)*n) = 1;
  else
    x(i*n+1:(i+1)*n) = 0;
  end
end
