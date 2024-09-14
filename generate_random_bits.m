#Function to generate random 0's and 1's by the program itself
function bits = generate_random_bits(num_bits)
  bits = randi([0 1], num_bits, 1)';
end
#

