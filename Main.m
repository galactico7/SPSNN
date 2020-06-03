melody_length = 500;
epoch = 10;

## Network size
n1 = 20;
nh = 1000;
nf = 20;

## melody
melody = zeros(1, melody_length+4);
melody(1:melody_length) = randi([1 n1], 1, melody_length);

## weight initialization
w1 = 0.08 + 0.15 * rand(n1, nh);  # excitatory weight matrix
w2 = 0.08 + 0.15 * rand(n1, nh);  # excitatory weight matrix
w3 = 0.08 + 0.15 * rand(n1, nh);  # excitatory weight matrix
w4 = 0.08 + 0.15 * rand(n1, nh);  # excitatory weight matrix
wh = 0.2 * ones(nh, nf);          # excitatory weight matrix

OperationEpoch = zeros(10, 3);
AccuracyEpoch = zeros(1, 10);

for epo = 1:epoch
  [w1, w2, w3, w4, wh, operation] = SequenceLearning(n1, nh, nf, melody, melody_length, w1, w2, w3, w4, wh);
  [accuracy, pred_seq] = SequenceTestSpike(n1, nh, nf, melody, melody_length, w1, w2, w3, w4, wh);
  OperationEpoch(epoch, :) = operation(1, :);
  AccuracyEpoch(1, epoch) = accuracy; 
end  
