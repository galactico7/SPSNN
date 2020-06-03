function x = ChainReset(x)
  MemorySize = LIFConstants.MemorySize;
  BufferSize = LIFConstants.BufferSize;
  
  if length(size(x)) == 3
    x(:, :, 1:BufferSize) = x(:, :, 1 + MemorySize:BufferSize + MemorySize); 
    x(:, :, 1 + BufferSize:end) = 0;
  
  elseif length(size(x)) == 2
    x(1:BufferSize, :) = x(1 + MemorySize:BufferSize + MemorySize, :);
    x(1 + BufferSize:end, :) = 0;
  end
end