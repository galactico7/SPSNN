function [accuracy, pred_seq] = SequenceTestSpike(n1, nh, nf, melody, melody_length, w1, w2, w3, w4, wh)
  ## setup parameters and state variables
  T = (melody_length +4) * 1 * 0.1;    # total time to simulate (sec), melody_length * 1_repeat * 100ms
  dt = 0.001;    # simulation time step (sec)
  Num = floor(T/dt) + 1; # total time step

  ## LIF properties
  Rm = LIFConstants.Rm;  
  tm = LIFConstants.tm;  
  td = LIFConstants.td;  
  ts = LIFConstants.ts; 
  usth = LIFConstants.usth; 
  LI = LIFConstants.LI;  
  sf = LIFConstants.sf;
  I_0 = LIFConstants.I_0;
  chaindelay = LIFConstants.chaindelay; 
  layerdelay = LIFConstants.layerdelay;   
  period = LIFConstants.period;   
  MemorySize = LIFConstants.MemorySize;
  BufferSize = LIFConstants.BufferSize;

  tusf_sp = zeros(Num, nf); # last layer firing time

  tI_syn1h = zeros(n1, nh, MemorySize + BufferSize);
  tI_syn2h = zeros(n1, nh, MemorySize + BufferSize);
  tI_syn3h = zeros(n1, nh, MemorySize + BufferSize);
  tI_syn4h = zeros(n1, nh, MemorySize + BufferSize);
  tI_synhf = zeros(nh, nf, MemorySize + BufferSize);

  ## membrane potential
  tus1 = zeros(MemorySize + BufferSize, n1);  # 1st layer somatic potential
  tush = zeros(MemorySize + BufferSize, nh);  # hidden layer somatic potential  
  tusf = zeros(MemorySize + BufferSize, nf);  # last layer somatic potential
  
  spikefind = 0;
  count = 1;
  pred_seq = zeros(1, melody_length+4); 
 
  ## iterate over each time step
  for ii = 1:Num-1        
    Iinj1 = zeros(1, n1);
    
    sp1 = zeros(n1, 1);
    sph = zeros(nh, 1);
    
    if ceil(ii/period) <= melody_length
      Iinj1(1, melody(ceil(ii/period))) = I_0;
    end
    
    i = mod(ii, MemorySize);
    if i == 0 i = MemorySize; end
    
    ## Buffer reset
    if i == 1
      tI_syn1h = ChainReset(tI_syn1h); tI_syn2h = ChainReset(tI_syn2h); tI_syn3h = ChainReset(tI_syn3h); tI_syn4h = ChainReset(tI_syn4h); tI_synhf = ChainReset(tI_synhf);
      tus1 = ChainReset(tus1); tush = ChainReset(tush); tusf = ChainReset(tusf);
    end
    
    ## EPSC integration
    tus1(i+1, :) = tus1(i, :) + (1 ./ tm) .* (-tus1(i, :) + Rm .* Iinj1(1, :)) .* dt;
     
    if(max(tus1(i, :))>=usth)
      l1 = find(tus1(i, :) >= usth);
      sp1(l1) = I_0;
      tus1(i+1, l1) = 0;
      
      tI_syn1h(:, : , i + layerdelay) = tI_syn1h(:, : , i + layerdelay) + sf * sp1.* w1;
      tI_syn2h(:, : , i + layerdelay + chaindelay) = tI_syn2h(:, : , i + layerdelay + chaindelay) + sf * sp1.* w2;
      tI_syn3h(:, : , i + layerdelay + 2 * chaindelay) = tI_syn3h(:, : , i + layerdelay + 2 * chaindelay) + sf *sp1.* w3;
      tI_syn4h(:, : , i + layerdelay + 3 * chaindelay) = tI_syn4h(:, : , i + layerdelay + 3 * chaindelay) + sf *sp1.* w4; 
    end                    
        
    tI_syn1h(:, :, i + layerdelay + 1) = tI_syn1h(:, :, i + layerdelay) - (1 ./ ts) .* tI_syn1h(:, :, i + layerdelay) .* dt;
    tI_syn2h(:, :, i + layerdelay + chaindelay + 1) = tI_syn2h(:, :, i + layerdelay + chaindelay) - (1 ./ ts) .* tI_syn2h(:, :, i + layerdelay + chaindelay) .* dt;
    tI_syn3h(:, :, i + layerdelay + 2 * chaindelay + 1) = tI_syn3h(:, :, i + layerdelay + 2 * chaindelay) - (1 ./ ts) .* tI_syn3h(:, :, i + layerdelay + 2 * chaindelay) .* dt;
    tI_syn4h(:, :, i + layerdelay + 3 * chaindelay + 1) = tI_syn4h(:, :, i + layerdelay + 3 * chaindelay) - (1 ./ ts) .* tI_syn4h(:, :, i + layerdelay + 3 * chaindelay) .* dt;
    
    ## EPSC integration
    tush(i+1, :) = tush(i, :) + (1 ./ tm) .* (-tush(i, :) + Rm .* (sum(tI_syn1h(:, :, i)) + sum(tI_syn2h(:, :, i)) + sum(tI_syn3h(:, :, i)) + sum(tI_syn4h(:, :, i)))) .* dt;

    if(max(tush(i, :)) >= usth)
      lh = find(tush(i, :)>=usth);
      sph(lh) = I_0;
      tush(i+1, lh) = 0;
      
      ## lateral inhibition   
      tush(i+1, :) -= 0.9*usth;        
      tush(i+1, find(tush(i+1, :)<0)) = 0;
      
      if i + layerdelay < Num    
        tI_synhf(:, : , i + layerdelay) = tI_synhf(:, : , i + layerdelay) + sf * sph.* wh;
      end
    end

    if i + layerdelay + 1 < Num   ## update the synaptic current
      tI_synhf(:, :, i + layerdelay + 1) = tI_synhf(:, :, i + layerdelay) - (1 ./ ts) .* tI_synhf(:, :, i + layerdelay) .* dt; 
    end
    
    ## EPSC integration
    tusf(i+1, :) = tusf(i, :) + (1 ./ tm) .* (-tusf(i, :) + Rm .* (sum(tI_synhf(:, :, i)))) .* dt;
          
    if(max(tusf(i, :)) >= usth)
      lf = find(tusf(i, :)>=usth);
      tusf(i+1, lf) = 0;
      tusf_sp(ii, lf) = ii;
     
      ## lateral inhibition
      tusf(i+1, :) -= LI * usth;           
      tusf(i+1, find(tusf(i+1, :)<0)) = 0;
    end
 
    if ii >= 4 * period && mod(ii, period) == 25
      spikefind = 1;
    end
    
    if spikefind == 1 && max(tusf_sp(ii, :)) != 0
      pred_seq(count) = min(find(tusf_sp(ii, :) != 0));
      count += 1;
      spikefind = 0;
    end  
  end 

  count2 = 0;
  for j6=5:melody_length
    if melody(j6) == pred_seq(j6-4)
      count2+=1;
    end
  end
  accuracy = count2 / (melody_length - 4);
end
