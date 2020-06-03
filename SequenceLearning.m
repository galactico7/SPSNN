function [w1, w2, w3, w4, wh, operation] = SequenceLearning(n1, nh, nf, melody, melody_length, w1, w2, w3, w4, wh)
  
  ## setup parameters and state variables
  T = (melody_length + 4) * 1 * 0.1;    # total time to simulate (sec), melody_length * 1_repeat * 100ms
  dt = 0.001;                           # simulation time step (sec)
  Num = floor(T/dt) + 1;                # total time step

  ## LIF properties
  Rm = LIFConstants.Rm;  
  tm = LIFConstants.tm;  
  td = LIFConstants.td;  
  ts = LIFConstants.ts; 
  usth = LIFConstants.usth; 
  udth1 = LIFConstants.udth1; 
  udth2 = LIFConstants.udth2; 
  LI = LIFConstants.LI;  
  sf = LIFConstants.sf;
  I_0 = LIFConstants.I_0;
  chaindelay = LIFConstants.chaindelay; 
  layerdelay = LIFConstants.layerdelay;   
  period = LIFConstants.period;  
  MemorySize = LIFConstants.MemorySize;
  BufferSize = LIFConstants.BufferSize;
  
  ## weight matrix
  maxw = 0.25;    # max weight
  maxw2 = 0.75;   # max weight
  dw = 0.03;      # delta weight  

  operation = zeros(1, 3);

  I_syn1h = zeros(n1, nh, MemorySize + BufferSize);
  I_syn2h = zeros(n1, nh, MemorySize + BufferSize);
  I_syn3h = zeros(n1, nh, MemorySize + BufferSize);
  I_syn4h = zeros(n1, nh, MemorySize + BufferSize);
  I_synhf = zeros(nh, nf, MemorySize + BufferSize);

  I_synd1h = zeros(n1, nh, MemorySize + BufferSize);
  I_synd2h = zeros(n1, nh, MemorySize + BufferSize);
  I_synd3h = zeros(n1, nh, MemorySize + BufferSize);
  I_synd4h = zeros(n1, nh, MemorySize + BufferSize);

  ## membrane potential
  us1 = zeros(MemorySize + BufferSize, n1);      # input layer somatic potential
  ush = zeros(MemorySize + BufferSize, nh);      # hidden layer somatic potential  
  usf = zeros(MemorySize + BufferSize, nf);      # output layer somatic potential

  ud1 = zeros(n1, nh, MemorySize + BufferSize);  # dendritic potential
  ud2 = zeros(n1, nh, MemorySize + BufferSize);  # dendritic potential
  ud3 = zeros(n1, nh, MemorySize + BufferSize);  # dendritic potential
  ud4 = zeros(n1, nh, MemorySize + BufferSize);  # dendritic potential
  udh = zeros(nh, nf, MemorySize + BufferSize);  # dendritic potential

  ## iterate over each time step
  for ii = 1:Num-1
    Iinj1 = zeros(1, n1);
    Iinj2 = zeros(1, nf);
    
    sp1 = zeros(n1, 1);
    sph = zeros(nh, 1);
    
    if ceil(ii/period) <= melody_length
      Iinj1(1, melody(ceil(ii/period))) = I_0;
    end
    
    if ii > 0.6 * period && ceil((ii + 0.4 * period)/period) <= melody_length
      Iinj2(1, melody(ceil((ii + 0.4 * period)/period))) = I_0;
    end    
  
    i = mod(ii, MemorySize);
    if i == 0 i = MemorySize; end
    
    ## Buffer reset
    if i == 1
      I_syn1h = ChainReset(I_syn1h); I_syn2h = ChainReset(I_syn2h); I_syn3h = ChainReset(I_syn3h); I_syn4h = ChainReset(I_syn4h); I_synhf = ChainReset(I_synhf);
      I_synd1h = ChainReset(I_synd1h); I_synd2h = ChainReset(I_synd2h); I_synd3h = ChainReset(I_synd3h); I_synd4h = ChainReset(I_synd4h);
      ud1 = ChainReset(ud1); ud2 = ChainReset(ud2); ud3 = ChainReset(ud3); ud4 = ChainReset(ud4); udh = ChainReset(udh);
      us1 = ChainReset(us1); ush = ChainReset(ush); usf = ChainReset(usf);   
    end
    
    us1(i+1, :) = us1(i, :) + (1 ./ tm) .* (-us1(i, :) + Rm .* Iinj1(1, :)) .* dt;    ## EPSC integration
    
    if(max(us1(i, :))>=usth)
      l1 = find(us1(i, :) >= usth);
      sp1(l1) = I_0;
      us1(i+1, l1) = 0;
      
      operation(1, 1) += 4;
      
      I_syn1h(:, : , i + layerdelay) = I_syn1h(:, : , i + layerdelay) + sf * sp1.* w1;
      I_syn2h(:, : , i + layerdelay + chaindelay) = I_syn2h(:, : , i + layerdelay + chaindelay) + sf * sp1.* w2;
      I_syn3h(:, : , i + layerdelay + 2 * chaindelay) = I_syn3h(:, : , i + layerdelay + 2 * chaindelay) + sf * sp1.* w3;
      I_syn4h(:, : , i + layerdelay + 3 * chaindelay) = I_syn4h(:, : , i + layerdelay + 3 * chaindelay) + sf * sp1.* w4;
      
      I_synd1h(:, : , i + layerdelay) = I_synd1h(:, : , i + layerdelay) + sf * sp1.* w1;
      I_synd2h(:, : , i + layerdelay + chaindelay) = I_synd2h(:, : , i + layerdelay + chaindelay) + sf * sp1.* w2;
      I_synd3h(:, : , i + layerdelay + 2 * chaindelay) = I_synd3h(:, : , i + layerdelay + 2 * chaindelay) + sf * sp1.* w3;
      I_synd4h(:, : , i + layerdelay + 3 * chaindelay) = I_synd4h(:, : , i + layerdelay + 3 * chaindelay) + sf * sp1.* w4;   
    end                    
        
    I_syn1h(:, :, (i + layerdelay + 1)) = I_syn1h(:, :, i + layerdelay) - (1 ./ ts) .* I_syn1h(:, :, i + layerdelay) .* dt;
    I_syn2h(:, :, (i + layerdelay + chaindelay + 1)) = I_syn2h(:, :, i + layerdelay + chaindelay) - (1 ./ ts) .* I_syn2h(:, :, i + layerdelay + chaindelay) .* dt;
    I_syn3h(:, :, (i + layerdelay + 2 * chaindelay + 1)) = I_syn3h(:, :, i + layerdelay + 2 * chaindelay) - (1 ./ ts) .* I_syn3h(:, :, i + layerdelay + 2 * chaindelay) .* dt;
    I_syn4h(:, :, (i + layerdelay + 3 * chaindelay + 1)) = I_syn4h(:, :, i + layerdelay + 3 * chaindelay) - (1 ./ ts) .* I_syn4h(:, :, i + layerdelay + 3 * chaindelay) .* dt;
    
    I_synd1h(:, : , (i + layerdelay + 1)) = I_synd1h(:, :, i + layerdelay) - (1 ./ ts) .* I_synd1h(:, :, i + layerdelay) .* dt;
    I_synd2h(:, : , (i + layerdelay + chaindelay + 1)) = I_synd2h(:, :, i + layerdelay + chaindelay) - (1 ./ ts) .* I_synd2h(:, :, i + layerdelay + chaindelay) .* dt;
    I_synd3h(:, : , (i + layerdelay + 2 * chaindelay + 1)) = I_synd3h(:, :, i + layerdelay + 2 * chaindelay) - (1 ./ ts) .* I_synd3h(:, :, i + layerdelay + 2 * chaindelay) .* dt;
    I_synd4h(:, : , (i + layerdelay + 3 * chaindelay + 1)) = I_synd4h(:, :, i + layerdelay + 3 * chaindelay) - (1 ./ ts) .* I_synd4h(:, :, i + layerdelay + 3 * chaindelay) .* dt;

    ## EPSC integration
    ud1(:, :, i+1) = ud1(:, :, i) + (1 ./ td) .* (-ud1(:, :, i) + Rm .* I_synd1h(:, :, i)) .* dt;
    ud2(:, :, i+1) = ud2(:, :, i) + (1 ./ td) .* (-ud2(:, :, i) + Rm .* I_synd2h(:, :, i)) .* dt;
    ud3(:, :, i+1) = ud3(:, :, i) + (1 ./ td) .* (-ud3(:, :, i) + Rm .* I_synd3h(:, :, i)) .* dt;
    ud4(:, :, i+1) = ud4(:, :, i) + (1 ./ td) .* (-ud4(:, :, i) + Rm .* I_synd4h(:, :, i)) .* dt;

    ush(i+1, :) = ush(i, :) + (1 ./ tm) .* (-ush(i, :) + Rm .* (sum(I_syn1h(:, :, i)) + sum(I_syn2h(:, :, i)) + sum(I_syn3h(:, :, i)) + sum(I_syn4h(:, :, i)))) .* dt;
    
    if(max(ush(i, :)) >= usth)     
      lh = find(ush(i, :)>=usth);
      sph(lh) = I_0;
      ush(i+1, lh) = 0;
      
      ## SynOps
      operation(1, 2) += 1;
      
      ## lateral inhibition
      ush(i+1, :) -= 0.9 * usth;           
      ush(i+1, find(ush(i+1, :)<0)) = 0;
      
      I_synhf(:, : , i + layerdelay) = I_synhf(:, : , i + layerdelay) + sf * sph.* wh;
      
      l = zeros(n1, nh);
      l(:, lh) = 1;
      l(find(ud1(:, :, i)<udth2)) = -l(find(ud1(:, :, i)<udth2)) * dw;
      l(find(ud1(:, :, i)>=udth2)) = l(find(ud1(:, :, i)>=udth2)) * dw;
      w1 += l;
      w1(find(w1>maxw)) = maxw;
      w1(find(w1<0)) = 0;    
      
      l = zeros(n1, nh);
      l(:, lh) = 1;
      l(find(ud2(:, :, i)<udth2)) = -l(find(ud2(:, :, i)<udth2)) * dw;
      l(find(ud2(:, :, i)>=udth2)) = l(find(ud2(:, :, i)>=udth2)) * dw;
      w2 += l;
      w2(find(w2>maxw)) = maxw;
      w2(find(w2<0)) = 0;      
      
      l = zeros(n1, nh);
      l(:, lh) = 1;
      l(find(ud3(:, :, i)<udth2)) = -l(find(ud3(:, :, i)<udth2)) * dw;
      l(find(ud3(:, :, i)>=udth2)) = l(find(ud3(:, :, i)>=udth2)) * dw;
      w3 += l;
      w3(find(w3>maxw)) = maxw;
      w3(find(w3<0)) = 0;

      l = zeros(n1, nh);
      l(:, lh) = 1;
      l(find(ud4(:, :, i)<udth2)) = -l(find(ud4(:, :, i)<udth2)) * dw;
      l(find(ud4(:, :, i)>=udth2)) = l(find(ud4(:, :, i)>=udth2)) * dw;
      w4 += l;
      w4(find(w4>maxw)) = maxw;
      w4(find(w4<0)) = 0;      
    end    

    ## update the synaptic current
    I_synhf(:, :, (i + layerdelay + 1)) = I_synhf(:, :, i + layerdelay) - (1 ./ ts) .* I_synhf(:, :, i + layerdelay) .* dt;
    
    ## EPSC integration
    udh(:, :, i+1) = udh(:, :, i) + (1 ./ td) .* (-udh(:, :, i) + Rm .* I_synhf(:, :, i)) .* dt;    
    usf(i+1, :) = usf(i, :) + (1 ./ tm) .* (-usf(i, :) + Rm .* (Iinj2(1, :))) .* dt + (1 ./ tm) .* (Rm .* (sum(I_synhf(:, :, i)))) .* dt;
        
    if(max(usf(i, :)) >= usth)
      lf = find(usf(i, :)>=usth);
      usf(i+1, lf) = 0;
      
      ## SynOps
      operation(1, 3) += 1;
     
      ## lateral inhibition
      usf(i+1, :) -= LI*usth;           
      usf(i+1, find(usf(i+1, :)<0)) = 0;

      l = zeros(nh, nf);
      l(:, lf) = 1;
      l(find(udh(:, :, i)<udth2)) = -l(find(udh(:, :, i)<udth2)) * dw * 0.6;
      l(find(udh(:, :, i)<udth1)) = -l(find(udh(:, :, i)<udth1)) * 0;
      l(find(udh(:, :, i)>=udth2)) = l(find(udh(:, :, i)>=udth2)) * dw;
      wh += l;

      wh(find(wh>maxw2)) = maxw2;
      wh(find(wh<0)) = 0;      
    end   
  end
end
