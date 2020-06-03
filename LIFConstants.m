classdef LIFConstants
   properties (Constant)
   Rm = 1e4;          # resistance (Ohm)
   tm = 20E-3;        # membrane potential time constant (s)
   td = 15E-3;        # dendritic potential time constant (s)
   ts = 15E-3;        # EPSC time constant (s)
   usth = 10E-3;      # spike threshold (V)
   udth1 = 0.05E-3;   # dendritic threshold for LTD (V)
   udth2 = 1E-3;      # Threshold for a bAP boost (V)
   LI = 0.9;          # lateral inhibition
   sf = 1.5;          # scaling factor
   I_0 = 1.62E-6;     # synaptic current
   chaindelay = 100;  # axonal delays in synaptic chain (ms)
   layerdelay = 20;   # axonal delays between layers (ms)
   period = 100;      # input period (ms)
   MemorySize = 500;  # synaptic chain size
   BufferSize = 321;  # synaptic chain buffer size
   end
end