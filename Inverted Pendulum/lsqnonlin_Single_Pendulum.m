[Bp1,Jp1] = runtracklsq;

function [Bp,Jp] = runtracklsq
mdl = 'EEE4022S_Single_Pendulum_lsqnonlin';
open_system(mdl)                             % Load the model
in = Simulink.SimulationInput(mdl);          % Create simulation input object
in = in.setModelParameter('StopTime','100'); % Stop time 100
Damp0 = [0.002,0.007];                       % Initial values

options = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt',...
   'Display','iter','StepTolerance',0.001,'OptimalityTolerance',0.001);

set_param(mdl,'FastRestart','on');           % Fast restart
Damp = lsqnonlin(@tracklsq,Damp0,[],[],options);
set_param(mdl,'FastRestart','off');

Bp = Damp(1);
Jp = Damp(2);

    function F = tracklsq(Damp)
      % Track the output of optsim to a signal of 1
      % Set the simulation input object parameters
      in = in.setVariable('B_p',Damp(1),'Workspace',mdl);
      in = in.setVariable('J_p',Damp(2),'Workspace',mdl);
      
      % Simulate
      out = sim(in);
      F = out.get('resnorm').Data;
      
    end
end