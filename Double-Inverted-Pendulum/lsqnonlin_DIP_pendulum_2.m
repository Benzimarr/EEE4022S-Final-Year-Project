[Bp21,Jp21] = runtracklsq;

function [Bp2,Jp2] = runtracklsq 
mdl = 'MediumPenModellingSIM';
open_system(mdl)                                % Load the model
in = Simulink.SimulationInput(mdl);             % Create simulation input object
in = in.setModelParameter('StopTime','60');     % Stop time 60
Damp0 = [0.00013,0.0008];                       % Initial gain values

options = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt',...
   'Display','iter','StepTolerance',0.001,'OptimalityTolerance',0.001);
% Optimize the Damping and Inertia
set_param(mdl,'FastRestart','on');           % Fast restart
Damp = lsqnonlin(@tracklsq,Damp0,[],[],options);
set_param(mdl,'FastRestart','off');
% Return the Damping and Inertia
Bp2 = Damp(1);
Jp2 = Damp(2);

    function F = tracklsq(Damp)
      % Track the output of optsim to a signal of 1
      % Set the simulation input object parameters
      in = in.setVariable('B_p2',Damp(1),'Workspace',mdl);
      in = in.setVariable('J_p2',Damp(2),'Workspace',mdl);
      
      % Simulate
      out = sim(in);
      F = out.get('resnorm').Data;
      
    end
end