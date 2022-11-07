[Bc1] = runtracklsq;

function [Bc] = runtracklsq
mdl = 'CartModelling';
open_system(mdl)                                % Load the model
in = Simulink.SimulationInput(mdl);             % Create simulation input object
in = in.setModelParameter('StopTime','30');     % Stop time 60
Damp0 = [1];                                    % Initial values

options = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt',...
   'Display','iter','StepTolerance',0.001,'OptimalityTolerance',0.001);
% Optimize the damping
set_param(mdl,'FastRestart','on');           % Fast restart
Damp = lsqnonlin(@tracklsq,Damp0,[],[],options);
set_param(mdl,'FastRestart','off');
% Return the Damping
Bc = Damp(1);

    function F = tracklsq(Damp)
      % Track the output of optsim to a signal of 1
      % Set the simulation input object parameters
      in = in.setVariable('B_c',Damp(1),'Workspace',mdl);
      
      % Simulate
      out = sim(in);
      F = out.get('resnorm').Data;
      
    end
end