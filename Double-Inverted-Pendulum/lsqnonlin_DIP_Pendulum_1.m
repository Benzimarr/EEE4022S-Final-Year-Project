[Bp11,Jp11] = runtracklsq; 

function [Bp1,Jp1] = runtracklsq
mdl = 'ShortPenModellingSIM';
open_system(mdl)                            % Load the model
in = Simulink.SimulationInput(mdl);         % Create simulation input object
in = in.setModelParameter('StopTime','60'); % Stop time 60
Damp0 = [0.0004,0.0009];                    % Initial values 

options = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt',...
   'Display','iter','StepTolerance',0.001,'OptimalityTolerance',0.001);
% Optimize the Damping and Inertia
set_param(mdl,'FastRestart','on');           % Fast restart
Damp = lsqnonlin(@tracklsq,Damp0,[],[],options);
set_param(mdl,'FastRestart','off');
% Return the Damping and Inertia
Bp1 = Damp(1);
Jp1 = Damp(2);

    function F = tracklsq(Damp)
      % Track the output of optsim to a signal of 1
      % Set the simulation input object parameters
      in = in.setVariable('B_p1',Damp(1),'Workspace',mdl);
      in = in.setVariable('J_p',Damp(2),'Workspace',mdl);
      
      % Simulate
      out = sim(in);
      F = out.get('resnorm').Data;
            
    end
end