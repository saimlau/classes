% In MATLAB, the main function (in this case, simulateTunedTrackSystem)
% must match the file name exactly. To that end, do not rename the file!
% Do not modify the file except in as noted in the two steps below:
% Step 1: Fill in the two lines where noted in tunedTrack_equationsOfMotion
% Step 2: Fill in the three lines where noted in simulateTunedTrackSystem

function simulateTunedTrackSystem
    % ODE parameters (end time, initial conditions)
    t_final = 5;
    xr_initial = 0;
    xrd_initial = 0.1;
    xt_initial = 0;
    
    % System parameters
    m = 50;
    km = 1000;
    zeta = 0.55;
    b = 2*sqrt(m*km)*zeta;
    %-- FILL IN THE THREE LINES BELOW, e.g., kt_compliant = 0.1*km; --%
    kt_compliant = km*1e-1; %-- FILL IN HERE --%
    kt_tuned = 3*km; %-- FILL IN HERE --%
    kt_stiff = km*1e2; %-- FILL IN HERE --%
    
    % Use Matlab's ode45 to integrate equations of motion
    % Syntax:
    %   [t_solution, X_solution] = ode45( <function handle to equations of motion>, [t_initial, t_final], initial_conditions)
    [time_compliant, X_compliant] = ode45(@(t,X) tunedTrack_equationsOfMotion(t, X, m, km, b, kt_compliant), [0, t_final], [xr_initial, xrd_initial, xt_initial]);
    [time_tuned, X_tuned] = ode45(@(t,X) tunedTrack_equationsOfMotion(t, X, m, km, b, kt_tuned), [0, t_final], [xr_initial, xrd_initial, xt_initial]);
    [time_stiff, X_stiff] = ode45(@(t,X) tunedTrack_equationsOfMotion(t, X, m, km, b, kt_stiff), [0, t_final], [xr_initial, xrd_initial, xt_initial]);
    
    % Plot results
    figure
    subplot(2,1,1)
    plot(time_compliant, X_compliant(:,1)), hold all
    plot(time_tuned, X_tuned(:,1))
    plot(time_stiff, X_stiff(:,1))
    plot([0,t_final], [0,0], 'k--')
    grid on
    ylabel('x_r')
    title('Displacement of runner')
    legend('Compliant track', 'Tuned track', 'Stiff track')
    
    subplot(2,1,2)
    plot(time_compliant, X_compliant(:,3)), hold all
    plot(time_tuned, X_tuned(:,3))
    plot(time_stiff, X_stiff(:,3))
    plot([0,t_final], [0,0], 'k--')
    grid on
    ylabel('x_t')
    title('Displacement of track')
    xlabel('t')
    hT_c_rigid = pi/sqrt(km/m*(1-b^2/(4*m*km)))
end


function Xd = tunedTrack_equationsOfMotion(t, X, m, km, b, kt)
    % This function implements the equations of motion for the 2-degree of
    % freedom tuned track system.
    %
    % Variables:
    %   xr: runner displacement
    %   xrd, xrdd: first and second time derivatives, respectively, of xr
    %      i.e., xrd = d(xr)/dt, and xrdd = d^2(xr)/dt^2
    %   xt: track displacement
    %   xtd: first time derivative of xt
    %      i.e., xtd = d(xt)/dt
    %
    % System state (X), and state time derivative (Xd)
    % X = [xr, d(xr)/dt, xt]
    % Xd = [d(xr)/dt, d^2(xr)/dt^2, d(xt)/dt]
    %
    % The equations of motion from dynamics are:
    %   m(xrdd) = km(xt - xr) + b(xtd - xrd)
    %   km(xt - xr) + b(xtd - xrd) = -kt(xt)
    
    % Get state values from X
    xr  = X(1);
    xrd = X(2);
    xt  = X(3);
    
    % Compute state derivatives
    %-- FILL IN THE TWO LINES BELOW --%
    %-- Remember to use the same variables names from above! --%
    xtd = xrd+km./b.*(xr-xt)-kt/b.*xt; %-- FILL IN HERE --%
    xrdd = -kt/m*xt; %-- FILL IN HERE --%
    
    % Return vector of state derivatives
    Xd(1,1) = xrd;
    Xd(2,1) = xrdd;
    Xd(3,1) = xtd;

end
