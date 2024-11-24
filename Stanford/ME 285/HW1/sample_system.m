function sample_system
% This is the main function for this file.  It solves the sample
% problem given in the assignment and should be used as a guide
% to solving the elastance based systemic circulation problem.

% We are given the following parameter values: 
pressure_conversion = 1333.2237; % dyne*sec/cm^5 to 1 mmHg
p_u = 50; % upstream pressure in mmHg
L_u = 0.001; % upstream inductance in mmHg s^2/ml
R_d = 1200/pressure_conversion; % downstream resistance in mmHg s/ml
v_d0 = 205+370+401; % downstream volume at zero pressure in ml

% The system is started at time 0.
initial_t = 0; % seconds

% The initial pressure has been chosen such that the valve will be
% closed at the start of the simulation.
initial_p_d = 120; % mmHg

% The initial value of the downstream volume, one of our state
% variables, must be determined from the pressure and the elastance
% at the initial time.  The function that returns the elastance at a
% given time is called "elastance" and is defined by the function below.
initial_v_d = v_d0+initial_p_d/elastance(initial_t);

% The flow through the valve is initially zero, because the pressure
% downstream of it is higher than its upstream pressure.
initial_Q_v = 0;

% Set the ending time of the simulation.
final_t = 2; % seconds

% Initialize variables to store the solution.  The solution of this
% problem can be completely defined in terms of the flow through the
% valve and the downstream volume.  Because the ODE
% solver requires a column array with all of the state
% variables, we will combine the flow and downstream volume into a column array called
% "y".
t_soln=initial_t; % time solution
Q_v_soln=initial_Q_v; % flow solution
v_d_soln=initial_v_d; % downstream volume solution

% The ODE solver is controlled by setting its options with a function
% called "odeset".  The command below says that a function will be
% used to test for when the valve opens and that the maximum time step
% size that will be allowed is 0.01 seconds.  "valve_opening_test" is
% called within an anonymous function.  For an explanation of how
% anonymous functions work, see the MATLAB advice file.
options = odeset('Events', ...
	@(t,y)valve_opening_test(t,y,v_d0,p_u),'Maxstep',0.01);

% The ODE solver we will use is "ode23", which implements a Runge-
% Kutta method.  Please note that "..." tells MATLAB to combine the
% current line of code with the next line of actual code.  Numerical
% approximations to the solutions of systems of ODE's are found by
% this and other solvers by using the current values of the state
% variables and expressions for the derivatives of these state
% variables in terms of the state variables themselves and time.

% The parameters that are returned by ode23 are an array of time
% values where the solution was computed called "t_ode_sol" here, an
% array of the y solution values at those times called "y_ode_sol", 
% and information about the event that occurred.  The first
% parameter is the name of the function that is called to find the
% time derivative of the current state variables, which is
% "dydt_valve_closed" here.  The second is an array that gives the
% starting and ending times of the solution.  In this case, we will
% not reach the ending time, because the valve will open, and the
% solver will detect that and stop the solution process.  The third
% parameter is an array containing the initial values of the state
% variables that we set above.  The last is the options set above.
[t_ode_soln,y_ode_soln,t_event,y_event,i_event] = ... 
    ode23(@(t,y)dydt_valve_closed(t,y,R_d,v_d0), ... 
    [initial_t final_t], ... % desired time limits of solution
    [initial_Q_v initial_v_d], ... % initial values of state variables
    options); % ODE solver options defined above

% Store the solution in arrays that contain the complete history of
% the solution so that you can eventually plot the entire solution
% instead of the portion only up to the valve opening event.  This is
% done by concatenating the solution obtained from ode23 with previous
% values.  For examples describing how concatenation works, see the
% MATLAB advice file.  In this case, the previous value
% is the initial value.  "t_ode_soln" and "y_ode_soln", which were
% returned by ode23, contain the initial values and the remainder of
% the solution before the valve event.  "y_ode_soln" contains a column
% for each state variable.  "t_ode_soln(2:end)" represents all entries
% in "t_ode_soln" besides the first.  To concatenate it with previous
% solution values and the value at the valve event, use the brackets
% and semicolons as follows.  The jth column of "y_ode_soln", without
% its initial value is "y_ode_soln(2:end,j)".
t_soln = [t_soln; t_ode_soln(2:end); t_event];
Q_v_soln = [Q_v_soln; y_ode_soln(2:end,1); y_event(1)];
v_d_soln = [v_d_soln; y_ode_soln(2:end,2); y_event(2)];

% The valve has opened, so now we call the ODE solver again, using the
% values at the time of the valve opening, the end values of the
% previous solutions, as initial values.  We are no longer looking for
% any valve-related event, so we can leave that out of the options,
% the return values from ode23, and the concatenation step.  The ode23
% outputs are named the same as above, and they will be replaced with
% the new solutions.  The concatenation step appends the solution from
% the step after the opening of the valve until the end of the
% simulation to the solution found above.
options = odeset('Maxstep',0.01);
[t_ode_soln,y_ode_soln] = ode23( ...
	@(t,y)dydt_valve_open(t,y,p_u,L_u,R_d,v_d0), ...
    [t_soln(end) final_t], [Q_v_soln(end) v_d_soln(end)], options);
t_soln = [t_soln; t_ode_soln(2:end)];
Q_v_soln = [Q_v_soln; y_ode_soln(2:end,1)];
v_d_soln = [v_d_soln; y_ode_soln(2:end,2)];

% Compute the downstream pressure p_d using the elastance and volume.
p_d_soln=elastance(t_soln).*(v_d_soln-v_d0);

% Plot the solution.  See the MATLAB advice file for more
% information on plotting.
figure; % creates a new figure
plot(t_soln,v_d_soln); % plots the downstream volume vs. time
title('Downstream volume vs. time'); % provides a title
xlabel('time (seconds)'); % labels the x axis
ylabel('volume (ml)'); % labels the y axis

% We now create a new figure and plot the downstream PV line for the
% solution, marking the endpoints.
figure; % creates a new figure
plot(v_d_soln,p_d_soln); % plots the downstream volume vs. time
title('Downstream pressure vs. volume'); % provides a title
xlabel('volume (ml)'); % labels the x axis
ylabel('pressure (mmHg)'); % labels the y axis
hold on; % tells the plot not to erase what has been plotted
plot(v_d_soln(end),p_d_soln(end),'rx','MarkerSize',10); % marks the
% end of the PV loop
plot(v_d_soln(1),p_d_soln(1),'ro','MarkerSize',10); % marks the
% beginning of the PV loop


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function E_d = elastance(t)
% This function returns the elastance at a given time or the array of
% elastance values for a given array of times.  The function of time
% given here is not meant to mimic any real life behavior.  It is
% demonstrating that the elastance value can be a function of time.

C_d = 0.7; % ml/mmHg
E = 1/C_d; % mmHg/ml, elastance is the inverse of capacitance.
E_d=E*(1+(t>1).*(t-1)); % This function is equal to E until t=1s and
% then varies linearly between E and 2*E at t=2s.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dydt = dydt_valve_closed(t_current,y_current,R_d,v_d0)
% This function is called by the ODE solver to evaluate the time
% derivatives of the state variables when the valve is closed.  The
% current time is t_current, and the current values of the state
% variables are contained in a vector called y_current.

% The function must return a column array containing the time
% derivatives of the state variables.  The first row contains the time
% derivative of the first state variable, Q_v, and the second row
% contains the time derivative of the second state variable, v_d.

% When the valve is closed, Q_v, which we have chosen to be the first
% entry in y, does not change with time.  It remains at zero, because 
% there is no flow through the valve.  Therefore, its time derivative
% is also zero. We could even avoid solving for this variable, but we
% have chosen to keep it in the system of ODE's for the sake of
% simplicity.
dydt(1,1)=0;

% The downstream volume v_d, which we have chosen to be the second entry
% in y, has a time derivative given by
% d(v_d)/dt=Q_v-E_d(t)*(v_d-v_d0)/R_d.  We can give the current
% time value to the "elastance" function to find E_d.  The current
% value of v_d is equal to y_current(2), and we can even get the value
% of Q_v from y_current(1), even though we know it will be zero.
dydt(2,1)=y_current(1)-elastance(t_current)*(y_current(2)-v_d0)/R_d;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dydt = dydt_valve_open(t_current,y_current,p_u,L_u,R_d,v_d0)
% This function is called by the ODE solver to evaluate the time
% derivatives of the state variables when the valve is open.  The
% current time is t_current, and the current values of the state
% variables are contained in a vector called y_current.

% The function must return a column array containing the time
% derivatives of the state variables.  The first row contains the time
% derivative of the first state variable, Q_v, and the second row
% contains the time derivative of the second state variable, v_d.

% When the valve is open, Q_v, which we have chosen to be the first
% entry in y, has a time derivative given by
% d(Q_v)/dt=1/L_u*(p_u-E_d(t)*(v_d-v_d0)).  We can give the current time
% value to the "elastance" function to find E_d.  The current value of
% v_d is equal to y_current(2), so we have all the necessary variables.
dydt(1,1)=1/L_u*(p_u-elastance(t_current)*(y_current(2)-v_d0));

% The downstream volume v_d, which we have chosen to be the second entry
% in y, has a time derivative given by
% d(v_d)/dt=Q_v-E_d(t)*(v_d-v_d0)/R_d.  We can give the current time
% value to the "elastance" function to find E_d.  The current value of
% v_d is equal to y_current(2), and we can even get the value of Q_v from
% y_current(1), so we're ready to compute d(v_d)/dt.
dydt(2,1)=y_current(1)-elastance(t_current)*(y_current(2)-v_d0)/R_d;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [value,stop,dir] = valve_opening_test(t_current,y_current, ...
	v_d0,p_u)
% This function is called by the ODE solver to test whether or not
% the "dp_valve", the pressure downstream of the valve minus the
% pressure upstream of the valve, changes sign during the solution
% process.  In this system, the valve opens when "dp_valve" is less
% than zero, changing from positive to negative.  When this function
% is called and the ODE solver sees that "dp_valve" is negative, the
% ODE solver will backtrack to figure out when it was equal to zero.

% Compute the current pressure difference given by
% E_d(t)*(v_d-v_d0)-p_u.  We can give the current time value to the
% "elastance" function to find E_d.  The current value of v_d is equal
% to y_current(2), so we have all the necessary variables.
dp_valve = elastance(t_current)*(y_current(2)-v_d0)-p_u;

% The first entry in the list of returned variables above, the list
% to the right of "function" above, is the function that the ODE
% solver is searching for a zero of.  That means that "value" needs
% to be set to the current value of "dp_valve".
value = dp_valve;

% The second returned variable tells the ODE solver whether or not to
% stop when a zero of the first returned variable has been found.
stop = 1;

% The third returned variable tells the ODE solver that it should
% consider only zero crossings in which the first returned variable
% has changed from positive to negative.  If we wanted to detect zero
% crossings in only the other direction, we could change this to +1;
dir = -1;