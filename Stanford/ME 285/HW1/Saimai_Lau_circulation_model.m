% Please change the name of this file and of the main function to include
% your name.  
function Saimai_Lau_circulation_model
close all,clc

% This quantity is needed to convert data from the literature in cgs units
% to the units used in this file.
pressure_conversion=1333.2237; % number of dyne*sec/cm^5 in 1 mmHg

% The following parameters were taken from the literature for this model.
% The units and the authors of the work from which they are taken are
% listed.
L_la=0.00005; % mmHg s^2/ml, Danielsen and Ottesen 2004
L_lv=0.000416; % mmHg s^2/ml, Danielsen and Ottesen 2004
v_lv0=2; % ml, Danielsen and Ottesen 2004
v_d0=205+370+401; % ml, Danielsen and Ottesen 2004
R_c=90/pressure_conversion; % mmHg s/ml, Westerhof, Elzinga, Sipkema 1971
R_d=1200/pressure_conversion; % mmHg s/ml, Westerhof, Elzinga, Sipkema 1971
C_d=8e-4*pressure_conversion; % ml/mmHg, Westerhof, Elzinga, Sipkema 1971

% The following parameters are given to approximate actual physiology.
p_la = 12; % mmHg, chosen to give the desired results
initial_v_lv = 128.35607011204; % ml
initial_v_d = 1053.96962417265; % ml

% Set the number of cycles for the first set of conditions.
num_cycles=8;

% Run the simulation once for the initial set of conditions and then for the
% reduced value of R_d.
for first_set_of_conditions=[true false]
	if ~first_set_of_conditions % i.e. if this is not the first simulation
		% Lower R_d for the second set of conditions.
		R_d=R_d*2/3;
		% Run the second set of conditions for 20 cycles.
		num_cycles=20;
	end

	% Start the simulation at time 0.
	initial_t=0;

	% Use the following order for the state variables.
	% y=[Q_m v_lv Q_a v_d]

	% Initialize variables to store the solution.
	Q_m_soln=0;
	v_lv_soln=initial_v_lv;
	Q_a_soln=0;
	v_d_soln=initial_v_d;
	t_soln=initial_t;

    final_t = 30;
    last_contraction_indx = 0;

	for cycle=1:num_cycles
%         if ~first_set_of_conditions
%             fprintf("looping\n")
%         end
		% The code inside this "for" loop will be executed num_cycles times.
		% Because of the shape of the elastance vs. time curve, each cardiac
		% cycle must begin with isovolumic contraction.

		% Please use a 'Maxstep' of 0.01 as in the sample system.
		
		% Solve the system for isovolumic contraction.
%         fprintf("Isovolumic contraction, t= %f.\n", t_soln(end))
        last_contraction_indx = length(t_soln)+1;
        if Q_a_soln(end)*R_c+1/C_d*(v_d_soln(end)-v_d0)-elastance(t_soln(end))*(v_lv_soln(end)-v_lv0) > 0
%         last_contraction_indx = length(t_soln)+1;
        options = odeset('Events', @(t,y)valve_a_opening_test(t,y,v_lv0,R_c,C_d,v_d0),'Maxstep',0.01);
		[t_ode_soln,y_ode_soln,t_event,y_event,i_event] = ... 
            ode45(@(t,y)dydt_m_closed_a_closed(t,y,p_la,L_la,R_d,v_d0,R_c,C_d,v_lv0,L_lv), ... 
            [t_soln(end) final_t], ... % desired time limits of solution
            [Q_m_soln(1) v_lv_soln(end) Q_a_soln(1) v_d_soln(end)], ... % initial values of state variables
            options); % ODE solver options defined above

        t_soln = [t_soln; t_ode_soln(2:end); t_event];
        Q_m_soln = [Q_m_soln; y_ode_soln(2:end,1); y_event(1)];
        v_lv_soln = [v_lv_soln; y_ode_soln(2:end,2); y_event(2)];
        Q_a_soln = [Q_a_soln; y_ode_soln(2:end,3); y_event(3)];
        v_d_soln = [v_d_soln; y_ode_soln(2:end,4); y_event(4)];
        end

		% Solve the system for ejection.
%         fprintf("Ejection, t= %f.\n", t_soln(end))
        options = odeset('Events', @(t,y)valve_a_closing_test(t,y,v_lv0,R_c,C_d,v_d0),'Maxstep',0.01);
		[t_ode_soln,y_ode_soln,t_event,y_event,i_event] = ... 
            ode45(@(t,y)dydt_m_closed_a_open(t,y,p_la,L_la,R_d,v_d0,R_c,C_d,v_lv0,L_lv), ... 
            [t_soln(end) final_t], ... % desired time limits of solution
            [Q_m_soln(1) v_lv_soln(end) Q_a_soln(end) v_d_soln(end)], ... % initial values of state variables
            options); % ODE solver options defined above

        t_soln = [t_soln; t_ode_soln(2:end); t_event];
        Q_m_soln = [Q_m_soln; y_ode_soln(2:end,1); y_event(1)];
        v_lv_soln = [v_lv_soln; y_ode_soln(2:end,2); y_event(2)];
        Q_a_soln = [Q_a_soln; y_ode_soln(2:end,3); y_event(3)];
        v_d_soln = [v_d_soln; y_ode_soln(2:end,4); y_event(4)];
		
		% Solve the system for isovolumic relaxation.
%         fprintf("Isovolumic relaxation, t= %f.\n", t_soln(end))
        options = odeset('Events', @(t,y)valve_m_opening_test(t,y,v_lv0,p_la),'Maxstep',0.01);
		[t_ode_soln,y_ode_soln,t_event,y_event,i_event] = ... 
            ode45(@(t,y)dydt_m_closed_a_closed(t,y,p_la,L_la,R_d,v_d0,R_c,C_d,v_lv0,L_lv), ... 
            [t_soln(end) final_t], ... % desired time limits of solution
            [Q_m_soln(1) v_lv_soln(end) Q_a_soln(1) v_d_soln(end)], ... % initial values of state variables
            options); % ODE solver options defined above

        t_soln = [t_soln; t_ode_soln(2:end); t_event];
        Q_m_soln = [Q_m_soln; y_ode_soln(2:end,1); y_event(1)];
        v_lv_soln = [v_lv_soln; y_ode_soln(2:end,2); y_event(2)];
        Q_a_soln = [Q_a_soln; y_ode_soln(2:end,3); y_event(3)];
        v_d_soln = [v_d_soln; y_ode_soln(2:end,4); y_event(4)];
		
		% Solve the system for filling.
%         fprintf("Filling, t= %f.\n", t_soln(end))
        options = odeset('Events', @(t,y)valve_m_closing_test(t,y,v_lv0,p_la),'Maxstep',0.01);
		[t_ode_soln,y_ode_soln,t_event,y_event,i_event] = ... 
            ode45(@(t,y)dydt_m_open_a_closed(t,y,p_la,L_la,R_d,v_d0,R_c,C_d,v_lv0,L_lv), ... 
            [t_soln(end) final_t], ... % desired time limits of solution
            [Q_m_soln(end) v_lv_soln(end) Q_a_soln(1) v_d_soln(end)], ... % initial values of state variables
            options); % ODE solver options defined above

        t_soln = [t_soln; t_ode_soln(2:end); t_event];
        Q_m_soln = [Q_m_soln; y_ode_soln(2:end,1); y_event(1)];
        v_lv_soln = [v_lv_soln; y_ode_soln(2:end,2); y_event(2)];
        Q_a_soln = [Q_a_soln; y_ode_soln(2:end,3); y_event(3)];
        v_d_soln = [v_d_soln; y_ode_soln(2:end,4); y_event(4)];
		
	end
	
	% Create the required plots and compute the requested pressure and work
    if first_set_of_conditions
        condit = "(Normal)";
    else
        condit = "(Injured)";
    end
	% outputs.
%     t_soln = t_soln(1:end-1);
    figure()
    plot(t_soln,Q_a_soln,'r-'), hold on
    plot(t_soln,Q_m_soln,'b-')
    legend('Q_a','Q_m',Location='southeast')
    title("Q vs. time "+condit); % provides a title
    xlabel('time (seconds)'); % labels the x axis
    ylabel('flow rate (ml/s)'); % labels the y axis

    figure()
    p_lv_soln = elastance(t_soln).*(v_lv_soln-v_lv0);
    p_d_soln = 1/C_d.*(v_d_soln-v_d0);
    plot(t_soln,p_lv_soln,'r-'), hold on
    plot(t_soln,p_d_soln,'b-')
    legend('P_{lv}','P_d')
    title("p vs. time "+condit); % provides a title
    xlabel('time (seconds)'); % labels the x axis
    ylabel('pressure (mmHg)'); % labels the y axis
%     xlim([0 1.2])
%     ylim([-1 300])

    figure()
    plot(v_lv_soln,p_lv_soln,'r-')
    title("PV loop for left ventricle "+condit); % provides a title
    xlabel('v_{lv} (ml)'); % labels the x axis
    ylabel('p_{lv} (mmHg)'); % labels the y axis
%     xlim([0 150])
%     ylim([0 250])
    
    work = 0.0;
    for i=length(v_lv_soln):-1:last_contraction_indx
        work = work - (p_lv_soln(i)+p_lv_soln(i-1))*(v_lv_soln(i)-v_lv_soln(i-1))/2;
    end
    fprintf("Maximum p_d "+condit+" = %f.\n", max(p_d_soln))
    fprintf("Minimum p_d "+condit+" = %f.\n", min(p_d_soln))
    fprintf("Work done by the last cycle "+condit+" = %f.\n", work)
    
end
end



% End of main function.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create functions that will be used by the ODE solver to evaluate the
% time derivatives of the state variables.

function dydt = dydt_m_closed_a_open(t_current,y_current,p_la,L_la,R_d,v_d0,R_c,C_d,v_lv0,L_lv)
dydt(1,1) = 0;
dydt(2,1) = y_current(1)-y_current(3);
dydt(3,1) = 1/L_lv*(elastance(t_current)*(y_current(2)-v_lv0)-y_current(3)*R_c-1/C_d*(y_current(4)-v_d0));
dydt(4,1) = y_current(3)-1/R_d/C_d*(y_current(4)-v_d0);
end

function dydt = dydt_m_open_a_open(t_current,y_current,p_la,L_la,R_d,v_d0,R_c,C_d,v_lv0,L_lv)
dydt(1,1) = 1/L_la*(p_la-elastance(t_current)*(y_current(2)-v_lv0));
dydt(2,1) = y_current(1)-y_current(3);
dydt(3,1) = 1/L_lv*(elastance(t_current)*(y_current(2)-v_lv0)-y_current(3)*R_c-1/C_d*(y_current(4)-v_d0));
dydt(4,1) = y_current(3)-1/R_d/C_d*(y_current(4)-v_d0);
end

function dydt = dydt_m_open_a_closed(t_current,y_current,p_la,L_la,R_d,v_d0,R_c,C_d,v_lv0,L_lv)
dydt(1,1) = 1/L_la*(p_la-elastance(t_current)*(y_current(2)-v_lv0));
dydt(2,1) = y_current(1)-y_current(3);
dydt(3,1) = 0;
dydt(4,1) = y_current(3)-1/R_d/C_d*(y_current(4)-v_d0);
end

function dydt = dydt_m_closed_a_closed(t_current,y_current,p_la,L_la,R_d,v_d0,R_c,C_d,v_lv0,L_lv)
dydt(1,1) = 0;
dydt(2,1) = y_current(1)-y_current(3);
dydt(3,1) = 0;
dydt(4,1) = y_current(3)-1/R_d/C_d*(y_current(4)-v_d0);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create functions that will be used by the ODE solver to test for the
% openings and closings of valves.

function [value,stop,dir] = valve_m_opening_test(t_current,y_current, v_lv0,p_la)
dp_m_valve = elastance(t_current)*(y_current(2)-v_lv0)-p_la;
value = dp_m_valve;
stop = 1;
dir = -1;
end

function [value,stop,dir] = valve_a_opening_test(t_current,y_current, v_lv0,R_c,C_d,v_d0)
dp_a_valve = y_current(3)*R_c+1/C_d*(y_current(4)-v_d0)-elastance(t_current)*(y_current(2)-v_lv0);
value = dp_a_valve;
stop = 1;
dir = -1;
end

function [value,stop,dir] = valve_m_closing_test(t_current,y_current, v_lv0,p_la)
dp_m_valve = p_la-elastance(t_current)*(y_current(2)-v_lv0);
value = dp_m_valve;
stop = 1;
dir = -1;
end

function [value,stop,dir] = valve_a_closing_test(t_current,y_current, v_lv0,R_c,C_d,v_d0)
dp_a_valve = elastance(t_current)*(y_current(2)-v_lv0)-(y_current(3)*R_c+1/C_d*(y_current(4)-v_d0));
value = dp_a_valve;
stop = 1;
dir = -1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function E = elastance(t)
% Return the elastance value and its derivative at a time or a column
% vector of times t from the normal data in the following paper:
% Single-Beat Estimation of End-Systolic Pressure-Volume Relation in
% Humans, A New Method With the Potential for Noninvasive Application
% Authors: Hideaki Senzaki, MD; Chen-Huan Chen, MD; David A. Kass, MD
% Circulation 1996; 94:2497-2506

% pi/2 was added to this table to cause sin()=1 for the first term.
Senzaki_table_4=[28.38975 pi/2;
    37.58583 .08367674;
    21.02345 -1.486758;
    7.665592 2.865675;
    4.809436 .1677238;
    4.181973 4.630239;
    1.940692 3.088379;
    .5870049 -.3053668;
    1.181256 4.410703;
    .84039 3.181538;
    .02259011 1.242886;
    .3071458 4.156753;
    .3226207 2.946186];
% Table 1 states that normal subjects had heart rates of 95 beats per
% minute.
hrt_rate=95;
% Table 1 states that normal subjects had end-systolic elastances of
% 2.0 mmHg/ml.
E_es=2.0;

% Compute the elastance values from the sine series.
period=60/hrt_rate;
t_sin=t/period*2*pi;
E=E_es*sum(10^-2*ones(size(t))*Senzaki_table_4(:,1)'.* ...
    sin(t_sin*[0:12]+ones(size(t))*Senzaki_table_4(:,2)'),2);
end

