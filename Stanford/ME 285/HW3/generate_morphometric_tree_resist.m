% This function initializes the creation of a morphometry-based tree from the rat
% pulmonary arterial data contained in the following paper.

% Z. L. Jiang, G. S. Kassab, and Y. C. Fung, "Diameter-defined Strahler system and
% connectivity matrix of the pulmonary arterial tree," J Appl Physiol, v. 76, p. 882-892,
% 1994.

% It is given the order of the initial element and returns the vascular resistance to of
% the arterial tree that is generated.  All calculations are performed with CGS
% (centimeters, grams, and seconds) units.

function generate_morphometric_tree(order)

% Create global matrices to store the morphometric data from Jiang et al. in CGS units.
global connectivity length diameter viscosity
connectivity=[...
    0.19 4.06 2.48 1.36 0 0 0 0 0 0 0;
    0 0.12 1.55 0.93 0 0 0 0 0 0 0;
    0 0 0.17 2.26 0.31 0.05 0.11 0.13 0 0 0;
    0 0 0 0.05 2.00 0.70 0.62 0.43 0.39 0 0;
    0 0 0 0 0.18 1.92 1.56 0.85 0.83 0.50 0;
    0 0 0 0 0 0.18 2.05 1.15 1.00 0.67 0;
    0 0 0 0 0 0 0.05 1.55 1.06 0.67 0;
    0 0 0 0 0 0 0 0.08 1.67 0.83 6;
    0 0 0 0 0 0 0 0 0.06 1.67 4;
    0 0 0 0 0 0 0 0 0 0.17 2;
    0 0 0 0 0 0 0 0 0 0 0;
    ];
length=[0.005 0.015 0.020 0.027 0.040 0.072 0.124 0.133 0.174 0.264 1.811]; % centimeters
diameter=[0.00133 0.00317 0.00444 0.00615 0.00881 0.0152 0.0266 0.0417 0.0602 0.0929 ...
    0.1639]; % centimeters
viscosity=0.04; % Poise

% Create and initialize global matrices to store the remainders of the number of child
% elements created and the total number of elements of each order.  Create a remainder
% for each entry in the connectivity matrix so that the generated vascular tree has the
% given connectivity matrix.
global remainder total_elements_created
remainder=zeros(11,11);
total_elements_created=zeros(11,1);
total_elements_created(order)=1;

% Call calculate_resistance to recursively calculate the vascular resistance, and display
% the resulting number of elements and total resistance.
resistance=calculate_resistance(order);
display(' ');
for display_order=1:11
    display(['Number of order ' num2str(display_order) ' elements: ' ...
        num2str(total_elements_created(display_order))]);
end
display(['Resistance: ' num2str(resistance) ' dyne*s/cm^5']); % CGS units
display(' ');

return