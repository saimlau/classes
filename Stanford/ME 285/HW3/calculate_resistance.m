% This function creates the children of a parent element of a given order using the data
% of Jiang et al.  It returns the vascular resistance at the inlet node of the parent.
% All calculations are performed with CGS (centimeters, grams, and seconds) units.  This
% method of calculating the impedance is recursive, because this function calls itself.
% Each call of this function is associated with one element and computes its total
% resistance based on the resistance of the downstream elements.

% Z. L. Jiang, G. S. Kassab, and Y. C. Fung, "Diameter-defined Strahler system and
% connectivity matrix of the pulmonary arterial tree," J Appl Physiol, v. 76, p. 882-892,
% 1994.

function resistance=calculate_resistance(parent_order)

% List the global variables so that they can be accessed.
global connectivity length diameter viscosity remainder total_elements_created

% From the study of Poiseuille flow, we have the following formula for
% the resistance of each segment.
%       resistance = 8 * dynamic viscosity * length of the segment / ( pi * radius^4 )
% The resistances of the segments and downstream elements must be combined to obtain the
% resistance at the inlet node of the parent element.  The resistance of each child
% element must be found to compute the resistance of the parent element.  This is done
% by calling the calculate_resistance function with the order of the child element as
% input.  We also use the following results from analagous resistive electrical circuits.
% The total resistance R_t of two resistors with resistance values R_1 and R_2 in series
% is their sum, R_t=R_1+R_2.  The total resistance of two resistors in parallel is
% R_t=1/(1/R_1+1/R_2).

% Calculate the number of child elements of each order by rounding the sum of the mean
% value in table 4 and the previous remainder to the nearest whole number.  Since the
% data on the number of pruned elements is omitted from the paper, their contribution to
% the morphometry-based tree is neglected.  Add the number of child elements to the
% running total of elements created.  Update the remainders based on the number of child
% elements created.
tem1 = connectivity(:,parent_order)+remainder(:,parent_order);
child_counts = round(tem1);
n = sum(child_counts);
total_elements_created = total_elements_created+child_counts;
remainder(:,parent_order) = tem1-child_counts;
l = length(parent_order);
r = diameter(parent_order)/2;
if child_counts(parent_order)>=2
    fprintf("%d same order child, more than 1.\n", child_counts(parent_order))
end

% If there are no children, the vessel is of order 1, and the recursive method has
% "bottomed out".  Therefore the resistance is calculated and returned to the calling
% function, and the rest of this function is not executed.
if sum(child_counts) == 0
    resistance = 8 * viscosity * l / ( pi * r^4 );
    return
end

% Divide the element evenly into segments, assuming that each segment bifurcates.  Each
% point where segments meet is called a node.  Connect the child elements to the parent
% based on the following rules:
%       No child elements are connected to the inlet node, which is either the inlet of
%           the entire tree or connected to the upstream element.
%       If the parent and child elements are both of order 1, change the number of
%           segments in the element to 2, and the child element is connected to the node
%           at the midpoint of the parent element.
%       Connect two child elements to the outlet node of the parent element, choosing
%           first elements of one order smaller than the parent and then decreasing in
%           order.  If vessels of the same order as the parent remain unconnected,
%           connect them at the outlet node.  This method accounts for the definition of
%           an element in the Strahler ordering system terminating where no downstream
%           elements are of the same order.
%       Connect the remainder of the child elements to the interior nodes of the parent
%           element with the higher order child elements closer to the downstream end of
%           the parent element and the lower order elements closer to the upstream end.
if n==1
    if parent_order==1
        segment_count = 2;
        l = l/segment_count;
        tem2 = 8*viscosity*l/(pi*r^4);
        resistance = tem2 + 1/(1/tem2+1/calculate_resistance(1));
        return
    else
        % segment_count = 1;
        resistance = 8*viscosity*l/(pi*r^4)+calculate_resistance(parent_order-1);
        return
    end
else
    tem3 = child_counts;
    tem3(parent_order) = 0;
    tem4 = 0.0;
    try
        idx1 = find(tem3,1,'last');
        tem3(idx1) = tem3(idx1)-1;
        tem4 = tem4+1/calculate_resistance(idx1);
        idx2 = find(tem3,1,'last');
        tem3(idx2) = tem3(idx2)-1;
        tem4 = tem4+1/calculate_resistance(idx2);
    catch
    end
    for i=1:child_counts(parent_order)
        tem4 = tem4+1/calculate_resistance(parent_order);
    end
    m = sum(tem3);
    segments_count = m+1;
    l = l/segments_count;
    tem2 = 8*viscosity*l/(pi*r^4);
    resistance = tem2+1/tem4;
    for i=1:m
        idx = find(tem3,1,'last');
        resistance = tem2+1/(1/resistance+1/calculate_resistance(idx));
        tem3(idx) = tem3(idx)-1;
    end
    if sum(tem3)>0
        fprintf("something wrong, unassigned child")
    end
end
return