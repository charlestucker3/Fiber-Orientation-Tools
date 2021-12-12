function Rstar = fldRstar(n, S)
%RSTAR = FLDRSTAR(N, S) returns the NxN matrix RSTAR used in calculating
%     fiber length distributions using the Phelps-Tucker model.  S times
%     the parent fiber length is the standard deviation of a normal
%     distribution for the child-fiber lengths.  Use S << 1 to get the
%     fiber to always break in the middle, and S >> 1 to distribute the
%     breaks evenly along the fiber length.  Values of S less than 0.005
%     will give NaN results, and are not allowed.
%     
%     Multiplying each column of Rstar by the breakage rate for the
%     corresponding fiber length times the time step gives a transition
%     matrix that will advance the fiber length distribution (as a column
%     vector) by one time step.
%
%     The calculation relies on the fiber length for any row or column i
%     being given by i*Delta_L.

% Confirm that the value of S is legal
if S < 0.005
    error('S must be >= 0.005')
end

Rstar = zeros(n);
for col = 2:n
    rows = 1:col-1; % Rows in upper triangle, for this column
    % Normal PDF for this column
    Rstar(rows, col) = exp(-0.5*((rows-(col/2))/(S*col)).^2) ...
                       / (S*col*sqrt(2/pi));
    % Normalize the column values
    Rstar(rows, col) = 2 * Rstar(rows,col) / sum(Rstar(rows,col));
    % Add a -1 on the diagonal, to include the breakage term
    Rstar(col,col) = -1;
end

% Allow the shortest length to have a breakage rate, if P_1 ~= 0.
Rstar(1,1) = -1;

return
    