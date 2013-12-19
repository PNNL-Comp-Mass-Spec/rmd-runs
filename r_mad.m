function [MAD]= r_mad(Z)

% DESCRIPTION
    % This function calculates the robust estimate of the scale parameter
    % sigma using the median absolute deviation.

% INPUT
    % Z  
        % A (n x q) matrix for which n is the number of observations (rows) 
        % and q is the number of variables (columns).

% OUTPUT
    % MAD
        % A row vector containing the scale parameter estimate (sigma) 
        % for each variable in Z (1 x q)

% OTHER FUNCTIONS REQUIRED
    % Statistics Toolbox
        % Obtain from MathWorks, Inc.


MAD = mad(Z,1) .* 1.4826;

