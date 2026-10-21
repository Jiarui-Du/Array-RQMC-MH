function sobol_base = array_rqmc_sobol_base()
%ARRAY_RQMC_SOBOL_BASE Return the paper experiment's Sobol point set.
%
% Returns an unscrambled three-dimensional Sobol object. The kernel scrambles
% it afresh at each transition: coordinate 1 matches sorted particles, while
% coordinates 2 and 3 drive the proposal and acceptance decision.


    %% Construct the shared base net configuration
    sobol_base = sobolset(3, ...
        'Skip', 0, ...
        'Leap', 0, ...
        'PointOrder', 'graycode');
end
