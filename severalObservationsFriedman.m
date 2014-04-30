% SEVERALOBSERVATIONSFRIEDMAN friedman's procedure in case of several
%    observations per experimental unit. 
%       
%    [SIG, P] = SEVERALOBSERVATIONSFRIEDMAN( X, SIGLEVEL ) returns the 
%       significance the friedman's procedure found. If there exists a 
%       significant difference at the level SIGLEVEL, P gives
%       an incidence matrix of all significantly different pairs of 
%       treatments.
%    
%    X is a three dimensional matrix [ x_{ijn} ]_{b x k x m} containing
%       b blocks of
%       k mutually independent treatments with
%       m observations
%
%    SIGLEVEL denotes the significance level, default 0.05
%    the block design needs to be complete.
%
%    According to "Practical Nonparametric Statistics"
%    Third Edition, W.J. Connor. Wiley Series in Probability
%    and Statistics, pages 369-384

function [sig, P] = severalObservationsFriedman( X, sigLevel )
    P = [];
    if( nargin < 2 )
        sigLevel = 0.05;
    end
    [b,k,m] = size(X); % get nr of blocks, treatments, and observations
    
    R = zeros(1,k);    % sums the rank of 

    v = 0;             % used to determine the variance of R
    
    % determine ranking for all blocks and sum up the sum of ranks
    for i = 1 : b
        % retrieve observations
        obs = reshape( X(i,:,:), k*m,1 );
        % get ranks of the observations
        r = tiedrank( obs );
        % get the sum of ranks
        R = R + [sum( reshape(r,m,k) )];
        % update the variance information
        v = v + sum(r.^2);
    end
    
    % determine the expected value of R
    Rexp = b*m*(m*k+1)/2;
    
    % determine the expected variance of R
    Rvar = m*(k-1)/(k*(m*k-1))* ( v - m*k*b*(m*k+1)^2/4 );
    
    % determine the test statistic
    T6 = sum( (k-1)/k* ( R - Rexp ).^2/Rvar );
    
    % determine the significance of the Friedman test
    sig = 1 - chi2cdf(T6, k-1 );
    
    % given the difference is significant, find the significant differences
    % post-hoc
    
    if( sig < sigLevel )
        % Differences in rank sums
        Rdiff = abs( repmat(R,k,1) - repmat(R',1,k) );

        % Test Statistic
        T = tinv( 1-sigLevel/2, m*b*k -k -b + 1 )* ...
            ( ( 2*k*b*(m*k-1)*Rvar ) / ( (k-1)*(m*b*k -k -b + 1) )* ...
            (1 - T6/( b*(m*k-1) ) ) )^(1/2);

        P = Rdiff > T;
    end
    
    
    
    