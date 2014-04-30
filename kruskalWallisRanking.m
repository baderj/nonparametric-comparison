%KRUSKALWALLISRANKING Ranking of algorithms
%   [KW,L] = KRUSKALWALLISRANKING( X, ALPHA ) ranks the samples
%   corresponding to the columns of X according to the 
%   number of other samples it is significantly better to.
%  
%   X is a matrix of size MAX(N_i) x K, where N_i denotes the
%       number of instances in sample i. The remaining elements
%       need to be set to NaN.
%   ALPHA denotes the significance level.
%
%   KW represents the ranking where
%       the first column goes from 1 to K
%       the second column gives the rank, zero meaning no other
%            algorithm is significantly better, one meaning one 
%            algorithm is better etc.
%       the third column gives the index of the corresponding 
%            algorithm in the data array.
%
%   For example
%       [KW,L] = KRUSKALWALLISRANKING( X, 0.1)
%   with
%   X = [83 71 101; 
%        91 70 100; 
%        94 NaN 91]
%   
%   Ranks the three samples of size 3,2 and 3 respectively 
%
%   
%
% According to "Practical Nonparametric Statistics"
% Third Edition, W.J. Connor. Wiley Series in Probability
% and Statistics, pages 288-290
%
% Uses the Conover-Inman procedure, Fisher's least significant 
% difference method performed on ranks.
%
% Tested using the example on page 291f in the same book


function [KW,L] = kruskalWallisRanking( X, alpha )
    % Number of Samples
    k = size(X,2);
    assert( size(X,1) > 1, 'Sample size needs to be at least 2');  
    % Sample size
    n = size(X,1) - sum( isnan(X) );
    assert( min(n) > 1, 'Sample size needs to be at least 2' );
      
    % Total number of observations
    N = sum(n);
    
    % Rank all observations together
    Ranks = rankWithTies( X, n );
    R = zeros(1,k);
    for i = 1 : k
        R(i) = sum( Ranks( find( Ranks(:,3) == i ) , 1 ) );
    end
        
    % Test Statistics
    S = sqrt( 1/(N-1)*( sum( Ranks(:,1).^2 ) - N*(N+1)^2/4 ));
    T = 1/S^2*( sum( R.^2./n) - N*(N+1)^2/4 );
    
    % Null Distribution - chi-squared approximation
    p = 1 - chi2cdf( T, k-1 );
    
    if( p < alpha )
        % There exists a difference, test all pairs
        Pairs = abs( repmat( R./n, k, 1) - repmat( [R./n]', 1, k ) );        
        t = tinv( 1-alpha/2, N-k);
        C = t*( S^2*(N-1-T)/(N-k) )^(0.5)* ...
            sqrt( ( repmat( 1./n, k, 1) + repmat( [1./n]', 1, k ) ) );
        Sig = Pairs > C;     
        Worse = ( repmat( R./n, k, 1) - repmat( [R./n]', 1, k ) < 0 ).*...
            Sig;
        Better = ( repmat( R./n, k, 1) - repmat( [R./n]', 1, k ) > 0 ).*...
            Sig; 
        
        % Build L       
        Ranking = sum( Worse, 2);
        [KW,tmp] = sortrows( [ Ranking,[1:k]'] );
        Better = Better(tmp,:);
        Better = Better(:,tmp);
        L = [];
        for i = 1 : k
            L = [L, {find( Better(i,:) > 0 )}];
        end
        KW = [[1:k]',KW];
    else
        % There is no difference
        KW = [[1:k];zeros(1,k);1:k]';
        L = [];
        for i = 1 : k
            L = [L, {[]}];
        end
    end
    
    function R = rankWithTies( X, n )
        N = sum(n);
        R = [];
        for i = 1 : size(X,2)
            R = [ R; X(1:n(i),i), i*ones( n(i), 1) ];
        end
        R = sortrows(R,-1);
        R = [ [1:N]', R ];
        % Correct for ties
        tiefrom = 1;
        tieto = 1;
        for i = 2 : N
            if( R(tiefrom, 2 ) == R( i, 2) )
                tieto = i;                
            else          
                if( tieto > tiefrom )
                    R( tiefrom:tieto, 1 ) = mean( R( tiefrom:tieto, 1) );                    
                end
                tiefrom = i;
                tieto = i;
            end
        end
        if( tieto > tiefrom )
            R( tiefrom:tieto, 1 ) = mean( R( tiefrom:tieto, 1) );                    
        end
        