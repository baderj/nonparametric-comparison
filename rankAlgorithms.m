%RANKALGORITHMS Rank the algorithm according to statistical significance
%   RANKALGORITHMS( H, N, ALPHA, FILE ) ranks the algorithms N according
%   to the number of competitors they are significantly better to at
%   the level ALPHA considering the indicator values H (Bigger values
%   lead to better, lower ranks). If FILE is specified, the results are 
%   written to FILE, otherwise to standard out.
%
%   H is a matrix of size NxK, where N denotes the number of runs and
%       K the number of algorithms. If the number of runs varies, then
%       the columns need to be stuffed with NaN.
%   N is a cellvector containing the names of the algorithms
%   ALPHA is a significance level, default 0.05
%   FILE is a path to the output file, default stdout
%
%   The output consists of a four tab seperated columns, where the
%   second entry denotes the number of algorithms dominated by A
%   where A is given by the third entry. The fourth column gives a list 
%   of algorithms dominated by A, refering to column 1
%
%   The statistical test is described in KRUSKALWALLISRANKING.
%
%   See also KRUSKALWALLISRANKING


function rankAlgorithms( H, N, alpha, file )
    Hmin = min( min( H ) );
    Hmax = max( max( H ) );
    Hspread = (Hmax - Hmin);

    if( nargin < 4 )
        fid = -1;
    else
        fid = fopen( file, 'w' );
    end
    if( nargin < 3 )
        alpha = 0.05;
    end
    assert( size(H,1) > 1 );
    assert( size(H,2) > 1 );
    [KW, L] = kruskalWallisRanking( H, alpha );
    
        
    for i = 1 : size(H,2)
        if( fid == -1 )           
            fprintf( '%d\t%d\t%s\t( %1.10f )\t',KW(i,1), KW(i,2), cell2mat( N(KW(i,3)) ), (mean( H(:,KW(i,3))) - Hmin) / Hspread );
            B = cell2mat( L(i) );
            for j = 1 : length(B)
                fprintf('%d', B(j) );
                if( j < length(B) )
                    fprintf(', ');
                end
            end
            fprintf('\n');
        else
            fprintf( fid,'%d\t%d\t%s\t( %1.10f )\t',KW(i,1), KW(i,2), cell2mat( N(KW(i,3)) ), (mean( H(:,KW(i,3))) - Hmin) / Hspread );
            B = cell2mat( L(i) );
            for j = 1 : length(B)
                fprintf(fid,'%d', B(j) );
                if( j < length(B) )
                    fprintf(fid,', ');
                end
            end
            fprintf(fid,'\n');
        end
    end
    
    if( fid ~= -1 )
        fclose(fid);
    end