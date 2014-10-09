function [nhood, ss] = strel_neighbors( idx, s, SE )

ndim = numel( s );

% Get strucuring element offset from center

n = getneighbors( SE );
nse = size( n, 1 ); 
oset = zeros( nse, ndim );
oset( :,1:size(n, 2) ) = n;

% Convert index to subscripts

subs = cell( ndim, 1 );
[subs{:}] = ind2sub( s, idx );

% Add offsets

ss = repmat( cell2mat(subs)', nse, 1) + oset;

% Clip to boundaries

ss( ~all(bsxfun(@le, ss, s) & ss>0, 2),: ) = [];

% Convert to indeces: 
%              x        +       y    *s(1) + (   z-1   )*s(1)*s(2) + ...
% nhood = (ss(:,1)-0)*1 + (ss(:,2)-1)*s(1) + (ss(:,3)-1)*s(1)*s(2) + ...

ssm = bsxfun( @minus, ss, [0 ones(1,ndim-1)] );
nhood = sum( bsxfun(@times, ssm, [1 cumprod(s(1:end-1))]), 2 );



