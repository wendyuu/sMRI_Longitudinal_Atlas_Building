function est = apply_logistic_model( img, model, mask, ti, varargin )

%  Parse varargin
for k=1:2:length( varargin ) 
    eval( [ varargin{k}, '=varargin{', int2str(k+1), '};' ] );
end

params = reshape( model, [], 2 );       

% Get mask indeces 
mask_idx = reshape( find(mask), 1, [] );

% Evaluate model
est = double( img );

for ii=mask_idx
    est(ii) = logistic_fcn( 1, params(ii,1), params(ii,2), ti );
end

% Scale
est(mask_idx) = est(mask_idx) .* ( ymax-ymin ) + ymin;
