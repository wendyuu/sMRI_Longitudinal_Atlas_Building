function im = imgresize( im, omega, m, varargin )

if( isempty(im) ), return; end

for k=1:2:length(varargin),
  eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

dim = numel(m);
sz  = size( im );

if( ~exist('method','var') )
    method = 'linear';
end

if( ~isequal(size(im), m) )
    switch( method )
        case {'nearest', 'linear', 'spline', 'cubic'}
            % MATLAB
            step = sz./m;
            
            switch( dim )
                case 2
                    [xi,yi] = meshgrid(1:step(2):sz(2),1:step(1):sz(1));
                    im = interp2( im, xi, yi, method );                    
                case 3
                    [xi,yi,zi] = meshgrid(1:step(2):sz(2),1:step(1):sz(1),1:step(3):sz(3));
                    im = interp3( im, xi, yi, zi, method );
            end
        otherwise
            % FAIR
            inter( 'reset', 'inter', method );
            xc = getCenteredGrid( omega, m );
            im = reshape( inter(im, omega, xc), m );
    end
end