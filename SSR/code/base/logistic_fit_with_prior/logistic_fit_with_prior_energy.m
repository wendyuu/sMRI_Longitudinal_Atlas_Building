function [E,grad] = logistic_fit_with_prior_energy( p, ti, yi, c, w, tip, yip, gammas )
% Logistic function fit energy.
%
% p       - function parameters
% ti      - data time-points
% yi      - data-points
% c       - upper asymptote
% w       - weighting for the prior
% tip     - prior time-points
% yip     - prior data-points
% gammas  - local weighting for the prior

b = p(1);
k = p(2);

vi = logistic_fcn( c, b, k, ti );
vip = logistic_fcn( c, b, k, tip ); % sampled at prior time-points

gammas = w*gammas/numel(gammas); % normalize gammas, so prior weighting doesn't depend on the number of sample points

E = double( 0.5*(sum( (yi-vi).^2 )  +  sum(gammas.*( yip-vip ).^2 )) );

if ( nargout>1 )
  grad = compute_logistic_fit_with_prior_gradient( c, b, k, ti, yi, vi, tip, yip, vip, gammas );
end

function grad = compute_logistic_fit_with_prior_gradient( c, b, k, ti, yi, vi, tip, yip, vip, gammas )

% gradB = sum( c*   ( (yi-vi) + gammas.*(yip-vi) ).*exp( -k*ti )    ./((1+b*exp(-k*ti)).^2) );
% gradK = sum( -c*b*( (yi-vi) + gammas.*(yip-vi) ).*exp( -k*ti ).*ti./((1+b*exp(-k*ti)).^2) );

gradB = sum(  c*  (yi-vi).*exp( -k*ti )     ./((1+b*exp(-k*ti)).^2) )   +   sum(  c*  gammas.*(yip-vip).*exp( -k*tip )     ./((1+b*exp(-k*tip)).^2) );
gradK = sum( -c*b*(yi-vi).*exp( -k*ti ).*ti ./((1+b*exp(-k*ti)).^2) )   +   sum( -c*b*gammas.*(yip-vip).*exp( -k*tip ).*tip./((1+b*exp(-k*tip)).^2) );

grad = double([gradB;gradK]);
