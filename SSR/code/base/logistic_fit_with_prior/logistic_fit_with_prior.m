function [bhat,khat] = logistic_fit_with_prior( ti, yi, c, w, tip, yip, gammas, varargin )
% Fit logistic function with prior.
%
%  ti        - timepoints
%  yi        - measurement points
%  c         - upper asymptote for logistic model l = c/(1+b*exp(-k*t))
%  w         - weighting for the prior
%  tip       - prior time-points
%  yip       - measurements given as prior
%  gammas    - local strength of prior
%
% Optional arguments:
%  b0,k0     - initial conditions

b0 = 0;
k0 = 0;

%  Parse varargin
for k=1:2:length( varargin ) 
    eval( [ varargin{k}, '=varargin{', int2str(k+1), '};' ] );
end

% optimizer options fmincon
options_fmincon = optimset('fmincon');
options_fmincon = optimset( options_fmincon, 'Display', 'off' );
options_fmincon = optimset( options_fmincon, 'FunValCheck', 'off' );
%options_fmincon = optimset( options_fmincon, 'DerivativeCheck', 'on' );
options_fmincon = optimset( options_fmincon, 'GradObj', 'on' );
options_fmincon = optimset( options_fmincon, 'algorithm', 'interior-point' );
options_fmincon = optimset( options_fmincon, 'MaxFunEvals', 100 );

% constraints; b needs to be greater zero
x0 = [b0,k0];

lb = [-Inf; -Inf];
ub = [ Inf; Inf];

bkHat = fmincon( @(x)logistic_fit_with_prior_energy( x, ti, yi, c, w, tip, yip, gammas ), x0, [], [], [], [], lb, ub, [], options_fmincon );

bhat = bkHat(1);
khat = bkHat(2);

