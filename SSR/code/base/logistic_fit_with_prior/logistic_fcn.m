function v = logistic_fcn( c, b, k, ti )
v = c./(1+b*exp(-k*ti) );

