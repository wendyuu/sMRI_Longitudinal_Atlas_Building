function options = parseArgs( defaults, vararg )

options = {};

% let's see if we can find options
for ii=1:length(vararg)
    % let's see if we can find options, else storse parameters
    if strcmp(vararg{ii},'options')
        options = vararg{ii+1};
    end
end

% fill in the default values into options
if(isempty(options))
    options = defaults;
else
    % let's see which ones were given and fill in the defaults for the remaining ones
    structNames = fieldnames( defaults );
  
    for iI=1:length(structNames)
        if(~isfield( options, structNames{iI} ))
            eval( sprintf('options.%s = defaults.%s;', structNames{iI}, structNames{iI}) );
        end
    end
end

% now, update options with any additional varargin parameters
if( ~isempty(vararg) )
    for ii=1:2:length(vararg)
        if( ~strcmp(vararg{ii},'options') )
            eval( sprintf('options.%s = vararg{%i+1};', vararg{ii}, ii) );
        end
    end   
end