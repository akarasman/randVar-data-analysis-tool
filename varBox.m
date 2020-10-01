%% varBox class: Used for generating random variables from data set varBox 
%% class implements a black box enclosing all randVar class functionalities. 
%% Performs transformation of random variables from one type to the other 
%% based on tolerance properties, also performs given type of analysis at 
%% parallel buffer-path joins in circuit which allows us to obtain different
%% interpratations of data.

classdef varBox
    properties
        log_tol = 5*10^-3 %% Limit for "a" parameter to transdorm mixed variable to lognorm 
        norm_tol = 5*10^-3 %% Limit for "a" parameter to transdorm mixed variable to normal 
        analysis_type = 'worst' %% Analysis type at path joints 
        sample_limit
    end
    methods
        
        %% varBox : class constructor, define an object by lognorm, norm %%
        %% tolerance and analysis type                                   %%
        
        function vb = varBox(arg1, arg2, arg3, arg4)
            
            %% Constructs a varBox set to given parameters
            
            vb.log_tol = arg1;
            vb.norm_tol = arg2;
            vb.analysis_type = arg3;
            vb.sample_limit = arg4;
        end
        
        %% makeVar: Obtain a random variable from a given dataset. varBox %%
        %% decides type of randVar object based on tolerance properties   %%
        
        function rv = makeVar(vb, data, valspace, acc)
            
            %% Obtain a random variable from a given dataset
            
            %% Initial "guess" is mixed %%
            rv = mixVar(data, valspace, min(acc, vb.sample_limit));
            
            %% Decide type of variable based on varBox tolerances %%
            if(1 - rv.a < vb.log_tol)
               rv = lognVar(rv);
            elseif(rv.a < vb.norm_tol)
               rv = normVar(rv);
            end
        end
        
        %% max: max analysis at buffer path joints. Multiplexes different %%
        %% functionallities for different implementations of analysis     %%
        
        function rv = max(vb, rv1, rv2)
            
            %% Computes max between two variables based on set type of analysis
            
            if(strcmp(vb.analysis_type, 'worst'))
                rv = worst(rv1,rv2); %% worst case analysis %%
            elseif(strcmp(vb.analysis_type, 'best'))
                rv = best(rv1,rv2); %% best case analysis %%
            elseif(strcmp(vb.analysis_type, 'lognmax'))
                rv = lognmax(rv1, rv2); %% max and logn fit analysis %%
            elseif(strcmp(vb.analysis_type, 'mixmax'))
                rv = mixmax(rv1, rv2); %% max and mix fit analysis %%
            elseif(strcmp(vb.analysis_type, 'normmax'))
                rv = normmax(rv1,rv2); %% max and norm fit analysis %%
            else
                rv = worst(rv1,rv2); %% worst case analysis as default %%
            end
        end
        
    end     
end