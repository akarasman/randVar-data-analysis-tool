%% Normal Random Variable class: Declares a random variable with a normal
%% distribution.Contains optimized methods for manipulating random
%% variables.
%%
%% All variables assumed independent.
%%
%% Author: Karasmanoglou Apostolos @ apostolis.kar98@gmail.com
%%

classdef normVar < randVar
    methods
        
        %% normVar: Class constructor, can define a normal random   %%
        %% variable from another random variable, a mean and a      %%
        %% variance or a data set for fitting. Domain must be given %%
        
        function nv = normVar(arg1, arg2, arg3, arg4)
            
            %% Construct new normaly distributed random variable or create a copy
            
            if(nargin == 1)
                
                %% Define normVar from a randVar object %%
                if(isa(arg1, 'randVar'))
                    
                    %% Copy properties %%
                    cpy_rv = arg1; %% arg1 is the random variable %%
                    nv.mu = cpy_rv.mu;
                    nv.var = cpy_rv.var;
                    nv.dev = cpy_rv.dev;
                    nv.valspace = cpy_rv.valspace;
                    nv.acc = cpy_rv.acc;
                    nv.mu3 = cpy_rv.mu3;
                    nv.skew = cpy_rv.skew;
                    
                    %% Calculate p.d.f for conversion from mixVar/randVar %%
                    %% to normVar                                         %%
                    if(isa(cpy_rv, 'mixVar'))
                        nv.pd = cpy_rv.npart;
                    else
                        xeval = linspace(nv.valspace(1), nv.valspace(2), nv.acc);
                        nv.pd = normpdf(xeval, nv.mu, nv.dev); %% Moment match approx %%
                        nv.pd = nv.pd/sum(nv.pd);
                    end
                    
                end
            elseif(nargin == 3)
                
                if(isa(arg1, 'double') && isa(arg2, 'double') && (length(arg1) == 1))
                    
                    %% Define a normVar from a mu and sigma^2 pair %%
                    
                    %% Copy all properties %%
                    nv.mu = arg1; %% arg1 is the mean %%
                    nv.var = arg2; %% arg2 is the variance %%
                    nv.valspace = arg3; %% arg3 is the domain %%
                    nv.dev = sqrt(nv.var);
                    xeval = linspace(nv.valspace(1), nv.valspace(2), nv.acc);
                    nv.pd = normpdf(xeval, nv.mu, nv.var);
                    
                    nv.skew = 0;
                else
                    
                    %% Define a normVar from a data set and a plot accuracy %%
                    
                    data = arg1; %% arg1 is the data set %%
                    acc = arg2; %% arg2 is the plot accuracy %%
                    nv.valspace = arg3; %% arg3 is the domain %%
                    nv.mu = mean(data);
                    nv.var = var(data);
                    nv.mu3 = 0;
                    nv.skew = 0;
                    nv.dev = sqrt(nv.var);
                    nv.acc = acc;
                    xeval = linspace(nv.valspace(1), nv.valspace(2), nv.acc);
                    nv.pd = normpdf(xeval, nv.mu, nv.var);
                end
            elseif(nargin == 4)
                
                %% Define a normVar variable from a mean a sigma^2 pair and a %%
                %% plotting accuracy && domain                                %%
                
                nv.mu = arg1; %% arg1 is the mean %%
                nv.var = arg2; %% arg2 is the variance %%
                nv.valspace = arg3; %% arg4 is the domain %%
                nv.acc = arg4; %% arg3 is the plot accuracy %%
                nv.dev = sqrt(nv.var);
                xeval = linspace(nv.valspace(1), nv.valspace(2), nv.acc);
                nv.pd = normpdf(xeval, nv.mu, nv.var);
                nv.pd = nv.pd/sum(nv.pd);
                nv.mu3 = 0;
                nv.skew = 0;
            end
        end
        
        %% print: prints relevant information about random variable to screen %%
        
        function print(rv)
            fprintf("Random Variable");
            fprintf("------------\n");
            fprintf("Type: Normal\n");
            fprintf("Mean: %.5f\n", rv.mu);
            fprintf("Variance: %.5f\n", rv.var);
            fprintf("Deviance: %.5f\n", rv.dev);
            fprintf("Valspace: " + mat2str(rv.valspace) + "\n");
            fprintf("Sample num : %d\n", rv.acc);
            fprintf("------------\n");
        end
        
        %% plus(overload+): Overloading sumation function, adds    %%
        %% any two independent random variables by adding their    %%
        %% cumulants to obtain a new normal distribution           %%
        
        function sumVar = plus(arg1, arg2)
            
            if(isa(arg1, 'randVar') && isa(arg2, 'randVar'))
                
                %% Primary addition functionality %%
                rv_left = arg1;
                rv_right = arg2;
                sumVar = randVar();
                
                
                %% Compute moments of distribution via addition and   %%
                %% efficiently                                        %%
                sumVar.mu = rv_left.mu + rv_right.mu;
                sumVar.var = rv_left.var + rv_right.var;
                sumVar.dev = sqrt(sumVar.var);
                sumVar.mu3 = rv_left.mu3 + rv_right.mu3;
                sumVar.skew = sumVar.mu3/sumVar.var;
                
                sumVar.acc = rv_left.acc + rv_right.acc - 1; %% Increase amount of samples %%
                
                %% Calculate new domain for distribution             %%
                sumVar.valspace(1) = rv_left.valspace(1) + rv_right.valspace(1);
                sumVar.valspace(2) = rv_left.valspace(2) + rv_right.valspace(2);
                
                %% Compute sum distribution %%
                if(isa(rv_left, 'normVar') && isa(rv_right, 'normVar'))
                    
                    %% If both variables are normal convolution unnecessary %%
                    sumVar = normVar(sumVar);
                    xeval = linspace(sumVar.valspace(1), sumVar.valspace(2), sumVar.acc);
                    sumVar.pd = normpdf(xeval, sumVar.mu, sumVar.dev);
                else
                    sumVar.pd = conv(rv_left.pd, rv_right.pd);
                end
            elseif(isa(arg1, 'randVar') && isa(arg2, 'double'))
                
                %% Secondary shift functionallity rv + shift %%
                rv = arg1;
                shift = arg2;
                
                rv.mu = rv.mu + shift; %% Shift distribution mean %%
                rv.valspace = rv.valspace + shift; %% Shift domain %%
                sumVar = rv;
            elseif(isa(arg2, 'randVar') && isa(arg1, 'double'))
                
                %% Secondary shift functionallity shift + rv %%
                rv = arg2;
                shift = arg1;
                
                rv.mu = rv.mu + shift; %% Shift distribution mean %%
                rv.valspace = rv.valspace + shift; %% Shift domain %%
                sumVar = rv;
            end
        end
        
        %% minus(overload-): Overloading minus function, subtracts any %%
        %% two independent random variables by subtracting their       %%
        %% cumulants to obtain a new normal variable                   %%
        
        function diffVar = minus(arg1, arg2)
            
            diffVar = randVar();
            
            if(isa(arg1, 'randVar') && isa(arg2, 'randVar'))
                
                %% Primary subtraction functionality %%
                rv_left = arg1;
                rv_right = arg2;
                
                %% Compute moments of distribution via addition and   %%
                %% subtraction efficiently                            %%
                diffVar.mu = rv_left.mu - rv_right.mu;
                diffVar.var = rv_left.var + rv_right.var;
                diffVar.dev = sqrt(diffVar.var);
                diffVar.mu3 = rv_left.mu3 - rv_right.mu3;
                diffVar.skew = diffVar.mu3/diffVar.var;
                
                diffVar.acc = rv_left.acc + rv_right.acc - 1; %% Increase amount of samples %%
                
                %% Calculate new domain for distribution             %%
                diffVar.valspace(1) = rv_left.valspace(1) - rv_right.valspace(2);
                diffVar.valspace(2) = rv_left.valspace(2) - rv_right.valspace(1);
                
                %% Compute sum distribution %%
                if(isa(rv_left, 'normVar') && isa(rv_right, 'normVar'))
                    diffVar = normVar(diffVar);
                    xeval = linspace(diffVar.valspace(1), diffVar.valspace(2), diffVar.acc);
                    diffVar.pd = normpdf(xeval, diffVar.mu, diffVar.dev);
                else
                    diffVar.pd = conv(rv_left.pd, rv_right.pd);
                end
            elseif(isa(arg1, 'randVar') && isa(arg2, 'double'))
                
                %% Secondary shift functionallity rv - shift %%
                rv = arg1;
                shift = arg2;
                
                rv.mu = rv.mu - shift; %% Shift distribution mean %%
                rv.valspace = rv.valspace - shift; %% Shift domain %%
                diffVar = rv;
            elseif(isa(arg1, 'randVar') && isa(arg2, 'double'))
                
                %% Secondary shift functionallity shift - rv %%
                rv = arg1;
                shift = arg2;
                
                rv.mu = rv.mu - shift; %% Shift distribution mean %%
                rv.valspace = rv.valspace - shift; %% Shift domain %%
                diffVar = rv;
            end
        end
        
        %% configure_valspace: Reset random variable domain %%
        
        function nv = configure_valspace(nv, new_valspace)
            
            %% Change function domain
            
            nv.valspace = new_valspace;
            xeval = linspace(nv.valspace(1), nv.valspace(2), nv.acc);
            nv.pd = normpdf(xeval, nv.mu, nv.var);
        end
        
        %% decimatepd: reduce amount of samples of pdf %%
        
        function nv = decimatepd(nv, nsample)
            nv.acc = nsample;
            nv.pd = normpdf(linspace(nv.valspace(1), nv.valspace(2), nv.acc), nv.mu, nv.dev);
            nv.pd = nv.pd/sum(nv.pd);
        end
    end
    
end