%% Log-Normal Random Variable class: Declares a random variable following a  
%% lognorm distribution. Contains optimized methods for manipulating random
%% variables.  
%%
%% All variables assumed independent.                                  
%%                                                                  
%% Author: Karasmanoglou Apostolos @ apostolis.kar98@gmail.com      
%%                                                                  

classdef lognVar < randVar
    properties
        kshift %% lognorm p.d.f shift value 
        lmu %% log-domain mean 
        lsigma %% log-domain sigma 
    end
    
    methods
        
        %% lognVar: Class constructor can define a lognorm random variable %%
        %% from another random variable, mean and sigma^2 pair or a data   %%
        %% set for fitting.                                                %%
        function lnv = lognVar(arg1, arg2, arg3, arg4, arg5)
            
            %% Construct new lognormally distributed random variable or create copy
            
            if(nargin == 1)
                
                %% Define a mixVar object from a randVar object %%
                if(isa(arg1, 'mixVar'))
                    
                    %% Copy all standard properties %%
                    cpy_rv = arg1;
                    lnv.mu = cpy_rv.mu;
                    lnv.var = cpy_rv.var;
                    lnv.dev = cpy_rv.dev;
                    lnv.valspace = cpy_rv.valspace;
                    
                    %% Compute log-scale mean and variance %%
                    lnv.lmu = log((lnv.mu^2)/sqrt(lnv.var+lnv.mu^2));
                    lnv.lsigma = sqrt(log(lnv.var/(lnv.mu^2)+1));
                    
                    %% Calculate shift to positive domain values %%
                    lnv.kshift = cpy_rv.valspace(1);
                    if(lnv.kshift > 0)
                        lnv.kshift = 0;
                    end
                    lnv.acc = cpy_rv.acc;
                    lnv.pd = cpy_rv.lnpart;
                    lnv.pd = lnv.pd / sum(lnv.pd);
                    lnv.mu3 = cpy_rv.mu3;
                    lnv.skew = cpy_rv.skew;
                    
                elseif(isa(arg1, 'lognVar'))
                    
                    %% Copy all standard properties %%
                    cpy_rv = arg1;
                    lnv.mu = cpy_rv.mu;
                    lnv.var = cpy_rv.var;
                    lnv.dev = cpy_rv.dev;
                    lnv.lmu = cpy_rv.lmu;
                    lnv.lsigma = cpy_rv.lsigma;
                    lnv.valspace = cpy_rv.valspace;
                    lnv.kshift = cpy_rv.kshift;
                    lnv.acc = cpy_rv.acc;
                    lnv.pd = cpy_rv.pd;
                    lnv.mu3 = cpy_rv.mu3;
                    lnv.skew = cpy_rv.skew;
                elseif(isa(arg1, 'randVar'))
                    
                    %% Copy all standard properties %%
                    cpy_rv = arg1;
                    lnv.mu = cpy_rv.mu;
                    lnv.var = cpy_rv.var;
                    lnv.dev = cpy_rv.dev;
                    lnv.valspace = cpy_rv.valspace;
                    lnv.kshift = cpy_rv.valspace(1);
                    if(lnv.kshift > 0)
                        lnv.kshift = 0;
                    end
                    
                    %% Calculate log domain mean and sigma %%
                    lnv.lmu = log(((lnv.mu - lnv.kshift)^2)/sqrt(lnv.var+(lnv.mu - lnv.kshift)^2));
                    lnv.lsigma = sqrt(log(lnv.var/((lnv.mu - lnv.kshift)^2)+1));
                    
                    lnv.acc = cpy_rv.acc;
                    lnv.mu3 = cpy_rv.mu3;
                    lnv.skew = cpy_rv.skew;
                    
                    %% Moment match to approximate p.d.f %%
                    if(lnv.mu3 >= 0)
                        xeval = linspace(lnv.valspace(1) - lnv.kshift, lnv.valspace(2) - lnv.kshift, lnv.acc);
                        lnv.pd = lognpdf(xeval, lnv.lmu, lnv.lsigma);
                        lnv.pd = lnv.pd/sum(lnv.pd);
                    else
                        xeval = 2*lnv.mu - lnv.kshift - linspace(lnv.valspace(1), lnv.valspace(2), lnv.acc);
                        lnv.pd = lognpdf(xeval, lnv.lmu, lnv.lsigma);
                        lnv.pd = lnv.pd/sum(lnv.pd);
                    end
                end
            elseif(nargin == 3)
                
                %% Define a Var from a data set and a plot accuracy %%
                
                data = arg1; %% arg1 is the data set %%
                lnv.acc = arg2; %% arg2 is the plot accuracy %%
                lnv.valspace = arg3; %% arg3 is the domain %%
                lnv.mu = mean(data);
                lnv.var = moment(data,2);
                lnv.dev = sqrt(lnv.var);
                lnv.mu3 = moment(data,3);
                lnv.kshift = lnv.valspace(1);
                if(lnv.kshift > 0)
                    lnv.kshift = 0;
                end
                
                %% Calculate log mean and sigma %%
                lnv.lmu = log(((lnv.mu - lnv.kshift)^2)/sqrt(lnv.var+(lnv.mu - lnv.kshift)^2));
                lnv.lsigma = sqrt(log(lnv.var/((lnv.mu - lnv.kshift)^2)+1));
                
                if(lnv.mu3 >= 0)
                    lnv.pd = lognpdf(linspace(lnv.valspace(1) - lnv.kshift, lnv.valspace(2) - lnv.kshift, lnv.acc), lnv.lmu, lnv.lsigma);
                    lnv.pd = lnv.pd/sum(lnv.pd);
                else
                    lnv.pd = lognpdf(2*lnv.mu - lnv.kshift - linspace(lnv.valspace(1), lnv.valspace(2), lnv.acc), lnv.lmu, lnv.lsigma);
                    lnv.pd = lnv.pd/sum(lnv.pd);
                end
                lnv.skew = lnv.mu3/lnv.var;
            elseif(nargin == 5)
                
                %% Define a lognVar variable from a mean a sigma^2 pair and a %%
                %% plotting accuracy                                          %%
                
                lnv.mu = arg1; %% arg1 is the mean %%
                lnv.var = arg2; %% arg2 is the variance %%
                lnv.valspace = arg3; %% arg4 is the domain %%
                lnv.acc = arg4; %% arg3 is the plot accuracy %%
                sign = arg5; %% arg5 is the skewness sign %%
                lnv.dev = sqrt(lnv.var);
                lnv.kshift = lnv.valspace(1);
                if(lnv.kshift > 0)
                    lnv.kshift = 0;
                end
                
                %% Calculate log mean and sigma %%
                lnv.lmu = log(((lnv.mu - lnv.kshift)^2)/sqrt(lnv.var+(lnv.mu - lnv.kshift)^2));
                lnv.lsigma = sqrt(log(lnv.var/((lnv.mu - lnv.kshift)^2)+1));
                
                %% Calculate 3rd central moment and skewness %%
                lnv.mu3 = sign*(((exp(lnv.lsigma^2) - 1)^2)*exp(3*lnv.lmu + 1.5*lnv.lsigma^2)*(exp(lnv.lsigma^2) + 2));
                lnv.skew = lnv.mu3/lnv.var;
                if(lnv.mu3 > 0)
                    xeval = linspace(lnv.valspace(1), lnv.valspace(2), lnv.acc) - lnv.kshift;
                    lnv.pd = lognpdf(xeval, lnv.lmu, lnv.lsigma);
                else
                    xeval = 2*lnv.mu - linspace(lnv.valspace(1), lnv.valspace(2));
                    lnv.pd = lognpdf(xeval, lnv.lmu, lnv.lsigma);
                end
                lnv.pd = lnv.pd/sum(lnv.pd);
            end
        end
        
        
        %% plus(overload+): Overloading sumation function, adds    %%
        %% any two independent random variables by adding their    %%
        %% cumulants to obtain a new distribution via moment match %%
        function sumVar = plus(rv1, rv2)
            if(isa(rv1, 'randVar') && isa(rv2, 'randVar'))
                
                sumVar = randVar();
                
                %% Compute moments of distribution via addition and   %%
                %% efficiently                                        %%
                sumVar.mu = rv1.mu + rv2.mu;
                sumVar.var = rv1.var + rv2.var;
                sumVar.dev = sqrt(sumVar.var);
                sumVar.mu3 = rv1.mu3 + rv2.mu3;
                sumVar.skew = sumVar.mu3/sumVar.var;
                
                sumVar.acc = rv1.acc + rv2.acc - 1; %% Increase amount of samples %%
                
                %% Calculate new domain for distribution             %%
                sumVar.valspace(1) = sumVar.mu - 6*sumVar.dev;
                sumVar.valspace(2) = sumVar.mu + 6*sumVar.dev;
                
                %% Compute sum distribution %%
                if(isa(rv1, 'lognVar') && isa(rv2, 'lognVar'))
                    sumVar = lognVar(sumVar);
                    sumVar.kshift = rv1.kshift + rv2.kshift;
                    size(sumVar.kshift)
                    if(sumVar.mu3 > 0)
                        sumVar.pd = lognpdf(linspace(sumVar.valspace(1), sumVar.valspace(2), sumVar.acc) - sumVar.kshift, sumVar.lmu, sumVar.lsigma);
                    else
                        sumVar.pd = lognpdf(2*sumVar.mu - linspace(sumVar.valspace(2), sumVar.valspace(1), sumVar.acc) - sumVar.kshift, sumVar.lmu, sumVar.lsigma);
                    end
                else
                    sumVar.pd = conv(rv1.pd, rv2.pd);
                end
            elseif(isa(rv1, 'randVar') && isa(rv2, 'double'))
                rv1.mu = rv1.mu + rv2; %% Shift distribution mean %%
                rv1.valspace = rv1.valspace + rv2; %% Shift domain %%
                sumVar = rv1;
            elseif(isa(rv2, 'randVar') && isa(rv1, 'double'))
                rv2.mu = rv2.mu + rv1; %% Shift distribution mean %%
                rv2.valspace = rv2.valspace + rv1; %% Shift domain %%
                sumVar = rv2;
            end
        end
        
        %% minus(overload+): Overloading minus function, adds any  %%
        %% two independent random variables by subtracting their   %%
        %% cumulants to obtain a new distribution via moment match %%
        function diffVar = minus(rv1, rv2)
            
            diffVar = randVar();
            
            if(isa(rv1, 'randVar') && isa(rv2, 'randVar'))
                %% Compute moments of distribution via addition and   %%
                %% subtraction efficiently                            %%
                diffVar.mu = rv1.mu - rv2.mu;
                diffVar.var = rv1.var + rv2.var;
                diffVar.dev = sqrt(diffVar.var);
                diffVar.mu3 = rv1.mu3 - rv2.mu3;
                diffVar.skew = diffVar.mu3/diffVar.var;
                diffVar.acc = rv1.acc + rv2.acc - 1; %% Increase amount of samples %%
                
                %% Calculate new domain for distribution             %%
                diffVar.valspace(1) = diffVar.mu - 6*diffVar.dev;
                diffVar.valspace(2) = diffVar.mu + 6*diffVar.dev;
                
                %% Compute difference distribution %%
                if(isa(rv1, 'lognVar') && isa(rv2, 'lognVar'))
                    diffVar = lognVar(diffVar);
                    diffVar.kshift = diffVar.valspace(1);
                    if(diffVar.kshift > 0)
                        diffVar.kshift = 0;
                    end
                    
                    if(diffVar.mu3 > 0)
                        diffVar.pd = lognpdf(linspace(diffVar.valspace(1), diffVar.valspace(2), diffVar.acc) - diffVar.kshift, diffVar.lmu, diffVar.lsigma);
                    else
                        diffVar.pd = lognpdf(2*diffVar.mu - linspace(diffVar.valspace(2), diffVar.valspace(1), diffVar.acc) - diffVar.kshift, diffVar.lmu, diffVar.lsigma);
                    end
                else
                    diffVar.pd = conv(rv1.pd, rv2.pd);
                end
            elseif(isa(rv1, 'randVar') && isa(rv2, 'double'))
                rv1.mu = rv1.mu - rv2; %% Shift distribution mean %%
                rv1.valspace = rv1.valspace - rv2; %% Shift domain %%
                diffVar = rv1;
            elseif(isa(rv2, 'randVar') && isa(rv1, 'double'))
                rv2.mu = rv2.mu - rv1; %% Shift distribution mean %%
                rv2.valspace = rv2.valspace - rv1; %% Shift domain %%
                diffVar = rv2;
            end
        end
        
        %% configure_valspace: Reset random variable domain %%
        
        function lnv = configure_valspace(lnv, new_valspace)
            
            %% Reset pdf domain
            
            lnv.valspace = new_valspace;
            lnv.kshift = lnv.valspace(1);
            if(lnv.kshift > 0)
                lnv.kshift = 0;
            end
            
            %% Compute log scale mean and sigma %%
            lnv.lmu = log(((lnv.mu - lnv.kshift)^2)/sqrt(lnv.var+(lnv.mu - lnv.kshift)^2));
            lnv.lsigma = sqrt(log(lnv.var/(lnv.mu - lnv.kshift)^2+1));
            
            %% Compute p.d.f at given domain %%
            if(lnv.mu3 >= 0)
                lnv.pd = lognpdf(linspace(lnv.valspace(1) - lnv.kshift, lnv.valspace(2) - lnv.kshift, lnv.acc), lnv.lmu, lnv.lsigma);
                lnv.pd = lnv.pd/sum(lnv.pd);
            else
                lnv.pd = lognpdf(2*lnv.mu - lnv.kshift - linspace(lnv.valspace(1), lnv.valspace(2), lnv.acc), lnv.lmu, lnv.lsigma);
                lnv.pd = lnv.pd/sum(lnv.pd);
            end
        end
        
        %% decimatepd: Reduce amount of samples of p.d.f %%
        
        function lnv = decimatepd(lnv, nsample)
            
            lnv.acc = nsample;
            if(lnv.mu3 >= 0)
                xeval = linspace(lnv.valspace(1) - lnv.kshift, lnv.valspace(2) - lnv.kshift, lnv.acc);
                lnv.pd = lognpdf(xeval, lnv.lmu, lnv.lsigma);
                lnv.pd = lnv.pd/sum(lnv.pd);
            else
                xeval = 2*lnv.mu - lnv.kshift - linspace(lnv.valspace(1), lnv.valspace(2), lnv.acc);
                lnv.pd = lognpdf(xeval, lnv.lmu, lnv.lsigma);
                lnv.pd = lnv.pd/sum(lnv.pd);
            end
        end
    end
end