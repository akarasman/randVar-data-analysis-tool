%% Mixed Random Variable class: Declares a random variable with a mixed
%% norm/lognorm distribution. Contains optimized methods for manipulating 
%% random variables.
%%
%% All variables assumed independent.                                
%%                                                                  
%% Author: Karasmanoglou Apostolos @ apostolis.kar98@gmail.com     
%%                                                                  

classdef mixVar < randVar
    properties
        a %% P.D.F mixing parameter 
        kshift %% lognorm part p.d.f shift value 
        npart %% normal p.d.f part 
        lnpart %% log-normal p.d.f part
    end
    
    methods
        
        %% mixVar: Class constructor can define a mixed random variable   %%
        %% from another random variable or a data set for fitting. Domain %%
        %% must be given                                                  %%
        
        function mv = mixVar(arg1, arg2, arg3)
            
            %% Construct new mixed normal/lognormal random variable or create a copy
            
            if(nargin == 1)
                
                %% Define a mixVar object from a randVar object %%
                if(isa(arg1, 'mixVar'))
                    
                    %% Copy all standard properties of mixVar variable %%
                    
                    cpy_rv = arg1; %% arg1 is the mixed variable to copy %%
                    mv.mu = cpy_rv.mu;
                    mv.var = cpy_rv.var;
                    mv.dev = cpy_rv.dev;
                    mv.valspace = cpy_rv.valspace;
                    mv.acc = cpy_rv.acc;
                    mv.pd = cpy_rv.pd;
                    mv.mu3 = cpy_rv.mu3;
                    mv.skew = cpy_rv.skew;
                    
                elseif(isa(arg1, 'randVar'))
                    
                    %% Copy all standard properties of randVar variable %%
                    %% and approximate p.d.f with mixed norm - lognorm  %%
                    %% using moment match method                        %%
                    
                    cpy_rv = arg1; %% arg1 is the random variable to copy %%
                    mv.mu = cpy_rv.mu;
                    mv.var = cpy_rv.var;
                    mv.dev = cpy_rv.dev;
                    mv.valspace = cpy_rv.valspace;
                    mv.acc = cpy_rv.acc;
                    mv.mu3 = cpy_rv.mu3;
                    mv.skew = cpy_rv.skew;
                    mv.kshift = mv.valspace(1);
                    if(mv.kshift > 0)
                        mv.kshift = 0;
                    end
                    
                    %% Perform transformation of mean/deviation for lognormal %%
                    lmu = log(((mv.mu - mv.kshift)^2)/sqrt(mv.var+(mv.mu - mv.kshift)^2));
                    lsigma = sqrt(log(mv.var/(mv.mu - mv.kshift)^2+1));
                    mv.a = min(abs(mv.mu3)/(((exp(lsigma^2) - 1)^2)*exp(3*lmu + 1.5*lsigma^2)*(exp(lsigma^2) + 2)),1);
                    
                    %% Compute normal part of estimation %%
                    mv.npart = normpdf(linspace(mv.valspace(1), mv.valspace(2), mv.acc), mv.mu, mv.dev);
                    
                    %% Compute lognormal part of estimation, flip tail if needed %%
                    if(mv.mu3 > 0)
                        mv.lnpart = lognpdf(linspace(mv.valspace(1), mv.valspace(2), mv.acc) - mv.kshift, lmu, lsigma);
                    else
                        mv.lnpart = lognpdf(2*mv.mu - linspace(mv.valspace(1), mv.valspace(2), mv.acc) - mv.kshift, lmu, lsigma);
                    end
                    
                    %% Mix normal and lognormal part to get estimated distribution %%
                    mv.pd = (1 - mv.a)*mv.npart + mv.a*mv.lnpart;
                    mv.pd = mv.pd/sum(mv.pd);
                    mv.npart = mv.npart/sum(mv.npart);
                    mv.lnpart = mv.lnpart/sum(mv.lnpart);
                end
            elseif(nargin == 3)
                
                %% Derive all properties from data set %%
                
                data = arg1; %% arg1 is the data set %%
                mv.acc = arg3; %% arg3 is the accuracy of the plot %%
                mv.valspace = arg2; %% arg2 is the domain  %%
                
                %% compute moments of distribution %%
                mv.mu = mean(data);
                mv.var = var(data);
                mv.dev = sqrt(mv.var);
                mv.mu3 = moment(data, 3);
                mv.skew = mv.mu3/mv.var;
                
                %% Use data fitting method %%
                [mv.pd, mv.kshift, mv.a, mv.npart, mv.lnpart] = mmfit(data, linspace(mv.valspace(1), mv.valspace(2), mv.acc));
                
                %% Compute normal/lognormal and total p.d.f %%
                mv.pd = mv.pd/sum(mv.pd);
                mv.npart = mv.npart/sum(mv.npart);
                mv.lnpart = mv.lnpart/sum(mv.lnpart);
            end
        end
        
        %% print: prints relevant information about random variable to screen %%
        
        function print(rv)
            disp("Random Variable");
            fprintf("------------\n");
            fprintf("Type: Mixed\n");
            fprintf("Mean: %.5f\n", rv.mu);
            fprintf("Variance: %.5f\n", rv.var);
            fprintf("Deviance: %.5f\n", rv.dev);
            fprintf("Skew: %.5f\n", rv.skew);
            fprintf("Valspace: " + mat2str(rv.valspace) + "\n");
            fprintf("Sample num : %d\n", rv.acc);
            fprintf("mix parameter(a): %d\n", rv.a);
            fprintf("kshift value: %d\n", rv.kshift);
            fprintf("------------\n");
        end
        
        %% plus(overload+): Overloading sumation function, adds    %%
        %% any two independent random variables by adding their    %%
        %% cumulants to obtain a new distribution via moment match %%
        
        function sumVar = plus(arg1, arg2)
            
            if(isa(arg1, 'randVar') && isa(arg2, 'randVar'))
                
                %% Primary addition functionality  %%
                rv_left = arg1;
                rv_right = arg2;
                sumVar = randVar();
                
                %% Compute moments of derived distribution %%
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
                if(isa(rv_left, 'mixVar') && isa(rv_right, 'mixVar'))
                    sumVar = mixVar(sumVar);
                    sumVar.kshift = rv_left.kshift +  rv_right.kshift;
                    xeval = linspace(sumVar.valspace(1), sumVar.valspace(2), sumVar.acc);
                    
                    %% Perform transformation of mean/deviation for lognormal %%
                    lmu = log((sumVar.mu^2)/sqrt(sumVar.var+sumVar.mu^2));
                    lsigma = sqrt(log(sumVar.var/(sumVar.mu^2)+1));
                    
                    
                    %% Compute mixing parameter %%
                    sumVar.a = abs(sumVar.mu3)/((exp(lsigma^2) - 1)*exp(2*lmu + lsigma^2)*(exp(lsigma^2) + 2));
                    
                    %% Compute normal part of estimation %%
                    sumVar.npart = normpdf(xeval, sumVar.mu, sumVar.dev);
                    
                    %% Compute lognormal part of estimation, flip tail if needed %%
                    if(sumVar.mu3 > 0)
                        sumVar.lnpart = lognpdf(xeval - sumVar.kshift, lmu, lsigma);
                    else
                        sumVar.lnpart = lognpdf(2*sumVar.mu - xeval - sumVar.kshift, lmu, lsigma);
                    end
                    
                    %% Mix normal and lognormal part to get estimated distribution %%
                    sumVar.pd = (1 - sumVar.a)*sumVar.npart + sumVar.a*sumVar.lnpart;
                    
                    sumVar.pd = sumVar.pd / sum(sumVar.pd);
                    sumVar.lnpart = sumVar.lnpart / sum(sumVar.lnpart);
                    sumVar.npart = sumVar.npart / sum(sumVar.npart);
                else
                    sumVar.pd = conv(rv_left.pd, rv_right.pd);
                end
            elseif(isa(arg1, 'randVar') && isa(arg2, 'double'))
                
                %% Secondary shift functionality rv + shift %%
                rv = arg1;
                shift = arg2;
                
                rv.mu = rv.mu + shift; %% Shift distribution mean %%
                rv.valspace = rv.valspace + rv; %% Shift domain %%
                sumVar = rv;
            elseif(isa(arg1, 'randVar') && isa(arg2, 'double'))
                
                %% Secondary shift functionality shift + rv %%
                rv = arg2;
                shift = arg1;
                
                rv.mu = rv.mu + shift; %% Shift distribution mean %%
                rv.valspace = rv.valspace + rv; %% Shift domain %%
                sumVar = rv;
            end
        end
        
        %% minus(overload-): Overloading minus function, subs any  %%
        %% two independent random variables by subtracting their   %%
        %% cumulants to obtain a new distribution via moment match %%
        %% Secondarilly shifts distribution by any double amount   %%
        
        function diffVar = minus(arg1, arg2)
            
            diffVar = randVar();
            
            if(isa(arg1, 'randVar') && isa(arg2, 'randVar'))
                
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
                if(isa(arg1, 'mixVar') && isa(arg2, 'mixVar'))
                    
                    %% Can create new mixVar without convolving distributions %%
                    diffVar = mixVar(diffVar);
                    xeval = linspace(diffVar.valspace(1), diffVar.valspace(2), diffVar.acc);
                    
                    %% Perform transformation of mean/deviation for lognormal %%
                    lmu = log((diffVar.mu^2)/sqrt(diffVar.var+diffVar.mu^2));
                    lsigma = sqrt(log(diffVar.var/(diffVar.mu^2)+1));
                    
                    %% Compute mixing parameter %%
                    diffVar.a = abs(diffVar.mu3)/((exp(lsigma^2) - 1)*exp(2*lmu + lsigma^2)*(exp(lsigma^2) + 2));
                    
                    %% Compute normal part of estimation %%
                    diffVar.npart = normpdf(xeval, diffVar.mu, diffVar.dev);
                    
                    %% Compute lognormal part of estimation, flip tail if needed %%
                    if(diffVar.mu3 > 0)
                        diffVar.lnpart = lognpdf(xeval, lmu, lsigma);
                    else
                        diffVar.lnpart = lognpdf(2*diffVar.mu - xeval, lmu, lsigma);
                    end
                    diffVar.pd = diffVar.a*diffVar.lnpart + (1 - diffVar.a)*diffVar.npart;
                    
                elseif(isa(arg1, 'randVar') && isa(arg2, 'double'))
                    
                    %% Secondary shift functionality rv - shift %%
                    rv = arg1;
                    shift = arg2;
                    
                    rv.mu = rv.mu - shift; %% Shift distribution mean %%
                    rv.valspace = rv.valspace - shift; %% Shift domain %%
                    diffVar = rv;
                elseif(isa(rv_right, 'randVar') && isa(rv_left, 'double'))
                    
                    %% Secondary shift functionality shift - rv %%
                    rv = arg2;
                    shift = arg1;
                    
                    rv.mu = rv.mu - shift; %% Shift distribution mean %%
                    rv.valspace = rv.valspace - shift; %% Shift domain %%
                    diffVar = rv;
                end
            end
        end
        
        %% decimatepd : P.D.F Decimation function. Helps curb the amount  %%
        %% of samples to a given integer                                  %%
        
        function mv = decimatepd(mv, nsamples)
            
            mv.acc = nsamples;
            
            %% Perform transformation of mean/deviation for lognormal %%
            lmu = log(((mv.mu - mv.kshift)^2)/sqrt(mv.var+(mv.mu - mv.kshift)^2));
            lsigma = sqrt(log(mv.var/(mv.mu - mv.kshift)^2+1));
            
            %% Compute normal part of estimation %%
            mv.npart = normpdf(linspace(mv.valspace(1), mv.valspace(2), mv.acc), mv.mu, mv.dev);
            
            %% Compute lognormal part of estimation, flip tail if needed %%
            if(mv.mu3 > 0)
                mv.lnpart = lognpdf(linspace(mv.valspace(1), mv.valspace(2), mv.acc) - mv.kshift, lmu, lsigma);
            else
                mv.lnpart = lognpdf(2*mv.mu - linspace(mv.valspace(1), mv.valspace(2), mv.acc) - mv.kshift, lmu, lsigma);
            end
            
            %% Mix normal and lognormal part to get estimated distribution %%
            mv.pd = (1 - mv.a)*mv.npart + mv.a*mv.lnpart;
            mv.pd = mv.pd/sum(mv.pd);
            mv.npart = mv.npart/sum(mv.npart);
            mv.lnpart = mv.lnpart/sum(mv.lnpart);
        end
        
        %% configure_valspace: Reset random variable domain %%
        
        function mv = configure_valspace(mv, new_valspace)
            
            %% Reset pdf domain
            
            mv.valspace = new_valspace;
            
            xeval = linspace(mv.valspace(1), mv.valspace(2), mv.acc);
            mv.npart = normpdf(xeval, mv.mu, mv.dev);
            mv.npart = mv.npart/sum(mv.npart);
            
            mv.kshift = mv.valspace(1);
            if(mv.kshift > 0)
                mv.kshift = 0;
            end
            
            %% Compute log scale mean and sigma %%
            lmu = log(((mv.mu - mv.kshift)^2)/sqrt(mv.var+(mv.mu - mv.kshift)^2));
            lsigma = sqrt(log(mv.var/(mv.mu - mv.kshift)^2+1));
            
            %% Compute p.d.f at given domain %%
            %% Compute lognormal part of estimation, flip tail if needed %%
            if(mv.mu3 > 0)
                mv.lnpart = lognpdf(linspace(mv.valspace(1), mv.valspace(2), mv.acc) - mv.kshift, lmu, lsigma);
                mv.lnpart = mv.lnpart/sum(mv.lnpart);
            else
                mv.lnpart = lognpdf(2*mv.mu - linspace(mv.valspace(1), mv.valspace(2), mv.acc) - mv.kshift, lmu, lsigma);
                mv.lnpart = mv.lnpart/sum(mv.lnpart);
            end
            
            mv.pd = mv.a*mv.lnpart + (1 - mv.a)*mv.npart;
            mv.pd = mv.pd / sum(mv.pd);
        end
    end
end