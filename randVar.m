

%% Random Variable class: Declares a random variable with unspecified type
%% Contains generalized methodsfor manipulating random variables.
%%
%% All variables assumed independent.
%%
%% Use not recommended if distribution is known to be normal, lognormal or
%% mixed normal-lognormal. For operation optimality use normVar, lognVar or
%% mixVar classes.
%%
%% Author: Karasmanoglou Apostolos @ apostolis.kar98@gmail.com
%%

classdef randVar
    properties
        pd %% Variable density function
        mu %% Mean
        var %% Variance
        dev %% Deviation
        mu3 %% 3rd central moment
        skew %% Skewness
        valspace %% P.D.F Domain of definition
        acc = 40 %% Plotting/Estimation evaluation sample number
    end
    methods
        
        %% randVar: Class Constructor method %%
        
        function rv = randVar(arg1, arg2)
            
            %% Define a new random variable or create a copy
            
            if(nargin == 1)
                
                %% Copy constructor functionality : Single argument is %%
                %% random variable that is to be copied                %%
                cpy_rv = arg1;
                
                %% Copy all random variable properties into new object %%
                rv.pd = cpy_rv.pd;
                rv.mu = cpy_rv.mu;
                rv.var = cpy_rv.var;
                rv.dev = cpy_rv.dev;
                rv.mu3 = cpy_rv.mu3;
                rv.skew = cpy_rv.skew;
                rv.valspace = cpy_rv.valspace;
                rv.acc = cpy_rv.acc;
                
            elseif(nargin == 2)
                
                %% Construct random variable from x->P(x) mapping. Every %%
                %% value x of xspace corresponds to a probability P(x)   %%
                %% in dist vector at the same index                      %%
                dist = arg2;
                xspace = arg1;
                
                %% Compute random variable properties %%
                rv.pd = dist;
                rv.mu = dist*xspace';
                rv.var = dist*((xspace' - rv.mu').^2);
                rv.dev = sqrt(rv.var');
                rv.mu3 = dist*(xspace' - rv.mu').^3;
                rv.skew = rv.mu3/rv.var;
                rv.valspace(1) = xspace(1);
                rv.valspace(2) = xspace(length(xspace));
                rv.acc = length(xspace);
            end
        end
        
        %% print: Prints all relevant information about a random  %%
        %% variable to screen.                                    %%
        
        function print(rv)
            
            %% prints random variable properties
            
            disp("Random Variable");
            fprintf("------------\n");
            fprintf("Type: Unspecified\n");
            fprintf("Mean: %.5f\n", rv.mu);
            fprintf("Variance: %.5f\n", rv.var);
            fprintf("Deviance: %.5f\n", rv.dev);
            fprintf("Skew: %.5f\n", rv.skew);
            fprintf("Valspace: " + mat2str(rv.valspace) + "\n");
            fprintf("Sample num : %d\n", rv.acc);
            fprintf("------------\n");
        end
        
        %% plotdist: Prints distribution of random variable to    %%
        %% figure #dest. Color/Print style information contained  %%
        %% in cstr                                                %%
        
        function plotdist(rv, dest, cstr)
            
            %% plots random variable distribution
            
            xplot = linspace(rv.valspace(1), rv.valspace(2), rv.acc); %% Generate pdf domain %%
            figure(dest); %%Open figure #dest %%
            
            if(rv.acc > 52)
                
                %% Large accuracy plots ploted continuously %%
                if(nargin == 3)
                    plot(xplot, rv.pd, cstr);
                else
                    plot(xplot, rv.pd);
                end
            else
                
                %% Small accuracy plots printed as bar graphs %%
                if(nargin == 3)
                    bar(xplot, rv.pd, cstr);
                else
                    bar(xplot, rv.pd);
                end
            end
        end
        
        %% plus(overload+): Overloading sumation function, adds     %%
        %% any two independent random variables convolving their    %%
        %% distributions to obtain the derived variable's p.d.f     %%
        %% Secondary functionality allows shifting a function by    %%
        %% addition of a double value representing the shift amount %%
        
        function sumVar = plus(arg1, arg2)
            
            %% Computes random variable sum of two random variable objects
            
            sumVar = randVar();
            
            if(isa(arg1, 'randVar') && isa(arg2, 'randVar'))
                
                %% Primary functionality. Convolve random variable %%
                %% distributions to obtain sum distribution        %%
                
                rv_left = arg1;
                rv_right = arg2;
                sumVar.pd = conv(rv_left.pd, rv_right.pd); %% Compute sum distribution %%
                
                %% Compute moments of distribution via addition %%
                sumVar.mu = rv_left.mu + rv_right.mu;
                sumVar.var = rv_left.var + rv_right.var;
                sumVar.dev = sqrt(sumVar.var);
                sumVar.mu3 = rv_left.mu3 + rv_right.mu3;
                sumVar.skew = sumVar.mu3/sumVar.var;
                
                sumVar.acc = rv_left.acc + rv_right.acc - 1; %% Increase amount of samples %%
                
                %% Calculate new domain for distribution             %%
                sumVar.valspace(1) = rv_left.valspace(1) + rv_right.valspace(1);
                sumVar.valspace(2) = rv_left.valspace(2) + rv_right.valspace(2);
            elseif(isa(arg1, 'randVar') && isa(arg2, 'double'))
                
                %% Secondary shift functionallity rv + shift %%
                sumVar = arg1;
                shift = arg2;
                
                sumVar.mu = sumVar.mu + shift; %% Shift distribution mean %%
                sumVar.valspace = sumVar.valspace + shift; %% Shift domain %%
            elseif(isa(arg2, 'randVar') && isa(arg1, 'double'))
                
                %% Secondary shift functionallity shift + rv %%
                sumVar  = arg2;
                shift = arg1;
                
                sumVar.mu = sumVar.mu + shift; %% Shift distribution mean %%
                sumVar.valspace = sumVar.valspace + shift; %% Shift domain %%
            end
        end
        
        %% minus(overload-): Overloading subtraction function,      %%
        %% (randVar, randVar) sums any two independent random       %%
        %% variables convolving their distributions to obtain the   %%
        %% derived variable's p.d.f. Secondary functionality allows %%
        %% shifting a function by addition of a double value        %%
        %% representing the shift amount.                           %%
        
        function diffVar = minus(arg1, arg2)
            
            %% Computes random variable difference of two random variable objects
            
            diffVar = randVar();
            
            if(isa(arg1, 'randVar') && isa(arg2, 'randVar'))
                
                %% Primary functionality. Convolve random variable %%
                %% distributions to obtain diff distribution       %%
                rv_left = arg1;
                rv_right = arg2;
                diffVar.pd = conv(rv_left.pd, rv_right.pd);
                
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
            elseif(isa(arg1, 'randVar') && isa(arv_right, 'double'))
                
                %% Secondary shift functionallity rv - shift %%
                rv = argv;
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
        
        %% similar: Calculate similarity between two variables  rv_left,rv_right %%
        %% in natural and log scale.                                    %%
        
        function [sim, logsim] = similar(rv_left, rv_right)
            
            %% Computes/plots similarity measure of two random variable distributions
            
            sim = rv_left.pd./rv_right.pd; %% Similarity function rv_left/rv_right %%
            logsim = log(sim); %% Log-Scale similarity function %%
            
            %% Plot similarity %%
            grid on;
            subplot(2, 1, 1);
            grid on;
            plot(sim, "-r");
            subplot(2, 1, 2);
            grid on;
            plot(logsim, "-b");
        end
        
        %% worst: Performs worst case analysis between two variables %%
        %% i.e selects that variable which is most likely to exhibit %%
        %% large values                                              %%
        function maxVar = worst(rv_left, rv_right)
            
            %% Selects random variable most likely to exhibit high values
            
            %% Calculate difference variable (convolution)%%
            diff = rv_right - rv_left;
            
            %% If values are entirely negative or entirely positive %%
            %% skip integral calculation                            %%
            if(diff.valspace(1) > 0)
                maxVar = rv_right;
                return;
            end
            if(diff.valspace(2) < 0)
                maxVar = rv_left;
                return;
            end
            
            %% Compute boundry of negative to positive samples space %%
            ii = 1;
            xval = linspace(diff.valspace(1), diff.valspace(2), diff.acc);
            while(xval(ii) < 0)
                ii = ii + 1;
            end
            
            %% Compute integral of difference distribution at negative & %%
            %% positive range                                            %%
            q1 = trapz(diff.pd(1:ii));
            q2 = 1 - q1;
            
            %% If positive/negative space integrals are approximately  %%
            %% equal decide max based on skew.                         %%
            if(abs(q1 - q2) < rv_left.tol)
                
                %% If skew's approximately equal decide distribution based %%
                %% on maximum variance                                     %%
                if(abs(rv_left.skew - rv_right.skew) < rv_left.tol)
                    if(rv_left.var < rv_right.var)
                        maxVar = rv_right.var;
                        return;
                    else
                        maxVar = rv_left.var;
                        return;
                    end
                end
                
                if(rv_left.skew < rv_right.skew)
                    maxVar = rv_right;
                    return;
                else
                    maxVar = rv_left;
                    return;
                end
            end
            
            %% Decide max from integral area at positive vs negative %%
            %% sample ranges.                                        %%
            if(q2 < q1)
                maxVar = rv_left;
            else
                maxVar = rv_right;
            end
        end
        
        %% gt(Overloading >): Returns logical value corresponding to the %%
        %% validity of the condition rv_left > rv_right in worst case terms  %%
        function isgt = gt(rv_left, rv_right)
            
            %% Returns true if right variable is more likely to exhibit larger values than the left
            
            %% Calculate difference variable (convolution)%%
            diff = rv_right - rv_left;
            
            %% If values are entirely negative or entirely positive %%
            %% skip integral calculation                            %%
            if(diff.valspace(1) > 0)
                isgt = false;
                return;
            end
            if(diff.valspace(2) < 0)
                isgt = true;
                return;
            end
            
            %% Compute boundry of negative to positive samples space %%
            ii = 1;
            xval = linspace(diff.valspace(1), diff.valspace(2), diff.acc);
            while(xval(ii) < 0)
                ii = ii + 1;
            end
            
            %% Compute integral of difference distribution at negative & %%
            %% positive range                                            %%
            q1 = trapz(diff.pd(1:ii));
            q2 = 1 - q1;
            
            %% If positive/negative space integrals are approximately  %%
            %% equal decide if condition holds based on skew.          %%
            if(abs(q1 - q2) < rv_left.tol)
                
                %% If skew's approximately equal decide distribution based %%
                %% on maximum variance                                     %%
                if(abs(rv_left.skew - rv_right.skew) < rv_left.tol)
                    if(rv_left.var < rv_right.var)
                        isgt = false;
                        return;
                    else
                        isgt = true;
                        return;
                    end
                end
                
                if(rv_left.skew < rv_right.skew)
                    isgt = false;
                    return;
                else
                    isgt = true;
                    return;
                end
            end
            
            %% Decide max from integral area at positive vs negative %%
            %% sample ranges.                                        %%
            if(q2 < q1)
                isgt = true;
            else
                isgt = false;
            end
        end
        
        function minVar = best(rv_left, rv_right)
            
            %% Selects random variable most likely to exhibit smaller values
            
            %% Calculate difference variable (convolution)%%
            diff = rv_right - rv_left;
            
            %% If values are entirely negative or entirely positive %%
            %% skip integral calculation                            %%
            if(diff.valspace(1) > 0)
                minVar = rv_left;
                return;
            end
            if(diff.valspace(2) < 0)
                minVar = rv_right;
                return;
            end
            
            %% Compute boundry of negative to positive samples space %%
            ii = 1;
            xval = linspace(diff.valspace(1), diff.valspace(2), diff.acc);
            while(xval(ii) < 0)
                ii = ii + 1;
            end
            
            %% Compute integral of difference distribution at negative & %%
            %% positive range                                            %%
            q1 = trapz(diff.pd(1:ii));
            q2 = 1 - q1;
            
            %% If positive/negative space integrals are approximately  %%
            %% equal decide min based on skew.                         %%
            if(abs(q1 - q2) < rv_left.tol)
                
                %% If skew's approximately equal decide distribution based %%
                %% on minimum variance                                     %%
                if(abs(rv_left.skew - rv_right.skew) < rv_left.tol)
                    if(rv_left.var < rv_right.var)
                        minVar = rv_left.var;
                        return;
                    else
                        minVar = rv_right.var;
                        return;
                    end
                end
                
                if(rv_left.skew < rv_right.skew)
                    minVar = rv_left;
                    return;
                else
                    minVar = rv_right;
                    return;
                end
            end
            
            %% Decide min from integral area at positive vs negative %%
            %% sample ranges.                                        %%
            if(q2 < q1)
                minVar = rv_right;
            else
                minVar = rv_left;
            end
        end
        
        function islt = lt(rv_left, rv_right)
            
            %% Returns true if right variable is more likely to exhibit larger values than the left
            
            %% Calculate difference variable (convolution)%%
            diff = rv_right - rv_left;
            
            %% If values are entirely negative or entirely positive %%
            %% skip integral calculation                            %%
            if(diff.valspace(1) > 0)
                islt = true;
                return;
            end
            if(diff.valspace(2) < 0)
                islt = false;
                return;
            end
            
            %% Compute boundry of negative to positive samples space %%
            ii = 1;
            xval = linspace(diff.valspace(1), diff.valspace(2), diff.acc);
            while(xval(ii) < 0)
                ii = ii + 1;
            end
            
            %% Compute integral of difference distribution at negative & %%
            %% positive range                                            %%
            q1 = trapz(diff.pd(1:ii));
            q2 = 1 - q1;
            
            %% If positive/negative space integrals are approximately  %%
            %% equal decide min based on skew.                         %%
            if(abs(q1 - q2) < rv_left.tol)
                
                %% If skew's approximately equal decide distribution based %%
                %% on minimum variance                                     %%
                if(abs(rv_left.skew - rv_right.skew) < rv_left.tol)
                    if(rv_left.var < rv_right.var)
                        islt = true;
                        return;
                    else
                        islt = false;
                        return;
                    end
                end
                
                if(rv_left.skew < rv_right.skew)
                    islt = true;
                    return;
                else
                    islt = false;
                    return;
                end
            end
            
            %% Decide min from integral area at positive vs negative %%
            %% sample ranges.                                        %%
            if(q2 < q1)
                islt = false;
            else
                islt = true;
            end
        end
        
        %% Computes new random variable Z = max{X,Y} attributes %%
        function rv = max(rv_left, rv_right)
            
            %% Computes max random variable distribution
            
            rv = randVar;
            
            %% Compute new random variable value space %%
            xspace = linspace(min(rv_left.valspace(1), rv_right.valspace(1)),max(rv_left.valspace(2), rv_right.valspace(2)), max(rv_left.acc, rv_right.acc));
            
            %% Compute and normalize pdf according to formula: %%
            %% fz = fx*Fy + Fx*fy                              %%
            rv.pd = cumsum(rv_left.pd).*rv_right.pd + cumsum(rv_right.pd).*rv_left.pd;
            rv.pd = rv.pd/sum(rv.pd); %% Normalize %%
            
            %% Compute attributes from x->P(x) relation %%
            rv.mu = rv.pd*xspace';
            rv.var = rv.pd*((xspace' - rv.mu').^2);
            rv.dev = sqrt(rv.var');
            rv.mu3 = rv.pd*(xspace' - rv.mu').^3;
            rv.skew = rv.mu3/rv.var;
            rv.valspace(1) = xspace(1);
            rv.valspace(2) = xspace(length(xspace));
            rv.acc = length(xspace);
        end
        
        %% Compute max{X,Y} P.D.F and approximate with mixed norm - %%
        %% lognorm distribution                                     %%
        
        function mv = mixmax(rv_left, rv_right)
            
            %% Computes mix approximation of max random variable distribution
            
            rv = max(rv_left, rv_right);
            mv = mixVar(rv);
        end
        
        %% Compute max{X,Y} P.D.F and approximate with normal  %%
        %% distribution                                        %%
        
        function nv = normmax(rv_left, rv_right)
            
            %% Computes norm approximation of max random variable
            
            rv = max(rv_left, rv_right);
            nv = normVar(rv);
        end
        
        %% Compute max{X,Y} P.D.F and approximate with lognnorm  %%
        %% distribution                                          %%
        
        function lnv = lognmax(rv_left, rv_right)
            
            %% Computes lognorm approximation of max random variable
            
            rv = max(rv_left, rv_right);
            lnv = lognVar(rv);
        end
        
        %% decimatepd : P.D.F Decimation function. Helps curb the amount  %%
        %% of samples to a given integer                                  %%
        
        function rv = decimatepd(rv, nsamples)
            
            %% Reduce sample amount of distribution
            
            resample_rate = ceil(rv.acc/nsamples); %% calculate decimation resample_rate %%
            rv.pd = decimate(rv.pd, resample_rate);
            rv.acc = length(rv.pd);
        end
        
        %% calcProb : Calculates the probability that the random variable %%
        %% rv has value >= left_bound and <= right_bound. Use -inf, inf   %%
        %% for lower/upper_bound when evaluating the probability of rv    %%
        %% being less than or greater than a certain value                %%
        
        function prob = calcProb(rv, left_bound, right_bound)
            
            %% Calculate the probability of left_bound <= rv <= right_bound
            
            h = (rv.valspace(2) - rv.valspace(1))/rv.acc;
            
            %% left and right side bounds may not exceed domain %%
            left = max(left_bound, rv.valspace(1));
            right = min(right_bound, rv.valspace(2));
            
            prob = 0;
            ii = 1;
            
            %% Calculate integral in a rectangular-rule way %%
            for x = left:h:right-h
                prob = prob + rv.pd(ii);
                ii = ii + 1;
            end
        end
    end
end