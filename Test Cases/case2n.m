%% SSTA Matlab tool case 2 %%

%% SSTA Matlab tool size varying version of case2 %%
function rv = case2n(vb, n)

    %% Generate data batch %%
    xsample = normrnd(1, 0.1, 1, 1000);
    
    %% Fit distibution/ randVar on data 
    rv_basic = mixVar(xsample, [0 2], 25);
    
    %% Compute total delay of buffer path %% 
    rv = rv_basic;
    for ii = 2:n
        rv = rv + rv_basic;
        rv = decimatepd(rv, 25);
        
        %% Enforce mixed type %%
        if(~isa(rv, 'normVar') && ~isa(rv, 'lognVar') && ~isa(rv, 'mixVar'))
            rv = mixVar(rv);
        end
    end
    
    %% plot results %%
    rv = rv.configure_valspace([rv.mu - 3*rv.dev rv.mu + 3*rv.dev]);
    plotdist(rv, 1);
end