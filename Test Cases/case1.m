
%% SSTA Matlab tool test case 1 %%

function case1(vb)

    clf;
    clc;
    close all;
    
    %% 3 different batches of data %%
    xsample_norm1 = normrnd(5,2,1,1000);
    xsample_logn = lognrnd(2,0.1,1,1000);
    xsample_norm2 = normrnd(10,3,1,1000);
    
    %% Generate random variable object from data using varBox %%
    x1 = vb.makeVar(xsample_norm1, [0 30], 25);
    x2 = vb.makeVar(xsample_logn, [0 30], 25);
    x4 = vb.makeVar(xsample_norm2, [0 30], 25);
    
    %% plot distributions of two input buffer delays %%
    hold on;
    plotdist(x1, 1);
    plotdist(x2, 1);
    
    %% x3 is delay at input of 3rd buffer %%
    x3 = vb.max(x1,x2);
    plotdist(x3,1);
    legend('x1','x2','x3');
    
    %% x5 is total delay %%
    x5 = x3 + x4;
    x5_norm = normVar(x5); %% Convert x5 to normal var %%
    
    %% Plot total delay %%
    plotdist(x5, 2);
    hold on;
    plotdist(x5_norm, 2);
    legend('Total delay randVar', 'Normal pdf estimation');
    
    %% Print problem variable data  %%
    disp(x1);
    disp(x2);
    disp(x3);
    disp(x4);
    disp(x5);
    disp(x5_norm);
    prob = x5.calcProb(0,25)
    
end
