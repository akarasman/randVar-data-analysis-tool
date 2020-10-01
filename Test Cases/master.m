
function master(case_num, n)
    clf;
    clc;
    close all;
    
    vb1 = varBox(0.1, 0.1, 'mixmax', 25);
    vb2 = varBox(0.05, 0.5, 'mixmax', 25);
    vb3 = varBox(1, 0, 'mixmax', 25);
    vb4 = varBox(0.5, 0.5, 'mixmax', 25);
    
    switch(case_num)
        case 1
            case1(vb1);
            input('-CONTINUE-\n');
            case1(vb2);
            input('-CONTINUE-\n');
            case1(vb3);
            input('-CONTINUE-\n');
            case1(vb4);
            input('-CONTINUE-\n');
        case 2
            case2(vb1);
            input('-CONTINUE-\n');
            case2(vb2);
            input('-CONTINUE-\n');
            case2(vb3);
            input('-CONTINUE-\n');
            case2(vb4);
            input('-CONTINUE-\n');
        case 3
            case3(vb1);
            input('-CONTINUE-\n');
            case3(vb2);
            input('-CONTINUE-\n');
            case3(vb3);
            input('-CONTINUE-\n');
            case3(vb4);
            input('-CONTINUE-\n');
        case 4
            case4(vb1);
            input('-CONTINUE-\n');
            case4(vb2);
            input('-CONTINUE-\n');
            case4(vb3);
            input('-CONTINUE-\n');
            case4(vb4);
            input('-CONTINUE-\n');
        case 5
            t = zeros(1,n);
            for ii = 1:1:n
                tic;
                case2n(vb1, 2^ii);
                t(ii) = toc;
            end
            
            clf;
            hold on;
            plot(log(100*t), '-b');
            plot(log(100*t), '*m');
    end
            
    
end