function [f,NC] = evaluate_objective(x,problem, alfa, Re, ~,sample_inputs)

% Function to evaluate the objective functions for the given input vector
% x. x has the decision variables

switch problem
    case 1
        f = [];
        %% Objective function one
        f(1) = 1 - exp(-4*x(1))*(sin(6*pi*x(1)))^6;
        sum = 0;
        for i = 2 : 6
            sum = sum + x(i)/4;
        end
        %% Intermediate function
        g_x = 1 + 9*(sum)^(0.25);
        %% Objective function one
        f(2) = g_x*(1 - ((f(1))/(g_x))^2);
    case 2
        f = [];
        %% Intermediate function
        g_x = 0;
        for i = 3 : 12
            g_x = g_x + (x(i) - 0.5)^2;
        end
        %% Objective function one
        f(1) = (1 + g_x)*cos(0.5*pi*x(1))*cos(0.5*pi*x(2));
        %% Objective function two
        f(2) = (1 + g_x)*cos(0.5*pi*x(1))*sin(0.5*pi*x(2));       
        %% Objective function three
        f(3) = (1 + g_x)*sin(0.5*pi*x(1));
    case 3
        No_convergence = zeros;
        F = zeros;
        samples = Sampling(sample_inputs(1),sample_inputs(2),...
                            sample_inputs(3),sample_inputs(4));
        s = length(samples);
        f = [];
        no = 1;
        Mach = samples;
            coord = x;
               
            pol = xfoil(coord, alfa, Re, Mach,...
                               'PCOP','PANE','oper iter 120',...
                               'oper/vpar n 9');
       for i=1:s
        a = pol.CL(:,i);
        b = pol.CD(:,i);
%         c = s;
        %Individual will be penalized by each time the result of the CFD
        %solver does not converge by the following criteria:
        %mu = 0; sigma = 0; This way all samples will be evaluated and the
        %ones with no result will impact negatively in the total fitness of
        %the indiviual as mean value and standard deviation.
        if a == 1000 || b == 1000
              F(i) = 0; %mean value & standard deviation penalization
%        %Is important to know at which values of Mach, the solver is not
         %able to converge for the specific geometry.
         No_convergence(no,1) = pol.Mach(:,i);
         no = no + 1;
        else
              F(i) = (pol.CL(:,i)./pol.CD(:,i));
              No_convergence(no,1) = 0;
              no = no + 1;
            %When using CD, to MINIMIZE it, multiply by *-1
        end
        end
        f(1) = mean(F)*-1;
        f(2) = std(F);
        NC = No_convergence;
end
