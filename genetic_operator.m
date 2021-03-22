function [f,NC]  = genetic_operator(parent_chromosome,pro,mu,mum,nPoints,...
                                alfa, Re, Mach, limit_values,sample_inputs)

% This function is utilized to produce offsprings from parent chromosomes.
% The genetic operators corssover and mutation which are carried out with
% slight modifications from the original design. For more information read
% the document enclosed.

%  Copyright (c) 2009, Aravind Seshadri
%  All rights reserved.
%
%  Redistribution and use in source and binary forms, with or without 
%  modification, are permitted provided that the following conditions are 
%  met:
%
%     * Redistributions of source code must retain the above copyright 
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright 
%       notice, this list of conditions and the following disclaimer in 
%       the documentation and/or other materials provided with the distribution
%      
%  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
%  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
%  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
%  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
%  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
%  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
%  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
%  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
%  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
%  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
%  POSSIBILITY OF SUCH DAMAGE.

[N,M] = size(parent_chromosome);
switch pro
    case 1
        M = 2;
        V = 6;
    case 2
        M = 3;
        V = 12;
    case 3
        M = 2;
        V = 12;
end
p = 1;
was_crossover = 0;
was_mutation = 0;
l_limit = 0;
u_limit = 1;
conv = 1; %Variable to count the non-convergence iterations and save them
          % inside NoC which and at the end save them inside NC to send
          % them to the Main structure for analysis.
for i = 1 : N
    if rand(1) < 0.9
        child_1 = [];
        child_2 = [];
        parent_1 = round(N*rand(1));
        if parent_1 < 1
            parent_1 = 1;
        end
        parent_2 = round(N*rand(1));
        if parent_2 < 1
            parent_2 = 1;
        end
        while isequal(parent_chromosome(parent_1,:),parent_chromosome(parent_2,:))
            parent_2 = round(N*rand(1));
            if parent_2 < 1
                parent_2 = 1;
            end
        end
        parent_1 = parent_chromosome(parent_1,:);
        parent_2 = parent_chromosome(parent_2,:);
        
%*******%Generation of Children%**************        
        for j = 1 : V
            % SBX (Simulated Binary Crossover)
            % Generate a random number
            %IMPORTANT for Airfoil parameters with NEGATIVE values:
            %In the case of airfoils with PARSEC parametrization is
            %necesary to generate children angles ALFA (negative) and by a
            %different way, using a general restriction, as well as for
            %(negative) curvature values
            if pro == 3
                if j == 1
                    u_limit = (limit_values.Rleup(1));
                    l_limit = (limit_values.Rleup(2));
                    u(j) = (u_limit-l_limit).*rand(1)+l_limit;
                    v(j) = rand(1);
                elseif j == 2
                    u_limit = (limit_values.Rlelo(1));
                    l_limit = (limit_values.Rlelo(2));
                    u(j) = (u_limit-l_limit).*rand(1)+l_limit;
                    v(j) = rand(1);
                elseif j == 3
                    u_limit = (limit_values.ALFAte(1));
                    l_limit = (limit_values.ALFAte(2));
                    u(j) = ((u_limit-l_limit).*rand(1)+l_limit);
                    v(j) = rand(1);
                elseif j == 4
                    u_limit = (limit_values.BETAte(1));
                    l_limit = (limit_values.BETAte(2));
                    u(j) = ((u_limit-l_limit).*rand(1)+l_limit);
                    v(j) = rand(1);
                elseif j == 5
                    u_limit = (limit_values.Zte(1));
                    l_limit = (limit_values.Zte(2));
                    u(j) = (u_limit-l_limit).*rand(1)+l_limit;
                    v(j) = rand(1);
                elseif j == 6
                    u_limit = (limit_values.DELTAte(1));
                    l_limit = (limit_values.DELTAte(2));
                    u(j) = (u_limit-l_limit).*rand(1)+l_limit;
                    v(j) = rand(1);                    
                elseif j == 7
                    u_limit = (limit_values.Xup(1));
                    l_limit = (limit_values.Xup(2));
                    u(j) = (u_limit-l_limit).*rand(1)+l_limit;
                    v(j) = rand(1);                    
                elseif j == 8
                    u_limit = (limit_values.Zup(1));
                    l_limit = (limit_values.Zup(2));
                    u(j) = (u_limit-l_limit).*rand(1)+l_limit;
                    v(j) = rand(1);                    
                elseif j == 9
                    u_limit = (limit_values.Zxxup(1));
                    l_limit = (limit_values.Zxxup(2));
                    u(j) = (u_limit-l_limit).*rand(1)+l_limit;
                    v(j) = rand(1);
                elseif j == 10
                    u_limit = (limit_values.Xlo(1));
                    l_limit = (limit_values.Xlo(2));
                    u(j) = (u_limit-l_limit).*rand(1)+l_limit;
                    v(j) = rand(1);                    
                elseif j == 11
                    u_limit = (limit_values.Zlo(1));
                    l_limit = (limit_values.Zlo(2));
                    u(j) = (u_limit-l_limit).*rand(1)+l_limit;
                    v(j) = rand(1);
                else
                    u_limit = (limit_values.Zxxlo(1));
                    l_limit = (limit_values.Zxxlo(2));
                    u(j) = (u_limit-l_limit).*rand(1)+l_limit;
                    v(j) = rand(1);
                end   
            if v(j) <= 0.5
                bq(j) = (2*u(j))^(1/(mu+1));
                %This will make sure the values of bq are always real
                %numbers
                test_for_real = isreal(bq);
                if test_for_real == 0
                    bq = real(bq);
                end
            else
                bq(j) = (1/(2*(1 - u(j))))^(1/(mu+1));
                %This will make sure the values of bq are always real
                %numbers
                test_for_real = isreal(bq);
                if test_for_real == 0
                     bq = real(bq);
                end                
            end
            child_1(j) = ...
                0.5*(((1 + bq(j))*parent_1(j)) + (1 - bq(j))*parent_2(j));
            child_2(j) = ...
                0.5*(((1 - bq(j))*parent_1(j)) + (1 + bq(j))*parent_2(j));
            %IMPORTANT: The children should be tested, in order to have the
            %variables inside the specified limits.
            if child_1(j) > u_limit
                child_1(j) = u_limit;
            elseif child_1(j) < l_limit
                child_1(j) = l_limit;
            end
            if child_2(j) > u_limit
                child_2(j) = u_limit;
            elseif child_2(j) < l_limit
                child_2(j) = l_limit;
            end
            l_limit = 0;
            u_limit = 1; 
%***********%For normal optimization problems%****************
            else
                u(j) = rand(1);
            
                if u(j) <= 0.5
                   bq(j) = (2*u(j))^(1/(mu+1));
                else
                   bq(j) = (1/(2*(1 - u(j))))^(1/(mu+1));
                end
                child_1(j) = ...
                0.5*(((1 + bq(j))*parent_1(j)) + (1 - bq(j))*parent_2(j));
                child_2(j) = ...
                0.5*(((1 - bq(j))*parent_1(j)) + (1 + bq(j))*parent_2(j));
                if child_1(j) > u_limit
                    child_1(j) = u_limit;
                elseif child_1(j) < l_limit
                    child_1(j) = l_limit;
                end
                if child_2(j) > u_limit
                    child_2(j) = u_limit;
                elseif child_2(j) < l_limit
                    child_2(j) = l_limit;
                end
            end
        end
%*******%Evaluation of children 1 & 2 %************
        if pro==3
           n = nPoints; 
           child_1 = transpose(child_1);
           child_2 = transpose(child_2);
           coord(:,:,1) = parsec(n,child_1);
           coord(:,:,2) = parsec(n,child_2);
           plot(coord(:,1,1),coord(:,2,1),'g')
           axis([0,1.1,-.2,0.4])
           title(i)
           xlabel('x/c')
           ylabel('y/c')
           plot(coord(:,1,2),coord(:,2,2),'g')
           axis([0,1.1,-.2,0.4])
           title(i)
           xlabel('x/c')
           ylabel('y/c')           
           child_1 = transpose(child_1);
           child_2 = transpose(child_2);
           [child_1(:,V + 1: M + V),NoC(:,conv)] = evaluate_objective(coord(:,:,1),pro,...
                                                        alfa, Re,...
                                                        Mach,sample_inputs);
                                                    conv = conv + 1;                                                  
           [child_2(:,V + 1: M + V),NoC(:,conv)] = evaluate_objective(coord(:,:,2),pro,...
                                                        alfa, Re,...
                                                        Mach,sample_inputs);
                                                    conv = conv + 1;
           was_crossover = 1;
           was_mutation = 0;
        
        else
            [child_1(:,V + 1: M + V),NC] = evaluate_objective(child_1,pro);
            [child_2(:,V + 1: M + V),NC] = evaluate_objective(child_2,pro);
            was_crossover = 1;
            was_mutation = 0;
        end
    else
        parent_3 = round(N*rand(1));
        if parent_3 < 1
            parent_3 = 1;
        end
        % Make sure that the mutation does not result in variables out of
        % the search space. For both the MOP's the range for decision space
        % is [0,1]. In case different variables have different decision
        % space each variable can be assigned a range.
        child_3 = parent_chromosome(parent_3,:);
        for j = 1 : V
           r(j) = rand(1);
           if r(j) < 0.5
               delta(j) = (2*r(j))^(1/(mum+1)) - 1;
           else
               delta(j) = 1 - (2*(1 - r(j)))^(1/(mum+1));
           end
           child_3(j) = child_3(j) + delta(j);
%*******%For Airfoils only using PARSEC%***********************************
            if pro == 3
                if j == 1
                    u_limit = (limit_values.Rleup(1));
                    l_limit = (limit_values.Rleup(2));
                elseif j == 2
                    u_limit = (limit_values.Rlelo(1));
                    l_limit = (limit_values.Rlelo(2));
                elseif j == 3
                    u_limit = (limit_values.ALFAte(1));
                    l_limit = (limit_values.ALFAte(2));
                elseif j == 4
                    u_limit = (limit_values.BETAte(1));
                    l_limit = (limit_values.BETAte(2));
                elseif j == 5
                    u_limit = (limit_values.Zte(1));
                    l_limit = (limit_values.Zte(2));
                elseif j == 6
                    u_limit = (limit_values.DELTAte(1));
                    l_limit = (limit_values.DELTAte(2));                    
                elseif j == 7
                    u_limit = (limit_values.Xup(1));
                    %% Why is l_limit deactivated?
                    l_limit = (limit_values.Xup(2));                    
                elseif j == 8
                    u_limit = (limit_values.Zup(1));
                    l_limit = (limit_values.Zup(2));                   
                elseif j == 9
                    u_limit = (limit_values.Zxxup(1));
                    l_limit = (limit_values.Zxxup(2));
                elseif j == 10
                    u_limit = (limit_values.Xlo(1));
                    l_limit = (limit_values.Xlo(2));                    
                elseif j == 11
                    u_limit = (limit_values.Zlo(1));
                    l_limit = (limit_values.Zlo(2));
                else
                    u_limit = (limit_values.Zxxlo(1));
                    l_limit = (limit_values.Zxxlo(2));
                end
            end
%*******%For 1 & 2 & 3 Optimization Problems%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           if child_3(j) > u_limit
               child_3(j) = u_limit;
           elseif child_3(j) < l_limit
               child_3(j) = l_limit;
           end
        end
        %*******%Evaluation of children%************
        if pro==3
           n = nPoints; 
           child_3 = transpose(child_3);
           coord(:,:,3) = parsec(n,child_3);
           plot(coord(:,1,3),coord(:,2,3),'g')
           axis([0,1.1,-.2,0.4])
           title(i)
           xlabel('x/c')
           ylabel('y/c')
           child_3 = transpose(child_3);
           [child_3(:,V + 1: M + V),NoC(:,conv)] = evaluate_objective(coord(:,:,3),pro,...
                                                        alfa, Re,...
                                                        Mach,sample_inputs);
                                                    conv = conv + 1;
                                    NoC(:,conv) = zeros(sample_inputs(1),1);
                                                    conv = conv + 1;
           was_crossover = 0;
           was_mutation = 1;
        else
        child_3(:,V + 1: M + V) = evaluate_objective(child_3,pro);
        was_mutation = 1;
        was_crossover = 0;
        end
    end
    if was_crossover
        child(p,:) = child_1;
        child(p+1,:) = child_2;
        was_cossover = 0;
        p = p + 2;
    elseif was_mutation
        child(p,:) = child_3(1,1 : M + V);
        was_mutation = 0;
        p = p + 1;
    end
end
f = child;
NC = NoC;