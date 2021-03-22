%% %Robust Design Optimization of Airfoils at Low Re Number
%MAIN
tic
close all
clear all
clc
%% %*****>General Values<*****
nProfiles = 10;  %Number of individuals
nPoints = 64;    %Number of points per profile. MULTIPLES of 4
alfa = 2;        %Angle of attack, if a vector use [a1,a2,...an]  
Re = 1e5;        %Reynolds number
Mach = 0.4735;   %Mach number
pop = nProfiles; %Name of the variable needed for the MOEA
gen = 3;         %Number of generations for MOEA
pro = 3;         %Number of design objectives, 1>pro<=3. For Airfoils Number 3
active = 1;      %variable that initiates the optimization tool
initial_parameters = zeros(nProfiles,12); %initial array of parameters for PARSEC
sample_size = 10; %number of samples, MULTIPLES OF 10
sampling_method = 2; %1: Simple Random Sample. 2: Latin Hypercube Sample
standard_deviation = 0.1; % Corresponds to a variation of 22 m/s of wind speed
mean_value = Mach; %Operation condition considered as uncertainty
sample_inputs = [sample_size,sampling_method,standard_deviation,mean_value];
parametrization = 1; %Use 1 for Parzec or 2 for Bezier-Parzec
%% %********************%Ranges for parameters%**************************%%
%% FOR PARSEC
if parametrization == 1;
Rleup_max= 0.0126;
Rleup_min= 0.0085*0.8;
Rlelo_max= 0.004;
Rlelo_min= 0.002*0.8;
ALFAte_max= 7.0 *pi/180*-1;
ALFAte_min= 10.0 *pi/180*-1;
BETAte_max= 14.0 *pi/180;
BETAte_min= 10.0 *pi/180;
Zte_max= -0.003;
Zte_min= -0.006;
DELTAZte_max= 0.0050;
DELTAZte_min= 0.0025;
Xup_max= 0.46;
Xup_min= 0.41*0.8;
Zup_max= 0.13;
Zup_min= 0.11*0.8;
Zxxup_max= -0.7;
Zxxup_min= -0.9;
Xlo_max= 0.26;
Xlo_min= 0.20;
Zlo_max= -0.015;
Zlo_min= -0.023*0.8;
Zxxlo_max= 0.20;
Zxxlo_min= 0.05;
end
          
limit_values.Rleup = [Rleup_max,Rleup_min];
limit_values.Rlelo = [Rlelo_max,Rlelo_min];
limit_values.ALFAte = [ALFAte_max,ALFAte_min];
limit_values.BETAte = [BETAte_max,BETAte_min];
limit_values.Zte = [Zte_max,Zte_min];
limit_values.DELTAte = [DELTAZte_max,DELTAZte_min];
limit_values.Xup = [Xup_max,Xup_min];
limit_values.Zup = [Zup_max,Zup_min];
limit_values.Zxxup = [Zxxup_max,Zxxup_min];
limit_values.Xlo = [Xlo_max,Xlo_min];
limit_values.Zlo = [Zlo_max,Zlo_min];
limit_values.Zxxlo = [Zxxlo_max,Zxxlo_min];
limit_values.param = parametrization;
while active
%% %*****>Initiate Profiles<*****

obj = MOEA_MAIN_nsga_2(pop, gen, pro, nPoints, alfa, Re, Mach, limit_values, active, ...
                        initial_parameters,sample_inputs);
chromosome = obj.chromosome;
No_convergence = obj.No_convergence;
No_convergence.initial = unique(No_convergence.initial(No_convergence.initial ~= 0));
No_convergence.generations = unique(No_convergence.generations(No_convergence.generations ~= 0));
active = 0;
%% %*****>Final Results<*****%
initial_parameters = chromosome(:, 1:12);
obj = profiles(nProfiles,nPoints,limit_values, active, initial_parameters);
coord = obj.coord;

for i = 1:nProfiles
    figure 
    plot(coord(:,1,i),coord(:,2,i),'b')
    axis([0,1.1,-.2,0.4])
    title(i)
    xlabel('x/c')
    ylabel('y/c')
end  
end
toc
timerVal = tic;
Minutes = toc / 60;
display(Minutes);
save Time.txt Minutes -ASCII
save Aerofoil_coord.txt coord -ASCII