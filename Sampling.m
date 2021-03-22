%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%SAMPLING METHODS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sample = Sampling(Nsample,Sample_method,standard_value,mean_value)

lower_limit = 0.2;
upper_limit = 0.7;
mean_v = mean_value; %mean value
standard = standard_value; %standard deviation
n = Nsample; %Sampling size

switch Sample_method
%% ********************************************************************* %%
%%%%%%%%%%%%%%%%%%%Normal Distribution Sample%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 1
            rng('shuffle');
            sample = standard.*randn(n,1) + mean_v;
            for i = 1:length(sample)
                while sample(i,1)>upper_limit || sample(i,1)<lower_limit
                    sample(i,1) = standard.*randn(1,1) + mean_v;

                end
            end
%% ********************************************************************* %%
%%%%%%%%%%%%%%%%%%Latin Hypercube Sample%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 2
%%%%%%USE ONLY IF THE MATLAB SESION HAS STATISTICS PACKAGE%%%%%%%%%%%%%%%%%        
%              sample = lhsnorm(mean_v,standard^2,n);
%              for i = 1:length(sample)
%                  while sample(i,1)>upper_limit || sample(i,1)<lower_limit
%                           sample(i,1) = lhsnorm(mean_v,standard^2,1);
% 
%                  end
%              end
%%%%%IN CASE THE USED MATLAB SESION DOES NOT HAVE STATISTICS PACKAGE%%%%%%%             
               ran=rand(n,1);
               sample=zeros(n,1);
               % method of Stein
               % for j=1: 1
                    idx=randperm(n);
                    P=(idx'-ran(:,1))/n;       % probability of the cdf
                    sample(:,1) = mean_v(1) + ltqnorm(P).* standard(1); % this can be replaced by any inverse distribution function
                        for i = 1:length(sample)
                            while sample(i,1)>upper_limit || sample(i,1)<lower_limit
                                  sample(i,1) = mean_v + randn(1,1).* standard;

                            end
                        end
% end
end

%% USE FOR PLOTING DIFFERENT SAMPLES

% standard_deviation = [std(sample_1),std(sample_2),std(sample_3)]
% mean_value = [mean(sample_1),mean(sample_2),mean(sample_3)]
% 
% % plot(sample_1(:,1),'. b')
% % hold on
% plot(sample_2(:,1),'o g')
% hold on
% plot(sample_3(:,1),'^ m')
    

end