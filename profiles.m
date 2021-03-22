classdef profiles
    
    properties
        n  %Number of individuals
        profi   %Array in wich the profiles will be saved
        %For Parsec Parametrization
        rba
        rbai
        abs
        bbs
        zbs
        dzbs
        xe
        ze
        zxxe
        xi
        zi
        zxxi
        coord
        points
        limit_values = []   %upper and lower limits for PARSEC parameters
    end
    methods      
        
        function obj = profiles(nProfiles, nPoints, limit_values, active, ...
                                initial_parameters)
            obj.n = nProfiles;
            obj.points = nPoints;
            rng('default');
            rng(1);
            
            for i=1:nProfiles
                if active == 1                             
                  obj.profi(:,i) = [(limit_values.Rleup(1) - limit_values.Rleup(2)).*rand(1) + limit_values.Rleup(2) ...
                                  (limit_values.Rlelo(1) - limit_values.Rlelo(2)).*rand(1) + limit_values.Rlelo(2) ...
                                  ((limit_values.ALFAte(1) - limit_values.ALFAte(2)).*rand(1) + limit_values.ALFAte(2)) ...
                                  (limit_values.BETAte(1) - limit_values.BETAte(2)).*rand(1) + limit_values.BETAte(2) ...
                                  (limit_values.Zte(1) - limit_values.Zte(2)).*rand(1) + limit_values.Zte(2) ...
                                  (limit_values.DELTAte(1) - limit_values.DELTAte(2)).*rand(1) + limit_values.DELTAte(2) ...
                                  (limit_values.Xup(1) - limit_values.Xup(2)).*rand(1) + limit_values.Xup(2) ...
                                  (limit_values.Zup(1) - limit_values.Zup(2)).*rand(1) + limit_values.Zup(2) ...
                                  (limit_values.Zxxup(1) - limit_values.Zxxup(2)).*rand(1) + limit_values.Zxxup(2) ...
                                  (limit_values.Xlo(1) - limit_values.Xlo(2)).*rand(1) + limit_values.Xlo(2) ...
                                  (limit_values.Zlo(1) - limit_values.Zlo(2)).*rand(1) + limit_values.Zlo(2) ...
                                  (limit_values.Zxxlo(1) - limit_values.Zxxlo(2)).*rand(1) + limit_values.Zxxlo(2) ...
                                   ];



                elseif active == 0
                    obj.profi(:,i) = transpose(initial_parameters(i,:));
                end
            
            end
            obj.coord = parametrization(nPoints, obj.profi, nProfiles, ...
                                        active, limit_values);
        end
        
    end
        
end

function coord = parametrization(npo, parameters, npr, ~, ...
                                 limit_values)
    
    if limit_values.param == 1;
         coord = zeros((npo*2-2),2,npr);
         for k=1:npr
             pars = parameters(:,k);
             %Calls the parametrization function to generate coordinates
             xy = parsec(npo, pars);
                  for i = 2:(npo*2-1)
                      coord(i-1,:,k) = xy(i,:);
                  end
                  
          end        
    elseif limit_values.param == 2;
%            coord = zeros((npo*4-3),2,npr);
           for k=1:npr
               pars = parameters(:,k);
               %Calls the parametrization function to generate coordinates
               xy = bezier_parsec(npo, pars);
                for i = 2:length(xy)-1
                      coord(i-1,:,k) = xy(i,:);
                end

           end
    end

end