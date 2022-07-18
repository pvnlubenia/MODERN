% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                                           %
%    odeToCRN                                                               %
%                                                                           %
%                                                                           %
% OUTPUT: Returns the chemical reaction network corresponding to a system   %
%    of ordinary differential equations (ODEs) with mass action kinetics.   %
%    If the system does not satisfy the Hars-Toth criterion, a message      %
%    appears stating what needs to be added in a flux. The output variable  %
%    'ode' allows the user to view the complete system of ODEs with all the %
%    species and fluxes listed in the 'species' and 'flux' fields,          %
%    respectively.                                                          %
%                                                                           %
% INPUT: ode: a structure, representing the system of ODEs (see README.txt  %
%    for details on how to fill out the structure)                          %                   %
%                                                                           %
% Note: It is assumed that the ODE system is a MASS ACTION SYSTEM.          %
%                                                                           %
% Reference: Chellaboina V, Bhat S, Haddad W, Bernstein D (2009) Modeling   %
%    and analysis of mass-action kinetics. IEEE Control Syst 29(4):60-78.   %
%    https://doi.org/10.1109/MCS.2009.932926                                %
%                                                                           %
% Created: 9 June 2022                                                      %
% Last Modified: 18 July 2022                                               %
%                                                                           %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



function ode = odeToCRN(ode)
    
    %
    % Step 1: Create a list of all species indicated in the reactions
    %

    % Initialize list of species
    ode.species = { };

    % Get all species from each equation
    for i = 1:numel(ode.equation)
        ode.species{end+1} = ode.equation(i).id.var;
    end

    % Count the number of species
    m = numel(ode.species);





    %
    % Step 2: Create a list of all fluxes indicated in the reactions
    %

    % Initialize list of fluxes
    ode.flux = { };

    % Get all fluxes from each kinetics
    for i = 1:numel(ode.reaction)
        ode.flux{end+1} = ode.reaction(i).id.var;
    end

    % Count the number of fluxes/reactions
    r = numel(ode.flux);





    %
    % Step 3: Hars-Toth Criterion
    %    Basically, it is saying that a chemical reaction network (CRN) exists for an ODE system if and only if a state variable's outgoing flux contains the state variable
    %

    % Go through each equation
    for i = 1:numel(ode.equation)
        
        % Go through each flux in the equation
        for j = 1:numel(ode.equation(i).flux)
            
            % If the flux has a negative coefficient
            if ode.equation(i).coeff{j} < 0
                
                % Get the reaction number of the flux in question
                % {[ode.kinetics(:).id].var} collects all reaction 'id.var' names
                % ode.equation(i).flux{j} is the flux in question
                collect_reaction_ids = [ode.reaction(:).id];
                num = find(ismember({collect_reaction_ids.var}, ode.equation(i).flux{j}));
                
                % Make sure the flux has the state variable of the equation
                % (e.g., if the equation is dX1/dt, the flux should have X1 in it)
                checker = ismember(ode.reaction(num).species, ode.equation(i).id.var);
                
                if checker == 0
                    disp(['Make sure outgoing flux ' ode.equation(i).flux{j} ' from state variable ' ode.equation(i).id.var ' contains ' ode.equation(i).id.var '.']);
                    return
                end
            end
        end
    end





    %
    % Step 4: Construct matrix A (the kinetic orders)
    %

    % Initialize matrix A
    A = [ ];

    % Go through each kinetics
    for i = 1:numel(ode.reaction)
        
        % Add a row of zeros
        A(end+1, :) = zeros(1, numel(ode.equation));
        
        % Fill out the kinetic order of all the species
        for j = 1:numel(ode.reaction(i).species)
            A(end, find(strcmp(ode.reaction(i).species{j}, ode.species))) = ode.reaction(i).kinetic{j};
        end
    end





    %
    % Step 5: Construct THE stoichiometric matrix N = (B - A)'
    %

    % Initialize matrix N
    N = [ ];

    % Go through each equation
    for i = 1:numel(ode.equation)
        
        % Add a row of zeros
        N(end+1, :) = zeros(1, r);
        
        % Fill out the stochiometric coefficients
        for j = 1:numel(ode.equation(i).flux)
            N(end, find(strcmp(ode.equation(i).flux{j}, ode.flux))) = ode.equation(i).coeff{j};
        end
    end





    %
    % Step 6: Solve for B using N and A
    %

    B = N' + A;





    %
    % Step 7: Construct the CRN
    %

    % Go through each row of A and B
    for i = 1:size(A, 1)
        
        % Locate the index of nonzero entries of A and B
        A_nnz = find(A(i, :));
        B_nnz = find(B(i, :));
        
        % If the row of A is 0
        if isempty(A_nnz)
            
            % If the row of B is also 0
            if isempty(B_nnz)
                disp('There should be no reaction 0 -> 0.');
            
            % Inflow reaction
            else
                for j = 1:numel(B_nnz)
                    
                    % Start the first species of the product complex
                    if j == 1
                        
                        % Do not show coefficient if it is 1
                        if B(i, B_nnz(j)) == 1
                            R = ['R' num2str(i) ': 0 -> ' ode.species{B_nnz(j)}];
                        else
                            R = ['R' num2str(i) ': 0 -> ' num2str(B(i, B_nnz(j))) ode.species{B_nnz(j)}];
                        end
                    
                    % Add the succeeding species into the product complex
                    else
                        
                        % Do not show coefficient if it is 1
                        if B(i, B_nnz(j)) == 1
                            R = [R '+' ode.species{B_nnz(j)}];
                        else
                            R = [R '+' num2str(B(i, B_nnz(j))) ode.species{B_nnz(j)}];
                        end
                    end
                end
                
                % Display the reaction
                disp(R);
            end
        
        % The row of A is not 0
        else
            for j = 1:numel(A_nnz)
                
                % Start the first species of the reactant complex
                if j == 1
                    
                    % Do not show coefficient if it is 1
                    if A(i, A_nnz(j)) == 1
                        R = ['R' num2str(i) ': ' ode.species{A_nnz(j)}];
                    else
                        R = ['R' num2str(i) ': ' num2str(A(i, A_nnz(j))) ode.species{A_nnz(j)}];
                    end
                
                % Add the succeeding species into the reactant complex
                else
                    
                    % Do not show coefficient if it is 1
                    if A(i, A_nnz(j)) == 1
                        R = [R '+' ode.species{A_nnz(j)}];
                    else
                       R = [R '+' num2str(A(i, A_nnz(j))) ode.species{A_nnz(j)}];
                    end
                end
            end
            
            % Outflow reaction
            if isempty(B_nnz)
                disp([R ' -> 0']);
            
            % Construct the product complex
            else
                for j = 1:numel(B_nnz)
                    
                    % Start the first species of the product complex
                    if j == 1
                        
                        % Do not show coefficient if it is 1
                        if B(i, B_nnz(j)) == 1
                            R = [R ' -> ' ode.species{B_nnz(j)}];
                        else
                            R = [R ' -> ' num2str(B(i, B_nnz(j))) ode.species{B_nnz(j)}];
                        end
                    
                    % Add the succeeding species into the product complex
                    else
                        
                        % Do not show coefficient if it is 1
                        if B(i, B_nnz(j)) == 1
                            R = [R '+' ode.species{B_nnz(j)}];
                        else
                            R = [R '+' num2str(B(i, B_nnz(j))) ode.species{B_nnz(j)}];
                        end
                    end
                end
                
                % Display the reaction
                disp(R);
            end
        end
    end

end