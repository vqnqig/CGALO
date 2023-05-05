
function [Elite_antlion_position,Convergence_curve]=CGALO(N,MaxFEs,lb,ub,dim,fobj)

% Initialize the positions of antlions and ants
antlion_position=initialization(N,dim,ub,lb);
ant_position=initialization(N,dim,ub,lb);
FEs=0;
delta = 0.1;
% Initialize variables to save the position of elite, sorted antlions,
% convergence curve, antlions fitness, and ants fitness
Sorted_antlions=zeros(N,dim);
Elite_antlion_position=zeros(1,dim);
Elite_antlion_fitness=inf;
Convergence_curve=[];
antlions_fitness=zeros(1,N);
ants_fitness=zeros(1,N);

% Calculate the fitness of initial antlions and sort them
for i=1:size(antlion_position,1)
    antlions_fitness(1,i)=fobj(antlion_position(i,:));
    FEs=FEs+1;
end

[sorted_antlion_fitness,sorted_indexes]=sort(antlions_fitness);

for newindex=1:N
    Sorted_antlions(newindex,:)=antlion_position(sorted_indexes(newindex),:);
end

Elite_antlion_position=Sorted_antlions(1,:);
Elite_antlion_fitness=sorted_antlion_fitness(1);

% Main loop start from the second iteration since the first iteration
% was dedicated to calculating the fitness of antlions
Current_iter=2;
while FEs<MaxFEs

    % This for loop simulate random walks
    for i=1:size(ant_position,1)
        % Select ant lions based on their fitness (the better anlion the higher chance of catching ant)
        Rolette_index=RouletteWheelSelection(1./sorted_antlion_fitness);
        if Rolette_index==-1
            Rolette_index=1;
        end

        % RA is the random walk around the selected antlion by rolette wheel
        RA=Random_walk_around_antlion(dim,MaxFEs,lb,ub, Sorted_antlions(Rolette_index,:),FEs);

        % RA is the random walk around the elite (best antlion so far)
        [RE]=Random_walk_around_antlion(dim,MaxFEs,lb,ub, Elite_antlion_position(1,:),FEs);

        ant_position(i,:)= (RA(FEs,:)+RE(FEs,:))/2; % Equation (2.13) in the paper
      
        eta = FEs/MaxFEs;
        ant_position(i,:) = ant_position(i,:) * (1 + delta * ( eta * randn() + (1-eta) * trnd(1) ));
    end



    for i=1:size(ant_position,1)

        % Boundar checking (bring back the antlions of ants inside search
        % space if they go beyoud the boundaries
        Flag4ub=ant_position(i,:)>ub;
        Flag4lb=ant_position(i,:)<lb;
        ant_position(i,:)=(ant_position(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        ants_fitness(1,i)=fobj(ant_position(i,:));
        FEs=FEs+1;
    end

    % Update antlion positions and fitnesses based of the ants (if an ant
    % becomes fitter than an antlion we assume it was cought by the antlion
    % and the antlion update goes to its position to build the trap)
    double_population=[Sorted_antlions;ant_position];
    double_fitness=[sorted_antlion_fitness ants_fitness];

    [double_fitness_sorted I]=sort(double_fitness);
    double_sorted_population=double_population(I,:);

    antlions_fitness=double_fitness_sorted(1:N);
    Sorted_antlions=double_sorted_population(1:N,:);

    % Update the position of elite if any antlinons becomes fitter than it
    if antlions_fitness(1)<Elite_antlion_fitness
        Elite_antlion_position=Sorted_antlions(1,:);
        Elite_antlion_fitness=antlions_fitness(1);
    end

    % Keep the elite in the population
    Sorted_antlions(1,:)=Elite_antlion_position;
    antlions_fitness(1)=Elite_antlion_fitness;

    % Update the convergence curve
    Convergence_curve(Current_iter)=Elite_antlion_fitness;
    Convergence_curve(1)=Convergence_curve(2);

    % Display the iteration and best optimum obtained so far
    %     if mod(Current_iter,50)==0
    %         display(['At iteration ', num2str(Current_iter), ' the elite fitness is ', num2str(Elite_antlion_fitness)]);
    %     end

    Current_iter=Current_iter+1;
end






