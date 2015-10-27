function [xq,fval,exitflag,output,population,score] = StartGA2(nvars,lb,ub,PopulationSize_Data, ...
  Generations_Data, StallGenLimit_Data)
%% This function for calling Genetic Algorithm for problem

%% Start with the default options
% Pearson correlation,      von Liebig approach
options = gaoptimset;
%% Modify options setting
options = gaoptimset(options,'PopulationSize', PopulationSize_Data);
options = gaoptimset(options,'Generations', Generations_Data);
options = gaoptimset(options,'StallGenLimit', StallGenLimit_Data);
options = gaoptimset(options,'Display', 'diagnose');
options = gaoptimset(options,'EliteCount', 15);
options = gaoptimset(options,'CrossoverFraction', 0.9);
%% options = gaoptimset(options,'PlotFcns', {  @gaplotbestf  }); %     @gaplotbestindiv
[xq,fval,exitflag,output,population,score] = ga(@CORPEAR2, nvars,[],[],[],[],lb,ub,[],[],options);
