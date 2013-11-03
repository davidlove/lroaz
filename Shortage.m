function [ shortageTotal, shortageByYear ] = ...
    Shortage( solution, varNames, numVars )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if nargin < 3
    numVars = length(varNames);
end

numYears = length(solution)/numVars;
solution = reshape(solution, numVars, numYears);

dummy = ~cellfun(@isempty,regexp(varNames, 'Dummy'));
shortageByYear = sum(solution(dummy,:));
shortageTotal = sum(shortageByYear);

end
