function compare_lro_deterministic( inputLocation, lroFile, periods1 )

% clear get_stage_vectors

if nargin < 3
    periods1 = 10;
    if nargin < 2
        lroFile = 'saved_variables_01_99.mat';
        if nargin < 1
            inputLocation = 'all_scenarios/5/';
        end
    end
end

% one_stage_az uses the code from lroaz_version15 to do a one-stage problem
% that takes a given input

% ConnectionsFile = [ inputLocation 'Connections.xlsx'];
% cellInputFile = { [inputLocation 'Inputs.xlsx'] };

% Problem Parameters
Period = 41;
% periods1 = 41;

% clear get_stage_vectors
[Q] = one_stage_az( inputLocation );

% clear get_stage_vectors
% [c,A_full,b_vec,LB,UB] = get_stage_vectors(1,1, ...
%     ConnectionsFile,cellInputFile,Period,periods1);
[Nc,~] = get_stage_vectors('names');
[~,Ap,~] = get_stage_vectors('base');
stageSize = size(Ap,2);

% [F fval] = linprog( c, [], [], A_full, b_vec, LB, UB );

% [q,D,d,l,u,B] = get_stage_vectors(2,1);
% B = zeros(size(B));

% [S,fval2] = linprog(q,[],[],D,d+B*F,l,u);

gammaPrime = [0.03:0.01:0.99];

lroData = load(lroFile);
x = lroData.x;
gp = lroData.gp;
gpIndex = false(size(gp));
for ii = 1:length(gammaPrime)
    gpIndex = gpIndex | (abs(gp-gammaPrime(ii))<1e-4);
end
F = x(:,gpIndex);

Q = reshape(Q,stageSize,Period);
F = reshape(F,stageSize,periods1,length(gammaPrime));
% S = reshape(S,stageSize,Period-periods1);

notZero = Q(:,periods1) ~= 0;
relChange = zeros(sum(notZero),length(gammaPrime));
for ii = 1:length(gammaPrime)
    relChange(:,ii) = (F(notZero,periods1,ii) - Q(notZero,periods1)) ...
    ./Q(notZero,periods1);
end
figure(1)
plot(find(Q(:,periods1)),relChange,'.')
xlabel('Decision variable number', 'FontSize',16)
ylabel('Relative Difference from Deterministic', 'FontSize',16)
title(['Relative Difference Between Deterministic and LRO Solution ' ...
    'for time period ' num2str(periods1)], 'FontSize',16)
disp(['Max percent change in solution = ' num2str(max(abs(relChange)))])

isZero = ~notZero & (F(:,periods1) >= 1e-3);
disp(['Solution increased from zero: ' Nc(isZero)])
% figure(2)
% plot(F(isZero,periods1),'.')
disp(['Max above zero = ' num2str(max(abs(F(isZero,periods1))))])

numListChanges = 50;

[~,sortIndex] = sort(abs(relChange(:,1)),'descend');
maxIndexRel = sortIndex(1:numListChanges);

maxIndex = zeros(1,numListChanges);
for ii=1:numListChanges
    temp = find(notZero,maxIndexRel(ii));
    maxIndex(ii) = temp(end);
end

disp(' ')
disp(['Decisions in stage ' num2str(periods1) ' with largest '...
 'change from deterministic to LRO'])
for ii = 1:numListChanges
    disp([Nc(maxIndex(ii)), num2str(relChange(maxIndexRel(ii))), ...
        ['(' num2str(maxIndex(ii)) ')']])
end
% disp([Nc(maxIndex)')

% relDiffCost = (c*F(:)+q*S(:))/([c,q]*Q(:));

% disp(['Relative increase in two-stage cost = ' num2str(relDiffCost)])

% Plot a single index to see how it changes


% Note: Storage variables are 103-107
plotIndex = 104;
newAmount = F(plotIndex,periods1,:);
newAmount = newAmount(:);
baseAmount = Q(plotIndex,periods1)*[1,1];
figure(3)
plot(gammaPrime,newAmount,'.', [0,1],baseAmount,'k-')
xlabel('\gamma''', 'FontSize',16)
ylabel('Value of Decision Variable', 'FontSize',16)
title(Nc(plotIndex), 'FontSize',16)
