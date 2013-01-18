function plot_p_dec_v_gammaprime(gp,numScen,pWorst)

% plot_p_dec_v_gammaprime plots the probability of a decrease in cost
% agains the value gammaprime output for LRLP-2.

lgp = length(gp);
lns = length(numScen);
N = sum(numScen);
pMLE = numScen/N;
pMLE = pMLE(:);
numAtoms = length(numScen);

prob_dec = zeros(size(gp));

% Nl stores the upper bound probabilty of decrease, sum(numscen with
% decrease)/numsamples
Nl = zeros(size(gp));

for ii=1:length(gp)
    dec_cost = ( pMLE > (N+1)/N*pWorst(:,ii) )';
    prob_dec(ii) = prob_dec_cost( numScen, dec_cost, gp(ii) );
    Nl(ii) = sum(dec_cost.*numScen)/N;
end

plot( gp,prob_dec,'k.', 'MarkerSize',10)
hold on
% plot( gp,prob_dec_test,'ro')
% plot(gp,prob_dec_primal,'ro' )
plot(gp,Nl,'r--', 'LineWidth',2)
hold off
xlabel( '\gamma''', 'FontSize',16 )
ylabel( 'Probability of Cost Decrease', 'FontSize',16)
legend( 'Prob. of Cost Decrease', 'N_L/N' )
title( [ num2str(numAtoms) ' scenarios, ' num2str(N) ' samples'], ...
    'FontSize',16)
% figure(2)
% plot( gp,gp_pworst,'bo', gp,gp,'k--' )
% xlabel( '\gamma''', 'FontSize',16 )
% ylabel( '\gamma'' of worst case distribution', 'FontSize',16 )
