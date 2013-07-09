function TestLRLP
%TESTLRLP Run tests on LRLP

simpleLP = InitializeSimpleTwoStageLP();

phi = PhiDivergence('burg');
lrlp = LRLP(simpleLP,phi,ones(1,simpleLP.numScenarios));
lrlp.Test()

phi = PhiDivergence('kl');
lrlp = LRLP(simpleLP,phi,ones(1,simpleLP.numScenarios));
lrlp.Test()

phi = PhiDivergence('chi2');
lrlp = LRLP(simpleLP,phi,ones(1,simpleLP.numScenarios));
lrlp.Test()

phi = PhiDivergence('hellinger');
lrlp = LRLP(simpleLP,phi,ones(1,simpleLP.numScenarios));
lrlp.Test()

end

