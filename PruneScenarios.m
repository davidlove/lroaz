function outLP = PruneScenarios( lp, scens )
%PruneScenarios Generate lp with only selected scenarios

n = length(scens);
outLP = LPModel('2 stage',n);

if max(scens) > lp.numScenarios
    error(['lp only has ' num2str(lp.numScenarios) ' scenarios, cannot '...
        'satisfy request for scenario ' num2str(max(scens))])
end

outLP.Setb(lp.b);
outLP.Setc(lp.c);
outLP.Setl(lp.l);
outLP.Setu(lp.u);
outLP.SetA(lp.A);

for ss = 1:n
    outLP.Setd(lp.Getd(scens(ss)), ss);
    outLP.Setq(lp.Getq(scens(ss)), ss);
    outLP.Setl2(lp.Getl2(scens(ss)), ss);
    outLP.Setu2(lp.Getu2(scens(ss)), ss);
    outLP.SetD(lp.GetD(scens(ss)), ss);
    outLP.SetB(lp.GetB(scens(ss)), ss);
end

end

