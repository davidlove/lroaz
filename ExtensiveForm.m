function [c,A,b,l,u] = ExtensiveForm(lp)

c = lp.c;
b = lp.b;
l = lp.l;
u = lp.u;
[i,j,v] = find(lp.A);
if isempty(i)
    [i,j] = size(lp.A);
    v = 0;
end
prob = 1./lp.numScenarios;

for ss=1:lp.numScenarios
    q = lp.Getq(ss);
    q = q(:)';
    c = [c, q];
    b = [b; lp.Getd(ss)];
    l = [l; lp.Getl2(ss)];
    u = [u; lp.Getu2(ss)];
    
    rows = max(i);
    cols = max(j);
    [i2,j2,v2] = find(-lp.GetB(ss));
    i = [i; i2+rows];
    j = [j; j2];
    v = [v; v2];
    [i3,j3,v3] = find(lp.GetD(ss));
    i = [i; i3+rows];
    j = [j; j3+cols];
    v = [v; v3];
end
A = sparse(i,j,v);

% x = linprog(c,[],[],A,b,l,u);