function lp = InitializeCSVLP(directory)

files = dir([directory 'B_*.csv']);
n = length(files);
lp = LPModel('2 stage',n);

scens = zeros(1,n);
for ii = 1:length(scens) 
    scens(ii) = sscanf(files(ii).name,'B_%d.csv'); 
end
scens = sort(scens);

lp.Setb(load([directory 'b.csv']));
lp.Setc(load([directory 'c.csv'])');
lp.Setl(load([directory 'l.csv']));
lp.Setu(load([directory 'u.csv']));
temp = load([directory 'A.csv']);
lp.SetA(sparse(temp(:,1),temp(:,2),temp(:,3)));
nc = max(temp(:,2));

for ss = scens
    lp.Setd(load([directory 'd_' num2str(ss) '.csv']), ss);
    lp.Setq(load([directory 'q_' num2str(ss) '.csv']), ss);
    lp.Setl2(load([directory 'l2_' num2str(ss) '.csv']), ss);
    lp.Setu2(load([directory 'u2_' num2str(ss) '.csv']), ss);
    temp = load([directory 'D_' num2str(ss) '.csv']);
    lp.SetD(sparse(temp(:,1),temp(:,2),temp(:,3)), ss);
    nr = max(temp(:,1));
    temp = load([directory 'B_' num2str(ss) '.csv']);
    lp.SetB(sparse(temp(:,1),temp(:,2),temp(:,3),nr,nc), ss);
end