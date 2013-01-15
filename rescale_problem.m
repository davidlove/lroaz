function scale = rescale_problem(firstStageCost,indivScens,numscen,Nbar)

solnCost = firstStageCost + sum(numscen/sum(numscen).*indivScens');
scale = 2*Nbar/solnCost;