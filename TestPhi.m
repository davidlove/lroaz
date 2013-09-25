function TestPhi()

% (-1,0) is the only ratio all conjugates take values in.
dv = 0.01;
vals = (-1+dv):dv:(0-dv);
dlimit = 1e-6;

rho = 3.2;
q = 1:5;

% LRO:
%  conjugate  = -log(-s)-1
%  conjugate' = -1/s
phi = PhiDivergence( 'lro' );
s = -exp(-vals-1);
assertElementsAlmostEqual( phi.Conjugate(s), vals, 'relative', 1e-6 )
s = 1./vals;
assertElementsAlmostEqual( phi.ConjugateDerivative(s), -vals, 'relative', 1e-6 )
assertElementsAlmostEqual( phi.SecondDerivativeAt1, 1, 'relative', 1e-6 )
if isfinite(phi.limit)
    slim = phi.limit + dlimit*[-1,1];
    testLimits = phi.Conjugate(slim);
    assertTrue( isreal(testLimits(1)) )
    assertTrue( isinf(testLimits(2)) )
end
assertEqual( phi.divergence, 'lro' )
assertEqual( phi.Contribution(0,1), Inf )
assertEqual( phi.Contribution(1,0), 0 )

% Burg Entropy:
%  conjugate  = -log(1-s)
%  conjugate' = 1/(1-s)
phi = PhiDivergence( 'burg' );
assertEqual( phi.Contribution(1,1), 0 )
s = 1 - exp(-vals);
assertElementsAlmostEqual( phi.Conjugate(s), vals, 'relative', 1e-6 )
s = 1 + 1./vals;
assertElementsAlmostEqual( phi.ConjugateDerivative(s), -vals, 'relative', 1e-6 )
assertElementsAlmostEqual( phi.SecondDerivativeAt1, 1, 'relative', 1e-6 )
if isfinite(phi.limit)
    slim = phi.limit + dlimit*[-1,1];
    testLimits = phi.Conjugate(slim);
    assertTrue( isreal(testLimits(1)) )
    assertTrue( isinf(testLimits(2)) )
end
assertEqual( phi.divergence, 'burg' )
assertEqual( phi.Contribution(0,1), Inf )
assertEqual( phi.Contribution(1,0), 1 )

% Kullback-Leibler:
%  conjugate  = exp(s) - 1
%  conjugate' = exp(s)
phi = PhiDivergence( 'kl' );
assertEqual( phi.Contribution(1,1), 0 )
s = log(1+vals);
assertElementsAlmostEqual( phi.Conjugate(s), vals, 'relative', 1e-6 )
s = log(-vals);
assertElementsAlmostEqual( phi.ConjugateDerivative(s), -vals, 'relative', 1e-6 )
assertElementsAlmostEqual( phi.SecondDerivativeAt1, 1, 'relative', 1e-6 )
if isfinite(phi.limit)
    slim = phi.limit + dlimit*[-1,1];
    testLimits = phi.Conjugate(slim);
    assertTrue( isreal(testLimits(1)) )
    assertTrue( isinf(testLimits(2)) )
end
phi.SetComputationLimit(q,rho)
assertElementsAlmostEqual( rho, min(q)/sum(q)*...
    phi.Contribution(phi.ConjugateDerivative(phi.computationLimit),1) )
assertEqual( phi.divergence, 'kl' )
assertEqual( phi.Contribution(0,1), 1 )
assertEqual( phi.Contribution(1,0), Inf )

% Chi^2:
%  conjugate  = 2 - 2sqrt(1-s)
%  conjugate' = 1/sqrt(1-s)
phi = PhiDivergence( 'chi2' );
assertEqual( phi.Contribution(1,1), 0 )
s = 1 - ((2-vals)/2).^2;
assertElementsAlmostEqual( phi.Conjugate(s), vals, 'relative', 1e-6 )
s = 1 - 1./(vals.^2);
assertElementsAlmostEqual( phi.ConjugateDerivative(s), -vals, 'relative', 1e-6 )
assertElementsAlmostEqual( phi.SecondDerivativeAt1, 2, 'relative', 1e-6 )
if isfinite(phi.limit)
    slim = phi.limit + dlimit*[-1,1];
    testLimits = phi.Conjugate(slim);
    assertTrue( isreal(testLimits(1)) )
    assertTrue( isinf(testLimits(2)) )
end
assertEqual( phi.divergence, 'chi2' )
assertEqual( phi.Contribution(0,1), Inf )
assertEqual( phi.Contribution(1,0), 1 )

% Modified Chi^2:
%  conjugate  = -1 for s < -2,  s+s^2/4 for s >= -2
%  conjugate' = 0(s < -2) or 1 + s/2 (s >= -2)
phi = PhiDivergence( 'mchi2' );
assertEqual( phi.Contribution(1,1), 0 )
s = -2 + 2*sqrt(1+vals);
assertElementsAlmostEqual( phi.Conjugate(s), vals, 'relative', 1e-6 )
s = 2*(-vals-1);
assertElementsAlmostEqual( phi.ConjugateDerivative(s), -vals, 'relative', 1e-6 )
assertElementsAlmostEqual( phi.SecondDerivativeAt1, 2, 'relative', 1e-6 )
s = -2.1:-3.03:-40;
assertEqual( phi.Conjugate(s), -ones(size(s)) )
if isfinite(phi.limit)
    slim = phi.limit + dlimit*[-1,1];
    testLimits = phi.Conjugate(slim);
    assertTrue( isreal(testLimits(1)) )
    assertTrue( isinf(testLimits(2)) )
end
phi.SetComputationLimit(q,rho)
assertElementsAlmostEqual( rho, min(q)/sum(q)*...
    phi.Contribution(phi.ConjugateDerivative(phi.computationLimit),1) )
assertEqual( phi.divergence, 'mchi2' )
assertEqual( phi.Contribution(0,1), 1 )
assertEqual( phi.Contribution(1,0), Inf )

% Hellinger Distance:
%  conjugate  = s/(1-s)
%  conjugate' = 1/(s-1)^2
phi = PhiDivergence( 'hellinger' );
assertEqual( phi.Contribution(1,1), 0 )
s = vals./(1+vals);
assertElementsAlmostEqual( phi.Conjugate(s), vals, 'relative', 1e-6 )
s = 1 - 1./sqrt(-vals);
assertElementsAlmostEqual( phi.ConjugateDerivative(s), -vals, 'relative', 1e-6 )
assertElementsAlmostEqual( phi.SecondDerivativeAt1, 1/2, 'relative', 1e-6 )
if isfinite(phi.limit)
    slim = phi.limit + dlimit*[-1,1];
    testLimits = phi.Conjugate(slim);
    assertTrue( isreal(testLimits(1)) )
    assertTrue( isinf(testLimits(2)) )
end
assertEqual( phi.divergence, 'hellinger' )
assertEqual( phi.Contribution(0,1), 1 )
assertEqual( phi.Contribution(1,0), 1 )

% Test error for unknown phi
assertExceptionThrown( @() PhiDivergence( 'fakey-super-made-up' ), ...
    'PhiDivergence:PhiDivergence:Unknown' )