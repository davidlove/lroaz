function TestSolution()

lp = InitializeSimpleTwoStageLP;
phi = PhiDivergence('burg');
s = Solution(lp,phi,'multi');

% Size errors for all variables
assertExceptionThrown( @() s.SetX(zeros(size(lp.A,2)+1,1)), ...
    'Solution:SetX:size')
assertExceptionThrown( @() s.SetLambda([1 1]), 'Solution:SetLambda:size' )
assertExceptionThrown( @() s.SetMu([1 1]), 'Solution:SetMu:size' )
assertExceptionThrown( @() s.SetTheta(zeros(lp.numScenarios-1,1),'master'), ...
    'Solution:SetTheta:size' )

% Sign error for lambda
assertExceptionThrown( @() s.SetLambda(-2), 'Solution:SetLambda:sign' )

% Only master and true for theta
assertExceptionThrown( @() s.SetTheta(rand(1,lp.numScenarios),'blah'), ...
    'Solution:SetTheta:type' )

% Assign legitimate variables
x = 0.5;
l = 0.1;
m = 2;
th = ones(1,lp.numScenarios);
s.SetX(x);
s.SetLambda(l);
s.SetMu(m);
s.SetTheta(th,'master');
assertEqual( s.Limit, phi.limit )

% Not allowed to set X or Lambda a second time
assertExceptionThrown( @() s.SetX(1), 'Solution:SetX:setagain' )
assertExceptionThrown( @() s.SetLambda(5), 'Solution:SetLambda:setagain' )

% Give trust region a non-logical
assertExceptionThrown( @() s.SetTrustRegionInterior( 2 ), ...
    'Solution:SetTrustRegionInterior:logical' )
% Give trust region interior status and check
s.SetTrustRegionInterior( true );
assertTrue( s.TrustRegionInterior )
% Try to assign trust region status again
assertExceptionThrown( @() s.SetTrustRegionInterior( false ), ...
    'Solution:SetTrustRegionInterior:setagain' )

% Bad scenario number tests
assertExceptionThrown( @() s.SetSecondStageValue( 0, 32 ), ...
    'Solution:SetSecondStageValue:badscen' )
assertExceptionThrown( @() s.SetSecondStageValue( lp.numScenarios+3, 32 ), ...
    'Solution:SetSecondStageValue:badscen' )
assertExceptionThrown( @() s.SetSecondStageDual( 0, 32, 'slope' ), ...
    'Solution:SetSecondStageDual:badscen' )
assertExceptionThrown( @() s.SetSecondStageDual( lp.numScenarios+2, 32, 'slope' ), ...
    'Solution:SetSecondStageDual:badscen' )

% Size error for duals
assertExceptionThrown( @() s.SetSecondStageDual( 1, ones(1,size(lp.A,2)+1), 'slope' ), ...
    'Solution:SetSecondStageDual:size' )
assertExceptionThrown( @() s.SetSecondStageDual( 1, [2 3], 'int' ), ...
    'Solution:SetSecondStageDual:size' )

% Assign second stage costs and dual values
v = zeros(lp.numScenarios,1);
ds = v;
di = v;
for scen=1:lp.numScenarios
    v(scen) = scen-1;
    ds(scen) = scen*2;
    di(scen) = scen/3;
    s.SetSecondStageValue( scen, v(scen) );
    s.SetSecondStageDual( scen, ds(scen), 'slope' );
    s.SetSecondStageDual( scen, di(scen), 'int' );
end
s.SetTheta(2*ones(1,lp.numScenarios),'true');

% Return second stage duals and values
s.SecondStageValues;
s.SecondStageSlope(1);
s.SecondStageIntercept(1);

% Check that s variables are being calculated correctly
assertElementsAlmostEqual( s.S(), (v-m)/l, 'relative', 1e-6 )

% Determine whether mu should be feasible
muShouldBeFeasible = max(s.S()) < s.Limit();

% Assert whether mu should be feasible
assertTrue( s.MuFeasible == muShouldBeFeasible )

% Feasibility of Mu should not change when Mu is updated
s.SetMu( max(s.SecondStageValues) + 1 );
assertTrue( s.MuFeasible == muShouldBeFeasible )

% Reset, then assign varaibles again
s.Reset()
s.SetX(-0.5)
s.SetLambda(lp.numScenarios)
s.SetMu(0)
s.SetTheta(4*ones(1,lp.numScenarios),'master')
s.SetTrustRegionInterior( false )

% Assign second stage costs and duals, this time with Mu feasible
for scen=1:lp.numScenarios
    s.SetSecondStageValue( scen, s.Mu + scen - 1 );
    s.SetSecondStageDual( scen, s.Mu - scen*2, 'slope' );
    s.SetSecondStageDual( scen, s.Mu - scen/3, 'int' );
end
s.SetTheta(8*ones(1,lp.numScenarios),'true');

assertTrue( s.MuFeasible )

% ------------------------------------------------------------------------

% Test single cut version

t = Solution( lp, phi, 'single' );

assertExceptionThrown( @() t.SetTheta(zeros(2,1),'master'), ...
    'Solution:SetTheta:size' )

t.SetTheta(1,'master')
