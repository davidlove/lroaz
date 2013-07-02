function simpleLP = InitializeSimpleTwoStageLP()

simpleLP = LPModel('2 stage',4);

% First stage data
simpleLP.Setb(0);
simpleLP.SetA(0);
simpleLP.Setl(-1);
simpleLP.Setu(1);
simpleLP.Setc(3.24e6);
% simpleLP.Setc(1);

% Second stage costs
simpleLP.Setq( 3.59e7 * ones(1,4), 1 );
simpleLP.Setq( 3.53e7 * ones(1,4), 2 );
simpleLP.Setq( 2.22e7 * ones(1,4), 3 );
simpleLP.Setq( 2.19e7 * ones(1,4), 4 );
% simpleLP.Setq( ones(1,4), 1 );
% simpleLP.Setq( ones(1,4), 2 );
% simpleLP.Setq( ones(1,4), 3 );
% simpleLP.Setq( ones(1,4), 4 );

% Second stage constant RHS
simpleLP.Setd( [2;0]  , 1 );
simpleLP.Setd( [-2;4] , 2 );
simpleLP.Setd( [1;-1] , 3 );
simpleLP.Setd( [-2;-3], 4 );

% Second stage technology matrices
simpleLP.SetB( [3;2] , 1 );
simpleLP.SetB( [3;2] , 2 );
simpleLP.SetB( [4;2] , 3 );
simpleLP.SetB( [3;-2], 4 );

% Second stage constraint matrix, lower and upper bounds are identical
for ii = 1:4
    simpleLP.Setl2( zeros(4,1), ii );
    simpleLP.Setu2( Inf(4,1), ii );
    simpleLP.SetD( [eye(2), -eye(2)], ii );
end
