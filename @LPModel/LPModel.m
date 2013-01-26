% LPModel constructs and stores all the matrices and vectors needed to
% formulate a 1- or 2-stage LP.  It can construct the LP from either matrix
% data entered through MATLAB, or using Alicia Forrester's code to read
% water models from excel files.

classdef LPModel < handle
    
%     First stage/single stage LP matrices
    properties %(GetAccess=public, SetAccess=private)
        c;
        A;
        b;
        l;
        u;
        variableNames;
    end
    
%     Second stage LP matrices
    properties %(Access=private)
        q
        D
        d
        l2
        u2
        B
    end
    
%     Properties of the LP Model
    properties %(GetAccess=public, SetAccess=immutable)
        numStages
        numScenarios
        timePeriods
        timeLag
        firstStagePeriods
        folderCellArray
    end
    
    methods
        
%         The constructor determines what type of LP model to work with:
%           1. One-stage explicit model
%           2. Two-stage explicit model
%           3. One-stage excel water model
%           4. Two-stage excel water model
%         And sets the immutable data accordingly
        function obj = LPModel(varargin)
            if nargin == 0
                error('Must supply initialization arguments')
            end            
            if ischar(varargin{1})
                switch varargin{1}
                    case '1 stage'
                        obj.numStages = 1;
                        obj.numScenarios = 0;
                        obj.timePeriods = 0;
                        obj.timeLag = 0;
                        obj.firstStagePeriods = 0;
                    case '2 stage'
                        obj.numStages = 2;
                        obj.numScenarios = varargin{2};
                        obj.timePeriods = 0;
                        obj.timeLag = 0;
                        obj.firstStagePeriods = 0;
                        obj.SetBlankSecondStage;
                    otherwise
                        obj.folderCellArray = { varargin{1} };
                        obj.numStages = 1;
                        obj.numScenarios = 0;
                        obj.timePeriods = varargin{2};
                        obj.timeLag = varargin{3};
                        obj.firstStagePeriods = obj.timePeriods;
                end
            elseif iscellstr(varargin{1})
                obj.folderCellArray = varargin{1};
                obj.numStages = 2;
                obj.numScenarios = length(obj.folderCellArray);
                obj.timePeriods = varargin{2};
                obj.timeLag = varargin{3};
                obj.firstStagePeriods = varargin{4};
                obj.SetBlankSecondStage;
            else
                error('Unrecognized initialization')
            end
        end
    end
    
    methods %(Access=private)
%         SetBlankSecondStage sets all second stage data to blank arrays of
%         the correct size
        function SetBlankSecondStage( obj )
            if obj.numStages == 2
                obj.q = cell( 1, obj.numScenarios );
                obj.D = cell( 1, obj.numScenarios );
                obj.d = cell( 1, obj.numScenarios );
                obj.l2 = cell( 1, obj.numScenarios );
                obj.u2 = cell( 1, obj.numScenarios );
                obj.B = cell( 1, obj.numScenarios );
            else
                error( 'Must be a two stage problem' )
            end
        end
        
        function IsValidScenario( obj, inScenario )
            if inScenario > 0 && inScenario <= obj.numScenarios
                tf = true;
            else
                error(['Input scenario number ' num2str(inScenario) ...
                    ' is wrong.  Scenarios must be numbered 1 to ' ...
                    num2str(obj.numScenarios)])
            end
        end
                
    end
    
    
%     Setting and Accessing Routines
    methods (Access=public)
        
%         Read in data for the first stage
        function Setc( obj, inc )
            obj.c = inc;
        end
        function SetA( obj, inA )
            obj.A = inA;
        end
        function Setb( obj, inb )
            obj.b = inb;
        end
        function Setl( obj, inl )
            obj.l = inl;
        end
        function Setu( obj, inu )
            obj.u = inu;
        end
        
%         Read in data for the second stage
        function Setq( obj, inq, inScenario )
            if obj.IsValidScenario( inScenario )
                obj.q{inScenario} = inq;
            end
        end
        function SetD( obj, inD, inScenario )
            if obj.IsValidScenario( inScenario )
                obj.D{inScenario} = inD;
            end
        end
        function Setd( obj, ind, inScenario )
            if obj.IsValidScenario( inScenario )
                obj.d{inScenario} = ind;
            end
        end
        function Setl2( obj, inl2, inScenario )
            if obj.IsValidScenario( inScenario )
                obj.l2{inScenario} = inl2;
            end
        end
        function Setu2( obj, inu2, inScenario )
            if obj.IsValidScenario( inScenario )
                obj.u2{inScenario} = inu2;
            end
        end
        function SetB( obj, inB, inScenario )
            if obj.IsValidScenario( inScenario )
                obj.B{inScenario} = inB;
            end
        end
        
%         Accessor routines for second stage data
        function outq = Getq( obj, inScenario )
            if obj.IsValidScenario( inScenario )
                outq = obj.q{inScenario};
            end
        end
        function outD = GetD( obj, inScenario )
            if obj.IsValidScenario( inScenario )
                outD = obj.D{inScenario};
            end
        end
        function outd = Getd( obj, inScenario )
            if obj.IsValidScenario( inScenario )
                outd = obj.d{inScenario};
            end
        end
        function outl2 = Getl2( obj, inScenario )
            if obj.IsValidScenario( inScenario )
                outl2 = obj.l2{inScenario};
            end
        end
        function outu2 = Getu2( obj, inScenario )
            if obj.IsValidScenario( inScenario )
                outu2 = obj.u2{inScenario};
            end
        end
        function outB = GetB( obj, inScenario )
            if obj.IsValidScenario( inScenario )
                outB = obj.B{inScenario};
            end
        end
        
    end
end
        