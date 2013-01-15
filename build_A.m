function [A,Ap,Nr,Nc,Zones,Cuser] = build_A(ConnectionsFile, InputFile)

%% Load Data
disp('Loading Data...');

[Con, text] = xlsread(ConnectionsFile,'Connections'); % System Connection Matrix
Name.user = text(3:end,1);
Name.source = text(1,4:end);

%% Arrange Connection Matrix by type and zone

disp('Arranging Connection Matrix...');

Zones = Con(2:end,1);
userTyp = Con(2:end,2);
Con = Con(:,3:end);
sourceN = size(Con,2);
userN = size(Con,1)-1;

%preallocate
gwN = 0; wtpN = 0; rchrgN = 0; wwtpN = 0; dmyN = 0; swN = 0; pN = 0;
npN = 0; rtrnN = 0;

gwsource(1:userN,1)= 0; wtpsource(1:userN,1)= 0; rchrgsource(1:userN,1)= 0;
wwtpsource(1:userN,1)= 0; dmysource(1:userN,1)= 0; swsource(1:userN,1)= 0;
psource(1:userN,1)= 0; npsource(1:userN,1)= 0; rtrnsource(1:userN,1)= 0;

gwID(1) = {[]}; wtpID(1) = {[]}; rchrgID(1) = {[]}; wwtpID(1) = {[]};
dmyID(1) = {[]}; swID(1) = {[]}; pID(1) = {[]}; npID(1) = {[]}; rtrnID(1) = {[]};

% Sort Connection Matrix by source type
for i = 1:sourceN
    if Con(1,i) == 1
        gwN = gwN + 1;
        gwsource(:,gwN) = Con(2:userN+1,i);
        gwID(gwN) = Name.source(i);
    elseif Con(1,i) == 2
        wwtpN = wwtpN + 1;
        wwtpsource(:,wwtpN) = Con(2:userN+1,i);
        wwtpID(wwtpN) = Name.source(i);
    elseif Con(1,i) == 3
        swN = swN + 1;
        swsource(:,swN) = Con(2:userN+1,i);
        swID(swN) = Name.source(i);
    elseif Con(1,i) == 4
        wtpN = wtpN + 1;
        wtpsource(:,wtpN) = Con(2:userN+1,i);
        wtpID(wtpN) = Name.source(i);
    elseif Con(1,i) == 5
        dmyN = dmyN + 1;
        dmysource(:,dmyN) = Con(2:userN+1,i);
        dmyID(dmyN) = Name.source(i);
    elseif Con(1,i) == 6
        rchrgN = rchrgN + 1;
        rchrgsource(:,rchrgN) = Con(2:userN+1,i);
        rchrgID(rchrgN) = Name.source(i);
    elseif Con(1,i) == 7 % potable main
        pN = pN + 1;
        psource(:,pN) = Con(2:userN+1,i);
        pID(pN) = Name.source(i);
    elseif Con(1,i) == 8 % non-potable main
        npN = npN + 1;
        npsource(:,npN) = Con(2:userN+1,i);
        npID(npN) = Name.source(i);
    elseif Con(1,i) == 9 % return main
        rtrnN = rtrnN + 1;
        rtrnsource(:,rtrnN) = Con(2:userN+1,i);
        rtrnID(rtrnN) = Name.source(i);
    end
end

S = [gwN,wwtpN,swN,wtpN,dmyN,rchrgN,pN,npN,rtrnN];
cumS = cumsum(S);
Con = [gwsource,wwtpsource,swsource,wtpsource,...
    dmysource,rchrgsource,psource,npsource,rtrnsource];
sourceID = [gwID,wwtpID,swID,wtpID,dmyID,rchrgID,pID,npID,rtrnID];
ST = gwN + swN;

empties = S==0; % identify the empty cells
S(empties) = [];
empties = find(cellfun(@isempty,sourceID));
sourceID(empties) = [];
Con(:,empties) = [];


% Arrange pressure zones

xa = 0;
xb = 0;
% xd = 0;
zuser = zeros(size(Con));
userID = cell(size(Name.user));
zones_new = zeros(size(Zones));
user_type = zeros(size(userTyp));

for xc = 1:max(Zones)
    xb = xb+1;
    for xe = 1:userN
        if Zones(xe) == xc
            xa = xa + 1;
            zuser(xa,:) =  Con(xe,:);
            userID(xa,1) = Name.user(xe);
            zones_new(xa,1) = Zones(xe);
            user_type(xa,1) = userTyp(xe);
%             if userTyp(xe) == 10
%                 xd = xd + 1;
%                 ret(xd,zones_new(xa,1)) = 0.98;
%                 retsourceID(xd,1) = userID(xa,1);
%             elseif userTyp(xe) == 12
%                 xd = xd + 1;
%                 ret(xd,zones_new(xa,1)) = 0.4;
%                 retsourceID(xd,1) = userID(xa,1);
%             end
        end
    end
end


Con = zuser; % resorted connection matrix
Zones = zones_new;
% ret = ret(:,2:end);


%% Check Return Matrix

Return_Matrix = xlsread(InputFile,'Returns'); % Use previous values
%% preallocate
disp('preallocating...');

A = zeros(userN*2+ST,sourceN*userN+userN+sourceN*2+1);
A_st = zeros(userN*2+ST,sourceN*userN+userN+sourceN*2+1);
Nr = cell(2*userN+ST,1);
Nr(1:userN,1) = userID';
Nc = cell(1,sourceN*userN+userN+sourceN*2+1);
% clc;

%% Build Equality Matrix A
disp('Building Matrix...');

% Indexing variables
xb = userN;
% xu = sourceN*userN; % length of source outflows
% clear xc
% xa = 0;
xi = 0;
xk = userN + ST;
xl = 0;
xn = 0;
xo = 0;
xp = 0;
Loss_percentage = 0.011;

% Centralized facilities
index = find(Zones==1 & user_type==2);
xq = length(index); % reclamation
index = find(Zones==1 & user_type==6);
xs = length(index); % recharge

for xd = 1:sourceN % For each source term
    
    if xd <= ST
        xb = xb+1; % start of storage rows in matrix
    end
    xu = sourceN*userN+xd; % storage column
    % xe = userN + ST;
    xg = sourceID(xd);
    
    % Determine zone of current source
    a1 = sourceID(xd);
    a1 = strcmp(userID,a1);
    a1 = a1==1;
    
    a1 = Zones(a1);
    
    source_type = find(xd <= cumS,1,'first');
    
    if source_type == 2
        if xd <= gwN + xq
            xn = xn + 1;
        else
            xn = 1;
        end
        
    elseif source_type == 4
        xo = xo + 1;
        
    elseif source_type == 6
        if xd <= gwN + wwtpN + swN + wtpN + dmyN + xs
            xp = xp +1;
        else
            xp = 1;
        end
    end
    
    for xr = 1:userN % For each user
        
        if Con(xr,xd) == 1
            xh = userID(xr); % user ID
            xa = (xd-1)*userN + xr; % current A matrix column
            xe = userN+ST+xr; % loss constraint row
            xf = Zones(xr); %user zone number
            xj = user_type(xr); %user type
            
            xl = xl + 1;
            
            %             Csource(xl) = xg; % Cost vector source label
            Cuser(xl) = xh; % Cost vector user label
%             Nc(1,xa) = xg; % source ID
%             David Love -- Changed column name to indicate which arc, not
%             which source!
            Nc(1,xa) = {[xg{1},'->',xh{1}]}; % source ID -> user ID
            %             Nc_type(1,xa) = source_type(xd); % source type number
            %             Nr_type(xb,1) = xj; % user type number
            %             Nr_zone(xb,1) = xf; % zone number
            
            if ~isempty(a1)
                index = find(Zones==a1 & user_type==source_type);
            end
            
            switch source_type
                % gw source conditions
                case 1                    
                    A(xr,xa) = 1; % user inflow
                    A(xe,xa) = -Loss_percentage;
                    Nr(xe,1) = strcat(xh,'-L');
                    A(xb,xa) = -1; % source outflow
                    Nr(xb,1) = xg; % Source ID
                    
                    % storage
                    A(xb,xu) = -1; % source end of year storage
                    A(length(Nr),length(A)) = -1; % total potable storage
                    Nc(1,xu) = strcat(xg,'-strg'); % source loss ID
                    A_st(xb,xu) = 1; % next period storage carry over
                    
                    % wwtp source conditions
                case 2
                    
                    if user_type(xr) == 8  % NP main
                        
                        if a1==1
                            A(index(xn),xa) = -1; % source outflow
                        else
                            A(index,xa) = -1; % source outflow
                        end
                        
                        A(xr,xa) = 0.95; % percent of total supply to NP node (5% lost to solids)
                        A(xe,xa) = -Loss_percentage;
                        Nr(xe,1) = strcat(xh,'-L');
                        A(index(xn),xu+sourceN) = -1;% source release to surface water
                        Nc(1,xu+sourceN) = strcat(xg,'-release');
                        
                    elseif user_type(xr) == 6 % recharge
                        A(index(xn),xa) = -1; % source outflow
                        A_st(xr,xa) = 0.95; % percent of total supply to NP node (5% lost to solids)
                        A_st(xe,xa) = -Loss_percentage;
                        Nr(xe,1) = strcat(xh,'-L');
                        A(index(xn),xu+sourceN) = -1;% source release to surface water
                        Nc(1,xu+sourceN) = strcat(xg,'-release');
                        
                    else    % release
                        
                        A(index(xn),xa) = -1; % source outflow
                        A(xr,xa) = 0.95; % percent of total supply to user (5% lost to solids)
                        A(xe,xa) = -Loss_percentage;
                        Nr(xe,1) = strcat(xh,'-L');
                        A(index(xn),xu+sourceN) = -1;% source release to surface water
                        Nc(1,xu+sourceN) = strcat(xg,'-release');
                    end
                    
                    % surface water conditions
                case 3
                    if xj ~= 6 %(recharge user)
                        A(xr,xa) = 1; % user inflow from sw
                        A(xe,xa) = -Loss_percentage;
                    else
                        A_st(xr,xa) = 1;
                        A_st(xe,xa) = -Loss_percentage;
                    end
                    A(xb,xa) = -1; % source outflow
                    
                    Nr(xe,1) = strcat(xh,'-L');
                    A(xb,xu+sourceN) = -1;% downstream flow
                    Nc(1,xu+sourceN) = strcat(xg,'-DSflow');
                    Nr(xb,1) = xg;
                    
                    % wtp conditions/ reservoir
                case 4                    
                    A(xr,xa) = 1; % user inflow
                    A(xe,xa) = -Loss_percentage;
                    Nr(xe,1) = strcat(xh,'-L');
                    
                    A(index(xo),xa) = -1; % source outflow
                    
                    % dummy source conditions
                case 5
                    A(xr,xa) = 1; % user inflow from dmy source
                    
                    % recharge conditions
                case 6
                    A(xr,xa) = 1; % user inflow
                    A(xe,xa) = -Loss_percentage;
                    Nr(xe,1) = strcat(xh,'-L');
                    Nc(1,xu) = strcat(xg,'-strg');
                    
                    
                    if a1==1
                        A(index(xp),xa) = -1; % source outflow
                        A(index(xp),xu) = -1; % source end of year storage
                        A_st(index(xp),xu) = 1; % next period storage carry over
                    else
                        A(index,xa) = -1; % source outflow
                        A(index,xu) = -1; % source end of year storage
                        A_st(index,xu) = 1; % next period storage carry over
                    end
                    
                    % potable zone source
                case 7
                    
                    A(xr,xa) = 1; % user inflow from zonal main
                    A(xe,xa) = -Loss_percentage;
                    Nr(xe,1) = strcat(xh,'-L');
                    
                    A(index,xa) = -1; % source outflow
                    if xj == 10 || xj == 12
                        xi = xi+1;
                        index = find(Zones==xf & user_type==9);
                        A(index,xa) = Return_Matrix(xi,xi);
                        A(index+xk,xa) = -Loss_percentage;
                        Nr(index+xk) = strcat(userID(index),'-L');
                    end
                    
                    % non-potable zone source
                case 8
                    
                    A(xr,xa) = 1; % user inflow from zonal main
                    A(xe,xa) = -Loss_percentage;
                    Nr(xe,1) = strcat(xh,'-L');
                    A(index,xa) = -1; % source outflow
                    
                    % return zone source
                case 9
                    
                    A(xr,xa) = 1; % user inflow from zonal main
                    A(xe,xa) = -Loss_percentage;
                    Nr(xe,1) = strcat(xh,'-L');
                    A(index,xa) = -1; % source outflow
            end
        end
    end
end

%% Add total loss variables to A matrix

A( 1:userN          ,sourceN*userN+2*sourceN+1:end-1) = -eye(userN);
A((1:userN)+userN+ST,sourceN*userN+2*sourceN+1:end-1) = eye(userN);
Nc(1,sourceN*userN+2*sourceN + (1:userN)) = strcat(userID(1:userN)','-loss');
% clc;

%% Reduce A matrix to non-zero constraints
disp('Reducing A matrix to non-zero constraints...');

empties = find(cellfun(@isempty,Nc)); % identify the empty cells
Nc(empties) = [];
A(:,empties) = [];
A_st(:,empties) = [];
Ap = vertcat(A,A_st);