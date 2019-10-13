% g = (randn(1,Nscatter)+1i*randn(1,Nscatter))/sqrt(Nscatter); random
% g is the channel gains/path gains, 1xNs (number of scatteres =
% clusters*rays)


Ncl = 2;
Nray = 4;
Nscatter = Nray*Ncl;

%rng(7);
g = randn(1000,8)+1i*randn(1000,8);
size(g)

txclang = [-13.0458,166.9542];
rxclang= [-25.4037,25.4037];

Nt = 64; % number of transmitting(tx) antenna
NtRF = 4; % number of rf chains tx
 
Nr = 16; % number of recieving(rx) antenna
NrRF = 4; % number of rf chains rx

% rng(1); % set seed for reproducibility if random
c = 3e8;
fc = 28e9;  % frequency 
lambda = c/fc;


txarray = phased.PartitionedArray(...
    'Array',phased.URA([sqrt(Nt) sqrt(Nt)],lambda/2),...
    'SubarraySelection',ones(NtRF,Nt),'SubarraySteering','Custom');
rxarray = phased.PartitionedArray(...
    'Array',phased.URA([sqrt(Nr) sqrt(Nr)],lambda/2),...
    'SubarraySelection',ones(NrRF,Nr),'SubarraySteering','Custom');

txpos = getElementPosition(txarray)/lambda;
rxpos = getElementPosition(rxarray)/lambda;
txang = [6, 4, -18, -8, 26, 29, 1, -4];
rxang = [43, -73, 78, 10, -6, -5, -40, 0];

H = scatteringchanmtx(txpos,rxpos,txang,rxang,g(2,:));
Hnext = scatteringchanmtx(txpos,rxpos,txang,rxang,g(3,:));

figure
F = diagbfweights(H);
F = F(1:NtRF,:);
pattern(txarray,fc,-90:90,-90:90,'Type','efield',...
    'ElementWeights',F','PropagationSpeed',c);

At = steervec(txpos,txang);
Ar = steervec(rxpos,rxang);

Ns = NtRF;
[Fbb,Frf] = helperOMPHybridPrecodingWeights(H,NtRF,Ns,At);

figure
pattern(txarray,fc,-90:90,-90:90,'Type','efield',...
    'ElementWeights',Frf'*Fbb','PropagationSpeed',c);



figure
Fnext = diagbfweights(Hnext);
Fnext = Fnext(1:NtRF,:);
pattern(txarray,fc,-90:90,-90:90,'Type','efield',...
    'ElementWeights',Fnext','PropagationSpeed',c);

At = steervec(txpos,txang);
Ar = steervec(rxpos,rxang);

Ns = NtRF;
[Fbb,Frf] = helperOMPHybridPrecodingWeights(Hnext,NtRF,Ns,At);

figure
pattern(txarray,fc,-90:90,-90:90,'Type','efield',...
    'ElementWeights',Frf'*Fbb','PropagationSpeed',c);