% MATLAB EXAMPLE WITH RANDOMIZATION TAKEN OUT

% SET UP FOR A 1D LINEAR ARRAY
% physical setup 

Nt = 64; % number of transmitting(tx) antenna
NtRF = 4; % number of rf chains tx

Nr = 16; % number of recieving(rx) antenna
NrRF = 1; % number of rf chains rx

% load H;
load txpos;


c = 3e8;
fc = 28e9;  % frequency 
lambda = c/fc;


txarray = phased.PartitionedArray(...
    'Array',phased.URA([sqrt(Nt) sqrt(Nt)],lambda/2),...
    'SubarraySelection',ones(NtRF,Nt),'SubarraySteering','Custom');

rxarray = phased.PartitionedArray(...
    'Array',phased.URA([sqrt(Nr) sqrt(Nr)],lambda/2),...
    'SubarraySelection',ones(NrRF,Nr),'SubarraySteering','Custom');

Ncl = 8;
Nray = 1;
Nscatter = Nray*Ncl;



Roptarray = zeros(9,2,10);
Rhybarray = zeros(9,2,10); 


    

%for w = 1:10
snr_param = -40:5:0;

Nsnr = numel(snr_param);

Ns_param = [1 2];
NNs = numel(Ns_param);

Ropt = zeros(Nsnr,NNs);
Rhyb = zeros(Nsnr,NNs);

Niter = 50;


w=1 ;
for m = 1:Nsnr

    snr = db2pow(snr_param(m));

    for n = 1:Niter    

        At = steervec(txpos,txang);
        Ar = steervec(rxpos,rxang);


        for k = 1:NNs

            Ns = Ns_param(k);

            % Compute optimal weights and its spectral efficiency

            [Fopt,Wopt] = helperOptimalHybridWeights(H(:,:,w),Ns,1/snr);
            Ropt(m,k) = Ropt(m,k)+helperComputeSpectralEfficiency(H(:,:,w),Fopt,Wopt,Ns,snr);

            % Compute hybrid weights and its spectral efficiency

            [Fbb,Frf,Wbb,Wrf] = helperOMPHybridWeights(H(:,:,w),NtRF,NrRF,Ns,At,Ar,1/snr);
            Rhyb(m,k) = Rhyb(m,k)+helperComputeSpectralEfficiency(H(:,:,w),Fbb*Frf,Wrf*Wbb,Ns,snr);

        end

    end

end

Ropt = Ropt/Niter;
Rhyb = Rhyb/Niter;

Roptarray(:,:,w) = Ropt;
Rhybarray(:,:,w) = Rhyb;
   
%end

% Ropt = Ropt/Niter;
% Rhyb = Rhyb/Niter;

% plot(snr_param,Ropt(:,1),'--sm', snr_param,Ropt(:,2),'--b',snr_param,Rhyb(:,1),'-k',snr_param,Rhyb(:,2),'-b');

% xlabel('SNR (dB)');
% ylabel('Spectral Efficiency (bits/s/Hz');
% 
% legend('Ns=1 optimal','Ns=2 optimal','Ns=1 hybrid', 'Ns=2 hybrid',...
%     'Location','best');
% 
% grid on;