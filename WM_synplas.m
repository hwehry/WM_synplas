

clear all;
tic

% simulation params
dt=0.01;
T = 1e3;
steps = T/dt;

% Single-cell parameters
Vt=20;  %threshold
Vre = 16;
Vri = 13;
tau_mE = 15; %ms for E
tau_mI = 10; %10 ms for I
tau_arp = 2; %ms for both E & I

% Network parameters
f = 0.10; % coding level
p = 5; % # of memories
c = 0.2; %connection prob
Ne=8000; %8000 excitatory
Ni=2000;  %250 inhibitory
Je0=23.10;%23.10; % mV mean external drives
Ji0=21.0;%21.0; %mV must add random noise to both
J0std = sqrt(1);

% Synaptic parameters
Jie=0.135;
Jei=0.25; %
Jii=0.20; %
Jb=0.10; % baseline level of EE strenght
Jp=0.45; %potentiated level of EE strength
gamma0=0.10; %fraction of potentiated synapses before learning
delay=0.1;% 0.1-1ms

%short-term synaptic dynamics params
U = 0.20; %baseline utilization factor
tau_f= 1500; %recovery time of utilization factor
tau_d= 200; %recovery time of synaptic resources

%selective stimulation
Tcue=350; %duration
Acue=1.15; %contrast factor

% reactivating signals
Treact = 250;
Areact = 1.05;
Tperiodic = 100; %duration of periodic reactivating signal
Pperiodic = 250; %period of periodic reactivating signal
Aperiodic = 1.075; %contrast factor
nonspecificInputs = rand(Ne,1) < 0.15;

%connection matrices
% selective E-E within the same selective population: Jp
% selective E-E between different selective populations: Jb
% non-selective E-E to selective and non-selective to non-selevective: Jp
% wp 0.1

%E to E
cEE=zeros(Ne);
cEE_t=rand(Ne);
cEE(cEE_t< c)=1;

% E to E J matrix
% Es : each selective population
% En : non-selective pop
EEs = f*Ne; % selective
EEn = Ne-p*EEs; % # non-selective
Jnn = rand(EEn);
Jnn(Jnn > gamma0) = Jb;
Jnn(Jnn < gamma0) = Jp;
Jns = ones(EEs*p,EEn).*Jb;
Jsn = Jns';
% selective EE
Jsmemcell = cellfun(@double,repmat({Jp*ones(EEs,EEs)},1,p),'Un',0);
Jsmem = blkdiag(Jsmemcell{:});
Jss = Jb.*(Jsmem<Jp)+Jsmem;
Jee = blkdiag(Jss,Jnn)+fliplr(blkdiag(Jns,Jsn));

%E to I
cIE=zeros(Ni,Ne);
cIE_t=rand(Ni,Ne);
cIE(cIE_t<c)=1;

%I to E
cEI=zeros(Ne,Ni);
cEI_t=rand(Ne,Ni);
cEI(cEI_t<c)=1;

%I to I
cII=zeros(Ni);
cII_t=rand(Ni);
cII(cII_t<c)=1;

% intial conditions
Ve=Vre + (Vt-Vre).*rand(Ne,1);
Vi=Vri + (Vt-Vri).*rand(Ni,1);
u = U*ones(Ne,1);
x = ones(Ne,1);

Dee = round((rand(Ne)*4+1)*100,0); %round to .01
Die = round((rand(Ni,Ne)*4+1)*100,0); %round to .01
Dei = round((rand(Ne,Ni)*4+1)*100,0); %round to .01
Dii = round((rand(Ni)*4+1)*100,0); %round to .01

delayidx = 5/.01+1;
udelay = zeros(Ne,delayidx);
udelay(:,1) = u;
xdelay = zeros(Ne,delayidx);
xdelay(:,1) = x;
edelay = zeros(Ne,delayidx);
idelay = zeros(Ni,delayidx);

% spiketime arrays
maxspk=100000;
spktime_e=zeros(maxspk,1);
spkindex_e=zeros(maxspk,1);
spktime_i=zeros(maxspk,1);
spkindex_i=zeros(maxspk,1);

% refractory period
steps_refrac = tau_arp/dt;
resetpenaltye = false(Ne,steps+steps_refrac);
resetpenaltyi = false(Ni,steps+steps_refrac);

counte=1;
counti=1;

storeu = zeros(steps,1);
storex = zeros(steps,1);
storev = zeros(steps,1);
storev(1,1) = Ve(1);

% time loop
for t=[dt:dt:T]
    idx=floor(t/dt);
    
    edelay(:,2:delayidx) = edelay(:,1:delayidx-1);
    idelay(:,2:delayidx) = idelay(:,1:delayidx-1);
    
    spke = Ve>=Vt;
    index_spke=find(spke);  %find spikers
    
    if (~isempty(index_spke))
        spktime_e(counte:counte+length(index_spke)-1)=t; %update arrays
        spkindex_e(counte:counte+length(index_spke)-1)=index_spke; %update arrays
        Ve(index_spke)=Vre;  %reset
        counte=counte+length(index_spke)+1;
    end
    
    spki = Vi>=Vt;
    index_spki=find(spki);
    
    if (~isempty(index_spki))
        spktime_i(counti:counti+length(index_spki)-1)=t;
        spkindex_i(counti:counti+length(index_spki)-1)=index_spki;
        Vi(index_spki)=Vri;
        counti=counti+length(index_spki)+1;
    end
    
    %update refractory penalty box
    timeout = idx+steps_refrac;
    %or use idx+1 and timeout
    resetpenaltye(:,idx:timeout-1) = repmat(spke,1,steps_refrac) | resetpenaltye(:,idx:timeout-1) ;
    resetpenaltyi(:,idx:timeout-1) = repmat(spki,1,steps_refrac) | resetpenaltyi(:,idx:timeout-1) ;
    
    
    edelay(:,1) = spke;
    idelay(:,1) = spki;
    
    
    
    % %     if T <= 350
    % vector of presynaptic 'calcium' and 'neurotransmitter'
    u = udelay(:,1);     
    udelay(:,2:delayidx) = udelay(:,1:delayidx-1);
    udelay(:,1) = u+dt/tau_f.*(U-u)+dt*U.*(1-u).*spke;
    
    x = xdelay(:,1);
    xdelay(:,2:delayidx) = xdelay(:,1:delayidx-1);
    xdelay(:,1) = x+dt/tau_d.*(1-x)-dt*u.*x.*spke;
    xdelay(xdelay(:,1)<0,1) = 0; %flatten to zero
    
    % E Cells Presynaptic
    uee = false(Ne,Ne);
    xee = false(Ne,Ne);
    kee = false(Ne,Ne);
    kie = false(Ni,Ne);
    %presynaptic neuron j
    for j = 1:Ne
        % postsynaptic neuron i
        for i = 1:Ne
            uee(i,j) = udelay(j,Dee(i,j)+1); % because zero delay is idx 1
            xee(i,j) = xdelay(j,Dee(i,j)+1);
            kee(i,j) = edelay(j,Dee(i,j)+1);
        end
        for i = 1:Ni
            kie(i,j) = edelay(j,Die(i,j)+1);
        end
    end
    
    %I Cells Presynaptic
    %presynaptic neuron j
    kei = false(Ne,Ni);
    kii = false(Ni,Ni);
    for j = 1:Ni
        % postsynaptic neuron i
        for i = 1:Ne
            kei(i,j) = idelay(j,Dei(i,j)+1);
        end
        for i = 1:Ni
            kii(i,j) = idelay(j,Dii(i,j)+1);
        end
    end
    
    Jhat = Jee.*uee.*xee.*cEE;
    Irecee = sum(Jhat.*kee,2);
    Irecie = sum(Jie.*cIE.*kie,2);
    Irecei = sum(Jei.*cEI.*kei,2);
    Irecii = sum(Jii.*cII.*kii,2);
    
    
    storeu(idx) = u(1);
    storex(idx) = x(1);
    
    
    Ve=Ve+(Irecee-Irecei).*~resetpenaltye(:,idx);  %spike interaction
    Vi=Vi+(Irecie-Irecii).*~resetpenaltyi(:,idx);
    
    Ve=Ve+dt/tau_mE*(-Ve+Je0+J0std*rand(Ne,1)).*~resetpenaltye(:,idx);  %E membrane integration
    Vi=Vi+dt/tau_mI*(-Vi+Ji0+J0std*rand(Ni,1)).*~resetpenaltyi(:,idx);  %I membrane integration
    
    storev(idx) = Ve(1);
    
    % selective stimulation
    % must add the reset here because Vre ~= 0 

    if t<350,
        Ancue = Acue.*~resetpenaltye(:,idx);
        Ancue(Ancue<1) = 1;
        Ve(1:EEs*1) = Ve(1:EEs*1).*Ancue(1:EEs*1);
    end
    % periodic reactivating signal
    if ( (t > 600 && t < 700) || (t > 850 && t < 950) )
        Anperiodic = Aperiodic.*~resetpenaltye(:,idx);
        Anperiodic(Anperiodic<1) = 1;
        Ve(nonspecificInputs) = Ve(nonspecificInputs).*Anperiodic;
    end
    
    if mod(100*t,10) == 10, progress = t/T; disp(progress); end 
end

toc

figure;

subplot(2,1,1), plot(spktime_e,spkindex_e,'.k', 'MarkerSize',8); xlabel('Time (ms)', 'fontsize', 16, 'fontweight', 'b'); ylabel('E cell index', 'fontsize', 16, 'fontweight', 'b')
subplot(2,1,2), plot(spktime_i,spkindex_i,'.k', 'MarkerSize',8); xlabel('Time (ms)', 'fontsize', 16, 'fontweight', 'b'); ylabel('I cell index', 'fontsize', 16, 'fontweight', 'b')






