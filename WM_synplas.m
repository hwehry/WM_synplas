% MATH 3370 final project: synaptic theory of working memory
% recurrent LIF network with modeled calcium and neurotransmitter dynamics
% Hillary Wehry, May 1, 2018


clear all;

% simulation params
dt=.01;
T = 3e3;
steps = T/dt;

% Single-cell parameters
Vt = 20;  %threshold
Vre = 16;
Vri = 13;
tau_mE = 15; %ms for E
tau_mI = 10; %10 ms for I
tau_arp = 2; %ms for both E & I

% Network parameters
f = 0.10; % coding level
p = 1; % # of memories
c = 0.2; %connection prob
Ne=8000; %8000 excitatory
Ni=2000;  %2000 inhibitory
Je0=23.1*ones(Ne,1);%20.5 %23.10; % mV mean external drives
Ji0=21.0*ones(Ni,1);%20 %21.0; %mV must add random noise to both
J0std = sqrt(1);%sqrt(1);

% Synaptic parameters
Jie=0.135; % e --> i
Jei=0.25;%0.25; %  i --> e
Jii=0.20; %  i --> i
Jb=0.10; % baseline level of EE strength
Jp=0.45; %potentiated level of EE strength
gamma0=0.10; %fraction of potentiated synapses "before learning"
delay=0.1;% synaptic delays: 0.1-1ms

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
%E to E
cEE=zeros(Ne);
cEE_r=rand(Ne);
cEE(cEE_r< c)=1;

%E to I
cIE=zeros(Ni,Ne);
cIE_r=rand(Ni,Ne);
cIE(cIE_r<c)=1;

%I to E
cEI=zeros(Ne,Ni);
cEI_r=rand(Ne,Ni);
cEI(cEI_r<c)=1;

%I to I
cII=zeros(Ni);
cII_r=rand(Ni);
cII(cII_r<c)=1;
clear cEE_r cIE_r cEI_r cII_r

% J matrix for E-E cells

% selective E-E within the same selective population: Jp
% selective E-E between different selective populations: Jb
% non-selective E-E to selective and non-selective to non-selevective: Jp
% w.p. 0.1, o.w. Jb

%code is written for one memory
% E-E J matrix
EEs = f*Ne; % each selective population
EEn = Ne-p*EEs; % # non-selective population
Jr = rand(EEn);
Jnn = zeros(EEn);
Jnn(Jr > gamma0) = Jb;
Jnn(Jr < gamma0) = Jp;
Jns = Jb.*ones(EEn,EEs*p);
Jr = rand(EEs*p,EEn);
Jsn = zeros(EEs*p,EEn);
Jsn(Jr > gamma0) = Jb;
Jsn(Jr < gamma0) = Jp;
% selective EE
Jsmemcell = cellfun(@double,repmat({Jp*ones(EEs,EEs)},1,p),'Un',0);
Jsmem = blkdiag(Jsmemcell{:});
Jss = Jb.*(Jsmem<Jp)+Jsmem;
Jee = blkdiag(Jss,Jnn)+fliplr(blkdiag(Jsn,Jns));


% intial conditions
Ve=Vre + (Vt-Vre).*rand(Ne,1);
Vi=Vri + (Vt-Vri).*rand(Ni,1);
u = U*ones(Ne,1);
x = ones(Ne,1);

% delay matrices 1-5ms
Dee = round((rand(Ne)*4+1)*1/dt,0); %round to .01
Die = round((rand(Ni,Ne)*4+1)*1/dt,0); %round to .01
Dei = round((rand(Ne,Ni)*4+1)*1/dt,0); %round to .01
Dii = round((rand(Ni)*4+1)*1/dt,0); %round to .01

delayidx = 5/dt+1;
udelay = ones(Ne,delayidx)*U;
% udelay(:,1) = u;
xdelay = ones(Ne,delayidx);
% xdelay(:,1) = x;
edelay = zeros(Ne,delayidx);
idelay = zeros(Ni,delayidx);

% spiketime arrays
maxspike=100000;
spiketime_e=zeros(maxspike,1);
spikeindex_e=zeros(maxspike,1);
spiketime_i=zeros(maxspike,1);
spikeindex_i=zeros(maxspike,1);

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
idx = 1;
for t=[dt:dt:T]
    
    edelay(:,2:delayidx) = edelay(:,1:delayidx-1); %shift delayed spike train
    idelay(:,2:delayidx) = idelay(:,1:delayidx-1);
    
    spikee = Ve>=Vt;
    index_spike=find(spikee);  %find spiking neurons
    
    if (~isempty(index_spike))
        spiketime_e(counte:counte+length(index_spike)-1)=t; %update
        spikeindex_e(counte:counte+length(index_spike)-1)=index_spike; %update
        Ve(index_spike)=Vre;  %reset
        counte=counte+length(index_spike)+1;
    end
    
    spikei = Vi>=Vt;
    index_spikei=find(spikei);
    
    if (~isempty(index_spikei))
        spiketime_i(counti:counti+length(index_spikei)-1)=t; %update
        spikeindex_i(counti:counti+length(index_spikei)-1)=index_spikei; %update
        Vi(index_spikei)=Vri; %reset
        counti=counti+length(index_spikei)+1;
    end
    
    %update refractory penalty box
    timeout = idx+steps_refrac;
    %or use idx+1 and timeout
    resetpenaltye(:,idx:timeout-1) = repmat(spikee,1,steps_refrac) | resetpenaltye(:,idx:timeout-1) ;
    resetpenaltyi(:,idx:timeout-1) = repmat(spikei,1,steps_refrac) | resetpenaltyi(:,idx:timeout-1) ;
    
    
    edelay(:,1) = spikee; 
    idelay(:,1) = spikei;
    
    % vector of presynaptic 'calcium' and 'neurotransmitter'
    % the functions multiplying the spike train must be evaulated first
    u = udelay(:,1);
    u = u +(dt/tau_f).*(U-u);
    udelay(:,2:delayidx) = udelay(:,1:delayidx-1);
    udelay(:,1) =  u + U.*(1-u).*spikee;
    
    x = xdelay(:,1);
    x = x + (dt/tau_d).*(1-x);
    xdelay(:,2:delayidx) = xdelay(:,1:delayidx-1);
    xdelay(:,1) = x - u.*x.*spikee;
    xdelay(xdelay(:,1)<0,1) = 0; %flatten to zero

    
    % E Cells Presynaptic
    uee = zeros(Ne,Ne);
    xee = zeros(Ne,Ne);
    kee = zeros(Ne,Ne);
    kie = zeros(Ni,Ne);
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
    kei = zeros(Ne,Ni);
    kii = zeros(Ni,Ni);
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
    
    Ve=Ve+(Irecee-Irecei).*~resetpenaltye(:,idx);  %spike interaction
    Vi=Vi+(Irecie-Irecii).*~resetpenaltyi(:,idx);
    
    %external current with time-dependent stimulation
    Iexte = Je0+J0std*randn(Ne,1);
    Iexti = Ji0+J0std*randn(Ni,1);
    if  t>1000 && t<1000+350
        Iexte =   Iexte+Acue.*[ones(EEs*1,1);zeros(Ne-EEs*1,1)]; %Iexte.*[Acue*ones(EEs*1,1);ones(Ne-EEs*1,1)];
    elseif ( (t > 2200 && t < 2300) || (t > 2550 && t < 2650) )
        Iexte =   Iexte+Aperiodic.*ones(Ne,1); %Iexte.*Aperiodic.*ones(Ne,1); 
    else
        Iexte = Iexte;
    end
    
    Ve=Ve+dt/tau_mE*(-Ve+Iexte).*~resetpenaltye(:,idx);  %E membrane integration
    Vi=Vi+dt/tau_mI*(-Vi+Iexti).*~resetpenaltyi(:,idx);  %I membrane integration
    
    % save variables
    storev(idx) = Ve(1);
    storeu(idx) = udelay(1,1);%mean(udelay(1:EEs,1));
    storex(idx) = xdelay(1,1);%mean(xdelay(1:EEs,1));
        
    idx = idx+1;
end



figure;

subplot(2,1,1), plot(spiketime_e,spikeindex_e,'.k', 'MarkerSize',8); xlabel('Time (ms)', 'fontsize', 16, 'fontweight', 'b'); ylabel('E cell index', 'fontsize', 16, 'fontweight', 'b')
subplot(2,1,2), plot(spiketime_i,spikeindex_i,'.k', 'MarkerSize',8); xlabel('Time (ms)', 'fontsize', 16, 'fontweight', 'b'); ylabel('I cell index', 'fontsize', 16, 'fontweight', 'b')

figure;
plot(storex(1:idx-1))
hold on
plot(storeu(1:idx-1))

figure;
plot(storev(1:idx-1))
% savefile = 'WM_synplasrun_10000_add';
% save(savefile,'spiketime_e','spikeindex_e','spiketime_i','spikeindex_i','storex','storeu');
