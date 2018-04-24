

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
tau_mI = 5; %10 ms for I
tau_arp = 2; %ms for both E & I

% Network parameters
f = 0.10; % coding level
p = 5; % # of memories
c = 0.2; %connection prob
Ne=1000; %8000 excitatory
Ni=250;  %250 inhibitory
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



%connection matrices
% selective E-E within the same selective population: Jp
% selective E-E between different selective populations: Jb
% non-selective E-E to selective and non-selective to non-selevective: Jp
% wp 0.1

%E to E
cEE=zeros(Ne);
cEE_t=rand(Ne);
cEE(cEE_t< p)=1;

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
cIE=zeros(Ne,Ni);
cIE_t=rand(Ne,Ni);
cIE(cIE_t<p)=1;

%I to E
cEI=zeros(Ni,Ne);
cEI_t=rand(Ni,Ne);
cEI(cEI_t<p)=1;

%I to I
cII=zeros(Ni);
cII_t=rand(Ni);
cII(cII_t<p)=1;

% intial conditions
Ve=Vre + (Vt-Vre).*rand(Ne,1);
Vi=Vri + (Vt-Vri).*rand(Ni,1);
u = zeros(Ne,1);
x = zeros(Ne,1);

spke=zeros(Ne,1);
kickei=zeros(Ne,1);
kickie=zeros(Ni,1);
kickii=zeros(Ni,1);
    
D = round((rand(Ne+Ni)*4+1)*100,0); %round to .01
delayidx = 5/.01+1;
udelay = zeros(Ne,delayidx);
xdelay = zeros(Ne,delayidx);
eedelay = zeros(Ne,delayidx);
eidelay = zeros(Ne,delayidx);
iedelay = zeros(Ni,delayidx);
iidelay = zeros(Ni,delayidx);

% spiketime arrays
maxspk=100000;
spktime_e=zeros(maxspk,1);
spkindex_e=zeros(maxspk,1);
spktime_i=zeros(maxspk,1);
spkindex_i=zeros(maxspk,1);

counte=1;
counti=1;

storeu = zeros(steps,1);
storex = zeros(steps,1);
storev = zeros(steps,1);




% time loop
for t=[dt:dt:T]
    
% zero the interactions from the last step

    
    % work out delay
    spke = eedelay(:,1);
    kickei = eidelay(:,1);
    kickie = iedelay(:,1);
    kickii = iidelay(:,1);
    
    eedelay(:,2:delayidx) = eedelay(:,1:delayidx-1);
    eidelay(:,2:delayidx) = eidelay(:,1:delayidx-1);
    iedelay(:,2:delayidx) = iedelay(:,1:delayidx-1);
    iidelay(:,2:delayidx) = iidelay(:,1:delayidx-1);

    
    index_spke=find(Ve>=Vt);  %find spikers
    spke = (Ve>=Vt);
    
    if (~isempty(index_spke))
    spktime_e(counte:counte+length(index_spke)-1)=t; %update arrays
    spkindex_e(counte:counte+length(index_spke)-1)=index_spke; %update arrays
    Ve(index_spke)=Vre;  %reset
    counte=counte+length(index_spke)+1;
    end
    
    index_spki=find(Vi>=Vt);
    
    if (~isempty(index_spki))
    spktime_i(counti:counti+length(index_spki)-1)=t;
    spkindex_i(counti:counti+length(index_spki)-1)=index_spki;
    Vi(index_spki)=Vri;
    counti=counti+length(index_spki)+1;
    end
    
    for j=1:length(index_spke) %update kick arrays
%      
%         
%         
%         %need to edit this part because it matters which e gives what spike
%        kickee_index=find(cEE(index_spke(j),:)>0);
%        kickee(kickee_index)=kickee(kickee_index)+1;
%        
       kickie_index=find(cIE(index_spke(j),:)>0);
       kickie(kickie_index)=kickie(kickie_index)+1;         
   end    
    
    for j=1:length(index_spki)
     
       kickei_index=find(cEI(index_spki(j),:)>0);
       kickei(kickei_index)=kickei(kickei_index)+1;
       
       kickii_index=find(cII(index_spki(j),:)>0);
       kickii(kickii_index)=kickii(kickii_index)+1;
        
    end 

    % work out delay
    eedelay(:,1) = spke; 
    eidelay(:,1) = kickei;
    iedelay(:,1) = kickie;
    iidelay(:,1) = kickii;
    
% %     if T <= 350
    % vector of presynaptic 'calcium' and 'neurotransmitter'
    u = udelay(:,1);
    udelay(:,2:delayidx) = udelay(:,1:delayidx-1);
    udelay(:,1) = u+dt/tau_f.*(U-u)+dt*U.*(1-u).*spke;
    
    x = xdelay(:,1);
    xdelay(:,2:delayidx) = xdelay(:,1:delayidx-1);
    xdelay(:,1) = x+dt/tau_d.*(1-x)-dt*u.*x.*spke;
    xdelay(xdelay(:,1)<0,1) = 0; %flatten to zero
    
    %presynaptic neuron j
    for j = 1:Ne
        % postsynaptic neuron i
        for i = 1:Ne
            uee(i,j) = udelay(j,D(i,j)); % because zero delay is idx 1
            xee(i,j) = xdelay(j,D(i,j));
            kee(i,j) = eedelay(j,D(i,j));
        end
    end
    Jhat = Jee.*uee.*xee.*cEE;
    Irec = sum(Jhat*kee,2);
    
    %column is pre
%     idx=floor(t/dt);
%     storeu(idx) = u(1);
%     storex(idx) = x(1);

       
    Ve=Ve+Jee*spke-Jei*kickei;  %spike interaction 
    Vi=Vi+Jie*kickie-Jii*kickii;  
    
    Ve=Ve+dt/tau_mE*(-Ve+Je0+J0std*rand(Ne,1));  %E membrane integration 
    Vi=Vi+dt/tau_mI*(-Vi+Ji0+J0std*rand(Ni,1));  %I membrane integration
   
    storev(idx) = Ve(1);

end 

toc

figure;

subplot(2,1,1), plot(spktime_e,spkindex_e,'.k', 'MarkerSize',8); xlabel('Time (ms)', 'fontsize', 16, 'fontweight', 'b'); ylabel('E cell index', 'fontsize', 16, 'fontweight', 'b')
subplot(2,1,2), plot(spktime_i,spkindex_i,'.k', 'MarkerSize',8); xlabel('Time (ms)', 'fontsize', 16, 'fontweight', 'b'); ylabel('I cell index', 'fontsize', 16, 'fontweight', 'b')






