% Faisal Baqai

%% 2 Neuron Test
maxt=1000;
N_Kick=[10,50,70,100,200,220,500]; %Presynaptic Neuron
%Vpre=zeros(1,1000);
V=zeros(1,1000);
u=V;
x=ones(1,1000);
J=.4;

tau_m=10;
tau_f=1500;
tau_d=200;
U=.2;
dt=1;
for i=2:1000
    
    if ismember (i,N_Kick)  % synaptic kick?
        kick=1;
    else    
        kick=0;
    end
    %I_ext=normrnd(.1,.2);
    I_rec=J*u(i-1)*x(i-1)*kick;%.05*kick;
    V(i)=V(i-1)+(dt/tau_m)*(-V(i-1)+I_rec);
    if V(i)>1
        V(i)=0;
    end
    u(i)=u(i-1)+dt*((U-u(i-1))/tau_f+U*(1-u(i-1))*kick);
    x(i)=x(i-1)+dt*((1-x(i-1))/tau_d-u(i-1)*x(i-1)*kick);
    if x(i)<0
        x(i)=0;
    end
end

%x=(x-min(x))/(max(x)-min(x));  % X normalization which makes things weird
figure
subplot(2,1,1)
plot(u);
subplot(2,1,2)
plot(x)

figure
subplot(2,1,1)
Raster(maxt,{N_Kick})
subplot(2,1,2)
plot(V)



%% Network architecture 

%D=rand(1000)*4+1;


% for i=1:1000
%     
%     if (i==N_kick)  % synaptic kick?
%         kick=1;
%     else    
%         kick=0;
%     end
%     
%     I_rec=
%     V(i)=V(i-1)+(dt/tau_m)*(V(i-1)+I_rec+I_ext);
%     u(i)=u(i-1)+dt*((U-u(i-1))/tau_f+U*(1-u(i-1))*kick);
%     x(i)=
% end
% maxt=1200;
% [ PopTrains,S] = LIF( maxt,1000,1,0,.08,0.25);