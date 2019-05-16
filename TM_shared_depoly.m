close all
clear all
% two flagella diffusive- motors and tubulin shared with active disassembly

% uncommenting a section below would implement the feedback mechanism. The
% section to uncomment starts with 'begin feedback' and ends with
% 'end feedback' comment line

% to implement severing uncomment section starting with comment 'begin severing' and
% ending with 'end severing'


% input: starts with comment 'parameters of assembly process' and ends with 'input of user ends'
% output: plot of length of both the flagella(len1, len2) vs time

% variable motor stores properties of the motors
% first row stores if in base or in flagella 1 or in flagella 2
% second one stores if ballistic or diffusive, 
% third one stores the position of the motor in the flagella and 0 if it is
% in basal body
% fourth row stores the amount of tubulin which the motor carries

%boolean values describing various states of the motors
% basal=0;%if motors are in basal body
% flagella=1;%if motors are in the flagellum 1
% flagella=2;%if motors are in the flagellum 2
% ballistic=1;%undergoing ballistic motion
% diffuse=2;%undergoing diffusing motion

tic

%%(input by user)
% paramters of assembly process
n=127;  % amount of time to simulate the process in minutes
D=1.7; % Diffusion constant (microm^2/s)
T=84; % total tubulin(microm)
gamma=2.5*10^-4; % constant gamma
gammakM = 0.2/60; %value of gamma * kon *M
da = 2/60; % Disassembly rate(um/s)
tot = 3/60;
v=2.5; % velocity of ballistic path(microm/s)
Lss = 10.5; % length of steady state of flagella
t_m = 10*60; % replenishment timescale for motors in seconds
t_t = 10*60; % replenishment timescale for tubulin in seconds
sevt=20; % time at which flagella is severed in minutes
dt=0.032; %time step dt in seconds
% input of user ends

% calculating other parameters from input parameters
kon = (Lss/v+Lss^2/(2*D))^(-1)*((0.5/tot)*gammakM*(T-2*Lss)-1); %value of k
gammaM = gammakM/kon; 
M=ceil(gammaM/gamma); %number of motors
Jfac = 1/gamma;
Jss = Jfac*0.5*gammakM/(1+kon*(Lss/v+Lss^2/(2*D)));
d1 = D/(Lss*Jss)*(tot-da);
mydelta = da+d1*(Lss*Jss/(D));
Nt = T-2*Lss;  % amount of free tubulin in the basal pool at steady state
gammaJ = 0.5*gammakM/(1+kon*Lss/v+0.5*kon*Lss^2/D);
gammaMd = gammaJ*Lss^2/(D);
gammaMb = 2*Lss*gammaJ/v;
Nm = (gammaM-gammaMb-gammaMd)/gamma; % amount of free motors in the basal pool at steady state

%program parameters
dx=v*dt; %distance moved per time step dt
time_total=0:dt:100*60*n/100; %in seconds
N=3000000/4*n/100*0.008/dt; % N=n/dt where n is in time of simulation

%parameters updated at each time step in simulation
motor=zeros(4, M); % stores properties of the motors
% first row stores if in base or flagella 1 or flagella 2
% second one stores if ballistic or diffusive, 
% third one stores the position of the motor in the flagella and 0 if it is
% in basal body
% fourth row stores the amount of tubulin which the motor carries
len1 = ones(1, N+1)*Lss; %initialise and store length vs time
len2 = ones(1, N+1)*Lss; %initialise and store length vs time
freeM=M; % number of free motors in basal pool 
freeT=T-len1(1)-len2(1); % amount of free tubulin in basal pool
time=0; % present discrete time
check = 0.1; %checks fraction of simulation completed

while time<N
    time=time+1;
    
    %prints the fraction of job completed
    if (time/N > check)
        'fraction done is'
        time/N
        check = check+0.1;
    end
    
    
    % checks if the time = sevt minutes, the flagella is severed at that point
    % to have 0 length
    
    %%begin severing
%     if time==3000000/4*sevt/100*0.008/dt
%         sever_amt = len1(time);
%         len1(time)=len1(time)-sever_amt;
%         find_b2 = find(motor(1,:)==1 & motor(3,:)>(len1(time)-sever_amt));
%         if ~isempty(find_b2)
%             motor(1:3, find_b2) = -1; %finding motors which are in ballistic and move them by a constant amount
%         end
%     end
    %%end severing
    
    %%(uncomment below to implement feedback)
    
    %%begin feedback
%     freeT = freeT+(Nt-freeT)/t_t*dt; %feedback for tubulins
%  
%     % feedback system for free motors greater than steady state in basal pool
% 
%     prob_out=dt*(freeM-Nm)/t_m;
%     if rand<prob_out
%         find_f_o = find(motor(1,:)==0);
%         if ~isempty(find_f_o)
%             change_f_o = randi(length(find_f_o));
%             motor(1:4, find_f_o(change_f_o))=-1;
%             freeM=freeM-1;
%         end
%     end
%     
%     % feedback system for free motors less than steady state in basal pool
% 
%     prob_in=dt*(Nm-freeM)/t_m;
%     if rand<prob_in
%         find_f = find(motor(1,:)==-1);
%         if ~isempty(find_f)
%             change_f = randi(length(find_f));
%             motor(1:4, find_f(change_f))=0;
%             freeM=freeM+1;
%         else
%             motor = cat(2, zeros(size(motor,1),1), motor);
%             freeM=freeM+1;
%         end
%     end
    %%end feedback
    
    %%uncomment above to implement feedback
    
    % disassembly proces in both the flagella depending on
    % concentration(active disassembly) ensuring that
    % lengths are non-negative
    if (len1(time) <= da*dt)
        freeT=freeT+len1(time);
        len1(time+1)=0;
    else
        f1 = find(motor(2,:)==2 & motor(1,:)==1 & motor(3,:)>=len1(time)-1 & motor(3,:)<=len1(time));
        len1(time+1)=len1(time)-da*dt-length(f1)*d1*dt;
        freeT=freeT+da*dt+length(f1)*d1*dt;
    end
    if (len2(time) <= da*dt)
        freeT=freeT+len2(time);
        len2(time+1)=0;
    else
        f2 = find(motor(2,:)==2 & motor(1,:)==2 & motor(3,:)>=len2(time)-1 & motor(3,:)<=len2(time));
        len2(time+1)=len2(time)-da*dt-length(f2)*d1*dt;
        freeT=freeT+da*dt+length(f2)*d1*dt;
    end
    
    % choosing which flagella to be injected into and injecting a
    % motor into the flagella with probability prob_motor
    
    if rand<0.5
        %flagella 1 motor injection
        prob_motor = kon*freeM*dt; %probability that motor is injected
        if rand<=prob_motor
            %randomly choosing a motor which is in basal body to move
            find_motor = find(~motor(1,:));
            change = randi(length(find_motor));
            mot = find_motor(change); %motor chosen to be inserted
            motor(1:2,mot) = 1; % motor number in mot inserted into flagellum
            motor(4,mot) = gamma*freeT; %amount of tubulin the motor carries
            freeT = freeT-gamma*freeT; 
            freeM=freeM-1; % decrease the free motors in basal pool by 1
        end
    else
         % flagella 2 motor injection
        prob_motor = kon*freeM*dt; %probability that motor is injected
        if rand<=prob_motor
            %randomly choosing a motor which is in basal body to move
            find_motor = find(~motor(1,:));
            change = randi(length(find_motor));
            mot = find_motor(change); %motor chosen to be inserted
            motor(1,mot) = 2; % motor number in mot inserted into flagellum
            motor(2,mot) = 1;
            motor(4,mot) = gamma*freeT; %amount of tubulin the motor carries
            freeT = freeT-gamma*freeT; 
            freeM=freeM-1; % decrease the free motors in basal pool by 1
        end
    end
    
    find_b = find(motor(2,:)==1);
    if ~isempty(find_b)
        motor(3, find_b) = motor(3, find_b) + dx; %finding motors which are in ballistic and move them by a constant amount
    end
    
    find_d = find(motor(2,:)==2); %finding motors which are diffusive
    if ~isempty(find_d)
        d=sqrt(2*D*dt)*randn(1,length(find_d));
        motor(3, find_d) = motor(3, find_d) + d;
    end   
    
    find_t = find(motor(3,:)>=len1(time) & motor(1,:)==1 & motor(2, :)==1); %finding motors which have reached the tip of flagella 1 and are ballistic
    if ~isempty(find_t) % in flagella 1
        for i=1:length(find_t)
            len1(time+1)=len1(time+1)+ motor(4, find_t(i));  %adding tubulin pool to the flagella 1 from motor
        end
        motor(2,find_t) = 2;
        motor(3,find_t) = len1(time+1);
    end
    
    find_td= find(motor(3,:)>=len1(time+1) & motor(2, :)==2 & motor(1,:)==1); %finding motors which have reached the tip of flagella 1 and are diffusive
    if ~isempty(find_td) % in flagella 1
        motor(3,find_td) =  2*len1(time+1) - motor(3,find_td); % reflecting boundary condition
    end
    
    find_t1 = find(motor(3,:)>=len2(time) & motor(2, :)==1 & motor(1,:)==2); %finding motors which have reached the tip of flagella 2 and are ballistic
    if ~isempty(find_t1) % in flagella 2
        for i=1:length(find_t1)
            len2(time+1)=len2(time+1)+ motor(4, find_t1(i)); %adding tubulin pool to the flagella 2 from motor
        end
        motor(2,find_t1) = 2;
        motor(3,find_t1) = len2(time+1);
    end
    
    find_td1= find(motor(3,:)>=len2(time+1) & motor(2, :)==2 & motor(1,:)==2); %finding motors which have reached the tip of flagella 2 and are diffusive
    if ~isempty(find_td1) % in flagella 2
        motor(3,find_td1) =  2*len2(time+1) - motor(3,find_td1); % reflecting boundary condition
    end
    
    %finding motors which have reached the basal body
    find_ba = find((motor(3,:)<=0 & motor(1,:)==1 & motor(2,:)==2)); % in flagella 1
    find_ba1= find((motor(3,:)<=0 & motor(1,:)==2 & motor(2,:)==2)); % in flagella 2
    if ~isempty(find_ba)
        motor(:,find_ba)=0;       
        freeM=freeM+length(find_ba);
    end
    if ~isempty(find_ba1)
        motor(:,find_ba1)=0;       
        freeM=freeM+length(find_ba1);
    end
    
end

% plot length vs time
figure
plot(time_total/60, len1(1:N+1),'Linewidth',2);
hold on
plot(time_total/60, len2(1:N+1),'Linewidth',2);
hold off
legend('Flagellum 1', 'Flagellum 2');
ylabel('Flagellum length(\mum)')
set(gca, 'fontsize', 18);
ylim([0,14]);    
toc