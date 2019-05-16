close all
clear all

% two flagella diffusive- motors not shared tubulin shared with constant disassembly

% to implement severing uncomment section starting with comment 'begin severing' and
% ending with 'end severing'

% input: starts with comment 'parameters of assembly process' and ends with 'input of user ends'
% output: plot of length of both the flagella(len, len1) vs time

% variable motor and motor1 stores properties of the motors in basal pool 1
% and pool 2 respectively
% first row stores if in base or flagella
% second one stores if ballistic or diffusive, 
% third one stores the position of the motor
% fourth row stores the amount of tubulin which the motor carries

% boolean values
% basal=0;%if motors are in basal body
% flagella=1;%if motors are in the flagella
% ballistic=1;%undergoing ballistic motion
% diffuse=2;%undergoing diffusing motion

tic

%%(input by user)
% paramters of assembly process
n=127;% amount of time to simulate the process in minutes
D=1.7; % Diffusion constant (microm^2/s) 
T=84; %total tubulin(microm)
gamma=2.5*10^-4; %constant gamma
gammakM = 0.03/60; %value of gamma * kon *M
da=0.5/60; %Disassembly rate(um/s)
v=2.5; %velocity of ballistic path(microm/s) 
Lss = 10.5; % length of steady state of flagella
sevt=20; % time at which flagella is severed in minutes
dt=0.032; %time step dt in seconds
% input of user ends

% calculating other parameters from input parameters
kon = (Lss/v+Lss^2/(2*D))^(-1)*(gammakM*(T-2*Lss)/da-1); %value of k
gammaM = gammakM/kon;
M=ceil(gammaM/gamma); %number of motors in each flagella

%program parameters
dx=v*dt; %distance moved per time step dt
time_total=0:dt:100*60*n/100; %in seconds
N=3000000/4*n/100*0.008/dt; % N=n/dt where n is in time of simulation

%parameters updated at each time step in simulation
motor=zeros(4, M); % stores properties of the motors in basal pool 1
motor1=zeros(4, M); % stores properties of the motors in basal pool 2
% first row stores if in base or flagella
% second one stores if ballistic or diffusive, 
% third one stores the position of the motor
% fourth row stores the amount of tubulin which the motor carries
len = ones(1, N+1)* Lss; %initialise and store length vs time
len1 = ones(1, N+1)* Lss; %initialise and store length vs time
freeM=M; % number of free motors in basal pool 1
freeT=T-len(1)-len1(1); % amount of free tubulin in basal pool
freeM1=M; % number of free motors in basal pool 2
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
    
    % checks if the time = sevt minutes, the flagella at that point is
    % severed to 1.5 micron 
    
    %%begin severing
%     if time==3000000/4*sevt/100*0.008/dt
%         sever_amt = len(time)-1.5;
%         find_b2 = find(motor(1,:)==1 & motor(3,:)>(len(time)-sever_amt));
%         len(time)=len(time)-sever_amt;
%         if ~isempty(find_b2)
%             motor(1:3, find_b2) = -1; %finding motors which are in severed part and removing them
%         end
%     end
    %%end severing

    % constant disassembly proces in both the flagella ensuring that
    % lengths are non-negative
    if (len(time) <= da*dt)
        freeT=freeT+len(time);
        len(time+1)=0;
    else
        len(time+1)=len(time)-da*dt;
        freeT=freeT+da*dt;
    end
    if (len1(time) <= da*dt)
        freeT=freeT+len1(time);
        len1(time+1)=0;
    else
        len1(time+1)=len1(time)-da*dt;
        freeT=freeT+da*dt;
    end    
    
    % injecting a motor into the flagella with probability prob_motor
    
    % flagella 1 motor injection
    prob_motor = kon*freeM*dt; %probability that motor is injected
    if rand<=prob_motor
        %randomly choosing a motor which is in basal body to move
        find_motor = find(~motor(1,:));
        change = randi(length(find_motor));
        mot = find_motor(change); %motor chosen to be inserted
        motor(1:2,mot) = 1; % motor number in mot inserted into flagellum
        motor(4, mot) = gamma*freeT; %amount of tubulin the motor carries
        freeT = freeT-gamma*freeT;
        freeM=freeM-1; % decrease the free motors in basal pool 1 by 1
    end
    % flagella 2 motor injection
    prob_motor1 = kon*freeM1*dt; %probability that motor is injected
    if rand<=prob_motor1
        %randomly choosing a motor which is in basal body to move
        find_motor1 = find(~motor1(1,:));
        change1 = randi(length(find_motor1));
        mot1 = find_motor1(change1); %motor chosen to be inserted
        motor1(1:2,mot1) = 1; % motor number in mot inserted into flagellum
        motor1(4, mot1) = gamma*freeT; %amount of tubulin the motor carries
        freeT = freeT-gamma*freeT;
        freeM1=freeM1-1; % decrease the free motors in basal pool 2 by 1
    end
    
    find_b = find(motor(2,:)==1);
    if ~isempty(find_b)
        motor(3, find_b) = motor(3, find_b) + dx; %finding motors in flagella 1 which are in ballistic and move them by a constant amount
    end
    find_b1 = find(motor1(2,:)==1);
    if ~isempty(find_b1)
        motor1(3, find_b1) = motor1(3, find_b1) + dx; %finding motors in flagella 2 which are in ballistic and move them by a constant amount
    end
    find_d = find(motor(2,:)==2); %finding motors which are diffusive(in flagella 1)
    if ~isempty(find_d)
        d=sqrt(2*D*dt)*randn(1,length(find_d));
        motor(3, find_d) = motor(3, find_d) + d;
    end   
    find_d1 = find(motor1(2,:)==2); %finding motors which are diffusive(in flagella 2)
    if ~isempty(find_d1)
        d1=sqrt(2*D*dt)*randn(1,length(find_d1));
        motor1(3, find_d1) = motor1(3, find_d1) + d1;
    end  
    find_t = find(motor(3,:)>=len(time) & motor(2, :)==1); %finding motors which have reached the tip of flagella 1 and are ballistic
    if ~isempty(find_t)
        for i=1:length(find_t)
            len(time+1)=len(time+1)+ motor(4, find_t(i)); %adding tubulin pool to the flagella 1 from motor
        end
        motor(2,find_t) = 2;
        motor(3,find_t) = len(time+1);
    end
    find_td= find(motor(3,:)>=len(time+1) & motor(2, :)==2); %finding motors which have reached the tip of flagella 1 and are diffusive
    if ~isempty(find_td)
        motor(3,find_td) =  2*len(time+1) - motor(3,find_td); % reflecting boundary condition
    end
    find_t1 = find(motor1(3,:)>=len1(time) & motor1(2, :)==1); %finding motors which have reached the tip of flagella 2 and are ballistic
    if ~isempty(find_t1)
        for i=1:length(find_t1)
            len1(time+1)=len1(time+1)+ motor1(4, find_t1(i)); %adding tubulin pool to the flagella 2 from motor
        end
        motor1(2,find_t1) = 2;
        motor1(3,find_t1) = len1(time+1);
    end
    find_td1= find(motor1(3,:)>=len1(time+1) & motor1(2, :)==2); %finding motors which have reached the tip of flagella 2 and are diffusive
    if ~isempty(find_td1)
        motor1(3,find_td1) =  2*len1(time+1) - motor1(3,find_td1); % reflecting boundary condition
    end
    find_ba = find(motor(3,:)<=0 & motor(1,:)==1 & motor(2,:)==2); %finding motors which have reached the basal body 1
    if ~isempty(find_ba)
        motor(:,find_ba)=0;       
        freeM=freeM+length(find_ba);
    end
    find_ba1 = find(motor1(3,:)<=0 & motor1(1,:)==1 & motor1(2,:)==2); %finding motors which have reached the basal body 2
    if ~isempty(find_ba1)
        motor1(:,find_ba1)=0;       
        freeM1=freeM1+length(find_ba1);
    end
end

% plot length vs time
figure
plot(time_total/60, len(1:N+1),'Linewidth',2);
hold on
plot(time_total/60, len1(1:N+1),'Linewidth',2);
hold off
legend('Flagellum 1', 'Flagellum 2');
ylabel('Flagellum length(\mum)')
set(gca, 'fontsize', 18);
ylim([0,14]);
toc
