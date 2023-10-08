% IKOCP Protocol
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%% Variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
InitialNodeEnergy=5;
InitialTemperatureValue=30;
RWP_Point_Count=100;
RWP_Point=zeros(RWP_Point_Count,3); % Index | X | Y

SimulationTime=600;
TotalRoundList=zeros(SimulationTime,21); % 1: RoundNumber | 2: Minimum distance between actors | 3: Maximum distance between actors | 4: Average distance between actors | 5: Minimum distance between nodes | 6: Maximum distance between nodes | 7: Average distance between nodes | 8: Average data of the actors | 9: Average data of the nodes | 10: Average remaining energy of the nodes | 11: Total consumption energy of the nodes | 12: Average consumption energy of the nodes | 13: QOD | 14: Live nodes count | 15: Dead nodes count | 16: Filled buffer count of the nodes | 17: Filled buffer average of the nodes | 18: Filled buffer percent of the nodes | 19: Filled buffer count of the BaseStation | 20: Filled buffer percent of the BaseStation | 21: Produced packet count
TotalRoundList_Index=0;
ProducedPacketCount=0;
PacketGeneration_Count=100;

Network_Length=1000;
Network_Width=1000;

ActorNumber=50;
ActorRange=100;
MaximumActorValidValue=40;
ActorList=zeros(ActorNumber,6); % ID | X | Y | Z | Data | RWP_Point_Index
Actor_ChangePosition_RoundNumber=100;

NodeNumber=2000;
NodeRange=50;
NodeBufferSize=200;
NodeList=zeros(NodeNumber,20); % ID | X | Y | Z | RemainingEnergy | Data | Status (has sendable data or no) | Dead RoundNumber | Is Clusterhead (0=no, 1=yes) | Sink_ID | Clusterhead_ID1 | Clusterhead_Percentage1 | Clusterhead_ID2 | Clusterhead_Percentage2 | Clusterhead_ID3 | Clusterhead_Percentage3 | Clusterhead_ID4 | Clusterhead_Percentage4 | Clusterhead_ID5 | Clusterhead_Percentage5
%NodeBuffer=zeros(NodeNumber,NodeBufferSize);
NodeBuffer_Index=zeros(NodeNumber,1);
NodeInitiatorSeqNo=zeros(NodeNumber,1);
NodeNeighbors=zeros(NodeNumber,NodeNumber);
NodeNeighbors_Index=zeros(NodeNumber,1);
PeriodTime_CH_Sink=5; % Period time to transmit data packets from CHs to sinks

Alpha1=0.15;
Beta1=0.05;
bias=0.1;

Distance_ShapeFactor=80;
RemainingEnergy_ShapeFactor=5;
NumberOfMembers_ShapeFactor=10;
HitRate_ShapeFactor=1;

Distance_Type='b';
RemainingEnergy_Type='b';
NumberOfMembers_Type='i';
HitRate_Type='b';

SinkBufferSize=10000;
%SinkBuffer=zeros(SinkNumber,SinkBufferSize);
%SinkBuffer_Index=zeros(SinkNumber,1);

BaseStation=[Network_Length/2 Network_Width+30]; % X | Y
BaseStationBufferSize=20000;
%BaseStationBuffer=zeros(BaseStationBufferSize,1);
BaseStationBuffer_Index=0;

PacketSize=3000;
NodeThresholdEnergy=5*(10^(-9))*PacketSize;
d0=87.7;





%%%%% Membership functions of the fuzzy decision
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Input fuzzy
%%% Distance
MaximumDistance=80;
U_Distance=[0 20 40 60 80];
% Near
mu_Distance_Near=[1.0 0.0 0.0 0.0 0.0];
% Medium
mu_Distance_Medium=[0.0 0.0 1.0 0.0 0.0];
% Far
mu_Distance_Far=[0.0 0.0 0.0 0.0 1.0];

%%% Remaining energy
MaximumRemainingEnergy=5;
U_RemainingEnergy=[1 2 3 4 5];
% Low
mu_RemainingEnergy_Low=[1.0 0.0 0.0 0.0 0.0];
% Medium
mu_RemainingEnergy_Medium=[0.0 0.0 1.0 0.0 0.0];
% High
mu_RemainingEnergy_High=[0.0 0.0 0.0 0.0 1.0];

%%% Number of members
MaximumNumberOfMembers=10;
U_NumberOfMembers=[1 3 5 7 10];
% Few
mu_NumberOfMembers_Few=[1.0 0.0 0.0 0.0 0.0];
% Medium
mu_NumberOfMembers_Medium=[0.0 0.0 1.0 0.0 0.0];
% Large
mu_NumberOfMembers_Large=[0.0 0.0 0.0 0.0 1.0];

%%%%% Output fuzzy
%%% Hit rate
MaximumHitRate=1;
U_HitRate=[0 0.25 0.5 0.75 1];
% Low
mu_HitRate_Low=[1.0 0.0 0.0 0.0 0.0];
% Medium
mu_HitRate_Medium=[0.0 0.0 1.0 0.0 0.0];
% High
mu_HitRate_High=[0.0 0.0 0.0 0.0 1.0];





%%%%% Fuzzy rules
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Rule 1
mu_ABCD=fuzzyand(mu_Distance_Near,mu_RemainingEnergy_Low,mu_NumberOfMembers_Few);
R1=rulemakem(mu_ABCD,mu_HitRate_Medium);
%%% Rule 2
mu_ABCD=fuzzyand(mu_Distance_Near,mu_RemainingEnergy_Low,mu_NumberOfMembers_Medium);
R2=rulemakem(mu_ABCD,mu_HitRate_Medium);
%%% Rule 3
mu_ABCD=fuzzyand(mu_Distance_Near,mu_RemainingEnergy_Low,mu_NumberOfMembers_Large);
R3=rulemakem(mu_ABCD,mu_HitRate_Low);
%%% Rule 4
mu_ABCD=fuzzyand(mu_Distance_Near,mu_RemainingEnergy_Medium,mu_NumberOfMembers_Few);
R4=rulemakem(mu_ABCD,mu_HitRate_High);
%%% Rule 5
mu_ABCD=fuzzyand(mu_Distance_Near,mu_RemainingEnergy_Medium,mu_NumberOfMembers_Medium);
R5=rulemakem(mu_ABCD,mu_HitRate_Medium);
%%% Rule 6
mu_ABCD=fuzzyand(mu_Distance_Near,mu_RemainingEnergy_Medium,mu_NumberOfMembers_Large);
R6=rulemakem(mu_ABCD,mu_HitRate_Medium);
%%% Rule 7
mu_ABCD=fuzzyand(mu_Distance_Near,mu_RemainingEnergy_High,mu_NumberOfMembers_Few);
R7=rulemakem(mu_ABCD,mu_HitRate_High);
%%% Rule 8
mu_ABCD=fuzzyand(mu_Distance_Near,mu_RemainingEnergy_High,mu_NumberOfMembers_Medium);
R8=rulemakem(mu_ABCD,mu_HitRate_High);
%%% Rule 9
mu_ABCD=fuzzyand(mu_Distance_Near,mu_RemainingEnergy_High,mu_NumberOfMembers_Large);
R9=rulemakem(mu_ABCD,mu_HitRate_Medium);
%%% Rule 10
mu_ABCD=fuzzyand(mu_Distance_Medium,mu_RemainingEnergy_Low,mu_NumberOfMembers_Few);
R10=rulemakem(mu_ABCD,mu_HitRate_Medium);
%%% Rule 11
mu_ABCD=fuzzyand(mu_Distance_Medium,mu_RemainingEnergy_Low,mu_NumberOfMembers_Medium);
R11=rulemakem(mu_ABCD,mu_HitRate_Low);
%%% Rule 12
mu_ABCD=fuzzyand(mu_Distance_Medium,mu_RemainingEnergy_Low,mu_NumberOfMembers_Large);
R12=rulemakem(mu_ABCD,mu_HitRate_Low);
%%% Rule 13
mu_ABCD=fuzzyand(mu_Distance_Medium,mu_RemainingEnergy_Medium,mu_NumberOfMembers_Few);
R13=rulemakem(mu_ABCD,mu_HitRate_High);
%%% Rule 14
mu_ABCD=fuzzyand(mu_Distance_Medium,mu_RemainingEnergy_Medium,mu_NumberOfMembers_Medium);
R14=rulemakem(mu_ABCD,mu_HitRate_Medium);
%%% Rule 15
mu_ABCD=fuzzyand(mu_Distance_Medium,mu_RemainingEnergy_Medium,mu_NumberOfMembers_Large);
R15=rulemakem(mu_ABCD,mu_HitRate_Low);
%%% Rule 16
mu_ABCD=fuzzyand(mu_Distance_Medium,mu_RemainingEnergy_High,mu_NumberOfMembers_Few);
R16=rulemakem(mu_ABCD,mu_HitRate_High);
%%% Rule 17
mu_ABCD=fuzzyand(mu_Distance_Medium,mu_RemainingEnergy_High,mu_NumberOfMembers_Medium);
R17=rulemakem(mu_ABCD,mu_HitRate_Medium);
%%% Rule 18
mu_ABCD=fuzzyand(mu_Distance_Medium,mu_RemainingEnergy_High,mu_NumberOfMembers_Large);
R18=rulemakem(mu_ABCD,mu_HitRate_Medium);
%%% Rule 19
mu_ABCD=fuzzyand(mu_Distance_Far,mu_RemainingEnergy_Low,mu_NumberOfMembers_Few);
R19=rulemakem(mu_ABCD,mu_HitRate_Low);
%%% Rule 20
mu_ABCD=fuzzyand(mu_Distance_Far,mu_RemainingEnergy_Low,mu_NumberOfMembers_Medium);
R20=rulemakem(mu_ABCD,mu_HitRate_Low);
%%% Rule 21
mu_ABCD=fuzzyand(mu_Distance_Far,mu_RemainingEnergy_Low,mu_NumberOfMembers_Large);
R21=rulemakem(mu_ABCD,mu_HitRate_Low);
%%% Rule 22
mu_ABCD=fuzzyand(mu_Distance_Far,mu_RemainingEnergy_Medium,mu_NumberOfMembers_Few);
R22=rulemakem(mu_ABCD,mu_HitRate_Medium);
%%% Rule 23
mu_ABCD=fuzzyand(mu_Distance_Far,mu_RemainingEnergy_Medium,mu_NumberOfMembers_Medium);
R23=rulemakem(mu_ABCD,mu_HitRate_Medium);
%%% Rule 24
mu_ABCD=fuzzyand(mu_Distance_Far,mu_RemainingEnergy_Medium,mu_NumberOfMembers_Large);
R24=rulemakem(mu_ABCD,mu_HitRate_Low);
%%% Rule 25
mu_ABCD=fuzzyand(mu_Distance_Far,mu_RemainingEnergy_High,mu_NumberOfMembers_Few);
R25=rulemakem(mu_ABCD,mu_HitRate_Medium);
%%% Rule 26
mu_ABCD=fuzzyand(mu_Distance_Far,mu_RemainingEnergy_High,mu_NumberOfMembers_Medium);
R26=rulemakem(mu_ABCD,mu_HitRate_Medium);
%%% Rule 27
mu_ABCD=fuzzyand(mu_Distance_Far,mu_RemainingEnergy_High,mu_NumberOfMembers_Large);
R27=rulemakem(mu_ABCD,mu_HitRate_Low);

%%% Aggregation of rules 
R=totalrule(R1,R2,R3,R4,R5,R6,R7,R8,R9,R10,R11,R12,R13,R14,R15,R16,R17,R18,R19,R20,R21,R22,R23,R24,R25,R26,R27);





%%%%% Initialize
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Number of regions
NumberOfRegions=MCDM_Sink(Network_Length,Network_Width,NodeNumber);
Sink=[];
for i=1:NumberOfRegions
    % X
    if mod(i,sqrt(NumberOfRegions))==0
        Sink(i,1)=(Network_Length*(2*sqrt(NumberOfRegions)-1))/(2*sqrt(NumberOfRegions));
    elseif mod(i,sqrt(NumberOfRegions))==1
        Sink(i,1)=Network_Length/(2*sqrt(NumberOfRegions));
    elseif mod(i,sqrt(NumberOfRegions))==2
        Sink(i,1)=(3*Network_Length)/(2*sqrt(NumberOfRegions));  
    elseif mod(i,sqrt(NumberOfRegions))==3
        Sink(i,1)=(5*Network_Length)/(2*sqrt(NumberOfRegions));  
    elseif mod(i,sqrt(NumberOfRegions))==4
        Sink(i,1)=(7*Network_Length)/(2*sqrt(NumberOfRegions));        
    end;
    
    % Y
    if ceil(i/sqrt(NumberOfRegions))==1
        Sink(i,2)=Network_Width/(2*sqrt(NumberOfRegions));        
    elseif ceil(i/sqrt(NumberOfRegions))==2
        Sink(i,2)=(3*Network_Width)/(2*sqrt(NumberOfRegions));        
    elseif ceil(i/sqrt(NumberOfRegions))==3
        Sink(i,2)=(5*Network_Width)/(2*sqrt(NumberOfRegions));          
    elseif ceil(i/sqrt(NumberOfRegions))==4
        Sink(i,2)=(7*Network_Width)/(2*sqrt(NumberOfRegions));          
    elseif ceil(i/sqrt(NumberOfRegions))==5
        Sink(i,2)=(Network_Width*(2*sqrt(NumberOfRegions)-1))/(2*sqrt(NumberOfRegions));        
    end;      
end;

%%% Set RWP points
for i=1:RWP_Point_Count
    RWP_Point(i,1)=i;     
end;

RWP_Point=RWP_Point_Set(RWP_Point,RWP_Point_Count,Network_Length,Network_Width);

%%% NodeBuffer
for i=1:NodeNumber
    for j=1:NodeBufferSize
        NodeBuffer(i,j).Initiator_ID=0;
        NodeBuffer(i,j).InitiatorSeqNo=0;
        NodeBuffer(i,j).PartialRoute=[];
        NodeBuffer(i,j).Data=0;
        NodeBuffer(i,j).StartSend_RoundNumber=0;
    end;
end;

%%% SinkBuffer
SinkNumber=NumberOfRegions;
for i=1:SinkNumber
    for j=1:SinkBufferSize
        SinkBuffer(i,j).Initiator_ID=0;
        SinkBuffer(i,j).InitiatorSeqNo=0;
        SinkBuffer(i,j).PartialRoute=[];
        SinkBuffer(i,j).Data=0;
        SinkBuffer(i,j).StartSend_RoundNumber=0;
    end;
end;
SinkBuffer_Index=zeros(SinkNumber,1);

%%% BaseStationBuffer
for i=1:BaseStationBufferSize
    BaseStationBuffer(i,1).Initiator_ID=0;
    BaseStationBuffer(i,1).InitiatorSeqNo=0;
    BaseStationBuffer(i,1).PartialRoute=[];
    BaseStationBuffer(i,1).Data=0; % Data
    BaseStationBuffer(i,1).StartSend_RoundNumber=0;
    BaseStationBuffer(i,1).FinishSend_RoundNumber=0;
    BaseStationBuffer(i,1).Total_RoundNumber=0;
end;

%%% Actor
ActorList=ActorList_Set(ActorList,ActorNumber,RWP_Point_Count);

for i=1:ActorNumber
    ActorList(i,1)=i;
    ActorList(i,2)=RWP_Point(ActorList(i,6),2);
    ActorList(i,3)=RWP_Point(ActorList(i,6),3);
end;

%%% NodeList
NodeList=NodeList_Set(NodeList,Network_Length,Network_Width);

for i=1:NodeNumber
    NodeList(i,1)=i;
    NodeList(i,5)=InitialNodeEnergy;
    NodeList(i,6)=InitialTemperatureValue;     
end;

%%% PacketGeneration_List
PacketGeneration_List=zeros(1,PacketGeneration_Count);
PacketGeneration_List_Index=1;
PacketGeneration_List=PacketGeneration_List_Set(PacketGeneration_List);

%%% PacketGenerationRate
PacketGenerationRate=1:1:SimulationTime; 
 




%%%%% Cycle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Continue_Flag=1;
RoundNumber=1;
while (Continue_Flag==1)&&(RoundNumber<=SimulationTime)
    
    SensedDataCount=0;    
    AcceptedDataCount=0;

    %%% NodeConsumptionEnergy
    NodeConsumptionEnergy=zeros(NodeNumber,1);
    for i=1:NodeNumber
        NodeConsumptionEnergy(i,1)=NodeList(i,5);
    end;     
    
    %%% Change current position of the actors
    if mod(RoundNumber,Actor_ChangePosition_RoundNumber)==0
        for i=1:ActorNumber        
            if ActorList(i,6)<RWP_Point_Count
                ActorList(i,6)=ActorList(i,6)+1;
            else
                ActorList(i,6)=1;
            end;

            ActorList(i,2)=RWP_Point(ActorList(i,6),2);
            ActorList(i,3)=RWP_Point(ActorList(i,6),3);       
        end;
    end;    

    %%% Update NodeNeighbors
    if RoundNumber==1
        NodeNeighbors=zeros(NodeNumber,NodeNumber);
        NodeNeighbors_Index=zeros(NodeNumber,1);    
        for i=1:NodeNumber
            for j=1:NodeNumber
                if i~=j
                    Distance=sqrt(((NodeList(i,2)-NodeList(j,2))^2)+((NodeList(i,3)-NodeList(j,3))^2));
                    if Distance<=NodeRange
                        NodeNeighbors_Index(i)=NodeNeighbors_Index(i)+1;                               
                        NodeNeighbors(i,NodeNeighbors_Index(i))=j;
                    end;
                end;
            end;
        end; 
    end;
    
    %%% Associate nodes to sinks
    if RoundNumber==1 
        for i=1:NodeNumber
            AssociateSink=zeros(SinkNumber,2); % Sink_ID | Distance

            for j=1:SinkNumber
                AssociateSink(j,1)=j;
                AssociateSink(j,2)=sqrt(((NodeList(i,2)-NodeList(j,1))^2)+((NodeList(i,3)-NodeList(j,2))^2));
            end;

            AssociateSink=sortrows(AssociateSink,2);
            NodeList(i,10)=AssociateSink(1,1);
        end;
    end;
    
    %%% Clustering
    if RoundNumber==1        
        % Select clusterheads
        for i=1:NodeNumber
            %\ Determine the value of p
            SA_i=(Network_Length*Network_Width)/NumberOfRegions;
            SA=Network_Length*Network_Width;
            
            Sink_ID=NodeList(i,10);
            NN_i=0;
            for j=1:NodeNumber
                if NodeList(j,10)==Sink_ID
                    NN_i=NN_i+1;
                end;
            end;
            
            NN=NodeNumber;
            
            P=(Alpha1*SA_i/SA)+(Beta1*NN_i/NN)+bias;
            
            %\ T(n)
            r=RoundNumber;
            T_n=(P*NodeList(i,5))/((1-(P*(mod(r,(1/P)))))*InitialNodeEnergy);
            RandomNumber=rand();
            if RandomNumber<T_n
                NodeList(i,9)=1;
            else
                NodeList(i,9)=0;
            end;
        end;
    
        % Construct Clusters
        %/ Find the nearest Clusterhead nodes
        for i=1:NodeNumber
            if NodeList(i,9)==0  
                % Decrease consumption energy of the non cluster head node to check the neighboring cluster heads with the broadcast process
                Distance=NodeRange;                        
                if Distance<=d0
                    NodeList(i,5)=NodeList(i,5)-(50*(10^(-9))*PacketSize+10*(10^(-12))*PacketSize*(Distance^2));
                else
                    NodeList(i,5)=NodeList(i,5)-(50*(10^(-9))*PacketSize+13*(10^(-16))*PacketSize*(Distance^4));
                end;    
                
                for j=1:NodeNeighbors_Index(i,1)
                    Neighbor_ID=NodeNeighbors(i,j);
                    
                    if NodeList(Neighbor_ID,9)==1
                        if NodeList(i,11)==0
                            NodeList(i,11)=Neighbor_ID;
                        elseif NodeList(i,13)==0
                            NodeList(i,13)=Neighbor_ID;
                        elseif NodeList(i,15)==0
                            NodeList(i,15)=Neighbor_ID;
                        elseif NodeList(i,17)==0
                            NodeList(i,17)=Neighbor_ID;
                        elseif NodeList(i,19)==0
                            NodeList(i,19)=Neighbor_ID;
                        end;
                        
                        % Decrease consumption energy of the cluster head node to receive the broadcast
                        NodeList(Neighbor_ID,5)=NodeList(Neighbor_ID,5)-(50*(10^(-9))*PacketSize);
                        
                        % Decrease consumption energy of the cluster head node to send the response
                        Distance=sqrt(((NodeList(Neighbor_ID,2)-NodeList(i,2))^2)+((NodeList(Neighbor_ID,3)-NodeList(i,3))^2));                        
                        if Distance<=d0
                            NodeList(Neighbor_ID,5)=NodeList(Neighbor_ID,5)-(50*(10^(-9))*PacketSize+10*(10^(-12))*PacketSize*(Distance^2));
                        else
                            NodeList(Neighbor_ID,5)=NodeList(Neighbor_ID,5)-(50*(10^(-9))*PacketSize+13*(10^(-16))*PacketSize*(Distance^4));
                        end; 
                        
                        % Decrease consumption energy of the non cluster head node to receive the response
                        NodeList(i,5)=NodeList(i,5)-(50*(10^(-9))*PacketSize);                        
                    end;
                end;
            end;
        end;
        
        %%% Associate with appropriate CHs
        for i=1:NodeNumber
            if NodeList(i,9)==0                
                if NodeList(i,11)~=0
                    CH_ID=NodeList(i,11);
                    
                    % Distance
                    Sink_ID=NodeList(CH_ID,10);
                    Distance=sqrt(((NodeList(CH_ID,2)-Sink(Sink_ID,1))^2)+((NodeList(CH_ID,3)-Sink(Sink_ID,2))^2));
                    if Distance<=MaximumDistance
                        Crisp=Distance; 
                    else
                        Crisp=MaximumDistance;
                    end;                            
                    mu_Distance=fuzzifysn(U_Distance,Crisp,Distance_Type,Distance_ShapeFactor-1);
            
                    % Remaining energy
                    RemainingEnergy=NodeList(CH_ID,5);
                    if RemainingEnergy<=MaximumRemainingEnergy
                        Crisp=RemainingEnergy; 
                    else
                        Crisp=MaximumRemainingEnergy;
                    end;
                    mu_RemainingEnergy=fuzzifysn(U_RemainingEnergy,Crisp,RemainingEnergy_Type,RemainingEnergy_ShapeFactor-1); 
            
                    % Number of members
                    NumberOfMembers=0;
                    for j=1:NodeNumber
                        if (NodeList(j,11)==CH_ID)||(NodeList(j,13)==CH_ID)||(NodeList(j,15)==CH_ID)||(NodeList(j,17)==CH_ID)||(NodeList(j,19)==CH_ID)
                            NumberOfMembers=NumberOfMembers+1;
                        end;                        
                    end;
                    if NumberOfMembers<=MaximumNumberOfMembers
                        Crisp=NumberOfMembers; 
                    else
                        Crisp=MaximumNumberOfMembers;
                    end;
                    mu_NumberOfMembers=fuzzifysn(U_NumberOfMembers,Crisp,NumberOfMembers_Type,NumberOfMembers_ShapeFactor-1); 
            
                    %%% Hit rate
                    mu_HitRate=ruleresp(R,fuzzyand(mu_Distance,mu_RemainingEnergy,mu_NumberOfMembers));
                    HitRate=defuzzyg(U_HitRate,mu_HitRate);            
                    NodeList(i,12)=HitRate;                     
                end;
                if NodeList(i,13)~=0
                    CH_ID=NodeList(i,13);
                    
                    % Distance
                    Sink_ID=NodeList(CH_ID,10);
                    Distance=sqrt(((NodeList(CH_ID,2)-Sink(Sink_ID,1))^2)+((NodeList(CH_ID,3)-Sink(Sink_ID,2))^2));
                    if Distance<=MaximumDistance
                        Crisp=Distance; 
                    else
                        Crisp=MaximumDistance;
                    end;                            
                    mu_Distance=fuzzifysn(U_Distance,Crisp,Distance_Type,Distance_ShapeFactor-1);
            
                    % Remaining energy
                    RemainingEnergy=NodeList(CH_ID,5);
                    if RemainingEnergy<=MaximumRemainingEnergy
                        Crisp=RemainingEnergy; 
                    else
                        Crisp=MaximumRemainingEnergy;
                    end;
                    mu_RemainingEnergy=fuzzifysn(U_RemainingEnergy,Crisp,RemainingEnergy_Type,RemainingEnergy_ShapeFactor-1); 
            
                    % Number of members
                    NumberOfMembers=0;
                    for j=1:NodeNumber
                        if (NodeList(j,11)==CH_ID)||(NodeList(j,13)==CH_ID)||(NodeList(j,15)==CH_ID)||(NodeList(j,17)==CH_ID)||(NodeList(j,19)==CH_ID)
                            NumberOfMembers=NumberOfMembers+1;
                        end;                        
                    end;
                    if NumberOfMembers<=MaximumNumberOfMembers
                        Crisp=NumberOfMembers; 
                    else
                        Crisp=MaximumNumberOfMembers;
                    end;
                    mu_NumberOfMembers=fuzzifysn(U_NumberOfMembers,Crisp,NumberOfMembers_Type,NumberOfMembers_ShapeFactor-1); 
            
                    %%% Hit rate
                    mu_HitRate=ruleresp(R,fuzzyand(mu_Distance,mu_RemainingEnergy,mu_NumberOfMembers));
                    HitRate=defuzzyg(U_HitRate,mu_HitRate);            
                    NodeList(i,14)=HitRate; 
                end;
                if NodeList(i,15)~=0
                    CH_ID=NodeList(i,15);
                    
                    % Distance
                    Sink_ID=NodeList(CH_ID,10);
                    Distance=sqrt(((NodeList(CH_ID,2)-Sink(Sink_ID,1))^2)+((NodeList(CH_ID,3)-Sink(Sink_ID,2))^2));
                    if Distance<=MaximumDistance
                        Crisp=Distance; 
                    else
                        Crisp=MaximumDistance;
                    end;                            
                    mu_Distance=fuzzifysn(U_Distance,Crisp,Distance_Type,Distance_ShapeFactor-1);
            
                    % Remaining energy
                    RemainingEnergy=NodeList(CH_ID,5);
                    if RemainingEnergy<=MaximumRemainingEnergy
                        Crisp=RemainingEnergy; 
                    else
                        Crisp=MaximumRemainingEnergy;
                    end;
                    mu_RemainingEnergy=fuzzifysn(U_RemainingEnergy,Crisp,RemainingEnergy_Type,RemainingEnergy_ShapeFactor-1); 
            
                    % Number of members
                    NumberOfMembers=0;
                    for j=1:NodeNumber
                        if (NodeList(j,11)==CH_ID)||(NodeList(j,13)==CH_ID)||(NodeList(j,15)==CH_ID)||(NodeList(j,17)==CH_ID)||(NodeList(j,19)==CH_ID)
                            NumberOfMembers=NumberOfMembers+1;
                        end;                        
                    end;
                    if NumberOfMembers<=MaximumNumberOfMembers
                        Crisp=NumberOfMembers; 
                    else
                        Crisp=MaximumNumberOfMembers;
                    end;
                    mu_NumberOfMembers=fuzzifysn(U_NumberOfMembers,Crisp,NumberOfMembers_Type,NumberOfMembers_ShapeFactor-1); 
            
                    %%% Hit rate
                    mu_HitRate=ruleresp(R,fuzzyand(mu_Distance,mu_RemainingEnergy,mu_NumberOfMembers));
                    HitRate=defuzzyg(U_HitRate,mu_HitRate);            
                    NodeList(i,16)=HitRate;  
                end;
                if NodeList(i,17)~=0
                    CH_ID=NodeList(i,17);
                    
                    % Distance
                    Sink_ID=NodeList(CH_ID,10);
                    Distance=sqrt(((NodeList(CH_ID,2)-Sink(Sink_ID,1))^2)+((NodeList(CH_ID,3)-Sink(Sink_ID,2))^2));
                    if Distance<=MaximumDistance
                        Crisp=Distance; 
                    else
                        Crisp=MaximumDistance;
                    end;                            
                    mu_Distance=fuzzifysn(U_Distance,Crisp,Distance_Type,Distance_ShapeFactor-1);
            
                    % Remaining energy
                    RemainingEnergy=NodeList(CH_ID,5);
                    if RemainingEnergy<=MaximumRemainingEnergy
                        Crisp=RemainingEnergy; 
                    else
                        Crisp=MaximumRemainingEnergy;
                    end;
                    mu_RemainingEnergy=fuzzifysn(U_RemainingEnergy,Crisp,RemainingEnergy_Type,RemainingEnergy_ShapeFactor-1); 
            
                    % Number of members
                    NumberOfMembers=0;
                    for j=1:NodeNumber
                        if (NodeList(j,11)==CH_ID)||(NodeList(j,13)==CH_ID)||(NodeList(j,15)==CH_ID)||(NodeList(j,17)==CH_ID)||(NodeList(j,19)==CH_ID)
                            NumberOfMembers=NumberOfMembers+1;
                        end;                        
                    end;
                    if NumberOfMembers<=MaximumNumberOfMembers
                        Crisp=NumberOfMembers; 
                    else
                        Crisp=MaximumNumberOfMembers;
                    end;
                    mu_NumberOfMembers=fuzzifysn(U_NumberOfMembers,Crisp,NumberOfMembers_Type,NumberOfMembers_ShapeFactor-1); 
            
                    %%% Hit rate
                    mu_HitRate=ruleresp(R,fuzzyand(mu_Distance,mu_RemainingEnergy,mu_NumberOfMembers));
                    HitRate=defuzzyg(U_HitRate,mu_HitRate);            
                    NodeList(i,18)=HitRate;  
                end;
                if NodeList(i,19)~=0
                    CH_ID=NodeList(i,19);
                    
                    % Distance
                    Sink_ID=NodeList(CH_ID,10);
                    Distance=sqrt(((NodeList(CH_ID,2)-Sink(Sink_ID,1))^2)+((NodeList(CH_ID,3)-Sink(Sink_ID,2))^2));
                    if Distance<=MaximumDistance
                        Crisp=Distance; 
                    else
                        Crisp=MaximumDistance;
                    end;                            
                    mu_Distance=fuzzifysn(U_Distance,Crisp,Distance_Type,Distance_ShapeFactor-1);
            
                    % Remaining energy
                    RemainingEnergy=NodeList(CH_ID,5);
                    if RemainingEnergy<=MaximumRemainingEnergy
                        Crisp=RemainingEnergy; 
                    else
                        Crisp=MaximumRemainingEnergy;
                    end;
                    mu_RemainingEnergy=fuzzifysn(U_RemainingEnergy,Crisp,RemainingEnergy_Type,RemainingEnergy_ShapeFactor-1); 
            
                    % Number of members
                    NumberOfMembers=0;
                    for j=1:NodeNumber
                        if (NodeList(j,11)==CH_ID)||(NodeList(j,13)==CH_ID)||(NodeList(j,15)==CH_ID)||(NodeList(j,17)==CH_ID)||(NodeList(j,19)==CH_ID)
                            NumberOfMembers=NumberOfMembers+1;
                        end;                        
                    end;
                    if NumberOfMembers<=MaximumNumberOfMembers
                        Crisp=NumberOfMembers; 
                    else
                        Crisp=MaximumNumberOfMembers;
                    end;
                    mu_NumberOfMembers=fuzzifysn(U_NumberOfMembers,Crisp,NumberOfMembers_Type,NumberOfMembers_ShapeFactor-1); 
            
                    %%% Hit rate
                    mu_HitRate=ruleresp(R,fuzzyand(mu_Distance,mu_RemainingEnergy,mu_NumberOfMembers));
                    HitRate=defuzzyg(U_HitRate,mu_HitRate);            
                    NodeList(i,20)=HitRate;  
                end;                
            end;
        end;        
    end;        

    %%% Generate and transmit data  
    if length(strfind(PacketGenerationRate,RoundNumber))>0
        for i=1:ActorNumber
            ActorList(i,5)=PacketGeneration_List(1,PacketGeneration_List_Index);
            
            % PacketGeneration_List_Index
            if PacketGeneration_List_Index<PacketGeneration_Count
                PacketGeneration_List_Index=PacketGeneration_List_Index+1;
            else
                PacketGeneration_List_Index=1;
            end;
                        
            for j=1:NodeNumber
                Distance=sqrt(((ActorList(i,2)-NodeList(j,2))^2)+((ActorList(i,3)-NodeList(j,3))^2));
                if (Distance<=ActorRange)&&(NodeList(j,5)>NodeThresholdEnergy)
                    if ActorList(i,5)<=MaximumActorValidValue
                        NodeList(j,6)=ActorList(i,5);
                        NodeList(j,7)=1;    
                        
                        % Update remaining energy of the node
                        NodeList(j,5)=NodeList(j,5)-(50*(10^(-9))*PacketSize);
                        
                        AcceptedDataCount=AcceptedDataCount+1;                       
                    end;
                    SensedDataCount=SensedDataCount+1;
                end;
            end;
        end; 
    end;
    
    %%% Save sensed data of the nodes on their buffer or transmit to CH
    for i=1:NodeNumber
        if NodeList(i,7)==1
            % Non clusterhead
            if NodeList(i,9)==0
                % Select an appropriate CH
                CH_List=zeros(5,3); % CH_ID | Hit rate | Percentage
                CH_List_Index=0;
            
                if NodeList(i,11)~=0
                    CH_List_Index=CH_List_Index+1;
                    
                    CH_List(CH_List_Index,1)=NodeList(i,11);
                    CH_List(CH_List_Index,2)=NodeList(i,12);
                end;
                if NodeList(i,13)~=0
                    CH_List_Index=CH_List_Index+1;
                    
                    CH_List(CH_List_Index,1)=NodeList(i,13);
                    CH_List(CH_List_Index,2)=NodeList(i,14);                    
                end;
                if NodeList(i,15)~=0
                    CH_List_Index=CH_List_Index+1;
                    
                    CH_List(CH_List_Index,1)=NodeList(i,15);
                    CH_List(CH_List_Index,2)=NodeList(i,16);                    
                end;
                if NodeList(i,17)~=0
                    CH_List_Index=CH_List_Index+1;
                    
                    CH_List(CH_List_Index,1)=NodeList(i,17);
                    CH_List(CH_List_Index,2)=NodeList(i,18);                    
                end;
                if NodeList(i,19)~=0
                    CH_List_Index=CH_List_Index+1;
                    
                    CH_List(CH_List_Index,1)=NodeList(i,19);
                    CH_List(CH_List_Index,2)=NodeList(i,20);                    
                end;
                
                S=sum(CH_List(:,2));
                for j=1:CH_List_Index
                    CH_List(j,3)=(CH_List(j,2)/S)*100;
                end;
                
                R=rand()*100;
                if R==0
                    R=1;
                end;
                
                Selected_CH_ID=NodeList(i,11);
                S=0;
                for j=1:CH_List_Index
                    S=S+CH_List(j,3);
                    if R<=S
                        Selected_CH_ID=CH_List(j,1);
                        break;
                    end;
                end;  
                
                if Selected_CH_ID~=0
                    if NodeBuffer_Index(Selected_CH_ID,1)<NodeBufferSize
                        NodeBuffer_Index(Selected_CH_ID,1)=NodeBuffer_Index(Selected_CH_ID,1)+1;
                        NodeInitiatorSeqNo(i,1)=NodeInitiatorSeqNo(i,1)+1;
                        ProducedPacketCount=ProducedPacketCount+1; 
                    
                        NodeBuffer(Selected_CH_ID,NodeBuffer_Index(Selected_CH_ID,1)).Initiator_ID=i;
                        NodeBuffer(Selected_CH_ID,NodeBuffer_Index(Selected_CH_ID,1)).InitiatorSeqNo=NodeInitiatorSeqNo(i,1);                                
                        NodeBuffer(Selected_CH_ID,NodeBuffer_Index(Selected_CH_ID,1)).PartialRoute=[i];
                        NodeBuffer(Selected_CH_ID,NodeBuffer_Index(Selected_CH_ID,1)).Data=NodeList(i,6);
                        NodeBuffer(Selected_CH_ID,NodeBuffer_Index(Selected_CH_ID,1)).StartSend_RoundNumber=RoundNumber;

                        NodeList(i,7)=0;                         
                    end;
                    
                    % Decrease consumption energy of the node
                    Distance=sqrt(((NodeList(i,2)-NodeList(Selected_CH_ID,2))^2)+((NodeList(i,3)-NodeList(Selected_CH_ID,3))^2));                        
                    if Distance<=d0
                        NodeList(i,5)=NodeList(i,5)-(50*(10^(-9))*PacketSize+10*(10^(-12))*PacketSize*(Distance^2));
                    else
                        NodeList(i,5)=NodeList(i,5)-(50*(10^(-9))*PacketSize+13*(10^(-16))*PacketSize*(Distance^4));
                    end;
                    NodeList(Selected_CH_ID,5)=NodeList(Selected_CH_ID,5)-(50*(10^(-9))*PacketSize);                    
                end;
                
            % Clusterhead
            else 
                if NodeBuffer_Index(i,1)<NodeBufferSize
                    NodeBuffer_Index(i,1)=NodeBuffer_Index(i,1)+1;
                    NodeInitiatorSeqNo(i,1)=NodeInitiatorSeqNo(i,1)+1;
                    ProducedPacketCount=ProducedPacketCount+1; 

                    NodeBuffer(i,NodeBuffer_Index(i,1)).Initiator_ID=i;
                    NodeBuffer(i,NodeBuffer_Index(i,1)).InitiatorSeqNo=NodeInitiatorSeqNo(i,1);                                
                    NodeBuffer(i,NodeBuffer_Index(i,1)).PartialRoute=[i];
                    NodeBuffer(i,NodeBuffer_Index(i,1)).Data=NodeList(i,6);
                    NodeBuffer(i,NodeBuffer_Index(i,1)).StartSend_RoundNumber=RoundNumber;

                    NodeList(i,7)=0;                          
                end;                 
            end;
        end;
    end;    

    %%% Transmit data packets of the nodes' buffer to sinks
    for i=1:NodeNumber
        if mod(RoundNumber,PeriodTime_CH_Sink)==0
            % Clusterhead        
            if NodeList(i,9)==1
                DataList=[];
                DataList_Index=0;
                
                if NodeBuffer_Index(i,1)>0
                    for j=1:NodeBuffer_Index(i,1)
                        DataList_Index=DataList_Index+1;
                        DataList(DataList_Index)=NodeBuffer(i,j).Data;
                    end;
                    
                    NodeBuffer_Index(i,1)=0;
                end;
                
                % Aggregate data
                VeryLow=[];
                VeryLow_Index=0;
                
                Low=[];
                Low_Index=0;                
                
                Medium=[];
                Medium_Index=0;
                
                High=[];
                High_Index=0;
                
                VeryHigh=[];
                VeryHigh_Index=0; 
                
                for j=1:DataList_Index
                    Data=DataList(j);
                    
                    if (Data<-10)
                        VeryLow_Index=VeryLow_Index+1;
                        VeryLow(VeryLow_Index)=Data;
                    elseif (Data>=-10)&&(Data<5)
                        Low_Index=Low_Index+1;
                        Low(Low_Index)=Data;                        
                    elseif (Data>=5)&&(Data<20)
                        Medium_Index=Medium_Index+1;
                        Medium(Medium_Index)=Data;                        
                    elseif (Data>=20)&&(Data<35)
                        High_Index=High_Index+1;
                        High(High_Index)=Data;                        
                    else
                        VeryHigh_Index=VeryHigh_Index+1;
                        VeryHigh(VeryHigh_Index)=Data;                        
                    end;
                end;
                
                SortList=[];
                SortList(1,1)=1;
                SortList(1,2)=VeryLow_Index;
                SortList(2,1)=2;
                SortList(2,2)=Low_Index; 
                SortList(3,1)=3;
                SortList(3,2)=Medium_Index; 
                SortList(4,1)=4;
                SortList(4,2)=High_Index; 
                SortList(5,1)=5;
                SortList(5,2)=VeryHigh_Index;  
                
                sortrows(SortList,2);
                No=SortList(5,1);
                if No==1
                    Data_Index=VeryLow_Index;
                    Data=VeryLow;
                elseif No==2
                    Data_Index=Low_Index;
                    Data=Low;                    
                elseif No==3
                    Data_Index=Medium_Index;
                    Data=Medium;                    
                elseif No==4
                    Data_Index=High_Index;
                    Data=High;                    
                else
                    Data_Index=VeryHigh_Index;
                    Data=VeryHigh;                    
                end;
                
                if Data_Index~=0
                    sort(Data);
                    Min=Data(1);
                    Avg=sum(Data)/Data_Index;
                    Max=Data(Data_Index);
                    
                    Sink_ID=NodeList(i,10);
                    
                    % Transmit Min to the sink
                    if SinkBuffer_Index(Sink_ID,1)<SinkBufferSize
                        SinkBuffer_Index(Sink_ID,1)=SinkBuffer_Index(Sink_ID,1)+1;

                        SinkBuffer(Sink_ID,SinkBuffer_Index(Sink_ID,1)).Initiator_ID=i;
                        SinkBuffer(Sink_ID,SinkBuffer_Index(Sink_ID,1)).InitiatorSeqNo=NodeInitiatorSeqNo(i,1);                                
                        SinkBuffer(Sink_ID,SinkBuffer_Index(Sink_ID,1)).PartialRoute=[i];
                        SinkBuffer(Sink_ID,SinkBuffer_Index(Sink_ID,1)).Data=Min;
                        SinkBuffer(Sink_ID,SinkBuffer_Index(Sink_ID,1)).StartSend_RoundNumber=RoundNumber; 
                        
                        % Decrease consumption energy of the node
                        Distance=sqrt(((NodeList(i,2)-Sink(Sink_ID,1))^2)+((NodeList(i,3)-Sink(Sink_ID,2))^2));                        
                        if Distance<=d0
                            NodeList(i,5)=NodeList(i,5)-(50*(10^(-9))*PacketSize+10*(10^(-12))*PacketSize*(Distance^2));
                        else
                            NodeList(i,5)=NodeList(i,5)-(50*(10^(-9))*PacketSize+13*(10^(-16))*PacketSize*(Distance^4));
                        end;                        
                    end;                      
                    
                    % Transmit Avg to the sink
                    if SinkBuffer_Index(Sink_ID,1)<SinkBufferSize
                        SinkBuffer_Index(Sink_ID,1)=SinkBuffer_Index(Sink_ID,1)+1;

                        SinkBuffer(Sink_ID,SinkBuffer_Index(Sink_ID,1)).Initiator_ID=i;
                        SinkBuffer(Sink_ID,SinkBuffer_Index(Sink_ID,1)).InitiatorSeqNo=NodeInitiatorSeqNo(i,1);                                
                        SinkBuffer(Sink_ID,SinkBuffer_Index(Sink_ID,1)).PartialRoute=[i];
                        SinkBuffer(Sink_ID,SinkBuffer_Index(Sink_ID,1)).Data=Avg;
                        SinkBuffer(Sink_ID,SinkBuffer_Index(Sink_ID,1)).StartSend_RoundNumber=RoundNumber; 
                        
                        % Decrease consumption energy of the node
                        Distance=sqrt(((NodeList(i,2)-Sink(Sink_ID,1))^2)+((NodeList(i,3)-Sink(Sink_ID,2))^2));                        
                        if Distance<=d0
                            NodeList(i,5)=NodeList(i,5)-(50*(10^(-9))*PacketSize+10*(10^(-12))*PacketSize*(Distance^2));
                        else
                            NodeList(i,5)=NodeList(i,5)-(50*(10^(-9))*PacketSize+13*(10^(-16))*PacketSize*(Distance^4));
                        end;                         
                    end; 
                    
                    % Transmit Max to the sink                    
                    if SinkBuffer_Index(Sink_ID,1)<SinkBufferSize
                        SinkBuffer_Index(Sink_ID,1)=SinkBuffer_Index(Sink_ID,1)+1;

                        SinkBuffer(Sink_ID,SinkBuffer_Index(Sink_ID,1)).Initiator_ID=i;
                        SinkBuffer(Sink_ID,SinkBuffer_Index(Sink_ID,1)).InitiatorSeqNo=NodeInitiatorSeqNo(i,1);                                
                        SinkBuffer(Sink_ID,SinkBuffer_Index(Sink_ID,1)).PartialRoute=[i];
                        SinkBuffer(Sink_ID,SinkBuffer_Index(Sink_ID,1)).Data=Max;
                        SinkBuffer(Sink_ID,SinkBuffer_Index(Sink_ID,1)).StartSend_RoundNumber=RoundNumber; 
                        
                        % Decrease consumption energy of the node
                        Distance=sqrt(((NodeList(i,2)-Sink(Sink_ID,1))^2)+((NodeList(i,3)-Sink(Sink_ID,2))^2));                        
                        if Distance<=d0
                            NodeList(i,5)=NodeList(i,5)-(50*(10^(-9))*PacketSize+10*(10^(-12))*PacketSize*(Distance^2));
                        else
                            NodeList(i,5)=NodeList(i,5)-(50*(10^(-9))*PacketSize+13*(10^(-16))*PacketSize*(Distance^4));
                        end;                         
                    end;                     
                end;
            end;
        end;
    end;
    
    %%% Transmit data packets from the sinks to the base station
    for i=1:SinkNumber
        for j=1:SinkBuffer_Index(i,1)
            if BaseStationBuffer_Index<BaseStationBufferSize                
                BaseStationBuffer_Index=BaseStationBuffer_Index+1;

                BaseStationBuffer(BaseStationBuffer_Index,1).Initiator_ID=SinkBuffer(i,j).Initiator_ID;
                BaseStationBuffer(BaseStationBuffer_Index,1).InitiatorSeqNo=SinkBuffer(i,j).InitiatorSeqNo;
                BaseStationBuffer(BaseStationBuffer_Index,1).PartialRoute=SinkBuffer(i,j).PartialRoute;
                BaseStationBuffer(BaseStationBuffer_Index,1).Data=SinkBuffer(i,j).Data;
                BaseStationBuffer(BaseStationBuffer_Index,1).StartSend_RoundNumber=SinkBuffer(i,j).StartSend_RoundNumber;                                                 
            end;        
        end;
        
        SinkBuffer_Index(i,1)=0;
    end;    
        
    %%% TotalRoundList
    TotalRoundList(RoundNumber,1)=RoundNumber;    
    % 
    Distance=zeros(1,ActorNumber^2);
    Distance_Index=0;
    for i=1:ActorNumber
        for j=1:ActorNumber
            if i<j
                Distance_Index=Distance_Index+1;
                
                Distance_Value=sqrt(((ActorList(i,2)-ActorList(j,2))^2)+((ActorList(i,3)-ActorList(j,3))^2));
                Distance(1,Distance_Index)=Distance_Value;
            end;
        end;
    end;
    TotalRoundList(RoundNumber,2)=min(Distance(1,1:Distance_Index)); 
    % 
    Distance=zeros(1,ActorNumber^2);
    Distance_Index=0;
    for i=1:ActorNumber
        for j=1:ActorNumber
            if i<j
                Distance_Index=Distance_Index+1;
                
                Distance_Value=sqrt(((ActorList(i,2)-ActorList(j,2))^2)+((ActorList(i,3)-ActorList(j,3))^2));
                Distance(1,Distance_Index)=Distance_Value;
            end;
        end;
    end;
    TotalRoundList(RoundNumber,3)=max(Distance(1,1:Distance_Index)); 
    % 
    Distance=zeros(1,ActorNumber^2);
    Distance_Index=0;
    for i=1:ActorNumber
        for j=1:ActorNumber
            if i<j
                Distance_Index=Distance_Index+1;
                
                Distance_Value=sqrt(((ActorList(i,2)-ActorList(j,2))^2)+((ActorList(i,3)-ActorList(j,3))^2));
                Distance(1,Distance_Index)=Distance_Value;
            end;
        end;
    end;
    TotalRoundList(RoundNumber,4)=round((sum(Distance(1,1:Distance_Index)))/ActorNumber);     
    % 
    Distance=zeros(1,NodeNumber^2);
    Distance_Index=0;
    for i=1:NodeNumber
        for j=1:NodeNumber
            if i<j
                Distance_Index=Distance_Index+1;
                
                Distance_Value=sqrt(((NodeList(i,2)-NodeList(j,2))^2)+((NodeList(i,3)-NodeList(j,3))^2));
                Distance(1,Distance_Index)=Distance_Value;
            end;
        end;
    end;
    TotalRoundList(RoundNumber,5)=min(Distance(1,1:Distance_Index)); 
    %
    Distance=zeros(1,NodeNumber^2);
    Distance_Index=0;
    for i=1:NodeNumber
        for j=1:NodeNumber
            if i<j
                Distance_Index=Distance_Index+1;
                
                Distance_Value=sqrt(((NodeList(i,2)-NodeList(j,2))^2)+((NodeList(i,3)-NodeList(j,3))^2));
                Distance(1,Distance_Index)=Distance_Value;
            end;
        end;
    end;
    TotalRoundList(RoundNumber,6)=max(Distance(1,1:Distance_Index)); 
    %
    Distance=zeros(1,NodeNumber^2);
    Distance_Index=0;
    for i=1:NodeNumber
        for j=1:NodeNumber
            if i<j
                Distance_Index=Distance_Index+1;
                
                Distance_Value=sqrt(((NodeList(i,2)-NodeList(j,2))^2)+((NodeList(i,3)-NodeList(j,3))^2));
                Distance(1,Distance_Index)=Distance_Value;
            end;
        end;
    end;
    TotalRoundList(RoundNumber,7)=round((sum(Distance(1,1:Distance_Index)))/NodeNumber);
    %
    Sum=0;
    for i=1:ActorNumber
        Sum=Sum+ActorList(i,5);
    end;
    TotalRoundList(RoundNumber,8)=Sum/ActorNumber; 
    %
    Sum=0;
    for i=1:NodeNumber
        Sum=Sum+NodeList(i,6);
    end;
    TotalRoundList(RoundNumber,9)=Sum/NodeNumber;            
    %
    Sum=0;
    for i=1:NodeNumber
        Sum=Sum+NodeList(i,5);
    end;
    if (Sum/NodeNumber)>NodeThresholdEnergy
        TotalRoundList(RoundNumber,10)=Sum/NodeNumber;    
    else
        TotalRoundList(RoundNumber,10)=NodeThresholdEnergy;
    end;   
    %
    for i=1:NodeNumber
        NodeConsumptionEnergy(i,1)=NodeConsumptionEnergy(i,1)-NodeList(i,5);
    end;    
    Sum=0;
    for i=1:NodeNumber
        Sum=Sum+NodeConsumptionEnergy(i,1);
    end;
    TotalRoundList(RoundNumber,11)=Sum;    
    %
    TotalRoundList(RoundNumber,12)=TotalRoundList(RoundNumber,11)/NodeNumber;
    %
    if SensedDataCount>0
        TotalRoundList(RoundNumber,13)=(AcceptedDataCount/SensedDataCount)*100;
    else
        TotalRoundList(RoundNumber,13)=0;
    end;
    %
    Count=0;
    for i=1:NodeNumber
        if NodeList(i,5)>NodeThresholdEnergy
            Count=Count+1;
        end;
    end;
    TotalRoundList(RoundNumber,14)=Count; 
    %
    TotalRoundList(RoundNumber,15)=NodeNumber-TotalRoundList(RoundNumber,14);        
    %
    Sum=0;
    for i=1:NodeNumber
        Sum=Sum+NodeBuffer_Index(i,1);
    end; 
    TotalRoundList(RoundNumber,16)=Sum; 
    %
    TotalRoundList(RoundNumber,17)=TotalRoundList(RoundNumber,16)/NodeNumber;    
    %
    TotalRoundList(RoundNumber,18)=(TotalRoundList(RoundNumber,16)/(NodeNumber*NodeBufferSize))*100;                 
    % 
    TotalRoundList(RoundNumber,19)=BaseStationBuffer_Index;     
    % 
    TotalRoundList(RoundNumber,20)=BaseStationBuffer_Index/BaseStationBufferSize;     
    % 
    TotalRoundList(RoundNumber,21)=ProducedPacketCount;
    
    %%% Check dead state of the nodes
    for i=1:NodeNumber
        if (NodeList(i,8)==0)&&(NodeList(i,5)<NodeThresholdEnergy)
            NodeList(i,8)=RoundNumber;    
        end;
    end; 
    
    %%% Check dead state of the nodes
    if TotalRoundList(RoundNumber,15)>=(NodeNumber/2)
        %Continue_Flag=0;
    end;    
    
    %%% Display and increase RoundNumber
    RoundNumber
    RoundNumber=RoundNumber+1;       
end; 





%%%%% Display summary tables and charts of the rounds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%%% Summary table
TotalRoundList=TotalRoundList(1:(RoundNumber-1),:);
TotalRoundList

%%% Charts
% Average remaining energy of the nodes
figure(1);
plot(TotalRoundList(:,1),TotalRoundList(:,10),'-b','LineWidth',2);
%title('Average Remaining Energy of the Nodes');
xlabel('Simulation time (min)');
ylabel('Average remaining energy of the nodes (J)');

% Total energy consumption of the nodes
figure(2);
plot(TotalRoundList(:,1),TotalRoundList(:,11),'-r','LineWidth',2);
%title('Total energy consumption of the nodes');
xlabel('Simulation time (min)');
ylabel('Total energy consumption of the nodes (J)'); 

% Average energy consumption of the nodes
figure(3);
plot(TotalRoundList(:,1),TotalRoundList(:,12),'-k','LineWidth',2);
%title('Average energy consumption of the nodes');
xlabel('Simulation time (min)');
ylabel('Average energy consumption of the nodes (J)');     

% Alive nodes count
figure(4);
plot(TotalRoundList(:,1),TotalRoundList(:,14),'-b','LineWidth',2);
%title('Number of Alive Nodes');
xlabel('Simulation time (min)');
ylabel('Number of alive nodes');    

% Filled buffer count of the nodes
figure(5);
plot(TotalRoundList(:,1),TotalRoundList(:,16),'-r','LineWidth',2);
%title('Number of packets on nodes queues');
xlabel('Simulation time (min)');
ylabel('Number of packets on nodes queues');   

% Filled buffer average of the nodes
figure(6);
plot(TotalRoundList(:,1),TotalRoundList(:,17),'-k','LineWidth',2);
%title('Average of packets on nodes queues');
xlabel('Simulation time (min)');
ylabel('Average of packets on nodes queues'); 

% Filled buffer percent of the nodes
figure(7);
plot(TotalRoundList(:,1),TotalRoundList(:,18),'-m','LineWidth',2);
%title('Percent of packets on nodes queues');
xlabel('Simulation time (min)');
ylabel('Percent of packets on nodes queues (%)'); 

% Number of packets delivered to the base station
figure(8);
plot(TotalRoundList(:,1),TotalRoundList(:,19),'-b','LineWidth',2);
%title('Number of packets delivered to the base station');
xlabel('Simulation time (min)');
ylabel('Number of packets delivered to the base station');     





%%%%% Set NodeThresholdEnergy of the deaded nodes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:NodeNumber
    if NodeList(i,5)<NodeThresholdEnergy
        NodeList(i,5)=NodeThresholdEnergy;
    end;
end;





%%%%% Display summary tables and charts of the BaseStation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Summary table
TotalBaseStationDataList=zeros(BaseStationBuffer_Index,8); % 1: PacketNumber | 2: Initiator_ID | 3: InitiatorSeqNo | 4: IntermediateNodesCount | 5: Data | 6: StartSend_RoundNumber
for i=1:BaseStationBuffer_Index
    TotalBaseStationDataList(i,1)=i;
    TotalBaseStationDataList(i,2)=BaseStationBuffer(i,1).Initiator_ID;
    TotalBaseStationDataList(i,3)=BaseStationBuffer(i,1).InitiatorSeqNo;
    TotalBaseStationDataList(i,4)=length(BaseStationBuffer(i,1).PartialRoute);
    TotalBaseStationDataList(i,5)=BaseStationBuffer(i,1).Data;
    TotalBaseStationDataList(i,6)=BaseStationBuffer(i,1).StartSend_RoundNumber;        
end;
TotalBaseStationDataList


            


%%%%% Display information of the nodes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

NodeList

 



%%%%% Display single outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

% Number of Produced Packets
ProducedPacketCount

% Number of Delivered Packets
DeliveredPacketCount=BaseStationBuffer_Index

% Packet Delivery Ratio
PacketDeliveryRatio=(DeliveredPacketCount/ProducedPacketCount)*100

% Number of Alive Nodes
Count=0;
for i=1:NodeNumber
    if NodeList(i,5)>NodeThresholdEnergy
        Count=Count+1;    
    end;
end; 
NumberOfAliveNodes=Count

% Average Remaining Energy
AverageRemainingEnergy=TotalRoundList(RoundNumber-1,10)

% Number of Lost Packets
LostPacketCount=ProducedPacketCount-DeliveredPacketCount

% Packet Delivery Latency
SumLatency=0;
for i=1:BaseStationBuffer_Index
    SumLatency=SumLatency+BaseStationBuffer(i,1).Total_RoundNumber;    
end; 
PacketDeliveryLatency=SumLatency/BaseStationBuffer_Index
