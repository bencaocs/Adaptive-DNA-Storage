%_______________________________________________________________________________________%
%  B Multi-Verse Optimizer (MVO) source codes demo version 1.0                            %
%                                                                                       %
%  Developed in MATLAB R2011b(7.13)                                                     %
%                                                                                       %
%  Author and programmer: Seyedali Mirjalili                                            %
%                                                                                       %
%         e-Mail: ali.mirjalili@gmail.com                                               %
%                 seyedali.mirjalili@griffithuni.edu.au                                 %
%                                                                                       %
%       Homepage: http://www.alimirjalili.com                                           %
%                                                                                       %
%   Main paper:                                                                         %
%                                                                                       %
%   S. Mirjalili, S. M. Mirjalili, A. Hatamlou                                          %
%   Multi-Verse Optimizer: a nature-inspired algorithm for global optimization          %
%   Neural Computing and Applications, in press,2015,                                   %
%   DOI: http://dx.doi.org/10.1007/s00521-015-1870-7                                    %
%                                                                                       %
%_______________________________________________________________________________________%

% You can simply define your cost in a seperate file and load its handle to fobj
% The initial parameters that you need are:
%__________________________________________
% fobj = @YourCostFunction
% dim = number of your variables
% Max_iteration = maximum number of generations
% SearchAgents_no = number of search agents
% lb=[lb1,lb2,...,lbn] where lbn is the lower bound of variable n
% ub=[ub1,ub2,...,ubn] where ubn is the upper bound of variable n
% If all the variables have equal lower bound you can just
% define lb and ub as two single number numbers

% To run MVO: [Best_score,Best_pos,cg_curve]=MVO(Universes_no,Max_iteration,lb,ub,dim,fobj)
%__________________________________________

function [newDNA]=EGNA_BMVO(N,Max_time,lb,ub,dim,Length,H_target,count)

%Two variables for saving the position and inflation rate (fitness) of the best universe
Best_universe=zeros(1,dim);
Best_universe_Inflation_rate=0;

%Initialize the positions of universes
%Universes=initialization(N,dim,ub,lb);%N 随机的一列数

Universes=NG_InitDNASet(Length,H_target,count);
EditDist = EditDistMatrix(Universes); 
%Minimum and maximum of Wormhole Existence Probability (min and max in
% Eq.(3.3) in the paper
WEP_Max=1;
WEP_Min=0.2;
n = size(Universes,1); 
Convergence_curve=zeros(1,Max_time);
FirstUniverses=[];
%Iteration(time) counter
Time=1;
%nDNA=[];
nDNA=NG_InitDNASet(Length,H_target,10000);
%Main loop
while Time<Max_time+1
      individual=[];
      DNA=[]; 
      temp = randperm(n);
      C1=(0.45*exp(1)/2+0.4*cos(2*Time))*((exp(-((8*Time)/(3*Max_time))^2)));
    for i=1:size(Universes,1) 
        individual(i).n=i;               %individual结构体
        individual(i).p=[];
        DNA(i).p = [temp(i)];
        DNA(i).n=[];
        A=find(EditDist(temp(i),:)>=H_target);%hammingdist  是DNAset的距离对称矩阵
        individual(i).n=sum(A); %我加的
        individual(i).p = A;%tenp(i)行如果hamming distance >H_target 存入individua(i).p
    end

%     a=0;
%     for i=1:size(DNASet,1) 
%         if a<length(individual(i).p)
%             a=length(individual(i).p);
%         end
%     end
%     b=size(individual,2);
%     index = 1;
%     del = [];
%     for k=1:b
%         if  a>length(individual(k).p)
%             del(index) = k;
%             index = index + 1;
%         end
% 
%     end
% 
%         individual(del)=[];
%     for j=1:size(individual,2)
%         Universes(j,:)=individual(j).p;
%     end
  
%         for j=1:a
%             Universes(i,j)=individual(i).p(j);
%         end
%         individual(i).n=sum(-A); %我加的
%         individual(i).p = A;%tenp(i)行如果hamming distance >H_target 存入individua(i).p
    
    %Eq. (3.3) in the paper
    WEP=WEP_Min+Time*((WEP_Max-WEP_Min)/Max_time);
    
    %Travelling Distance Rate (Formula): Eq. (3.4) in the paper
    TDR=1-((Time)^(1/6)/(Max_time)^(1/6));
    
    %Inflation rates (I) (fitness values)
    Inflation_rates=zeros(1,size(Universes,1));
    
    for i=1:size(Universes,1)
        
        %Boundary checking (to bring back the universes inside search
        % space if they go beyoud the boundaries
        Flag4ub=Universes(i,:)>ub;
        Flag4lb=Universes(i,:)<lb;
        Universes(i,:)=(Universes(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        
        %Calculate the inflation rate (fitness) of universes
        if Time==1
            Inflation_rates(1,i)=individual(i).n;%通胀率=下面函数的结果
        else
            Inflation_rates(1,i)=EditDistSUM(FirstUniverses,Universes(i,:));
        end
        
        
        %fobj 
        %function o = F17(x)
        %    o=(x(2)-(x(1)^2)*5.1/(4*(pi^2))+5/pi*x(1)-6)^2+10*(1-1/(8*pi))*cos(x(1))+10;
        
        
        %Elitism
        if Inflation_rates(1,i)>Best_universe_Inflation_rate
            Best_universe_Inflation_rate=Inflation_rates(1,i);%找最小的Inflation_rate = best_universe
            Best_universe=Universes(i,:);%最好的宇宙是最小的Inflation_rate的i时的Universes
        end
%         if i==size(Universes,1)&&Time==1
%         
%         end
    end
       
        if i==size(Universes,1)&&Time==1
         %if size(Universes,1)~=1
          del = [];
         B=[];
         
         P=ismember(Universes,Best_universe,'rows');
         a=find(P==1);%best_universe的位置
         B=individual(temp(a)).p;
         newDNA=Best_universe;
         index = 1;
         b=find(temp==a);
         FirstUniverses=cat(1,Universes(1:a,:),Universes(a+1:n,:));
       for i=1:length(B)
        DisH=Editconstraint(Universes(B(i),:),newDNA);
        len1=find(DisH<H_target);
       if isempty(len1) 
           newDNA=[Universes(B(i),:);newDNA];
            elseif length(len1)==1  
                newDNA(len1,:)=[];
                newDNA=[Universes(B(i),:);newDNA];
       end
       end
        
       %else
           %newDNA=Best_universe;
        % end
        end
%          for  k= 1:length(B)
%                if HammingDist(B(k),a) < H_target 
%                  del(index) = k;
%                  index = index + 1;
%                else
%                  DNASeat= [B(k);] 
%                end  
%          end
%          individual(a).p(del)=[]; 
%          for j=1:size(individual,2)
%          DNA(j).n=individual(j).p;
%          end
    [sorted_Inflation_rates,sorted_indexes]=sort(Inflation_rates);%对Inflation_rates排序 [a,b]=sort(A) a是排序结果，b储存下标 
    
    for newindex=1:size(Universes,1)
        Sorted_universes(newindex,:)=Universes(sorted_indexes(newindex),:);
    end%按照inflation_rate 大小把Universes排序存到sorted_universes

    %Normaized inflation rates (NI in Eq. (3.1) in the paper)标准膨胀率
    normalized_sorted_Inflation_rates=normr(sorted_Inflation_rates);%归一化 使该矩阵每一列的平均值为0,标准差为1
    %normr是使一个矩阵的行或列标准化，即标准化使m矩阵的列，使其长度为1
    
    Universes(1,:)= Sorted_universes(1,:);

    %Update the Position of universes
    for i=2:size(Universes,1)%Starting from 2 since the firt one is the elite
        Back_hole_index=i;
        for j=1:size(Universes,2)
            r1=rand();
            if r1<normalized_sorted_Inflation_rates(i)
                White_hole_index=RouletteWheelSelection(sorted_Inflation_rates);% for maximization problem -sorted_Inflation_rates should be written as sorted_Inflation_rates
                if White_hole_index==-1
                    White_hole_index=1;
                end
                %Eq. (3.1) in the paper
                Universes(Back_hole_index,j)=Sorted_universes(White_hole_index,j);
            end
            
            if (size(lb,2)==1)
                %Eq. (3.2) in the paper if the boundaries are all the same
                r2=rand();
                if r2<WEP
                    r3=rand();
                    if r3<0.5
                        Universes(i,j)=Best_universe(1,j)+round(C1*TDR*((ub-lb)*rand+lb));
                        if Universes(i,j)>ub||Universes(i,j)<lb
                            Universes(i,j)=round(rand()*3);
                        end
                    end
                    if r3>0.5
                        Universes(i,j)=Best_universe(1,j)-round(C1*TDR*((ub-lb)*rand+lb));
                        if Universes(i,j)>ub||Universes(i,j)<lb
                            Universes(i,j)=round(rand()*3);
                        end
                    end
                end
            end
%             if i>= dim*9/10
%               temp=Best_universe(1,1:(dim/2));
%               Best_universe(1,1:(dim/2))=Universes(i,1:(dim/2));
%               Universes(i,1:(dim/2))=temp;
%             end
             if i>= dim*3/10                
               Universes(i,:)=Universes(i,:).*brownian(dim)';
             end


        end
                
               %nDNA=InitDNASet(Length,H_target,500);
            if   GCContent(Universes(i,:))~=1||norun(Universes(i,:))==1
                nDNA=[];
                nDNA=NG_InitDNASet(Length,H_target,10000);
                %nDNA=empty_new(nDNA,Length,H_target,count);
                Universes(i,:)=nDNA(1,:);
            end
    end
    for i=1:size(Universes,1)
        DisH=Editconstraint(Universes(i,:),newDNA);
        len1=find(DisH<H_target);
       if isempty(len1) 
           newDNA=[Universes(i,:);newDNA];
            elseif length(len1)==1  
                newDNA(len1,:)=[];
                newDNA=[Universes(i,:);newDNA];
       end
    end
%     if Time==1
%      for i=1:size(newDNA,1)
%          MFEvalue(i,:)=MFE_fun(fliplr(DNAComplement(newDNA(i,:))),DNAComplement(newDNA(i,:)));
%      end
%     
%      MFEAve=sum(MFEvalue)/size(newDNA,1);
%     end
%      px = size(newDNA,1);
%      index2 = 1;
%      temp=[];
%      for q = 1:px
%          if  MFE_fun(fliplr(DNAComplement(newDNA(q,:))),DNAComplement(newDNA(q,:)))>=MFEAve%||  GCContent(DNASet(q,:))~=1
%              temp(index2) =q;
%              index2 = index2 + 1;
%          end
%      end
%      if length(temp)~=0
%          newDNA(temp,:)=[];
%      end
    %Update the convergence curve
    Convergence_curve(Time)=Best_universe_Inflation_rate;
    
    %Print the best universe details after every 50 iterations
    if mod(Time,50)==0
        display(['At iteration ', num2str(Time), ' the best universes fitness is ', num2str(length(newDNA))]);
    end
    Time=Time+1;
end



