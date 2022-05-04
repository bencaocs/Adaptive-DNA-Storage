function result = Aco(DNASet,len,H_target)
t_min=0.01;
n = size(DNASet,1);                           
m = 20;                                      %���ϵ�����
alpha = 1;                                  %��Ϣ��Ũ����Ҫ�̶�����
beta = 1;                                   %����������Ҫ�̶�����
Q = 1.5;                                      %��ϵ��
rho = 0.01;                                  %��Ϣ�ػӷ�����
Tau = ones(n,1);                            %��Ϣ�ؾ���
iter = 1;                                   %����������ֵ
iter_max =500 ;
if m>n
    m=n;
end
HammingDist = DistHammingMatrix(DNASet); 
%  
 
BestMaxSet=[];
while  iter <= iter_max%iter����������ֵ
      individual=[];
      DNA=[]; 
      temp = randperm(n);
    for i=1:m  
        individual(i).n=i;               %individual�ṹ��
        individual(i).p=[];
        DNA(i).p = [temp(i)];
        A=find(HammingDist(temp(i),:)>=H_target);%hammingdist  ��DNAset�ľ���Գƾ���
       
        individual(i).p = A;%tenp(i)�����hamming distance >H_target ����individua(i).p
    end
     for i = 1:m                      %m����������
         while 1
            if length([individual(i).p]) ==0
                break;
            end
             del = [];
         for  k= 1:length([individual(i).p])
             P(k) = Tau(individual(i).p(k))^alpha;%��Ϣ��Ũ����Ҫ�̶�����alpha= 1;%��Ϣ�ؾ���Tau = ones(n,1);
         end
         P = P / sum(P);%��������̶�
         Pc = cumsum(P);%cumsum����ͨ�����ڼ���һ����������ۼӵ���һλ��
         target_index = find(Pc >= rand);%���������(0-1)��Pc
         target =individual(i).p(target_index(1));%�� ����������ĵ�һ��λ��
         individual(i).p(target_index(1)) = [];%ɾ��         individual(i).p()һֱ��ʾ�������   
         Pc=[];
         P=[];
         DNA(i).p = [DNA(i).p target];
         index = 1;
         for  k= 1:length([individual(i).p])
               if HammingDist(individual(i).p(k),target) < H_target 
                 del(index) = k;
                 index = index + 1;
               end  
         end
         individual(i).p(del)=[];
         end
       
         if length([DNA(i).p]) > length(BestMaxSet)
             BestMaxSet = [DNA(i).p];
         end
         for j = 1:length([DNA(i).p])
           Tau(DNA(i).p(j)) = Tau(DNA(i).p(j)) - t_min*length(BestMaxSet)/length([DNA(i).p]);%��Ϣ�ؾ���Tau = ones(n,1);
         end
         Tau(find(Tau(:,1)>=6))=6;
        Tau(find(Tau(:,1)<=0.01))=0.01;
     end
     Length =  zeros(m,1);
     for i = 1:m
        Length(i) = length([DNA(i).p]);
     end
     [LocalMax_Length,LocalMax_index] = max(Length);%[x,y] x�����ֵ y������±�
     LocalMaxSet = DNA(LocalMax_index).p;
     Delta_Tau = zeros(n,1);
     for i=1:LocalMax_Length
         Delta_Tau(LocalMaxSet(i)) = Delta_Tau(LocalMaxSet(i)) + Q/(length(BestMaxSet)-LocalMax_Length+1);     
     end  %��ϵ��Q=1.5
     Tau = (1-rho)*Tau + Delta_Tau;%��Ϣ�ػӷ�����rho=0.01
     Tau(find(Tau(:,1)>=6))=6;
     Tau(find(Tau(:,1)<=0.01))=0.01;
     iter = iter + 1;
end
px = length(BestMaxSet);
for i=1:px
    result(i,:) = DNASet(BestMaxSet(i),:);
end
end
