function result = Aco(DNASet,len,H_target)
t_min=0.01;
n = size(DNASet,1);                           
m = 20;                                      %蚂蚁的数量
alpha = 1;                                  %信息素浓度重要程度因子
beta = 1;                                   %启发函数重要程度因子
Q = 1.5;                                      %常系数
rho = 0.01;                                  %信息素挥发因子
Tau = ones(n,1);                            %信息素矩阵
iter = 1;                                   %迭代次数初值
iter_max =500 ;
if m>n
    m=n;
end
HammingDist = DistHammingMatrix(DNASet); 
%  
 
BestMaxSet=[];
while  iter <= iter_max%iter迭代次数初值
      individual=[];
      DNA=[]; 
      temp = randperm(n);
    for i=1:m  
        individual(i).n=i;               %individual结构体
        individual(i).p=[];
        DNA(i).p = [temp(i)];
        A=find(HammingDist(temp(i),:)>=H_target);%hammingdist  是DNAset的距离对称矩阵
       
        individual(i).p = A;%tenp(i)行如果hamming distance >H_target 存入individua(i).p
    end
     for i = 1:m                      %m是蚂蚁数量
         while 1
            if length([individual(i).p]) ==0
                break;
            end
             del = [];
         for  k= 1:length([individual(i).p])
             P(k) = Tau(individual(i).p(k))^alpha;%信息素浓度重要程度因子alpha= 1;%信息素矩阵Tau = ones(n,1);
         end
         P = P / sum(P);%诡异的轮盘赌
         Pc = cumsum(P);%cumsum函数通常用于计算一个数组各行累加到下一位置
         target_index = find(Pc >= rand);%大于随机数(0-1)的Pc
         target =individual(i).p(target_index(1));%存 大于随机数的第一个位置
         individual(i).p(target_index(1)) = [];%删除         individual(i).p()一直表示编码距离   
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
           Tau(DNA(i).p(j)) = Tau(DNA(i).p(j)) - t_min*length(BestMaxSet)/length([DNA(i).p]);%信息素矩阵Tau = ones(n,1);
         end
         Tau(find(Tau(:,1)>=6))=6;
        Tau(find(Tau(:,1)<=0.01))=0.01;
     end
     Length =  zeros(m,1);
     for i = 1:m
        Length(i) = length([DNA(i).p]);
     end
     [LocalMax_Length,LocalMax_index] = max(Length);%[x,y] x是最大值 y是最大下标
     LocalMaxSet = DNA(LocalMax_index).p;
     Delta_Tau = zeros(n,1);
     for i=1:LocalMax_Length
         Delta_Tau(LocalMaxSet(i)) = Delta_Tau(LocalMaxSet(i)) + Q/(length(BestMaxSet)-LocalMax_Length+1);     
     end  %常系数Q=1.5
     Tau = (1-rho)*Tau + Delta_Tau;%信息素挥发因子rho=0.01
     Tau(find(Tau(:,1)>=6))=6;
     Tau(find(Tau(:,1)<=0.01))=0.01;
     iter = iter + 1;
end
px = length(BestMaxSet);
for i=1:px
    result(i,:) = DNASet(BestMaxSet(i),:);
end
end
