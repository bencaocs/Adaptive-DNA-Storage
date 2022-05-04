function [FE] = FreeEnergy(DNA)
%依Bioinformatics toolbox计算DNAs(整数矩阵)的[GC含量值,Tm值,Hairpin值],DNA长度不能小于8
global SaltValue PrimerconcValue 
SaltValue=1; % Shin2005 盐浓度1M  Bioinformatics Toolbox 默认 0.05
PrimerconcValue=10^(-8); %Shin2005 10nM (10^(-8))DNA分子浓度. 1nM是纳摩尔，1纳摩尔/升=10^(-9)摩尔/升；默认50*10^(-6)
DCopy = DNA;
m=size(DCopy,1);% DNA序列数量
GC=zeros(m,1);
FE=GC;
for i=1:m
    oligoprops=oligoprop(DCopy(i,:),'Salt', SaltValue,'Primerconc', PrimerconcValue);
   FE(i)=round(oligoprops.Thermo(4,3)*10000)/10000; % Tm值，保留小数点后四位
end
FE = sum(FE)/7;