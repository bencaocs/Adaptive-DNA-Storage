function [Tm,GC] = GCTmBioBox(DNAs )
%��Bioinformatics toolbox����DNAs(��������)��[GC����ֵ,Tmֵ,Hairpinֵ],DNA���Ȳ���С��8
global SaltValue PrimerconcValue 
SaltValue=1; % Shin2005 ��Ũ��1M  Bioinformatics Toolbox Ĭ�� 0.05
PrimerconcValue=10^(-8); %Shin2005 10nM (10^(-8))DNA����Ũ��. 1nM����Ħ����1��Ħ��/��=10^(-9)Ħ��/����Ĭ��50*10^(-6)
m=size(DNAs,1);% DNA��������
GC=zeros(m,1);
Tm=GC;
for i=1:m 
    oligoprops=oligoprop(DNAs(i,:),'Salt', SaltValue,'Primerconc', PrimerconcValue);
%     oligoprops=oligoprop(DNAs(i,:));
    GC(i)=round(oligoprops.GC);            %%%% GC-Content��ȡ��
%     Tm(i)=oligoprops.Tm(5); %%%% Tmֵ
   Tm(i)=round(oligoprops.Tm(5)*10000)/10000; %%%% Tmֵ������С�������λ
end



