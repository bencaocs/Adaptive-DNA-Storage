

Param1=['GTGATCACT'];%[2,0,2,3,0,1,3,1,0]
Param2=['AGTGATCAC'];
b=MFE_fun(Param1,Param2);
c=MFE_fun(DNAComplement(newDNA(1,:)),DNAnum2let(newDNA(1,:)));

Param3=[2,0,3,2,3,2,3];%GTAGAGA²¹ CATCTCT

[Tm,GC]=GCTmBioBox(Param2);
DNAc=DNAComplement(Param3)