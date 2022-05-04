function MFE = MFE_fun(Param1,Param2)
ExeFileName='pairfold.exe';
ExeFilePath=fullfile('.\MultiRNAFold\',ExeFileName);
% Param1=[' ATGCATCGATCAG'];%第一个参数，一定要有' '
% Param2=[' TACGTAGCTAGTC'];

P1=[' ',Param1];%第一个参数，一定要有' '
P2=[' ',Param2];
Cmd=[ExeFilePath,P1,P2];

[a,b]=system(Cmd);
MFE=str2num(b);
end