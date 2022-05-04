function nDNA = empty_new(a,Length,H_target,count)
 if isempty(a)
        nDNA=InitDNASet(Length,H_target,count);
        empty_new(nDNA,Length,H_target,count);
 end