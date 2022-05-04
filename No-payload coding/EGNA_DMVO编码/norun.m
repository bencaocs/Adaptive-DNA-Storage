function norun = norun(DNASetsingle)
 for p = 1:(size(DNASetsingle,2)-1)
   if DNASetsingle(p)==DNASetsingle(p+1)
       norun=1;
   break;
   else
       norun=0;
   end
 
 end