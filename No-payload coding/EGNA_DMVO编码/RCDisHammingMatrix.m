function RCHammingDistMatrix= RCDisHammingMatrix(IntDNAMatrix)
 IntDNA = IntDNAMatrix;
 RIntDNAMatrix = fliplr(3-IntDNAMatrix);
 RCHammingDistMatrix = pdist2(IntDNA,RIntDNAMatrix,'hamming')*size(IntDNAMatrix,2);
end
