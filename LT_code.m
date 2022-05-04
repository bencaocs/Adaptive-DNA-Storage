% Copyright (c) 2021 by Jae-Won Kim and Jaeho Jeong,
% from Coding and Cryptography Lab (CCL), Department of Electrical and Computer Engineering,
% Seoul National University, South Korea.
% Email address: jaehoj@ccl.snu.ac.kr
% All Rights Reserved.

% Code for encoding method of DNA Fountain enables a robust and efficient storage architecture,
% Erlich and Zielinski, Science (2017).
% input file: image_restart.txt
% output file: Encoded_oligo_%d_%d.txt, LT_K, LT_N

clear;

%%%%%%%%%%%%%%%%%LT parameter + seed parameter%%%%%%%%%%%%%%%%%%%
LT_N = 17200;
LT_K = 15350;
LT_L = 256; % number of bits in payload
LT_seed_bit = 4*8;
LT_c = 0.025;
LT_delta = 0.001;
LT_s = LT_c * sqrt(LT_K) * log(LT_K/LT_delta);


Max_seed_num = power(2,LT_seed_bit)-1;
LT_generator = zeros(LT_K+1,LT_N); % seed at the very last row. messages that were used in encoding for each columns.
LFSR_initial_seed = 42;
RNG_input_seed = 10;

soliton = zeros(1,LT_K);
for i=1:LT_K
    if(i==1)
        soliton(i) = 1/LT_K;
    end
    
    if(i~=1)
        soliton(i) = 1/i/(i-1);
    end
end

tau = zeros(1,LT_K);
KSratio = round(LT_K / LT_s);
for i=1:LT_K
    if(i == KSratio)
        tau(i) = LT_s * log(LT_s/LT_delta) / LT_K;
    end
    
    if(i<KSratio)
        tau(i) = LT_s / LT_K / i;
    end
end

LT_Z = 0;
for i=1:LT_K
    LT_Z = LT_Z + soliton(i) + tau(i);
end

Robust_soliton = zeros(1,LT_K);
for i=1:LT_K
    Robust_soliton(i) = (soliton(i) + tau(i)) / LT_Z;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sum_test = zeros(1,LT_K);

for i=1:LT_K
    if(i == 1)
        sum_test(i) = Robust_soliton(1);
    else
        sum_test(i) = sum_test(i-1) + Robust_soliton(i);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

total_degree_sum = 0;

%%%%%%%%%%%%%%%%%%%%%%RS parameter%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RS_field_size = 256;
RS_parity_num = 2;    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%% Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rng(RNG_input_seed);
Binary_Data_input = zeros(LT_K,LT_L);

input_filename = sprintf('DNAstorage_start.txt');
[FP] = fopen(input_filename,'r');

input_data_save = fscanf(FP,'%d');
fclose(FP);

input_file_bit_size = size(input_data_save,1);

for i=1:LT_K
    for j=1:LT_L
        if((i-1)*LT_L+j<input_file_bit_size + 1)
            Binary_Data_input(i,j)=input_data_save((i-1)*LT_L+j);
        else
            Binary_Data_input(i,j)=randi(2)-1;
        end
            
    end
end


Binary_Data_output = zeros(LT_N,LT_seed_bit+LT_L+8*RS_parity_num);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%% encoding trial %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try_encoding = 0;
success_encoding = 0;
encoding_max_number = Max_seed_num;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%% homopolymer / GC content%%%%%%%%%%%%%%%%%%%%

Max_run_length = 3;
Min_GC_content = 0.45;
MAX_GC_content = 0.55;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







while((try_encoding < encoding_max_number) && (success_encoding<LT_N))
    total_degree_sum_test = 0;
    success_GC = 1;
    success_run = 1;
    
    if(try_encoding == 0)
        temp1_seed = de2bi(LFSR_initial_seed,LT_seed_bit);
        Galois_register = temp1_seed;
    end
    
    temp2_seed = Galois_register;
    for j=1:LT_seed_bit
        if(j~=25 && j~=26 && j~=30 && j~=32)
           Galois_register(j) = temp2_seed(j+1); 
        end

        if(j==32)
           Galois_register(j) = temp2_seed(1);
        end

        if(j==25 || j==26 || j==30)
           Galois_register(j) = mod(temp2_seed(j+1) + temp2_seed(1),2);
        end

    end
    
    %Here, seed is stored at Galois_register. Then, we need to make LT code one symbol.
    
    Galois_register_dec_temp = bi2de(Galois_register);
    rng(Galois_register_dec_temp);
    degree_prop = rand;
    LT_generator_temp = zeros(LT_K+1,1);
    
    
    for j=1:LT_K
        if(j==1)
            if(degree_prop < Robust_soliton(1))
                symbolselection = datasample([1:LT_K],1,'Replace',false);
                LT_generator_temp(symbolselection(1),1)=1;
                total_degree_sum_test = total_degree_sum_test + 1;
            end
        end
        
        if(j~=1)
            if((sum_test(j-1)<=degree_prop) && (degree_prop < sum_test(j)))
                symbolselection = datasample([1:LT_K],j,'Replace',false);
                total_degree_sum_test = total_degree_sum_test + j;
                for m=1:j
                    LT_generator_temp(symbolselection(m),1)=1;
                end
            end
        end
        
        if(degree_prop==1)
            LT_generator_temp(:,1)=1;
            total_degree_sum_test = total_degree_sum_test + LT_K;
        end
    end
    
    LT_generator_temp(LT_K+1,1) = Galois_register_dec_temp;
    
    %Here, LT code one symbol is generated.
    
    One_strand_bit_temp = zeros(1,LT_seed_bit+LT_L+8*RS_parity_num);
    Galois_seed_temp = fliplr(Galois_register);
    One_strand_bit_temp(1:LT_seed_bit) = Galois_seed_temp;
    LT_code_one_symbol_encoding_temp = LT_generator_temp(1:LT_K);
    LT_code_one_symbol_temp = mod(transpose(LT_code_one_symbol_encoding_temp)*Binary_Data_input,2);
    One_strand_bit_temp(LT_seed_bit+1:LT_seed_bit+LT_L) = LT_code_one_symbol_temp;
    
    % Seed and LT code is done. RS encoding is needed.
    % For every 8 bits, correspond to GF(256) and encode them.
    
    
    RS_input_temp = zeros(1,(LT_seed_bit+LT_L)/8);
    eightbit_temp = zeros(1,8);
    
    for i=1:(LT_seed_bit+LT_L)/8
        for j=1:8
            eightbit_temp(j)=One_strand_bit_temp(8*i-j+1);
        end
        RS_input_temp(i) = bi2de(eightbit_temp);
    end
    msg_temp = gf(RS_input_temp,8);
    RS_code_temp = rsenc(msg_temp,(LT_seed_bit+LT_L+8*RS_parity_num)/8,(LT_seed_bit+LT_L)/8);
    
    RS_code_parity_binary_temp = zeros(1,8*RS_parity_num);
    eightbit_temp2 = zeros(1,8);
    galois_test = RS_code_temp.x;
    
    for i=(LT_seed_bit+LT_L)/8+1:(LT_seed_bit+LT_L+8*RS_parity_num)/8
        eightbit_temp2 = de2bi(galois_test(i),8);
        One_strand_bit_temp((i-1)*8+1:8*i)=fliplr(eightbit_temp2);
    end
    
    % Here, strand after RS is stored in One_strand_bit_temp.
    % Changing into ACGT is needed.
    
    for i=1:(LT_seed_bit+LT_L+8*RS_parity_num)/2
        if((One_strand_bit_temp(2*i-1)==0) && (One_strand_bit_temp(2*i)==0))
            One_strand_ACGT_temp(i) = 'A';
        end
        
        if((One_strand_bit_temp(2*i-1)==0) && (One_strand_bit_temp(2*i)==1))
            One_strand_ACGT_temp(i) = 'C';
        end
        
        if((One_strand_bit_temp(2*i-1)==1) && (One_strand_bit_temp(2*i)==0))
            One_strand_ACGT_temp(i) = 'G';
        end
        
        if((One_strand_bit_temp(2*i-1)==1) && (One_strand_bit_temp(2*i)==1))
            One_strand_ACGT_temp(i) = 'T';
        end 
    end
    
    %One_strand_ACGT_temp is generated. size is (1, x)
    %length is (LT_seed_bit+LT_L+8*RS_parity_num)/2.
    
    for i=1:(LT_seed_bit+LT_L+8*RS_parity_num)/2-Max_run_length
        for j=1:Max_run_length
            if(One_strand_ACGT_temp(i)~=One_strand_ACGT_temp(i+j))
                break;
            end
            
            if(j==Max_run_length)
                success_run = 0;
            end
        end
        
        if(success_run == 0)
            break;
        end
    end  
    if(success_run == 0)
        display(['success GC=0 ',One_strand_ACGT_temp,'itar',success_encoding]);
        try_encoding = try_encoding + 1;
        continue;
        
    end
    GC_count = 0;
%     for i=1:(LT_seed_bit+LT_L+8*RS_parity_num)/2
%         if((One_strand_ACGT_temp(i)=='G') || (One_strand_ACGT_temp(i)=='C'))
%             GC_count = GC_count + 1;
%         end
%     end
    for j=1:((LT_seed_bit+LT_L+8*RS_parity_num)/2)-70+1
        GC_count=0;
        for i=j:70+j-1
            if((One_strand_ACGT_temp(i)=='G') || (One_strand_ACGT_temp(i)=='C'))
                GC_count = GC_count + 1;
            end
        end
         GC_ratio = GC_count / 70;
        if((GC_ratio < Min_GC_content) || (GC_ratio > MAX_GC_content))
            success_GC = 0;
            display(['success GC=0 ',One_strand_ACGT_temp,'itar',success_encoding]);
            break;
        else
            success_GC = 1;
            display('success_encoding GC=1 ',num2str(success_encoding));
        end
    end
%     GC_ratio = GC_count / ((LT_seed_bit+LT_L+8*RS_parity_num)/2);
%     
%     if((GC_ratio < Min_GC_content) || (GC_ratio > MAX_GC_content))
%         success_GC = 0;
%     end
    
    if((success_run == 1) && (success_GC == 1))
        success_encoding = success_encoding + 1;
        
        LT_generator(:,success_encoding) = LT_generator_temp;
        Binary_Data_output(success_encoding,:) = One_strand_bit_temp;
        ACGT_Data_output(success_encoding,:) = One_strand_ACGT_temp;
        total_degree_sum = total_degree_sum + total_degree_sum_test;     
    end
    
    try_encoding = try_encoding + 1;
    
end

output_file_name = sprintf('Encoded_oligo_%d_%d.txt', LT_K, LT_N);
[FP_out] = fopen(output_file_name, 'wt');
DNA_len = (LT_seed_bit + LT_L + 8*RS_parity_num)/2;
for i=1:LT_N
    for j=1:DNA_len
        if(Binary_Data_output(i,2*j-1)==0 && Binary_Data_output(i,2*j)==0)
            fprintf(FP_out, 'A');
        elseif(Binary_Data_output(i,2*j-1)==0 && Binary_Data_output(i,2*j)==1)
            fprintf(FP_out, 'C');
        elseif(Binary_Data_output(i,2*j-1)==1 && Binary_Data_output(i,2*j)==0)
            fprintf(FP_out, 'G');
        elseif(Binary_Data_output(i,2*j-1)==1 && Binary_Data_output(i,2*j)==1)
            fprintf(FP_out, 'T');
        end
    end
    fprintf(FP_out, '\n');
end

