function [seq1,seq2,seq3]=isonepf( A)
% A= find(ispf_cxutemp(c,:,u))  % your array of values
A(end+1)=2   % adds new endpoint to very end of A so code picks up end of last group of consecutive values
I_1=find(diff(A)~=1);  % finds where sequences of consecutive numbers end
[m,n]=size(I_1);   % finds dimensions of I_1 i.e. how many sequences of consecutive numbers you have
startpoint=1;    % sets start index at the first value in your array
seq=cell(1,n)  % had to preallocate because without, it only saved last iteration of the for loop below
seq2=[]    ;       % used n because this array is a row vector
seq3=[];
for i=1:n
    End_Idx=I_1(i);   %set end index
    seq{i}=A(startpoint:End_Idx);  %finds sequences of consecutive numbers and assigns to cell array
    startpoint=End_Idx+1;   %update start index for the next consecutive sequence
end


    seq1=cell2mat(seq(1,1));
    if size(seq,2)==2
    seq2=cell2mat(seq(1,2));
    elseif size(seq,2)==3
    seq3=cell2mat(seq(1,3));
    end
end