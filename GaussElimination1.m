function  [result] = GaussElimination1(A,B)
%  inputs:
%        A:系数矩阵，为n*n维方阵
%        B :载荷矩阵，为n*1维矩阵,输入B时要求时单列
%       说白了，就是AX=B    
%  outputs：
%         result:计算结果向量,为n*1为矩阵

% 1、判断输入的矩阵维度是不是满足要求 n*n
[A_Row,A_Column] = size(A);
[B_Row,~] = size(B);

% 2、初始化result向量
result = zeros(B_Row,1); %和B同行数
%3、判断输入的维度是不是满足要求
if (A_Row ~= A_Column) ||(A_Row ~= B_Row) % A要行列数等，还要和B同行数，有可能唯一解
    print('输入错误！');
else       
    for L =1:A_Row-1        %需要进行(A_Row-1)次调整
        if A(L,L)==0
            print('主对角线元素不能为0');
            break;
        else
            for k = L+1:A_Row    % 循环计算第l+1行到最后一行
                L_lk = A(k,L)*A(L,L); % 每次更新 ，为图示的C(K,L)/C(L,L)
                for j =L+1:A_Column  % 更新每一行从第i个元素开始后的所有元素
                    A(k,j)=  A(k,j)-L_lk*A(L,j);
                end
                %需要对B也进行改动
                B(k)=B(k) - L_lk *B(L);
            end %for
        end% 第二个if
    end% for
    
    %倒序求解过程，获得结果,倒三角形
    result(B_Row) = B(B_Row)/A(A_Row,A_Column);
    for k=A_Row-1:-1:1
        sumTemp =0;
        for j = k+1:A_Row
            sumTemp =sumTemp+A(k,j)*result(j);
        end
        result(k) =( B(k)-sumTemp)/A(k,k);
    end
end % if