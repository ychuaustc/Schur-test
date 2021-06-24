function  [result] = GaussElimination1(A,B)
%  inputs:
%        A:ϵ������Ϊn*nά����
%        B :�غɾ���Ϊn*1ά����,����BʱҪ��ʱ����
%       ˵���ˣ�����AX=B    
%  outputs��
%         result:����������,Ϊn*1Ϊ����

% 1���ж�����ľ���ά���ǲ�������Ҫ�� n*n
[A_Row,A_Column] = size(A);
[B_Row,~] = size(B);

% 2����ʼ��result����
result = zeros(B_Row,1); %��Bͬ����
%3���ж������ά���ǲ�������Ҫ��
if (A_Row ~= A_Column) ||(A_Row ~= B_Row) % AҪ�������ȣ���Ҫ��Bͬ�������п���Ψһ��
    print('�������');
else       
    for L =1:A_Row-1        %��Ҫ����(A_Row-1)�ε���
        if A(L,L)==0
            print('���Խ���Ԫ�ز���Ϊ0');
            break;
        else
            for k = L+1:A_Row    % ѭ�������l+1�е����һ��
                L_lk = A(k,L)*A(L,L); % ÿ�θ��� ��Ϊͼʾ��C(K,L)/C(L,L)
                for j =L+1:A_Column  % ����ÿһ�дӵ�i��Ԫ�ؿ�ʼ�������Ԫ��
                    A(k,j)=  A(k,j)-L_lk*A(L,j);
                end
                %��Ҫ��BҲ���иĶ�
                B(k)=B(k) - L_lk *B(L);
            end %for
        end% �ڶ���if
    end% for
    
    %���������̣���ý��,��������
    result(B_Row) = B(B_Row)/A(A_Row,A_Column);
    for k=A_Row-1:-1:1
        sumTemp =0;
        for j = k+1:A_Row
            sumTemp =sumTemp+A(k,j)*result(j);
        end
        result(k) =( B(k)-sumTemp)/A(k,k);
    end
end % if