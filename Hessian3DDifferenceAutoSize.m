function [Dxx, Dyy, Dzz, Dxy, Dxz, Dyz] = Hessian3DDifferenceAutoSize(Volume,Sigma,window_size) 
 
% �ò�ͬ�Ĵ��ڽ��о��,�ٰ��մ���ֵȡ���ֵ 
 
[l,m,n] = size(Volume); 
 
Dxx = zeros(l,m,n); 
Dyy = zeros(l,m,n); 
Dzz = zeros(l,m,n); 
 
Dxy = zeros(l,m,n); 
Dxz = zeros(l,m,n); 
Dyz = zeros(l,m,n); 
 
%% 
% ���о��ֿ��ܵĴ��ڴ�С1 3 5 7 9 1 1 13 1 5 17 ����Ÿ����ڽ��о�� 
% ���ǵ�����ά��ϸ������Ч��,ʵ����ֻ����1-9 
 
% window_size = 3 
kernel_3_xx =[1,-2,1]; 
kernel_3_yy =[1;-2;1]; 
kernel_3_zz =ones(1,1,3); 
kernel_3_zz(:,:,2)=-2;  %�õ�zά���ϵ�1 -2 1 
 
Dxx3 = imfilter(Volume,kernel_3_xx,'corr','replicate'); 
Dyy3 = imfilter(Volume,kernel_3_yy,'corr','replicate'); 
Dzz3 = imfilter(Volume,kernel_3_zz,'corr','replicate'); 
 
kernel_3_xy =[-1,0,1;0,0,0;1,0,-1]; 
kernel_3_xz = zeros(3,1,3); 
kernel_3_xz(:,:,1)=[-1,0,1]; 
kernel_3_xz(:,:,3)=[1,0,-1]; 
kernel_3_yz =zeros(1,3,3); 
kernel_3_yz(:,:,1)=[-1,0,1]; 
kernel_3_yz(:,:,3)=[1,0,-1]; 
 
Dxy3 = 0.25*imfilter(Volume,kernel_3_xy,'corr','replicate'); 
Dxz3 = 0.25*imfilter(Volume,kernel_3_xz,'corr','replicate'); 
Dyz3 = 0.25*imfilter(Volume,kernel_3_yz,'corr','replicate'); 
 
 
% window_size = 5 
kernel_5_xx =[1,0,-2,0,1]; 
kernel_5_yy =[1;0;-2;0;1]; 
kernel_5_zz =ones(1,1,5); 
kernel_5_zz(:,:,2)=0; 
kernel_5_zz(:,:,3)=-2; 
kernel_5_zz(:,:,4)=0;  %�õ�zά���ϵ�1 0 -2 0 1 
 
Dxx5 = imfilter(Volume,kernel_5_xx,'corr','replicate'); 
Dyy5 = imfilter(Volume,kernel_5_yy,'corr','replicate'); 
Dzz5 = imfilter(Volume,kernel_5_zz,'corr','replicate'); 
 
kernel_5_xy =[-1,0,0,0,1;0,0,0,0,0;0,0,0,0,0;0,0,0,0,0;1,0,0,0,-1]; 
kernel_5_xz = zeros(5,1,5); 
kernel_5_xz(:,:,1)=[-1,0,0,0,1]; 
kernel_5_xz(:,:,5)=[1,0,0,0,-1]; 
kernel_5_yz =zeros(1,5,5); 
kernel_5_yz(:,:,1)=[-1,0,0,0,1]; 
kernel_5_yz(:,:,5)=[1,0,0,0,-1]; 
 
Dxy5 = 0.25*imfilter(Volume,kernel_5_xy,'corr','replicate'); 
Dxz5 = 0.25*imfilter(Volume,kernel_5_xz,'corr','replicate'); 
Dyz5 = 0.25*imfilter(Volume,kernel_5_yz,'corr','replicate'); 
 
 
% window_size = 7 
kernel_7_xx =[1,0,0,-2,0,0,1]; 
kernel_7_yy =[1;0;0;-2;0;0;1]; 
kernel_7_zz =zeros(1,1,7); 
kernel_7_zz(:,:,1)=1; 
kernel_7_zz(:,:,4)=-2; 
kernel_7_zz(:,:,7)=-1;  %�õ�zά���ϵ�1 0 0 -2 0 0 1 
 
Dxx7 = imfilter(Volume,kernel_7_xx,'corr','replicate'); 
Dyy7 = imfilter(Volume,kernel_7_yy,'corr','replicate'); 
Dzz7 = imfilter(Volume,kernel_7_zz,'corr','replicate'); 
 
kernel_7_xy =[-1,0,0,0,0,0,1; 
    0,0,0,0,0,0,0; 
    0,0,0,0,0,0,0; 
    0,0,0,0,0,0,0; 
    0,0,0,0,0,0,0; 
    0,0,0,0,0,0,0 
    1,0,0,0,0,0,-1]; 
kernel_7_xz = zeros(7,1,7); 
kernel_7_xz(:,:,1)=[-1,0,0,0,0,0,1]; 
kernel_7_xz(:,:,7)=[1,0,0,0,0,0,-1]; 
kernel_7_yz =zeros(1,7,7); 
kernel_7_yz(:,:,1)=[-1,0,0,0,0,0,1]; 
kernel_7_yz(:,:,7)=[1,0,0,0,0,0,-1]; 
 
Dxy7 = 0.25*imfilter(Volume,kernel_7_xy,'corr','replicate'); 
Dxz7 = 0.25*imfilter(Volume,kernel_7_xz,'corr','replicate'); 
Dyz7 = 0.25*imfilter(Volume,kernel_7_yz,'corr','replicate'); 
 
 
% window_size = 9 
kernel_9_xx =[1,0,0,0,-2,0,0,0,1]; 
kernel_9_yy =[1;0;0;0;-2;0;0;0;1]; 
kernel_9_zz =zeros(1,1,9); 
kernel_9_zz(:,:,1)=1; 
kernel_9_zz(:,:,5)=-2; 
kernel_9_zz(:,:,9)=-1;  %�õ�zά���ϵ�1 0 0 0 -2 0 0 0 1 
 
Dxx9 = imfilter(Volume,kernel_9_xx,'corr','replicate'); 
Dyy9 = imfilter(Volume,kernel_9_yy,'corr','replicate'); 
Dzz9 = imfilter(Volume,kernel_9_zz,'corr','replicate'); 
 
kernel_9_xy =zeros(9,9,1); 
kernel_9_xy(1,:,:)=[-1,0,0,0,0,0,0,0,1]; 
kernel_9_xy(9,:,:)=[1,0,0,0,0,0,0,0,-1]; 
kernel_9_xz = zeros(9,1,9); 
kernel_9_xz(:,:,1)=[-1,0,0,0,0,0,0,0,1]; 
kernel_9_xz(:,:,9)=[1,0,0,0,0,0,0,0,-1]; 
kernel_9_yz =zeros(1,9,9); 
kernel_9_yz(:,:,1)=[-1,0,0,0,0,0,0,0,1]; 
kernel_9_yz(:,:,9)=[1,0,0,0,0,0,0,0,-1]; 
 
Dxy9 = 0.25*imfilter(Volume,kernel_9_xy,'corr','replicate'); 
Dxz9 = 0.25*imfilter(Volume,kernel_9_xz,'corr','replicate'); 
Dyz9 = 0.25*imfilter(Volume,kernel_9_yz,'corr','replicate'); 
 
 
%% ���մ���ֵȡ���ֵ ��������ά��ϸ����Ϊ���ֱ��Ϊ9   
 
names = ['Dxx','Dyy','Dzz','Dxy','Dxz','Dyz']; 
 
for i = 1:8     
    for j = 1:6 
        if(i<4) 
            a = strcat(names(j*3-2:j*3),'(window_size==',num2str(i),')= ',names(j*3-2:j*3),num2str(2*i+1),'(window_size==',num2str(i),');');       
        else 
            a = strcat(names(j*3-2:j*3),'(window_size==',num2str(i),')= ',names(j*3-2:j*3),num2str(9),'(window_size==',num2str(i),');'); 
        end 
        eval(a); 
    end 
end 