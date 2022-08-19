im_baby=imread("..\\images\\baby.jpg");
sz=size(im_baby);
mask=cast(ones(sz),'uint8');
nbhd=3;
for i=1:sz(1)
    for j=1:sz(2)
        if (rem(i,2)==0 && rem(j,2) == 1)
            mask(i,j,1)=0;
            mask(i,j,2)=0;
            mask(i,j,3)=0;
        elseif(rem(j,2)==0 && rem(i,2)==1)
            mask(i,j,1)=0;
            mask(i,j,2)=0;
            mask(i,j,3)=0;
        end   
    end
end
im_noisy=cast(im_baby.*mask,'double');
imshow(cast(im_noisy,'uint8'));
result=im_noisy;
for t_step=1:20
    [gx1,gy1]=imgradientxy(im_noisy(:,:,1));
    [gx2,gy2]=imgradientxy(im_noisy(:,:,2));
    [gx3,gy3]=imgradientxy(im_noisy(:,:,3));
    im_noisy_pad=padarray(im_noisy,[nbhd,nbhd]);
    for i=1:sz(1)
        for j=1:sz(2)
            if mask(i,j,1) == 0
                delta_I1=[gx1(i,j);gy1(i,j)];
                delta_I2=[gx2(i,j);gy2(i,j)];
                delta_I3=[gx3(i,j);gy3(i,j)];
                G=(delta_I1*transpose(delta_I1))+(delta_I2*transpose(delta_I2))+(delta_I3*transpose(delta_I3));
                G = imgaussfilt(G,1);
                if(i==1 && j==1)
                    disp(G);
                end
                [V1,D1] = eig(G);
                [D,perm] = sort(diag(D1), 'ascend');
                V = V1(:, perm);
                lambda1 = D(1);
                lambda2 = D(2);
                theta1= V(:,1);
                theta2= V(:,2);
                f1=1/(1+lambda1+lambda2);
                f2=1/(sqrt(1+lambda1+lambda2));
                T=f2*(theta2*transpose(theta2))+f1*(theta1*transpose(theta1));
                T_inv=inv(T);
                [A,B]=meshgrid(-nbhd:nbhd,-nbhd:nbhd);
                gauss=(1/(4*pi))*exp(-((A.*A)*T(1,1)+(A.*B)*(T(1,2)+T(2,1))+B.*B*T(2,2))/(4));
                S = sum(gauss,'all');
                gauss=gauss./S;
                im_subarr=im_noisy_pad(i:i+2*nbhd,j:j+2*nbhd,:);
                for k=1:3
                    conv_res=conv2(im_subarr(:,:,k),gauss,'same');
                    result(i,j,k)=conv_res(1+nbhd,1+nbhd);
                end
            end
        end
    end
    im_noisy=result;
end
figure;
imshow(im_baby);
figure;
imshow(cast(result,'uint8'));

            