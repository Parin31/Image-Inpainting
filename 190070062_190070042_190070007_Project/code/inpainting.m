img = imread("..\\images\\parrot.PNG");
img_og=img;
inp_mask = imread("..\\images\\mask.PNG");
inp_mask=inp_mask(1:360,1:359,:);
inp_mask_bw=rgb2gray(inp_mask);
img=cast(img,'double');
mask = imbinarize(inp_mask_bw, 0.5);
nbhd = 20;
[m,n,no_of_channels] = size(img);
result = img;

for t_step=1:20
    [gx1,gy1]=imgradientxy(img(:,:,1));
    [gx2,gy2]=imgradientxy(img(:,:,2));
    [gx3,gy3]=imgradientxy(img(:,:,3));
    padded_img = padarray(img,[nbhd,nbhd]);
    %iterating over the non-padded part of the image
    for i=1:m
        for j=1:n
            if mask(i,j) == 1
                delta_I1=[gx1(i,j);gy1(i,j)];
                delta_I2=[gx2(i,j);gy2(i,j)];
                delta_I3=[gx3(i,j);gy3(i,j)];
                G=(delta_I1*transpose(delta_I1))+(delta_I2*transpose(delta_I2))+(delta_I3*transpose(delta_I3));
                G_sigma = imgaussfilt(G,1);
                [V1,D1] = eig(G_sigma);
                [D,perm] = sort(diag(D1), 'ascend');
                V = V1(:, perm);
                lambda1 = D(1);
                lambda2 = D(2);
                theta1= V(:,1);
                theta2= V(:,2);
                f1=1/(1+lambda1+lambda2);
                f2=1/(sqrt(1+lambda1+lambda2));
                T=f2*(theta2*transpose(theta2))+f1*(theta1*transpose(theta1));
                [A,B]=meshgrid(-nbhd:nbhd,-nbhd:nbhd);
                gauss=(1/(4*pi))*exp(-((A.*A)*T(1,1)+(A.*B)*(T(1,2)+T(2,1))+B.*B*T(2,2))/(4));
                S = sum(gauss,'all');
                gauss=gauss./S;
                im_subarr=padded_img(i:i+2*nbhd,j:j+2*nbhd,:);
                for k=1:3
                    conv_res=conv2(im_subarr(:,:,k),gauss,'same');
                    result(i,j,k) = conv_res(1+nbhd,1+nbhd);
                end
            end
        end
    end
    img=result;
end
img=cast(img,'uint8');
result=cast(result,'uint8');
figure;
imshow(img_og);
figure;
imshow(mask);
figure;
imshow(img);  