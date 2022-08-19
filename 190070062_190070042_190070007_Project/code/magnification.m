im_face=(imread("..\images\\jerry_32.png"));
figure;
imshow(im_face);
mag_img1 = imresize(im_face,2,"bilinear");
figure;
imshow(mag_img1);
result=im2double(mag_img1);
sz= size(mag_img1);

for t_step=1:3
    [gx1,gy1]=imgradientxy(result(:,:,1));
    [gx2,gy2]=imgradientxy(result(:,:,2));
    [gx3,gy3]=imgradientxy(result(:,:,3));

    for i=2:sz(1)-1
        for j=2:sz(2)-1
            
            delta_I1=[gx1(i,j);gy1(i,j)];
            delta_I2=[gx2(i,j);gy2(i,j)];
            delta_I3=[gx3(i,j);gy3(i,j)];
            G=(delta_I1*transpose(delta_I1))+(delta_I2*transpose(delta_I2))+(delta_I3*transpose(delta_I3));
%             filter = fspecial('gaussian',2,1);
%             G = imfilter(G,filter);
            [V1,D1] = eig(G);
            [D,P] = sort(diag(D1), 'ascend');
            V = V1(:, P);
            lambda1 = D(1);
            lambda2 = D(2);
            theta1= V(:,1);
            theta2= V(:,2);
            f1=1/(1+lambda1+lambda2);
            f2=1/(sqrt(1+lambda1+lambda2));
            T=f2*(theta2*transpose(theta2))+f1*(theta1*transpose(theta1));
            
           
            Gauss = zeros(3,3);
            sum =0;
            for a= 1:3
                for b = 1:3
                    x = [a-2;b-2];
                    Gauss(a,b) = 1/exp((transpose(x)*(inv(T))*x)/4*t_step);
                        sum = sum + Gauss(a,b);
                end
            end
            Gauss = Gauss/sum;
            subimg = result(i-1:i+1,j-1:j+1,:);
            for k = 1:3
                C = conv2(subimg(:,:,k),Gauss, 'same');
                result(i,j,k) = C(2,2);
            end
            
            
        end
    end
end

figure;

imshow(result,[]);
mag_img2 = imresize(mag_img1,2,"bilinear");
figure;
imshow(mag_img2);
result=im2double(mag_img2);
sz= size(mag_img2);

for t_step=1:3
    [gx1,gy1]=imgradientxy(result(:,:,1));
    [gx2,gy2]=imgradientxy(result(:,:,2));
    [gx3,gy3]=imgradientxy(result(:,:,3));

    for i=2:sz(1)-1
        for j=2:sz(2)-1
            
            delta_I1=[gx1(i,j);gy1(i,j)];
            delta_I2=[gx2(i,j);gy2(i,j)];
            delta_I3=[gx3(i,j);gy3(i,j)];
            G=(delta_I1*transpose(delta_I1))+(delta_I2*transpose(delta_I2))+(delta_I3*transpose(delta_I3));
%             filter = fspecial('gaussian',2,1);
%             G = imfilter(G,filter);
            [V1,D1] = eig(G);
            [D,P] = sort(diag(D1), 'ascend');
            V = V1(:, P);
            lambda1 = D(1);
            lambda2 = D(2);
            theta1= V(:,1);
            theta2= V(:,2);
            f1=1/(1+lambda1+lambda2);
            f2=1/(sqrt(1+lambda1+lambda2));
            T=f2*(theta2*transpose(theta2))+f1*(theta1*transpose(theta1));
            
           
            Gauss = zeros(3,3);
            sum =0;
            for a= 1:3
                for b = 1:3
                    x = [a-2;b-2];
                    Gauss(a,b) = 1/exp((transpose(x)*(inv(T))*x)/4*t_step);
                        sum = sum + Gauss(a,b);
                end
            end
            Gauss = Gauss/sum;
            subimg = result(i-1:i+1,j-1:j+1,:);
            for k = 1:3
                C = conv2(subimg(:,:,k),Gauss, 'same');
                result(i,j,k) = C(2,2);
            end
            
            
        end
    end
end

figure;

imshow(result,[]);