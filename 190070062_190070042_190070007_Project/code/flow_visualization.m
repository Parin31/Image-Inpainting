im_face=imread("..\\images\\face.jpeg");
figure;
imshow(im_face);
sz=size(im_face);
size(im_face(:,:,1));
result=im_face;
G=zeros(sz);
F1=[0;3];
F2=[3;0];
t=100;
for k=1:sz(3)
    for t_step=1:t
        [gx,gy]=imgradientxy(result(:,:,k),'sobel');
        [gxx,gxy]=imgradientxy(gx,'sobel');
        [gyx,gyy]=imgradientxy(gy,'sobel');
        for i=1:sz(1)
            for j=1:sz(2)
                H = [gxx(i,j),gxy(i,j);gyx(i,j),gyy(i,j)];
                if j>sz(2)/2
                    T=(F1*transpose(F1))./norm(F1);
                else
                    T=(F2*transpose(F2))./norm(F2);
                end
                result(i,j,k)=result(i,j,k)+0.01*trace(T*H);
            end
        end
    end
end
figure;
imshow(uint8(result));

            