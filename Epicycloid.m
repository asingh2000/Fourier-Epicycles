% Handling the input image
Input='Pi-symbol.svg.png'
img=imread(Input);
% img=imrotate(img,-90,'bilinear');
% imshow(img);

% Binarize the image
BI=im2bw(img);
BI=~BI;
imshow(BI);

% Getting the pixels of boundary figure
dim=size(BI);
col=round(dim(2)/2);
row=min(find(BI(:,col)));
boundary=bwtraceboundary(BI,[row,col],'N');
center=round((min(boundary)+max(boundary))./2);
pts= [boundary(:,2)-center(2) -boundary(:,1)+center(1)];
x=pts(:,1);
y=pts(:,2);
% plot(x,y);

% Obtaining the Fourier Coefficients
m=50;
z=x+1i.*y;
N=numel(z);
Cn=zeros(2*m+1,1);
for j=1:2*m+1
    cf=0;
    for k=1:N
        cf=cf+z(k)*exp(-2i*(pi/N)*(j-m-1)*k);
    end
    Cn(j)=cf/N;
end

% Plotting the Epicycloid
Output='Contour.avi'
w=VideoWriter(Output);
w.FrameRate=5;
open(w);

t=linspace(0,2*pi,100);
r=zeros(2*m+1,1);
for ind=1:length(t)
    T=t(ind);
    for k=1:2*m+1
        r(k)=Cn(k)*exp(1i*(k-m-1)*T);
    end
    a1=fliplr(r(1:m));
    a2=r(m+2:2*m+1);
    a3=[a2;a1];
    a3=a3(:).';
    r=[r(m+1) a3];
    r=cumsum(r);
    R=abs(r(2:2*m+1)-r(1:2*m));
    
    % Plotting circles
    for cntr=1:2*m
        x0=real(r(cntr));
        y0=imag(r(cntr));
        rad=R(cntr);
%         plot(x0+rad*cos(t),y0+rad*sin(t));
%         axis([min(x)-40 max(x)+40 min(y)-25 max(y)+25]);
%         hold on
    end
    
    % Plotting Radius
    plot(real(r),imag(r),'--.');
%     hold off
    
    % Plotting contour
    if ind~=1
        plot([real(rnew),real(r(2*m+1))],[imag(rnew),imag(r(2*m+1))],'g','Linewidth',4);
        axis([min(x)-40 max(x)+40 min(y)-25 max(y)+25]);
        writeVideo(w,getframe(gcf));
        hold on
    end
    rnew=r(2*m+1);
end
close(w);





