function [ finalImage ] = foregroundStars( inputImage,displayImages,slideshow,smoothImage)

%DislpayImages: Boolean, if true displays images of each calculation, if
%   false displays only final image
%slideshow: Boolean, if true waits for keyboard input before displaying the
%   next image, if false displays images as soon as they are calculated
%SmoothImage, intager from [0,3]
%   >1 adds randomness over star replacements
%   ==2 smooths images using average from 4 pixels (up, down, left, right)
%   ==3 smooths image using average from 8 surrounding pixels


maximaSearchRadius=2;
minimumStarBrightness=0;
differenceFromAve=8;
galaxyAveBrightness=20;
%Finding info about the image
[sizeX,sizeY]=size(inputImage);
[maximum]=max(inputImage(:));
[galCentreX,galCentreY]=find(inputImage==maximum);
minimum=min(inputImage(:));

if displayImages
    if slideshow
        pause
    end
    figure('Name','Original image'), image(inputImage);
end

%Getting Log of image and limiting it from 1 to 64, done twice to look
%nicer
inputImage=inputImage-minimum;
b=log(inputImage+1.5);
b=b-min(b(:));
b=b/max(b(:))*64;

b=log(b+1.5);
b=b-min(b(:));
b=b/max(b(:))*64;


if displayImages
    if slideshow
        pause
    end
    figure('Name','Log of image between 0 and 64'), image(b);
end

finalImage=b;

%Finding the average brightness in a ring around the galaxy in increasing
%distance

total(1:ceil(2*sqrt((sizeX-galCentreX)^2+(sizeY-galCentreY)^2)))=0;
number=total;
for x=1:sizeX
    for y=1:sizeY
        total(1+ceil(sqrt((x-galCentreX)^2+(y-galCentreY)^2)))=total(1+ceil(sqrt((x-galCentreX)^2+(y-galCentreY)^2)))+b(x,y);
        number(1+ceil(sqrt((x-galCentreX)^2+(y-galCentreY)^2)))=number(1+ceil(sqrt((x-galCentreX)^2+(y-galCentreY)^2)))+1;
    end
end
aveRadBright=total(:)./number(:);

if displayImages
    if slideshow
        pause
    end
    figure('Name','Average Radial Brightness around Galaxy'), plot(aveRadBright);
end

%Estimate the radius of the galaxy 
galRadius=endOfGalaxy(aveRadBright,galaxyAveBrightness);


%Average out the image of the galaxy, this makes the image smoother, and
%makes it easier to find local peaks (stars)
averageBright(1:sizeX,1:sizeY)=0;
for i=2:sizeX-1
    for j=2:sizeY-1
        averageBright(i,j)=(b(i,j)+b(i+1,j)+b(i+1,j-1)+b(i+1,j+1)+b(i-1,j)+b(i-1,j-1)+b(i-1,j+1)+b(i,j+1)+b(i,j-1))/9;
    end
end

%Look for any local maxima, it must be the brighteset point in
% "maximaSearchRadius" around each point.
%it must also be brighter (by "differenceFromAve") than the average brightness that distance from
%the galaxy. And it must not be inside the galaxy's centre ("galRadius")
localMaxima(1:sizeX,1:sizeY)=0;
for i=1+maximaSearchRadius:sizeX-maximaSearchRadius
    for j=1+maximaSearchRadius:sizeY-maximaSearchRadius
        if isLocalMax(averageBright,i,j,maximaSearchRadius,minimumStarBrightness);
            if ((b(i,j)-aveRadBright(1+ceil(sqrt((i-galCentreX)^2+(j-galCentreY)^2))))>differenceFromAve)&&(((sqrt((i-galCentreX)^2+(j-galCentreY)^2)))>galRadius)
                localMaxima(i,j)=64;
            end
        end
    end
end


if displayImages
    if slideshow
        pause
    end
    figure('Name','Local Maxima'), image(localMaxima)
    colormap(hot)
end

if displayImages
    c=imfuse(localMaxima,b);
        if slideshow
            pause
        end
    figure('Name','Star Position'), imshow(c,'InitialMagnification',200);
end

%Find the radius of each star, determined by the average brightness around
%that area.
rad=b*0;
for i=1:sizeX
   for j=1:sizeY
        if localMaxima(i,j)==64
            rad(i,j)=findStarRad(b,aveRadBright,i,j,galCentreX,galCentreY);
        end
   end
end

%Makes an array of the image where 0 means that that pixel is not in any stars radius
newRad=spreadRadius(rad);

%Replaces all pixels of stars with the average brightness of that radius around the galaxy
newImage=replaceStars(b,newRad,aveRadBright,galCentreX,galCentreY);

%if smoothImage>1 makes star replacements look more natural
if slideshow
    for i=0:3
        tempImage=makeRepNicer(newImage,newRad,i);
        pause
        figure('Name',strcat('Smoothed image Number ',int2str(i))), image(tempImage);
    end
end
   newImage=makeRepNicer(newImage,newRad,smoothImage);

if slideshow
    pause
end
figure('Name','Finished Product'), image(newImage);

finalImage=newImage;
end



function result=isLocalMax(table,x,y,r,minBright)
result=true;
for i=-r:r
    for j=-r:r
        if table(x,y)<table(x+i,y+j)||table(x,y)<minBright
            result=false;
        end
    end
end
end



function galRad=endOfGalaxy(aveRadBrightness,galBrightness)
    [length,~]=size(aveRadBrightness(:));
    found=false;
    i=1;
    while (found==0)&&(i<=length)
        if aveRadBrightness(i)<galBrightness
            found=true;
            galRad=i-1;
        end
        i=i+1;
    end
end

function starRad=findStarRad(table,aveRadBrightness,sX,sY,galX,galY) 
    distanceFromGal=ceil(sqrt((sX-galX)^2+(sY-galY)^2));
    aveBright=64;
    currentRadius=1;
    [sizeX,sizeY]=size(table);
    while aveBright>5+(aveRadBrightness(distanceFromGal))
        total=0;
             for theta=1:360
                if (sX-currentRadius>=1&&sY-currentRadius>=1&&sX+currentRadius<=sizeX&&sY+currentRadius<=sizeY)
                     total=total+table(floor(sX+currentRadius*cosd(theta)),floor(sY+currentRadius*sind(theta)));
                else
                   total=total+total/theta;
                 end
             end
         aveBright=total/360;
         currentRadius=currentRadius+1;
    end
    starRad=currentRadius;
end


function result=spreadRadius(rad)
[sizeX,sizeY]=size(rad);
result=rad;
    for x=1:sizeX
        for y=1:sizeY
            if rad(x,y)>0
                radius=rad(x,y);
                for tempRad=1:radius
                    for i=-tempRad:tempRad
                        for j=-tempRad:tempRad
                            if (i^2+j^2<=tempRad^2)
                                if((x+i<=sizeX)&&(x+i>=1)&&(y+j<=sizeY)&&(y+j>=1))
                                    result(x+i,y+j)=rad(x,y);
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end


function newTable=replaceStars(image,stars,aveRadbrightness,galCentreX,galCentreY)
    newTable=image;
    [sizeX,sizeY]=size(image);
    for i=1:sizeX
        for j=1:sizeY
            if(stars(i,j)>0)
                newTable(i,j)=(aveRadbrightness(1+ceil(sqrt((i-galCentreX)^2+(j-galCentreY)^2))));
            end
        end
    end
end


function newIm=makeRepNicer(image,star,smoothImage)
if smoothImage>0
    [sizeX,sizeY]=size(image);

    for i=1:sizeX+2
        for j=1:sizeY+2
             randomM(i,j)=rand();  
        end
    end


    for i=2:sizeX+1
        for j=2:sizeY+1
            newRand(i-1,j-1)=(randomM(i,j)+randomM(i+1,j)+randomM(i+1,j-1)+randomM(i+1,j+1)+randomM(i-1,j)+randomM(i-1,j-1)+randomM(i-1,j+1)+randomM(i,j+1)+randomM(i,j-1))/9;
        end
    end

        newRand=newRand-min(newRand(:));
        newRand=newRand/max(newRand(:))*64;

    for i=1:sizeX
        for j=1:sizeY
            if (star(i,j)>0)
                image(i,j)=(0.9*image(i,j)+0.07*newRand(i,j));
            end
        end
    end


    for i=2:sizeX-1
        for j=2:sizeY-1
            if smoothImage==2
                newIm(i-1,j-1)=(image(i,j)+image(i+1,j)+image(i-1,j)+image(i,j+1)+image(i,j-1))/4;
            elseif smoothImage==3
                newIm(i-1,j-1)=(image(i,j)+image(i+1,j)+image(i+1,j-1)+image(i+1,j+1)+image(i-1,j)+image(i-1,j-1)+image(i-1,j+1)+image(i,j+1)+image(i,j-1))/9;
            else
                newIm=image;
            end
        end
    end
else
    newIm=image;
end
    
end






