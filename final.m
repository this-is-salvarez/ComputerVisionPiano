clear; close all;
vidReader = VideoReader('piano1.mp4');
result = VideoWriter('result.mp4');
open(result);
vidFrame = readFrame(vidReader);
[userinputx, userinputy] = transform(vidFrame);
itransformed = transformImage(userinputx, userinputy, vidFrame);
[bx, by] = findBlackKeys(itransformed); %coordinates of the black key centroids
[wx, wy] = findWhiteKeys(itransformed);%coordinates of the white key centroids

%Assigning notes

notes = ["C", "C#", "D", "D#", "E", "F","F#", "G", "G#", "A", "A#", "B"];
notesW =["C", "D", "E", "F", "G", "A", "B"];
notesB =["C#", "D#", "F#", "G#", "A#"];

newB = [notesB notesB notesB notesB "Null"; bx; by];
newW = [notesW notesW notesW notesW "C"; wx; wy];
counter=0;
while (hasFrame(vidReader ) && counter < 5)
    vidFrame = readFrame(vidReader);
    currentCroppedFrame = transformImage(userinputx, userinputy, vidFrame);
    temp = itransformed - currentCroppedFrame;
    %Iblur1 = imgaussfilt(temp,2);
    S2= strel('line',30, 0);
    %S3= strel('disk',5);
    temp1 = imopen(temp,S2);
    %temp1 = imerode(temp1,S3);
    figure, imshow(currentCroppedFrame);
    
    [L, num] = bwlabel(temp1);
    f = regionprops(L);
    max=0;
    count=0;
 
    for i=1:num
      area = f(i).Area;
    
        if area > max
            max=area;
            count=i;
        end
    end
    
    DISS  = [notesB notesB notesB notesB "Null" notesW notesW notesW notesW "C"]
    
    if max > 20
        counter=counter+1;
        rectangle('Position',f(count).BoundingBox,'EdgeColor','r', 'LineWidth',3);
        centFinger=f(count).Centroid;
        vec1=[f(count).Centroid(1) f(count).Centroid(2)]

        for i=1:length(newB(1,:))
            vec2=[str2num(newB(2,i)) str2num(newB(3, i))];
            DISS(2,i)=sqrt(sum((vec1-vec2).^2,2))
        end

        for i=1:length(newW(1,:))
            vec3=[str2num(newW(2,i)) str2num(newW(3, i))];
            DISS(2,i+length(newB(1,:)))=sqrt(sum((vec1-vec3).^2,2))
        end
        values= str2double(DISS(2,:));
        currentMin=min(values);
        for i=1:length(values)
            if(DISS(2,i)==num2str(currentMin))
                disp("The matching key is " + DISS(:, i));
            end
        end
    end
   
    writeVideo(result, currentCroppedFrame); 
end

close(result);

function keyboard = findKeyboard(vidFrame)
vidFrame = imrotate(vidFrame, 90);
Igrey = rgb2gray(vidFrame);
mask1 = vidFrame>210;
mask = uint8(mask1);
I2 = mask .* Igrey;
I2bw = uint8(imbinarize(I2)*255);
Igrey2 = rgb2gray(I2bw);
%Morphological analysis
S1= strel('disk',10, 0);
out1 = imdilate(Igrey2,S1); 
%figure, imshow(out1);
%Do connected component analysis
L = bwlabel(out1); 
%figure, imshow( label2rgb(L) );
% compute blob features & find large elements with area of the piano
f = regionprops(L);
areas = cat(2, f(:).Area);
indeces = find(areas>200000); % Find indices of nonzero elements
fXL = f(indeces);
m = numel(fXL); % get number of elements in fXL
%Crop image
keyboard= imcrop(vidFrame, [fXL.BoundingBox]);
end

function [ginputx, ginputy] = transform(vidFrame)
%Use thresholding to get the focus of the keyboard
keyboard = findKeyboard(vidFrame);
figure, imshow(keyboard);
%Do affine transformation
%[ginputx,ginputy] = ginput(4);
ginputx=[1.0e+03 *0.0360; 1.0e+03 *1.4540; 1.0e+03 *0.044; 1.0e+03 *1.4480];
ginputy=[104; 84; 308; 268];
end

function transformedImage = transformImage(xinput, yinput, vidFrame)
%Use thresholding to get the focus of the keyboard
keybpard = findKeyboard(vidFrame);
x1 = min(xinput);
x2 = max(xinput);
y1 = min(yinput);
y2 = max (yinput);

fpoints = [x1, y1; x2, y1; x1, y2; x2, y2];
transformation = fitgeotrans([xinput,yinput], fpoints, 'projective');
Ioutput = imwarp(keybpard, transformation);
%figure, imshow(Ioutput(:,:,:)); 
%Final crop the image
Ifinal = imcrop(Ioutput, [x1+xinput(1)/2, y1+yinput(1)/2, x2-x1, y2-y1]);
%figure, imshow(Ifinal(:,:,:));
Ifinalbw = rgb2gray(Ifinal);
mask1 = Ifinalbw>180;
mask = uint8(mask1);
I2 = mask .* Ifinalbw;
transformedImage = uint8(imbinarize(I2)*255);
%figure, imshow(transformedImage);
end

function [blackx, blacky] = findBlackKeys(itransformed)
S1= strel('disk',4);
I2 = imopen(itransformed,S1); 
%black Keys
blackKeys = imcomplement(I2);
S2= strel('line',5, 0);
blackKeys = imopen(blackKeys,S2);
%figure, imshow(blackKeys);
[Lblack,numblack] = bwlabel(blackKeys,8);% L is output matrix, num is number of blobs
fblack = regionprops(Lblack);
areasblackKeys = cat(2, fblack(:).Area);
indeces = find(areasblackKeys>5000); % Find indices of nonzero elements
fXLblack = fblack(indeces);
blackkeys = numel(fXLblack); % get number of elements in fXL
f=regionprops(Lblack);
for i = 1 : length(f)
    blackx(i) = f(i).Centroid(1); 
    blacky(i) = f(i).Centroid(2);
% line([blackx(i)-5 blackx(i)+5], [blacky(i) blacky(i)], 'Color', 'r');
% line([blackx(i) blackx(i)], [blacky(i)-5 blacky(i)+5], 'Color', 'r');
% rectangle('Position', fXLblack(i).BoundingBox, 'EdgeColor','r', 'linewidth', 3);
end
end

function [whitex, whitey] = findWhiteKeys(itransformed)
S1= strel('disk',4);
I2 = imopen(itransformed,S1); 
%white keys
S1= strel('line',3, 0);
whiteKeys = imerode(I2,S1);
%figure, imshow(whiteKeys);

[L_white, num_white] = bwlabel(whiteKeys);

f=regionprops(L_white);
for i = 1: length(f)
   whitex(i) = f(i).Centroid(1); 
   whitey(i) = f(i).Centroid(2);
%    line([whitex(i)-5 whitex(i)+5], [whitey(i) whitey(i)], 'Color', 'r');
%    line([whitex(i) whitex(i)], [whitey(i)-5 whitey(i)+5], 'Color', 'r');
%    rectangle('Position',f(i).BoundingBox,'EdgeColor','r', 'LineWidth', 1);
end

end