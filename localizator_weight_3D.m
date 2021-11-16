function [output, imgdim] = localizator_weight_3D(program_status, gain, camera_offset, rwin, filtersigma,  framestart, frameend, xmin, xmax, ymin, ymax,...
    imgfile, show_prev, pvalue, sigpx,noperframe,fittype,deflation)
% clear

%%
load('darkimage_cam1_wo_corr.mat');
load('gain_det1.mat');
% ROI=[897,897,1152,1152];   center of camera
L=1013;
T=789;
R=1352;
B=1252;
dark_img=imgseq{1,1}(T:B,L:R);
var_img=imgseq{2,1}(T:B,L:R);
gain_img=sCMOSgain(T:B,L:R);
% var_img=var_img(1:256,1:256);
% gain_img=gain_img(1:256,1:256);
% dark_img=dark_img(1:256,1:256);
gain_img(gain_img<1)=1;
Rs = 8;                         % size of segmented images
alpha = 0.1;                    % weight factor of noise contribution only
Pixelsize = 0.096;              % the pixel size on sample plane, unit is micron
Lambda = 0.581;                 % emission wavelength of the sample, unit is micron
NA = 1.1;                       % numerical aperture of the objective lens
iterationN = 5;                % number of iterations
sz=256;                         % image size
%%



t1=tic;
%show process status
set(program_status, 'String', 'getting image file information'); %show progress
drawnow update; %show progress
% t1=tic;
%get file information
imgfilesets = size(imgfile, 1); %get number of image files to process
for s = 1 : imgfilesets %go through all image files
    set(program_status, 'String', ['getting image file information: ' imgfile{s}]); %show progress
    drawnow update; %show progress
    IMGinfo{s} = imfinfo(imgfile{s}); %get information of image set s
    zdim(s) = length(IMGinfo{s}); %add number of frames in set s to zdim
    %     IMGinfo = imfinfo(imgfile{s}); %get information of image set s
    %     zdim(s) = length(IMGinfo); %add number of frames in set s to zdim
end

IMGinfo_temp = IMGinfo{1};
xdim = IMGinfo_temp(1).Height; %get x dimension
ydim = IMGinfo_temp(1).Width; %get y dimension
% 
% xdim = IMGinfo(1).Height; %get x dimension
% ydim = IMGinfo(1).Width; %get y dimension
imgdim = [xdim ydim sum(zdim)]; %store image dimensions in imgdim
% clear IMGinfo;
%set parameters
if frameend == -1 %all frames will be processed
    frameend = sum(zdim); %get number of frames of all image files
end
if xmax == -1%full frame will be processed
    xmax = xdim; %end of roi in x direction in pixel
    ymax = ydim; %end of roi in y direction in pixel
end
xmax = min(xmax,xdim);
xmin = max(xmin,1);
ymax = min(ymax,ydim);
ymin = max(ymin,1);

% rbox = 2;
%
% %compute square of theoretical PSF width
% sig2px = sigpx * sigpx;
%
% %compute intensity distribution given by theoretical gaussian psf
% numpx = 10;
% offsetpx = numpx + 1;
% for x0px = -numpx : numpx
%     for y0px = -numpx : numpx
%         PSF(x0px + offsetpx, y0px + offsetpx) = exp(-(x0px.^2 + y0px.^2) / 2 / sig2px);
%     end
% end
%
% %compute full intensity given by gaussian psf
% ginf = sum(sum(PSF));
%
% %radius of inner box
% ibox = rbox - 1;
%
% %get width of inner box
% innerbox = 2 * ibox + 1;
%
% %get number of pixel in inner box
% innerbox2 = innerbox * innerbox;
%
% %get width of full box
% fullbox = 2 * rbox + 1;
%
% %get number of pixel in full box
% fullbox2 = fullbox * fullbox;
%
% %get number of pixel in verge of full box
% outerbox2 = fullbox2 - innerbox2;
%
% %circular sum over inner box
% innerfilter = fspecial('disk', innerbox / 2) * innerbox2;
%
% %circular sum over full box
% flboxfilter = fspecial('disk', fullbox / 2) * fullbox2;
%
% %filter image
% INNERPSF = imfilter(PSF, innerfilter, 'replicate');
%
% %filter image
% FLBOXPSF = imfilter(PSF, flboxfilter, 'replicate');
%
% %get image averaged over verge of full box
% OUTERPSF = FLBOXPSF - INNERPSF;
%
% %compute differential image
% DIFFPSF = INNERPSF / innerbox2 - OUTERPSF / outerbox2;
%
% %get value at center of psf
% gdbox = max(max(DIFFPSF));

% gdbox = 0.364985147079997;
% ginf = 15.674212206566814;
% %compute theoretical difference intensity for clean procedure
% diffthresh = threshold / ginf * gdbox;
% %assignin('base','diffthresh',diffthresh);
%
% %compute theoretical difference intensity for clean procedure
% diffupthresh = upthreshold / ginf * gdbox;

%use dark state offset file as background
% dso = mean2(double(imread(dsofile))); %compute mean dark state offset

% %initialize image sum
% IMGSUM(1 : xdim, 1 : ydim) = 0; %initialization


%open new figure window for preview
if show_prev == 1 %preview enabled
    h = figure('Name', 'Current Image Frame'); %open a new figure window
end

% disp(['The time to initialization:' num2str(toc(t1)) 's']);


%hpxtemp = zeros(20,3);
output = zeros((frameend-framestart+1)*noperframe,19);
frames = 0; %initialization
noframe = frameend-framestart+1;
moleculescum = 0;
% rwin = 6; %fit window area is (2*rbox+1)^2 pixel
fitwinarea = (2*rwin+1)^2;
mem = feature('memstats');
% maxfitwin = floor(0.5*mem/8/fitwinarea)
maxfitwin = 20000;
disp(['Time for Initialization: ' num2str(toc(t1))]); %for debugging purpose
% t3 = tic;
tseek = 0;
% tread = 0;
for i = 1:ceil((frameend-framestart+1)*noperframe/maxfitwin)
    esmol = min(maxfitwin,(frameend-framestart+1)*noperframe);
    fitwindow = zeros(2*rwin+1,2*rwin+1,esmol);
    % fitwindow = zeros(2*rwin+1,2*rwin+1,maxfitwin);
    hpx = zeros(min(maxfitwin,(frameend-framestart+1)*noperframe),3);
    %fitwindowtemp = zeros(7,7);
    IMG = zeros(xdim,ydim);
    molecules = 0; %initialization
    for f = framestart : min(framestart+floor(maxfitwin/noperframe)-1,frameend)
        %    tic;
        set(program_status, 'String', ['analyzing frame number ' num2str(f) ' of ' num2str(frameend)]); %show progress
        drawnow update; %show progress
        %     disp(['t211:' num2str(toc)])
        %disp(['frame:' num2str(f)]);
        
        %     tic;
        s = 0; %initialization
        sf = 0; %initialization
        while sf < f %set frame count smaller than frame count
            s = s + 1; %increase set count
            sf = sf + zdim(s); %add frames of set s to frame count
        end
        sf = sf - zdim(s); %go back one frame set
        %     disp(['t212:' num2str(toc)])
        
        %         t2 = tic;
         IMG = (double(imread(imgfile{s}, f - sf,'info',IMGinfo{s})) - camera_offset)/gain; %store dso corrected image in matrix IMG
%          IMG=IMG(1:256,1:256);
         IMG=(IMG-dark_img)./gain_img;
         IMG(IMG<0)=1e-6;
%          [IMG] = reducenoise(Rs,IMG,var_img,gain_img,sz,Pixelsize,NA,Lambda,alpha,iterationN);

         
         IMG=medfilt2(IMG);
         IMG = max(IMG, 0); 
         %store dso corrected image in matrix IMG
      %  IMG = double(imread(imgfile{s}, f - sf))/gain; %store dso corrected image in matrix IMG
        %         tread = tread+toc(t2);
        %     disp(['t2131:' num2str(toc)])
        
        %         IMG = (IMG )/gain;
        %         assignin('base','IMG',IMG);
        %IMG = IMG / gain; %convert pixel intensity count to photon count
        %      tic;
        %IMG(IMG < 0) = 0; %set negative intensity counts to zero
        t2 = tic;
        %         seeker_robust_Acc_sub_exoverlap
        %          [hpxtemp fitwindowtemp m] = seeker_robust_Acc_test(IMG, rwin, filtersigma, xmin, xmax, ymin, ymax, pvalue,noperframe);
        if deflation ==1
            [hpxtemp fitwindowtemp m] = seeker_robust_Acc_test(IMG, rwin, filtersigma, xmin, xmax, ymin, ymax, pvalue,noperframe);
        elseif deflation ==2
            [hpxtemp fitwindowtemp m] = seeker_robust_Acc_sub_exoverlap_2(IMG, rwin, filtersigma, xmin, xmax, ymin, ymax, pvalue,noperframe);
        end
        tseek = tseek+toc(t2);
        %     assignin('base','hpx',hpxtemp);
        %     timer1 = timer1 + timer1temp;
        %     timer2 = timer2 + timer2temp;
        %     tic;
        
        if m ~= 0
            hpxtemp(:,3) = f;
            
            frames = frames+1;
            hpx(molecules+1:molecules+m,:) = hpxtemp;
            fitwindow(:,:,molecules+1:molecules+m) = fitwindowtemp;
            molecules = molecules+m;
        end;
        
        %     disp(['t23:' num2str(toc)]);
     %%   
        if show_prev == 1 %show preview if selected
            %             R = ceil(3*sigpx);   % cutoff radius of the gaussian kernel
            %             M = zeros(2*R+1); % KJ
            %             for i = -R:R,
            %                 for j = -R:R,
            %                     M(i+R+1,j+R+1) = exp(-(i*i+j*j)/2/sigpx/sigpx);
            %                 end
            %             end
            %             M = M/sum(M(:));
            %              imagesc(filter2(M,IMG(xmin : xmax, ymin : ymax)))
            %                 colormap(gray)
            %                 axis equal
            %                 axis off
            
            
            figure(h);
            imagesc(IMG(xmin : xmax, ymin : ymax)); %display current frame
            %             dipshow(IMG(xmin : xmax, ymin : ymax));
            colormap(gray); %set colormap
            axis equal; %set axis to same size
            axis off; %do not display axis
            hold on; %hold image frame
            if hpxtemp ~= 0
                for m = 1 : size(hpxtemp,1) %draw box around highpx found by the_seeker
                    plot(hpxtemp(m,2)- rwin - ymin + 1, hpxtemp(m,1)- rwin - xmin + 1 : hpxtemp(m,1)+ rwin - xmin + 1,...
                        hpxtemp(m,2) + rwin - ymin + 1, hpxtemp(m,1) - rwin - xmin + 1 : hpxtemp(m,1) + rwin - xmin + 1,...
                        hpxtemp(m,2) - rwin - ymin + 1 : hpxtemp(m,2) + rwin - ymin + 1, hpxtemp(m,1) - rwin - xmin + 1,...
                        hpxtemp(m,2) - rwin - ymin + 1 : hpxtemp(m,2) + rwin - ymin + 1, hpxtemp(m,1) + rwin - xmin + 1,...
                        'Marker', '.', 'Color', 'r');
                end
            end
            hold off; %clear for next image frame
            pause(0.3); %wait 0.5 seconds
        end
        
    end
    disp(['seeker:' num2str(tseek)]);
    %      disp(['readIm:' num2str(tread)]);
    
    hpx = hpx(1:molecules,:);
%         fitwindow = fitwindow(:,:,1:molecules);
%     assignin('base','fitwindow',fitwindow);
    
    set(program_status, 'String', ['fitting frames ' num2str(framestart) ' : ' num2str(min(framestart+floor(maxfitwin/noperframe)-1,frameend))]); %show progress
    drawnow update; %show progress
    t2 = tic;
    if molecules ~=0
        %         [X Y N BG S CRLBx CRLBy CRLBn CRLBb CRLBs LogL] = mGPUgaussMLE(dip_image(fitwindow(:,:,1:molecules)),sigpx,20,fittype);
        
        switch fittype
            case 1
                % [X Y N BG CRLBx CRLBy CRLBn CRLBb CRLBs LogL] = [P CRLB LL];
                [P CRLB LL]=GPUgaussMLEv2(permute(single(fitwindow(:,:,1:molecules)),[2 1 3]),sigpx,40,fittype);
                output(moleculescum+1:moleculescum+molecules,1:4) = [P(:,2)+hpx(:,1)-rwin P(:,1)+hpx(:,2)-rwin P(:,3) P(:,4)];
                output(moleculescum+1:moleculescum+molecules,7:10) = sqrt(CRLB);
                output(moleculescum+1:moleculescum+molecules,13) = LL;
                output(moleculescum+1:moleculescum+molecules,14) = hpx(:,3);
                output(moleculescum+1:moleculescum+molecules,17:18) = [P(:,2) P(:,1)];
            case 2
                [P CRLB LL]=GPUgaussMLEv2(permute(single(fitwindow(:,:,1:molecules)),[2 1 3]),1,40,fittype);
                output(moleculescum+1:moleculescum+molecules,1:5) = [P(:,2)+hpx(:,1)-rwin P(:,1)+hpx(:,2)-rwin P(:,3) P(:,4) P(:,5)];
                output(moleculescum+1:moleculescum+molecules,7:11) = sqrt(CRLB);
                output(moleculescum+1:moleculescum+molecules,13) = LL;
                output(moleculescum+1:moleculescum+molecules,14) = hpx(:,3);
                output(moleculescum+1:moleculescum+molecules,17:18) = [P(:,2) P(:,1)];
            case 4
                [P CRLB LL]=GPUgaussMLEv2(permute(single(fitwindow(:,:,1:molecules)),[2 1 3]),1,40,fittype);
                output(moleculescum+1:moleculescum+molecules,1:6) = [P(:,2)+hpx(:,1)-rwin P(:,1)+hpx(:,2)-rwin P(:,3) P(:,4) P(:,5) P(:,6)];
                output(moleculescum+1:moleculescum+molecules,7:12) = sqrt(CRLB);
                output(moleculescum+1:moleculescum+molecules,13) = LL;
                output(moleculescum+1:moleculescum+molecules,14) = hpx(:,3);
                output(moleculescum+1:moleculescum+molecules,17:18) = [P(:,2) P(:,1)];
        end
        %         output = [output;outputtemp];
    end;
    framestart = framestart+floor(maxfitwin/noperframe);
    moleculescum= molecules+moleculescum;
    clear fitwindow
    disp(['fitIm:' num2str(toc(t2))]);
end
set(program_status, 'String', 'Weight Molecules'); %show progress
drawnow update; %show progress
output = output(1:moleculescum,:);
%exculde molecules with X,Y<2 or X,Y>5 or likelihood<-1 or error fit CRLBx=0
% assignin('base','outputnoweith',output);
output = output(output(:,17)>1.5&output(:,17)<2*rwin+1-1.5&output(:,18)>1.5&output(:,18)<2*rwin+1-1.5&output(:,7)~=0,:);
% disp(['Avarage process time for each frame is no initial:' num2str(toc(t3)) 's/tframe']);
disp(['Avarage time for each frame without weighting is:' num2str(toc(t1)/noframe) 's/tframe']);
try
    output = weightlocation(output);
catch
    disp('All positions are at the same time... go back!');
end
% [Y,I]=sort(output(:,12));
% output = output(I,:);
output(1,19) = xdim;
output(2,19) = ydim;
disp(['number of frames with molecules found: ' num2str(frames)]); %for debugging purpose
disp(['total number of molecules found: ' num2str(moleculescum)]); %for debugging purpose
disp(['Avarage time for each frame is:' num2str(toc(t1)/noframe) 's/tframe']);
% assignin('base','output',output);














