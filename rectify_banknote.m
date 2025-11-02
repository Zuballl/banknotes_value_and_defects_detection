function [rectImg, rectMask, success] = rectify_banknote(imgColor, templateData)
%RECTIFY_BANKNOTE Try to remove perspective distortion of a photographed banknote.
%   [rectImg, rectMask, success] = rectify_banknote(imgColor, templateData)
%   Attempts to find the largest connected component corresponding to the
%   banknote, extract its contour, approximate four corner points and
%   compute a projective transform to warp the banknote to the canonical
%   template size (choosing the closest template by aspect ratio).

    success = false;
    rectImg = imgColor;
    rectMask = true(size(rgb2gray(imgColor)));

    try
        gray = rgb2gray(imgColor);
        % Smooth and detect edges
        I = imgaussfilt(gray, 2);
        bw = edge(I, 'Canny', [0.04 0.2]);
        bw = imclose(bw, strel('disk', 3));
        bw = imfill(bw, 'holes');
        bw = bwareaopen(bw, 500);

        cc = bwconncomp(bw);
        if cc.NumObjects == 0
            return;
        end
        stats = regionprops(cc, 'Area', 'PixelIdxList', 'BoundingBox');
        [~, idxMax] = max([stats.Area]);

        mask = false(size(bw));
        mask(stats(idxMax).PixelIdxList) = true;

        boundaries = bwboundaries(mask, 'noholes');
        srcPts = [];
        if ~isempty(boundaries)
            boundary = boundaries{1}; % rows, cols
            try
                k = convhull(boundary(:,2), boundary(:,1));
                hullPts = [boundary(k,2), boundary(k,1)];
            catch
                hullPts = [boundary(:,2), boundary(:,1)];
            end

            if size(hullPts,1) >= 4
                % cluster hull points to 4 corner groups (works well for noisy contours)
                try
                    idxk = kmeans(hullPts, 4, 'Replicates', 5);
                    centers = zeros(4,2);
                    for ii = 1:4
                        centers(ii,:) = mean(hullPts(idxk==ii,:),1);
                    end
                    % order points consistently (clockwise)
                    c = mean(centers,1);
                    ang = atan2(centers(:,2)-c(2), centers(:,1)-c(1));
                    [~, order] = sort(ang);
                    srcPts = centers(order,:);
                catch
                    srcPts = [];
                end
            end
        end

        % fallback -> bounding box
        if isempty(srcPts)
            bb = stats(idxMax).BoundingBox; % [x y w h]
            srcPts = [bb(1), bb(2); bb(1)+bb(3), bb(2); bb(1)+bb(3), bb(2)+bb(4); bb(1), bb(2)+bb(4)];
        end

        % choose template size by matching aspect ratio
        nTemp = length(templateData);
        targetAspects = ones(nTemp,1);
        for k = 1:nTemp
            try
                tmask = templateData(k).Mask;
                sbb = regionprops(tmask, 'BoundingBox');
                if isempty(sbb)
                    bbox = [1 1 size(tmask,2) size(tmask,1)];
                else
                    bbox = sbb(1).BoundingBox;
                end
                targetAspects(k) = bbox(3) / bbox(4);
            catch
                targetAspects(k) = 1;
            end
        end
        srcAspect = (max(srcPts(:,1))-min(srcPts(:,1))) / (max(srcPts(:,2))-min(srcPts(:,2)) + eps);
        [~, chosenIdx] = min(abs(targetAspects - srcAspect));

        tmplImg = templateData(chosenIdx).Image;
        hT = size(tmplImg,1); wT = size(tmplImg,2);
        dstPts = [1,1; wT,1; wT,hT; 1,hT];

        % Fit projective transform and warp
        tform = fitgeotrans(srcPts, dstPts, 'projective');
        outputView = imref2d([hT, wT]);
        rectImg = imwarp(imgColor, tform, 'OutputView', outputView);
        rectMask = imwarp(mask, tform, 'OutputView', outputView, 'Interp', 'nearest');
        success = true;
    catch ME
        % If anything goes wrong, return original image and success=false
        warning('rectify_banknote:failed', 'Rectification failed: %s', ME.message);
        success = false;
        rectImg = imgColor;
        rectMask = true(size(rgb2gray(imgColor)));
    end
end
