 function [det_pr, detTotProb] = detect_object_feature()

    sum_matched_ref = 0;
    sum_matched_rcvd = 0;
    
    origVid = 'parks.ts';
    srcVid = 'parks_ref.ts';
    rcvdVid = 'parks_copied.ts';
    origDir = 'orig_frames/';
    srcDir = 'src_frames/';     
    rcvdDir = 'rcvd_frames/';
    ffcmd = '/Users/sabur/Downloads/ffmpeg-2.8/ffmpeg';

    %fclose('all');
    nf = length(dir('orig_frames/*.png'));
    
    if(nf > 0)
        fprintf('\n Original video frames exist');
    else
        rmorig = sprintf('%s*.png', origDir);
        delete(rmorig);
        %origExe = sprintf('%s -i %s -r 24.04 %sorig_%%04d.png', ffcmd, origVid, origDir);
        origExe = sprintf('%s -i %s -r 25 %sorig_%%04d.png', ffcmd, origVid, origDir);
        system(origExe);
    end
    
    rmsrc = sprintf('%s*.png', srcDir);
    delete(rmsrc);
    rmrcvd = sprintf('%s*.png', rcvdDir);
    delete(rmrcvd);

    
    %srcExe = sprintf('%s -i %s -r 24.04 %ssrc_%%04d.png', ffcmd, srcVid, srcDir);
    srcExe = sprintf('%s -i %s -r 25 %ssrc_%%04d.png', ffcmd, srcVid, srcDir);
    system(srcExe);
    %rcvdExe = sprintf('%s -i %s -r 24.04 %srcvd_%%04d.png', ffcmd, rcvdVid, rcvdDir);
    rcvdExe = sprintf('%s -i %s -r 25 %srcvd_%%04d.png', ffcmd, rcvdVid, rcvdDir);
    system(rcvdExe);

    num_files = length(dir('orig_frames/*.png'));
    num_files1 = length(dir('src_frames/*.png'));
    num_files2 = length(dir('rcvd_frames/*.png'));
    
    display(num_files);
    display(num_files1);
    display(num_files2);
    
    if((num_files ~= num_files1) || (num_files1 ~= num_files2))
        display('Frame mismatch in reference and received video');
        return;
    else

        prob_vector = zeros(1, num_files);

        for i=1:num_files
            forig = sprintf('orig_frames/orig_%04d.png', i);
            orig = imread(forig);
            origImage = rgb2gray(orig);
            fsrc = sprintf('src_frames/src_%04d.png', i);
            ref = imread(fsrc);
            refImage = rgb2gray(ref);
            frcvd = sprintf('rcvd_frames/rcvd_%04d.png', i);
            rcv = imread(frcvd);
            rcvdImage = rgb2gray(rcv);

   
            %% Step 2: Detect Feature Points
            % Detect feature points in both images.
            origPoints = detectSURFFeatures(origImage);
            refPoints = detectSURFFeatures(refImage);
            rcvdPoints = detectSURFFeatures(rcvdImage);


            %% Step 3: Extract Feature Descriptors
            % Extract feature descriptors at the interest points in both images.
            [origFeatures, origPoints] = extractFeatures(origImage, origPoints);
            [refFeatures, refPoints] = extractFeatures(refImage, refPoints);
            [rcvdfFeatures, rcvdPoints] = extractFeatures(rcvdImage, rcvdPoints);


            %% Step 4: Find Putative Point Matches
            % Match the features using their descriptors. 
            %boxPairs = matchFeatures(refFeatures, rcvdfFeatures);
            [indrefpairs matchedref] =  matchFeatures(origFeatures, refFeatures);
            [indrcvdpairs matchedrcvd] =  matchFeatures(origFeatures, rcvdfFeatures);
            %[indices matches] = matchFeatures(refFeatures, rcvdfFeatures);
            %x=indref(:,1);
            %y=indrcvd(:,1);
            %d = intersect(x,y);
            
            src = size(matchedref);
            num_image_features = src(1);
            sum_matched_ref = sum_matched_ref + num_image_features;

            mtf = size(matchedrcvd);
            num_matched_features = mtf(1);
            sum_matched_rcvd = sum_matched_rcvd + num_matched_features;
            
            if(num_image_features  ~= 0)
                percentage_of_features_detected = num_matched_features/num_image_features;
            else
                percentage_of_features_detected = 1;
            end
            display(percentage_of_features_detected);
            prob_vector(i) = percentage_of_features_detected;

        end

        det_pr = mean(prob_vector);
        detTotProb = sum_matched_rcvd/sum_matched_ref;
        
        
    end
    
