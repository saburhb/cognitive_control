
function [psnrmean] = parks_psnr(badpackets_ref, badpackets)

    if (isequal(badpackets, badpackets_ref))
        fprintf('\n ###### D2D 0 same as reference case \t');
    else
        fprintf('\n ###### D2D 0 is NOT same as reference case \t');
    end

    for j=1:20
        badpackets(j) = 0;
    end
    for j=1:20
        badpackets_ref(j) = 0;
    end

    m = 0;
    n = size(badpackets);
    num_pkts = n(2);
    display(num_pkts);
    
    nr = size(badpackets_ref);
    num_pkts_ref = nr(2);
    display(num_pkts_ref);
    
    bytesPerPacket=188;
    count_bad = 0;
    
 
    delete('parks_copied.ts');
    delete('parks_ref.ts');
    copyfile('parks.ts', 'parks_copied.ts', 'f');
    copyfile('parks.ts', 'parks_ref.ts', 'f');
    
    
    %%%%%%%%%%%%%%%%%%%%% COURRUPT REFERENCE VIDEO FILE %%%%%%%%%%%%%%%%%%%%%%
    fp = fopen('parks_ref.ts', 'rb+', 'n');

    for index=1:num_pkts_ref
        if (badpackets_ref(1,index) == 1) %It indicates a bad packet
            count_bad = count_bad + 1;

            offset = (index-1)*bytesPerPacket;

            % READ before overwrite
            fseek(fp, offset, 'bof');
            data=fread(fp, bytesPerPacket);
            %fprintf('%c',dec2hex(data)); 
            %fprintf('\n');

            overwrite=repmat(0,bytesPerPacket,1);
            fseek(fp, offset, 'bof');
            fwrite(fp, data .* overwrite);

        end
    end
    fclose(fp);
 
    fprintf('\n REFERENCE FILE CORRUPTION DONE \n');
    
    
     %%%%%%%%%%%%%%%%%%%%% COURRUPT RECVD VIDEO FILE %%%%%%%%%%%%%%%%%%%%%%
    fp = fopen('parks_copied.ts', 'rb+', 'n');

    for index=1:num_pkts
        if (badpackets(1,index) == 1) %It indicates a bad packet
            count_bad = count_bad + 1;

            offset = (index-1)*bytesPerPacket;

            % READ before overwrite
            fseek(fp, offset, 'bof');
            data=fread(fp, bytesPerPacket);
            %fprintf('%c',dec2hex(data)); 
            %fprintf('\n');

            overwrite=repmat(0,bytesPerPacket,1);
            fseek(fp, offset, 'bof');
            fwrite(fp, data .* overwrite);

        end
    end
    fclose(fp);
 
    fprintf('\n RCVD FILE CORRUPTION DONE \n');
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%  Now Convert TS files to raw Video formats  %%%%%%%%%%%%%%%%%%

    srcVid='parks.ts';
    rcvVid='parks_copied.ts';

    srcRawVid='raw_parks_src.yuv';
    rcvRawVid='raw_parks_rcv.yuv';

    ffcmd = '/Users/sabur/Downloads/ffmpeg-2.8/ffmpeg';

    % stream_src.ts needs to converted once to yuv
    if exist(srcRawVid, 'file')
        % Do not convert src TS again again
        fprintf('\n src YUV exists\n');
    else
        srcExe = sprintf('%s -i %s -vcodec rawvideo -pix_fmt yuv420p %s', ffcmd, srcVid, srcRawVid);
        system(srcExe);
    end

    %%%%% Convert the new rcvd video %%%%%%
    delete(rcvRawVid);
    rcvExe = sprintf('%s -i %s -vcodec rawvideo -pix_fmt yuv420p %s', ffcmd, rcvVid, rcvRawVid);
    [status, result] = system(rcvExe);

    fc = fopen('decodelog.txt', 'w+');
    fprintf(fc, '%s', result);
    fclose(fc);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%% Now Calculate MSE of rcvd raw video w.r.t. src raw video %%%%%%%

    fwidth = 0.5;
    fheight= 0.5;

    %width = 480;
    %height= 360;
    
    width = 640;
    height= 360;

    %get Filesize and Framenumber
    filep = dir(srcRawVid); 
    fileBytes = filep.bytes; %Filesize1
    clear filep
    framenumber1 = fileBytes/(width*height*(1+2*fheight*fwidth)); %Framenumber1

    filep2 = dir(rcvRawVid); 
    fileBytes2 = filep2.bytes; %Filesize2
    clear filep2
    framenumber2 = fileBytes2/(width*height*(1+2*fheight*fwidth)); %Framenumber2

    display(framenumber1);
    display(framenumber2);

    if mod(framenumber1,1) ~= 0 | mod(framenumber2,1) ~= 0 | framenumber1~=framenumber2
            display('Error: wrong resolution, format, filesize or different video lengths');
            %psnrmean = 0;
    else

        mserr = zeros(1, framenumber1);
        frame_array=zeros(1,framenumber1);

        for cntf = 1:framenumber1
            YUV1 = loadFileYUVData(width,height,cntf,srcRawVid,fheight,fwidth);
            YUV2 = loadFileYUVData(width,height,cntf,rcvRawVid,fheight,fwidth);
            mserr(cntf) = sum((double(YUV1(:))-double(YUV2(:))).^2)/length(YUV1(:));
            frame_array(cntf) = cntf;
        end

            m = mean(mserr);
            display(m);

            avgPSNR= 10 * log10((255^2)/m);
            display(avgPSNR);
            psnrmean = avgPSNR;
    end


end
