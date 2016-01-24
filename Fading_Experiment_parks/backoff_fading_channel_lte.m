function backoff_fading_channel_lte()
  
    tic;
    
    num_quantize = 100;
    quantum_prob = (1/num_quantize);
    num_grids=(num_quantize * num_quantize);

    % Create the matrix for channel fading
    Cf_u = zeros(1,num_quantize);
    Cf_i = zeros(1,num_quantize);

    cum_steady_state_prob = zeros(1,num_quantize);
    T_mat_u = zeros(num_quantize,num_quantize);
    T_mat_i = zeros(num_quantize,num_quantize);
            
    ps = zeros(1,num_quantize);
    ps(num_quantize) = 5;
    Cf_u(num_quantize) = 2.2; %ps(num_quantize); %%% Channel fading for UE %%%%
    Cf_i(num_quantize) = 2.2; %ps(num_quantize);  %%% Channel fading for interfering node %%%%

    fvid = fopen('vid_out.txt','w+');

    %%%%%%%% Create Rayleigh Distribution for fading and its quantization %%%%

    k=0;
    pi = 3.14;
    b = sqrt(2/pi);

    for i=1:(num_quantize-1)
        k=k+quantum_prob;
        y = raylinv(k,b);
        ps(i) = y;
        Cf_u(i) = y;
        Cf_i(i) = y;
    end

    display(ps);


    %%%%%%% Get Cumulative transition probabilities %%%%%%%%
    [T_mat_u, T_mat_i] = trans_prob_matrix();


    %%%%%%%%%%%% physical layer properties %%%%%%%%%%%%%%%
    K = 1000; 
    d1 = 300; 
    P_u = 500; 
    N = 0.02;                
    d2=2.5; 
    d=200; 
    P_d = 4; %14; 
    threshold = 15; %% in dB
 
    fprintf(fvid, '\n ******************************************************** \n');
    
    fprintf(fvid, '\nUE Power : %.3f \n', P_u);
    fprintf(fvid, 'D2D Power : %.3f \n', P_d);
    
    %%%%%%%%%%%%%%%% Get the decode Maps %%%%%%%%%%%%%%%%

    decodable_mat_ue = zeros(num_quantize,num_quantize);
    decodable_mat_ue_ref = zeros(num_quantize,num_quantize);    
    decodable_mat_d2d = zeros(num_quantize,num_quantize);    

    num_decodable_ue = 0;
    num_decodable_d2d = 0;
    num_decodable_ue_ref = 0;

    %%%%%%%%%%%%% Decode Map for Reference video with no D2D %%%%%%%%%%%%% 
    for i=1:num_quantize
        for j=1:num_quantize
            s = (P_u * Cf_u(i) * K/(d1 * d1))/N ;  %% modified
            sinr = 10 * log10(s);
            if sinr >= threshold
                decodable_mat_ue_ref(i,j) = 1;
                num_decodable_ue_ref = num_decodable_ue_ref + 1;
            end
        end
    end

    fprintf('#### \n Percentage of packet loss for reference video = %.3f ####\n', (num_grids-num_decodable_ue_ref)*100/num_grids);
    fprintf(fvid, '#### Percentage of packet loss for reference video = %.3f ####\n', (num_grids-num_decodable_ue_ref)*100/num_grids);


    %%%%%%%%%%%%%% UE Decode Map for video frames %%%%%%%%%%%%%%%%%
    for i=1:num_quantize
        for j=1:num_quantize
            s = (P_u * Cf_u(i) * K/(d1 * d1))/ (N + (P_d * Cf_i(j) * K/(d1 * d1)) );
            sinr = 10 * log10(s);
            if sinr >= threshold
                decodable_mat_ue(i,j) = 1;
                num_decodable_ue  = num_decodable_ue + 1;
            end
        end
    end

    fprintf('#### \n Percentage of packet loss for video packets = %.3f ####\n', (num_grids-num_decodable_ue)*100/num_grids);
    fprintf(fvid, '#### Percentage of packet loss for video packets = %.3f ####\n', (num_grids-num_decodable_ue)*100/num_grids);

    
    %%%%%%%%%%%%%% D2D Decode Map for d2d frames %%%%%%%%%%%%%%%%%
    for i=1:num_quantize
        for j=1:num_quantize
            s = (P_d * Cf_i(i) * K /(d2 * d2))/ (N + (P_u * Cf_u(j) * K /(d * d)));
            sinr = 10 * log10(s);
            if sinr >= threshold
                decodable_mat_d2d(i,j) = 1;
                num_decodable_d2d  = num_decodable_d2d + 1;
            end
        end
    end

    fprintf('#### \n Percentage of d2d packet loss = %.3f ####\n', (num_grids-num_decodable_d2d)*100/num_grids);
    fprintf(fvid, '#### Percentage of d2d packet loss = %.3f ####\n', (num_grids-num_decodable_d2d)*100/num_grids);

    fprintf('\n*****************************************************\n');
    
    TOTAL_TS_PACKETS = 8154;
    
    %%%%%%%%%%%%%%%% Now apply protocol and calulate metrics on recvd video
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ptrans = [0 0.1 0.2 0.3 0.4 0.6 0.8 1];       
    %policy=1  ;  %% Policy : 0 - base ; Policy: 1 - our model 
    %MAX_ITERATIONS = 10;
    numset = 10; %10;
    repn = 8; 
    it = 1;
        
    
    cumAvgPsnr =zeros(1,repn);
    cumAvgThrptUe =zeros(1,repn);
    cumAvgThrptD2d =zeros(1,repn);
    cumAvgDetProb =zeros(1,repn);
    cumTotDetProb =zeros(1,repn);
    
    cumAvgPsnr1 =zeros(1,repn);
    cumAvgThrptUe1 =zeros(1,repn);
    cumAvgThrptD2d1 =zeros(1,repn);
    cumAvgDetProb1 =zeros(1,repn);
    cumTotDetProb1 =zeros(1,repn);
    
    cumAvgPsnr2 =zeros(1,repn);
    cumAvgThrptUe2 =zeros(1,repn);
    cumAvgThrptD2d2 =zeros(1,repn);
    cumAvgDetProb2 =zeros(1,repn);
    cumTotDetProb2 =zeros(1,repn);
    
    cumAvgPsnr3 =zeros(1,repn);
    cumAvgThrptUe3 =zeros(1,repn);
    cumAvgThrptD2d3 =zeros(1,repn);
    cumAvgDetProb3 =zeros(1,repn);
    cumTotDetProb3 =zeros(1,repn);
    
    cumAvgPsnr4 =zeros(1,repn);
    cumAvgThrptUe4 =zeros(1,repn);
    cumAvgThrptD2d4 =zeros(1,repn);
    cumAvgDetProb4 =zeros(1,repn);
    cumTotDetProb4 =zeros(1,repn);
    
    fprintf(fvid, '\n\n ################# START %d set of simulations ################# \n', it);
    
    for it = 1:numset
        
        fprintf(fvid, ' ********* Simulation iteration %d *********: \n', it);
        fprintf(fvid, ' *******************************************: \n');
        
        
        fprintf('\n\n The state transition is as below \n');

        row_cum_trans_prob_u = zeros(1,num_quantize);
        row_cum_trans_prob_i = zeros(1,num_quantize);
        channel_state = zeros(2, TOTAL_TS_PACKETS);

        [ue_state, in_state] =  getSteadyState(); 
        row_cum_trans_prob_u = T_mat_u(ue_state, :);
        row_cum_trans_prob_i = T_mat_i(in_state, :);

        %%%%%%%%%%%% Run the Markov chain over its %%%%%%%%%
        num_lines = 0;
        fp = fopen('bunny_tsframes.txt', 'r+');
        tline = fgetl(fp);

        while (ischar(tline)) 
            num_lines = num_lines + 1;  

            c = strsplit(tline, ',');
            k = str2num(char(c(1))); 
            frame_type = str2num(char(c(2)));
            pkt_idx = k;

            if(frame_type ~= 0)
                [next_state_u, next_state_i] = getNextState();
                channel_state(1, k) = next_state_u;
                channel_state(2, k) = next_state_i;
                row_cum_trans_prob_u = T_mat_u(next_state_u, :);
                row_cum_trans_prob_i = T_mat_u(next_state_i, :);
                
                %%%%%%%%%% save the same channel propagation to use for different
                %%%%%%%%%% Tx probabilities in the protocol %%%%%%%%%%%%%
               
            end
            
            tline = fgetl(fp);
        end
        fclose(fp);
        
        var = 0;
        for t=1:3
            
            if(var == 0)
                policy = 0;
                Delta_Max = 0;
            elseif(var == 1)
                policy = 1;
                Delta_Max = 0;
            elseif(var == 2)
                policy = 2;
                Delta_Max = 5;
            elseif(var == 3)
                policy = 2;
                Delta_Max = 7;
            else
                policy = 2;
                Delta_Max = 10;
            end
            
            fprintf(fvid, '\n ******* START iteration %d for policy %d ******** \n', it, policy);
             
            ind = 1;
            retry=0;
            while (ind < (repn+1))   %%%%% for loop with different packet loss percentage %%%%

                
                try
                    tx_prob = ptrans(ind); %% transmission probability of D2D  
                    %Delta_Max =  3; %ceil(1/tx_prob);

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%%%%%%%%%%%% Now run a Markov chain %%%%%%%%%%%%%
                    badpackets = zeros(1,TOTAL_TS_PACKETS);
                    badpackets_ref = zeros(1,TOTAL_TS_PACKETS);
                    pkt_idx = 1;
                    num_bad_pkts = 0;
                    num_bad_pkts_ref = 0;
                    num_bad_pkts_d2d = 0;
                    num_total_pkts_d2d = 0;
                    num_total_pkts = 0;
                    num_lines = 0;
                    lte_prev = 0;
                    lte_curr = 0;
                    curr_frame = 0;
                    prev_frame = 0;
                    failure_flag = 0;
                    delta_slot = Delta_Max ;


                    fprintf(fvid, '\n ******* Simulation for specific Tx probability = %0.2f ******** \n', tx_prob); 

                    %%%%%%%%%%%% Run the Markov chain over its %%%%%%%%%
                    fp = fopen('bunny_tsframes.txt', 'r+');
                    tline = fgetl(fp);

                    while (ischar(tline)) 
                        num_lines = num_lines + 1;  

                        c = strsplit(tline, ',');
                        k = str2num(char(c(1))); 
                        frame_type = str2num(char(c(2)));
                        pkt_idx = k;

                        if(frame_type ~= 0)
                            curr_frame = frame_type;
                            num_total_pkts = num_total_pkts + 1;

                            fprintf('..... Processing packet : %d , of frame type: %d \n ', num_total_pkts, frame_type);

                            next_state_u = channel_state(1, k);
                            next_state_i = channel_state(2, k);

                            %%%%%%%% Transmission policy for UE and D2D %%%%%%%%%%
                            if(policy == 0)
                                policy0(ind, next_state_u, next_state_i);
                                
                            elseif(policy == 1) %%% flat transmission control
                                if(curr_frame == 1) %% I frame
                                    d2d_idle(next_state_u, next_state_i);
                                else
                                    policy0(ind, next_state_u, next_state_i);
                                end

                            else
                                if(failure_flag == 1)
                                    if(curr_frame == prev_frame) %% same frame continues
                                        if(curr_frame == 1) %% I frame
                                            d2d_idle(next_state_u, next_state_i);
                                        else                %% P frame
                                            if(delta_slot > 0)
                                                delta_slot = delta_slot-1;
                                                d2d_idle(next_state_u, next_state_i);
                                            else
                                                failure_flag = 0;
                                                delta_slot = Delta_Max;
                                                policy0(ind, next_state_u, next_state_i);
                                            end
                                        end

                                    else   %%%%% If the current packet is from different frame
                                        failure_flag = 0;
                                        policy0(ind, next_state_u, next_state_i);
                                    end
                                else
                                    policy0(ind, next_state_u, next_state_i);
                                end
                            end %% end of backoff policy  

                        end  %% end of functionality for video frames %%

                        lte_prev = lte_curr;
                        if(lte_prev == 0)
                            fprintf('>>>>>>>>>>>>>>>> LTE packet# %d failed >>>>>>>>>>>>>\n', num_total_pkts);
                            failure_flag = 1; %% failure flag is only set to 1 on failure
                        end
                        prev_frame = curr_frame;

                        tline = fgetl(fp);

                    end %% end of file activity for all TS packets

                    fclose(fp);

                    %%%%%%%%% for all tx_prob call
                    %%%%%%%%% calculateMetrics %%%%%%%%%%%%%
                    getCumAvg(var);
                    calculateMetrics(ind);
                    setCumAvg(var);
                    
                    retry = 0;
                    ind = ind + 1;
                catch err
                    fprintf('\n BAD ITERATION ....... RETRY \n ');
                    fprintf(fvid, '\n BAD ITERATION ....... RETRY \n ');
                    %exception = MException.last;
                    msgString = getReport(err, 'extended');
                    fprintf('\n Exception Message: %s \n ', msgString);
                    fprintf(fvid, '\n Exception Message: %s \n ', msgString);
                    
                    retry = retry + 1;
                    if(retry == 3)
                        ind = ind + 1;
                        retry = 0;
                    end
                end

            end %% END of one set of simulations for all Tx probabilities %%%
            fprintf(fvid, '\n');
            fprintf(fvid, '%%%%%% END of iteration %d for policy%d %%%%%% \n', it, policy);
            fprintf(fvid, '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n');

            fprintf(fvid, '\n cumAvgThrptD2d%d :  \n', var);
            fprintf(fvid, '%d\t', cAvgD2dThrpt);
            fprintf('\n');
            fprintf(fvid, '\n cumAvgPsnr%d :  \n', var);
            fprintf(fvid, '%d\t', cAvgPsnr);
            fprintf('\n');
            fprintf(fvid, '\n cumAvgDetProb%d :  \n', var);
            fprintf(fvid, '%d\t', cAvgDetProb);
            fprintf('\n');
            fprintf(fvid, '\n cumTotDetProb%d :  \n', var);
            fprintf(fvid, '%d\t', cTotDetProb);
            fprintf('\n');
            
            var = var + 1;
            
        end
        fprintf(fvid, '\n');
        fprintf(fvid, '%%%%%% END of iteration %d for all Policies %%%%%% \n', it);
        fprintf(fvid, '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n');
        %%%% end of cases for base, different deltas %%%%%

    end %%%% end for n set of simulations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    fprintf(fvid, '\n\n ############ END of %d set of simulations ############ \n', it);

    
    cAvgPsnr=zeros(1,repn);
    cAvgUethrpt=zeros(1,repn);
    cAvgD2dThrpt=zeros(1,repn);
    cAvgDetProb=zeros(1,repn);
    cTotDetProb =zeros(1,repn);

    for q=1:3
        var = q -1;
        getCumAvg(var);
        fprintf(fvid, '\n cumAvgThrptD2d%d :  \n', var);
        fprintf(fvid, '%d\t', cAvgD2dThrpt);
        fprintf('\n');
        fprintf(fvid, '\n cumAvgPsnr%d :  \n', var);
        fprintf(fvid, '%d\t', cAvgPsnr);
        fprintf('\n');
        fprintf(fvid, '\n cumAvgDetProb%d :  \n', var);
        fprintf(fvid, '%d\t', cAvgDetProb);
        fprintf('\n');
        fprintf(fvid, '\n cumTotDetProb%d :  \n', var);
        fprintf(fvid, '%d\t', cTotDetProb);
        fprintf('\n');


        figure(1);
        hold on;
        plot(cAvgD2dThrpt, cAvgPsnr);
        xlabel('D2D Throughput (%)');
        ylabel('Video PSNR');
        title('PSNR vs D2D throughput');
        legend('base case', 'fixed transmission control', 'backoff control with ?=5');


        figure(2);
        hold on;
        plot(cAvgD2dThrpt, cAvgDetProb);
        xlabel('D2D Throughput (%)');
        ylabel('Object Detection probability');
        title('Detection Probability vs D2D throughput');
        legend('base case', 'fixed transmission control', 'backoff control with ?=5');
        
        figure(3);
        hold on;
        plot(cAvgD2dThrpt, cTotDetProb);
        xlabel('D2D Throughput (%)');
        ylabel('Total Detection probability');
        title('Total Detection Probability vs D2D throughput');
        legend('base case', 'fixed transmission control', 'backoff control with ?=5');
        
    end
    
    fclose(fvid);
    
    %%%%%%%%%%%%% END of execution of main function %%%%%%%%%%%%%%%%
  
    
    %%%%%%%%%%%%%%%%%%%%%%% OTHER FUNCTION DEFINITIONS %%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    function getCumAvg(var)
        if(var == 0)
            cAvgPsnr = cumAvgPsnr;
            cAvgUethrpt = cumAvgThrptUe;
            cAvgD2dThrpt = cumAvgThrptD2d;
            cAvgDetProb = cumAvgDetProb;
            cTotDetProb = cumTotDetProb;
        elseif(var == 1)
            cAvgPsnr = cumAvgPsnr1;
            cAvgUethrpt = cumAvgThrptUe1;
            cAvgD2dThrpt = cumAvgThrptD2d1;
            cAvgDetProb = cumAvgDetProb1;
            cTotDetProb = cumTotDetProb1;
        elseif(var == 2)
            cAvgPsnr = cumAvgPsnr2;
            cAvgUethrpt = cumAvgThrptUe2;
            cAvgD2dThrpt = cumAvgThrptD2d2;
            cAvgDetProb = cumAvgDetProb2;
            cTotDetProb = cumTotDetProb2;
        elseif(var == 3)
            cAvgPsnr = cumAvgPsnr3;
            cAvgUethrpt = cumAvgThrptUe3;
            cAvgD2dThrpt = cumAvgThrptD2d3;
            cAvgDetProb = cumAvgDetProb3;
            cTotDetProb = cumTotDetProb3;
        else
            cAvgPsnr = cumAvgPsnr4;
            cAvgUethrpt = cumAvgThrptUe4;
            cAvgD2dThrpt = cumAvgThrptD2d4;
            cAvgDetProb = cumAvgDetProb4;
            cTotDetProb = cumTotDetProb4;
        end
    end


    function setCumAvg(var)
        if(var == 0)
            cumAvgPsnr = cAvgPsnr;
            cumAvgThrptUe = cAvgUethrpt;
            cumAvgThrptD2d = cAvgD2dThrpt;
            cumAvgDetProb = cAvgDetProb;
            cumTotDetProb = cTotDetProb;
        elseif(var == 1)
            cumAvgPsnr1 = cAvgPsnr;
            cumAvgThrptUe1 = cAvgUethrpt;
            cumAvgThrptD2d1 = cAvgD2dThrpt;
            cumAvgDetProb1 = cAvgDetProb;
            cumTotDetProb1 = cTotDetProb;
        elseif(var == 2)
            cumAvgPsnr2 = cAvgPsnr;
            cumAvgThrptUe2 = cAvgUethrpt;
            cumAvgThrptD2d2 = cAvgD2dThrpt;
            cumAvgDetProb2 = cAvgDetProb;
            cumTotDetProb2 = cTotDetProb;
        elseif(var == 3)
            cumAvgPsnr3 = cAvgPsnr;
            cumAvgThrptUe3 = cAvgUethrpt;
            cumAvgThrptD2d3 = cAvgD2dThrpt;
            cumAvgDetProb3 = cAvgDetProb;
            cumTotDetProb3 = cTotDetProb;
        else
            cumAvgPsnr4 = cAvgPsnr;
            cumAvgThrptUe4 = cAvgUethrpt;
            cumAvgThrptD2d4 = cAvgD2dThrpt;
            cumAvgDetProb4 = cAvgDetProb;
            cumTotDetProb4 = cTotDetProb;
        end
    end
    
    function calculateMetrics(ind)                
        %%%%%%%%%%%%%%%%%%%% Calculate metrics %%%%%%%%%%%%%%%%%%%%%
        fprintf(fvid, '\n Number of badpackets for Tx_prob: %.3f = %d\t', tx_prob, sum(badpackets));
        fprintf(fvid, '\n Number of refernce badpackets for Tx_prob: %.3f = %d\t', tx_prob, sum(badpackets_ref));
        if (isequal(badpackets, badpackets_ref))
            fprintf(fvid, '\n D2D 0 same as reference case \t');
        else
            fprintf(fvid, '\n D2D 0 is NOT same as reference case \t');
        end
        
        t_ue = (num_total_pkts - num_bad_pkts) * 100 /num_total_pkts;
        cAvgUethrpt(ind) = ((it-1)*cAvgUethrpt(ind) + t_ue)/it;
        fprintf(fvid, '\n\t UE Throughput percentage for TX_prob = %.3f at it : %d = %.3f \n', tx_prob, it, t_ue);
        fprintf(fvid, '\t Cumulative average UE Throughput percentage for TX_prob = %.3f = %.3f \n', tx_prob, cAvgUethrpt(ind));
        

        t_dd = (num_total_pkts - num_bad_pkts_d2d) * 100 /num_total_pkts;
        cAvgD2dThrpt(ind) = ((it-1)*cAvgD2dThrpt(ind) + t_dd)/it;
        fprintf(fvid, '\t D2D Throughput percentage for Tx_prob: %.3f at it : %d = %.3f \n', tx_prob, it, t_dd);
        fprintf(fvid, '\t Cumulative average D2D Throughput percentage for Tx_prob: %.3f = %.3f \n', tx_prob, cAvgD2dThrpt(ind));
        
       
        fprintf('Now Calculate the PSNR w.r.t damaged packets \n');
        [m_psnr] = bunny_psnr(badpackets_ref, badpackets);
        fprintf(fvid, '\t PSNR obtained for Tx_prob: %.3f at iteration: %d = %.3f \n', tx_prob, it, m_psnr);

        if(cAvgPsnr(ind) == 0)
            cAvgPsnr(ind) = m_psnr;
        else
            if(m_psnr ~= 0)
               cAvgPsnr(ind) = ((it-1)*cAvgPsnr(ind) + m_psnr)/it;
               fprintf(fvid, '\t Average cumulative PSNR for Tx_prob: %.3f = %.3f \n', tx_prob, cAvgPsnr(ind));
            end
        end
        
        
        if(tx_prob == 0)
            detection_prob = 1;
            detTotProb = 1;
        else
            [detection_prob, detTotProb] = detect_object_feature();
        end
        fprintf(fvid, '\t Avg. object detection probability for Tx_prob: %.3f at iteration: %d = %.3f \n', tx_prob, it, detection_prob);
        fprintf(fvid, '\t Total object detection probability for Tx_prob: %.3f at iteration: %d = %.3f \n', tx_prob, it, detTotProb);
        
        if(detection_prob ~= 0)
           cAvgDetProb(ind) = ((it-1)*cAvgDetProb(ind) + detection_prob)/it;
           fprintf(fvid, '\t Average cumulative Detection Probability for Tx_prob: %.3f = %.3f \n', tx_prob, cAvgDetProb(ind));
        end
        
        if(detTotProb ~= 0)
           cTotDetProb(ind) = ((it-1)*cTotDetProb(ind) + detTotProb)/it;
           fprintf(fvid, '\t Total cumulative Detection Probability for Tx_prob: %.3f = %.3f \n', tx_prob, cTotDetProb(ind));
        end
          
        
    end
    
    
    
    function policy0(ind, next_state_u, next_state_i)
        tx_prob = ptrans(ind);
        if(rand() < tx_prob)    %% d2d may decides to transmit
            if(decodable_mat_ue(next_state_u, next_state_i) == 0)
                badpackets(pkt_idx) = 1;
                num_bad_pkts = num_bad_pkts + 1;
                lte_curr = 0;
            else
                lte_curr = 1;
            end
            
            if(decodable_mat_d2d(next_state_u, next_state_i) == 0)
                num_bad_pkts_d2d = num_bad_pkts_d2d + 1;
            end
            
            %%%% Reference video transmission
            if(decodable_mat_ue_ref(next_state_u, next_state_i) == 0)
                badpackets_ref(pkt_idx) = 1;
                num_bad_pkts_ref = num_bad_pkts_ref + 1;
            end
            
        else
           d2d_idle(next_state_u, next_state_i);
        end
        
        
        
    end %% end of function policy0(ind, next_state_u, next_state_i)

    
    function d2d_idle(next_state_u, next_state_i)
        %%next_state_i = 1;
        %%row_cum_trans_prob_i = T_mat_i(next_state_i, :);
        fprintf('D2D remains idle for packet# %d for frame = %d due to LTE failure\n', num_total_pkts, curr_frame);
        num_bad_pkts_d2d = num_bad_pkts_d2d + 1;
        
        if(decodable_mat_ue_ref(next_state_u, next_state_i) == 0)
            badpackets(pkt_idx) = 1;
            num_bad_pkts = num_bad_pkts + 1;
            lte_curr = 0;
        else
            lte_curr = 1;
        end
        
        %%%% Reference video transmission
        if(decodable_mat_ue_ref(next_state_u, next_state_i) == 0)
            badpackets_ref(pkt_idx) = 1;
            num_bad_pkts_ref = num_bad_pkts_ref + 1;
        end    
    end


    function [ue_state, in_state] =  getSteadyState()
        %%%%%% Get the Initial state %%%%%%%
  
         a = quantum_prob;
         for st=1:num_quantize
             cum_steady_state_prob(st) = a;
             a = a+quantum_prob;
         end

         r=rand();
         for ix=1:(num_quantize-1)
             if(r <= quantum_prob)
                 ue_state = 1;
                 break;
             else
                 cur=cum_steady_state_prob(ix);
                 nxt=cum_steady_state_prob(ix+1);
                 if((r > cur) && (r <= nxt))
                     ue_state = ix+1;
                     break;
                 end
             end
         end
         
         r=rand();
         for ix=1:(num_quantize-1)
             if(r <= quantum_prob)
                 in_state = 1;
                 break;
             else
                 cur=cum_steady_state_prob(ix);
                 nxt=cum_steady_state_prob(ix+1);
                 if((r > cur) && (r <= nxt))
                     in_state = ix+1;
                     break;
                 end
             end
         end

         fprintf ('The first state is : (%d, %d) -->', ue_state, in_state);
    end
    
    
    function [next_state_u, next_state_i] = getNextState()
        r=rand();
        for in=1:(num_quantize-1)
             if(r <= row_cum_trans_prob_u(1))
                 next_state_u = 1;
                 break;
             else
                 cur=row_cum_trans_prob_u(in);
                 nxt=row_cum_trans_prob_u(in+1);
                 if((r > cur) && (r <= nxt))
                     next_state_u = in+1;
                     break;
                 end
             end
         end %% end of for %

         r=rand();
         for in=1:(num_quantize-1)
             if(r <= row_cum_trans_prob_i(1))
                 next_state_i = 1;
                 break;
             else
                 cur=row_cum_trans_prob_i(in);
                 nxt=row_cum_trans_prob_i(in+1);
                 if((r > cur) && (r <= nxt))
                     next_state_i = in+1;
                     break;
                 end
             end
         end %% end of for %
    end
    
    
    
    toc;
    
end
