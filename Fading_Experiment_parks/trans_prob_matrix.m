function [T_mat_u, T_mat_i] = trans_prob_matrix()

num_quantize = 100;
quantum_prob = (1/num_quantize);

ps = zeros(1,num_quantize);
ps(num_quantize) = 5;

pi = 3.14;
b = sqrt(2/pi);
m=0;
for i=1:(num_quantize-1)
    m=m+quantum_prob;
    y = raylinv(m,b);
    ps(i) = y;
end

display(ps);

%k=0;
 %%%%%%%%% START for loop %%%%%%%%%
%for d = 0:quantum_prob:(1-quantum_prob)
    
    trans_mat = zeros(num_quantize,num_quantize);

    %%%%%%%%%%%%%% Creating joint probability Distribution %%%%%%%%%%%%%%

    f = @(x) ((pi/2) .* x .* (exp((-pi) .* (x.^2) ./4)));
    pxxx = integral(f, 1.9535, 5);
    %display(pxxx);
    
    
    k=0.6;
    z = (pi./2) .* (k./(1-(k.^2))) ;
    %display(z);


    %%%% Joint pdf %%%%%%%
    fun = @(x,y) ((pi.^2)/4) .* ((x .* y)./(1-(k.^2)))  .* (exp((-(pi)/4) .* ((x.^2 + y.^2) ./ (1-(k.^2))))) .* besseli(0,(z .* (x.*y)));

    for i=1:num_quantize
        for j=1:num_quantize
            if(i==1) 
                xmin=0;
            else
                xmin=ps(i-1);
            end
            if(j==1)
                ymin=0;
            else
                ymin=ps(j-1);
            end
            xmax = ps(i);
            ymax = ps(j);

            Pxy = integral2(fun, xmin, xmax, ymin, ymax, 'AbsTol',eps,'RelTol',100*eps);

            %display(Pxy);

            Px = quantum_prob; %%% Quantized probability %%%

            Pygx = Pxy/Px;

            %display(Pygx);

            trans_mat(i,j) = Pygx;

        end
    end


    cum_trans_prob = zeros(num_quantize,num_quantize);

    for i=1:num_quantize
        row_sum=0;
        for j=1:num_quantize
            row_sum = row_sum + trans_mat(i,j);
            cum_trans_prob(i,j) = row_sum;
        end
        %display(row_sum);
    end

    %{
    new_mat1 = zeros(num_quantize,num_quantize);
    new_mat2 = zeros(num_quantize,num_quantize);
    new_mat1 = trans_mat * trans_mat;
    new_mat2 = new_mat1 * new_mat1;
    
    %%%%%%%%% END for loop %%%%%%%%%
    
    hold on;
    plot(trans_mat(10, :));
    %}
    
    T_mat_u = cum_trans_prob;
    T_mat_i = cum_trans_prob;
    
    
%end








        