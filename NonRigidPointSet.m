function [NewPoints,G,W]=NonRigidPointSet(X,Y,D,M,N,sig2,w,W,G,alpha)
    NewPoints=ones(N,D);
    
    iter=0;

    gifname = 'nonrigid.gif';

    while (iter<1000) & (abs(sig2)>10*eps) 
        
        iter=iter+1;
              
        %E-Step:

        % Fill P with all cominations of (x_n - sRy_m+t)^2.
        T = bsxfun(@plus, Y,G*W);
        A = sum(X .* X, 2);
        B_m = -2*X*T';
        C = sum(T .* T, 2);
        P = bsxfun(@plus, A, B_m);
        P = bsxfun(@plus, C', P);
        
     
        % Transform every element p=exp(-1/(2*sig2) * p);
        P = exp(-P/(2*sig2));
        % The denominator is specific to each column. Sum over the rows.
        % Add constant term
        denom = sum(P,2) + (2*pi*sig2)^(D/2)*w/(1-w)*M/N;
        assert(length(denom) == N);
        % Divide each column by the denominator for that column.
        P = bsxfun(@rdivide, P, denom);
        P=P';
        
        E = 0;
        for n=1:N
            tmp = 0;
            for m=1:M
                tmp = tmp + P(n,m)/ sum(P(:,m));
            end
            E = E - log(tmp/M);
        end

        
        %M-Step:
        diagonal=diag(P*ones(N,1));
        diagonal1=diag(P'*ones(M,1));
        
        idP = diag(1 ./ (P*ones(N, 1)));
        A =  G + alpha * sig2 * idP;
        B = idP* P*X - Y;
        
%         A=(diagonal*G+alpha*sig2);%*((diagonal)^-1));
%         B=P*X-diagonal*Y;
        rcond(A)
        W=A\B;
         
        
        N_P=ones(M,1)'*P*ones(N,1);
        T=Y+G*W;

        sig2=1/(N_P*D)*(trace(X'*diagonal1*X)-2*trace((P*X)'*T)+trace(T'*diagonal*T));

        NewPoints=T;
        
           if (1)
%             Nrm=normalise(X);
            plot(X(:,1),X(:,2),'or'); hold on;
            plot(NewPoints(:,1),NewPoints(:,2),'*g'); hold off;
            title( sprintf('k = %d', iter) );
            axis off;
    
            % Copy the new image in the gif.
            frame = getframe(1);
            new_image = frame2im(frame);
            [ind_image,colormap] = rgb2ind(new_image,256);
            if iter == 1
                imwrite(ind_image,colormap,gifname,'gif','LoopCount',Inf,'DelayTime',1);
            else
                imwrite(ind_image,colormap,gifname,'gif','WriteMode','append','DelayTime',1);
            end
        end
        
    end
end