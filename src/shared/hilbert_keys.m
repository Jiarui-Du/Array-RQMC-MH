function keys = hilbert_keys(U,bits)
%HILBERT_KEYS Encode each point independently, so unchanged keys can be cached.
% In d=1 retain the exact coordinate: the historical order is ordinary sort,
% not a finite-bit cell order. In d>1 use a uint64 scalar for d*bits<=64,
% otherwise big-endian uint32 words (the final partial word is right-aligned).
%
% U is an N-by-d array in [0,1]. Each output row belongs to the same input row;
% call hilbert_order_from_keys to obtain a permutation. Integer packing preserves
% key bits without floating-point conversion.
% Omitted/empty bits selects min(20,floor(64/d)); explicit depths support audits.

    %% Validate coordinates and retain exact one-dimensional keys
    if nargin<2 || isempty(bits),bits=resolve_hilbert_bits(size(U,2));end
    validateattributes(U,{'double'},{'2d','real','finite','nonempty','>=',0,'<=',1});
    validateattributes(bits,{'numeric'},{'scalar','integer','>=',1,'<=',30});
    [n,d]=size(U);
    if d==1,keys=U;return;end

    %% Quantize coordinates and apply Hilbert transpose operations
    X=uint32(floor(min(U,1-eps)*2^bits));
    Q=bitshift(uint32(1),bits-1);
    while Q>1
        P=Q-1;
        for i=1:d
            mask=bitand(X(:,i),Q)~=0;
            % Full-column bit operations avoid repeated logical gathers.
            exchange=bitand(bitxor(X(:,1),X(:,i)),P.*uint32(~mask));
            X(:,1)=bitxor(X(:,1),bitor(exchange,P.*uint32(mask)));
            X(:,i)=bitxor(X(:,i),exchange);
        end
        Q=bitshift(Q,-1);
    end

    %% Apply Gray encoding and the final orientation correction
    for i=2:d,X(:,i)=bitxor(X(:,i),X(:,i-1));end
    t=zeros(n,1,'uint32');Q=bitshift(uint32(1),bits-1);
    while Q>1
        t=bitxor(t,(Q-1).*uint32(bitand(X(:,d),Q)~=0));
        Q=bitshift(Q,-1);
    end
    X=bitxor(X,t);

    %% Pack short keys into one unsigned 64-bit integer
    if bits*d<=64
        % Build in 32-bit temporaries to reduce per-bit memory traffic, then
        % merge exactly once. No full key ever passes through double.
        high=zeros(n,1,'uint32');low=zeros(n,1,'uint32');k=0;
        for b=bits:-1:1
            for i=1:d
                k=k+1;
                if k<=32
                    high=bitor(bitshift(high,1),bitget(X(:,i),b));
                else
                    low=bitor(bitshift(low,1),bitget(X(:,i),b));
                end
            end
        end
        if bits*d<=32
            keys=uint64(high);
        else
            keys=bitor(bitshift(uint64(high),bits*d-32),uint64(low));
        end
        return
    end

    %% Pack longer keys into lexicographically ordered words
    keys=zeros(n,ceil(bits*d/32),'uint32');k=0;
    for b=bits:-1:1
        for i=1:d
            k=k+1;word=ceil(k/32);
            keys(:,word)=bitor(bitshift(keys(:,word),1),uint32(bitget(X(:,i),b)));
        end
    end
end
