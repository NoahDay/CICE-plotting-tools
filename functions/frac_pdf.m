function output = frac_pdf(data,fstd)
% Make global fraction PDFS for ITD and FSD
if nargin < 2 || isempty(fstd)
    fstd = 0;
end
if fstd == 0 % 3 dim
    [~,~,depth] = size(data);
    idx = isnan(data(:,:,1)); % Land mask
    % ITD
    for i = 1:depth
        temp_data = data(:,:,i);
        num_pdf(i) = sum(sum(temp_data(~idx)));
    end
    if sum(num_pdf < eps)
        output = zeros(depth);
    else
        output = num_pdf/sum(num_pdf); % normalize
    end
elseif fstd == 1 % Plot p(ITD|FSD)
    % [LON, LAT, NFSD, NCAT) 
    [~,~,nfsd,ncat] = size(data);
    idx = isnan(data(:,:,1,1)); % Land mask
    for j = 1:nfsd
        % ITD
        for i = 1:ncat
            temp_data = data(:,:,i,:);
            num_pdf(i) = sum(sum(temp_data(~idx)));
        end
        if sum(num_pdf < eps)
            output(j,:) = zeros(1,ncat);
        else
            output(j,:) = num_pdf/sum(num_pdf); % normalize
        end
    end
    
    
elseif fstd == 2 % Plot p(FSD|ITD)
    [~,~,nfsd,ncat] = size(data);
    idx = isnan(data(:,:,1,1)); % Land mask
    for j = 1:ncat
        % ITD
        for i = 1:nfsd
            temp_data = data(:,:,i,:);
            num_pdf(i) = sum(sum(temp_data(~idx)));
        end
        if sum(num_pdf < eps)
            output(j,:) = zeros(1,nfsd);
        else
            output(j,:) = num_pdf/sum(num_pdf); % normalize
        end
    end
end

