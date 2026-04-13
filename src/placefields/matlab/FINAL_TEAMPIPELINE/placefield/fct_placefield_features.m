% Author(s): Marti Geoffrey
% Epsztein Lab 2019

% Based on the firing rate map 'fr_tr' (row: time bin (lap), col: spatial bin), a labeled vector of place
% fields per lap ('ispf_tx') and a labeled vector of corresponding mean place field('ismeanpf_x'),
% this function returns a set of features for each mean place field as a structure 'pf':
%       ->  ratein               = mean firing rate inside the place field
%       ->  rateout              = mean firing rate outside the place field
%       ->  outovin              = rateout/ratein ratio
%       ->  outovinc             = corrected rateout/ratein ratio
%       ->  size                 = place field width in spatial bins
%       ->  size2                = another measure of place field width in spatial bins
%       ->  height               = place field height
%       ->  ifrmax               = index of spatial bin with hightest firing rate for the place field of interest
%       ->  devcenter_index      = index to quantify how much a place field peak is shifted from its center
%       ->  lap_size_t           = vector of place fields width per lap
%       ->  lap_ifrmax_t         = index of spatial bin with place field peak for each lap
%       ->  lap_dev              = standard deviation of the place field peaks per lap
%       ->  lap_maxdev           = max deviation of the place field peaks per lap
%       ->  lap_devwithmean      = deviation of the place field peaks per lap compared to the mean place field
%       ->  lap_size             = mean width of place fields per lap
%       ->  lap_varsize          = variability of place field width per lap
%       ->  lap_diffsizewithmean = difference between mean place field width and mean width per lap
%       ->  lap_pct              = percentage of laps with detected fields
% Each feature vector is calculated for all place fields detected, and ordered from the place field with the highest peak (first component)
% to the lowest (last component)


function pf = fct_placefield_features(fr_tx, ispf_tx, ismeanpf_x)


meanfr_x = nanmean(fr_tx, 1);
nb_pf = max(ismeanpf_x);
[nb_tbin, nb_xbin] = size(ispf_tx);

N = 5; 
pf.ispf = false(1, N);
pf.ratein = NaN(1, N);
pf.rateout = NaN(1, N);
pf.outovin = NaN(1, N);
pf.outovinc = NaN(1, N);
pf.size = NaN(1, N);
pf.size2 = NaN(1, N);
pf.height = NaN(1, N);
pf.ifrmax = NaN(1, N);
pf.devcenter_index = NaN(1, N);
pf.lap_dev = NaN(1, N);
pf.lap_maxdev = NaN(1, N);
pf.lap_devwithmean = NaN(1, N);
pf.lap_size = NaN(1, N);
pf.lap_varsize = NaN(1, N);
pf.lap_diffsizewithmean = NaN(1, N);
pf.lap_pct = NaN(1, N);

pf.lap_size_t = NaN(nb_tbin, N);
pf.lap_ifrmax_t = NaN(nb_tbin, N);

if nb_pf > 0
    pf.ispf(1:nb_pf) = true;
end

for p = 1:nb_pf
    ind_in = find(ismeanpf_x == p);
    ind_out = find(ismeanpf_x == 0);
    pf.ratein(p) = nanmean(meanfr_x(ind_in));
    pf.rateout(p) = nanmean(meanfr_x(ind_out));
    pf.outovin(p) = pf.rateout(p) / pf.ratein(p);
    baseline = mean(meanfr_x(meanfr_x <= prctile(meanfr_x, 10)));
    meanfr_x_corr = meanfr_x - baseline;
    ratein = nanmean(meanfr_x_corr(ind_in));
    rateout = nanmean(meanfr_x_corr(ind_out));
    pf.outovinc(p) = rateout / ratein;
    
    pf.size(p) = length(ind_in);
    pf.size2(p) = length(find(meanfr_x(ind_in) > 0.60*max(meanfr_x)));
    pf.height(p) = nanmax(meanfr_x(ind_in)) - nanmin(meanfr_x(ind_in));
    
    [~, imax] = nanmax(meanfr_x(ind_in));
    pf.ifrmax(p) = imax + ind_in(1) - 1;
    icenter = floor(length(ind_in)/2);
    pf.devcenter_index(p) = abs(imax - icenter) / (pf.size(p) / 2);
    
    
    ispf_tx_tmp = ispf_tx == p;
    [istart, istop] = fct_find_startstop(ispf_tx_tmp);
    pf.lap_size_t(:, p) = istop - istart;
    
    idx1 = repmat((1:nb_tbin)', 1, nb_xbin);
    idx2 = ispf_tx_tmp + 1;
    maxfr = accumarray([idx1(:) idx2(:)], fr_tx(:), [nb_tbin, 2], @nanmax, NaN);
    
    %%%% do not work with me!!!!!
    [~, imax] = nanmax(fr_tx == maxfr(:, 2), [], 2);
    
    
    %{
    %I replace by : 
    crap = nan(length(fr_tx(:,1)),1);
    for yy=1:length(fr_tx(:,1))
        yy
        tp = find(fr_tx(yy,:)==maxfr(yy,2));
        crap(yy)=tp;
    end
    %}
   
    imax(isnan(maxfr(:, 2))) = NaN;
    pf.lap_ifrmax_t(:, p) = imax;
    
    pf.lap_dev(p) = nanstd(pf.lap_ifrmax_t(:, p));
    pf.lap_maxdev(p) = nanmax(pf.lap_ifrmax_t(:, p)) - nanmin(pf.lap_ifrmax_t(:, p));
    pf.lap_devwithmean(p) = sqrt(nansum(((pf.ifrmax(p) - pf.lap_ifrmax_t(:, p)).^2) / (nb_tbin)));
    
    pf.lap_size(p) = nanmean(pf.lap_size_t(:, p));
    pf.lap_varsize(p) = nanstd(pf.lap_size_t(:, p));
    pf.lap_diffsizewithmean(p) = pf.size(p) - pf.lap_size(p);
    
    pf.lap_pct(p) = sum(sum(ispf_tx_tmp, 2) > 0) / nb_tbin;
end