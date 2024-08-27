clc;clear all;close all;

%% define file names and info
wys = 1979:2023;
SNOTEL_data_filename_PNNL = '/Volumes/Pruina_External_Elements/SnowData/SNOTEL/SNOTEL_all_obs_PNNLbcqc_WY1979-2021.nc';
new_filename = '/Volumes/Pruina_External_Elements/SnowData/SNOTEL/SNOTEL_all_obs_PNNLbcqc_WY1979-2023_with_snowh.nc';
site_ids = ncread(SNOTEL_data_filename_PNNL,'site_id');

%for F to K conversin:
func_F_to_K =@(tempF) (tempF-32).*5/9 + 273.15;

%% add snow depth using the raw data:
raw_data_dir = '/Volumes/Pruina_External_Elements/WPO_SnowModeling/Observations/SNOTEL/raw/1979_2023_Daily_All_Data/';
total_datelist = [datenum([1978 10 1]):datenum([2023 9 30])]';
total_datevecs = datevec(total_datelist);
ndates = length(total_datelist);
nsites = length(site_ids);

%initialize outputs of raw snotel data:
store_swe_mm_all_sites=nan(ndates,nsites);
store_snowh_mm_all_sites=nan(ndates,nsites);
store_tmean_K_all_sites=nan(ndates,nsites);
store_tmax_K_all_sites=nan(ndates,nsites);
store_tmin_K_all_sites=nan(ndates,nsites);
store_prcp_mm_all_sites=nan(ndates,nsites);

%get raw snotel data
for s=1:length(site_ids)
    s
    current_ID = site_ids(s);
    infilename = sprintf('STA%d-1979-2023.csv',current_ID);
    raw_data = readtable([raw_data_dir,infilename]);
    snowh_mm_data_all = nan(length(total_datelist),1);

    if height(raw_data)>100 %only consider sites with at least 100 data points:
        raw_data.Properties.VariableNames = ["dates","swe_in","prcp_in","T2_F_start","T2mean_F","T2max_F","T2min_F","snowh_in","QC_snowh","QA_snowh"];
        dates = datenum(raw_data.dates);

        %convert SWE to mm:
        swe_in = raw_data.swe_in;
        swe_mm = swe_in.*25.4;

        %convert SNOWH to mm:
        snowh_in = raw_data.snowh_in;
        snowh_mm = snowh_in.*25.4;

        %convert prcp to mm:
        prcp_in = raw_data.prcp_in;
        prcp_mm = prcp_in.*25.4;

        %convert temp to K:
        T2mean_F = raw_data.T2mean_F;
        T2max_F = raw_data.T2max_F;
        T2min_F = raw_data.T2min_F;

        T2mean_K = func_F_to_K(T2mean_F);
        T2max_K = func_F_to_K(T2max_F);
        T2min_K = func_F_to_K(T2min_F);

        [c,ia,ib] = intersect(total_datelist,dates);

        store_swe_mm_all_sites(ia,s) = swe_mm(ib);
        store_snowh_mm_all_sites(ia,s) = snowh_mm(ib);
        store_tmean_K_all_sites(ia,s) = T2mean_K(ib);
        store_tmax_K_all_sites(ia,s) = T2max_K(ib);
        store_tmin_K_all_sites(ia,s) = T2min_K(ib);
        store_prcp_mm_all_sites(ia,s) = prcp_mm(ib);
    end
end

%% reorganize swe data to be 3d:
start_idx = 1;
swe_raw=[];
for wyi=1:length(wys)
    current_wy = wys(wyi);
    current_dates = datenum([current_wy-1 10 1]):datenum([current_wy 9 30]);
    ndays = length(current_dates);
    day_idx_end = start_idx+ndays-1;
    annual_data = store_swe_mm_all_sites(start_idx:day_idx_end,:);
    if ndays ==365
        annual_data = [annual_data;nan(1,length(annual_data))];
    end
    swe_raw = cat(3,swe_raw,annual_data');

    start_idx = day_idx_end+1;
end
swe_QC = swe_raw;

%% reorganize snowh data to be 3d:
start_idx = 1;
snowh_raw=[];
for wyi=1:length(wys)
    current_wy = wys(wyi);
    current_dates = datenum([current_wy-1 10 1]):datenum([current_wy 9 30]);
    ndays = length(current_dates);
    day_idx_end = start_idx+ndays-1;
    annual_data = store_snowh_mm_all_sites(start_idx:day_idx_end,:);
    if ndays ==365
        annual_data = [annual_data;nan(1,length(annual_data))];
    end
    snowh_raw = cat(3,snowh_raw,annual_data');

    start_idx = day_idx_end+1;
end
snowh_QC = snowh_raw;

%% reorganize tmean data to be 3d:
start_idx = 1;
tmean_raw=[];
for wyi=1:length(wys)
    current_wy = wys(wyi);
    current_dates = datenum([current_wy-1 10 1]):datenum([current_wy 9 30]);
    ndays = length(current_dates);
    day_idx_end = start_idx+ndays-1;
    annual_data = store_tmean_K_all_sites(start_idx:day_idx_end,:);
    if ndays ==365
        annual_data = [annual_data;nan(1,length(annual_data))];
    end
    tmean_raw = cat(3,tmean_raw,annual_data');

    start_idx = day_idx_end+1;
end
tmean_QC = tmean_raw;

%% reorganize tmax data to be 3d:
start_idx = 1;
tmax_raw=[];
for wyi=1:length(wys)
    current_wy = wys(wyi);
    current_dates = datenum([current_wy-1 10 1]):datenum([current_wy 9 30]);
    ndays = length(current_dates);
    day_idx_end = start_idx+ndays-1;
    annual_data = store_tmax_K_all_sites(start_idx:day_idx_end,:);
    if ndays ==365
        annual_data = [annual_data;nan(1,length(annual_data))];
    end
    tmax_raw = cat(3,tmax_raw,annual_data');

    start_idx = day_idx_end+1;
end
tmax_QC = tmax_raw;

%% reorganize tmin data to be 3d:
start_idx = 1;
tmin_raw=[];
for wyi=1:length(wys)
    current_wy = wys(wyi);
    current_dates = datenum([current_wy-1 10 1]):datenum([current_wy 9 30]);
    ndays = length(current_dates);
    day_idx_end = start_idx+ndays-1;
    annual_data = store_tmin_K_all_sites(start_idx:day_idx_end,:);
    if ndays ==365
        annual_data = [annual_data;nan(1,length(annual_data))];
    end
    tmin_raw = cat(3,tmin_raw,annual_data');

    start_idx = day_idx_end+1;
end
tmin_QC = tmin_raw;

%% reorganize prcp data to be 3d:
start_idx = 1;
prcp_raw=[];
for wyi=1:length(wys)
    current_wy = wys(wyi);
    current_dates = datenum([current_wy-1 10 1]):datenum([current_wy 9 30]);
    ndays = length(current_dates);
    day_idx_end = start_idx+ndays-1;
    annual_data = store_prcp_mm_all_sites(start_idx:day_idx_end,:);
    if ndays ==365
        annual_data = [annual_data;nan(1,length(annual_data))];
    end
    prcp_raw = cat(3,prcp_raw,annual_data');

    start_idx = day_idx_end+1;
end
prcp_QC = prcp_raw;

%% QC process from Yan et al. (2018) as done for the PNNL product
%ensure temperature values are reasonable:
%w/in +/- 40C
idx_too_hot = find(tmin_QC>313.15 | tmean_QC>313.15 | tmax_QC>313.15);
idx_too_cold = find(tmin_QC<233.15 | tmean_QC<233.15 | tmax_QC<233.15);
idx_bad_temp = unique([idx_too_hot;idx_too_cold]);
tmean_QC(idx_bad_temp) = NaN;
tmax_QC(idx_bad_temp) = NaN;
tmin_QC(idx_bad_temp) = NaN;

%stage 1 and 2 are from Serreze et al. (1999) to identify missing, erroneous and outlies in precip, air temp & SWE:
%(2a) ensure non-negative values:
idx_neg_SWE=find(swe_QC<0);
swe_QC(idx_neg_SWE)= NaN;

idx_neg_prcp=find(prcp_QC<0);
prcp_QC(idx_neg_prcp)= NaN;

nyears = length(wys);
for yi=1:nyears
    %(1a)Screen data from the total year if missing values were inserted for precip or SWE in first 15 days of a water year
    current_PRCP_days1_15 = squeeze(prcp_QC(:,1:15,yi));
    current_SWE_days1_15 = squeeze(swe_QC(:,1:15,yi));

    MISSING_PRCP = isnan(current_PRCP_days1_15);
    MISSING_SWE = isnan(current_SWE_days1_15);

    MISSING_PRCP = sum(MISSING_PRCP,2);
    MISSING_SWE=sum(MISSING_SWE,2);

    idx_bad_prcp = find(MISSING_PRCP>0);
    idx_bad_SWE = find(MISSING_SWE>0);

    prcp_QC(idx_bad_prcp,:,yi) = NaN;
    swe_QC(idx_bad_SWE,:,yi) = NaN;

    %(1b) or if accum PRE on oct 1 is >5 in (127mm)
    current_PRCP_day1 = squeeze(prcp_QC(:,1,yi));
    idx_bad_prcp = find(current_PRCP_day1>127);
    if length(idx_bad_prcp) >0
        prcp_QC(idx_bad_prcp,:,yi) = NaN;
    end
end

%(2) - still using cumulative prcp
S = size(prcp_QC);
for si=1:S(1)
    for yi=1:S(3)
        current_prcp = squeeze(prcp_QC(si,:,yi));
        current_swe = squeeze(swe_QC(si,:,yi));
        for di=2:length(current_prcp)-1
            prcp_N = current_prcp(di);
            prcp_Nminus1 = current_prcp(di-1);
            prcp_Nplus1 = current_prcp(di+1);
            %(2a) if precip on day N is higher or lower than N-1 or N+1 when N-1=N+1, set N to N-1:
            if prcp_Nplus1==prcp_Nminus1 && prcp_N~=prcp_Nminus1
                current_prcp(di) = prcp_Nminus1;
            end

            %(2b) if N+1 - N-1 > 0; but N<N-1 | N+1
            if ( (prcp_Nplus1-prcp_Nminus1) > 0) && ( (prcp_N<prcp_Nminus1) || (prcp_N>prcp_Nplus1) )
                if (prcp_Nplus1-prcp_Nminus1) < 12.7 %if the difference is small (0.5in, set to mean of N+1 & N-1)
                    prcp_N = (prcp_Nplus1+prcp_Nminus1)/2;
                else
                    prcp_N = NaN; %if difference is too large, screen the data point
                end
            end
            current_prcp(di) = prcp_N;

            %ensure prcp increments are reasonable (no too large / spikey)
            prcp_N = current_prcp(di);
            prcp_Nminus1 = current_prcp(di-1);

            swe_N = current_swe(di);
            swe_Nminus1 = current_swe(di-1);
            swe_Nplus1 = current_swe(di+1);
            dSWE1 = swe_N-swe_Nminus1;
            dSWE2 = swe_Nplus1-swe_N;

            if prcp_N-prcp_Nminus1 > 254
                current_prcp(di)=NaN;
            end
            if dSWE1>254
                current_swe(di) = NaN;
            end
            if prcp_N-prcp_Nminus1 < -12.7 %if prcp shows a significant decrease
                current_prcp(di)=NaN;
            elseif prcp_N-prcp_Nminus1 < 0
                current_prcp(di)=prcp_Nminus1;
            end

            %check to make sure the daily increment values for SWE are acceptable
            if dSWE1 < -254 %screen large decrease in SWE
                current_swe(di) = NaN;
            end
            if (abs(dSWE1) > 63.5 && abs(dSWE2)>63.5) && (dSWE1*dSWE2 < 0) %screen spike in SWE
                current_swe(di) = NaN;
            end
        end
        prcp_QC(si,:,yi) = current_prcp;
        swe_QC(si,:,yi) = current_swe;
    end
end

%convert QC prcp from cumulative prcp to mm/day:
S = size(prcp_QC);
prcp_QC_daily = nan(S);
for yi=1:S(3)
    WY = wys(yi);
    wy_dates = datenum([WY-1 10 1]):datenum([WY 9 30]);
    assert(S(3) == length(wys),'dim mistmatch')
    ndays = length(wy_dates);
    for si=1:S(1)
        current_prcp = squeeze(prcp_QC(si,:,yi));
        for di=1:ndays
            if di<ndays
                current_prcp_N = current_prcp(di);
                current_prcp_Nplus1 = current_prcp(di+1);
                daily_prcp = current_prcp_Nplus1-current_prcp_N;
                prcp_QC_daily(si,di,yi) = daily_prcp;
            elseif di==ndays && yi<S(3)
                current_prcp_N = current_prcp(di);
                current_prcp_Nplus1 = squeeze(prcp_QC(si,1,yi+1));
                daily_prcp = current_prcp_Nplus1;
                prcp_QC_daily(si,di,yi) = daily_prcp;
            elseif di==ndays && yi==S(3)
                prcp_QC_daily(si,di,yi) = NaN; %last day of record will be recorded as NaN
            end
        end

    end
end

%do the same for the raw prcp data, converting to mm/day
S = size(prcp_raw);
prcp_raw_daily = nan(S);
for yi=1:S(3)
    WY = wys(yi);
    wy_dates = datenum([WY-1 10 1]):datenum([WY 9 30]);
    assert(S(3) == length(wys),'dim mistmatch')
    ndays = length(wy_dates);
    for si=1:S(1)
        current_prcp = squeeze(prcp_raw(si,:,yi));
        for di=1:ndays
            if di<ndays
                current_prcp_N = current_prcp(di);
                current_prcp_Nplus1 = current_prcp(di+1);
                daily_prcp = current_prcp_Nplus1-current_prcp_N;
                prcp_raw_daily(si,di,yi) = daily_prcp;
            elseif di==ndays && yi<S(3)
                current_prcp_N = current_prcp(di);
                current_prcp_Nplus1 = squeeze(prcp_raw(si,1,yi+1));
                daily_prcp = current_prcp_Nplus1;
                prcp_raw_daily(si,di,yi) = daily_prcp;
            elseif di==ndays && yi==S(3)
                prcp_raw_daily(si,di,yi) = NaN; %last day of record will be recorded as NaN
            end
        end
    end
end

%screen instances of 0 change in temp from day N to N+1
%for tmean
for yi=1:S(3)
    for si=1:S(1)
        current_temp = squeeze(tmean_QC(si,:,yi));
        for di=1:length(current_temp)-1
            temp_N = current_temp(di);
            temp_Nplus1 = current_temp(di+1);
            if temp_N == temp_Nplus1
                tmean_QC(si,di,yi) = NaN;
                tmean_QC(si,di+1,yi) = NaN;
            end
        end
    end
end

%for tmin
for yi=1:S(3)
    for si=1:S(1)
        current_temp = squeeze(tmin_QC(si,:,yi));
        for di=1:length(current_temp)-1
            temp_N = current_temp(di);
            temp_Nplus1 = current_temp(di+1);
            if temp_N == temp_Nplus1
                tmin_QC(si,di,yi) = NaN;
                tmin_QC(si,di+1,yi) = NaN;
            end
        end
    end
end

%for tmax
for yi=1:S(3)
    for si=1:S(1)
        current_temp = squeeze(tmax_QC(si,:,yi));
        for di=1:length(current_temp)-1
            temp_N = current_temp(di);
            temp_Nplus1 = current_temp(di+1);
            if temp_N == temp_Nplus1
                tmax_QC(si,di,yi) = NaN;
                tmax_QC(si,di+1,yi) = NaN;
            end
        end
    end
end

%% apply statistical measure filter for SWE and PREC based on Yan et al. (2018):
%(2)
S = size(prcp_QC_daily);
for si=1:S(1)
    current_prcp_total_timeseries = squeeze(prcp_QC_daily(si,:,:));
    current_swe_total_timeseries = squeeze(swe_QC(si,:,:));

    %perform Yan et al. (2018) check on multiyear mean peak SWE vs. prcp (>20% difference are screened)
    BAD_SWE=0;
    BAD_PRCP=1;
    store_data=nan(S(3),2);
    for yi=1:S(3)
        [PeakSWE,i_PeakSWE] = max(current_swe);
        if ~isnan(PeakSWE)
            total_prcp = cumsum(current_prcp(1:i_PeakSWE));
            store_data(yi,1) = PeakSWE;
            store_data(yi,2) = total_prcp(end);
        end
    end
    peakSWE_mean = nanmean(store_data(:,1));
    PRCP_mean = nanmean(store_data(:,2));
    if peakSWE_mean > (PRCP_mean+PRCP_mean*0.2)
        BAD_SWE = 1;
        BAD_PRCP = 1;
    end

    for yi=1:S(3)
        [si,yi]
        WY=wys(yi);
        WY_dates = datenum([WY-1 10 1]):datenum([WY 9 30]);
        WY_datevecs = datevec(WY_dates);
        current_prcp = squeeze(prcp_QC_daily(si,:,yi));
        current_swe = squeeze(swe_QC(si,:,yi));
        current_tmean = squeeze(tmean_QC(si,:,yi));
        current_tmax = squeeze(tmax_QC(si,:,yi));
        current_tmin = squeeze(tmin_QC(si,:,yi));

        %initialize flags
        current_prcp_flag = nan(size(current_prcp));
        current_swe_flag = nan(size(current_swe));
        current_tmean_flag= nan(size(current_tmean));
        current_tmin_flag= nan(size(current_tmin));
        current_tmax_flag= nan(size(current_tmax));

        idx=find(isnan(current_prcp));
        current_prcp_flag(idx) = 1;
        idx=find(isnan(current_swe));
        current_swe_flag(idx) = 1;
        idx=find(isnan(current_tmean));
        current_tmean_flag(idx) = 1;
        idx=find(isnan(current_tmin));
        current_tmin_flag(idx) = 1;
        idx=find(isnan(current_tmax));
        current_tmax_flag(idx) = 1;

        %check if peak SWE for this WY is 5% greater than prec, if so screen prcp & swe for total WY:
        [PeakSWE,i_PeakSWE] = max(current_swe);
        if ~isnan(PeakSWE)
            total_prcp = cumsum(current_prcp(1:i_PeakSWE));
            total_prcp = total_prcp(end);
            if PeakSWE > (total_prcp+0.05*total_prcp)
                current_prcp_flag(1:length(current_prcp_flag)) = 1; %1 flag = screen data
                current_swe_flag(1:length(current_swe_flag)) = 1; %1 flag = screen data
            end
        end

        %if this site's mean peak SWE>20% of accum season prcp, screen all data
        if BAD_SWE==1 && BAD_PRCP==1
            current_swe_flag = ones(size(current_swe)); %1 flag = screen data
            current_prcp_flag= ones(size(current_swe)); %1 flag = screen data
        end

        for di=2:length(current_prcp)-1
            %define multiyear 60-day climateology for mean and std calcs for screening surrounding day di:
            current_month = WY_datevecs(di,2);
            current_midday = round(eomday(WY,current_month)/2);
            if current_month<10
                current_month_idx = find(WY_dates == datenum([WY,current_month,current_midday]));
            else
                current_month_idx = find(WY_dates == datenum([WY-1,current_month,current_midday]));
            end
            start_idx = current_month_idx-30;
            end_idx = current_month_idx+30;
            idxs = start_idx:end_idx;
            i = find(idxs<=0);
            if ~isempty(i)
                idxs(i) = idxs(i) + 366;
            end
            i = find(idxs>366);
            if ~isempty(i)
                idxs(i) = idxs(i) - 366;
            end
            assert(min(idxs)>0 & max(idxs)<367,'invalid range')

            prcp_window = squeeze(prcp_QC_daily(si,idxs,:));
            prcp_window_vec = prcp_window(:);
            mean_prcp = nanmean(prcp_window_vec);
            std_prcp = nanstd(prcp_window_vec);

            tmean_window = squeeze(tmean_QC(si,idxs,:));
            tmean_window_vec = tmean_window(:);
            mean_tmean = nanmean(tmean_window_vec);
            std_tmean = nanstd(tmean_window_vec);

            tmin_window = squeeze(tmin_QC(si,idxs,:));
            tmin_window_vec = tmin_window(:);
            mean_tmin = nanmean(tmin_window_vec);
            std_tmin = nanstd(tmin_window_vec);

            tmax_window = squeeze(tmax_QC(si,idxs,:));
            tmax_window_vec = tmax_window(:);
            mean_tmax = nanmean(tmax_window_vec);
            std_tmax = nanstd(tmax_window_vec);

            swe_window = squeeze(swe_QC(si,idxs,:));
            swe_window_diff = diff(swe_window,[],1);
            swe_window_vec = swe_window_diff(:);
            idx_accum = find(swe_window_vec>0);
            idx_ablation = find(swe_window_vec<0);
            mean_swe_accum = nanmean(swe_window_vec(idx_accum));
            std_swe_accum = nanstd(swe_window_vec(idx_accum));
            mean_swe_ablat = nanmean(swe_window_vec(idx_ablation));
            std_swe_ablat = nanstd(swe_window_vec(idx_ablation));
            assert(mean_swe_accum>0 | isnan(mean_swe_accum),'not accumulation')
            assert(mean_swe_ablat<0 | isnan(mean_swe_ablat),'not ablation')

            %run statistical filter for prcp & SWE based on 5 stdev increment flag:
            %PRCP:
            prcp_N = current_prcp(di);
            %ensure prcp increments <5 stds
            if prcp_N > mean_prcp+std_prcp*5 || isnan(prcp_N)
                current_prcp_flag(di)=1; %large prcp (5 st. dev above window avg)
            elseif prcp_N > mean_prcp+std_prcp*3
                current_prcp_flag(di)=2; %large prcp (3 st. dev above window avg / to be used later)
            end

            %SWE:
            swe_N = current_swe(di);
            swe_Nminus1 = current_swe(di-1);
            delta_SWE = swe_N-swe_Nminus1;
            if delta_SWE > (mean_swe_accum+std_swe_accum*5)
                current_swe_flag(di)=1; %large accumulation
            elseif delta_SWE< (mean_swe_ablat-std_swe_ablat*5)
                current_swe_flag(di)=2; %large ablation
            end

            %temp:
            tmean_N = current_tmean(di);
            tmax_N = current_tmax(di);
            tmin_N = current_tmin(di);
            %linear interpolate NaN values based on Sun et al. (2019)
            if isnan(tmean_N)
                current_tmean(di) = (current_tmean(di-1) + current_tmean(di+1))/2;
                tmean_N = current_tmean(di);
            end

            if isnan(tmin_N)
                current_tmin(di) = (current_tmin(di-1) + current_tmin(di+1))/2;
                tmin_N = current_tmin(di);
            end

            if isnan(tmax_N)
                current_tmax(di) = (current_tmax(di-1) + current_tmax(di+1))/2;
                tmax_N = current_tmax(di);
            end

            %ensure tmean < 3 stds below avg
            if tmean_N > mean_tmean+std_tmean*3
                current_tmean_flag(di)=1; %large tmean
            elseif tmax_N > mean_tmax+std_tmax*3
                current_tmax_flag(di)=1; %large tmax
            elseif tmin_N > mean_tmin+std_tmin*3
                current_tmin_flag(di)=1; %large tmin
            end

            %determine screening from flags:
            %remove SWE & temp flags if SWE flag = 2 & tmean flag =1
            %(large ablation corresponding with high temp)
            idx=find(current_swe_flag==2 & current_tmean_flag==1);
            current_swe_flag(idx) = NaN;
            current_tmean_flag(idx) = NaN;
            current_tmax_flag(idx) = NaN;
            current_tmin_flag(idx) = NaN;

            %remove SWE & prcp flags if SWE flag = 1 & prcp flag=2
            %(large SWE increase corresponding with large prcp event)
            idx=find(current_prcp_flag>=1 & current_swe_flag==1);
            current_prcp_flag(idx)= NaN;
            current_swe_flag(idx) = NaN;

            %apply screened data across total WY if NaN in first 15 days, if
            %not, then apply NaNs from 1st day of invalid data:
            idx=find(~isnan(current_swe_flag));
            if ~isempty(idx) %if there are flags...
                if idx(1)>15 %if the first flag is after day 15, screen from the 1st flag day to end of WY
                    current_swe_flag(idx(1):end) = 1;
                else
                    current_swe_flag(1:end) = 1; %if 1st flag occurs before day 15, screen total WY
                end
            end
            idx=find(~isnan(current_prcp_flag));
            if ~isempty(idx)
                if idx(1)>15
                    current_prcp_flag(idx(1):end) = 1; %if the first flag is after day 15, screen from the 1st flag day to end of WY
                else
                    current_prcp_flag(1:end) = 1; %if 1st flag occurs before day 15, screen total WY
                end
            end
            idx=find(~isnan(current_tmean_flag));
            if ~isempty(idx)
                if idx(1)>15
                    current_tmean_flag(idx(1):end) = 1; %if the first flag is after day 15, screen from the 1st flag day to end of WY
                else
                    current_tmean_flag(1:end) = 1;%if 1st flag occurs before day 15, screen total WY
                end
            end
            idx=find(~isnan(current_tmax_flag));
            if ~isempty(idx)
                if idx(1)>15
                    current_tmax_flag(idx(1):end) = 1;%if the first flag is after day 15, screen from the 1st flag day to end of WY
                else
                    current_tmax_flag(1:end) = 1;%if 1st flag occurs before day 15, screen total WY
                end
            end
            idx=find(~isnan(current_tmin_flag));
            if ~isempty(idx)
                if idx(1)>15
                    current_tmin_flag(idx(1):end) = 1;%if the first flag is after day 15, screen from the 1st flag day to end of WY
                else
                    current_tmin_flag(1:end) = 1; %if 1st flag occurs before day 15, screen total WY
                end
            end

            idx_bad_SWE = find(~isnan(current_swe_flag));
            current_swe(idx_bad_SWE) = NaN;
            idx_bad_prcp = find(~isnan(current_prcp_flag));
            current_prcp(idx_bad_prcp) = NaN;
            idx_bad_tmean = find(~isnan(current_tmean_flag));
            current_tmean(idx_bad_tmean) = NaN;
            idx_bad_tmax = find(~isnan(current_tmax_flag));
            current_tmax(idx_bad_tmax) = NaN;
            idx_bad_tmin = find(~isnan(current_tmin_flag));
            current_tmin(idx_bad_tmin) = NaN;
        end
        %apply screened data across total WY if NaN in first 15 days, if
        %not, then apply NaNs from 1st day of invalid data:
        idx=find(~isnan(current_swe_flag));
        if ~isempty(idx) %if there are flags...
            if idx(1)>15 %if the first flag is after day 15, screen from the 1st flag day to end of WY
                current_swe_flag(idx(1):end) = 1;
            else
                current_swe_flag(1:end) = 1; %if 1st flag occurs before day 15, screen total WY
            end
        end
        idx=find(~isnan(current_prcp_flag));
        if ~isempty(idx)
            if idx(1)>15
                current_prcp_flag(idx(1):end) = 1; %if the first flag is after day 15, screen from the 1st flag day to end of WY
            else
                current_prcp_flag(1:end) = 1; %if 1st flag occurs before day 15, screen total WY
            end
        end
        idx=find(~isnan(current_tmean_flag));
        if ~isempty(idx)
            if idx(1)>15
                current_tmean_flag(idx(1):end) = 1; %if the first flag is after day 15, screen from the 1st flag day to end of WY
            else
                current_tmean_flag(1:end) = 1;%if 1st flag occurs before day 15, screen total WY
            end
        end
        idx=find(~isnan(current_tmax_flag));
        if ~isempty(idx)
            if idx(1)>15
                current_tmax_flag(idx(1):end) = 1;%if the first flag is after day 15, screen from the 1st flag day to end of WY
            else
                current_tmax_flag(1:end) = 1;%if 1st flag occurs before day 15, screen total WY
            end
        end
        idx=find(~isnan(current_tmin_flag));
        if ~isempty(idx)
            if idx(1)>15
                current_tmin_flag(idx(1):end) = 1;%if the first flag is after day 15, screen from the 1st flag day to end of WY
            else
                current_tmin_flag(1:end) = 1; %if 1st flag occurs before day 15, screen total WY
            end
        end

        idx_bad_SWE = find(~isnan(current_swe_flag));
        current_swe(idx_bad_SWE) = NaN;
        idx_bad_prcp = find(~isnan(current_prcp_flag));
        current_prcp(idx_bad_prcp) = NaN;
        idx_bad_tmean = find(~isnan(current_tmean_flag));
        current_tmean(idx_bad_tmean) = NaN;
        idx_bad_tmax = find(~isnan(current_tmax_flag));
        current_tmax(idx_bad_tmax) = NaN;
        idx_bad_tmin = find(~isnan(current_tmin_flag));
        current_tmin(idx_bad_tmin) = NaN;

        prcp_QC_daily(si,:,yi) = current_prcp;
        swe_QC(si,:,yi) = current_swe;
        tmean_QC(si,:,yi) = current_tmean;
        tmax_QC(si,:,yi) = current_tmax;
        tmin_QC(si,:,yi) = current_tmin;
        if ~isnan(max(current_tmax)) && ~isnan(min(current_tmin))
            assert(max(current_tmax)<314 & min(current_tmin)>231,'temp are out of set bounds')
        end
    end
end

%% apply bias correction using Sun et al. (2019) / Livneh et al. (2014)
%https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2014WR015442
%count # of times the BC step is used:
count=0;
for si=1:S(1)
    for yi=1:S(3)
        current_prcp = squeeze(prcp_QC_daily(si,:,yi));
        current_swe = squeeze(swe_QC(si,:,yi));
        current_tmean = squeeze(tmean_QC(si,:,yi));
        current_tmax = squeeze(tmax_QC(si,:,yi));
        current_tmin = squeeze(tmin_QC(si,:,yi));

        %apply bias correction for air temperature from Sun et al. (2019) - correction formula is for [C]:
        current_tmean = 1.03.*(current_tmean-273.15) - 0.9; %B.C. in C
        current_tmax = 1.03.*(current_tmax-273.15) - 0.9; %B.C. in C
        current_tmin = 1.03.*(current_tmin-273.15) - 0.9; %B.C. in C

        current_tmean = current_tmean+273.15; %convert back to K
        current_tmax = current_tmax+273.15; %convert back to K
        current_tmin = current_tmin+273.15; %convert back to K

        %apply bias correction for PRCP undercatch from Sun et al. (2019):
        window_length = 7;%days from Livneh et al. (2014)
        start_idx=0;
        end_idx=0;
        while end_idx < length(current_prcp)-window_length
            
            start_idx=end_idx+1;
            end_idx=start_idx+window_length-1;
            assert(start_idx<length(current_prcp)-(window_length),'out of bounds')

            prcp_tmp = current_prcp(start_idx:end_idx); 
            swe_tmp = current_swe(start_idx:(end_idx+1));%b/c prcp on day n is realized as SWE_N+1 - SWE_N
            delta_SWE = swe_tmp(end) - swe_tmp(1);
            total_prcp = cumsum(prcp_tmp);
            total_prcp = total_prcp(end)-total_prcp(1);

            if yi==16 && si==225 && start_idx>=51 && end_idx<=51
                STORE_EXAMPLE_swe_tmp = swe_tmp;
                STORE_EXAMPLE_PRCP_TMP = prcp_tmp;
            end
            if delta_SWE>total_prcp
                for di=1:length(prcp_tmp)
                    day_prcp = prcp_tmp(di);
                    day_delta_SWE = swe_tmp(di+1)-swe_tmp(di);
                    if day_delta_SWE>=0
                        prcp_tmp(di) = day_delta_SWE;
                    else
                        prcp_tmp(di) = day_prcp;
                    end
                end
                total_prcp_new = sum(prcp_tmp);
                assert(total_prcp_new - delta_SWE > -0.00001,'prcp is less than delta SWE after bias correction')
                current_prcp(start_idx:end_idx) = prcp_tmp;
            elseif total_prcp>=delta_SWE
                for di=1:length(prcp_tmp)
                    day_prcp = prcp_tmp(di);
                    day_delta_SWE = swe_tmp(di+1)-swe_tmp(di);
                    if day_prcp>0
                        prcp_tmp(di) = max([day_delta_SWE,day_prcp]);
                    end
                end
                total_prcp_new = sum(prcp_tmp);
                assert(total_prcp_new - delta_SWE > -0.00001,'prcp is less than delta SWE after bias correction')
                current_prcp(start_idx:end_idx) = prcp_tmp;
            end
        end
        %store data:
        [2,si,yi]
        prcp_QC_daily(si,:,yi) = current_prcp;
        swe_QC(si,:,yi) = current_swe;
        tmean_QC(si,:,yi) = current_tmean;
        tmax_QC(si,:,yi) = current_tmax;
        tmin_QC(si,:,yi) = current_tmin;
        if ~isnan(max(current_tmax)) && ~isnan(min(current_tmin))
            assert(max(current_tmax)<314 & min(current_tmin)>231,'temp are out of set bounds')
        end
    end
end

%% QC snow depth:
% %insert nan values where swe obs are larger than snowh:
idx_nan = find(swe_QC > snowh_QC);
snowh_QC(idx_nan) = NaN;

%  snowh = 0 when swe = 0
idx= find(swe_QC == 0);
snowh_QC(idx) = 0;

%screen very high or low snowh values:
idx = find(snowh_QC<0 | snowh_QC>50000);
snowh_QC(idx) = NaN;

for s=1:length(site_ids)
    for y=1:length(wys)
        current_snowh = squeeze(snowh_QC(s,:,y));
        current_swe = squeeze(swe_QC(s,:,y));

        delta_snowh = current_snowh(2:end) - current_snowh(1:end-1);
        delta_swe = current_swe(2:end) - current_swe(1:end-1);

        %if swe decreases, so should snowh
        idx=find(delta_swe < 0 & delta_snowh>=0);
        snowh_QC(s,idx+1,y) = NaN;

        %screen 1m spike in a day
        idx=find(delta_snowh>1000);
        if length(idx) > 0
            spike_val = snowh_QC(s,idx+1,y);
            snowh_QC(s,idx+1,y) = NaN;

            %if spike never comes down - remove all data affected by spike:
            new_current_snowh = squeeze(snowh_QC(s,:,y));

            if (max(spike_val)*1.05 > max(current_snowh) && max(spike_val)*0.95 < max(current_snowh))
                idx = find(new_current_snowh > 0.9.*max(spike_val) & new_current_snowh<1.1.*max(spike_val));
                if length(idx) > 0
                    snowh_QC(s,idx,y) = NaN;
                end
            end
        end
        %second, ensure peak snowh is higher than peak swe, if not toss out the entire year of data
        if max(current_snowh) < max(current_swe)
            snowh_QC(s,:,y) = NaN;
        end
        %finally, if snowh decreases swe must decrease by at leat 1% of snowh
        idx = find(delta_snowh <= -100 & delta_swe> delta_snowh/100);
        snowh_QC(s,idx+1,y) = NaN;
    end
end

%finally before writing - perform an outlier analysis for SWE vs snow depth
%(a) compute linear relationship for swe vs snowh:
SWE_vec = swe_QC(:);
snowh_QC_vec = snowh_QC(:);
idx=find(isnan(SWE_vec) | isnan(snowh_QC_vec));
SWE_vec(idx)= [];
snowh_QC_vec(idx)=[];
p = polyfit(SWE_vec,snowh_QC_vec,1);
lm = @(swe) p(1).*swe + p(2);

SWE_vec = swe_QC(:);
yhat = lm(SWE_vec);

snowh_QC_vec = snowh_QC(:);
residule = abs(yhat - snowh_QC_vec);
idx = find(residule>2000);
snowh_QC(idx) = NaN;

%% write data
if exist(new_filename,'file')>0
    CMD_rm = ['rm ',new_filename];
    system(CMD_rm);
end
cmd_CP = ['cp ',SNOTEL_data_filename_PNNL,' ',new_filename];
system(cmd_CP);

%remove variables to overwrite:
variables_to_write = {'PRCP_pnnl','PRCP_raw','SWE_pnnl','SWE_raw','T2_raw','T2_pnnl','T2max_raw','T2max_pnnl','T2min_pnnl','T2min_raw','day','month','year','nyear'};
nvars = length(variables_to_write);
for i=1:nvars
    var_name = variables_to_write{i};
    %remove from outfile
    CMD_ncks = ['/opt/local/bin/ncks -C -O -x -v ',var_name,' ',new_filename,' ',new_filename];
    system(CMD_ncks);
end

%write new variables:
nccreate(new_filename,'WY','Datatype','double',...
    'Dimensions', {'nyear',length(wys)});
ncwrite(new_filename,'WY', wys);

nccreate(new_filename,'snowh_raw','Datatype','double',...
    'Dimensions', {'nsite',length(site_ids),'nday',366,'nyear',length(wys)});
ncwrite(new_filename,'snowh_raw', snowh_raw);
ncwriteatt(new_filename,'snowh_raw','units','mm');

nccreate(new_filename,'snowh_QC','Datatype','double',...
    'Dimensions', {'nsite',length(site_ids),'nday',366,'nyear',length(wys)});
ncwrite(new_filename,'snowh_QC', snowh_QC);
ncwriteatt(new_filename,'snowh_QC','units','mm');

nccreate(new_filename,'swe_raw','Datatype','double',...
    'Dimensions', {'nsite',length(site_ids),'nday',366,'nyear',length(wys)});
ncwrite(new_filename,'swe_raw', swe_raw);
ncwriteatt(new_filename,'swe_raw','units','mm');

nccreate(new_filename,'swe_QC','Datatype','double',...
    'Dimensions', {'nsite',length(site_ids),'nday',366,'nyear',length(wys)});
ncwrite(new_filename,'swe_QC', swe_QC);
ncwriteatt(new_filename,'swe_QC','units','mm');

nccreate(new_filename,'T2mean_raw','Datatype','double',...
    'Dimensions', {'nsite',length(site_ids),'nday',366,'nyear',length(wys)});
ncwrite(new_filename,'T2mean_raw', tmean_raw);
ncwriteatt(new_filename,'T2mean_raw','units','K');

nccreate(new_filename,'T2mean_QCBC','Datatype','double',...
    'Dimensions', {'nsite',length(site_ids),'nday',366,'nyear',length(wys)});
ncwrite(new_filename,'T2mean_QCBC', tmean_QC);
ncwriteatt(new_filename,'T2mean_QCBC','units','K');

nccreate(new_filename,'T2max_raw','Datatype','double',...
    'Dimensions', {'nsite',length(site_ids),'nday',366,'nyear',length(wys)});
ncwrite(new_filename,'T2max_raw', tmax_raw);
ncwriteatt(new_filename,'T2max_raw','units','K');

nccreate(new_filename,'T2max_QCBC','Datatype','double',...
    'Dimensions', {'nsite',length(site_ids),'nday',366,'nyear',length(wys)});
ncwrite(new_filename,'T2max_QCBC', tmax_QC);
ncwriteatt(new_filename,'T2max_QCBC','units','K');

nccreate(new_filename,'T2min_raw','Datatype','double',...
    'Dimensions', {'nsite',length(site_ids),'nday',366,'nyear',length(wys)});
ncwrite(new_filename,'T2min_raw', tmin_raw);
ncwriteatt(new_filename,'T2min_raw','units','K');

nccreate(new_filename,'T2min_QCBC','Datatype','double',...
    'Dimensions', {'nsite',length(site_ids),'nday',366,'nyear',length(wys)});
ncwrite(new_filename,'T2min_QCBC', tmin_QC);
ncwriteatt(new_filename,'T2min_QCBC','units','K');

nccreate(new_filename,'PRCP_raw','Datatype','double',...
    'Dimensions', {'nsite',length(site_ids),'nday',366,'nyear',length(wys)});
ncwrite(new_filename,'PRCP_raw', prcp_raw_daily);
ncwriteatt(new_filename,'PRCP_raw','units','mm/day');

nccreate(new_filename,'prcp_QCBC','Datatype','double',...
    'Dimensions', {'nsite',length(site_ids),'nday',366,'nyear',length(wys)});
ncwrite(new_filename,'prcp_QCBC', prcp_QC_daily);
ncwriteatt(new_filename,'prcp_QCBC','units','mm/day');



