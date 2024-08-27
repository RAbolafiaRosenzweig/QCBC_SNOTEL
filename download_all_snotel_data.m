clc;clear all;close all;

%% adjustable input data

%get SITE IDs and corresponding state IDs for each site you want data for
SNOTEL_wus_data_filename = '/Volumes/Pruina_External_Elements/SnowData/SNOTEL/SNOTEL_all_obs_PNNLbcqc_WY1979-2021.nc';
site_ids = double(ncread(SNOTEL_wus_data_filename,'site_id'));
states = ncread(SNOTEL_wus_data_filename,'state');

%define outdir:
outdir = '/Volumes/Pruina_External_Elements/WPO_SnowModeling/Observations/SNOTEL/raw/1979_2023_Daily_All_Data/';
if exist(outdir,'dir') ==0
    CMD_mkdir = ['mkdir -p ',outdir];
    system(CMD_mkdir);
end

%define WYs to get data for
start_WY = 1979;
end_WY = 2023;

%% do not change  --get data--
for si = 1:length(site_ids)
    site_ID = site_ids(si);
    state = states(si,:);
    cmd_get = sprintf('curl -k -L "https://wcc.sc.egov.usda.gov/reportGenerator/view_csv/customMultiTimeSeriesGroupByStationReport/daily/start_of_period/%d:%s:SNTL/%d-10-01,%d-09-30/stationId,name,WTEQ::value,PREC::value,TOBS::value,TAVG::value,TMAX::value,TMIN::value,SNWD::value,SNWD::qcFlag,SNWD::qaFlag?fitToScreen=false" -o %s/STA%d-%d-%d.csv',site_ID,state,start_WY,end_WY,outdir,site_ID,start_WY,end_WY);
    system(cmd_get)
    %%%%%other values can include: WTEQ::value,PREC::value,TOBS::value,TAVG::value,TMAX::value,TMIN::value
end

