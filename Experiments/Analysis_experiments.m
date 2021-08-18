% Analysis experiments
% To run by section

%% Table 3 Simulations experiments
clear;close all;clc;
load('../simulations_baseline.mat')
disp('       EY        EK       EB        SY        SI        RGY        RIY       RBY')
disp([mean(EY) mean(EK) mean(EB) mean(SDY) mean(SDI) mean(cor.GY) mean(cor.IY) mean(cor.BY)])
load('Simulations/simulations_fc_86.mat')
disp('       EY        EK       EB        SY        SI        RGY        RIY       RBY')
disp([mean(EY) mean(EK) mean(EB) mean(SDY) mean(SDI) mean(cor.GY) mean(cor.IY) mean(cor.BY)])
load('Simulations/simulations_fc_99.mat')
disp('       EY        EK       EB        SY        SI        RGY        RIY       RBY')
disp([mean(EY) mean(EK) mean(EB) mean(SDY) mean(SDI) mean(cor.GY) mean(cor.IY) mean(cor.BY)])


%% Table 4 Simulations Golden Rule
clear;close all;clc;
load('../simulations_baseline.mat')
baseline.y = mean(EY);
baseline.k = mean(EK);
baseline.l = mean(EL);
baseline.b = mean(EB); 
baseline.b_y = mean(EB)/mean(EY);
baseline.g = mean(EG);
baseline.g_y = mean(EG)/mean(EY);
baseline.c = mean(EC);
baseline.u = mean(EU);

baseline.d = mean(ED);
baseline.spread = mean(ES);

baseline.sdy = mean(SDY);
baseline.sdi = mean(SDI);
baseline.sdb = mean(SDB);
baseline.sdg = mean(SDG);
baseline.sdc = mean(SDC);
baseline.sds = mean(SDS);

baseline.gy = mean(cor.GY);
baseline.by = mean(cor.BY);
baseline.iy = mean(cor.IY);

clearvars -except baseline
load('Simulations/simulations_gr_1.mat')

gr_1.y = mean(EY);
gr_1.k = mean(EK);
gr_1.l = mean(EL);
gr_1.b = mean(EB);
gr_1.b_y = mean(EB)/mean(EY);
gr_1.g = mean(EG);
gr_1.g_y = mean(EG)/mean(EY);
gr_1.c = mean(EC);
gr_1.u = mean(EU);

gr_1.d = mean(ED);
gr_1.spread = mean(ES);

gr_1.sdy = mean(SDY);
gr_1.sdi = mean(SDI);
gr_1.sdb = mean(SDB);
gr_1.sdg = mean(SDG);
gr_1.sdc = mean(SDC);
gr_1.sds = mean(SDS);

gr_1.gy = mean(cor.GY);
gr_1.by = mean(cor.BY);
gr_1.iy = mean(cor.IY);

%%%%%%%%%%%%%%%%%%%%%%
gr_1per.y = mean(EY)/baseline.y;
gr_1per.k = mean(EK)/baseline.k;
gr_1per.l = mean(EL)/baseline.l;
gr_1per.b = mean(EB)/baseline.b;
gr_1per.b_y = gr_1.b_y/baseline.b_y;
gr_1per.g = mean(EG)/baseline.g;
gr_1per.g_y = 1-gr_1.g_y/baseline.g_y;
gr_1per.c = mean(EC)/baseline.c;
gr_1per.u = 1 + 1-mean(EU)/baseline.u;

gr_1per.d = mean(ED)/baseline.d;
gr_1per.spread = mean(ES)/baseline.spread;

gr_1per.sdy = 1-mean(SDY)/baseline.sdy;
gr_1per.sdi = 1-mean(SDI)/baseline.sdi;
gr_1per.sdb = 1-mean(SDB)/baseline.sdb;
gr_1per.sdg = 1-mean(SDG)/baseline.sdg;
gr_1per.sdc = 1-mean(SDC)/baseline.sdc;
gr_1per.sds = 1-mean(SDS)/baseline.sds;

gr_1per.gy = mean(cor.GY)/baseline.gy;
gr_1per.by = 1-mean(cor.BY)/baseline.by;
gr_1per.iy = 1-mean(cor.IY)/baseline.iy;

clearvars -except baseline gr_1 gr_1per 
load('Simulations/simulations_gr_0.5.mat')

gr_05.y = mean(EY);
gr_05.k = mean(EK);
gr_05.l = mean(EL);
gr_05.b = mean(EB);
gr_05.b_y = mean(EB)/mean(EY);
gr_05.g = mean(EG);
gr_05.g_y = mean(EG)/mean(EY);
gr_05.c = mean(EC);
gr_05.u = mean(EU);

gr_05.d = mean(ED);
gr_05.spread = mean(ES);

gr_05.sdy = mean(SDY);
gr_05.sdi = mean(SDI);
gr_05.sdb = mean(SDB);
gr_05.sdg = mean(SDG);
gr_05.sdc = mean(SDC);
gr_05.sds = mean(SDS);

gr_05.gy = mean(cor.GY);
gr_05.by = mean(cor.BY);
gr_05.iy = mean(cor.IY);

%%%%%%%%%%%%%%%%%%%%%%
gr_05per.y = mean(EY)/baseline.y;
gr_05per.k = mean(EK)/baseline.k;
gr_05per.l = mean(EL)/baseline.l;
gr_05per.b = mean(EB)/baseline.b;
gr_05per.b_y = gr_05.b_y/baseline.b_y;
gr_05per.g = mean(EG)/baseline.g;
gr_05per.g_y = 1-gr_05.g_y/baseline.g_y;
gr_05per.c = mean(EC)/baseline.c;
gr_05per.u = 1 + 1-mean(EU)/baseline.u;

gr_05per.d = mean(ED)/baseline.d;
gr_05per.spread = mean(ES)/baseline.spread;

gr_05per.sdy = 1-mean(SDY)/baseline.sdy;
gr_05per.sdi = 1-mean(SDI)/baseline.sdi;
gr_05per.sdb = 1-mean(SDB)/baseline.sdb;
gr_05per.sdg = 1-mean(SDG)/baseline.sdg;
gr_05per.sdc = 1-mean(SDC)/baseline.sdc;
gr_05per.sds = 1-mean(SDS)/baseline.sds;

gr_05per.gy = mean(cor.GY)/baseline.gy;
gr_05per.by = 1-mean(cor.BY)/baseline.by;
gr_05per.iy = 1-mean(cor.IY)/baseline.iy;

clearvars -except baseline gr_1 gr_1per gr_05 gr_05per 
load('Simulations/simulations_gr_0.01.mat')

gr_0.y = mean(EY);
gr_0.k = mean(EK);
gr_0.l = mean(EL);
gr_0.b = mean(EB);
gr_0.b_y = mean(EB)/mean(EY);
gr_0.g = mean(EG);
gr_0.g_y = mean(EG)/mean(EY);
gr_0.c = mean(EC);
gr_0.u = mean(EU);

gr_0.d = mean(ED);
gr_0.spread = mean(ES);

gr_0.sdy = mean(SDY);
gr_0.sdi = mean(SDI);
gr_0.sdb = mean(SDB);
gr_0.sdg = mean(SDG);
gr_0.sdc = mean(SDC);
gr_0.sds = mean(SDS);

gr_0.gy = mean(cor.GY);
gr_0.by = mean(cor.BY);
gr_0.iy = mean(cor.IY);

%%%%%%%%%%%%%%%%%%%%%%
gr_0per.y = mean(EY)/baseline.y;
gr_0per.k = mean(EK)/baseline.k;
gr_0per.l = mean(EL)/baseline.l;
gr_0per.b = mean(EB)/baseline.b;
gr_0per.b_y = gr_0.b_y/baseline.b_y;
gr_0per.g = mean(EG)/baseline.g;
gr_0per.g_y = 1-gr_0.g_y/baseline.g_y;
gr_0per.c = mean(EC)/baseline.c;
gr_0per.u = 1 + 1-mean(EU)/baseline.u;

gr_0per.d = mean(ED)/baseline.d;
gr_0per.spread = mean(ES)/baseline.spread;

gr_0per.sdy = 1-mean(SDY)/baseline.sdy;
gr_0per.sdi = 1-mean(SDI)/baseline.sdi;
gr_0per.sdb = 1-mean(SDB)/baseline.sdb;
gr_0per.sdg = 1-mean(SDG)/baseline.sdg;
gr_0per.sdc = 1-mean(SDC)/baseline.sdc;
gr_0per.sds = 1-mean(SDS)/baseline.sds;

gr_0per.gy = mean(cor.GY)/baseline.gy;
gr_0per.by = 1-mean(cor.BY)/baseline.by;
gr_0per.iy = 1-mean(cor.IY)/baseline.iy;

clearvars -except baseline gr_1 gr_1per gr_05 gr_05per gr_0 gr_0per 

disp(baseline)
disp(gr_1)
disp(gr_05)
disp(gr_0)

%% Table 5 debt ceiling
load('Simulations/simulations_dc_0.3.mat')

dc_3.y = mean(EY);
dc_3.k = mean(EK);
dc_3.l = mean(EL);
dc_3.b = mean(EB);
dc_3.b_y = mean(EB)/mean(EY);
dc_3.g = mean(EG);
dc_3.g_y = mean(EG)/mean(EY);
dc_3.c = mean(EC);
dc_3.u = mean(EU);

dc_3.d = mean(ED);
dc_3.spread = mean(ES);

dc_3.sdy = mean(SDY);
dc_3.sdi = mean(SDI);
dc_3.sdb = mean(SDB);
dc_3.sdg = mean(SDG);
dc_3.sdc = mean(SDC);
dc_3.sds = mean(SDS);

dc_3.gy = mean(cor.GY);
dc_3.by = mean(cor.BY);
dc_3.iy = mean(cor.IY);

%%%%%%%%%%%%%%%%%%%%%%
dc_3per.y = mean(EY)/baseline.y;
dc_3per.k = mean(EK)/baseline.k;
dc_3per.l = mean(EL)/baseline.l;
dc_3per.b = mean(EB)/baseline.b;
dc_3per.b_y = dc_3.b_y/baseline.b_y;
dc_3per.g = mean(EG)/baseline.g;
dc_3per.g_y = 1-dc_3.g_y/baseline.g_y;
dc_3per.c = mean(EC)/baseline.c;
dc_3per.u = 1 + 1-mean(EU)/baseline.u;

dc_3per.d = mean(ED)/baseline.d;
dc_3per.spread = mean(ES)/baseline.spread;

dc_3per.sdy = 1-mean(SDY)/baseline.sdy;
dc_3per.sdi = 1-mean(SDI)/baseline.sdi;
dc_3per.sdb = 1-mean(SDB)/baseline.sdb;
dc_3per.sdg = 1-mean(SDG)/baseline.sdg;
dc_3per.sdc = 1-mean(SDC)/baseline.sdc;
dc_3per.sds = 1-mean(SDS)/baseline.sds;

dc_3per.gy = mean(cor.GY)/baseline.gy;
dc_3per.by = 1-mean(cor.BY)/baseline.by;
dc_3per.iy = 1-mean(cor.IY)/baseline.iy;


clearvars -except baseline gr_1 gr_1per gr_05 gr_05per gr_0 gr_0per ...
    dc_3 dc_3per
load('Simulations/simulations_dc_0.4.mat')

dc_4.y = mean(EY);
dc_4.k = mean(EK);
dc_4.l = mean(EL);
dc_4.b = mean(EB);
dc_4.b_y = mean(EB)/mean(EY);
dc_4.g = mean(EG);
dc_4.g_y = mean(EG)/mean(EY);
dc_4.c = mean(EC);
dc_4.u = mean(EU);

dc_4.d = mean(ED);
dc_4.spread = mean(ES);

dc_4.sdy = mean(SDY);
dc_4.sdi = mean(SDI);
dc_4.sdb = mean(SDB);
dc_4.sdg = mean(SDG);
dc_4.sdc = mean(SDC);
dc_4.sds = mean(SDS);

dc_4.gy = mean(cor.GY);
dc_4.by = mean(cor.BY);
dc_4.iy = mean(cor.IY);

%%%%%%%%%%%%%%%%%%%%%%
dc_4per.y = mean(EY)/baseline.y;
dc_4per.k = mean(EK)/baseline.k;
dc_4per.l = mean(EL)/baseline.l;
dc_4per.b = mean(EB)/baseline.b;
dc_4per.b_y = dc_4.b_y/baseline.b_y;
dc_4per.g = mean(EG)/baseline.g;
dc_4per.g_y = 1-dc_4.g_y/baseline.g_y;
dc_4per.c = mean(EC)/baseline.c;
dc_4per.u = 1 + 1-mean(EU)/baseline.u;

dc_4per.d = mean(ED)/baseline.d;
dc_4per.spread = mean(ES)/baseline.spread;

dc_4per.sdy = 1-mean(SDY)/baseline.sdy;
dc_4per.sdi = 1-mean(SDI)/baseline.sdi;
dc_4per.sdb = 1-mean(SDB)/baseline.sdb;
dc_4per.sdg = 1-mean(SDG)/baseline.sdg;
dc_4per.sdc = 1-mean(SDC)/baseline.sdc;
dc_4per.sds = 1-mean(SDS)/baseline.sds;

dc_4per.gy = mean(cor.GY)/baseline.gy;
dc_4per.by = 1-mean(cor.BY)/baseline.by;
dc_4per.iy = 1-mean(cor.IY)/baseline.iy;


clearvars -except baseline gr_1 gr_1per gr_05 gr_05per gr_0 gr_0per ...
    dc_3 dc_3per dc_4 dc_4per
load('Simulations/simulations_dc_0.5.mat')
dc_5.y = mean(EY);
dc_5.k = mean(EK);
dc_5.l = mean(EL);
dc_5.b = mean(EB);
dc_5.b_y = mean(EB)/mean(EY);
dc_5.g = mean(EG);
dc_5.g_y = mean(EG)/mean(EY);
dc_5.c = mean(EC);
dc_5.u = mean(EU);

dc_5.d = mean(ED);
dc_5.spread = mean(ES);

dc_5.sdy = mean(SDY);
dc_5.sdi = mean(SDI);
dc_5.sdb = mean(SDB);
dc_5.sdg = mean(SDG);
dc_5.sdc = mean(SDC);
dc_5.sds = mean(SDS);

dc_5.gy = mean(cor.GY);
dc_5.by = mean(cor.BY);
dc_5.iy = mean(cor.IY);

%%%%%%%%%%%%%%%%%%%%%%
dc_5per.y = mean(EY)/baseline.y;
dc_5per.k = mean(EK)/baseline.k;
dc_5per.l = mean(EL)/baseline.l;
dc_5per.b = mean(EB)/baseline.b;
dc_5per.b_y = dc_5.b_y/baseline.b_y;
dc_5per.g = mean(EG)/baseline.g;
dc_5per.g_y = 1-dc_5.g_y/baseline.g_y;
dc_5per.c = mean(EC)/baseline.c;
dc_5per.u = 1 + 1-mean(EU)/baseline.u;

dc_5per.d = mean(ED)/baseline.d;
dc_5per.spread = mean(ES)/baseline.spread;

dc_5per.sdy = 1-mean(SDY)/baseline.sdy;
dc_5per.sdi = 1-mean(SDI)/baseline.sdi;
dc_5per.sdb = 1-mean(SDB)/baseline.sdb;
dc_5per.sdg = 1-mean(SDG)/baseline.sdg;
dc_5per.sdc = 1-mean(SDC)/baseline.sdc;
dc_5per.sds = 1-mean(SDS)/baseline.sds;

dc_5per.gy = mean(cor.GY)/baseline.gy;
dc_5per.by = 1-mean(cor.BY)/baseline.by;
dc_5per.iy = 1-mean(cor.IY)/baseline.iy;


clearvars -except baseline gr_1 gr_1per gr_05 gr_05per gr_0 gr_0per ...
    dc_3 dc_3per dc_4 dc_4per dc_5 dc_5per 

disp(baseline)
disp(dc_3)
disp(dc_4)
disp(dc_5)