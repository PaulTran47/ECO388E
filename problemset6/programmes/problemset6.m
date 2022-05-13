%==========================================================================
%% 2 percentage signs represent sections of code;
% 1 percentage sign represents comments for code or commented out code;

% Answers to question parts that don't involve code can be found at the
% bottom of the programme, in the section ``Questions asked in problemset x
% that don't involve code".

% Text answers to question parts that involve code will be between the
% sub-section label:
%=======
% ANSWER
%=======
% Answer here
%===========
% END ANSWER
%===========

% Comments that are important will be between the sub-section label:
%=====
% NOTE
%=====
% Important note here
%=========
% END NOTE
%=========
% ECO388E Problem Set 6
% Paul Le Tran, plt377
% 12 May, 2022
%==========================================================================

%==========================================================================
%% Setting up workspace
clear all;
close all;
clc;

home_dir = 'path\to\programmes';
data_dir = 'path\to\data';

% Setting text interpreter to latex
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaulttextinterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');

% Setting random seed for reproductibility
rng(47);

cd(home_dir);
%==========================================================================

%==========================================================================
%% Importing data
cd(data_dir);
% Importing training data
opts = delimitedTextImportOptions("NumVariables", 529);
% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";
% Specify column names and types
opts.VariableNames = ["WAP001", "WAP002", "WAP003", "WAP004", "WAP005", "WAP006", "WAP007", "WAP008", "WAP009", "WAP010", "WAP011", "WAP012", "WAP013", "WAP014", "WAP015", "WAP016", "WAP017", "WAP018", "WAP019", "WAP020", "WAP021", "WAP022", "WAP023", "WAP024", "WAP025", "WAP026", "WAP027", "WAP028", "WAP029", "WAP030", "WAP031", "WAP032", "WAP033", "WAP034", "WAP035", "WAP036", "WAP037", "WAP038", "WAP039", "WAP040", "WAP041", "WAP042", "WAP043", "WAP044", "WAP045", "WAP046", "WAP047", "WAP048", "WAP049", "WAP050", "WAP051", "WAP052", "WAP053", "WAP054", "WAP055", "WAP056", "WAP057", "WAP058", "WAP059", "WAP060", "WAP061", "WAP062", "WAP063", "WAP064", "WAP065", "WAP066", "WAP067", "WAP068", "WAP069", "WAP070", "WAP071", "WAP072", "WAP073", "WAP074", "WAP075", "WAP076", "WAP077", "WAP078", "WAP079", "WAP080", "WAP081", "WAP082", "WAP083", "WAP084", "WAP085", "WAP086", "WAP087", "WAP088", "WAP089", "WAP090", "WAP091", "WAP092", "WAP093", "WAP094", "WAP095", "WAP096", "WAP097", "WAP098", "WAP099", "WAP100", "WAP101", "WAP102", "WAP103", "WAP104", "WAP105", "WAP106", "WAP107", "WAP108", "WAP109", "WAP110", "WAP111", "WAP112", "WAP113", "WAP114", "WAP115", "WAP116", "WAP117", "WAP118", "WAP119", "WAP120", "WAP121", "WAP122", "WAP123", "WAP124", "WAP125", "WAP126", "WAP127", "WAP128", "WAP129", "WAP130", "WAP131", "WAP132", "WAP133", "WAP134", "WAP135", "WAP136", "WAP137", "WAP138", "WAP139", "WAP140", "WAP141", "WAP142", "WAP143", "WAP144", "WAP145", "WAP146", "WAP147", "WAP148", "WAP149", "WAP150", "WAP151", "WAP152", "WAP153", "WAP154", "WAP155", "WAP156", "WAP157", "WAP158", "WAP159", "WAP160", "WAP161", "WAP162", "WAP163", "WAP164", "WAP165", "WAP166", "WAP167", "WAP168", "WAP169", "WAP170", "WAP171", "WAP172", "WAP173", "WAP174", "WAP175", "WAP176", "WAP177", "WAP178", "WAP179", "WAP180", "WAP181", "WAP182", "WAP183", "WAP184", "WAP185", "WAP186", "WAP187", "WAP188", "WAP189", "WAP190", "WAP191", "WAP192", "WAP193", "WAP194", "WAP195", "WAP196", "WAP197", "WAP198", "WAP199", "WAP200", "WAP201", "WAP202", "WAP203", "WAP204", "WAP205", "WAP206", "WAP207", "WAP208", "WAP209", "WAP210", "WAP211", "WAP212", "WAP213", "WAP214", "WAP215", "WAP216", "WAP217", "WAP218", "WAP219", "WAP220", "WAP221", "WAP222", "WAP223", "WAP224", "WAP225", "WAP226", "WAP227", "WAP228", "WAP229", "WAP230", "WAP231", "WAP232", "WAP233", "WAP234", "WAP235", "WAP236", "WAP237", "WAP238", "WAP239", "WAP240", "WAP241", "WAP242", "WAP243", "WAP244", "WAP245", "WAP246", "WAP247", "WAP248", "WAP249", "WAP250", "WAP251", "WAP252", "WAP253", "WAP254", "WAP255", "WAP256", "WAP257", "WAP258", "WAP259", "WAP260", "WAP261", "WAP262", "WAP263", "WAP264", "WAP265", "WAP266", "WAP267", "WAP268", "WAP269", "WAP270", "WAP271", "WAP272", "WAP273", "WAP274", "WAP275", "WAP276", "WAP277", "WAP278", "WAP279", "WAP280", "WAP281", "WAP282", "WAP283", "WAP284", "WAP285", "WAP286", "WAP287", "WAP288", "WAP289", "WAP290", "WAP291", "WAP292", "WAP293", "WAP294", "WAP295", "WAP296", "WAP297", "WAP298", "WAP299", "WAP300", "WAP301", "WAP302", "WAP303", "WAP304", "WAP305", "WAP306", "WAP307", "WAP308", "WAP309", "WAP310", "WAP311", "WAP312", "WAP313", "WAP314", "WAP315", "WAP316", "WAP317", "WAP318", "WAP319", "WAP320", "WAP321", "WAP322", "WAP323", "WAP324", "WAP325", "WAP326", "WAP327", "WAP328", "WAP329", "WAP330", "WAP331", "WAP332", "WAP333", "WAP334", "WAP335", "WAP336", "WAP337", "WAP338", "WAP339", "WAP340", "WAP341", "WAP342", "WAP343", "WAP344", "WAP345", "WAP346", "WAP347", "WAP348", "WAP349", "WAP350", "WAP351", "WAP352", "WAP353", "WAP354", "WAP355", "WAP356", "WAP357", "WAP358", "WAP359", "WAP360", "WAP361", "WAP362", "WAP363", "WAP364", "WAP365", "WAP366", "WAP367", "WAP368", "WAP369", "WAP370", "WAP371", "WAP372", "WAP373", "WAP374", "WAP375", "WAP376", "WAP377", "WAP378", "WAP379", "WAP380", "WAP381", "WAP382", "WAP383", "WAP384", "WAP385", "WAP386", "WAP387", "WAP388", "WAP389", "WAP390", "WAP391", "WAP392", "WAP393", "WAP394", "WAP395", "WAP396", "WAP397", "WAP398", "WAP399", "WAP400", "WAP401", "WAP402", "WAP403", "WAP404", "WAP405", "WAP406", "WAP407", "WAP408", "WAP409", "WAP410", "WAP411", "WAP412", "WAP413", "WAP414", "WAP415", "WAP416", "WAP417", "WAP418", "WAP419", "WAP420", "WAP421", "WAP422", "WAP423", "WAP424", "WAP425", "WAP426", "WAP427", "WAP428", "WAP429", "WAP430", "WAP431", "WAP432", "WAP433", "WAP434", "WAP435", "WAP436", "WAP437", "WAP438", "WAP439", "WAP440", "WAP441", "WAP442", "WAP443", "WAP444", "WAP445", "WAP446", "WAP447", "WAP448", "WAP449", "WAP450", "WAP451", "WAP452", "WAP453", "WAP454", "WAP455", "WAP456", "WAP457", "WAP458", "WAP459", "WAP460", "WAP461", "WAP462", "WAP463", "WAP464", "WAP465", "WAP466", "WAP467", "WAP468", "WAP469", "WAP470", "WAP471", "WAP472", "WAP473", "WAP474", "WAP475", "WAP476", "WAP477", "WAP478", "WAP479", "WAP480", "WAP481", "WAP482", "WAP483", "WAP484", "WAP485", "WAP486", "WAP487", "WAP488", "WAP489", "WAP490", "WAP491", "WAP492", "WAP493", "WAP494", "WAP495", "WAP496", "WAP497", "WAP498", "WAP499", "WAP500", "WAP501", "WAP502", "WAP503", "WAP504", "WAP505", "WAP506", "WAP507", "WAP508", "WAP509", "WAP510", "WAP511", "WAP512", "WAP513", "WAP514", "WAP515", "WAP516", "WAP517", "WAP518", "WAP519", "WAP520", "LONGITUDE", "LATITUDE", "FLOOR", "BUILDINGID", "SPACEID", "RELATIVEPOSITION", "USERID", "PHONEID", "TIMESTAMP"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];
% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
% Import the data
train_table = readtable(append(data_dir, '\trainingData.csv'), opts);
clear opts;

% Importing validation data
opts = delimitedTextImportOptions("NumVariables", 529);
% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";
% Specify column names and types
opts.VariableNames = ["WAP001", "WAP002", "WAP003", "WAP004", "WAP005", "WAP006", "WAP007", "WAP008", "WAP009", "WAP010", "WAP011", "WAP012", "WAP013", "WAP014", "WAP015", "WAP016", "WAP017", "WAP018", "WAP019", "WAP020", "WAP021", "WAP022", "WAP023", "WAP024", "WAP025", "WAP026", "WAP027", "WAP028", "WAP029", "WAP030", "WAP031", "WAP032", "WAP033", "WAP034", "WAP035", "WAP036", "WAP037", "WAP038", "WAP039", "WAP040", "WAP041", "WAP042", "WAP043", "WAP044", "WAP045", "WAP046", "WAP047", "WAP048", "WAP049", "WAP050", "WAP051", "WAP052", "WAP053", "WAP054", "WAP055", "WAP056", "WAP057", "WAP058", "WAP059", "WAP060", "WAP061", "WAP062", "WAP063", "WAP064", "WAP065", "WAP066", "WAP067", "WAP068", "WAP069", "WAP070", "WAP071", "WAP072", "WAP073", "WAP074", "WAP075", "WAP076", "WAP077", "WAP078", "WAP079", "WAP080", "WAP081", "WAP082", "WAP083", "WAP084", "WAP085", "WAP086", "WAP087", "WAP088", "WAP089", "WAP090", "WAP091", "WAP092", "WAP093", "WAP094", "WAP095", "WAP096", "WAP097", "WAP098", "WAP099", "WAP100", "WAP101", "WAP102", "WAP103", "WAP104", "WAP105", "WAP106", "WAP107", "WAP108", "WAP109", "WAP110", "WAP111", "WAP112", "WAP113", "WAP114", "WAP115", "WAP116", "WAP117", "WAP118", "WAP119", "WAP120", "WAP121", "WAP122", "WAP123", "WAP124", "WAP125", "WAP126", "WAP127", "WAP128", "WAP129", "WAP130", "WAP131", "WAP132", "WAP133", "WAP134", "WAP135", "WAP136", "WAP137", "WAP138", "WAP139", "WAP140", "WAP141", "WAP142", "WAP143", "WAP144", "WAP145", "WAP146", "WAP147", "WAP148", "WAP149", "WAP150", "WAP151", "WAP152", "WAP153", "WAP154", "WAP155", "WAP156", "WAP157", "WAP158", "WAP159", "WAP160", "WAP161", "WAP162", "WAP163", "WAP164", "WAP165", "WAP166", "WAP167", "WAP168", "WAP169", "WAP170", "WAP171", "WAP172", "WAP173", "WAP174", "WAP175", "WAP176", "WAP177", "WAP178", "WAP179", "WAP180", "WAP181", "WAP182", "WAP183", "WAP184", "WAP185", "WAP186", "WAP187", "WAP188", "WAP189", "WAP190", "WAP191", "WAP192", "WAP193", "WAP194", "WAP195", "WAP196", "WAP197", "WAP198", "WAP199", "WAP200", "WAP201", "WAP202", "WAP203", "WAP204", "WAP205", "WAP206", "WAP207", "WAP208", "WAP209", "WAP210", "WAP211", "WAP212", "WAP213", "WAP214", "WAP215", "WAP216", "WAP217", "WAP218", "WAP219", "WAP220", "WAP221", "WAP222", "WAP223", "WAP224", "WAP225", "WAP226", "WAP227", "WAP228", "WAP229", "WAP230", "WAP231", "WAP232", "WAP233", "WAP234", "WAP235", "WAP236", "WAP237", "WAP238", "WAP239", "WAP240", "WAP241", "WAP242", "WAP243", "WAP244", "WAP245", "WAP246", "WAP247", "WAP248", "WAP249", "WAP250", "WAP251", "WAP252", "WAP253", "WAP254", "WAP255", "WAP256", "WAP257", "WAP258", "WAP259", "WAP260", "WAP261", "WAP262", "WAP263", "WAP264", "WAP265", "WAP266", "WAP267", "WAP268", "WAP269", "WAP270", "WAP271", "WAP272", "WAP273", "WAP274", "WAP275", "WAP276", "WAP277", "WAP278", "WAP279", "WAP280", "WAP281", "WAP282", "WAP283", "WAP284", "WAP285", "WAP286", "WAP287", "WAP288", "WAP289", "WAP290", "WAP291", "WAP292", "WAP293", "WAP294", "WAP295", "WAP296", "WAP297", "WAP298", "WAP299", "WAP300", "WAP301", "WAP302", "WAP303", "WAP304", "WAP305", "WAP306", "WAP307", "WAP308", "WAP309", "WAP310", "WAP311", "WAP312", "WAP313", "WAP314", "WAP315", "WAP316", "WAP317", "WAP318", "WAP319", "WAP320", "WAP321", "WAP322", "WAP323", "WAP324", "WAP325", "WAP326", "WAP327", "WAP328", "WAP329", "WAP330", "WAP331", "WAP332", "WAP333", "WAP334", "WAP335", "WAP336", "WAP337", "WAP338", "WAP339", "WAP340", "WAP341", "WAP342", "WAP343", "WAP344", "WAP345", "WAP346", "WAP347", "WAP348", "WAP349", "WAP350", "WAP351", "WAP352", "WAP353", "WAP354", "WAP355", "WAP356", "WAP357", "WAP358", "WAP359", "WAP360", "WAP361", "WAP362", "WAP363", "WAP364", "WAP365", "WAP366", "WAP367", "WAP368", "WAP369", "WAP370", "WAP371", "WAP372", "WAP373", "WAP374", "WAP375", "WAP376", "WAP377", "WAP378", "WAP379", "WAP380", "WAP381", "WAP382", "WAP383", "WAP384", "WAP385", "WAP386", "WAP387", "WAP388", "WAP389", "WAP390", "WAP391", "WAP392", "WAP393", "WAP394", "WAP395", "WAP396", "WAP397", "WAP398", "WAP399", "WAP400", "WAP401", "WAP402", "WAP403", "WAP404", "WAP405", "WAP406", "WAP407", "WAP408", "WAP409", "WAP410", "WAP411", "WAP412", "WAP413", "WAP414", "WAP415", "WAP416", "WAP417", "WAP418", "WAP419", "WAP420", "WAP421", "WAP422", "WAP423", "WAP424", "WAP425", "WAP426", "WAP427", "WAP428", "WAP429", "WAP430", "WAP431", "WAP432", "WAP433", "WAP434", "WAP435", "WAP436", "WAP437", "WAP438", "WAP439", "WAP440", "WAP441", "WAP442", "WAP443", "WAP444", "WAP445", "WAP446", "WAP447", "WAP448", "WAP449", "WAP450", "WAP451", "WAP452", "WAP453", "WAP454", "WAP455", "WAP456", "WAP457", "WAP458", "WAP459", "WAP460", "WAP461", "WAP462", "WAP463", "WAP464", "WAP465", "WAP466", "WAP467", "WAP468", "WAP469", "WAP470", "WAP471", "WAP472", "WAP473", "WAP474", "WAP475", "WAP476", "WAP477", "WAP478", "WAP479", "WAP480", "WAP481", "WAP482", "WAP483", "WAP484", "WAP485", "WAP486", "WAP487", "WAP488", "WAP489", "WAP490", "WAP491", "WAP492", "WAP493", "WAP494", "WAP495", "WAP496", "WAP497", "WAP498", "WAP499", "WAP500", "WAP501", "WAP502", "WAP503", "WAP504", "WAP505", "WAP506", "WAP507", "WAP508", "WAP509", "WAP510", "WAP511", "WAP512", "WAP513", "WAP514", "WAP515", "WAP516", "WAP517", "WAP518", "WAP519", "WAP520", "LONGITUDE", "LATITUDE", "FLOOR", "BUILDINGID", "SPACEID", "RELATIVEPOSITION", "USERID", "PHONEID", "TIMESTAMP"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];
% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
% Import the data
eval_table = readtable(append(data_dir, '\validationData.csv'), opts);
clear opts;

% Creating table of X variables from training and testing data
regressors_train = train_table(:, 1:520);
regressors_eval = eval_table(:, 1:520);

% Creating regressands from training and testing data
long_train = table2array(train_table(:, 'LONGITUDE'))./1000;
lat_train = table2array(train_table(:, 'LATITUDE'))./100000;
floor_train = table2array(train_table(:, 'FLOOR'));
bldg_train = table2array(train_table(:, 'BUILDINGID'));
long_eval = table2array(eval_table(:, 'LONGITUDE'))./1000;
lat_eval = table2array(eval_table(:, 'LATITUDE'))./100000;
floor_eval = table2array(eval_table(:, 'FLOOR'));
bldg_eval = table2array(eval_table(:, 'BUILDINGID'));

% Returning to home directory
cd(home_dir);
%==========================================================================

%==========================================================================
%% Training and evaluation ensemble of boosted learner trees for both classification (building ID and floor number) and regression (latitude and longitude) problems
% Creating template tree for classification trees
t = ClassificationTree.template('MinLeaf', 1);
%===================================================
% Variable to predict: Floor number (classification)
%===================================================
% Fitting ensemble of learner trees
tic;
disp('Training algorithm for floor');
floor_mdl = fitcensemble(regressors_train, floor_train, 'Method', 'Bag', 'NumLearningCycles', 40, 'Learners', t);
toc;

% Evaluating model
mse_floor = loss(floor_mdl, regressors_eval, floor_eval, 'mode', 'ensemble');
fprintf('Mean-square testing error = %f\n', mse_floor);

%==================================================
% Variable to predict: Building ID (classification)
%==================================================
% Fitting ensemble of learner trees
tic;
disp('Training algorithm for building ID');
bldg_mdl = fitcensemble(regressors_train, bldg_train, 'Method', 'Bag', 'NumLearningCycles', 19, 'Learners', t);
toc;

% Evaluating model
mse_bldg = loss(bldg_mdl, regressors_eval, bldg_eval, 'mode', 'ensemble');
fprintf('Mean-square testing error = %f\n', mse_bldg);

% Creating template tree for regression trees
t = RegressionTree.template('MinLeaf', 5);
%===========================================
% Variable to predict: Latitude (regression)
%===========================================
% Fitting ensemble of learner trees
tic;
disp('Training algorithm for latitude');
lat_mdl = fitrensemble(regressors_train, lat_train, 'Method', 'Bag', 'NumLearningCycles', 13, 'Learners', t);
toc;

% Evaluating model
mse_lat = loss(lat_mdl, regressors_eval, lat_eval, 'mode', 'ensemble');
fprintf('Mean-square testing error = %f\n', mse_lat);

%============================================
% Variable to predict: Longitude (regression)
%============================================
% Fitting ensemble of learner trees
tic;
disp('Training algorithm for longitude');
long_mdl = fitrensemble(regressors_train, long_train, 'Method', 'Bag', 'NumLearningCycles', 42, 'Learners', t);
toc;

% Evaluating model
mse_long = loss(long_mdl, regressors_eval, long_eval, 'mode', 'ensemble');
fprintf('Mean-square testing error = %f\n', mse_long);
%==========================================================================

%==========================================================================
%% Examining learning capability of algorithm by plotting actual values in training data with model fitted values for training data
%===================================================
% Variable to predict: Floor number (classification)
%===================================================
plot(floor_train, 'b', 'LineWidth', 2);
hold on;
plot(predict(floor_mdl, regressors_train), 'r.-', 'LineWidth', 1,'MarkerSize', 15);

% Only showing first 100 points. One has the option to pan to view more
% Observe first hundred points, pan to view more
xlim([0 100]);
legend({'Actual', 'Predicted'}, 'Location', 'Best');
xlabel('Training data point');
ylabel('Floor number');
hold off;
saveas(gcf, 'path\to\graphics\floor_1_plot.png');
close(gcf);

%==================================================
% Variable to predict: Building ID (classification)
%==================================================
plot(bldg_train, 'b', 'LineWidth', 2);
hold on;
plot(predict(bldg_mdl, regressors_train), 'r.-', 'LineWidth', 1,'MarkerSize', 15);

% Only showing first 100 points. One has the option to pan to view more
% Observe first hundred points, pan to view more
xlim([0 100]);
legend({'Actual', 'Predicted'}, 'Location', 'Best');
xlabel('Training data point');
ylabel('Building ID');
hold off;
saveas(gcf, 'path\to\graphics\bldg_1_plot.png');
close(gcf);

%===========================================
% Variable to predict: Latitude (regression)
%===========================================
plot(lat_train, 'b', 'LineWidth', 2);
hold on;
plot(predict(lat_mdl, regressors_train), 'r.-', 'LineWidth', 1,'MarkerSize', 15);

% Only showing first 100 points. One has the option to pan to view more
% Observe first hundred points, pan to view more
xlim([0 100]);
legend({'Actual', 'Predicted'}, 'Location', 'Best');
xlabel('Training data point');
ylabel('Latitude');
hold off;
saveas(gcf, 'path\to\graphics\lat_1_plot.png');
close(gcf);

%============================================
% Variable to predict: Longitude (regression)
%============================================
plot(long_train, 'b', 'LineWidth', 2);
hold on;
plot(predict(long_mdl, regressors_train), 'r.-', 'LineWidth', 1,'MarkerSize', 15);

% Only showing first 100 points. One has the option to pan to view more
% Observe first hundred points, pan to view more
xlim([0 100]);
legend({'Actual', 'Predicted'}, 'Location', 'Best');
xlabel('Training data point');
ylabel('Longitude');
hold off;
saveas(gcf, 'path\to\graphics\long_1_plot.png');
close(gcf);
%==========================================================================

%==========================================================================
%% Plotting MSE of training data, validation data, and out-of-bag with respect to number of boosted learner trees
%===================================================
% Variable to predict: Floor number (classification)
%===================================================
% Creating training error function
training_error = resubLoss(floor_mdl, 'mode', 'cumulative');
% Creating testing error function
test_error = loss(floor_mdl, regressors_eval, floor_eval, 'mode', 'cumulative');
% Creating out-of-bag error function
oob_error = oobLoss(floor_mdl, 'mode', 'cumulative');

plot(training_error);
hold on;
plot(test_error, 'r');
plot(oob_error, 'k--');
title({'Classification error (mean squared error),', 'predicting floor number'});
legend({'Training set error', 'Test set error', 'Out-of-bag error'}, 'Location', 'Best');
xlabel('Number of trees');
ylabel('Mean squared error');
set(gcf,'Position', [249 634 1009 420]);
hold off;
saveas(gcf, 'path\to\graphics\floor_2_plot.png');
close(gcf);

%==================================================
% Variable to predict: Building ID (classification)
%==================================================
% Creating training loss function
training_error = resubLoss(bldg_mdl, 'mode', 'cumulative');
% Creating testing loss function
test_error = loss(bldg_mdl, regressors_eval, bldg_eval, 'mode', 'cumulative');
% Creating out-of-bag error function
oob_error = oobLoss(bldg_mdl, 'mode', 'cumulative');

plot(training_error);
hold on;
plot(test_error, 'r');
plot(oob_error, 'k--');
title({'Classification error (mean squared error),', 'predicting building ID'});
legend({'Training set loss', 'Test set loss', 'Out-of-bag error'}, 'Location', 'Best');
xlabel('Number of trees');
ylabel('Mean squared error');
set(gcf,'Position', [249 634 1009 420]);
hold off;
saveas(gcf, 'path\to\graphics\bldg_2_plot.png');
close(gcf);

%===========================================
% Variable to predict: Latitude (regression)
%===========================================
% Creating training loss function
training_error = resubLoss(lat_mdl, 'mode', 'cumulative');
% Creating testing loss function
test_error = loss(lat_mdl, regressors_eval, lat_eval, 'mode', 'cumulative');
% Creating out-of-bag error function
oob_error = oobLoss(lat_mdl, 'mode', 'cumulative');

plot(training_error);
hold on;
plot(test_error, 'r');
plot(oob_error, 'k--');
title({'Regression error (mean squared error),', 'predicting latitude'});
legend({'Training set error', 'Test set error', 'Out-of-bag error'}, 'Location', 'Best');
xlabel('Number of trees');
ylabel('Mean squared error');
set(gcf,'Position', [249 634 1009 420]);
hold off;
saveas(gcf, 'path\to\graphics\lat_2_plot.png');
close(gcf);

%============================================
% Variable to predict: Longitude (regression)
%============================================
% Creating training loss function
training_error = resubLoss(long_mdl, 'mode', 'cumulative');
% Creating testing loss function
test_error = loss(long_mdl, regressors_eval, long_eval, 'mode', 'cumulative');
% Creating out-of-bag error function
oob_error = oobLoss(long_mdl, 'mode', 'cumulative');

plot(training_error);
hold on;
plot(test_error, 'r');
plot(oob_error, 'k--');
title({'Regression error (mean squared error),', 'predicting longitude'});
legend({'Training set error', 'Test set error', 'Out-of-bag error'}, 'Location', 'Best');
xlabel('Number of trees');
ylabel('Mean squared error');
set(gcf,'Position', [249 634 1009 420]);
hold off;
saveas(gcf, 'path\to\graphics\long_2_plot.png');
close(gcf);
%==========================================================================