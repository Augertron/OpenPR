libname = 'openpr_plsi';

names = ['plsiread','int_readdata'; 'plsitrain','int_plsitrain'];

gateway_path = get_absolute_file_path('builder_gateway_plsi.sce');

files = (listfiles(['*.h'; '*.cpp']))';

//if ~MSDOS then
//	hfiles = (listfiles('*.h'))';
//	files = [hfiles, files];
//end

tbx_build_gateway(libname, names, files, gateway_path);

clear libname names gateway_path files tbx_build_gateway;



