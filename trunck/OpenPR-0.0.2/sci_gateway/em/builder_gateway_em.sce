gateway_path = get_absolute_file_path('builder_gateway_em.sce');

cur_path = pwd();
chdir(gateway_path);

ilib_name  = 'openpr_em';

// objects files (but do not give mexfiles here) 
files = ['_ml.h', 'ml.h', 'transformation.h', 'ml.cpp', 'mlem.cpp', 'ml_inner_functions.cpp', 'transformation.cpp', 'emtrain.cpp', 'empredict.cpp'];

// table of (scilab_name,interface-name or mexfile-name, type) 
table =['emtrain',  'emtrain',  'cmex';
		'empredict','empredict','cmex'];

if ~MSDOS then
//	files = ['_ml.h', files];
	libs = [];
	opencv_version = unix_g('pkg-config --modversion opencv');
	if( length(opencv_version) == 0 | ( strtod( strsubst(opencv_version, '.', '')) < 200 ) )
		error(gettext("OpenCV (version 2.0.0) is needed for compiling OpenPR."));
	end;
	cflags = unix_g('pkg-config --cflags opencv');
	ldflags = unix_g('pkg-config --libs opencv');
else
    other_lib_path = SCI+"/contrib/OpenPR-0.0.2/opencv/lib/win";
    libs  = [other_lib_path+"/cxcore200"]; 		
    ldflags = "-LIBPATH:"""+SCI+"/contrib/OpenPR-0.0.2/opencv/lib/win"" ";
    cflags = "-I"""+SCI+"/contrib/OpenPR-0.0.2/opencv/include"" ";
end

ilib_mex_build(ilib_name,table,files,libs,'',ldflags,cflags);

clear gateway_path ilib_name files table other_lib_path libs ldflags cflags ilib_mex_build;

chdir(cur_path);

