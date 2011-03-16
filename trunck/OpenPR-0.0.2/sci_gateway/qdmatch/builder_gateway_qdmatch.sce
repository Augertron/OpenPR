libname = 'openpr_qdmatch';

names = ['qdmatch', 'int_qdmatch'];

gateway_path = get_absolute_file_path('builder_gateway_qdmatch.sce');

files = [(listfiles([gateway_path+'*.h'; gateway_path+'*.cpp']))', gateway_path+'common.c', gateway_path+'NyxMat.cc'];
//files = ['common.c', cppfiles];

if ~MSDOS then 
//	hfiles = (listfiles('*.h'))';
//	files = [hfiles, 'NyxMat.cc', files];
	libs = [];
	opencv_version = unix_g('pkg-config --modversion opencv');
	if( length(opencv_version) == 0 | ( strtod( strsubst(opencv_version, '.', '')) < 200 ) )
		error(gettext("OpenCV (version 2.0.0) is needed for compiling OpenPR."));
	end;
	cflags = unix_g('pkg-config --cflags opencv');
	ldflags = unix_g('pkg-config --libs opencv');
else
	other_lib_path = SCI+"/contrib/OpenPR-0.0.2/opencv/lib/win";
	libs = [other_lib_path+"/cxcore200", other_lib_path+"/cv200", other_lib_path+"/highgui200", other_lib_path+"/cvaux200"];
	ldflags = "-LIBPATH:"""+SCI+"/contrib/OpenPR-0.0.2/opencv/lib/win"" ";
	cflags = "-I"""+SCI+"/contrib/OpenPR-0.0.2/opencv/include"" ";
end

tbx_build_gateway(libname, names, files, gateway_path, libs, ldflags, cflags);

clear libname names cppfiles files gateway_path other_lib_path libs ldflags cflags tbx_builde_gateway;

