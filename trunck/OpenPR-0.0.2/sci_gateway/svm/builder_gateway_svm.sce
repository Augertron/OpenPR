lib_name = 'openpr_svm';

table = ['readsparse', 'read_sparse', 'cmex';
		 'svmtrain',   'svmtrain',    'cmex';
		 'svmpredict', 'svmpredict', 'cmex'];
		 
files = ['svm.h', 'svm_model_scilab.h', 'svm.cpp', 'svm_model_scilab.c', 'svmtrain.c', 'svmpredict.c', 'read_sparse.c'];


libs = [];

if ~MSDOS then
	if part(getenv('OSTYPE','no'),1:6)=='darwin' then
    	cflags = "-D__SCILAB__"
	    fflags = ""; 
	    ldflags= ""; 
	    cc = "g++";
    else 
		ldflags = "";
		cflags = "-lstdc++ -D__SCILAB__";
		fflags = "";
		cc = "";
	end
else
	ldflags = "";
	cflags = "-D__SCILAB__";
	fflags = "";
	cc = "";
end

ilib_mex_build(lib_name, table, files, libs, 'Makelib', ldflags, cflags, fflags, cc);

clear lib_name table files libs ldflags cflags fflags cc ilib_mex_build;

