sci_gateway_dir = get_absolute_file_path('builder_gateway.sce');

tbx_builder_gateway_lang('svm', sci_gateway_dir);
tbx_builder_gateway_lang('plsi', sci_gateway_dir);
tbx_builder_gateway_lang('qdmatch', sci_gateway_dir);
tbx_builder_gateway_lang('knearest', sci_gateway_dir);
tbx_builder_gateway_lang('kmeans', sci_gateway_dir);
tbx_builder_gateway_lang('nbayes', sci_gateway_dir);
tbx_builder_gateway_lang('em', sci_gateway_dir);


tbx_build_gateway_loader(['svm', 'plsi', 'qdmatch', 'knearest', 'kmeans', 'nbayes', 'em'], sci_gateway_dir);

clear tbx_builder_gateway_lang tbx_build_gateway_loader;
clear sci_gateway_dir;


