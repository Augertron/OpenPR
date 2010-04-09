
help_dir = get_absolute_file_path('builder_help.sce');

tbx_builder_help_lang("en_US", help_dir);

helpdoc_path = help_dir+'en_US';
//xmltopdf(helpdoc_path, 'OpenPR Manual', 'en_US');
//xmltohtml(helpdoc_path, 'OpenPR Manual', 'en_US');

clear help_dir;


