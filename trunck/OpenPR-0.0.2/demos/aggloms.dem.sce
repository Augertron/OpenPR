mode(-1);
lines(0);

my_handle = scf(0);
clf(my_handle,"reset");
demo_viewCode("aggloms.dem.sce");

// DEMO START

Data = read(SCI+'/contrib/OpenPR-0.0.2/etc/data/aggloms_data',2,1810);
Label = read(SCI+'/contrib/OpenPR-0.0.2/etc/data/aggloms_label',1,1810);

sigma = 0.4; // kernel bandwidth
ite_num = 60; // iteration times

[cluster_centers, cluster_id]=aggloms(Data', sigma, ite_num);

//-------------------- draw clusters --------------------------------

cluster_count = size(cluster_centers,1);
for (i=1:cluster_count)
    my_color = rand(1,3);
    plot(Data(1,find(cluster_id==i)), Data(2,find(cluster_id==i)),'o','marker','sq','markersize',6,'markforegroun',my_color,'markbackgro',my_color);
end
plot(cluster_centers(:,1), cluster_centers(:,2), 'gs','marker','o','markersize',10,'markbackgro','g');

// DEMO END
