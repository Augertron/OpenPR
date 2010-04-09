mode(-1);
lines(0);

my_handle = scf(0);
clf(my_handle,"reset");
demo_viewCode("ahclustering.dem.sce");

// DEMO START

samples = rand(2,30);
cluster_num = 3;
dist_type = 'avg';
[centers, labels] = ahclustering(samples, cluster_num, dist_type);
//show figures
//scf(0);
plot(samples(1,:), samples(2,:), 'b.', 'MarkerSize', 3);
scf(1);
clusters = unique(labels);
plot(samples(1,find(labels==clusters(1))),samples(2,find(labels==clusters(1))),'g.','MarkerSize',3);
set(gca(),"auto_clear","off");
plot(samples(1,find(labels==clusters(2))),samples(2,find(labels==clusters(2))),'y.','MarkerSize',3);
set(gca(),"auto_clear","off");
plot(samples(1,find(labels==clusters(3))),samples(2,find(labels==clusters(3))),'m.','MarkerSize',3);
set(gca(),"auto_clear","off");
plot(centers(1,:), centers(2,:), 'r.', 'MarkerSize', 4);

// DEMO END
