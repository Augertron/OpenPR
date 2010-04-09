//
openpr_demo_path = get_absolute_file_path('openpr.dem.gateway.sce');

subdemolist = ["confusion matrix to normalized mutual information", "confmatrix2ni_mi.dem.sce";..
			   "wpca", "wpca.dem.sce";..
			   "kmeans", "kmeans.dem.sce";..			   
			   "agglomerative hierarchical clustering", "ahclustering.dem.sce";..
			   "basic leader-follower clustering", "leader_follower.dem.sce";..
			   "agglomerative mean-shift clustering", "aggloms.dem.sce"];   //"svmtrain", "svmtrain.dem.sce";..

subdemolist(:, 2) = openpr_demo_path + subdemolist(:, 2);
