//
mode(-1);
lines(0);

my_handle = scf(100001);
clf(my_handle,"reset");
demo_viewCode("kmeans.dem.sce");

// DEMO START

my_plot_desc          = "kmeans";
my_handle.figure_name = my_plot_desc;

x = [rand(100, 2)+0.5*ones(100, 2);rand(100, 2)-0.5*ones(100, 2)];
idx = kmeans(x, 2);
plot(x(idx==1, 1), x(idx==1, 2), 'r.', 'MarkerSize', 5);
set(gca(),"auto_clear","off")
plot(x(idx==2, 1), x(idx==2, 2), 'b.', 'MarkerSize', 5);

// DEMO END
