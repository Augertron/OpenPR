mode(-1);
lines(0);

my_handle = scf(1);
clf(my_handle,"reset");
demo_viewCode("leader_follower.dem.sce");

// DEMO START

samples = [rand(2,10), -1*rand(2,10)];
theta = 1;
[centers, labels] = leader_follower(samples, theta);
scf(1);
plot(samples(1,:), samples(2,:), 'b.', 'MarkerSize', 3);
scf(2);
style = ['g.', 'c.', 'y.', 'k.', 'm.'];
ul = unique(labels)
for i = 1:length(ul),
  plot(samples(1,find(labels==ul(i))), samples(2,find(labels==ul(i))), style(i), 'MarkerSize', 3);
  set(gca(), "auto_clear", "off");
end
plot(centers(1,:), centers(2,:), 'r.', 'MarkerSize', 4);

// DEMO END
