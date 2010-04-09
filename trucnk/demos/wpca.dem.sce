//
mode(-1);
lines(0);

my_handle = scf(100001);
clf(my_handle,"reset");
demo_viewCode("wpca.dem.sce");

// DEMO START

my_plot_desc          = "wpca";
my_handle.figure_name = my_plot_desc;

rotation = [7 -cos(3.14/4);sin(3.14/4) 1];
Train_Patterns = rand(2,100);
Train_Patterns = rotation* Train_Patterns;                  
weight = ones(100,1) ;                                    //equal weight
plot2d(Train_Patterns(1,:),Train_Patterns(2,:),style=-4); //plot the data points
[eig_vec,m_vec,eig_val] = wpca(Train_Patterns,weight);
plot2d(m_vec(1),m_vec(2),style=-5);                       //plot the mean vector
//plot the first principal component
x=-1:0.1:7;
y=(eig_vec(2,1)/eig_vec(1,1))*(x-m_vec(1))+m_vec(2);
plot2d(x,y,style=2);
legends(["data";"data center";"eigen vector"],[-4,-5,2], with_box=%f, opt="?")

// DEMO END
