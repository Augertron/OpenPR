//
mode(-1);
lines(0);

// DEMO START
data_path = SCI+'/contrib/OpenPR-0.0.2/etc/data/heart_scale';

[label_vector, instance_vector] = readsparse(data_path);
model = svmtrain(label_vector, instance_vector, '-c 1 -g 0.07');

disp("svm model:")
disp(model)

// DEMO END
