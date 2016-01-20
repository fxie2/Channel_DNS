function save_data(t)
s = num2str(t);
s = strcat('timestep_',s);
mkdir('..\data',s);
save('p.mat','p');
save('u.mat','u');
save('v.mat','v');
save('w.mat','w');
save('t.mat','t');
end