function print_SCARA_kin(l_mat, R_mat)
% Function that prints a schematic representation
% of the SCARA robot seen from two different point of view,
% given the relative transformations between links

%% set the links lengths
l_links = [0.5-0.005, 0.5, 0.4, 0.3, 0.2];

%% get the extremes of each link

% allocate the space
start_points = zeros(5,3);
end_points = zeros(5,3);

% first link
end_points(1,:) = [0,0,l_links(1)];

% second link
end_points(2,:) = l_mat(:,1);
start_points(2,:) = l_mat(:,1) - R_mat(:,1:3)*[l_links(2);0;0];

% third link
end_points(3,:) = l_mat(:,2);
start_points(3,:) = l_mat(:,2) - R_mat(:,4:6)*[l_links(3);0;0];

% fourth link
end_points(4,:) = l_mat(:,3);
start_points(4,:) = l_mat(:,3)- R_mat(:,7:9)*[0;0;-l_links(4)];

% fifth link
end_points(5,:) = l_mat(:,4);
start_points(5,:) = l_mat(:,4)- R_mat(:,10:12)*[l_links(5);0;0];

%% 2D print

figure
subplot(1,2,1)
hold on
grid on
title('top view')
for i = 1 : 5
    x = [start_points(i,1),end_points(i,1)];
    y = [start_points(i,2),end_points(i,2)];
    plot(x,y, 'LineWidth',4)
end
xlim([-1.2,1.2])
ylim([-1.1,1.1])
legend('link 0','link 1', 'link 2', 'link 3', 'link 4')
xlabel('X axis')
ylabel('Y axis')
subplot(1,2,2)
hold on
grid on
title('lateral view')
for i = 1 : 5
    z = [start_points(i,3),end_points(i,3)];
    x = [start_points(i,1),end_points(i,1)];
    plot(x,z,'LineWidth',4)
end
xlim([-1.2,1.2])
ylim([0,0.8])
legend('link 0','link 1', 'link 2', 'link 3', 'link 4')
xlabel('X axis')
ylabel('Z axis')