close all
%folder_index = 3:10; % everything
%folder_index = [4 10]; %6Torr h10 h15 8.6g/s
%folder_index = [3 5]; %50Torr h10 8.6g/s
%folder_index = [6 9]; %50Torr h3 8.6g/s 0.32g/s
%folder_index = [5 7]; %50Torr h10 8.6g/s 0.32g/s
%folder_index = [5 10]; %50mTorr 6Torr h10 8.6g/s
folder_index = [8 4]; %50mTorr 6Torr h15 8.6g/s

test_data = strings(8,length(folder_index));

figure
grid on 
box on
hold on

x_bound = [-250 250];
aspect_ratio = 0.2;
xlim(x_bound)
ylim(x_bound*aspect_ratio)
pbaspect([1 aspect_ratio 1])
ylabel('Vertical position (mm)')
xlabel('Horizontal position (mm)')

for i = 1:length(folder_index)
% Get a list of all files and folders in this folder.
files = dir('./');
% Get a logical vector that tells which is a directory.
dirFlags = [files.isdir];
% Extract only those that are directories.
subFolders = files(dirFlags);

k = 2 + folder_index(i);
fprintf("%s\n",subFolders(k).name)

%extracting data from d
test_data(:,i) = string(split(subFolders(k).name,'_'));

A = load(sprintf("./%s/%s_X-Slice.csv",subFolders(k).name,subFolders(k).name));

x_crater = A(:,1);
y_crater = A(:,3);

[x_crater, index_sort] = sort(x_crater);

% if i == 2 %6Torr 8.6g/s
%     x_crater = x_crater-5;
% end

y_crater = y_crater(index_sort);
y_crater = y_crater - y_crater(1);

plot(x_crater,y_crater)
end

%legend(test_data(5,:))%6Torr 8.6g/s
%legend('repetion 1','repetion 2') %50Torr h10 8.6g/s
%legend(test_data(5,:)) %50Torr h10 h15 h3 8.6g/s
%legend('8.60 g/s','0.32 g/s') %50Torr h3 8.6g/s 0.32g/s
%legend('8.60 g/s','0.32 g/s') %50Torr h10 8.6g/s 0.32g/s
legend('50 mTorr','6 Torr') %50mTorr 6Torr h10 8.6g/s

x0=100;
y0=100;
width=1200;
height=400;
set(gcf,'position',[x0,y0,width,height])
set(gca,'fontname','times')

