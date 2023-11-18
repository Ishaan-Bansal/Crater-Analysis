function plotCrossSections(folder_index, legend_labels, legend_title, tit)
    dir_name = pwd;
    cd X-Slices\

    % test_data = strings(8, length(folder_index));

    figure
    grid on
    box on
    hold on

    x_bound = [-250 250];
    aspect_ratio = 0.2;
    xlim(x_bound)
    ylim(x_bound * aspect_ratio)
    pbaspect([1 aspect_ratio 1])
    ylabel('Vertical position (mm)')
    xlabel('Horizontal position (mm)')
    for i = 1:length(folder_index)
        % Get a list of all files and folders in this folder.
        files = dir('./');
        % Get a logical vector that tells which is a directory.
        dirFlags = [files.isdir];
        % Extract only those that are files.
        subFiles = files(~dirFlags);

        k = folder_index(i);
        fprintf("%s\n", subFiles(k).name)

        % Extracting data from d
        % test_data(:, i) = string(split(subFiles(k).name, '_'));

        A = load(subFiles(k).name);

        x_crater = A(:, 1);
        y_crater = A(:, 3);

        [x_crater, index_sort] = sort(x_crater);

        y_crater = y_crater(index_sort);
        y_crater = y_crater - y_crater(1);

        plot(x_crater, y_crater)
    end
    
    cd(dir_name);
    leg = legend(legend_labels);
    title(leg, legend_title)

    title(tit)

    x0 = 100;
    y0 = 100;
    width = 1200;
    height = 400;
    set(gcf, 'position', [x0, y0, width, height])
    set(gca, 'fontname', 'times')
    saveas(gcf, [tit + '.svg']);
    saveas(gcf, [tit + '.fig']);
end
