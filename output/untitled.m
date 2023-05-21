clear;
close;
AidDir = uigetdir(); 	% 通过交互的方式选择一个文件夹

if AidDir == 0 			% 用户取消选择
    fprintf('Please Select a New Folder!\n');
else
	cd(AidDir)
	RawFile = dir('**/*.dat'); %主要是这个结构，可以提取所有文件
	AllFile = RawFile([RawFile.isdir]==0);
    if isempty(fieldnames(AllFile))
    	fprintf('There are no files in this folder!\n');
    else	% 当前文件夹下有文件，反馈文件数量
    	fprintf('Number of Files: %i \n',size(AllFile,1));
    end
end

pic_index = 1;%记录图像编号
filename={AllFile.name}'
for i=1:size(AllFile,1)
    data = importdata(filename{i}).data;
    sizes = size(data);
    plot(data(1,:));
    ylim([-1, 1])
    title(num2str(i))
    %抓取当前的figure，保存为rgb图像后，再转化为索引图像
    [A,map] = rgb2ind(frame2im(getframe),256);
%     if pic_index == 1
%       imwrite(A,map,'test.gif','gif','Loopcount',inf,'DelayTime',0.2);
%     else
%       imwrite(A,map,'test.gif','gif','WriteMode','append','DelayTime',0.2);
%     end
    pic_index = pic_index + 1;
    pause(0.1)
end

data = importdata(filename{1}).data;
imagesc(data)
colorbar