function [edges2] = risingEdge2(time, data)
threshold = -.05;
offsetData = [data(2:end); NaN];
edges2 = find(data < threshold & offsetData > threshold);
% plot(time, data);
% hold on;
% plot(time(risingEdge, threshold, 'rx');
end