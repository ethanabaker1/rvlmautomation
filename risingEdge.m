function [edges] = risingEdge(time, data)
threshold = 5;
offsetData = [data(2:end); NaN];
edges = find(data < threshold & offsetData > threshold);
% plot(time, data);
% hold on;
% plot(time(risingEdge, threshold, 'rx');
end