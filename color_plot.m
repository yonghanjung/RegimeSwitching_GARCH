function color_plot(data_vector, color_vector)
styles={'ro','g.','bx','kd'};
hold off;
for i=unique(color_vector)
    thisIdx=find(color_vector==i);
    thisY=data_vector(color_vector==i);
    thisStyle=styles{mod(i-1,numel(styles))+1}; 
    plot(thisIdx,thisY,thisStyle);
    hold on;
end
hold off;
