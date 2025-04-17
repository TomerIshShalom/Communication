function marker=num2marker(num)
Markers={"o","+","*",".","x","_","|","square","diamond","^","v",">","<","pentagram","hexagram"};
marker=Markers{ mode(num-1,numel(Markers))+1 };