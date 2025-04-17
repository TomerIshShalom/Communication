function LineStyle=num2linestyle(num)
LineStyles={"-","--","-.",":"};
LineStyle=LineStyles{ mod(num-1,numel(LineStyles))+1};