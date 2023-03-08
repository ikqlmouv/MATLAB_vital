[text_area_1,text_area_2,axes_1] = setup();

%Add code here

function [text_area_1,text_area_2,axes_1] = setup(comp)
   
    fig = uifigure('Name' , 'Display' , 'Position' , [500 200 1200 650]);

    uilabel(fig , "Position", [20 483 200 45] , 'Text', 'Heart rate','FontSize', 36 , 'HorizontalAlignment', 'center');
    uilabel(fig , "Position", [20 183 200 45] , 'Text', 'Breath rate','FontSize', 36, 'HorizontalAlignment', 'center') ;
    
    text_area_1 = uitextarea(fig);
    text_area_1.Position = [50 415 150 50];
    text_area_1.Value = '0';
    text_area_1.HorizontalAlignment = 'center';
    text_area_1.FontSize = 32;

    text_area_2 = uitextarea(fig);
    text_area_2.Position = [50 120 150 50];
    text_area_2.Value = '0';
    text_area_2.HorizontalAlignment = 'center';
    text_area_2.FontSize = 32;

    axes_1 = uiaxes(fig,'Position',[250 75 900 500]);
    axes_1.XLabel.String = 'time (s)';
    axes_1.XLabel.FontSize = 24;
    axes_1.YLabel.String = 'Phase';
    axes_1.YLabel.FontSize = 24;
    
end
