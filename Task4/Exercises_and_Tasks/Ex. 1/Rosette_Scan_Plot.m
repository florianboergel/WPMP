function [xi,yi,vlosi] = Rosette_Scan_Plot(x,y,vlos,varargin)

coloraxis = [4 16];
colorbarstring = '$$v\rm_{LOS}\,\,\,[m\,\,s^{-1}]$$';
figuretitle = [];
data_availability = 0;
x_fullscan = [];
y_fullscan = [];
meshgridvector = -30:2:30;

if ~isempty(varargin)
    for i=1:length(varargin)
        if strcmp(varargin{i},'coloraxis')
            coloraxis=varargin{i+1};
        end
        if strcmp(varargin{i},'colorbarstring')
            colorbarstring=varargin{i+1};
        end
        if strcmp(varargin{i},'title')
            figuretitle=varargin{i+1};
        end
        if strcmp(varargin{i},'availability')
            data_availability=varargin{i+1};
        end
        if strcmp(varargin{i},'rosettebackground')
            x_fullscan=varargin{i+1};
            y_fullscan=varargin{i+2};
        end
        if strcmp(varargin{i},'meshgridvector')
            meshgridvector=varargin{i+1};
        end
    end
end
        
figure
[xi,yi] = meshgrid(meshgridvector,meshgridvector);
vlosi = griddata(x,y,vlos,xi,yi);
surf(xi,yi,0*yi,vlosi);
shading interp
hold on
    if ~isempty(x_fullscan)
        plot3(x_fullscan,y_fullscan,-ones(size(x_fullscan)),'k--')
    end
plot3(x,y,-ones(size(x)),'ko','linewidth',2)
xlabel('\ity \rm[m]','fontsize',14)
ylabel('\itz \rm[m]','fontsize',14)
    if ~isempty(figuretitle)
        title(figuretitle,'fontsize',14);
    else
    end
h = colorbar('peer',gca);
set(get(h,'ylabel'),'String',colorbarstring,'Interpreter','Latex','fontsize',16,'FontName','Arial');
caxis(coloraxis)
axis equal
set(gca,'fontsize',14)
view([-180 -90])
zlim([-2 2])
    if data_availability>0
        datastring=sprintf('%s%d%s','Data availability: ',round(data_availability),'%');
        text(15,-25,-3,datastring,'fontsize',14,'BackgroundColor','w','EdgeColor','k')
    else
    end

end