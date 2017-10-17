function f = f_anatplot(templates,RSN_names)
load DSI_enhanced
MNI_coord=tal2mni(talairach);
clear talairach
load BrainImageVertFaces

[N,K] = size(templates);
levels = 3;

% labels are not ordered! therefore order back
load Human_66 Order
x = 1:N;
x_Order = x(Order);
[~,reOrder] = sort(x_Order,'ascend');

% design the markers
[x,y,z] = sphere; % unit size sphere (made up of points along circles in different locations in 3D space)
cc = 3; % diameter of the circles

persps = [0 90; 90 0; -125 45]; % perspective - 1: top, 2: right, 3: "sideways"
for persp=1
    f = figure;
    f.Units = 'centimeters';
    f.PaperPositionMode = 'auto';
    for k=1:K % use three points of view
        subplot(1,K,k)
        title(RSN_names{k})
        %f = figure;
        %f.Units = 'centimeters';
        %f.PaperPositionMode = 'auto';
        activity = templates(reOrder,k);
        activity = activity/max(activity);
        qtemp = zeros(size(activity));
        for i=1:size(activity,2)
            limits = linspace(min(activity(:,i)),max(activity(:,i)),levels+1);
            partition = limits(2:end-1);
            qtemp(:,i) = quantiz(activity(:,i),partition);
        end
        
        acth = 0;
        
        %plots a mesh surface for the cortex
        patch('vertices',cortexvert,'faces',cortexface,'EdgeColor','none',...
            'FaceColor',[0.9 0.9 0.9],'FaceAlpha',0.1);
        hold on
        
        % not sure where the CLim comes from, but the colors definitely correspond
        % to the values in the PC vector; I suppose the CTickValues (if existent?)
        % could be changed
        for j = 1:N
            n = find(roi_lbls==j);
            for i = 1:length(n)
                act = activity(j);
                if abs(act) > acth
                    clrs = get(gca,'DefaultAxesColorOrder');
                    if act<0
                        clr = clrs(1,:); % blue
                    else
                        clr = clrs(2,:); % orange
                    end
                    surf(cc*x+MNI_coord(n(i),1),cc*y+MNI_coord(n(i),2),cc*z+MNI_coord(n(i),3),...
                        'FaceColor',clr,'EdgeColor','none','FaceAlpha',0.8*abs(act));
                end
            end
        end
        
        %set(gca,'XTick',[],'YTick',[])
        
        %set(gca,'ZLim',[10 78],'XLim',[-90 90])
        axis equal
        box off
        axis off
        %axis equal
        %colorbar
        %set(gca,'CLim',[act_rng(1),act_rng(end)]);
        
        view(persps(persp,:))
        pbaspect([1 1 1]);
        pbaspect([1 1 1]);
        %rotate3d;
        % view([0,90]) % top
        % view([90 -90]) % ventral
        % view([90 45]) % R sideways
        % view([-180 0]) % L sideways
        % view([45 20]) % perspective
        % view([90 0]) % front
        material dull;
        lighting phong;
        camlight;
        rotate3d;
    end
end
set(gcf,'color','white')
end