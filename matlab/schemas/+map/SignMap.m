%{
# Retinotopy Sign map
-> map.RetMap
---
sign_map                 : mediumblob      # sign map of the retinotopy
parameters               : mediumblob      # parameters of the sign map extraction
%}

classdef SignMap < dj.Imported
    methods(Access=protected)
        function makeTuples(obj,key) %create ref map
            insert( obj, key );
        end
    end
    
    methods
        function extractSign(self, key, varargin)
            
            params.manual = 1;
            params.pexp = 1.5;
            
            params = ne7.mat.getParams(params,varargin);
            
            % define functions & colors
            cm = rgb2hsv(parula(100));
            cm = cm(:,1);
            
            % set/get parameters
            if ~exists(self & key)
                sign_params.init_gauss = [0.1 0.1 200];
                sign_params.grad_gauss = [0.1 0.1 20]; % initial map gauss filter in sd parameter
                sign_params.diff_gauss = [12 0.1 20]; % gradient diff map gauss filter in sd parameter
                sign_params.diff_open = [0 0 20]; % gradient diff imopen param
                MAP = [];
            else
                [MAP, sign_params] = fetch1(self & key,'sign_map', 'parameters');
            end
            
            % initialize
            MAP = [];
            running = true;
            control = [];
            current_param= [];
            current_field = [];
            
            % fetch horizontal & vertical maps
            if ~exists(map.OptImageBar & (map.RetMapScan & key) & 'axis="horizontal"'); return;end
            [H, A1, Ves] = fetch1(map.OptImageBar & (map.RetMapScan & key) & 'axis="horizontal"','ang','amp','vessels');
            [Hor(:,:,1),Hor(:,:,2),Hor(:,:,3)] = plot(map.OptImageBar & (map.RetMapScan & key) & 'axis="horizontal"','exp',params.pexp,'sigma',sign_params.init_gauss(1));
            if ~exists(map.OptImageBar & (map.RetMapScan & key) & 'axis="vertical"'); return;end
            [V, A2] = fetch1(map.OptImageBar & (map.RetMapScan & key) & 'axis="vertical"','ang','amp');
            [Ver(:,:,1),Ver(:,:,2),Ver(:,:,3)] = plot(map.OptImageBar & (map.RetMapScan & key) & 'axis="vertical"','exp',params.pexp,'sigma',sign_params.init_gauss(1));
            Amp = ne7.mat.normalize(self.replaceNaNs(A1+A2));
            Ves = ne7.mat.normalize(self.replaceNaNs(single(Ves)));
            Ves = ones(size(Ves));
            Hor(:,:,3) = ones(size(Hor,1),size(Hor,2));
            Ver(:,:,3) = ones(size(Ver,1),size(Ver,2));
            Hor(:,:,2) = ones(size(Hor,1),size(Hor,2));
            Ver(:,:,2) = ones(size(Ver,1),size(Ver,2));
            V = wrapTo2Pi(V);
            H = wrapTo2Pi(H);
            H = self.replaceNaNs(H);
            V = self.replaceNaNs(V);
            createMAP
            
            if params.manual
                % plot
                f = figure;
                screensize = get( groot, 'Screensize' );
                set(f,'KeyPressFcn',@EvalEvent,'position',[0.3*screensize(3) 0.3*screensize(4) 0.7*screensize(4) 0.7*screensize(4)],...
                    'windowscrollWheelFcn',@scroll,'WindowButtonMotionFcn', @mouseMove)
                subplot(2,2,1)
                imagesc(hsv2rgb(Hor)); axis image; axis off; title('Horizontal Retinotopy')
                subplot(2,2,2)
                imagesc(hsv2rgb(Ver)); axis image; axis off; title('Vertical Retinotopy')
                plotMAP

                % Add the slider and slider label text to the figure
                fields = fieldnames(sign_params);
                positions = repmat([0.6,nan,0.3,0.05],length(fields),1);
                positions(:,2) = linspace(0.4,0.1,length(fields));
                for ifield = 1:length(fields)
                    uicontrol('Parent',f,'Style','text','Units', 'normal','Position',positions(ifield,:)+[0 0.05 0 0],...
                        'String',fields{ifield});
                    control(ifield)=uicontrol('Parent',f,'Style','slider','Units', 'normal','Position',positions(ifield,:),...
                        'value',sign_params.(fields{ifield})(1), ...
                        'min',sign_params.(fields{ifield})(2), 'max',sign_params.(fields{ifield})(3),...
                        'callback',@update_val,'KeyPressFcn',@EvalEvent, 'TooltipString',fields{ifield});
                end

                % wait until done
                while running
                    try if ~ishandle(f);MAP = [];break;end;catch;break;end
                    pause(0.1);
                end
            end
            
            % insert/update tuple
            if ~isempty(MAP)
                if exists(self & key)
                    update(self & key,'sign_map',MAP)
                    update(self & key,'parameters',sign_params)
                else
                    key.sign_map = MAP;
                    key.parameters = sign_params;
                    makeTuples(self,key)
                end
            end
            
            function scroll(~,event)
                new_value = (-event.VerticalScrollCount)*0.1 + get(control(current_field),'Value');
                if ~isempty(current_param)
                    if new_value>sign_params.(current_param)(2) && new_value<sign_params.(current_param)(3)
                        sign_params.(current_param)(1) = new_value;
                        set(control(current_field),'Value',new_value)
                        plotMAP
                    end
                end
            end
            
            function mouseMove(~, ~)
                C = get (gcf, 'CurrentPoint');
                pos = get (gcf,'position');
                X = C(1,1)/pos(3);
                Y = C(1,2)/pos(4);
                for ifield = 1:length(fields)
                    if X>positions(ifield,1) && X<positions(ifield,1)+positions(ifield,3) && Y>positions(ifield,2) && Y<positions(ifield,2)+positions(ifield,4)
                        current_field = ifield;
                        current_param = get(control(ifield),'TooltipString');
                    end
                end
            end
            
            function EvalEvent(~, event)
                switch event.Key
                    case 'escape'
                        yes = questdlg('Do you wish to abort segmenting this image?','Finish segmentation', 'yes','no','yes');
                        if strcmpi('yes', yes)
                            running = false;
                            MAP = [];
                            close(gcf)
                        end
                    case 'return'
                        yes = questdlg('Ready to commit?','Finish segmentation', 'yes','no','yes');
                        if strcmpi('yes', yes)
                            running = false;
                            close(gcf)
                        end
                    otherwise
                        disp '---'
                        disp(['key: "' event.Key '" not assigned!'])
                end
            end
            
            function update_val(source,~)
                val = source.Value;
                if strcmp(source.TooltipString,'diff_open'); val = round(val);end
                sign_params.(source.TooltipString)(1) = val;
                plotMAP
            end
            
            function createMAP
                % calculate gradients
                [~,dH] = imgradient(imgaussfilt((H),sign_params.init_gauss(1)));
                [~,dV] = imgradient(imgaussfilt((V),sign_params.init_gauss(1)));
                
                % filter gradients
                dH = imgaussfilt(dH,sign_params.grad_gauss(1));
                dV = imgaussfilt(dV,sign_params.grad_gauss(1));
                grad_diff = imopen(imgaussfilt(sind(dH  - dV),sign_params.diff_gauss(1)),...
                    strel('disk',round(sign_params.diff_open(1))));
                
                % calculate maps
                MAP = round(ne7.mat.normalize(grad_diff)*length(cm));
                MAP(MAP>length(cm))= length(cm);
                MAP(MAP<1)= 1;
                MAP = cm(MAP);
                MAP(:,:,2) = Amp;
                MAP(:,:,3) = Ves;
            end
            
            function plotMAP
                createMAP
                subplot(2,2,3)
                cla
                imagesc(hsv2rgb(MAP)); axis image; axis off; title('Visual Field Sign Map'); colormap jet
            end
        end
        
        function masks = extractMasks(self, varargin)
            
            sign_params.sign_thr = 1; % threshold as standard deviations of signmap
            sign_params.step = .01; % mask expand step parameter
            sign_params.min_area_size = 250;  % min area size in microns^2
            sign_params.final_erode = 10;
            sign_params.final_dilate = 10;
            
            sign_params = ne7.mat.getParams(sign_params,varargin);
            
            % get signmap 
            sign_map = fetch1(self,'sign_map');
            grad_diff = sign_map(:,:,1);
            A = sign_map(:,:,2);
            
            % convert size threshold to pixels^2
            sign_params.min_area_size = (sign_params.min_area_size/...
                mean(fetchn(map.OptImageBar & (map.RetMapScan & self),'pxpitch')))^2;
 
            % filter & expand masks
            filtered_grad_diff = zeros(size(grad_diff));
            imStd = std(grad_diff(:));
            filtered_grad_diff(grad_diff>imStd*sign_params.sign_thr)= 1;
            filtered_grad_diff(grad_diff<-imStd*sign_params.sign_thr)= -1;
            SE = strel('arbitrary',[0 1 0;1 1 1;0 1 0]); %%% parameter
            SEMAX = strel('arbitrary',[0 0 1 0 0;0 1 1 1 0;1 1 1 1 1;0 1 1 1 0;0 0 1 0 0]); %%% parameter
            [all_areas, n] = bwlabel(filtered_grad_diff);
            for iThr = 1:-sign_params.step:0.05
                n = max(all_areas(:));
                thr = imStd*iThr;
                for iArea = 1:n
                    one_area = imfill(imdilate(all_areas==iArea,SE),'holes');
                    all_areas(one_area & all_areas==0 & abs(grad_diff)>thr) = iArea;
                end
                lmax = imdilate(imbinarize(all_areas),SEMAX);
                one_area = logical(abs(grad_diff)>thr & lmax == 0);
                stats = regionprops(one_area,'area');
                un = unique(one_area(:));
                for i = un(un>0)'
                    if stats(i).Area<sign_params.min_area_size
                        one_area(one_area==i) = 0;
                    end
                end
                one_area = bwlabel(one_area);
                if any(one_area(:))
                    all_areas(logical(one_area)) = one_area(logical(one_area))+n;
                end
            end
            
            % select areas that have at least half their area with amplitude more than the threshold
            area_map = zeros(size(all_areas));
            idx = 0;
            thr = mean(A(:));% max(A(:)) - 2 * std(A(:)); %%% parameter
            for iarea = 1:n
                if ~(sum(all_areas(:)==iarea)/sum(all_areas(:)==iarea & A(:)>thr)>2)
                    idx = idx+1;
                    mask = imdilate(imerode(all_areas==iarea,strel('disk',sign_params.final_erode)),strel('disk',sign_params.final_dilate));
                    area_map(mask) = idx;
                end
            end
            
            if ~nargout
                % plot unique areas
                im(:,:,1) = ne7.mat.normalize(area_map);
                im(:,:,2) = area_map>0;
                im(:,:,3) = sign_map(:,:,3);
                image(hsv2rgb(im)); axis image; axis off; title('Areas')
            else
                masks = area_map;
            end
        end
        
        function im = replaceNaNs(~,im)
            idx = find(isnan(im));
            while ~isempty(idx)
                [x,y] = ind2sub(size(im),idx);
                im(idx) = interp2(im,x,y,'spline',0);
                idx = find(isnan(im));
            end
        end
    end
    
end

