function out=refitSBE63calibrationcertificate(in,modeltype)
% function out=refitSBE63calibrationcertificate(in,modeltype)
%
% refit original Sea-Bird multipoint calibration certificates to match,
% e.g., the GEOMAR modification of the Uchida et al. 2008 model after
% Bittig et al. 2017
%
% Note that the Sea-Bird multipoint calibration's O2 reference data are
% given im mL(S.T.P.) O2 L-1.
%
% inputs: structure 'in' with fields
% temp            - 24x1 matrix of temperatures (in °C)
% phasedelay      - 24x1 matrix of phasedelay readings (in us)
% O2ref           - 24x1 matrix of O2 reference (in mL L-1)
% serial_number   - optode serial number
% batch_number (optional)         - note on calibration
% calibrationdate (optional)      - date of calibration certificate
%
%
% modeltype (optional) - functional model for refit 
%                        ('SVU' or 'Bittig' (default))
% 
%
% Henry Bittig, LOV/GEOMAR
% 12.05.2017
%
% reference:
% Bittig et al. (2017). Oxygen Optode Sensors: Principle, Characterization,
% Calibration and Application in the Ocean. Front. Mar. Sci.
% http://dx.doi.org/10.3389/fmars.2017.00429 
%
% requires Matlab O2 conversion functions following SCOR WG 142
% recommendations: 
% Bittig et al. (2016). SCOR WG 142: Quality Control Procedures for Oxygen
% and Other Biogeochemical Sensors on Floats and Gliders. Recommendations
% on the conversion between oxygen quantities for Bio-Argo floats and other
% autonomous sensor platforms. http://dx.doi.org/10.13155/45915

% usage example:
% out=refitSBE63calibrationcertificate(sbe63_0800_factory);
% out=refitSBE63calibrationcertificate(sbe63_0800_factory,'Bittig');

%% check inputs
if ~isempty(setdiff({'temp';'phasedelay';'O2ref';'serial_number'},fieldnames(in)));
    disp(['     Mandatory field(s) ' strjoin(setdiff({'temp';'phasedelay';'O2ref';'serial_number'},fieldnames(in)),', ') ' are missing!'])
    return
end
if nargin<2, modeltype='bittig'; else modeltype=lower(modeltype); end

%% store metadata
if isfield(in,'batch_number'), out.batch_number=in.batch_number; batchno=in.batch_number;
else batchno=''; end
out.model_number='SBE63';
out.serial_number=in.serial_number;
if isfield(in,'calibrationdate'), out.calibrationdate=in.calibrationdate; end
out.refitdate=datestr(now);
out.modeltype=modeltype;

%% start to reorganize input data
phase=in.phasedelay(:); % / us
temp=in.temp(:); % / °C
O2conc=in.O2ref(:).*44.6596; % / convert to umol L-1; freshwater
pO2=O2ctoO2p(O2conc,temp,0,0); % and convert to pO2 / hPa
% and keep data in output
out.temp=temp(:);
out.phasedelay=phase(:); % / us
out.phase=phase(:)*1e-6*3840*360; % convert phase delay / us to phase shift / °
out.pO2=pO2(:);
out.O2conc=O2conc(:);

%% refit 24x1 foil calibration data according to specified modeltype 
% (individual multi-point data)
options=statset('MaxIter',1200); % optimization options
if strcmp(modeltype,'bittig')
    option=1; % fit pO2 / hPa
    % initial guess
    %foilcoef0=[5e-3 8e-5 1e-7 1e-1 -4e-5 -1e-2 1e-3]'; % 7 parameter Bittig model
    foilcoef0=[2e-3 8e-5 1e-7 1e-1 -4e-5 -1e-2 1e-3]'; % 7 parameter Bittig model
    % do optimisation (of pO2)
    foilcoef=nlinfit([temp phase],pO2,@optodefunUchidaBittig,foilcoef0,options);
    % then statistics for actual pO2 fit
    pO2fit=optodefunUchidaBittig(foilcoef,[temp phase]); % / hPa
    O2concfit=O2ptoO2c(pO2fit,temp,0,0); % and converted back to O2conc / umol L-1
elseif strcmp(modeltype,'svu')
    option=2; % fit O2concentration / umol L-1
    % initial guess
    %foilcoef0=[8e-3 -1e-4 -1e-6 0 0 1/70 0]'; % 7 parameter Uchida model
    foilcoef0=[5e-3 8e-5 1e-7 1e-1 -4e-5 -1e-2 1e-3]'; % 7 parameter Uchida model
    % do optimisation (of freshwater concentration)
    foilcoef=nlinfit([temp phase],O2conc,@optodefunUchidaSVU,foilcoef0,options);
    % then statistics for (freshwater) oxygen concentration / umol L-1
    O2concfit=optodefunUchidaSVU(foilcoef,[temp phase]); % / umol L-1
else
    disp('Functional model not recognized. Please use ''Bittig'' for Bittig et al. 2017 model (pO2) or ''SVU'' for Uchida et al. 2008 model (O2concentration).')
    return
end
out.foilcoef=foilcoef;
out.r=O2concfit-O2conc; % residuals / umol L-1
out.rmse=sqrt(sum(out.r.^2)./length(out.r)); % rmse / umol L-1

%%{
%% add visual check
note=char({'SBE63 multi-point certificate refit',['' num2str(batchno,'%.4u') ' -> SBE63 SN ' num2str(in.serial_number,'%.4u')]});
% prepare plot of results
tlim=[min(temp)-2 max(temp)+2]; %freshwater; only T and Phase
if option==12
    plim=[min(phase)-2e-6 max(phase)+2e-6];
else
    plim=[min(phase)-2 max(phase)+2];
end
% fitted surface
[plotT,plotP]=meshgrid([tlim(1):(tlim(2)-tlim(1))/15:tlim(2)],[plim(1):(plim(2)-plim(1))/15:plim(2)]);
if option==1
    plotO=O2ptoO2c(optodefunUchidaBittig(foilcoef(:,1),[plotT(:),plotP(:)]),plotT(:),0,0);
elseif option==2
    plotO=optodefunUchidaSVU(foilcoef(:,1),[plotT(:),plotP(:)]);
end
plotO=reshape(plotO,size(plotT));
plotO(plotO>max(O2conc+60))=NaN;
plotO(plotO<min(O2conc-60))=NaN;

% plot surface
fha=figure;
set(fha,'color',[1 1 1])
set(fha,'Position',[068 350 (1.7+0.6)*560 1*420])
subplot(1,4,4) % statistics/info plot
subplot(1,4,3) % difference plot
subplot(1,4,[1 2]) % surface plot

brightness=[.5 .5 .5];
mhandle=mesh(plotP,plotT,plotO,'EdgeColor',brightness); % fitted surface
set(gca,'XDir','reverse','YDir','reverse')
xlim(plim)
ylim(tlim)
hold on
shandle=scatter3(phase,temp,O2conc,72,[0 0 1],'filled','MarkerEdgeColor',[0 0 0]); % reference samples
set(shandle,'CData',out.r); % rednblue
for i = 1:length(O2conc)    % and sample-fit surface distance line for each
    plot3([phase(i) phase(i)],[temp(i) temp(i)],[O2conc(i) O2concfit(i)],'-','Linewidth',1,'Color',brightness)
end
hold off
view(-75,30)
xlabel('TCPhase / °')
ylabel('T / °C')
zlabel('O_2 / \mumol L^{-1}')
lhandle=legend({'fitted surface';'reference samples'},'Location','NorthEast');
hidden off
cbah=colorbar;
set(get(cbah,'Title'),'String','\Delta / \mumol L^{-1}');
caxis(get(shandle,'Parent'),[-5 5])

% add plain difference plot
subplot(1,4,3)
scah2(1)=scatter(O2conc,out.r,72,out.r,'filled','MarkerEdgeColor',[0 0 0]);
grid on
xlabel('O_2 / \mumol L^{-1}')
ylabel('\Delta = optode fit - reference / \mumol L^{-1}')
caxis([-5 5])
ylim([-5 5])

% info plot
subplot(1,4,4)
if option==1
    functionstring={'Bittig et al. 2017 model (pO2)'};
elseif option==2
    functionstring={'Uchida et al. 2008 model (O2 conc)'};
end
betastring{1}='new coefficients:';
for i=1:size(foilcoef,1)
	eval(['betastring{i+1}= ''c'  num2str(i) ' = ' num2str(foilcoef(i,1),'%6.4e') ''';']);
end
set(gca,'Visible','Off')
thandle=text(0,.5,[note, functionstring ,{'','',['RMSE = ' num2str(out.rmse,'%6.2f') ' \mumol L^{-1}  (N = ' num2str(length(find(~isnan(O2concfit)))) ')'],''}, betastring ]);

try cmocean('balance'), catch me, end % try to set diverging rednblue colormap from cmocean

% finally set renderer of figure to painters (not zbuffer)
set(fha,'Renderer','painters')
%%}