function out=refit4330batchcalibrationcertificate(in,modeltype,refittype)
% function out=refit4330batchcalibrationcertificate(in,modeltype,refittype)
%
% refit original Aanderaa batch calibration certificates to match, e.g.,
% the Uchida et al. 2008 model ("SVU") or the GEOMAR modification after
% Bittig et al. 2017
%
% inputs: structure 'in' with fields
% temp            - 63x1 matrix of temperatures (in °C)
% phase           - 63x1 matrix of phase readings (in °)
% O2ref           - 63x1 matrix of O2 reference (in umol L-1)
% airsaturated    - 1x3 matrix of phase, temperature, and air pressure for
%                   air saturated water calibration point (2-point
%                   calibration)
% zerosolution    - 1x2 matrix of phase and temperature for zero solution
%                   (2-point calibration)
% serial_number   - optode serial number
% batch_number (optional)         - number of foil batch
% batchcalibrationdate (optional) - date of batch calibration
% calibrationdate (optional)      - date of calibration certificate's /
%                                   2-point calibration
%
%
% modeltype (optional) - functional model for refit 
%                        ('SVU' or 'Bittig' (default))
% refittype (optional) - type of refit equation (see reference for details):
%                        'a' (not recommended), 'b' or 'e' (default)
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
% out=refit4330batchcalibrationcertificate(aadi4330_1463_b1206E);
% out=refit4330batchcalibrationcertificate(aadi4330_1463_b1206E,'Bittig');
% out=refit4330batchcalibrationcertificate(aadi4330_1463_b1206E,[],'c');

%% check inputs
if ~isempty(setdiff({'temp';'phase';'O2ref';'airsaturated';'zerosolution';'serial_number'},fieldnames(in)));
    disp(['     Mandatory field(s) ' strjoin(setdiff({'temp';'phase';'O2ref';'airsaturated';'zerosolution';'serial_number'},fieldnames(in)),', ') ' are missing!'])
    return
end
if nargin<2 || isempty(modeltype), modeltype='bittig'; else modeltype=lower(modeltype); end
if nargin<3 || isempty(refittype), refittype='e'; else refittype=lower(refittype); end


%% store metadata
if isfield(in,'batch_number'), out.batch_number=in.batch_number; batchno=in.batch_number;
else batchno=''; end
if isfield(in,'batchcalibrationdate'), out.batchcalibrationdate=in.batchcalibrationdate; end
out.model_number=4330;
out.serial_number=in.serial_number;
if isfield(in,'calibrationdate'), out.calibrationdate=in.calibrationdate; end
out.refitdate=datestr(now);
out.modeltype=modeltype;

%% start to reorganize input data
phase=in.phase(:); % / °
temp=in.temp(:); % / °C
O2conc=in.O2ref(:); % / umol L-1; freshwater
pO2=O2ctoO2p(O2conc,temp,0,0); % and convert to pO2 / hPa

%% refit 63x1 foil calibration data according to specified modeltype 
% (batch data)
options=statset('MaxIter',1200); % optimization options
if strcmp(modeltype,'bittig')
    option=1; % fit pO2 / hPa
    % initial guess
    %foilcoef0=[5e-3 8e-5 1e-7 1e-1 -4e-5 -1e-2 1e-3]'; % 7 parameter Bittig model
    foilcoef0=[2e-3 8e-5 1e-7 1e-1 -4e-5 -1e-2 1e-3]'; % 7 parameter Bittig model
    % do optimisation (of pO2)
    foilcoef=nlinfit([temp phase],pO2,@optodefunUchidaBittig,foilcoef0,options);
elseif strcmp(modeltype,'svu')
    option=2; % fit O2concentration / umol L-1
    % initial guess
    %foilcoef0=[8e-3 -1e-4 -1e-6 0 0 1/70 0]'; % 7 parameter Uchida model
    foilcoef0=[5e-3 8e-5 1e-7 1e-1 -4e-5 -1e-2 1e-3]'; % 7 parameter Uchida model
    % do optimisation (of freshwater concentration)
    foilcoef=nlinfit([temp phase],O2conc,@optodefunUchidaSVU,foilcoef0,options);
else
    disp('Functional model not recognized. Please use ''Bittig'' for Bittig et al. 2017 model (pO2) or ''SVU'' for Uchida et al. 2008 model (O2concentration).')
    return
end

%% add individual optode 2-point calibration data
% get equilibrium pO2 partial pressure for 100 % air saturated water
pO2air=O2stoO2p(100,in.airsaturated(2),0,0,in.airsaturated(3)); % / hPa

%% adjustments for two-point calibration:
switch refittype
    case 'a'
%% follow Bittig et al. 2017: refit a 
% traditional Aanderaa method: slope and offset on phase (PhaseCoefs)
% estimated accuracy: 18 hPa (=legacy processing; not recommended)
out.estimated_accuracy_hPa=18; % / hPa

% slope on phase, offset on phase
% get scaling of batch phase for current optode (based on two-point calibration)
if option==1
    slopeoffset=nlinfit([in.zerosolution(2) in.zerosolution(1);in.airsaturated(2) in.airsaturated(1)],[0;pO2air],@(beta,X)optodefunUchidaBittig(foilcoef,[X(:,1) X(:,2).*beta(1)+beta(2)]),[1;0]);
elseif option==2
    slopeoffset=nlinfit([in.zerosolution(2) in.zerosolution(1);in.airsaturated(2) in.airsaturated(1)],[0;O2ptoO2c(pO2air,in.airsaturated(2),0,0)],@(beta,X)optodefunUchidaSVU(foilcoef,[X(:,1) X(:,2).*beta(1)+beta(2)]),[1;0]);
end
% scale reference phase
phase=(phase-slopeoffset(2))./slopeoffset(1); % / °

    case 'b'
%% follow Bittig et al. 2017: refit b 
% estimated accuracy: 9 hPa
out.estimated_accuracy_hPa=9; % / hPa

% slope on pO2/O2, offset on phase
% get scaling of batch equilibrium pO2 for current optode (based on two-point calibration)
% start with phase offset at zero O2
if option==1
    offset=fzero(@(beta)optodefunUchidaBittig(foilcoef,[in.zerosolution(2) in.zerosolution(1)+beta]),0); % / °
elseif option==2
    offset=fzero(@(beta)optodefunUchidaSVU(foilcoef,[in.zerosolution(2) in.zerosolution(1)+beta]),0); % / °
end
% scale reference phase for zero O2 phase offset
phase=phase-offset; % / °
% and get updated foil batch coefficients (including zero O2 phase offset)
if option==1
    % do optimisation (of pO2)
    foilcoef=nlinfit([temp phase],pO2,@optodefunUchidaBittig,foilcoef0,options);
elseif option==2
    % do optimisation (of freshwater concentration)
    foilcoef=nlinfit([temp phase],O2conc,@optodefunUchidaSVU,foilcoef0,options);
end
% get pO2 partial pressure for 100 % air saturated water optode measurement
% according to (zero-point updated) foil batch coefficients
if option==1
    % gives pO2 directly
    pO2air_batch=optodefunUchidaBittig(foilcoef,[in.airsaturated(2),in.airsaturated(1)]); % / hPa
elseif option==2
    % needs to be converted to pO2
    pO2air_batch=O2ctoO2p(optodefunUchidaSVU(foilcoef,[in.airsaturated(2),in.airsaturated(1)]),in.airsaturated(2),0,0); % / hPa
end
slope=pO2air./pO2air_batch;

% scale batch foil data for 100 % calibration
pO2=pO2.*slope;
O2conc=O2ptoO2c(pO2,temp,0,0); % / umol L-1; freshwater

    case 'e'
%% follow Bittig et al. 2017: refit e
% estimated accuracy: 7 hPa
out.estimated_accuracy_hPa=7; % / hPa

% offset on pO2/O2, offset on phase
% get scaling of batch phase for current optode (based on two-point calibration)
if option==1
    offsetoffset=nlinfit([in.zerosolution(2) in.zerosolution(1);in.airsaturated(2) in.airsaturated(1)],[0;pO2air],@(beta,X)optodefunUchidaBittig(foilcoef,[X(:,1) X(:,2)+beta(1)])-beta(2),[0;0]);
    % scale batch foil data for 100 % calibration
    pO2=pO2-offsetoffset(2); % / hPa
    O2conc=O2ptoO2c(pO2,temp,0,0); % / umol L-1; freshwater
elseif option==2
    offsetoffset=nlinfit([in.zerosolution(2) in.zerosolution(1);in.airsaturated(2) in.airsaturated(1)],[0;pO2air],@(beta,X)O2ctoO2p(optodefunUchidaSVU(foilcoef,[X(:,1) X(:,2)+beta(1)]),X(:,1),0,0)-beta(2),[0;0]);
    % scale batch foil data for 100 % calibration
    pO2=pO2-offsetoffset(2); % / hPa
    O2conc=O2ptoO2c(pO2,temp,0,0); % / umol L-1; freshwater
end
% scale reference phase
phase=phase-offsetoffset(1); % / °

    otherwise
disp(['Did not recognize refittype ''' refittype '''. No 2-point adjustment of batch calibration coefficients performed.'])
return
end


%% recalculate foil coefficients based on two-point adjusted reference data
% batch foil data adjusted to individual optode
out.temp=temp(:);
out.phase=phase(:);
out.pO2=pO2(:);
out.O2conc=O2conc(:);
% and re-do fit for single optode calibration
if strcmp(modeltype,'bittig')
    % do optimisation (of pO2)
    foilcoef=nlinfit([temp phase],pO2,@optodefunUchidaBittig,foilcoef0,options);
    % then statistics for actual pO2 fit
    pO2fit=optodefunUchidaBittig(foilcoef,[temp phase]); % / hPa
    O2concfit=O2ptoO2c(pO2fit,temp,0,0); % and converted back to O2conc / umol L-1
elseif strcmp(modeltype,'svu')
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

%% check result
if option==1
    checkvalue=optodefunUchidaBittig(foilcoef,[in.zerosolution(2) in.zerosolution(1);in.airsaturated(2) in.airsaturated(1)])-[0;pO2air]; % / hPa
elseif option==2
    checkvalue=O2ctoO2p(optodefunUchidaSVU(foilcoef,[in.zerosolution(2) in.zerosolution(1);in.airsaturated(2) in.airsaturated(1)]),[in.zerosolution(2);in.airsaturated(2)],0,0)-[0;pO2air]; % / hPa
end
disp(['Refit residuals: ' num2str(checkvalue(1),'%+.2f') ' hPa at 0 % O2 / ' num2str(in.zerosolution(2),'%.1f') ' °C and ' num2str(checkvalue(2),'%+.2f') ' hPa at 100 % O2 / ' num2str(in.airsaturated(2),'%.1f') ' °C'])

%%{
%% add visual check
note=char({'Aanderaa batch certificate refit',['batch no. ' num2str(batchno,'%.4u') ' -> 4330 SN ' num2str(in.serial_number,'%.4u')]});
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