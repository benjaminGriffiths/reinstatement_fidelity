function varargout = csl_getRes(sCon)
	

% xSPM      - structure containing SPM, distribution & filtering details
% .swd      - SPM working directory - directory containing current SPM.mat
% .title    - title for comparison (string)
% .Z        - minimum of n Statistics {filtered on u and k}
% .n        - number of conjoint tests        
% .STAT     - distribution {Z, T, X, F or P}     
% .df       - degrees of freedom [df{interest}, df{residual}]
% .STATstr  - description string     
% .Ic       - indices of contrasts (in SPM.xCon)
% .Im       - indices of masking contrasts (in xCon)
% .pm       - p-value for masking (uncorrected)
% .Ex       - flag for exclusive or inclusive masking
% .u        - height threshold
% .k        - extent threshold {voxels}
% .XYZ      - location of voxels {voxel coords}
% .XYZmm    - location of voxels {mm}
% .S        - search Volume {voxels}
% .R        - search Volume {resels}
% .FWHM     - smoothness {voxels}     
% .M        - voxels -> mm matrix
% .iM       - mm -> voxels matrix
% .VOX      - voxel dimensions {mm} - column vector
% .DIM      - image dimensions {voxels} - column vector
% .Vspm     - Mapped statistic image(s)
% .Ps       - list of P values for voxels at SPM.xVol.XYZ (used by FDR)
%- sCon = struct( ...
%-     'spmmat', '', ...    % full path to SPM.mat file
%-     'Ic', [], ...        % no of contrast (or contrasts for conjunction)
%-     'mask', 0,...        % whether to mask with another contrast
%-     'masking_Ic', [],... % no of contrast to mask with
%-     'mask_ucp', [],...   % masking contrast uncorrected p
%-     'exclusive', [],...  % whether masking is inclusive or exclusive
%-     'title', '',...      % if empty results in default contrast title
%-     'thresh', 0.01,...  % threshold (corrected or uncorrected, as above)
%-     'Mcp', 'none',...  % Mutiple comp : FWE|FDR|none
%-     'extent_thresh', 0); % extent threshold
%- 
%- 
%- This function is just a fragment  of spm_results_ui.m, with a liitle change a line 
%- to call csl_getSPM2 with sCon as input.
%====================================================================
%- From spm_results_ui 
%- KHERIF Ferath csl, Cambridge, UK.

%-Initialise 
%-----------------------------------------------------------------------
%SPMid      = spm('FnBanner',mfilename,SCCSid);
[Finter,Fgraph,CmdLine] = spm('FnUIsetup','Stats: Results');
FS         = spm('FontSizes');

% clear satfig if it exists
%-----------------------------------------------------------------------
hSat       = findobj('tag','Satellite');
spm_figure('clear',hSat);

%-Get thresholded xSPM data and parameters of design
%=======================================================================
[SPM,xSPM] = csl_getSPM2(sCon);
M          = SPM.xVol.M;
DIM        = SPM.xVol.DIM;

% ensure pwd = swd so that relative filenames are valid
%-----------------------------------------------------------------------
cd(SPM.swd)

%-Setup Results User Interface; Display MIP, design matrix & parameters
%=======================================================================
spm('FigName',['SPM{',xSPM.STAT,'}: Results'],Finter,CmdLine);


%-Setup results GUI
%-----------------------------------------------------------------------
spm_figure('Clear',Finter)
hReg      = spm_results_ui('SetupGUI',M,DIM,xSPM,Finter);

%-Setup design interrogation menu
%-----------------------------------------------------------------------
hDesRepUI = spm_DesRep('DesRepUI',SPM);
figure(Finter)

%-Setup Maximium intensity projection (MIP) & register
%-----------------------------------------------------------------------
hMIPax = axes('Parent',Fgraph,'Position',[0.05 0.60 0.55 0.36],'Visible','off');
hMIPax = spm_mip_ui(xSPM.Z,xSPM.XYZmm,M,DIM,hMIPax);
spm_XYZreg('XReg',hReg,hMIPax,'spm_mip_ui');
if xSPM.STAT == 'P'
	str = xSPM.STATstr;
else
	str = ['SPM\{',xSPM.STATstr,'\}'];
end
text(240,260,str,...
	'Interpreter','TeX',...
	'FontSize',FS(14),'Fontweight','Bold',...
	'Parent',hMIPax)


%-Print comparison title
%-----------------------------------------------------------------------
hTitAx = axes('Parent',Fgraph,...
		'Position',[0.02 0.95 0.96 0.02],...
		'Visible','off');

text(0.5,0,xSPM.title,'Parent',hTitAx,...
	'HorizontalAlignment','center',...
	'VerticalAlignment','baseline',...
	'FontWeight','Bold','FontSize',FS(14))


%-Print SPMresults: Results directory & thresholding info
%-----------------------------------------------------------------------
hResAx = axes('Parent',Fgraph,...
		'Position',[0.05 0.55 0.45 0.05],...
		'DefaultTextVerticalAlignment','baseline',...
		'DefaultTextFontSize',FS(9),...
		'DefaultTextColor',[1,1,1]*.7,...
		'Units','points',...
		'Visible','off');
AxPos = get(hResAx,'Position'); set(hResAx,'YLim',[0,AxPos(4)])
h     = text(0,24,'SPMresults:','Parent',hResAx,...
	'FontWeight','Bold','FontSize',FS(14));
text(get(h,'Extent')*[0;0;1;0],24,spm_str_manip(SPM.swd,'a30'),'Parent',hResAx)
text(0,12,sprintf('Height threshold %c = %0.2f',xSPM.STAT,xSPM.u),'Parent',hResAx)
text(0,00,sprintf('Extent threshold k = %0.0f voxels',xSPM.k), 'Parent',hResAx)


%-Plot design matrix
%-----------------------------------------------------------------------
hDesMtx   = axes('Parent',Fgraph,'Position',[0.65 0.55 0.25 0.25]);
hDesMtxIm = image((SPM.xX.nKX + 1)*32);
xlabel('Design matrix')
set(hDesMtxIm,'ButtonDownFcn','spm_DesRep(''SurfDesMtx_CB'')',...
	'UserData',struct(...
		'X',		SPM.xX.xKXs.X,...
		'fnames',	{reshape({SPM.xY.VY.fname},size(SPM.xY.VY))},...
		'Xnames',	{SPM.xX.name}))

%-Plot contrasts
%-----------------------------------------------------------------------
nPar   = size(SPM.xX.X,2);
xx     = [repmat([0:nPar-1],2,1);repmat([1:nPar],2,1)];
nCon   = length(xSPM.Ic);
xCon   = SPM.xCon;
if nCon
	dy     = 0.15/max(nCon,2);
	hConAx = axes('Position',[0.65 (0.80 + dy*.1) 0.25 dy*(nCon-.1)],...
		'Tag','ConGrphAx','Visible','off');
	title('contrast(s)')
	htxt   = get(hConAx,'title'); 
	set(htxt,'Visible','on','HandleVisibility','on')
end

for ii = nCon:-1:1
    axes('Position',[0.65 (0.80 + dy*(nCon - ii +.1)) 0.25 dy*.9])
    if xCon(xSPM.Ic(ii)).STAT == 'T' & size(xCon(xSPM.Ic(ii)).c,2) == 1

	%-Single vector contrast for SPM{t} - bar
	%---------------------------------------------------------------
	yy = [zeros(1,nPar);repmat(xCon(xSPM.Ic(ii)).c',2,1);zeros(1,nPar)];
	h  = patch(xx,yy,[1,1,1]*.5);
	set(gca,'Tag','ConGrphAx',...
		'Box','off','TickDir','out',...
		'XTick',spm_DesRep('ScanTick',nPar,10) - 0.5,'XTickLabel','',...
		'XLim',	[0,nPar],...
		'YTick',[-1,0,+1],'YTickLabel','',...
		'YLim',[min(xCon(xSPM.Ic(ii)).c),max(xCon(xSPM.Ic(ii)).c)] +...
		       [-1 +1] * max(abs(xCon(xSPM.Ic(ii)).c))/10	)

    else

	%-F-contrast - image
	%---------------------------------------------------------------
	h = image((xCon(xSPM.Ic(ii)).c'/max(abs(xCon(xSPM.Ic(ii)).c(:)))+1)*32);
	set(gca,'Tag','ConGrphAx',...
		'Box','on','TickDir','out',...
		'XTick',spm_DesRep('ScanTick',nPar,10),'XTickLabel','',...
		'XLim',	[0,nPar]+0.5,...
		'YTick',[0:size(SPM.xCon(xSPM.Ic(ii)).c,2)]+0.5,....
		'YTickLabel','',...
		'YLim',	[0,size(xCon(xSPM.Ic(ii)).c,2)]+0.5	)

    end
    ylabel(num2str(xSPM.Ic(ii)))
    set(h,'ButtonDownFcn','spm_DesRep(''SurfCon_CB'')',...
    	'UserData',	struct(	'i',		xSPM.Ic(ii),...
    				'h',		htxt,...
    				'xCon',		xCon(xSPM.Ic(ii))))
end


%-Store handles of results section Graphics window objects
%-----------------------------------------------------------------------
H  = get(Fgraph,'Children');
H  = findobj(H,'flat','HandleVisibility','on');
H  = findobj(H);
Hv = get(H,'Visible');
set(hResAx,'Tag','PermRes','UserData',struct('H',H,'Hv',{Hv}))


%-Finished results setup
%-----------------------------------------------------------------------
varargout = {hReg,xSPM,SPM};
