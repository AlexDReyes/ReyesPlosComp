function AudCTX_analysis(baseName) %e.g. basename of "myFile1" is "myFile"

cThold=-0.27;%minimum synaptic current to evoke action potentials; used to calculate diameter of effective synaptic current
  fname = [baseName + num2str(1)];
  load(fname,'ncells','R_length','delT','repetition');
  xArray=zeros(1,repetition);dArray=xArray;fdArray=xArray;xCMarray=xArray;yCMarray=xArray;synArray=xArray;totLarray=xArray;totLarray_syn=xArray;
    histoSum=zeros(1,round(sqrt(2*ncells^2)));histoSum_syn=zeros(1,round(sqrt(2*ncells^2)));
    for iii=1:repetition
     fname = [baseName + num2str(iii)];
        [xArray(iii),fdArray(iii),synArray(iii),totLarray(iii),totLarray_syn(iii)]=seeCMPmap(fname);%figure(20);view([-25 40]);figure(21);view([-25 40])
    end
    meanArea=mean(xArray(1:repetition));  % = mean number of evoked spikes
        stdArea=std(xArray(1:repetition));
    meanfD=mean(fdArray(1:repetition)); % diameter of fitted circle to spikes (black trace in Fig. 3)
        stdfD=std(fdArray(1:repetition));
    meanSdiam=mean(synArray(1:repetition)); % diameter of fitted circle to net synaptic current that exceeded rheobase (cThold above)
        stdSdiam=std(synArray(1:repetition));
    meanL=mean(totLarray(1:repetition));   %mean length of projection of spikes to tonotopic axis
        stdL=std(totLarray(1:repetition));
    meanL_syn=mean(totLarray_syn(1:repetition)); %mean length of projection of synaptic input to tonotopic axis
        stdL_syn=std(totLarray_syn(1:repetition));
        
    fName="AnalysisFile"; %file for output
    save(fName,"meanArea","stdArea","meanfD","stdfD","meanSdiam","stdSdiam","meanL","stdL","meanL_syn","stdL_syn");

function [area,fDiam,synDiam,totL_spk,totL_syn]=seeCMPmap(fname)

center=round(ncells/2);
ms = delT*(1:R_length);
nxn=ncells^2;
mgEtoP = zeros(ncells,ncells);mgItoP= zeros(ncells,ncells);mFreqE=zeros(ncells,ncells);mFreqI=mFreqE;%mFreqE=NaN(ncells,ncells);mFreqI=mFreqE;
load(fname,'cmpE','cmpI','cmpgEtoP','cmpgItoP');
get_Conts(); %used to calculate maximum conductance evoked in E cells (Figs 1, 2)

%--------------------------------
synContour=mgetgEIprofiles();
[dum1,synDiam,dum3,dum4]=decode(synContour,1);
synContour(find(synContour~=1))=0;
[totL_syn,histL_syn]=calc_Proj(synContour);
[area,fDiam,xCM,yCM]=findSpikes(cmpE,cmpI,1,R_length);
[totL_spk,histL_spk]=calc_Proj(mFreqE);
clear cmpE cmpI cmpgEtoP cmpgItoP

    function [tL,hist]=calc_Proj(S)
        vDiag=[ncells ncells];vLength=norm(vDiag);hist=zeros(1,round(vLength));
        for ii=1:ncells
            for jj=1:ncells
                if (S(ii,jj)>0)
                    vPoint=[ii jj];
                    proj=dotprod(vDiag,vPoint')/vLength;
                    hist(round(proj))=hist(round(proj))+1;
                end
            end
        end
        tHist=hist;
        tHist(find(hist>0))=1;
        tL=sum(tHist);
    end


function [area,fDiam,xCM,yCM]=decode(S,showPlot)
    area=sum(sum(S));
   "area:" + area;
    [fDiam,xCM,yCM]=findCoords();

    function [fDiam,xCM,yCM]=findCoords()
        xCoord=zeros(1,2);yCoord=xCoord;
        nZ=find(S==1);
        for i=1:numel(nZ)
            yCoord(i)= mod(fix((nZ(i)-1)),ncells)+1;  %finds x and y values of activated cells  ***COLUMN
            xCoord(i) = floor(fix((nZ(i)-1)/ncells))+2;
        end
        xCoord(find(xCoord==0))=[];
        yCoord(find(yCoord==0))=[];
        xCM = sum(xCoord)/numel(xCoord);
        yCM = sum(yCoord)/numel(yCoord);

        k = boundary(xCoord',yCoord');

        th = linspace(0,2*pi,20)';R=1.1111111111;
        [xc,yc,Re] =circfit(xCoord(k),yCoord(k));
        fDiam=2*Re;
        "diameter:" + fDiam;
        xe = Re*cos(th)+xc; ye = Re*sin(th)+yc;
        if (showPlot==2)
            figure(3);clf;hold on; title('Network activity');ylabel('cell #');ax=gca; ax.FontSize=14;
            plot(xCoord,yCoord,'.');axis square;axis([0 ncells 0 ncells]);
            plot(xCoord(k),yCoord(k));
            plot(xe,ye,'k'); legend('active cells','outermost boundary','fitted circle');

        end
            function   [xc,yc,R,a] = circfit(x,y)  %fits circle to outermost points of active cells
               x=x(:); y=y(:);
               a=[x y ones(size(x))]\[-(x.^2+y.^2)];
               xc = -.5*a(1);
               yc = -.5*a(2);
               R  =  sqrt((a(1)^2+a(2)^2)/4-a(3));
            end
    end
 
end


  function [area,fDiam,xCM,yCM]=findSpikes(mFileE,mFileI,startt,endd)
        for i=1:nxn
                cellPosX = mod(fix((i-1)),ncells)+1;  %finds x and y values of activated inhcells  ***COLUMN
                cellPosY = floor(fix((i-1)/ncells));
            if ((cellPosX>0)&&(cellPosY>0))
                if (find(mFileE(i,startt:endd))>0)
                    mFreqE(cellPosX,cellPosY)=1;
                end
                if (find(mFileI(i,startt:endd))>0)
                    mFreqI(cellPosX,cellPosY)=1;
                end
           end
        end
        
      [nSpikes,fDiam,xCM,yCM]=decode(mFreqE,2);
      area=nSpikes;
  end

    function get_Conts()
    
    for i=1:nxn
        cellPosX = mod(fix((i-1)),ncells)+1;  %finds x and y values of activated inhcells  ***COLUMN
        cellPosY = floor(fix((i-1)/ncells));%mod((j-1),ncells)+1; %***ROW
        if ((cellPosX>0)&&(cellPosY>0))
            mgEtoP(cellPosX,cellPosY) = max(cmpgEtoP(i,:))/repetition;
            mgItoP(cellPosX,cellPosY) = max(cmpgItoP(i,:))/repetition;
        end
    end
    tM=thin_it(mgEtoP);
    figure(1);clf;mesh(tM,'facecolor','none','EdgeColor','k');title('Excitatory synaptic conductance to E cells');xlabel('cell #');ylabel('cell #');zlabel('maximum conductance');ax=gca; ax.FontSize=14;
    c=tM;
    tM=thin_it(mgItoP);figure(2);clf;mesh(tM,'facecolor','none','EdgeColor','k');title('Inhibitory synaptic conductance to E cells');xlabel('cell #');ylabel('cell #');zlabel('maximum conductance');ax=gca; ax.FontSize=14;
    clear tM;clear c;
    clear mgEtoP;clear mgItoP;
    end

function tFile=thin_it(mFile)
tFile = mFile(1:1:ncells,1:1:ncells);
end

function synCont=mgetgEIprofiles()
    synCont=zeros(ncells,ncells);synContE=synCont;synContI=synContE;synContIb=synContE;delE=zeros(1,R_length);
    sweeps=0;
    PSTH=getgEIprofiles(center,cmpgEtoP);
    sweeps=0;
    Espread=getgEIprofiles(center,cmpgEtoP);
    Ispread=getgEIprofiles(center,cmpgItoP);
   for xx=1:ncells
       for yy=1:ncells
           synContE(xx,yy) = max(Espread(xx,yy,:));
           synContI(xx,yy) = max(Ispread(xx,yy,:));
       end
   end
    synCont= synContE*(-45)+synContI*(-45+80);
    synCont(find(synCont<cThold))=1;
   
    function twoDMat=getgEIprofiles(center,cmpEI)
        x_length=round(ncells/2);ms=0.1*(1:R_length);
        fProfE=zeros(x_length);fProfI=fProfE;
        rast=zeros(1,R_length);rast=rast*NaN;
        PSTH=zeros(1,R_length);tPSTH=zeros(1,R_length);
        twoDMat=oneD2D(cmpEI);

        function tMat=oneD2D(cmpEI)
            tMat=zeros(ncells,ncells,R_length);
            for i=1:nxn
            cellPosX = mod(fix((i-1)),ncells)+1;  %finds x and y values of activated inhcells  ***COLUMN
            cellPosY = floor(fix((i-1)/ncells));
            if ((cellPosX>0)&&(cellPosY>0))
               tMat(cellPosX,cellPosY,:)=cmpEI(i,:);
            end
            end
        end
    end
end
end

end
