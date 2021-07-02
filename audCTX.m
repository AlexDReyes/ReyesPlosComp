
function audCTX() 

    %-----simulation parameters---------
    ncells =100; % number of neurons in network= ncells x ncells
    repetition = 1; %number of sweeps
    R_length=500; delT=0.1; %record length and time increment; total duration = delT*R_length
    
    %-----external excitatory (thalamic) input to E and I cells-------
    inSigmaA=10; %sigma of external input to cluster of neurons on left; in units of cells (multiply by 7.4 for equivalent in microns)
    inSigmaB=0;  %sigma of external input to cluster of neurons on right;   Set to 0 to toggle off, 
    sepX=0;    %separation between 2 Gaussians
    
    thalNum=100;    %maximum number of thalamic inputs to each neuron

    %-----local inhibitory (fast spiking cell) input to E -------
    inSigmaI = 10; %sigma of external input to cluster of INHIBITORY neurons; in units of cells (multiply by 7.4 for equivalent in microns)
    inhMode=0;  %inhMode=0, normal; inhMode=1, independent; see Readme for details
    
    %-------Network  architecture--------------
     RW=0; %if RW=0,reads already created network file; RW=1 creates new ncells x ncells network; see Readme for details    
    
     %------Start simulations-----------
     for rep=1:repetition
       fname=['myFile' num2str(rep)]; %file name for output
       hname ='networkFileName';% file name of network; contains info about network connections
      mult(thalNum,inSigmaA,inSigmaB,inSigmaI,inhMode,fname,RW,sepX,rep); %runs the program
     end
     
    AudCTX_analysis("myFile"); %runs analysis program; comment to toggle off


%-------------gateway to main function----------------
    function mult(nIn,inSigA,inSigB,inSigI,inhMode,fname,RW,sepX,rep)
        if ((RW==1)&&(rep>1))
            RW=0;   %only makes network once...takes a long time
        end
        fname
        GoAud1_0(nIn,ncells,inSigA,inSigB,inSigI,hname,fname,inhMode,RW,R_length,delT,repetition,sepX);
    end

end

    
%-------------main function--------------
function GoAud1_0(nInput,ncells,inSigA,inSigB,inSigI,iName,oName,inhMode,RW,R_length,delT,repetition,sepX)
STPon=1;STPthal=1; %1= synaptic depression/facilitation is on; 0=off

%------connection architecture from Levy and Reyes, 2012---------
eSigma=14.3;eiSigma=10.125;ieSigma=9.855;%ieSigmaB=4.455;eiSigmaB=14.715; %These are in "cell units" (multiply by 7.4 for equivalent microns)

nIn=nInput;RorW=RW;fname=oName;
inhcells = round(sqrt(0.1*ncells^2));
dscalar=1;    
inDataName = iName;
totalE = ncells * ncells; totalI = inhcells*inhcells;
nInb=5*nIn; inMatrix = zeros(nInb,R_length);inMatrixI=inMatrix;inMatrixIB=inMatrix;
vTrace=zeros(1,R_length);gTrace=vTrace;vETrace=zeros(1,R_length);gETrace=vTrace;gITrace=vTrace;sEtrace=vTrace;

%-------Synaptic Parameters from Levy and Reyes, 2011 and 2012---------------

    EPSGtoPamp = 0.000145;EPSGtoIamp =0.00047;thalToPamp=0.0007;thalToIamp=0.0027;%amplitudes       
    IPSGtoPamp=0.0005; 
    
    EPSGtoPalpha=0.375;EPSGtoIalpha=0.55; IPSGtoPalpha=0.25;thalToPalpha=0.9;thalToIalpha= 3;   %alpha functions
    
    %Depression/facilitation paramters
    UPtoP=0.55;tauFacPtoP=0;tauRecPtoP=0.45;DfPtoP=0.9;RtPtoP=zeros(totalE,2);RtPtoP(:,1)=1-UPtoP;
    UPtoI=0.55;tauFacPtoI=0;tauRecPtoI=0.6;DfPtoI=0.7;  RtPtoI=zeros(totalE,2);RtPtoI(:,1)=1-UPtoI;
    UItoP=0.555;tauFacItoP=0;tauRecItoP=0.375;DfItoP=0.5;  RtItoP=zeros(totalE,2);RtItoP(:,1)=1-UItoP;
    UTtoP=0.25;tauFacTtoP=0;tauRecTtoP=0.5;DfTtoP=0.755;
    UTtoI=0.55;tauFacTtoI=0;tauRecTtoI=0.6;DfTtoI=0.55;

    %parameters for adaptive exponential LIFs
    CP=0.099;gLP=0.0055;ELP=-65;VTP=-45;slpTP=0.8; aP=-0.00095;tauWP=60;bP=0.25;VrP=-48;   
    CFS=0.067;gLFS=0.0084;ELFS=-67;VTFS=-47;slpTFS=3;aFS=0.00;tauWFS=16;bFS=0.08;VrFS=-50;
    
    if (inhMode==1)  %IPSPs to E cells are converted thalamic inputs, not from I cells and can be varied  "independently" from excitatory input
        IPSGtoPamp=thalToPamp;
        IPSGtoPalpha=thalToPalpha;
    end

%----------estimate number of connected cells for center cell
nh=round(ncells/2);
for x=1:ncells
    for y=1:ncells
        distance = (sqrt(((y-nh))^2 + ((x-nh))^2));
        gFunct1(x,y) = 0.2*exp(-((distance)^2)/(2*eSigma^2));
        gFunct2(x,y) = 0.5*exp(-((distance)^2)/(2*ieSigma^2));
    end
end
nEconn=round(sum(sum(gFunct1)));nIconn=round(sum(sum(gFunct2)));
nEIconn=round(max(nEconn,nIconn)+0.1*max(nEconn,nIconn));  %maximum + 10%

EPSGtoPMax=EPSGtoPalpha*exp(-1)*delT;EPSGtoIMax=EPSGtoIalpha*exp(-1)*delT;IPSGtoPMax=IPSGtoPalpha*exp(-1)*delT;
%EPSGtoIMaxB=EPSGtoIalphaB*exp(-1)*delT;IPSGtoPMaxB=IPSGtoPalphaB*exp(-1)*delT;


thalToPMax=thalToPalpha*exp(-1)*delT;
thalToIMax=thalToIalpha*exp(-1)*delT;
%thalToIMaxB=thalToIalphaB*exp(-1)*delT;
 if (RorW==0)       %read from files
   load(inDataName,'PtoPMatrix','PtoIMatrix','ItoPMatrix','inhList');
       PtoPMatrix = int32(zeros(ncells,ncells,nEconn))-1;
       PtoIMatrix=int32(zeros(ncells,ncells,nEIconn))-1;
       
        L=1;
        for i=1:ncells
            for j=1:ncells            %the first z value contains the identifier for the neuron; converted to X and y coord by 'mod' and 'fix' below
                PtoPMatrix(j,i,1) = L;      %for next step; j is index;  j is converted to X and y coord by 'mod' and 'fix' below;I is pre
                PtoIMatrix(j,i,1) = L;
                PtoIMatrixB(j,i,1) = L;
                L=L+1;
            end         %j loop
        end             %i loop
 end
if (RorW==1)
    PtoPMatrix = int32(zeros(ncells,ncells,nEconn))-1;     %identifies the postsynaptic targets of a pyr cell in the network
    PtoIMatrix=int32(zeros(ncells,ncells,nEIconn))-1;		%pyr to inhibi;  identifies the inhibitory target neurons of a pyramidal cell
    ItoPMatrix=int32(zeros(ncells,ncells,nEIconn))-1;		%inhibitory to pyr; identifies the excitatory targets of an interneuron
    inhList = int32(zeros(1,2*totalI))-1;
end
inputList = zeros(ncells,ncells,2); 
gEtoP = single(zeros(totalE,5));gEtoP(:,4)=EPSGtoPalpha;gEtoP(:,5)=EPSGtoPMax;    %conductances: 1st element=xOld;2nd = yOld;3rd = amp; 4th=alpha;5th=max val
gItoP = single(zeros(totalE,5));gItoP(:,4)=IPSGtoPalpha;gItoP(:,5)=IPSGtoPMax;

gEtoI = single(zeros(totalE,5));gEtoI(:,4)=EPSGtoIalpha;gEtoI(:,5)=EPSGtoIMax;
gthalToP = single(zeros(totalE,5));gthalToP(:,4)=thalToPalpha;gthalToP(:,5)=thalToPMax;
gthalToI = single(zeros(totalE,5));gthalToI(:,4)=thalToIalpha;gthalToI(:,5)=thalToIMax;
%--------randomize amplitudes------------
for i= 1:totalE
    gEtoP(i,3)=PSCran(EPSGtoPamp);
    gItoP(i,3)=PSCran(IPSGtoPamp);
    gEtoI(i,3) = PSCran(EPSGtoIamp);
    gthalToP(i,3)=PSCran(thalToPamp);
    gthalToI(i,3)=PSCran(thalToIamp);
end

spkEtoP = single(zeros(totalE,1));  %presynaptic spikes
spkItoP = single(zeros(totalE,1));
spkEtoI = single(zeros(totalE,1));
spkThalToP=zeros(totalE,R_length);
spkThalToI=zeros(totalE,R_length);
EspkTimes=zeros(totalE,R_length);         %will store spike times; should grow automatically if needed.
IspkTimes=zeros(totalE,R_length);
Eidx=ones(totalE,1);                %stores indices
Iidx=ones(totalE,1);
gEwave=zeros(totalE,R_length);gEwaveR=zeros(totalE,R_length);
gIwave=zeros(totalE,R_length);


vE = single(zeros(totalE,5));
vE(:,3)=ELP;         %resting potential
vI = single(zeros(totalE,5));
vI(:,3)=ELFS;

tic
     prime_it(ncells,inhcells,R_length);%'primed'
     getNetwork(RorW);%'gotNetwork'
     makeInputs(nIn,inSigA,inSigB,inSigI,sepX);%'made inputs'
     run();
toc
cmpE=EspkTimes;cmpI = IspkTimes;cmpgEtoP = gEwave;cmpgItoP = gIwave;cmpIB=0*cmpI;cmpgItoPB = 0*cmpgItoP;cmpgEtoPr=0*cmpgEtoP; %necessary to rename files for analysis program
save(fname,'cmpE','cmpI','cmpgEtoP','cmpgItoP','ncells','R_length','delT','repetition');


 

 %------------------prepare simulation------------------//f
 function prime_it(ncells,inhcells,R_length)
    totalE = ncells * ncells; totalI = inhcells*ncells;
    if (RorW==1)
        PtoPMatrix = 0*PtoPMatrix-1;     %identifies the postsynaptic targets of a pyr cell in the network
        PtoIMatrix=0*PtoIMatrix-1;		%pyr to inhibi;  identifies the inhibitory target neurons of a pyramidal cell
        ItoPMatrix=0*ItoPMatrix-1;		%inhibitory to pyr; identifies the excitatory targets of an interneuron
        inhList = 0*inhList-1;
    end
    inputList = 0*inputList;
  end
 %------------------------------------------------------//
 
 %---------input or make new network-------------------//
function getNetwork(RorW)       %***NOTE: 1ST INDEX IS ROW, 2ND IS COLUMN e.g. A(2,3) points to 2nd row, 3rd column
 if (RorW==1)       %make new network
     getNet();
 end
 %---------start getNet--------
     function getNet()
     labelMatrices();
     ItoPConnections(ncells);
save(inDataName,'PtoPMatrix','PtoIMatrix','ItoPMatrix','inhList','-v7.3');


     %----------Label matrices---------
     function labelMatrices()				%assigns wave number to all excitatory cells in the network e.g. PtoPmatrix[24][2] = cell76
        L=1;
        incr=round(ncells/inhcells);
        for i=1:ncells
            for j=1:ncells            %the first z value contains the identifier for the neuron; converted to X and y coord by 'mod' and 'fix' below
                PtoPMatrix(j,i,1) = L;      %for next step; j is index;  j is converted to X and y coord by 'mod' and 'fix' below;I is pre
                PtoIMatrix(j,i,1) = L;
                ItoPMatrix(j,i,1) = L;               
                L=L+1;
            end         %j loop
        end             %i loop        
        z=PtoPMatrix(:,:,1);
        z=z(incr:incr:ncells-incr,incr:incr:ncells-incr);
        inhList=z(find(z~=0));
     end
    %-----------------------------------
     function cProb = calcProb(distance,IorP)
        connected=0;
        if (IorP==0)			%pyramidal
            sqrxSigma = 2*eSigma^2;
            prob = (0.077)*exp(-(distance*distance)/(sqrxSigma));
        end
        if (IorP==1)			%inhibitory to PYR
            sqrxSigma = 2*ieSigma^2;
            
            prob = (0.39)*exp(-(distance*distance)/(sqrxSigma));
            
        end
        if (IorP==2)			%PYR to inhibitory
            sqrxSigma = 2*eiSigma^2;
            prob = (0.27)*exp(-(distance*distance)/(sqrxSigma));
        end
        On = abs(unifrnd(0,1));
        if (On<=prob)
            connected=1;
        end
        cProb=connected;
     end

    %------------start ItoPConnections-------------
    function ItoPConnections(ncells)			
    dscalar=1;
    idx=1;
    stop=numel(inhList);
    while ((idx<=stop)&&(inhList(idx)>-1))    %look at all active inhibitory cells
      cellPosY=mod(inhList(idx),ncells);
      if (cellPosY==0)
          cellPosY=ncells;
      end
      cellPosX=find(PtoPMatrix(cellPosY,:,1)==inhList(idx));

        idx=idx+1;
        Z=2;Zb=2;
        counts=0;
        for X=1:ncells      %go through each potential target pyramidal cell and find which cells the inhib (inhList) projects to
            for Y=1:ncells
              distance = (sqrt(double(dscalar*(Y-cellPosY))^2 + (dscalar*(X-cellPosX))^2));		
                connected = calcProb(distance,1);
                if ((connected==1)&& (cellPosY>0)&& (cellPosX>0))					%if the output of function is 1,i.e. connected=1
                    ItoPMatrix(cellPosY,cellPosX,Z) = PtoPMatrix(Y,X,1);
                    counts=counts+1;
                    Z=Z+1;
                end
                ItoPMatrix(cellPosY,cellPosX,numel(ItoPMatrix(cellPosY,cellPosX,:)))=-1;
            end
        end 		%row of ncells
    end
 

    end
    %------------End ItoPConnections----------------

    end
 %--------end getNet-------
 
end
%----------end GetNetwork---------------------//

%----------start Make_inputs-----------
function makeInputs(nIn,inSigA,inSigB,inSigI,sepX)       %puts all functions necessary for inputs into one function

nInb = 5*nIn;
discX=round(ncells/2);discY = round(ncells/2);  %center of disc;round((1/4)*ncells)
    DiscInput(inSigA,discX-sepX,discY-sepX,1);
    if (inSigB>0)
     DiscInput(inSigB,discX+sepX,discY+sepX,1);
    end
    if (inSigI>0)
        if (inhMode==0)
            if (sepX==0)
                DiscInput(inSigI,discX,discY,2);
            end
            if ((sepX>0)&&(inSigB>0)&&(inSigA>0))
              DiscInput(inSigI,discX-sepX,discY-sepX,2);
              DiscInput(inSigI,discX+sepX,discY+sepX,2);
            end
        end
        if (inhMode==1)
          DiscInput(inSigI,discX+sepX,discY+sepX,2);
        end
    end



getThalInB(nIn);
%----Start DisscInput------
function DiscInput(inSigma,discX,discY,EorI)       %calculates disc input to network sheet
    for i=1:ncells
        for j=1:ncells
            distance = sqrt((i-discX)^2 + (j-discY)^2);
             inputList(i,j,EorI)=CalcProbThal(distance,inSigma)+inputList(i,j,EorI);    %1 is for excitatory
        end
    end
    %----------CalcProbThal------
    function prob=CalcProbThal(distance,inSigma)				
        sqrxSigma = 2*inSigma^2;
        prob = exp(-(distance*distance)/(sqrxSigma));
    end
    %----------end CalcProbThal-------
    
end
%------end DiscInput------

%-------start getthalInB-------
function getThalInB(nIn)            %gets inputs to cells in the first layer
freq=100; %average firing rate of thalamic cell
spikes=zeros(1,R_length);
inMatrix=inMatrix*0;inMatrixI=inMatrixI*0;inMatrixIB=inMatrixIB*0;tSpikes=zeros(1,R_length);
%TC inputs to excitatory cells------------//
nInb = nIn*5;  nxn=ncells^2;   dIndex=1;
for i=1:ncells
    for j=1:ncells
        spikes=spikes*0;
        nPreIn = (nIn*inputList(i,j,1));
        nPreInB=nPreIn;
        while (nPreInB<0)       %this is not working!! should be set to 0 so it enters loop
          nPreInB=normrnd(nPreIn,0.25*nPreIn);
        end
        fPreIn=(freq*inputList(i,j,1));
        fPreInB=fPreIn;
         while (fPreInB<0)
          fPreInB=normrnd(fPreIn,0.25*fPreIn);
        end

        idx=1;n=1;
        while (n<nPreInB)
            tSpikes=inputTrainsB(fPreInB,0);
            spikes=spikes+tSpikes(1:R_length);
            n=n+1;
        end
        spkThalToP(dIndex,:)=spkThalToP(dIndex,:)+spikes;
        
        
           %TC inputs to inhibitory cells---------------//
        spikes=spikes*0;
        tSpikes=tSpikes*0;
        nPreIn = round((nIn*inputList(i,j,2)));
        nPreInB=nPreIn;
         while (nPreInB<0)
          nPreInB=normrnd(nPreIn,0.25*nPreIn);
         end
         fPreIn=(freq*inputList(i,j,2));
        fPreInB=fPreIn;
         while (fPreInB<0)
          fPreInB=normrnd(fPreIn,0.25*fPreIn);
        end

        idx=1;n=1;
        while (n<nPreInB)
            tSpikes=inputTrainsB(fPreInB,1);
            spkThalToI(dIndex,:)=spkThalToI(dIndex,:)+tSpikes(1:R_length);
            n=n+1;
        end    
        dIndex=dIndex+1;
    end    
    
end
    function amp = PSCran(inAmp)

    amp=-1;
        while (amp<=0)
           amp = normrnd(inAmp,0.25*inAmp);       
        end
    end
%------start inputTrains-------
function Rtrain=inputTrainsB(xFreq,cType)
    Rtrain=zeros(1,R_length);
    rendd = R_length-100;%R_length-500;				%ends at 100 before R_length
    ISI = 1000/(delT*xFreq);			%ISI in counts
    refract=5/delT; %refractory period
    gnoiseSTD = 0.5*ISI;			%adds jitter; 10% of ISI
    AnoiseSTD=0.25;     %adds jitter to amplitude
    incrOld=-1;
    while (fix(incrOld)<=0)
        incrOld=normrnd(400,100);
    end
    ArandScale=0;
    while (ArandScale<=0);
        ArandScale=normrnd(1,AnoiseSTD); 
    end
    if (cType==0)   %pyr
     tauRecT=tauRecTtoP;   
     tauFacT=tauFacTtoP;
     UT=UTtoP;
     Tscale = 1/(UT*(1-UT));
     RtT=1-UT;
    end
    if (cType==1)   %FS
     tauRecT=tauRecTtoI;   
     tauFacT=tauFacTtoI;
     UT=UTtoI;
     Tscale = 1/(UT*(1-UT));
     RtT=1-UT;
    end
    if (cType==2)   %nonFS
        tauRecT=tauRecTtoIB;
        tauFacT=tauFacTtoIB;
        UT=UTtoIB;
        Tscale = 1/(UT*(1-UT));
        RtT=1-UT;
    end
    UTn=0;
    RTUT=ArandScale;
    Rtrain(fix(incrOld))=RTUT;%baseRTUoldP;
    latency = 0;
    incr=0;
    while (incr<=rendd)
        tmp=0;tstop=0;
        while (tmp<=tstop)
            nse=0;
            while(nse<=refract+incrOld)
                nse=incrOld+ISI+normrnd(0,gnoiseSTD);
            end

            tmp= nse;
        end
        incr = fix(tmp+latency);
        if (STPthal==1)
            tDiff=delT*(incr-incrOld);           
              [RTUT,RtT,UTn]=calc_height(tDiff,RtT,UTn,tauRecT,tauFacT,UT);RTUT=RTUT*Tscale*ArandScale;
            
            if (numel(find((RTUT>0) | (RTUT<=0)))==0)
                RTUT=0;
            end

        end
        incrOld=incr;	
        if (incr<rendd)
            Rtrain(incr)=RTUT;    %each row represents an input train
        end
    end
   end

end
%-------end getthalInB----------
end

%-----------end Make_inputs---------


%---------update the conductances---------
    function gIn=updateG(gIn,spkIn,idxx)
        totalcells = ncells^2;
        for k=1:totalE
            xOld=gIn(k,1);yOld = gIn(k,2);alphaIn=gIn(k,4);
            delX=yOld*delT;
            delY=(-alphaIn^2*xOld-2*alphaIn*yOld + alphaIn^2*(spkIn(k,idxx)))*delT;
            gIn(k,1) = xOld+delX;
            gIn(k,2)=yOld+delY;          
        end
    end
%-----------end update-----------------------

%-----------------Start calc_ISI--------------
function [RTU,Rn,Un]=calc_height(ISI,RnOld,UnOld,tauRec,tauFac,U)
ISI=ISI/1000;
    Un=UnOld*exp(-ISI/tauFac)+U*(1-UnOld*exp(-ISI/tauFac));
    Rn=RnOld*(1-Un)*exp(-ISI/tauRec)+1-exp(-ISI/tauRec);
    RTU=Rn*Un;
end


%------------start run-------------//
function run()
 totalcells = ncells*ncells;    %total number of cells in the network
 gthalToP(:,1)=0;gthalToI(:,1)=0;gthalToIB(:,1)=0;
 gthalToP(:,2)=0;gthalToI(:,2)=0;gthalToIB(:,2)=0;
 
for i=3:R_length        %1st column is used for cell identification
    gEtoP=updateG(gEtoP,spkEtoP,1);
    gEtoI=updateG(gEtoI,spkEtoI,1);
    
    if (inhMode==0)
        gItoP=updateG(gItoP,spkItoP,1);
    end
    if (inhMode==1)
        gItoP=updateG(gItoP,spkThalToI,1);
    end

    gthalToP=updateG(gthalToP,spkThalToP,i);
    gthalToI=updateG(gthalToI,spkThalToI,i);
    
    if (inhMode==1)
        gItoP=gthalToI;
    end    
    idx=1;
    i;
    ctr=1;
    spkEtoP(:,1)=0;spkItoP(:,1)=0;spkEtoI(:,1)=0;spkItoPB(:,1)=0;spkEtoIB(:,1)=0;
    for j=1:totalcells
      cellPosY=mod(j,ncells);
      if (cellPosY==0)
          cellPosY=ncells;
      end
      cellPosX=find(PtoPMatrix(cellPosY,:,1)==j);
      gEtoPtotal=gEtoP(j,1)*gEtoP(j,3)/gEtoP(j,5)+gthalToP(j,1)*gthalToP(j,3)/gthalToP(j,5);
      gEtoItotal=gEtoI(j,1)*gEtoI(j,3)/gEtoI(j,5)+gthalToI(j,1)*gthalToI(j,3)/gthalToI(j,5);      
      gItoPtotal=(gItoP(j,3)/gItoP(j,5))*gItoP(j,1);
      
     if ((cellPosX>0)&&(cellPosY>0))
            gEwave(j,i)=gEtoPtotal; %stores conductance waves
            gEwaveR(j,i)=gEtoP(j,1)*gEtoP(j,3)/gEtoP(j,5);
            gIwave(j,i)=gItoPtotal;            
            gItoPtotal=gIwave(j,i);            
            [vE(j,3),vE(j,2)] = adEXP(vE(j,3),vE(j,2),(gEtoPtotal),gItoPtotal,CP,gLP,ELP,VTP,slpTP,aP,tauWP,bP,VrP);  
            if (vE(j,3)>=0)      
                EspkTimes(j,Eidx(j,1))=i;
                RTUPtoP=1;RTUPtoI=1;%RTUPtoIB=1;
                if ((STPon==1)&&(Eidx(j,1)>1))
                    PtoPscale=1/(UPtoP*(1-UPtoP));PtoIscale=1/(UPtoI*(1-UPtoI));%PtoIBscale=1/(UPtoIB*(1-UPtoIB));
                     tDiff=delT*(i-EspkTimes(j,Eidx(j,1)-1));
                     [RTUPtoP,RtPtoP(j,1),RtPtoP(j,2)]=calc_height(tDiff,RtPtoP(j,1),RtPtoP(j,2),tauRecPtoP,tauFacPtoP,UPtoP);
                     RTUPtoP=RTUPtoP*PtoPscale;
                     [RTUPtoI,RtPtoI(j,1),RtPtoI(j,2)]=calc_height(tDiff,RtPtoI(j,1),RtPtoI(j,2),tauRecPtoI,tauFacPtoI,UPtoI);
                     RTUPtoI=RTUPtoI*PtoIscale;
                end
                k=2;
                while (PtoPMatrix(cellPosY,cellPosX,k)~=-1)      %change conductance in postsynaptic pyrs
                    target = PtoPMatrix(cellPosY,cellPosX,k);%+1;
                    k=k+1;
                      spkEtoP(target,1) = spkEtoP(target,1) + RTUPtoP;
                end
                k=2;
                while (PtoIMatrix(cellPosY,cellPosX,k)~=-1)      %change conductance in postsynaptic inhibs
                    target = PtoIMatrix(cellPosY,cellPosX,k);%+1;
                        spkEtoI(target,1) = spkEtoI(target,1) + RTUPtoI;
                    k=k+1;
                end
                
            Eidx(j,1)=Eidx(j,1)+1;
          end
            if find(inhList==j)
                [vI(j,3),vI(j,2)] = adEXP(vI(j,3),vI(j,2),(gEtoItotal),0,CFS,gLFS,ELFS,VTFS,slpTFS,aFS,tauWFS,bFS,VrFS);
                if ((cellPosY==21)&&(cellPosX==18))
                end

                if (vI(j,3)>=0)       %threshold       use 25 for Vinh; -45 for LIF
                    IspkTimes(j,Iidx(j,1))=i;                    
                    RTUItoP=1;
                    if ((STPon==1)&&(Iidx(j,1)>1))
                        ItoPscale=1/(UItoP*(1-UItoP)); 
                        tDiff=delT*(i-IspkTimes(j,Iidx(j,1)-1));                           
                        [RTUItoP,RtItoP(j,1),RtItoP(j,2)]=calc_height(tDiff,RtItoP(j,1),RtItoP(j,2),tauRecItoP,tauFacItoP,UItoP);
                        RTUItoP=RTUItoP*ItoPscale;
                    end
                    k=2;
                     while (ItoPMatrix(cellPosY,cellPosX,k)~=-1)      %change inhib conductance in postsynaptic pyrs
                        target = ItoPMatrix(cellPosY,cellPosX,k);%+1;
                          spkItoP(target,1) = spkItoP(target,1) +RTUItoP;
                  k=k+1;
                     end
                    Iidx(j,1)=Iidx(j,1)+1;
                end
                
            end  %  end if
            
        end

    end     %end j loop
    
end         %end i loop


function [lV,lw]=adEXP(V,w,gE,gI,C,gL,EL,VT,slpT,a,tauW,b,Vr)
 
  delV=0;delW=0;
  if ((V<0)||(V>=0))
      if (V<0)
        delV= (delT/C)*( -gL*(V-EL)+gL*slpT*exp(((V-VT)/slpT)) + (gI*(-80-V))+(gE*(0-V))-w);  
        lV=V+delV;
      delW=(delT/tauW)*(a*(V-EL)-w);
      lw=w+delW;

        if (lV>=0)
          lV=0;
          lw=w+b;
        end

      end
     if (V>=0)
         lV=Vr;
         lw=w;
     end
  else
      lw=0;
      lV=EL;
  end
end    
end
%----------------end run-------------------//

end
    
    function amp = PSCran(inAmp)
        amp = abs(normrnd(inAmp,0.25*inAmp));       
    end





