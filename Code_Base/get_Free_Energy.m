function fes = get_Free_Energy(T, masked_residues)
%% FesCalc_Block_Function - Computes 1D free energy profile with charge masking
% INPUTS:
%   T               - Temperature in Kelvin (scalar)
%   masked_residues - List of residue indices whose charges should be set to zero (vector, can be empty [])
% OUTPUTS:
%   fes             - [N x 2] array: column 1 = # structured blocks, column 2 = Free Energy (kJ/mol)

%% Input Parameters
pdb = char('1UBQ'); % Input PDB file name, Make sure to run cmapCalcElec_Block.m and generate necessary input files
stridefile = char('struct.txt'); % File containing raw output from the STRIDE server
ene = -0.091; %vdW internation energy (J/mol) per native contact
DS = -14/1000; % Entropic cost (J/mol.K) per residue
DCp = -0.36/1000; % Heat capacity change (J/mol.K) per native contact
IS = 0.1; % Ionic strength in Molar units

%% Contact Maps obtained through PDB

aa=pdb(1,:);
eval(['load contactmapmatElecB',aa,'.dat;']);
eval(['cmapmask=contactmapmatElecB',aa,';']);
eval(['load contdistElecB',aa,'.dat;']);
eval(['contdist=contdistElecB',aa,';']);
eval(['load BlockSize',aa,'.dat;']);
eval(['bs=BlockSize',aa,';']);
eval(['load BlockDet',aa,'.dat;']);
eval(['BlockDet=BlockDet',aa,';']);
nres=length(cmapmask);
nres1 = length(unique(BlockDet(:,1)));
Mw=nres1*110;
wi=10;             % for number of calculating windows in case of 'NaN's in the free energy projections

fes_array = cell(length(ene), 1);
Fpath_array = cell(length(ene), 1);
fes2D_array = cell(length(ene), 1);
ResProb_array = cell(length(ene), 1);

disr = []; % non-helical/non-strand/non-310 helix residues (see lines 46-62)
ppos=[]; % proline residues (see lines 46-62)
C=0;

aa2={'GLY','ALA','VAL','LEU','ILE','MET','PHE','TYR','TRP','SER','ASP','ASN','THR','GLU','GLN','HIS','LYS','ARG','PRO','CYS'};
aacode={'G','A','V','L','I','M','F','Y','W','S','D','N','T','E','Q','H','K','R','P','C'};

hh9=fopen(stridefile(1,:),'rt');
lkk=fgetl(hh9);
lengthprot = 1;
while lkk > 0
    if length(lkk)> 25 && strcmp(lkk(1:3),'ASG')
        if strcmp(lkk(6:8),'PRO')
            ppos(end+1)=str2double(lkk(17:20));
        elseif strcmp(lkk(6:8),'GLY') || (strcmp(lkk(25),'H')==0 && strcmp(lkk(25),'E')==0 && strcmp(lkk(25),'G')==0)
            disr(end+1)=str2double(lkk(17:20));
        end
        [~,x2]=find(strcmp(lkk(6:8),aa2));
        protseq(1,lengthprot)=aacode(x2);
        lengthprot = lengthprot+1;
    end
    lkk=fgetl(hh9);
end
fclose(hh9);

%% Apply charge masking for specified residues
for i = 1:length(masked_residues)
    res_idx = masked_residues(i);
    if res_idx >= 1 && res_idx <= size(contdist,1)
        % Zero out the charges for this residue in contdist
        mask = (contdist(:,1) == res_idx) | (contdist(:,2) == res_idx);
        contdist(mask, 5) = 0;
    end
end

%% Constants
R=0.008314;
Tref = 385;
zval=exp(DS./R);
zvalc=exp((DS-(6.0606/1000))./R); %% Excess entropic cost: dDS = -6.0606 J/mol.K (DOI: 10.1021/acs.jpcb.6b00658 )
zjj=zval.*ones(nres1,1);
zjj(disr)=zvalc;
zjj(ppos)=1;

zvec = ones(BlockDet(end,2),1);
for i=1:length(BlockDet)
    zvec(BlockDet(i,2)) = zvec(BlockDet(i,2)) * zjj(BlockDet(i,1));
end

%% loops for each ene value
for iei = 1:length(ene)
    %% calculate Zfin
    fesmat = zeros(nres,length(T));
    Zfin=zeros(length(T));
    for kk=1:length(T)
        for ll=1:length(C)
            %% Elec. energy involving IS contribution
            ISfac=5.66*sqrt(IS/T(kk))*sqrt(80/29);
            emapmask=zeros(nres);
            for i=1:nres
                for iin=i:nres
                    x1=find(contdist(:,1)==i & contdist(:,2)==iin);
                    emapmask(i,iin)=sum(contdist(x1,5).*exp(-ISfac.*contdist(x1,3)));
                end
            end
            pepval=zeros(nchoosek(nres+1,2)+2*nchoosek(nres+1,4),7);
            k=1;
            %% Generating Combinations for SSA
            for i=1:nres
                for iin=1:nres-i+1
                    sw=0;
                    stabEtemp=0; eneEtemp=0; nconttemp=0;
                    wii=floor(iin/wi);
                    swgen=zeros(wi,1);
                    wistart=i;
                    for wiii=1:wi-1
                        nconttemp=sum(sum(cmapmask(wistart:i+wii*wiii,i:i+iin-1)));
                        stabEtemp=(nconttemp*(ene+DCp*(T(kk)-Tref)-T(kk)*DCp*log(T(kk)/Tref)));
                        eneEtemp=sum(sum(emapmask(wistart:i+wii*wiii,i:i+iin-1)));
                        swgen(wiii,1)=exp(-(stabEtemp+eneEtemp)/(R*T(kk)))*prod(zvec(wistart:i+wii*wiii));
                        wistart=i+wii*wiii+1;
                    end
                    nconttemp=sum(sum(cmapmask(wistart:i+iin-1,i:i+iin-1)));
                    stabEtemp=(nconttemp*(ene+DCp*(T(kk)-Tref)-T(kk)*DCp*log(T(kk)/Tref)));
                    eneEtemp=sum(sum(emapmask(wistart:i+iin-1,i:i+iin-1)));
                    swgen(wi,1)=exp(-(stabEtemp+eneEtemp)/(R*T(kk)))*prod(zvec(wistart:i+iin-1));
                    sw=prod(swgen);
                    
                    pepval(k,:)=[(iin) sw i iin 0 0 1];
                    k=k+1;
                end
                
            end
            
            %% Generating Combinations for DSA
            for i=1:nres
                for iin=1:nres-i+1
                    for j=i+iin+1:nres
                        for jin=1:nres-j+1
                            
                            %for island 1
                            sw1=0;
                            stabEtemp=0; eneEtemp=0; nconttemp=0;
                            wii=floor(iin/wi);
                            swgen=zeros(wi,1);
                            wistart=i;
                            for wiii=1:wi-1
                                nconttemp=sum(sum(cmapmask(wistart:i+wii*wiii,i:i+iin-1)));
                                stabEtemp=(nconttemp*(ene+DCp*(T(kk)-Tref)-T(kk)*DCp*log(T(kk)/Tref)));
                                eneEtemp=sum(sum(emapmask(wistart:i+wii*wiii,i:i+iin-1)));
                                swgen(wiii,1)=exp(-(stabEtemp+eneEtemp)/(R*T(kk)))*prod(zvec(wistart:i+wii*wiii));
                                wistart=i+wii*wiii+1;
                                stabEtemp=0; eneEtemp=0; nconttemp=0;
                            end
                            nconttemp=sum(sum(cmapmask(wistart:i+iin-1,i:i+iin-1)));
                            stabEtemp=(nconttemp*(ene+DCp*(T(kk)-Tref)-T(kk)*DCp*log(T(kk)/Tref)));
                            eneEtemp=sum(sum(emapmask(wistart:i+iin-1,i:i+iin-1)));
                            swgen(wi,1)=exp(-(stabEtemp+eneEtemp)/(R*T(kk)))*prod(zvec(wistart:i+iin-1));
                            sw1=prod(swgen);
                            
                            %for island 2
                            sw2=0;
                            stabEtemp=0; eneEtemp=0; nconttemp=0;
                            wii=floor(jin/wi);
                            swgen=zeros(wi,1);
                            wistart=j;
                            for wiii=1:wi-1
                                nconttemp=sum(sum(cmapmask(wistart:j+wii*wiii,j:j+jin-1)));
                                stabEtemp=(nconttemp*(ene+DCp*(T(kk)-Tref)-T(kk)*DCp*log(T(kk)/Tref)));
                                eneEtemp=sum(sum(emapmask(wistart:j+wii*wiii,j:j+jin-1)));
                                swgen(wiii,1)=exp(-(stabEtemp+eneEtemp)/(R*T(kk)))*prod(zvec(wistart:j+wii*wiii));
                                wistart=j+wii*wiii+1;
                            end
                            nconttemp=sum(sum(cmapmask(wistart:j+jin-1,j:j+jin-1)));
                            stabEtemp=(nconttemp*(ene+DCp*(T(kk)-Tref)-T(kk)*DCp*log(T(kk)/Tref)));
                            eneEtemp=sum(sum(emapmask(wistart:j+jin-1,j:j+jin-1)));
                            swgen(wi,1)=exp(-(stabEtemp+eneEtemp)/(R*T(kk)))*prod(zvec(wistart:j+jin-1));
                            sw2=prod(swgen);
                            
                            sw=sw1.*sw2;
                            
                            pepval(k,:)=[(iin+jin) sw i iin j jin 2];
                            k=k+1;
                        end
                    end
                end
                
            end
            
            %% Generating Combinations for DSAw/L
            for i=1:nres
                for iin=1:nres-i+1
                    for j=i+iin+1:nres
                        for jin=1:nres-j+1
                            stabE=0; eneE=0; sw=0;
                            vv = [(i:i+iin-1) (j:j+jin-1)];
                            ncont=sum(sum(cmapmask(vv,vv)));
                            stabE=(ncont*(ene+DCp*(T(kk)-Tref)-T(kk)*DCp*log(T(kk)/Tref)));
                            eneE=sum(sum(emapmask(vv,vv)));
                            if sum(sum(cmapmask(i:i+iin-1,j:j+jin-1)))~=0 || sum(sum(emapmask(i:i+iin-1,j:j+jin-1)))~=0
                                %from island 1
                                sw1=0;
                                stabEtemp=0; eneEtemp=0; nconttemp=0;
                                wii=floor(iin/wi);
                                swgen=zeros(wi,1);
                                wistart=i;
                                for wiii=1:wi-1
                                    nconttemp=sum(sum(cmapmask(wistart:i+wii*wiii,vv)));
                                    stabEtemp=(nconttemp*(ene+DCp*(T(kk)-Tref)-T(kk)*DCp*log(T(kk)/Tref)));
                                    eneEtemp=sum(sum(emapmask(wistart:i+wii*wiii,vv)));
                                    swgen(wiii,1)=exp(-(stabEtemp+eneEtemp)/(R*T(kk)))*prod(zvec(wistart:i+wii*wiii));
                                    wistart=i+wii*wiii+1;
                                end
                                nconttemp=sum(sum(cmapmask(wistart:i+iin-1,vv)));
                                stabEtemp=(nconttemp*(ene+DCp*(T(kk)-Tref)-T(kk)*DCp*log(T(kk)/Tref)));
                                eneEtemp=sum(sum(emapmask(wistart:i+iin-1,vv)));
                                swgen(wi,1)=exp(-(stabEtemp+eneEtemp)/(R*T(kk)))*prod(zvec(wistart:i+iin-1));
                                sw1=prod(swgen);
                                
                                %from island 2
                                sw2=0;
                                stabEtemp=0; eneEtemp=0; nconttemp=0;
                                wii=floor(jin/wi);
                                swgen=zeros(wi,1);
                                wistart=j;
                                for wiii=1:wi-1
                                    nconttemp=sum(sum(cmapmask(wistart:j+wii*wiii,vv)));
                                    stabEtemp=(nconttemp*(ene+DCp*(T(kk)-Tref)-T(kk)*DCp*log(T(kk)/Tref)));
                                    eneEtemp=sum(sum(emapmask(wistart:j+wii*wiii,vv)));
                                    swgen(wiii,1)=exp(-(stabEtemp+eneEtemp)/(R*T(kk)))*prod(zvec(wistart:j+wii*wiii));
                                    wistart=j+wii*wiii+1;
                                end
                                nconttemp=sum(sum(cmapmask(wistart:j+jin-1,vv)));
                                stabEtemp=(nconttemp*(ene+DCp*(T(kk)-Tref)-T(kk)*DCp*log(T(kk)/Tref)));
                                eneEtemp=sum(sum(emapmask(wistart:j+jin-1,vv)));
                                swgen(wi,1)=exp(-(stabEtemp+eneEtemp)/(R*T(kk)))*prod(zvec(wistart:j+jin-1));
                                sw2=prod(swgen);
                                
                                sw=sw1*sw2*zvalc^(j-(i+iin));
                            else
                                sw = 0;
                            end
                            pepval(k,:)=[(iin+jin) sw i iin j jin 3];
                            k=k+1;
                        end
                    end
                end
              
            end
            
            Zfin(kk) = sum(pepval(:,2));
            
            pepval = pepval(pepval(:,2)~=0,:);
            
            %% 1D free energy profile
            fes = zeros(nres,1);
            for i=1:nres
                fes(i) = sum(pepval(pepval(:,1)==i,2));
            end
            fes = fes./sum(fes);
            fes = -R*T(kk)*log(fes);
            
            fesmat(:,kk)=fes;
        end
    end
    fes_array{iei} = fes;
end
fes = [(1:nres)' fes];

end