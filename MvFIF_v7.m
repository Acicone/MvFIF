function [IMF,stats] = MvFIF_v7(f,options,M)

%
%  function IMF = MvFIF_v7(f,options,M)
%
% It generates the decomposition of the multivariate signal f :
%
%  f = IMF(1,:) + IMF(2,:) + ... + IMF(K, :)
%
% where the last row in the matrix IMF is the trend and the other rows
% are actual IMFs
%
%                                Inputs
%
%   f         Signal to be decomposed. The algorithm assumes that the data
%             are arranged as a matrix whose columns correspond to different channels
%             (along the rows we have the evolution over time)
%
%   options    Structure, generated using function Settings_FIF_v3, containing
%              all the parameters needed in the various algorithms
%
%   M         Mask length values for each Inner Loop
%
%                               Output
%
%   IMF       Cell of D entris 8as many as the channels in f) each entry contains a matrix whose row i is the i-th IMF. The last row
%              contains the remainder-trend.
%
%   stats     Statistics regarding the IMFs
%               logM      Mask length values used for each IMF
%               posF      position of the first minimum in the filter DFT which is forced to become zero
%               valF      filter DFT first minimum value before the downward shift
%               inStepN   Number of inners steps performed by the method
%
%   See also SETTINGS_FIF_V3, GETMASK_V1_1, MAXMINS_v3_7.
%
%  Please cite: 
%
%  A. Cicone, H. Zhou. "Numerical Analysis for Iterative Filtering with 
%  New Efficient Implementations Based on FFT". Numerische Mathematik, 147 (1), pages 1-28, 2021.
%  doi: 10.1007/s00211-020-01165-5
%  ArXiv http://arxiv.org/abs/1802.01359
%
%  A. Cicone. "Iterative Filtering as a direct method for the decomposition 
%  of nonstationary signals". Numerical Algorithms, Volume 373, 2020,  112248. 
%  doi: 10.1007/s11075-019-00838-z
%  ArXiv http://arxiv.org/abs/1811.03536%
%
%  A. Cicone, "Multivariate fast iterative Filtering for the decomposition 
%  of nonstationary signals"
%  ArXiv https://arxiv.org/abs/1902.04860
%



%% deal with the input

if nargin < 1,  help MvFIF_v7; return; end
if nargin < 2, options = Settings_FIF_v3; end
if nargin < 3, M = []; end

FigCol = 'ckmygr'; % Plot Colors

[N,D] = size(f);
% if D>N
%     f = f.';
%     [N,D] = size(f);
% end

IMF =cell(1,D);
for i=1:D
    IMF{i} = zeros(N,options.NIMFs);
end

nameFile=sprintf('%1.0d',sum(round(clock*1000)));

load('prefixed_double_filter','MM');


%% Create a signal without zero regions and compute the number of extrema
g=f';
f_pp=real(acos(dot(normc(g(:,1:end-1)),normc(g(:,2:end)))));
f_pp(f_pp<=10^-18)=[];
if isempty(f_pp)
    disp('Signal too small')
    IMF=[];
    stats=[];
    return
end
maxmins_pp=Maxmins_v3_7(f_pp);
if isempty(maxmins_pp)
    IMF=f;
    stats=[];
    return
end
diffMaxmins_pp=diff(maxmins_pp);
N_pp=length(f_pp);
k_pp = length(maxmins_pp);

countIMFs=0;

while countIMFs < options.NIMFs && min(k_pp)>=options.ExtPoints
    countIMFs=countIMFs+1;
    
    SD=1;
    h=f;
    
    if isempty(M) || length(M)<countIMFs
        
        if isa(options.alpha,'char')
            if strcmp(options.alpha,'ave') % Using an average mask length
                m = 2*round(N_pp/k_pp*options.Xi);
            elseif strcmp(options.alpha,'Almost_min') % Using an almost min mask length
                if 2*round(options.Xi*prctile(diffMaxmins_pp,30))<2*round(N_pp/k_pp*options.Xi)
                    m = 2*round(options.Xi*prctile(diffMaxmins_pp,30));
                else
                    m = 2*round(N_pp/k_pp*options.Xi);
                end
            else
                disp(' Value of alpha not recognized')
                return
            end
        else % using a fixed value alpha
            m = 2*round(options.Xi*prctile(diffMaxmins_pp,options.alpha));
        end
        if countIMFs>1
            if m<=stats(countIMFs-1).logM
                if options.verbose>0
                    fprintf('Warning mask length is decreasing at step %1d. ',countIMFs)
                end
                if options.MonotoneMaskLength==true
                    m=ceil(stats(countIMFs-1).logM*1.1);
                    if options.verbose>0
                        fprintf('The old mask length is %1d whereas the new one is forced to be %1d.\n',stats(countIMFs-1).logM,ceil(stats(countIMFs-1).logM*1.1))
                    end
                else
                    if options.verbose>0
                        fprintf('The old mask length is %1d whereas the new one is %1d.\n',stats(countIMFs-1).logM,m)
                    end
                end
            end 
        end
    else
        m=M(countIMFs);
    end
    
    inStepN=0;
    if options.verbose>0
        fprintf('\n IMF # %1.0d   -   # Extreme points %5.0d\n',countIMFs,min(k_pp))
        fprintf('\n  step #            SD             Mask length\n\n')
    end
    
    stats(countIMFs).logM=m;
    a = get_mask_v1_1(MM,m,options.verbose);
    ExtendSig=1==0;
    if N < length(a) % we need to extend the signal
        ExtendSig=1==1;
        Nxs=ceil(length(a)/N);
        N_old=N;
        if rem(Nxs,2)==0
            Nxs=Nxs+1;
        end
        h_n=[];
        for ii=1:Nxs
            h_n=[h_n; h];
        end
        h=h_n;
        N=Nxs*N;
    end
    
    Nza=N-length(a);
    if rem(Nza,2)==0
        a = [zeros(1,Nza/2) a zeros(1,Nza/2)];
        fftA=real(fft([a((length(a)-1)/2+1:end) a(1:(length(a)-1)/2)]));
        % figure,plot(circshift(a,(length(a)-1)/2+1)-ifft(real(fft(circshift(a,(length(a)-1)/2+1)))),'r')
    else
        a = [zeros(1,(Nza-1)/2) a zeros(1,(Nza-1)/2+1)];
        %csA=circshift(a,(length(a))/2+1);
        fftA=real(fft([a((length(a))/2:end) a(1:(length(a))/2-1)]));
        % figure,plot(circshift(a,(length(a))/2+1)-ifft(real(fft(circshift(a,(length(a))/2+1)))),'r')
    end
    if options.plots>0 %&& rem(inStepN,5)==0
        if gcf > 30
            close all
        end
        figN=figure;
        set(figN,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    end
    fftH=fft(h);
    %fft_h_new=fftH;
        
    %% we compensate for the aliasing effect in the DFT of the filter
    % we look for the first minimum in fftA
    
    posF=find((diff(fftA)>0)~=0,1,'first');
    
    stats(countIMFs).posF = posF;
    stats(countIMFs).valF = fftA(posF);
    
    fftA=fftA-fftA(posF);
    fftA(fftA<0)=0;
    
    %plot([ifftA(end/2+2:end) ifftA(1:end/2+1)])
    if options.plots>=1
        figMask=figure;
        figRem=figure;
        set(figMask,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    end
    
    while SD>options.delta && inStepN < options.MaxInner
        inStepN=inStepN+options.NumSteps;
                
        fft_h_old=(1-fftA').^(inStepN-1).*fftH;
        fft_h_new=(1-fftA').^inStepN.*fftH;
        
        %%%%%%%%%%%%%%%% Updating stopping criterium %%%%%%%%%%%%%%%%%
        
        for i=1:D
            SD(i)=norm(fft_h_new(:,i)-fft_h_old(:,i))^2/norm(fft_h_old(:,i))^2;            
        end
        SD=max(SD);        
        
        if options.verbose>0
            fprintf('    %2.0d      %1.14f          %2.0d\n',inStepN,SD,m)
        end
    end
    
    h=ifft(fft_h_new);
    
    if ExtendSig % we reduce the signal
        N=N_old;
        h=h(N*(Nxs-1)/2+1:N*((Nxs-1)/2+1),:);
    end
    if inStepN >= options.MaxInner && options.verbose>0
        disp('Max # of inner steps reached')
        %return
    end
    stats(countIMFs).inStepN=inStepN;
    for i=1:D
        IMF{i}(:,countIMFs) = h(:,i);
    end
    f=f-h;
    
    %% Create a signal without zero regions and compute the number of extrema
    g=f';
    f_pp=real(acos(dot(normc(g(:,1:end-1)),normc(g(:,2:end)))));
    f_pp(f_pp<=10^-18)=[];
    if isempty(f_pp)
        break
    end    
    maxmins_pp=Maxmins_v3_7(f_pp);
    if isempty(maxmins_pp)
        break
    end
    diffMaxmins_pp=diff(maxmins_pp);
    N_pp=length(f_pp);
    k_pp = length(maxmins_pp);
    
    if options.saveInter==1
        save([nameFile '_intermediate_MvFIF_v7.mat'],'IMF','f','stats','-v7.3');
    end
end %end of while
for i=1:D
    IMF{i} = [IMF{i}(:,1:countIMFs) f(:,i)];
end

if options.saveEnd == 1
    save([ 'Final_' nameFile '_MvFIF_v7.mat'],'IMF','stats','-v7.3');
end

end


%% Auxiliar functions


function a=get_mask_v1_1(y,k,verbose)
%
% Rescale the mask y so that its length becomes 2*k+1.
% k could be an integer or not an integer.
% y is the area under the curve for each bar

n=length(y);
m=(n-1)/2;

if k<=m % The prefixed filter contains enough points
    
    if mod(k,1)==0     % if the mask_length is an integer
        
        a=zeros(1,2*k+1);
        
        for i=1:2*k+1
            s=(i-1)*(2*m+1)/(2*k+1)+1;
            t=i*(2*m+1)/(2*k+1);
            
            %s1=s-floor(s);
            s2=ceil(s)-s;
            
            t1=t-floor(t);
            %t2=ceil(t)-t;
            
            if floor(t)<1
                disp('Ops')
            end
            a(i)=sum(y(ceil(s):floor(t)))+s2*y(ceil(s))+t1*y(floor(t));
        end
        
    else   % if the mask length is not an integer
        new_k=floor(k);
        extra = k-new_k;
        c=(2*m+1)/(2*new_k+1+2*extra);
        
        a=zeros(1,2*new_k+3);
        
        t=extra*c+1;
        t1=t-floor(t);
        %t2=ceil(t)-t;
        if k<0
            disp('Ops')
            a=[];
            return
        end
        a(1)=sum(y(1:floor(t)))+t1*y(floor(t));
        
        for i=2:2*new_k+2
            s=extra*c+(i-2)*c+1;
            t=extra*c+(i-1)*c;
            %s1=s-floor(s);
            s2=ceil(s)-s;
            
            t1=t-floor(t);
            
            
            a(i)=sum(y(ceil(s):floor(t)))+s2*y(ceil(s))+t1*y(floor(t));
        end
        t2=ceil(t)-t;
        
        a(2*new_k+3)=sum(y(ceil(t):n))+t2*y(ceil(t));
    end
else % We need a filter with more points than MM, we use interpolation
    dx=0.01;
    % we assume that MM has a dx = 0.01, if m = 6200 it correspond to a
    % filter of length 62*2 in the physical space
    f=y/dx; % function we need to interpolate
    dy=m*dx/k;
    b=interp1(0:m,f(m+1:2*m+1),0:m/k:m);
    if size(b,1)>size(b,2)
        b=b.';
    end
    if size(b,1)>1
        fprintf('\n\nError!')
        disp('The provided mask is not a vector!!')
        a=[];
        return
    end
    a=[fliplr(b(2:end)) b]*dy;
    if abs(norm(a,1)-1)>10^-14
        if verbose>0
        fprintf('\n\n Warning!\n\n')
        fprintf(' Area under the mask equals %2.20f\n',norm(a,1))
        fprintf(' it should be equal to 1\n We rescale it using its norm 1\n\n')
        end
        a=a/norm(a,1);
    end
end

end


function varargout = Maxmins_v3_7(f)

% v3_7 Based on version 3_6
% Minor revision: we extend f as df in order to make possible to compute
% the relative derivative like df(i)/abs(f(i)) for every i=1:N
%
% v3_6 Based on version 3_5
% Minor revision: we consider relative differences to avoid having problems with small or big values in a signal
%
% v3_5 Based on version 3_4
% Minor revision: we consider only periodical extension
%
% v3_4 Based on version 3
% Minor revisions: 1) added for constant extention the checking for Mins and
%                     Maxs emptiness
%                  2) completed the code for the periodical case
%
% v3 is Based on Version 2.
% Modified the way zero-derivative regions are handled.
%
% Identify the maxima and minima of a signal f

tol=10^-15;

N = length(f);
Maxs = zeros(1,N);
Mins = zeros(1,N);
df = diff(f);


h = 1;

    while h<N && abs(df(h)/f(h)) <= tol        
        h=h+1;
    end    
    if h==N
        if nargout<=1
            varargout{1}=[];
        elseif nargout==2
            varargout{1}=[];
            varargout{2}=[];
        end
        return
    end

cmaxs=0;
cmins=0;

c = 0;

N_old=N;

df=diff([f f(2:h+1)]);
f=[f f(2:h)];
N=N+h;


last_df=[];
for i=h:N-2
    if   df(i)*df(i+1)/abs(f(i))^2 <= tol && df(i)*df(i+1)/abs(f(i))^2 >= -tol
        if df(i)/abs(f(i)) < -tol
            last_df=-1;
            posc = i;
        elseif df(i)/abs(f(i)) > tol
            last_df=+1;
            posc = i;
        end
        c = c + 1;
        if df(i+1)/abs(f(i)) < -tol
            if last_df==+1
                cmaxs=cmaxs+1;
                Maxs(cmaxs)=mod(posc+floor((c-1)/2)+1,N_old);
            end
            c=0;
        end
        if df(i+1)/abs(f(i)) > tol
            if last_df==-1
                cmins=cmins+1;
                Mins(cmins)=mod(posc+floor((c-1)/2)+1,N_old);
            end
            c=0;
        end
        
    end
    if   df(i)*df(i+1)/abs(f(i))^2 < -tol
        if df(i)/abs(f(i)) < -tol && df(i+1)/abs(f(i)) > tol
            cmins=cmins+1;
            Mins(cmins)=mod(i+1,N_old);
            if Mins(cmins)==0
                Mins(cmins)=1;
            end
            last_df=-1;
        elseif df(i)/abs(f(i)) > tol && df(i+1)/abs(f(i)) < -tol
            cmaxs=cmaxs+1;
            Maxs(cmaxs)=mod(i+1,N_old);
            if Maxs(cmaxs)==0
                Maxs(cmaxs)=1;
            end
            last_df=+1;
        end
    end
end
if c > 0    
    %         % we deal with the boundary
    %         df_0=f(N)-f(1);
    %         if df_0==0
    %             if Initial_df < 0
    %                 if last_df==+1
    %                     cmaxs=cmaxs+1;
    %                     Maxs(cmaxs)=mod(posc+floor((c+cIn-1)/2)+1,N);
    %                 end
    %             elseif Initial_df > 0
    %                 if last_df==-1
    %                     cmins=cmins+1;
    %                     Mins(cmins)=mod(posc+floor((c+cIn-1)/2)+1,N);
    %                 end
    %             end
    %         else
    %             disp('Code missing!')
    %         end 
    if cmins>0 && Mins(cmins)==0 
        Mins(cmins)=N;
    end
    if cmaxs>0 && Maxs(cmaxs)==0
        Maxs(cmaxs)=N;
    end
end

Maxs=Maxs(1:cmaxs);
Mins=Mins(1:cmins);
maxmins=sort([Maxs Mins]);
%     disp('Code to be completed')
%     if isempty(maxmins)
%         maxmins = 1;
%     else
%         if maxmins(1)~=1 && maxmins(end)~=N
%             if (f(maxmins(end)) > f(maxmins(end)+1) && f(maxmins(1)) > f(maxmins(1)-1)) || (f(maxmins(end)) < f(maxmins(end)+1) && f(maxmins(1)) < f(maxmins(1)-1))
%                 maxmins=[1 maxmins];
%             end
%         end
%     end

if sum(Maxs==0)>0
    Maxs(Maxs==0)=1;
end
if sum(Mins==0)>0
    Mins(Mins==0)=1;
end
if sum(maxmins==0)>0
    maxmins(maxmins==0)=1;
end

if nargout<=1
    varargout{1}=maxmins;
elseif nargout==2
    varargout{1}=Maxs;
    varargout{2}=Mins;
end

end