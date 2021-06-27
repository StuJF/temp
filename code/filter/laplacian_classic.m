function LAP=laplacian_classic(de_chanlocs,num_near)
% laplacian spatial filter
% Input:
%   num_near: number of neigbour electrode,defalut is 4
%   de_chanlocs: channel location file, either use standard eeglab channel
%       location file or specify channel position in a matrix of nchanXndim
%       nchan is number of channel,ndim is dimension of distance,e.g.
%       2->x,y;3->x,y,z
if nargin<2
    num_near=4;
end
if nargin<1
    error('no chanlocs input');
end

if isstruct(de_chanlocs)
    num_chan=length(de_chanlocs);
    
    dist=zeros(num_chan);
    for n=1:num_chan %num_chan-1
        for m=(n+1):num_chan
            xi=de_chanlocs(n).X;
            yi=de_chanlocs(n).Y;
            zi=de_chanlocs(n).Z;
            xj=de_chanlocs(m).X;
            yj=de_chanlocs(m).Y;
            zj=de_chanlocs(m).Z;
            dist(n,m)=sqrt((xi-xj)^2+(yi-yj)^2+(zi-zj)^2);
        end
    end
    dist=dist+dist';
else
    try
        num_chan=size(de_chanlocs,1);
        dist=zeros(num_chan);
        for n=1:num_chan
            for m=(n+1):num_chan
                dist(n,m)=sqrt(sum((de_chanlocs(n,:)-de_chanlocs(m,:)).^2));
            end
        end
        dist=dist+dist';
    catch
        error('undefined chanloc information')
    end
end


if num_near>=num_chan
    error('number of near electrode > total number of electrode');
else
    near=zeros(num_chan);%确定每个电极的邻近电极，1表示对i来说j为它的邻近电极
    [~,ind]=sort(dist,2);
    nearby=ind(:,2:num_near+1);
    for n=1:num_chan
        near(n,nearby(n,:))=1;
    end
    distc=dist+eye(num_chan);%对角线元素取1只是为了后续计算方便，实际值应为0
end

%G=zeros(num_chan);
S=sum(near./distc,2);

G=bsxfun(@rdivide,near./distc,S);

LAP=eye(num_chan)-G;



end