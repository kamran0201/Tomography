classdef CSHelperCode
    properties
        transpose=0;
        theta=[];
        numberOfAngles=0;
        n=0;
        m=0;
        size=0;
    end
    methods
       function obj = CSHelperCode(m,n,theta)
            obj.m=m;
            obj.n=n;
            obj.size=round(sqrt(n));
            obj.theta=theta;
            obj.numberOfAngles=numel(theta);
       end
       function res = mtimes(obj,beta)
            angles=obj.theta;
            if obj.transpose == 0
                dctCoeffMat=reshape(beta,obj.size,obj.size);
                x=idct2(dctCoeffMat);
                Ax=radon(x,angles);
                res=reshape(Ax,obj.m*obj.numberOfAngles,1);
            else
                projection=reshape(beta,obj.m,obj.numberOfAngles);
                Atx=dct2(iradon(projection,angles,'linear','Ram-Lak',1,obj.size));
                res=reshape(Atx,obj.n,1);
            end
       end
       function  res = ctranspose(obj)
           obj.transpose=xor(obj.transpose,1);
           res=obj;
       end
    end
end