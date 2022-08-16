classdef CoupledCSHelperCode
    properties
        transpose=0;
        coupleSize=1;
        theta={};
        noOfAngle=0;
        n=0;
        m=0;
        size=0;
        mapMtx=[];
    end
    methods
       function obj = CoupledCSHelperCode(m,n,coupleSize,theta)
            obj.m=m;
            obj.noOfAngle=numel(theta{1});
            obj.n=n;
            obj.size=round(sqrt(n));
            obj.theta=theta;
            obj.coupleSize=coupleSize;
       end
       function res = mtimes(obj,beta)
            if obj.transpose == 0
                predictY=[];
                betaCell = getBetaMtxFromVec(obj,beta);
                for i=1:obj.coupleSize
                    betaSum=zeros(obj.size,obj.size);
                    for j=1:i
                        betaSum=betaSum+betaCell{j};
                    end
                    angles=obj.theta{i};
                    dctCoffMat=betaSum;
                    x=idct2(dctCoffMat);
                    Ax=radon(x,angles);
                    predictY=horzcat(predictY,Ax);
                end
                res=reshape(predictY,obj.m*obj.noOfAngle*obj.coupleSize,1);
            else
                yCell = getYMtxFromVec(obj,beta);
                predictBeta=[];
                for i=1:obj.coupleSize
                    pSum=zeros(obj.size,obj.size);
                    for j=i:obj.coupleSize
                        angles=obj.theta{j};
                        projectionMatrix=yCell{j};
                        Atx=dct2(iradon(projectionMatrix,angles,'linear','Ram-Lak',1,obj.size));
                        pSum=pSum+Atx;
                    end
                    predictBeta=horzcat(predictBeta,Atx);
                end
                res=reshape(predictBeta,obj.n*obj.coupleSize,1);
            end
       end
       function  res = ctranspose(obj)
           obj.transpose=xor(obj.transpose,1);
           res=obj;
       end
       function betaCell = getBetaMtxFromVec(obj,vecBeta)
             betaCell=cell(obj.coupleSize,1);
             vsize=obj.n;
             for i=1:obj.coupleSize
                offset=(i-1)*vsize;
                tbeta=vecBeta(offset+1:offset+vsize,1);
                dctCoffMat=reshape(tbeta,obj.size,obj.size);
                betaCell{i}=dctCoffMat;
             end
       end
       function yCell = getYMtxFromVec(obj,vecY)
             yCell=cell(obj.coupleSize,1);
             vsize=obj.m*obj.noOfAngle;
             for i=1:obj.coupleSize
                offset=(i-1)*vsize;
                y=vecY(offset+1:offset+vsize,1);
                proj=reshape(y,obj.m,obj.noOfAngle);
                yCell{i}=proj;
             end
       end
    end    
end