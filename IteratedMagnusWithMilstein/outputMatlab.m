function outputMatlab(fileName,...
                      T,dtEulerRef,dtEuler,dtMagnus,M,kappas,nx,nv,...
                      ctimeEulerRefTotal,ctimeEulerTotal,ctimeMagnus1Total,ctimeMagnus2Total,ctimeMagnus3Total,...
                      errDisExactEulerRef,errDisExactEuler,errDisExactMagnus1,errDisExactMagnus2,errDisExactMagnus3,...
                      errDisEulerRefEuler,errDisEulerRefMagnus1,errDisEulerRefMagnus2,errDisEulerRefMagnus3,...
                      errAbsAvgExactEulerRef,errAbsAvgExactEuler,errAbsAvgExactMagnus1,errAbsAvgExactMagnus2,errAbsAvgExactMagnus3,...
                      errAbsAvgEulerRefEuler,errAbsAvgEulerRefMagnus1,errAbsAvgEulerRefMagnus2,errAbsAvgEulerRefMagnus3,...
                      hasExact,useExact,useEulerRef,useEuler,useMagnus1,useMagnus2,useMagnus3)
root=[pwd, '/' ,'Results'];
matRoot=[root,'/','Mat','/',fileName];
delDir(matRoot)
mkDir(matRoot)
Result=struct;
if useEulerRef
    Result.EulerRef=struct;
    Result.EulerRef=addMethod(Result.EulerRef,ctimeEulerRefTotal,dtEulerRef);
end
if useEuler
    Result.Euler=struct;
    Result.Euler=addMethod(Result.Euler,ctimeEulerTotal,dtEuler);
end
if useMagnus1
    Result.Magnus1=struct;
    Result.Magnus1=addMethod(Result.Magnus1,ctimeMagnus1Total,dtMagnus);
end
if useMagnus2
    Result.Magnus2=struct;
    Result.Magnus2=addMethod(Result.Magnus2,ctimeMagnus2Total,dtMagnus);
end
if useMagnus3
    Result.Magnus3=struct;
    Result.Magnus3=addMethod(Result.Magnus3,ctimeMagnus3Total,dtMagnus);
end
if useExact && hasExact
    if useEulerRef
        Result.EulerRef=addAbsError(Result.EulerRef,'Exact',errAbsAvgExactEulerRef{end});
        Result.EulerRef=addRelError(Result.EulerRef,'Exact',errDisExactEulerRef{end});
    end
    if useEuler
        Result.Euler=addAbsError(Result.Euler,'Exact',errAbsAvgExactEuler{end});
        Result.Euler=addRelError(Result.Euler,'Exact',errDisExactEuler{end});
    end
    if useMagnus1
        Result.Magnus1=addAbsError(Result.Magnus1,'Exact',errAbsAvgExactMagnus1{end});
        Result.Magnus1=addRelError(Result.Magnus1,'Exact',errDisExactMagnus1{end});
    end
    if useMagnus2
        Result.Magnus2=addAbsError(Result.Magnus2,'Exact',errAbsAvgExactMagnus2{end});
        Result.Magnus2=addRelError(Result.Magnus2,'Exact',errDisExactMagnus2{end});
    end
    if useMagnus3
        Result.Magnus3=addAbsError(Result.Magnus3,'Exact',errAbsAvgExactMagnus3{end});
        Result.Magnus3=addRelError(Result.Magnus3,'Exact',errDisExactMagnus3{end});
    end
end
if useEulerRef 
    if useEuler
        Result.Euler=addAbsError(Result.Euler,'EulerRef',errAbsAvgEulerRefEuler{end});
        Result.Euler=addRelError(Result.Euler,'EulerRef',errDisEulerRefEuler{end});
    end
    if useMagnus1
        Result.Magnus1=addAbsError(Result.Magnus1,'EulerRef',errAbsAvgEulerRefMagnus1{end});
        Result.Magnus1=addRelError(Result.Magnus1,'EulerRef',errDisEulerRefMagnus1{end});
    end
    if useMagnus2
        Result.Magnus2=addAbsError(Result.Magnus2,'EulerRef',errAbsAvgEulerRefMagnus2{end});
        Result.Magnus2=addRelError(Result.Magnus2,'EulerRef',errDisEulerRefMagnus2{end});
    end
    if useMagnus3
        Result.Magnus3=addAbsError(Result.Magnus3,'EulerRef',errAbsAvgEulerRefMagnus3{end});
        Result.Magnus3=addRelError(Result.Magnus3,'EulerRef',errDisEulerRefMagnus3{end});
    end
end
save([matRoot,'/',fileName,'.mat'],'Result')
    function sct=addMethod(sct,ctime,dt)
        sct.ctime=ctime;
        sct.dt=dt;
        sct.kappas=kappas;
        sct.nx=nx;
        sct.nv=nv;
        sct.T=T;
        sct.M=M;
    end
    function sct=addAbsError(sct,ref,absError)
        sct.(ref).absError=absError;
    end
    function sct=addRelError(sct,ref,relError)
        temp=zeros(1,length(relError));
        for it1 = 1:1:length(temp)
            temp(it1)=mean(relError{it1},1);
        end
        sct.(ref).relError=temp;
    end
end
function delFile(file)
    if exist(file)
        delete(file);
    end
end
function delDir(dir)
    if exist(dir)==7
        rmdir(dir,'s');
    end
end
function cleanDir(mdir,except)
    except{end+1}='.';
    except{end+1}='..';
    for d = dir(mdir).'
      if ~any(strcmp(d.name,except))
          if d.isdir
              rmdir([d.folder,'/',d.name],'s');
          else
              delete([d.folder,'/',d.name]);
          end
      end
    end
end
function mkDir(dir)
    if exist(dir)==0
        mkdir(dir);
    end
end