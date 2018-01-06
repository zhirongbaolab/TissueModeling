function saveObjFile(name, v,f)
fid = fopen(name,'w');
    
    for i=1:size(v,1)
        fprintf(fid,'v %f %f %f\n',v(i,1),v(i,2),v(i,3));
    end
    
    fprintf(fid,'g foo\n');
    
    for i=1:size(f,1);
        fprintf(fid,'f %d %d %d\n',f(i,1),f(i,2),f(i,3));
    end
    fprintf(fid,'g\n');
    
    fclose(fid);
end