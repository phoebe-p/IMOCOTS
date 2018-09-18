function parsave(fname, mat)
mat = sparse(mat);
save(strcat('results/', fname), 'mat')
end