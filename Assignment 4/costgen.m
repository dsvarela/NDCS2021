function [H,h]=costgen(T,S,dim,LTI)

H = S'*S + eye(dim.N*dim.nu);
h = 2*(S'*T*LTI.x0)';
 
end
