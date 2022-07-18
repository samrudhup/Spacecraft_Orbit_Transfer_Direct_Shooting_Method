function grd = OrbitTransferGrd(Z)
% computes the gradient

output = OrbitTransferObj_Jac(Z);
grd    = output;

end

