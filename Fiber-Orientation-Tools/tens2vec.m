function Vec = tens2vec(Tens)
%VEC = TENS2VEC(TENS) Convert a symmetric 2nd order tensor TENS, stored in 
%     3x3 matrix form, to VEC in 6x1 column vector (contracted) form.
Vec = [Tens(1,1), Tens(2,2), Tens(3,3), Tens(2,3), Tens(3,1), Tens(1,2)]';