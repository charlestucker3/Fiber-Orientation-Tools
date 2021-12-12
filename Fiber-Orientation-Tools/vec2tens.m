function Tens = vec2tens(Vec)
% TENS = VEC2TENS(VEC) converts a 6x1 contracted symmetry tensor VEC to 
%       a 3x3 matrix form TENS.
Tens = [Vec(1), Vec(6), Vec(5);
        Vec(6), Vec(2), Vec(4);
        Vec(5), Vec(4), Vec(3)];