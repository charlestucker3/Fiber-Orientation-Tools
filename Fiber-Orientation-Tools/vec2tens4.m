function A4 = vec2tens4(v4)
%A4 = VEC2TENS4(V4) converts a columnt vector of the independent components of
%     a fully symmetric fourth-order tensor to matrix (contracted) form.
%     For 3-D tensors, V4 is 15x1 and A4 is 6x6.  For 2-D tensors, V4 is
%     5x1 and A4 is 3x3.  The components in V4 are arranged as in
%     TENS2VEC4.
%
%  See also: TENS2VEC4

% Determine whether v4 is 2-D or 3-D (by its size),
% and expand accordingly
[rows, cols] = size(v4);
if (rows == 15 && cols == 1) || (rows == 1 && cols == 15) 
    % 3-D tensor
    A4 = [v4( 1), v4( 6), v4( 5), v4( 7), v4( 8), v4( 9); ...
          v4( 6), v4( 2), v4( 4), v4(10), v4(11), v4(12); ...
          v4( 5), v4( 4), v4( 3), v4(13), v4(14), v4(15); ...
          v4( 7), v4(10), v4(13), v4( 4), v4(15), v4(11); ...
          v4( 8), v4(11), v4(14), v4(15), v4( 5), v4( 7); ...
          v4( 9), v4(12), v4(15), v4(11), v4( 7), v4( 6)];
  
elseif (rows == 5 && cols == 1) || (rows == 1 && cols == 5) 
    % 2-D tensor
    A4 = [v4( 1), v4( 3), v4( 4); ...
          v4( 3), v4( 2), v4( 5); ...
          v4( 4), v4( 5), v4( 3)];
else
    error('v4 must be 15x1 or 5x1')

end


  
end