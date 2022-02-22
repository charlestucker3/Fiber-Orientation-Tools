function A4 = vec2tens4(v4)
%A4 = VEC2TENS4(V4) converts a 15x1 vector of the independent components of
%     a fully symmetric fourth-order tensor to 6x6 matrix (contracted)
%     form.  The components in V4 are arranged as in TENS2VEC4.
%
%  See also: TENS2VEC4

A4 = [v4( 1), v4( 6), v4( 5), v4( 7), v4( 8), v4( 9); ...
      v4( 6), v4( 2), v4( 4), v4(10), v4(11), v4(12); ...
      v4( 5), v4( 4), v4( 3), v4(13), v4(14), v4(15); ...
      v4( 7), v4(10), v4(13), v4( 4), v4(15), v4(11); ...
      v4( 8), v4(11), v4(14), v4(15), v4( 5), v4( 7); ...
      v4( 9), v4(12), v4(15), v4(11), v4( 7), v4( 6)];
  
end