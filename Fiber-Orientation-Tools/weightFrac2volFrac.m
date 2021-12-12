function vf = weightFrac2volFrac(wf, rhof, rhom)
%VF = WEIGHTFRAC2VOLFRAC(WF, RHOF, RHOM) returns the volume fraction of fibers VF in a 
%     two-phase composite with fiber weight fraction WF, fiber density RHOF
%     and matrix density RHOM.  0 <= WF <= 1.

vf = wf ./ (wf + (rhof/rhom)*(1-wf));

return