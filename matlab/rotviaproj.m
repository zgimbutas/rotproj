function mpout=rotviaproj(nterms,mpole,beta,alpha)
%ROTVIAPROJ Rotate the complex spherical harmonics expansion about the y-axis.
%
%  MPOUT = rotviaproj(NTERMS,MPOLE,BETA) rotates the complex spherical 
%  harmonics expansion of degree NTERMS about the y-axis by degree BETA. 
%
%  MPOUT = rotviaproj(NTERMS,MPOLE,BETA,ALPHA) rotates the complex spherical 
%  harmonics expansion of degree NTERMS about the x-axis by degree ALPHA 
%  and about the y-axis by degree BETA.
%
%  After rotation, the expansion pole is moved to location (beta, alpha) 
%  in spherical coordinates (theta, phi).  
%
%  Both MPOUT and MPOLE are (NTERMS+1)-by-(2*NTERMS+1) complex matrices.
%
%
%  Fast and stable algorithm for applying rotation operator about
%  the y-axis determined by angle beta.
%
%  The method is based on computing the induced potential and
%  its theta-derivative on the rotated equator
%  for each order (first index). The coefficients of  the rotated
%  expansion can then be obtained by FFT and projection.
%
%  There is some loss in speed over using recurrence relations 
%  but it is stable to all orders whereas the recurrence schemes 
%  are not.
%
%  Our definition of complex spherical harmonics is
%
%  Ynm(theta,phi)= sqrt( 2n+1) sqrt((n-m)!/(n+m)!) 
%                  Pnm(cos theta) e^(im phi), 
%  Yn,-m(theta,phi) = sqrt( 2n+1) sqrt((n-m)!/(n+m)!) 
%                  Pnm(cos theta) e^(-im phi),   for m >= 0.
%       
%  Note that we do not include the Condon-Shortley phase (-1)^m, if m<0.
%

if( nargin == 0 ), mpout=rotviaproj_test(); return; end;

m1=nterms;
m2=nterms;

mpout = zeros(nterms+1,2*nterms+1)+1i*zeros(nterms+1,2*nterms+1);

lmp = nterms;
lmpn = nterms;

if( nargin == 3 ),
mex_id_ = 'rotviaprojf90(i double[x], i int[x], i int[x], i int[x], io dcomplex[], i int[x], io dcomplex[], i int[x])';
[mpole, mpout] = rotproj_r2012a(mex_id_, beta, nterms, m1, m2, mpole, lmp, mpout, lmpn, 1, 1, 1, 1, 1, 1);
end

if( nargin >= 4 ),
mex_id_ = 'rotviaproj2f90(i double[x], i double[x], i int[x], i int[x], i int[x], io dcomplex[], i int[x], io dcomplex[], i int[x])';
[mpole, mpout] = rotproj_r2012a(mex_id_, beta, alpha, nterms, m1, m2, mpole, lmp, mpout, lmpn, 1, 1, 1, 1, 1, 1, 1);
end

function mpout=rotviaproj_test()

nterms = 2;
beta = pi/2;
nrot = 2;

mpole = zeros(nterms+1,2*nterms+1)+1i*zeros(nterms+1,2*nterms+1);

j=nterms+1;
mpole(1,j) = 1;
mpole(2,j+(-1:1)) = 1; 
mpole(3,j+(-2:2)) = 1; 

mpout = rotviaproj(nterms,mpole,beta);

