function signAA = huberAbs(A,mu)


   signA1  = sign(A).*(abs(A)>mu);
    signA2 =(abs(A)<mu).*A/mu;
     signAA = signA1+signA2;
