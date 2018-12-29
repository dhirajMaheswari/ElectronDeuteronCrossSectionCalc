      subroutine struct_bodek_p(gp0,gp2,fm2,f2,f1)        
*******************************************************
* Structure function f2, f2 of protons within Bodek   *
* parameterization                                    *
*******************************************************
      PM  = 0.938279   
      R   = 0.18
*      gp0 = gp2/2.0/pm/x                                              
*      FM2 = PM**2+2.*PM*GP0-GP2                                           
      W2H = 0.                                                            
      IF(FM2.LT.PM**2) RETURN                                           
      WI  = SQRT(FM2)                                                      
      W2H = GP_H(GP0,GP2)*B(WI,GP2)/GP0  
      W1H = (1.+GP0**2/GP2)/(1.+R)*W2H
      f2  = gp0 * w2h                           
      f1  = pm  * w1h                                   
      RETURN                                                            
      END                                                               
C     .........................                                         
      FUNCTION GP_H(Q0,Q2)                                                
      PM = 0.938279                                        
      XX = Q2/(2.*PM*Q0)                                         
      GI = 2.*PM*Q0                                              
      WW = (GI+1.642)/(Q2+0.376)                                 
      T  = (1.-1./WW)                                             
      WP = 0.256*T**3+2.178*T**4+0.898*T**5-6.716*T**6+3.756*T**7
      GP_H=WP*WW*Q2/(2.*PM*Q0) 
      RETURN                                                            
      END                                                               
C---------------------------------                                      





C..................................             
      subroutine struct_bodek_n(gp0,gp2,fm2,f2,f1)         
*******************************************************
* Structure function f2, f2 of neutrons within Bodek  *
* parameterization                                    *
*******************************************************
      common/ths_ths/ths0
      PM  = 0.938279     
      R   = 0.18            
*      gp0 = gp2/2.0/pm/x                                             
*      FM2=PM**2+2.*PM*GP0-GP2                                           
      W2NT=0.                                                            
      IF(FM2.LT.PM**2) RETURN                                           
      WI=SQRT(FM2)                                                      
      W2NT =GP_N(GP0,GP2)*B(WI,GP2)/GP0       
      W1NT = (1.+GP0**2/GP2)/(1.+R)*W2NT
      f2  = gp0 * w2NT                           
      f1  = pm  * w1NT       
*      write(15,*)"bodek",ths0,f1,f2,gp0,gp2,wi
      RETURN                                                            
      END                                                               
C     .........................                                         
      FUNCTION GP_N(Q0,Q2)                                                
      PM = 0.938279                                        
      XX = Q2/(2.*PM*Q0)                                         
      GI = 2.*PM*Q0                                              
      WW = (GI+1.642)/(Q2+0.376)                                 
      T  = (1.-1./WW)                                             
      WN = 0.064*T**3+0.225*T**4+4.106*T**5-7.079*T**6+3.055*T**7    
      GP_N=WN*WW*Q2/(2.*PM*Q0) 
      RETURN                                                            
      END                                                               
C---------------------------------                                      

 
C  *******************************                                       
C  *    BODEK PARAMETERIZATION    *                                       
C  *******************************                                       
      FUNCTION B(WM,QSQ)                                                 
      DIMENSION C(24)                                                    
      INTEGER LSPIN(4)                                                   
      DATA PMSQ/0.880324/,PM2/1.876512/,PM/0.938256/                     
      DATA NRES/4/,NBKG/5/                                               
      DATA LSPIN/1,2,3,2/                                                
      DATA C/1.0741163,0.75531124,3.3506491,1.7447015,3.5102405,1.040004 
     *,1.2299128,0.10625394,0.48132786,1.5101467,0.081661975,0.65587179, 
     *1.7176216,0.12551987,0.7473379,1.953819,0.19891522,-0.17498537,    
     *0.0096701919,-0.035256748,3.5185207,-0.59993696,4.7615828,0.411675  
     *89/                                                                
      B=0.                                                               
      IF(WM.LE.0.939)RETURN                                               
      WSQ=WM**2                                                          
      OMEGA=1.+(WSQ-PMSQ)/QSQ                                            
      X=1./OMEGA                                                         
      XPX=C(22)+C(23)*(X-C(24))**2                                       
      PIEMSQ=(C(1)-PM)**2                                                
********************************************************
*     added part
********************************************************
      B1 = 0.0
      IF(WM.EQ.C(1))GOTO 11 
******************************************************** 
      B1=AMAX1(0.,(WM-C(1)))/(WM-C(1))*C(2)       ! 0/0                  
********************************************************
11    EB1=C(3)*(WM-C(1))                                                 
      IF(EB1.GT.25.)GO TO 1                                              
      B1=B1*(1.0-EXP(-EB1))                                              
*********************************************************
*     added part
*********************************************************  
      B2 = 0.0
      IF(WM.EQ.C(4))GOTO 12  
*********************************************************
1     B2=AMAX1(0.,(WM-C(4)))/(WM-C(4))*(1.-C(2))   ! 0/0                 
*********************************************************
12    EB2=C(5)*(WSQ-C(4)**2)                                             
      IF(EB2.GT.25.) GO TO 2                                             
      B2=B2*(1.-EXP(-EB2))                                               
2     CONTINUE                                                           
      BBKG=B1+B2                                                         
      BRES=C(2)+B2                                                       
      RESSUM=0.                                                          
      DO 30 I=1,NRES                                                     
      INDEX=(I-1)*3+1+NBKG                                               
      RAM=C(INDEX)                                                       
      IF(I.EQ.1)RAM=C(INDEX)+C(18)*QSQ+C(19)*QSQ**2                      
      RMA=C(INDEX+1)                                                     
      IF(I.EQ.3)RMA=RMA*(1.+C(20)/(1.+C(21)*QSQ))                        
      RWD=C(INDEX+2)                                                     
      QSTARN=SQRT(AMAX1(0.,((WSQ+PMSQ-PIEMSQ)/(2.*WM))**2-PMSQ))         
      QSTARO=SQRT(AMAX1(0.,((RMA**2-PMSQ+PIEMSQ)/(2.*RMA))**2-PIEMSQ))   
      IF(QSTARO.LE.1.E-10)GO TO 40                                       
      TERM=6.08974*QSTARN                                                
      TERMO=6.08974*QSTARO                                               
      J=2*LSPIN(I)                                                       
      K=J+1                                                              
      GAMRES=RWD*(TERM/TERMO)**K*(1.+TERMO**J)/(1.+TERM**J)              
      GAMRES=GAMRES/2.                                                   
      BRWIG=GAMRES/((WM-RMA)**2+GAMRES**2)/3.1415926                     
      RES=RAM*BRWIG/PM2                                                  
      GO TO 30                                                           
40    RES=0.                                                             
30    RESSUM=RESSUM+RES                                                  
      B=BBKG*(1.+(1.-BBKG)*XPX)+RESSUM*(1.-BRES)                         
      RETURN                                                             
      END                                                                

