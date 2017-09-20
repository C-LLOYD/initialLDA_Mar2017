from numpy import *
from numpy.fft import *

def auto(t,u,w,atw,Tc,Fs,meanrm,locnor=False):
#
# for further information see
# http://ldvproc.nambis.de/programs/pyLDV.html
#
# Inputs:
#
# t: array of real values, random arrival times (single block
#    or whole LDV data set)
# u: array of real or complex values, velocities (single block
#    or whole LDV data set)
# w: array of real values, individual weights (single block
#    or whole LDV data set), e.g. transit times, used only,
#    if atw==False)
# atw: boolean, forward-backward inter-arrival time weighting
#      (if atw==True) or other weighting using w (if atw==False)
#      (atw==True derives the weights from the array t and
#      overrides w)
# Tc: real value, estimated correlation time (correlation is
#     estimated between -Tc/2 and Tc/2)
# Fs: real value, maximum sampling frequency (inverse temporal
#     resolution of the correlation)
# meanrm: boolean, if meanrm==True, mean is calculated and
#         removed from the block, this also enables Bessel's
#         correction
#
# Options:
#
# locnor: boolean, if locnor==True, local normalization is used
# If used without specifying locnor, local normalization is
# disabled
#
# Outputs:
#
# tau: array of real values, time lags of the correlation
#      (between -Tc/2 and Tc/2)
# R: array of complex values, correlation values (R(tau), with
#    dimension of u^2)
# f: array of real values, frequencies of the power density
#    spectrum (between -Fs/2 and Fs/2)
# S: array of complex values, power spectral density (S(f), with
#    dimension of u^2/Fs)
#
  tx=t-t[0]
  dt=1.0/float(Fs)
  K=int(round(Tc*Fs))
  JT=int(ceil(tx[len(tx)-1]*Fs))+1
  T=JT*dt
  J=JT+K//2
  if atw:
    wfw=append(tx[1:len(tx)]-tx[0:len(tx)-1],[0])
    for i in range(0,len(wfw)):
      if (wfw[i]>5*(tx[len(tx)-1]-tx[0])/(len(tx)-1.0)) or (wfw[i]<0):
        wfw[i]=0.0
    wbw=append([0],tx[1:len(tx)]-tx[0:len(tx)-1])
    for i in range(0,len(wbw)):
      if (wbw[i]>5*(tx[len(tx)-1]-tx[0])/(len(tx)-1.0)) or (wbw[i]<0):
        wbw[i]=0.0
  else:
    wbw=w
  if meanrm:
    ux=u-sum(array(wbw)*array(u))/float(sum(array(wbw)))
  else:
    ux=u
  if atw:
    df=float(Fs)/J
    fp=(((arange(0,J)+(J//2))%J)-(J//2))*df
    U=zeros(J)+0j
    Ubw=zeros(J)+0j
    Wbw=zeros(J)+0j
    W=zeros(J)+0j
    if locnor:
      Q1=zeros(J)+0j
      Q2=zeros(J)+0j
      Qbw=zeros(J)+0j
    for k in range(0,len(ux)):
      E=exp(-2j*pi*fp*round(tx[k]/dt)*dt)
      Ufw=wfw[k]*ux[k]*E
      U+=conj(Ubw)*Ufw
      Ubw+=wbw[k]*ux[k]*E
      Wfw=wfw[k]*E
      W+=conj(Wbw)*Wfw
      if locnor:
        Qfw=wfw[k]*abs(ux[k])**2*E
        Q1+=conj(Qbw)*Wfw
        Q2+=conj(Wbw)*Qfw
        Qbw+=wbw[k]*abs(ux[k])**2*E
      Wbw+=wbw[k]*E
    Su=T*(U+conj(U))/(W[0]+conj(W[0]))
    Ru=ifft(Su)*Fs
    if (not locnor) or meanrm:
      Sw=T*(W+conj(W))/(W[0]+conj(W[0]))
      Rw=ifft(Sw)*Fs
    if locnor:
      Sq1=T*(Q1+conj(Q2))/(W[0]+conj(W[0]))
      Rq1=ifft(Sq1)*Fs
      Sq2=T*(Q2+conj(Q1))/(W[0]+conj(W[0]))
      Rq2=ifft(Sq2)*Fs
  else:
    up=zeros(J)+0j
    wp=zeros(J)
    if locnor:
      qp=zeros(J)
    for k in range(0,len(ux)):
      j=int(round(tx[k]/dt))
      if j<J:
        up[j]+=wbw[k]*ux[k]
        wp[j]+=wbw[k]
        if locnor:
          qp[j]+=wbw[k]*abs(ux[k])**2
    Up=fft(up)
    Wp=fft(wp)
    Su=T*(conj(Up)*Up-sum(abs(wbw*ux)**2))/(conj(Wp[0])*Wp[0]-sum(abs(wbw)**2))
    Ru=ifft(Su)*Fs
    if (not locnor) or meanrm:
      Sw=T*(conj(Wp)*Wp-sum(abs(wbw)**2))/(conj(Wp[0])*Wp[0]-sum(abs(wbw)**2))
      Rw=ifft(Sw)*Fs
    if locnor:
      Qp=fft(qp)
      Sq1=T*(conj(Qp)*Wp-sum(abs(wbw*ux)**2))/(conj(Wp[0])*Wp[0]-sum(abs(wbw)**2))
      Sq2=T*(conj(Wp)*Qp-sum(abs(wbw*ux)**2))/(conj(Wp[0])*Wp[0]-sum(abs(wbw)**2))
      Rq1=ifft(Sq1)*Fs
      Rq2=ifft(Sq2)*Fs
  df=float(Fs)/K
  tau=(((arange(0,K)+(K//2))%K)-(K//2))*dt
  f=(((arange(0,K)+(K//2))%K)-(K//2))*df
  R=zeros(K)+0j
  if locnor:
    V=sum(array(wbw)*abs(array(ux))**2)/float(sum(array(wbw)))
  if meanrm:
    #Bessel's correction with correlation and weighting
    ru=zeros(K)+0j
    rw=zeros(K)+0j
    for k in range(-K//2,K-K//2):
      ru[k%K]=Ru[k%J]
      rw[k%K]=Rw[k%J]
    if atw:
      VW=W[0]+conj(W[0])
      VM=(sum(ru)/float(JT)+sum(array(wfw)*array(wbw)*abs(array(ux))**2)/VW)/(1-sum(real(rw))/float(JT))
    else:
      VW=conj(Wp[0])*Wp[0]-sum(abs(wbw)**2)
      VM=(sum(ru)/float(JT)+sum(abs(array(wbw)*array(ux))**2)/VW)/(1-sum(real(rw))/float(JT))
  for k in range(-K//2,K-K//2):
    if (tau[k%K]>-T) and (tau[k%K]<T):
      if locnor:
        if meanrm:
          #Bessel's correction with correlation and weighting
          if real((Rq1[k%J]+VM*Rw[k%J])*(Rq2[k%J]+VM*Rw[k%J]))>0:
            R[k%K]=(V+VM)*(Ru[k%J]+VM*Rw[k%J])/sqrt(real((Rq1[k%J]+VM*Rw[k%J])*(Rq2[k%J]+VM*Rw[k%J])))
          else:
            R[k%K]=0
        else:
          if real(Rq1[k%J]*Rq2[k%J])>0:
            R[k%K]=V*Ru[k%J]/sqrt(real(Rq1[k%J]*Rq2[k%J]))
          else:
            R[k%K]=0
      else:
        if real(Rw[k%J])>0:
          if meanrm:
            #Bessel's correction with correlation and weighting
            R[k%K]=Ru[k%J]/real(Rw[k%J])+VM
          else:
            R[k%K]=Ru[k%J]/real(Rw[k%J])
        else:
          R[k%K]=0
    else:
      R[k%K]=0
  S=dt*fft(R)
  return tau,R,f,S

from numpy import *
from numpy.fft import *

def fuzzy(t,u,w,atw,Tc,Fs,meanrm,locnor=False,fuzzy=False):
#
# for further information see
# http://ldvproc.nambis.de/programs/pyLDV.html
#
# Inputs:
#
# t: array of real values, random arrival times (single block
#    or whole LDV data set)
# u: array of real or complex values, velocities (single block
#    or whole LDV data set)
# w: array of real values, individual weights (single block
#    or whole LDV data set), e.g. transit times, used only,
#    if atw==False)
# atw: boolean, forward-backward inter-arrival time weighting
#      (if atw==True) or other weighting using w (if atw==False)
#      (atw==True derives the weights from the array t and
#      overrides w)
# Tc: real value, estimated correlation time (correlation is
#     estimated between -Tc/2 and Tc/2)
# Fs: real value, maximum sampling frequency (inverse temporal
#     resolution of the correlation)
# meanrm: boolean, if meanrm==True, mean is calculated and
#         removed from the block, this also enables Bessel's
#         correction
#
# Options:
#
# locnor: boolean, if locnor==True, local normalization is used
# fuzzy: boolean, if fuzzy==True, fuzzy quantization of arrival
#        times is used
# If used without specifying locnor and fuzzy, local normalization
# and fuzzy quantization of arrival times are disabled
#
# Outputs:
#
# tau: array of real values, time lags of the correlation
#      (between -Tc/2 and Tc/2)
# R: array of complex values, correlation values (R(tau), with
#    dimension of u^2)
# f: array of real values, frequencies of the power density
#    spectrum (between -Fs/2 and Fs/2)
# S: array of complex values, power spectral density (S(f), with
#    dimension of u^2/Fs)
#
  tx=t-t[0]
  dt=1.0/float(Fs)
  K=int(round(Tc*Fs))
  JT=int(ceil(tx[len(tx)-1]*Fs))
  T=JT*dt
  J=JT+K//2+1
  if locnor:
    J+=(J+1)%2 #for loc. nor. allways odd
  if atw:
    wfw=append(tx[1:len(tx)]-tx[0:len(tx)-1],[0])
    for i in range(0,len(wfw)):
      if (wfw[i]>5*(tx[len(tx)-1]-tx[0])/(len(tx)-1.0)) or (wfw[i]<0):
        wfw[i]=0.0
    wbw=append([0],tx[1:len(tx)]-tx[0:len(tx)-1])
    for i in range(0,len(wbw)):
      if (wbw[i]>5*(tx[len(tx)-1]-tx[0])/(len(tx)-1.0)) or (wbw[i]<0):
        wbw[i]=0.0
  else:
    wfw=w
    wbw=w
  if meanrm:
    ux=u-sum(array(wbw)*array(u))/float(sum(array(wbw)))
  else:
    ux=u
  df=float(Fs)/J
  fp=(((arange(0,J)+(J//2))%J)-(J//2))*df
  U=zeros(J)+0j
  Ubw=zeros(J)+0j
  Wbw=zeros(J)+0j
  W=zeros(J)+0j
  if locnor:
    Q1=zeros(J)+0j
    Q2=zeros(J)+0j
    Qbw=zeros(J)+0j
  for k in range(0,len(ux)):
    E=exp(-2j*pi*fp*tx[k])
    Ufw=wfw[k]*ux[k]*E
    U+=conj(Ubw)*Ufw
    Ubw+=wbw[k]*ux[k]*E
    Wfw=wfw[k]*E
    W+=conj(Wbw)*Wfw
    if locnor:
      Qfw=wfw[k]*abs(ux[k])**2*E
      Q1+=conj(Qbw)*Wfw
      Q2+=conj(Wbw)*Qfw
      Qbw+=wbw[k]*abs(ux[k])**2*E
    Wbw+=wbw[k]*E
  if fuzzy:
    U*=sinc(fp/Fs)**2
    W*=sinc(fp/Fs)**2
    if locnor:
      Q1*=sinc(fp/Fs)**2
      Q2*=sinc(fp/Fs)**2
  Su=T*(U+conj(U))/(W[0]+conj(W[0]))
  Ru=ifft(Su)*Fs
  if (not locnor) or meanrm:
    Sw=T*(W+conj(W))/(W[0]+conj(W[0]))
    Rw=ifft(Sw)*Fs
  if locnor:
    Sq1=T*(Q1+conj(Q2))/(W[0]+conj(W[0]))
    Rq1=ifft(Sq1)*Fs
    Sq2=T*(Q2+conj(Q1))/(W[0]+conj(W[0]))
    Rq2=ifft(Sq2)*Fs
  df=float(Fs)/K
  tau=(((arange(0,K)+(K//2))%K)-(K//2))*dt
  f=(((arange(0,K)+(K//2))%K)-(K//2))*df
  R=zeros(K)+0j
  if locnor:
    V=sum(array(wbw)*abs(array(ux))**2)/float(sum(array(wbw)))
  if meanrm:
    #Bessel's correction with correlation and weighting
    ru=zeros(K)+0j
    rw=zeros(K)+0j
    for k in range(-K//2,K-K//2):
      ru[k%K]=Ru[k%J]
      rw[k%K]=Rw[k%J]
    VW=W[0]+conj(W[0])
    VM=(sum(ru)/float(JT)+sum(array(wfw)*array(wbw)*abs(array(ux))**2)/VW)/(1-sum(real(rw))/float(JT))
  for k in range(-K//2,K-K//2):
    if (tau[k%K]>-T) and (tau[k%K]<T):
      if locnor:
        if meanrm:
          #Bessel's correction with correlation and weighting
          if real((Rq1[k%J]+VM*Rw[k%J])*(Rq2[k%J]+VM*Rw[k%J]))>0:
            R[k%K]=(V+VM)*(Ru[k%J]+VM*Rw[k%J])/sqrt(real((Rq1[k%J]+VM*Rw[k%J])*(Rq2[k%J]+VM*Rw[k%J])))
          else:
            R[k%K]=0
        else:
          if real(Rq1[k%J]*Rq2[k%J])>0:
            R[k%K]=V*Ru[k%J]/sqrt(real(Rq1[k%J]*Rq2[k%J]))
          else:
            R[k%K]=0
      else:
        if real(Rw[k%J])>0:
          if meanrm:
            #Bessel's correction with correlation and weighting
            R[k%K]=Ru[k%J]/real(Rw[k%J])+VM
          else:
            R[k%K]=Ru[k%J]/real(Rw[k%J])
        else:
          R[k%K]=0
    else:
      R[k%K]=0
  S=dt*fft(R)
  return tau,R,f,S
