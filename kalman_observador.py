import pandas as pd
import matplotlib.pyplot as plt
import mplcyberpunk
import numpy as np
import math
import random
from tqdm import tqdm

#plt.style.use("cyberpunk")

def generate_FEM_sqr(escala, Pi, Ps):
    FEM_DTC=np.zeros((2*Ps+1,3),dtype=float)
    for t in range(0,int(1*Ps/3)+1):
        FEM_DTC[t,0]=6/Ps*t-1

    for t in range(round((1*Ps/3)),int(3*Ps/3)+1):
        FEM_DTC[t,0]=1

    for t in range(round((3*Ps/3)),round((4*Ps/3))):
        FEM_DTC[t,0]=-(6/(Ps))*t+7

    for t in range(round((4*Ps/3)),round((6*Ps/3))):
        FEM_DTC[t,0]=-1

    
    desloca = 90 #graus de deslocamento
    FEM_DTC[:,1]=FEM_DTC[:,0]
    for t in range(0,2*Ps):
        a=t+round((1/(180/desloca)*Ps))
        if a>2*Ps-1:
            FEM_DTC[t,1]=FEM_DTC[a-2*(Ps),0]
        else:
            FEM_DTC[t,1]=FEM_DTC[a,0]
            
    FEM_DTC[:,0]=FEM_DTC[:,1]
    for t in range(0,round((2*Pi)/escala)):
        a=t-round((2*Pi/3)/escala)
        if a<0:
            FEM_DTC[t,1]=FEM_DTC[int(2*Pi/escala-round((2*Pi/3)/escala))+t,0]
        else:
            FEM_DTC[t,1]=FEM_DTC[a,0]

    for t in range(0,round((2*Pi)/escala)):
        a=t+round((2*Pi/3)/escala)
        if a>2*Pi/escala-1:
            FEM_DTC[t,2]=FEM_DTC[int(a-2*Pi/escala),0]
        else:
            FEM_DTC[t,2]=FEM_DTC[a,0]

    return FEM_DTC

sign = lambda x: math.copysign(1, x)

#Parametros da simulação
escala=0.0001
PI=3.1415
PS=int(PI/escala)
theta=np.arange(0,2*PI,escala)
t1=np.arange(0,2*PS,1)


FEM_sin=np.zeros((2*PS,3),dtype=float)
DFEM_sin=np.zeros((2*PS,3),dtype=float)
for i in range(0,2*PS):
    FEM_sin[i,1]=math.sin(theta[i] + math.pi/3.0)
    FEM_sin[i,2]=math.sin(theta[i] - (2.0/3.0)*math.pi + math.pi/3.0)
    FEM_sin[i,0]=math.sin(theta[i] + (2.0/3.0)*math.pi + math.pi/3.0)
    
    DFEM_sin[i,1]=math.cos(theta[i] + math.pi/3.0)
    DFEM_sin[i,2]=math.cos(theta[i] - (2.0/3.0)*math.pi + math.pi/3.0)
    DFEM_sin[i,0]=math.cos(theta[i] + (2.0/3.0)*math.pi + math.pi/3.0)
    
FEM_sqr = generate_FEM_sqr(escala,PI,PS)
DFEM_sqr=np.zeros((2*PS,3),dtype=float)
DFEM_sqr[:,0] = np.diff(FEM_sqr[:,0])
DFEM_sqr[:,1] = np.diff(FEM_sqr[:,1])
DFEM_sqr[:,2] = np.diff(FEM_sqr[:,2])

#Parametros.
R           = 5.75
L           = 0.55*pow(10,-3)
J           = 0.87*pow(10,-3)
BV          = 0.362*pow(10,-3)
Ke          = 0.31
p           = 4
Pp          = p/2
M           = 0
Tmax        = 200
Rmax        = 3500*0.1047
rvolt       = Tmax/Rmax
cmax        = 120
Tl          = 0.0667

E           = Tmax
Nbits_pwm   = 8
P           = pow(2,Nbits_pwm)
Vmax        = Tmax

EAD         = 3
Nbits_ad    = 12
Q_ad        = pow(2,Nbits_ad)
offsetCS    = 0.5
K_ad        = EAD/Rmax
Kkal        = EAD/cmax

kp          = 40
ki          = 3
kd          = 6

kp2         = 0.1
ki2         = 0.0012
kd2         = 0.001

kp3         = 0.1
ki3         = 0.0012
kd3         = 0.001

ui          = 0
ui2         = 0
ui3         = 0

ud          = 0
ud2         = 0
ud3         = 0

ui_kal          = 0
ui2_kal         = 0
ui3_kal         = 0

ud_kal          = 0
ud2_kal         = 0
ud3_kal         = 0

clarke      = np.array([[math.sqrt(1/2), math.sqrt(1/2),   math.sqrt(1/2)],
                        [1,              -1/2,             -1/2],
                        [0,              math.sqrt(3)/2,   -math.sqrt(3)/2]])

clarke     = np.multiply(clarke,math.sqrt(2/3))
inversor   = np.array([[2/3,       -1/3,        -1/3],
                       [-1/3,      2/3,         -1/3],
                       [-1/3,      -1/3,        2/3]])

inv_clarke = np.linalg.inv(clarke)

Ai         = np.array([[-R/(L-M),    0,           0,           -Pp*Ke/(L-M),         0],
                       [0,           -R/(L-M),    0,           -Pp*Ke/(L-M),         0],
                       [0,           0,           -R/(L-M),    -Pp*Ke/(L-M),         0],
                       [Pp*Ke/J,     Pp*Ke/J,     Pp*Ke/J,     -BV/J,                0],
                       [0,           0,           0,           Pp,                   0]])


B          = np.array([[1/(L-M),      0,           0,           0],
                       [0,            1/(L-M),     0,           0],
                       [0,            0,           1/(L-M),     0],
                       [0,            0,           0,           -1/J],
                       [0,            0,           0,           0]])

I5         = np.array([[1, 0, 0, 0, 0],
                       [0, 1, 0, 0, 0],
                       [0, 0, 1, 0, 0],
                       [0, 0, 0, 1, 0],
                       [0, 0, 0, 0, 1]])


I3         = np.array([[1, 0, 0],
                       [0, 1, 0],
                       [0, 0, 1]])



FS        = 10000
Ts        = 1/FS
Tc        = Ts/P

ciclos    = int(FS/10)
alin      = int(FS/4)
alto      = 100
alto2     = 150
N         = FS*4

r         = np.zeros(N)
ri        = np.zeros(N)
N2        = int(N/4-ciclos)
N3        = int((N/2))-ciclos
N4        = int(N*2/3)-ciclos

t         = np.zeros(N)
tc        = np.zeros(N*P+P)

x         = np.zeros((5,N*P+P))
yr        = np.zeros((3,N*P+P))
ym        = np.zeros((3,N))
ym_vel    = np.zeros(N)
backe     = np.zeros((3,N*P+P))
torque_kal = np.zeros(N*P+P)
torque_e_kal = np.zeros(N*P+P)

PWM       = np.zeros((3,N*P+P))
Vn        = np.zeros((4,N*P+P))
V         = np.zeros((3,N))

e         = np.zeros(N)
e2        = np.zeros(N)
e3        = np.zeros(N)

e_kal         = np.zeros(N)
e2_kal        = np.zeros(N)
e3_kal        = np.zeros(N)

iab       = np.zeros((3,N))


for i in range(alin,N):
    if i<ciclos+alin:
        ri[i]=(alto/ciclos)*(i-alin)
    elif i>=ciclos+alin and i < N2-ciclos:
        ri[i]=alto
    elif i>=N2-ciclos and i <N2:
        ri[i]=((alto2-alto)/ciclos)*(i-(N2-ciclos)) + alto
    elif i>=N2 and i < N3-ciclos:
        ri[i]=alto2
    elif i>=N3-ciclos and i< N3+3*ciclos:
        ri[i]=(-alto2/(4*ciclos))*i+ciclos-alto2-25
    elif i>= N3+4*ciclos and i < N4-ciclos:
        ri[i] = ((alto)/(N4-ciclos-(N3+4*ciclos)))*(i-(N3+4*ciclos))
    elif i>=N4-ciclos and i < N4+ciclos:
        ri[i] = alto
    elif i>=N4+ciclos and i < N4+3*ciclos:
        ri[i] = -(alto)/((N4+3*ciclos)-(N4+ciclos))*i + ciclos+433

for i in range(1,N):
    r[i]=0.98*r[i-1]+0.02*ri[i]


C_kal = np.array([[1, 0, 0, 0, 0],
                  [0, 1, 0, 0, 0],
                  [0, 0, 1, 0, 0]],dtype=np.float64)

x_pri = np.zeros((5,N),dtype=np.float64)
y_k = np.zeros((N*P+P,5),dtype=np.float64)
x_pos = np.zeros((5,N),dtype=np.float64)


Qn = I5*np.multiply([[0.1], [0.1], [0.1], [0.22], [0.22]],1) #rever estes valores para o meu sistema
Rn = I3*np.multiply([[0.1],[0.1],[0.1]],1) # Este também

Qn=Qn*0.022
Rn=Rn*2.2

ang = 73
pi=math.pi


theta0_kal_est = np.float64(ang*pi/180)
angulo_kal_est = np.zeros((N),np.float64)
angulo_kal_est[0] = (np.float64(theta0_kal_est))
thetae_kal_est = (round(((theta0_kal_est*(Pp))%(2*pi))/escala))

theta0_kal = ang*pi/180
angulo_kal = np.zeros((N*P+P),np.float64)
angulo_kal[0] = (np.float64(theta0_kal_est))
thetae_kal = round(((theta0_kal*(Pp))%(2*pi))/escala)


V_pos = 0
#P_pos = np.zeros((5,5),dtype=float)
P_pos = I5
kc=0

backe[:,0]=([0,0,0])


pwm = np.zeros((3,N*P+P),dtype=np.float64)
vg = []


for k in tqdm(range(0,alin),ascii=True, desc="Alinhamento: "):
    t[k]=k*Ts
    
    
    Ax     = np.array([[1,                             1,                            1,                          FEM_sin[thetae_kal,0],         1],
                       [1,                             1,                            1,                          FEM_sin[thetae_kal,1],         1],
                       [1,                             1,                            1,                          FEM_sin[thetae_kal,2],         1],
                       [FEM_sin[thetae_kal,0],         FEM_sin[thetae_kal,1],        FEM_sin[thetae_kal,2],      1,                             1],
                       [1,                             1,                            1,                          1,                             1]],dtype=np.float64)

    Ac = Ai*Ax
    Ad = I5+Ac*Tc+(Ac**2)/2*Tc**2
    Bd = ((I5*Tc)+(Ac/2)*(Tc**2))@B

    Ax_kal = np.array([[1,                                 1,                                1,                              FEM_sin[thetae_kal_est,0],         1],
                       [1,                                 1,                                1,                              FEM_sin[thetae_kal_est,1],         1],
                       [1,                                 1,                                1,                              FEM_sin[thetae_kal_est,2],         1],
                       [FEM_sin[thetae_kal_est,0],         FEM_sin[thetae_kal_est,1],        FEM_sin[thetae_kal_est,2],      1,                                 1],
                       [1,                                 1,                                1,                              1,                                 1]],dtype=np.float64)

    Ac_kal = Ai*Ax_kal
    Ad_kal = I5+Ac_kal*Ts+((Ac_kal**2)/2)*Ts**2
    Bd_kal = ((I5*Ts)+(Ac_kal/2)*(Ts**2))@B
    
    wr = x_pos[3,k]
    [ia, ib, ic] = ym[:,k]
    aux = (Pp*Ke*(DFEM_sin[thetae_kal_est,0]*ia + DFEM_sin[thetae_kal_est,1]*ib + DFEM_sin[thetae_kal_est,2]*ic))/J
    
    Fd = np.array([[-R/(L-M),                               0,                                         0,                                         -Pp*Ke*FEM_sin[thetae_kal_est,0]/(L-M),           -Pp*Ke*wr*DFEM_sin[thetae_kal_est,0]/(L-M)],
                   [0,                                      -R/(L-M),                                  0,                                         -Pp*Ke*FEM_sin[thetae_kal_est,1]/(L-M),           -Pp*Ke*wr*DFEM_sin[thetae_kal_est,1]/(L-M)],
                   [0,                                      0,                                         -R/(L-M),                                  -Pp*Ke*FEM_sin[thetae_kal_est,2]/(L-M),           -Pp*Ke*wr*DFEM_sin[thetae_kal_est,2]/(L-M)],
                   [Pp*Ke*FEM_sin[thetae_kal_est,0]/J,      Pp*Ke*FEM_sin[thetae_kal_est,1]/J,         Pp*Ke*FEM_sin[thetae_kal_est,2]/J,         -BV/J,                                            aux],
                   [0,                                      0,                                         0,                                         Pp,                                               0]],dtype=np.float64)
    
    
    #dFd = I5 + Fd*Ts
    dFd = I5 + Fd*Ts + ((Fd**2)/2)*(Ts**2)
    
    if V_pos < Vmax/8:
        V_pos = V_pos+0.01
    
    V[0:3,k] = inversor@[V_pos,V_pos,0]+np.multiply(backe[:,k],1/3)
    
    d = [abs(V_pos)*P/E,abs(V_pos)*P/E,0]

    #Etapa de superdiscretização
    for cc in range(0,P):
        kc = k*P+cc
        tc[kc] = kc*Tc
        
        
        for aux in range(0,3):
            if cc < d[aux]:
                PWM[aux,kc] = sign(V_pos)*E
            else:
                PWM[aux,kc] = 0
                
        torque_m   = 0.02*x[3,kc]+Tl
        Vn[0:3,kc] = (inversor@PWM[:,kc]) + np.multiply(backe[:,k],1/3)
        Vn[3,kc]   = torque_m
        
        
        x[:,kc+1] = Ad@x[:,kc]+Bd@Vn[:,kc]
        yr[:,kc+1] = C_kal@x[:,kc+1] #+ (Rn@[1,1,1])*random.random()
        
        theta_kal = theta0_kal + Pp*x[3,kc+1]*Tc
        theta0_kal = theta_kal%(Pp*2*pi)
        angulo_kal[kc] =((theta_kal/Pp)%(2*pi))
        thetae_kal = round((theta0_kal%(2*pi))/escala)
        
        if thetae_kal >= 62830:
            thetae_kal = 62829
            
            
        
        torque_e_kal[kc] = ((FEM_sqr[thetae_kal,:]*x[0:3,kc+1]).sum()*Ke)
        torque_kal[kc]   = (torque_e_kal[kc]-torque_m -BV*x[3,kc])
    
    ym[:,k+1] = [((((round(yr[0,kc+1])*Kkal+offsetCS)*Q_ad/EAD)/Q_ad)*EAD-offsetCS)/Kkal,((((round(yr[1,kc+1])*Kkal+offsetCS)*Q_ad/EAD)/Q_ad)*EAD-offsetCS)/Kkal,((((round(yr[2,kc+1])*Kkal+offsetCS)*Q_ad/EAD)/Q_ad)*EAD-offsetCS)/Kkal]
    backe[:,k+1]  = (np.multiply(FEM_sqr[thetae_kal,:],Ke*x[3,kc+1]))
    
    #x_pri[:,kc+1]  = x_pri[:,kc] + (Ad_kal@x_pri[:,kc] + Bd_kal@[Vn[0,kc],Vn[1,kc],Vn[2,kc],0.02*x_pos[3,kc+1]])*Tc
    x_pri[:,k+1]  = Ad_kal@x_pri[:,k]+Bd_kal@[V[0,k],V[1,k],0,0.02*x_pos[3,k]+Tl]
    P_pri         = dFd@P_pos@dFd.transpose() + Qn

    yk            = C_kal@x_pri[:,k+1]
    y_kal         = ym[:,k+1] - yk
    Sk            = C_kal@P_pri@C_kal.transpose() + Rn
    K_kal         = P_pri@C_kal.transpose()@np.linalg.inv(Sk)
    x_pos[:,k+1]  = x_pri[:,k+1] + (K_kal@y_kal)
    P_pos         = (I5 -(K_kal)@(C_kal))@(P_pri)
    
    theta_kal_est = theta0_kal_est + Pp*x_pos[3,k+1]*Ts
    theta0_kal_est = theta_kal_est%(Pp*2*pi)
    angulo_kal_est[k] = (theta_kal_est/Pp)%(2*pi)
    thetae_kal_est = round((theta0_kal_est%(2*pi))/escala)

    if thetae_kal_est >= 62830:
        thetae_kal_est = 62829
    
    
    ym_vel[k+1] = ((((round(x[3,kc+1])*K_ad+offsetCS)*Q_ad/EAD)/Q_ad)*EAD-offsetCS)/K_ad

## Simulação da etapa de aceleração e frenagem do motor.

Qn = I5*np.multiply([[0.1], [0.1], [0.1], [0.1], [0.1]],1) #rever estes valores para o meu sistema
Rn = I3*np.multiply([[0.1],[0.1],[0.1]],1) # Este também

Qn=Qn*0.099
Rn=Rn*56

for k in tqdm(range(alin,N-1),ascii=True,desc="Regime: "):
    t[k]=k*Ts
    
    
    Ax     = np.array([[1,                             1,                            1,                          FEM_sin[thetae_kal,0],         1],
                       [1,                             1,                            1,                          FEM_sin[thetae_kal,1],         1],
                       [1,                             1,                            1,                          FEM_sin[thetae_kal,2],         1],
                       [FEM_sin[thetae_kal,0],         FEM_sin[thetae_kal,1],        FEM_sin[thetae_kal,2],      1,                             1],
                       [1,                             1,                            1,                          1,                             1]],dtype=np.float64)

    Ac = Ai*Ax
    Ad = I5+Ac*Tc+((Ac**2)/2)*(Tc**2)
    Bd = ((I5*Tc)+(Ac/2)*(Tc**2))@B

    Ax_kal = np.array([[1,                                 1,                                1,                              FEM_sin[thetae_kal_est,0],         1],
                       [1,                                 1,                                1,                              FEM_sin[thetae_kal_est,1],         1],
                       [1,                                 1,                                1,                              FEM_sin[thetae_kal_est,2],         1],
                       [FEM_sin[thetae_kal_est,0],         FEM_sin[thetae_kal_est,1],        FEM_sin[thetae_kal_est,2],      1,                                 1],
                       [1,                                 1,                                1,                              1,                                 1]],dtype=np.float64)

    Ac_kal = Ai*Ax_kal
    Ad_kal = I5+Ac_kal*Ts+((Ac_kal**2)/2)*(Ts**2)
    Bd_kal = ((I5*Ts)+(Ac_kal/2)*(Ts**2))@B
    
    wr = x_pos[3,k]
    [ia, ib, ic] = ym[:,k]
    aux = (Pp*Ke*(DFEM_sin[thetae_kal_est,0]*ia + DFEM_sin[thetae_kal_est,1]*ib + DFEM_sin[thetae_kal_est,2]*ic))/J
    
    Fd = np.array([[-R/(L-M),                               0,                                         0,                                         -Pp*Ke*FEM_sin[thetae_kal_est,0]/(L-M),           -Pp*Ke*wr*DFEM_sin[thetae_kal_est,0]/(L-M)],
                   [0,                                      -R/(L-M),                                  0,                                         -Pp*Ke*FEM_sin[thetae_kal_est,1]/(L-M),           -Pp*Ke*wr*DFEM_sin[thetae_kal_est,1]/(L-M)],
                   [0,                                      0,                                         -R/(L-M),                                  -Pp*Ke*FEM_sin[thetae_kal_est,2]/(L-M),           -Pp*Ke*wr*DFEM_sin[thetae_kal_est,2]/(L-M)],
                   [Pp*Ke*FEM_sin[thetae_kal_est,0]/J,      Pp*Ke*FEM_sin[thetae_kal_est,1]/J,         Pp*Ke*FEM_sin[thetae_kal_est,2]/J,         -BV/J,                                            aux],
                   [0,                                      0,                                         0,                                         Pp,                                               0]],dtype=np.float64)
    
    
    dFd = I5 + Fd*Ts + ((Fd**2)/2)*(Ts**2)
    
    e[k] = r[k] - ym_vel[k]
    
    up = e[k]*kp
    ui = e[k]*ki + ui
    ud = 0.9*ud + 0.1*(kd*(e[k] - e[k-1]))
    
    if up > Vmax:
        up=Vmax
    if up < -Vmax:
        up=-Vmax
    
    if ui > Vmax:
        ui=Vmax
    if ui < -Vmax:
        ui=-Vmax
        
    if ud > Vmax:
        ud=Vmax
    if ud < -Vmax:
        ud=-Vmax
        
    pid = up+ui+ud
    if pid > Vmax:
        pid = Vmax
    if pid < -Vmax:
        pid = -Vmax
    
    alpha_beta = clarke@[ia, ib, ic]
    # Verificar se utiliza angulo mecanico ou eletrico nesta etapa.
    park = [[1, 0,                     0],
            [0, math.cos(theta0_kal),  math.sin(theta0_kal)],
            [0, -math.sin(theta0_kal), math.cos(theta0_kal)]]
    
    inv_park = np.array([[1,           0,                            0],
                [0,           math.cos(theta0_kal),        -math.sin(theta0_kal)],
                [0,           math.sin(theta0_kal),         math.cos(theta0_kal)]],dtype=np.float64)
    
    
    park_kal = [[1, 0,                     0],
            [0, math.cos(theta0_kal_est),  math.sin(theta0_kal_est)],
            [0, -math.sin(theta0_kal_est), math.cos(theta0_kal_est)]]
    
    inv_park_kal = np.array([[1,           0,                            0],
                [0,           math.cos(theta0_kal_est),        -math.sin(theta0_kal_est)],
                [0,           math.sin(theta0_kal_est),         math.cos(theta0_kal_est)]],dtype=np.float64)
    
    idq   = park@alpha_beta
    iq    = idq[2]
    
    e2[k] = pid - iq
    up2   = e2[k]*kp2
    ui2   = ki2*e2[k] + ui2
    ud2   = 0.9*ud2 + 0.1*kd2*(e2[k] - e2[k-1])
    
    pid2 = up2+ui2+ud2
    if pid2 >Vmax:
        pid2=Vmax
    if pid2 <-Vmax:
        pid2 = -Vmax
    
    i_d = idq[1]
    e3[k] = 0 -i_d
    up3   = e3[k]*kp3
    ui3   = ki3*e3[k] + ui3
    ud3   = 0.9*ud3 + 0.1*kd3*(e3[k] - e3[k-1])
    
    pid3 = up3+ui3+ud3
    if pid3 > Vmax:
        pid3=Vmax
    if pid3 <-Vmax:
        pid3 = -Vmax
    
    q = pid2
    d = pid3
    
    iab[:,k] = inv_park@[idq[0], d, q]
    
    V[:,k] = inv_clarke@iab[:,k]
    
    V_PWM = [abs(V[0,k])*P/E, abs(V[1,k])*P/E, abs(V[2,k])*P/E]
    
    
    
    ## Etapa de kalman
    
    e_kal[k] = r[k] - x_pos[3,k]
    
    up_kal = e_kal[k]*kp
    ui_kal = e_kal[k]*ki + ui_kal
    ud_kal = 0.9*ud_kal + 0.1*(kd*(e_kal[k] - e_kal[k-1]))
    
    if up_kal > Vmax:
        up_kal=Vmax
    if up_kal < -Vmax:
        up_kal=-Vmax
    
    if ui_kal > Vmax:
        ui_kal=Vmax
    if ui_kal < -Vmax:
        ui_kal=-Vmax
        
    if ud_kal > Vmax:
        ud_kal=Vmax
    if ud_kal < -Vmax:
        ud_kal=-Vmax
        
    pid = up_kal+ui_kal+ud_kal
    if pid > Vmax:
        pid = Vmax
    if pid < -Vmax:
        pid = -Vmax
    
    
    idq   = park_kal@alpha_beta
    iq    = idq[2]
    
    e2_kal[k] = pid - iq
    up2_kal   = e2_kal[k]*kp2
    ui2_kal   = ki2*e2_kal[k] + ui2_kal
    ud2_kal   = 0.9*ud2_kal + 0.1*kd2*(e2_kal[k] - e2_kal[k-1])
    
    pid2 = up2_kal+ui2_kal+ud2_kal
    if pid2 >Vmax:
        pid2=Vmax
    if pid2 <-Vmax:
        pid2 = -Vmax
    
    i_d = idq[1]
    e3_kal[k] = 0 -i_d
    up3_kal   = e3_kal[k]*kp3
    ui3_kal   = ki3*e3_kal[k] + ui3_kal
    ud3_kal   = 0.9*ud3_kal + 0.1*kd3*(e3_kal[k] - e3_kal[k-1])
    
    pid3 = up3_kal+ui3_kal+ud3_kal
    if pid3 > Vmax:
        pid3=Vmax
    if pid3 <-Vmax:
        pid3 = -Vmax
    
    q = pid2
    d = pid3
    
    V_kal  = inv_clarke@(inv_park_kal@[idq[0], d, q])
    V_kal  = (inversor@V_kal) + np.multiply(backe[:,k],1/3)
    
    #Etapa de superdiscretização
    for cc in range(0,P):
        kc = k*P+cc
        tc[kc] = kc*Tc
        
        
        for aux in range(0,3):
            if cc < V_PWM[aux]:
                PWM[aux,kc] = sign(V[aux,k])*E
            else:
                PWM[aux,kc] = 0
                
        torque_m   = 0.02*x[3,kc]+Tl
        Vn[0:3,kc] = (inversor@PWM[:,kc]) + np.multiply(backe[:,k],1/3)
        Vn[3,kc]   = torque_m
        
        
        x[:,kc+1] = Ad@x[:,kc]+Bd@Vn[:,kc]
        yr[:,kc+1] = C_kal@x[:,kc+1] #+ (Rn@[1,1,1])*random.random()
        
        theta_kal = theta0_kal + Pp*x[3,kc+1]*Tc
        theta0_kal = theta_kal%(Pp*2*pi)
        angulo_kal[kc] =((theta_kal/Pp)%(2*pi))
        thetae_kal = round((theta0_kal%(2*pi))/escala)
        
        if thetae_kal >= 62830:
            thetae_kal = 62829
            

        torque_e_kal[kc] = ((FEM_sqr[thetae_kal,:]*x[0:3,kc+1]).sum()*Ke)
        torque_kal[kc]   = (torque_e_kal[kc]-torque_m -BV*x[3,kc])
    
    ym[:,k+1] = [((((round(yr[0,kc+1])*Kkal+offsetCS)*Q_ad/EAD)/Q_ad)*EAD-offsetCS)/Kkal,((((round(yr[1,kc+1])*Kkal+offsetCS)*Q_ad/EAD)/Q_ad)*EAD-offsetCS)/Kkal,((((round(yr[2,kc+1])*Kkal+offsetCS)*Q_ad/EAD)/Q_ad)*EAD-offsetCS)/Kkal]
    backe[:,k+1]  = (np.multiply(FEM_sqr[thetae_kal,:],Ke*x[3,kc+1]))
    
    
    #x_pri[:,k+1]  = x_pri[:,k] + (Ad_kal@x_pri[:,k] + Bd_kal@[V[0,k],V[1,k],V[2,k],0.02*x_pos[3,k+1]+Tl])*Ts
    x_pri[:,k+1]  = Ad_kal@x_pri[:,k]+(Bd_kal@[V_kal[0],V_kal[1],V_kal[2],0.02*x_pos[3,k]+Tl])
    P_pri         = dFd@P_pos@dFd.transpose() + Qn

    yk            = C_kal@x_pri[:,k+1]
    y_kal         = ym[:,k+1] - yk
    Sk            = C_kal@P_pri@C_kal.transpose() + Rn
    K_kal         = P_pri@C_kal.transpose()@np.linalg.inv(Sk)
    x_pos[:,k+1]  = x_pri[:,k+1] + (K_kal@y_kal)
    P_pos         = (I5 -(K_kal)@(C_kal))@(P_pri)
    
    theta_kal_est = theta0_kal_est + Pp*x_pos[3,k+1]*Ts
    theta0_kal_est = theta_kal_est%(Pp*2*pi)
    angulo_kal_est[k] = (theta_kal_est/Pp)%(2*pi)
    thetae_kal_est = round((theta0_kal_est%(2*pi))/escala)

    if thetae_kal_est >= 62830:
        thetae_kal_est = 62829
    
    ym_vel[k+1] = ((((round(x[3,kc])*K_ad+offsetCS)*Q_ad/EAD)/Q_ad)*EAD-offsetCS)/K_ad

plt.plot(tc[0:N*P-P],angulo_kal[0:N*P-P])
plt.plot(t[0:N-1],angulo_kal_est[0:N-1],'--')
plt.grid(True)
#mplcyberpunk.add_glow_effects()
plt.legend(['Real','Estimado'])
plt.savefig('resultados/angulo_ob.pdf')
plt.show()
plt.clf()

plt.plot(t[0:N-1],ym_vel[0:N-1],'--')
plt.plot(t[0:N-1],x_pos[3,0:N-1],'--')
plt.plot(t[0:N-1],r[0:N-1],'--')
plt.grid(True)
plt.legend(['Real','Estimado','Ref'])
#mplcyberpunk.add_glow_effects()
plt.savefig('resultados/velocidade_ob.pdf')
plt.show()
plt.clf()

plt.plot(t[0:N-1],x_pos[3,0:N-1] - ym_vel[0:N-1])
plt.grid(True)
#mplcyberpunk.add_glow_effects()
plt.savefig('resultados/erro_vel_ob.pdf')
plt.show()
plt.clf()