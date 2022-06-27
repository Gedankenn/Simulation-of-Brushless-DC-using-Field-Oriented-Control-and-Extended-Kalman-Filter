# Simulation-of-Brushless-DC-using-Field-Oriented-Control-and-Extended-Kalman-Filter
Simulation of a brushless dc motor using FOC method and Extended Kalman Filter as estimation method for the velocity of the rotor.
 This code is a python simulation of a BLDC motor using Extended Kalman Filter for speed and position estimation.
This code as based on the works of:
* XIA, C. L. Permanent Magnet Brushless DC Motor Drives and Controls. [S. l.: s. n.], 2012.
ISBN 9781118188330. DOI: 10.1002/9781118188347.
* XIA, K. et al. Model predictive control method of torque ripple reduction for BLDC Motor.
IEEE Transactions on Magnetics, v. 56, n. 1, p. 1–6, 2020. ISSN 19410069. DOI:
10.1109/TMAG.2019.2950953.

* O’ROURKE, C. J. et al. A Geometric Interpretation of Reference Frames and Transformations:
Dq0, Clarke, and Park. IEEE Transactions on Energy Conversion, v. 34, n. 4, p. 2070–2083,2019. ISSN 15580059. DOI: 10.1109/TEC.2019.2941175.
* SCHMITZ, C. Projeto E Otimização De Motores Bldc De Imãs Permanentes Superficiais.
Dissertação de Mestrado, UFSC, 2017.
* TERZIC, B.; JADRIC, M. Design and implementation of the extended Kalman filter for the
speed and rotor position estimation of brushless dc motor. IEEE Transactions on Industrial
Electronics, IEEE, v. 48, n. 6, p. 1065–1073, 2001. ISSN 02780046. DOI:
10.1109/41.969385.

* LIU, Y.; ZHU, Z. Q.; HOWE, D. Instantaneous torque estimation in sensorless
direct-torque-controlled brushless DC motors. IEEE Transactions on Industry Applications,
v. 42, n. 5, p. 1275–1283, 2006. ISSN 00939994. DOI: 10.1109/TIA.2006.880854

The mathematical model implemented as its follows:<br>
$
\frac{d}{dt}
\left(\begin{array}{c}
I_a\\ 
I_b \\
I_c \\
\omega_r \\
\theta_e
\end{array}\right)
$=
$
\left(\begin{array}{ccccc} 
-\frac{R}{L-M} & 0 & 0 & -\frac{P_p\phi_mF_a(\theta_e)}{L-M} & 0\\
0 & -\frac{R}{L-M} & 0 & -\frac{P_p\phi_mF_b(\theta_e)}{L-M} & 0\\
0 & 0 & -\frac{R}{L-M} & -\frac{P_p\phi_mF_c(\theta_e)}{L-M} & 0\\
0 & 0 & -\frac{R}{L-M} & -\frac{P_p\phi_mF_c(\theta_e)}{L-M} & 0\\
\frac{P_p\phi_mF_A(\theta_e)}{J} & \frac{P_p\phi_mF_B(\theta_e)}{J} & \frac{P_p\phi_mF_C(\theta_e)}{J} & -\frac{B}{J} & 0 \\
0 & 0 & 0 & P_p & 0
\end{array}\right)
\left(\begin{array}{c} 
I_a \\
I_b \\
I_c \\
\omega_r \\
\theta_e
\end{array}\right)
+
\left(\begin{array}{cccc} 
\frac{1}{L-M} & 0 & 0 & 0\\
0 & \frac{1}{L-M} & 0 & 0\\
0 & 0 & \frac{1}{L-M} & 0\\
0 & 0 & 0 & -\frac{1}{J} \\
0 & 0 & 0 & 0
\end{array}\right)
\left(\begin{array}{c} 
V_an \\
V_bn \\
V_cn \\
T_c
\end{array}\right)
$
<br>$
\left(\begin{array}{c} 
I_a \\
I_b \\
I_c
\end{array}\right)
\left(\begin{array}{ccccc} 
1 & 0 & 0 & 0 & 0\\
0 & 1 & 0 & 0 & 0\\
0 & 0 & 1 & 0 & 0
\end{array}\right)
$
$
\left(\begin{array}{c} 
I_a \\
I_b \\
I_c \\
\omega_r \\
\theta_e
\end{array}\right)
$

$
\dot{x} = f(x,t) + v(t)
$

$
dF(x,t) = \frac{\partial f(x,t)}{\partial x}
$

$
dF(x,t) =
\left(\begin{array}{ccccc} 
-\frac{R}{L-M} & 0 & 0 & -\frac{P_p\phi_mF_a(\theta_e)}{L-M} & -\frac{P_p\phi_m\omega_r}{L-M}\frac{dF_a(theta_e)}{d\theta_e}\\
0 & -\frac{R}{L-M} & 0 & -\frac{P_p\phi_mF_b(\theta_e)}{L-M} & -\frac{P_p\phi_m\omega_r}{L-M}\frac{dF_b(theta_e)}{d\theta_e}\\
0 & 0 & -\frac{R}{L-M} & -\frac{P_p\phi_mF_c(\theta_e)}{L-M} & -\frac{P_p\phi_m\omega_r}{L-M}\frac{dF_c(theta_e)}{d\theta_e}\\
\frac{P_p\phi_mF_a(\theta_e)}{J} & \frac{P_p\phi_mF_b(\theta_e)}{J} & \frac{P_p\phi_mF_c(\theta_e}{J} & -\frac{B}{J} &
\frac{P_p\phi_m}{J}(\frac{dF_a(\theta_e)}{d\theta_e}i_a + \frac{dF_b(\theta_e)}{d\theta_e}i_b + \frac{dF_c(\theta_e)}{d\theta_e}i_c) \\
0 & 0 & 0 & P_p & 0
\end{array}\right)
$
<br>
## The discrete equation for F is:
$F_d = \phi = I + F.Ts$<br>
$dF_d = \phi = I + dF.Ts$<br>
$\hat{x}_{k|k-1} = \hat{x}_{k-1|k-1} + f(\hat{x}_{k-1|k-1},u_{k-1})Tc$

## For The Extended Kalman Filter the following equations are used:
### Estimation
<br>
$\hat{x}_{k|k+1} = F_d\hat{x}_{k|k} + B_du_{k}\\
P_{k|k+1} = dF_dP_{k|k}dF_d^T_s + Q $
<br><br>

### Prediction
<br>
$ \hat{y}_{k|k-1} = C\hat{x}_{k|k-1}\\
\tilde{y}_k = y_k - \tilde{y}_{k|k-1}\\
S_k = CP_{k|k-1}C^T + R\\
K_k = P_{k|k-1}C^TS_k^-1\\
\hat{x}_{k|k} = x_{k|k-1} + K_k\tilde{y}_k\\
P_{k|k} = (I -K_kC)P_{k|k-1} $
<br><br>

### Matrizes de covariancia
<br>

$R = diag[\sigma^2 \sigma^2 \sigma^2]$<br>
$Q = diag[q_{11}, q_{22}, q_{33}, q_{44}, q_{55}]$<br>
<br>
onde:<br>
<ul>
    <li> 

$\sigma$- 
    <li> 

$q_{11}, q_{22}, q_{33}$ - 
    <li>

$q_{44}, q_{55}$ - 
