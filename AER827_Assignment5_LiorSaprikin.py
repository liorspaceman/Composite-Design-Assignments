# -*- coding: utf-8 -*-
"""
Created on Sun Jan 24 13:32:37 2021

@author: liors
"""
import numpy as np
import math
from fractions import Fraction as frac
#********************************************************************************************************************************************************************************
# 1.1 Choose material to set property values for Modulus of elasticity (Ex,Ey,Es), Poisson's Ratio (NUx), and layer thickness (Ho)
material = int (input("Choose a material (input #): \n 1 = T300/N5208 \n 2 = E-Glass/Epoxy \n 3 = Kevlar/Epoxy \n 4 = AS/H3501 \n 5 = AS4/PEEK"))
if material == 1:
    mat_name = "T300/N5208"
    Ex=181.0;Ey=10.3;Es=7.17;NUx=0.28;Xt=1500.0;Xc=1500.0;Yt=40.0;Yc=246.0;Sc=68.0;Ho=0.125;   
elif material == 2:
    mat_name = "E-Glass/Epoxy"
    Ex=38.6;Ey=8.27;Es=4.14;NUx=0.26;Xt=1062.0;Xc=610.0;Yt=31.0;Yc=118.0;Sc=72.0;Ho=0.125;   
elif material == 3:
    mat_name = "Kevlar/Epoxy"
    Ex=76;Ey=5.5;Es=2.3;NUx=0.34;Xt=1400.0;Xc=235.0;Yt=12.0;Yc=53.0;Sc=34.0;Ho=0.125;  
elif material == 4:
    mat_name = "AS/H3501"
    Ex=138.0;Ey=8.96;Es=7.1;NUx=0.3;Xt=1447.0;Xc=1447.0;Yt=51.7;Yc=206.0;Sc=93.0;Ho=0.125;  
elif material == 5:
    mat_name = "AS4/PEEK"
    Ex=134.0;Ey=8.9;Es=5.1;NUx=0.28;Xt=2130.0;Xc=1100.0;Yt=80.0;Yc=200.0;Sc=160.0;Ho=0.125;  

# 1.2 Calculate On-Axis Compliance [Sxys] and Stiffness [Qxys] 
NUy=NUx*(Ey/Ex); Sxx = 1/Ex; Sxy=-NUy/Ey; Syx=-NUx/Ex; Syy=1/Ey; Sss=1/Es;
Sxys = np.array([[Sxx,Sxy,0],[Syx,Syy,0],[0,0,Sss]])
S_TPa = np.array ([[Sxx*1000,Sxy*1000,0],[Syx*1000,Syy*1000,0],[0,0,Sss*1000]])
Qxys = np.linalg.inv(Sxys)
#Qxx = Qxys[0,0];
#Qxy = Qxys[0,1];
#Qyy = Qxys[1,1];
#Qss = Qxys[2,2];
Qxx=134.7;
Qyy=8.96;
Qss=5.1;
Qxy=2.51;
# 1.3 Find U from Q
U1Q = 1/8*(3*Qxx + 3*Qyy + 2*Qxy + 4*Qss);
U2Q = 1/2*(Qxx - Qyy);
U3Q = 1/8*(Qxx + Qyy - 2*Qxy - 4*Qss);
U4Q = 1/8*(Qxx + Qyy + 6*Qxy - 4*Qss);
U5Q = 1/8*(Qxx + Qyy - 2*Qxy + 4*Qss);
    
# 1.4 Find U from S
U1S = 1/8*(3*Sxx + 3*Syy + 2*Sxy + Sss);
U2S = 1/2*(Sxx - Syy);
U3S = 1/8*(Sxx + Syy - 2*Sxy - Sss);
U4S = 1/8*(Sxx + Syy + 6*Sxy - Sss);
U5S = 1/8*(Sxx + Syy - 2*Sxy + Sss);

#**********************************************************************************************************************************************************************************
# 2.1 Build laminate code
laminatecode = []
# Input number of layers ex: 4*3*2(symmetrical)=24 layers
ntotal = int(input("Enter total # of layers in the laminate : \nie:[45/0/90/-45]3S=4*3*2=24 "))
# Option to add core with thickness (2Zc)
TwoZc = float(input("What is core thickness (mm)? \n(0 mm for no core)"))
Zc=TwoZc/2;
# Is it symmetrical ex: yes
sym = input("Is the laminate symmetrical? (yes/no)")
if sym == "yes":
    sym_int = 2.0;
    sym_code = "S";
elif sym == "no":
    sym_int = 1.0;
    sym_code = " ";
else:
    print("Invalid Answer. Rerun code and Retry")
# How many pattern repeats ex: 3
repeat_int = float(input("How many repeats?"))
# so laminatecode/(2(if symmetrical)*n())
n = int(ntotal/(sym_int*repeat_int));
total_thickness = (ntotal*Ho)+TwoZc;
Zc_star = TwoZc/total_thickness;
# Based on above n will be the simplified orientation code ex:[45/0/90/-45]
# iterating till the range
for i in range(0, ntotal):
   print("Enter layer #{}".format(i+1) + "'s orientation ie. 0,+45,-45,90")
   elm = int(input())
   laminatecode.append(elm) # adding the element
plythickness_matrix = np.full((ntotal), Ho)
laminate_short = []
for i in range (0,n):
    laminate_short.append(laminatecode[i])
    

#**********************************************************************************************************************************************************************************
# 3.1 Input Applied Stress Vectors
Force = 'Stress';
Axis = 'In-Plane';
unit = 'N/m or Pa*m';
print("\nPlease input the applied Bending Moments in Nm")
P = float(input("Total Load (N)?"))
b = float(input("Beam Width (m) ?"))
L = float(input("Beam Length (m) ?"))
#print("\nPlease input the applied stress vectors in N/m")
N1 = 0.5*P/b
N2 = 0
N6 = 0
#N1 = float(input("N1?"))
#N2 = float(input("N2?"))
#N6 = float(input("N6?"))
Force1 = N1; Force2 = N2; Force3 = N6; 
N126 = np.array([[N1],[N2],[N6]])

# 3.2 Find Angle Fraction 
numply0 = laminatecode.count(0)
v0=numply0/ntotal;
numply90 = laminatecode.count(90)
v90=numply90/ntotal;
numply45 = laminatecode.count(45)
v45=numply45/ntotal;
numplyn45 = laminatecode.count(-45)
vn45=numplyn45/ntotal;
numply30 = laminatecode.count(30)
v30=numply30/ntotal;
numplyn30 = laminatecode.count(-30)
vn30=numplyn30/ntotal;
print("v0:",frac(numply0,ntotal)," = ",v0,"\nv90:",frac(numply90,ntotal)," = ",v90,"\nv45:",frac(numply45,ntotal)," = ",v45,"\nv-45:",frac(numplyn45,ntotal)," = ",vn45,"\nv30:",frac(numply30,ntotal)," = ",v30,"\nv-30:",frac(numplyn30,ntotal)," = ",vn30)


# 3.3 Calc Vstars
V1star = v0*int(math.cos(2*math.radians(0)))+v90*int(math.cos(2*math.radians(90)))+v45*int(math.cos(2*math.radians(45)))+vn45*int(math.cos(2*math.radians(-45)));
print(V1star)
V2star =v0*int(math.cos(4*math.radians(0)))+v90*int(math.cos(4*math.radians(90)))+v45*int(math.cos(4*math.radians(45)))+vn45*int(math.cos(4*math.radians(-45)));
print(V2star)
V3star = v0*int(math.sin(2*math.radians(0)))+v90*int(math.sin(2*math.radians(90)))+v45*int(math.sin(2*math.radians(45)))+vn45*int(math.sin(2*math.radians(-45)));
print(V3star)
V4star = v0*int(math.sin(4*math.radians(0)))+v90*int(math.sin(4*math.radians(90)))+v45*int(math.sin(4*math.radians(45)))+vn45*int(math.sin(4*math.radians(-45)));
print(V4star)

# 3.4 Calculating Table 4.4, [A], [a]
t44_1UU = np.array([[1],[U2Q],[U3Q]]);
Table44 = np.array([[U1Q, V1star, V2star],[U1Q, -1*V1star, V2star],[U4Q, 0, -1*V2star],[U5Q, 0, -1*V2star],[0, 0.5*V3star, V4star],[0, 0.5*V3star, -1*V4star]])
Ah_value = [0,0,0,0,0,0]
for i in range(len(Table44)):
    for j in range(len(t44_1UU[0])):
        for k in range(len(t44_1UU)):
            Ah_value[i] += Table44[i][k] * t44_1UU[k][j]
        for r in Ah_value:
            print(r)
Ah_matrix =np.array([[Ah_value[0], Ah_value[2], Ah_value[4]],[Ah_value[2], Ah_value[1], Ah_value[5]],[Ah_value[4], Ah_value[5], Ah_value[3]]])
A_matrix = np.dot(total_thickness/7.66708,Ah_matrix)
a_matrix = np.dot(1000,np.linalg.inv(A_matrix))

# 3.5 Calculate off-axis Strain
EPSnot126big =  [0,0,0]
for i in range(len(a_matrix)):
    for j in range(len(N126[0])):
        for k in range(len(N126)):
            EPSnot126big[i] += a_matrix[i][k] * N126[k][j]
        for r in EPSnot126big:
            print(r)
        EPSnot126 = np.dot(1/100000,EPSnot126big)    
        EPSnot1 = EPSnot126[0];
        EPSnot2 = EPSnot126[1];
        EPSnot6 = EPSnot126[2];

#***********************************************************************************************************************************************
# 4 Input Bending Moments
# 4.1 Input Applied Stress Vectors
Force = 'Bending Moment';
Axis = 'In-Plane';
unit = 'N/m or Pa*m';

M1 = -1*P*L/(4*b);
M2 = 0;
M6 = 0;
M126 = np.array([[M1],[M2],[M6]])

# 4.2 V values from Eq. 5.43
Zc_star = float(TwoZc/total_thickness);
h_star = (total_thickness**3)/12*(1-Zc_star**3);
V1_term_set = {}
V2_term_set = {}
V3_term_set = {}
V4_term_set = {}
for i in range (6,ntotal):
    if i <= n-1:
        V1_term = 2/3*int(math.cos(2*math.radians(laminatecode[i])))*((((n-i)*Ho)+Zc)**3-(((n-i-1)*Ho)+Zc)**3)
        V2_term = 2/3*int(math.cos(4*math.radians(laminatecode[i])))*((((n-i)*Ho)+Zc)**3-(((n-i-1)*Ho)+Zc)**3)
        V3_term = 2/3*int(math.sin(2*math.radians(laminatecode[i])))*((((n-i)*Ho)+Zc)**3-(((n-i-1)*Ho)+Zc)**3)
        V4_term = 2/3*int(math.sin(4*math.radians(laminatecode[i])))*((((n-i)*Ho)+Zc)**3-(((n-i-1)*Ho)+Zc)**3)
    else:
        V1_term = 2/3*int(math.cos(2*math.radians(laminatecode[i])))*((((i-n+1)*Ho)+Zc)**3-(((i-n)*Ho)+Zc)**3)
        V2_term = 2/3*int(math.cos(4*math.radians(laminatecode[i])))*((((i-n+1)*Ho)+Zc)**3-(((i-n)*Ho)+Zc)**3)
        V3_term = 2/3*int(math.sin(2*math.radians(laminatecode[i])))*((((i-n+1)*Ho)+Zc)**3-(((i-n)*Ho)+Zc)**3)
        V4_term = 2/3*int(math.sin(4*math.radians(laminatecode[i])))*((((i-n+1)*Ho)+Zc)**3-(((i-n)*Ho)+Zc)**3)
    V1_term_set[i] = V1_term;
    V2_term_set[i] = V2_term;
    V3_term_set[i] = V3_term;
    V4_term_set[i] = V4_term;
print(V1_term_set)
V1_values = V1_term_set.values()
V1 = sum(V1_values)
print(V1)
V2_values = V2_term_set.values()
V2 = sum(V2_values)
print(V2)
V3_values = V3_term_set.values()
V3 = sum(V3_values)
print(V3)
V4_values = V4_term_set.values()
V4 = sum(V4_values)
print(V4)

# 4.3 Calculate Table 5.3 to find [D]
t53_hUU = [[h_star],[U2Q],[U3Q]];
Table53 = [[U1Q, V1, V2],[U1Q, -1*V1, V2],[U4Q, 0, -1*V2],[U5Q, 0, -1*V2],[0, 0.5*V3, V4],[0, 0.5*V3, -1*V4]]
D_value = [0,0,0,0,0,0]
for i in range(len(Table53)):
    for j in range(len(t53_hUU[0])):
        for k in range(len(t53_hUU)):
            D_value[i] += Table53[i][k] * t53_hUU[k][j]
        for r in D_value:
            print(r)
D_matrix =np.array([[D_value[0], D_value[2], D_value[4]],[D_value[2], D_value[1], D_value[5]],[D_value[4], D_value[5], D_value[3]]])
d_matrix = np.dot(1000,np.linalg.inv(D_matrix))

# 4.4 Use Table5.2 to calculate [k126] curvature values
k126 =  [0,0,0]
for i in range(len(d_matrix)):
    for j in range(len(M126[0])):
        for k in range(len(M126)):
            k126[i] += d_matrix[i][k] * M126[k][j]
        for r in k126:
            print(r)
        #126 = np.dot(1/100000,EPSnot126big)    
        k1 = k126[0];
        k2 = k126[1];
        k6 = k126[2];

# 4.5 In-plane strain EPS126 due to curvature [k126] for each ply 
EPS126_top_set = {}
EPS126_bottom_set = {}
EPSxys_set = {}
EPSxys_top_set = {}
EPSxys_bottom_set = {}
SIGxys_top_set = {}
SIGxys_bottom_set = {}
for i in range (0,ntotal):
    if i <= n-1:
        EPS126_top = np.dot(1/1000,np.array([EPSnot1+((n-i)*Ho+Zc)*k1, EPSnot2+((n-i)*Ho+Zc)*k2, EPSnot6+((n-i)*Ho+Zc)*k6]))
        EPS126_bottom = np.dot(1/1000,np.array([EPSnot1+((n-i-1)*Ho+Zc)*k1, EPSnot2+((n-i-1)*Ho+Zc)*k2, EPSnot6+((n-i-1)*Ho+Zc)*k6]))
    else:
        EPS126_top = np.dot(-1/1000, np.array([EPSnot1+((i-n+1)*Ho+Zc)*k1, EPSnot2+((i-n+1)*Ho+Zc)*k2, EPSnot6+((i-n+1)*Ho+Zc)*k6]))
        EPS126_bottom = np.dot(-1/1000, np.array([EPSnot1+((i-n)*Ho+Zc)*k1, EPSnot2+((i-n)*Ho+Zc)*k2, EPSnot6+((i-n)*Ho+Zc)*k6]))
    EPS126_top_set[i] = np.dot(1/1000,EPS126_top);
    EPS126_bottom_set[i] = np.dot(1/1000,EPS126_bottom);

# 4.6 Calculating On-Axis strain [EPSxys] for each layer using off-axis strain [EPS126] and Table 3.5

    for z in range(0,ntotal):
        theta = laminatecode[z];
        # Top strains and stresses
        p_top = 0.5*(EPS126_top[0]+EPS126_top[1]);
        q_top = 0.5*(EPS126_top[0]-EPS126_top[1]);
        r1_top = 0.5*EPS126_top[2];
        pqr_top = [[p_top],[q_top],[r1_top]]
        Table35_top = np.array([[1, int(math.cos(2*math.radians(theta))), int(math.sin(2*math.radians(theta)))],[1, int(-1*math.cos(2*math.radians(theta))), int(-1*math.sin(2*math.radians(theta)))],[0, int(-2*math.sin(2*math.radians(theta))), 2*int(math.cos(2*math.radians(theta)))]])
        EPSxys_top = np.dot(Table35_top, pqr_top)
        EPSxys_top_set[z] = np.dot(1/1000,EPSxys_top);
        if z < 6:
            SIGxys_top =  np.dot(Qxys,EPSxys_top)
        else:
            SIGxys_top =  np.dot(Qxys,-1*EPSxys_top)
        SIGxys_top_set[z] = SIGxys_top;
        
        # Bottom strains ans stresses
        p_bottom = 0.5*(EPS126_bottom[0]+EPS126_bottom[1]);
        q_bottom = 0.5*(EPS126_bottom[0]-EPS126_bottom[1]);
        r1_bottom = 0.5*EPS126_bottom[2];
        pqr_bottom = [[p_bottom],[q_bottom],[r1_bottom]]
        Table35_bottom = np.array([[1, int(math.cos(2*math.radians(theta))), int(math.sin(2*math.radians(theta)))],[1, int(-1*math.cos(2*math.radians(theta))), int(-1*math.sin(2*math.radians(theta)))],[0, int(-2*math.sin(2*math.radians(theta))), 2*int(math.cos(2*math.radians(theta)))]])
        EPSxys_bottom = np.dot(Table35_bottom,pqr_bottom)
        EPSxys_bottom_set[z] = np.dot(1/1000,EPSxys_bottom);
        if z < 6:
            SIGxys_bottom =  np.dot(Qxys,EPSxys_bottom)
        else:
            SIGxys_bottom =  np.dot(Qxys,-1*EPSxys_bottom)
        SIGxys_bottom_set[z] = SIGxys_bottom;
        
# 4.7 Additional Calculations 
d11 = d_matrix[0,0];
EPSx_top_max = EPSxys_top_set[11][0]
delta = (d11*P*L**3)/(4.8*b);

#***********************************************************************************************************************************************
# 5 Max Stress Criteria
R_top_set = {}
R_bottom_set = {}
R_min_set = {}
for i in range (0,ntotal):
    SIGx_top = (SIGxys_top_set[i])[0]
    SIGx_bottom = (SIGxys_bottom_set[i])[0]
    SIGy_top = (SIGxys_top_set[i])[1]
    SIGy_bottom = (SIGxys_bottom_set[i])[1]
    SIGs_top = (SIGxys_top_set[i])[2]
    SIGs_bottom = (SIGxys_bottom_set[i])[2]
    # Solving Rs in X direction
    if SIGx_top > 0:
        Rx_top = Xt/SIGx_top
    else:
        Rx_top = -1*Xc/SIGx_top
    if SIGx_bottom > 0:
        Rx_bottom = Xt/SIGx_bottom
    else:
        Rx_bottom = -1*Xc/SIGx_bottom
    # Solving Rs in Y direction
    if SIGy_top > 0:
        Ry_top = Yt/SIGy_top
    else:
        Ry_top = -1*Yc/SIGy_top
    if SIGy_bottom > 0:
        Ry_bottom = Yt/SIGy_bottom
    else:
        Ry_bottom = -1*Yc/SIGy_bottom
    # Solving Rs in S direction
    Rs_top = Sc/abs(SIGs_top)
    Rs_bottom = Sc/abs(SIGs_bottom)
    R_top_set[i] = [Rx_top,Ry_top,Rs_top]
    R_bottom_set[i] = [Rx_bottom,Ry_bottom,Rs_bottom]
    R_min_set[i] = min(Rx_top,Ry_top,Rs_top,Rx_bottom,Ry_bottom,Rs_bottom)

Fail_ply = 11
R_min_maxstress = R_min_set[Fail_ply]
maxstress_N = R_min_maxstress*N126
maxstress_M = R_min_maxstress*M126

#***********************************************************************************************************************************************
# 6 Tsai-Wu Failure Criteria 
# 6.1 Calc Tsai-Wu parameters
Fxx = 1/(Xt*Xc)
Fx = 1/Xt - 1/Xc
Fyy = 1/(Yt*Yc)
Fy = 1/Yt - 1/Yc
Fss = 1/(Sc**2)
Fxystar = -1/2
Fxy = Fxystar*math.sqrt(Fxx*Fyy)

# 6.1 Solve Rs using quadratic equation, Eq7.8, new stresses for each ply
R_top_tsai_set = {}
R_bottom_tsai_set = {}
R_min_tsai_set = {}
for i in range (0,ntotal):
    X_top = (SIGxys_top_set[i])[0]
    X_bottom = (SIGxys_bottom_set[i])[0]
    Y_top = (SIGxys_top_set[i])[1]
    Y_bottom = (SIGxys_bottom_set[i])[1]
    S_top = (SIGxys_top_set[i])[2]
    S_bottom = (SIGxys_bottom_set[i])[2]
    R_top_pos = (-1*(Fx*X_top+Fy*Y_top) + math.sqrt(((Fx*X_top+Fy*Y_top)**2)-(4*(Fxx*X_top**2 + 2*Fxy*X_top*Y_top + Fyy*Y_top**2 + Fss*S_top**2)*(-1))))/(2*(Fxx*X_top**2 + 2*Fxy*X_top*Y_top + Fyy*Y_top**2 + Fss*S_top**2))
    R_bottom_pos = (-1*(Fx*X_bottom+Fy*Y_bottom) + math.sqrt(((Fx*X_bottom+Fy*Y_bottom)**2)-(4*(Fxx*X_bottom**2 + 2*Fxy*X_bottom*Y_bottom + Fyy*Y_bottom**2 + Fss*S_bottom**2)*(-1))))/(2*(Fxx*X_bottom**2 + 2*Fxy*X_bottom*Y_bottom + Fyy*Y_bottom**2 + Fss*S_bottom**2))
    R_top_neg = (-1*(Fx*X_top+Fy*Y_top) - math.sqrt(((Fx*X_top+Fy*Y_top)**2)-(4*(Fxx*X_top**2 + 2*Fxy*X_top*Y_top + Fyy*Y_top**2 + Fss*S_top**2)*(-1))))/(2*(Fxx*X_top**2 + 2*Fxy*X_top*Y_top + Fyy*Y_top**2 + Fss*S_top**2))
    R_bottom_neg = (-1*(Fx*X_bottom+Fy*Y_bottom) - math.sqrt(((Fx*X_bottom+Fy*Y_bottom)**2)-(4*(Fxx*X_bottom**2 + 2*Fxy*X_bottom*Y_bottom + Fyy*Y_bottom**2 + Fss*S_bottom**2)*(-1))))/(2*(Fxx*X_bottom**2 + 2*Fxy*X_bottom*Y_bottom + Fyy*Y_bottom**2 + Fss*S_bottom**2))
    R_top_tsai_set[i] = [R_top_pos,R_top_neg]
    R_bottom_tsai_set[i] = [R_bottom_pos,R_bottom_neg]
    R_min_tsai_set[i] = min(R_top_pos,abs(R_top_neg),R_bottom_pos,abs(R_bottom_neg))
R_min_tsai = R_min_tsai_set[Fail_ply]
tsai_N = R_min_tsai*N126
tsai_M = R_min_tsai*M126

#***********************************************************************************************************************************************
# 7 2D Hashin Failure Criteria
R_top_hashin_set = {}
R_bottom_hashin_set = {}
R_min_hashin_set = {}
for i in range (0,ntotal):
    X_top = (SIGxys_top_set[i])[0]
    X_bottom = (SIGxys_bottom_set[i])[0]
    Y_top = (SIGxys_top_set[i])[1]
    Y_bottom = (SIGxys_bottom_set[i])[1]
    S_top = (SIGxys_top_set[i])[2]
    S_bottom = (SIGxys_bottom_set[i])[2]
    # 7.1 Fiber Mode
    if X_top > 0:
        Rf_top = math.sqrt(1/((X_top/Xt)**2 + (S_top/Sc)**2))
    else:
        Rf_top = -1*Xc/X_top
    if X_bottom > 0:
        Rf_bottom = math.sqrt(1/((X_bottom/Xt)**2 + (S_bottom/Sc)**2))
    else:
        Rf_bottom = -1*Xc/X_bottom
    # 7.2 Matrix Mode
    if Y_top > 0:
        Rm_top = math.sqrt(1/((Y_top/Yt)**2 + (S_top/Sc)**2))
    else:
        a_top = (Y_top/(2*Sc))**2 + (S_top/Sc)**2
        b_top = ((Yc/(2*Sc))**2 - 1)*(Y_top/Yc)
        c_top = -1
        Rm_top = (-b_top + math.sqrt(b_top**2 - 4*a_top*c_top))/(2*a_top)
    if Y_bottom > 0:
        Rm_bottom = math.sqrt(1/((Y_bottom/Yt)**2 + (S_bottom/Sc)**2))
    else:
        a_bottom = (Y_bottom/(2*Sc))**2 + (S_bottom/Sc)**2
        b_bottom = ((Yc/(2*Sc))**2 - 1)*(Y_bottom/Yc)
        c_bottom = -1
        Rm_bottom = (-b_bottom + math.sqrt(b_bottom**2 - 4*a_bottom*c_bottom))/(2*a_bottom)
    R_top_hashin_set[i] = [Rf_top,Rm_top]
    R_bottom_hashin_set[i] = [Rf_bottom,Rm_bottom]
    R_min_hashin_set[i] = min(Rf_top,Rm_top,Rf_bottom,Rm_bottom)
R_min_hashin = R_min_hashin_set[Fail_ply]
hashin_N = R_min_hashin*N126
hashin_M = R_min_hashin*M126
#**********************************************************************************************************************************************************************************
# 8.1 Print Out Material Data
print("\n************************Results************************")
print("Your material is: "+ str(mat_name))
print("\n****Material Properties****\nModulus Parameters:\nEx =",float(Ex),"GPa, Ey =",float(Ey),"GPa, Es =",float(Es),"GPa, Nu/x =",float(NUx))
print("\nStrength Parameters:\nX =",float(Xt),"MPa, X' =",float(Xc),"MPa, Y =",float(Yt),"MPa, Y' =",float(Yc),"MPa, S =",float(Sc),"MPa")

# 8.2 Print Out Laminate Data
print("\n*******Laminate Data******\nThe laminate code is: \n",laminate_short,int(repeat_int),str(sym_code))
print("\nComplete Angles: ",laminatecode)
print("With individual ply thickness",plythickness_matrix," mm")
print("\nCore thickness: 2Zc = ",TwoZc," mm")
print("The total laminate thickness is " + str(total_thickness) + " mm.")

# 8.3 Print out Non-normailized Matrices [A] and [a]
print("\n********Matrices*********\nNOTE: All resulting strain/stress arrays are in order 1,2,6 for Off-Axis & x,y,s for On-Axis\nYou Input "+str(Axis)," "+str(Force)," with values:\n"+str(Force1)," "+str(unit),"\n"+str(Force2)," "+str(unit),"\n"+str(Force3)," "+str(unit),"\nBending Moments: M1 = ",M1," Nm; M2 = ",M2," Nm; M6 = ",M6," Nm")
print("\nThe volume fraction of the ply orientations are:\nv0:",frac(numply0,ntotal)," = ",v0,"\nv90:",frac(numply90,ntotal)," = ",v90,"\nv45:",frac(numply45,ntotal)," = ",v45,"\nv-45:",frac(numplyn45,ntotal)," = ",vn45,"\nv30:",frac(numply30,ntotal)," = ",v30,"\nv-30:",frac(numplyn30,ntotal)," = ",vn30)
print("\nThe Invariant U values are\nU1=",U1Q," GPa, U2=",U2Q," GPa, U3=",U3Q," GPa, U4=",U4Q," GPa, U5=",U5Q," GPa")
print("\nThe Normalized Geometric Parameters V* are:\nV1*=",V1star,",V2*=",V2star,", V3*=",V3star,", V4*=",V4star)
print("\nA) Non-normalized In-plane Modulus [A]:\n",A_matrix," MPa*m")
print("\nB) Non-normalized In-plane Compliance [a]:\n",a_matrix,"1/GPa*m")
print("\nC) Non-normalized Flexural Modulus [D]:\n",D_matrix," MPa*m^2")
print("\nD) Non-normalized Flexural Compliance [d]:\n",d_matrix,"1/GPa*m^2")
print("\nE) In-plane Off-Axis Strain [EPS°126]:\n",EPSnot126,"mm/mm\n\nCurvature [k126]",np.dot(1/1000,k126)," m^-1")
print("\nF) In-Plane Off-Axis Strain [EPS126] due to curvature 'k126' for each layer, beginning at ply 0(initial layer):\nTop of ply:\n")
for key in EPS126_top_set.keys():
    print("Ply : {} , EPS126_top : {}".format(key,EPS126_top_set[key])," mm/mm")
print("\nBottom of ply:\n")
for key in EPS126_bottom_set.keys():
    print("Ply : {} , EPS126_bottom : {}".format(key,EPS126_bottom_set[key])," mm/mm")
print("\nG) On-Axis Strain [EPSxys] for each layer, beginning at ply 0(initial layer):\nTop of ply:\n")
for key in EPSxys_top_set.keys():
    print("Ply : {} , EPSxys_top : {}".format(key,EPSxys_top_set[key])," mm/mm")
print("\nBottom of ply:\n")
for key in EPSxys_bottom_set.keys():
    print("Ply : {} , EPSxys_bottom : {}".format(key,EPSxys_bottom_set[key])," mm/mm")
print("\nH) On-Axis Stress [SIGxys] for each layer, beginning at ply 0(initial layer):\nTop of ply:\n")
for key in SIGxys_top_set.keys():
    print("Ply : {} , SIGxys_top : {}".format(key,SIGxys_top_set[key])," MPa")
print("\nBottom of ply:\n")
for key in SIGxys_bottom_set.keys():
    print("Ply : {} , SIGxys_bottom : {}".format(key,SIGxys_bottom_set[key])," MPa")
print("\nI) R values for Max Stress Criterion:\nTop of ply:\n")
for key in R_top_set.keys():
    print("Ply : {} , R_top : {}".format(key,R_top_set[key]))
print("\nBottom of ply:\n")
for key in R_bottom_set.keys():
    print("Ply : {} , R_bottom : {}".format(key,R_bottom_set[key]))
print("The first ply to fail Max Stress Criterion is: ",Fail_ply," with R= ",R_min_maxstress, " due to Fiber Compression\nWith N vector:",maxstress_N," N/m or Pa*m\nAnd M vector:",maxstress_M," Nm")
print("\nJ) R values for Tsai-Wu Criterion:\nTop of ply:\n")
for key in R_top_tsai_set.keys():
    print("Ply : {} , R_top_tsai-wu : {}".format(key,R_top_tsai_set[key]))
print("\nBottom of ply:\n")
for key in R_bottom_tsai_set.keys():
    print("Ply : {} , R_bottom_tsai-wu : {}".format(key,R_bottom_tsai_set[key]))
print("The first ply to fail Tsai-Wu Criterion is: ",Fail_ply," with R= ",R_min_tsai," \nWith N vector:",tsai_N," N/m or Pa*m\nAnd M vector:",tsai_M," Nm")
print("\nK) R values for 2D Hashin Criterion:\nTop of ply:\n")
for key in R_top_hashin_set.keys():
    print("Ply : {} , R_top_hashin : {}".format(key,R_top_hashin_set[key]))
print("\nBottom of ply:\n")
for key in R_bottom_hashin_set.keys():
    print("Ply : {} , R_bottom_hashin- : {}".format(key,R_bottom_hashin_set[key]))
print("The first ply to fail 2D Hashin Criterion is: ",Fail_ply," with R= ",R_min_hashin," due to Fiber Compressive Mode\nWith N vector:",hashin_N," N/m or Pa*m\nAnd M vector:",hashin_M," Nm")


