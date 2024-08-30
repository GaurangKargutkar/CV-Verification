import numpy as np      #Gaurang Kargutkar #234103417
import math as mt


def mat_mult(A,B):
    r1=A.shape[0]
    c1=A.shape[1]
    r2=B.shape[0]
    c2=B.shape[1]
    sol=np.zeros([r1,c2])
    for i in range(r1):
            for k in range(c2):
                sum=0
                for j in range(r2):
                    sum+=A[i][j]*B[j][k]
                sol[i][k]=(sum)  
    return sol


def ABBD_formation(Q_bar):
    A=np.zeros([3,3])
    B=np.zeros([3,3])
    D=np.zeros([3,3])
    zk_minus_1=-tk*(n/2)   
    for i in range(n):
        A+=Q_bar[i]*tk
        zk=zk_minus_1+tk   
        B+=0.5*Q_bar[i]*(zk**2-zk_minus_1**2)
        D+=(1/3)*Q_bar[i]*(zk**3-zk_minus_1**3)
        zk_minus_1+=tk
    ABBD=np.zeros([6,6])
    for i in range(3):
        for j in range(3):
            ABBD[i][j]=A[i][j]
        for j in range(3,6):
            ABBD[i][j]=B[i][j-3] 
    for i in range(3,6):
        for j in range(3):
            ABBD[i][j]=B[i-3][j]
        for j in range(3,6):
            ABBD[i][j]=D[i-3][j-3] 
    return ABBD

def N_M_calculation (alpha_xy,beta_xy,N_M):     # Function to calculate Total N_M vector
    N_T=np.zeros([3,1])
    N_H=np.zeros([3,1])
    M_T=np.zeros([3,1])
    M_H=np.zeros([3,1])
    zk_minus_1=-tk*(n/2)
    for i in range(n):
        N_T+=mat_mult(Q_bar[i],alpha_xy[i])*tk*delta_T
        N_H+=mat_mult(Q_bar[i],beta_xy[i])*tk*delta_c
        zk=zk_minus_1+tk   
        M_T+=0.5*(zk**2-zk_minus_1**2)*delta_T*mat_mult(Q_bar[i],alpha_xy[i])
        M_H+=0.5*(zk**2-zk_minus_1**2)*delta_c*mat_mult(Q_bar[i],beta_xy[i])
        zk_minus_1+=tk
    N_M_total=np.zeros([6,1])
    for i in range(3):
        N_M_total[i][0]=N_T[i][0]+N_H[i][0]+N_M[i][0]
        N_M_total[i+3][0]=M_T[i][0]+M_H[i][0]+N_M[i+3][0]    
    return (N_M_total)

def Strain_calculation(e_k):            #Function to find Strain at each lamina
    epsilon_xy_0=np.zeros([3,1])
    k=np.zeros([3,1])
    for i in range(3):
        epsilon_xy_0[i][0]=e_k[0][i]
        k[i][0]=e_k[0][i+3]
    zk_minus_1=-tk*(n/2)
    epsilon_xy=np.zeros([n,3,1])
    for i in range(n):
        if zk_minus_1!=0:
            epsilon_xy[i]=epsilon_xy_0+zk_minus_1*k
        else:
            zk_minus_1+=tk
            epsilon_xy[i]=epsilon_xy_0+zk_minus_1*k
        zk_minus_1+=tk
    return epsilon_xy


#Code Starts from Here
N_M=np.array([[100000],[0],[0],[0],[0],[0]])   # {Nx,Ny,Nxy,Mx,My,Mxy}
su=np.array([[1062*10**6],[610*10**6],[31*10**6],[118*10**6],[72*10**6]])   #{Ultimate Tensile , Compressive and Shear Stress} 
delta_T=50  # Delta T
delta_c=0
tk=0.125*10**(-3)       # Lamina Thickness
E1=38.6*10**9
E2=8.27*10**9
v12=0.28
G12=4.14*10**9
v21=E2*v12/E1
Q11 =E1/(1-v12*v21)
Q22 = E2/(1-v12*v21)
Q66=G12
Q12 =Q21 = (v12*E2)/(1-v12*v21)
Q=np.array([[Q11,Q12,0],[Q21,Q22,0],[0,0,Q66]])
#print(Q)
Theta=np.array(np.radians([0,45,-45,90,90,-45,45,0]))
n=len(Theta) #No. of Layers of Lamina
symm=1
for i in range(n//2) :
    if Theta[i] != Theta[n-1-i]:
        symm = 0
        break
    
Q_bar=np.zeros([n,3,3])
T=np.zeros([n,3,3])
#Calculating Q bar for eac lamina
for i in range(n):
    Q_bar[i][0][0]=Q[0][0]*(np.cos(Theta[i]))**4+2*(Q[0][1]+2*Q[2][2])*((np.sin(Theta[i]))*(np.cos(Theta[i])))**2+Q[1][1]*np.sin(Theta[i])**4
    Q_bar[i][1][1]=Q[0][0]*(np.sin(Theta[i]))**4+2*(Q[0][1]+2*Q[2][2])*((np.sin(Theta[i]))*(np.cos(Theta[i])))**2+Q[1][1]*np.cos(Theta[i])**4
    Q_bar[i][0][1]=(Q[0][0]+Q[1][1]-4*Q[2][2])*((np.sin(Theta[i]))*(np.cos(Theta[i])))**2+Q[0][1]*(np.sin(Theta[i])**4+np.cos(Theta[i])**4)
    Q_bar[i][1][0]=(Q[0][0]+Q[1][1]-4*Q[2][2])*((np.sin(Theta[i]))*(np.cos(Theta[i])))**2+Q[0][1]*(np.sin(Theta[i])**4+np.cos(Theta[i])**4)
    Q_bar[i][2][2]=(Q[0][0]+Q[1][1]-2*Q[0][1]-2*Q[2][2])*((np.cos(Theta[i]))*(np.sin(Theta[i])))**2+Q[2][2]*(np.sin(Theta[i])**4+np.cos(Theta[i])**4)
    Q_bar[i][0][2]=(Q[0][0]-Q[0][1]-2*Q[2][2])*(np.sin(Theta[i]))*(np.cos(Theta[i]))**3-(Q[1][1]-Q[0][1]-2*Q[2][2])*(np.cos(Theta[i]))*(np.sin(Theta[i]))**3
    Q_bar[i][2][0]=(Q[0][0]-Q[0][1]-2*Q[2][2])*(np.sin(Theta[i]))*(np.cos(Theta[i]))**3-(Q[1][1]-Q[0][1]-2*Q[2][2])*(np.cos(Theta[i]))*(np.sin(Theta[i]))**3
    Q_bar[i][1][2]=(Q[0][0]-Q[0][1]-2*Q[2][2])*(np.cos(Theta[i]))*(np.sin(Theta[i]))**3-(Q[1][1]-Q[0][1]-2*Q[2][2])*(np.sin(Theta[i]))*(np.cos(Theta[i]))**3
    Q_bar[i][2][1]=(Q[0][0]-Q[0][1]-2*Q[2][2])*(np.cos(Theta[i]))*(np.sin(Theta[i]))**3-(Q[1][1]-Q[0][1]-2*Q[2][2])*(np.sin(Theta[i]))*(np.cos(Theta[i]))**3
#print(Q_bar) 

# Rotation matrix T[] for each lamina
for i in range(n):
    T[i][0][0]=(np.cos(Theta[i]))**2
    T[i][0][1]=(np.sin(Theta[i]))**2
    T[i][0][2]=2*(np.sin(Theta[i]))*np.cos(Theta[i])
    T[i][1][0]=(np.sin(Theta[i]))**2
    T[i][1][1]=(np.cos(Theta[i]))**2
    T[i][1][2]=-2*(np.sin(Theta[i]))*np.cos(Theta[i])
    T[i][2][0]=-1*(np.sin(Theta[i]))*np.cos(Theta[i])
    T[i][2][1]=(np.sin(Theta[i]))*np.cos(Theta[i])
    T[i][2][2]=(np.cos(Theta[i]))**2-(np.sin(Theta[i]))**2
No_of_Present_Lamina=n

iteration=0
if symm==1: #if Laminate is symmetric
    print(" This is a Symmetric Laminate")
    while No_of_Present_Lamina>0:
        ABBD = ABBD_formation(Q_bar)
        #print(ABBD)
        ABBD_inv=np.linalg.inv(ABBD)
        #print(ABBD_inv)
        e_k=np.zeros([6,1])

        # Hygro-Thermal load
        alpha_12=np.array([[8.6*10**(-6)],[22.1*10**(-6)],[0]])   #insert the values of CTEs in material axes
        beta_12=np.array([[1],[2],[0]])   #random values of beta taken
        beta_xy=np.zeros([n,3,1])
        alpha_xy=np.zeros([n,3,1])
        for i in range(n):
            alpha_xy[i]=mat_mult(np.linalg.inv(T[i]),alpha_12)
            beta_xy[i]=mat_mult(np.linalg.inv(T[i]),beta_12)
        
        #Total N_M Matrix is calculated
        N_M_total = N_M_calculation(alpha_xy,beta_xy,N_M)

        #Strain values at mid surface layer
        e_k=ABBD_inv*N_M_total
        
        #Strain and stress Calculation at Each Lamina layer
        epsilon_xy=Strain_calculation(e_k)
        sigma_xy=np.zeros([n,3,1])
        for i in range(n):
            sigma_xy[i]=mat_mult(Q_bar[i],epsilon_xy[i]) 
        sigma_12=np.zeros([n,3,1])

        SR=np.zeros([n,4,1])
        Tsai_Hill=np.zeros(n)
        Failure_ply=[]
        Strength_ratio = 0
        Strength_ratio_new=0
        for i in range(n):
            sigma_12[i]=mat_mult(T[i],sigma_xy[i])
            if sigma_12[i][0][0]>=0:
                SR[i][0][0]=sigma_12[i][0][0]/su[0]
            else:
                SR[i][0][0]=sigma_12[i][0][0]/su[1]

            if sigma_12[i][1][0]>=0:
                SR[i][1][0]= sigma_12[i][1][0]/su[0]
            else:
                SR[i][1][0]=sigma_12[i][1][0]/su[1]

            if sigma_12[i][1][0]>=0:
                SR[i][2][0]=sigma_12[i][1][0]/su[2]
            else:
                SR[i][2][0]=sigma_12[i][1][0]/su[3]

            SR[i][3][0]=sigma_12[i][2][0]/su[4]

            Tsai_Hill[i]=SR[i][0][0]**2-SR[i][0][0]*SR[i][1][0]+SR[i][2][0]**2+SR[i][3][0]**2
            if Tsai_Hill[i]>=1:
                Failure_ply.append(i)
                print("Lamina ",i+1,"Fails")
                No_of_Present_Lamina-=1
                Strength_ratio = max(Strength_ratio,max(SR[i]))
            Strength_ratio_new = max(Strength_ratio_new,max(SR[i]))
       
        if No_of_Present_Lamina == n-2 and Strength_ratio!=0:
            FPF = N_M/Strength_ratio
        if No_of_Present_Lamina== 0:
            LPF =N_M/Strength_ratio
            print("FPF=",FPF[0]," N/m")
            print("LPF",LPF[0]," N/m")
        if Strength_ratio==0:
            N_M=N_M/Strength_ratio_new
        else:
            N_M=N_M/Strength_ratio
        
        for k in Failure_ply: 
            for i in range(3):
                for j in range(3):
                    Q_bar[k][i][j]=0    #Complete Degradation 
        iteration+=1
        print("Iteration no ", iteration ," is completed",)
else:       #if Laminate is Unsymmetric
    print(" This is a UnSymmetric Laminate")
    while No_of_Present_Lamina>0:
        ABBD = ABBD_formation(Q_bar)
        #print(ABBD)
        ABBD_inv=np.linalg.inv(ABBD)
        #print(ABBD_inv)
        e_k=np.zeros([6,1])
        #Thermal load
        alpha_12=np.array([[8.6*10**(-6)],[22.1*10**(-6)],[0]])   #insert the values of CTEs in material axes
        alpha_xy=np.zeros([n,3,1])
        for i in range(n):
            alpha_xy[i]=mat_mult(np.linalg.inv(T[i]),alpha_12)
        #Hygroscopic load
        beta_12=np.array([[1],[2],[0]])   #random values of beta taken
        beta_xy=np.zeros([n,3,1])
        for i in range(n):
            beta_xy[i]=mat_mult(np.linalg.inv(T[i]),beta_12)
        #Total N_M Matrix is calculated
        N_M_total = N_M_calculation(alpha_xy,beta_xy,N_M)

        #Strain values at mid surface layer
        e_k=ABBD_inv*N_M_total
        
        #Strain Calculation at Each Lamina layer
        epsilon_xy=Strain_calculation(e_k)
        sigma_xy=np.zeros([n,3,1])
        for i in range(n):
            sigma_xy[i]=mat_mult(Q_bar[i],epsilon_xy[i]) 
        sigma_12=np.zeros([n,3,1])

        SR=np.zeros([n,4,1])
        Tsai_Hill=np.zeros(n)
        Failure_ply=[]
        Strength_ratio = 0
        Strength_ratio_new=0
        for i in range(n):
            sigma_12[i]=mat_mult(T[i],sigma_xy[i])
            if sigma_12[i][0][0]>=0:
                SR[i][0][0]=sigma_12[i][0][0]/su[0]
            else:
                SR[i][0][0]=sigma_12[i][0][0]/su[1]

            if sigma_12[i][1][0]>=0:
                SR[i][1][0]=sigma_12[i][1][0]/su[0]
            else:
                SR[i][1][0]=sigma_12[i][1][0]/su[1]

            if sigma_12[i][1][0]>=0:
                SR[i][2][0]=sigma_12[i][1][0]/su[2]
            else:
                SR[i][2][0]=sigma_12[i][1][0]/su[3]

            SR[i][3][0]=sigma_12[i][2][0]/su[4]

            Tsai_Hill[i]=SR[i][0][0]**2-SR[i][0][0]*SR[i][1][0]+SR[i][2][0]**2+SR[i][3][0]**2
            if Tsai_Hill[i]>=1:
                Failure_ply.append(i)
                print("Lamina ",i+1,"Fails")
                No_of_Present_Lamina-=1
                Strength_ratio = max(Strength_ratio,max(SR[i]))
            else:
                print("Lamina ",i+1," is Safe ")
            Strength_ratio_new = max(Strength_ratio_new,max(SR[i]))
       
        if No_of_Present_Lamina == n-1 and Strength_ratio!=0:
            FPF = N_M/Strength_ratio
        if No_of_Present_Lamina== 0:
            LPF =N_M/Strength_ratio
            print("FPF =",FPF[0]," N/m")
            print("LPF = ",LPF[0]," N/m")
        if Strength_ratio==0:
            N_M=N_M/Strength_ratio_new
        else:
            N_M=N_M/Strength_ratio
        
        for k in Failure_ply: 
            for i in range(3):
                for j in range(3):
                    Q_bar[k][i][j]=0    #Complete Degradation 
        iteration+=1
        print("Iteration no ", iteration ," is completed",)

        



    
    


 