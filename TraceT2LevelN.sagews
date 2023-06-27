#note that in all these computations N is odd because N needs to be coprime to 2 and 4 for these formulas to be correct
def omega(n): #this is the number of prime divisors of n
    omegan=len(list(factor(n)))
    return(omegan)

def mu2(t,N): #this computes the value of mu for T2 and level N
    mutwo= solve_mod(x^2-t*x+2==(0),N) #Finds all solutions in Z/NZ
    mu2total=0
    for i in range (0,len(mutwo)):
        if gcd(mutwo[i][0],N)==1: #Finds which solutions are elements of (Z/NZ)*
            mu2total=mu2total+1
    return(mu2total)

def mu4(t,N):#this computes the values of mu relavant to TrT4 and level N
    mufour= solve_mod(x^2-t*x+4==(0),N) #Finds all solutions
    mu4total=0
    for i in range (0,len(mufour)):
        if gcd(mufour[i][0],N)==1: #Finds which solutions are elements of (Z/NZ)*
            mu4total=mu4total+1
    return(mu4total)
    
def A2T2(N,k): #here k is half of the weight, this is the A2 component given in the book "Traces of Hecke Operators" for T2
    if k>1:
        return(0)
    if k==1:
        Atwo=2
        if N%2==1:
            Atwo=3
    return(Atwo)

def A2T4(N,k): #here k is half of the weight, this is the A2 component given in the book "Traces of Hecke Operators" for T4
    if k>1:
        return(0)
    if k==1:
        AT4=4
        if N%2==1:
            AT4=7
    return(AT4)

def P2k(k,t,m): #here k is half of the weight
    Ptot=0
    for j in range (0,k):
        #this is the combinatorial formula given in Chiriac and Jorza's paper "The trace of T2 takes no repeated values"
        Ptot= (-1)^j*binomial(2*k-2-j,j)*m^j*t^(2*k-2-2*j)+Ptot 
    return(Ptot)

#for TrT2 and TrT4, these are from the version of the Eichler-Selberg trace formula found in the book "Traces of Hecke Operators"
def TrT2(k,N): #k is half of the weight and N is the level
    Term2=-1/4*P2k(k,2,2)*(mu2(2,N)+mu2(-2,N))-1/2*P2k(k,1,2)*(mu2(1,N)+mu2(-1,N))-1/2*P2k(k,0,2)*mu2(0,N)
    Term3=-2^(omega(N))
    Trace2=A2T2(N,k)+Term2+Term3
    return(Trace2)

def TrT4(k,N):
    Term2=(2*k-1)/12*Gamma0(N).index()*4^(k-1)
    Term3=-1/2*(P2k(k,3,4)*(mu4(-3,N)+mu4(3,N))+P2k(k,2,4)*(4/3*(mu4(-2,N)+mu4(2,N)))+P2k(k,1,4)*(mu4(-1,N)+mu4(1,4))+P2k(k,0,4)*3/2*mu4(0,N))
    Term4=-2^(omega(N))-2^(2*k-2)*Gamma0(N).ncusps() #The number of cusps in level N showing up here is interesting
    Trace4=A2T4(N,k)+Term2+Term3+Term4
    return(Trace4)

def EigenSum(k,N): #This is sum(a(i)a(j))
    A=TrT2(k,N)^2
    B=TrT4(k,N)
    dk=ModularForms(N,2*k).cuspidal_subspace().dimension()
    Eigens=1/2*(A-B)-2^(2*k-2)*dk
    return(Eigens)








