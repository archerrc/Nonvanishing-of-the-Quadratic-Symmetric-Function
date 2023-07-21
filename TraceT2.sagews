#Note that this is only for level 1. Code that is valid for any odd level is given in the file named TraceT2LevelN
def P2k(k,t,m): #here k is half of the weight
    Ptot=0
    for j in range (0,k):
        #this is the combinatorial formula given in Chiriac and Jorza's paper "The trace of T2 takes no repeated values"
        Ptot= (-1)^j*binomial(2*k-2-j,j)*m^j*t^(2*k-2-2*j)+Ptot 
    return(Ptot)

#for TrT2 and TrT4, these are from the version of the Eichler-Selberg trace formula found in Zagier's appendix to Lang's "Introduction to modular forms"
def TrT2(k): #k is half of the weight
    Trace2=-1/2*P2k(k,0,2)-P2k(k,1,2)-1/2*P2k(k,2,2)-1
    return(Trace2)

def TrT4(k):
    Trace4=-3/4*P2k(k,0,4)-2*P2k(k,1,4)-4/3*P2k(k,2,4)-P2k(k,3,4)+1/12*P2k(k,4,4)-2^(2*k-2)-1
    return(Trace4)

def EigenSum(k): #This is sum(a(i)a(j))
    A=TrT2(k)^2
    B=TrT4(k)
    dk=floor(k/6)
    if k%6==1 & k!=1:
        dk=dk-1
    Eigens=1/2*(A-B)-2^(2*k-2)*dk
    return(Eigens)

for i in range (1,68):
    if EigenSum(i)==0:
            print('zero at weight', 2*i)
    if i==67:
            print('Done')










